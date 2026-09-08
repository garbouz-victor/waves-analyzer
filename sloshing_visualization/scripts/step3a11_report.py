"""Read-only analysis of STEP3A11 temporary I/O measurements, never PDE."""
import csv
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

import numpy as np
from step3a11_io_attribution import (ROOT, PHASES, read, output, quantiles, wilson,
    operation_classification, context_classification, policy_guard, check_bindings)


def summary(batch):
    records = batch["records"]; stalls = [r for r in records if r["stall"]]
    keys = [k+"_s" for k in PHASES]+["total_append_s", "total_io_s", "encode_hash_s", "non_fsync_s"]
    return {"n":len(records), "stall_count":len(stalls), "stall_fraction":len(stalls)/len(records) if records else None,
        "wilson95":wilson(len(stalls),len(records)), "statistics":{k:quantiles([r[k] for r in records]) for k in keys},
        "stall_records":stalls, "stall_spacing":np.diff([r["repetition"] for r in stalls]).tolist(),
        "modulo_exploratory":{str(m):{str(k):sum(r["repetition"]%m == k for r in stalls) for k in range(m)} for m in (2,4,5,8,10,16)}}


def correlation(x, y):
    pairs = [(a,b) for a,b in zip(x,y) if a is not None and b is not None]
    if len(pairs)<3: return None
    a,b = np.asarray(pairs,dtype=float).T
    if np.std(a)==0 or np.std(b)==0: return None
    return float(np.corrcoef(a,b)[0,1])


def kernel_analysis(records):
    result = {}
    for key in ("Dirty", "Writeback", "WritebackTmp", "Buffers", "Cached"):
        groups = {name:[r["kernel_before_fsync"]["meminfo"][key] for r in records
            if bool(r["stall"]) == flag and r.get("kernel_before_fsync") and r["kernel_before_fsync"]["meminfo"][key] is not None]
            for name,flag in (("stall",True),("nonstall",False))}
        result[key] = {k:{"n":len(v), **quantiles(v)} for k,v in groups.items()}
        result[key]["units"] = "kB from proc/meminfo; descriptive, not causal"
    return result


def table(headers, rows):
    def fmt(x):
        if x is None: return "N/A"
        if isinstance(x, float): return f"{x:.12g}"
        return str(x).replace("|", "/").replace("\n", " ")
    return "\n".join(["| "+" | ".join(headers)+" |", "| "+" | ".join("---" for _ in headers)+" |"]+
                     ["| "+" | ".join(fmt(x) for x in row)+" |" for row in rows])+"\n"


def figures(hist, batches):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out = ROOT/"analysis"; out.mkdir(parents=True,exist_ok=True)
    def finish(fig,name):
        fig.tight_layout(); fig.savefig(out/name); plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4)); ax.plot([r["step"] for r in hist["records"]],[r["archival_s"] for r in hist["records"]],".-")
    ax.axhline(1.,ls="--",color="red",label="1 s descriptive stall");ax.set(yscale="log",xlabel="accepted step",ylabel="ordinary archival wall [s]");ax.legend()
    finish(fig,"historical_ordinary_archival_latency.pdf")
    bind=batches["bind"]["records"]; x=[r["repetition"] for r in bind]
    fig,ax=plt.subplots(figsize=(10,4));bottom=np.zeros(len(bind))
    for key,label in (("encode_hash_s","prepare/hash/encode"),("open_s","open"),("write_s","write"),("flush_s","flush"),("fsync_s","fsync"),("close_s","close")):
        vals=np.array([r[key] for r in bind]);ax.bar(x,vals,bottom=bottom,label=label,width=.95);bottom+=vals
    ax.axhline(1.,ls="--",color="red");ax.set(xlabel="bind burst repetition",ylabel="phase wall [s]");ax.legend(ncol=3)
    finish(fig,"bind_append_phase_breakdown.pdf")
    fig,ax=plt.subplots(figsize=(8,4))
    for i,name in enumerate(("bind","host","container_tmp","docker_volume")):
        vals=[r["fsync_s"] for r in batches[name]["records"]]
        ax.scatter(np.full(len(vals),i),vals,s=12,alpha=.4)
        ax.plot([i-.2,i+.2],[quantiles(vals)["p95"]]*2,color="red")
    ax.set_xticks(range(4),["container bind","host","container /tmp","Docker volume"]);ax.set(yscale="log",ylabel="fsync wall [s]")
    finish(fig,"fsync_context_comparison.pdf")
    fig,ax=plt.subplots(figsize=(8,4))
    for name in ("bind","paced"):
        rec=batches[name]["records"];ax.plot([r["repetition"] for r in rec],[r["fsync_s"] for r in rec],".-",label=name)
    ax.set(yscale="log",xlabel="repetition (sleep excluded)",ylabel="fsync wall [s]");ax.legend()
    finish(fig,"burst_vs_paced_fsync.pdf")
    valid=[r for r in bind if r.get("kernel_before_fsync") and r["kernel_before_fsync"]["meminfo"]["Dirty"] is not None]
    if valid:
        fig,ax=plt.subplots(figsize=(7,4));ax.scatter([r["kernel_before_fsync"]["meminfo"]["Dirty"] for r in valid],[r["fsync_s"] for r in valid],c=[int(r["stall"]) for r in valid])
        ax.set(yscale="log",xlabel="Dirty before fsync [kB]",ylabel="fsync wall [s]");finish(fig,"fsync_vs_dirty_pages.pdf")
    fig,axes=plt.subplots(1,2,figsize=(10,4));last_size=bind[0]["file_size_before"] if bind else 0;since=[]
    for r in bind:
        since.append(r["file_size_before"]-last_size)
        if r["stall"]:last_size=r["file_size_after"]
    for ax,xv,label in ((axes[0],[r["file_size_before"] for r in bind],"journal bytes before append"),
                        (axes[1],since,"bytes since previous stall / initial size")):
        ax.scatter(xv,[r["fsync_s"] for r in bind],c=[int(r["stall"]) for r in bind]);ax.set(yscale="log",xlabel=label,ylabel="fsync wall [s]")
    finish(fig,"fsync_vs_file_size.pdf")


def generate():
    policy=policy_guard();check_bindings(read(ROOT/"preflight/prefix_integrity.json")["files"])
    hist=read(ROOT/"historical/ordinary_archive.json");suite=read(ROOT/"benchmarks/suite_status.json")
    names={"bind":"bind_burst", "host":"host_burst", "container_tmp":"container_tmp_burst", "docker_volume":"docker_volume_burst", "paced":"bind_paced"}
    batches={k:read(ROOT/"benchmarks"/(v+".json")) for k,v in names.items() if (ROOT/"benchmarks"/(v+".json")).exists()}
    if len(batches)!=5 or suite["status"]!="complete":
        output(ROOT/"summary.json", {"status":"incomplete diagnostics", "suite":suite, "operation":"JOURNAL FSYNC NOT YET ISOLATED", "storage_path":"STORAGE PATH UNRESOLVED", "overall":"MODEL NOT YET VALIDATED"})
        raise RuntimeError("Incomplete suite: no positive attribution; inspect retained batch records")
    summaries={k:summary(v) for k,v in batches.items()}
    operation=operation_classification(batches["bind"]["records"])
    operation["primary_sample"]="bind burst, no pooling of paced or controls"
    operation["complete_N"] = len(batches["bind"]["records"])==120
    metrics={k:{"p95":s["statistics"]["fsync_s"]["p95"],"stall_fraction":s["stall_fraction"]} for k,s in summaries.items() if k!="paced"}
    op_pass=operation["classification"]=="FSYNC_DOMINATED" and operation["complete_N"] and suite["actual_cap_pass"]
    storage=context_classification(metrics,op_pass)
    overlap=max(hist["wilson95"][0],summaries["bind"]["wilson95"][0]) <= min(hist["wilson95"][1],summaries["bind"]["wilson95"][1])
    operation.update(historical_wilson95=hist["wilson95"],benchmark_wilson95=summaries["bind"]["wilson95"],wilson_intervals_overlap=overlap)
    kernel={k:kernel_analysis(v["records"]) for k,v in batches.items()}
    cadence={"historical_correlations":{k:correlation([r[k] for r in hist["records"]],[r["archival_s"] for r in hist["records"]])
        for k in ("step","offset_before","line_bytes","block","dt","core_s","previous_step_core_s")},
        "bind_file_size_latency_correlation":correlation([r["file_size_before"] for r in batches["bind"]["records"]],[r["fsync_s"] for r in batches["bind"]["records"]]),
        "bind_stall_spacing":summaries["bind"]["stall_spacing"], "historical_stall_spacing":hist["stall_step_spacing"],
        "byte_boundary_remainders":{str(m):[{"step":r["step"],"offset_mod":r["offset_before"]%m,"line_crosses_boundary":r["offset_before"]//m != r["file_size_after"]//m} for r in hist["records"] if r["stall"]] for m in (4096,65536,1048576)},
        "interpretation":"exploratory only; no causal or periodicity significance claim"}
    controls={m:read(ROOT/"benchmarks"/("control_"+m+".json")) for m in ("fsync","fdatasync","flush_only")}
    atomic={k:read(ROOT/"benchmarks"/("atomic_"+k+".json")) for k in ("bind","container_tmp")}
    npz=read(ROOT/"benchmarks/archive_npz.json")
    trace=read(ROOT/"benchmarks/strace_sequence.json"); trace_stats={}
    if (ROOT/"benchmarks/strace.log").exists():
        for call in ("fsync","fdatasync","write","openat","close","rename"):
            vals=[]
            for line in (ROOT/"benchmarks/strace.log").read_text().splitlines():
                if call+"(" in line:
                    match=re.search(r"<([0-9.]+)>$",line)
                    if match:vals.append(float(match.group(1)))
            trace_stats[call]={"n":len(vals),**quantiles(vals)}
    operation_verdict="JOURNAL FSYNC LATENCY IDENTIFIED" if op_pass else "JOURNAL FSYNC NOT YET ISOLATED"
    path_verdict={"BIND_MOUNT_SPECIFIC_EVIDENCE":"BIND-MOUNT-SPECIFIC EVIDENCE", "HOST_FILESYSTEM_DURABILITY_EVIDENCE":"HOST-FILESYSTEM DURABILITY EVIDENCE", "DOCKER_STORAGE_PATH_EVIDENCE":"DOCKER-STORAGE-PATH EVIDENCE", "STORAGE_PATH_UNRESOLVED":"STORAGE PATH UNRESOLVED"}[storage["classification"]]
    final=read(ROOT/"final_integrity.json")
    tests={}
    for path in (ROOT/"tests").glob("*.xml"):
        tests[path.name]=ET.parse(path).getroot().find("testsuite").attrib
    result={"operation_verdict":operation_verdict,"storage_path_verdict":path_verdict,"suite":suite,"contexts":summaries,
        "operation":operation,"storage":storage,"historical":{k:v for k,v in hist.items() if k!="records"},
        "kernel":kernel,"cadence":cadence,"strace":trace_stats,"final_integrity":final,"tests":tests,
        "remote_CI":"not observed; no push", "full_CHNS_transfer":"incomplete; STEP3A10 authorization remains false",
        "interval128_executed":False,"accepted_steps":127,"overall":"MODEL NOT YET VALIDATED"}
    for name,value in (("classification/operation.json",operation),("classification/storage_path.json",storage),
                       ("classification/summary.json",{"operation":operation_verdict,"storage_path":path_verdict}),
                       ("analysis/kernel_state.json",kernel),("analysis/cadence.json",cadence),
                       ("analysis/strace.json",trace_stats),("summary.json",result)):
        output(ROOT/name,value)
    with (ROOT/"summary.csv").open("x",newline="") as stream:
        writer=csv.writer(stream);writer.writerow(["context","n","stalls","stall_fraction","fsync_p50","fsync_p95","fsync_max"])
        for k,s in summaries.items():writer.writerow([k,s["n"],s["stall_count"],s["stall_fraction"],*[s["statistics"]["fsync_s"][q] for q in ("p50","p95","max")]])
    figures(hist,batches)
    bind=summaries["bind"];stalls=bind["stall_records"]
    phase_rows=[]
    for label,keys in (("prepare/hash",["copy_prepare_s","record_hash_s"]),("encode",["json_encode_s"]),*[ (k,[k+"_s"]) for k in ("open","write","flush","fsync","close")]):
        q=quantiles([sum(r[k] for k in keys) for r in batches["bind"]["records"]])
        share=float(np.median([sum(r[k] for k in keys)/r["total_append_s"] for r in stalls])) if stalls else None
        phase_rows.append([label,*[q[k] for k in ("p50","p90","p95","max")],share])
    context_rows=[["historical ordinary archival","historical bind",hist["n"],hist["stall_fraction"],"not phase-timed","not phase-timed","not phase-timed"]]
    for k in ("bind","host","container_tmp","docker_volume"):
        s=summaries[k];context_rows.append([k,batches[k]["environment"]["mount"]["filesystem"],s["n"],s["stall_fraction"],*[s["statistics"]["fsync_s"][q] for q in ("p50","p95","max")]])
    sections=[]
    def section(i,title,text): sections.append(f"## {i}. {title}\n\n{text}\n")
    section(1,"Goal","Attribute ordinary journal latency without a PDE solve, cost reauthorization or storage change.")
    section(2,"Inherited frozen scientific state","Prefix 127/1068 at 6.148107277665776e-10 s. Independent guard and all 72 frozen archive members PASS. The historical execution HEAD is preserved; this administrative execution base is "+policy["base_HEAD"]+".")
    section(3,"Why STEP3A10 cost audit failed","Frozen forecast 16354.593801770276 s exceeded 14400 s; authorization=false remains unchanged. Ordinary archival tails were not confined to sparse events.")
    section(4,"Exact ordinary archival code path","RateTrajectory.advance times scalar Journal.append only on F0_C0_A0. Copy/previous hash → record hash → open(ab) → JSON encode → write → flush → fsync → close → memory publication. Timing journal append is outside both archival_s and total_s. Field/checkpoint NPZ and expanded-audit files are absent inside ordinary archival_s.")
    section(5,"Journal durability semantics","fsync is intentional in the restart model. Atomic_json fsyncs the temporary file then os.replace; archive fsyncs NPZ, publishes metadata, then renames the directory. Neither fsyncs the parent directory. Atomic rename is not by itself a guarantee of persistence after host/VM crash. No durability semantics were changed.")
    section(6,"Attribution policy and predeclared rules","See policy.json and STEP3A11_DESIGN.md. Primary operation decision uses bind burst ONLY: >=5 stalls, median fsync share>=0.90, non-fsync p95<=0.10 s, no other phase median share>=0.20. All repetitions/factors frozen before measurements; no pooling, outlier deletion or retuning. Probe overhead is included conservatively in non-fsync and total append, not the syscall timer.")
    section(7,"Environment and mount topology",table(["context","filesystem","N","stall fraction","fsync p50 [s]","p95 [s]","max [s]"],context_rows)+"\nDocker context: "+policy["docker_context"]+"; image "+policy["image_digest"]+". Exact mount sources/options/device IDs, Python versions and complete mountinfo are archived. Filesystem claims use observations, not wrapper comments.")
    hist_rows=[[r[k] for k in ("step","time","archival_s","total_s","core_s","offset_before","line_bytes")] for r in hist["records"] if r["stall"]]
    section(8,"Historical ordinary stalls",f"N={hist['n']}, stalls={hist['stall_count']}, fraction={hist['stall_fraction']}; statistics: `{json.dumps(hist['statistics'])}`. All stall events follow.\n\n"+table(["step","physical time [s]","archival [s]","total [s]","core [s]","byte offset","line bytes"],hist_rows))
    section(9,"Instrumented append equivalence","Ten representative real payloads independently matched frozen Journal.append byte-for-byte, with valid Journal.read hash chains, in each of the four contexts. Initial journal bytes and payload family are identical. Setup copying/loading is outside append timing. Secondary atomic JSON also matches exact bytes; NPZ arrays/metadata are validated after each archive.")
    section(10,"Bind burst results",f"N={bind['n']}, stalls={bind['stall_count']}; total append `{json.dumps(bind['statistics']['total_append_s'])}`; non-fsync `{json.dumps(bind['statistics']['non_fsync_s'])}`.")
    section(11,"Bind paced results",f"60 appends after the first 60 historical ordinary core durations; no delay cap. Sleep total {sum(r['core_s'] for r in policy['paced_delays'])} s is outside append statistics. Stall fraction {summaries['paced']['stall_fraction']}, fsync statistics `{json.dumps(summaries['paced']['statistics']['fsync_s'])}`.")
    for i,k,title in ((12,"host","Host results"),(13,"container_tmp","Container tmp results"),(14,"docker_volume","Docker volume results")):
        section(i,title,f"{summaries[k]['n']} exact burst appends; fsync `{json.dumps(summaries[k]['statistics']['fsync_s'])}`; stalls {summaries[k]['stall_count']}. Only temporary benchmark storage was used.")
    events=[]
    for r in stalls:
        mem=r["kernel_before_fsync"]["meminfo"] if r.get("kernel_before_fsync") else {}
        events.append([r["repetition"],r["total_append_s"],r["fsync_s"],r["phase_fractions"]["fsync"],r["file_size_before"],mem.get("Dirty"),mem.get("Writeback"),r["total_append_cpu_s"],r["total_append_s"]])
    section(15,"Phase attribution",table(["phase","p50 [s]","p90 [s]","p95 [s]","max [s]","median stall share"],phase_rows)+"\nAll primary bind stall events:\n\n"+table(["rep","total [s]","fsync [s]","fsync fraction","file size before","Dirty kB","Writeback kB","CPU [s]","wall [s]"],events)+f"\nFrozen checks: `{json.dumps(operation['checks'])}`. Wilson historical {hist['wilson95']}, bind {bind['wilson95']}; overlap={overlap}. Frequency comparison is descriptive only.")
    control_rows=[]
    for mode,batch in controls.items():
        values=[r["total_append_s"] for r in batch["records"]];q=quantiles(values)
        control_rows.append([mode,len(values),q["p50"],q["p95"],q["max"],sum(v>=1 for v in values)])
    section(16,"Sync-mode controls",table(["mode","N","total p50 [s]","p95 [s]","max [s]","stalls"],control_rows)+"\nTemporary controls only; no recommendation to replace production fsync. Optional strace context/status: "+str(trace.get("context",trace.get("status")))+"; syscall statistics `"+json.dumps(trace_stats)+"`.")
    section(17,"Atomic-json secondary results",table(["context","N","temp-file fsync p50","p95","max"],[[k,len(v["records"]),*[quantiles([r["fsync_s"] for r in v["records"]])[q] for q in ("p50","p95","max")]] for k,v in atomic.items()])+"\nSeparate distribution; not an ordinary scalar append operation.")
    section(18,"NPZ/archive secondary results",f"N={len(npz['records'])}; NPZ fsync `{json.dumps(quantiles([r['fsync_s'] for r in npz['records']]))}`; metadata fsync `{json.dumps(quantiles([r['metadata_atomic']['fsync_s'] for r in npz['records']]))}`; compression/write `{json.dumps(quantiles([r['npz_compress_write_s'] for r in npz['records']]))}`. These costs cannot cause historical F0_C0_A0 latency because that path does not archive NPZ.")
    section(19,"Dirty/writeback correlation","Full stall/nonstall medians, p90 and ranges are in analysis/kernel_state.json. Bind observations: `"+json.dumps({k:kernel['bind'][k] for k in ('Dirty','Writeback')})+"`. Disk device mapping is used only if directly observed; null means unavailable. Correlation is not causation.")
    section(20,"File-size/cadence analysis","Historical step spacings: `"+json.dumps(hist['stall_step_spacing'])+"`; bind repetition spacings: `"+json.dumps(bind['stall_spacing'])+"`. Correlations: `"+json.dumps(cadence['historical_correlations'])+"`. Full offsets, modulo exploratory counts and byte-boundary remainders are retained. No periodicity significance or deterministic size-trigger mechanism is inferred from this sample.")
    section(21,"Operation attribution verdict",operation_verdict+". Phase measurement identifies the present diagnostic operation only; matching historical tail frequency does not retrospectively instrument old syscalls.")
    section(22,"Storage-path attribution verdict",path_verdict+". Frozen criteria: `"+json.dumps(storage)+"`. This is evidence about observed paths, not proof of an OS or Docker defect.")
    section(23,"What is NOT proven","No new PDE result, no full-horizon coupling evidence, no cost authorization, no proof of crash durability under a changed scheme, and no kernel/VM bug established. Host/container controls also differ in execution environment; the path classification is empirical evidence, not elimination of all possible confounders.")
    next_work=("A separate durable Docker-volume journal/trajectory-storage experiment with reproducible export, byte/hash equivalence and forced-crash/restart validation; keep current fsync semantics initially." if storage['classification']=='BIND_MOUNT_SPECIFIC_EVIDENCE' else "A separate repeated native-host versus VM durable-append experiment with externally observed writeback/device timing; do not change production durability until localization is stronger.")
    section(24,"Candidate next storage experiments",next_work+" Not implemented. Removing fsync is not a harmless optimization.")
    section(25,"Prefix final integrity",f"Final guard `{json.dumps(final)}`. All inherited STEP3A9 prefix and STEP3A10 files remain byte-identical. Accepted=127, no row/field/checkpoint128, no new primary session, no FAILED/COMPLETE. interval128 NOT RUN.")
    section(26,"Tests / CI","Only pure I/O/policy tests; no numerical FEM solve in tests. Results: `"+json.dumps(tests)+f"`. Diagnostic wall {suite['wall_s']} s / 3600 s; all batches/cleanup included, pytest/report rendering excluded. Remote CI not observed; no push.")
    section(27,"Exact verdict",operation_verdict+"; "+path_verdict+". Full CHNS transfer remains incomplete; STEP3A10 authorization=false unchanged. MODEL NOT YET VALIDATED.")
    Path("STEP3A11_REPORT.md").write_text("# STEP3A11 — Journal / fsync latency attribution\n\n"+"\n".join(sections))
    print(json.dumps({"operation":operation_verdict,"storage_path":path_verdict,"wall_s":suite["wall_s"]},indent=2))


if __name__ == "__main__": generate()
