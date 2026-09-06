"""Extended measured-event storyboard, independent of video frame signatures."""

from pathlib import Path
import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from PIL import Image

from .panels import MainFigure
from .render import frame_metrics,write_json


def storyboard_events(times,meta):
    items=[dict(r) for r in meta["key_frames"]]
    for event in meta["modal"]["extrema"]:
        if event["time"]<2:continue
        i=int(np.argmin(abs(times-event["time"])))
        if any(r["index"]==i for r in items):continue
        items.append({"name":f"late_{i}","phase":"Damped modal "+event["kind"],"event_time_s":event["time"],
                      "index":i,"selected_time_s":float(times[i])})
    return sorted(items,key=lambda r:r["index"])


def storyboard(cache,meta,output):
    output=Path(output)
    with h5py.File(cache,"r") as h,PdfPages(output/"key_phases.pdf") as pdf:
        rows=[];fig=MainFigure(h,meta)
        for item in storyboard_events(h["times"][:],meta):
            i=item["index"];fig.update(i)
            values=[float(h[key][i]) for key in ("kinetic_energy","potential_energy","total_energy")]
            fig.phase.set_text(item["phase"]+f" | K={values[0]:.4e}, P={values[1]:.4e}, E={values[2]:.4e} m^4/s²")
            pdf.savefig(fig.fig,dpi=100)
            rows.append(frame_metrics(h,meta,item))
        plt.close(fig.fig)
    write_json(output/"storyboard_summary.json",{"pages":len(rows),"frames":rows,"source":meta["sources"]["fine"]})
    # Six required representative previews; the PDF additionally covers later extrema.
    sheet=Image.new("RGB",(1920,1620),"white")
    for j,item in enumerate(meta["key_frames"]):
        with Image.open(output/"frames_preview"/(item["name"]+".png")) as im:
            resample=Image.Resampling.LANCZOS if hasattr(Image,"Resampling") else Image.LANCZOS
            sheet.paste(im.resize((960,540),resample),((j%2)*960,(j//2)*540))
    sheet.save(output/"key_phases.png")
