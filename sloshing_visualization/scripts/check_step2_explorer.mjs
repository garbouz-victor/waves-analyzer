// Optional real-browser smoke check. No npm dependencies; Node >=22 + Chrome.
// Temporary Chrome profile only; no user's browser profile is accessed.
import {spawn} from 'node:child_process';
import {mkdtemp,readFile,writeFile,mkdir} from 'node:fs/promises';
import {tmpdir} from 'node:os';
import {resolve,join} from 'node:path';
import {pathToFileURL} from 'node:url';
import assert from 'node:assert/strict';

const out=resolve(process.argv[2]||'output/step2');
const profile=await mkdtemp(join(tmpdir(),'sloshing-step2-browser-'));
const child=spawn(process.env.CHROME_BIN||'google-chrome',[
  '--headless=new','--disable-gpu','--no-first-run','--no-default-browser-check',
  '--disable-background-networking','--remote-debugging-port=0',`--user-data-dir=${profile}`,'about:blank'
],{stdio:['ignore','ignore','pipe']});
let ws;
try{
  const url=await new Promise((ok,fail)=>{
    let log='';const timeout=setTimeout(()=>fail(new Error('Chrome startup timeout: '+log)),20000);
    child.on('error',fail);child.on('exit',code=>fail(new Error('Chrome exited '+code+': '+log)));
    child.stderr.on('data',chunk=>{log+=chunk;const match=log.match(/DevTools listening on (ws:\/\/[^\s]+)/);
      if(match){clearTimeout(timeout);ok(match[1]);}});
  });
  ws=new WebSocket(url);await new Promise((ok,fail)=>{ws.onopen=ok;ws.onerror=fail;});
  let id=0;const pending=new Map(),exceptions=[];
  ws.onmessage=event=>{const m=JSON.parse(event.data);
    if(m.method==='Runtime.exceptionThrown')exceptions.push(m.params.exceptionDetails);
    if(pending.has(m.id)){const [ok,fail]=pending.get(m.id);pending.delete(m.id);m.error?fail(new Error(JSON.stringify(m.error))):ok(m.result);}};
  const send=(method,params={},sessionId)=>new Promise((ok,fail)=>{const n=++id;pending.set(n,[ok,fail]);ws.send(JSON.stringify({id:n,method,params,sessionId}));});
  const {targetId}=await send('Target.createTarget',{url:'about:blank'});
  const {sessionId}=await send('Target.attachToTarget',{targetId,flatten:true});
  const call=(m,p={})=>send(m,p,sessionId);
  await call('Runtime.enable');await call('Page.enable');
  await call('Emulation.setDeviceMetricsOverride',{width:1600,height:1150,deviceScaleFactor:1,mobile:false});
  const evaluate=async expression=>{const r=await call('Runtime.evaluate',{expression,returnByValue:true,awaitPromise:true});
    if(r.exceptionDetails)throw new Error(JSON.stringify(r.exceptionDetails));return r.result.value;};
  await call('Page.navigate',{url:pathToFileURL(join(out,'explorer.html')).href});
  let ready=false;
  for(let k=0;k<120;k++){ready=await evaluate('window.step2Ready === true');if(ready)break;await new Promise(r=>setTimeout(r,500));}
  assert(ready,'Explorer must finish first Plotly rendering');
  const state=()=>evaluate(`(()=>{const p=document.getElementById('plot');return {
    clim:[p.data[0].zmin,p.data[0].zmax],topEdge:p.data[0].y.at(-1),
    time:document.getElementById('status').textContent,range:p.layout.yaxis.range,
    warning:document.getElementById('warning').textContent,arrows:p.data[1].visible,tracers:p.data[3].visible,
    shapes:p.layout.shapes.length};})()`);
  const initial=await state();assert.equal(initial.topEdge,0);
  const change=async(id,value)=>{await evaluate(`document.getElementById('${id}').${typeof value==='boolean'?'checked':'value'}=${JSON.stringify(value)};document.getElementById('${id}').dispatchEvent(new Event('${id==='time'?'input':'change'}'));window.step2RenderPromise.then(()=>true)`);};
  await change('time','8');const moving=await state();assert.deepEqual(moving.clim,initial.clim);assert.match(moving.time,/0.400/);
  await change('vorticity',true);const vort=await state();assert.equal(vort.clim[0],-vort.clim[1]);
  await change('time','90');assert.deepEqual((await state()).clim,vort.clim);
  await change('contact',false);const unshaded=await state();assert.match(unshaded.warning,/not physically resolved/);assert(unshaded.shapes>=5);
  await change('arrows',false);await change('tracers',false);assert.equal((await state()).arrows,false);assert.equal((await state()).tracers,false);
  await change('view','full');const full=await state();assert.deepEqual(full.range,[-10,.04]);assert.equal(full.topEdge,0);
  await change('view','bulk');assert.deepEqual(await evaluate("document.getElementById('plot').layout.xaxis.range"),[-.9,.9]);
  await change('view','near');await change('contact',true);await change('arrows',true);await change('tracers',true);await change('speed',true);await change('time','8');
  await evaluate("document.getElementById('play').click()");await new Promise(r=>setTimeout(r,700));await evaluate("document.getElementById('play').click();window.step2RenderPromise.then(()=>true)");
  assert.notEqual((await state()).time,moving.time,'Play advances real saved frames');
  await change('time','8');
  const png=await call('Page.captureScreenshot',{format:'png',captureBeyondViewport:false});
  await mkdir(join(out,'video_inspection'),{recursive:true});
  await writeFile(join(out,'video_inspection/explorer_smoke.png'),Buffer.from(png.data,'base64'));
  assert.deepEqual(exceptions,[],'No browser runtime exceptions');
  const generated=JSON.parse(await readFile(join(out,'explorer_summary.json'),'utf8'));
  const result={status:'passed',html_sha256:generated.sha256,checks:['offline load','time slider','fixed speed scale','fixed omega scale',
    'contact warning persists','arrow/tracer toggles','full/bulk/near views','reference-domain heatmap edges','play/pause'],initial,moving,vort,full,browser_exceptions:exceptions};
  await writeFile(join(out,'explorer_inspection.json'),JSON.stringify(result,null,2)+'\n');
  console.log(JSON.stringify({status:result.status,checks:result.checks,html_sha256:generated.sha256},null,2));
}finally{
  if(ws)ws.close();child.kill('SIGTERM');
}
