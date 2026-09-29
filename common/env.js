/* ---------- 环境效应 ---------- */
(function(){
const ENV=DATA.params.ENV, IS_RICE=ENV.crop==="rice";
const f0=v=>v==null?"–":Math.abs(v)>=100?Math.round(v).toLocaleString():Math.abs(v)>=10?v.toFixed(0):v.toFixed(1);
const pctTxt=p=>`${p[1]>0?"+":""}${p[1].toFixed(0)}%`;
const ci=p=>`<small>（${p[0].toFixed(0)}~${p[2].toFixed(0)}）</small>`;
const cls=v=>v<0?"d":v>0?"u":"";
const med=a=>{a=[...a].sort((x,y)=>x-y);return a[Math.floor(a.length/2)]};

/* 汇总 */
function renderEnvTiles(){
const diffs=k=>SEASONS.map(({s})=>s.env.diff[k].pct[1]);
const benefit=SEASONS.map(({s})=>s.env.scen.FP.DMG[1]-s.env.scen.CRF.DMG[1]);
const tiles=[
  [`${med(diffs("NH3")).toFixed(0)}%`,"氨挥发"],
  [`${med(diffs("N2O")).toFixed(0)}%`,"N₂O（直接 + 间接）"],
  [`${med(diffs("LEACH")).toFixed(0)}%`,"淋溶 + 径流"],
  [`${med(diffs("GHG")).toFixed(0)}%`,"化肥相关温室气体"],
  [`${Math.round(med(benefit))}`,"社会效益 USD/ha·季"],
];
document.getElementById("envTiles").innerHTML=tiles.map(([b,t])=>`<div><b>${b}</b><span>${t}（中位数）</span></div>`).join("");
}
renderEnvTiles();

/* 全表 */
const tb=document.querySelector("#envAll tbody");
function renderEnvTable(){
  const q=document.getElementById("eq").value.trim().toLowerCase(); let n=0; tb.innerHTML="";
  SEASONS.forEach(({z,s})=>{
    if(q && !(z.country_name+z.name+s.name).toLowerCase().includes(q)) return; n++;
    const e=s.env, d=e.diff, fp=e.scen.FP, cr=e.scen.CRF;
    const tr=document.createElement("tr");
    tr.innerHTML=`<td><span class="dot" style="background:${CROPC[z.crop]}"></span>${z.country_name}</td><td>${z.name}</td><td>${s.name}</td>
      <td class="n">${s.FN} → ${s.total}</td>
      <td class="n ${cls(d.NH3.pct[1])}">${pctTxt(d.NH3.pct)}</td><td class="n ${cls(d.N2O.pct[1])}">${pctTxt(d.N2O.pct)}</td><td class="n ${cls(d.LEACH.pct[1])}">${pctTxt(d.LEACH.pct)}</td>
      <td class="n">${f0(fp.Nr[1])} → <b>${f0(cr.Nr[1])}</b></td>
      <td class="n">${e.decomp.rate>0?"+":""}${e.decomp.rate.toFixed(1)} / ${e.decomp.product.toFixed(1)}</td>
      <td class="n">${f0(fp.CF[1])} → ${f0(cr.CF[1])}</td>
      <td class="n"><b>${Math.round(fp.DMG[1]-cr.DMG[1])}</b></td>`;
    tr.onclick=()=>selectZone(z.id,s.idx,true);
    tb.appendChild(tr);
  });
  document.getElementById("eqCount").textContent=`${n} / ${SEASONS.length}`;
}
document.getElementById("eq").addEventListener("input",renderEnvTable);
renderEnvTable();
window.renderEnvAll=()=>{ renderEnvTiles(); renderEnvTable(); };

/* 参数表 */
const R=ENV.RR, rr=k=>`${Math.round(R[k][0]*100)}%（${Math.round(R[k][1]*100)}–${Math.round(R[k][2]*100)}）`;
const rows=[
  ["控释尿素减 NH₃（同氮量）",rr("NH3")],["控释尿素减 N₂O",rr("N2O")],["控释尿素减淋溶径流",rr("LEACH")],
  ...(IS_RICE?[["控释尿素减 CH₄（单列）",rr("CH4")]]:[]),
  ["尿素 NH₃ 挥发比例（中性土）",`${Math.round(ENV.NH3_BASE[0]*100)}%，酸性土 ×0.7，碱性土 ×1.4`],
  ["直接 N₂O EF1",IS_RICE?"淹水稻田 0.4%":"湿润/灌溉 1.6%，干旱 0.5%"],
  ["间接 N₂O EF4 / EF5","湿润 1.4%、干旱 0.5% / 1.1%"],
  ["淋溶径流比例",IS_RICE?"稻田 10%（5–18%）":"湿润/灌溉 20%，干旱 5%"],
  ["化肥生产（kg CO₂e/kg N）",`中国 ${ENV.FERT_PROD.china}，其他 ${ENV.FERT_PROD.other}；包膜附加 ${ENV.COATING_EXTRA}`],
  ["损害成本 NH₃ / 水体（USD/kg N）",`${ENV.DAMAGE.NH3[0]}（${ENV.DAMAGE.NH3[1]}–${ENV.DAMAGE.NH3[2]}） / ${ENV.DAMAGE.LEACH[0]}（${ENV.DAMAGE.LEACH[1]}–${ENV.DAMAGE.LEACH[2]}）`],
  ["碳价（USD/t CO₂e）",`${ENV.CARBON_PRICE[0]}（${ENV.CARBON_PRICE[1]}–${ENV.CARBON_PRICE[2]}）`],
  ["GWP₁₀₀ N₂O / CH₄","273 / 27（IPCC AR6）"],
  ["Monte Carlo 次数",ENV.N_MC],
];
document.getElementById("envPtbl").innerHTML=`<tr><th>参数</th><th class="n">取值（范围）</th></tr>`+rows.map(([a,b])=>`<tr><td>${a}</td><td class="n" style="white-space:normal">${b}</td></tr>`).join("");

/* 详情面板中的环境卡片 */
const box=document.createElement("div"); box.className="chart-box"; box.id="envBox";
document.querySelector(".charts").appendChild(box);
let area=10;
function renderEnv(s){
  const e=s.env, S=e.scen, D=e.diff;
  const line=(label,k,unit,digits)=>{
    const d=D[k]?D[k].pct:null;
    return `<tr><td>${label}<small> ${unit}</small></td><td>${f0(S.FP[k][1])}</td><td>${f0(S.OPT[k][1])}</td><td><b>${f0(S.CRF[k][1])}</b></td>`+
      (d?`<td class="${cls(d[1])}">${pctTxt(d)} ${ci(d)}</td>`:`<td></td>`)+`</tr>`;
  };
  const tot=e.decomp.rate+e.decomp.product, W=Math.max(Math.abs(e.decomp.rate)+Math.abs(e.decomp.product),1e-6);
  const benefit=[S.FP.DMG[0]-S.CRF.DMG[2], S.FP.DMG[1]-S.CRF.DMG[1], S.FP.DMG[2]-S.CRF.DMG[0]];
  box.innerHTML=`
    <h3>环境效应（每公顷每季）</h3>
    <p class="note" style="margin-top:2px">包膜尿素占 ${Math.round(e.c_crf*100)}% · 土壤 ${({acid:"酸性",neutral:"中性",alkaline:"碱性"})[e.ph]} · ${e.climate==="wet"?"湿润/灌溉":"干旱"}气候 · 化肥产地 ${e.prod==="china"?"中国":"其他"}</p>
    <div class="scroll" style="border:0;background:transparent">
    <table class="envtbl">
      <tr><th></th><th>农户常规</th><th>分次施尿素</th><th>控释掺混</th><th>控释比常规</th></tr>
      <tr><td>施氮量<small> kg N</small></td><td>${S.FP.N[1]}</td><td>${S.OPT.N[1]}</td><td><b>${S.CRF.N[1]}</b></td><td></td></tr>
      ${line("NH₃ 挥发","NH3","kg N",1)}
      ${line("N₂O（直接+间接）","N2O","kg N",2)}
      ${line("淋溶 + 径流","LEACH","kg N",1)}
      ${line("Nr 损失合计","Nr","kg N",1)}
      ${line("温室气体（化肥生产 + N₂O）","GHG","kg CO₂e",0)}
      ${IS_RICE?line("CH₄（证据较弱，单列）","CH4","kg CH₄",0):""}
      ${line("氮足迹","NF","kg Nr/t",1)}
      ${line("碳足迹","CF","kg CO₂e/t",0)}
      ${line("社会损害成本","DMG","USD",0)}
    </table></div>
    <p class="note" style="margin-top:8px">Nr 损失变化 <b class="num">${tot>0?"+":""}${tot.toFixed(1)}</b> kg N/ha：减量 <span class="num">${e.decomp.rate>0?"+":""}${e.decomp.rate.toFixed(1)}</span>，控释 <span class="num">${e.decomp.product.toFixed(1)}</span></p>
    <div class="decomp" title="减量效应 / 控释效应">
      <div style="width:${100*Math.abs(e.decomp.rate)/W}%;background:${e.decomp.rate<=0?"var(--urea)":"var(--bad)"}"></div>
      <div style="width:${100*Math.abs(e.decomp.product)/W}%;background:var(--cr1)"></div>
    </div>
    <div class="chart-legend"><span><i style="background:var(--urea)"></i>减量效应</span><span><i style="background:var(--cr1)"></i>控释效应</span></div>
    <div class="area-row"><label for="envArea">推广面积</label><input id="envArea" type="number" min="0" step="1" value="${area}"><span class="muted">万公顷 · 每季合计</span></div>
    <div class="area-out" id="areaOut"></div>`;
  const upd=()=>{
    const A=Math.max(0,+document.getElementById("envArea").value||0); area=A; const ha=A*1e4;
    const dd=k=>(S.FP[k][1]-S.CRF[k][1])*ha;
    document.getElementById("areaOut").innerHTML=[
      ["少施氮",(S.FP.N[1]-S.CRF.N[1])*ha/1000,"t N"],
      ["少挥发 NH₃",dd("NH3")/1000,"t N"],
      ["少排 N₂O",dd("N2O")/1000,"t N"],
      ["少淋溶径流",dd("LEACH")/1000,"t N"],
      ["少排温室气体",dd("GHG")/1e7,"万 t CO₂e"],
      ["社会效益",benefit[1]*ha/1e6,`百万 USD（${(benefit[0]*ha/1e6).toFixed(0)}~${(benefit[2]*ha/1e6).toFixed(0)}）`],
    ].map(([a,v,u])=>`<div>${a}<b>${v.toLocaleString(undefined,{maximumFractionDigits:v<10?1:0})}</b>${u}</div>`).join("");
  };
  document.getElementById("envArea").oninput=upd; upd();
}
const _rd=renderDetail;
renderDetail=function(){ _rd(); renderEnv(selZone.seasons[selSeason]); };
const _zt=zoneTip;
zoneTip=function(z){
  const extra=z.seasons.map(s=>{const d=s.env.diff;return `<div style="font-size:12px"><span class="muted">${s.name} 比常规：</span>NH₃ ${pctTxt(d.NH3.pct)} · N₂O ${pctTxt(d.N2O.pct)} · 淋溶 ${pctTxt(d.LEACH.pct)}</div>`}).join("");
  return _zt(z).replace(`<div class="muted" style="margin-top:6px;font-size:12px">`,`<div class="s">${extra}</div><div class="muted" style="margin-top:6px;font-size:12px">`);
};
})();
