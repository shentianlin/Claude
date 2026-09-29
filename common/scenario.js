/* ---------- 产品方案情景：线性型（现有）/ S 型单档 / S 型两档 ---------- */
(function(){
const SC=DATA.params.SCEN, KEYS=["lin","s1","s2"];
const FIELDS=["rec","rec_str","total","simple","simple_str","simple_total","ideal","kg_cr","kg_urea","kg_product","RE","unreleased","sens","robust","env"];
// 页面数据默认就是线性型方案：先把它存成 alt.lin，切换时整体替换字段
SEASONS.forEach(({s})=>{
  const lin={beta:SC.lin_beta, series_comps:s.series.comps, series_pool:s.series.pool};
  FIELDS.forEach(f=>lin[f]=s[f]);
  s.alt=Object.assign({lin}, s.alt);
});
const NOTES={
  lin:`现有线性型产品（β ${SC.lin_beta}）：尿素 + 1–2 个控释期（CR 档）。`,
  s1:`假设的 S 型产品（β ${SC.stype_beta}，有滞后期；控释期同样按 25 °C 静水 80% 释放天数标定）：尿素 + 1 个 S 型档位。尚无实测样品，仅供评估。`,
  s2:`假设的 S 型产品（β ${SC.stype_beta}）：尿素 + 至多 2 个 S 型档位；取整后若只剩一档，则与单档相同。尚无实测样品，仅供评估。`,
};
const segEl=document.getElementById("scenSeg");
segEl.innerHTML=KEYS.map(k=>`<button type="button" data-k="${k}">${SC.labels[k]}</button>`).join("");
segEl.querySelectorAll("button").forEach(b=>b.addEventListener("click",()=>applyScenario(b.dataset.k,true)));

window.applyScenario=function(k, rerender){
  if(!KEYS.includes(k)) k="lin";
  SCEN=k;
  SEASONS.forEach(({s})=>{
    const a=s.alt[k];
    FIELDS.forEach(f=>{ if(f in a) s[f]=a[f]; });
    s.series.comps=a.series_comps; s.series.pool=a.series_pool;
  });
  segEl.querySelectorAll("button").forEach(b=>b.setAttribute("aria-pressed", String(b.dataset.k===k)));
  document.getElementById("scenNote").textContent=NOTES[k];
  try{ localStorage.setItem("crf-scenario",k); }catch(e){}
  try{ history.replaceState(null,"",k==="lin"?location.pathname+location.search:"#"+k); }catch(e){}
  if(!rerender) return;
  renderFindings(); renderTable(); recomputeUsed();
  exB=SEASONS[0].s.alt[k].beta;
  seg2Beta();
  drawC3();
  if(window.renderEnvAll) renderEnvAll();
  if(selZone) renderDetail();
};
function seg2Beta(){ seg($("#exB"),[1.0,1.3,2.0,2.5],exB,v=>v.toFixed(1),v=>exB=v); }

// 详情面板：标出当前方案
const _rp=renderPanel;
renderPanel=function(z,s){
  _rp(z,s);
  const tag=document.createElement("div");
  tag.style.marginBottom="8px";
  tag.innerHTML=`<span class="pill ${SCEN==="lin"?"":"warn"}">${SC.labels[SCEN]}</span>`+
    (SCEN==="s2"&&s.alt.s2.same_as_s1?`<span class="pill">两档取整后与单档相同</span>`:"")+
    (SCEN!=="lin"?`<span class="pill good">比线性型省 ${Math.round(100*(1-s.total/s.alt.lin.total))}%</span>`:"");
  const panel=document.getElementById("dPanel"); panel.insertBefore(tag, panel.firstChild);
};

// 初始方案：链接中的 #s1 / #s2，其次是本机上次的选择
let init="lin";
try{ init=localStorage.getItem("crf-scenario")||"lin"; }catch(e){}
const h=(location.hash||"").slice(1); if(KEYS.includes(h)) init=h;
applyScenario(init,false);
if(init!=="lin"){ exB=SEASONS[0].s.alt[init].beta; }
seg2Beta(); recomputeUsed(); renderFindings(); if(window.renderEnvAll) renderEnvAll();
})();
