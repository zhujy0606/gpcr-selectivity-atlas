const state={summary:null,receptors:[],pairs:[],seeds:[],compounds:[]};
const $=selector=>document.querySelector(selector);
const esc=value=>String(value??'—').replace(/[&<>"']/g,char=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[char]));
const fmt=(value,digits=2)=>value===null||value===undefined||value===''?'—':typeof value==='number'?value.toLocaleString('en-US',{maximumFractionDigits:digits}):esc(value);
const includes=(values,query)=>values.join(' ').toLowerCase().includes(query.trim().toLowerCase());
const badge=(text,kind='')=>`<span class="badge ${kind}">${esc(text)}</span>`;
const boolBadge=(value,yes='是',no='否')=>`<span class="badge ${value?'good':''}">${value?yes:no}</span>`;
const table=(headers,rows,empty='无匹配记录')=>`<table><thead><tr>${headers.map(header=>`<th>${header}</th>`).join('')}</tr></thead><tbody>${rows.length?rows.join(''):`<tr><td class="empty-cell" colspan="${headers.length}">${empty}</td></tr>`}</tbody></table>`;
const evidenceLabel=compound=>compound.evidence.strict_final_selected_111?badge('严格精选111','final'):compound.evidence.mmgbsa_dual_endpoint_improved?badge('双端MM/GBSA改善141','good'):compound.evidence.mmgbsa_baseline_complete?badge('基线完整323'):compound.evidence.cross_domain_shortlist?badge('跨域短名单438','warn'):badge('robust最终候选904');

async function load(){
  const names=['summary','receptors','pairs','seeds','compounds'];
  const data=await Promise.all(names.map(name=>fetch(`data/${name}.json`).then(response=>{
    if(!response.ok)throw new Error(`${name}: ${response.status}`);
    return response.json();
  })));
  [state.summary,state.receptors,state.pairs,state.seeds,state.compounds]=data;
  renderReceptors();
  renderPairs();
  $('#footer-version').textContent=`${state.summary.version} · ${state.summary.build_date}`;
}

function renderReceptors(){
  const query=$('#receptor-search').value;
  const items=state.receptors.filter(receptor=>includes([receptor.uniprot,receptor.name,receptor.subfamily],query));
  $('#receptor-count').textContent=`${items.length} / ${state.receptors.length} 个受体`;
  $('#receptor-table').innerHTML=table(
    ['UniProt','受体名称','亚家族','dMaSIF资产','核心资产','最近表面邻居','全局距离','进入163对','详情'],
    items.map(receptor=>`<tr class="clickable" data-receptor="${receptor.uniprot}">
      <td class="mono"><b>${receptor.uniprot}</b></td>
      <td>${esc(receptor.name||'未注释')}</td>
      <td>${esc(receptor.subfamily||'—')}</td>
      <td><b>${receptor.surface_asset_count}</b><br><small>${(receptor.surface_asset_bytes/1048576).toFixed(1)} MB</small></td>
      <td>${boolBadge(receptor.core_assets_complete,'完整','缺失')}</td>
      <td class="mono">${esc(receptor.nearest_surface_neighbor?.uniprot)}</td>
      <td>${fmt(receptor.nearest_surface_neighbor?.distance,4)}</td>
      <td>${receptor.selected_pair_count}</td>
      <td><span class="row-action">查看资产表 →</span></td>
    </tr>`)
  );
  document.querySelectorAll('[data-receptor]').forEach(row=>{
    row.onclick=()=>showReceptor(state.receptors.find(receptor=>receptor.uniprot===row.dataset.receptor));
  });
}

function showReceptor(receptor){
  const assetRows=receptor.assets.map(asset=>`<tr>
    <td class="mono">${esc(asset.filename)}</td>
    <td>${asset.kind==='npy'?'NumPy数组':asset.kind==='vtk'?'VTK表面':'其他'}</td>
    <td>${asset.kind==='npy'&&asset.filename.includes('coords')?'表面点坐标':asset.filename.includes('features_emb1')?'dMaSIF embedding 1':asset.filename.includes('features_emb2')?'dMaSIF embedding 2':asset.filename.includes('emb1')?'embedding 1 可视化':'embedding 2 可视化'}</td>
    <td>${asset.size_bytes?(asset.size_bytes/1048576).toFixed(2)+' MB':'清单已登记'}</td>
    <td>${badge('冻结资产','good')}</td>
  </tr>`);
  showDialog(`
    <div class="detail-head receptor-detail-head">
      <p class="eyebrow">MODULE 01 · RECEPTOR SURFACE RECORD</p>
      <h2>${esc(receptor.name||receptor.uniprot)}</h2>
      <p class="mono">${receptor.uniprot}</p>
    </div>
    <div class="pair-summary-strip receptor-summary-strip">
      <span><small>dMaSIF资产</small><b>${receptor.surface_asset_count}</b></span>
      <span><small>冻结大小</small><b>${(receptor.surface_asset_bytes/1048576).toFixed(1)} MB</b></span>
      <span><small>最近邻</small><b class="mono">${esc(receptor.nearest_surface_neighbor?.uniprot)}</b></span>
      <span><small>全局距离</small><b>${fmt(receptor.nearest_surface_neighbor?.distance,4)}</b></span>
      <span><small>进入受体对</small><b>${receptor.selected_pair_count}</b></span>
    </div>
    <div class="evidence-section-head">
      <div><p class="eyebrow">SYSTEMATIC ASSET TABLE</p><h3>dMaSIF文件清单</h3></div>
      <span>${receptor.core_assets_complete?'核心资产完整':'核心资产存在缺失'}</span>
    </div>
    <div class="data-table evidence-table">${table(['文件名','格式','数据内容','大小','状态'],assetRows)}</div>
    <aside class="notice amber detail-note"><strong>公开状态</strong><span>文件级清单已公开；大型原始二进制资产仍标记为 ${esc(receptor.public_raw_asset_status)}，待外部公开档案固定永久下载地址。</span></aside>
  `);
}

function pairMolecules(pairId){return state.compounds.filter(compound=>compound.pair_id===pairId)}
function pairSeeds(pairId){return state.seeds.filter(seed=>seed.pair_id===pairId)}

function renderPairs(){
  const query=$('#pair-search').value;
  const items=state.pairs.filter(pair=>includes([
    pair.pair_id,pair.receptor_a.uniprot,pair.receptor_a.name,pair.receptor_b.uniprot,pair.receptor_b.name,
    ...pair.input_seed_zinc_ids,...pair.hotspots.flatMap(hotspot=>[hotspot.bw,hotspot.residues])
  ],query));
  $('#pair-count').textContent=`${items.length} / ${state.pairs.length} 对受体`;
  $('#pair-table').innerHTML=table(
    ['排名','受体A','受体B','dMaSIF/MaSIF距离','Top 3差异热点','输入种子任务','robust生成分子','严格精选','关联表'],
    items.map(pair=>`<tr class="clickable" data-pair="${pair.pair_id}">
      <td><b>#${pair.rank}</b></td>
      <td><b>${esc(pair.receptor_a.name)}</b><br><span class="mono">${pair.receptor_a.uniprot}</span></td>
      <td><b>${esc(pair.receptor_b.name)}</b><br><span class="mono">${pair.receptor_b.uniprot}</span></td>
      <td><b>${fmt(pair.surface_distance,3)}</b></td>
      <td>${pair.hotspots.map(hotspot=>badge(hotspot.bw,'good')).join(' ')}</td>
      <td><b>${pairSeeds(pair.pair_id).length}</b><br><small>${pair.input_seed_zinc_ids.length}个ZINC</small></td>
      <td><b class="molecule-count">${pairMolecules(pair.pair_id).length}</b></td>
      <td>${pair.strict_final_selected_count}</td>
      <td><span class="row-action">打开三张表 →</span></td>
    </tr>`)
  );
  document.querySelectorAll('[data-pair]').forEach(row=>{
    row.onclick=()=>showPair(state.pairs.find(pair=>pair.pair_id===row.dataset.pair));
  });
  const linked=state.pairs.reduce((total,pair)=>total+pairMolecules(pair.pair_id).length,0);
  const orphan=state.compounds.filter(compound=>!state.pairs.some(pair=>pair.pair_id===compound.pair_id)).length;
  $('#pair-total-check').textContent=linked===904&&orphan===0
    ?'已核对：163对受体下关联分子合计904个，0个游离记录。'
    :`当前关联分子${linked}个，未匹配受体对${orphan}个；请查看数据版本。`;
}

function hotspotPanel(pair){
  const rows=pair.hotspots.map(hotspot=>`<tr>
    <td><b>#${hotspot.rank}</b></td>
    <td>${badge(hotspot.bw,'good')}</td>
    <td class="mono">${esc(hotspot.residues)}</td>
    <td><b>${fmt(hotspot.fingerprint_difference,3)}</b></td>
    <td class="mono">${esc(hotspot.hotspot_id)}</td>
  </tr>`);
  return `<div class="table-explainer"><b>Top 3局部差异热点</b><span>按局部表面指纹差异排序；BW为Ballesteros–Weinstein通用编号。</span></div><div class="data-table evidence-table">${table(['热点排名','BW位点','两受体残基','局部指纹差异Δ','热点ID'],rows)}</div>`;
}

function seedPanel(pair,seeds){
  const rows=seeds.map(seed=>`<tr>
    <td><span class="seed-chip"><small>INPUT SEED</small><b class="mono">${esc(seed.seed_zinc_id)}</b></span></td>
    <td>${esc(seed.target_name||seed.target_uniprot)} → ${esc(seed.offtarget_name||seed.offtarget_uniprot)}<br><small class="mono">${esc(seed.task)}</small></td>
    <td>${badge(seed.hotspot_bw,'good')}</td>
    <td>${fmt(seed.fast.target,3)}</td>
    <td>${fmt(seed.fast.offtarget,3)}</td>
    <td><b>${fmt(seed.fast.dd,3)}</b></td>
    <td>${fmt(seed.detail.target_median,3)}</td>
    <td>${fmt(seed.detail.offtarget_median,3)}</td>
    <td><b>${fmt(seed.detail.dd_median,3)}</b></td>
    <td>${fmt(seed.detail.dd_worst,3)}</td>
    <td>${fmt(seed.detail.dd_sd,3)}</td>
    <td>${boolBadge(seed.detail.target_pose_stable&&seed.detail.offtarget_pose_stable,'双端稳定','未双稳')}</td>
    <td><b>${seed.generated_compound_count}</b></td>
  </tr>`);
  return `<div class="table-explainer"><b>PocketXMol输入种子</b><span>每一行是一条有向种子任务；保留fast mode及detail mode三次重复汇总。</span></div><div class="data-table evidence-table wide-table">${table(['输入种子ZINC','选择性方向','热点','fast目标','fast脱靶','fast DD','detail目标中位数','detail脱靶中位数','detail DD中位数','worst DD','DD SD','姿势稳定','生成分子数'],rows,'该受体对在904集合中没有关联的输入种子任务')}</div>`;
}

function compoundPanel(pair,compounds){
  const rows=compounds.map(compound=>{
    const structure=compound.structure_download;
    const structureText=structure.complex_pdb_count?`1 SDF + ${structure.complex_pdb_count} PDB`:'1 SDF';
    return `<tr>
      <td><b class="mono">${compound.compound_id}</b><span class="truncate mono molecule-smiles" title="${esc(compound.canonical_smiles)}">${esc(compound.canonical_smiles)}</span></td>
      <td>${esc(compound.target_name||compound.target_uniprot)} → ${esc(compound.offtarget_name||compound.offtarget_uniprot)}</td>
      <td><span class="seed-chip compact-seed"><small>INPUT SEED</small><b class="mono">${esc(compound.seed_zinc_id)}</b></span></td>
      <td>${fmt(compound.similarity_to_seed,3)}</td>
      <td>${fmt(compound.properties.MW,1)}</td>
      <td>${fmt(compound.properties.cLogP,2)}</td>
      <td>${fmt(compound.properties.QED,2)}</td>
      <td>${fmt(compound.properties.SA,2)}</td>
      <td><b>${fmt(compound.docking.dd_median,3)}</b></td>
      <td>${fmt(compound.docking.dd_worst,3)}</td>
      <td><b>${fmt(compound.docking.dd_change_vs_seed,3)}</b></td>
      <td>${boolBadge(compound.docking.both_pose_stable,'双端稳定','未双稳')}</td>
      <td>${evidenceLabel(compound)}</td>
      <td><a class="structure-link" href="${esc(structure.bundle_url)}" download title="下载${compound.compound_id}结构包">${structureText} ↓</a></td>
    </tr>`;
  });
  return `<div class="table-explainer"><b>robust选择性生成分子</b><span>这里展示该受体对在904个最终候选集合中的全部分子；结构包含配体SDF，438子集另含计算复合物PDB。</span></div><div class="data-table evidence-table wide-table molecule-table">${table(['生成分子 / SMILES','选择性方向','输入种子ZINC','与种子相似度','MW','cLogP','QED','SA','detail DD','worst DD','ΔDD vs seed','姿势','下游证据','结构下载'],rows,'该受体对在904集合中没有robust生成分子')}</div>`;
}

function showPair(pair){
  const seeds=pairSeeds(pair.pair_id);
  const compounds=pairMolecules(pair.pair_id);
  showDialog(`
    <div class="detail-head pair-detail-head">
      <p class="eyebrow">MODULE 02 · RECEPTOR PAIR · RANK ${pair.rank}</p>
      <h2>${esc(pair.receptor_a.name)} / ${esc(pair.receptor_b.name)}</h2>
      <p class="mono">${pair.pair_id}</p>
    </div>
    <div class="pair-summary-strip">
      <span><small>dMaSIF/MaSIF距离</small><b>${fmt(pair.surface_distance,3)}</b></span>
      <span><small>Top差异热点</small><b>${esc(pair.hotspots[0]?.bw)}</b></span>
      <span><small>输入种子任务</small><b>${seeds.length}</b></span>
      <span><small>robust生成分子</small><b>${compounds.length}</b></span>
      <span><small>严格精选111子集</small><b>${pair.strict_final_selected_count}</b></span>
    </div>
    <div class="evidence-tabs" role="tablist" aria-label="受体对关联数据表">
      <button class="evidence-tab active" role="tab" aria-selected="true" data-panel-target="hotspots">Top 3热点 <b>${pair.hotspots.length}</b></button>
      <button class="evidence-tab" role="tab" aria-selected="false" data-panel-target="seeds">输入种子 <b>${seeds.length}</b></button>
      <button class="evidence-tab" role="tab" aria-selected="false" data-panel-target="molecules">robust生成分子 <b>${compounds.length}</b></button>
    </div>
    <section class="evidence-panel active" data-panel="hotspots">${hotspotPanel(pair)}</section>
    <section class="evidence-panel" data-panel="seeds" hidden>${seedPanel(pair,seeds)}</section>
    <section class="evidence-panel" data-panel="molecules" hidden>${compoundPanel(pair,compounds)}</section>
  `);
  document.querySelectorAll('[data-panel-target]').forEach(button=>{
    button.onclick=()=>activateEvidencePanel(button.dataset.panelTarget);
  });
}

function activateEvidencePanel(panelName){
  document.querySelectorAll('[data-panel-target]').forEach(button=>{
    const active=button.dataset.panelTarget===panelName;
    button.classList.toggle('active',active);
    button.setAttribute('aria-selected',String(active));
  });
  document.querySelectorAll('[data-panel]').forEach(panel=>{
    const active=panel.dataset.panel===panelName;
    panel.classList.toggle('active',active);
    panel.hidden=!active;
  });
}

function showDialog(html){
  $('#detail-content').innerHTML=html;
  $('#detail-dialog').showModal();
}

function route(){
  const legacyMap={surface:'receptors',seeds:'pairs',compounds:'pairs',structures:'overview',downloads:'overview'};
  const requested=location.hash.slice(1)||'overview';
  const id=legacyMap[requested]||requested;
  const valid=['overview','receptors','pairs'].includes(id)?id:'overview';
  document.querySelectorAll('.page').forEach(page=>page.classList.toggle('active',page.id===valid));
  document.querySelectorAll('nav a').forEach(link=>link.classList.toggle('active',link.getAttribute('href')===`#${valid}`));
  $('#nav').classList.remove('open');
  $('#menu-button').setAttribute('aria-expanded','false');
  window.scrollTo(0,0);
}

window.addEventListener('hashchange',route);
$('#menu-button').onclick=()=>{
  const open=$('#nav').classList.toggle('open');
  $('#menu-button').setAttribute('aria-expanded',String(open));
};
$('#dialog-close').onclick=()=>$('#detail-dialog').close();
$('#detail-dialog').addEventListener('click',event=>{if(event.target===$('#detail-dialog'))$('#detail-dialog').close()});
$('#receptor-search').addEventListener('input',renderReceptors);
$('#pair-search').addEventListener('input',renderPairs);
route();
load().catch(error=>{
  document.querySelectorAll('.data-table').forEach(element=>{
    if(!element.innerHTML)element.innerHTML=`<div class="load-error"><b>数据载入失败</b><span>${esc(error.message)}。请通过网站地址或本地HTTP服务器访问，不要直接双击HTML文件。</span></div>`;
  });
  $('#pair-total-check').textContent='数据尚未载入，无法核对904条关联。';
  console.error(error);
});
