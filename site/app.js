const state={summary:null,receptors:[],pairs:[],seeds:[],compounds:[],detailModeCompounds:[]};
const $=selector=>document.querySelector(selector);
const esc=value=>String(value??'—').replace(/[&<>"']/g,char=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[char]));
const fmt=(value,digits=2)=>value===null||value===undefined||value===''?'—':typeof value==='number'?value.toLocaleString('en-US',{maximumFractionDigits:digits}):esc(value);
const includes=(values,query)=>values.join(' ').toLowerCase().includes(query.trim().toLowerCase());
const badge=(text,kind='')=>`<span class="badge ${kind}">${esc(text)}</span>`;
const boolBadge=(value,yes='Yes',no='No')=>`<span class="badge ${value?'good':''}">${value?yes:no}</span>`;
const table=(headers,rows,empty='No matching records')=>`<table><thead><tr>${headers.map(header=>`<th>${header}</th>`).join('')}</tr></thead><tbody>${rows.length?rows.join(''):`<tr><td class="empty-cell" colspan="${headers.length}">${empty}</td></tr>`}</tbody></table>`;

async function load(){
  const names=['summary','receptors','pairs','seeds','compounds','detail_mode_compounds'];
  const data=await Promise.all(names.map(name=>fetch(`data/${name}.json`).then(response=>{
    if(!response.ok)throw new Error(`${name}: ${response.status}`);
    return response.json();
  })));
  [state.summary,state.receptors,state.pairs,state.seeds,state.compounds,state.detailModeCompounds]=data;
  renderReceptors();
  renderPairs();
  $('#footer-version').textContent=`${state.summary.version} · ${state.summary.build_date}`;
}

function renderReceptors(){
  const query=$('#receptor-search').value;
  const items=state.receptors.filter(receptor=>includes([receptor.uniprot,receptor.name,receptor.subfamily],query));
  $('#receptor-count').textContent=`${items.length} / ${state.receptors.length} receptors`;
  $('#receptor-table').innerHTML=table(
    ['UniProt','Receptor','Subfamily','dMaSIF assets','Core assets','Nearest surface neighbor','Global distance','Selected pairs','Record'],
    items.map(receptor=>`<tr class="clickable" data-receptor="${receptor.uniprot}">
      <td class="mono"><b>${receptor.uniprot}</b></td>
      <td>${esc(receptor.name||'Unannotated')}</td>
      <td>${esc(receptor.subfamily||'—')}</td>
      <td><b>${receptor.surface_asset_count}</b><br><small>${(receptor.surface_asset_bytes/1048576).toFixed(1)} MB</small></td>
      <td>${boolBadge(receptor.core_assets_complete,'Complete','Missing')}</td>
      <td class="mono">${esc(receptor.nearest_surface_neighbor?.uniprot)}</td>
      <td>${fmt(receptor.nearest_surface_neighbor?.distance,4)}</td>
      <td>${receptor.selected_pair_count}</td>
      <td><span class="row-action">View asset inventory →</span></td>
    </tr>`)
  );
  document.querySelectorAll('[data-receptor]').forEach(row=>{
    row.onclick=()=>showReceptor(state.receptors.find(receptor=>receptor.uniprot===row.dataset.receptor));
  });
}

function showReceptor(receptor){
  const assetRows=receptor.assets.map(asset=>`<tr>
    <td class="mono">${esc(asset.filename)}</td>
    <td>${asset.kind==='npy'?'NumPy array':asset.kind==='vtk'?'VTK surface':'Other'}</td>
    <td>${asset.kind==='npy'&&asset.filename.includes('coords')?'Surface point coordinates':asset.filename.includes('features_emb1')?'dMaSIF embedding 1':asset.filename.includes('features_emb2')?'dMaSIF embedding 2':asset.filename.includes('emb1')?'Embedding 1 visualization':'Embedding 2 visualization'}</td>
    <td>${asset.size_bytes?(asset.size_bytes/1048576).toFixed(2)+' MB':'Manifest entry'}</td>
    <td>${badge('Frozen asset','good')}</td>
  </tr>`);
  showDialog(`
    <div class="detail-head receptor-detail-head">
      <p class="eyebrow">MODULE 01 · RECEPTOR SURFACE RECORD</p>
      <h2>${esc(receptor.name||receptor.uniprot)}</h2>
      <p class="mono">${receptor.uniprot}</p>
    </div>
    <div class="pair-summary-strip receptor-summary-strip">
      <span><small>dMaSIF assets</small><b>${receptor.surface_asset_count}</b></span>
      <span><small>Frozen size</small><b>${(receptor.surface_asset_bytes/1048576).toFixed(1)} MB</b></span>
      <span><small>Nearest neighbor</small><b class="mono">${esc(receptor.nearest_surface_neighbor?.uniprot)}</b></span>
      <span><small>Global distance</small><b>${fmt(receptor.nearest_surface_neighbor?.distance,4)}</b></span>
      <span><small>Selected pairs</small><b>${receptor.selected_pair_count}</b></span>
    </div>
    <div class="evidence-section-head">
      <div><p class="eyebrow">SYSTEMATIC ASSET TABLE</p><h3>dMaSIF file inventory</h3></div>
      <span>${receptor.core_assets_complete?'Core assets complete':'Core assets missing'}</span>
    </div>
    <div class="data-table evidence-table">${table(['File','Format','Contents','Size','Status'],assetRows)}</div>
    <aside class="notice amber detail-note"><strong>Availability</strong><span>The file-level inventory is public. Large binary assets remain pending deposition in an external archive with permanent download links.</span></aside>
  `);
}

function pairMolecules(pairId){return state.compounds.filter(compound=>compound.pair_id===pairId)}
function pairSeeds(pairId){return state.seeds.filter(seed=>seed.pair_id===pairId)}
function pairDetailModeCompounds(pairId){return state.detailModeCompounds.filter(compound=>compound.pair_id===pairId)}

function renderPairs(){
  const query=$('#pair-search').value;
  const items=state.pairs.filter(pair=>includes([
    pair.pair_id,pair.receptor_a.uniprot,pair.receptor_a.name,pair.receptor_b.uniprot,pair.receptor_b.name,
    ...pair.input_seed_zinc_ids,...pair.detail_mode_selected_zinc_ids,...pair.hotspots.flatMap(hotspot=>[hotspot.bw,hotspot.residues])
  ],query));
  $('#pair-count').textContent=`${items.length} / ${state.pairs.length} receptor pairs`;
  $('#pair-table').innerHTML=table(
    ['Rank','Receptor A','Receptor B','dMaSIF/MaSIF distance','Top 3 differential hotspots','Input seed records','Pocketxmol-generated compounds','Record'],
    items.map(pair=>`<tr class="clickable" data-pair="${pair.pair_id}">
      <td><b>#${pair.rank}</b></td>
      <td><b>${esc(pair.receptor_a.name)}</b><br><span class="mono">${pair.receptor_a.uniprot}</span></td>
      <td><b>${esc(pair.receptor_b.name)}</b><br><span class="mono">${pair.receptor_b.uniprot}</span></td>
      <td><b>${fmt(pair.surface_distance,3)}</b></td>
      <td>${pair.hotspots.map(hotspot=>badge(hotspot.bw,'good')).join(' ')}</td>
      <td><b>${pairSeeds(pair.pair_id).length}</b><br><small>${pair.input_seed_zinc_ids.length} ZINC IDs</small></td>
      <td><b class="molecule-count">${pairMolecules(pair.pair_id).length}</b>${pairMolecules(pair.pair_id).length===0?`<br><small>${pairDetailModeCompounds(pair.pair_id).length} detail-mode selected</small>`:''}</td>
      <td><span class="row-action">Open evidence tables →</span></td>
    </tr>`)
  );
  document.querySelectorAll('[data-pair]').forEach(row=>{
    row.onclick=()=>showPair(state.pairs.find(pair=>pair.pair_id===row.dataset.pair));
  });
  const linked=state.pairs.reduce((total,pair)=>total+pairMolecules(pair.pair_id).length,0);
  const orphan=state.compounds.filter(compound=>!state.pairs.some(pair=>pair.pair_id===compound.pair_id)).length;
  $('#pair-total-check').textContent=linked===904&&orphan===0
    ?'Validated: all 904 compounds are linked to the 163 receptor pairs; no orphan records.'
    :`${linked} compounds are linked and ${orphan} records are unmatched. Check the current data release.`;
}

function hotspotPanel(pair){
  const rows=pair.hotspots.map(hotspot=>`<tr>
    <td><b>#${hotspot.rank}</b></td>
    <td>${badge(hotspot.bw,'good')}</td>
    <td class="mono">${esc(hotspot.residues)}</td>
    <td><b>${fmt(hotspot.fingerprint_difference,3)}</b></td>
    <td class="mono">${esc(hotspot.hotspot_id)}</td>
  </tr>`);
  return `<div class="table-explainer"><b>Top 3 differential hotspots</b><span>Ranked by local surface-fingerprint difference. BW denotes Ballesteros–Weinstein numbering.</span></div><div class="data-table evidence-table">${table(['Hotspot rank','BW position','Paired residues','Local fingerprint difference Δ','Hotspot ID'],rows)}</div>`;
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
    <td>${boolBadge(seed.detail.target_pose_stable&&seed.detail.offtarget_pose_stable,'Stable at both receptors','Not stable at both')}</td>
    <td><b>${seed.generated_compound_count}</b></td>
  </tr>`);
  return `<div class="table-explainer"><b>Pocketxmol input seeds</b><span>Each row represents a directed seed task with fast-mode results and the three-repeat detail-mode summary. ΔE = E<sub>target</sub> − E<sub>off-target</sub>; more negative values indicate a stronger predicted preference for the target receptor.</span></div><div class="data-table evidence-table wide-table">${table(['Input seed ZINC','Selectivity direction','Hotspot','Fast E<sub>target</sub> (kcal/mol)','Fast E<sub>off-target</sub> (kcal/mol)','Fast ΔE (kcal/mol)','Detail E<sub>target</sub> median (kcal/mol)','Detail E<sub>off-target</sub> median (kcal/mol)','Detail ΔE median (kcal/mol)','Worst ΔE (kcal/mol)','ΔE SD (kcal/mol)','Pose stability','Generated compounds'],rows,'No input seed records are linked to this receptor pair in the 904-compound collection.')}</div>`;
}

function detailModeCompoundPanel(pair,compounds){
  const rows=compounds.map(compound=>`<tr>
    <td><span class="seed-chip"><small>DETAIL-MODE COMPOUND</small><b class="mono">${esc(compound.zinc_id)}</b></span></td>
    <td>${esc(compound.target_name||compound.target_uniprot)} → ${esc(compound.offtarget_name||compound.offtarget_uniprot)}<br><small class="mono">${esc(compound.detail_record_id)}</small></td>
    <td>${badge(compound.hotspot_bw,'good')}</td>
    <td>${fmt(compound.fast.rank,0)}</td>
    <td>${fmt(compound.fast.target,3)}</td>
    <td>${fmt(compound.fast.offtarget,3)}</td>
    <td><b>${fmt(compound.fast.dd,3)}</b></td>
    <td>${fmt(compound.detail.target_median,3)}</td>
    <td>${fmt(compound.detail.offtarget_median,3)}</td>
    <td><b>${fmt(compound.detail.dd_median,3)}</b></td>
    <td>${fmt(compound.detail.dd_best,3)}</td>
    <td>${fmt(compound.detail.dd_worst,3)}</td>
    <td>${fmt(compound.detail.dd_sd,3)}</td>
    <td>${boolBadge(compound.high_confidence,'High confidence','Selected with limitations')}</td>
  </tr>`);
  return `<div class="table-explainer"><b>Detail-mode selected compounds</b><span>Frozen seed compounds retained after fast-mode screening and three-repeat detail-mode evaluation for receptor pairs with no compounds in the 904 Pocketxmol collection. These records are not Pocketxmol-generated compounds. ΔE = E<sub>target</sub> − E<sub>off-target</sub>; more negative values indicate a stronger predicted preference for the target receptor.</span></div><div class="data-table evidence-table wide-table">${table(['ZINC compound','Selectivity direction','Hotspot','Fast rank','Fast E<sub>target</sub> (kcal/mol)','Fast E<sub>off-target</sub> (kcal/mol)','Fast ΔE (kcal/mol)','Detail E<sub>target</sub> median (kcal/mol)','Detail E<sub>off-target</sub> median (kcal/mol)','Detail ΔE median (kcal/mol)','Best ΔE (kcal/mol)','Worst ΔE (kcal/mol)','ΔE SD (kcal/mol)','Evidence status'],rows,'No detail-mode selected compounds are linked to this receptor pair.')}</div>`;
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
      <td><a class="structure-link" href="${esc(structure.bundle_url)}" download title="Download the ${compound.compound_id} structure bundle">${structureText} ↓</a></td>
    </tr>`;
  });
  return `<div class="table-explainer"><b>Pocketxmol-generated compounds</b><span>All records assigned to this receptor pair within the 904-compound collection. ΔE = E<sub>target</sub> − E<sub>off-target</sub>; more negative values indicate a stronger predicted preference for the target receptor. Downloads include a ligand SDF and computational complex PDB files when available.</span></div><div class="data-table evidence-table wide-table molecule-table">${table(['Generated compound / SMILES','Selectivity direction','Input seed ZINC','Similarity to seed','MW','cLogP','QED','SA','Detail ΔE (kcal/mol)','Worst ΔE (kcal/mol)','ΔE change vs seed (kcal/mol)','Structure download'],rows,'No Pocketxmol-generated compounds are linked to this receptor pair in the 904-compound collection.')}</div>`;
}

function showPair(pair){
  const seeds=pairSeeds(pair.pair_id);
  const compounds=pairMolecules(pair.pair_id);
  const detailModeCompounds=pairDetailModeCompounds(pair.pair_id);
  showDialog(`
    <div class="detail-head pair-detail-head">
      <p class="eyebrow">MODULE 02 · RECEPTOR PAIR · RANK ${pair.rank}</p>
      <h2>${esc(pair.receptor_a.name)} / ${esc(pair.receptor_b.name)}</h2>
      <p class="mono">${pair.pair_id}</p>
    </div>
    <div class="pair-summary-strip">
      <span><small>dMaSIF/MaSIF distance</small><b>${fmt(pair.surface_distance,3)}</b></span>
      <span><small>Top differential hotspot</small><b>${esc(pair.hotspots[0]?.bw)}</b></span>
      <span><small>Input seed records</small><b>${seeds.length}</b></span>
      <span><small>Pocketxmol-generated compounds</small><b>${compounds.length}</b></span>
      ${compounds.length===0?`<span><small>Detail-mode selected compounds</small><b>${detailModeCompounds.length}</b></span>`:''}
    </div>
    <div class="evidence-tabs" role="tablist" aria-label="Receptor-pair evidence tables">
      <button class="evidence-tab active" role="tab" aria-selected="true" data-panel-target="hotspots">Top 3 hotspots <b>${pair.hotspots.length}</b></button>
      <button class="evidence-tab" role="tab" aria-selected="false" data-panel-target="seeds">Input seeds <b>${seeds.length}</b></button>
      <button class="evidence-tab" role="tab" aria-selected="false" data-panel-target="molecules">Pocketxmol compounds <b>${compounds.length}</b></button>
      ${compounds.length===0?`<button class="evidence-tab" role="tab" aria-selected="false" data-panel-target="detail-mode">Detail-mode compounds <b>${detailModeCompounds.length}</b></button>`:''}
    </div>
    <section class="evidence-panel active" data-panel="hotspots">${hotspotPanel(pair)}</section>
    <section class="evidence-panel" data-panel="seeds" hidden>${seedPanel(pair,seeds)}</section>
    <section class="evidence-panel" data-panel="molecules" hidden>${compoundPanel(pair,compounds)}</section>
    ${compounds.length===0?`<section class="evidence-panel" data-panel="detail-mode" hidden>${detailModeCompoundPanel(pair,detailModeCompounds)}</section>`:''}
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
    if(!element.innerHTML)element.innerHTML=`<div class="load-error"><b>Data could not be loaded</b><span>${esc(error.message)}. Open the published site or serve the files through a local HTTP server.</span></div>`;
  });
  $('#pair-total-check').textContent='Data have not loaded; the 904 compound-to-pair relationships could not be validated.';
  console.error(error);
});
