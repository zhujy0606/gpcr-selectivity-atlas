# 数据字典

## 稳定标识

- `uniprot`：受体的UniProt accession。
- `pair_id`：两个UniProt排序后以`__`连接，表示无向受体对。
- `task`：目标/脱靶方向明确的Pocketxmol任务。
- `seed_record_id`：`task::ZINC_ID`，同一亲本化学分子在不同方向任务中保留独立记录。
- `compound_id`：冻结的A6PX预筛ID。
- `hotspot_id`：`pair_id::H1`至`H3`。

## 关键语义

- `surface_distance`：dMaSIF表面指纹距离。不是序列相似度、活性态结构RMSD或实验亲和力。
- `fingerprint_difference`：局部BW位点的表面指纹差异强度。
- `detail_dd_median`：网页公开标签为 `Detail ΔE`，定义为三次 detail-mode 中 `E_target − E_off-target` 的中位数，单位为 kcal/mol；越负表示计算上越偏向目标受体。
- `detail_dd_worst`：网页公开标签为 `Worst ΔE`，表示三次重复中最不利的 ΔE，单位为 kcal/mol，用于稳健性评估。
- `detail_dd_sd`：三次重复 ΔE 的标准差，单位为 kcal/mol。
- `dd_change_vs_seed`：生成分子相对于输入种子的 ΔE 变化，单位为 kcal/mol。内部字段名为兼容既有冻结数据而保留；公开网页统一使用 ΔE 表述。
- `final_candidate_904`：属于网站公开展示的904个Pocketxmol生成分子；这些分子均满足采样稳健且 detail-mode ΔE 优于对应输入种子。
- `seed_zinc_id`：生成该分子时输入Pocketxmol的亲本种子ZINC号。
- `structure_download.bundle_url`：以A6PX分子ID命名的公开ZIP结构包。所有904个包都含1个生成分子配体SDF；在有计算复合物证据时，包内还包含目标/脱靶复合物PDB。
- `structure_download.complex_pdb_count`：该分子结构包内的计算复合物PDB数量；0表示该结构包没有可公开的计算复合物PDB，不代表配体结构文件丢失。
- `structure_download.complex_pdbs`：ZIP内每个复合物PDB的受体角色（target/offtarget）、UniProt、姿势簇、文件名、大小及SHA-256。

## 结构文件

结构包内的SDF和PDB均来自冻结计算证据并经过哈希核验。PDB是计算受体–配体复合物，不是PDB数据库中的实验解析结构。904个Pocketxmol生成分子均有配体SDF；在具有计算复合物证据时，结构包同时提供对应PDB。

## 空值

`null`表示未计算、不可用或无法从冻结来源可靠连接。它不表示数值为零，也不表示阴性结果。

## MM-GBSA终点重评分（2026-11-12网页更新）

`compounds.json`中的`mmgbsa`关联全部904个候选物；`seeds.json`中的`mmgbsa`关联178条母体任务。仅以ZINC号连接会混淆不同任务和方向，因此使用`seed_record_id = task::ZINC_ID`并核对目标/非目标受体。父记录`record_id`对应候选物的`parent_id`。

| 字段 | 含义 |
| --- | --- |
| `status` | 候选物为`paired_complete`；母体为`complete` |
| `target_energy_kcal_mol` | 目标受体代表构象的MM-GBSA结合能估计 |
| `offtarget_energy_median_kcal_mol` | 非目标受体1–3个代表构象的结合能估计的中位数 |
| `differential_kcal_mol` | D = 目标能量 − 非目标构象能量中位数 |
| `parent_differential_kcal_mol` | 实际母体种子的D |
| `delta_vs_parent_kcal_mol` | δ = D候选物 − D母体，负值为更有利 |
| `improved_vs_parent` | δ<0；447条为true，457条为false，不代表实验活性或选择性 |
| `target_pose_count` / `offtarget_pose_count` | 分别用于该分子目标/非目标聚合的构象数 |
| `target_pose_stable` / `offtarget_pose_stable` | 该分子原对接姿态稳定性标记，非MD稳定性；未据此排除记录 |
| `parent_target_pose_stable` / `parent_offtarget_pose_stable` | 对应实际母体的原对接姿态稳定性 |

下载：`mmgbsa_candidates_904.csv`、`mmgbsa_parents_178.csv`和`mmgbsa_pose_endpoints_3733.csv`（UTF-8 BOM，便于Windows Excel打开）。构象表包含分子、任务、方向、受体、簇序号及结合能估计，可重建D和δ；不公开原始服务器路径。常规`compounds_904.csv`和`seeds.csv`同步加入`mmgbsa_`前缀字段。协议、输入哈希及范围见`site/data/mmgbsa_protocol.json`。

采用ff14SB/GAFF2/AM1-BCC、两阶段约束最小化和GB（igb=5、0.150 M）。每个代表构象仅使用一个最小化快照，无熵校正和膜MD采样；打分模型不同，但初始结构来自既有对接。该网站更新未收录于固定Zenodo v0.2.0 DOI。原有结构ZIP未改变，不能据其文件数推断本次数值重评分覆盖数。

## 904与438的集合关系

- `pareto.shortlisted`：是否属于原冻结`cross_domain_shortlist.csv`中的438个主要计算优选候选物。全部属于904个候选集合；网站以淡紫色行及“Final Pareto shortlist”标签标注。
- `pareto.target_pose_eligible`：是否通过靶受体姿态收敛门，共749个。
- `pareto.seed_relative_eligible`：是否满足原母体相对筛选条件。
- `pareto.docking_front`、`pareto.druglike_front`、`pareto.admet_front`：分别属于对接、类药性、ADMET领域的Pareto前沿；前沿在相同任务和姿态证据层级内计算。
- `pareto.domain_front_count`：所属前沿数；最终438个中，225个为2，213个为3。
- `pair.pareto_shortlist_count`：该无向受体对中的优选候选数；总计438，分布于108个无向受体对。原Pareto汇总的129个pair值采用有向/原始pair键，不用于本字段。
- 下载`pareto_shortlist_438.csv`包括438个优选项的分子、对接、MM-GBSA及领域成员标记；完整904表保留相应`pareto_`前缀字段。原`final_candidate_904`和`final_candidates`兼容字段继续指完整904集合，主要交付优选集合以`pareto.shortlisted`及`counts.pareto_shortlist`为准。

高亮不代表MM-GBSA改善。438个Pareto优选集合与447个MM-GBSA改善集合按不同规则定义，均不能替代实验活性或选择性证据。集合定义、原规则及来源哈希见`site/data/cohort_definitions.json`。

## 预测性质

网页展示MW、cLogP、TPSA、HBD、HBA、可旋转键、QED、SA及AMES、BBB、口服生物利用度、主要CYP抑制、ClinTox、DILI、HIA、PAMPA、P-gp、hERG、Caco-2、清除率、半衰期、溶解度和VDss。所有ADMET均为模型预测。
