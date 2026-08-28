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
- `detail_dd_median`：三次detail-mode中，目标得分减脱靶得分的中位差；越负表示计算上的目标选择性越有利。
- `detail_dd_worst`：三次重复中最不利的DD，用于稳健性评估。
- `final_candidate_904`：属于网站公开展示的904个Pocketxmol生成分子；这些分子均满足采样稳健且detail-mode对接DD优于对应输入种子。
- `seed_zinc_id`：生成该分子时输入Pocketxmol的亲本种子ZINC号。
- `structure_download.bundle_url`：以A6PX分子ID命名的公开ZIP结构包。所有904个包都含1个生成分子配体SDF；在有计算复合物证据时，包内还包含目标/脱靶复合物PDB。
- `structure_download.complex_pdb_count`：该分子结构包内的计算复合物PDB数量；0表示该结构包没有可公开的计算复合物PDB，不代表配体结构文件丢失。
- `structure_download.complex_pdbs`：ZIP内每个复合物PDB的受体角色（target/offtarget）、UniProt、姿势簇、文件名、大小及SHA-256。

## 结构文件

结构包内的SDF和PDB均来自冻结计算证据并经过哈希核验。PDB是计算受体–配体复合物，不是PDB数据库中的实验解析结构。904个Pocketxmol生成分子均有配体SDF；在具有计算复合物证据时，结构包同时提供对应PDB。

## 空值

`null`表示未计算、不可用或无法从冻结来源可靠连接。它不表示数值为零，也不表示阴性结果。

## 预测性质

网页展示MW、cLogP、TPSA、HBD、HBA、可旋转键、QED、SA及AMES、BBB、口服生物利用度、主要CYP抑制、ClinTox、DILI、HIA、PAMPA、P-gp、hERG、Caco-2、清除率、半衰期、溶解度和VDss。所有ADMET均为模型预测。
