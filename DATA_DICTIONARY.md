# 数据字典

## 稳定标识

- `uniprot`：受体的UniProt accession。
- `pair_id`：两个UniProt排序后以`__`连接，表示无向受体对。
- `task`：目标/脱靶方向明确的PocketXMol任务。
- `seed_record_id`：`task::ZINC_ID`，同一亲本化学分子在不同方向任务中保留独立记录。
- `compound_id`：冻结的A6PX预筛ID。
- `hotspot_id`：`pair_id::H1`至`H3`。

## 关键语义

- `surface_distance`：dMaSIF表面指纹距离。不是序列相似度、活性态结构RMSD或实验亲和力。
- `fingerprint_difference`：局部BW位点的表面指纹差异强度。
- `detail_dd_median`：三次detail-mode中，目标得分减脱靶得分的中位差；越负表示计算上的目标选择性越有利。
- `detail_dd_worst`：三次重复中最不利的DD，用于稳健性评估。
- `mmgbsa_baseline_complete`：候选及其正确亲本种子均有完整、有效、可配对的MM/GBSA端点。
- `mmgbsa_dual_endpoint_improved`：median MM/GBSA和worst-pose MM/GBSA均优于配对种子。
- `final_selected`：进入冻结的最终111个候选。

## 空值

`null`表示未计算、不可用或无法从冻结来源可靠连接。它不表示数值为零，也不表示阴性结果。

## 预测性质

网页展示MW、cLogP、TPSA、HBD、HBA、可旋转键、QED、SA及AMES、BBB、口服生物利用度、主要CYP抑制、ClinTox、DILI、HIA、PAMPA、P-gp、hERG、Caco-2、清除率、半衰期、溶解度和VDss。所有ADMET均为模型预测。
