# GPCR Selectivity Atlas

这是一个由项目冻结证据自动生成的静态公开数据库。首版整合287个受体的dMaSIF资产清单、286×286全局表面距离矩阵、163组受体对及489个局部热点、178条“有向任务×PocketXMol输入种子”记录、138个唯一ZINC输入种子，以及904个PocketXMol最终候选及其对接、理化、ADMET、姿势与MM/GBSA状态。111个分子作为更严格的终态精选子集单独标注。每个候选均提供配体SDF结构包；进入MM/GBSA的438个候选还提供合计1,444个目标/脱靶计算复合物PDB。

- 公开网站：https://zhujy0606.github.io/gpcr-selectivity-atlas/
- GitHub仓库：https://github.com/zhujy0606/gpcr-selectivity-atlas
- 作者：朱景一（山东大学），ORCID：[0009-0003-8404-0455](https://orcid.org/0009-0003-8404-0455)

## 构建

在项目根目录运行：

```bash
./.venv/bin/python public_database/build.py
./.venv/bin/python -m unittest public_database/tests/test_public_build.py
```

## 本地浏览

```bash
cd public_database/site
python3 -m http.server 8080
```

然后访问`http://localhost:8080/`。该站点不依赖后端服务或第三方JavaScript库，可直接部署到GitHub Pages、Cloudflare Pages或任何静态网站服务器。

## 推荐公开架构

本项目采用GitHub Pages托管网页，Zenodo托管约1.1 GB dMaSIF原始资产并提供DOI。发布步骤见`PUBLISHING.md`，GitHub Pages工作流已配置在`.github/workflows/deploy-pages.yml`。

## 发布边界

当前站点已经通过GitHub Pages公开。候选详情中的结构下载包按A6PX ID组织，ZIP内含配体SDF以及可用的计算复合物PDB；每个包和内部结构文件均记录SHA-256。约1.1 GB的dMaSIF二进制资产只发布元数据清单；应在确认再分发条款后上传Zenodo，再回填公开URL和SHA-256。活性结构页面目前只有36条结构参考，系统性活性态结构相似度尚未计算。
