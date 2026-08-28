# GPCR Selectivity Atlas

这是一个由项目冻结证据自动生成的静态公开数据库。公众页面收敛为两个表格化模块：

1. **287受体 dMaSIF 数据**：逐受体汇总表面坐标、embedding、VTK资产、完整性、最近表面邻居和全局距离，并可展开文件级清单。
2. **163对受体与选择性分子**：逐受体对汇总dMaSIF/MaSIF距离；点击任一受体对，可在三张关联表中查看Top 3差异热点、Pocketxmol输入种子和该受体对对应的Pocketxmol生成分子。

第二模块完整纳入904个Pocketxmol生成分子及其输入种子ZINC号、对接与理化性质。每个分子均提供配体SDF结构包；有复合物结构证据的分子还提供目标/脱靶计算复合物PDB。

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

当前站点已经通过GitHub Pages公开。受体对的Pocketxmol生成分子表中，结构下载包按A6PX ID组织，ZIP内含配体SDF以及可用的计算复合物PDB；每个包和内部结构文件均记录SHA-256。约1.1 GB的dMaSIF二进制资产当前只发布元数据清单；应在确认再分发条款后上传Zenodo，再回填公开URL和SHA-256。
