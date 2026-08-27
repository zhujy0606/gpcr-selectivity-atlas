# 公开发布方案

## 平台决定

- 网页与小型CSV/JSON：GitHub Pages；
- 约1.1 GB dMaSIF原始资产：Zenodo；
- 可选高速镜像：Cloudflare R2，不作为学术主档案。

网站仅约4.2 MB，适合静态托管。大型资产不进入Git仓库；Zenodo记录负责永久版本、DOI、作者、许可证和推荐引用。网页中的每个原始资产下载链接最终指向Zenodo版本记录。

## GitHub Pages发布

`public_database/`作为独立公开仓库根目录。工作流`.github/workflows/deploy-pages.yml`会在`main`分支更新时发布`site/`。

- 公开仓库：https://github.com/zhujy0606/gpcr-selectivity-atlas
- 公开网站：https://zhujy0606.github.io/gpcr-selectivity-atlas/

发布步骤：

1. [x] 创建公开GitHub仓库`zhujy0606/gpcr-selectivity-atlas`；
2. [x] 将`public_database/`内容作为仓库根目录推送到`main`；
3. [x] 将Pages Source设为GitHub Actions；
4. [x] 完成Deploy GPCR Selectivity Atlas工作流；
5. [x] 将生成的网址写回README；
6. [ ] Zenodo发布后回填DOI和大型资产链接。

## Zenodo发布

建议一个Zenodo Dataset记录，文件控制在三个以内：

1. `gpcr_selectivity_atlas_v0.1.0_public_tables.zip`；
2. `gpcr_selectivity_atlas_v0.1.0_dmasif_assets.zip`；
3. `SHA256SUMS.txt`。

发布前必须确认：

- 正式标题、作者顺序、单位和ORCID；
- 数据许可证；
- PocketXMol、ZINC、GPCRdb、PDB、UniProt和ADMET-AI衍生数据的再分发边界；
- 版本号与发布日期；
- 说明所有结果属于计算优先级排序证据。

Zenodo DOI产生后，回填网站下载链接、`CITATION.cff`和README。

## 更新策略

- GitHub仓库：网页代码与小型规范数据的持续版本；
- Zenodo：每次正式数据冻结创建新版本；
- DOI引用固定版本，网页默认指向最新版本；
- 每次构建保留源文件SHA-256和导出测试结果。
