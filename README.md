# GPCR Selectivity Atlas

这是一个由项目冻结证据自动生成的静态公开数据库。首版整合287个受体的dMaSIF资产清单、286×286全局表面距离矩阵、163组受体对及489个局部热点、178条“有向任务×亲本种子”记录、138个唯一ZINC亲本种子、904个PocketXMol化合物及其对接、理化、ADMET、姿势与MM/GBSA状态。

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

当前站点是可部署MVP，尚未自动上传到公开互联网。约1.1 GB的dMaSIF二进制资产只发布元数据清单；应在确认再分发条款后上传Zenodo或对象存储，再回填公开URL和SHA-256。活性结构页面目前只有36条结构参考，系统性活性态结构相似度尚未计算。
