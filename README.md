# 仿真项目

用于集装箱与储能箱静力分析的 Codex Skill。按当前装配配置和试验工况识别框架，整理可清理件和待复核接口，并验收另存的 Parasolid X_T。

## 集装箱框架简化

通过 [Skill 入口](container-frame-simplification/SKILL.md) 了解输入和流程。保留目录结构安装后，可使用 `$container-frame-simplification` 调用。

- 读取同配置 STEP，以几何、装配位置和连接关系做匿名识别。
- 分别记录零件用途、清理建议和工程批准。输出有歧义位置的分组截图。
- 逐体匹配批准清理件，保护活动梁连接，并独立重开新版 X_T 核对留件。

依赖、命令和模板见 [运行接口](container-frame-simplification/references/workflow.md)。包内包括通用脚本、规则、输入示例和快速回归。

## 使用许可

本仓库按 [PolyForm Noncommercial License 1.0.0](LICENSE) 提供。非商业使用、修改及分发适用该许可的具体条件，超出许可范围的商业使用需另行获得授权。分发时保留 [版权声明](NOTICE) 和许可文本。
