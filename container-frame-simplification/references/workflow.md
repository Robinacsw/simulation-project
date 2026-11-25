# 运行接口与数据约定

## 输入与工具

工作目录由本次任务指定，每个命令写到独立的新子目录。使用可导入 OCP、NumPy、SciPy 的 Python 3.12。截图还使用 Pillow、Numba 和 Windows 微软雅黑字体。报告或指令中的模型位置、版本名、配置和实体号均从此次输入读取。

直接读取带层级的毫米制 STEP/STP。原生 SolidWorks 装配以任务指定版本打开，核对配置、引用、抑制和隐藏后导出同配置 STEP，再进入 `scan`。已有原生清单时，以独立叶实例数核对并传给 `--expected-leaves`。格式桥接后保留源文件和配置的关联。SpaceClaim 几何写出及回读使用已验证的 2025 R2 API V252，`--spaceclaim` 传其可执行文件路径。当前工程以已有 X_T 为基线时，`deliver --geometry` 指定该版 X_T。

`scan` 自动建立箱体局部坐标。如果原始几何无法确定方向，测量并传 `--datum`，文件包含毫米制 `origin` 和三个按长、宽、高排列的单位向量 `axes`。未知单位、缩放或镜像变换、缺几何和不可确认的关键接口分别记录，补齐后再裁定。遇到 CAD 软件不可用时可以交付已经验收的识别表和待复核件。

## 命令示例

下面变量都由当前任务填写。`$reviewProfile` 是确认前工况，`$approvedProfile` 只更新已确认的工况许可。`$stepPath` 和 `$geometryPath` 应有可追溯的相同装配坐标。`$caseDir` 是本次工作目录。

```powershell
& $cadPython -B $entry scan --step $stepPath --configuration $configuration --profile $reviewProfile --out "$caseDir/recognition"
& $cadPython -B $entry review --analysis "$caseDir/recognition" --rule ROOF_HOLDER --out "$caseDir/review"
& $cadPython -B $entry replay --analysis "$caseDir/recognition" --profile $approvedProfile --out "$caseDir/approved"
& $cadPython -B $entry deliver --analysis "$caseDir/approved" --geometry $geometryPath --approval $approvalFile --protection $protectionFile --spaceclaim $spaceClaimExe --name $newXtName --out "$caseDir/delivery"
```

`$entry` 指向当前 Skill 的 `scripts/frame_simplify.py`。`scan --mode G` 只跑形状位置，默认 `G+C` 量测关系。已有 G+C 包可用 `replay --mode G` 以相同几何做对照。G 包升为 G+C 时重新运行量测。已知待复核规则可重复传 `--rule`，用 `逐件识别与去留.csv` 中的判据编号定位。截图为源 STEP 的 CAD 曲面渲染，标注全部候选位置及代表性局部。以当前 X_T 为去留依据时先完成对应对账。

## 工况、批准和保护文件

- `profile.json` 按 [工况样例](../assets/profile.example.json) 写 `cases` 和 `scope`。每个工况记录 `name`、`entrance` 和 `support`，其中支持描述包含限位、夹紧或角位组合。`fixture_boundary` 可取 `discrete_castings_and_floor`、`other` 或 `unconfirmed`。底蒙板旁垫板、顶盖连接板、下方固定座的专项布尔许可见样例，按任务确认后设置。
- `approval.json` 记录当前几何的 `source_sha256`、工程依据 `decision`、本次批准的 `rules`。若只批准某些实例，另填同一封存结果的 `instances`。范围和工况来自项目资料与用户确认。
- `protection.json` 可含 `active` 和 `other` 两组实测实体签名，用于复核已知活动梁连接和其他关键几何。每条签名包含 `tight_box`、`center`、`volume`、`faces`，使用源装配世界坐标和毫米单位。可从可信的 SpaceClaim 几何回读按项目确认范围抽取。批量清理已有工程时应提供项目已掌握的保护件签名。

## 输出与回读

- 匿名结果、输入与程序指纹、接口证据和中性编号由运行目录记录。`trace_private.json` 提供原始标签回查，`逐件识别与去留.csv` 汇集判据、建议动作、证据编号和待核条件。
- `review` 输出整箱位置、局部几何、可用邻件及实测间隙，`截图覆盖实例.csv` 记录每组覆盖的实例。候选几何缺失时在 `截图索引.json` 中记原因。
- `deliver` 先复制并独立读取当前几何。按紧包围盒、质心、体积、曲面面积和面数建立唯一对应。本批编号来自实际导入顺序，重新写出前以指纹及样本量测核对。位姿不同、同形重合、候选争用、存在未对应的曲面体时，在 `清理前后对账.csv` 记录待复核。另存 X_T 后重开并逐体配对，检查体数量、位置、活动梁保护及关键几何签名。`X_T几何验收.json` 保存结果。
- 现有匹配数值容差为包围盒及质心 0.01 mm、体积绝对阈值 0.005 mm³ 或相对阈值 0.005%。STEP 到 X_T 的对应允许 1 个面的分面差，同时要求两边实测表面积满足绝对阈值 0.005 mm² 或相对阈值 0.005%。同一 X_T 导出前后的留件验收要求面数一致。几何容差用于已经核实的对应与回读，工程焊缝按实际界面另行定义。
