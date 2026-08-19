# PyMOL 抗体结构分析工具集

本仓库目前包含两个互补工具：

1. `pymol_antibody_surface_hydrophobic.py`：分析抗体表面 RSA、侧链 SASA 和疏水候选斑块；
2. `pymol_ab_interface_analyzer.py`：分析抗原–抗体界面残基、BSA 及几何支持的相互作用候选。

抗原–抗体界面分析器的完整中文文档见：

- [PyMOL 抗原–抗体结合界面分析器](docs/ANTIBODY_ANTIGEN_INTERFACE_ANALYZER_CN.md)

下面首先是表面疏水分析器的使用说明。

## 表面疏水分析器

一个面向抗体、纳米抗体和其他蛋白结构的 PyMOL 脚本，用于：

- 按残基计算溶剂可接触表面积（SASA）
- 根据 Tien et al. 2013 的 MaxASA 标准计算相对溶剂可及性（RSA）
- 筛选暴露在溶剂表面的疏水氨基酸
- 计算疏水侧链 SASA
- 统计邻近表面疏水残基，辅助识别疏水斑块
- 在 PyMOL 中建立选择集、着色并添加残基标签
- 导出 CSV 报告，便于抗体可开发性分析和突变位点筛选

> 该工具用于生成候选位点，不能单独决定最终突变方案。抗体改造还需要结合 CDR、抗原界面、VH–VL 界面、保守性、结构稳定性和实验数据。

## 文件

```text
pymol_antibody_surface_hydrophobic.py
README.md
```

## 1. 基本原理

### SASA

SASA（solvent-accessible surface area）表示残基能够被溶剂探针接触的面积，单位为 Å²。

### RSA

不同氨基酸大小不同，因此脚本使用 RSA 对 SASA 进行归一化：

```text
RSA = 残基 SASA / 该氨基酸的最大参考 SASA
```

脚本采用 Tien et al. 2013 的理论 Gly-X-Gly MaxASA 数值。

### 表面疏水候选判定

默认情况下，一个残基必须同时满足以下条件：

```text
属于 A/V/I/L/M/F/W/Y
并且 RSA >= 0.25
并且侧链 SASA >= 15 Å²
```

才会进入表面疏水候选选择集。

## 2. 安装

### 方法一：当前 PyMOL 会话手动加载

将脚本保存到本地，例如：

```text
C:\Users\Administrator\Desktop\pymol_antibody_surface_hydrophobic.py
```

在 PyMOL 底部命令行输入：

```pml
run C:/Users/Administrator/Desktop/pymol_antibody_surface_hydrophobic.py
```

同一次 PyMOL 会话只需要执行一次 `run`。

> Windows 路径推荐使用正斜杠 `/`。不要把 `run ...` 输入到 PyMOL 文本编辑器的文件名位置。

### 方法二：启动时自动加载

在 Windows 用户目录创建：

```text
C:\Users\Administrator\pymolrc.pml
```

文件中写入：

```pml
run C:/Users/Administrator/Desktop/pymol_antibody_surface_hydrophobic.py
```

重启 PyMOL 后，可直接使用 `ab_hydro` 或 `surface_hydrophobics`。

也可以在 PyMOL 中使用：

```text
File -> Edit pymolrc
```

加入相同的 `run` 命令。

## 3. 快速开始

加载结构：

```pml
load antibody.pdb, antibody
```

分析默认 H、L 链：

```pml
ab_hydro antibody
```

如果当前只启用了一个结构对象，也可以直接输入：

```pml
ab_hydro
```

如果抗体链名是 A 和 B：

```pml
ab_hydro antibody, chains=A+B
```

## 4. 完整命令

```pml
surface_hydrophobics antibody and chain H+L, context=antibody, rsa_cutoff=0.25, sidechain_sasa_cutoff=15, prefix=ab_hydro
```

### 抗体–抗原复合物

假设完整复合物对象为 `complex`，抗体是 H、L 链：

```pml
surface_hydrophobics complex and chain H+L, context=complex, rsa_cutoff=0.25, sidechain_sasa_cutoff=15, prefix=ab_hydro
```

这里：

- 第一个选择只决定报告哪些残基
- `context=complex` 决定 SASA 计算时保留哪些结构
- 抗原遮挡的抗体界面不会被误判为自由溶剂暴露面

如果需要让糖链或配体参与表面遮挡计算：

```pml
surface_hydrophobics complex and chain H+L, context=complex, include_hetatm=1, prefix=ab_hydro
```

水分子始终会被排除。

## 5. 参数

| 参数 | 含义 | 默认值 |
|---|---|---:|
| `selection` | 需要分析并写入报告的蛋白残基 | `polymer.protein` |
| `context` | SASA 计算的完整结构环境；为空时等于 selection | 空 |
| `rsa_cutoff` | 表面暴露 RSA 阈值 | 0.25 |
| `sidechain_sasa_cutoff` | 疏水侧链最小暴露面积，Å² | 15.0 |
| `hydrophobic` | 疏水残基单字母集合 | `AVILMFWY` |
| `patch_distance` | 统计疏水邻居的侧链中心距离，Å | 7.5 |
| `dot_density` | PyMOL 表面采样密度，1–4 | 3 |
| `state` | 分析的坐标状态 | 1 |
| `prefix` | 输出对象、选择集和默认 CSV 的前缀 | `hydro` |
| `csv_file` | CSV 完整输出路径；为空则写入当前工作目录 | 空 |
| `include_hetatm` | 是否保留非蛋白原子作为遮挡环境 | 0 |
| `quiet` | 是否减少候选残基终端输出 | 0 |

如需把 Cys 和 Pro 也纳入疏水集合：

```pml
surface_hydrophobics antibody, context=antibody, hydrophobic=AVILMFWYCP, prefix=ab_hydro
```

## 6. 输出结果

运行后会生成以下 PyMOL 对象或选择集。

### `<prefix>_view`

用于显示的结构副本。原始对象不会被修改。

脚本将每个残基的：

```text
RSA × 100
```

写入该副本的 B-factor 字段。

### `<prefix>_all_exposed`

所有达到 RSA 阈值的残基，包括亲水、带电和疏水残基。

### `<prefix>_exposed_hydrophobic`

最重要的选择集：同时通过疏水类型、RSA 和侧链 SASA 筛选的候选残基。

### CSV 报告

默认文件名：

```text
<prefix>_rsa.csv
```

例如：

```text
ab_hydro_rsa.csv
```

默认写入 PyMOL 当前工作目录。查看当前目录：

```pml
pwd
```

明确保存到桌面：

```pml
surface_hydrophobics antibody, context=antibody, prefix=ab_hydro, csv_file=C:/Users/Administrator/Desktop/ab_hydro_rsa.csv
```

## 7. 颜色解释

脚本使用两套显示逻辑。

### 蓝–白–红：RSA 暴露程度

```text
蓝色：RSA 低，残基大部分埋藏
白色：中等暴露
红色：RSA 高，高度暴露
```

这些颜色不表示疏水强弱。

即使蓝色区域出现在可见分子表面，也只表示对应残基整体暴露比例较低；该残基仍可能有少量原子贡献可见表面。

### 橙色：表面疏水候选

橙色表示该残基通过了以下筛选：

```text
疏水残基类型 + RSA 阈值 + 侧链 SASA 阈值
```

所有候选统一为橙色，橙色深浅不表示疏水强度。

如果希望显示更直观，可将普通表面改为灰色：

```pml
color gray80, ab_hydro_view
color orange, ab_hydro_exposed_hydrophobic
show surface, ab_hydro_view
show sticks, ab_hydro_exposed_hydrophobic
```

## 8. CSV 字段

| 字段 | 含义 |
|---|---|
| `chain` | 链名 |
| `segi` | segment 标识 |
| `resi` | 残基编号，保留插入码 |
| `resn` | 三字母氨基酸名称 |
| `aa` | 单字母氨基酸名称 |
| `sasa_A2` | 残基总 SASA |
| `rsa` | 相对溶剂可及性 |
| `sidechain_sasa_A2` | 侧链 SASA |
| `exposed` | 是否达到 RSA 阈值 |
| `hydrophobic_candidate` | 是否为表面疏水候选 |
| `hydrophobic_neighbors_within_7.5A` | 7.5 Å 内其他候选疏水残基数量 |

建议首先筛选：

```text
hydrophobic_candidate = 1
```

然后优先检查：

- RSA 较大
- 侧链 SASA 较大
- 疏水邻居数较多
- 多个候选在空间上形成连续斑块

## 9. 查看候选残基

只显示候选：

```pml
hide sticks, all
show sticks, ab_hydro_exposed_hydrophobic
color orange, ab_hydro_exposed_hydrophobic
zoom ab_hydro_exposed_hydrophobic
```

查看点击残基的 RSA：

```pml
iterate (byres pk1) and name CA, print(chain, resi, resn, b/100.0)
```

最后一列为 RSA。

## 10. 调整筛选严格程度

默认推荐：

```text
RSA >= 0.25
侧链 SASA >= 15 Å²
```

更严格：

```pml
surface_hydrophobics antibody, context=antibody, rsa_cutoff=0.35, sidechain_sasa_cutoff=25, prefix=ab_strict
```

更宽松：

```pml
surface_hydrophobics antibody, context=antibody, rsa_cutoff=0.20, sidechain_sasa_cutoff=10, prefix=ab_relaxed
```

建议同时比较默认、严格和宽松结果，不要把单一阈值当作绝对分类。

## 11. 抗体改造建议

优先检查：

1. 多个候选聚集形成的疏水斑块
2. RSA 和侧链 SASA 同时较高的残基
3. 位于框架区、远离抗原界面的候选
4. 暴露的 Phe、Trp、Leu、Ile、Val、Met 等残基

突变前应排除或谨慎处理：

- CDR 关键残基
- 抗原直接接触残基
- VH–VL 界面残基
- 维持疏水核心的残基
- 二硫键相关 Cys
- 关键保守位点
- 参与局部氢键或芳香堆积的残基

候选位点还应结合序列保守性、结构稳定性预测、静电势、聚集实验和表达数据评估。

## 12. 常见问题

### `cmd.label() got an unexpected keyword argument 'space'`

旧版 PyMOL 的 `cmd.label` 不支持 `space=` 参数。本仓库版本不使用该参数，兼容性更好。

### `SyntaxError: expected 'except' or 'finally' block`

通常是手动修改标签代码时破坏了 `try` 内部缩进。建议重新下载仓库中的完整脚本，不要逐行修补。

### `UnicodeDecodeError: 'gbk' codec...`

这是部分 Windows PyMOL 内置文本编辑器读取 UTF-8 文件时的问题。不要使用内置编辑器打开脚本；直接在底部命令行使用 `run` 加载。

### `ab_hydro` 不是已知命令

说明脚本尚未加载。先执行：

```pml
run C:/路径/pymol_antibody_surface_hydrophobic.py
```

### 当前启用了多个对象

`ab_hydro` 无参数时只支持自动识别一个启用对象。请明确指定：

```pml
ab_hydro antibody
```

### 找不到 H、L 链

检查实际链名并指定：

```pml
ab_hydro antibody, chains=A+B
```

### 找不到 CSV

执行：

```pml
pwd
```

或者显式指定 `csv_file`。只有脚本完整运行结束后才会写出 CSV。

## 13. 方法限制

- PyMOL `get_area` 使用表面点采样，数值是近似结果
- 结果依赖结构完整性、缺失侧链、构象状态和复合物装配方式
- 单一静态结构不能描述抗体全部构象波动
- 蛋白表面疏水风险不等同于单残基疏水性
- 如需发表级定量结果，建议使用 FreeSASA、DSSP 或其他独立方法交叉验证

## 14. 参考资料

- Tien MZ, Meyer AG, Sydykova DK, Spielman SJ, Wilke CO. Maximum Allowed Solvent Accessibilities of Residues in Proteins. PLoS ONE. 2013;8(11):e80635. https://doi.org/10.1371/journal.pone.0080635
- PyMOL `get_area` command reference: https://pymol.org/pymol-command-ref.html
