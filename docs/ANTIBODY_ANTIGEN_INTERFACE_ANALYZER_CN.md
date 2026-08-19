# PyMOL 抗原–抗体结合界面分析器

文件：`pymol_ab_interface_analyzer.py`

当前版本：0.1.0

## 1. 它解决什么问题

这个脚本在 PyMOL 中接收两个互不重叠的蛋白选择，例如“抗体”和“抗原”，然后同时完成：

1. 识别结合界面残基；
2. 计算复合物埋藏表面积（BSA）和每个残基的 ΔSASA；
3. 识别具有几何支持的相互作用候选；
4. 在 PyMOL 中建立选择、棒状结构和虚线对象；
5. 导出相互作用 CSV、界面残基 CSV 和完整 JSON 审计报告。

它是命令式插件：加载一次脚本后，可以像 PyMOL 自带命令一样执行 `ab_interface`。

## 2. 最重要的科学边界

本工具分析的是**静态结构中的几何接触和化学类型候选**，不是实验测得的作用力，也不计算结合自由能。

因此报告使用以下措辞：

- `distance_defined`：由公开距离阈值定义；
- `geometry_supported`：距离和方向都符合几何规则；
- `putative`：候选作用，例如没有显式氢时的氢键；
- `charge_state_assumed`：依赖标准质子化状态假设；
- `geometry_warning`：结构冲突警告。

不能仅凭一个晶体结构和距离阈值断言某个残基“贡献了多少 kcal/mol”。若要判断热点残基或能量贡献，应进一步进行 Rosetta InterfaceAnalyzer、分子动力学、自由能计算或实验突变验证。

## 3. 方法依据

实现综合参考了以下成熟项目：

- [BindCraft](https://github.com/martinpacesa/BindCraft)：界面残基采用任意原子 4.0 Å 接触；其更深层界面评分由 Rosetta InterfaceAnalyzer 完成。本脚本采用前者，不冒充 Rosetta 能量结果。
- [PDBe Arpeggio](https://github.com/PDBeurope/arpeggio)：强调按原子类型分类接触、保留规则和结构质量警告。
- [GetContacts](https://github.com/getcontacts/getcontacts)：盐桥、氢键、π 相互作用等几何阈值的主要参考。
- [PRODIGY](https://github.com/haddocking/prodigy)：蛋白–蛋白接触及埋藏表面积的验证背景。
- [RCSB PDB 1MLC](https://www.rcsb.org/structure/1MLC)：公开 Fab–溶菌酶复合物，用于端到端 PyMOL 集成验证。

## 4. 安装和加载

### 临时加载

每次重新打开 PyMOL 后执行：

```pml
run C:/你的路径/pymol_ab_interface_analyzer.py
```

加载成功后可输入：

```pml
ab_interface_rules
```

如果能打印规则和阈值，说明命令已经注册。

### 自动加载

在 PyMOL 的启动脚本 `pymolrc.py` 中加入：

```python
run C:/你的路径/pymol_ab_interface_analyzer.py
```

以后打开 PyMOL 就可以直接使用 `ab_interface`NaN`complex`，抗体链为 H、L，抗原链为 A：

```pml
ab_interface (complex and chain H+L), (complex and chain A), prefix=case1
```

约定上 partner1 写抗体、partner2 写抗原，但数学计算是对称的。

指定结果目录：

```pml
ab_interface (complex and chain H+L), (complex and chain A), prefix=case1, output_dir=C:/Users/Administrator/Desktop/interface_results
```

关闭残基标签：

```pml
ab_interface (complex and chain H+L), (complex and chain A), prefix=case1, labels=0
```

更改少数常用阈值：

```pml
ab_interface (complex and chain H+L), (complex and chain A), prefix=case1, contact_cutoff=4.0, bsa_cutoff=1.0, hbond_cutoff=3.5, salt_cutoff=4.0, hydrophobic_cutoff=4.5
```

两个选择必须非空、互不重叠，并建议位于同一个复合物对象中且使用不同链 ID。

## 6. 公开示例：1MLC

1MLC 含有两套 Fab–溶菌酶复合物。第一套中：

- 抗体轻链：A；
- 抗体重链：B；
- 溶菌酶抗原：E。

在联网的 PyMOL 中执行：

```pml
fetch 1mlc, async=0
run C:/你的路径/pymol_ab_interface_analyzer.py
ab_interface (1mlc and chain A+B), (1mlc and chain E), prefix=test_1mlc, output_dir=C:/Users/Administrator/Desktop/interface_results
```

不要把 A+B+C+D 全部当成一个抗体再与 E+F 同时计算，因为那是两个晶体学复合物副本，会让结果难以解释。

## 7. 界面残基如何定义

脚本同时报告两种定义，并提供二者并集用于 PyMOL 显示。

### 接触界面

某个残基只要有任意重原子与另一侧任意重原子的距离 ≤4.0 Å，就属于接触界面。这与 BindCraft 的热点/界面残基距离定义一致。

### BSA 界面

脚本计算：

```text
BSA = (SASA_partner1 + SASA_partner2 - SASA_complex) / 2
```

并计算每个残基从孤立 partner 到复合物状态损失的 SASA，即 ΔSASA。默认 ΔSASA ≥1.0 Å² 的残基属于 BSA 界面。

PyMOL SASA 设置为：

- `dot_solvent = 1`
- `dot_density = 4`
- 不计氢原子和溶剂原子

接触界面回答“是否足够靠近”，BSA 界面回答“结合后是否发生可检测的表面埋藏”。两者不完全相同是正常现象。

## 8. 相互作用规则

| 类型 | 默认规则 | 解释限制 |
|---|---|---|
| close_contact | 任意重原子距离 ≤4.0 Å | 仅代表接近，不等同于特定作用力 |
| hydrogen_bond | PyMOL donor/acceptor 类型；供受体 ≤3.5 Å；mode=1 方向检查；angle=63° | 无显式氢时属于候选氢键 |
| salt_bridge | Lys/Arg 与 Asp/Glu 指定带电侧链原子 ≤4.0 Å | 假设标准质子化；默认排除 His 和链末端 |
| hydrophobic_contact | A/V/I/L/M/F/W/Y 的侧链 C/S 原子 ≤4.5 Å | 是疏水接触候选，不是疏水能量 |
| pi_stacking | 芳环质心 ≤7.0 Å；法向夹角 ≤30°；psi ≤45° | 基于环平面几何 |
| t_stacking | 芳环质心 ≤5.0 Å；法向接近 90°±30°；psi ≤45° | 基于环平面几何 |
| cation_pi | Lys/Arg 阳离子中心至芳环质心 ≤6.0 Å；环轴夹角 ≤60° | 依赖质子化状态 |
| disulfide | Cys SG–SG 1.8–2.3 Å | 距离支持；未读取化学键级 |
| vdw_contact | 原子距离与范德华半径和之差绝对值 ≤0.5 Å | 使用 PyMOL 半径，缺失时用元素默认值 |
| steric_clash | 原子距离小于范德华半径和 0.6 Å 以上 | 应检查分辨率、altloc、缺失原子和建模质量 |

芳香环包括 Phe、Tyr、His 和 Trp。Trp 使用吲哚的六元环部分定义稳定平面。

## 9. PyMOL 中会出现什么

以 `prefix=case1` 为例：

- `case1_p1_interface`：partner1 界面残基；
- `case1_p2_interface`：partner2 界面残基；
- `case1_interface`：双方界面并集；
- `case1_hydrogen_bond`、`case1_salt_bridge` 等：各类虚线对象；
- `case1`：包含上述对象的 PyMOL group。

默认配色：

- 氢键：青色；
- 盐桥：洋红；
- 疏水接触：橙色；
- π–π：紫色；
- T 型堆积：紫罗兰；
- 阳离子–π：海蓝；
- 二硫键：黄色；
- 空间冲突：红色。

π 和阳离子–π 的 CSV 距离是质心距离；PyMOL 虚线连接的是用于定位的代表原子，因此虚线本身是示意，不应拿其屏幕测量值代替 CSV 中的质心距离。

## 10. 导出文件

默认保存在 PyMOL 当前工作目录。为了避免找不到文件，推荐总是明确给出 `output_dir`。

### prefix_interactions.csv

每行是一种“相互作用类型–残基对”的最佳几何记录，包含：

- 类型；
- 可信度/限定词；
- 双方残基和原子；
- 距离；
- 额外角度或范德华间隙；
- 完整判据。

同一个残基对可能同时属于氢键和盐桥，这是正常的不同层级描述。

### prefix_residues.csv

每行是一个界面残基，包含：

- 属于 partner1 还是 partner2；
- 残基编号；
- 是否属于 4.0 Å 接触界面；
- 是否属于 ΔSASA 界面；
- ΔSASA；
- 最短接触距离；
- 原子接触数；
- 参与的相互作用类型。

### prefix_summary.json

这是审计主文件，保存：

- 输入选择；
- 原子数和界面残基数；
- BSA；
- 各相互作用数量；
- 本次实际使用的全部阈值；
- 假设和警告；
- 方法参考链接；
- 脚本版本。

复现实验时应把 JSON 与原始 PDB/mmCIF 一起保存。

## 11. 结构准备建议

分析前至少检查：

1. 抗体、抗原链 ID 是否正确；
2. 是否错误地同时保留了多个 biological assembly 副本；
3. 是否存在 alternate locations；
4. 关键侧链原子是否缺失；
5. 结构分辨率和 wwPDB validation report；
6. 是否需要补氢和指定 pH/质子化状态；
7. 是否有结晶接触被误认为生物学界面。

补氢不会自动解决所有问题。His、Asp/Glu、链末端和配体周围的质子化状态可能需要专业工具及人工判断。

## 12. 当前明确不做的内容

0.1.0 不自动判定：

- 水桥；
- 金属配位；
- 卤键；
- 配体–蛋白相互作用；
- 非标准氨基酸完整化学类型；
- 抗体 CDR 编号；
- 结合自由能；
- 热点残基能量贡献。

水桥高度依赖结晶水、分辨率和水分子占有率；金属配位和非标准残基需要更完整的化学感知。对于这些需求，建议用 PDBe Arpeggio 或专门的化学信息学后端交叉验证，而不是在纯 PyMOL 脚本里给出过度确定的结论。

## 13. 结果怎么用于抗体改造

建议按照以下顺序判断，而不是看到接触就直接突变：

1. 从 `residues.csv` 找到抗体侧界面残基；
2. 区分核心界面（ΔSASA 较大、原子接触多）与边缘界面；
3. 查看是否参与盐桥、氢键或芳香相互作用；
4. 检查是否位于 CDR、框架区或可能影响结构稳定性的位置；
5. 对候选突变进行结构建模和冲突检查；
6. 用 Rosetta/MD/自由能方法和实验验证排序。

界面残基不是天然的“可突变位点”。一个表面疏水残基如果同时是结合热点，直接改成亲水残基可能降低聚集倾向，也可能显著损害亲和力。

## 14. 常见报错

### partner1 selection is empty

链名或对象名写错。先执行：

```pml
select test_partner1, complex and chain H+L
count_atoms test_partner1
```

### partner1 and partner2 overlap

两个选择包含了同一批原子。重新限定链或对象。

### Atom identifiers are not unique

通常表示两个不同对象使用了重复链号和残基号。最稳妥的做法是把真正的一套复合物放入同一个对象，并给每条链唯一 chain ID；必要时也可设置不同 segi。

### No putative hydrogen bonds were found

检查侧链是否完整、是否有氢、原子类型是否被 PyMOL 正确识别，以及结构是否经过合理准备。不要为了得到更多氢键而随意放宽阈值。

### CSV 在哪里

命令中明确写：

```pml
output_dir=C:/Users/Administrator/Desktop/interface_results
```

运行结束时，PyMOL 控制台也会打印三个文件的绝对路径。

## 15. 验证策略

仓库自动检查包括：

1. Python 3.8、3.11、3.13 语法编译；
2. 空间网格、角度、π–π、T 型堆积和阳离子–π单元测试；
3. 在 Ubuntu 开源 PyMOL 中下载并分析公开 1MLC；
4. 检查界面残基数量、BSA 合理范围、接触和氢键候选；
5. 检查 CSV/JSON 是否确实生成且 JSON 可重新读取。

自动测试通过只能证明实现行为与当前规则一致，不能替代对每个输入结构的质量审查或实验验证。
