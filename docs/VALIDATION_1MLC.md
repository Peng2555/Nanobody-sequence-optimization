# 1MLC 集成验证记录

## 验证对象

- PDB：1MLC
- 结构：D44.1 Fab–鸡卵清溶菌酶复合物
- 方法：X 射线晶体学
- 分辨率：2.50 Å
- 抗体：轻链 A + 重链 B
- 抗原：溶菌酶链 E
- 运行环境：Ubuntu 24.04、Python 3.12、开源 PyMOL 3.2.0a0
- 脚本版本：0.1.0
- 首次完整通过的工作流：[Validate PyMOL interface analyzer #23](https://github.com/Peng2555/Nanobody-sequence-optimization/actions/runs/32238265620)

1MLC 的不对称单元中有两套几乎相同的 Fab–溶菌酶复合物。验证只使用第一套 A+B:E，避免把两个复合物副本混入同一次界面计算。

## 验证命令

```python
result = ab_interface(
    "1mlc and chain A+B",
    "1mlc and chain E",
    prefix="validation_1mlc",
    output_dir=".",
    labels=0,
    quiet=1,
)
```

## 实测结果

| 指标 | 结果 |
|---|---:|
| BSA | 686.6 Å² |
| 抗体侧 4.0 Å 接触残基 | 19 |
| 抗体侧 ΔSASA 界面残基 | 22 |
| 抗体侧两种定义并集 | 22 |
| 抗原侧 4.0 Å 接触残基 | 16 |
| 抗原侧 ΔSASA 界面残基 | 19 |
| 抗原侧两种定义并集 | 19 |
| close contact 残基对 | 33 |
| 氢键候选残基对 | 10 |
| 盐桥残基对 | 3 |
| 疏水接触残基对 | 3 |
| 阳离子–π残基对 | 3 |
| 范德华接触残基对 | 26 |

本结构在默认规则下没有报告 π–π、T 型堆积、跨界面二硫键或空间冲突。

## 自动检查内容

集成测试不仅检查脚本能否导入，还真正启动 PyMOL、下载结构、计算 SASA/BSA、识别相互作用并导出文件。当前断言包括：

- 抗体和抗原两侧都有合理数量的界面残基；
- BSA 位于预先设定的宽松合理范围 300–2000 Å²；
- 至少存在 8 个 close-contact 残基对；
- 至少识别到 1 个候选氢键；
- 两个 CSV 和一个 JSON 均成功生成；
- JSON 可以重新读取且版本字段一致。

每次 PR 更新还会在 Python 3.8、3.11、3.13 上执行语法编译和纯几何单元测试。

## 解释边界

RCSB 对 1MLC 的说明指出，界面中存在 3 个埋藏水分子，它们参与氢键并改善表面互补性。本工具 0.1.0 明确排除水分子和水桥，因此这里的 10 个氢键候选只代表直接、由 PyMOL 供体/受体类型和重原子方向支持的候选，不能与包含水桥的文献氢键总数直接比较。

同理，BSA 是 PyMOL 点采样结果；不同软件的原子半径、探针、氢原子处理和采样精度不同，数值允许存在方法学差异。验证的意义是证明实现可运行、结果数量级和内部关系合理，而不是宣称所有算法会给出完全相同的数字。

## 可审计输出

工作流保存以下文件作为 GitHub Actions artifact：

- `validation_1mlc_interactions.csv`
- `validation_1mlc_residues.csv`
- `validation_1mlc_summary.json`

JSON 中记录了本次实际阈值、警告、假设和方法参考。
