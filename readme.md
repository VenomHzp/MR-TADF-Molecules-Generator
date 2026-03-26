# 🧬 MR-TADF Generator: 基于图论与严苛化学约束的分子演化引擎


> **面向下一代高色纯度发光材料的高通量虚拟发现平台**
> 本项目致力于解决传统的一维 SMILES 随机变异在多重共振（MR-TADF）刚性稠环体系生成中极易产生“化学乱码”的痛点，创新性地提出了一套基于 **2D 拓扑网络外科手术** 的定向分子演化算法。

---

## 💡 为什么开发这个项目？

多重共振热活化延迟荧光（MR-TADF）材料（如 DABNA 衍生物）因其极窄的发射光谱（FWHM）和高发光效率，在 OLED 领域备受瞩目。
然而，MR-TADF 分子的核心在于其**精密稠合的刚性共振网络**。传统的分子生成模型（如基于 SMILES 的 RNN/VAE 或随机变异算法）在处理此类分子时面临灾难性问题：
1. 极易破坏 Kekule 芳香性，生成无法被量化软件识别的“化学垃圾”。
2. 缺乏对 B/N/O/S 杂原子推拉电子交替网络的化学直觉，产生大量无物理意义的死位点（如 N-N、B-O 直接相连）。

**本项目采用底层 `RDKit.Chem.RWMol` 节点级图编辑技术，赋予了生成算法顶尖化学家的“直觉护城河”。**

---

## 🚀 核心技术与四大演化法则

本引擎内置了四大基于真实《Nature》/《JACS》顶刊合成路线的“虚拟化学手术”法则。所有的变异均采用**“事前安全验证 (Pre-validation)”**机制，保证 100% 的生成合法率，极致压榨算力。

### 🔪 法则 1：精准氮掺杂 (Targeted N-Doping)
模拟氮杂芳环合成，精细调控 EA 与 LUMO 能级。
* 🔒 **外围 C-H 保护锁**：自动避开桥头稠合碳，防止底层多环芳香网络崩溃。
* 🔒 **绝对杂原子隔离 (Absolute Heteroatom Isolation)**：强制审查拓扑邻居，确保新生成的 N 周围 100% 为碳网络，维持交替共振发光机制。

### 🌉 法则 2：邻接桥连平面化 (Topological Planarization)
模拟真实合成中的硼化/胺化锁环反应，引入 `O=S=O`, `P=O`, `C=O`, `N-Ph`, `O/S/Se` 等复杂前沿桥连基团。
* 🔒 **纯图论海湾探测器**：彻底弃用不可靠的 3D 构象测距，使用**精确 5 步最短拓扑寻路**完美锁定分子边缘的“海湾 (Bay Region)”。
* 🔒 **同环排斥锁 (Fused Ring Lock)**：基于底层 RingInfo，精准识别并拦截已被骨架锁死的“V型发散碳”，彻底根除假阳性。

### 🧩 法则 3：导向边缘稠环扩增 (Directed $\pi$-Extension)
在不破坏内部共振网络的前提下扩大共轭面，实现靶向波长红移。
* 采用**节点级无损缝合技术**，摒弃极易报错的 SMARTS 反应，支持苯环、呋喃、苯并呋喃、N-苯基吲哚的完美接枝。
* 支持非对称片段的**正反双向随机缝合**，极大丰富生成空间的构象多样性。

### 🔗 法则 4：海湾脱氢偶联 (Scholl-type Cyclization)
模拟 Scholl 氧化交叉偶联反应，将悬垂芳香环强制闭合成 5 元刚性稠环。
* 🔒 **全局原子消耗锁 (Consumption Lock)**：解决旋转基团导致的拓扑对称性 Bug，杜绝同一个碳原子被错误匹配两次引发的五价碳崩溃。
* 极限提升分子刚性，压榨重组能 ($\lambda$)。

---

## ⚙️ 安装与依赖

本项目极为轻量化，只需安装 RDKit 即可运行全部功能：
```bash
# 推荐使用 conda 环境
conda create -n mr_gen python=3.9
conda activate mr_gen
conda install -c conda-forge rdkit
````

-----

## 🏃 快速开始 (Quick Start)

你可以直接运行自动化控制脚本，一键生成万级别的未知 MR 分子：

```bash
python mr_generator_main.py
```

### 代码调用示例

```python
from mr_generator_main import run_dft_ml_generator

# 1. 准备纯净的“裸骨架”作为创世种子
seeds = [
   C1=CC(=CC=C1)N1C2C(B3C4=CC=CC=C4N(C4C=CC=CC=4)C4C=CC=C1C=43)=CC=CC=2
]

# 2. 启动蒙特卡洛定向演化引擎
novel_cores = run_dft_ml_generator(
    seed_smiles_list=seeds, 
    num_iterations=5000,    # 迭代次数
    max_heavy_atoms=55      # 体积膨胀限制
)
```

### 📊 自动落盘与输出

引擎运行结束后，会自动在当前目录下生成：

1.  `AI_Generated_MR_Cores.csv`：包含所有全新的唯一 SMILES 及其重原子数，完美对接 Pandas。
2.  `/generator_results/Top20_Showcase.png`：随机抽取 20 个生成骨架的高清 2D 渲染图谱，供直观审核。

-----

## 🔮 未来工作 (Future Work)

本项目生成的结构 100% 满足高阶量化计算的拓扑要求。设计输出的 CSV 数据库可无缝对接 DFT-ML 工作流：

  - [ ] 结合 `GFN2-xTB` 进行自动化大通量基态构象预优化。
  - [ ] 调用 `PySCF` / `gpu4pyscf` 执行 sTD-DFT 计算。
  - [ ] 自动化提取 $\Delta E_{ST}$ 和振子强度 ($f$)，构建发光性能的机器学习回归模型 (ML Regression)。

-----

## 📝 贡献与许可 (License)

本项目为学术研究开发，欢迎提出 Issue 或 Pull Request 探讨底层算法逻辑！
遵循 [MIT License](https://www.google.com/search?q=LICENSE) 开源协议。

*Designed with ❤️ for Materials Discovery.*