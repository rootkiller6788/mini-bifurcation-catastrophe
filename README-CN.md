# 分岔与突变理论（Bifurcation & Catastrophe Theory）

**从零开始、零依赖的 C 语言实现**，涵盖非线性动力系统的分岔与突变理论。每个子模块对应 MIT、Caltech、Stanford 等顶尖大学的课程，将教科书中的方程转化为可运行的 C 代码，实现理论与实践的桥接。基于 Kuznetsov (2004)、Thom (1975)、Zeeman (1977)。

## 子模块总览

| 子模块 | 主题 | 参考课程 |
|--------|--------|-------------|
| [mini-saddle-node-hopf-pitchfork](mini-saddle-node-hopf-pitchfork/) | 鞍结分岔、Hopf 分岔、叉形分岔；中心流形约化；余维-2 分析；特征值穿越检测 | MIT 18.385J, Caltech CDS140 |
| [mini-normal-form-theory](mini-normal-form-theory/) | 规范型计算、中心流形近似约化、共振分析、普适开折 | MIT 18.385J, Cornell MAE 5790 |
| [mini-numerical-continuation](mini-numerical-continuation/) | 伪弧长延拓（Keller 方法）、预测-校正、分支切换、Moore-Penrose 延拓 | MIT 18.085, Caltech CDS140 |
| [mini-global-bifurcations-homoclinic](mini-global-bifurcations-homoclinic/) | 同宿/异宿轨道、Shilnikov 定理、Melnikov 方法、Lin 方法、八字形轨道 | MIT 18.385J, Princeton MAE 546 |
| [mini-bifurcation-in-pde-patterns](mini-bifurcation-in-pde-patterns/) | Turing 不稳定性、反应扩散 PDE 分岔、线性稳定性分析、斑图形成 | MIT 18.385J, Stanford MATH 220 |
| [mini-bifurcation-control-delay](mini-bifurcation-control-delay/) | 延迟微分方程（DDE）、Pyragas 延迟反馈控制（DFC）、Lambert W 稳定性、延迟系统 Hopf 分岔 | MIT 6.832, ETH 227-0216 |
| [mini-thom-elementary-catastrophes](mini-thom-elementary-catastrophes/) | Thom 七种初等突变：折叠、尖点、燕尾、蝴蝶、双曲/椭圆/抛物脐点 | MIT 18.385J, Harvard Math 118 |
| [mini-zeeman-applications](mini-zeeman-applications/) | Zeeman 突变机构、尖点几何、生物/社会/经济突变模型 | MIT 18.385J, Cambridge Part III |

## 设计理念

- **零外部依赖** — 纯 C（C99/C11），仅使用 `libc` 和 `libm`
- **模块自包含** — 每个目录自带 `Makefile`、`include/`、`src/`、`examples/`、`demos/`、`tests/`
- **理论到代码的映射** — 每个子模块包含 `docs/` 目录，内有课程对齐说明
- **实用演示程序** — 分岔检测器、延拓求解器、突变曲面可视化工具等

## 构建方式

每个子模块相互独立。进入子模块目录后运行：

```bash
cd mini-bifurcation-control-delay
make all    # 构建全部
make test   # 运行测试
```

需要 **GCC** 和 **GNU Make**。

## 项目结构

```
mini-bifurcation-catastrophe/
├── mini-saddle-node-hopf-pitchfork/  # 鞍结分岔、Hopf 分岔、叉形分岔；中心流形
├── mini-normal-form-theory/          # 规范型计算与普适开折
├── mini-numerical-continuation/      # 数值延拓与分支追踪
├── mini-global-bifurcations-homoclinic/ # 全局分岔：同宿轨道、Melnikov 方法
├── mini-bifurcation-in-pde-patterns/ # PDE 系统中的分岔分析与斑图形成
├── mini-bifurcation-control-delay/   # 延迟微分方程下的分岔控制
├── mini-thom-elementary-catastrophes/ # Thom 七种初等突变
└── mini-zeeman-applications/         # Zeeman 突变机构与现实应用
```

## 许可证

MIT
