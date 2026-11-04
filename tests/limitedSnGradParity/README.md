# limitedSnGrad 的普通版 / ADR 不一致

在同一个高阶 L3 保存状态，普通版和 ADR 的 U/p/T/nuTilda、rho/rhoE、nut/alphat、
thermo:mu/alpha 内部场完全一致，但残差相差约 15.96%（L2 相对差），最大绝对差为
2703.214。质量残差完全一致。重复评估及坐标不变的 movePoints 均不消除差异。

仅将有限非正交修正改成 uncorrected 后，两版残差最大差降至 2.22e-12。
更关键的是：**原 ADR 的 limited 残差与普通版 uncorrected 残差相等至舍入精度**。
AD 方向导数与 ADR 自身的有限差分明显比与普通版有限差分更吻合。
这是求导所对应的离散方程不一致，不可用梯度验收容差掩盖。

## 根因及修复

已安装 OpenFOAM v2506 的 limitedSnGrad.H 声明顺序为：

```cpp
tmp<snGradScheme<Type>> correctedScheme_;
scalar limitCoeff_;
```

构造 correctedScheme_ 时调用 lookupCorrectedScheme，先写 limitCoeff_；
但 C++ 按声明顺序构造成员，AD scalar 的构造发生在这次赋值之后，覆盖已读系数。
普通 double 没有对应的默认构造清零操作，所以普通版保留设定系数，而 ADR 丢失系数。

HiSA 新增独立运行时方案 hisaLimited：系数成员先声明并初始化，再执行读取方案的
辅助函数。其公式保持 OpenFOAM v2506 原式，保留原文件 GPL 许可及版权。
使用新类型名，避免修改已安装 OpenFOAM 文件或改变现有类型的 ABI / 注册表。

## 使用

需要编译包含新实现的 HiSA 普通版与 ADR 库，并同时加载；fvSchemes 中替换：

```cpp
laplacian(muEff,U)    Gauss linear hisaLimited 0.33;
laplacian(alphaEff,e) Gauss linear hisaLimited 0.33;
// 原 default 中的 limited 0.5 同样替换为 hisaLimited 0.5。
```

snGradSchemes 中原有的 limited 也须替换；不能改动 limitedLinear 等其他方案名称。
wingAero/prepare.py 已为新测试应用这一配置，历史案例不改写。

## 验证流程

highOrderNewton/check_parity.py 对同一保存状态比较两版内部场及残差，并做无位移几何刷新。
highOrderNewton/check_derivatives.py 比较 4 类状态扰动和 4 类残差投影，扫描 6 个步长（温度额外细化至 1e-10）。
AD 输出已含 normalizeStates，有限差分相应扰动 x+h*scale*v；记录与 native/ADR
各自有限差分的误差。随机块投影不是完整雅可比验证，也不是稳态气动梯度验证。

before_parity.json、uncorrected_parity.json、before_derivatives.json 保存修复前证据。
后续修复结果需同时满足：普通版结果保持一致、native/ADR 残差一致、方向导数检查。
该初始化问题位于 OpenFOAM 依赖层，不是 wVanLeer 对流限制器的错误；修复主要保证
伴随方程与原流场离散一致，本身不保证普通版高阶流场停滞消失。

修复复核通过：普通版残差变化为 0，native/ADR 残差相对差为 1.035e-16；
16 个块投影中各自最优步长对应的最大相对误差为 8.823e-6（门槛 1e-4）。
温度扰动最初的最小步长 1e-8 不足以通过，细化后通过；完整扫描保留在 after_derivatives.json。
本轮未改变验收门槛。普通及 ADR 隔离库均编译、加载通过。

当前普通版 / ADR 修复库路径：`dafoam_case/validation/hisa_limited_fix/build/lib`（工作区隔离构建，未覆盖系统安装）。
