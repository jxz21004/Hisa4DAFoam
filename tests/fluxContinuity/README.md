# AUSMPlusUp 黏性功的流向切换不连续

在 L3 2000 步保存状态附近做冻结状态残差测试，没有推进流场。生产源码此前将
`Up = U_f * magSf` 用于黏性功，`U_f` 根据面 Mach 数选左/右重构速度。
对流项的这一选择乘有质量流量，质量流量过零时通量仍连续；黏性功没有这个因子。
因此左右速度不同时，流向切换会引起有限的能量残差跳变。

内部面 4701（owner 1567，neighbour 1667）初始 Mach=-1.05825e-5，左右重构速度
差的模为 1.25321 m/s。只将 owner 压力提高约 1.469685 Pa（原压力 29424.4146 Pa）
即可跨零。在二分定位的零点两侧扰动压力：

| 压力扰动 ±Pa | 原实现：能量残差最大差值 | 仅改黏性面速度为线性插值 |
|---:|---:|---:|
| 1e-4 | 0.203174 | 2.56169e-3 |
| 1e-6 | 0.204452 | 2.56170e-5 |
| 1e-8 | 0.204465 | 2.56051e-7 |

这些数值是相同 DAFoam 残差定义下两次函数评估之差，不是收敛残差或梯度误差。
原实现的 `Up` 跳变模始终为 0.0271193；线性插值下为零（此次只扰动压力）。
质量及对流能量通量差值随扰动缩小，对照前后完全相同。重复评估基准状态的残差
最大差为 2.27e-12，不足以解释 0.204 的跳变。

修复仅使 AUSMPlusUp 黏性功使用 `linearInterpolate(U_)*magSf`，保持对流迎风项。
这也与已有 LaxFriedrichs 黏性功的面速度处理一致。普通版与 ADR 版共用此源码。
源码核对表明 HiSA 推进及 DAFoam 实际残差都调用该通量，黏性功公式一致；这不等于
已验证全部源项、边界、离散导数的一致性。既有压力修正的四项连续性测试仍通过。

## 重现

`prepare_probe.py BASE_BUILD NEW_BUILD` 复制**修复前**的隔离完整构建，加入诊断代码。
以普通 OpenFOAM 环境执行 `wmake -j2 libso NEW_BUILD/src`。诊断构建保留原 ADR 库；
环境开关只用于普通版残差对照，不用来验证伴随。安装库不修改。

复制保存状态案例的 `0/constant/system/final_state.npz` 到新的诊断目录，设置
`LD_LIBRARY_PATH` 指向诊断库后运行：

```sh
python tests/fluxContinuity/probe_state.py --case CASE --library-dir NEW_BUILD/lib
```

脚本逐面寻找附近可由压力扰动触发的流向反转，二分定位，再比较两种黏性面速度。
结果见 `continuity_report.json`。生产修复须另行编译普通版和 ADR 库，使用新案例副本：

```sh
python tests/fluxContinuity/verify_fix.py --case FIXED_CASE \
  --library-dir FIXED_BUILD/lib --reference continuity_report.json
```

该检查验证生产普通版残差与诊断线性插值对照一致，且差值随扰动缩小。
尚未验证多 MPI rank、壁面/运动网格的专门回归、AD 导数或修复后的长期收敛。
发现这个不连续问题不等于证明它是当前全局停滞的唯一原因，也未证明修复后梯度准确。

正式修复已通过上述冻结状态复核，结果见 `verified_fix.json`；普通版及 ADR 编译/链接通过，实际加载库的 SHA256 已记录。
