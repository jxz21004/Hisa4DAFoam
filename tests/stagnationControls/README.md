# L3 高阶重构与伪时间推进对照

基准为 6c1d633 的 AUSMPlusUp 连续压力修正、连续黏性功及 GMRES 提前退出。
三项追加测试各执行 2000 步：

1. 一阶：仅将 reconstruct(rho/U/T) 的 wVanLeer 换为零斜率重构。
2. 高阶低 CFL：仅将 pseudoCoNumMax 从 100000 改成 5。
3. 一阶初始化高阶：由第一项的第 2000 步初始化原高阶配置，再推进 2000 步。
   删除重启 pseudoCoField 和 uniform 缓存，CFL 初值恢复为 1。

前两项从相同均匀来流开始。GMRES rtol=0.05、atol=1e-12、nKrylov=40、maxIter=60；
网格、边界、湍流、黏性项均不变。每个测试实际加载的 native/ADR 库哈希见 comparison.json。

| 算例 | 最终/均匀初始残差 L2 比 | 最终原始残差最大值 | 最终 HiSA 归一化指标 |
|---|---:|---:|---:|
| 原高阶基准 | 2.575e-2 | 103058 | 1.890e-2 |
| 一阶 | 1.444e-6 | 17.375 | 2.010e-7 |
| 高阶、CFL 上限 5 | 4.412e-3 | 26899 | 7.421e-4 |
| 一阶初始化后恢复高阶 | 未评估对应重启初场比值 | 17567 | 3.994e-4 |

L2 比使用 DAFoam getResiduals 的固定定义，在各自离散方程下评估初始及最终状态；
HiSA 逐步归一化指标不是该相对残差。跨一阶/高阶的结果属于不同离散方程，不能把
一阶解当作高阶精度验证。第一项后期仍有振荡，原始残差最大项为 SA（17.375），
动量仍为 11.177，因此没有达到既定稳态验收门槛。其他高阶案例也未通过。
未启动梯度测试。

## 源码检查结论

- hisaModuleTimestepping.C 中 CFL 更新基于 first residual / current residual，
  再限制增减倍率及范围。残差停滞时提高 CFL 上限不一定增加实际 CFL。这是当前策略，
  本轮没有证据将其直接判为 bug。初始化改变也会改变此自适应控制的参考残差。
- wVanLeer 在 r>0 分支计算非线性限制系数，并做 min/max 限制；逐分量重构。
  分支和限幅可能导致导数不光滑，但本轮没有确认其公式存在有限跳变或实现错误。
- 一阶明显改善、高阶恢复后再度停滞，支持继续检查高阶离散与推进的配合。
  一阶也增加了数值耗散，不能仅据这个对照把高阶限制器宣布为唯一根因。

## 重现及隔离

prepare.py 复制已有修复库的构建树，并仅在隔离构建中注册 controlFirstOrder：
limiter 恒为 0，经已有权重转换得到内部面 owner/neighbour 状态；边界保持原有处理。
生产源码不增加该测试专用选项。普通版和 ADR 都需编译：wmake -j2 libso BUILD/src。
然后使用 wingAero/run_suite.py WORK/first_order --levels L3，以及 WORK/high_order_cfl5。
一阶计算完成后运行 prepare_restart.py，再使用 wingAero/run_case.py 的 --restart。

均匀初场的额外 --probe 评估不推进迭代，使用原 primalInitCondition 构造状态；
必须使用对应重构配置及对应库。summarize.py 生成数值对照和曲线。

结论：没有获得通过验收的高阶稳态解；没有据此修改生产限制器或宣称新 bug。
建议下一步对高阶残差做方向导数一致性检查，再验证带阻尼/线搜索的更新。
