# HiSA 高阶质量通量的 SA 雅可比依赖

sa_hisa_second_ring_velocity.patch 为 DAFoam 的可压缩 SA 残差添加第二圈 U 依赖。
HiSA 的质量通量由高阶左右重构计算，SA 对流使用该通量，因而会依赖第二圈速度。
现有依赖表遗漏 U 会让完整着色有限差分雅可比遗漏非零项，并可能造成颜色混叠。

补丁根据本机 DAFoam 源码生成；相同差分步长的列检查已确认遗漏项，证据见
../../tests/highOrderNewton/connectivity_report.json。
这是待应用、待重编译验证的依赖补丁，不会因构建 HiSA 而自动生效。
本次没有修改 /home/dafoamuser/dafoam 下的已安装源码或库。

在单独 DAFoam checkout 中可先执行 git apply --check，再应用并构建普通版和 ADR 版。
修改依赖表后必须重新生成着色文件；旧着色缓存不可直接沿用。
验收须比较修复后的完整组装列与独立差分，再验证方向乘积和实际 Newton 下降。
增加保守依赖不改变流场方程，但可能增加稀疏存储和着色工作量。
