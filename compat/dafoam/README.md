# HiSA 高阶质量通量的 SA 雅可比依赖

`sa_hisa_second_ring_velocity.patch` 只修改 DAFoam 的
`DAStateInfoHisaFoam.C`：在读取共用模型依赖后，向 `nuTildaRes` 的第二圈
加入 `U`（若已存在则不重复添加）。不修改共用 `DASpalartAllmaras.C`。
HiSA 的质量通量由高阶左右重构计算，SA 对流使用该通量，因而会依赖第二圈速度。

2026-11-04 的早期诊断通过修改共用 SA 表确认了原因；其历史结果保存在
`tests/highOrderNewton/connectivity_fixed_report.json`。该历史构建不代表本文件
当前补丁的修改位置。现在的补丁将影响范围限制在 HiSA 专属适配中。

用户克隆的相邻 `dafoam` 仓库包含正式修改、独立构建安装脚本和测试入口，见
`../dafoam/tests/hisaConnectivity/README.md`（相对于两个仓库的共同父目录）。
从该仓库运行 `python tools/build_hisa_local.py --hisa-lib ABSOLUTE_HISA_LIB_DIR`，
再加载 `build/hisa-local/activate.sh` 即可使用工作区安装，不替换系统安装。
该方式适用于新版本源码，也会重新构建匹配的 Python 扩展。

本目录旧版 `build_isolated.py` 用于复制当前系统安装的源码和缓存后仅重编译
普通/ADR C++ 库；它也已改为应用 HiSA 专属补丁。不要用此方式混用不同版本的
Python 扩展与 C++ 库。

修改依赖后必须重新生成着色。增加正确依赖可能增加稀疏存储和着色工作量，
不改变流场残差方程，也不能单独解决原生 HiSA 非线性收敛停滞。
