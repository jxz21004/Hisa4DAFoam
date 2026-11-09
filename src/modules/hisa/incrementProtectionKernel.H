// GPL-3.0-or-later. CFL3D update.F scalar increment damping.
#ifndef hisaIncrementProtectionKernel_H
#define hisaIncrementProtectionKernel_H
namespace hisaIncrementProtection {
template<class S> S coefficient(const S& r) {
    return r >= S(-0.2) ? S(1) : S(1)/(S(0.6)-S(2)*r);
}
}
#endif
