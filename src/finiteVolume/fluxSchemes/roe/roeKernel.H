// GPL-3.0-or-later. Calorically perfect gas, unit-normal Roe flux and frozen blocks.
#ifndef hisaRoeKernel_H
#define hisaRoeKernel_H
#include <cmath>
namespace hisaRoe {
template<class S> struct State { S q[5]; };
template<class S> struct Average { S u[3], H, a, gamma; };
template<class S> S pressure(const State<S>& w, const S& g) {
    S m2=0; for(int i=0;i<3;++i)m2+=w.q[i+1]*w.q[i+1];
    return (g-1)*(w.q[4]-m2/(2*w.q[0]));
}
template<class S> State<S> conservative(const S& rho,const S* u,const S& p,const S& g) {
    State<S> w; w.q[0]=rho; S u2=0;
    for(int i=0;i<3;++i){w.q[i+1]=rho*u[i];u2+=u[i]*u[i];}
    w.q[4]=p/(g-1)+rho*u2/2;return w;
}
template<class S> State<S> physical(const State<S>& w,const S* n,const S& g) {
    State<S> f;S un=0;for(int i=0;i<3;++i)un+=w.q[i+1]*n[i]/w.q[0];
    S p=pressure(w,g);f.q[0]=w.q[0]*un;
    for(int i=0;i<3;++i)f.q[i+1]=w.q[i+1]*un+p*n[i];
    f.q[4]=(w.q[4]+p)*un;return f;
}
template<class S> Average<S> average(const State<S>& l,const State<S>& r,const S& g) {
    using std::sqrt;Average<S> b;b.gamma=g;S sl=sqrt(l.q[0]),sr=sqrt(r.q[0]),u2=0;
    for(int i=0;i<3;++i){b.u[i]=(l.q[i+1]/sl+r.q[i+1]/sr)/(sl+sr);u2+=b.u[i]*b.u[i];}
    b.H=((l.q[4]+pressure(l,g))/sl+(r.q[4]+pressure(r,g))/sr)/(sl+sr);
    b.a=sqrt((g-1)*(b.H-u2/2));return b;
}
template<class S> S entropyAbs(const S& x,const S& delta) {
    S ax=x>=0 ? x : -x;
    return delta>0 && ax<delta ? (x*x/delta+delta)/2 : ax;
}
template<class S> State<S> dissipate(const Average<S>& b,const S* n,const State<S>& v,const S& entropy,const bool cfl3d=false) {
    S un=0,u2=0,udm=0,dmn=0;
    for(int i=0;i<3;++i){un+=b.u[i]*n[i];u2+=b.u[i]*b.u[i];udm+=b.u[i]*v.q[i+1];dmn+=n[i]*v.q[i+1];}
    S dp=(b.gamma-1)*(v.q[4]-udm+u2*v.q[0]/2),dn=dmn-un*v.q[0];
    S am=(dp-b.a*dn)/(2*b.a*b.a),ap=(dp+b.a*dn)/(2*b.a*b.a),ac=v.q[0]-dp/(b.a*b.a);
    // CFL3D fhat.F: epsilon = entropy*(a+|u|+|v|+|w|), delta=2*epsilon.
    // Keep the scale active for AD; only the configuration switch is passive.
    S delta=entropy*b.a;
    if(cfl3d){S speed=b.a;for(int i=0;i<3;++i)speed+=b.u[i]>=0 ? b.u[i] : -b.u[i];delta=2*entropy*speed;}
    S lm=entropyAbs(S(un-b.a),delta),lp=entropyAbs(S(un+b.a),delta),lc=entropyAbs(un,cfl3d ? delta : S(0));
    State<S> d;d.q[0]=lm*am+lc*ac+lp*ap;S ushear=0;
    for(int i=0;i<3;++i){S shear=v.q[i+1]-b.u[i]*v.q[0]-n[i]*dn;ushear+=b.u[i]*shear;
        d.q[i+1]=lm*am*(b.u[i]-b.a*n[i])+lp*ap*(b.u[i]+b.a*n[i])+lc*(ac*b.u[i]+shear);}
    d.q[4]=lm*am*(b.H-b.a*un)+lp*ap*(b.H+b.a*un)+lc*(ac*u2/2+ushear);return d;
}
template<class S> State<S> flux(const State<S>& l,const State<S>& r,const S* n,const S& g,const S& entropy,const bool cfl3d=false) {
    State<S> v;for(int i=0;i<5;++i)v.q[i]=r.q[i]-l.q[i];
    auto d=dissipate(average(l,r,g),n,v,entropy,cfl3d);auto fl=physical(l,n,g),fr=physical(r,n,g);
    for(int i=0;i<5;++i) { fl.q[i]=(fl.q[i]+fr.q[i]-d.q[i])/2; }
    return fl;
}
template<class S> State<S> physicalAction(const State<S>& w,const State<S>& v,const S* n,const S& g) {
    S un=0,u2=0,udm=0,dmn=0,u[3];
    for(int i=0;i<3;++i){u[i]=w.q[i+1]/w.q[0];un+=u[i]*n[i];u2+=u[i]*u[i];udm+=u[i]*v.q[i+1];dmn+=v.q[i+1]*n[i];}
    S dp=(g-1)*(v.q[4]-udm+u2*v.q[0]/2),du=(dmn-un*v.q[0])/w.q[0];
    State<S> f;f.q[0]=dmn;for(int i=0;i<3;++i)f.q[i+1]=v.q[i+1]*un+w.q[i+1]*du+dp*n[i];
    f.q[4]=(v.q[4]+dp)*un+(w.q[4]+pressure(w,g))*du;return f;
}
template<class S> void frozenBlocks(const State<S>& l,const State<S>& r,const S* n,const S& g,const S& entropy,S L[5][5],S R[5][5],const bool cfl3d=false) {
    auto b=average(l,r,g);
    for(int j=0;j<5;++j){State<S> v;for(int i=0;i<5;++i)v.q[i]=S(i==j ? 1 : 0);
        auto d=dissipate(b,n,v,entropy,cfl3d),al=physicalAction(l,v,n,g),ar=physicalAction(r,v,n,g);
        for(int i=0;i<5;++i){L[i][j]=(al.q[i]+d.q[i])/2;R[i][j]=(ar.q[i]-d.q[i])/2;}}
}
// Geometrically corrected U-MUSCL: linear-exact before limiting, including off-midpoint faces.
template<class S> S reconstruct(const S& q,const S& other,const S& gradToFace,const S& gradToNeighbour,const S& kappa) {
    return q+gradToFace+kappa*(other-q-gradToNeighbour)/2;
}
}
#endif
