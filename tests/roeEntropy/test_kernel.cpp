#include "roeKernel.H"
#include <cassert>
#include <cmath>
#include <iostream>
using namespace hisaRoe;
void check(double a,double b){assert(std::abs(a-b)<1e-10*(1+std::abs(a)+std::abs(b)));}
int main(){
 double n[3]={1,0,0},nr[3]={-1,0,0},u[3]={0,0,0},g=1.4;
 auto l=conservative(1.,u,1.,g),r=conservative(2.,u,1.,g);
 auto av=average(l,r,g);double eps=.1*av.a;
 auto f=flux(l,r,n,g,.1,true);check(f.q[0],-.5*eps);check(f.q[4],0.);check(f.q[1],1.);
 auto rev=flux(r,l,nr,g,.1,true);for(int i=0;i<5;i++)check(f.q[i],-rev.q[i]);
 auto same=flux(l,l,n,g,.1,true),phys=physical(l,n,g);for(int i=0;i<5;i++)check(same.q[i],phys.q[i]);
 auto off=flux(l,r,n,g,0.,true),old=flux(l,r,n,g,0.);for(int i=0;i<5;i++)check(off.q[i],old.q[i]);
 double h=1e-6,delta=2*eps;check(entropyAbs(delta,delta),delta);
 assert(std::abs((entropyAbs(delta+h,delta)-entropyAbs(delta-h,delta))/(2*h)-1)<1e-5);
 std::cout<<"CFL3D entropy contact, conservation, consistency, zero-switch, C1 checks passed\n";
}
