#include "incrementProtectionKernel.H"
#include <cassert>
#include <cmath>
#include <iostream>
using namespace hisaIncrementProtection;
int main(){
 for(double r:{-.2,-.1,0.,.5,100.})assert(coefficient(r)==1.);
 for(double r:{-1e100,-100.,-2.,-1.,-.21}) {
  double a=coefficient(r);assert(a>0 && a<1);
  assert(r*a>=-.5 && r*a<-.2);
 }
 assert(std::abs(coefficient(-1.)-1./2.6)<1e-15);
 assert(std::abs(coefficient(-.2-1e-10)-1)<1e-8);
 // A pure momentum increment has zero linear pressure change at rest,
 // but the nonlinear internal energy can still become negative.
 double alpha=1;int n=0;
 while(2.5-.5*alpha*alpha*100<1.25){alpha*=.5;++n;}
 assert(n>0 && 2.5-.5*alpha*alpha*100>=1.25);
 std::cout<<"PASS CFL3D coefficient, identity, continuity, loss bounds, nonlinear backtrack example\n";
}
