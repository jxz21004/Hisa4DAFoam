#include "roeKernel.H"
#include <codi.hpp>
#include <cassert>
#include <cmath>
#include <iostream>
using hisaRoe::venkatIncrement;
static void close(double a,double b,double tol=1e-9) {
    assert(std::abs(a-b)<=tol*(1+std::abs(a)+std::abs(b)));
}
int main() {
    close(venkatIncrement(0.,0.,.01),0.);
    close(venkatIncrement(1.,1.,0.),.75);
    close(venkatIncrement(1.,2.,0.),1.);
    close(venkatIncrement(10.,0.,.01),10.*.01/200.01);
    for(double d: {0.,1e-12,.03,1.,10.}) {
        const double a=.3,e=.01;
        close(venkatIncrement(-d,-a,e),-venkatIncrement(d,a,e));
        close(venkatIncrement(d*220,a*220,e*220*220),220*venkatIncrement(d,a,e));
        // All three active inputs, including geometric smoothing, retain derivatives.
        for(int j=0;j<3;++j) {
            codi::RealForward x=d,y=a,z=e;
            if(j==0)x.setGradient(1);if(j==1)y.setGradient(1);if(j==2)z.setGradient(1);
            codi::RealForward out=venkatIncrement(x,y,z);
            double h=1e-6,plus[3]={d,a,e},minus[3]={d,a,e};
            plus[j]+=h;minus[j]-=h;
            double fd=(venkatIncrement(plus[0],plus[1],plus[2])-venkatIncrement(minus[0],minus[1],minus[2]))/(2*h);
            close(out.getGradient(),fd,2e-7);
        }
    }
    // Distinct upper/lower bounds meet with unit increment derivative at zero.
    const double h=1e-8;
    close(venkatIncrement(h,.4,.01)/h,1.,1e-7);
    close(venkatIncrement(-h,-.2,.01)/(-h),1.,1e-7);
    std::cout<<"Venkat values, sign/scale invariance, zero join and CoDi derivatives passed\n";
}
