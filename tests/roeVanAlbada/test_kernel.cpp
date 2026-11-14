#include "roeKernel.H"
#include <codi.hpp>
#include <cassert>
#include <cmath>
#include <iostream>
using hisaRoe::vanAlbadaTypeIncrement;
void close(double a,double b,double tol=1e-9){assert(std::abs(a-b)<tol*(1+std::abs(a)+std::abs(b)));}
int main(){
    close(vanAlbadaTypeIncrement(1.,0.,0.,.01),1.);
    for(double a: {-2.,-.1,0.,.1,2.})for(double b: {-2.,-.1,0.,.1,2.}){
        double d=.3,e=.02;
        double f=vanAlbadaTypeIncrement(1.,a,b,e);
        assert(f>0 && f<=1+1e-15);
        close(f,vanAlbadaTypeIncrement(1.,b,a,e));
        close(vanAlbadaTypeIncrement(-d,-a,-b,e),-vanAlbadaTypeIncrement(d,a,b,e));
        close(vanAlbadaTypeIncrement(d*220,a*220,b*220,e*220*220),220*vanAlbadaTypeIncrement(d,a,b,e));
        for(int j=0;j<4;++j){
            codi::RealForward x=d,y=a,z=b,t=e;
            if(j==0)x.setGradient(1);if(j==1)y.setGradient(1);if(j==2)z.setGradient(1);if(j==3)t.setGradient(1);
            codi::RealForward out=vanAlbadaTypeIncrement(x,y,z,t);
            double h=1e-6,p[4]={d,a,b,e},m[4]={d,a,b,e};p[j]+=h;m[j]-=h;
            double fd=(vanAlbadaTypeIncrement(p[0],p[1],p[2],p[3])-vanAlbadaTypeIncrement(m[0],m[1],m[2],m[3]))/(2*h);
            close(out.getGradient(),fd,1e-7);
        }
    }
    // Linear-exact at non-midpoint faces and for all supported kappa values.
    for(double k:{0.,1./3.,.5,1.}){
        double q=3.,jump=.7,gface=.21;
        double unlimited=hisaRoe::reconstruct(q,q+jump,gface,jump,k)-q;
        close(q+vanAlbadaTypeIncrement(unlimited,jump,jump,.01),q+gface);
    }
    // Squared type intentionally retains symmetric extrema, unlike TVD clipping.
    close(vanAlbadaTypeIncrement(.3,1.,-1.,.01),.3);
    std::cout<<"van Albada type bounds, scaling, linear preservation and CoDi derivatives passed\n";
}
