#ifndef HisaCubicLineSearch_H
#define HisaCubicLineSearch_H
#include <algorithm>
#include <cmath>
#include <functional>
namespace hisaNK {
struct TrialMerit { double value=0, turbulence=0; bool valid=false; };
struct SearchResult { bool accepted=false, minimumStep=false, prelimited=false; double alpha=0; };
// ADflow LSCubic policy, retaining the adapter's physical admissibility checks.
// Invalid states are never accepted, including at the minimum step.
inline SearchResult cubicSearch(double f0, double slope, double turb0,
    double minStep, bool relaxed, const std::function<TrialMerit(double)>& evaluate)
{
    SearchResult out;
    double a=1,prevA=0,prevF=0;
    TrialMerit t=evaluate(a);
    if(relaxed && (!t.valid || t.turbulence>2*turb0)) {
        out.prelimited=true;
        a=t.valid ? std::max(turb0/t.turbulence,.1) : .5;
        for(int i=0;i<10;++i) {
            t=evaluate(a);
            // Match ADflow's preliminary branch (no lambda in its Armijo RHS).
            if(t.valid && (t.value<=f0+.01*slope || a<minStep)) {
                out.accepted=true;out.alpha=a;
                out.minimumStep=t.value>f0+.01*slope;return out;
            }
            a*=.5;
        }
        return out;
    }
    for(int i=0;i<30;++i) {
        if(t.valid && (t.value<=f0+.01*a*slope || (relaxed && a<=minStep))) {
            out.accepted=true;out.alpha=a;
            out.minimumStep=t.value>f0+.01*a*slope;return out;
        }
        if(!t.valid) {a*=.5;prevA=0;}
        else {
            double candidate=a/2;
            if(prevA==0) {
                double den=2*(t.value-f0-a*slope);
                if(den!=0)candidate=-slope*a*a/den;
            } else {
                double t1=t.value-f0-a*slope,t2=prevF-f0-prevA*slope;
                double c=(t1/(a*a)-t2/(prevA*prevA))/(a-prevA);
                double b=(-prevA*t1/(a*a)+a*t2/(prevA*prevA))/(a-prevA);
                candidate=c!=0?(-b+std::sqrt(std::max(b*b-3*c*slope,0.)))/(3*c)
                    :(b!=0?-slope/(2*b):a/2);
            }
            prevA=a;prevF=t.value;
            if(!std::isfinite(candidate))candidate=a/2;
            a=std::min(.5*a,std::max(.1*a,candidate));
        }
        if(a<1e-8 && !relaxed)break;
        t=evaluate(a);
    }
    return out;
}
}
#endif
