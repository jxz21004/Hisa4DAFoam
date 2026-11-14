#include "HisaCubicLineSearch.H"
#include <cassert>
#include <iostream>
using namespace hisaNK;
int main() {
 auto full=cubicSearch(.5,-1,1,.01,true,[](double a){return TrialMerit{.5*(1-a)*(1-a),1,true};});
 assert(full.accepted && full.alpha==1 && !full.minimumStep);
 // A locally non-descending merit: relaxed policy takes a small step;
 // strict policy refuses it. This exercises the former stagnation path.
 auto uphill=[](double a){return TrialMerit{.5+a,1,true};};
 auto relaxed=cubicSearch(.5,-1,1,.01,true,uphill);
 auto strict=cubicSearch(.5,-1,1,.01,false,uphill);
 assert(relaxed.accepted && relaxed.minimumStep && relaxed.alpha<=.01);
 assert(!strict.accepted);
 auto invalid=cubicSearch(.5,-1,1,.01,true,[](double){return TrialMerit{};});
 assert(!invalid.accepted);
 auto turb=cubicSearch(.5,-1,1,.01,true,[](double a){return TrialMerit{.1,1+3*a,true};});
 assert(turb.accepted && turb.prelimited && turb.alpha==.25);
 auto recovery=cubicSearch(.5,-1,1,.01,true,[](double a){return TrialMerit{.1,1,a<=.25};});
 assert(recovery.accepted && recovery.prelimited && recovery.alpha==.25);
 std::cout<<"cubic search policy checks passed\n";
}
