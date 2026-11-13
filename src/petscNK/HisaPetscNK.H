#ifndef HisaPetscNK_H
#define HisaPetscNK_H
// Passive primal NK core. The adapter owns the physical residual and state.
#include <petscksp.h>
#include <vector>
#include <functional>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <algorithm>
namespace hisaNK {
using Array=std::vector<double>;
inline void check(PetscErrorCode e) { if(e) throw std::runtime_error("PETSc NK error "+std::to_string(e)); }
inline double dot(const Array&a,const Array&b) { double v=0,g=0; for(size_t i=0;i<a.size();++i)v+=a[i]*b[i]; MPI_Allreduce(&v,&g,1,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD); return g; }
inline double norm(const Array&a) {return std::sqrt(dot(a,a));}
struct Options {int subspace=60,maxLinear=60,overlap=1,fill=2,lag=20; double rtol=.3,stagnationTol=.99; bool ew=true;};
struct Result {bool accepted=false; double alpha=0,linear=1; int iterations=0,reason=0,retries=0; bool stagnated=false; double firstLinear=0;};
class Solver {
    Options opt_; Array weights_,base_,res_; Mat A_=nullptr,P_=nullptr; KSP ksp_=nullptr;
    std::function<Array(const Array&)> residual_;
    std::function<bool(const Array&)> admissible_;
    std::function<Mat()> matrix_;
    int count_=0; double previousNorm_=0,previousEta_=.3;
    static PetscErrorCode multiply(Mat A,Vec x,Vec y) {
        void* c=nullptr; PetscErrorCode e=MatShellGetContext(A,&c); if(e)return e;
        Solver& s=*static_cast<Solver*>(c); const PetscScalar* xp=nullptr; PetscScalar* yp=nullptr;
        try {
            check(VecGetArrayRead(x,&xp)); Array v(s.base_.size()); for(size_t i=0;i<v.size();++i)v[i]=PetscRealPart(xp[i]);
            check(VecRestoreArrayRead(x,&xp)); double nv=norm(v); Array out(v.size(),0);
            if(nv>0) {double h=std::sqrt(std::numeric_limits<double>::epsilon())*std::sqrt(1+norm(s.base_))/nv;
                Array trial=s.base_; for(size_t i=0;i<v.size();++i)trial[i]+=h*v[i]; out=s.evaluate(trial);
                for(size_t i=0;i<v.size();++i)out[i]=(out[i]-s.res_[i])/h;}
            check(VecGetArray(y,&yp));for(size_t i=0;i<out.size();++i)yp[i]=out[i];check(VecRestoreArray(y,&yp));
        } catch(...) {return PETSC_ERR_LIB;} return 0;
    }
    Array evaluate(const Array& z) {Array r=residual_(z);for(size_t i=0;i<r.size();++i)r[i]*=weights_[i];return r;}
    void resetPC() {if(ksp_)KSPDestroy(&ksp_);if(P_)MatDestroy(&P_);}
public:
    Solver(const Options&o,const Array&w,std::function<Array(const Array&)>r,std::function<bool(const Array&)>a,std::function<Mat()>p)
    :opt_(o),weights_(w),residual_(r),admissible_(a),matrix_(p) {
        if(o.subspace<1||o.maxLinear<1||o.lag<1||o.overlap<0||o.fill<0||o.rtol<=0||o.rtol>=1||o.stagnationTol<=0||o.stagnationTol>1)throw std::runtime_error("Invalid NK options");
        check(MatCreateShell(PETSC_COMM_WORLD,w.size(),w.size(),PETSC_DETERMINE,PETSC_DETERMINE,this,&A_));
        check(MatShellSetOperation(A_,MATOP_MULT,(void(*)(void))multiply));
    }
    ~Solver(){resetPC();if(A_)MatDestroy(&A_);}
    void restart(){resetPC();count_=0;previousNorm_=0;}
    Result step(Array& z, bool retry=false) {
        const double savedNorm=previousNorm_,savedEta=previousEta_;
        Result result; base_=z;res_=evaluate(z);double rn=norm(res_);
        if(!std::isfinite(rn)||rn==0)return result;
        if(!ksp_||count_%opt_.lag==0){resetPC();P_=matrix_();Vec w=nullptr;check(VecCreateMPIWithArray(PETSC_COMM_WORLD,1,weights_.size(),PETSC_DECIDE,weights_.data(),&w));check(MatDiagonalScale(P_,w,nullptr));check(VecDestroy(&w));
            check(KSPCreate(PETSC_COMM_WORLD,&ksp_));check(KSPSetType(ksp_,KSPGMRES));check(KSPGMRESSetRestart(ksp_,opt_.subspace));check(KSPSetPCSide(ksp_,PC_RIGHT));check(KSPSetNormType(ksp_,KSP_NORM_UNPRECONDITIONED));check(KSPSetOperators(ksp_,A_,P_));
            PC pc;check(KSPGetPC(ksp_,&pc));check(PCSetType(pc,PCASM));check(PCASMSetOverlap(pc,opt_.overlap));check(KSPSetUp(ksp_));PetscInt n;KSP* subs;check(PCASMGetSubKSP(pc,&n,nullptr,&subs));
            for(PetscInt i=0;i<n;++i){PC sub;check(KSPSetType(subs[i],KSPPREONLY));check(KSPGetPC(subs[i],&sub));check(PCSetType(sub,PCILU));check(PCFactorSetLevels(sub,opt_.fill));check(PCFactorSetMatOrderingType(sub,MATORDERINGRCM));}
        }
        double eta=opt_.rtol;if(opt_.ew&&previousNorm_>0){double power=(1+std::sqrt(5.))/2;eta=std::pow(rn/previousNorm_,power);double safe=std::pow(previousEta_,power);if(safe>.1)eta=std::max(eta,safe);eta=std::min(eta,.8);}
        previousNorm_=rn;previousEta_=eta;check(KSPSetTolerances(ksp_,eta,1e-12,10,opt_.maxLinear));
        Vec rhs=nullptr,dx=nullptr,jd=nullptr;Array b=res_;for(double&v:b)v=-v;
        check(VecCreateMPIWithArray(PETSC_COMM_WORLD,1,b.size(),PETSC_DECIDE,b.data(),&rhs));check(VecDuplicate(rhs,&dx));check(VecDuplicate(rhs,&jd));check(KSPSolve(ksp_,rhs,dx));
        PetscInt its;KSPConvergedReason reason;check(KSPGetIterationNumber(ksp_,&its));check(KSPGetConvergedReason(ksp_,&reason));result.iterations=its;result.reason=reason;
        check(MatMult(A_,dx,jd));const PetscScalar *p;Array d(z.size()),j(z.size());check(VecGetArrayRead(dx,&p));for(size_t i=0;i<z.size();++i)d[i]=PetscRealPart(p[i]);check(VecRestoreArrayRead(dx,&p));check(VecGetArrayRead(jd,&p));for(size_t i=0;i<z.size();++i)j[i]=PetscRealPart(p[i]);check(VecRestoreArrayRead(jd,&p));
        Array defect=j;for(size_t i=0;i<z.size();++i)defect[i]+=res_[i];result.linear=norm(defect)/rn;
        VecDestroy(&jd);VecDestroy(&dx);VecDestroy(&rhs);++count_;
        // Rebuild once at the same base state before accepting a useless direction.
        if(!std::isfinite(result.linear)||result.linear>=opt_.stagnationTol) {
            residual_(base_);
            if(!retry) {
                resetPC(); --count_;
                previousNorm_=savedNorm;previousEta_=savedEta;
                Result second=step(z,true);
                second.retries=1;second.firstLinear=result.linear;
                second.iterations+=result.iterations;
                return second;
            }
            result.stagnated=true;
            return result;
        }
        // An iteration limit may still yield a useful direction; breakdown may not.
        if(reason<0&&reason!=KSP_DIVERGED_ITS){residual_(base_);return result;}
        double f0=.5*rn*rn,slope=-std::abs(dot(res_,j)),alpha=1,prevA=0,prevF=0;if(slope==0)slope=-1;
        for(int trialNo=0;trialNo<30;++trialNo){Array trial=base_;for(size_t i=0;i<z.size();++i)trial[i]+=alpha*d[i];
            if(!admissible_(trial)){alpha*=.5;prevA=0;continue;}
            Array r=evaluate(trial);double f=.5*dot(r,r);if(!std::isfinite(f)){alpha*=.5;prevA=0;continue;}
            // Do not accept a worsening minimum step in the production adapter.
            if(f<=f0+.01*alpha*slope){z=trial;result.accepted=true;result.alpha=alpha;return result;}
            double candidate=alpha/2;
            if(prevA==0){double den=2*(f-f0-alpha*slope);if(den!=0)candidate=-slope*alpha*alpha/den;}
            else {double t1=f-f0-alpha*slope,t2=prevF-f0-prevA*slope;double a=(t1/(alpha*alpha)-t2/(prevA*prevA))/(alpha-prevA);double bb=(-prevA*t1/(alpha*alpha)+alpha*t2/(prevA*prevA))/(alpha-prevA);candidate=a!=0?(-bb+std::sqrt(std::max(bb*bb-3*a*slope,0.)))/(3*a):(bb!=0?-slope/(2*bb):alpha/2);}
            prevA=alpha;prevF=f;if(!std::isfinite(candidate))candidate=alpha/2;alpha=std::min(.5*alpha,std::max(.1*alpha,candidate));if(alpha<1e-8)break;
        }
        residual_(base_);return result;
    }
};
}
#endif
