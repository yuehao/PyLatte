//
// Created by Yue Hao on 8/12/16.
//

#include "element.h"
#include "bunch.h"
#include "mathfunc.h"

CTPS approx_ps2(const psp& ps){//momentum in s
    return 1.0+ps[dE_];
}
std::vector<double> approx_ps2(const bunch& ps){//momentum in s
    return 1.0+ps[dE_];
}

T exact_ps2(CPS6D<T>& ps){//momentum in s
    if (GLOVAR.exactHamiltonian)
        return (1.0+ps[dE_])*(1.0+ps[dE_])-ps[px_]*ps[px_]-ps[py_]*ps[py_];
    else
        return (1.0+ps[dE_])*(1.0+ps[dE_]);
}

template<typename T>
void PassDrift(const double& L, CPS6D<T>& ps){
    
    
    //  T u=L/(1+ps[dE_]);
    //  ps[x_]+=ps[px_]*u;
    //  ps[y_]+=ps[py_]*u;
    //  ps[ct_]+=(ps[px_]*ps[px_]+ps[py_]*ps[py_])*u/2.0/(1+ps[dE_]);
    
    T u=L/sqrt_ps2(ps);
    ps[x_]+=ps[px_]*u;
    ps[y_]+=ps[py_]*u;
    ps[ct_]+=((1.0+ps[dE_])*u-L);
    
}


template<typename T>
void PassMatrix(CMatrix Map, CPS6D<T>& ps){
    int order=Map.row();
    CPS6D<T> temp;
    for (int i=0;i<order;i++){
        for (int j=0;j<order;j++) temp[i]+=Map(i,j)*ps[j];
    }
    ps=temp;
}
template<typename T>
void SyncRadiate(const double& cur, const double& l,  const T& bx, const T& by, CPS6D<T>& x){
    //	T rl=(1+x[x_]*cur+(x[px_]*x[px_]+x[py_]*x[py_])/2.0/(1+x[dE_])/(1+x[dE_]))*l;
    T u=sqrt_ps2(x);
    T rl=(x[x_]*cur+(1.0+x[dE_])/u)*l;
    //T b2=(Kickreal+cur)*(Kickreal+cur)+Kickimag*Kickimag;
    T absn=sqrt(pow(1.0+x[x_]*cur,2.0)+pow(x[px_]/u,2.0)+pow(x[py_]/u,2.0));
    T ex=x[px_]/u/absn;
    T ey=x[py_]/u/absn;
    T ez=(1+x[x_]*cur)/absn;
    T b2=pow(by*ez,2.0)+pow(bx*ez,2.0)+pow(bx*ey-by*ex,2.0);
    
    x[dE_]-=GLOVAR.Particle_In_Use->cgamma/2.0/M_PI*pow(GLOVAR.CurrentEnergy, 3.0)*pow(1.0+x[dE_],2.0)*b2*rl;
    T newps=sqrt_ps2(x);
    x[px_]=x[px_]/u*newps;
    x[py_]=x[py_]/u*newps;
    
}

template<typename T>
void PassMpole(CPS6D<T>& x, const vector<double>& bn, const vector<double>& an, const double& cur, double l){
    GLOVAR.numberMpoleKick++;
    if (l==0) {
        //x[px_]+=cur*x[dE_];
        x[px_]+=cur*(sqrt_ps2(x)-1);
        x[ct_]+=cur*x[x_];
        l=1;
    }
    else if (cur!=0) {
        //	x[px_]+=(cur*x[dE_]-cur*cur*x[x_])*l;
        x[px_]+=(cur*(sqrt_ps2(x)-1)-cur*cur*x[x_])*l;
        x[ct_]+=cur*l*x[x_];
    }
    if (bn.size()==0 && an.size()==0) return;
    
    
    
    
    T Kickreal, Kickimag, xaccureal, xaccuimag;
    Kickreal=bn[0];
    Kickimag=an[0];
    xaccureal=1;
    xaccuimag=0;
    for (int i=1;i<bn.size();i++){
        T temp=xaccureal*x[x_]-xaccuimag*x[y_];
        xaccuimag=xaccuimag*x[x_]+xaccureal*x[y_];
        xaccureal=temp;
        
        Kickreal+=(xaccureal*bn[i]-xaccuimag*an[i]);
        Kickimag+=(xaccureal*an[i]+xaccuimag*bn[i]);
    }
    
    if (GLOVAR.Radiation) {
        SyncRadiate(cur, l, Kickimag, Kickreal+cur, x);
    }
    x[px_]-=Kickreal*l;
    x[py_]+=Kickimag*l;
    return;
}
template<typename T>
void PassBeamBeam(CPS6D<T>& x, const vector<double>& bs, const vector<double>& betaw, const CParticle& weakp, const double& nstr, const double& cutoff, const double& nstep, const int& i, const double& charge2){
    double binsize=2*cutoff*bs[Z_]/nstep;
    
    
    
    double center=-cutoff*bs[Z_]+(i+0.5)*binsize;
    
    //double CP=center/2.0;
    double ps=bs[X_]*sqrt(1+center*center/4/betaw[0]/betaw[0]);
    if (i==0){
        T CP=(center-x[ct_])/2.0;
        x[x_]+=x[px_]*CP;
        x[y_]+=x[py_]*CP;
    }
    T factor=nstr*binsize*exp(-center*center/2.0/bs[Z_]/bs[Z_])*weakp.cradius/weakp.gamma(GLOVAR.CurrentEnergy)/sqrt(2*M_PI)/bs[Z_];
    T rr=x[x_]*x[x_]+x[y_]*x[y_]; T dpx,dpy;
    if (rr > 0){
        
        dpx=charge2*2.0*factor*x[x_]*(1-exp(-rr/2.0/ps/ps))/rr;
        dpy=charge2*2.0*factor*x[y_]*(1-exp(-rr/2.0/ps/ps))/rr;
    }
    else{
        
        dpx=charge2*2.0*factor*x[x_]/2.0/ps/ps;
        dpy=charge2*2.0*factor*x[y_]/2.0/ps/ps;
    }
    if (GLOVAR.Radiation) SyncRadiate(0,binsize,dpx/binsize,dpy/binsize,x);
    x[px_]+=dpx;
    x[py_]+=dpy;
    x[x_]+=x[px_]*binsize;
    x[y_]+=x[py_]*binsize;
    
    if (i==nstep-1){
        T CP=(center+binsize-x[ct_])/2.0;
        x[x_]-=x[px_]*CP;
        x[y_]-=x[py_]*CP;
    }
    
    
    
}