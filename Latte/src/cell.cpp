/*
 *  CCell.cpp
 *  ERLtrac
 *
 *  Created by Yue Hao on 4/27/11.
 *  Copyright 2011 Brookhaven National Laboratory. All rights reserved.
 *
 */

#include "cell.h"
#include "global.h"
#include "mathfunctions.h"
#include "line.h"

#include <cstdlib>
#include <complex>
using namespace std;
extern GLOBAL_VARIABLE GLOVAR;

template<typename T>
T sqrt_ps2(CPS6D<T>& ps){//momentum in s
  if (GLOVAR.exactHamiltonian)
    return sqrt((1.0+ps[dE_])*(1.0+ps[dE_])-ps[px_]*ps[px_]-ps[py_]*ps[py_]);
  else
    return 1.0+ps[dE_];
}
template<typename T>
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
/**********************************CCELL and its sons******************************************/

CCell::CCell(){
	l=0.0;
	s.resize(3);
	nstep=1;
	npass=0;//accumulated pass of this element
	eletype=NOTYPE;
	aperture=0;
	
	pos_error.resize(3);
	pos_error_rms.resize(3);
	tilt_error.resize(2);
	tilt_error_rms.resize(2);
	rot_angle_design=0;
	rot_angle_rms=0;
	rot_angle=0;
		
}
void CCell::SetParameter(const CElement& ele){
	 map<string,CVariable>::const_iterator  it;
	int checkflag=1;
	for (int i=0; i<requiredparameter.size(); i++) {
		checkflag*=ele.varlist.count(requiredparameter[i]);
		if (checkflag==0) {
			cout << "The element "<<ele.name<<"does not have required parameter "<<requiredparameter[i]<<endl;
			exit(-1);
		}
	}
	name=ele.name;
	eletype=ele.ct;
	for ( it = ele.varlist.begin(); it != ele.varlist.end(); it++) {
		switch (it->second.type) {
			case 1:
				SetP(it->first, it->second.value);
				break;
			case 2:
				SetP(it->first, it->second.svalue);
				break;
				
			case 3:
				SetP(it->first, it->second.vlist);
				break;
				
			case 4:
				SetP(it->first, it->second.slist);
				break;
				
			default:
				break;
		}
	}
	
}

void CCell::SetError(){
	pos_error[X_]+=pos_error_rms[X_]*Random(NORMAL,0,1);
	pos_error[Y_]+=pos_error_rms[Y_]*Random(NORMAL,0,1); 
	tilt_error[X_]+=tilt_error_rms[X_]*Random(NORMAL,0,1);
	tilt_error[Y_]+=tilt_error_rms[Y_]*Random(NORMAL,0,1);
	rot_angle=rot_angle_design+rot_angle_rms*Random(NORMAL,0,1);
}

void CCell::checkaperture(CPS6D<double>& x){
	if (this->aperture>0 && pow(x[x_],2.0)+pow(x[y_], 2.0)>pow(aperture, 2.0)) {
		cout << "The beam loss happened at "<<this->name<<", which is "<<this->passl[npass]<<endl;
		exit(-1);
	}
}

/***********Matrix Element*********************/
CLinearMap::CLinearMap(){
	dim=6;
	themap.resize(dim,dim);
	themap.setI();
	chrom.resize(2);
}
int CLinearMap::SetP(const string& str, const double& value){
	if (str=="L") this->l=value;
	else if (str=="R11") this->themap(0,0)=value;
	else if (str=="R12") this->themap(0,1)=value;
	else if (str=="R13") this->themap(0,2)=value;
	else if (str=="R14") this->themap(0,3)=value;
	else if (str=="R15") this->themap(0,4)=value;
	else if (str=="R16") this->themap(0,5)=value;
	
	else if (str=="R21") this->themap(1,0)=value;
	else if (str=="R22") this->themap(1,1)=value;
	else if (str=="R23") this->themap(1,2)=value;
	else if (str=="R24") this->themap(1,3)=value;
	else if (str=="R25") this->themap(1,4)=value;
	else if (str=="R26") this->themap(1,5)=value;
	
	else if (str=="R31") this->themap(2,0)=value;
	else if (str=="R32") this->themap(2,1)=value;
	else if (str=="R33") this->themap(2,2)=value;
	else if (str=="R34") this->themap(2,3)=value;
	else if (str=="R35") this->themap(2,4)=value;
	else if (str=="R36") this->themap(2,5)=value;
	
	else if (str=="R41") this->themap(3,0)=value;
	else if (str=="R42") this->themap(3,1)=value;
	else if (str=="R43") this->themap(3,2)=value;
	else if (str=="R44") this->themap(3,3)=value;
	else if (str=="R45") this->themap(3,4)=value;
	else if (str=="R46") this->themap(3,5)=value;
	
	else if (str=="R51") this->themap(4,0)=value;
	else if (str=="R52") this->themap(4,1)=value;
	else if (str=="R53") this->themap(4,2)=value;
	else if (str=="R54") this->themap(4,3)=value;
	else if (str=="R55") this->themap(4,4)=value;
	else if (str=="R56") this->themap(4,5)=value;
	
	else if (str=="R61") this->themap(5,0)=value;
	else if (str=="R62") this->themap(5,1)=value;
	else if (str=="R63") this->themap(5,2)=value;
	else if (str=="R64") this->themap(5,3)=value;
	else if (str=="R65") this->themap(5,4)=value;
	else if (str=="R66") this->themap(5,5)=value;
	
	else if (str=="BETX1") this->tstart.beta[X_]=value;
	else if (str=="BETX2") this->tend.beta[X_]=value;
	else if (str=="BETY1") this->tstart.beta[Y_]=value;
	else if (str=="BETY2") this->tend.beta[Y_]=value;
	
	else if (str=="ALFX1") this->tstart.alpha[X_]=value;
	else if (str=="ALFX2") this->tend.alpha[X_]=value;
	else if (str=="ALFY1") this->tstart.alpha[Y_]=value;
	else if (str=="ALFY2") this->tend.alpha[Y_]=value;
	
	else if (str=="DPX1") this->tstart.dp[X_]=value;
	else if (str=="DPX2") this->tend.dp[X_]=value;
	else if (str=="DPY1") this->tstart.dp[Y_]=value;
	else if (str=="DPY2") this->tend.dp[Y_]=value;
	
	else if (str=="DX1") this->tstart.d[X_]=value;
	else if (str=="DX2") this->tend.d[X_]=value;
	else if (str=="DY1") this->tstart.d[Y_]=value;
	else if (str=="DY2") this->tend.d[Y_]=value;
	
	else if (str=="PHIX1") this->tstart.tune[X_]=value;
	else if (str=="PHIX2") this->tend.tune[X_]=value;
	else if (str=="PHIY1") this->tstart.tune[Y_]=value;
	else if (str=="PHIY2") this->tend.tune[Y_]=value;
	
	else if (str=="CHROMX") this->chrom[X_]=value;
	else if (str=="CHROMY") this->chrom[Y_]=value;
	else {
		return 0;
	}
	return 1;
}
int CLinearMap::SetP(const string& str, const string& svalue) {return 0;}
int CLinearMap::SetP(const string& str, const vector<double>& lvalue) {return 0;}
int CLinearMap::SetP(const string& str, const vector<string>& slvalue) {return 0;}
void CLinearMap::Initialize(){
	this->tstart.gamma[X_]=(pow(this->tstart.alpha[X_],2.0)+1)/this->tstart.beta[X_];
	this->tstart.gamma[Y_]=(pow(this->tstart.alpha[Y_],2.0)+1)/this->tstart.beta[Y_];
	this->tend.gamma[X_]=(pow(this->tend.alpha[X_],2.0)+1)/this->tend.beta[X_];
	this->tend.gamma[Y_]=(pow(this->tend.alpha[Y_],2.0)+1)/this->tend.beta[Y_];
	tstart.inimatrix();
	tend.inimatrix();
	CMatrix t=transfermap(tstart, tend);
	themap*=t;
}

int CLinearMap::Pass(CPS6D<double>& x){
	PassMatrix(this->themap, x);
	this->passl.push_back(this->l);
	this->coor.push_back(x.getps());
	return 1;
}
int CLinearMap::Pass(CPS6D<CTPS>& x){
	PassMatrix(this->themap, x);
	this->passl.push_back(this->l);
	this->x.push_back(x.getps());
	this->x.push_back(x.getmap());
	return 1;
}
/***********CDrift Element*********************/
CDrift::CDrift(){
	requiredparameter.push_back("L");
}

int CDrift::Pass(CPS6D<double>& x){
	PassDrift(l, x);
	this->passl.push_back(this->l);
	this->coor.push_back(x.getps());
	return 1;
}
int CDrift::Pass(CPS6D<CTPS>& x){
	PassDrift(l, x);
	this->passl.push_back(this->l);
	this->x.push_back(x.getps());
	this->x.push_back(x.getmap());
	return 1;
}
int CDrift::SetP(const string& str, const double& value){
	if (str=="L") this->l=value;
	else if (str=="APERTURE") this->aperture=value;
	else {
		return 0;
	}
	return 1;
}
int CDrift::SetP(const string& str, const string& svalue) {return 0;}
int CDrift::SetP(const string& str, const vector<double>& lvalue) {return 0;}
int CDrift::SetP(const string& str, const vector<string>& slvalue) {return 0;}
void CDrift::Initialize(){
	if (this->aperture==0) {
		this->aperture=GLOVAR.Aperture;
	}
}
/***********CBeamBeam Element*********************/
CBeamBeam::CBeamBeam(){
	this->beamsize.resize(3);
	this->betaw.resize(2);
	this->nstrong=0;
	this->fcoll=0;
	this->cutoff=5.0;
	requiredparameter.push_back("NSTRONG");
	requiredparameter.push_back("FREQCOLL");
	requiredparameter.push_back("BETAW");
	requiredparameter.push_back("BEAMSIZE");
	
}



int CBeamBeam::Pass(CPS6D<double>& x){
  for (int i=0;i<this->nstep;i++){
    PassBeamBeam(x, this->beamsize, this->betaw, *(GLOVAR.Particle_In_Use), this->nstrong, this->cutoff, this->nstep, i,this->strbeam.charge*GLOVAR.Particle_In_Use->charge );
    this->passl.push_back(this->l/this->nstep);
    this->x.push_back(x.getps());
    
  }
  return 1;
}
int CBeamBeam::Pass(CPS6D<CTPS>& x){
  for (int i=0;i<this->nstep;i++){
    PassBeamBeam(x, this->beamsize, this->betaw, *(GLOVAR.Particle_In_Use), this->nstrong, this->cutoff, this->nstep, i, this->strbeam.charge*GLOVAR.Particle_In_Use->charge );
    this->passl.push_back(this->l/this->nstep);
    this->x.push_back(x.getps());
    this->x.push_back(x.getmap());
  }
  return 1;
}

int CBeamBeam::SetP(const string& str, const double& value){
	if (str=="NSTEP") this->nstep=value;
	else if (str=="NSTRONG") this->nstrong=value;
	else if (str=="FREQCOLL") this->fcoll=value;
	else if (str=="TILTERROR") this->rot_angle_rms=value;
	else if (str=="APERTURE") this->aperture=value;
	else return 0;
	return 1;
}
int CBeamBeam::SetP(const string& str, const string& svalue) {return 0;}
int CBeamBeam::SetP(const string& str, const vector<double>& lvalue) {
	if (str=="BEAMSIZE" && lvalue.size()==3) this->beamsize=lvalue;
	else if (str=="BETAW" && lvalue.size()==2) this->betaw=lvalue;
	else return 0;
	return 1;
}
int CBeamBeam::SetP(const string& str, const vector<string>& slvalue) {return 0;}
void CBeamBeam::Initialize(){
	SetError();
	if (this->aperture==0) {
		this->aperture=GLOVAR.Aperture;
	}
}
/*******************CCAVITY*******************************/
CCavity::CCavity(){
	this->voltage=0;
	this->freq=0;
	this->harmonic=0;
	this->phase=0;
	this->method=1;
	requiredparameter.push_back("L");
	requiredparameter.push_back("FREQ");
	requiredparameter.push_back("VOLT");
	if (GLOVAR.CavityAutoPhasing) autophasingstatus=1; else autophasingstatus=0;
	
}
int CCavity::SetP(const string& str, const double& value){
	if (str=="VOLT") this->voltage=value;
	else if (str=="FREQ") this->freq=value;
	else if (str=="HARMON") this->harmonic=value;
	else if (str=="PHASE") this->phase=value;
	else if (str=="L") this->l=value;
	else if (str=="APERTURE") this->aperture=value;
	else if (str=="NSTEP") this->nstep=value;
	else return 0;
	return 1;
}
int CCavity::SetP(const string& str, const string& svalue) {return 0;}
int CCavity::SetP(const string& str, const vector<double>& lvalue) {
	int len=lvalue.size();
	if (str=="DISPLACEMENT" && (len==2 || len==3)) {
		this->pos_error_rms=lvalue;
		this->pos_error_rms.resize(3);
	}
	else return 0;
	return 1;
}
int CCavity::SetP(const string& str, const vector<string>& slvalue) {return 0;}
void CCavity::Initialize(){
	SetError();
	if (this->aperture==0) {
		this->aperture=GLOVAR.Aperture;
	}
}

int CCavity::Pass(CPS6D<double>& x){
	if (autophasingstatus) {
		this->phase-=(Two_Pi*freq*(this->passl.back()-this->l/2.0));
		autophasingstatus=0;
	}
	if (method==1) {
		double oldde=x[dE_];
		double gmogmp=GLOVAR.CurrentEnergy/voltage*this->l;
		if (GLOVAR.IdealERLCavity!=0) {
			x[dE_]+=GLOVAR.IdealERLCavity*voltage/(GLOVAR.CurrentEnergy)*sin(Two_Pi*freq*x[ct_]+phase);
		}
		else {
			x[dE_]+=voltage/(GLOVAR.CurrentEnergy)*sin(Two_Pi*freq*(x[ct_]+this->passl.back()-this->l/2.0)+phase);
		}
		double lnfactor=log((1.0+x[dE_])/(1.0+oldde));
		double postemp=(1-lnfactor/2.0)*x[x_]+gmogmp*lnfactor*x[px_];
		x[px_]=-lnfactor/gmogmp/4.0*x[x_]+(1+lnfactor/2.0)*x[px_];
		x[x_]=postemp;
		
		postemp=(1-lnfactor/2.0)*x[y_]+gmogmp*lnfactor*x[py_];
		x[py_]=-lnfactor/gmogmp/4.0*x[y_]+(1+lnfactor/2.0)*x[py_];
		x[y_]=postemp;
		
		
	}
	return 1;
}
int CCavity::Pass(CPS6D<CTPS>& x){
	if (autophasingstatus) {
		this->phase-=(Two_Pi*freq*(this->passl.back()-this->l/2.0));
		autophasingstatus=0;
	}
	if (method==1) {
		CTPS oldde=x[dE_];
		double gmogmp=GLOVAR.CurrentEnergy/voltage*this->l;
		if (GLOVAR.IdealERLCavity!=0) {
			x[dE_]+=GLOVAR.IdealERLCavity*voltage/(GLOVAR.CurrentEnergy)*sin(Two_Pi*freq*x[ct_]+phase);
		}
		else {
			x[dE_]+=voltage/(GLOVAR.CurrentEnergy)*sin(Two_Pi*freq*(x[ct_]+this->passl.back()-this->l/2.0)+phase);
		}
		CTPS lnfactor=log((1.0+x[dE_])/(1.0+oldde));
		CTPS postemp=(1-lnfactor/2.0)*x[x_]+gmogmp*lnfactor*x[px_];
		x[px_]=-lnfactor/gmogmp/4.0*x[x_]+(1+lnfactor/2.0)*x[px_];
		x[x_]=postemp;
		
		postemp=(1-lnfactor/2.0)*x[y_]+gmogmp*lnfactor*x[py_];
		x[py_]=-lnfactor/gmogmp/4.0*x[y_]+(1+lnfactor/2.0)*x[py_];
		x[y_]=postemp;//*/
		
		
	}
	return 1;
}

/********************CMpole*******************************/
CMpole::CMpole(){
	this->bend_angle=0;
	this->curvature=0;
	this->design_order=0;
	this->max_order=0;
	usekvalue=0;
	design_bstr.resize(1);
	design_astr.resize(1);
	errrms_bstr.resize(1);
	errrms_astr.resize(1);
	edgeangle.resize(2);
	requiredparameter.push_back("L");
	
}

int CMpole::SetP(const string& str, const double& value) {
	if (str=="L") this->l=value;
	else if (str=="ANGLE") {
		this->bend_angle=value;
		this->design_order=0;
	}
	else if (str=="TILT") this->rot_angle_design=value;
	else if (str=="TILTERROR") this->rot_angle_rms=value;
	else if (str=="DESIGN_ORDER") this->design_order=(int)value;
	else if (str=="MAX_ORDER") this->max_order=(int)value;
	else if (str=="APERTURE") this->aperture=value;
	else if (str=="NSTEP") this->nstep=value;
	else if (str=="B0X") { //Bx over BRHO
		if (design_bstr.size()<1) {
			design_bstr.resize(1);
			max_order=0;
			design_order=0;
			design_astr.resize(1);
			errrms_bstr.resize(1);
			errrms_astr.resize(1);
		}
		design_bstr[0]=value;
	}
	else if (str=="B0Y") { //By over BRHO
		if (design_bstr.size()<1) {
			design_bstr.resize(1);
			max_order=0;
			design_order=0;
			design_astr.resize(1);
			errrms_bstr.resize(1);
			errrms_astr.resize(1);
		}
		design_astr[0]=value;
	}
	else if (str=="K1") {
		if (design_bstr.size()<2) {
			design_bstr.resize(2);
			max_order=1;
			design_order=1;
			design_astr.resize(2);
			errrms_bstr.resize(2);
			errrms_astr.resize(2);
		}
		design_bstr[1]=value;
	}
	else if (str=="K1S") {
		if (design_astr.size()<2) {
			design_astr.resize(2);
			
			max_order=1;
			design_order=1;
			design_bstr.resize(2);
			errrms_bstr.resize(2);
			errrms_astr.resize(2);
		}
		design_astr[1]=value;
	}
	else if (str=="K2") {
		if (design_bstr.size()<3) {
			design_bstr.resize(3);
			
			max_order=2;
			design_order=2;
			design_astr.resize(3);
			errrms_bstr.resize(3);
			errrms_astr.resize(3);
		}
		design_astr[1]=value;
	}
	else if (str=="K2S") {
		if (design_astr.size()<3) {
			design_astr.resize(3);
			
			max_order=2;
			design_order=2;
			design_bstr.resize(3);
			errrms_bstr.resize(3);
			errrms_astr.resize(3);
		}
		design_astr[2]=value/2.0;
	}
	else if (str=="K3") {
		if (design_bstr.size()<4) {
			design_bstr.resize(4);
			
			max_order=3;
			design_order=3;
			design_astr.resize(4);
			errrms_bstr.resize(4);
			errrms_astr.resize(4);
		}
		design_bstr[3]=value/6.0;
	}
	else if (str=="K3S") {
		if (design_astr.size()<4) {
			design_astr.resize(4);
			
			max_order=3;
			design_order=3;
			design_bstr.resize(4);
			errrms_bstr.resize(4);
			errrms_astr.resize(4);
		}
		design_astr[3]=value/6.0;
	}
	else if (str=="E1") edgeangle[0]=value;
	else if (str=="E2") edgeangle[1]=value;
	else if (str=="USEKVALUE") usekvalue=(int)value;
	else return 0;
	
	return 1;
	
}
int CMpole::SetP(const string& str, const string& svalue){return 0;}
int CMpole::SetP(const string& str, const vector<double>& lvalue){
	int len=lvalue.size();
	if (str=="DISPLACEMENT_RMS" && (len==2 || len==3)) {
		this->pos_error_rms=lvalue;
		this->pos_error_rms.resize(3);
	}
	else if (str=="DISPLACEMENT" && (len==2 || len==3)) {
		this->pos_error=lvalue;
		this->pos_error.resize(3);
	}
	else if (str=="ANGLE_DISPLACEMENT_RMS" && (len==2)) this->tilt_error_rms=lvalue;
	else if (str=="BN" || str=="KNL") {
	  this->design_bstr=lvalue; 
	  if (len >this->max_order+1) this->max_order=len-1;
	}
	else if (str=="AN" || str=="KSL") {
	  this->design_astr=lvalue; 
	  if (len >this->max_order+1) this->max_order=len-1;
	}
	else if (str=="BERROR" || str=="KNLERROR") {
	  this->errrms_bstr=lvalue;
	  if (len >this->max_order+1) this->max_order=len-1;
	}
	else if (str=="AERROR" || str=="KSLERROR") {
	  this->errrms_astr=lvalue;
	  if (len >this->max_order+1) this->max_order=len-1;
	}
	else return 0;
	//if (len >this->max_order+1) this->max_order=len-1;
	//this->design_bstr.resize(max_order+1);
	//this->design_astr.resize(max_order+1);
	//this->errrms_bstr.resize(max_order+1);
	//this->errrms_astr.resize(max_order+1);
	
	return 1;
}	
int CMpole::SetP(const string& str, const vector<string>& slvalue){return 0;}


int CMpole::Pass(CPS6D<double>& x){
	if (this->l!=0) 
		for (int i=0;i<this->nstep;i++)
			for (int j=0;j<GLOVAR.SI_index.size();j++){
				if (j%2==0) PassDrift(this->l/nstep*GLOVAR.SI_index[j],x);
				else PassMpole(x, this->bstr, this->astr, this->curvature,this->l/nstep*GLOVAR.SI_index[j]);
			}
	else  PassMpole(x, this->bstr, this->astr, this->curvature,0);
	return 1;
}
int CMpole::Pass(CPS6D<CTPS>& ps){
	if (this->l!=0) 
		for (int i=0;i<this->nstep;i++)
			for (int j=0;j<GLOVAR.SI_index.size();j++){
				if (j%2==0) PassDrift(this->l/nstep*GLOVAR.SI_index[j],ps);
				else PassMpole(ps, this->bstr, this->astr, this->curvature,this->l/nstep*GLOVAR.SI_index[j]);
			}
	else  PassMpole(ps, this->bstr, this->astr, this->curvature,0.0);
	return 1;
}
void CMpole::Initialize(){
	if(this->l!=0) this->curvature=this->bend_angle/this->l;
	else curvature=this->bend_angle;
	if (this->nstep<GLOVAR.nstep) {
		this->nstep=GLOVAR.nstep;
	}
	
	errran_bstr.resize(max_order+1);
	errran_astr.resize(max_order+1);
	design_bstr.resize(max_order+1);
	design_astr.resize(max_order+1);
	errrms_bstr.resize(max_order+1);
	errrms_astr.resize(max_order+1);
	bstr.resize(max_order+1);
	astr.resize(max_order+1);
	if (usekvalue) {
		for (int i=2; i<design_bstr.size(); i++) {
			design_bstr[i]/=factorial(i);
			design_astr[i]/=factorial(i);
			errrms_bstr[i]/=factorial(i);
			errrms_astr[i]/=factorial(i);
			
		}
	}
	for (int i=0; i<=max_order; i++) {
		errran_bstr[i]=errrms_bstr[i]*Random(NORMAL, 0.0, 1.0);
		errran_astr[i]=errrms_astr[i]*Random(NORMAL, 0.0, 1.0);
		bstr[i]=design_bstr[i]+errran_bstr[i];
		astr[i]=design_astr[i]+errran_astr[i];
	}
	
	SetError();
	if (this->aperture==0) {
		this->aperture=GLOVAR.Aperture;
	}
}
/*********************Spreader and Merger**************/
CSplitter::CSplitter()
{
	npass=0;
	espraccept=0.2;
	refenergy=GLOVAR.RefEnergy;
	refbend_angle=0;
}

int CSplitter::SetP(const string& str, const double& value){
	if (str=="L") this->l=value;
	else if (str=="REFENERGY") this->refenergy=value;
	else if (str=="ANGLE") {
		this->refbend_angle=value;
	}
	else if (str=="APERTURE") this->aperture=value;
	else if (str=="TILT") this->rot_angle_design=value;
	else return 0;
	
	return 1;
}
int CSplitter::SetP(const string& str, const string& svalue){return 0;}
int CSplitter::SetP(const string& str, const vector<double>& lvalue){
	int len=lvalue.size();
	if (str=="ENERGIES") energy=lvalue;
	else if (str=="DISPLACEMENT" && (len==2 || len==3)) {
		this->pos_error_rms=lvalue;
		this->pos_error_rms.resize(3);
	}
	else return 0;
	return 1;
}
int CSplitter::SetP(const string& str, const vector<string>& slvalue){
	if (str=="PATHES") passname=slvalue;
	else return 0;
	return 1;
}
int CSplitter::Pass(CPS6D<double>& x){
	vector<double> bn;
	bn.resize(1);
	int np=whichpass(GLOVAR.CurrentEnergy*(1+x[dE_]));
	if (np<0) {
		cout << "The energy does not fit the Spreader/Combiner:" <<this->name <<". Particle lost."<<endl;
		return 0;
	}
	double ratio;
	if (this->eletype==SPREADER) {
		ratio=GLOVAR.CurrentEnergy/energy[np];
		x[dE_]=(1+x[dE_])/ratio-1.0;
		x[px_]/=ratio;
		x[py_]/=ratio;
		PassMpole(x, bn, bn, curvature[np], arcl[np]);
	}
	else if (this->eletype==MERGER){
		PassMpole(x, bn, bn, curvature[np], arcl[np]);
		ratio=energy[np]/GLOVAR.RefEnergy;
		x[dE_]=(1+x[dE_])/ratio-1.0;
		x[px_]/=ratio;
		x[py_]/=ratio;
	} 
	
	//if (this->eletype==SPREADER)  return passes[np]->Pass(x,"DUMMY");
	return 1;
}
int CSplitter::Pass(CPS6D<CTPS>& x){
	vector<double> bn;
	bn.resize(1);
	int np=whichpass(GLOVAR.CurrentEnergy*(1+x[dE_].cst()));
	if (np<0) {
		cout << "The energy does not fit the Spreader/Combiner:" <<this->name <<". Particle lost."<<endl;
		return 0;
	}
	double ratio;
	if (this->eletype==SPREADER) {
		ratio=GLOVAR.CurrentEnergy/energy[np];
		x[dE_]=(1+x[dE_])/ratio-1.0;
		x[px_]/=ratio;
		x[py_]/=ratio;
		PassMpole(x, bn, bn, curvature[np], arcl[np]);
	}
	else if (this->eletype==MERGER){
		PassMpole(x, bn, bn, curvature[np], arcl[np]);
		ratio=energy[np]/GLOVAR.RefEnergy;
		x[dE_]=(1+x[dE_])/ratio-1.0;
		x[px_]/=ratio;
		x[py_]/=ratio;
	} 
	
	//if (this->eletype==SPREADER)  return passes[np]->Pass(x);
	return 1;
}
void CSplitter::Initialize(){
	if (energy.size()!=passname.size()) {
		cout << "The energylist and passlist do not have same size in Spreader/Merger:" <<this->name<<endl;
		exit(-1);
	}
	npass=energy.size();
	curvature.resize(npass);
	arcl.resize(npass);
	passes.resize(npass);
	for (int i=0;i<npass;i++){
		if (l>0) curvature[i]=tan(refbend_angle)/l*refenergy/energy[i]; else curvature[i]=refbend_angle*refenergy/energy[i]; 
		if (refbend_angle!=0 && l>0) arcl[i]=atan(tan(refbend_angle)*refenergy/energy[i])/curvature[i]; else arcl[i]=l;
		passes[i]=new CBeamLine;
		passes[i]->setbeamline(passname[i]);
	}
	SetError();
	if (this->aperture==0) {
		this->aperture=GLOVAR.Aperture;
	}
}
int CSplitter::whichpass(const double& en){
	
	if (en<=0) return -1;
	if (npass==1) {
		return 0;
	}
	for (int i=0; i<npass; i++) {
		if (i==0 && en<(energy[i]+energy[i+1])/2.0) return 0;	
		else if (i==npass-1) return npass-1;
		else if (en<(energy[i]+energy[i+1])/2.0 && en>(energy[i]+energy[i-1])/2.0) return i;
	}
	return -1;
}

/****************************CMARKER************************************************/

int CMarker::Pass(CPS6D<double>& x){
	if (isbpm) {
		position[X_]=x[x_];
		position[Y_]=x[y_];
	}
	if (name=="EXIT") {
		return 0;
	}
	if (name=="THEERLUSER" && GLOVAR.IdealERLCavity==1) {
		GLOVAR.IdealERLCavity=-1;
	}
	return 1;
}
int CMarker::Pass(CPS6D<CTPS>& ps){
	{
		if (isbpm) {
			position[X_]=ps[x_].cst();
			position[Y_]=ps[y_].cst();
		}
		if (name=="EXIT") {
			return 0;
		}
		if (name=="THEERLUSER" && GLOVAR.IdealERLCavity==1) {
			GLOVAR.IdealERLCavity=-1;
		}
		return 1;
	}
	
}

int CMarker::SetP(const string& str, const vector<double>& lvalue){
	int len=lvalue.size();
	if (str=="DISPLACEMENT" && (len==2 || len==3)) {
		this->pos_error_rms=lvalue;
		this->pos_error_rms.resize(3);
	}
	else return 0;
	return 1;
}

int CMarker::SetP(const string& str, const double& value){
	if (str=="TILT") this->rot_angle_design=value;
	else if (str=="NOISE") this->noise=value;
	else if (str=="APERTURE") this->aperture=value;
	else return 0;
	
	return 1;
}

void CMarker::Initialize(){
	SetError();
	if (this->aperture==0) {
		this->aperture=GLOVAR.Aperture;
	}
}
/****************************CCorrector************************************************/

int CCorrector::Pass(CPS6D<double>& x){
	PassDrift(this->l/2.0, x);
	x[px_]+=this->kick[X_];
	x[py_]+=this->kick[Y_];
	PassDrift(this->l/2.0, x);
	return 1;
}
int CCorrector::Pass(CPS6D<CTPS>& x){
	PassDrift(this->l/2.0, x);
	x[px_]+=this->kick[X_];
	x[py_]+=this->kick[Y_];
	PassDrift(this->l/2.0, x);
	return 1;
}
int CCorrector::SetP(const string& str, const double& value){
	if (str=="L") this->l=value;
	else return 0;
	
	return 1;
}

int CCorrector::SetP(const string& str, const vector<double>& lvalue){
	int len=lvalue.size();
	if (str=="KICKANGLE" && (len==2)) {
		this->kick=lvalue;
	}
	else return 0;
	return 1;
}
void CCorrector::Initialize(){
	if (this->aperture==0) {
		this->aperture=GLOVAR.Aperture;
	}
}

