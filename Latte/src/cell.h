/*
 *  cell.h
 *  ERLtrac
 *
 *  Created by Yue Hao on 4/27/11.
 *  Copyright 2011 Brookhaven National Laboratory. All rights reserved.
 *
 */

#ifndef LATTICE_H
#define LATTICE_H

#include "matrix.h"
#include "constants.h"
#include "psc.h"
#include "twiss4D.h"

#include <vector>
#include <string>

using namespace std;

class CBeamLine;

class CCell
{
public:
	CCell* next;
	double l;   //length
	vector<double> passl;
	vector<double> s;    //location x, y, z
	vector<double> sp;    //location x, y, z
	string name; //Element Name
	CellType eletype; //Element Type from enum CellType
	int nstep; //Integration Step
	int npass;
	/* the optics functions beta alpha gamma dispersion dispersion'  The length should be 2*nstep*/
	vector<twiss4D> twiss;
	vector<string> requiredparameter;
	
	double aperture;
	
	vector<CMatrix> linmap; //Linear Map after this Cell
	vector<vector<double> >coor;
	
	vector<double> pos_error_rms; //RMS Position Error , x ,  y, z
	vector<double> pos_error; //Real Position Error , x, y ,z
	
	vector<double> tilt_error_rms;  //tilt error in x, y
	vector<double> tilt_error;  //real tilt error in x, y
	
	double rot_angle; // Real Rotation angle between x and y with respect to its own axis
	double rot_angle_rms;  //Rotation error rms
	double rot_angle_design; //Design Rotation
	
	CCell();
	virtual ~CCell(){}
	
	
    void SetError();
	
	virtual int Pass(CPS6D<double>& x)=0;
	virtual int Pass(CPS6D<CTPS>& ps)=0;
	virtual int SetP(const string& str, const double& value)=0;
	virtual int SetP(const string& str, const string& svalue)=0;
	virtual int SetP(const string& str, const vector<double>& lvalue)=0;
	virtual int SetP(const string& str, const vector<string>& slvalue)=0;
	virtual void Initialize()=0;
	void SetParameter(const CElement& ele);
	void checkaperture(CPS6D<double>& x);
	
};


class CBeamBeam:public CCell {
private:
	vector<double> beamsize;
	vector<double> betaw;
	CParticle strbeam;
	double nstrong;
	double fcoll;  //Collision frequency to calculate luminosity based on weak beam;
	double cutoff;
	
public:
	CBeamBeam();
	~CBeamBeam(){}
	
	virtual int Pass(CPS6D<double>& x);
	virtual int Pass(CPS6D<CTPS>& ps);
	virtual int SetP(const string& str, const double& value);
	virtual int SetP(const string& str, const string& svalue);
	virtual int SetP(const string& str, const vector<double>& lvalue);
	virtual int SetP(const string& str, const vector<string>& slvalue);
	virtual void Initialize();
	
	
};

class CDrift:public CCell{
public: 
	CDrift();
	~CDrift(){}
	virtual int Pass(CPS6D<double>& x);
	virtual int Pass(CPS6D<CTPS>& ps);
	virtual int SetP(const string& str, const double& value);
	virtual int SetP(const string& str, const string& svalue);
	virtual int SetP(const string& str, const vector<double>& lvalue);
	virtual int SetP(const string& str, const vector<string>& slvalue);
	virtual void Initialize();
	
};

class CCavity:public CCell{
private:
	int autophasingstatus;
public:
	double voltage;
	double phase;
	double freq;
	double harmonic;
	int method;  //1 uniform acceleration tube
	
	vector<double> homfreq;
	vector<double> homvoltage;
	vector<double> extQ;
	vector<double> homphase;
	CCavity();
	~CCavity(){}
	virtual int Pass(CPS6D<double>& x);
	virtual int Pass(CPS6D<CTPS>& ps);
	virtual int SetP(const string& str, const double& value);
	virtual int SetP(const string& str, const string& svalue);
	virtual int SetP(const string& str, const vector<double>& lvalue);
	virtual int SetP(const string& str, const vector<string>& slvalue);
	virtual void Initialize();
	
};

class CLinearMap:public CCell{
public:	
	int dim; 
	CMatrix themap;
	
	CLinearMap();
	~CLinearMap(){}
	
	twiss4D tstart,tend;
	vector<double> chrom;
	
	virtual int Pass(CPS6D<double>& x);
	virtual int Pass(CPS6D<CTPS>& ps);
	virtual int SetP(const string& str, const double& value);
	virtual int SetP(const string& str, const string& svalue);
	virtual int SetP(const string& str, const vector<double>& lvalue);
	virtual int SetP(const string& str, const vector<string>& slvalue);
	virtual void Initialize();
	
};


class CMpole:public CCell{
private:
	vector<double> bstr,astr,design_bstr,design_astr,errrms_bstr,errrms_astr,errran_bstr,errran_astr;
	double bend_angle;
	double curvature;
	int design_order;  //Dipole:  0  Quadrupole: 1 ...
	int max_order;  
	int usekvalue;  //set the bn and an use Kvalue or not, bn=kn/(n-1)!
	vector<double> edgeangle;
	
public:
	CMpole();
	~CMpole(){}
	virtual int Pass(CPS6D<double>& x);
	virtual int Pass(CPS6D<CTPS>& ps);
	virtual int SetP(const string& str, const double& value);
	virtual int SetP(const string& str, const string& svalue);
	virtual int SetP(const string& str, const vector<double>& lvalue);
	virtual int SetP(const string& str, const vector<string>& slvalue);
	virtual void Initialize();
	friend class CBeamLine;
};


class CSplitter:public CCell{
private:
	int npass;
	double espraccept;  //energy spread acceptances
	double refenergy;
	double refbend_angle;
	vector<double> energy;
	vector<string> passname;
	vector<double> arcl;
	vector<double> curvature;
	vector<CBeamLine*> passes;  //Pointer to Energy recovery path
	int whichpass(const double& en);
public:
	CSplitter();
	~CSplitter(){}
	virtual int Pass(CPS6D<double>& x);
	virtual int Pass(CPS6D<CTPS>& ps);
	virtual int SetP(const string& str, const double& value);
	virtual int SetP(const string& str, const string& svalue);
	virtual int SetP(const string& str, const vector<double>& lvalue);
	virtual int SetP(const string& str, const vector<string>& slvalue);
	virtual void Initialize();
	
};

class CMarker:public CCell{
public:
	int isbpm;
	double noise;
	vector<double> position;
	CMarker(){isbpm=0; position.resize(2);}
	~CMarker(){}
	virtual int Pass(CPS6D<double>& x);
	virtual int Pass(CPS6D<CTPS>& ps);
	virtual int SetP(const string& str, const double& value);
	virtual int SetP(const string& str, const string& svalue) {return 0;}
	virtual int SetP(const string& str, const vector<double>& lvalue) ;
	virtual int SetP(const string& str, const vector<string>& slvalue) {return 0;}
	virtual void Initialize();
	
};

class CCorrector:public CCell{
public:
	vector<double> kick;
	CCorrector(){kick.resize(2);}
	~CCorrector(){}
	virtual int Pass(CPS6D<double>& x);
	virtual int Pass(CPS6D<CTPS>& ps);
	virtual int SetP(const string& str, const double& value);
	virtual int SetP(const string& str, const string& svalue) {return 1;}
	virtual int SetP(const string& str, const vector<double>& lvalue);
	virtual int SetP(const string& str, const vector<string>& slvalue) {return 1;}
	virtual void Initialize();
	
};

#endif
