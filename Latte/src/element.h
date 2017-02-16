//
// Created by Yue Hao on 8/12/16.
//

#ifndef LATTE_ELEMENT_H
#define LATTE_ELEMENT_H
#include <armadillo>
#include <vector>
#include "constants.h"
#include "twiss4D.h"
#include <string>
#include "psp.h"
#include "tpsa.h"

class element {
    element* next;
    double l;   //length
    std::vector<double> passl;
    std::string name; //Element Name
    CellType eletype; //Element Type from enum CellType
    int nstep; //Integration Step
    int npass;
    /* the optics functions beta alpha gamma dispersion dispersion'  The length should be 2*nstep*/
    std::vector<twiss4D> twiss;
    static std::vector<std::string> required_parameter;
    
    double aperture;
    
    std::vector<double> pos_error_rms; //RMS Position Error , x ,  y, z
    std::vector<double> pos_error; //Real Position Error , x, y ,z
    
    std::vector<double> tilt_error_rms;  //tilt error in x, y
    std::vector<double> tilt_error;  //real tilt error in x, y
    
    double rot_angle; // Real Rotation angle between x and y with respect to its own axis
    double rot_angle_rms;  //Rotation error rms
    double rot_angle_design; //Design Rotation
    
    element();
    virtual ~element(){}
    
    
    void SetError();
    
    virtual int Pass(std::vector<double>& )=0;
    //virtual int Pass(CPS6D<CTPS>& ps)=0;
    virtual int SetP(const std::string& str, const double& value)=0;
    virtual int SetP(const std::string& str, const std::string& svalue)=0;
    virtual int SetP(const std::string& str, const std::vector<double>& lvalue)=0;
    virtual int SetP(const std::string& str, const std::vector<std::string>& slvalue)=0;
    virtual void Initialize()=0;
    //void SetParameter(const CElement& ele);
    //void checkaperture(CPS6D<double>& x);
};


#endif //LATTE_ELEMENT_H
