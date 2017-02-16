/*
 *  twiss4D.cpp
 *  ERLtrac
 *
 *  Created by Yue Hao on 6/17/11.
 *  Copyright 2011 Brookhaven National Laboratory. All rights reserved.
 *
 */

#include "twiss4D.h"

#include <iostream>


using namespace std;
twiss4D gettwiss(arma::mat tmap, const twiss4D& ori){
	twiss4D result;
	tmap=tmap*ori.side_matrix;
	
	result.d[X_]=tmap(x_,dE_);
	result.dp[X_]=tmap(px_,dE_);
	result.d[Y_]=tmap(y_,dE_);
	result.dp[Y_]=tmap(py_,dE_);
	
	result.beta[X_]=pow(tmap(x_,x_),2.0)+pow(tmap(x_,px_),2.0);
	double phase=atan2(tmap(x_,px_), tmap(x_,x_));
	result.tune[X_]=ori.tune[X_]+phase/2/M_PI;
	result.alpha[X_]=-(cos(phase)*tmap(px_,x_)+sin(phase)*tmap(px_,px_))*sqrt(result.beta[X_]);
	result.gamma[X_]=(1+result.alpha[X_]*result.alpha[X_])/result.beta[X_];
	
	result.beta[Y_]=pow(tmap(y_,y_),2.0)+pow(tmap(y_,py_),2.0);
	phase=atan2(tmap(y_,py_), tmap(y_,y_));
	result.tune[Y_]=ori.tune[Y_]+phase/2/M_PI;
	result.alpha[Y_]=-(cos(phase)*tmap(py_,y_)+sin(phase)*tmap(py_,py_))*sqrt(result.beta[Y_]);
	result.gamma[Y_]=(1+result.alpha[Y_]*result.alpha[Y_])/result.beta[Y_];
	return result;
}

void twiss4D::getoneturntwiss(const arma::mat& tmap){//No coupling only
	
	double phase=acos((tmap(x_,x_)+tmap(px_,px_))/2.0);
	tune[X_]=0;tune[Y_]=0;
	beta[X_]=tmap(x_,px_)/sin(phase);
	if (beta[X_]<0) {
		beta[X_]*=(-1.0);
		phase=-phase;
	}
	alpha[X_]=(tmap(x_,x_)-tmap(px_,px_))/2.0/sin(phase);
	gamma[X_]=(1+alpha[X_]*alpha[X_])/beta[X_];
	d[X_]=((1-cos(phase)+alpha[X_]*sin(phase))*tmap(x_,dE_)+beta[X_]*sin(phase)*tmap(px_,dE_))/2.0/(1-cos(phase));
	dp[X_]=(-tmap(x_,dE_)*gamma[X_]*sin(phase)+tmap(px_,dE_)*(1-cos(phase)-alpha[X_]*sin(phase)))/2.0/(1-cos(phase));
	
	phase=acos((tmap(y_,y_)+tmap(py_,py_))/2.0);
	
	beta[Y_]=tmap(y_,py_)/sin(phase);
	if (beta[Y_]<0) {
		beta[Y_]*=(-1.0);
		phase=-phase;
	}
	alpha[Y_]=(tmap(y_,y_)-tmap(py_,py_))/2.0/sin(phase);
	gamma[Y_]=(1+alpha[Y_]*alpha[Y_])/beta[Y_];
	d[Y_]=((1-cos(phase)+alpha[Y_]*sin(phase))*tmap(y_,dE_)+beta[Y_]*sin(phase)*tmap(py_,dE_))/2.0/(1-cos(phase));
	dp[Y_]=(-tmap(y_,dE_)*gamma[Y_]*sin(phase)+tmap(py_,dE_)*(1-cos(phase)-alpha[Y_]*sin(phase)))/2.0/(1-cos(phase));

	this->inimatrix();
}


void twiss4D::inimatrix(){
	side_matrix.eye(6,6);
	side_matrix(x_,x_)=sqrt(beta[X_]);
	side_matrix(px_,x_)=-alpha[X_]/sqrt(beta[X_]);
	side_matrix(px_,px_)=1/sqrt(beta[X_]);
    side_matrix(y_,y_)=sqrt(beta[Y_]);
    side_matrix(py_,y_)=-alpha[Y_]/sqrt(beta[Y_]);
    side_matrix(py_,py_)=1/sqrt(beta[Y_]);
    side_matrix(x_,dE_)=d[X_];
    side_matrix(px_,dE_)=dp[X_];
    side_matrix(y_,dE_)=d[Y_];
    side_matrix(py_,dE_)=dp[Y_];
}

arma::mat transfermap(const twiss4D& ts, const twiss4D& te){
    arma::mat result,rot;
	double phix=te.tune[X_]-ts.tune[X_];
	double phiy=te.tune[Y_]-ts.tune[Y_];
	rot.eye(4,4);
	
	rot(x_,x_)=cos(phix); rot(x_,px_)=sin(phix);
	rot(px_,x_)=-sin(phix); rot(px_,px_)=cos(phix);
	
	rot(y_,y_)=cos(phiy); rot(y_,py_)=sin(phiy);
	rot(py_,y_)=-sin(phiy); rot(py_,py_)=cos(phiy);
	
	result.eye(6,6);
    result.submat(0,0,arma::size(4,4))=te.side_matrix.submat(0,0,arma::size(4,4))*rot*(ts.side_matrix.submat(0,0,arma::size(4,4)).i());
	
	return result;
}

