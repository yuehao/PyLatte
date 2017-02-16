/*
 *  twiss4D.h
 *  ERLtrac
 *
 *  Created by Yue Hao on 6/17/11.
 *  Copyright 2011 Brookhaven National Laboratory. All rights reserved.
 *
 */

#ifndef TWISS4D_H
#define TWISS4D_H
#include "constants.h"
#include <armadillo>
class twiss4D {
private:
	arma::mat side_matrix;
public:
	twiss4D():side_matrix(arma::mat(6,6)){
		beta[X_]=1;beta[Y_]=1;
		alpha[X_]=0;alpha[Y_]=0;
		gamma[X_]=1;gamma[Y_]=1;
		tune[X_]=0;tune[Y_]=0;
		d[X_]=0;d[Y_]=0;
		dp[X_]=0;dp[Y_]=0;
	}
	friend twiss4D gettwiss(arma::mat tmap, const twiss4D& ori);
	friend arma::mat transfermap(const twiss4D& ts, const twiss4D& te);
	void inimatrix();
	void getoneturntwiss(const arma::mat& tmap);
	double beta[2], alpha[2], gamma[2], tune[2];
	double d[2], dp[2];
	
	
};




#endif