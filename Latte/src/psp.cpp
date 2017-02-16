//
// Created by Yue Hao on 8/24/16.
//

#include "psp.h"
#include "constants.h"
psp::psp():_dim(Phase_Dimension), _nvar(Phase_Dimension){
    for (int i=0;i<_nvar;i++){
        CTPS::Initialize(_dim);
        CTPS varTPS;
        varTPS.assign(0.0, i+1);
        _pslist.push_back(varTPS);
    }
}

std::vector<double> psp::getcoor() const {
    std::vector<double> result(_nvar);
    for (int i=0; i<_nvar; i++) {
        result[i]=_pslist[i].cst();
    }
    return result;
}

arma::mat psp::getmap() const {
    arma::mat map(_nvar, _dim);
    for (int i=0; i<_nvar; i++)
        for (int j=0; j<_dim; j++){
            map(i,j)=_pslist[i].element(j+1);
    }
    return map;
}

void psp::clear() {
    for (int i=0; i<_nvar; i++){
        _pslist[i]=CTPS(0.0);
    }
    return;
}