//
// Created by Yue Hao on 8/12/16.
//

#include "bunch.h"
#include "constants.h"
#include "mathfunc.h"

void bunch::generate_distribution(const std::vector<double> sizes, random_generator& rdg) {
    if (n_particle==1){
        for (int i=0;i<dim;i++){
            coor[i][0]=0;
        }
        return;
    }
    for (int i=0;i<dim;i++){
        coor[i]=rdg.random_Gaussian_list(n_particle,0.0, sizes[i], 5.0);
    }
    
}
void bunch::statistics(std::vector<double> &ave, std::vector<double> &moment) {
    ave.resize(dim);
    moment.resize(dim*dim);
    for (int i =0;i<dim;i++){
        array_statistics(coor[i],ave[i], moment[i*dim]);
    }
    for (int i =1;i<dim;i++)
        for (int j=0; j<i; j++){
            array_correlation(coor[i], coor[j], ave[i], ave[j], moment[dim*i+j]);
            moment[dim*j+i]=moment[dim*i+j];
    }
}