//
// Created by Yue Hao on 8/25/16.
//

#include "../random_gen.h"
#include "../mathfunc.h"
#include "omp.h"
#include <iostream>

int main(){
    random_generator rd1;
    rd1.seed();
    double ave, sig, corr;
    double start_s = omp_get_wtime();
    std::vector<double> gl=rd1.random_Gaussian_list(10000000, 0, 1.0, 5.0);
    std::vector<double> gl2=rd1.random_Gaussian_list(10000000, 0, 1.0, 5.0);
    //for (int i=0;i<gl.size();i++){
    //    std::cout<<gl[i]<<'\t';
    //}
    array_statistics(gl, ave, sig);
    array_correlation(gl,gl2,0,0,corr);
    double end_s=omp_get_wtime();
    //std::cout<<test_c;
    std::cout<<(end_s-start_s)<<std::endl;
    std::cout<<std::endl;
    std::cout<<ave<<'\t'<<sig<<'\t'<<corr<<std::endl;
}
