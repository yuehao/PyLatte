//
// Created by Yue Hao on 8/18/16.
//

#include "random_gen.h"

void random_generator::seed() {
    std::random_device rd;
    mt_generator=std::mt19937(rd());
}
void random_generator::seed(const int& seed) {
    mt_generator.seed(seed);
}
std::vector<double> random_generator::random_Gaussian_list(const unsigned long& nc, const double &mean, const double &size, const double& cutoff) {
    std::vector<double> list(nc);
    std::normal_distribution<> gaudis(0.0,1.0);
    
    #pragma omp parallel for schedule(static)
    for (int i=0;i<nc;i++){
        double temp=gaudis(mt_generator);
        while (fabs(temp=gaudis(mt_generator)) > cutoff) {
            temp=gaudis(mt_generator);
        }
        list[i]=temp;
    }
    double sum=0;
    #pragma omp parallel for reduction(+:sum)
    for (int i=0;i<nc;i++){
        sum+=list[i];
    }
    double ave=sum/nc;
    #pragma omp parallel for schedule(static)
    for (int i=0;i<nc;i++){
        list[i]-=ave;
    }
    sum=0;
    #pragma omp parallel for reduction(+:sum)
    for (int i=0;i<nc;i++){
        sum+=list[i]*list[i];
    }
    double rmssize=sqrt(sum/nc);
    #pragma omp parallel for schedule(static)
    for (int i=0;i<nc;i++){
        list[i]/=rmssize;
    }//*/
    return list;
}
