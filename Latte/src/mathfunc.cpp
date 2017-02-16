//
// Created by Yue Hao on 8/17/16.
//

#include "mathfunc.h"
void array_statistics(const std::vector<double>& arr, double& ave, double& sqrsize){
    unsigned long size=arr.size();
    if (size==1) {
        ave = arr[0];
        sqrsize=0;
        return;
    }
    double sum=0;
    #pragma omp parallel for reduction(+:sum)
    for (int i=0;i<size;i++){
        sum+=arr[i];
    }
    ave=sum/size;
    sum=0;
    #pragma omp parallel for reduction(+:sum)
    for (int i=0;i<size;i++){
        sum+=(arr[i]-ave)*(arr[i]-ave);
    }
    sqrsize=sum/size;
    return;
    
}

void array_correlation(std::vector<double>& arr1, const std::vector<double>& arr2,
                       const double& ave1, const double& ave2, double& corr){
    double sum=0;
    if (arr1.size()==1) {
        corr=0;
        return;
    }
    #pragma omp parallel for reduction(+:sum)
    for (int i =0; i< arr1.size();i++){
        sum+=arr1[i]*arr2[i];
    }
    corr=sum/arr1.size();
    return;
}



