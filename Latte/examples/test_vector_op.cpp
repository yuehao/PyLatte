//
// Created by Yue Hao on 8/31/16.
//

#include "../mathfunc.h"
#include <iostream>
#include "omp.h"
int main(){
    std::vector<double> v1,v2;
    
    unsigned long size=1000000;
    v1.resize(size);
    v2.resize(size);
    #pragma omp parallel
    {
        if (omp_get_thread_num()==0) std::cout<<"# of thread: "<<omp_get_num_threads()<<std::endl;
    }
    
    
    double start_s = omp_get_wtime();
    #pragma omp parallel for
    for (int i=0;i<size;i++) {v1[i]=(i+1);}
    #pragma omp parallel for
    for (int i=0;i<size;i++) {v2[i]=(i);}
    double end_s=omp_get_wtime();
    //std::cout<<test_c;
    std::cout<<(end_s-start_s)<<std::endl;
    
    
    std::vector<double> v3(size), v4(size);
    start_s = omp_get_wtime();
    for (int j=0;j<1000;j++) {
        
        #pragma omp parallel for
        for (int i = 0; i < size; i++) {
            v3[i] = (v1[i] + 5.0) ;
        }
        
    }
    end_s=omp_get_wtime();
    //std::cout<<test_c;
    std::cout<<(end_s-start_s)<<std::endl;
    start_s = omp_get_wtime();
    for (int j=0;j<1000;j++) v3=(v1+=5.0);
    end_s=omp_get_wtime();
    //std::cout<<test_c;
    std::cout<<(end_s-start_s)<<std::endl;
    
    
    
    
}