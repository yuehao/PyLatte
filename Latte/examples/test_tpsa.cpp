//
// Created by Yue Hao on 8/25/16.
//

#include "../src/tpsa.h"
#include "omp.h"
#include <iostream>

int main(){
    #pragma omp parallel
    {
        if (omp_get_thread_num()==0) std::cout<<"# of thread: "<<omp_get_num_threads()<<std::endl;
    }
    double start_s = omp_get_wtime();
    CTPS::Initialize(6,20);
    double end_s=omp_get_wtime();
    std::cout<<"Time for initialization: "<<(end_s-start_s)<<std::endl;
    start_s = omp_get_wtime();
    CTPS test_a, test_b, test_c, test_d;
    test_a.assign(1.3,1);
    test_b.assign(1.0,2);
    test_c.assign(1.0,3);
    for (int i=0;i<10000;i++){
        test_d=test_a-test_b+test_c+test_a-test_b+test_c+test_a-test_b+test_c+test_a-test_b+test_c+test_a-test_b+test_c+test_a-test_b+test_c+test_a-test_b+test_c+test_a-test_b+test_c+test_a-test_b+test_c
        +test_a-test_b+test_c+test_a-test_b+test_c+test_a-test_b+test_c+test_a-test_b+test_c+test_a-test_b+test_c+test_a-test_b+test_c+test_a-test_b+test_c+test_a-test_b+test_c+test_a-test_b+test_c+test_a-test_b+test_c+test_a-test_b+test_c+test_a-test_b+test_c+test_a-test_b+test_c;;
    }
    end_s=omp_get_wtime();
    
    
    
    std::cout<<"Time for calculation: "<<(end_s-start_s)<<std::endl;
    //test_d.print(3);
}
