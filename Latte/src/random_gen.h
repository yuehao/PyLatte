//
// Created by Yue Hao on 8/18/16.
//

#ifndef LATTE_RANDOM_LIST_H
#define LATTE_RANDOM_LIST_H
#include <random>
#include <vector>

class random_generator {
private:
    std::mt19937 mt_generator;
public:
    random_generator(){
        seed();
    }
    random_generator(const int& seeding){
        seed(seeding);
    }
    void seed();
    void seed(const int&);
    std::vector<double> random_Gaussian_list(const unsigned long& nc, const double& mean, const double& size, const double& cutoff);
    
};


#endif //LATTE_RANDOM_LIST_H
