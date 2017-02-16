//
// Created by Yue Hao on 8/12/16.
//

#ifndef LATTE_BUNCH_H
#define LATTE_BUNCH_H
#include <vector>
#include <cstdlib>
#include "random_gen.h"
#include "constants.h"

class bunch {
private:
    std::vector<std::vector<double> > coor;
public:
    static int dim;
    unsigned long n_particle;
    
    bunch(const unsigned long & np): n_particle(np), dim(Phase_Dimension){
        std::vector<double> temp(np);
        for (int i=0; i<dim; i++){
            coor.push_back(temp);
        }
        if (np<1) exit(-1);
    };
    std::vector<double>& operator[](const int &i) { return coor[i]; }
    const std::vector<double>& operator[](const int &i) const { return coor[i]; }
    void generate_distribution(const std::vector<double> sizes, random_generator& rdg);
    void statistics(std::vector<double>& ave, std::vector<double>& moment);
    
    
};


#endif //LATTE_BUNCH_H
