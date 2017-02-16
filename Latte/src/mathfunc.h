//
// Created by Yue Hao on 8/17/16.
//

#ifndef LATTE_MATHFUNC_H
#define LATTE_MATHFUNC_H

#include <vector>
#include "assert.h"
#include <iostream>
void array_statistics(const std::vector<double>&, double& ave, double& sqrsize);
void array_correlation(std::vector<double>& arr1, const std::vector<double>& arr2,
                       const double& ave1, const double& ave2, double& corr);


template <typename T>
std::vector<T>& operator+=(std::vector<T>& a, const T& c) {
    #pragma omp parallel for schedule(static)
    for (int i=0; i < a.size(); i++ ){
        a[i]+=c;
    }
    return a;
}
template <typename T>
std::vector<T>& operator+=(std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size()==b.size());
    #pragma omp parallel for schedule(static)
    for (int i=0; i < a.size(); i++ ){
        a[i]+=b[i];
    }
    return a;
}
template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b){std::vector<T> result(a); return result+=b;}
template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const T& c){std::vector<T> result(a); return result+=c;}
template <typename T>
std::vector<T> operator+(const T& c, const std::vector<T>& a){std::vector<T> result(a); return result+=c;}

template <typename T>
std::vector<T>& operator*=(std::vector<T>& a, const T& c) {
    #pragma omp parallel for schedule(static)
    for (int i=0; i < a.size(); i++ ){
        a[i]*=c;
    }
    return a;
}
template <typename T>
std::vector<T>& operator*=(std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size()==b.size());
    #pragma omp parallel for schedule(static)
    for (int i=0; i < a.size(); i++ ){
        a[i]*=b[i];
    }
    return a;
}
template <typename T>
std::vector<T> operator*(const std::vector<T>& a, const std::vector<T>& b){std::vector<T> result=a; return result*=b;}
template <typename T>
std::vector<T> operator*(const std::vector<T>& a, const T& c){std::vector<T> result=a; return result*=c;}
template <typename T>
std::vector<T> operator*(const T& c, const std::vector<T>& a){std::vector<T> result=a; return result*=c;}


template <typename T>
std::vector<T>& operator-=(std::vector<T>& a, const T& c) {
    #pragma omp parallel for schedule(static)
    for (int i=0; i < a.size(); i++ ){
        a[i]-=c;
    }
    return a;
}
template <typename T>
std::vector<T>& operator-=(std::vector<T>& a, const std::vector<T>& b) {
    assert(a.size()==b.size());
    #pragma omp parallel for schedule(static)
    for (int i=0; i < a.size(); i++ ){
        a[i]-=b[i];
    }
    return a;
}
template <typename T>
std::vector<T> operator-(const std::vector<T>& a){std::vector<T> result=a; return result*=(-1.0);}
template <typename T>
std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b){std::vector<T> result=a; return result-=b;}
template <typename T>
std::vector<T> operator-(const std::vector<T>& a, const T& c){std::vector<T> result=a; return result-=c;}
template <typename T>
std::vector<T> operator-(const T& c, const std::vector<T>& a){std::vector<T> result=(-a); return result+=c;}




#endif //LATTE_MATHFUNC_H
