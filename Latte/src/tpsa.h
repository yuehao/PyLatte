/*
 *  ctpsa.h
 *  ERLtrac
 *
 *  Created by Yue Hao on 4/8/11.
 *  Copyright 2011 Brookhaven National Laboratory. All rights reserved.
 *
 */


#ifndef TPS
#define TPS

#include <vector>
#include <iostream>

//using namespace std;
const int MAX_TPS_ORDERS=6;
const int DEFAULT_TPS_DIMS=6;

class CPolyMap{
private:
    int dim;
    int max_order;
    std::vector<int> decomposite(const int& i);
public:
    CPolyMap();
    std::vector<std::vector<int> > map;
    std::vector<int> getindexmap(const int& i) {
        return map[i];
    }
    void setindexmap();
    
    
};

class CTPS{
private:
    int degree;
    unsigned long terms;//Total terms
    std::vector<double> map;
    
    
    void redegree(const int&);
    enum errors {
        DivZero, OverFlow, DiffDim, NegValue
    };
    void ErrMsg(const CTPS::errors& err, const std::string& fromfunc) const;

public:
    static int Maximum_TPS_Degree;
    static int TPS_Dim;
    static CPolyMap polymap;
    static void Initialize(const int& dim, const int& max_order=MAX_TPS_ORDERS){
        CTPS::TPS_Dim=dim;
        CTPS::Maximum_TPS_Degree=max_order;
        polymap=CPolyMap();
    }
    CTPS();
    CTPS(const double& a);
    CTPS(const CTPS &);
    
    
    void assign(const double&); //A constant
    void assign(const double &, const int& );// A Variable
    
    unsigned long findindex(const std::vector<int>& );
    std::vector<int> findpower(const unsigned long & n);
    
    inline const int get_dim() const {return TPS_Dim;}
    inline const unsigned long get_degree() const {return degree;}
    inline const unsigned long get_terms() const {return terms;}
    
    const double element(const int& ind) const;
    const double element(std::vector<int> indmap) const;
    
    double evaluate(const std::vector<double>& inivalue) const;
    
    void print(const int & max_print_degree) const;
    
    //operator double() const;
    CTPS& operator=(const CTPS &);
    CTPS& operator+=(const CTPS &);
    CTPS& operator-=(const CTPS &);
    CTPS& operator*=(const CTPS &);
    CTPS& operator/=(const CTPS &);
    
    inline const double cst() const {return map[0];}
    
    friend CTPS inv(const CTPS &);
    friend CTPS exp(const CTPS &);
    friend CTPS log(const CTPS &);
    friend CTPS sqrt(const CTPS &);
    friend CTPS pow(const CTPS &, const double &);
    friend CTPS sin(const CTPS &);
    friend CTPS cos(const CTPS &);
    friend CTPS tan(const CTPS &);
    friend CTPS sinh(const CTPS &);
    friend CTPS cosh(const CTPS &);
    friend std::ostream& operator<<(std::ostream& output, const CTPS& A);
    
};



inline const CTPS operator+(const CTPS & M) {return M;}
inline const CTPS operator-(const CTPS & M) {return CTPS(M)*=(-1.0);}

inline const CTPS operator+(const CTPS & M, const CTPS & N) {if (M.get_degree()>N.get_degree()) return CTPS(M)+=N;else return CTPS(N)+=M;}
inline const CTPS operator-(const CTPS & M, const CTPS & N) {if (M.get_degree()>N.get_degree()) return CTPS(M)-=N;else return CTPS(-N)+=M;}

inline const CTPS operator*(const CTPS & M, const CTPS & N) {if (M.get_degree()>N.get_degree()) return CTPS(M)*=N;else return CTPS(N)*=M;}
inline const CTPS operator/(const CTPS & M, const CTPS & N) {return CTPS(M)/=N;}

inline const bool operator==(const CTPS & M, const CTPS & N) {return (M.cst()==N.cst());}
inline const bool operator<=(const CTPS & M, const CTPS & N) {return (M.cst()<=N.cst());}
inline const bool operator>=(const CTPS & M, const CTPS & N) {return (M.cst()>=N.cst());}
inline const bool operator!=(const CTPS & M, const CTPS & N) {return (M.cst()!=N.cst());}
inline const bool operator<(const CTPS & M, const CTPS & N) {return (M.cst()<N.cst());}
inline const bool operator>(const CTPS & M, const CTPS & N) {return (M.cst()>N.cst());}


unsigned long doublefactorial(const int&);
unsigned long factorial(const int&);
unsigned long binomial(const int& n, const int& m);
#endif
