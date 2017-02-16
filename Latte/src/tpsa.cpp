/*
 *  ctpsa.cpp
 *  ERLtrac
 *
 *  Created by Yue Hao on 4/8/11.
 *  Copyright 2011 Brookhaven National Laboratory. All rights reserved.
 *
 */

#include "tpsa.h"
#include <cmath>
#include <iostream>

using namespace std;

unsigned long factorial(const int & n) {
    if (n<0) {
        return 0;
    }
    switch (n) {
        case 0:
            return 1;
        case 1:
            return 1;
        case 2:
            return 2;
        case 3:
            return 6;
        case 4:
            return 24;
        case 5:
            return 120;
        case 6:
            return 720;
        case 7:
            return 5040;
        case 8:
            return 40320;
        case 9:
            return 362880;
        case 10:
            return 3628800;
        case 11:
            return 39916800;
        case 12:
            return 479001600;
        
        
        default:
            return n*factorial(n-1);
    }
}
unsigned long doublefactorial(const int& n){
    if (n<0) {
        return 0;
    }
    switch (n) {
        case 0:
            return 1;
        case 1:
            return 1;
        case 2:
            return 2;
        case 3:
            return 3;
        case 4:
            return 8;
        case 5:
            return 15;
        case 6:
            return 48;
        case 7:
            return 105;
        case 8:
            return 384;
        case 9:
            return 945;
        case 10:
            return 3840;
        case 11:
            return 10395;
        case 12:
            return 46080;
        case 13:
            return 135135;
        case 14:
            return 645120;
        case 15:
            return 2027025;
        case 16:
            return 10321920;
        case 17:
            return 34459425;
        case 18:
            return 185794560;
        
        default:
            return n*doublefactorial(n-2);
    }
    
}
unsigned long binomial(const int& n, const int& m) {
    if (n<=0 || m>n || m<0) {
        return 0;
    }
    int ml;
    if (m>n/2) ml=n-m;
    else ml=m;
    switch (ml) {
        case 0:
            return 1;
        case 1:
            return (unsigned long)n;
        case 2:
            return (unsigned long)n*(n-1)/2;
            
        case 3:
            return (unsigned long)n*(n-1)*(n-2)/6;
            
        
        default:
            break;
    }
    if (n<=12) {
        return factorial(n)/factorial(m)/factorial(n-m);
    }
    else {
        return (n*binomial(n-1,ml-1)/ml);
    }
}


CPolyMap::CPolyMap():dim(CTPS::TPS_Dim),max_order(CTPS::Maximum_TPS_Degree) {
    setindexmap();
}
vector<int> CPolyMap::decomposite(const int& n){
    vector<int> result((unsigned long)dim+1);
    int itemp=n+1;
    for (int i=dim; i>0; i--) {
        int k=i-1;
        while (binomial(k, i)<itemp) {
            k++;
        }
        itemp-=binomial(k-1, i);
        result[dim-i]=k-i;
    }
    for (int i=dim;i>0;i--) {
        result[i]=result[i-1]-result[i];
    }
    return result;
}
void CPolyMap::setindexmap(){
    unsigned long totallength=binomial(max_order+dim, dim);
    map.resize(totallength);
    #pragma omp parallel for schedule(dynamic)
    for (int i=0; i<totallength; i++) {
        map[i]=this->decomposite(i);
    }
    return;
}


void CTPS::ErrMsg(const CTPS::errors& err, const string& fromfunc) const{
    switch (err) {
        case OverFlow:
            cout << "The function "<<fromfunc<<" in class CTPS, INDEX OverFlow";
            break;
        case DivZero:
            cout << "The function "<<fromfunc<<" in class CTPS, Divided by zero";
            break;
        case DiffDim:
            cout << "The function "<<fromfunc<<" in class CTPS, has different dimension";
            break;
        case NegValue:
            cout << "The function "<<fromfunc<<" in class CTPS, has negative parameter";
            break;
        default:
            break;
    }
    exit(-1);
    return;
}


int CTPS::Maximum_TPS_Degree=MAX_TPS_ORDERS;
int CTPS::TPS_Dim=DEFAULT_TPS_DIMS;
CPolyMap CTPS::polymap=CPolyMap();

CTPS::CTPS(){
    //int dim=CTPS::TPS_Dim;
    //degree=1;
    //terms=(unsigned long)dim+1;
    //map.reserve(binomial(dim+CTPS::Maximum_TPS_Degree, dim));
    this->assign(0.0);
    
    
    
}
CTPS::CTPS(const double& a){
    //int dim=TPS_Dim;
    //degree=1;
    //terms=(unsigned long)dim+1;
    //map.reserve(binomial(dim+CTPS::Maximum_TPS_Degree, dim));
    this->assign(a);
    
    
}

CTPS::CTPS(const CTPS &M){
    int dim=TPS_Dim;
    this->degree=M.degree;
    this->terms=M.terms;
    //map.reserve(binomial(dim+CTPS::Maximum_TPS_Degree, dim));
    this->map=M.map;
    
    
}
//CTPS::operator double() const{
//  return this->map[0];
//}
unsigned long CTPS::findindex(const vector<int>& indexmap){
    int dim=TPS_Dim;
    vector<int> sum((unsigned long)dim+1);
    sum[0]=indexmap[0];
    
    for (int i=1; i<=dim; i++) {
        sum[i]=sum[i-1]-indexmap[i];
    }
    unsigned long result=0;
    for (int i=dim; i>0; i--) {
        if (sum[dim-i]==0) {
            break;
        }
        result+=binomial(sum[dim-i]-1+i, i);
    }
    return result;
}
vector<int> CTPS::findpower(const unsigned long &n) {
    return this->polymap.getindexmap(n);
}

void CTPS::redegree(const int& degree){
    this->degree=degree;
    if (degree>CTPS::Maximum_TPS_Degree) this->degree=CTPS::Maximum_TPS_Degree;
    terms=binomial(TPS_Dim+degree, degree);
    map.resize(terms);
}


void CTPS::assign(const double& a, const int& n_var){
    
    if (n_var<=this->TPS_Dim && n_var>0) {
        this->degree=1;
        this->terms=(unsigned long)CTPS::TPS_Dim+1;
        map.clear();
        //map.reserve(binomial(CTPS::TPS_Dim+CTPS::Maximum_TPS_Degree, CTPS::TPS_Dim));
        map.assign(this->terms, 0.0);
        this->map[n_var]=1;
        this->map[0]=a;
    }
    else this->ErrMsg(OverFlow, "assign (,)");
}
void CTPS::assign(const double& a){
    this->degree=0;
    this->terms=1;
    map.clear();
    map.reserve(binomial(CTPS::TPS_Dim+CTPS::Maximum_TPS_Degree, CTPS::TPS_Dim));
    map.assign(this->terms, 0.0);
    this->map[0]=a;
}

const double CTPS::element(const int& ind) const{
    if (ind<0 || ind >=terms) {
        ErrMsg(OverFlow, "element (int)");
        exit(0);
    }
    return map[ind];
}

const double CTPS::element(vector<int> ind) const{
    int dim=TPS_Dim;
    if (ind.size()!=dim+1) {
        ErrMsg(OverFlow, "element(vector<int>)");
        exit(0);
    }
    for (int i=1; i<=dim; i++) {
        ind[i]=ind[i-1]-ind[i];
    }
    if (ind[dim]!=0 || ind[0]>this->degree) {
        ErrMsg(OverFlow, "element(vector<int>)");
        exit(0);
    }
    int result=0;
    for (int i=dim; i>0; i--) {
        if (ind[dim-i]==0) {
            break;
        }
        result+=binomial(ind[dim-i]-1+i, i);
    }
    return map[result];
}

double CTPS::evaluate(const vector<double>& inivalue) const{
    if (inivalue.size()!=this->TPS_Dim) {
        ErrMsg(DiffDim, "evaluate (const vector<double>&)");
    }
    double sum=map[0];
    #pragma omp parallel for
    for (int i=1; i<this->terms;i++){
        vector<int> temp=this->polymap.getindexmap(i);
        double product=1.0;
        for (int j=0; j<this->TPS_Dim; j++) {
            product=product*std::pow(inivalue[j], temp[j+1]);
        }
        sum+=product*map[i];
    }
    return sum;
}

void CTPS::print(const int &max_print_degree) const {
    int current_order=0;
    int acc_return=0;
    int just_endl=0;
    std::cout<<"Order 0:"<<endl;
    for (int i=0;i<this->terms;i++){
        
        vector<int> temp=this->polymap.getindexmap(i);
        if (temp[0] > current_order){
            if (just_endl==0) std::cout<<endl;
            std::cout<<endl<<"Order "<<temp[0]<<":"<<endl;
            current_order=temp[0];
            acc_return=0;
        }
        std::cout<<"(";
        for (int j=1;j<temp.size();j++) {
            std::cout<<temp[j];
            if (j<temp.size()-1) std::cout<<',';
        }
        std::cout<<"): ";
        std::cout.width(10);
        std::cout<<std::left<<this->map[i];
        acc_return++;
        
        if (acc_return % 3 == 0) {std::cout << endl; just_endl=1;}
        else {std::cout << '\t'; just_endl=0;}
           
        
    }
    if (just_endl==0) std::cout<<endl;
}


CTPS& CTPS::operator=(const CTPS &M){
    if (this != &M){
        this->degree=M.degree;
        this->terms=M.terms;
        this->map=M.map;
    }
    return *this;
}
CTPS& CTPS::operator+=(const CTPS &M){
    
    if (this->degree<M.degree) {
        this->redegree(M.degree);
    }
    //#pragma omp parallel for schedule(static)
    for (int i=0; i<M.terms; i++) {
        this->map[i]+=M.map[i];
    }
    return *this;
    
}
CTPS& CTPS::operator-=(const CTPS &M){
    
    if (this->degree<M.degree) {
        this->redegree(M.degree);
    }
    //#pragma omp parallel for schedule(static)
    for (int i=0; i<M.terms; i++) {
        this->map[i]-=M.map[i];
    }
    return *this;
}
CTPS& CTPS::operator*=(const CTPS &M){
    
    if (M.get_degree()==0){
        
        //#pragma omp parallel for schedule(static,4096)
        for (int i=0; i<this->terms; i++){
            this->map[i]*=M.map[0];
        }
        return *this;
    }
    
    CTPS temp(*this);
    this->map.clear();
    (*this).redegree(min(CTPS::Maximum_TPS_Degree, this->degree+M.degree));
    //#pragma omp parallel for schedule(dynamic)
    for (int i=0; i<temp.map.size(); i++ ) {
        if (temp.map[i]==0) {
            continue;
        }
        vector<int> vthis=polymap.getindexmap(i);
        //#pragma omp parallel for schedule(dynamic)
        for (int j=0; j<M.map.size();j++) {
            if (M.map[j]==0) {
                continue;
            }
            vector<int> vm=polymap.getindexmap(j);
            if (vthis[0]+vm[0]>CTPS::Maximum_TPS_Degree) break;
            vector<int> indexmap(vthis.size());
            for (int k=0; k<vthis.size(); k++) {
                indexmap[k] = vthis[k] + vm[k];
            }
            unsigned long target_ind=this->findindex(indexmap);
            //#pragma omp atomic
            this->map[target_ind]+=temp.map[i]*M.map[j];
        }
    }
    return *this;
}

CTPS& CTPS::operator/=(const CTPS &M){
    
    if (M.cst()==0) {
        this->ErrMsg(DivZero, "operator /= (CTPS)");
        exit(0);
    }
    if (M.get_degree()==0){
        #pragma omp parallel for schedule(static)
        for (int i=0; i<this->terms; i++){
            this->map[i]/=M.map[0];
        }
        return *this;
    }
    return ((*this)*=inv(M));
    
    
}


CTPS inv(const CTPS & M){
    CTPS temp(M), sum(0.0), temp2;
    double a0=M.cst();
    if (a0==0) {
        M.ErrMsg(CTPS::DivZero, "inv()");
        exit(0);
    }
    temp=temp-a0;
    for (int i=1; i<=CTPS::Maximum_TPS_Degree; i++) {
        double index=pow(-1.0, i)/pow(a0, i+1);
        if (i==1) {
            temp2=temp;
        }
        else {
            temp2=temp2*temp;
        }
        sum+=(temp2*index);
        
    }
    sum+=(1/a0);
    return sum;
    
}
CTPS exp(const CTPS & M){
    CTPS temp(M), sum(0.0), temp2;
    double a0=M.cst();
    temp=temp-a0;
    for (int i=1; i<=CTPS::Maximum_TPS_Degree; i++) {
        double index=1.0/factorial(i);
        if (i==1) {
            temp2=temp;
        }
        else {
            temp2=temp2*temp;
        }
        sum+=(temp2*index);
        
    }
    sum+=1.0;
    sum*=exp(a0);
    return sum;
}

CTPS log(const CTPS & M){
    CTPS temp(M), sum(0.0), temp2;
    double a0=M.cst();
    if (a0<=0) {
        M.ErrMsg(CTPS::NegValue, "log()");
        exit(0);
    }
    temp=temp-a0;
    for (int i=1; i<=CTPS::Maximum_TPS_Degree; i++) {
        double index=pow(-1.0, i+1)/i/pow(a0, i);
        if (i==1) {
            temp2=temp;
        }
        else {
            temp2=temp2*temp;
        }
        sum+=(temp2*index);
        
    }
    sum+=log(a0);
    return sum;
}
CTPS sqrt(const CTPS & M){
    CTPS temp(M), sum(0.0), temp2;
    double a0=M.cst();
    if (a0<0) {
        M.ErrMsg(CTPS::NegValue, "sqrt()");
        exit(0);
    }
    temp=temp-a0;
    for (int i=1; i<=CTPS::Maximum_TPS_Degree; i++) {
        double index=pow(-1.0, i+1)*doublefactorial(2*i-3)/pow(a0, i-0.5)/doublefactorial(2*i);
        
        if (i==1) {
            temp2=temp;
            index=1/2.0/sqrt(a0);
        }
        else {
            temp2=temp2*temp;
        }
        sum+=(temp2*index);
        
    }
    sum+=sqrt(a0);
    return sum;
}
CTPS pow(const CTPS & M, const double& b){
    if (b==1.0) return M;
    if (b==0.0) {
        CTPS temp(1.0);
        return temp;
    }
    CTPS temp(M), sum(0.0), temp2;
    double a0=M.cst();
    temp=temp-a0;
    double index=b, factor=pow(a0,b);
    double f0=factor;
    for (int i=1; i<=CTPS::Maximum_TPS_Degree; i++) {
        factor=factor/a0*index/(i*1.0);
        index--;
        if (i==1) {
            temp2=temp;
        }
        else {
            temp2=temp2*temp;
        }
        sum+=(temp2*factor);
        if (index==0.0) break;
    }
    sum+=f0;
    return sum;
}

CTPS sin(const CTPS & M){
    CTPS temp(M), sum(0.0), temp2;
    double a0=M.cst();
    temp=temp-a0;
    for (int i=1; i<=CTPS::Maximum_TPS_Degree; i++) {
        double index;
        if (i%2==1) {
            index=cos(a0)*pow(-1.0, (i-1)/2)/factorial(i);
        }
        else {
            index=sin(a0)*pow(-1.0, i/2)/factorial(i);
        }
        
        if (i==1) {
            temp2=temp;
        }
        else {
            temp2=temp2*temp;
        }
        sum+=(temp2*index);
        
    }
    sum+=sin(a0);
    return sum;
}
CTPS cos(const CTPS & M){
    CTPS temp(M), sum(0.0), temp2;
    double a0=M.cst();
    temp=temp-a0;
    for (int i=1; i<=CTPS::Maximum_TPS_Degree; i++) {
        double index;
        if (i%2==1) {
            index=sin(a0)*pow(-1.0, (i+1)/2)/factorial(i);
        }
        else {
            index=cos(a0)*pow(-1.0, i/2)/factorial(i);
        }
        
        if (i==1) {
            temp2=temp;
        }
        else {
            temp2=temp2*temp;
        }
        sum+=(temp2*index);
        
    }
    sum+=cos(a0);
    return sum;
}
CTPS tan(const CTPS & M){
    return sin(M)/cos(M);
}
CTPS sinh(const CTPS & M){
    CTPS temp(M), sum(0.0), temp2;
    double a0=M.cst();
    temp=temp-a0;
    for (int i=1; i<=CTPS::Maximum_TPS_Degree; i++) {
        double index;
        if (i%2==1) {
            index=cosh(a0)/factorial(i);
        }
        else {
            index=sinh(a0)/factorial(i);
        }
        
        if (i==1) {
            temp2=temp;
        }
        else {
            temp2=temp2*temp;
        }
        sum+=(temp2*index);
        
    }
    sum+=sinh(a0);
    return sum;
}
CTPS cosh(const CTPS & M){
    CTPS temp(M), sum(0.0), temp2;
    double a0=M.cst();
    temp=temp-a0;
    for (int i=1; i<=CTPS::Maximum_TPS_Degree; i++) {
        double index;
        if (i%2==1) {
            index=sinh(a0)/factorial(i);
        }
        else {
            index=cosh(a0)/factorial(i);
        }
        
        if (i==1) {
            temp2=temp;
        }
        else {
            temp2=temp2*temp;
        }
        sum+=(temp2*index);
        
    }
    sum+=cosh(a0);
    return sum;
}
std::ostream& operator<<(std::ostream& output, const CTPS& A){
    output.precision(8);
    
    for (int i=0;i<A.terms;i++) {
        output.width(16);
        output <<left<<A.map[i];
    }
    output << std::endl;
    return output;
}
