from libcpp.vector cimport vector
from libcpp cimport bool
cdef extern from "tpsalib.h":
    cdef cppclass CTPS:
        CTPS() except +
        CTPS(const double & ) except +
        CTPS(const CTPS &) except +

        @staticmethod
        void Initialize(const int&, const int&)

        unsigned long findindex(const vector[int]& ) const
        vector[int] findpower(const unsigned long & n) const

        const int get_dim() const
        const unsigned long get_degree() const
        const unsigned long get_terms() const

        const double element(const int& ) const
        const double element(const vector[int] ) const

        double evaluate(const vector[double] ) const

        void print_by_order(const int &) const

        CTPS& operator=(const CTPS &)
        CTPS& add_to(const CTPS &)
        CTPS& minus_to(const CTPS &)
        CTPS& time_to(const CTPS &)
        CTPS& divide_to(const CTPS &)

        const double cst() const

    CTPS inv(const CTPS &)
    CTPS exp(const CTPS &)
    CTPS log(const CTPS &)
    CTPS sqrt(const CTPS &)
    CTPS pow(const CTPS &, const double &)
    CTPS sin(const CTPS &)
    CTPS cos(const CTPS &)
    CTPS tan(const CTPS &)
    CTPS sinh(const CTPS &)
    CTPS cosh(const CTPS &)
    
    const CTPS operator+(const CTPS & M)
    const CTPS operator-(const CTPS & M)

    const CTPS operator+(const CTPS & M, const CTPS & N)
    const CTPS operator-(const CTPS & M, const CTPS & N)

    const CTPS operator*(const CTPS & M, const CTPS & N)
    const CTPS operator/(const CTPS & M, const CTPS & N)

    const bool operator==(const CTPS & M, const CTPS & N)
    const bool operator<=(const CTPS & M, const CTPS & N)
    const bool operator>=(const CTPS & M, const CTPS & N)
    const bool operator!=(const CTPS & M, const CTPS & N)
    const bool operator<(const CTPS & M, const CTPS & N)
    const bool operator>(const CTPS & M, const CTPS & N)

cdef class PyTPSA:
    cdef CTPS ctps
    def __cinit__(self, a=None):
        if a is None:
            self.ctps = new CTPS()
        elif isinstance(a,float) or isinstance(a,int):
            self.ctps = new CTPS(1.0 * a)
        elif isinstance(a, PyTPSA):
            self.ctps = new CTPS(a.ctps)

    @classmethod
    def initialize(cls, int dim, int order):
        cls.ctps.Initialize(dim, order)

    def get_dim(self):
        self.ctps.get_dim()
    def get_degree(self):
        self.ctps.get_degree()
    def get_term(self):
        self.ctps.get_terms()



