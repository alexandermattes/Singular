#ifndef LATTICE_HPP
#define LATTICE_HPP

#include <polys/monomials/p_polys.h>


// Algorithms based on chapter 2.6 of 
// "A course in computational algebraic numbertheory" by Henry Cohen
class lattice {
    
    private:
        
        //array of basisvectors
        bigintmat ** basis;
        
        // 2 dimensional triangular array for gram matrix
        // only used when the lattice is defined by a gram matrix
        // access via gram_matrix_rawset(int i,int j, number n) 
        // and gram_matrix_view(int i,int j)
        number ** gram_matrix_content;
        
        //size of basis
        int n;
        
        //length of basisvectors
        int m;
        
        coeffs coef;
                
        //constant for LLL
        number c;
        
        //array of basisvectors of LLL-basis
        bigintmat ** b; 
        
        //for LLL, see Cohen 2.6.3
        bigintmat ** b_star;
        
        //array of B_i, see Cohen 2.6.3
        number * B; 
        
        //for LLL, see Cohen 2.6.3
        bigintmat * H; 
        
        //for LLL, see Cohen 2.6.3
        bigintmat * my; 
        
        //array for integral LLL, see Cohen 2.6.7
        number * d;
        
        //for integral LLL, see Cohen 2.6.7
        bigintmat * lambda;
        
        //rank of the lattice
        int rank;
        
        //if true transformation matrix H will be computed
        bool trans_matrix;
        
        //true if basisvectors should be considered independent
        bool independentVectors;
        
        //true if gram matrix is integral
        bool integral;
        
        //true if the lattice is only defined by the gram matrix of the basis
        bool only_gram_matrix_given;
        
         //triangular matrix for enumeration
        bigintmat * Q;
        
        //testing element for enumeration
        bigintmat * x;
        
        //array bound for x_i in enumeration
        number * bound;
        
        //coef for output
        coeffs fieldcoef;
            
        //modify gram_matrix_content
        inline void gram_matrix_rawset(int i,int j, number n);
        
        //read from gram_matrix_content
        inline number gram_matrix_view(int i,int j);
        
        inline void delete_LLL_computations();
        
        inline void print_all_in_LLL();
        
        inline void cancel_ratios_in_LLL();
        
        //for LLL, see Cohen 2.6.3
        inline void RED(int k, int l);
        
        //for Integral LLL, see Cohen 2.6.7
        inline void REDI(int k, int l);
        
        //for LLL, see Cohen 2.6.3
        inline void SWAP(int k, int k_max);
        
        //for Integral LLL, see Cohen 2.6.7
        inline void SWAPI(int k, int k_max);
        
        //for MLLL, see Cohen 2.6.8
        inline void SWAPG(int k, int k_max);
        
        //for Integral MLLL, see Cohen 2.6.7 and 2.6.8
        inline void SWAPG_integral(int k, int k_max);
        
        //for LLL, see Cohen 2.6.3
        inline bool gram_schmidt(int k);
        
        //for Integral LLL, see Cohen 2.6.7
        inline bool gram_schmidt_integral(int k);
        
        //for MLLL, see Cohen 2.6.8
        inline void gram_schmidt_MLLL(int k);
        
        //for Integral MLLL, see Cohen 2.6.7 and 2.6.8
        inline void gram_schmidt_MLLL_integral(int k);
        
        //to get the next lattice element lower than the given bound starting with x an return the norm of the element
        inline number enumerate_get_next();
        
        //algorithm 4.4 from the script of Prof. Fieker of "Algorithmic Number Theory"
        inline bool quadratic_supplement();
        
        inline void increase_x(int index);
        
        //checks that x is lower than the given bound
        inline number check_bound(int index);

    public:
        
        //constructor creates new lattice spanned by columns of inputmatrix
        //if use_as_gram_matrix=true => use input as gram_matrix instead
        lattice(bigintmat* inputmatrix, bool use_as_gram_matrix=false); 
          
        //destructor
        ~lattice();
        
        
        void Write();
        char* String();
//         void Print();
        
        //LLL with c=3/4 and auto for other flags
        bool LLL();
        
        //Cohen Chapter 2.6
        bool LLL(number& c, coeffs c_coef, bool trans_matrix=true, bool integral=false, bool independentVectors=false);
        
        //return the number of rows
        inline int get_dim() {return n;};
        
        //return the basis from which the lattic was generated
        bigintmat * get_basis();
        
        //return a LLL reduced basis, after LLL was computed or NULL
        bigintmat * get_reduced_basis();
        
        //return the transformation matrix to get from basis to a LLL reduced basis if LLL was computed
        bigintmat * get_transformation_matrix();
        
        //return the gram matrix of the lattice
        bigintmat * get_gram_matrix();
        
        //return an element of the lattice which is a linear combination of basis with x
        bigintmat * get_lattice_element(bigintmat * x);
        
        //all enumeration assume the last nonzero element is greater 0, it starts at the lowest index
        //return a matrix with all lattice elements lower than a sorted by their length but at most 10^6 elements
        bigintmat * enumerate_all(number a);
        
        // return next element in lattice lower than a with basis representation close to x where the last nonzero element is greater 0 
        bigintmat * enumerate_next(number a, bigintmat * x);
        
        // return next element in lattice lower than a close to 0 or to x if it was defined
        bigintmat * enumerate_next(number a);
        
        // return next element in lattice with basis representation close to x if a was given
        bigintmat * enumerate_next(bigintmat * x);
        
        //return next element with given bound until none will found it return NULL
        bigintmat * enumerate_next();
                
};

//NOTE: Most of the following procedures should be moved to somewhere else


//NOTE: could be moved to bigintmat
number scalarproduct(bigintmat * a, bigintmat * b);

//NOTE: should be moved to somewhere else
number round(number r, coeffs coef);

//NOTE: could be moved to bigintmat
void bimnlNormalize(bigintmat * m);


//return r1 and latticeNF defined by mapping all elements to the conjugates in the number field defined by poly
int minkowski(bigintmat * basiselements, number * poly,int deg, coeffs coef, int precision, lattice * latticeNF);

//test if a is real
bool IsReal(number a, coeffs coef);

// test if the imaginary part of a is greater zero
bool ImagGreaterZero(number a, coeffs coef);

// return the squareroot of a computed by the heron algorithm
number squareroot(number a, coeffs coef, int iteration);// iteration in Heron algorithm

//get nice polynom for field over Q
poly get_nice_poly(poly polynom_in);

//test if an element from a lattice generated from an order is has a primitive root
bool is_primitive(bigintmat * element,int r1, int precision, poly out, const ring polyring);


//poly to number array and vice versa

//maps poly to an array of coefficients
int poly2numbers(poly gls, number * &pcoeffs, ring polyring, coeffs coef);

//return a polynomial over polyring from an array of number which represent the coefficents 
poly numbers2poly(number * univpol, int deg, coeffs coef, ring polyring);


//T2-norm
//return the T2-norm of an element mapped with minkowski
number t2norm(bigintmat * elt);

//return T2-norm of an polynomial definded by pol with precision
number t2norm(number * pol, int deg, coeffs coef, int precision);

//return T2-norm of polynom with precision
number t2norm(poly polynom, ring polyring, coeffs coef, int precision);

//
number elementary_symmetric_function(number * roots, int deg, int si, int lower_bound, coeffs coef);


#endif
