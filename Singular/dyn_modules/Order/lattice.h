#ifndef LATTICE_HPP
#define LATTICE_HPP

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
        
//         int precision;
    
        //modify gram_matrix_content
        inline void gram_matrix_rawset(int i,int j, number n);
        
        //read from gram_matrix_content
        inline number gram_matrix_view(int i,int j);
        
        inline void delete_LLL_computations();
        
        inline void DEBUG_LLL();
        
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
        
//         inline void INSERT(int k, int i);
//         inline bool gram_matrix(int k);
                
        inline number enumerate_get_next();
        
        inline bool quadratic_supplement();
        
        inline void increase_x(int index);
        
        inline number check_bound(int index);

    public:
        
        //constructor creates new lattice spanned by columns of inputmatrix
        //if use_as_gram_matrix=true => use input as gram_matrix instead
        lattice(bigintmat* inputmatrix, bool use_as_gram_matrix=false); 
          
        //destructor
        ~lattice();
        
        //LLL with c=3/4 and auto for other flags
        bool LLL();
        
        //Cohen Chapter 2.6
        bool LLL(number& c, coeffs c_coef, bool trans_matrix=true, bool integral=false, bool independentVectors=false);
        
        bigintmat * get_basis();
        
        bigintmat * get_reduced_basis();
        
        bigintmat * get_transformation_matrix();
        
        bigintmat * get_gram_matrix();
         
        bigintmat * get_lattice_element(bigintmat * x);
        
        bigintmat * enumerate_all(number a);
        
        bigintmat * enumerate_next(number a, bigintmat * x);
        
        bigintmat * enumerate_next(number a);
        
        bigintmat * enumerate_next(bigintmat * x);
        
        bigintmat * enumerate_next();
        
//         void set_precision(int a);
        
};


//NOTE: could be moved to bigintmat
inline number scalarproduct(bigintmat * a, bigintmat * b);


//minkowski
bigintmat * minkowksi(bigintmat ** elementarray,int size_elementarray, number * poly,int deg, coeffs coef, int precision);
bool IsReal(number a, coeffs coef);
bool ImagGreaterZero(number a, coeffs coef);
number squareroot(number a, coeffs coef, int iteration);// iteration in Heron algorithm


#endif