#ifndef LATTICE_HPP
#define LATTICE_HPP


//"A course in computational algebraic numbertheory" by Henry Cohen
class lattice {
    
    private:
        
        //basis of lattice, spanned by columns
        bigintmat *basis;
        
        coeffs coef;
        
        //Array of B_i, see Cohen
        number *B; 
        
        //constant for LLL
        number c; 
        
        //LLL-basis
        bigintmat *b; 
        
        bigintmat *b_star; //for LLL
        bigintmat *H; //for LLL
        bigintmat *my; //for LLL
        
//         int precision;
        bool relat;
//         bool Gram, indep, relat;// integral;
        
   
        
        //for LLL, see Cohen
        inline void RED(int k, int l);
        
        //for LLL, see Cohen
        inline void SWAP(int k, int k_max);
        
//         inline void INSERT(int k, int i);
      
//         //for MLLL
//         inline void SWAPG(int k, int k_max);
        
        
        
        inline bool gram_schmidt(int k);
        inline void gram_schmidt_MLLL(int k);
        inline bool gram_matrix(int k);

        
    public:
        
        //constructor creates new lattice spanned by columns of b
        lattice(bigintmat* basis); 
          
        //destructor
        ~lattice();
        
        bool LLL();
        void set_c(number a);
        bigintmat* get_basis();
        bigintmat* get_reduced_basis();

        
//         void is_Gram(bool a){Gram = a;};
//         void is_independent(bool a){indep = a;};
//         void is_integral(bool a){integral = a;};
//         void set_rel(bool a);
//         void set_lattice(bigintmat L);
// 
//         void set_precision(int a);
// 
//         get methods
// 
//         bigintmat get_relation();
// 
//         bool enumerate(number a, vector<bigintmat> / pointer? liste);
// 
//         bool MLLL();
//         bool L1();
//         bool L2();
        
};

#endif