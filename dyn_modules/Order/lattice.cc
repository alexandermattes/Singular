#include <libpolys/coeffs/bigintmat.h>
#include "lattice.h"
#include "kernel/febase.h"  // for Print, WerrorS
#include "libpolys/coeffs/numbers.h" 
#include "libpolys/coeffs/coeffs.h"
#include "Singular/ipid.h"
#include <iostream>
#include <vector>
#include <utility>
#include "kernel/numeric/mpr_numeric.h"
#include "libpolys/coeffs/gnumpc.cc"




///////////////////////////////////////
//         Debugging Stuff          ///
///////////////////////////////////////

// #define DEBUG_PRINTS 1 //remove this line to disable debugging
#ifdef DEBUG_PRINTS

  //DEBUG_BLOCK(true / false); to enable/disable debugging in this block
# define DEBUG_BLOCK(x) bool debug_block = x; 
  //standard if DEBUG_BLOCK( ); not used
  static bool debug_block = true;
  //use like: DEBUG_PRINT(("Status %d",i))
# define DEBUG_PRINT(x) do {if(debug_block) {Print x ;}} while(0)
# define DEBUG_CMD(x) do {if(debug_block) {x;}} while(0)
# define DEBUG_VAR(x) do {if(debug_block) {std::cout<<#x<<": "<<x<<std::endl;}}  while(0)
# define DEBUG_N(x) do {if(debug_block) {Print(#x);Print(": ");n_Print(x,coef);Print("\n");}} while(0)
# define DEBUG_BIM(x) do {if(debug_block) {Print(#x);Print(": ");x->Print();Print("\n");}} while(0)
#else
# define DEBUG_BLOCK(x) do {} while (0)
# define DEBUG_PRINT(x) do {} while (0)
# define DEBUG_CMD(x)   do {} while (0)
# define DEBUG_VAR(x)   do {} while (0)
# define DEBUG_N(x)     do {} while (0)
# define DEBUG_BIM(x)   do {} while (0)
#endif


///////////////////////////////////////
//     constructors/destructor      ///
///////////////////////////////////////
 
lattice::lattice(bigintmat* inputmatrix, bool use_as_gram_matrix){
    DEBUG_BLOCK(true);
    DEBUG_PRINT(("Creating new lattice\n"));
    
    out_coef = inputmatrix->basecoeffs();
//     DEBUG_PRINT(("always set coef = out_coef"));
    coef = out_coef;
            
    //NOTE: Add transformation from rings to fields here
        //         exact <-> rounded
    //     if(nCoeff_is_Ring_Z(coef)) {
    //         integral = true;        
    //     }
    
    if(use_as_gram_matrix) {
        DEBUG_PRINT(("Input is gram matrix\n"));
        
        only_gram_matrix_given = true;
        
        if(inputmatrix->cols() != inputmatrix->rows()) {
                Werror("gram matrix not square");
        }
        
        n = inputmatrix->cols();
        m = n;
        
        
        gram_matrix_content = new number*[n];
        for(int i = 0; i < n; i++) {
            gram_matrix_content[i] = new number[i+1];
            for(int j = 0; j<=i; j++) {
                gram_matrix_content[i][j] = NULL;
            }
        } 
        for(int i = 1; i <= n; i++) {
            for(int j = 1; j <= i; j++) {
                gram_matrix_rawset(i,j,inputmatrix->get(i,j));
            }       
        }      
        
        basis = NULL;
    } else {
        DEBUG_PRINT(("Input is basismatrix\n"));
        
        only_gram_matrix_given = false;
        
        n = inputmatrix->cols();
        m = inputmatrix->rows();       

        basis = new bigintmat*[n+1]; //basis[0] is not used
        basis[0] = NULL;
        for(int i=1; i<=n; i++) {
            basis[i] = new bigintmat(m,1,coef);
            inputmatrix->getcol(i,basis[i]);
        }
        
        gram_matrix_content = NULL;
    }
    
    c = NULL;
    B = NULL;
    b = NULL;
    b_star = NULL;
    H = NULL;
    my = NULL;
    Q = NULL;
    d = NULL;
    lambda = NULL;
    rank = -1;
    independentVectors = true; 
    integral = true;
    trans_matrix = true;
    
    //for enumeration
    x = NULL;
    bound = NULL;
    
    DEBUG_PRINT(("Done\n"));
}

lattice::~lattice() {
    DEBUG_BLOCK(true);
    DEBUG_PRINT(("Deleting lattice..."));
    
    if(basis != NULL) {
        for(int i=0; i<=n; i++) {
            delete basis[i];
        }
        delete[] basis;
    }
    
    if(b != NULL) {
        for(int i=0; i<=n; i++) {
            delete b[i];
        }
        delete[] b;
    }
    
    if(b_star != NULL) {
        for(int i=0; i<=n; i++) {
            delete b_star[i];
        }
        delete[] b_star;
    }
    
    delete H;
    
    delete my;
    
    delete[] B;
    
    n_Delete(&c,coef);
    
    delete Q;
        
    //NOTE: add other member variables later
    
    DEBUG_PRINT(("Done\n"));
}


///////////////////////////////////////
//               LLL                ///
///////////////////////////////////////

bool lattice::LLL(){
    // c = 3/4
    number three = n_Init(3, coef);
    number four  = n_Init(4, coef);
    number default_c = n_Div(three,four,coef);
    n_Delete(&three,coef);
    n_Delete(&four,coef);
    bool integral = false;
    if(nCoeff_is_Ring_Z(coef)) {
        integral = true;
    }
    return lattice::LLL(default_c,integral,true,false);
}


bool lattice::LLL(number& c, bool trans_matrix, bool integral, bool independentVectors){    
    DEBUG_PRINT(("Start LLL\n"));  
    
    if(c==NULL) {
        number three = n_Init(3, coef);
        number four  = n_Init(4, coef);
        number default_c = n_Div(three,four,coef);
        n_Delete(&three,coef);
        n_Delete(&four,coef);
        c = default_c;
    }
    
    this->trans_matrix = trans_matrix;
    this->integral = integral;
    this->independentVectors = independentVectors;
    
    DEBUG_PRINT(("Delete old Results\n"));
    
    if(!only_gram_matrix_given){
        if(b == NULL) {
            b = new bigintmat*[n+1]; //b[0] is not used
            b[0] = NULL;
            for(int j=1; j<=n; j++) {
                b[j] = bimCopy(basis[j]);
            }
        } else {
            b[0] = NULL;
            for(int j=1; j<=n; j++) {
                delete b[j];
                b[j] = bimCopy(basis[j]);
            }
        }
        
        if(b_star == NULL) {
            b_star = new bigintmat*[n+1]; //b_star[0] is not used
            for(int j=0; j<=n; j++) {
                b_star[j] = NULL;
            }
        } else {
            for(int j=0; j<=n; j++) {
                delete b_star[j];
                b_star[j] = NULL;
            }
        }
    }
    
    
    if(B == NULL) {
        B = new number[n+1]; //B[0] is not used   
    }
    for(int j=0; j<=n; j++) {
        B[j] = n_Init(0,coef);
    }
    
    if(d == NULL) {
        d = new number[n+1]; //B[0] is not used   
    }
    for(int j=0; j<=n; j++) {
        d[j] = n_Init(0,coef);
    }
    
    delete H;
    if(trans_matrix) {
        H = new bigintmat(n,n,coef);
    } else {
        H = NULL;
    }
    
    delete my;
    my = new bigintmat(m,n,coef);
    
    delete lambda;
    lambda = new bigintmat(m,n,coef);
    
    delete[] d;
    if(integral) {
        d = new number[n+1];
    } else {
        d = NULL;
    }
    
    this->c = c; //NOTE: make copy and check coeffs
    
    
    
    if (n == 0) { //NOTE: necessary?
        return true;
    }
    if (n == 1) {
        return false;
    }
    
    
    
    DEBUG_PRINT(("Initialize\n"));
    int k = 2;
    int k_max = 1;
    
    if(integral) {
        d[0] = n_Init(1,coef);
        if(only_gram_matrix_given) {
            d[1] = n_Copy(gram_matrix_view(1,1),coef);
        } else {
            d[1] = scalarproduct(b[1],b[1]);
        }
    } else {
        if(only_gram_matrix_given) {
            B[1] = n_Copy(gram_matrix_view(1,1),coef);
        } else {
            b_star[1] = bimCopy(b[1]);
            B[1] = scalarproduct(b[1],b[1]);
        }
    }
    
    if(trans_matrix) {
        H->one();
    }
    
    do{
        DEBUG_PRINT(("Incremental Gram-Schmidt\n"));
        DEBUG_VAR(k);
        DEBUG_VAR(k_max);
        if(k > k_max){
            k_max = k;
            
            if(integral) {
                if(independentVectors) {
                    if(gram_schmidt_integral(k)) {
                        return true; //"did not form a basis"
                    }                //NOTE: better error handling?
                } else {
                    gram_schmidt_MLLL_integral(k);
                }
            } else {
                if(independentVectors) {
                    if(gram_schmidt(k)){
                        return true; //"did not form a basis"
                    }                //NOTE: better error handling?
                } else {
                    gram_schmidt_MLLL(k);
                }
            }
        }
        
        while(true){
            DEBUG_PRINT(("Test LLL condition\n"));
            if(integral) {
                REDI(k,k-1);
                                
                number leftside = n_Mult(d[k], d[k-2], coef);
                number d_kminusone_squared = n_Mult(d[k-1], d[k-1], coef);
                number c_times_stuff = n_Mult(c,d_kminusone_squared,coef);
                number lambda_k_kminusone_squared = n_Mult(lambda->view(k,k-1), lambda->view(k,k-1), coef);
                number rightside = n_Sub(c_times_stuff, lambda_k_kminusone_squared, coef);
                if(n_Greater(rightside,leftside,coef)){
                    if(independentVectors) {
                        SWAPI(k,k_max);
                    } else {
                        SWAPG_integral(k,k_max);
                    }
                    if(k>2){
                        k--;
                    }
                } else {
                    for(int l=k-2; l>0; l--){
                        REDI(k,l);
                    }
                    k++;
                    break;
                }
            } else {
                RED(k,k-1);
                // if((B[k] < (c- my*my)*B[k-1]))
                if(n_Greater(n_Mult(n_Sub(c, n_Mult(my->view(k,k-1), my->view(k,k-1), coef), coef), B[k-1], coef), B[k], coef)){
                    if(independentVectors) {
                        SWAP(k,k_max);
                    } else {
                        SWAPG(k,k_max);
                    }
                    if(k>2){
                        k--;
                    }
                } else {
                    for(int l=k-2; l>0; l--){
                        RED(k,l);
                    }
                    k++;
                    break;
                }
            }
        }
    } while(k <= n);
    
    if(!only_gram_matrix_given) {
        rank = n;
        for(int i=1; b[i]->isZero() && i<=n; i++) {
            rank--;
        }
    } else {
        //NOTE:...
    }
    
        
        
    DEBUG_PRINT(("End of LLL\n"));
    return false;
}

void lattice::RED(int k, int l){
    DEBUG_PRINT(("Start RED with k=%d and l=%d\n",k,l));

    number n_1div2    = n_Div(n_Init( 1,coef),n_Init(2,coef),coef);
    number n_neg1div2 = n_Div(n_Init(-1,coef),n_Init(2,coef),coef);
    number my_kl = my->get(k,l);
    DEBUG_N(my_kl);
    
    // if(|my_kl| > 1/2)
    if(n_Greater (my_kl,n_1div2,coef) || n_Greater (n_neg1div2,my_kl,coef)) { 
        
        number my_klplus1div2;
        if(n_Greater (my_kl,n_Init(0,coef),coef)) {
            my_klplus1div2 = n_Add(my_kl, n_1div2, coef);
        } else {
            my_klplus1div2 = n_Add(my_kl, n_neg1div2, coef);
        }
                
        number q = n_IntDiv(n_GetNumerator(my_klplus1div2,coef),n_GetDenom(my_klplus1div2,coef),coef); 
        
        DEBUG_N(q);
        
        if(only_gram_matrix_given){
            for(int i=1; i<=n; i++){
                gram_matrix_rawset(i,k,n_Sub(gram_matrix_view(i,k),n_Mult(gram_matrix_view(i,l),q,coef),coef));
            }
        } else {
            b[k]->sub(bimMult(b[l],q,coef));
        }

        if(trans_matrix) {
            H->addcol(k,l,n_Mult(q,n_Init(-1,coef),coef),coef);
        }
        
        my->set(k,l,n_Sub(my->view(k,l),q,coef),coef);
        
        for(int i=1;i<=l-1;i++){
            my->set(k,i,n_Sub(my->view(k,i), n_Mult(q, my->view(l,i),coef), coef), coef);
        }
    }
    DEBUG_PRINT(("End of RED\n"));
}

void lattice::REDI(int k, int l){
    DEBUG_PRINT(("Start REDI with k=%d and l=%d\n",k,l));
    
    
    number two_lambda = n_Mult(n_Init(2,coef),lambda->view(k,l),coef);
    number minus_d_l = n_Mult(d[l],n_Init(-1,coef),coef);
    if(n_Greater(two_lambda,d[l],coef)||n_Greater(minus_d_l,two_lambda,coef)) {
        
        number dividend = n_Add(two_lambda,d[l],coef);
        number divisor = n_Mult(n_Init(2,coef),d[l],coef);
        DEBUG_N(dividend);
        DEBUG_N(divisor);
        
        number q = n_IntDiv(dividend,divisor,coef);
        
        DEBUG_N(q);
        
        if(trans_matrix) {
            H->addcol(k,l,n_Mult(q,n_Init(-1,coef),coef),coef);
        }
        
        if(only_gram_matrix_given){
            for(int i=1; i<=n; i++){
                gram_matrix_rawset(i,k,n_Sub(gram_matrix_view(i,k),n_Mult(gram_matrix_view(i,l),q,coef),coef));
            }
        } else {
            b[k]->sub(bimMult(b[l],q,coef));
        }
        
        lambda->set(k,l,n_Sub(lambda->view(k,l), n_Mult(q, d[l],coef),coef),coef);
        
        for(int i=1;i<=l-1;i++){
            lambda->set(k,i,n_Sub(lambda->view(k,i), n_Mult(q, lambda->view(l,i),coef), coef), coef);
        }
    }
    
    DEBUG_PRINT(("End of REDI\n"));
}

void lattice::SWAP(int k, int k_max){
    DEBUG_PRINT(("Start SWAP with k=%d and k_max=%d\n",k,k_max));   
    
    
    if(only_gram_matrix_given) {
        for(int i=1; i<=n; i++){
            number temp = n_Copy(gram_matrix_view(i,k),coef);
            gram_matrix_rawset(i,k,n_Copy(gram_matrix_view(i,k-1),coef));
            gram_matrix_rawset(i,k-1,temp);
        }
    } else {
        bigintmat * temp = b[k];
        b[k] = b[k-1];
        b[k-1] = temp;
    }
    
    if(trans_matrix) {
        H->swap(k,k-1);
    }
    
    for(int j = 1; j <= k-2; j++){
        number my_kj = my->get(k,j);
        my->set(k,j,my->view(k-1,j),coef);
        my->set(k-1,j,my_kj,coef);
    }
    
    number my_ = my->get(k,k-1);
    
    number B_ = n_Add(B[k], n_Mult(n_Mult(my_,my_,coef), B[k-1], coef), coef); 
    
    my->set(k,k-1,n_Div(n_Mult(my_, B[k-1], coef), B_, coef), coef);
    
    
    if(!only_gram_matrix_given) {
        bigintmat * b_ = bimCopy(b_star[k-1]);
        
        b_star[k-1] = bimAdd(b_star[k], bimMult(b_, my_, coef));
        
        b_star[k] = bimSub(bimMult(b_, n_Div(B[k], B_, coef),coef), bimMult( b_star[k], my->view(k,k-1), coef));

        delete b_;
    }
    
    
    
    B[k] = n_Div(n_Mult(B[k], B[k-1], coef), B_, coef);
    
    B[k-1] = n_Copy(B_, coef);
    for(int i = k+1; i <= k_max; i++){
        number t = my->get(i,k);
        my->set(i,k,n_Sub(my->get(i,k-1), n_Mult(my_, t, coef), coef), coef);
        my->set(i,k-1, n_Add(t, n_Mult(my->view(k,k-1), my->view(i,k), coef), coef), coef);
    }
    DEBUG_PRINT(("End of SWAP\n"));
}

void lattice::SWAPI(int k, int k_max){
    DEBUG_PRINT(("Start SWAPI with k=%d and k_max=%d\n",k,k_max));   
    
    if(trans_matrix) {
        H->swap(k,k-1);
    }
    
    if(only_gram_matrix_given) {
        for(int i=1; i<=n; i++){
            number temp = n_Copy(gram_matrix_view(i,k),coef);
            gram_matrix_rawset(i,k,n_Copy(gram_matrix_view(i,k-1),coef));
            gram_matrix_rawset(i,k-1,temp);
        }
    } else {
        bigintmat * temp = b[k];
        b[k] = b[k-1];
        b[k-1] = temp;
    }
    
    for(int j = 1; j <= k-2; j++){
        number lambda_kj = lambda->get(k,j);
        lambda->set(k,j,lambda->view(k-1,j),coef);
        lambda->set(k-1,j,lambda_kj,coef);
    }
    
    number lambda_ = lambda->get(k,k-1);
    
    number B = n_Div(n_Add(n_Mult(d[k-2],d[k],coef),n_Mult(lambda_,lambda_,coef),coef),d[k-1],coef);
    
    for(int i = k+1; i <= k_max; i++){
        number t = lambda->get(i,k);
        
        lambda->set(i,k,n_Div(n_Sub(n_Mult(d[k],lambda->view(i,k-1),coef),n_Mult(lambda_,t,coef),coef),d[k-1],coef),coef);
        
        lambda->set(i,k-1,n_Div(n_Add(n_Mult(B,t,coef),n_Mult(lambda_,lambda->view(i,k),coef),coef),d[k],coef),coef);
        
    }
    
    d[k-1] = B;
    
    DEBUG_PRINT(("End of SWAPI\n"));
}

void lattice::SWAPG(int k, int k_max){
    DEBUG_PRINT(("Start SWAPG with k=%d and k_max=%d\n",k,k_max));   
    
    if(only_gram_matrix_given) {
        for(int i=1; i<=n; i++){
            number temp = n_Copy(gram_matrix_view(i,k),coef);
            gram_matrix_rawset(i,k,n_Copy(gram_matrix_view(i,k-1),coef));
            gram_matrix_rawset(i,k-1,temp);
        }
    } else {
        bigintmat * temp = b[k];
        b[k] = b[k-1];
        b[k-1] = temp;
    }
    
    if(trans_matrix) {
        H->swap(k,k-1);
    }
    
    for(int j = 1; j <= k-2; j++){
        number my_kj = my->get(k,j);
        my->set(k,j,my->get(k-1,j),coef);
        my->set(k-1,j,my_kj,coef);
    }
    
    number my_ = my->get(k,k-1);
    
    number B_ = n_Add(B[k], n_Mult(n_Mult(my_,my_,coef), B[k-1], coef), coef); 
    
    if(n_IsZero(B[k],coef)) {
        if(n_IsZero(my_,coef)) {
            number tempnumber = B[k];
            B[k] = B[k-1];
            B[k-1] = tempnumber;
            
            if(!only_gram_matrix_given) {
                bigintmat * temp = b_star[k];
                b_star[k] = b_star[k-1];
                b_star[k-1] = temp;
            }
            
            for(int i = k+1; i <= k_max; i++){
                number tempnumber = my->get(i,k);
                my->set(i,k,my->get(i,k-1), coef);
                my->set(i,k-1,tempnumber, coef);
            }
        } else {
            B[k-1] = n_Copy(B_, coef); //delete B[k-1] ?
            
            if(!only_gram_matrix_given) {
                b_star[k-1]->skalmult(my_, coef);
            }
            
            my->set(k,k-1, n_Div(n_Init(1,coef),my_, coef), coef);
            
            for(int i = k+1; i <= k_max; i++){
                my->set(i,k-1,n_Div(my->view(i,k-1), my_, coef), coef);
            }
        }
    } else {
        number t = n_Div(B[k-1], B_, coef);
        
        my->set(k,k-1,n_Mult(my_,t,coef),coef);
        
        if(!only_gram_matrix_given) {
            bigintmat * b_ = b_star[k-1];
                    
            b_star[k-1] = bimAdd(b_star[k], bimMult(b_, my_, coef));
            
            b_star[k] = bimSub(bimMult(b_, n_Div(B[k], B_, coef),coef), bimMult( b_star[k], my->view(k,k-1), coef));
            
            delete b_;
        }
        
        B[k] = n_Mult(B[k],t,coef); // n_InpMult
        
        B[k-1] = n_Copy(B_, coef);
        
        for(int i = k+1; i <= k_max; i++){
            t = my->get(i,k);
            my->set(i,k,n_Sub(my->view(i,k-1),n_Mult(my_,t,coef), coef), coef);
            my->set(i,k-1,n_Add(t,n_Mult(my->view(k,k-1),my->view(i,k),coef),coef),coef);
        }
    }
    DEBUG_PRINT(("End of SWAPG\n"));
}

void lattice::SWAPG_integral(int k, int k_max){
    DEBUG_PRINT(("Start SWAPI_integral with k=%d and k_max=%d\n",k,k_max));   
    
    Werror("integral MLLL not implemented yet\n");
    
    DEBUG_PRINT(("End of SWAPI_integral\n"));
}

bool lattice::gram_schmidt(int k) {
    DEBUG_PRINT(("Start gram_schmidt(%d)\n",k));
    
    if(only_gram_matrix_given) { //see Cohen 2.6.3 Remark (2)
        for(int j=1; j<k; j++) {
            number sum = n_Init(0,coef);
            for(int i=1; i<k; i++) {
                n_InpAdd(sum,n_Mult(n_Mult(my->view(j,i),my->view(k,i),coef),B[i],coef),coef);
            }
            my->set(k,j,n_Mult(n_Sub(gram_matrix_view(k,j),sum,coef),B[j],coef),coef);
            n_Delete(&sum,coef);
        }
        
        number sum = n_Init(0,coef);
        for(int i=1; i<k; i++) {
            n_InpAdd(sum,n_Mult(n_Mult(my->view(k,i),my->view(k,i),coef),B[i],coef),coef);
        }
        B[k] = n_Sub(gram_matrix_view(k,k),sum,coef);
        n_Delete(&sum,coef);
    } else {
        delete b_star[k];
        b_star[k] = bimCopy(b[k]);  
        
        for(int j=1; j<k; j++) {
            my->set(k,j,n_Div(scalarproduct(b[k],b_star[j]),B[j],coef),coef);
            
            b_star[k]->sub(bimMult(b_star[j],my->view(k,j),coef));
        }
        
        B[k] = scalarproduct(b_star[k],b_star[k]);
            
        
    }
    if(n_IsZero(B[k],coef)){
        Werror("did not form a basis\n");
        DEBUG_PRINT(("End of gram_schmidt(%d)\n",k));
        return true;
    } else {
        DEBUG_PRINT(("End of gram_schmidt(%d)\n",k));
        return false;
    }
}

bool lattice::gram_schmidt_integral(int k) {
    DEBUG_PRINT(("Start gram_schmidt_integral(%d)\n",k));
        
    for(int j=1; j<=k; j++) {
        
        number u = scalarproduct(b[k],b[j]);
        
        for(int i=1; i<j; i++) {
            u = n_Div(n_Sub(n_Mult(d[i],u,coef),n_Mult(lambda->view(k,i),lambda->view(j,i),coef),coef),d[i-1],coef);
        }
        
        if(j<k) {
             lambda->set(k,j,u,coef);
        } else {
            d[k] = n_Copy(u,coef);
        }
    }
       
    if(n_IsZero(d[k],coef)){
        Werror("did not form a basis\n");
        DEBUG_PRINT(("End of gram_schmidt_integral(%d)\n",k));
        return true;
    } else {
        DEBUG_PRINT(("End of gram_schmidt_integral(%d)\n",k));
        return false;
    }
}

void lattice::gram_schmidt_MLLL(int k) {
    DEBUG_PRINT(("Start gram_schmidt_MLLL(%d)\n",k));
    
    if(only_gram_matrix_given) { //see Cohen 2.6.3 Remark (2)  
        for(int j=1; j<k; j++) {
            if(n_IsZero(B[j],coef)) {
                my->set(k,j,n_Init(0,coef));
            } else {
                number sum = n_Init(0,coef);
                for(int i=1; i<k; i++) {
                    n_InpAdd(sum,n_Mult(n_Mult(my->view(j,i),my->view(k,i),coef),B[i],coef),coef);
                }
                my->set(k,j,n_Mult(n_Sub(gram_matrix_view(k,j),sum,coef),B[j],coef),coef);
                n_Delete(&sum,coef);
            }
        }
        
        number sum = n_Init(0,coef);
        for(int i=1; i<k; i++) {
            n_InpAdd(sum,n_Mult(n_Mult(my->view(k,i),my->view(k,i),coef),B[i],coef),coef);
        }
        B[k] = n_Sub(gram_matrix_view(k,k),sum,coef);
        n_Delete(&sum,coef);
    } else {
        for(int j=1; j<k; j++) {
            if(n_IsZero(B[j],coef)) {
                my->set(k,j,n_Init(0,coef));
            } else {
                my->set(k,j,n_Div(scalarproduct(b[k],b_star[j]),B[j],coef),coef);
            }
        }
        
        delete b_star[k];
        b_star[k] = bimCopy(b[k]);
        for(int j=1; j<k; j++) {
            b_star[k]->sub(bimMult(b_star[j],my->view(k,j),coef));
        }
        
        B[k] = scalarproduct(b_star[k],b_star[k]);
    }
    DEBUG_PRINT(("End of gram_schmidt_MLLL(%d)\n",k));
}

void lattice::gram_schmidt_MLLL_integral(int k) {
    DEBUG_PRINT(("Start gram_schmidt_MLLL_integral(%d)\n",k));
    
    Werror("integral MLLL not implemented yet\n");
    
    DEBUG_PRINT(("End of gram_schmidt_MLLL_integral(%d)\n",k));
}

// bool lattice::gram_matrix(int k){ Page 89 Remark 2
//     number* a = new number[k];
//     for(int j = 1; j<k;j++){
//         a[j] = n_Init(0,coef);
//         for(int i =1; i<=b->rows(); i++){
//             a[j] = n_Add(a[j],n_Mult(b->view(i,k),b->view(i,j),coef),coef);//a[j] += b->view(i,k) * b->view(i,j);
//         }
//         for(int i =1; i<=j-1; i++){
//             a[j] = n_Add(a[j],n_Mult(b->view(i,j),a[i],coef),coef);//a[j] += my->view(j,i) * a[i];
//         }
//         my->set(k,j,n_Div(a[j],B[j],coef),coef);
//     }
//     B[k]=n_Init(0,coef);
//     for(int i =1; i<=b->rows(); i++){
//         B[k]=n_Add(B[k],n_Mult(b->view(i,k),b->view(i,k),coef),coef);//B[k] += b->view(i,k) * b->view(i,k);
//     }
//     for(int i =1; i<=k-1; i++){
//         B[k] = n_Add(B[k],n_Mult(my->view(k,i),a[i],coef),coef);//B[k] += my->view(k,i) * a[i];
//     }
//     if(B[k] == 0){
//         Werror("did not form a basis\n");
//         return false;
//     }
//     return true;
// }


///////////////////////////////////////
//             Enumerate            ///
///////////////////////////////////////
//Public
bigintmat * lattice::enumerate_all(number a){
    //Quadratic Supplement
    //DEBUG_BLOCK(true);
    DEBUG_PRINT(("Start enumerate_all\n"));
    DEBUG_PRINT(("check input\n"));
    if(!n_GreaterZero(a,out_coef)){
        if(n_IsZero(a,out_coef) && n==m){
            return new bigintmat(m,1,out_coef);
        } else {
            DEBUG_PRINT(("negative input\n"));
            return NULL;
        }
    }
    if( Q == NULL){
        if(quadratic_supplement()){
            return NULL;
        }
    }
    DEBUG_PRINT(("Q defined\n"));
    //Q->Print();PrintS("\n");
    
    //usefull numbers
    number zero = n_Init(0,coef);
    number minusOne = n_Init(-1,coef);
    
    //backtracking for elements
    //DEBUG_BLOCK(true);
    DEBUG_PRINT(("Start backtracking\n"));
    DEBUG_PRINT(("Initialize vector and other variables\n"));
    std::vector<std::pair<number,bigintmat*> > elementsvector;
    elementsvector.push_back( std::make_pair(zero, new bigintmat(m,1,out_coef)));
    if( x != NULL){
        delete x;
        x=NULL;
    }
    x = new bigintmat(m,1,coef);
    //x->Print();PrintS("\n");
    DEBUG_PRINT(("map a\n"));
    if(bound != NULL){
        delete bound;
        bound = NULL;
    }
    bound = new number[n+1];
    nMapFunc f = n_SetMap(out_coef, coef);
    bound[1] = f(a, out_coef, coef);//map a to coef
    //n_Print(bound[1],coef);PrintS("\n");
    DEBUG_PRINT(("set bound\n"));
    for(int i = 2; i<n+1; i++){
        bound[i] = n_Copy(bound[1],coef);
        //n_Print(bound[i],coef);PrintS("\n");
    }
    DEBUG_PRINT(("find element\n"));
    //bigintmat* elements = enumerate_next(a);
    increase_x(1);
    number check = enumerate_get_next();
    while(!n_Equal(minusOne,check,coef)){
        //append x to elements
        DEBUG_PRINT(("new element to list\n"));
        //elements->appendCol(bimChangeout_coeff(x,out_coef));
        check = n_Sub(bound[1],check,coef);
        check = n_Sub(bound[n],check,coef);
        elementsvector.push_back(std::make_pair(n_Copy(check,coef),bimCopy(x)));
        //n_Print(elementsvector[elementsvector.size()-1].first,coef);PrintS("\n");
        for(unsigned i=1; i<elementsvector.size();i++){
            if(n_Greater(elementsvector[i].first,check,coef)){
                elementsvector.pop_back();
                elementsvector.insert(elementsvector.begin()+i,std::make_pair(n_Copy(check,coef),bimCopy(x)));
                DEBUG_VAR(elementsvector.size());
                break;
            }
        }
        if(elementsvector.size() >= 1000000){
            elementsvector.pop_back();
        }
        increase_x(1);
        check = enumerate_get_next();
    }
    DEBUG_PRINT(("generate bigintmat for return\n"));
    bigintmat* elements = new bigintmat(m,1,out_coef);
    if(b!=NULL){
        for(unsigned i=1; i<elementsvector.size();i++){
            elements->appendCol(bimChangeCoeff(elementsvector[i].second,out_coef));
        }
    } else {
        for(unsigned i=1; i<elementsvector.size();i++){
           elements->appendCol(bimChangeCoeff(elementsvector[i].second,out_coef));
        }
    }
    delete bound;
    bound = NULL;
    delete x;
    x = NULL;
    return elements;
}

bigintmat * lattice::enumerate_next(number a, bigintmat * in){//next element x with norm(x)<a and x=in possible
    //DEBUG_BLOCK(true);
    DEBUG_PRINT(("Start enumerate_next number and bigintmat\n"));
    if (in == NULL || in->rows()!=n || in->cols()!=1){
        DEBUG_PRINT(("Dimension error of input\n"));
        return NULL;
    }
    
    if(!n_GreaterZero(a,out_coef) || n_IsZero(a,out_coef) ){
        DEBUG_PRINT(("negative input\n"));
        return NULL;
    }
    
    DEBUG_PRINT(("check quadratic\n"));
    
    if( Q == NULL){
        if(!quadratic_supplement()){
            return NULL;
        }
    }
    DEBUG_PRINT(("Q defined\n"));
    //Q->Print();PrintS("\n");
    
    //usefull numbers
    number check;
    
    //backtracking for elements
    //DEBUG_BLOCK(true);
    DEBUG_PRINT(("Start backtracking\n"));
    DEBUG_PRINT(("Initialize variables\n"));
    if( x != NULL){
        delete x;
        x=NULL;
    }
    x = bimChangeCoeff(in,coef);
    if( bound != NULL){
        delete bound;
        bound=NULL;
    }
    bound = new number[n+1];
    DEBUG_PRINT(("set bound\n"));
    nMapFunc f = n_SetMap(out_coef, coef);
    bound[n] = f(a, out_coef, coef);//map a to coef
    for(int j = n; j>1; j--){
        check = check_bound(j);
        bound[j-1] = n_Sub(bound[j],check,coef);
        //n_Delete(&check, coef);
    }
    number minusOne = n_Init(-1,coef);
    DEBUG_PRINT(("find element\n"));
    number norm = enumerate_get_next();
    DEBUG_PRINT(("generate bigintmat for return\n"));
    if(n_Equal(minusOne,norm,coef)){
        return NULL;
    }
    bigintmat * out = bimChangeCoeff(x,out_coef);
    return out;
}

bigintmat * lattice::enumerate_next(number a){
    //DEBUG_BLOCK(true);
    DEBUG_PRINT(("Start enumerate_next number\n"));
    bigintmat * in =new bigintmat(m,1,out_coef);
    if(x == NULL){
        in->set(1,1,n_Init(1,out_coef),out_coef);
    } else {
        in = bimChangeCoeff(x,out_coef);
    }
    return enumerate_next(a,in);
}

bigintmat * lattice::enumerate_next(bigintmat * in){
    //DEBUG_BLOCK(true);
    DEBUG_PRINT(("Start enumerate_next bigintmat\n"));
    if(bound == NULL){
        Werror("no bound for elements given\n");
        return NULL;
    }
    if (in == NULL || in->rows()!=n || in->cols()!=1){
        DEBUG_PRINT(("Dimension error of input\n"));
        return NULL;
    }
    number a = bound[n];
    return enumerate_next(a,in);
}

bigintmat * lattice::enumerate_next(){
    //DEBUG_BLOCK(true);
    DEBUG_PRINT(("enumerate_next\n"));
    if(Q == NULL){
        Werror("not initialized\n");
        return NULL;
    }
    if(bound == NULL || x == NULL){
        return NULL;
    }
    increase_x(1);
    number minusOne = n_Init(-1,coef);
    number one = n_Init(1,coef);
    DEBUG_PRINT(("find element\n"));
    number norm = enumerate_get_next();
    DEBUG_PRINT(("generate bigintmat for return\n"));
    if(n_Equal(minusOne,norm,coef)){
        if(!n_Equal(minusOne, x->view(1,1),coef)){
            x->rawset(1,1, n_Add(one,x->view(1,1),coef),coef);
        }
        return NULL;
    }
    bigintmat * out = bimChangeCoeff(x,out_coef);
    return out;
}

//Private
number lattice::enumerate_get_next(){
    //DEBUG_BLOCK(true);
    DEBUG_PRINT(("enumerate_get_next\n"));
    number one = n_Init(1,coef);
    number zero = n_Init(0,coef);
    number minusOne = n_Init(-1,coef);
    int index =1;
    //x->Print();PrintS("\n");
    //DEBUG_PRINT(("first time changing x\n"));
    //increase_x(1);
    DEBUG_PRINT(("actual backtracking\n"));
    while (index <= m) {
        DEBUG_PRINT(("update check\n"));
        number check = check_bound(index);
        DEBUG_PRINT(("check check\n"));
        if (n_Greater(check,bound[index],coef)){
            DEBUG_PRINT(("element to great\n"));
            if(!(n_GreaterZero(x->view(index,1),coef) || n_IsZero(x->view(index,1),coef))){
                bound[index] = n_Init(0,coef);
                x->rawset(index,1,n_Init(0,coef),coef);
                index++;
                if(index<= m){
                    increase_x(index);
                }
            } else {
                if(index == n){
                    return minusOne;
                }
                x->rawset(index,1,minusOne,coef);
            }
        } else if(index == 1){
            DEBUG_PRINT(("possible new element\n"));
            if(n_IsZero(x->view(n,1),coef)){
                int j=n-1;
                while(n_IsZero(x->view(j,1),coef)){
                    j--;
                }
                if(n_GreaterZero(x->view(j,1),coef)){
                    return check;
                } else {
                    index = j+1;
                    for( j=1;j<index;j++){
                        x->set(j,1,zero,coef);
                    }
                    x->set(index,1,one,coef);
                }
            } else {
                return check;
            }
        } else {
            DEBUG_PRINT(("reduce index\n"));
            index--;
            bound[index] = n_Sub(bound[index+1],check,coef);
        }
    }
    return minusOne;
}

bool lattice::quadratic_supplement(){
    //DEBUG_BLOCK(true);
    delete Q;
    Q = NULL;
    if(n != m) {  //NOTE: rank?
        return true;
    }
    
    Q = get_gram_matrix();
    
    number zero = n_Init(0,coef);
    number mult;
    
    DEBUG_PRINT(("Begin Quadratic Suplement\n"));
    for(int i = 1; i<Q->cols();i++){
        if(n_IsZero( Q->view(i,i), coef)){
            DEBUG_PRINT(("matrix not positive definite\n"));
            delete Q;
            Q = NULL;
            return true;
        }
        for( int j=i+1; j<=Q->cols();j++){
            Q->rawset(j,i,Q->get(i,j),coef);
            Q->rawset(i,j,n_Div(Q->get(i,j),Q->view(i,i),coef),coef);
        }
        for(int m=i+1; m<=Q->rows();m++){
            for(int n=i+1; n<=Q->cols();n++){
                mult = n_Mult(Q->view(m,i),Q->view(i,n),coef);
                Q->rawset(m,n,n_Sub(Q->get(m,n),mult,coef),coef);
            }
        }
    }
    
    DEBUG_PRINT(("Set Zeros\n"));
    for(int i = 2; i<=Q->cols();i++){
        for(int j = 1; j<i;j++){
            Q->rawset(i,j,zero,coef);
        }
    }
    DEBUG_PRINT(("Test: matrix positive definite\n"));
    for(int i=1; i<=Q->cols();i++){
        if(!n_GreaterZero( Q->view(i,i), coef)){
            DEBUG_PRINT(("matrix not positive definite\n"));
            delete Q;
            Q = NULL;
            return true;
        }
    }
    return false;
}

void lattice::increase_x(int index){
    number one = n_Init(1,coef);
    if (n_GreaterZero(x->view(index,1),coef) || n_IsZero(x->view(index,1),coef)){
        x->rawset(index,1, n_Add(one,x->view(index,1),coef),coef); //x_i=x_i+1
    } else {
        x->rawset(index,1, n_Sub(x->view(index,1),one,coef),coef);//x_i=x_i-1
    }
}

number lattice::check_bound(int index){
    number check = n_Init(0,coef);
    for(int i=index + 1;i<=Q->cols();i++){
        number mult = n_Mult(x->view(i,1),Q->view(index,i),coef);
        n_InpAdd(check,mult,coef);
        n_Delete(&mult,coef);
    }
    n_InpAdd(check, x->view(index,1), coef);
    n_InpMult(check, check, coef);
    n_InpMult(check, Q->get(index,index), coef);
    return check;
}


///////////////////////////////////////
//               Getter             ///
///////////////////////////////////////

bigintmat * lattice::get_basis() {
    if(only_gram_matrix_given) {
        Werror("lattice is only defined by its gram matrix");
    }
    bigintmat * r = new bigintmat(m,n,coef);
    for(int j=1; j<=n; j++) {
        r->setcol(j,basis[j]);
    }
    return r;
}

bigintmat * lattice::get_reduced_basis() {
    if(b == NULL) {
        Werror("reduced basis needs to be calculated via LLL first");
    }
    bigintmat * r = new bigintmat(m,n,coef);
    for(int j=1; j<=n; j++) {
        r->setcol(j,b[j]);
    }
    return r;
}

bigintmat * lattice::get_transformation_matrix() {
    if(H == NULL) {
        Werror("transformation matrix needs to be calculated via LLL first");
    }
    return bimCopy(H);
}

bigintmat * lattice::get_gram_matrix() {
    bigintmat * r = new bigintmat(n,n,coef);
    if(only_gram_matrix_given) {
        for(int i = 1; i <= n; i++) {
            for(int j = 1; j <= n; j++) {
                r->set(i,j,gram_matrix_view(i,j));
            }       
        }   
    } else {
        for(int i = 1; i <= n; i++) {
            for(int j = 1; j <= n; j++) {
                r->set(i,j,scalarproduct(basis[i],basis[j]));
            }
        }
    }
    return r;
}


///////////////////////////////////////
//               Setter             ///
///////////////////////////////////////

// void lattice::set_c(number a){
//     if (n_Greater (n_Mult(a,n_Init(4,coef),coef),n_Init(1,coef),coef) && n_Greater (n_Init(1,coef),a,coef)) {//(1<4*a && a<1){
//         c = n_Copy(a, coef);
//     } else if(n_IsOne(a, coef)){
//         c = n_Copy(a, coef);
//         Werror("probably not in polynomial time\n");
//     } else {
//         Werror("not a possible value\n");
//     }
// }


///////////////////////////////////////
//               Other             ///
///////////////////////////////////////

void lattice::gram_matrix_rawset(int i, int j, number n){
    if(i<j) {
        gram_matrix_rawset(j, i, n);
        return;
    }
    n_Delete(&gram_matrix_content[i-1][j-1], coef);
    gram_matrix_content[i-1][j-1] = n;    
}

number lattice::gram_matrix_view(int i, int j){
    if(i<j) {
        return gram_matrix_view(j, i);
    }
    return gram_matrix_content[i-1][j-1];
}

number scalarproduct(bigintmat * a, bigintmat * b) {
    if(a->cols()!=1) {
        Werror("a->cols()!=1 in scalarproduct(a,b)\n");
        return NULL;
    }
    if(b->cols()!=1) {
        Werror("b->cols()!=1 in scalarproduct(a,b)\n");
        return NULL;
    }
    if(a->rows()!=b->rows()) {
        Werror("b->cols()!=1 in scalarproduct(a,b)\n");
        return NULL;
    }
    if(a->basecoeffs()!=b->basecoeffs()) {
        Werror("a->basecoeffs()!=b->basecoeffs() in scalarproduct(a,b)\n");
        return NULL;
    }

    coeffs coef = a->basecoeffs();
    number p = n_Init(0,coef);
    for(int i = 1; i <= b->rows(); i++){
        number prod = n_Mult(a->view(i,1), b->view(i,1), coef);
        n_InpAdd(p, prod, coef);
        n_Delete(&prod,coef);
    }
    return p;
}


///////////////////////////////////////
//         Minkowski map            ///
///////////////////////////////////////

bigintmat * minkowksi(bigintmat ** elementarray,int size_elementarray, number * poly,int deg, coeffs coef, int precision){
    DEBUG_BLOCK(true);
    DEBUG_PRINT(("Begin Minkowski map\n"));
    DEBUG_PRINT(("Input check\n"));
    if(elementarray == NULL || poly == NULL || coef != elementarray[0]->basecoeffs()){
        WerrorS("wrong input!\n");
        return NULL;
    }
    
    for(int i=1;i<size_elementarray;i++){
        if(coef != elementarray[0]->basecoeffs()){
            WerrorS("wrong input!\n");
            return NULL;
        }
    }
    
    //char = 0
    if ( !(nCoeff_is_Ring_Z(coef) || nCoeff_is_R(coef) || nCoeff_is_Q(coef) ||
            nCoeff_is_long_R(coef) || nCoeff_is_long_C(coef)) ){
        WerrorS("Ground field not implemented!\n");
        return NULL;
    }
    
    if(deg<2){
        WerrorS("degree of polynomial to small\n");
        return NULL;
    }
    //check and define precision for Q
    if(precision<6){
        precision = 6;
    }
    if ( !(nCoeff_is_R(coef) || nCoeff_is_long_R(coef) || nCoeff_is_long_C(coef)) ){
        setGMPFloatDigits( precision+6,precision+6);
    }
    
    DEBUG_PRINT(("find roots\n"));
    ring CurrOldRing = rCopy(currRing);//need to change currRing, because rootContainer uses the coef of it
    char* n[] = {(char*)"s"};
    ring newring = rDefault(coef, 1, n);
    rChangeCurrRing(newring);
    DEBUG_PRINT(("initialize rootContainer\n"));
    rootContainer * rootcont= new rootContainer();
    rootcont->fillContainer( poly, NULL, 1, deg, rootContainer::onepoly, 1 );///
    rootcont->solver( precision+12);
    int number_roots = rootcont ->getAnzRoots();
    if(number_roots != deg){
        WerrorS("something went wrong: \n\tnot all roots found\n");
        return NULL;
    }
    LongComplexInfo paramComp;
    paramComp.float_len = si_min (precision+6, 32767);
    paramComp.float_len2 = si_min (precision+8, 32767);
    paramComp.par_name=(const char*)"i";
    
    coeffs comp = nInitChar(n_long_C, &paramComp);
    DEBUG_PRINT(("find r1, r2 and save roots in array\n"));
    number* roots = new number[deg+1];
    number* complexroots = new number[deg+1];
    int r1 = 0;
    int r2 = 0;
    for(int j=0; j<deg; j++){
        number a = n_Copy((number)(rootcont->getRoot(j)),comp);
        if( IsReal(a,comp)){
            roots[r1] = n_Copy(a,comp);
            r1++;
        }else if(ImagGreaterZero(a, comp)){
            complexroots[r2] = n_Copy(a,comp);
            r2++;
        }
        n_Delete(&a,comp);
    }
    rChangeCurrRing(CurrOldRing);
    rDelete(newring);
    DEBUG_VAR(r1);
    DEBUG_VAR(r2);
    for(int j=0;j<r2;j++){
        roots[r1+j]= n_Copy(complexroots[j],comp);
        n_Delete(&complexroots[j],comp);
    }
    delete complexroots;
    DEBUG_PRINT(("map elementarray to complex\n"));
    bigintmat ** elements = new bigintmat*[size_elementarray];
    for(int i=0;i<size_elementarray;i++){
        elements[i] = bimChangeCoeff(elementarray[i],comp);
    }
    
    DEBUG_PRINT(("generate output matrix\n"));
    DEBUG_PRINT(("real part\n"));
    bigintmat * complexmat = new bigintmat(r1+2*r2,size_elementarray,comp);
    for(int i=1; i<= r1; i++){
        number pot = n_Init(1,comp);
        for(int l=1; l<= deg; l++){
            for(int j=0; j<size_elementarray;j++){
                number mult = n_Mult(pot,elements[j]->view(l,1),comp);
                complexmat->rawset(i,j+1,n_Add(complexmat->view(i,j+1),mult,comp),comp);
                n_Delete(&mult,comp);
            }
            pot = n_Mult(pot, roots[i-1],comp);
        }
        n_Delete(&pot,comp);
    }
    DEBUG_PRINT(("imaginary part\n"));
    number sqrt2 = n_Init(1,comp);
    if(r2>0){
        number two = n_Init(2,comp);
        number sqrt2 = squareroot(two,comp,precision+10);
        n_Delete(&two,comp);
        for(int i=1; i<= r2; i++){
            number pot = n_Init(1,comp);
            for(int l=1; l<= deg; l++){
                for(int j=0; j<size_elementarray;j++){
                    number mult = n_Mult(pot,elements[j]->view(l,1),comp);
                    complexmat->set(r1+2*i,j+1,n_Add(complexmat->view(r1+2*i,j+1),mult,comp),comp);
                    n_Delete(&mult,comp);
                }
                pot = n_Mult(pot, roots[r1+i-1],comp);
            }
            n_Delete(&pot,comp);
            for(int j=1;j<=size_elementarray;j++){
                complexmat->set(r1+2*i,j,n_Mult(complexmat->view(r1+2*i,j),sqrt2,comp),comp);
                complexmat->set(r1+2*i-1,j,ngcRePart(complexmat->view(r1+2*i,j),comp),comp);
                complexmat->set(r1+2*i,j,ngcImPart(complexmat->view(r1+2*i,j),comp),comp);
            }
        }
    }
    DEBUG_PRINT(("delete Variables\n"));
    for(int i=0;i<size_elementarray;i++){
        delete elements[i];
    }
    delete elements;
    for(int i=0;i<r1+r2;i++){
        n_Delete(&roots[i],comp);
    }
    delete roots;
    
    DEBUG_PRINT(("to real\n"));
    LongComplexInfo paramReal;
    paramReal.float_len = si_min (precision, 32767);
    paramReal.float_len2 = si_min (precision, 32767);
    paramComp.par_name=(const char*)"i";
    coeffs real = nInitChar(n_long_R, &paramReal);
    //setGMPFloatDigits( precision, precision);
    bigintmat * realmat = bimChangeCoeff(complexmat,real);
    delete complexmat;
    return realmat;
}

bool IsReal(number a, coeffs coef){ //Im(a)==0
    number imag = ngcImPart(a, coef);
    bool out = n_IsZero(imag,coef);
    n_Delete(&imag,coef);
    return out;
}

bool ImagGreaterZero(number a, coeffs coef){ //Im(a)>0
    number imag = ngcImPart(a, coef);
    bool out = n_GreaterZero(imag,coef);
    n_Delete(&imag,coef);
    return out;
}

number squareroot(number a, coeffs coef, int prec){
    if(n_IsZero(a,coef)){
        return a;
    }
    if(!n_GreaterZero(a,coef)){
        return n_Init(0,coef);
    }
    number two = n_Init(2,coef);
    number xn = n_Copy(a,coef);
    number xn1,xn2;
    for(int i=0;i<prec;i++){
        xn1 = n_Div(a,xn,coef);
        xn2 = n_Add(xn,xn1,coef);
        n_Delete(&xn,coef);
        xn = n_Div(xn2,two,coef);
        n_Delete(&xn1,coef);
        n_Delete(&xn2,coef);
    }
    return xn;
}
