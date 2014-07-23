#include <libpolys/coeffs/bigintmat.h>
#include "lattice.h"
#include "reporter/reporter.h"  // for Print, WerrorS
#include "libpolys/coeffs/numbers.h" 
#include "libpolys/coeffs/coeffs.h"
#include "Singular/ipid.h"
#include <iostream>
#include <vector>
#include <utility>
#include "kernel/numeric/mpr_numeric.h"
#include "libpolys/coeffs/gnumpc.cc"
#include "nforder.h"

#include <polys/monomials/p_polys.h>
#include <polys/clapsing.h>


#include <Singular/ipshell.h>




///////////////////////////////////////
//         Debugging Stuff          ///
///////////////////////////////////////

 #define DEBUG_PRINTS 0 //remove this line to disable debugging
#ifdef DEBUG_PRINTS

  //DEBUG_BLOCK(true / false); to enable/disable debugging in this block
# define DEBUG_BLOCK(x) bool debug_block = x; 
  //standard if DEBUG_BLOCK( ); not used
  static bool debug_block = false;
  //use like: DEBUG_PRINT(("Status %d",i))
# define DEBUG_PRINT(x) do {if(debug_block) {Print x ;}} while(0)
# define DEBUG_CMD(x) do {if(debug_block) {x;}} while(0)
# define DEBUG_VAR(x) do {if(debug_block) {std::cout<<#x<<": "<<x<<std::endl;}}  while(0)
# define DEBUG_N(x) do {if(debug_block) {Print(#x);Print(": ");n_Print(x,coef);Print("\n");}} while(0)
# define DEBUG_BIM(x) do {if(debug_block) {Print(#x);Print(": ");x->Print();Print("\n");}} while(0)
# define DEBUG_LLL(x) do {if(debug_block) {print_all_in_LLL();}} while(0)
#else
# define DEBUG_BLOCK(x) do {} while (0)
# define DEBUG_PRINT(x) do {} while (0)
# define DEBUG_CMD(x)   do {} while (0)
# define DEBUG_VAR(x)   do {} while (0)
# define DEBUG_N(x)     do {} while (0)
# define DEBUG_BIM(x)   do {} while (0)
# define DEBUG_LLL(x)   do {cancel_ratios_in_LLL();} while (0)
#endif


///////////////////////////////////////
//     constructors/destructor      ///
///////////////////////////////////////
 
lattice::lattice(bigintmat* inputmatrix, bool use_as_gram_matrix){
    DEBUG_BLOCK(true);
    DEBUG_PRINT(("Creating new lattice\n"));
    
    coef = inputmatrix->basecoeffs();
    
    
    if(nCoeff_is_Ring_Z(coef)) {
        fieldcoef = nInitChar(n_Q,NULL);  
    } else {
        fieldcoef = coef;
    }
            
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
    
    c      = NULL;
    B      = NULL;
    b      = NULL;
    b_star = NULL;
    H      = NULL;
    my     = NULL;
    Q      = NULL;
    d      = NULL;
    lambda = NULL;
    rank   = -1;
    
    independentVectors = true; 
    integral           = true;
    trans_matrix       = true;
    
    //for enumeration
    x     = NULL;
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
        basis = NULL;
    }
    
    if(gram_matrix_content != NULL) {
        for(int i = 0; i < n; i++) {
            delete[] gram_matrix_content[i];
        } 
        delete[] gram_matrix_content;
        gram_matrix_content = NULL;
    }
    
    delete_LLL_computations();
    
    delete Q;
    Q = NULL;
    
    delete x;
    x = NULL;
    
    if(bound != NULL) {
        for(int i=0; i<=n; i++) {
            n_Delete(&bound[i],fieldcoef);
        }
        delete[] bound;
        bound = NULL;
    }
        
    DEBUG_PRINT(("Done\n"));
}


///////////////////////////////////////
//      Write / String / Print      ///
///////////////////////////////////////

void lattice::Write() {
    StringAppend("Lattice:\nof dimension %d\n", n);
    
    if(only_gram_matrix_given) {
        StringAppend("defined by gram matrix:\n");
        bigintmat * gm = get_gram_matrix();
        gm->Write();
        delete gm;
        StringAppend("\n");
    } else {
        StringAppend("defined by basis:\n");
        for(int i=1; i<=n; i++) {
            StringAppend("%d: ", i);
            basis[i]->Write();
            StringAppendS("\n");
        }
    }
}

char * lattice::String() {
    StringSetS("");
    Write();
    return StringEndS();
}
// void lattice::Print() {
//     char * s = String();
//     PrintS(s);
//     PrintS("\n");
//     omFree(s);
// }


///////////////////////////////////////
//               LLL                ///
///////////////////////////////////////

void lattice::delete_LLL_computations(){
    n_Delete(&c,coef);
    
    if(b != NULL) {
        for(int i=0; i<=n; i++) {
            delete b[i];
        }
        delete[] b;
        b = NULL;
    }
    
    if(b_star != NULL) {
        for(int i=0; i<=n; i++) {
            delete b_star[i];
        }
        delete[] b_star;
        b_star = NULL;
    }
    
    if(B != NULL) {
        for(int i=0; i<=n; i++) {
            n_Delete(&B[i],coef);
        }
        delete[] B;
        B = NULL;
    }
    
    delete H;
    H = NULL;
    
    delete my;
    my = NULL;
     
    if(d != NULL) {
        for(int i=0; i<=n; i++) {
            n_Delete(&d[i],coef);
        }
        delete[] d;
        d = NULL;
    }
    
    delete lambda;
    lambda = NULL; 
}

void lattice::print_all_in_LLL(){
//    DEBUG_BLOCK(true);
    
    if(b != NULL) {
        for(int i=1; i<=n;i++){
            if(b[i] != NULL) {
                std::cout<<"b["<<i<<"]: ";
                b[i]->Print();
                Print("\n");
            }
        }
    }
    
    if(b_star != NULL) {
        for(int i=1; i<=n;i++){
            if(b_star[i] != NULL) {
                std::cout<<"b_star["<<i<<"]: ";
                b_star[i]->Print();
                Print("\n");
            }
        }
    }
    
    if(B != NULL) {
        for(int i=1; i<=n;i++){
            if(B[i] != NULL) {
                std::cout<<"B["<<i<<"]: ";
                n_Print(B[i],coef);
                Print("\n");
            }
        }  
    }
    
    if(d != NULL) {
        for(int i=1; i<=n;i++){
            if(d[i] != NULL) {
                std::cout<<"d["<<i<<"]: ";
                n_Print(d[i],coef);
                Print("\n");
            }
        }  
    }
    
    if(my != NULL) {
        Print("my: ");
        my->Print();
        Print("\n");
    }
    
    if(lambda != NULL){
        Print("lambda: ");
        lambda->Print();
        Print("\n");
    }
    
    if(H != NULL){
        Print("H: ");
        H->Print();
        Print("\n");
    }
    
//     getchar();
}


//NOTE: use nlNormalize; simplyfy x
void lattice::cancel_ratios_in_LLL(){

    
    if(b != NULL) {
        for(int i=1; i<=n;i++){
            if(b[i] != NULL) {
                b[i]->String();
//                 bimnlNormalize(b[i]);
            }
        }
    }
    
    if(b_star != NULL) {
        for(int i=1; i<=n;i++){
            if(b_star[i] != NULL) {
                b_star[i]->String();
            }
        }
    }
    
    if(B != NULL) {
        for(int i=1; i<=n;i++){
            if(B[i] != NULL) {
                nlNormalize(B[i],coef);
            }
        }  
    }
    
    if(d != NULL) {
        for(int i=1; i<=n;i++){
            if(d[i] != NULL) {
                nlNormalize(d[i],coef);
            }
        }  
    }
    
    if(my != NULL)
        my->String();
    
    if(lambda != NULL)
        lambda->String();
    
    if(H != NULL)
        H->String();
}

bool lattice::LLL(){
    // c = 3/4
    number three = n_Init(3, fieldcoef);
    number four  = n_Init(4, fieldcoef);
    number default_c = n_Div(three,four,fieldcoef);
    n_Delete(&three,fieldcoef);
    n_Delete(&four,fieldcoef);
    bool integral = false;
    if(nCoeff_is_Ring_Z(coef)) {
        integral = true;
    }
    return lattice::LLL(default_c,fieldcoef,integral,true,false);
}


bool lattice::LLL(number& c, coeffs c_coef, bool trans_matrix, bool integral, bool independentVectors){
    DEBUG_PRINT(("Start LLL\n"));  
    
    DEBUG_PRINT(("Delete old results\n"));
    delete_LLL_computations();
    
    if(c==NULL) {
        DEBUG_PRINT(("c==NULL => set c=3/4\n"));  
        number three = n_Init(3, coef);
        number four  = n_Init(4, coef);
        number default_c = n_Div(three,four,coef);
        n_Delete(&three,coef);
        n_Delete(&four,coef);
        this->c = default_c;
    } else {
        if(c_coef == fieldcoef) {
            this->c = n_Copy(c,c_coef);
        } else {
            nMapFunc f = n_SetMap(c_coef, fieldcoef);
            this->c = f(c, c_coef, fieldcoef);
        }
    }
     
    this->trans_matrix       = trans_matrix;
    this->integral           = integral;
    this->independentVectors = independentVectors;
      
    DEBUG_PRINT(("Create new arrays and matrices\n"));
    if(!only_gram_matrix_given){
        b = new bigintmat*[n+1]; //b[0] is not used
        b[0] = NULL;
        for(int j=1; j<=n; j++) {
            b[j] = bimCopy(basis[j]);
        }
        
        b_star = new bigintmat*[n+1]; //b_star[0] is not used
        for(int j=0; j<=n; j++) {
            b_star[j] = NULL;
        }        
    }
        
    if(trans_matrix) {
        H = new bigintmat(n,n,coef);
    }
    
    if(!integral) {
        B = new number[n+1]; //B[0] is not used
        B[0] = NULL;
        for(int j=1; j<=n; j++) {
            B[j] = n_Init(0,coef);
        }
    }
    
    if(!integral) {
        my = new bigintmat(m,n,coef);
    }
    
    if(integral) {
        d = new number[n+1];
        for(int j=0; j<=n; j++) {
            d[j] = n_Init(0,coef);
        }
    }
    
    if(integral) {
        lambda = new bigintmat(m,n,coef);
    }
        
    
    DEBUG_PRINT(("Initialize\n"));
    
    DEBUG_LLL();
    
    int k = 2;
    int k_max = 1;
    
    if(integral) {
        n_Delete(&d[0],coef);
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
    
    
    if (n == 0) { //NOTE: necessary?
        return true;
    }
    if (n == 1) {
        return false;
    }
    
    
    do{
        DEBUG_PRINT(("Incremental Gram-Schmidt\n"));
        
        DEBUG_LLL();
        
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
        
       
        bool LLL_condition = false;
        do{
            DEBUG_PRINT(("Test LLL condition\n"));
            
            DEBUG_LLL();
            
            if(integral) {
                REDI(k,k-1);
            } else {
                RED(k,k-1);
            }         
            
            DEBUG_LLL();
            
            if(integral) {
                number leftside = n_Mult(d[k], d[k-2], coef);
                number d_kminusone_squared = n_Mult(d[k-1], d[k-1], coef);
                number c_times_stuff = n_Mult(this->c,d_kminusone_squared,coef);
                number lambda_k_kminusone_squared = n_Mult(lambda->view(k,k-1), lambda->view(k,k-1), coef);
                number rightside = n_Sub(c_times_stuff, lambda_k_kminusone_squared, coef);
                LLL_condition = ! n_Greater(rightside,leftside,coef);
                n_Delete(&leftside,coef);
                n_Delete(&d_kminusone_squared,coef);
                n_Delete(&c_times_stuff,coef);
                n_Delete(&lambda_k_kminusone_squared,coef);
                n_Delete(&rightside,coef);
            } else {
                number my_squared = n_Mult(my->view(k,k-1), my->view(k,k-1), coef);
                number difference = n_Sub(this->c,my_squared, coef);
                number product = n_Mult(difference, B[k-1], coef);
                
                LLL_condition = ! n_Greater(product, B[k], coef);
                
                n_Delete(&my_squared,coef);
                n_Delete(&difference,coef);
                n_Delete(&product,coef);
            }
            
            DEBUG_VAR(LLL_condition);
            
            if(!LLL_condition) {
                if(integral) {
                    if(independentVectors) {
                        SWAPI(k,k_max);
                    } else {
                        SWAPG_integral(k,k_max);
                    }
                } else {
                    if(independentVectors) {
                        SWAP(k,k_max);
                    } else {
                        SWAPG(k,k_max);
                    }
                }
                if(k>2){
                    k--;
                }
            }
        } while(!LLL_condition);
        
        for(int l=k-2; l>0; l--){
            if(integral) {
                REDI(k,l);
            } else {
                RED(k,l);
            }
            DEBUG_LLL();
        }
        k++;
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

    number n_1    = n_Init( 1,coef);
    number n_neg1 = n_Init(-1,coef);
    number n_2    = n_Init( 2,coef);
    
    number n_1div2    = n_Div(n_1,n_2,coef);
    number n_neg1div2 = n_Div(n_neg1,n_2,coef);
    
    number my_kl = my->get(k,l);
    
    DEBUG_N(my_kl);
    
    bool abs_my_bigger_than_half = n_Greater(my_kl,n_1div2,coef) || n_Greater(n_neg1div2,my_kl,coef);

    if(abs_my_bigger_than_half) { 
                
        number q = round(my_kl,coef);
        DEBUG_N(q);
        
        if(only_gram_matrix_given){
            for(int i=1; i<=n; i++){
                number prod = n_Mult(gram_matrix_view(i,l),q,coef);
                number diff = n_Sub(gram_matrix_view(i,k),prod,coef);
                n_Delete(&prod,coef);
                gram_matrix_rawset(i,k,diff);
            }
        } else {
            bigintmat * prod = bimMult(b[l],q,coef);
            b[k]->sub(prod);
            delete prod;
        }

        if(trans_matrix) {
            number prod = n_Mult(q,n_neg1,coef);
            H->addcol(k,l,prod,coef);
            n_Delete(&prod,coef);
        }
        
        {
            number diff = n_Sub(my->view(k,l),q,coef);
            my->set(k,l,diff,coef);
            n_Delete(&diff,coef);
        }
        
        
        for(int i=1;i<=l-1;i++){
            number prod = n_Mult(q, my->view(l,i),coef);
            number diff = n_Sub(my->view(k,i), prod, coef);
            my->set(k,i,diff,coef);
            n_Delete(&prod,coef);
            n_Delete(&diff,coef);
        }
        
        n_Delete(&q,coef);
    }
    
    n_Delete(&n_1,coef);
    n_Delete(&n_neg1,coef);
    n_Delete(&n_2,coef);
    n_Delete(&n_1div2,coef);
    n_Delete(&n_neg1div2,coef);
    n_Delete(&my_kl,coef);
    
    DEBUG_PRINT(("End of RED\n"));
}

void lattice::REDI(int k, int l){
    DEBUG_PRINT(("Start REDI with k=%d and l=%d\n",k,l));
    
    number n_neg1 = n_Init(-1,coef);
    number n_2 = n_Init( 2,coef);
    
    number two_lambda = n_Mult(n_2,lambda->view(k,l),coef);
    number minus_d_l = n_Mult(n_neg1,d[l],coef);
    
    if(n_Greater(two_lambda,d[l],coef)||n_Greater(minus_d_l,two_lambda,coef)) {
        
        number q;
        {
            number dividend = n_Add(two_lambda,d[l],coef);
            number divisor = n_Mult(n_2,d[l],coef);
            DEBUG_N(dividend);
            DEBUG_N(divisor);
            
            number mod = n_IntMod(dividend,divisor,coef);
            number difference = n_Sub(dividend,mod,coef);
            
            q = n_IntDiv(difference,divisor,coef);
            
            n_Delete(&dividend,coef);
            n_Delete(&divisor,coef);
            n_Delete(&mod,coef);
            n_Delete(&difference,coef);
        }
        // number q = n_IntDiv(dividend,divisor,coef);
        // Doesn't work without error
        DEBUG_N(q);
        
        if(trans_matrix) {
            number prod = n_Mult(q,n_neg1,coef);
            H->addcol(k,l,prod,coef);
            n_Delete(&prod,coef);
        }
        
        if(only_gram_matrix_given){
            for(int i=1; i<=n; i++){
                number prod = n_Mult(gram_matrix_view(i,l),q,coef);
                gram_matrix_rawset(i,k,n_Sub(gram_matrix_view(i,k),prod,coef));
                n_Delete(&prod,coef);
            }
        } else {
            bigintmat * prod = bimMult(b[l],q,coef);
            b[k]->sub(prod);
            delete prod;
        }
        
        {
            number prod = n_Mult(q, d[l],coef);
            number diff = n_Sub(lambda->view(k,l), prod,coef);
            lambda->set(k,l,diff,coef);
            n_Delete(&prod,coef);
            n_Delete(&diff,coef);
        }
        
        for(int i=1;i<=l-1;i++){
            number prod = n_Mult(q, lambda->view(l,i),coef);
            number diff = n_Sub(lambda->view(k,i), prod, coef);
            lambda->set(k,i,diff,coef);
            n_Delete(&prod,coef);
            n_Delete(&diff,coef);
        }
    }
    
    n_Delete(&n_neg1,coef);
    n_Delete(&n_2,coef);
    
    n_Delete(&two_lambda,coef);
    n_Delete(&minus_d_l,coef);
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
        DEBUG_VAR(j);
        number my_kj = my->get(k,j);
        DEBUG_N(my_kj);
        my->set(k,j,my->view(k-1,j),coef);
        my->rawset(k-1,j,my_kj,coef);
    }
    
    number my_ = my->get(k,k-1);
    
    number B_;
    {
        number my_squared = n_Mult(my_,my_,coef);
        number prod = n_Mult(my_squared, B[k-1], coef);
        B_ = n_Add(B[k], prod, coef);
        n_Delete(&my_squared,coef);
        n_Delete(&prod,coef);
    }
    DEBUG_N(B_);
    
    
    {
        number prod = n_Mult(my_, B[k-1], coef);
        number div = n_Div(prod, B_, coef);
        my->set(k,k-1,div,coef);
        n_Delete(&prod,coef);
        n_Delete(&div,coef);
    }
    
    if(!only_gram_matrix_given) {
        bigintmat * b_ = bimCopy(b_star[k-1]);
        DEBUG_BIM(b_);
        
        {
            bigintmat * prod = bimMult(b_, my_, coef);
            delete b_star[k-1];
            b_star[k-1] = bimAdd(b_star[k], prod);
            delete prod;
        }
        
        {
            number quot = n_Div(B[k], B_, coef);
            bigintmat * prod1 = bimMult(b_, quot,coef);
            bigintmat * prod2 = bimMult( b_star[k], my->view(k,k-1),coef);
            
            delete b_star[k];
            b_star[k] = bimSub(prod1, prod2);
            
            n_Delete(&quot,coef);
            delete prod1;
            delete prod2;
        }
        
        delete b_;
    }
    
    
    {
        number prod = n_Mult(B[k], B[k-1], coef);
        n_Delete(&B[k],coef);
        B[k] = n_Div(prod, B_, coef);
        n_Delete(&prod,coef);
    }
    DEBUG_N(B[k]);
    
    n_Delete(&B[k-1],coef);
    B[k-1] = n_Copy(B_, coef);
    for(int i = k+1; i <= k_max; i++){
        DEBUG_VAR(i);
        number t = my->get(i,k);
        number prod1 = n_Mult(my_, t, coef);
        number diff = n_Sub(my->view(i,k-1), prod1, coef);
        my->set(i,k,diff,coef);
        number prod2 = n_Mult(my->view(k,k-1), my->view(i,k), coef);
        number sum = n_Add(t, prod2, coef);
        my->set(i,k-1, sum, coef);
        
        n_Delete(&t,coef);
        n_Delete(&prod1,coef);
        n_Delete(&diff,coef);
        n_Delete(&prod2,coef);
        n_Delete(&sum,coef);
    }
    
    n_Delete(&my_,coef);
    n_Delete(&B_,coef);
    
    
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
        lambda->rawset(k-1,j,lambda_kj,coef);
    }
    
    number lambda_ = lambda->get(k,k-1);
    
    
    number B;
    {
        number prod = n_Mult(d[k-2],d[k],coef);
        number lambda_squared = n_Mult(lambda_,lambda_,coef);
        number sum = n_Add(prod,lambda_squared,coef);
        B = n_Div(sum,d[k-1],coef);
        n_Delete(&prod,coef);
        n_Delete(&lambda_squared,coef);
        n_Delete(&sum,coef);
    }
    
    for(int i = k+1; i <= k_max; i++){
        number t = lambda->get(i,k);
        
        {
            number prod1 = n_Mult(d[k],lambda->view(i,k-1),coef);
            number prod2 = n_Mult(lambda_,t,coef);
            number diff = n_Sub(prod1,prod2,coef);
            number quot = n_Div(diff,d[k-1],coef);
            lambda->set(i,k,quot,coef);
            n_Delete(&prod1,coef);
            n_Delete(&prod2,coef);
            n_Delete(&diff,coef);
            n_Delete(&quot,coef);
        }
        
        {
            number prod1 = n_Mult(B,t,coef);
            number prod2 = n_Mult(lambda_,lambda->view(i,k),coef);
            number sum = n_Add(prod1,prod2,coef);
            number quot = n_Div(sum,d[k],coef);
            lambda->set(i,k-1,quot,coef);
            n_Delete(&prod1,coef);
            n_Delete(&prod2,coef);
            n_Delete(&sum,coef);
            n_Delete(&quot,coef);
        }
        
        n_Delete(&t,coef);
    }
    
    n_Delete(&d[k-1],coef);
    d[k-1] = B;

    n_Delete(&lambda_,coef);
    
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
                number prod1 = n_Mult(my->view(j,i),my->view(k,i),coef);
                number prod2 = n_Mult(prod1,B[i],coef);
                n_InpAdd(sum,prod2,coef);
                n_Delete(&prod1,coef);
                n_Delete(&prod2,coef);
            }
            number diff = n_Sub(gram_matrix_view(k,j),sum,coef);
            my->rawset(k,j,n_Mult(diff,B[j],coef),coef);
            n_Delete(&diff,coef);
            n_Delete(&sum,coef);
        }
        
        number sum = n_Init(0,coef);
        for(int i=1; i<k; i++) {
            number prod1 = n_Mult(my->view(k,i),my->view(k,i),coef);
            number prod2 = n_Mult(prod1,B[i],coef);
            n_InpAdd(sum,prod2,coef);
            n_Delete(&prod1,coef);
            n_Delete(&prod2,coef);
        }
        B[k] = n_Sub(gram_matrix_view(k,k),sum,coef);
        n_Delete(&sum,coef);
    } else {
        delete b_star[k];
        b_star[k] = bimCopy(b[k]);  
        
        for(int j=1; j<k; j++) {
            DEBUG_VAR(j);
            {
                number scalprd = scalarproduct(b[k],b_star[j]);
                my->rawset(k,j,n_Div(scalprd,B[j],coef),coef);
                n_Delete(&scalprd,coef);
            }
            
            {
                bigintmat * prod = bimMult(b_star[j],my->view(k,j),coef);
                b_star[k]->sub(prod);
                delete prod;
            }
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
            number prod1 = n_Mult(d[i],u,coef);
            number prod2 = n_Mult(lambda->view(k,i),lambda->view(j,i),coef);
            number diff = n_Sub(prod1,prod2,coef);
            
            n_Delete(&u,coef);
            u = n_Div(diff,d[i-1],coef);
            
            n_Delete(&prod1,coef);
            n_Delete(&prod2,coef);
            n_Delete(&diff,coef);
        }
        
        if(j<k) {
            lambda->rawset(k,j,u,coef);
        } else {
            d[k] = u;
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
                my->rawset(k,j,n_Init(0,coef));
            } else {
                number sum = n_Init(0,coef);
                
                for(int i=1; i<k; i++) {
                    number prod1 = n_Mult(my->view(j,i),my->view(k,i),coef);
                    number prod2 = n_Mult(prod1,B[i],coef);
                    n_InpAdd(sum,prod2,coef);
                    n_Delete(&prod1,coef);
                    n_Delete(&prod2,coef);
                }
                
                {
                    number diff = n_Sub(gram_matrix_view(k,j),sum,coef);
                    my->rawset(k,j,n_Mult(diff,B[j],coef),coef);
                    n_Delete(&diff,coef);
                }
                
                n_Delete(&sum,coef);
            }
        }
        
        {
            number sum = n_Init(0,coef);
            for(int i=1; i<k; i++) {
                number prod1 = n_Mult(my->view(k,i),my->view(k,i),coef);
                number prod2 = n_Mult(prod1,B[i],coef);
                n_InpAdd(sum,prod2,coef);
                n_Delete(&prod1,coef);
                n_Delete(&prod2,coef);                
            }
            B[k] = n_Sub(gram_matrix_view(k,k),sum,coef);
            n_Delete(&sum,coef);
        }
    } else {
        for(int j=1; j<k; j++) {
            if(n_IsZero(B[j],coef)) {
                my->rawset(k,j,n_Init(0,coef));
            } else {
                number scalprd = scalarproduct(b[k],b_star[j]);
                my->rawset(k,j,n_Div(scalprd,B[j],coef),coef);
                n_Delete(&scalprd,coef);
            }
        }
        
        delete b_star[k];
        b_star[k] = bimCopy(b[k]);
        for(int j=1; j<k; j++) {
            bigintmat * prod = bimMult(b_star[j],my->view(k,j),coef);
            b_star[k]->sub(prod);
            delete prod;
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
    DEBUG_N(a);
    DEBUG_PRINT(("check input\n"));
    if(!n_GreaterZero(a,coef)){
        if(n_IsZero(a,coef) && n==m){
            return new bigintmat(m,1,coef);
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
    number minusOne = n_Init(-1,fieldcoef);
    
    //backtracking for elements
    //DEBUG_BLOCK(true);
    DEBUG_PRINT(("Start backtracking\n"));
    DEBUG_PRINT(("Initialize vector and other variables\n"));
    std::vector<std::pair<number,bigintmat*> > elementsvector;
    elementsvector.push_back( std::make_pair(n_Init(0,fieldcoef), new bigintmat(m,1,coef)));
    if( x != NULL){
        delete x;
        x=NULL;
    }
    x = new bigintmat(n,1,fieldcoef);
    //x->Print();PrintS("\n");
    DEBUG_PRINT(("map a\n"));
    if(bound != NULL){
        for(int i=1;i<=n;i++){
            n_Delete(&bound[i],fieldcoef);
        }
        delete[] bound;
        bound = NULL;
    }
    bound = new number[n+1];
    nMapFunc f = n_SetMap(coef, fieldcoef);
    bound[1] = f(a, coef, fieldcoef);////////////
    //map a to fieldcoef
    DEBUG_PRINT(("set bound\n"));
    for(int i = 2; i<n+1; i++){
        bound[i] = n_Copy(bound[1],fieldcoef);
        //n_Print(bound[i],fieldcoef);PrintS("\n");
    }
    DEBUG_PRINT(("find element\n"));
    //bigintmat* elements = enumerate_next(a);
    increase_x(1);
    number check = enumerate_get_next();
    while(!n_Equal(minusOne,check,fieldcoef)){
        //append x to elements
        DEBUG_PRINT(("new element to list\n"));
        //elements->appendCol(bimChangecoeff(x,coef));
        check = n_Sub(bound[1],check,fieldcoef);
        check = n_Sub(bound[n],check,fieldcoef);
        elementsvector.push_back(std::make_pair(n_Copy(check,fieldcoef),bimCopy(x)));
        //n_Print(elementsvector[elementsvector.size()-1].first,fieldcoef);PrintS("\n");
        for(unsigned i=1; i<elementsvector.size();i++){
            if(n_Greater(elementsvector[i].first,check,fieldcoef)){
                elementsvector.pop_back();
                elementsvector.insert(elementsvector.begin()+i,std::make_pair(n_Copy(check,fieldcoef),bimCopy(x)));
                DEBUG_VAR(elementsvector.size());
                break;
            }
        }
        if(elementsvector.size() >= 1000000){
            elementsvector.pop_back();
        }
        increase_x(1);
        n_Delete(&check,fieldcoef);
        check = enumerate_get_next();
        DEBUG_PRINT(("got it\n"));
    }
    DEBUG_PRINT(("generate bigintmat for return\n"));
    bigintmat* elements = new bigintmat(m,1,coef);
    
    
    for(unsigned i=1; i<elementsvector.size();i++){
        //elements->appendCol(bimChangeCoeff(elementsvector[i].second,coef));
        bigintmat * temp = bimChangeCoeff(elementsvector[i].second, coef);
        elements->appendCol(temp);
        delete temp;
    }
    for(int i=1;i<=n;i++){
        n_Delete(&bound[i],fieldcoef);
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
    if (in == NULL || in->rows() * in->cols()!=n){
        DEBUG_BIM(in);
        DEBUG_PRINT(("Dimension error of input\n"));
        return NULL;
    }
    
    if(!n_GreaterZero(a,coef)){
        DEBUG_PRINT(("negative input\n"));
        return NULL;
    }
    
    if( bound != NULL){
        for(int i=1;i<=n;i++){
            n_Delete(&bound[i],fieldcoef);
        }
        delete[] bound;
        bound=NULL;
    }
    
    DEBUG_PRINT(("check quadratic\n"));
    
    if( Q == NULL){
        if(quadratic_supplement()){
            return NULL;
        }
    }
    DEBUG_PRINT(("Q defined\n"));
    
    //usefull numbers
    number minusOne = n_Init(-1,fieldcoef);
    
    //backtracking for elements
    //DEBUG_BLOCK(true);
    DEBUG_PRINT(("Start backtracking\n"));
    DEBUG_PRINT(("Initialize variables\n"));
    nMapFunc f = n_SetMap(coef, fieldcoef);
    delete x;
    x = new bigintmat(n,1,fieldcoef);
    for(int i=0;i<n;i++){
        x->rawset(i,f(in->view(i),coef,fieldcoef),fieldcoef);
    }
    x = bimChangeCoeff(in,fieldcoef);
    bound = new number[n+1];
    DEBUG_PRINT(("set bound\n"));
    bound[n] = f(a, coef, fieldcoef);//map a to fieldcoef
    for(int j = n; j>1; j--){
        number check = check_bound(j);
        bound[j-1] = n_Sub(bound[j],check,fieldcoef);
        n_Delete(&check, fieldcoef);
    }
    DEBUG_PRINT(("find element\n"));
    number norm = enumerate_get_next();
    DEBUG_PRINT(("generate bigintmat for return\n"));
    bigintmat * out;
    if(n_Equal(minusOne,norm,fieldcoef)){
        out = NULL;
    } else {
        out = bimChangeCoeff(x, coef);
    }
    n_Delete(&minusOne, fieldcoef);
    n_Delete(&norm, fieldcoef);
    return out;
}

bigintmat * lattice::enumerate_next(number a){
    //DEBUG_BLOCK(true);
    DEBUG_PRINT(("Start enumerate_next number\n"));
    bigintmat * in =new bigintmat(m,1,coef);
    if(x == NULL){
        in->rawset(1,1,n_Init(1,coef),coef);
    } else {
        in = bimChangeCoeff(x,coef);
    }
    bigintmat * out = enumerate_next(a,in);
    delete in;
    return out;
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
    nMapFunc f = n_SetMap(fieldcoef, coef);
    number a = f(bound[n],fieldcoef,coef);
    DEBUG_PRINT(("enumerate_next bigintmat\n"));
    bigintmat * out = enumerate_next(a,in);
    n_Delete(&a,coef);
    return out;
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
    number minusOne = n_Init(-1,fieldcoef);
    DEBUG_PRINT(("find element\n"));
    number norm = enumerate_get_next();
    DEBUG_PRINT(("generate bigintmat for return\n"));
    if(n_Equal(minusOne,norm,fieldcoef)){
        DEBUG_PRINT(("no element found\n"));
        number one = n_Init(1,fieldcoef);
        if(n_GreaterZero(x->view(1,1),fieldcoef)){
            x->rawset(1,1,n_Add(x->view(1,1),minusOne,fieldcoef),fieldcoef);
        } else {
            if(!n_Equal(minusOne, x->view(1,1),fieldcoef)){
                x->rawset(1,1,n_Add(x->view(1,1),one,fieldcoef),fieldcoef);
            }
        }
        n_Delete(&minusOne,fieldcoef);
        n_Delete(&one,fieldcoef);
        n_Delete(&norm,fieldcoef);
        return NULL;
    }
    n_Delete(&minusOne,fieldcoef);
    n_Delete(&norm,fieldcoef);
    bigintmat * out = bimChangeCoeff(x, coef);
    return out;
}

//Private
number lattice::enumerate_get_next(){
    DEBUG_BLOCK(false);
    DEBUG_PRINT(("enumerate_get_next\n"));
    int index =1;
    //x->Print();PrintS("\n");
    //DEBUG_PRINT(("first time changing x\n"));
    //increase_x(1);
    DEBUG_PRINT(("actual backtracking\n"));
    number check;
    while (index <= m) {
        DEBUG_PRINT(("update check\n"));
        check = check_bound(index);
        //DEBUG_PRINT(("check check\n"));
        if (n_Greater(check,bound[index],fieldcoef)){
            //DEBUG_PRINT(("element to great\n"));
            if(!(n_GreaterZero(x->view(index,1),fieldcoef) || n_IsZero(x->view(index,1),fieldcoef))){
                n_Delete(&bound[index],fieldcoef);
                bound[index] = n_Init(0,fieldcoef);
                x->rawset(index,1,n_Init(0,fieldcoef),fieldcoef);
                index++;
                if(index<= m){
                    increase_x(index);
                }
            } else {
                if(index == n){
                    return n_Init(-1,fieldcoef);
                }
                x->rawset(index,1,n_Init(-1,fieldcoef),fieldcoef);
            }
        } else if(index == 1){
            DEBUG_PRINT(("possible new element\n"));
            if(n_IsZero(x->view(n,1),fieldcoef)){
                int j=n-1;
                while(j>=1 && n_IsZero(x->view(j,1),fieldcoef)){
                    j--;
                }//DEBUG_VAR(j);
                if(j==0){
                    return check;
                }
                if(n_GreaterZero(x->view(j,1),fieldcoef)){
                    return check;
                } else {
                    index = j+1;
                    for( j=1;j<index;j++){
                        x->rawset(j,1,n_Init(0,fieldcoef),fieldcoef);
                    }
                    x->rawset(index,1,n_Init(1,fieldcoef),fieldcoef);
                }
            } else {
                return check;
            }
        } else {
            //DEBUG_PRINT(("reduce index\n"));
            index--;
            n_Delete(&bound[index],fieldcoef);
            bound[index] = n_Sub(bound[index+1],check,fieldcoef);
        }
        n_Delete(&check,fieldcoef);
    }
    return n_Init(-1,fieldcoef);
}

bool lattice::quadratic_supplement(){
    //DEBUG_BLOCK(true);
    delete Q;
    Q = NULL;
    if(n != m) {  //NOTE: rank?
        return true;
    }
    Q = get_gram_matrix();
    
    number zero = n_Init(0,fieldcoef);
    
    DEBUG_PRINT(("Begin Quadratic Suplement\n"));
    for(int i = 1; i<Q->cols();i++){
        if(n_IsZero( Q->view(i,i), fieldcoef)){
            DEBUG_PRINT(("matrix not positive definite\n"));
            delete Q;
            Q = NULL;
            n_Delete(&zero,fieldcoef);
            return true;
        }
        for( int j=i+1; j<=Q->cols();j++){
            Q->set(j,i,Q->view(i,j),fieldcoef);
            Q->rawset(i,j,n_Div(Q->view(i,j),Q->view(i,i),fieldcoef),fieldcoef);
        }
        for(int m=i+1; m<=Q->rows();m++){
            for(int n=i+1; n<=Q->cols();n++){
                number mult = n_Mult(Q->view(m,i),Q->view(i,n),fieldcoef);
                Q->rawset(m,n,n_Sub(Q->view(m,n),mult,fieldcoef),fieldcoef);
                n_Delete(&mult,fieldcoef);
            }
        }
    }
    
    DEBUG_PRINT(("Set Zeros\n"));
    for(int i = 2; i<=Q->cols();i++){
        for(int j = 1; j<i;j++){
            Q->set(i,j,zero,fieldcoef);
        }
    }
    n_Delete(&zero,fieldcoef);
    DEBUG_PRINT(("Test: matrix positive definite\n"));
    for(int i=1; i<=Q->cols();i++){
        if(!n_GreaterZero( Q->view(i,i), fieldcoef)){
            DEBUG_PRINT(("matrix not positive definite\n"));
            delete Q;
            Q = NULL;
            return true;
        }
    }
    return false;
}

void lattice::increase_x(int index){
    number one = n_Init(1,fieldcoef);
    if (n_GreaterZero(x->view(index,1),fieldcoef) || n_IsZero(x->view(index,1),fieldcoef)){
        x->rawset(index,1, n_Add(one,x->view(index,1),fieldcoef),fieldcoef); //x_i=x_i+1
    } else {
        x->rawset(index,1, n_Sub(x->view(index,1),one,fieldcoef),fieldcoef);//x_i=x_i-1
    }
    n_Delete(&one,fieldcoef);
}

number lattice::check_bound(int index){
    //DEBUG_BLOCK(true);DEBUG_PRINT(("check bound\n"));DEBUG_VAR(index);
    number check = n_Init(0,fieldcoef);
    for(int i=index + 1;i<=Q->cols();i++){
        number mult = n_Mult(x->view(i,1),Q->view(index,i),fieldcoef);
        n_InpAdd(check,mult,fieldcoef);
        n_Delete(&mult,fieldcoef);
        delete mult;
    }
    n_InpAdd(check, x->view(index,1), fieldcoef);
    n_InpMult(check, check, fieldcoef);
    n_InpMult(check, Q->view(index,index), fieldcoef);
    return check;
}


///////////////////////////////////////
//               Getter             ///
///////////////////////////////////////

bigintmat * lattice::get_basis() {
    if(only_gram_matrix_given) {
        Werror("lattice is only defined by its gram matrix");
        return NULL;
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
                r->rawset(i,j,scalarproduct(basis[i],basis[j]));
            }
        }
    }
    return r;
}

bigintmat * lattice::get_lattice_element(bigintmat * x){
    bigintmat * out = new bigintmat(m,1,coef);
    if(x==NULL){
        return NULL;
    }
    if(x->rows() * x->cols()!=m){
        Werror("dimension error of input\n");
        return NULL;
    }
    nMapFunc f = n_SetMap(x->basecoeffs(),coef);
    for(int i=1;i<=n;i++){
        number a = f(x->view(i-1),x->basecoeffs(),coef);
        bigintmat * temp = bimMult(basis[i],a,coef);
        n_Delete(&a,coef);
        out->add(temp);
        delete temp;
    }
    return out;
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


///////////////////////////////////////
//        None-Member-Funtions      ///
///////////////////////////////////////

number scalarproduct(bigintmat * a, bigintmat * b) {
    DEBUG_BLOCK(false);
    DEBUG_PRINT(("Start scalarproduct\n"));  
    
    if(a->cols()!=1) {
        Werror("a->cols()!=1 in scalarproduct(a,b)\n");
        return NULL;
    }
    if(b->cols()!=1) {
        Werror("b->cols()!=1 in scalarproduct(a,b)\n");
        return NULL;
    }
    if(a->rows()!=b->rows()) {
        Werror("a->rows()!=b->rows() in scalarproduct(a,b)\n");
        return NULL;
    }
    if(a->basecoeffs()!=b->basecoeffs()) {
        Werror("a->basecoeffs()!=b->basecoeffs() in scalarproduct(a,b)\n");
        return NULL;
    }

    coeffs coef = a->basecoeffs();
    number p = n_Init(0,coef);
    for(int i = 1; i <= b->rows(); i++){
        DEBUG_VAR(i);
        number prod = n_Mult(a->view(i,1), b->view(i,1), coef);
        DEBUG_N(prod);
        n_InpAdd(p, prod, coef);
        n_Delete(&prod,coef);
    }
    DEBUG_N(p);
    DEBUG_PRINT(("End scalarproduct\n"));  
    return p;
}

//based on /Singular/LIB/atkins.lib
//NOTE: rounds .5 up 
number round(number r,coeffs coef) {
    DEBUG_PRINT(("Start round\n"));
    DEBUG_N(r);
    
    if((r == NULL) || (coef == NULL)) {
        Werror("input of round is NULL");
        DEBUG_PRINT(("End round\n"));
        return NULL;
    }
    if(getCoeffType(coef)==n_Z){
        DEBUG_PRINT(("getCoeffType(coef)==n_Z\n"));
        DEBUG_PRINT(("End round\n"));;  
        return n_Copy(r,coef);
    }
    if(getCoeffType(coef)==n_Q){
        DEBUG_PRINT(("getCoeffType(coef)==n_Q\n"));
        number n_0    = n_Init( 0,coef);
        number n_1    = n_Init( 1,coef);
        number n_neg1 = n_Init(-1,coef);
        number n_2    = n_Init( 2,coef);
        
        number n_1div2    = n_Div(n_1,n_2,coef);
        number n_neg1div2 = n_Div(n_neg1,n_2,coef);
        
        number r_;
        
        if(n_Greater(r,n_0,coef)) {
            r_ = n_Add(r, n_1div2, coef);
        } else {
            r_ = n_Add(r, n_neg1div2, coef);
        }

        
        number numerator = n_GetNumerator(r_,coef);
        DEBUG_N(numerator);
        number denominator = n_GetDenom(r_,coef);
        DEBUG_N(denominator);
        
        number q = n_IntDiv(numerator,denominator,coef);
        DEBUG_N(q);
        
        n_Delete(&n_0,coef);
        n_Delete(&n_1,coef);
        n_Delete(&n_neg1,coef);
        n_Delete(&n_2,coef);
        n_Delete(&n_1div2,coef);
        n_Delete(&n_neg1div2,coef);
        n_Delete(&r_,coef);
        n_Delete(&numerator,coef);
        n_Delete(&denominator,coef);        
        
        DEBUG_PRINT(("End round\n"));
        return q;
    }
    if((getCoeffType(coef)==n_long_R) || (getCoeffType(coef)==n_R)){
        DEBUG_PRINT(("(getCoeffType(coef)==n_long_R) || (getCoeffType(coef)==n_R)\n"));  
        number n_0    = n_Init( 0,coef);
        number n_1    = n_Init( 1,coef);
        number n_neg1 = n_Init(-1,coef);
        number n_2    = n_Init( 2,coef);
        
        number n_1div2    = n_Div(n_1,n_2,coef);
        
        number a;
        if(n_Greater(r,n_0,coef)) {
            a = n_Copy(r,coef);
        } else {
            a = n_Mult(n_neg1,r,coef);
        }
        DEBUG_N(a);
        
        number v = n_Div(r,a,coef);
        DEBUG_N(v);
        
        number d = n_Init(10,coef);
        DEBUG_N(d);
        
        int e = 0;
        while(1){
            e=e+1;
            DEBUG_VAR(e);
            
            number pow;
            n_Power(d,e,&pow,coef);
            number diff = n_Sub(a,pow,coef);
            n_Delete(&pow,coef);
            if(n_Greater(n_0,diff,coef)){
                n_Delete(&diff,coef);
                e=e-1;
                break;
            }
            n_Delete(&diff,coef);
        }

        number b = n_Copy(a,coef);
        DEBUG_N(b);
        
        for(int k=0;k<=e;k++){
            DEBUG_VAR(k);
            
            number s;
            n_Power(d,e-k,&s,coef);
            s = n_InpNeg(s,coef);
            while(1)
            {
                n_InpAdd(b,s,coef);
                if(n_Greater(n_0,b,coef)){
                    s = n_InpNeg(s,coef);
                    n_InpAdd(b,s,coef);
                    break;
                }
            }
            DEBUG_N(b);
            n_Delete(&s,coef);
        }
        number diff = n_Sub(a,b,coef);
        number diffmin1 = n_Sub(diff,n_neg1,coef);
        
        number result;
        
        if(n_Greater(n_1div2,b,coef)){
            result = n_Mult(v,diff,coef);
        }else{
            result = n_Mult(v,diffmin1,coef);
        }
        
        n_Delete(&n_0,coef);
        n_Delete(&n_1,coef);
        n_Delete(&n_neg1,coef);
        n_Delete(&n_2,coef);
        n_Delete(&n_1div2,coef);
        n_Delete(&a,coef);
        n_Delete(&v,coef);
        n_Delete(&d,coef);
        n_Delete(&b,coef);
        
        DEBUG_N(result);
        DEBUG_PRINT(("End round\n"));
        return result;
    }
    
    if(getCoeffType(coef)==n_long_C){
        DEBUG_PRINT(("getCoeffType(coef)==n_long_C\n"));
        
        number r_real = ngcRePart(r,coef);    
        number r_imag = ngcImPart(r,coef);
        
        
        if(n_IsZero(r_imag,coef)) {
            r = r_real;
        } else {
            number tmp = n_Sub(r,r_real,coef);
            number i_unit = n_Div(tmp,r_imag,coef);
            n_Delete(&tmp,coef);
            
            if(n_IsZero(r_real,coef)) {
                number rounded_imag = round(r_imag,coef);
                n_InpMult(rounded_imag,i_unit,coef);
                
                n_Delete(&r_real,coef);
                n_Delete(&r_imag,coef);
                n_Delete(&i_unit,coef);
                
                return rounded_imag;
            } else {
                number rounded_real = round(r_real,coef);
                number rounded_imag = round(r_imag,coef);
                n_InpMult(rounded_imag,i_unit,coef);
                number result = n_Add(rounded_real,rounded_imag,coef);
                
                n_Delete(&rounded_real,coef);
                n_Delete(&rounded_imag,coef);
                n_Delete(&r_real,coef);
                n_Delete(&r_imag,coef);
                n_Delete(&i_unit,coef);
                
                return result;
            }
        }

            
        
        number n_0    = n_Init( 0,coef);
        number n_1    = n_Init( 1,coef);
        number n_neg1 = n_Init(-1,coef);
        number n_2    = n_Init( 2,coef);
        
        number n_1div2    = n_Div(n_1,n_2,coef);
        
        number a;
        if(n_Greater(r,n_0,coef)) {
            a = n_Copy(r,coef);
        } else {
            a = n_Mult(n_neg1,r,coef);
        }
        DEBUG_N(a);
        
        number v = n_Div(r,a,coef);
        DEBUG_N(v);
        
        number d = n_Init(10,coef);
        DEBUG_N(d);
        
        int e = 0;
        while(1){
            e=e+1;
            DEBUG_VAR(e);
            
            number pow;
            n_Power(d,e,&pow,coef);
            number diff = n_Sub(a,pow,coef);
            n_Delete(&pow,coef);
            if(n_Greater(n_0,diff,coef)){
                n_Delete(&diff,coef);
                e=e-1;
                break;
            }
            n_Delete(&diff,coef);
        }

        number b = n_Copy(a,coef);
        DEBUG_N(b);
        
        for(int k=0;k<=e;k++){
            DEBUG_VAR(k);
            
            number s;
            n_Power(d,e-k,&s,coef);
            s = n_InpNeg(s,coef);
            while(1)
            {
                n_InpAdd(b,s,coef);
                if(n_Greater(n_0,b,coef)){
                    s = n_InpNeg(s,coef);
                    n_InpAdd(b,s,coef);
                    break;
                }
            }
            DEBUG_N(b);
            n_Delete(&s,coef);
        }
        number diff = n_Sub(a,b,coef);
        number diffmin1 = n_Sub(diff,n_neg1,coef);
        
        number result;
        
        if(n_Greater(n_1div2,b,coef)){
            result = n_Mult(v,diff,coef);
        }else{
            result = n_Mult(v,diffmin1,coef);
        }
        
        n_Delete(&n_0,coef);
        n_Delete(&n_1,coef);
        n_Delete(&n_neg1,coef);
        n_Delete(&n_2,coef);
        n_Delete(&n_1div2,coef);
        n_Delete(&a,coef);
        n_Delete(&v,coef);
        n_Delete(&d,coef);
        n_Delete(&b,coef);
        
        DEBUG_N(result);
        DEBUG_PRINT(("End round\n"));
        return result;
    }
    
    Werror("round is not defined on this ring");
    return NULL;
}

void bimnlNormalize(bigintmat * m){
    for(int i=1; i<=m->rows(); i++){
        for(int j=1; j<=m->cols(); j++){
            number a = m->view(i,j);
            coeffs coef = m->basecoeffs();
            nlNormalize(a,coef);
        }
    }
}


///////////////////////////////////////
//         Minkowski map            ///
///////////////////////////////////////

int minkowski(bigintmat * basiselements, number * poly,int deg, coeffs coef, int precision, lattice * latticeNF){
//    DEBUG_BLOCK(true);
    DEBUG_PRINT(("Begin Minkowski map\n"));
    DEBUG_PRINT(("Input check\n"));
//     DEBUG_PRINT(("clear latticeNF\n"));
//     delete latticeNF;
//     latticeNF = NULL;
    
    
    if(basiselements == NULL || poly == NULL || !nCoeffs_are_equal(coef,basiselements->basecoeffs())){
        WerrorS("wrong input!\n");
        return -1;
    }
    
    //char = 0
    if ( !(nCoeff_is_Ring_Z(coef) || nCoeff_is_R(coef) || nCoeff_is_Q(coef) ||
            nCoeff_is_long_R(coef) || nCoeff_is_long_C(coef)) ){
        WerrorS("Ground field not implemented!\n");
        return -1;
    }
    
    if(deg<2){
        WerrorS("degree of polynomial to small\n");
        return -1;
    }
    //check and define precision for Q
    if(precision<6){
        precision = 6;
    }
    if(precision>32767){
        precision = 32767;
    }
    if ( !(nCoeff_is_R(coef) || nCoeff_is_long_R(coef) || nCoeff_is_long_C(coef)) ){
        setGMPFloatDigits( precision+6,precision+6);
    }
    
    DEBUG_PRINT(("find roots\n"));
    ring CurrOldRing = rCopy(currRing);//need to change currRing, because rootContainer uses the coef of it
    char* n[] = {(char*)"x"};///
    ring newring = rDefault(coef, 1, n);
    rChangeCurrRing(newring);
    DEBUG_PRINT(("initialize rootContainer\n"));
    number * rootpoly = new number[deg+1];//doesn't need to be deleted since rootContaine delete it...
    for(int i=0;i<=deg;i++){
        DEBUG_VAR(i);
        rootpoly[i] = n_Copy(poly[i],coef);
    }
    rootContainer * rootcont= new rootContainer();
    DEBUG_PRINT(("rootcont->fillContainer\n"));
    rootcont->fillContainer(rootpoly, NULL, 1, deg, rootContainer::onepoly, 1 );
    DEBUG_PRINT(("irootcont->solver(1)\n"));
    rootcont->solver(1);
    int number_roots = rootcont ->getAnzRoots();
    if(number_roots != deg){
        WerrorS("something went wrong: \n\tnot all roots found\n");
        return -1;
    }
    LongComplexInfo paramComp = {si_min (precision+6, 32767),si_min (precision+8, 32767),(const char*)"i"};
    
    coeffs comp = nInitChar(n_long_C, &paramComp);
    
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
    DEBUG_PRINT(("delete some variables\n"));
    delete rootcont;
    rChangeCurrRing(CurrOldRing);
    rDelete(newring);
    DEBUG_VAR(r1);
    DEBUG_VAR(r2);
    for(int j=0;j<r2;j++){
        roots[r1+j]= n_Copy(complexroots[j],comp);
        n_Delete(&complexroots[j],comp);
    }
    delete[] complexroots;
    DEBUG_PRINT(("map elementarray to complex\n"));
    bigintmat * elements = bimChangeCoeff(basiselements,comp);
    DEBUG_PRINT(("generate output matrix\n"));
    DEBUG_PRINT(("real part\n"));
    bigintmat * complexmat = new bigintmat(r1+2*r2, elements->cols(), comp);
    for(int i=1; i<= r1; i++){
        for(int l=deg; l>0; l--){
            for(int j=1; j<=elements->cols();j++){
                DEBUG_VAR(i);
                DEBUG_VAR(j);
                DEBUG_VAR(l);
                number val = complexmat->get(i,j);
//                 number val = complexmat->get(i,j+1); //NOTE ???
                 DEBUG_PRINT(("asdas\n"));
                n_InpMult(val,roots[i-1],comp);
                 DEBUG_PRINT(("tzu\n"));
                complexmat->rawset(i,j,n_Add(val,elements->view(j,l),comp),comp);
//                 complexmat->rawset(i,j+1,n_Add(val,elements->view(j,l),comp),comp); //NOTE ???
                n_Delete(&val,comp);
            }
        }
    }
    DEBUG_PRINT(("imaginary part\n"));
    if(r2>0){
        number two = n_Init(2,comp);
        number sqrt2 = squareroot(two,comp,precision+10);
        n_Delete(&two,comp);
        for(int i=1; i<= r2; i++){
            for(int l=deg; l>0; l--){
                for(int j=1; j<=elements->cols();j++){
                    DEBUG_VAR(i);
                    DEBUG_VAR(j);
                    DEBUG_VAR(l);
                    number val = complexmat->get(r1+2*i,j);
//                     number val = complexmat->get(r1+2*i,j+1); //NOTE
                    n_InpMult(val,roots[i-1],comp);
                    complexmat->rawset(r1+2*i,j,n_Add(val,elements->view(j,l),comp),comp);
//                     complexmat->rawset(r1+2*i,j+1,n_Add(val,elements->view(j,l),comp),comp); //NOTE ???
                    n_Delete(&val,comp);
                }
            }
            for(int j=1;j<=elements->cols();j++){
                complexmat->rawset(r1+2*i,j,n_Mult(complexmat->view(r1+2*i,j),sqrt2,comp),comp);
                complexmat->rawset(r1+2*i-1,j,ngcRePart(complexmat->view(r1+2*i,j),comp),comp);
                complexmat->rawset(r1+2*i,j,ngcImPart(complexmat->view(r1+2*i,j),comp),comp);
            }
        }
        n_Delete(&two,comp);
        n_Delete(&sqrt2,comp);
    }
    DEBUG_PRINT(("delete Variables\n"));
    delete elements;
    for(int i=0;i<r1+r2;i++){
        n_Delete(&roots[i],comp);
    }
    delete[] roots;
    
    DEBUG_PRINT(("to real\n"));
    LongComplexInfo paramReal = {si_min (precision, 32767), si_min (precision, 32767), (const char*)"i"};
    coeffs real = nInitChar(n_long_R, &paramReal);
    bigintmat * realmat = bimChangeCoeff(complexmat,real);
    delete complexmat;
    nKillChar(comp);
    
    lattice * l= new lattice(realmat);
    *latticeNF = *l;
    
    delete realmat;
    return r1;
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
    n_Delete(&two,coef);
    return xn;
}


///////////////////////////////////////
//       Get nice Polynomial        ///
///////////////////////////////////////

poly get_nice_poly(poly polynom){
    DEBUG_BLOCK(true);
    DEBUG_PRINT(("Start get_nice_poly\n"));
    if ( !(nCoeff_is_Ring_Z(currRing->cf) || nCoeff_is_Q(currRing->cf) ) ){
        WerrorS("Ground ring not implemented!\n");
        return NULL;
    }

    //primes<1000
    
    //Use EquationOrder from nforder.lib
    //works because of iiLibCmd(omStrDup("nforder.lib"), TRUE,TRUE,TRUE) in singular.cc
    
    idhdl EquationOrder=ggetid("EquationOrder");
    
    if(EquationOrder == NULL){
        WerrorS("EquationOrder not found\n");
        return NULL;
    }
    
    poly polynomCopy = p_Copy(polynom,currRing);
    
    leftv arg = new sleftv();
    arg->rtyp = POLY_CMD;
    arg->data = (void*) polynomCopy;
    
    if(iiMake_proc(EquationOrder,NULL,arg)){
        WerrorS("Error: Can't use EquationOrder\n");
        return NULL;
    }
   
    coeffs R = (coeffs) iiRETURNEXPR.Data();
    nforder *  maxord = (nforder*) R->data;
    
    iiRETURNEXPR.CleanUp();
    
    DEBUG_PRINT(("\n"));
    maxord->Print();
    DEBUG_PRINT(("\n"));
    
    
    
    
    coeffs coef ;
     if(nCoeff_is_Ring_Z(coef)) {
         DEBUG_PRINT(("nCoeffis_Ring_Z\n"));
         coef = nInitChar(n_Q,NULL);  
     } else {
         DEBUG_PRINT(("!nCoeffis_Ring_Z\n"));
         coef = currRing->cf;
     }
    
    {
        number disc = maxord->getDisc();
        DEBUG_N(disc);
        n_Delete(&disc,coef);
    }
    
    int deg = (int) p_Totaldegree(polynom,currRing);
    DEBUG_VAR(deg);
    
    number * poly_in = poly2numbers(polynom,currRing,coef);
    
    
    int primes_1000[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997};
    static int size_primes_1000 = 168;
    
    for(int i=0;i<size_primes_1000;i++){
//         number d = maxord->getDisc();
//         if()
        DEBUG_PRINT(("Maximize for primes_1000[%d]=%d\n",i,primes_1000[i]));
        nforder * temp = maxord;
        number p = n_Init(primes_1000[i],coef);
        maxord = pmaximal(temp, p);
        n_Delete(&p,coef);
//        nforder_delete(temp);
    }
    //we assume that maxord is maximal
    
    {
        number disc = maxord->getDisc();
        DEBUG_N(disc);
        n_Delete(&disc,coef);
    }
    
   
   
    bigintmat * basis = bimChangeCoeff(maxord->viewBasis(),coef);
//     nMapFunc f = n_SetMap(currRing->cf,coef);
//     number div = maxord->getDiv();
//     number temp2 = f(div,currRing->cf,coef);
//     basis->skaldiv(temp2);
//     n_Delete(&div,currRing->cf);
//     n_Delete(&temp2,coef);
    

    DEBUG_PRINT(("maxord->getDiv()\n"));
    number div = maxord->getDiv();
//     /D*EBUG_PRINT(("basis->skaldiv(div)\n"));
//     */basis->skaldiv(div);
    
    int precision = 42;//magic number
    lattice * latticeNF =  (lattice *)omAlloc(sizeof(lattice));
    
    DEBUG_BIM(basis);
    for(int i=0; i<=deg; i++) {
        DEBUG_N(poly_in[i]);
    }
    DEBUG_VAR(precision);
    
    
    
    
    DEBUG_VAR(precision);
    int r1 = minkowski(basis,poly_in,deg,coef,precision,latticeNF);
    number c = NULL;
    while(latticeNF->LLL(c,NULL,false,false,true) && precision < 32767){
//         delete latticeNF;
        precision = precision +5;
        DEBUG_VAR(precision);
        r1 = minkowski(basis,poly_in,deg,coef,precision,latticeNF);
    }
    
    
    bigintmat * LLLbasis = latticeNF->get_reduced_basis();
    DEBUG_BIM(LLLbasis);
    
    coeffs real = LLLbasis->basecoeffs();
    
    number bound = t2norm(polynom,currRing,real,precision);
    DEBUG_PRINT(("bound: "));
    DEBUG_CMD(n_Print(bound,real));
    DEBUG_PRINT(("\n"));
    
    poly reduced_poly = (poly) omAlloc(sizeof(poly));
    ///check for primitive element in LLLbasis
    for(int j=1;j<=deg;j++){
        DEBUG_VAR(j);
        bigintmat * testelement = new bigintmat(LLLbasis->rows(),1,real); 
        LLLbasis->getcol(j,testelement);
        
        DEBUG_BIM(testelement);
        
        number testnorm = t2norm(testelement);
        DEBUG_PRINT(("testnorm: "));
        DEBUG_CMD(n_Print(testnorm,real));
        DEBUG_PRINT(("\n"));
        
        if(n_Greater(testnorm,bound,real)){
            DEBUG_PRINT(("n_Greater(testnorm,bound,real)\n"));
            delete testelement;
            break;
        }
        if(is_primitive(testelement, r1, precision, reduced_poly, currRing)){
            DEBUG_PRINT(("is_primitive(testelement, r1, precision, reduced_poly, currRing)\n"));
            n_Delete(&bound,real);
            delete latticeNF;
            delete LLLbasis;
            return reduced_poly; //return NULL; //correct   NOTE: ????
        }
        delete testelement;
    }
    delete LLLbasis;
    
    /*//enumerate all elements less bound and check
    bigintmat * elementbasis = latticeNF->enumerate_next(bound);
    bigintmat * testelement = latticeNF->get_lattice_element(elementbasis);
    while(elementbasis != NULL){
        if(is_primitive(testelement, r1, precision, reduced_poly, currRing)){
            n_Delete(&bound,real);
            delete latticeNF;
            delete LLLbasis;
            return NULL;//correct
        }
        delete elementbasis;
        elementbasis = latticeNF->enumerate_next();
        delete testelement;
        testelement = latticeNF->get_lattice_element(elementbasis);
    }
    delete elementbasis;
    delete testelement;
    //*/
    
    //None found, so delete unused stuff
    n_Delete(&bound,real);
    delete latticeNF;
    return NULL;//none found
}

bool is_primitive(bigintmat * element,int r1, int precision, poly out, const ring polyring){
//     return true;
    
    //check if all numbers in minkowski element are diffrent
    coeffs coef = element->basecoeffs();
    for(int i=0;i<element->rows()-1;i++){
        for(int j=i+1;j<element->rows();j++){
            if(n_Equal(element->view(i),element->view(j),coef)){
                return false;
            }
        }
    }
    
    //generate roots from minkowski element
    LongComplexInfo paramComp = {si_min (precision, 32767),si_min (precision, 32767),(const char*)"i"};
    
    if(r1>element->rows()){
        return false;
    }
    int deg = element->rows();
    coeffs comp = nInitChar(n_long_C, &paramComp);
    nMapFunc f = n_SetMap(coef, comp);
    
    number * roots = new number[deg];
    for(int i = 0; i<r1;i++){
        roots[i] = f(element->view(i),coef,comp);
    }
    if(r1<deg){
        number imagunit=ngcInit(1,comp);
        number two = n_Init(2,comp);
        number sqrt2 = squareroot(two,comp,precision+10);
        for(int i=r1;i<deg;i=i+2){
            number minkowreal = f(element->view(i),coef,comp);
            number minkowimag = f(element->view(i+1),coef,comp);
            number realpart = n_Div(minkowreal,sqrt2,comp);
            number impart = n_Div(minkowimag,sqrt2,comp);
            n_InpMult(impart,imagunit,comp);
            roots[i] = n_Add(realpart,impart,comp);
            roots[i+1] = n_Sub(realpart,impart,comp);
            n_Delete(&minkowreal,comp);
            n_Delete(&minkowimag,comp);
            n_Delete(&realpart,comp);
            n_Delete(&impart,comp);
        }
        n_Delete(&imagunit,comp);
        n_Delete(&two,comp);
        n_Delete(&sqrt2,comp);
    }
    number * polycoef = new number[deg+1];
    
    nMapFunc g = n_SetMap(comp, coef);
    
    DEBUG_VAR(deg);
    for(int i=0;i<deg;i++){
        DEBUG_VAR(i);
        number temp = elementary_symmetric_function(roots,deg,i+1,0,comp);
        number temp2 = ngcRePart(temp,comp);
        n_Delete(&temp,comp);
        temp = round(temp2,comp);
        polycoef[i] = g(temp,comp,coef);
        DEBUG_N(polycoef[i]);
        n_Delete(&temp,comp);
        n_Delete(&temp2,comp);
    }
    polycoef[deg] = n_Init(1,coef);
    poly testpoly = numbers2poly(polycoef,deg,coef,polyring);
    
    int var = p_IsUnivariate(testpoly, polyring);
    poly diff_testpoly = p_Diff(testpoly, var, polyring);
    
  /*  
    
    
//     poly test = singclap_gcd(p_Copy(testpoly,polyring), p_Copy(diff_testpoly,polyring), polyring);
    
    
    
//     int err;
// 
//     // calling a kernel function via the interpreter interface
//     sleftv zu; memset(&zu,0,sizeof(zu));
//     sleftv arg; memset(&arg,0,sizeof(zu));
//     arg.rtyp=STRING_CMD;
//     arg.data=omStrDup("9+3");
//     err=iiExprArith1(&zu,&arg,TYPEOF_CMD);
// 
//     printf("interpreter returns %d\n",err);
//     if (err)
//         errorreported = 0; // reset error handling
//     else
//         printf("typeof returned type %d, >>%s<<\n",zu.Typ(),zu.Data());
// 
//     // clean up zu:
//     zu.CleanUp();
// 
//     return 0;
    
    
    
    idhdl myGCD=ggetid("myGCD");
    
    if(myGCD == NULL){
        WerrorS("myGCD not found\n");
        return NULL;
    }
    
    poly testpoly_copy = p_Copy(testpoly,polyring);
    
    leftv arg = new sleftv();
    arg->rtyp = POLY_CMD;
    arg->data = (void*) testpoly_copy;
    
//     leftv arg2 = new sleftv();
//     arg2->rtyp = POLY_CMD;
//     arg2->data = (void*) p_Copy(diff_testpoly,polyring);
//     
//     arg->next = arg2;
    
    if(iiMake_proc(myGCD,NULL,arg)){
        WerrorS("Error: Can't use myGCD\n");
        return NULL;
    }
   
    printf("iiMake_proc: myGCD type %d, >>%s<<\n", iiRETURNEXPR.Typ(), (char *)iiRETURNEXPR.Data());
    
    poly test = (poly) iiRETURNEXPR.Data();
    
    
    iiRETURNEXPR.CleanUp();
    
    
    
    
    */
    
    
//     number * pcoeffs = poly2numbers(test,polyring, coef);
//     int degg = (int) p_Totaldegree(test,polyring);
//     for(int i=0; i<=degg; i++) {
//         n_Print(pcoeffs[i],coef);
//         PrintS("\n");
//     }
    
//     DEBUG_PRINT(("trz\n"));
//     DEBUG_VAR(p_Deg(test,polyring));
//     if(p_Deg(test,polyring)<=1){//test if degree less than 1
//         //Delete things
//         DEBUG_PRINT(("p_Deg(test,polyring)<=1\n"));
        *out = *numbers2poly(polycoef,deg,coef,polyring);
        return true;
//     }
//     DEBUG_PRINT(("sdfs\n"));
//     //delete things
//     return false;   
}

number * poly2numbers(poly gls, ring polyring, coeffs coef){
    DEBUG_PRINT(("Start poly2numbers\n"));
    if(gls == NULL){
        WerrorS("No Input!");
        return NULL;
    }
     

    DEBUG_PRINT(("degree\n"));
    int deg = (int) p_Totaldegree(gls,polyring);
    DEBUG_VAR(deg);
    
    
    number * pcoeffs = (number *) omAlloc0((deg+1) * sizeof(number));
    
    
    DEBUG_PRINT(("univ\n"));
    int vpos = p_IsUnivariate(gls, polyring);
    if(vpos <= 0){
        WerrorS("not univariate");
        return NULL;
    }
    poly piter = gls;
    
    nMapFunc f = n_SetMap(polyring->cf,coef);
    DEBUG_PRINT(("iterate\n"));
    for (int i= deg; i >= 0; i-- ) {
        if ( piter && p_Totaldegree(piter,polyring) == i ) {
            number temp = n_Copy( p_GetCoeff( piter, polyring), polyring->cf );
            pcoeffs[i] = f(temp,polyring->cf,coef);
            n_Delete(&temp,polyring->cf);
            pIter( piter );
        } else {
            pcoeffs[i]= n_Init(0,coef);
        }
        DEBUG_VAR(i);
        DEBUG_N(pcoeffs[i]);
    }
    DEBUG_PRINT(("finished\n"));
    return pcoeffs;
}

poly numbers2poly(number * univpol, int deg, coeffs coef, ring polyring){
    DEBUG_PRINT(("Start numbers2poly\n"));
    poly result= NULL;
    poly ppos;
    nMapFunc f = n_SetMap(coef,polyring->cf);
    DEBUG_PRINT(("set coefficients\n"));
    for ( long i= deg; i >= 0; i-- ){
        if ( !n_IsZero(univpol[i],coef) ){
            poly p= p_One(polyring);
            //pSetExp( p, var+1, i);
            p_SetExp( p, 1, i, polyring);
            p_SetCoeff( p, f(univpol[i],coef,polyring->cf),polyring);
            p_Setm( p ,polyring);
            if (result) {
                ppos->next=p;
                ppos=ppos->next;
            } else {
                result=p;
                ppos=p;
            }
        }
    }
    if (result!=NULL) p_Setm( result, polyring );

    return result;
}



///////////////////////////////////////
//              T2-norm             ///
///////////////////////////////////////
number t2norm(bigintmat * elt){
    int n = elt->rows();
    coeffs coef = elt->basecoeffs();
    number sum = n_Init(0,coef);
    for(int i=1;i<=n;i++){
        number mult = n_Mult(elt->view(i,1),elt->view(i,1),coef);
        n_InpAdd(sum,elt->view(i,1),coef);
        n_Delete(&mult,coef);
    }
    return sum;
}

number t2norm(number * pol, int deg, coeffs coef, int precision){
    DEBUG_PRINT(("t2 norm numbers\n"));
    DEBUG_PRINT(("Input check\n"));
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
    if(precision>32767){
        precision = 32767;
    }
    if ( !(nCoeff_is_R(coef) || nCoeff_is_long_R(coef) || nCoeff_is_long_C(coef)) ){
        setGMPFloatDigits( precision+6,precision+6);
    }
    ring CurrOldRing = rCopy(currRing);//need to change currRing, because rootContainer uses the coef of it
    char* n[] = {(char*)"x"};///
    ring newring = rDefault(coef, 1, n);
    rChangeCurrRing(newring);
    
    DEBUG_PRINT(("initialize rootContainer\n"));
    number * rootpoly = new number[deg+1];//doesn't need to be deleted since rootContaine delete it...
    for(int i=0;i<=deg;i++){
        rootpoly[i] = n_Copy(pol[i],coef);
    }
    rootContainer * rootcont= new rootContainer();
    rootcont->fillContainer( rootpoly, NULL, 1, deg, rootContainer::onepoly, 1 );
    DEBUG_PRINT(("rootcont->solver(1);\n"));
    rootcont->solver(1);
    int number_roots = rootcont ->getAnzRoots();
    if(number_roots != deg){
        WerrorS("something went wrong: \n\tnot all roots found\n");
        return NULL;
    }
    LongComplexInfo paramComp = {si_min (precision+6, 32767),si_min (precision+8, 32767),(const char*)"i"};
    
    
    coeffs comp = nInitChar(n_long_C, &paramComp);
    number sumcomp = n_Init(0,comp);
    for(int j=0; j<deg; j++){
        number a = n_Copy((number)(rootcont->getRoot(j)),comp);
        number re = ngcRePart(a, comp);
        n_InpMult(re,re,comp);
        number im = ngcImPart(a, comp);
        n_InpMult(im,im,comp);
        n_Delete(&a,comp);
        a = n_Add(re,im,comp);
        n_InpAdd(sumcomp,a,comp);
        n_Delete(&a,comp);
    }
    delete rootcont;
    rChangeCurrRing(CurrOldRing);
    rDelete(newring);
    LongComplexInfo paramReal = {si_min (precision, 32767), si_min (precision, 32767), (const char*)"i"};
    coeffs real = nInitChar(n_long_R, &paramReal);
    nMapFunc f = n_SetMap(comp, real);
    number norm_real = f(sumcomp,comp,real);
    
    nMapFunc g = n_SetMap(real, coef);
    number norm = g(norm_real,real,coef);
    n_Delete(&sumcomp,comp);
    n_Delete(&norm_real,real);
    nKillChar(comp);
    nKillChar(real);
    return norm;
}

number t2norm(poly polynom, ring polyring, coeffs coef, int precision){
    DEBUG_PRINT(("t2 norm poly\n"));
    
    int deg = (int) p_Totaldegree(polynom,polyring);
    number * univpol = poly2numbers(polynom,polyring,coef);

    if(deg == -1){
        WerrorS("Not an univariate polynomial");
        return NULL;
    }
    number norm = t2norm(univpol,deg,coef,precision);
    for(int i=0;i<deg;i++){
        n_Delete(&univpol[i],coef);
    }
    delete[] univpol;
    return norm;
}

number elementary_symmetric_function(number * roots, int deg, int si, int lower_bound, coeffs coef){
    if (si < 1 || si > deg){
        WerrorS("it should 1<=si<=deg");
        return n_Init(0,coef);
    }
    if (lower_bound<0 || lower_bound>=deg){
        WerrorS("it should 0<=lower_bound<deg");
        return n_Init(0,coef);
    }
    
    if (si == 1){
        number sum = n_Init(0,coef);
        for(int i=lower_bound;i<deg;i++){
            n_InpAdd(sum,roots[i],coef);
        }
        return sum;
    } else {
        number sum = n_Init(0,coef);
        for(int i=lower_bound;i<deg+1-si;i++){
            number temp = elementary_symmetric_function(roots,deg,si-1,i+1,coef);
            n_InpMult(temp,roots[i],coef);
            n_InpAdd(sum,temp,coef);
            n_Delete(&temp,coef);
        }
        return sum;
    }
}


