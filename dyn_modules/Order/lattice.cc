#include <libpolys/coeffs/bigintmat.h>
#include "lattice.h"
#include "kernel/febase.h"  // for Print, WerrorS
#include "libpolys/coeffs/numbers.h" 
#include "libpolys/coeffs/coeffs.h"
#include "Singular/ipid.h"
#include <iostream>



///////////////////////////////////////
//         Debugging Stuff          ///
///////////////////////////////////////

#define DEBUG_PRINTS 1 //remove this line to disable debugging
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
 
lattice::lattice(bigintmat* basis){
    DEBUG_BLOCK(false);
    DEBUG_PRINT(("Creating new lattice..."));
    this->basis = bimCopy(basis);
    coef = basis->basecoeffs();
    c = n_Div(n_Init(3, coef),n_Init(4, coef),coef);// 3/4;
    B = (number*) omAlloc((basis->cols()+1) * sizeof(number));
    b = NULL;
    b_star = NULL;
    H = NULL;
    my = NULL;
    relat = false; //otherwise error?
    DEBUG_PRINT(("Done\n"));
}

lattice::~lattice() {
    DEBUG_BLOCK(false);
    DEBUG_PRINT(("Deleting lattice..."));
    delete basis;
    delete b;
    delete b_star;
    delete H;
    delete my;
    omFree(B);
    DEBUG_PRINT(("Done\n"));
}


///////////////////////////////////////
//               LLL                ///
///////////////////////////////////////

bool lattice::LLL(){
    DEBUG_BLOCK(true);
    DEBUG_PRINT(("Start LLL\n"));
    
    delete b;
    delete b_star;
    delete H;
    delete my;

    
    b      = bimCopy(basis);
    b_star = new bigintmat(b->rows(),b->cols(),coef);
    H      = new bigintmat(b->rows(),b->cols(),coef);
    my     = new bigintmat(b->rows(),b->cols(),coef);
    
    for(int i=1; i<=b->cols(); i++) {
            B[i] = n_Init(0,coef);
    }
    
    if (b->cols() == 0) {
        return true;
    }
    if (b->cols() == 1) {
        return false;
    }
    
    DEBUG_BIM(b);
    
    ///Initialize
    DEBUG_PRINT(("Initialize\n"));
    int k = 2;
    int k_max = 1;
    
    // b*_1 <- b_1
    for(int i=1; i<=b->rows(); i++) {
        b_star->set(i,1,b->view(i,1),coef);
    }
    
    // B_1 <- b_1 * b_1
    B[1]= n_Init(0,coef);
    for(int i = 1; i <= b->rows(); i++){
         n_InpAdd(B[1], n_Mult(b_star->view(i,1), b_star->view(i,1), coef), coef);
    }
    
    // H <- I_n
    H->one();
    
    do{
        ///Incremental Gram-Schmidt
        DEBUG_PRINT(("Incremental Gram-Schmidt\n"));
        DEBUG_VAR(k);
        DEBUG_VAR(k_max);
        DEBUG_BIM(b);
        DEBUG_BIM(b_star);
        DEBUG_BIM(my);
        DEBUG_N(B[1]);
        DEBUG_N(B[2]);
        DEBUG_N(B[3]);
        if(k > k_max){
            k_max = k;
            if(gram_schmidt(k)){
                return true;
            }
        }
        ///Test LLL condition
        DEBUG_PRINT(("Test LLL condition\n"));
        while(true){            
            RED(k,k-1);
            
            n_Print(B[3],coef); //<-most important line in the world
            Print("\n");
            
            // if((B[k] < (c- my*my)*B[k-1]))
            if(n_Greater(n_Mult(n_Sub(c, n_Mult(my->view(k,k-1), my->view(k,k-1), coef), coef), B[k-1], coef), B[k], coef)){
                SWAP(k,k_max);
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
    }while(k <= b->cols());
    DEBUG_VAR(k);
    DEBUG_BIM(b);
    DEBUG_BIM(b_star);
    DEBUG_BIM(my);
    DEBUG_N(B[1]);
    DEBUG_N(B[2]);
    DEBUG_N(B[3]);
    DEBUG_PRINT(("End of LLL\n"));
    return false;
}

void lattice::RED(int k, int l){
    DEBUG_PRINT(("Start RED with k=%d and l=%d\n",k,l));
    number n_1div2    = n_Div(n_Init( 1,coef),n_Init(2,coef),coef);
    number n_neg1div2 = n_Div(n_Init(-1,coef),n_Init(2,coef),coef);
    number n_kl = my->view(k,l);
    
    // if(|my_kl| > 1/2)
    if (n_Greater (n_kl,n_1div2,coef) || n_Greater (n_neg1div2,n_kl,coef)) { 
        number n_klplus1div2 = n_Add(n_kl, n_1div2, coef);
        
        //round(my_kl);
        number q = n_IntDiv(n_GetNumerator(n_klplus1div2,coef),n_GetDenom(n_klplus1div2,coef),coef); 
        DEBUG_N(q);
        b->addcol(k,l,n_Neg(q,coef),coef);
        
        if(relat) {
            H->addcol(k,l,n_Neg(q,coef),coef);
        }
        
        my->set(k,l,n_Sub(my->view(k,l),q,coef),coef);
        
        //my_ki <- my_ki - q*my_li
        for(int i=1;i<=l-1;i++){
            my->set(k,i,n_Sub(my->view(k,i), n_Mult(q, my->view(l,i),coef), coef), coef);
        }
    }
    DEBUG_PRINT(("End of RED\n"));
}

void lattice::SWAP(int k, int k_max){
    DEBUG_PRINT(("Start SWAP with k=%d and k_max=%d\n",k,k_max));   
    b->swap(k,k-1);
    
    if(relat) {
        H->swap(k,k-1);
    }
    
    for(int j = 1; j <= k-2; j++){
        DEBUG_PRINT(("j=%d\n",j));
        number my_kj = my->get(k,j);
        my->set(k,j,my->view(k-1,j),coef);
        my->set(k-1,j,my_kj,coef);
    }
    
    number my_ = my->view(k,k-1);
    
    //B_ <- B[k] + my_*my_*B[k-1];
    number B_ = n_Add(B[k], n_Mult(n_Mult(my_,my_,coef), B[k-1], coef), coef); 
    
    my->set(k,k-1,n_Div(n_Mult(my_, B[k-1], coef), B_, coef), coef);
    
    bigintmat *b_ = new bigintmat(b->rows(),1,coef);
    b_star->getcol(k-1, b_);
    
    bigintmat *b_star_k = new bigintmat(b->rows(),1,coef);
    b_star->getcol(k, b_star_k);
    
    b_star->setrow(k-1, bimAdd(b_star_k, bimMult(b_, my_, coef)));
    
    b_star->setrow(k, bimSub(bimMult(b_, n_Div(B[k], B_, coef),coef), bimMult( b_star_k, my->view(k,k-1), coef)));
    
    delete b_;
    delete b_star_k;
    
    // B[k] <- (B[k]* B[k-1])/B_;
    B[k] = n_Div( n_Mult( B[k], B[k-1], coef), B_, coef);
    
    B[k-1] = n_Copy(B_, coef);
    for(int i = k+1; i <= k_max; i++){
        DEBUG_PRINT(("i=%d\n",i));
        number t = my->get(i,k);
        my->set(i,k,n_Sub(my->get(i,k-1), n_Mult(my_, t, coef), coef), coef);
        my->set(i,k-1, n_Add(t, n_Mult(my->view(k,k-1), my->view(i,k), coef), coef), coef);
    }
    DEBUG_PRINT(("End of SWAP\n"));
}

bool lattice::gram_schmidt(int k) {
    DEBUG_PRINT(("Start gram_schmidt(%d)\n",k));
    
    // b*_k <- b_k
    for(int i=1; i<=b->rows(); i++) {
        b_star->set(i,k,b->view(i,k),coef);
    }
    
    
    for(int j=1; j<k; j++){
        // my_kj <- b_k * b*_j / B_j
        number my_ = n_Init(0,coef);
        for(int i = 1; i <= b->rows(); i++){
            // my_ += b*_ij * b_ik;
            n_InpAdd(my_,n_Mult(b_star->view(i,j), b->view(i,k),coef),coef); 
        }
        my_ = n_Div(my_,B[j],coef);
        my->set(k,j,my_,coef);
        
        // b*_k <- b*_k - my_kj * b*j
        b_star->addcol(k,j,n_Neg(my->get(k,j),coef),coef);
//         for(int i=1; i<=b->rows(); i++) {
//             number b_star_ = n_Sub(b_star->view(i,k), n_Mult(my->view(k,j),b_star->view(i,j),coef),coef);
//             b_star->set(i,k,b_star_,coef);
//         }
    }
    
    // B_k <- b*_k * b*_k
    B[k] = n_Init(0,coef);
    for(int i = 1; i <= b_star->rows(); i++){
         n_InpAdd(B[k],n_Mult(b_star->view(i,k),b_star->view(i,k),coef),coef);
    }
    
    if ( n_IsZero(B[k],coef)){
        Werror("did not form a basis\n");
        DEBUG_PRINT(("End of gram_schmidt(%d)\n",k));
        return true;
    } else {
        DEBUG_PRINT(("End of gram_schmidt(%d)\n",k));
        return false;
    }
}

bool lattice::gram_matrix(int k){
    number* a = new number[k];
    for(int j = 1; j<k;j++){
        a[j] = n_Init(0,coef);
        for(int i =1; i<=b->rows(); i++){
            a[j] = n_Add(a[j],n_Mult(b->view(i,k),b->view(i,j),coef),coef);//a[j] += b->view(i,k) * b->view(i,j);
        }
        for(int i =1; i<=j-1; i++){
            a[j] = n_Add(a[j],n_Mult(b->view(i,j),a[i],coef),coef);//a[j] += my->view(j,i) * a[i];
        }
        my->set(k,j,n_Div(a[j],B[j],coef),coef);
    }
    B[k]=n_Init(0,coef);
    for(int i =1; i<=b->rows(); i++){
        B[k]=n_Add(B[k],n_Mult(b->view(i,k),b->view(i,k),coef),coef);//B[k] += b->view(i,k) * b->view(i,k);
    }
    for(int i =1; i<=k-1; i++){
        B[k] = n_Add(B[k],n_Mult(my->view(k,i),a[i],coef),coef);//B[k] += my->view(k,i) * a[i];
    }
    if(B[k] == 0){
        Werror("did not form a basis\n");
        return false;
    }
    return true;
}


///////////////////////////////////////
//               MLLL               ///
///////////////////////////////////////

// bool lattice::MLLL(){
//     /// Algorithm 2.6.8.
//     
//     if (b->cols() == 0) return true;
//     if (b->cols() == 1) {
//         return false;
//     }
//     ///Initialize
//     int k = 2;
//     int k_max = 1;
//     b_star = bigintmat(b->rows(),b->cols(),coef);
//     
//     bigintmat b_;
//     b->getcol(1,&b_);
//     b_star->setcol(1,&b_);
//     delete[] &b_;
//     for(int i = 1; i <= b->rows(); i++){
//          B[1] = n_Add(B[1], n_Mult(b_star->view(i,1), b_star->view(i,1), coef), coef);//B[1] += b_star->view(i,1) * b_star->view(i,1);
//     }
//     if (relat){
//         do{
//             ///Incemental Gram-Schmidt
//             if(k > k_max){
//                 //not for Gram matrix
//                 k_max = k;
//                 gram_schmidt_MLLL(k);
//             }
//             ///Test LLL condition
//             while(1){
//                 RED_rel(k,k-1);
//                 if (n_Greater ( n_Mult(n_Sub(c, n_Mult(my->view(k,k-1), my->view(k,k-1), coef), coef), B[k-1], coef), B[k], coef)){//(B[k] < (c- my*my)*B[k-1]){
//                     SWAPG_rel(k,k_max);
//                     if (k>2) k--;
//                 } else {
//                     for( int l= k-2; l>0; l--){
//                         RED_rel(k,l);
//                     }
//                     k++;
//                     break;
//                 }
//             }
//             ///Finished ?
//         }while(k<=b->cols());
//     } else {
//         do{
//             ///Incemental Gram-Schmidt
//             if(k > k_max){
//                 //not for Gram matrix
//                 k_max = k;
//                 gram_schmidt_MLLL(k);
//             }
//             ///Test LLL condition
//             while(1){
//                 RED(k,k-1);
//                 if (n_Greater ( n_Mult(n_Sub(c, n_Mult(my->view(k,k-1), my->view(k,k-1), coef), coef), B[k-1], coef), B[k], coef)){//(B[k] < (c- my*my)*B[k-1]){
//                     SWAPG(k,k_max);
//                     if (k>2) {
//                         k--;
//                     }
//                 } else {
//                     for( int l= k-1; l>0; l--){
//                         RED(k,l);
//                         k++;
//                         break;
//                     }
//                 }
//             }
//             ///Finished ?
//         }while(k<=b->cols());
//     }
//     b_star = NULL;
//     
//     
//     return false;
// }
// 
// void lattice::gram_schmidt_MLLL(int k){
//     number my_;
//     for(int j = 1; j < k; j++){
//         if (! n_IsZero(B[j],coef)){
//             my_ = n_Init(0,coef);
//             for( int i = 1; i <= b->rows(); i++){
//                 my_ = n_Add(my_,n_Mult(b_star->view(i,j),b->view(i,k),coef),coef);//my_ += b_star->view(i,j) * b->view(i,k);
//             }
//             my_ = n_Div(my_,B[j],coef);
//             my->set(k,j,my_,coef);
//         } else {
//             b_star->set(k,j,0,coef);
//         }
//     }
//     bigintmat b_;
//     b->getcol(k,&b_);
//     b_star->setcol(k,&b_);
//     for(int j = 1; j <= k-1; j++){
//         b_star->addcol(k,j,n_Mult(n_Init(-1,coef),my->view(k,j),coef),coef);
//     }
//     B[k] = n_Init(0,coef);
//     for(int i = 1; i < b_star->rows(); i++){
//         B[k] = n_Add( n_Mult(b_star->view(i,k), b_star->view(i,k), coef), B[k], coef);//B[k] += b_star->view(i,k) * b_star->view(i,k);
//     }
// }

// void lattice::SWAPG(int k, int k_max){
//     number my_, B_, t;// my_k, my_i;
//     bigintmat b_, bk;
//     b->swap(k,k-1);
//     for(int j = 1; j <= k-2; j++){
//         my_ = my->get(k,j);
//         my->set(k,j,my->view(k-1,j),coef);
//         my->set(k-1,j,my_,coef);
//     }
//     my_=my->get(k,k-1);
//     B_ = n_Add( B[k], n_Mult( n_Mult( my_, my_, coef), B[k-1], coef), coef);//B[k] + my_ * my_ * B[k-1];
//     if(n_IsZero(B_,coef)){
//         t = B[k];
//         B[k] = B[k-1];
//         B[k-1] = t;
//         b_star->swap(k,k-1);
//         for(int i = k+1; i <= k_max; i++){
//             t = my->get(i,k);
//             my->set(i,k,my->view(i,k-1),coef);
//             my->set(i,k-1,t,coef);
//         }
//     } else if (n_IsZero(B[k],coef)){
//         B[k-1] = B_;
//         b_star->colskalmult(k-1,my_,coef);
//         my->set(k,k-1,n_Div(n_Init(1,coef),my_,coef),coef);
//         for(int i = k+1; i <= k_max;i++){
//             my->set(i,k-1,n_Div(my->view(i,k-1),my_,coef),coef);
//         }
//     } else {
//         t = n_Div(B[k-1],B_,coef);
//         my->set(k,k-1,n_Mult(my_,t,coef),coef);
//         b_star->getcol(k-1, &b_);
//         b_star->colskalmult(k-1,my_,coef);
//         b_star->addcol(k-1,k,n_Init(1,coef),coef);
//         b_star->colskalmult(k,my->view(k,k-1),coef);
//         b_star->getcol(k-1, &bk);
//         b_star->setcol(k,bimSub(bimMult( &b_, n_Div(B[k],B_,coef), coef),&bk));
//         B[k] = n_Mult(B[k],t,coef);
//         B[k-1] = B_;
//         for(int i=k+1; i<=k_max;i++){
//             t=my->get(i,k);
//             my->set(i,k,n_Sub(my->view(i,k-1),n_Mult(my_,t,coef),coef),coef);
//             my->set(i,k-1,n_Add(t,n_Mult(my->view(k,k-1),my->view(i,k),coef),coef),coef);
//         }
//     }
// }
// 
// void lattice::INSERT(int k, int i){
//     if(k > i){
//         for(int j = k; j > i; j--){
//             b->swap(j,j-1);
//         }
//     }
// }



///////////////////////////////////////
//             Enumerate            ///
///////////////////////////////////////

// bool lattice::enumerate(number a, bigintmat* liste){
//     //Quadratic Supplement
//     if(b->cols() != b->rows()) return true;
//     bigintmat Q;
//     Q.copy(b);
//     for(int i = 1; i<Q.cols();i++){
//         for( int j=i+1; j<Q.cols();j++){
//             Q.set(j,i,Q.view(i,j),coef);
//             Q.set(i,j,nExactDiv(Q.view(i,j),Q.view(i,i)),coef);
//         }
//         for(int m=i+1; m<Q.cols();m++){
//             for(int n=i+1; n<Q.cols();n++){
//                 Q.set(m,n,nSub(Q.view(m,n),nMult(Q.view(m,i),Q.view(i,n))),coef);
//             }
//         }
//     }
//     for(int i = 2; i<Q.cols();i++){
//         for(int j = 1; j<i;j++){
//             Q.set(i,j,nInit(0),coef);
//         }
//     }
//     
//     //backtracking for elements
//     int elements = 0;
//     int index = 1;
//     number check;
//     number x[Q.cols()+1] = new number;//////////
//     number grenze[Q.cols()+1] = new number;
//     for(int i = 1; i<=Q.cols(); i++){
//         x[i] = n_Init(0,coef);
//         grenze[i] = a;
//     }
//     while (index <= Q.cols()) {
//         //update check
//         check = n_Init(0,coef);
//         for(int i=index + 1;i<=Q.cols();i++){
//             check = n_Add(check,n_Mult(x[i],Q.view(index,i),coef),coef);
//         }
//         check = n_Add(x[index],check,coef);
//         check = n_Mult(Q.view(index,index),n_Mult(check,check,coef),coef);
//         //check check
//         if (n_Greater(check,grenze[index],coef)){
//             grenze[index] = n_Init(0,coef);
//             for(int i = 1; i<=index;i++){ //perhaps only x[index] = n_Init(0,coef);
//                 x[i] = n_Init(0,coef);
//             }
//             index++;
//             if (n_GreaterZero(x[index], coef)){
//                 x[index] = n_Mult(n_Init(-1,coef),x[index],coef); //x_i=-x_i
//             } else {
//                 x[index] = n_Add(n_Mult(n_Init(-1,coef),x[index],coef),N_Init(1,coef),coef);//x_i=-x_i+1
//             }
//         } else if(index == 1){
//             //hinzufügen zu liste
//             elements++;
//             liste = (bigintmat*) omAlloc(elements + 1);
//             liste[elements]->copy(x);
//             x[1] = n_Add(x[1],N_Init(1,coef),coef);
//         } else {
//             index--;
//             grenze[index] = n_Sub(grenze[index+1],check,coef);
//         }
//     }
//     return false;
// }


///////////////////////////////////////
//               Getter             ///
///////////////////////////////////////

bigintmat* lattice::get_basis() {
    return bimCopy(basis);
}

bigintmat* lattice::get_reduced_basis() {
    return bimCopy(b);
}


///////////////////////////////////////
//               Setter             ///
///////////////////////////////////////

void lattice::set_c(number a){
    if (n_Greater (n_Mult(a,n_Init(4,coef),coef),n_Init(1,coef),coef) && n_Greater (n_Init(1,coef),a,coef)) {//(1<4*a && a<1){
        c = n_Copy(a, coef);
    } else if(n_IsOne(a, coef)){
        c = n_Copy(a, coef);
        Werror("probably not in polynomial time\n");
    } else {
        Werror("not a possible value\n");
    }
}


