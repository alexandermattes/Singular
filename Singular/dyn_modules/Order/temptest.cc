#include <libpolys/coeffs/bigintmat.h>
#include "temptest.h"
#include <reporter/reporter.h>
#include "libpolys/coeffs/numbers.h" 
#include "libpolys/coeffs/coeffs.h"
#include "Singular/ipid.h"
#include"reporter/reporter.h"  // for Print, WerrorS
#include "lattice.h"
#include <iostream>
#include <math.h> 
#include <vector>
#include <utility>
#include <stdlib.h>

using namespace std;

//Temporary testfunction to play arround with new functions
//NOTE: remove later
number temp_test(number a) {    
    number c = nCopy(a);
    
    PrintS("\n\n\n");
    test_Enumerate();
    PrintS("\n\n\n");
    test_Minkowski();
    PrintS("\n\n\n");
    test_LLL();
    
    return c;    
}


void test_LLL() {
    coeffs coef = nInitChar(n_Q,NULL);
    PrintS("Test LLL\n");
    
    //Vector containing test matrices
    std::vector<bigintmat*> vec;
    
    {
        int ar[] = {1,-1,3,1,0,5,1,2,6};
        int n = 3;
        int m = 3;
        bigintmat * temp = new bigintmat(n,m,coef);
        for(int i=0; i<(n*m); i++) {
            temp->set(i,n_Init(ar[i],coef), coef);
        }
        vec.push_back(temp);
    }
    
    {
        int ar[] = {1,9,1,2,1,8,8,3,7,4,5,1,2,6,7,1};
        int n = 4;
        int m = 4;
        bigintmat * temp = new bigintmat(n,m,coef);
        for(int i=0; i<(n*m); i++) {
            temp->set(i,n_Init(ar[i],coef), coef);
        }
        vec.push_back(temp);
    }
    
//     {
//         int ar[] = {0};
//         int n = 0;
//         int m = 0;
//         bigintmat * temp = new bigintmat(n,m,coef);
//         for(int i=0; i<(n*m); i++) {
//             temp->set(i,n_Init(ar[i],coef), coef);
//         }
//         vec.push_back(temp);
//     }
    
//     {
//         int ar[] = {1};
//         int n = 1;
//         int m = 1;
//         bigintmat * temp = new bigintmat(n,m,coef);
//         for(int i=0; i<(n*m); i++) {
//             temp->set(i,n_Init(ar[i],coef), coef);
//         }
//         vec.push_back(temp);
//     }
        
    {
        int ar[] = {1,-1,3,1,0,5,1,2,6};
        int n = 3;
        int m = 3;
        bigintmat * temp = new bigintmat(n,m,coef);
        for(int i=0; i<(n*m); i++) {
            temp->set(i,n_Init(ar[i],coef), coef);
        }
        vec.push_back(temp);
    }
    
    {
        int ar[] = {1,2,3,31,41,51,101,201,301};
        int n = 3;
        int m = 3;
        bigintmat * temp = new bigintmat(n,m,coef);
        for(int i=0; i<(n*m); i++) {
            temp->set(i,n_Init(ar[i],coef), coef);
        }
        vec.push_back(temp);
    }
    
    {
        int ar[] = {10000,4545545,4557454,465445241};
        int n = 2;
        int m = 2;
        bigintmat * temp = new bigintmat(n,m,coef);
        for(int i=0; i<(n*m); i++) {
            temp->set(i,n_Init(ar[i],coef), coef);
        }
        vec.push_back(temp);
    }

//     //Random matrices
//     srand(1234);
//     for(int t = 0; t < 100; t++) {
//         int matrix_size = rand() % 10 + 1;
//         bigintmat * matrix = new bigintmat(matrix_size,matrix_size,coef);
//         matrix->one();
//         for(int ops = 0; ops < matrix_size*5; ops++) {
//             matrix->addcol(rand() % matrix_size + 1
//                           ,rand() % matrix_size + 1
//                           ,n_Init(rand() % 5 + 1,coef),coef);
//         }
//         vec.push_back(matrix);
//     }
    
    
    for(unsigned i= 0; i<vec.size(); i++) {
        PrintS("\nInput:\n");
// // //         cout << orthogonality_defect(vec[i]) << '\n';
        vec[i]->Print();
        PrintS("\n");
        lattice * l = new lattice(vec[i]);
        number c = NULL;
        l->LLL(c,true,false,false);
        bigintmat * reduced = l->get_reduced_basis();
        bigintmat * H = l->get_transformation_matrix();
        PrintS("Output:\n");
//         cout << orthogonality_defect(reduced) << '\n';
        reduced->Print();
        assume(bimSub(reduced,bimMult(vec[i],H))->isZero());
        number H_det = H->det();
        assume(n_IsOne(H_det,coef) || n_IsOne(n_InpNeg(H_det,coef),coef));

//         bigintmat * gram_matrix = l->get_gram_matrix();
//         lattice * l_gram = new lattice(vec[i],true);
//         l_gram->LLL(c,true,false,true);
//         bigintmat * H_gram = l_gram->get_transformation_matrix();
//         assume(bimSub(H_gram,H)->isZero());
    }
    PrintS("\n");
}


void test_Enumerate() {   
    PrintS("Test Enumerate\n");
    
    coeffs coef = nInitChar(n_Q,NULL);
    number c = n_Init(8,coef);
    
    //cout << "CoeffType of currRing: " << getCoeffType(currRing->cf) << '\n';
    //cout << "n_Greater(n_Init(3,coef),n_Init(3,coef),coef): " << n_Greater(n_Init(3,coef),n_Init(3,coef),coef) << '\n';
    
    bigintmat *m = new bigintmat(3,3,coef);
    for(int i=1; i<=3; i++) {
        m->set(i,i,n_Init(1,coef), coef);
    }
    m->set(1,2,n_Init(1,coef), coef);
    m->set(2,3,n_Init(1,coef), coef);
    m->set(3,1,n_Init(1,coef), coef);
    
    m->Print();
    PrintS("\n");
    
    lattice* l = new lattice(m);
    //l->LLL();
    bigintmat* enumer = NULL;// = new bigintmat(3,,coef);
    enumer = l->enumerate_all(c);
    //bigintmat *reduced = l->get_reduced_basis();
    //reduced->Print();
    if(enumer !=NULL){
        enumer->transpose()->Print();
        PrintS("\n");
    }
    PrintS("new number\n");
    enumer = l->enumerate_next( c);
    if(enumer !=NULL){
        enumer->Print();
        PrintS("\n");
    }
    PrintS("new number and bigintmat\n");
    bigintmat* x= new bigintmat(3,1,coef);
    x->set(2,1,n_Init(-4,coef), coef);
    x->set(3,1,n_Init(1,coef), coef);
    enumer = l->enumerate_next( c, x);
    if(enumer !=NULL){
        enumer->Print();
        PrintS("\n");
    }
    for(int i=1;i<4;i++){
        enumer = l->enumerate_next();
        if(enumer !=NULL){
            enumer->Print();
        }
        PrintS("\n");
    }
    PrintS("new bigintmat\n");
    x->set(3,1,n_Init(8,coef), coef);
    enumer = l->enumerate_next(x);
    if(enumer !=NULL){
        enumer->Print();
    }
    PrintS("\n");
    
        
    delete l;   
}


void test_Minkowski() {
    coeffs coef = nInitChar(n_Q,NULL);

    PrintS("Test Minkowski\n");
    PrintS("Elements\n");
    bigintmat ** elementarray = new bigintmat*[4];
    elementarray[0] = new bigintmat(4,1,coef);
    elementarray[1] = new bigintmat(4,1,coef);
    elementarray[2] = new bigintmat(4,1,coef);
    elementarray[3] = new bigintmat(4,1,coef);
    elementarray[0]->rawset(2,1,n_Init(1,coef),coef);
    
    elementarray[1]->rawset(1,1,n_Init(3,coef),coef);
    elementarray[1]->rawset(2,1,n_Init(7,coef),coef);
    elementarray[1]->rawset(3,1,n_Init(13,coef),coef);
    
    elementarray[2]->rawset(1,1,n_Init(-1,coef),coef);
    elementarray[2]->rawset(2,1,n_Init(42,coef),coef);
    elementarray[2]->rawset(3,1,n_Init(-5,coef),coef);
    elementarray[2]->rawset(4,1,n_Init(21,coef),coef);
    
    elementarray[3]->rawset(1,1,n_Init(2,coef),coef);
    elementarray[3]->rawset(2,1,n_Init(-1,coef),coef);
    elementarray[3]->rawset(3,1,n_Init(-5,coef),coef);
    elementarray[3]->rawset(4,1,n_Init(1,coef),coef);
    /// poly with real and imag roots
    PrintS("polynomial\n");
    number * poly = new number[5];//(number *)omAlloc( (5) * sizeof( number ) );//new number[5];
    //poly[0] = n_Init(6,coef);
    //poly[1] = n_Init(0,coef);
    //poly[2] = n_Init(5,coef);//positiv imaginär, negatic reelle wurzeln
    //poly[3] = n_Init(0,coef);
    //poly[4] = n_Init(1,coef);
    poly[0] = n_Init(-1,coef);
    poly[1] = n_Init(0,coef);
    poly[2] = n_Init(0,coef);
    poly[3] = n_Init(3,coef);
    poly[4] = n_Init(1,coef);
    int prec = 42;
    //coeffs rea = nInitChar(n_long_R,NULL);
    //setGMPFloatDigits( prec, prec);
    //number abc = n_Init(1,rea);
    //abc = n_Div(abc,n_Init(3333,rea),rea);
    //n_Print(abc,rea);
    //PrintS("\n");
    bigintmat * gitter = NULL;PrintS("Call function\n");
    gitter = minkowksi(elementarray,4,poly,4,coef,prec);
    if(gitter !=NULL){
        gitter->Print();
        PrintS("\n");
        //cout << "CoeffType of gitter: " << getCoeffType(gitter->basecoeffs()) << '\n';
    }
    //*/
}




void memorytest() {
    coeffs coef = nInitChar(n_Q,0);
    
    int ar[] = {1,-1,3,1,0,5,1,2,6};
    bigintmat * m = new bigintmat(3,3,coef);
    for(int i=0; i<9; i++) {
        m->set(i,n_Init(ar[i],coef), coef);
    }
    
    for(int i=0; i<100000000; i++) {
        lattice * l = new lattice(m);
        l->LLL();
//         bigintmat * reduced = l->get_reduced_basis();
        delete l;
    }
}

int average_length(bigintmat * m) {
    coeffs coef = nInitChar(n_Q,0);
    int total = 0;
    for(int j=1; j<=m->cols(); j++) {
        int col = 0;
        for(int i=1; i<=m->rows(); i++) {
            number m_ij = m->view(i,j);
            col += n_Int(m_ij,coef)*n_Int(m_ij,coef);
        }
        total += (int) sqrt ((double)col);
    }
    return total / m->cols();
}

double orthogonality_defect(bigintmat * m) {
    coeffs coef = m->basecoeffs();
    double total = 1;
    for(int j=1; j<=m->cols(); j++) {
        double col = 0;
        for(int i=1; i<=m->rows(); i++) {
            number m_ij = m->view(i,j);
            col += n_Int(m_ij,coef)*n_Int(m_ij,coef);
        }
        total *= sqrt(col);
    }
    number m_det = m->det();
    total /= ((double)n_Int(m_det,coef));
    return total>=0 ? total : -total;
}