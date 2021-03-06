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
    
//     PrintS("\n\n\n");
//     test_Enumerate();
//     PrintS("\n\n\n");
//     test_Enumerate(currRing->cf);
//     PrintS("\n\n\n");
//     test_LLL();
        
    test_Poly();
//     coeffs coef = nInitChar(n_R, NULL);
//     number t = n_Div(n_Init(110,coef),n_Init(100,coef),coef);
//     n_Print(t,coef);
//     PrintS("\n");
//     number r = round(t,coef);
//     PrintS("\n\n\n");
//     n_Print(r,coef);
//     PrintS("\n\n\n");
    
//     test_Minkowski();
    return c;    
}


void test_LLL() {
//     coeffs coef = nInitChar(n_Q,NULL);
//     int precision = 6;
//     LongComplexInfo paramReal;
//     paramReal.float_len = si_min (precision, 32767);
//     paramReal.float_len2 = si_min (precision, 32767);
//     paramReal.par_name=(const char*)"i";
//     coeffs coef = nInitChar(n_Z, NULL);
//     coeffs coef = nInitChar(n_long_R, &paramReal);
    coeffs coef = nInitChar(n_R, NULL);
    
    PrintS("Test LLL\n");    
    
    
    //Vector containing test matrices
    std::vector<bigintmat*> vec;
    
//     {
//         int ar[] = {1,-1,3,1,0,5,1,2,6};
//         int n = 3;
//         int m = 3;
//         bigintmat * temp = new bigintmat(n,m,coef);
//         for(int i=0; i<(n*m); i++) {
//             temp->set(i,n_Init(ar[i],coef), coef);
//         }
//         vec.push_back(temp);
//     }
    
//     {
//         int ar[] = {1,9,1,2,1,8,8,3,7,4,5,1,2,6,7,1};
//         int n = 4;
//         int m = 4;
//         bigintmat * temp = new bigintmat(n,m,coef);
//         for(int i=0; i<(n*m); i++) {
//             temp->set(i,n_Init(ar[i],coef), coef);
//         }
//         vec.push_back(temp);
//     }
    
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
        
//     {
//         int ar[] = {1,-1,3,1,0,5,1,2,6};
//         int n = 3;
//         int m = 3;
//         bigintmat * temp = new bigintmat(n,m,coef);
//         for(int i=0; i<(n*m); i++) {
//             temp->set(i,n_Init(ar[i],coef), coef);
//         }
//         vec.push_back(temp);
//     }
    
//     {
//         int ar[] = {1,2,3,31,41,51,101,201,301};
//         int n = 3;
//         int m = 3;
//         bigintmat * temp = new bigintmat(n,m,coef);
//         for(int i=0; i<(n*m); i++) {
//             temp->set(i,n_Init(ar[i],coef), coef);
//         }
//         vec.push_back(temp);
//     }
    
//     {
//         int ar[] = {10000,4545545,4557454,465445241};
//         int n = 2;
//         int m = 2;
//         bigintmat * temp = new bigintmat(n,m,coef);
//         for(int i=0; i<(n*m); i++) {
//             temp->set(i,n_Init(ar[i],coef), coef);
//         }
//         vec.push_back(temp);
//     }

    //Random matrices
    srand(1234);
    for(int t = 0; t < 1; t++) {
        int matrix_size = rand() % 10 + 1;
        bigintmat * matrix = new bigintmat(matrix_size,matrix_size,coef);
        matrix->one();
        for(int ops = 0; ops < matrix_size*5; ops++) {
            matrix->addcol(rand() % matrix_size + 1
                          ,rand() % matrix_size + 1
                          ,n_Init(rand() % 5 + 1,coef),coef);
        }
        vec.push_back(matrix);
    }
    
    
    for(unsigned i= 0; i<vec.size(); i++) {
        cout << "\nMatrix #" << i;
        PrintS("\nInputmatrix:\n");
        vec[i]->Print();
        PrintS("\n");
        
        lattice * l = new lattice(vec[i]);
        number c = NULL;
        l->LLL(c,coef,true,true,true);
        bigintmat * reduced = l->get_reduced_basis();
        bigintmat * H = l->get_transformation_matrix();
        
        PrintS("Output:\n");
        reduced->Print();
        PrintS("\n");
        
        bigintmat * input_times_H = bimMult(vec[i],H);
        bigintmat * errormatrix = bimSub(reduced,input_times_H);
        bool transformation_is_correct = errormatrix->isZero();
        
        assume(transformation_is_correct);
        if(transformation_is_correct) {
            PrintS("transformation is correct\n");
        } else {
            PrintS("transformation is NOT correct:\n");
            PrintS("H =\n");
            H->Print();
            PrintS("\nerrormatrix =\n");
            errormatrix->Print();
            PrintS("\n");
            getchar();
        }
        
//         number H_det = H->det();
//         bool H_is_regular = n_IsOne(H_det,coef) || n_IsOne(n_Neg(H_det,coef),coef);
//         
//         assume(H_is_regular);
//         if(H_is_regular) {
//             PrintS("H is regular\n");
//         } else {
//             PrintS("H is not regular:\n");
//             PrintS("H =\n");
//             H->Print();
//             PrintS("\n");
//             getchar();
//         }

//         bigintmat * gram_matrix = l->get_gram_matrix();
//         lattice * l_gram = new lattice(vec[i],true);
//         l_gram->LLL(c,true,false,true);
//         bigintmat * H_gram = l_gram->get_transformation_matrix();
//         PrintS("H_gram :\n");
//         H_gram->Print();
//         PrintS("\n");
// //         assume(bimSub(H_gram,H)->isZero());

        delete reduced;
        delete H;
        delete input_times_H;
        delete errormatrix;
        delete l;
    }
    
    
    
    for(unsigned i= 0; i<vec.size(); i++) {
        delete vec[i];
    }
    
    PrintS("\n");
}


void test_Enumerate(coeffs coef) {   
    PrintS("Test Enumerate\n");
    
    //Random matrices
    srand(1234);
    
    number c = n_Init(1000,coef);
    
    for(int t = 0; t < 1; t++) {
        int matrix_size = rand() % 50 + 1;
        if(matrix_size<2 ) matrix_size=2;
        bigintmat * matrix = new bigintmat(matrix_size,matrix_size,coef);
        for(int col =1; col <=matrix_size;col++){
            for(int row =1; row <=matrix_size;row++){
                number nom = n_Init(rand() % 4096 -2047,coef);
                number den = n_Init(rand() % 1024 + 1,coef);
                matrix->rawset(row,col,n_Div(nom,den,coef),coef);
                n_Delete(&nom,coef);
                n_Delete(&den,coef);
            }
        }
        
        cout << "\ni = " << t;
        PrintS("\nInput:\n");
        PrintS("\n");
        
        lattice * l = new lattice(matrix);
        bigintmat * x = new bigintmat(matrix->cols(),1,coef);
        x->rawset(matrix->cols()-1,1,n_Init(-1,coef),coef);
        bigintmat * elem = l->enumerate_next(c,x);
        
        PrintS("Output:\n");
        if(elem!=NULL){
            //elem->Print();PrintS("\n");
        } else {
            PrintS("NULL\n");
        }
        delete x;
        delete elem;
        
        PrintS("del lattice\n");
        delete l;
        delete matrix;
    }//*/
    n_Delete(&c,coef);
    PrintS("\n");
}


void test_Minkowski() {
    coeffs coef = nInitChar(n_Q,NULL);

    PrintS("Test Minkowski\n");
//     PrintS("Elements\n");
//     bigintmat ** elementarray = new bigintmat*[4];
//     elementarray[0] = new bigintmat(4,1,coef);
//     elementarray[1] = new bigintmat(4,1,coef);
//     elementarray[2] = new bigintmat(4,1,coef);
//     elementarray[3] = new bigintmat(4,1,coef);
//     elementarray[0]->rawset(2,1,n_Init(1,coef),coef);
//     
//     elementarray[1]->rawset(1,1,n_Init(3,coef),coef);
//     elementarray[1]->rawset(2,1,n_Init(7,coef),coef);
//     elementarray[1]->rawset(3,1,n_Init(13,coef),coef);
//     
//     elementarray[2]->rawset(1,1,n_Init(-1,coef),coef);
//     elementarray[2]->rawset(2,1,n_Init(42,coef),coef);
//     elementarray[2]->rawset(3,1,n_Init(-5,coef),coef);
//     elementarray[2]->rawset(4,1,n_Init(21,coef),coef);
//     
//     elementarray[3]->rawset(1,1,n_Init(2,coef),coef);
//     elementarray[3]->rawset(2,1,n_Init(-1,coef),coef);
//     elementarray[3]->rawset(3,1,n_Init(-5,coef),coef);
//     elementarray[3]->rawset(4,1,n_Init(1,coef),coef);
    
    
    int ar[] = {1,-1,3,1,0,5,1,2,6};
    int n = 3;
    int m = 3;
    bigintmat * mat = new bigintmat(n,m,coef);
    for(int i=0; i<(n*m); i++) {
        mat->set(i,n_Init(ar[i],coef), coef);
    }
    
    
    /// poly with real and imag roots
    PrintS("polynomial\n");
    number * poly = new number[4];//(number *)omAlloc( (5) * sizeof( number ) );//new number[5];
    //poly[0] = n_Init(6,coef);
    //poly[1] = n_Init(0,coef);
    //poly[2] = n_Init(5,coef);//positiv imaginär, negatic reelle wurzeln
    //poly[3] = n_Init(0,coef);
    //poly[4] = n_Init(1,coef);
    poly[0] = n_Init(-1,coef);
    poly[1] = n_Init(0,coef);
    poly[2] = n_Init(0,coef);
    poly[3] = n_Init(3,coef);
//     poly[4] = n_Init(1,coef);
    int prec = 42;
    //coeffs rea = nInitChar(n_long_R,NULL);
    //setGMPFloatDigits( prec, prec);
    //number abc = n_Init(1,rea);
    //abc = n_Div(abc,n_Init(3333,rea),rea);
    //n_Print(abc,rea);
    //PrintS("\n");
    lattice * gitter = (lattice *)omAlloc(sizeof(lattice));
    minkowski(mat,poly,3,coef,prec,gitter);

    if(gitter != NULL){
        bigintmat * b = gitter->get_basis();
        b->Print();
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



void test_Poly() {
    ring r = currRing;
    coeffs coef = r->cf;
    int deg = 3;
    
    number * univpol = new number[3+1]; 
//     for(int j=0; j<=deg; j++) {
//         univpol[j] = n_Init(1,coef);
//     }
    univpol[0] = n_Init(169,coef);
    univpol[1] = n_Init(27,coef);
    univpol[2] = n_Init(-13,coef);
    univpol[3] = n_Init(1,coef);
    
    poly f = numbers2poly(univpol, deg, coef, r);
    
//     number * pcoeffs = poly2numbers(f,r, coef);
//     int degf = (int) p_Totaldegree(f,r);
//     for(int i=0; i<=degf; i++) {
//         n_Print(pcoeffs[i],coef);
//         PrintS("\n");
//     }
    
    poly g = get_nice_poly(f);
    
    number * pcoeffs = NULL;
    int degg = poly2numbers(g,pcoeffs,r,coef);
    for(int i=0; i<=degg; i++) {
        n_Print(pcoeffs[i],coef);
        PrintS("\n");
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
