#include <libpolys/coeffs/bigintmat.h>
#include "temptest.h"
#include <reporter/reporter.h>
#include "libpolys/coeffs/numbers.h" 
#include "libpolys/coeffs/coeffs.h"
#include "Singular/ipid.h"
#include "kernel/febase.h"
#include "lattice.h"
#include <iostream>



using namespace std;

//Temporary testfunction to play arround with new functions
//NOTE: remove later
bigintmat* temp_test(bigintmat& a) {
    PrintS("This is a Test\n");
    bigintmat* bim = new bigintmat(&a);
    bim->Print();
    
    cout << "TEST " << getCoeffType(bim->basecoeffs()) << '\n';
//     lattice* l = new lattice(bim);
//     l->LLL();
//     delete l;
//     l->MLLL();
//     bigintmat n = l->get_lattice();
    bigintmat* c = bimAdd(bim,bim);
    delete bim;
    return c;    
}

number temp_test2(number a) {    
    number c = nCopy(a);
    coeffs coef = currRing->cf;
    
    cout << "CoeffType of currRing: " << getCoeffType(currRing->cf) << '\n';
    
//     bigintmat *m = new bigintmat(4,4,coef);
//     for(int i=0; i<=15; i++) {
//         m->set(i,n_Init(i,coef), coef);
//     }
    
    int ar[] = {1,-1,3,1,0,5,1,2,6};
    bigintmat *m = new bigintmat(3,3,coef);
    for(int i=0; i<9; i++) {
        m->set(i,n_Init(ar[i],coef), coef);
    }

//     int ar[] = {10000,4545545,4557454,465445241};
//     bigintmat *m = new bigintmat(2,2,coef);
//     for(int i=0; i<4; i++) {
//         m->set(i,n_Init(ar[i],coef), coef);
//     }
    
    PrintS("Input: ");
    m->Print();
    PrintS("\n");
    
    lattice* l = new lattice(m);
    l->LLL();
    bigintmat *reduced = l->get_reduced_basis();
    
    PrintS("Output: ");
    reduced->Print();
    PrintS("\n");
    
    
    delete l;
    return c;    
}

