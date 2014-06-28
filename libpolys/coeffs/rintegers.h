#ifndef RINTEGERS_H
#define RINTEGERS_H
/****************************************
*  Computer Algebra System SINGULAR     *
****************************************/
/*
* ABSTRACT: integers (as in ZZ)
*/

/*-----------------------------------------------------------------*/
/**
**  'SR_INT' is the type of those integers small enough to fit into  29  bits.
**  Therefor the value range of this small integers is: $-2^{28}...2^{28}-1$.
**
**  Small integers are represented by an immediate integer handle, containing
**  the value instead of pointing  to  it,  which  has  the  following  form:
**
**      +-------+-------+-------+-------+- - - -+-------+-------+-------+
**      | guard | sign  | bit   | bit   |       | bit   | tag   | tag   |
**      | bit   | bit   | 27    | 26    |       | 0     | 0     | 1     |
**      +-------+-------+-------+-------+- - - -+-------+-------+-------+
**
**  Immediate integers handles carry the tag 'SR_INT', i.e. the last bit is 1.
**  This distuingishes immediate integers from other handles which  point  to
**  structures aligned on 4 byte boundaries and therefor have last bit  zero.
**  (The second bit is reserved as tag to allow extensions of  this  scheme.)
**  Using immediates as pointers and dereferencing them gives address errors.
**
**  To aid overflow check the most significant two bits must always be equal,
**  that is to say that the sign bit of immediate integers has a  guard  bit.
**
**  The macros 'INT_TO_SR' and 'SR_TO_INT' should be used to convert  between
**  a small integer value and its representation as immediate integer handle.
**/

#define SR_HDL(A) ((long)(A))

#define SR_INT    1L
#define INT_TO_SR(INT)  ((number) (((long)INT << 2) + SR_INT))
#define SR_TO_INT(SR)   (((long)SR) >> 2)
#define n_Z_IS_SMALL(A)     (SR_HDL(A) & SR_INT)
#define INT_IS_SMALL(A) ( ((A << 1) >> 1) == A )

#define MP_SMALL 1


#ifdef HAVE_RINGS
#include <coeffs/coeffs.h>

BOOLEAN nrzInitChar    (coeffs r,  void * parameter);
number  nrzCopy        (number a, const coeffs r);
int     nrzSize        (number a, const coeffs r);
void    nrzDelete      (number *a, const coeffs r);
BOOLEAN nrzGreaterZero (number k, const coeffs r);
number  nrzMult        (number a, number b, const coeffs r);
number  nrzInit        (long i, const coeffs r);
int     nrzInt         (number &n, const coeffs r);
number  nrzAdd         (number a, number b, const coeffs r);
number  nrzSub         (number a, number b, const coeffs r);
void    nrzPower       (number a, int i, number * result, const coeffs r);
BOOLEAN nrzIsZero      (number a, const coeffs r);
BOOLEAN nrzIsOne       (number a, const coeffs r);
BOOLEAN nrzIsMOne      (number a, const coeffs r);
BOOLEAN nrzIsUnit      (number a, const coeffs r);
number  nrzGetUnit     (number a, const coeffs r);
number  nrzDiv         (number a, number b, const coeffs r);
number  nrzIntDiv      (number a, number b, const coeffs r);
number  nrzIntMod      (number a, number b, const coeffs r);
number  nrzNeg         (number c, const coeffs r);
number  nrzInvers      (number c, const coeffs r);
BOOLEAN nrzGreater     (number a, number b, const coeffs r);
BOOLEAN nrzDivBy       (number a, number b, const coeffs r);
int     nrzDivComp     (number a, number b, const coeffs r);
BOOLEAN nrzEqual       (number a, number b, const coeffs r);
number  nrzLcm         (number a,number b, const coeffs r);
number  nrzGcd         (number a,number b, const coeffs r);
number  nrzExtGcd      (number a, number b, number *s, number *t, const coeffs r);
number  nrzQuotRem     (number a, number b, number *q, const coeffs r);
nMapFunc nrzSetMap     (const coeffs src, const coeffs dst);
void    nrzWrite       (number &a, const coeffs r);
const char *  nrzRead  (const char *s, number *a, const coeffs r);
char *  nrzName        (number n, const coeffs r);
void    nrzCoeffWrite  (const coeffs r, BOOLEAN details);
#ifdef LDEBUG
BOOLEAN nrzDBTest      (number a, const char *f, const int l, const coeffs r);
#endif
void    nrzSetExp(int c, coeffs r);
void    nrzInitExp(int c, coeffs r);
void    nrzDelete(number *a, const coeffs r);
coeffs  nrzQuot1(number c, const coeffs r);

number nrzMapQ(number from, const coeffs src, const coeffs dst);
#endif
#endif
