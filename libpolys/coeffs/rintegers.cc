/****************************************
*  Computer Algebra System SINGULAR     *
****************************************/
/*
* ABSTRACT: numbers (integers)
*/
#ifdef HAVE_CONFIG_H
#include "libpolysconfig.h"
#endif /* HAVE_CONFIG_H */
#include <misc/auxiliary.h>
#include <factory/factory.h>

#ifdef HAVE_RINGS

#include <string.h>
#include <misc/mylimits.h>
#include <coeffs/coeffs.h>
#include <reporter/reporter.h>
#include <omalloc/omalloc.h>
#include <coeffs/numbers.h>
#include <coeffs/longrat.h>
#include <coeffs/mpr_complex.h>
#include <coeffs/rintegers.h>
#include "si_gmp.h"

/// Our Type!
static const n_coeffType ID = n_Z;

omBin gmp_nrz_bin = omGetSpecBin(sizeof(mpz_t));


//make sure that a small number is an immediate integer
//bascially coped from longrat.cc nlShort3
//TODO: is there any point in checking 0 first???
//TODO: it is not clear that this works in 32/64 bit everywhere.
//      too many hacks.

#define CF_DEBUG 0
static inline number nrz_short(number x) 
{
#if CF_DEBUG
  StringAppendS("short(");
  nrzWrite(x, NULL);
#endif
  if (mpz_cmp_ui((int_number) x,(long)0)==0)
  {
    mpz_clear((int_number)x);
    omFreeBin(x, gmp_nrz_bin);
#if CF_DEBUG
    StringAppendS(")=0");
#endif
    return INT_TO_SR(0);
  }
  if (mpz_size1((int_number)x)<=MP_SMALL)
  {
    int ui=mpz_get_si((int_number)x);
    if ((((ui<<3)>>3)==ui)
    && (mpz_cmp_si((int_number)x,(long)ui)==0))
    {
      mpz_clear((int_number)x);
      omFreeBin(x, gmp_nrz_bin);
#if CF_DEBUG
    StringAppendS(")=imm");
#endif
      return INT_TO_SR(ui);
    }
  }
#if CF_DEBUG
  StringAppendS(")");
#endif
  return x;
}


/*
 * Multiply two numbers
 * check for 0, 1, -1 maybe
 */
#if CF_DEBUG
number _nrzMult(number, number, const coeffs);
number nrzMult(number a, number b, const coeffs R)
{
  StringSetS("Mult: ");
  nrzWrite(a, R);
  StringAppendS(" by ");
  nrzWrite(b, R);
  number c = _nrzMult(a, b, R);
  StringAppendS(" is ");
  nrzWrite(c, R);
  char * s = StringEndS();
  Print("%s\n", s);
  omFree(s);
  return c;
}
number _nrzMult (number a, number b, const coeffs R)
#else
number nrzMult (number a, number b, const coeffs R)
#endif
{
  if (IS_SMALL(a) && IS_SMALL(b)) {
  //from longrat.cc 
    if (SR_TO_INT(b)==0)
      return b;
    long r=(long)((unsigned long)(SR_HDL(a)-1L))*((unsigned long)(SR_HDL(b)>>1));
    if ((r/(SR_HDL(b)>>1))==(SR_HDL(a)-1L))
    {
      number u=((number) ((r>>1)+SR_INT));
    //  if (((((long)SR_HDL(u))<<1)>>1)==SR_HDL(u)) return (u);
      return nrzInit(SR_HDL(u)>>2, R);
    }
    int_number erg = (int_number) omAllocBin(gmp_nrz_bin);
    mpz_init(erg);
    mpz_set_si(erg, SR_TO_INT(a));
    mpz_mul_si(erg, erg, SR_TO_INT(b));
    return (number) erg;
  } else if (IS_SMALL(a)) {
    int_number erg = (int_number) omAllocBin(gmp_nrz_bin);
    mpz_init_set(erg, (int_number) b);
    mpz_mul_si(erg, erg, SR_TO_INT(a));
    return (number) erg;
  } else if (IS_SMALL(b)) {
    int_number erg = (int_number) omAllocBin(gmp_nrz_bin);
    mpz_init_set(erg, (int_number) a);
    mpz_mul_si(erg, erg, SR_TO_INT(b));
    return (number) erg;
  } else {
    int_number erg = (int_number) omAllocBin(gmp_nrz_bin);
    mpz_init(erg);
    mpz_mul(erg, (int_number) a, (int_number) b);
    return (number) erg;
  }
}


static int int_gcd(int a, int b) 
{
  int r;
  a = ABS(a);
  b = ABS(b);
  if (!a) return b;
  if (!b) return a;
  do {
    r = a % b;
    a = b;
    b = r;
  } while (b);
  return ABS(a); // % in c doeas not imply a signn
                 // it woould be unlikely to see a negative here
                 // but who knows
}

/*
 * Give the smallest non unit k, such that a * x = k = b * y has a solution
 */
number nrzLcm (number a, number b, const coeffs R)
{
  int_number erg;
  if (IS_SMALL(a) && IS_SMALL(b)) {
    int g = int_gcd(SR_TO_INT(a), SR_TO_INT(b));
    return nrzMult(a, INT_TO_SR(SR_TO_INT(b)/g), R);
  } else if (IS_SMALL(a)) {
    erg = (int_number) omAllocBin(gmp_nrz_bin);
    mpz_init_set(erg, (int_number) b);
    unsigned long g = mpz_gcd_ui(NULL, erg, (unsigned long) ABS(SR_TO_INT(a)));
    mpz_mul_si(erg, erg, SR_TO_INT(a)/g);
  } else if (IS_SMALL(b)) {
    erg = (int_number) omAllocBin(gmp_nrz_bin);
    mpz_init_set(erg, (int_number) a);
    unsigned long g = mpz_gcd_ui(NULL, erg, (unsigned long) ABS(SR_TO_INT(b)));
    mpz_mul_si(erg, erg, SR_TO_INT(a)/g);
  } else {
    erg = (int_number) omAllocBin(gmp_nrz_bin);
    mpz_init(erg);
    mpz_lcm(erg, (int_number) a, (int_number) b);
  }
  return (number) erg;
}

/*
 * Give the largest non unit k, such that a = x * k, b = y * k has
 * a solution.
 */
number nrzGcd (number a,number b,const coeffs R)
{
  if (IS_SMALL(a) && IS_SMALL(b)) {
    int g = int_gcd(SR_TO_INT(a), SR_TO_INT(b));
    return INT_TO_SR(g);
  } else if (IS_SMALL(a)) { 
    if (a==INT_TO_SR(0))
      return nrzCopy(b, R);
    unsigned long g = mpz_gcd_ui(NULL, (int_number)b, (unsigned long) ABS(SR_TO_INT(a)));
    return INT_TO_SR((int) g);
  } else if (IS_SMALL(b)) {
    if (b==INT_TO_SR(0))
      return nrzCopy(a, R);
    unsigned long g = mpz_gcd_ui(NULL, (int_number)a, (unsigned long) ABS(SR_TO_INT(b)));
    return INT_TO_SR((int) g);
  } else {
    int_number erg = (int_number) omAllocBin(gmp_nrz_bin);
    mpz_init(erg);
    mpz_gcd(erg, (int_number) a, (int_number) b);
    return (number) erg;
  }
}

/*
 * Give the largest non unit k, such that a = x * k, b = y * k has
 * a solution and r, s, s.t. k = s*a + t*b
 */
static int int_extgcd(int a, int b, int * u, int* x, int * v, int* y)
{
  int q, r;
  if (!a) {
    *u = 0;
    *v = 1;
    *x = -1;
    *y = 0;
    return b;
  }
  if (!b) {
    *u = 1;
    *v = 0;
    *x = 0;
    *y = 1;
    return a;
  }
  *u=1;
  *v=0;
  *x=0;
  *y=1;
  do {
    q = a/b;
    r = a%b;
    assume (q*b+r == a);
    a = b;
    b = r;

    r = -(*v)*q+(*u);
    (*u) =(*v);
    (*v) = r;

    r = -(*y)*q+(*x);
    (*x) = (*y);
    (*y) = r;
  } while (b);

  return a;
}

number  nrzExtGcd (number a, number b, number *s, number *t, const coeffs)
{
  if (IS_SMALL(a) && IS_SMALL(b)) {
    int u, v, x, y;
    int g = int_extgcd(SR_TO_INT(a), SR_TO_INT(b), &u, &v, &x, &y);
    *s = INT_TO_SR(u);
    *t = INT_TO_SR(v);
    return INT_TO_SR(g);
  } else {
    mpz_t aa, bb;
    if (IS_SMALL(a)) {
      mpz_init_set_si(aa, SR_TO_INT(a));
    } else {
      mpz_init_set(aa, (int_number) a);
    }
    if (IS_SMALL(b)) {
      mpz_init_set_si(bb, SR_TO_INT(b));
    } else {
      mpz_init_set(bb, (int_number) b);
    }
    int_number erg = (int_number) omAllocBin(gmp_nrz_bin);
    int_number bs = (int_number) omAllocBin(gmp_nrz_bin);
    int_number bt = (int_number) omAllocBin(gmp_nrz_bin);
    mpz_init(erg);
    mpz_init(bs);
    mpz_init(bt);
    mpz_gcdext(erg, bs, bt, (int_number) a, (int_number) b);
    *s = nrz_short((number) bs);
    *t = nrz_short((number) bt);
    return nrz_short((number) erg);
  }
}
#if CF_DEBUG
number _nrzXExtGcd(number, number, number *, number *, number *, number *, const coeffs);
number nrzXExtGcd(number a, number b, number *x, number * y, number * u, number * v, const coeffs R)
{
  char * s;
  StringSetS("XExtGcd: ");
  nrzWrite(a, R);
  StringAppendS(" by ");
  nrzWrite(b, R);
  number c = _nrzXExtGcd(a, b, x, y, u, v, R);
  StringAppendS(" is ");
  nrzWrite(c, R);
  StringAppendS("[[");
  nrzWrite(*x, R);
  StringAppendS(", ");
  nrzWrite(*y, R);
  StringAppendS("], ");
  nrzWrite(*u, R);
  StringAppendS(", ");
  nrzWrite(*v, R);
  s=StringEndS();
  Print("%s]]\n", s);
  omFree(s);
  return c;
}
number  _nrzXExtGcd (number a, number b, number *s, number *t, number *u, number *v, const coeffs )
#else
number  nrzXExtGcd (number a, number b, number *s, number *t, number *u, number *v, const coeffs )
#endif
{
  if (IS_SMALL(a) && IS_SMALL(b)) {
    int uu, vv, x, y;
    int g = int_extgcd(SR_TO_INT(a), SR_TO_INT(b), &uu, &vv, &x, &y);
    *s = INT_TO_SR(uu);
    *t = INT_TO_SR(vv);
    *u = INT_TO_SR(x);
    *v = INT_TO_SR(y);
    return INT_TO_SR(g);
  } else {
    mpz_t aa, bb;
    if (IS_SMALL(a)) {
      mpz_init_set_si(aa, SR_TO_INT(a));
    } else {
      mpz_init_set(aa, (int_number) a);
    }
    if (IS_SMALL(b)) {
      mpz_init_set_si(bb, SR_TO_INT(b));
    } else {
      mpz_init_set(bb, (int_number) b);
    }
    int_number erg = (int_number) omAllocBin(gmp_nrz_bin);
    int_number bs = (int_number) omAllocBin(gmp_nrz_bin);
    int_number bt = (int_number) omAllocBin(gmp_nrz_bin);
    mpz_init(erg);
    mpz_init(bs);
    mpz_init(bt);

    mpz_gcdext(erg, bs, bt, aa, bb);

    int_number bu = (int_number) omAllocBin(gmp_nrz_bin);
    int_number bv = (int_number) omAllocBin(gmp_nrz_bin);

    mpz_init_set(bu, (int_number) bb);
    mpz_init_set(bv, (int_number) aa);

    mpz_clear(aa);
    mpz_clear(bb);
    assume(mpz_cmp_si(erg, 0));

    mpz_div(bu, bu, erg);
    mpz_div(bv, bv, erg);

    mpz_mul_si(bu, bu, -1);
    *u = nrz_short((number) bu);
    *v = nrz_short((number) bv);

    *s = nrz_short((number) bs);
    *t = nrz_short((number) bt);
    return nrz_short((number) erg);
  }
}
#if CF_DEBUG
number _nrzQuotRem(number, number, number *, const coeffs);
number nrzQuotRem(number a, number b, number * r, const coeffs R)
{
  StringSetS("QuotRem: ");
  nrzWrite(a, R);
  StringAppendS(" by ");
  nrzWrite(b, R);
  number c = _nrzQuotRem(a, b, r, R);
  StringAppendS(" is ");
  nrzWrite(c, R);
  if (r) {
    StringAppendS("+R(");
    nrzWrite(*r, R);
    StringAppendS(")");
  }
  char * s = StringEndS();
  Print("%s\n", s);
  omFree(s);
  return c;
}
number _nrzQuotRem (number a, number b, number * r, const coeffs )
#else
number nrzQuotRem (number a, number b, number * r, const coeffs )
#endif
{
  assume(SR_TO_INT(b));
  if (IS_SMALL(a) && IS_SMALL(b)) {
    if (r)
      *r = INT_TO_SR(SR_TO_INT(a) % SR_TO_INT(b));
    return INT_TO_SR(SR_TO_INT(a)/SR_TO_INT(b));
  } else if (IS_SMALL(a)) {
    //a is small, b is not, so q=0, r=a
    if (r)
      *r = a;
    return INT_TO_SR(0);
  } else if (IS_SMALL(b)) {
    unsigned long rr;
    int_number qq = (int_number) omAllocBin(gmp_nrz_bin);
    mpz_init(qq);
    mpz_t rrr;
    mpz_init(rrr);
    rr = mpz_divmod_ui(qq, rrr, (int_number) a, (unsigned long)ABS(SR_TO_INT(b)));
    mpz_clear(rrr);

    if (r)
      *r = INT_TO_SR(rr);
    if (SR_TO_INT(b)<0) {
      mpz_mul_si(qq, qq, -1);
    }
    return nrz_short((number)qq);
  }
  int_number qq = (int_number) omAllocBin(gmp_nrz_bin),
             rr = (int_number) omAllocBin(gmp_nrz_bin);
  mpz_init(qq);
  mpz_init(rr);
  mpz_divmod(qq, rr, (int_number)a, (int_number)b);
  if (r) 
    *r = (number) rr;
  else {
    mpz_clear(rr);
  }
  return (number) qq;
}


void nrzPower (number a, int i, number * result, const coeffs)
{
  int_number erg = (int_number) omAllocBin(gmp_nrz_bin);
  mpz_init(erg);
  mpz_t aa;
  if (IS_SMALL(a))
    mpz_init_set_si(aa, SR_TO_INT(a));
  else
    mpz_init_set(aa, (int_number) a);
  mpz_pow_ui(erg, aa, i);
  *result = nrz_short((number) erg);
}

/*
 * create a number from int
 * TODO: do not create an mpz if not necessary
 */
number nrzInit (long i, const coeffs)
{
  int_number erg = (int_number) omAllocBin(gmp_nrz_bin);
  mpz_init_set_si(erg, i);
  return nrz_short((number) erg);
}

void nrzDelete(number *a, const coeffs)
{
  if (*a == NULL) return;
  if (IS_SMALL(*a)) {
    *a = NULL;
    return;
  }
  mpz_clear((int_number) *a);
  omFreeBin((ADDRESS) *a, gmp_nrz_bin);
  *a = NULL;
}

number nrzCopy(number a, const coeffs)
{
  if (IS_SMALL(a)) return a;
  int_number erg = (int_number) omAllocBin(gmp_nrz_bin);
  mpz_init_set(erg, (int_number) a);
  return (number) erg;
}

int nrzSize(number a, const coeffs)
{
  if (a == NULL) return 0;
  return sizeof(mpz_t);
}

/*
 * convert a number to int
 */
int nrzInt(number &n, const coeffs)
{
  if (IS_SMALL(n)) return SR_TO_INT(n);
  return (int) mpz_get_si( (int_number)n);
}
#if CF_DEBUG
number _nrzAdd(number, number, const coeffs);
number nrzAdd(number a, number b, const coeffs R)
{
  StringSetS("Add: ");
  nrzWrite(a, R);
  StringAppendS(" to ");
  nrzWrite(b, R);
  number c = _nrzAdd(a, b, R);
  StringAppendS(" is ");
  nrzWrite(c, R);
  char * s = StringEndS();
  Print("%s\n", s);
  omFree(s);
  return c;
}
number _nrzAdd (number a, number b, const coeffs R)
#else
number nrzAdd (number a, number b, const coeffs R)
#endif
{
  if (IS_SMALL(a) && IS_SMALL(b)) {
    int c = SR_TO_INT(a) + SR_TO_INT(b);
    if (INT_IS_SMALL(c))
      return INT_TO_SR(c);
    int_number erg = (int_number) omAllocBin(gmp_nrz_bin);
    mpz_init_set_si(erg, c);
    return (number) erg;
  } else if (IS_SMALL(a)) { 
    int_number erg = (int_number) omAllocBin(gmp_nrz_bin);
    mpz_init(erg);
    if (SR_TO_INT(a)>0) 
      mpz_add_ui(erg, (int_number) b, (unsigned long)SR_TO_INT(a));
    else
      mpz_sub_ui(erg, (int_number) b, (unsigned long)-(SR_TO_INT(a)));
    return nrz_short((number) erg);
  } else if (IS_SMALL(b)) {
    int_number erg = (int_number) omAllocBin(gmp_nrz_bin);
    mpz_init(erg);
    if (SR_TO_INT(b)>0) 
      mpz_add_ui(erg, (int_number) a, (unsigned long)SR_TO_INT(b));
    else
      mpz_sub_ui(erg, (int_number) a, (unsigned long)-(SR_TO_INT(b)));
    return nrz_short((number) erg);
  } else {
    int_number erg = (int_number) omAllocBin(gmp_nrz_bin);
    mpz_init(erg);
    mpz_add(erg, (int_number) a, (int_number) b);
    return nrz_short((number) erg);
  }
}

number nrzSub (number a, number b, const coeffs)
{
  if (IS_SMALL(a) && IS_SMALL(b)) {
    int c = SR_TO_INT(a) - SR_TO_INT(b);
    if (INT_IS_SMALL(c))
      return INT_TO_SR(c);
    int_number erg = (int_number) omAllocBin(gmp_nrz_bin);
    mpz_init_set_si(erg, c);
    return (number) erg;
  } else if (IS_SMALL(a)) { 
    int_number erg = (int_number) omAllocBin(gmp_nrz_bin);
    mpz_init(erg);
   
    if (SR_TO_INT(a)>0)
      mpz_ui_sub(erg, (unsigned long)SR_TO_INT(a), (int_number) b);
    else {
      mpz_add_ui(erg, (int_number) b, (unsigned long)-SR_TO_INT(a));
      mpz_neg(erg, erg);
    }
    return nrz_short((number) erg);
  } else if (IS_SMALL(b)) {
    int_number erg = (int_number) omAllocBin(gmp_nrz_bin);
    mpz_init(erg);
    if (INT_TO_SR(b)>0)
      mpz_sub_ui(erg, (int_number) a, (unsigned long)SR_TO_INT(b));
    else
      mpz_add_ui(erg, (int_number) a, (unsigned long)-SR_TO_INT(b));
    return nrz_short((number) erg);
  } else {
    int_number erg = (int_number) omAllocBin(gmp_nrz_bin);
    mpz_init(erg);
    mpz_sub(erg, (int_number) a, (int_number) b);
    return nrz_short((number) erg);
  }

}

number  nrzGetUnit (number n, const coeffs r)
{
  if (IS_SMALL(n)) {
    if (SR_TO_INT(n)<0)
      return INT_TO_SR(-1);
    else
      return INT_TO_SR(1);
  }
  if (nrzGreaterZero(n, r))
    return INT_TO_SR(1);
  else
    return INT_TO_SR(-1);
}

BOOLEAN nrzIsUnit (number a, const coeffs)
{
  return IS_SMALL(a) && ABS(SR_TO_INT(a))==1;
}

BOOLEAN nrzIsZero (number  a, const coeffs)
{
  return IS_SMALL(a) && a==INT_TO_SR(0);
}

BOOLEAN nrzIsOne (number a, const coeffs)
{
  return IS_SMALL(a) && a==INT_TO_SR(1);
}

BOOLEAN nrzIsMOne (number a, const coeffs)
{
  return IS_SMALL(a) && a==INT_TO_SR(-1);
}

BOOLEAN nrzEqual (number a,number b, const coeffs)
{
  if (IS_SMALL(a) && IS_SMALL(b))
    return a==b;
  else if (IS_SMALL(a) || IS_SMALL(b))
    return FALSE;
  else
    return 0 == mpz_cmp((int_number) a, (int_number) b);
}

BOOLEAN nrzGreater (number a,number b, const coeffs)
{
  if (IS_SMALL(a) && IS_SMALL(b))
    return a>b;
  else if (IS_SMALL(a))
    return FALSE;
  else if (IS_SMALL(b))
    return TRUE;
  return 0 < mpz_cmp((int_number) a, (int_number) b);
}

BOOLEAN nrzGreaterZero (number k, const coeffs C)
{
  return nrzGreater(k, INT_TO_SR(0), C);
}

int nrzDivComp(number a, number b, const coeffs r)
{
  if (nrzDivBy(a, b, r))
  {
    if (nrzDivBy(b, a, r)) return 2;
    return -1;
  }
  if (nrzDivBy(b, a, r)) return 1;
  return 0;
}

BOOLEAN nrzDivBy (number a,number b, const coeffs)
{
  if (IS_SMALL(a) && IS_SMALL(b)) {
    return SR_TO_INT(a) %SR_TO_INT(b) ==0;
  } else if (IS_SMALL(a)) {
    return a==INT_TO_SR(0);
  } else if (IS_SMALL(b)) {
    return mpz_divisible_ui_p((int_number)a, (unsigned long)ABS(SR_TO_INT(b))) != 0;
  } else     
    return mpz_divisible_p((int_number) a, (int_number) b) != 0;
}

number nrzDiv (number a,number b, const coeffs R)
{
  assume(SR_TO_INT(b));
  if (IS_SMALL(a) && IS_SMALL(b)) {
    if (SR_TO_INT(a) % SR_TO_INT(b)) {
      WerrorS("1:Division by non divisible element.");
      WerrorS("Result is without remainder.");
    }
    return INT_TO_SR(SR_TO_INT(a)/SR_TO_INT(b));
  } else if (IS_SMALL(a)) {
    if (SR_TO_INT(a)) {
      WerrorS("2:Division by non divisible element.");
      WerrorS("Result is without remainder.");
    }
    return INT_TO_SR(0);
  } else if (IS_SMALL(b)) {
    int_number erg = (int_number) omAllocBin(gmp_nrz_bin);
    mpz_t r;
    mpz_init(r);
    mpz_init(erg);
    if (mpz_divmod_ui(erg, r, (int_number) a, (unsigned long)ABS(SR_TO_INT(b)))) {
      WerrorS("3:Division by non divisible element.");
      WerrorS("Result is without remainder.");
    }
    mpz_clear(r);
    if (SR_TO_INT(b)<0)
      mpz_neg(erg, erg);
    return nrz_short((number) erg);
  }
  int_number erg = (int_number) omAllocBin(gmp_nrz_bin);
  mpz_init(erg);
  mpz_t r;
  mpz_init(r);
  mpz_tdiv_qr(erg, r, (int_number) a, (int_number) b);
#if CF_DEBUG
  StringSetS("division of");
  nrzWrite(a, R);
  StringAppendS(" by ");
  nrzWrite(b, R);
  StringAppendS(" is ");
  number du;
  nrzWrite(du = (number)erg, R);
  StringAppendS(" rest ");
  nrzWrite(du = (number)r, R);
  char * s = StringEndS();
  Print("%s\n", s);
  omFree(s);
#endif

  if (mpz_cmp_si(r, 0)!=0)
  {
    WerrorS("4:Division by non divisible element.");
    WerrorS("Result is without remainder.");
  }
  mpz_clear(r);
  return nrz_short((number) erg);
}

number nrzIntDiv (number a,number b, const coeffs)
{
  assume(SR_TO_INT(b));
  mpz_t aa, bb;
  if (IS_SMALL(a))
    mpz_init_set_si(aa, SR_TO_INT(a));
  else
    mpz_init_set(aa, (int_number) a);
  if (IS_SMALL(b))
    mpz_init_set_si(bb, SR_TO_INT(b));
  else
    mpz_init_set(bb, (int_number) b);
  int_number erg = (int_number) omAllocBin(gmp_nrz_bin);
  mpz_init(erg);
  mpz_tdiv_q(erg, (int_number) aa, (int_number) bb);
  mpz_clear(aa);
  mpz_clear(bb);
  return (number) erg;
}

number nrzIntMod (number a,number b, const coeffs)
{
  mpz_t aa, bb;
  assume(SR_TO_INT(b));
  if (IS_SMALL(a))
    mpz_init_set_si(aa, SR_TO_INT(a));
  else
    mpz_init_set(aa, (int_number) a);
  if (IS_SMALL(b))
    mpz_init_set_si(bb, SR_TO_INT(b));
  else
    mpz_init_set(bb, (int_number) b);

  mpz_t erg;
  mpz_init(erg);
  int_number r = (int_number) omAllocBin(gmp_nrz_bin);
  mpz_init(r);
  mpz_tdiv_qr(erg, r, (int_number) aa, (int_number) bb);
  mpz_clear(erg);
  mpz_clear(aa);
  mpz_clear(bb);
  return (number) r;
}

number  nrzInvers (number c, const coeffs r)
{
  if (!nrzIsUnit((number) c, r))
  {
    WerrorS("Non invertible element.");
    return (number)0; //TODO
  }
  return c; // has to be 1 or -1....
}

number nrzNeg (number c, const coeffs)
{
// nNeg inplace !!!
  if (IS_SMALL(c))
    return INT_TO_SR(-SR_TO_INT(c));
  mpz_mul_si((int_number) c, (int_number) c, -1);
  return c;
}

number nrzMapMachineInt(number from, const coeffs /*src*/, const coeffs /*dst*/)
{
  int_number erg = (int_number) omAllocBin(gmp_nrz_bin);
  mpz_init_set_ui(erg, (NATNUMBER) from);
  return nrz_short((number) erg);
}

number nrzMapZp(number from, const coeffs /*src*/, const coeffs /*dst*/)
{
  int_number erg = (int_number) omAllocBin(gmp_nrz_bin);
  mpz_init_set_si(erg, (long) from);
  return nrz_short((number) erg);
}

number nrzMapQ(number from, const coeffs src, const coeffs /*dst*/)
{
  int_number erg = (int_number) omAllocBin(gmp_nrz_bin);
  mpz_init(erg);
  nlGMP(from, (number) erg, src);
  return nrz_short((number) erg);
}

number nrzModNMap(number from, const coeffs src, const coeffs /*dst*/)
{
  int_number erg = (int_number) omAllocBin(gmp_nrz_bin);
  mpz_init_set(erg, (int_number) from);
  return nrz_short((number) erg);
}

nMapFunc nrzSetMap(const coeffs src, const coeffs /*dst*/)
{
  /* dst = currRing */
  if (nCoeff_is_Ring_ModN(src) || nCoeff_is_Ring_PtoM(src))
    return nrzModNMap;

  if (nCoeff_is_Ring_Z(src))
  {
    return ndCopyMap; //nrzCopyMap;
  }
  if (nCoeff_is_Ring_2toM(src))
  {
    return nrzMapMachineInt;
  }
  if (nCoeff_is_Zp(src))
  {
    return nrzMapZp;
  }
  if (nCoeff_is_Q(src))
  {
    return nrzMapQ;
  }
  return NULL;      // default
}


/*
 * set the exponent (allocate and init tables) (TODO)
 */

void nrzSetExp(int, coeffs)
{
}

void nrzInitExp(int, coeffs)
{
}

#ifdef LDEBUG
BOOLEAN nrzDBTest (number, const char *, const int, const coeffs)
{
  return TRUE;//TODO
}
#endif

void nrzWrite (number &a, const coeffs)
{
  char *s,*z;
  if (a==NULL)
  {
    StringAppendS("o");
  }
  else
  {
    if (IS_SMALL(a)) {
      StringAppend("%d", SR_TO_INT(a));
    } else {
      int l=mpz_sizeinbase((int_number) a, 10) + 2;
      s=(char*)omAlloc(l);
      z=mpz_get_str(s,10,(int_number) a);
      StringAppendS(z);
      omFreeSize((ADDRESS)s,l);
    }
  }
}

/*2
* extracts a long integer from s, returns the rest    (COPY FROM longrat0.cc)
*/
static const char * nlEatLongC(char *s, mpz_ptr i)
{
  const char * start=s;

  if (*s<'0' || *s>'9')
  {
    mpz_set_si(i,1);
    return s;
  }
  while (*s >= '0' && *s <= '9') s++;
  if (*s=='\0')
  {
    mpz_set_str(i,start,10);
  }
  else
  {
    char c=*s;
    *s='\0';
    mpz_set_str(i,start,10);
    *s=c;
  }
  return s;
}

const char * nrzRead (const char *s, number *a, const coeffs)
{
  int_number z = (int_number) omAllocBin(gmp_nrz_bin);
  {
    mpz_init(z);
    s = nlEatLongC((char *) s, z);
  }
  *a = nrz_short((number) z);
  return s;
}

void    nrzCoeffWrite  (const coeffs, BOOLEAN /*details*/)
{
  PrintS("//   coeff. ring is : Integers\n");
}

static char* nrzCoeffString(const coeffs)
{
  return omStrDup("integer");
}

static CanonicalForm nrzConvSingNFactoryN( number n, BOOLEAN setChar, const coeffs /*r*/ )
{
  if (setChar) setCharacteristic( 0 );

  CanonicalForm term;
  if ( IS_SMALL(n))
  {
    term = SR_TO_INT(n);
  }
  else
  {
    mpz_t dummy;
    mpz_init_set( dummy,n->z );
    term = make_cf( dummy );
  }
  return term;
}

static number nrzConvFactoryNSingN( const CanonicalForm n, const coeffs r)
{
  if (n.isImm())
  {
    return nrzInit(n.intval(),r);
  }
  else
  {
    WerrorS("not implemented yet");
    number z;
    return z;
  }
}


BOOLEAN nrzInitChar(coeffs r,  void *)
{
  assume( getCoeffType(r) == ID );
  r->nCoeffIsEqual = ndCoeffIsEqual;
  r->cfCoeffString = nrzCoeffString;
  r->cfKillChar = ndKillChar;
  r->cfMult  = nrzMult;
  r->cfSub   = nrzSub;
  r->cfAdd   = nrzAdd;
  r->cfDiv   = nrzDiv;
  r->cfIntDiv= nrzDiv;
  r->cfIntMod= nrzIntMod;
  r->cfExactDiv= nrzDiv;
  r->cfInit = nrzInit;
  r->cfSize  = nrzSize;
  r->cfInt  = nrzInt;
  //#ifdef HAVE_RINGS
  r->cfDivComp = nrzDivComp; // only for ring stuff
  r->cfIsUnit = nrzIsUnit; // only for ring stuff
  r->cfGetUnit = nrzGetUnit; // only for ring stuff
  r->cfExtGcd = nrzExtGcd; // only for ring stuff
  r->cfXExtGcd = nrzXExtGcd; // only for ring stuff
  r->cfQuotRem = nrzQuotRem;
  r->cfDivBy = nrzDivBy; // only for ring stuff
  r->cfInit_bigint = nrzMapQ;
  //#endif
  r->cfNeg   = nrzNeg;
  r->cfInvers= nrzInvers;
  r->cfCopy  = nrzCopy;
  r->cfWriteLong = nrzWrite;
  r->cfRead = nrzRead;
  r->cfGreater = nrzGreater;
  r->cfEqual = nrzEqual;
  r->cfIsZero = nrzIsZero;
  r->cfIsOne = nrzIsOne;
  r->cfIsMOne = nrzIsMOne;
  r->cfGreaterZero = nrzGreaterZero;
  r->cfPower = nrzPower;
  r->cfGcd  = nrzGcd;
  r->cfLcm  = nrzLcm;
  r->cfDelete= nrzDelete;
  r->cfSetMap = nrzSetMap;
  r->cfCoeffWrite = nrzCoeffWrite;
  r->convSingNFactoryN = nrzConvSingNFactoryN;
  r->convFactoryNSingN = nrzConvFactoryNSingN;

  // debug stuff

#ifdef LDEBUG
  r->cfDBTest=nrzDBTest;
#endif

  r->nNULL = 0;
  r->ch = 0;
  r->has_simple_Alloc=FALSE;
  r->has_simple_Inverse=FALSE;
  return FALSE;
}

#endif
