//usr/bin/env root.exe -l ${0}\(\""${0}"\",\""${*}"\"\); exit $?
// ============================================================================
#include <cmath>
#include <iostream>
#include "Ostap/ValueWithError.h"
// ============================================================================

typedef Ostap::Math::ValueWithError VE ;

namespace Examples
{
  void ex_ve_001 ()
  {

    VE a  ( 1 , 1 ) ;
    VE b  ( 2 , 2 ) ;
    
    std::cout
      << " a          = " << a          << "\n"
      << " b          = " << b          << "\n"
      << "a+b         = " << (  a + b ) << "\n"
      << "a-b         = " << (  a - b ) << "\n"
      << "a*b         = " << (  a * b ) << "\n"
      << "a/b         = " << (  a / b ) << "\n"
      //
      << "a+a         = " << (  a + a ) << "\n"
      << "b-b         = " << (  b - b ) << "\n"
      << "a/a         = " << (  a / a ) << "\n"
      << "b*b         = " << (  b * b ) << "\n"
      //
      << "a/(a+b)     = "     <<  a.frac  ( b ) << "\n"
      << "(a-b)/(a+b) = " <<  a.asym  ( b ) << "\n"
      << "a/(a+a)     = "     <<  a.frac  ( a ) << "\n"
      << "(b-b)/(b+b) = " <<  b.asym  ( b ) << "\n"
      //
      << "a*2         = "  <<  a*2.           << "\n"
      << "2*a         = "  <<  2.*a           << "\n"
      << "a/2         = "  <<  a/2.           << "\n"
      << "2/a         = "  <<  2./a           << "\n"
      //
      << "pow(a,2)    = " << pow  ( a , 2  )  << "\n"
      << "pow(2,a)    = " << pow  ( 2 , a  )  << "\n"
      << "pow(a,b)    = " << pow  ( a , b  )  << "\n"
      << "pow(a,a)    = " << pow  ( a , a  )  << "\n"
      //
      << "mean(a,b)   = " << mean ( a , b  )  << "\n"
      << "chi2(a,b)   = " << chi2 ( a , b  )  << "\n"
      << "mean(a,a)   = " << mean ( a , a  )  << "\n"
      << "chi2(a,a)   = " << chi2 ( a , a  )  << "\n"
      //
      << "-a          = " << -a           << "\n"
      << "abs(a-b)    = " << abs( a - b ) << "\n"
      //
      << "exp(a)      = " << exp(a)       << "\n"
      << "exp2(a)     = " << exp2(a)      << "\n"
      << "expm1(a)    = " << expm1(a)     << "\n"
      << "log(b)      = " << log(b)       << "\n"
      << "log2(b)     = " << log2(b)      << "\n"
      << "log10(b)    = " << log10(b)     << "\n"
      << "log1p(b)    = " << log1p(b)     << "\n"
      //
      << "sqrt(b)     = " << sqrt(b)      << "\n"
      << "cbrt(b)     = " << cbrt(b)      << "\n"
      //
      << "sin(a)      = " << sin(a)       << "\n"
      << "cos(a)      = " << cos(a)       << "\n"
      << "tan(a)      = " << tan(a)       << "\n"
      //
      << "sinh(a)     = " << sinh(a)       << "\n"
      << "cosh(a)     = " << cosh(a)       << "\n"
      << "tanh(a)     = " << tanh(a)       << "\n"
      << "sech(a)     = " << sech(a)       << "\n"
      //
      << "erf(a)      = " << erf(a)       << "\n"
      << "erfc(a)     = " << erfc(a)      << "\n"
      << "erfi(a)     = " << erfi(a)      << "\n"
      << "erfcx(a)    = " << erfcx(a)     << "\n"
      //
      << "acos(a/2)   = " << asin(a/2.)     << "\n"
      << "acos(a/2)   = " << acos(a/2.)     << "\n"
      << "atan(a/2)   = " << atan(a/2.)     << "\n"
      //
      << "acosh(a)    = " << asinh(a)      << "\n"
      << "acosh(a)    = " << acosh(a)      << "\n"
      << "atanh(a/2)  = " << atan(a/2.)     << "\n"
      //
      << "atan2(a,b)  = " << atan2(a,b)    << "\n"
      //
      << "tgamma(b)   = " << tgamma(b)     << "\n"
      << "lgamma(b)   = " << lgamma(b)     << "\n"
      << "igamma(b)   = " << igamma(b)     << "\n"
      //
      << std::endl ;   
  }
}
