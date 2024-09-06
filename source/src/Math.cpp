// $Id$
// ============================================================================
// Include files
// ============================================================================
// STD& STL
// ============================================================================
#include <cstdlib>
#include <cmath>
#include <limits>
// ============================================================================
// GSL 
// ============================================================================
#include "gsl/gsl_sys.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Math.h"
// ============================================================================
/** @file 
 *  Implementation file for functions from namespace LHcb::Math 
 *  @author Vanya BELYAEV ibelyaev@physics.syr.edu
 *  @date   2007-11-27
 */
// ============================================================================
/*  compare two double numbers with the relative
 *  precision   'epsilon'
 *
 *  Essentially it is a wrapper to gsl_fcmp function 
 *  from GSL library
 *
 *  See D.E.Knuth, "Seminumerical Algorithms", section 4.2.2
 *  
 *  @author Vanya BELYAEV ibelyaev@physics.syr.edu
 *  @date 2007-11-27
 */
// ============================================================================
namespace 
{
  static_assert ( std::numeric_limits<long>::is_specialized     && 
                  std::numeric_limits<long>::is_integer         && 
                  std::numeric_limits<long>::is_signed           ,                  
                  "std::numeric_limits<long> is not appropriate" ) ;
  //
  const long   s_long_min  = std::numeric_limits<long>::min () ;
  const long   s_long_max  = std::numeric_limits<long>::max () ;
  const double sd_long_min = float ( std::numeric_limits<long>::min () ) ;
  const double sd_long_max = float ( std::numeric_limits<long>::max () ) ;
}
// ============================================================================
/* compare two double numbers with relative precision 'epsilon'
 *
 *  Essentially it is a wrapper to gsl_fcmp function from GSL library
 *  See D.E.Knuth, "Seminumerical Algorithms", section 4.2.2
 *
 *  @param value1  (INPUT) the first value 
 *  @param value2  (INPUT) the second value 
 *  @param epsilon (INPUT) the (relative) precision 
 *  @author Vanya BELYAEV ibelyaev@physics.syr.edu
 *  @date 2007-11-27
 */
// ============================================================================
bool Ostap::Math::knuth_equal_to_double
( const double value1  , 
  const double value2  , 
  const double epsilon ) 
{
  return 
    !epsilon    ? 0 == gsl_fcmp ( value1 , value2 , 1.0e-6  ) : 
    0 < epsilon ? 0 == gsl_fcmp ( value1 , value2 , epsilon ) : 
    0 == gsl_fcmp ( value1 , value2 , -epsilon ) ;
}
// ============================================================================
namespace 
{
  // ==========================================================================
  /// check the specialization 
  static_assert ( std::numeric_limits<long long> ::is_specialized     , 
                  "std::numeric_limits<long long> is not specialized" ) ;
  static_assert ( std::numeric_limits<long> ::is_specialized     , 
                  "std::numeric_limits<long> is not specialized" ) ;
  static_assert ( std::numeric_limits<int>  ::is_specialized     , 
                  "std::numeric_limits<int>  is not specialized" ) ;
  static_assert ( std::numeric_limits<short> ::is_specialized     , 
                  "std::numeric_limits<short> is not specialized" ) ;
  // ==========================================================================
  const double s_MAX_L  =  0.1 + double ( std::numeric_limits<long>::max      () ) ;
  const double s_MIN_L  = -0.1 - double ( std::numeric_limits<long>::max      () ) ;
  const double s_MAX_LL =  0.1 + double ( std::numeric_limits<long long>::max () ) ;
  const double s_MIN_LL = -0.1 - double ( std::numeric_limits<long long>::max () ) ;
  const double s_MAX_I  =  0.1 + 1.0 * std::numeric_limits<int>::max       () ;
  const double s_MIN_I  = -0.1 - 1.0 * std::numeric_limits<int>::max       () ;
  const double s_MAX_S  =  0.1 + 1.0 * std::numeric_limits<short>::max     () ;
  const double s_MIN_S  = -0.1 - 1.0 * std::numeric_limits<short>::max     () ;
  // ==========================================================================
  static_assert ( std::numeric_limits<unsigned long long> ::is_specialized     , 
                  "std::numeric_limits<unsigned long long> is not specialized" ) ;
  static_assert ( std::numeric_limits<unsigned long>  ::is_specialized     , 
                  "std::numeric_limits<long> is not specialized" ) ;
  static_assert ( std::numeric_limits<unsigned int>   ::is_specialized     , 
                  "std::numeric_limits<unsigned int>  is not specialized" ) ;
  static_assert ( std::numeric_limits<unsigned short> ::is_specialized     , 
                  "std::numeric_limits<unisgned short> is not specialized" ) ;
  // ==========================================================================
  const double s_MAX_UL  =  0.1 + double ( std::numeric_limits<unsigned long>::max      () ) ;
  const double s_MAX_ULL =  0.1 + double ( std::numeric_limits<unsigned long long>::max () ) ;
  const double s_MAX_UI  =  0.1 + 1.0 * std::numeric_limits<unsigned int>::max       () ;
  const double s_MAX_US  =  0.1 + 1.0 * std::numeric_limits<unsigned short>::max     () ;
  // ==========================================================================
}
// ============================================================================
/*  is the value actually long ?
 *  @author Vanya BELYAEV Ivan.Belyaev       
 *  @date 2011-07-18
 */
// ============================================================================
bool Ostap::Math::islong ( const double x ) 
{
  return 
    x <= s_MIN_L  ? false :
    x >= s_MAX_L  ? false :
    lomont_compare_double ( x                   , 
                            std::lround   ( x ) , 
                            mULPS_double        ) ;
}
// ============================================================================
/*  is the value actually int ?
 *  @author Vanya BELYAEV Ivan.Belyaev       
 *  @date 2011-07-18
 */
// ============================================================================
bool Ostap::Math::isint ( const double x ) 
{
  return 
    x <= s_MIN_I  ? false :
    x >= s_MAX_I  ? false :
    lomont_compare_double ( x , round ( x ) , mULPS_double ) ;
}
// ============================================================================
/*  is the value actually int ?
 *  @author Vanya BELYAEV Ivan.Belyaev       
 *  @date 2011-07-18
 */
// ============================================================================
bool Ostap::Math::isint ( const float x ) 
{
  return 
    x <= s_MIN_I  ? false :
    x >= s_MAX_I  ? false :
    lomont_compare_double ( x , round ( x ) , mULPS_float ) ;
}
// ============================================================================
/*  is the value actually short int ?
 */
// ============================================================================
bool Ostap::Math::isshort ( const double x ) 
{
  return 
    x <= s_MIN_S  ? false :
    x >= s_MAX_S  ? false :
    lomont_compare_double ( x , round ( x ) , mULPS_double ) ;
}
// ============================================================================
/*  is the value actually long long ?
 *  @author Vanya BELYAEV Ivan.Belyaev       
 *  @date 2011-07-18
 */
// ============================================================================
bool Ostap::Math::islonglong ( const double x ) 
{
  const long double x_ = x  ;
  return 
    x <= s_MIN_LL  ? false :
    x >= s_MAX_LL  ? false :
    lomont_compare_double ( x , std::llround ( x_ ) , mULPS_double ) ;
}
// ============================================================================
/*  is the value actually unsigned int ?
 *  @author Vanya BELYAEV Ivan.Belyaev       
 *  @date 2011-07-18
 */
// ============================================================================
bool Ostap::Math::isuint ( const double x ) 
{
  const long double x_ = x  ;
  return 
    x <= -0.1      ? false :
    x >= s_MAX_UI  ? false :
    lomont_compare_double ( x , std::llround ( x ) , mULPS_double ) ;
}
// ============================================================================
/*  is the value actually unsigned short int ?
 */
// ============================================================================
bool Ostap::Math::isushort ( const double x ) 
{
  return 
    x <= -0.1      ? false :
    x >= s_MAX_US  ? false :
    lomont_compare_double ( x , round ( x ) , mULPS_double ) ;
}
// ============================================================================
/*  is the value actually unsigned long ?
 *  @author Vanya BELYAEV Ivan.Belyaev       
 *  @date 2011-07-18
 */
// ============================================================================
bool Ostap::Math::isulong ( const double x ) 
{
  const long double x_ = x  ;
  return 
    x <= -0.1      ? false :
    x >= s_MAX_UL  ? false :
    lomont_compare_double ( x , std::llround ( x ) , mULPS_double ) ;
}
// ============================================================================
/*  is the value actually unisgned long long ?
 *  @author Vanya BELYAEV Ivan.Belyaev       
 *  @date 2011-07-18
 */
// ============================================================================
bool Ostap::Math::isulonglong ( const double x ) 
{
  const long double x_ = x  ;
  return 
    x <= -0.1      ? false :
    x >= s_MAX_ULL ? false :
    lomont_compare_double ( x , std::llround ( x_ ) , mULPS_double ) ;
}
// ============================================================================
/* check if the double value is actually equal to the integer value  
 *  @param val value to be compared with the integer 
 *  @param ref the reference integer number 
 *  @param mULPS the precision 
 *  @see Ostap::Math::lomont_compare_double 
 *  @see LHCb::Math::mULPS_double
 *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
 *  @date 2008-09-17
 */        
// ============================================================================
bool Ostap::Math::equal_to_int 
( const double       val   , 
  const int          ref   , 
  const unsigned int mULPS ) 
{
  const double tmp = ref ;
  return lomont_compare_double ( val , tmp , mULPS ) ;
}
// ============================================================================
/*  check if the double value is actually equal to the unsigned integer value  
 *  @param val value to be compared with the unsigned integer 
 *  @param ref the reference unsigned integer number 
 *  @param mULPS the precision 
 *  @see Ostap::Math::lomont_compare_double 
 *  @see LHCb::Math::mULPS_double
 *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
 *  @date 2008-09-17
 */        
// ============================================================================
bool Ostap::Math::equal_to_uint 
( const double       val   , 
  const unsigned int ref   , 
  const unsigned int mULPS ) 
{
  const double tmp = ref ;
  return lomont_compare_double ( val , tmp , mULPS ) ;
}
// ============================================================================
/*  get mantissa and exponent 
 *  similar to std::frexp, but radix=10)
 *  @param x  INPUT  value 
 *  @param e  UPDATE exponent 
 *  @return  pair of mantissa (0.1<=m<1) and exponent 
 *  @author Vanya BELYAEV Ivan.Belyaev       
 *  @date 2015-07-21
 */
// ============================================================================
std::pair<double,int> 
Ostap::Math::frexp10 ( const double x ) 
{
  //
  if ( 0 == x )  { return std::pair<double,int>( x , 0 ) ; }
  //
  long double xa = 0 < x ?  x : -x ;
  //
  int q  = (int) std::floor ( std::log10 ( xa ) ) ;
  if      ( 0 < q ) { xa /= std::pow ( 10.0L , (unsigned long)           q   ) ; }
  else if ( 0 > q ) { xa *= std::pow ( 10.0L , (unsigned long) std::abs( q ) ) ; }
  //
  if ( 1 <= xa ) { xa /= 10 ; ++q ; }
  //
  return 0 < x ? 
    std::pair<double,int>(  xa , q ) : 
    std::pair<double,int>( -xa , q ) ;
}
// ============================================================================
/*  get mantissa and binary exponent 
 *  similar to std::frexp
 *  @param x  INPUT  value 
 *  @return   pair of mantissa (0.5<=m<1) and exponent 
 *  @author Vanya BELYAEV Ivan.Belyaev       
 *  @date 2015-07-21
 */
// ============================================================================
std::pair<double,int>
Ostap::Math::frexp2 ( const double x ) 
{
  int    e = 0 ;
  double m = std::frexp ( x , &e ) ;
  return std::make_pair ( m , e ) ;
}
// ============================================================================
/*  get mantissa and exponent 
 *  similar to std::frexp, but radix=10)
 *  @param x  INPUT  value 
 *  @param e  UPDATE exponent 
 *  @return  mantissa  (0.1<=m<1)
 *  @author Vanya BELYAEV Ivan.Belyaev       
 *  @date 2011-07-18
 */
// ============================================================================
double Ostap::Math::frexp10 ( const double  x , int& e ) 
{
  //
  const std::pair<double,int> r = frexp10 ( x ) ;
  e    = r.second ;
  return r.first  ;
}
// ============================================================================
/*  get mantissa and exponent 
 *  similar to std::frexp, but radix=10)
 *  @param x  INPUT  value 
 *  @param e  UPDATE exponent 
 *  @return  mantissa  (0.1<=m<1)
 *  @author Vanya BELYAEV Ivan.Belyaev       
 *  @date 2011-07-18
 */
// ============================================================================
float Ostap::Math::frexp10 ( const float  x , int& e ) 
{
  const double xx = x ;
  return frexp10 ( xx , e ) ;
}
// ============================================================================
/*  Multiplies a floating point value num by the number 10(radix)
 *  raised to the exp power. 
 *  similar to std::lxep, but radix is 10 
 *  @param num  input value 
 *  @param exp  the 10-base exponent 
 *  @return     propertly scaled input value 
 *  @attention  it is not very efficient, but OK for our purpose
 */
// ============================================================================
double Ostap::Math::ldexp10
( const double value ,
  const short  expo  )
{
  if ( 0 == value )  { return 0 ; }
  //
  long double xa = 0 < value ? value : -value ;
  //
  const int qq  = (int) std::floor ( std::log10 ( xa ) ) ;
  //
  static const int     s_maxe10 =  std::numeric_limits<double>::max_exponent10 + 3 ;
  static const double  s_max    =  std::numeric_limits<double>::max () ;
  //
  if ( qq + expo > s_maxe10 ) { return 0 <= value ? s_max : -s_max ; }
  //
  xa *= std::pow ( 10.0L , expo ) ;
  const long double xx= 0 <= value ? xa  : -xa ;
  //
  return s_max <= xa ? ( 0 <= value ? s_max : -s_max ) : xx ;
}
// ============================================================================
/*  round to N-significant digits 
 *  @param x  INPUT  input value 
 *  @param n  INPUT  number of significant digits 
 *  @return rounded value 
 *  @author Vanya BELYAEV Ivan.Belyaev       
 *  @date 2015-07-21
 */
// ============================================================================
double Ostap::Math::round_N ( const double x , const unsigned short n ) 
{
  //
  if ( 0 == n ) { return 0 ; } // none of correct digits is required 
  //
  const std::pair<double,int> r = frexp10 ( x ) ;
  //
  const long    e  = r.second -  1 ;
  const double  m  = r.first  * 10 ;
  //
  const long    ni = n - 1 ;
  //
  const double  f1 = std::pow ( 10.0 , ni )  ;
  const double  f2 =
    ni < e ?     std::pow ( 10.0 , (unsigned long) ( e  - ni ) ) : 
    ni > e ? 1.0/std::pow ( 10.0 , (unsigned long) ( ni - e  ) ) : 1 ;
//
return round ( m * f1 ) * f2 ;
}
// ============================================================================
/*  round to N-significant digits 
 *  @param x  INPUT  input value 
 *  @param n  INPUT  number of significnat digits 
 *  @return rounded value 
 *  @author Vanya BELYAEV Ivan.Belyaev       
 *  @date 2015-07-21
 */
// ============================================================================
float Ostap::Math::round_N ( const float x , const unsigned short n ) 
{
  const double xd = x ;
  return round_N ( xd , n ) ;
}
// ========================================================================
/** round to nearest integer, rounds half integers to nearest even integer 
 *  It is just a simple wrapper around boost::numeric::converter 
 *  @author Vanya BELYAEV Ivan.BElyaev
 */
// ========================================================================
long Ostap::Math::round ( const double x ) 
{
  return 
    x <= sd_long_min ? s_long_min :
    x >= sd_long_max ? s_long_max : long ( std::lround ( x )  ) ;
}
// ============================================================================


// ============================================================================
//                                                                      The END 
// ============================================================================
