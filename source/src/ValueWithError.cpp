// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <cmath>
#include <string>
#include <sstream>
#include <climits>
// ============================================================================
// Ostap 
// ============================================================================
#include "Ostap/StatusCode.h"
#include "Ostap/Math.h"
#include "Ostap/MoreMath.h"
#include "Ostap/Clenshaw.h"
#include "Ostap/Combine.h"
#include "Ostap/Interpolation.h"
#include "Ostap/ValueWithError.h"
// ============================================================================
// Local
// ============================================================================
#include "Exception.h"
#include "format.h"
// ============================================================================
/** @file
 *  Implementation file for class Gaudi::Math::ValueWithError
 *  @date 2009-06-03
 *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
 */
// ============================================================================
// local namespace to hide the details
// ============================================================================
namespace
{
  // ==========================================================================
  // const unsigned int _maxULPs = 10000 ;
  // ==========================================================================
  /// equality of doubles 
  const Ostap::Math::Equal_To<double> s_equal{} ;
  /// equality of doubles
  inline bool _equal ( const double value1 ,
                       const double value2 )
  { return value1 == value2 || s_equal ( value1 , value2 ) ; }
  // ==========================================================================
  // check if the double value close to zero
  const Ostap::Math::Zero<double> s_zero{} ;
  // check if the double value close to zero
  inline bool _zero  ( const double value ) { return s_zero ( value) ; }
  // check if the double value close to zero
  inline bool _zero  ( const Ostap::Math::ValueWithError& a ) 
  { return s_zero ( a.value() ) && s_zero ( a.cov2() ) ; }  
  // ==========================================================================
  // check if the double value close to one
  inline bool _one   ( const double value ) { return _equal ( value , 1 ) ; }
  // check if the double value close to one
  inline bool _one  ( const Ostap::Math::ValueWithError& a ) 
  { return _one ( a.value() ) && _zero ( a.cov2() ) ; }  
  // ==========================================================================
  inline bool _is_long ( const double value ) 
  { return Ostap::Math::islong ( value ) ; }
  // ==========================================================================
  /// precomputed value of ln(2) 
  const double s_ln2    = std::log( double ( 2 ) ) ;
  /// precomputed value of ln(2) squared 
  const double s_ln2_sq = s_ln2 * s_ln2 ;  
  // precomputed value of 1/ln(10) 
  const double s_ln10_i =  1 / std::log ( double(10) ) ;
  // precomputed value of 1/ln(2) 
  const double s_ln2_i  =  1 / std::log ( double( 2) ) ;
  // ==========================================================================
}
// ============================================================================
// constructor from the value and covariance
// ============================================================================
Ostap::Math::ValueWithError::ValueWithError
( const double value      ,
  const double covariance )
  : m_value ( value      )
  , m_cov2  ( covariance )
{
  if ( _zero ( m_cov2 ) ) { m_cov2 = 0 ; }
}
// ============================================================================
// constructor from the (value,error)-pair
// ============================================================================
Ostap::Math::ValueWithError::ValueWithError
( const std::pair<double,double>& value )
  : m_value ( value.first )
  , m_cov2  ( 0 )
{
  setError ( value.second ) ;
}
// ============================================================================
// set the error
// ============================================================================
void Ostap::Math::ValueWithError::setError ( const double e )
{
  if    ( _zero ( e ) ) { m_cov2 = 0 ; }
  else 
  {
    m_cov2  = e * e ;
    //
    if ( 0 > e ) { m_cov2 = -m_cov2 ; }
    //
    if ( _zero ( m_cov2 ) ) { m_cov2 = 0 ; }
  }
}
// ============================================================================
// get the error
// ============================================================================
double Ostap::Math::ValueWithError::error      () const
{ 
  return 
    _zero ( m_cov2 ) ? 0.0 : 
    0 <=    m_cov2   ? std::sqrt ( m_cov2 ) : -std::sqrt ( -m_cov2 ) ; 
}
// ============================================================================
// +=
// ============================================================================
Ostap::Math::ValueWithError&
Ostap::Math::ValueWithError::operator+=
( const Ostap::Math::ValueWithError& right )                             // +=
{
  //
  if ( &right == this ) 
  {
    m_value *= 2  ;
    m_cov2  *= 4  ;
    return  *this ;
  }
  //
  m_value += right.m_value ;
  if ( 0 < right.m_cov2 ) { m_cov2  += right.m_cov2  ; }
  //
  if ( _zero ( m_cov2 ) ) { m_cov2 = 0 ; }
  //
  return *this ;
}
// ============================================================================
// -=
// ============================================================================
Ostap::Math::ValueWithError&
Ostap::Math::ValueWithError::operator-=
( const Ostap::Math::ValueWithError& right )                              // -=
{
  //
  if ( &right == this ) 
  {
    m_value = 0   ;
    m_cov2  = 0   ;
    return  *this ;
  }
  //
  m_value -= right.m_value ;
  if ( 0 < right.m_cov2 ) { m_cov2  += right.m_cov2  ; }
  //
  if ( _zero ( m_cov2 ) ) { m_cov2 = 0 ; }
  //
  return *this ;
}
// ============================================================================
// *=
// ============================================================================
Ostap::Math::ValueWithError&
Ostap::Math::ValueWithError::operator*=
( const Ostap::Math::ValueWithError& right )                              // *=
{
  if ( &right == this ) 
  {
    const double a = value() ;
    m_value  =     a * a ;
    m_cov2  *= 4 * a * a ;
    //
    if ( _zero ( m_cov2 ) ) { m_cov2 = 0 ; }
    //
    return  *this ;
  }
  //
  const double _a2 =       m_value *       m_value ;
  const double _b2 = right.m_value * right.m_value ;
  m_cov2  *= _b2                 ;
  if ( 0 < right.m_cov2 ) { m_cov2  += _a2 * right.m_cov2  ; }
  m_value *=      right.m_value ;
  //
  if ( _zero ( m_cov2 ) ) { m_cov2 = 0 ; }
  //
  return *this ;
}
// ============================================================================
// /=
// ============================================================================
Ostap::Math::ValueWithError&
Ostap::Math::ValueWithError::operator/=
( const Ostap::Math::ValueWithError& right )                              // /=
{
  if ( &right == this ) 
  {
    m_value  =  1 ;
    m_cov2   =  0 ;
    return  *this ;
  }
  //
  const double _a2 =       m_value *       m_value ;
  const double _b2 = right.m_value * right.m_value ;
  const double _b4 = _b2 * _b2 ;
  //
  m_cov2  /= _b2 ;
  if ( 0 < right.m_cov2 ) { m_cov2  += ( _a2 / _b4 ) * right.m_cov2 ; }
  m_value /= right.m_value ;
  //
  if ( _zero ( m_cov2 ) ) { m_cov2 = 0 ; }
  //
  return *this ;
}
// ============================================================================
// *=
// ============================================================================
Ostap::Math::ValueWithError&
Ostap::Math::ValueWithError::operator*= ( const double v )                // *=
{
  m_value *= v     ;
  m_cov2  *= (v*v) ;
  //
  if ( _zero ( m_cov2 ) ) { m_cov2 = 0 ; }
  //
  return *this ;
}
// ============================================================================
// /=
// ============================================================================
Ostap::Math::ValueWithError&
Ostap::Math::ValueWithError::operator/= ( const double v )                // /=
{
  m_value /= v     ;
  m_cov2  /= (v*v) ;
  //
  if ( _zero ( m_cov2 ) ) { m_cov2 = 0 ; }
  //
  return *this ;
}
// ============================================================================
// +=
// ============================================================================
Ostap::Math::ValueWithError&
Ostap::Math::ValueWithError::operator+=( const double right )             // +=
{
  m_value += right ;
  return *this ;
}
// ============================================================================
// -=
// ============================================================================
Ostap::Math::ValueWithError&
Ostap::Math::ValueWithError::operator-= ( const double right )            // -=
{
  m_value -= right ;
  return *this ;
}
// ============================================================================
// unary-
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::operator-() const                        // unary-
{ return ValueWithError( -value() , cov2() ) ; }
// ============================================================================
// printout
// ============================================================================
std::ostream&
Ostap::Math::ValueWithError::fillStream ( std::ostream& s ) const
{ return s << "( " << m_value << " +- " << error() << " )" ; }
// ============================================================================
// printout using format
// ============================================================================
std::ostream&
Ostap::Math::ValueWithError::fillStream
( std::ostream&      s   ,
  const std::string& fmt ) const
{ return s << Ostap::format ( fmt , value() , error() ) ; }
// ============================================================================
// conversion to string
// ============================================================================
std::string Ostap::Math::ValueWithError::toString   () const
{
  std::ostringstream s ;
  fillStream ( s ) ;
  return s.str() ;
}
// ============================================================================
// conversion to the string using format
// ============================================================================
std::string Ostap::Math::ValueWithError::toString
( const std::string& fmt ) const
{
  std::ostringstream s ;
  fillStream ( s , fmt ) ;
  return s.str() ;
}
// ============================================================================
// evaluate the mean of a and b
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::mean
( const Ostap::Math::ValueWithError& b ) const
{
  //
  if ( &b == this ) { return *this ; } // self-mean ???
  //
  if      ( 0 >=   cov2 () && 0 >= b.cov2 () ) 
  { return 0.5 * ( value() + b.value() ) ; }
  else if ( 0 >=   cov2 ()                   ) { return *this ; }
  else if ( 0 >= b.cov2 ()                   ) { return b     ; }
  //
  double _cov2 = 1.0/( 1.0/cov2() + 1.0/b.cov2() ) ;
  if ( _zero ( _cov2 ) ) { _cov2 = 0 ; }
  //
  return Ostap::Math::ValueWithError
    ( _cov2 * ( value()/cov2() + b.value()/b.cov2() ) ,  _cov2 ) ;
}
// =============================================================================
// evaluate chi2
// =============================================================================
double Ostap::Math::ValueWithError::chi2
( const Ostap::Math::ValueWithError& b ) const
{
  //
  if ( _equal ( value () , b.value() ) ) { return 0 ; } // RETURN
  //
  const double s_cov2 = cov2() + b.cov2() ;
  if ( 0 >= s_cov2 )                     { return -1 ; } // RETURN
  //
  const double diff = value() - b.value() ;
  return diff*diff/s_cov2 ;
}
// =============================================================================
// evaluate chi2
// =============================================================================
double Ostap::Math::ValueWithError::chi2 ( const double b ) const
{
  //
  if ( _equal ( value() , b ) ) { return  0 ; } // RETURN
  //
  if ( 0 >= cov2 ()           ) { return -1 ; } // RETURN
  const double diff = value() - b ;
  return diff*diff/cov2() ;
}
// =============================================================================
/*  get Kullback-Liebler divergency 
 *  return the divergency for valid arguments, -1 otherwise
 */
// =============================================================================
double Ostap::Math::ValueWithError::kullback
( const Ostap::Math::ValueWithError& b ) const
{
  //
  if ( 0 >= cov2() || 0 >= b.cov2 () ) { return -1 ; }
  //
  const double c1 =   cov2 () ;
  const double c2 = b.cov2 () ;
  //
  return ( c1 - c2 ) * ( 1.0 / c2 - 1.0 / c1 ) + chi2 ( b ) ;  
}
// =============================================================================
/*  get (squared) Hellinger distance
 *  @see https://en.wikipedia.org/wiki/Hellinger_distance
 *  @return heilinger distance for valid arguments, -1 otherwise
 */
// =============================================================================
double Ostap::Math::ValueWithError::hellinger2
( const ValueWithError& right ) const 
{
  const bool n1 = 0 >=       cov2() || s_zero (       cov2() ) ;
  const bool n2 = 0 >= right.cov2() || s_zero ( right.cov2() ) ;
  //
  if       ( n1  && n2  ) { return -1 ; }
  else if  ( n1  || n2  ) { return  1 ; }
  //
  const double m1 =       value() ;
  const double m2 = right.value() ;
  const double dm = m1 - m2 ;
  //
  const double sq1 =       cov2() ;
  const double sq2 = right.cov2() ;
  //
  return 1 - std::sqrt ( 2.0 *std::sqrt( sq1 * sq2 ) / ( sq1 + sq2 ) ) * std::exp ( -0.25 * dm * dm / ( sq1 + sq2 ) ) ;
}
// =============================================================================
// evaluate residual: signed sqrt(chi2)
// =============================================================================
double Ostap::Math::ValueWithError::residual
( const Ostap::Math::ValueWithError& b ) const
{
  //
  if ( _equal ( value () , b.value() ) ) { return     0 ; } // RETURN
  //
  const double s_cov2 = cov2() + b.cov2() ;
  if ( 0 >= s_cov2 )                     { return -1000 ; } // RETURN
  //
  const double diff = value() - b.value() ;
  //
  return diff / std::sqrt ( s_cov2 ) ;
}
// =============================================================================
// evaluate residual: signed sqrt(chi2)
// =============================================================================
double Ostap::Math::ValueWithError::residual
( const double b ) const
{
  //
  if ( _equal ( value() , b ) ) { return     0 ; } // RETURN
  //
  if ( 0 >= cov2 () )           { return -1000 ; } // RETURN
  //
  const double diff = value() - b ;
  //
  return diff / error () ;
}
// ============================================================================
/*  evaluate the "fraction" \f$  \frac{a}{a+b} \f$
 *  @param  b the parameter "b" for the fraction
 *  @return a/(a+b)
 */
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::frac
( const Ostap::Math::ValueWithError& b ) const
{
  //
  if ( &b == this ) { return fraction ( *this , b , 1 ) ; }
  //
  const double r  = value() / ( value() + b.value() ) ;
  //
  const double s  = value() + b.value() ;
  const double s2 = s  * s  ;
  const double s4 = s2 * s2 ;
  const double c2 =
    std::fabs (   cov2 () ) * b.value () * b.value () +
    std::fabs ( b.cov2 () ) *   value () *   value () ;
  //
  return ValueWithError
    ( r , 0 <= cov2() && 0 <= b.cov2() ? c2/s4 : -1.0 * c2 / s4 ) ;
  //
}
// ============================================================================
/*  evaluate the "fraction" \f$  \frac{a}{a+b} \f$
 *  @param  b the parameter "b" for the fraction
 *  @return a/(a+b)
 */
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::frac ( const double  b ) const
{ return frac ( ValueWithError ( b ) ) ; }
// ============================================================================
/*  evaluate the "asymmetry" \f$  \frac{a-b}{a+b} \f$
 *  @param  b the parameter "b" for the fraction
 *  @return (a-b)/(a+b)
 */
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::asym
( const Ostap::Math::ValueWithError& b ) const
{
  //
  if ( &b == this ) { return asymmetry ( *this , b , 1 ) ; }
  //
  const double r  = ( value() - b.value() ) / ( value() + b.value() ) ;
  //
  const double s  = value() + b.value() ;
  const double s2 = s  * s  ;
  const double s4 = s2 * s2 ;
  //
  const double c2 =
    4 * std::fabs (   cov2 () ) * b.value () * b.value () +
    4 * std::fabs ( b.cov2 () ) *   value () *   value () ;
  //
  return ValueWithError
    ( r , 0 <= cov2() && 0 <= b.cov2() ? c2/s4 : -1.0 * c2 / s4 ) ;
  //
}
// ============================================================================
/*  evaluate the "asymmetry" \f$  \frac{a-b}{a+b} \f$
 *  @param  b the parameter "b" for the fraction
 *  @return (a-b)/(a+b)
 */
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::asym ( const double  b ) const
{ return asym ( ValueWithError ( b ) ) ; }
// =============================================================================
// check for NaN
// =============================================================================
bool Ostap::Math::ValueWithError::isnan    () const
{
  return std::isnan    ( m_value ) || std::isnan    ( m_cov2  )  ;
}
// =============================================================================
// check for finiteness
// =============================================================================
bool Ostap::Math::ValueWithError::isfinite () const
{
  return std::isfinite ( m_value ) && std::isfinite ( m_cov2  )  ;
}
// =============================================================================
// check for finiteness
// =============================================================================
bool Ostap::Math::ValueWithError::isnormal () const
{ return std::isnormal ( m_value ) && std::isfinite ( m_cov2  )  ; }
// =============================================================================
// check for finiteness
// =============================================================================
bool Ostap::Math::ValueWithError::isinf () const
{ return std::isinf ( m_value ) || std::isinf ( m_cov2  )  ; }
// ============================================================================
// check for goodness: finite values and non-negative covariance 
// ============================================================================
bool Ostap::Math::ValueWithError::isgood   () const 
{ return isfinite () &&  ( 0 <= m_cov2 || _zero ( m_cov2 ) ) ; }
// =============================================================================
// for easy pythonization
// =============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__add__
( const Ostap::Math::ValueWithError& right ) const
{
  //
  if ( &right == this ) { return right * 2.0 ; }
  //
  ValueWithError tmp ( *this ) ;
  return tmp += right ;
}
// =============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__sub__
( const Ostap::Math::ValueWithError& right ) const
{
  //
  if ( &right == this ) { return  ValueWithError(0,0) ; }
  //
  ValueWithError tmp ( *this ) ;
  return tmp -= right ;
}
// =============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__mul__
( const Ostap::Math::ValueWithError& right ) const
{
  //
  if ( &right == this ) { return  pow ( *this , 2 ) ; }
  //
  ValueWithError tmp ( *this ) ;
  return tmp *= right ;
}
// =============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__div__
( const Ostap::Math::ValueWithError& right ) const
{
  //
  if ( &right == this ) { return  ValueWithError ( 1 , 0 ) ; }
  //
  ValueWithError tmp ( *this ) ;
  return tmp /= right ;
}
// =============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__add__ ( const double right ) const
{
  ValueWithError tmp ( *this ) ;
  return tmp += right ;
}
// =============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__sub__ ( const double right ) const
{
  ValueWithError tmp ( *this ) ;
  return tmp -= right ;
}
// =============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__mul__ ( const double right ) const
{
  ValueWithError tmp ( *this ) ;
  return tmp *= right ;
}
// =============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__div__ ( const double right ) const
{
  ValueWithError tmp ( *this ) ;
  return tmp /= right ;
}
// =============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__rsub__ ( const double right ) const
{ return ValueWithError ( right - value() , cov2()  ) ; }
// =============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__rdiv__ ( const double right ) const
{  
  ValueWithError tmp ( right ) ;
  return tmp /= (*this) ;
}
// ============================================================================
// abs(a)
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__abs__ () const
{ return Ostap::Math::abs ( *this ) ; }
// ============================================================================
// me**e
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__pow__  ( const int             e ) const
{ return pow ( *this , e ) ; }
// ============================================================================
// me**e
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__pow__  ( const double          e ) const
{ return pow ( *this , e ) ; }
// ============================================================================
// me**e
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__pow__
( const Ostap::Math::ValueWithError&  e ) const
{ return pow ( *this , e ) ; }
// ============================================================================
// e**me
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__rpow__  ( const int             e ) const
{ return pow ( e , *this ) ; }
// ============================================================================
// e**me
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__rpow__ ( const double          e ) const
{ return pow ( e , *this ) ; }
// ============================================================================
// -me
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__neg__() const
{ return Ostap::Math::ValueWithError ( -value() , cov2() ) ; }
// ============================================================================
// +me (no-effect)
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__pos__() const { return *this ; }
// ============================================================================
// exp(me)
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__exp__   () const { return exp   ( *this ) ; }
// ============================================================================
// exp2(me)
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__exp2__  () const { return exp2  ( *this ) ; }
// ============================================================================
// expm1(me)
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__expm1__ () const { return expm1 ( *this ) ; }
// ============================================================================
// log(me)
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__log__   () const { return log   ( *this ) ; }
// ============================================================================
// log2(me)
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__log2__  () const { return log2  ( *this ) ; }
// ============================================================================
// log10(me)
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__log10__ () const { return log10 ( *this ) ; }
// ============================================================================
// log1p(me)
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__log1p__ () const { return log1p ( *this ) ; }
// ============================================================================
// sqrt(me)
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__sqrt__  () const { return sqrt ( *this ) ; }
// ============================================================================
// cbrt(me)
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__cbrt__  () const { return cbrt ( *this ) ; }
// ============================================================================
// sin(me)
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__sin__   () const { return sin  ( *this ) ; }
// ============================================================================
// cos(me)
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__cos__   () const { return cos  ( *this ) ; }
// ============================================================================
// tan(me)
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__tan__   () const { return tan  ( *this ) ; }
// ============================================================================
// sinh(me)
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__sinh__  () const { return sinh ( *this ) ; }
// ============================================================================
// cosh(me)
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__cosh__  () const { return cosh ( *this ) ; }
// ============================================================================
// tanh(me)
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__tanh__  () const { return tanh ( *this ) ; }
// ============================================================================
// erf(me)
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__erf__   () const { return erf   ( *this ) ; }
// ============================================================================
// erfc(me)
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__erfc__  () const { return erfc  ( *this ) ; }
// ============================================================================
// asin(me)
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__asin__  () const { return asin  ( *this ) ; }
// ============================================================================
// acos(me)
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__acos__  () const { return acos  ( *this ) ; }
// ============================================================================
// atan(me)
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__atan__  () const { return atan  ( *this ) ; }
// ============================================================================
// asinh(me)
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__asinh__ () const { return asinh ( *this ) ; }
// ============================================================================
// acosh(me)
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__acosh__ () const { return acosh ( *this ) ; }
// ============================================================================
// atanh(me)
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__atanh__ () const { return atanh ( *this ) ; }
// ============================================================================
// tgamma(me)
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__tgamma__ () const { return tgamma ( *this ) ; }
// ============================================================================
// lgamma(me)
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ValueWithError::__lgamma__ () const { return lgamma ( *this ) ; }
// ============================================================================


// ============================================================================
/* Does this object represent natural number?
 *  - non-negative integer value 
 *  - cov2 == value  or cov2 == 0 
 */
// =============================================================================
bool Ostap::Math::natural_number
( const Ostap::Math::ValueWithError& v ) 
{
  return 
    0 <= v.value() && 0<= v.cov2() 
    && _is_long ( v.value() ) 
    && ( _zero ( v.cov2 () ) || _equal ( v.value() , v.cov2() ) ) ;
}
// ============================================================================
/** Does this object represent natural entry in histogram
 *  - non-negative integer value 
 *  - cov2 == value  or ( 0 == value && 1 == cov2 )
 */
// =============================================================================
bool Ostap::Math::natural_entry 
( const Ostap::Math::ValueWithError& v ) 
{
  return 
    0 <= v.value() && 0<= v.cov2() 
    && _is_long ( v.value() ) 
    && ( _equal ( v.value() , v.cov2() ) ||
         ( _zero ( v.value() ) && _equal ( 1 , v.cov2() ) ) ) ;
}
// ============================================================================
/*  evaluate the mean of a and b 
 *  taking into account correlation coefficient <code>rho</code>
 *  @param a (INPUT) the first argument 
 *  @param b (INPUT) the second argument 
 *  @param rho (INPUT) correlation coefficient \f$-1\le\rhi\le 1\f$
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::mean
( const Ostap::Math::ValueWithError& a   , 
  const Ostap::Math::ValueWithError& b   , 
  const double                       rho ) { return combine ( a , b , rho ) ; }
// ============================================================================
/*  evaluate abs(a)
 *  @param a (INPUT) the value
 *  @return the absolute value
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::abs
( const Ostap::Math::ValueWithError& a )
{ return ValueWithError ( std::fabs ( a.value() ) , a.cov2() ) ; }
// ============================================================================
/* evaluate the binomial efficiency for Bernulli scheme with
 *  @param n (INPUT) number of 'success'
 *  @param N (INPUT) total number
 *  @return the binomial efficiency
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::binomEff
( const size_t n ,
  const size_t N )
{
  if       ( n >  N ) { return binomEff       ( N , n ) ; }
  else if  ( 0 == N ) { return ValueWithError ( 1 , 1 ) ; }
  //
  const long n1 = 0 == n ? 1 :     n ;
  const long n2 = n == N ? 1 : N - n ;
  //
  const double eff = double ( n       ) / N         ;
  const double c2  = double ( n1 * n2 ) / N / N / N ;
  //
  return Ostap::Math::ValueWithError  ( eff , c2 ) ;
}
// ============================================================================
/*  evaluate the binomial efficiency interval using Wilson's prescription
 *  @param n (INPUT) number of 'success'
 *  @param N (INPUT) total number
 *  @return the binomial efficiency
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::wilsonEff
( const size_t n ,
  const size_t N )
{
  //
  if      ( n >  N ) { return wilsonEff      ( N , n ) ; }
  else if ( 0 == N ) { return ValueWithError ( 1 , 1 ) ; }
  //
  const long n1       = 0 == n ? 1 :     n ;
  const long n2       = n == N ? 1 : N - n ;
  //
  const double p      = double ( n1 ) / N ;
  const double q      = double ( n2 ) / N ;
  //
  const double kappa  =             1 ; // "1*sigma"
  const double kappa2 = kappa * kappa ;
  //
  const double nK     = N + kappa2 ;
  const double eff    = ( n + 0.5 * kappa2 ) / nK ;
  //
  const double prefix = kappa2 * N / ( nK * nK ) ;
  const double c2     = prefix * ( q * p + 0.25 * kappa2 / N ) ;
  //
  return Ostap::Math::ValueWithError  ( eff , c2 ) ;
}
// ============================================================================
/*  evaluate the binomial efficiency interval using Agresti-Coull's prescription
 *  @param n (INPUT) number of 'success'
 *  @param N (INPUT) total number
 *  @return the binomial efficiency
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::agrestiCoullEff
( const size_t n ,
  const size_t N )
{
  //
  if      ( n >  N ) { return wilsonEff      ( N , n ) ; }
  else if ( 0 == N ) { return ValueWithError ( 1 , 1 ) ; }
  //
  const double kappa  =             1 ; // "1*sigma"
  const double kappa2 = kappa * kappa ;
  //
  const double n1 = n + 0.5 * kappa2 ;
  const double n2 = N +       kappa2 ;
  //
  const double p  = n1/n2 ;
  const double q  = 1 - p ;
  //
  const double eff = p ;
  const double c2  = kappa2 * p * q / n2 ;
  //
  return Ostap::Math::ValueWithError  ( eff , c2 ) ;
}
// ============================================================================
/*  Simple evaluation of efficiency from statistically independend
 *  "exclusive" samples "accepted" and "rejected"
 *  \f$ \varepsilon = \frac{1}{ 1 + \frac{N_{rejected}}{N_accepted}}\f$ 
 *  @param accepted  (IN) accepted sample 
 *  @param rejected  (IN) rejected sample 
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::exclusiveEff
( const Ostap::Math::ValueWithError& accepted , 
  const Ostap::Math::ValueWithError& rejected )
{ return Ostap::Math::binomEff2 ( accepted , rejected ) ; }
// ============================================================================
/*  evaluate the binomial efficiency for Bernulli scheme with weights 
 *  @param nAccepted (INPUT) number of accepted (weighted) events 
 *  @param nRejected (INPUT) number of rejected (weighted) events 
 *  @return the binomial efficiency 
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::binomEff2
( const ValueWithError& nAccepted , 
  const ValueWithError& nRejected ) 
{
  const double vA = nAccepted.value() ;
  const double vR = nRejected.value() ;
  //
  const bool zeroA = _zero ( vA      ) ;
  const bool zeroR = _zero ( vR      ) ;
  //
  if ( zeroA && zeroR ) { return ValueWithError ( 1 , -1 ) ; }
  //
  const double vB  = vA + vR ;
  const bool zeroB = _zero ( vB ) ;
  //
  if ( zeroB          ) { return ValueWithError ( 0 , -1 ) ; }
  //
  double cov2   =  vA * vA * nRejected.cov2() ;
  cov2         +=  vR * vR * nAccepted.cov2() ;
  cov2         /=  vB * vB   ;
  //
  return ValueWithError ( vA / vB , cov2 ) ;
}
// ============================================================================
/*  calculate the ratio of weighted to unweighted sample with uncertainties
 *  \f[ R = \frac{N_w}{N}  = \frac{ \sum_1^{N} w_i }{N} \f] 
 *  using jackknife method:
 *  \f[ \sigma^2(R) = \left( \sum_1^N w_i^2 - NR^2 \right) / (N-1)^2 \f] 
 *  - thanks to Wouter Hulsbergen 
 *  @see http://en.wikipedia.org/wiki/Jackknife_%28statistics%29
 *  The result has proper behaviour : 
 *  uncertainty in R goes to zero if 
 *  dispersion on weights go to zero.
 *  @param   nWeighted (input) statistic of weighted sample 
 *  @param   n         (input) size      of origial sample 
 *  @return  ratio R with the proper uncertaities 
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::effJackknife 
( const ValueWithError& nWeighted , 
  const unsigned long   n         ) 
{
  //
  if      ( 0 == n ) { return ValueWithError (-1,-1)                   ; }
  else if ( 1 == n ) { return ValueWithError ( nWeighted.value() , 0 ) ; }
  //
  const unsigned long n1 = n - 1 ;
  //
  const double r  = nWeighted.value() / n ;
  //
  double c2 = nWeighted.cov2 () - r*r*n ;
  //
  c2 /= n1 ;
  c2 /= n1 ;
  //
  return ValueWithError ( r , c2 ) ;
}
// ============================================================================
/*  Simple evaluation of efficiency using Zech's prescription 
 *  @param accepted  (IN) accepted sub-sample 
 *  @param total     (IN) total     sample 
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::zechEff
( const Ostap::Math::ValueWithError& accepted , 
  const Ostap::Math::ValueWithError& total    ) 
{
  //
  const double e   =        accepted.value () / total.value () ;
  const double v2  = total.value() * total.value() ;
  const double t1  =           total.cov2  () / v2 ;
  const double t2  =        accepted.cov2  () / v2 ;
  //
  const double c2  = e * e * t1  + ( 1 - 2 * e ) * t2 ;
  //
  return Ostap::Math::ValueWithError ( e , c2 ) ;
}
// ============================================================================
/*  evaluate pow(a,b)
 *  @param a (INPUT) the base
 *  @param b (INPUT) the exponent
 *  @return the <c>a</c> rased to power <c>b</b>
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::pow
( const Ostap::Math::ValueWithError& a ,
  const int                          b )
{
  //
  if      ( 0 == b         ) { return 1 ; }          // RETURN
  else if ( 1 == b         ) { return a ; }          // RETURN
  //
  else if ( 0 >= a.cov2 () || _zero ( a.cov2() ) )
  { return std::pow ( a.value() , b ) ;  }               // RETURN
  //
  const double v  =     std::pow ( a.value () , b     ) ;
  const double e1 = b * std::pow ( a.value () , b - 1 ) ;
  //
  return Ostap::Math::ValueWithError ( v , e1 * e1 * a.cov2 () ) ;
  //
}
// ============================================================================
/*  evaluate pow(a,b)
 *  @param a (INPUT) the base
 *  @param b (INPUT) the exponent
 *  @return the <c>a</c> raised to power <c>b</b>
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::pow
( const Ostap::Math::ValueWithError& a ,
  const double                       b )
{
  //
  if      ( _zero ( b )    ) { return 1 ; }         // RETURN
  else if ( _one  ( b )    ) { return a ; }         // RETURN
  else if ( 0 >= a.cov2 () || _zero ( a.cov2() ) )
  { return std::pow ( a.value() , b ) ; }           // RETURN
  //
  const double v  =     std::pow ( a.value () , b     ) ;
  const double e1 = b * std::pow ( a.value () , b - 1 ) ;
  //
  return Ostap::Math::ValueWithError ( v , e1 * e1 * a.cov2 () ) ;
}
// ============================================================================
/*  evaluate pow(a,b)
 *  @param a (INPUT) the base
 *  @param b (INPUT) the exponent
 *  @return the <c>a</c> raised to power <c>b</b>
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::pow
( const int                          a ,
  const Ostap::Math::ValueWithError& b )
{
  if      ( 0 == a && 0 < b.value()       ) { return 0 ; }    // RETURN
  else if ( 1 == a && _zero ( b.cov2 () ) ) { return 1 ; }    // RETURN
  else if ( 0 >= b.cov2() || _zero ( b.cov2() ) )
  { return std::pow ( double ( a ) , b.value() ) ; }    // RETURN
  //
  const double v  =     std::pow ( double ( a ) , b.value() ) ;
  const double e2 = v * std::log ( double ( a ) ) ;
  //
  return Ostap::Math::ValueWithError ( v , e2 * e2 * b.cov2 () ) ;
}
// ============================================================================
/*  evaluate pow(a,b)
 *  @param a (INPUT) the base
 *  @param b (INPUT) the exponent
 *  @return the <c>a</c> raised to power <c>b</b>
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::pow
( const double                       a ,
  const Ostap::Math::ValueWithError& b )
{
  if      ( _zero ( a ) && 0 < b.value() ) { return 0 ; }    // RETURN
  else if ( 0 >= b.cov2() || _zero ( b.cov2() ) )
  { return std::pow ( a , b.value() ) ; }    // RETURN
  //
  const double v  =     std::pow ( a , b.value() ) ;
  const double e2 = v * std::log ( a ) ;
  //
  return Ostap::Math::ValueWithError ( v , e2 * e2 * b.cov2 () ) ;
  //
}
// ============================================================================
/*  evaluate pow(a,b)
 *  @param a (INPUT) the base
 *  @param b (INPUT) the exponent
 *  @return the <c>a</c> raised to power <c>b</b>
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::pow
( const Ostap::Math::ValueWithError& a ,
  const Ostap::Math::ValueWithError& b )
{
  //
  if ( &a == &b ) 
  {
    if      ( 0 >= a.cov2 () || _zero ( a.cov2() ) )
    { return std::pow ( a.value() , a.value() ) ; }
    //
    const double v2 = std::pow ( a.value() , a.value() ) ;
    const double v3 = std::log ( a.value() ) + 1 ;
    //
    return Ostap::Math::ValueWithError
      ( v2 , v2 * v2 * v3 * v3 * a.cov2 () ) ;
  }
  //
  if      ( 0 >= a.cov2 () || _zero ( a.cov2() ) )
  { return pow ( a.value() , b         ) ; }
  else if ( 0 >= b.cov2 () || _zero ( b.cov2() ) )
  { return pow ( a         , b.value() ) ; }
  //
  const double v  = std::pow ( a.value () , b.value ()     ) ;
  const double v1 = std::pow ( a.value () , b.value () - 1 ) ;
  //
  const double e1 = v1 *            b.value ()   ;
  const double e2 = v  * std::log ( a.value () ) ;
  //
  return Ostap::Math::ValueWithError
    ( v , e1 * e1 * a.cov2 () + e2 * e2 * b.cov2 () ) ;
}
// ============================================================================
/*  evaluate exp(b)
 *  @param b (INPUT) the exponent
 *  @return the <c>e</c> raised to power <c>b</b>
 *  @warning invalid and small covariances are ignored
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::exp
( const Ostap::Math::ValueWithError& b )
{
  if ( 0 >= b.cov2 () || _zero ( b.cov2() ) )
  { return std::exp ( b.value() ) ; }
  //
  const double v = std::exp ( b.value() ) ;
  return Ostap::Math::ValueWithError ( v , v * v * b.cov2 () ) ;
}
// ============================================================================
/*  evaluate exp2(b)
 *  @param b (INPUT) the exponent
 *  @return the <c>e</c> raised to power <c>b</b>
 *  @warning invalid and small covariances are ignored
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::exp2
( const Ostap::Math::ValueWithError& b )
{
  if ( 0 >= b.cov2 () || _zero ( b.cov2() ) )
  { return std::exp2 ( b.value() ) ; }
  //
  const double v = std::exp2 ( b.value() ) ;
  return Ostap::Math::ValueWithError ( v , v * v * b.cov2 () * s_ln2_sq ) ;
}
// ============================================================================
/*  evaluate expm1(b)
 *  @param b (INPUT) the exponent
 *  @return  expm1
 *  @warning invalid and small covariances are ignored
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::expm1
( const Ostap::Math::ValueWithError& b )
{
  if ( 0 >= b.cov2 () || _zero ( b.cov2() ) )
  { return std::expm1 ( b.value() ) ; }
  //
  const double v  = std::expm1 ( b.value() ) ;
  const double d1 = ( v + 1 )  ;
  const double d2 = d1 * d1    ;
  const double e2 = d2 * b.cov2() ;
  //
  return Ostap::Math::ValueWithError ( v , e2 ) ;
}
// ============================================================================
/*  evaluate tgamma(b)
 *  @param b (INPUT) the exponent
 *  @return  tgamma
 *  @warning invalid and small covariances are ignored
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::tgamma
( const Ostap::Math::ValueWithError& b )
{
  if ( 0 >= b.cov2 () || _zero ( b.cov2() ) )
  { return std::tgamma ( b.value() ) ; }
  //
  const double bv = b.value() ;
  const double v  =      std::tgamma ( bv ) ;
  //
  // Gamma'/Gamma:
  const double p  = Ostap::Math::psi ( bv ) ;
  const double e1 = v * p * b.error() ;
  //
  return Ostap::Math::ValueWithError ( v , e1 * e1  ) ;
}
// ============================================================================
/*  evaluate lgamma(b)
 *  @param b (INPUT) the exponent
 *  @return  lgamma
 *  @warning invalid and small covariances are ignored
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::lgamma
( const Ostap::Math::ValueWithError& b )
{
  if ( 0 >= b.cov2 () || _zero ( b.cov2() ) )
  { return std::lgamma ( b.value() ) ; }
  //
  const double bv = b.value() ;
  const double v  = std::lgamma  ( bv ) ;
  //
  const double d1 = Ostap::Math::psi ( bv ) ;
  const double d2 = d1 * d1 ;
  const double e2 = d2 * b.cov2() ;
  //
  return Ostap::Math::ValueWithError ( v , e2 ) ;
}
// ============================================================================
/*  evaluate igamma(b)
 *  @param b (INPUT) the exponent
 *  @return  1/Gamma(b)
 *  @warning invalid and small covariances are ignored
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::igamma
( const Ostap::Math::ValueWithError& b )
{
  if ( 0 >= b.cov2 () || _zero ( b.cov2() ) )
  { return Ostap::Math::igamma ( b.value() ) ; }
  //
  const double bv = b.value() ;
  const double v  = Ostap::Math::igamma  ( bv ) ;
  //
  const double d1 = - Ostap::Math::psi( bv ) * v ;
  const double d2 = d1 * d1 ;
  const double e2 = d2 * b.cov2() ;
  //
  return Ostap::Math::ValueWithError ( v , e2 ) ;
}
// ============================================================================
/*  evaluate Pochhammer symbol 
 *  \f[ (x)^n = x ( x + 1) ( x + 1 ) ... ( x + n - 1 ) = \Pi^{k-1}_{k=0} (x + k) \f] 
 *  @see https://en.wikipedia.org/wiki/Falling_and_rising_factorials
 *  @param x (INPUT) the parameter 
 *  @param n (INPUT) the parameter 
 *  @return  pochhammer  symbol 
 *  @warning invalid and small covariances are ignored 
 *  @see Ostap::Math::rising_factorial
 *  @see Ostap::Math::falling_factorial
 *  @see Ostap::Math::pochhammer 
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::pochhammer 
( const Ostap::Math::ValueWithError& x , 
  const unsigned short               n )
{
  if      ( 0 == n )   { return 1 ; }  // simple case
  ///
  if ( 0 >= x.cov2 () || _zero ( x.cov2() ) ) 
  { return pochhammer ( x.value() , n ) ; }
  ///
  std::pair<double,double> r = pochhammer_with_derivative ( x , n ) ;
  //
  const double v  = r.first  ;
  const double d  = r.second ;
  const double e2 = x.cov2() * d * d ;
  //
  return Ostap::Math::ValueWithError ( v , e2 ) ;
}
// ============================================================================
/*  evaluate <code>hypot(x,y)</code>
 *  \f$ \sqrt( x^2 + y^2 ) \f$
 *   @param x (INPUT) the first parameter
 *   @param y (INPUT) the second parameter
 *   @param c (INPUT) the correlation coefficient  (-1<=c<=1)
 *   @return the valueof <code>hypot</code> function
 */
// ============================================================================
Ostap::Math::ValueWithError  Ostap::Math::hypot
( const Ostap::Math::ValueWithError& x , 
  const Ostap::Math::ValueWithError& y , 
  const double                       c ) 
{
  const bool x0 = 0 >= x.cov2() || _zero ( x.cov2() ) ;
  const bool y0 = 0 >= y.cov2() || _zero ( y.cov2() ) ;
  //
  const double r  = std::hypot ( x.value() , y.value() ) ;
  if ( x0 && y0 ) { return r ; }             // RETURN
  //
  double e2 = 0 ;
  if ( !_zero ( r ) ) 
  {
    e2 += x0 ? 0.0 : x.cov2() * x.value() * x.value() ;
    e2 += y0 ? 0.0 : y.cov2() * y.value() * y.value() ;
    e2 += x0 || y0 || _zero ( c ) ? 0.0 :  
      2 * std::max ( std::min ( c , 1.0 ) , -1.0 ) 
      * x.value() 
      * y.value() 
      * std::max ( 0.0 , x.error() ) 
      * std::max ( 0.0 , y.error() ) ;
    e2 /= r * r ;
  }
  else 
  {
    e2 += x0 ? 0.0 : x.cov2() ; 
    e2 += y0 ? 0.0 : y.cov2() ;    
    e2 += x0 || y0 || _zero ( c ) ? 0.0 :
      2 * std::max ( std::min ( c , 1.0 ) , -1.0 ) 
      * std::max ( 0.0 , x.error() ) 
      * std::max ( 0.0 , y.error() ) ;
  }
  //
  return ValueWithError ( r , e2 ) ;
}
// ============================================================================
/* evaluate fma(x,y,z) = x*y+x 
 *  @param y    (INPUT) the parameter 
 *  @param x    (INPUT) the parameter 
 *  @param z    (INPUT) the parameter 
 *  @param cxy  (INPUT) the correlation coefficient   -1<=c_xy<=1 
 *  @param cxz  (INPUT) the correlation coefficient   -1<=c_xz<=1 
 *  @param cyz  (INPUT) the correlation coefficient   -1<=c_yz<=1 
 *  @return  fma(x,y,z)
 *  @warning invalid and small covariances are ignored
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::fma 
( const Ostap::Math::ValueWithError& x   ,
  const Ostap::Math::ValueWithError& y   , 
  const Ostap::Math::ValueWithError& z   , 
  const double                       cxy ,
  const double                       cxz ,
  const double                       cyz ) 
{
  //
  const bool x0 = 0 >= x.cov2() || _zero ( x.cov2() ) ;
  const bool y0 = 0 >= y.cov2() || _zero ( y.cov2() ) ;
  const bool z0 = 0 >= z.cov2() || _zero ( z.cov2() ) ;
  //
  const double xv = x.value() ;
  const double yv = y.value() ;
  const double zv = z.value() ;
  //
  const double r = std::fma ( xv , yv , zv ) ;
  if ( x0 && y0 && z0 ) { return r ; }                              // RETURN 
  //
  
  double e2 = 0.0 ;
  e2 +=  x0 ? 0.0 : yv * yv * x.cov2 () ;
  e2 +=  y0 ? 0.0 : xv * xv * y.cov2 () ;
  e2 +=  z0 ? 0.0 :           z.cov2 () ;
  // correlations 
  e2 +=  x0 || y0 || s_zero ( cxy ) ? 0.0 : 
    2 * std::max ( std::min ( cxy , 1.0 ) , -1.0 ) 
    * yv * x.error () 
    * xv * y.error () ;
  e2 +=  x0 || z0 || s_zero ( cxz ) ? 0.0 : 
    2 * std::max ( std::min ( cxz , 1.0 ) , -1.0 ) 
    * yv * x.error () 
    *      z.error () ;
  e2 +=  y0 || z0 || s_zero ( cyz ) ? 0.0 : 
    2 * std::max ( std::min ( cyz , 1.0 ) , -1.0 ) 
    * xv * y.error () 
    *      z.error () ;
  //
  return Ostap::Math::ValueWithError ( r , e2 ) ;
}
// ============================================================================
/*  evaluate log(b)
 *  @param b (INPUT) the parameter
 *  @return logarithm
 *  @warning invalid and small covariances are ignored
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::log
( const Ostap::Math::ValueWithError& b )
{
  if ( 0 >= b.cov2 () || _zero ( b.cov2() ) ) { return std::log ( b.value() ) ; }
  //
  const double v  = std::log ( b.value () ) ;
  const double e1 = 1.0 /      b.value ()   ;
  //
  return Ostap::Math::ValueWithError ( v , e1 * e1 * b.cov2 () ) ;
}
// ============================================================================
/*  make a sum two elements taking into account the 
 *  correlation coefficient  
 *  @param a  (input) the first value 
 *  @param b  (input) the second value 
 *  @param c  (input) the correlation coefficient
 *  @return a+b 
 *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
 *  @date 2012-11-09
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::sum 
( const Ostap::Math::ValueWithError& a , 
  const Ostap::Math::ValueWithError& b , 
  const double                       c ) 
{
  //
  // the simplest case, ignore correlation  
  if ( &a == &b ) { return a*2.0 ; }     // RETURN
  //
  // few more trivial cases 
  //
  if      ( _zero ( a )  ) { return b ; }
  else if ( _zero ( b )  ) { return a ; }
  //
  // the second trivial case, no correlation  
  if      ( _zero ( c )  ) { return a + b ; } 
  //
  if      ( 0 > a.cov2() ) { return sum ( a.value() , b         , c ) ; }
  else if ( 0 > b.cov2() ) { return sum ( a         , b.value() , c ) ; }
  //
  const double v = a.value() + b.value();
  //
  // adjust the correlation coefficient 
  const double r   = std::max ( -1.0 , std::min ( 1.0 , c ) ) ;
  //
  const double ac2 = std::max ( a.cov2 () , 0.0 ) ;
  const double bc2 = std::max ( b.cov2 () , 0.0 ) ;
  //
  if ( _zero ( ac2 ) ) { return ValueWithError ( v , bc2 ) ; }  // RETURN
  if ( _zero ( bc2 ) ) { return ValueWithError ( v , ac2 ) ; }  // RETURN 
  //
  return ValueWithError ( v , ac2 + bc2 + 2 * r * std::sqrt ( ac2 * bc2 ) ) ;    
}
// ============================================================================
/*  make a sum two elements taking into account the 
 *  correlation coefficient  
 *  @param a  (input) the first value 
 *  @param b  (input) the second value 
 *  @param c  (input) the correlation coefficient
 *  @return a+b 
 *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
 *  @date 2012-11-09
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::sum2
( const Ostap::Math::ValueWithError& a , 
  const Ostap::Math::ValueWithError& b , 
  const double                       c ) 
{ return sum ( a , b , c ) ; }
// ============================================================================
/*  make a difference  two elements taking into acocunt the 
 *  correlation coefficient  
 *  @param a  (input) the first value 
 *  @param b  (input) the second value 
 *  @param c  (input) the correlation coefficient
 *  @return a-b 
 *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
 *  @date 2012-11-09
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::subtract
( const Ostap::Math::ValueWithError& a , 
  const Ostap::Math::ValueWithError& b , 
  const double                       c ) 
{
  //
  // the simplest case, ignore correlation  
  if ( &a == &b ) { return  ValueWithError ( 0 , 0 )  ; }     // RETURN
  //
  // few more trivial cases 
  //
  if      ( _zero ( a ) ) { return -b ; }
  else if ( _zero ( b ) ) { return  a ; }
  //
  const double v = a.value() - b.value();
  //
  if      ( 0 > a.cov2() ) { return subtract ( a.value() , b         , c ) ; }
  else if ( 0 > b.cov2() ) { return subtract ( a         , b.value() , c ) ; }
  //
  // the second trivial case, no correlation  
  if ( _zero ( c ) ) { return a - b ; } 
  //
  // adjust the correlation coefficient 
  const double r   = std::max ( -1.0 , std::min ( 1.0 , c ) ) ;
  //
  const double ac2 = std::max ( a.cov2 () , 0.0 ) ;
  const double bc2 = std::max ( b.cov2 () , 0.0 ) ;
  //
  if ( _zero ( ac2 ) ) { return ValueWithError ( v , bc2 ) ; }  // RETURN
  if ( _zero ( bc2 ) ) { return ValueWithError ( v , ac2 ) ; }  // RETURN 
  //
  return ValueWithError ( v , ac2 + bc2 - 2 * r * std::sqrt ( ac2 * bc2 ) ) ;    
}
// ============================================================================
/*  make a multiplication of two elements taking into acocunt the 
 *  correlation coefficient  
 *  @param a  (input) the first value 
 *  @param b  (input) the second value 
 *  @param c  (input) the correlation coefficient
 *  @return a*b 
 *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
 *  @date 2012-11-09
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::multiply
( const Ostap::Math::ValueWithError& a , 
  const Ostap::Math::ValueWithError& b , 
  const double                       c ) 
{
  //
  // the simplest case, ignore correlation  
  //
  if ( &a == &b ) 
  {
    const  double v = a.value() * a.value() ;
    return ValueWithError ( v , 4 * v * a.cov2()  ) ;   // RETURN
  }
  //
  // few more trivial cases 
  //
  if      ( _zero ( a ) ) { return ValueWithError ( 0 , 0 ) ; }
  else if ( _zero ( b ) ) { return ValueWithError ( 0 , 0 ) ; }
  else if ( _one  ( a ) ) { return b ; }
  else if ( _one  ( b ) ) { return a ; }
  //
  // ignore negative uncertainties 
  //
  if      ( 0 > a.cov2() ) { return multiply ( a.value() , b         , c ) ; }
  else if ( 0 > b.cov2() ) { return multiply ( a         , b.value() , c ) ; }
  //
  // the second trivial case, no correlation  
  //
  if ( _zero ( c ) ) { return a * b ; }                         // RETURN
  //
  const double  v   = a.value () * b.value () ;
  const double av2  = a.value () * a.value () ;
  const double bv2  = b.value () * b.value () ;
  //
  // adjust the correlation coefficient 
  const double r   = std::max ( -1.0 , std::min ( 1.0 , c ) ) ;
  //
  const double ac2 = std::max ( a.cov2 () , 0.0 ) ;
  const double bc2 = std::max ( b.cov2 () , 0.0 ) ;
  //
  if ( _zero ( ac2 ) ) { return ValueWithError ( v , av2 * bc2 ) ; }  // RETURN
  if ( _zero ( bc2 ) ) { return ValueWithError ( v , bv2 * ac2 ) ; }  // RETURN 
  //
  return ValueWithError ( v ,
                          bv2 * ac2 + 
                          av2 * bc2 + 
                          2 * v * r * std::sqrt ( ac2 * bc2 ) ) ;
}
// ============================================================================
/*  make a division of two elements taking into account the 
 *  correlation coefficient  
 *  @param a  (input) the first value 
 *  @param b  (input) the second value 
 *  @param c  (input) the correlation coefficient
 *  @return a/b 
 *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
 *  @date 2012-11-09
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::divide
( const Ostap::Math::ValueWithError& a , 
  const Ostap::Math::ValueWithError& b , 
  const double                       c ) 
{
  //
  // the simplest case, ignore correlation  
  //
  if ( &a == &b ) { return ValueWithError ( 1, 0 ) ; }
  //
  // few more trivial cases 
  //
  if      ( _zero ( a ) ) { return ValueWithError ( 0 , 0 ) ; }
  else if ( _one  ( a ) ) { return 1./b ; }
  else if ( _one  ( b ) ) { return a    ; }
  //
  // ignore negative uncertainties 
  //
  if      ( 0 > a.cov2() ) { return divide ( a.value() , b         , c ) ; }
  else if ( 0 > b.cov2() ) { return divide ( a         , b.value() , c ) ; }
  //
  // the second trivial case, no correlation  
  //
  if ( _zero ( c ) ) { return a / b ; }                         // RETURN
  //
  const double  v   = a.value () / b.value () ;
  const double av2  = a.value () * a.value () ;
  const double bv2  = b.value () * b.value () ;  
  //
  // adjust the correlation coefficient 
  const double r   = std::max ( -1.0 , std::min ( 1.0 , c ) ) ;
  //
  const double ac2 = std::max ( a.cov2 () , 0.0 ) ;
  const double bc2 = std::max ( b.cov2 () , 0.0 ) ;
  //
  const double ac2_n = ac2 / bv2 ;
  const double bc2_n = bc2 / bv2 ;
  const double av2_n = av2 / bv2 ;
  //
  if ( _zero ( ac2 ) ) { return ValueWithError ( v ,  av2_n * bc2_n ) ; }  // RETURN
  if ( _zero ( bc2 ) ) { return ValueWithError ( v ,          ac2_n ) ; }  // RETURN 
  //
  return ValueWithError ( v ,
                          ac2_n +  
                          av2_n * bc2_n - 
                          2 * v * r * std::sqrt ( ac2_n * bc2_n ) ) ;
}
// ===========================================================================
/*  calculate "fraction" of two elements (a/(a+b)) taking into account the 
 *  correlation coefficient  
 *  @param a  (input) the first value 
 *  @param b  (input) the second value 
 *  @param c  (input) the correlation coefficient
 *  @return a/(a+b)
 *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
 *  @date 2012-11-09
 */
// ===========================================================================
Ostap::Math::ValueWithError  
Ostap::Math::fraction
( const Ostap::Math::ValueWithError& a , 
  const Ostap::Math::ValueWithError& b , 
  const double                       c ) 
{
  //
  const double av = std::abs ( a.value() ) ;
  const double bv = std::abs ( b.value() ) ;
  return 
    av > bv ? 
    1.0       / ( 1.0 + divide ( b , a , c ) ) : 
    1.0 - 1.0 / ( 1.0 + divide ( a , b , c ) ) ;
}
// ===========================================================================
/* calculate "asymmetry" of two elements $\frac{a-b}{a+b}$
 *  taking into account the correlation coefficient  
 *  @param a  (input) the first value 
 *  @param b  (input) the second value 
 *  @param c  (input) the correlation coefficient
 *  @return (a-b)/(a+b)
 *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
 *  @date 2012-11-09
 */
// ===========================================================================
Ostap::Math::ValueWithError Ostap::Math::asymmetry
( const Ostap::Math::ValueWithError& a , 
  const Ostap::Math::ValueWithError& b , 
  const double                       c ) 
{
  //
  const double av = std::abs ( a.value() ) ;
  const double bv = std::abs ( b.value() ) ;
  if ( av > bv ) 
  {
    const ValueWithError d = divide ( b , a , c )  ;
    return divide ( 1.0 - d , 1.0 + d , -1.0 ) ;
  }
  const ValueWithError d = divide ( a , b , c )  ;
  return divide ( d - 1.0 , d + 1.0 , 1.0 ) ;
}
// ============================================================================
/*  evaluate log2(b)
 *  @param b (INPUT) the parameter
 *  @return logarithm
 *  @warning invalid and small covariances are ignored
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::log2
( const Ostap::Math::ValueWithError& b )
{
  if ( 0 >= b.cov2 () || _zero ( b.cov2() ) )
  { return std::log2 ( b.value() ) ; }
  //
  const double v  = std::log2 ( b.value() ) ;
  ///
  const double e1 = s_ln2_i / b.value() ;
  //
  return Ostap::Math::ValueWithError ( v , e1 * e1 * b.cov2 () ) ;
}
// ============================================================================
/*  evaluate log10(b)
 *  @param b (INPUT) the parameter
 *  @return logarithm
 *  @warning invalid and small covariances are ignored
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::log10
( const Ostap::Math::ValueWithError& b )
{
  if ( 0 >= b.cov2 () || _zero ( b.cov2() ) )
  { return std::log10 ( b.value() ) ; }
  //
  const double v  = std::log10 ( b.value() ) ;
  ///
  const double e1 = s_ln10_i / b.value() ;
  //
  return Ostap::Math::ValueWithError ( v , e1 * e1 * b.cov2 () ) ;
}
// ============================================================================
/*  evaluate log1p(b)
 *  @param b (INPUT) the parameter
 *  @return  log1p(b)
 *  @warning invalid and small covariances are ignored
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::log1p
( const Ostap::Math::ValueWithError& b )
{
  if ( 0 >= b.cov2 () || _zero ( b.cov2() ) )
  { return std::log1p ( b.value() ) ; }
  //
  const double bv = b.value() ;
  const double v  = std::log1p ( bv ) ;
  //
  const double d1 = 1 / ( 1 + bv ) ;
  const double e2 = d1 * d1 * b.cov2() ;
  //
  return Ostap::Math::ValueWithError ( v , e2 ) ;
}
// ============================================================================
/*  evaluate sqrt(b)
 *  @param b (INPUT) the parameter
 *  @warning invalid and small covariances are ignored
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::sqrt
( const Ostap::Math::ValueWithError& b )
{
  if ( 0 >= b.cov2 () || _zero ( b.cov2() ) )
  { return std::sqrt ( b.value() ) ; }
  //
  const double v  = std::sqrt ( b.value() ) ;
  ///
  const double e2 = 0.25 * b.cov2() / b.value() ;
  //
  return Ostap::Math::ValueWithError ( v , e2 ) ;
}
// ============================================================================
/*  evaluate "signed-sqrt" (a)
 *  @param a (INPUT) the value
 *  @return the signed-sqrt value
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::signed_sqrt
( const Ostap::Math::ValueWithError& a )
{ 
  return 
    0 >= a.cov2 () || _zero ( a.cov2() ) ? 
    Ostap::Math::ValueWithError ( Ostap::Math::signed_sqrt ( a.value() ) ) :
    0 <= a.value() ? 
    Ostap::Math::sqrt ( a ) :
    Ostap::Math::sqrt ( Ostap::Math::ValueWithError( -a.value() , a.cov2 () ) ) ;
}
// ============================================================================
/*  evaluate cbrt(b)
 *  @param b (INPUT) the parameter
 *  @warning invalid and small covariances are ignored
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::cbrt
( const Ostap::Math::ValueWithError& b )
{
  if ( 0 >= b.cov2 () || _zero ( b.cov2() ) )
  { return std::cbrt ( b.value() ) ; }
  //
  const double v  = std::cbrt ( b.value() ) ;
  //
  const double e2 = b.cov2() / ( v * b.value() ) ;
  //
  return Ostap::Math::ValueWithError ( v , e2 / 9.0 ) ;
}
// ============================================================================
/*  evaluate sin(b)
 *  @param b (INPUT) the parameter
 *  @warning invalid and small covariances are ignored
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::sin
( const Ostap::Math::ValueWithError& b )
{
  if ( 0 >= b.cov2 () || _zero ( b.cov2() ) )
  { return std::sin ( b.value() ) ; }
  //
  const double v  = std::sin ( b.value() ) ;
  const double d2 = std::max ( 1 - v*v , 0.0 ) ;
  //
  const double e2 = std::min ( d2 * b.cov2() , 1.0 ) ;
  //
  return Ostap::Math::ValueWithError ( v , e2 ) ;
}
// ============================================================================
/*  evaluate cos(b)
 *  @param b (INPUT) the parameter
 *  @warning invalid and small covariances are ignored
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::cos
( const Ostap::Math::ValueWithError& b )
{
  if ( 0 >= b.cov2 () || _zero ( b.cov2() ) )
  { return std::cos ( b.value() ) ; }
  //
  const double v  = std::cos ( b.value() ) ;
  const double d2 = std::max ( 1 - v*v , 0.0 ) ;
  //
  const double e2 = std::min ( d2 * b.cov2() , 1.0 ) ;
  //
  return Ostap::Math::ValueWithError ( v , e2 ) ;
}
// ============================================================================
/*  evaluate tan (b)
 *  @param b (INPUT) the parameter
 *  @warning invalid and small covariances are ignored
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::tan
( const Ostap::Math::ValueWithError& b )
{
  if ( 0 >= b.cov2 () || _zero ( b.cov2() ) )
  { return std::tan ( b.value() ) ; }
  //
  const double v  = std::tan ( b.value() ) ;
  const double d  = 1 + v * v ;
  //
  const double e2 = d * d  * b.cov2() ;
  //
  return Ostap::Math::ValueWithError ( v , e2 ) ;
}
// ============================================================================
/*  evaluate sinh(b)
 *  @param b (INPUT) the parameter
 *  @warning invalid and small covariances are ignored
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::sinh
( const Ostap::Math::ValueWithError& b )
{
  if ( 0 >= b.cov2 () || _zero ( b.cov2() ) )
  { return std::sinh ( b.value() ) ; }
  //
  const double v  = std::sinh ( b.value() ) ;
  const double d2 =  1 + v * v ;
  //
  const double e2 = d2 * b.cov2() ;
  //
  return Ostap::Math::ValueWithError ( v , e2 ) ;
}
// ============================================================================
/*  evaluate cosh(b)
 *  @param b (INPUT) the parameter
 *  @warning invalid and small covariances are ignored
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::cosh
( const Ostap::Math::ValueWithError& b )
{
  if ( 0 >= b.cov2 () || _zero ( b.cov2() ) )
  { return std::cosh ( b.value() ) ; }
  //
  const double v  = std::cosh ( b.value() ) ;
  const double d2 = v * v - 1 ;
  //
  const double e2 = std::max ( d2 * b.cov2() , 0.0 ) ;
  //
  return Ostap::Math::ValueWithError ( v , e2 ) ;
}
// ============================================================================
/*  evaluate tanh (b)
 *  @param b (INPUT) the parameter
 *  @warning invalid and small covariances are ignored
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::tanh
( const Ostap::Math::ValueWithError& b )
{
  if ( 0 >= b.cov2 () || _zero ( b.cov2() ) )
  { return std::tanh ( b.value() ) ; }
  //
  const double v  = std::tanh ( b.value() ) ;
  const double d  = 1 - v * v ;
  //
  const double e2 = std::min ( d * d  * b.cov2() , 1.0 ) ;
  //
  return Ostap::Math::ValueWithError ( v , e2 ) ;
}
// ============================================================================
/*  evaluate sech (b)
 *  @param b (INPUT) the parameter
 *  @warning invalid and small covariances are ignored
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::sech
( const Ostap::Math::ValueWithError& b )
{
  if ( 0 >= b.cov2 () || _zero ( b.cov2() ) )
  { return Ostap::Math::sech ( b.value() ) ; }
  //
  const double v  = Ostap::Math::sech ( b.value() ) ;
  const double d  = -v *    std::tanh ( b.value() ) ;
  //
  const double e2 = std::min ( d * d  * b.cov2() , 1.0 ) ;
  return Ostap::Math::ValueWithError ( v , e2 ) ;
}
// ============================================================================
/*  evaluate erf(b)
 *  @param b (INPUT) the parameter 
 *  @return  erf(b)
 *  @warning invalid and small covariances are ignored 
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::erf 
( const Ostap::Math::ValueWithError& b )
{
  if ( 0 >= b.cov2 () || _zero ( b.cov2() ) )
  { return std::erf( b.value() ) ; }
  //
  const double bv = b.value() ;
  const double v  = std::erf ( bv ) ;
  //
  static const double factor  = 4.0 / M_PI ;
  //
  const double d2 = factor * std::exp ( - bv * bv ) ;
  const double e2 = d2     * b.cov2() ;
  //
  return Ostap::Math::ValueWithError ( v , e2 ) ;
}
// ============================================================================
/*  evaluate erfc(b)
 *  @param b (INPUT) the parameter 
 *  @return  erfc(b)
 *  @warning invalid and small covariances are ignored 
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::erfc 
( const Ostap::Math::ValueWithError& b )
{
  if ( 0 >= b.cov2 () || _zero ( b.cov2() ) )
  { return std::erfc( b.value() ) ; }
  //
  const double bv = b.value() ;
  const double v  = std::erfc ( bv ) ;
  //
  static const double factor  = 4.0 / M_PI ;
  //
  const double d2 = factor * std::exp ( - bv * bv ) ;
  const double e2 = d2     * b.cov2() ;
  //
  return Ostap::Math::ValueWithError ( v , e2 ) ;
}
// ============================================================================
/*  evaluate erfi(b)
 *  @param b (INPUT) the parameter 
 *  @return  erfc(b)
 *  @warning invalid and small covariances are ignored 
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::erfi
( const Ostap::Math::ValueWithError& b )
{
  if ( 0 >= b.cov2 () || _zero ( b.cov2() ) )
  { return Ostap::Math::erfi( b.value() ) ; }
  //
  const double bv = b.value() ;
  const double v  = Ostap::Math::erfi ( bv ) ;
  //
  static const double factor  = 2.0 / std::sqrt ( M_PI ) ;
  //
  const double d2 = factor * std::exp ( bv * bv ) ;
  const double e2 = d2     * b.cov2() ;
  //
  return Ostap::Math::ValueWithError ( v , e2 ) ;
}
// ============================================================================
/*  evaluate erfcx(b)
 *  @param b (INPUT) the parameter 
 *  @return  erfc(b)
 *  @warning invalid and small covariances are ignored 
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::erfcx 
( const Ostap::Math::ValueWithError& b )
{
  if ( 0 >= b.cov2 () || _zero ( b.cov2() ) )
  { return Ostap::Math::erfcx ( b.value() ) ; }
  //
  const double bv = b.value() ;
  const double v  = Ostap::Math::erfcx ( bv ) ;
  //
  static const double factor  = 2.0 / std::sqrt ( M_PI ) ;
  //
  // derivative 
  const double d  = 2 * bv * v - factor ; //  
  const double e2 = d * d  * b.cov2()   ;
  //
  return Ostap::Math::ValueWithError ( v , e2 ) ;
}
// ============================================================================
/*  evaluate probit(b)
 *  @param b (INPUT) the parameter 
 *  @return  erfc(b)
 *  @warning invalid and small covariances are ignored 
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::probit 
( const Ostap::Math::ValueWithError& b )
{
  if ( 0 >= b.cov2 () || _zero ( b.cov2() ) )
  { return Ostap::Math::probit ( b.value() ) ; }
  //
  const double bv = b.value() ;
  const double v  = Ostap::Math::probit ( bv ) ;
  //
  static const double factor  = std::sqrt ( 2 * M_PI ) ;
  //
  // derivative 
  const double d  = factor * std::exp ( 0.5 * v * v );
  const double e2 = d * d  * b.cov2()   ;
  //
  return Ostap::Math::ValueWithError ( v , e2 ) ;
}
// ============================================================================
/*  evaluate asin(b)
 *  @param b (INPUT) the parameter 
 *  @return  asin(b)
 *  @warning invalid and small covariances are ignored 
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::asin
( const Ostap::Math::ValueWithError& b )
{
  if ( 0 >= b.cov2 () || _zero ( b.cov2() ) )
  { return std::asin( b.value() ) ; }
  //
  const double bv = b.value() ;
  const double v  = std::asin ( bv ) ;
  //
  const double b2 = bv * bv ;
  if ( _one  ( b2 ) ) { return Ostap::Math::ValueWithError ( v , -1 ) ; }
  //
  const double d2 = 1.0 / ( 1 - b2 ) ;
  const double e2 = d2     * b.cov2() ;
  //
  return Ostap::Math::ValueWithError ( v , e2 ) ;
}
// ============================================================================
/*  evaluate acos(b)
 *  @param b (INPUT) the parameter 
 *  @return  acos(b)
 *  @warning invalid and small covariances are ignored 
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::acos
( const Ostap::Math::ValueWithError& b )
{
  if ( 0 >= b.cov2 () || _zero ( b.cov2() ) )
  { return std::acos ( b.value() ) ; }
  //
  const double bv = b.value() ;
  const double v  = std::acos( bv ) ;
  //
  const double b2 = bv * bv ;
  if ( _one  ( b2 ) ) { return Ostap::Math::ValueWithError ( v , -1 ) ; }
  //
  const double d2 = 1.0 / ( 1 - b2 ) ;
  const double e2 = d2     * b.cov2() ;
  //
  return Ostap::Math::ValueWithError ( v , e2 ) ;
}
// ============================================================================
/*  evaluate atan(b)
 *  @param b (INPUT) the parameter 
 *  @return  atan(b)
 *  @warning invalid and small covariances are ignored 
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::atan
( const Ostap::Math::ValueWithError& b )
{
  if ( 0 >= b.cov2 () || _zero ( b.cov2() ) )
  { return std::atan ( b.value() ) ; }
  //
  const double bv = b.value() ;
  const double v  = std::atan( bv ) ;
  //
  const double b2 = bv * bv ;
  const double d1 = 1.0 / ( 1 + b2 ) ;
  const double d2 = d1 * d1 ;
  const double e2 = d2      * b.cov2() ;
  //
  return Ostap::Math::ValueWithError ( v , e2 ) ;
}
// ============================================================================
/*  evalute atan2(y,x)
 *  @param y    (INPUT) the parameter 
 *  @param x    (INPUT) the parameter 
 *  @param corr (INPUT) the correlation coefficient: -1<=corr<=1 
 *  @return  atan2(y,x)
 *  @warning invalid and small covariances are ignored 
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::atan2
( const Ostap::Math::ValueWithError& y    ,
  const Ostap::Math::ValueWithError& x    ,
  const double                       corr ) 
{
  //
  const double yv = y.value () ;
  const double xv = x.value () ;
  //
  const double v  = std::atan2 ( yv , xv ) ;
  //  
  const double y2 = yv * yv ;
  const double x2 = xv * xv ;
  const double r2 = x2 + y2 ;
  //
  const bool x_err =  0 < x.cov2() && !s_zero ( x.cov2() ) ;
  const bool y_err =  0 < y.cov2() && !s_zero ( y.cov2() ) ;
  
  if ( s_zero ( r2 ) ) 
  {
    // no reliable error estimates is possible
    if      ( x_err || y_err )
    { return  Ostap::Math::ValueWithError ( v , M_PI * M_PI ) ; }
    else  
    { return  Ostap::Math::ValueWithError ( v               ) ; }
  }  
  //
  const double      r4 = r2 * r2 ;
  const double dphidx2 = y2 / r4  ;
  const double dphidy2 = x2 / r4  ;
  //
  const double cor = std::max ( std::min ( corr , 1.0 ) , -1.0 ) ; 
  //
  const double e2 = 
    ( x_err ? dphidx2 * x.cov2() : 0.0 ) + 
    ( y_err ? dphidy2 * y.cov2() : 0.0 ) + 
    ( x_err && y_err && !s_zero ( cor ) ? 
      - 2 * cor * xv * yv * x.error() * y.error() / r4 : 0.0 ) ;
  //
  return Ostap::Math::ValueWithError ( v , e2 ) ;
}
// ============================================================================
/*  evaluate asinh(b)
 *  @param b (INPUT) the parameter 
 *  @return  asinh(b)
 *  @warning invalid and small covariances are ignored 
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::asinh
( const Ostap::Math::ValueWithError& b )
{
  if ( 0 >= b.cov2 () || _zero ( b.cov2() ) )
  { return std::asinh( b.value() ) ; }
  //
  const double bv = b.value() ;
  const double v  = std::asinh ( bv ) ;
  //
  const double b2 = bv * bv ;
  const double d2 = 1.0 / ( 1 + b2 ) ;
  const double e2 = d2     * b.cov2() ;
  //
  return Ostap::Math::ValueWithError ( v , e2 ) ;
}
// ============================================================================
/*  evaluate acosh(b)
 *  @param b (INPUT) the parameter 
 *  @return  acosh(b)
 *  @warning invalid and small covariances are ignored 
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::acosh
( const Ostap::Math::ValueWithError& b )
{
  if ( 0 >= b.cov2 () || _zero ( b.cov2() ) )
  { return std::acosh( b.value() ) ; }
  //
  const double bv = b.value() ;
  const double v  = std::acosh ( bv ) ;
  //
  const double b2 = bv * bv ;
  const double d2 = 1.0 / ( b2 - 1 ) ;
  const double e2 = d2     * b.cov2() ;
  //
  return Ostap::Math::ValueWithError ( v , e2 ) ;
}
// ============================================================================
/*  evaluate atanh(b)
 *  @param b (INPUT) the parameter 
 *  @return  acosh(b)
 *  @warning invalid and small covariances are ignored 
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::atanh
( const Ostap::Math::ValueWithError& b )
{
  if ( 0 >= b.cov2 () || _zero ( b.cov2() ) )
  { return std::atanh( b.value() ) ; }
  //
  const double bv = b.value() ;
  const double v  = std::atanh( bv ) ;
  //
  const double b2 = bv * bv ;
  const double d1 = 1.0 / ( 1 - b2 ) ;
  const double d2 = d1 * d1 ;
  const double e2 = d2      * b.cov2() ;
  //
  return Ostap::Math::ValueWithError ( v , e2 ) ;
}
// ============================================================================
/*  simple linear interpolation 
 *  @param x  the value to evaluate the function 
 *  @param x0 the abscissa for the first  point
 *  @param y0 the function value for the first  point
 *  @param x1 the abscissa for the second point
 *  @param y1 the function value for the second point
 *  @return linear interpolation at point x
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::interpolate_1D
( const double                       x  , 
  const double                       x0 ,
  const Ostap::Math::ValueWithError& y0 , 
  const double                       x1 ,
  const Ostap::Math::ValueWithError& y1 ) 
{
  //
  const double c0 = ( x - x1 ) / ( x0 - x1 ) ;
  const double c1 = ( x - x0 ) / ( x1 - x0 ) ;
  //
  return c0 * y0 + c1 * y1 ;
}
// ============================================================================
/*  simple (bi)linear interpolation 
 *  @param x  the x-coordiate to evaluate the function 
 *  @param y  the y-coordiate to evaluate the function 
 *  @param x0 the x-coordinate for the first  pair of points
 *  @param x1 the x-coordinate for the second pair of points
 *  @param y0 the y-coordinate for the first  pair of points
 *  @param y1 the y-coordinate for the second pair of points
 *  @param v00 the function value 
 *  @param v01 the function value 
 *  @param v10 the function value 
 *  @param v11 the function value 
 *  @return bilinear interpolation at point (x,y)
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::Math::interpolate_2D 
( const double                       x   , 
  const double                       y   , 
  const double                       x0  ,
  const double                       x1  ,
  const double                       y0  ,
  const double                       y1  ,
  const Ostap::Math::ValueWithError& v00 , 
  const Ostap::Math::ValueWithError& v01 , 
  const Ostap::Math::ValueWithError& v10 , 
  const Ostap::Math::ValueWithError& v11 ) 
{
  //
  const double c00 =  ( x - x1 ) * ( y - y1 ) / ( x0 - x1 ) / ( y0 - y1 ) ;
  const double c01 =  ( x - x1 ) * ( y - y0 ) / ( x0 - x1 ) / ( y1 - y0 ) ;
  const double c10 =  ( x - x0 ) * ( y - y1 ) / ( x1 - x0 ) / ( y0 - y1 ) ;
  const double c11 =  ( x - x0 ) * ( y - y0 ) / ( x1 - x0 ) / ( y1 - y0 ) ;
  //
  return c00 * v00 + c01 * v01 + c10 * v10 + c11 * v11  ;
} 
// ============================================================================
/*  get the sum of the vector 
 *  @param vct the vector
 *  @param ini the intial value 
 *  @return the sum over the vector 
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::sum 
( const std::vector<Ostap::Math::ValueWithError>& vct , 
  Ostap::Math::ValueWithError                     ini ) 
{
  //
  for ( std::vector<Ostap::Math::ValueWithError>::const_iterator iv = 
          vct.begin() ; vct.end() != iv ; ++iv ) { ini += (*iv) ; }
  //
  return ini ;
  //
}
// ============================================================================
/*  get the sum of the vector 
 *  @param vct the vector
 *  @param ini the intial value 
 *  @return the sum over the vector 
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::accumulate  
( const std::vector<Ostap::Math::ValueWithError>& vct , 
  Ostap::Math::ValueWithError                     ini ) 
{ return sum ( vct , ini ) ; }
// ============================================================================
/*  get the sum of absolute values for the vector 
 *  @param vct the vector
 *  @return the sum over the vector 
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::abssum 
( const std::vector<Ostap::Math::ValueWithError>& vct )
{
  //
  ValueWithError val ;
  for ( std::vector<Ostap::Math::ValueWithError>::const_iterator iv = 
          vct.begin() ; vct.end() != iv ; ++iv ) 
  { val += abs (*iv) ; }
  //
  return val ;
  //
}
// ============================================================================
/*  evaluate polynomial
 *  \f$f(x) = a_0 + a_1x + a_2x^2 + ... + a_{n-1}x^{n-1} + a_nx^n\f$
 *  such as \f$f(0) = a_0 \f$      
 *  using Horner rule
 *  @param poly  INPUT the coefficients
 *  @param x     INPUT argument 
 *  @return value of polynomial
 *  @warning invalid and small covariances are ignored 
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::horner_a0
( const std::vector<double>&         poly ,
  const Ostap::Math::ValueWithError& x    ) 
{
  if ( poly.empty() ) { return Ostap::Math::ValueWithError(0,0) ; }
  //
  const std::pair<long double,long double> r = 
    Ostap::Math::Clenshaw::monomial_sum 
    ( poly.rbegin() , poly.rend() , x.value() ) ;
  //
  if ( 0 >= x.cov2() || _zero( x.cov2() ) ) { return r.first  ; }
  //
  return Ostap::Math::ValueWithError( r.first , r.second * r.second * x.cov2() );
}
// ============================================================================
/*  evaluate polynomial
 *  \f$f(x) = a_0x^n + a_1x^{n-1}+ ... + a_{n-1}x + a_n\f$, 
 *   such as \f$f(0) = a_n \f$      
 *  using Horner rule
 *  @param poly  INPUT the coefficients 
 *  @param x     INPUT argument 
 *  @return value of polynomial
 *  @warning invalid and small covariances are ignored 
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::horner_aN
( const std::vector<double>&         poly ,
  const Ostap::Math::ValueWithError& x    ) 
{
  if ( poly.empty() ) { return Ostap::Math::ValueWithError(0,0) ; }
  //
  const std::pair<double,double> r = 
    Ostap::Math::Clenshaw::monomial_sum 
    ( poly.begin() , poly.end() , x.value() ) ;
  //
  if ( 0 >= x.cov2() || _zero( x.cov2() ) ) { return r.first  ; }
  return Ostap::Math::ValueWithError( r.first , r.second * r.second * x.cov2() ) ;
}
// ============================================================================
// Utiilties 
// ============================================================================
/*  simple interpolation 
 *  @param values     INPUT  vector of yi 
 *  @param abscissas  INPUT  vector of xi
 *  @param x          INPUT  the point where the function to be evaluated 
 *  @param correlated INPUT  correlated uncertaties in yi?
 */
// ===========================================================================
Ostap::Math::ValueWithError 
Ostap::Math::interpolate 
( const std::vector<Ostap::Math::ValueWithError>& y_i        ,
  const std::vector<double>&                      x_i        ,
  const double                                    x          , 
  const bool                                      correlated ) 
{
  // simple  cases 
  if      ( x_i.empty() || y_i.empty() ) { return Ostap::Math::ValueWithError() ; }
  else if ( 1 == x_i.size ()           ) { return y_i.front() ; }  
  //
  typedef Ostap::Math::ValueWithError VE ;
  if ( !correlated ) 
  { return Ostap::Math::Interpolation::lagrange
      ( x_i.begin() ,  
        x_i.end  () , 
        y_i.begin() ,  
        y_i.end  () , 
        x           , 
        Ostap::Math::ValueWithError()  , 
        [] ( double    x ) { return x ; } ,
        [] ( const VE& y ) { return y ; } ) ;
  }
  //
  std::vector<double> _y ( x_i.size() , 0.0 ) ;
  //
  //
  std::transform ( y_i.begin () , 
                   y_i.begin () +  std::min ( x_i.size() , y_i.size() )  , 
                   _y .begin () , [] ( const VE& y ) { return y.value() ; } ) ;                 
  const double r_0 = 
    Ostap::Math::Interpolation::neville 
    ( x_i.begin() ,  
      x_i.end  () , 
      _y .begin() ,  
      x           , 
      [] ( double x ) { return x         ; } ) ;
  //
  //
  std::transform ( y_i.begin () , 
                   y_i.begin () +  std::min ( x_i.size() , y_i.size() )  , 
                   _y .begin () , 
                   [] ( const VE& y ) { return y.value() + y.error() ; } );                 
  std::fill ( _y.begin() + y_i.size() , _y.end ( ) , 0.0 )  ;
  const double r_plus = 
    Ostap::Math::Interpolation::neville 
    ( x_i.begin() ,  
      x_i.end  () , 
      _y .begin() ,  
      x           , 
      [] ( double x ) { return x         ; } ) ;
  //
  //
  std::transform ( y_i.begin () , 
                   y_i.begin () +  std::min ( x_i.size() , y_i.size() )  , 
                   _y .begin () ,
                   [] ( const VE& y ) { return y.value() - y.error() ; } );                 
  std::fill ( _y.begin() + y_i.size() , _y.end ( ) , 0.0 )  ;
  const double r_minus = 
    Ostap::Math::Interpolation::neville 
    ( x_i.begin() ,  
      x_i.end  () , 
      _y .begin() ,  
      x           , 
      [] ( double x ) { return x         ; } ) ;
  //
  // get an estimate for the error 
  const double e = std::max ( std::abs ( r_plus  - r_0     ) ,
                              std::abs ( r_0     - r_minus ) ) ;
  return Ostap::Math::ValueWithError ( r_0 , e * e ) ;
}
// ============================================================================
/*  simple interpolation 
 *  - if vector of y is larger  than vector of x, extra values are ignored 
 *  - if vector of y is shorter than vector of x, missinge entries assumed to be zero 
 *  @param values     INPUT  vector of yi 
 *  @param abscissas  INPUT  vector of xi
 *  @param x          INPUT  the point where the function to be evaluated 
 */
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::interpolate 
( const std::vector<double>&         y_i , 
  const std::vector<double>&         x_i ,
  const Ostap::Math::ValueWithError& x   ) 
{
  //
  // simple  cases 
  if      ( x_i.empty() || y_i.empty() ) { return 0           ; }
  else if ( 1 == x_i.size ()           ) { return y_i.front() ; }  
  //
  std::vector<double> _y ( x_i.size() , 0.0 ) ;
  std::vector<double> _d ( x_i.size() , 0.0 ) ;
  //
  std::copy ( y_i.begin () , 
              y_i.begin () + std::min ( x_i.size() , y_i.size() ) , _y .begin () ) ;
  //
  std::pair<double,double> r = 
    Ostap::Math::Interpolation::neville2
    ( x_i.begin() ,  
      x_i.end  () , 
      _y .begin() ,  
      _d .begin() ,  
      x.value  () , 
      [] ( double x ) { return x  ; } ) ;
  //
  if (  0 >= x.cov2() || _zero ( x.cov2() ) ) { return r.first ; }
  return ValueWithError ( r.first , r.second * r.second * x.cov2() ) ;
}
// ========================================================================
/*  simple interpolation 
 *  - if vector of y is larger  than vector of x, extra values are ignored 
 *  - if vector of y is shorter than vector of x, missing entries assumed to be zero 
 *  @param values     INPUT  vector of yi 
 *  @param abscissas  INPUT  vector of xi
 *  @param x          INPUT  the point where the function to be evaluated 
 *  @param correlated INPUT  correlated uncertaties in yi?
 */
// ========================================================================
Ostap::Math::ValueWithError 
Ostap::Math::interpolate 
( const std::vector<Ostap::Math::ValueWithError>& y_i        , 
  const std::vector<double>&                      x_i        ,
  const Ostap::Math::ValueWithError&              x          , 
  const bool                                      correlated ) 
{
  // calculate value (with uncertainty with respect to y)
  ValueWithError rl  =  interpolate (  y_i , x_i , x.value() , correlated ) ;
  if ( 0 >= x.cov2() || _zero ( x.cov2() ) ){ return rl ; } // RETURN 
  //
  // calculate the uncertainty with respect to x 
  std::pair<double,double> rn = 
    Ostap::Math::Interpolation::neville2 
    ( x_i.begin() , x_i.end() , 
      y_i.begin() , y_i.end() , 
      x , 
      [] ( double x ) { return x  ; } , 
      [] ( double x ) { return x  ; } ) ;
  //
  const double d = rn.second ;
  rl.setCov2 ( std::max( rl.cov2() , 0.0 ) + x.cov2() * d * d  ) ;
  return rl ;  
}

// ============================================================================
/*  evaluate standard Gauss PDF 
 *  @param x the value 
 *  @return valeu of the standard gaussian PDF  
 */
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::gauss_pdf  ( const Ostap::Math::ValueWithError& x )
{
  if ( 0 >= x.cov2() || s_zero ( x.cov2() ) ) { return gauss_pdf ( x.value() ) ; }
  // the value  
  const double v =   gauss_pdf  ( x.value() ) ;
  // derivative  
  const double d = - v * x.value() ;
  // result 
  return Ostap::Math::ValueWithError ( v  , d  * d * x.cov2() ) ;
}
// ============================================================================
/* evaluate standard Gauss CDF 
 *  @param x the value 
 *  @return value of the standard gaussian CDF  
 */
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::gauss_cdf  ( const Ostap::Math::ValueWithError& x )
{
  if ( 0 >= x.cov2() || s_zero ( x.cov2() ) ) { return gauss_cdf ( x.value() ) ; }
  // the value  
  const double v = gauss_cdf ( x.value() ) ;
  // derivative  
  const double d = gauss_pdf ( x.value() )  ;
  // result 
  return Ostap::Math::ValueWithError ( v  , d  * d * x.cov2() ) ;
}
// ========================================================================    
 


// =============================================================================
// The END
// =============================================================================


