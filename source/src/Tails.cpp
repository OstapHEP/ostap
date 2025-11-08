// ============================================================================
// Include files 
// ============================================================================
// STD & STL
// ============================================================================
#include <cmath>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Hash.h"
#include "Ostap/Tails.h"
#include "Ostap/MoreMath.h"
#include "Ostap/StatusCode.h"
// ============================================================================
//  Local
// ============================================================================
#include "local_math.h"
#include "local_hash.h"
#include "status_codes.h" // the cache 
// ============================================================================
/** \f$ n \rightarrow N \f$ transformation
 *  @param  n (input) n-paeameter (external) 
 *  @return transformed N-parameter (internal)
 */
// ============================================================================
double Ostap::Math::Tail::N ( const double n )
{ return std::hypot ( 1.0 , n ) ; }
// ============================================================================
/*  Tail parameters 
 *  @param alpha alpha-parameter
 *  @param n     n-parameter
 */
// ============================================================================
Ostap::Math::Tail::Tail
( const double alpha  ,
  const double n      )
  : m_alpha ( std::abs ( alpha ) )
  , m_n     ( std::abs ( n     ) )
{
  Ostap::Assert ( m_alpha ,
                  "Invalid parameter `alpha` : must be non-zero!" ,
                  "Ostap::Math::Tail" ,
                  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
}
// ======================================================================
bool Ostap::Math::Tail::setAlpha  ( const double value )
{
  //
  const double value_ = std::abs ( value );
  if ( s_equal ( value_ , m_alpha ) ) { return false ; }
  m_alpha = value_  ;
  Ostap::Assert ( m_alpha ,
                  "Invalid parameter `alpha` : must be non-zero!" ,
                  "Ostap::Math::Tail" ,
                  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  return true ;
}
// ============================================================================
bool Ostap::Math::Tail::setN      ( const double value )
{
  double value_ = std::abs ( value  ) ;
  value_        = s_zero   ( value_ ) ? 0.0 : value_ ;
  if ( s_equal ( value_ , m_n     ) ) { return false ; }
  //
  m_n = value_ ;
  return true ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Tail::tag () const 
{
  static const std::string s_name = "Tail" ;
  return Ostap::Utils::hash_combiner ( s_name , m_alpha , m_n ) ; 
}
// ============================================================================

// ============================================================================
// Left Tail 
// ============================================================================
Ostap::Math::LeftTail::LeftTail
( const double alpha ,
  const double n     )
  : Ostap::Math::Tail ( alpha , n )
{}
// ============================================================================
Ostap::Math::LeftTail::LeftTail
( const Ostap::Math::Tail& tail ) 
  : Ostap::Math::Tail ( tail )
{}
// ============================================================================
/*  evaluate the (left) tail function
 *  @param x  x point 
 *  @param x0 normalization point 
 *  @param F    function value \f$ f(x_0) \ff at normalization point 
 *  @param dFoF value of log-dervative \f$ \frac{f^\prime(x_0)}{f(x_0)}\f$  at normalization point 
 */
// ============================================================================
double Ostap::Math::LeftTail::evaluate
( const double x      ,
  const double x0     ,
  const double F      ,
  const double dFoF   ) const
{
  if ( x0 < x || !F || s_zero ( F ) ) { return 0 ; }
  //}
  const double delta = x - x0 ;
  //
  Ostap::Assert ( 0 < F && 0 < dFoF ,
                  "Invalid parameter `dFoF`: log-derivative must be positive!" ,
                  "Ostap::Math::LeftTail::evaluate" ,
                  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  const double nn = N  () ;
  const double yy = nn / ( nn - dFoF * delta ) ;
  //
  return F * std::pow ( yy , nn ) ;
}
// ===========================================================================
/*  get the integral of power-law function
 *  @paral low  low  interal edge 
 *  @paral high high interal edge
 *  @param F    function value \f$ f(x_0) \f$ at normalization point 
 *  @param dFoF value of log-dervative \f$ \frac{f^\prime(x_0)}{f(x_0)}\f$  at normalization point 
 */
// ===========================================================================
double Ostap::Math::LeftTail::integral
( const double low  ,
  const double high ,
  const double x0   , 
  const double F    ,
  const double dFoF ) const
{
  //
  if      ( !F || s_zero ( F )     ) { return 0 ; }
  else if ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high <  low            ) { return - integral ( high , low , x0 , F , dFoF ) ; }
  else if ( x0   <= low            ) { return 0 ; }
  else if ( x0   <  high           ) { return   integral ( low ,  x0  , x0 , F , dFoF ) ; }
  //
  Ostap::Assert ( 0 < F && 0 < dFoF ,
                  "Invalid parameter `dFoF`: log-derivative must be positive!" ,
                  "Ostap::Math::LeftTail::integral" ,
                  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  // Cavalieri's integral:
  //
  const double nn = N () ;
  const double a  = -        dFoF / nn ;
  const double b  = 1 + x0 * dFoF / nn ;
  //
  return F * Ostap::Math::cavalieri ( -nn , low , high , a , b ) ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::LeftTail::tag () const 
{
  static const std::string s_name = "LefTail" ;
  return Ostap::Utils::hash_combiner ( s_name , Ostap::Math::Tail::tag() ) ; 
}
// ============================================================================

// ============================================================================
// RightTail 
// ============================================================================
Ostap::Math::RightTail::RightTail
( const double alpha ,
  const double n     )
  : Ostap::Math::Tail ( alpha , n )
{}
// ============================================================================
Ostap::Math::RightTail::RightTail
( const Ostap::Math::Tail& tail ) 
  : Ostap::Math::Tail ( tail )
{}
// ============================================================================
/*  evaluate the (right) tail function
 *  @param x  x point 
 *  @param x0 normalization point 
 *  @param F    function value \f$ f(x_0) \ff at normalization point 
 *  @param dFoF value of log-dervative \f$ \frac{f^\prime(x_0)}{f(x_0)}\f$  at normalization point 
 */
// ============================================================================
double Ostap::Math::RightTail::evaluate
( const double x      ,
  const double x0     ,
  const double F      ,
  const double dFoF   ) const

{
  if ( x < x0 || !F || s_zero ( F ) ) { return 0 ; }
  //}
  const double delta = x - x0 ;
  //
  Ostap::Assert ( 0 < F && 0 > dFoF ,
                  "Invalid parameter `dFoF`: log-derivative must be negative!" ,
                  "Ostap::Math::RightTail::evaluate" ,
                  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  const double nn = N  () ;
  const double yy = nn / ( nn - dFoF * delta ) ;
  //
  return F * std::pow ( yy , nn ) ;
}
// ============================================================================
/*  get the integral of power-law function
 *  @paral low  low  interal edge 
 *  @paral high high interal edge
 *  @param F    function value \f$ f(x_0) \f$ at normalization point 
 *  @param dFoF value of log-dervative \f$ \frac{f^\prime(x_0)}{f(x_0)}\f$  at normalization point 
 */
// ===========================================================================
double Ostap::Math::RightTail::integral
( const double low  ,
  const double high ,
  const double x0   , 
  const double F    ,
  const double dFoF ) const
{
  //
  if      ( !F || s_zero ( F )     ) { return 0 ; }
  else if ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high <  low            ) { return - integral ( high , low  , x0 , F , dFoF ) ; }
  else if ( high <= x0             ) { return 0 ; }
  else if ( low  <  x0             ) { return   integral ( x0   , high , x0 , F , dFoF ) ; }
  //
  Ostap::Assert ( 0 < F && 0 > dFoF ,
                  "Invalid parameter `dFoF`: log-derivative must be negative!" ,
                  "Ostap::Math::RightTail::integral" ,
                  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  // Cavalieri's integral:
  //
  const double nn = N () ;
  const double a  = -        dFoF / nn ;
  const double b  = 1 + x0 * dFoF / nn ;
  //
  return F * Ostap::Math::cavalieri ( -nn , low , high , a , b ) ;
}

// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::RightTail::tag () const 
{
  static const std::string s_name = "RightTail" ;
  return Ostap::Utils::hash_combiner ( s_name , Ostap::Math::Tail::tag() ) ; 
}
// ============================================================================


// ============================================================================
//                                                                      The END 
// ============================================================================

