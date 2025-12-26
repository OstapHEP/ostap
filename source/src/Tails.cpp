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
Ostap::Math::AlphaTail::AlphaTail
( const double alpha )
: m_alpha ( std::abs ( alpha ) )
{
  Ostap::Assert ( m_alpha ,
                  "Invalid parameter `alpha` : must be non-zero!" ,
                  "Ostap::Math::AlphaTail" ,
                  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
}
// ======================================================================
bool Ostap::Math::AlphaTail::setAlpha  ( const double value )
{
  //
  const double value_ = std::abs ( value );
  if ( s_equal ( value_ , m_alpha ) ) { return false ; }
  m_alpha = value_  ;
  Ostap::Assert ( m_alpha ,
                  "Invalid parameter `alpha` : must be non-zero!" ,
                  "Ostap::Math::AlphaTail" ,
                  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  return true ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::AlphaTail::tag () const 
{
  static const std::string s_name = "AlphaTail" ;
  return Ostap::Utils::hash_combiner ( s_name , m_alpha ) ; 
}
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
  : AlphaTail ( alpha ) 
  , m_n (                        std::abs ( n )   )
  , m_N ( Ostap::Math::Tail::N ( std::abs ( n ) ) )
{}
// ============================================================================
bool Ostap::Math::Tail::setN      ( const double value )
{
  double value_ = std::abs ( value  ) ;
  value_        = s_zero   ( value_ ) ? 0.0 : value_ ;
  if ( s_equal ( value_ , m_n     ) ) { return false ; }
  //
  m_n = value_    ;
  m_N = N ( m_n ) ;
  return true ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Tail::tag () const 
{
  static const std::string s_name = "Tail" ;
  return Ostap::Utils::hash_combiner ( s_name , Ostap::Math::AlphaTail::tag ()  , m_n ) ; 
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
// ===========================================================================
/*  get the integral of power-law function from negatove infinity to <code>high</code>
 *  @paral high high interal edge
 *  @param F    function value \f$ f(x_0) \f$ at normalization point 
 *  @param dFoF value of log-dervative \f$ \frac{f^\prime(x_0)}{f(x_0)}\f$  at normalization point 
 */
// ===========================================================================
double Ostap::Math::LeftTail::integral
( const double high ,
  const double x0   , 
  const double F    ,
  const double dFoF ) const
{
  //
  if      ( !F || s_zero ( F ) ) { return 0 ; }
  else if ( x0   <  high       ) { return integral ( x0  , x0 , F , dFoF ) ; }
  //
  Ostap::Assert ( 0 < F && 0 < dFoF ,
                  "Invalid parameter `dFoF`: log-derivative must be positive!" ,
                  "Ostap::Math::LeftTail::integral" ,
                  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  const double nn = N () ;
  /// infinite ?
  if ( nn < 1 || s_equal ( nn , 1.0 ) ) { return s_POSINF  ; }  // ATTENTION! 
  /// get the integral
  const double kappa = dFoF ; 
  const double beta  = nn - kappa * ( high - x0 ) ;
  //
  return - F * beta * std::pow ( nn / beta , nn ) / ( kappa * ( 1 - nn ) ) ;
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
/*  get the integral of power-law function from low to positive infinity 
 *  @paral low  low  interal edge 
 *  @param F    function value \f$ f(x_0) \f$ at normalization point 
 *  @param dFoF value of log-dervative \f$ \frac{f^\prime(x_0)}{f(x_0)}\f$  at normalization point 
 */
// ===========================================================================
double Ostap::Math::RightTail::integral
( const double low  ,
  const double x0   , 
  const double F    ,
  const double dFoF ) const
{
  //
  if      ( !F || s_zero ( F ) ) { return 0 ; }
  else if ( low  <  x0         ) { return integral ( x0 , x0 , F , dFoF ) ; }
  //
  Ostap::Assert ( 0 < F && 0 > dFoF ,
                  "Invalid parameter `dFoF`: log-derivative must be negative!" ,
                  "Ostap::Math::RightTail::integral" ,
                  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  // Cavalieri's integral:
  //
  const double nn = N () ;
  // infinite ?
  if ( nn < 1 || s_equal ( nn , 1.0 ) ) { return s_POSINF  ; }  // ATTENTION! 
  
  /// get the integral
  const double kappa = dFoF ; 
  const double beta  = nn - kappa * ( low - x0 ) ;
  //
  return F * beta * std::pow ( nn / beta , nn ) / ( kappa * ( 1 - nn ) ) ;
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
// ExpTail 
// ============================================================================
Ostap::Math::LeftExpTail::LeftExpTail
( const double alpha )
  : Ostap::Math::AlphaTail ( alpha )
{}
// ============================================================================
Ostap::Math::LeftExpTail::LeftExpTail
( const Ostap::Math::AlphaTail& tail ) 
  : Ostap::Math::AlphaTail ( tail )
{}

// ============================================================================
// ExpTail 
// ============================================================================
Ostap::Math::RightExpTail::RightExpTail
( const double alpha )
  : Ostap::Math::AlphaTail ( alpha )
{}
// ============================================================================
Ostap::Math::RightExpTail::RightExpTail
( const Ostap::Math::AlphaTail& tail ) 
  : Ostap::Math::AlphaTail ( tail )
{}

// ============================================================================
/*  evaluate the (left) tail function
 *  @param x  x point 
 *  @param x0 normalization point 
 *  @param F    function value \f$ f(x_0) \ff at normalization point 
 *  @param dFoF value of log-dervative \f$ \frac{f^\prime(x_0)}{f(x_0)}\f$  at normalization point 
 */
// ============================================================================
double Ostap::Math::LeftExpTail::evaluate
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
                  "Ostap::Math::LeftExpTail::evaluate" ,
                  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  const double kappa = dFoF ;
  //
  return F * std::exp ( kappa * delta  ) ;
}
// ============================================================================
/*  evaluate the (right) tail function
 *  @param x  x point 
 *  @param x0 normalization point 
 *  @param F    function value \f$ f(x_0) \ff at normalization point 
 *  @param dFoF value of log-dervative \f$ \frac{f^\prime(x_0)}{f(x_0)}\f$  at normalization point 
 */
// ============================================================================
double Ostap::Math::RightExpTail::evaluate
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
  const double kappa = dFoF ; 
  //
  return F * std::exp ( kappa * delta ) ;
}


// ===========================================================================
/*  get the integral of power-law function
 *  @paral low  low  interal edge 
 *  @paral high high interal edge
 *  @param F    function value \f$ f(x_0) \f$ at normalization point 
 *  @param dFoF value of log-dervative \f$ \frac{f^\prime(x_0)}{f(x_0)}\f$  at normalization point 
 */
// ===========================================================================
double Ostap::Math::LeftExpTail::integral
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
                  "Ostap::Math::LeftExpTail::integral" ,
                  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  const double kappa = dFoF ; 
  //
  const double B = std::exp ( kappa * ( high - x0 ) ) / kappa ; 
  const double A = std::exp ( kappa * ( low  - x0 ) ) / kappa ;  
  //
  return F * ( B - A ) ; 
}
// ============================================================================
/*  get the integral of power-law function
 *  @paral low  low  interal edge 
 *  @paral high high interal edge
 *  @param F    function value \f$ f(x_0) \f$ at normalization point 
 *  @param dFoF value of log-dervative \f$ \frac{f^\prime(x_0)}{f(x_0)}\f$  at normalization point 
 */
// ===========================================================================
double Ostap::Math::RightExpTail::integral
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
  // 
  const double kappa = dFoF ; 
  //
  const double B = std::exp ( kappa * ( high - x0 ) ) / kappa ; 
  const double A = std::exp ( kappa * ( low  - x0 ) ) / kappa ;  
  //
  return F * ( B - A ) ; 
  // 
}


// ===========================================================================
/*  get the integral of power-law function
 *  @paral high high interal edge
 *  @param F    function value \f$ f(x_0) \f$ at normalization point 
 *  @param dFoF value of log-dervative \f$ \frac{f^\prime(x_0)}{f(x_0)}\f$  at normalization point 
 */
// ===========================================================================
double Ostap::Math::LeftExpTail::integral
( const double high ,
  const double x0   , 
  const double F    ,
  const double dFoF ) const
{
  //
  if      ( !F || s_zero ( F )     ) { return 0 ; }
  else if ( high > x0              ) { return integral ( x0 , x0 , F , dFoF ) ; }
  //
  Ostap::Assert ( 0 < F && 0 < dFoF ,
                  "Invalid parameter `dFoF`: log-derivative must be positive!" ,
                  "Ostap::Math::LeftExpTail::integral" ,
                  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  const double kappa = dFoF ; 
  //
  const double B = std::exp ( kappa * ( high - x0 ) ) / kappa ; 
  //
  return F * B ; 
}
// ============================================================================
/*  get the integral of power-law function
 *  @paral low  low  interal edge 
 *  @param F    function value \f$ f(x_0) \f$ at normalization point 
 *  @param dFoF value of log-dervative \f$ \frac{f^\prime(x_0)}{f(x_0)}\f$  at normalization point 
 */
// ===========================================================================
double Ostap::Math::RightExpTail::integral
( const double low  ,
  const double x0   , 
  const double F    ,
  const double dFoF ) const
{
  //
  if      ( !F || s_zero ( F )     ) { return 0 ; }
  else if ( low < x0               ) { return integral ( x0 , x0 , F , dFoF ); }
  //
  Ostap::Assert ( 0 < F && 0 > dFoF ,
                  "Invalid parameter `dFoF`: log-derivative must be negative!" ,
                  "Ostap::Math::RightTail::integral" ,
                  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  // 
  const double kappa = dFoF ; 
  //
  const double A = std::exp ( kappa * ( low - x0 ) ) / kappa ;  
  //
  return -F * A ; 
  // 
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::LeftExpTail::tag () const 
{
  static const std::string s_name = "LeftExtTail" ;
  return Ostap::Utils::hash_combiner ( s_name , Ostap::Math::AlphaTail::tag() ) ; 
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::RightExpTail::tag () const 
{
  static const std::string s_name = "RightExtTail" ;
  return Ostap::Utils::hash_combiner ( s_name , Ostap::Math::AlphaTail::tag() ) ; 
}
//
/* get alpha-parameter for the (left) tail of Needham function
 *  @see Ostap::Math::Needham
 *  @param sigma (INOUT) sigma-parameter
 *  @param c0    (INPUT) c0-parameter
 *  @param c1    (INPUT) c1-parameterter 
 *  @param c2    (INPUT) c2-parameter
 *  @param amin  (INPUT) a_min parammeter
 */
// ============================================================================ 
double Ostap::Math::needham_alpha 
( const double sigma , 
  const double c0    , 
  const double c1    , 
  const double c2    ,
  const double amin  ) 
{
  Ostap::Assert ( 0 < c1 , 
                  "Invalid parameter c1: must be positive!" ,
                  "Ostap::Math::needham_alpha" , 
                  INVALID_PARAMETER , __FILE__ , __LINE__) ;

  const double sc1 = std::abs ( sigma / c1 ) ;
  /// avoid overflows (1)
  if ( 1 <= sc1 )
    {
      const double q = std::pow ( sc1 , c2 ) ;
      const double a = c0 * q / ( 1 + q ) ;
      return std::hypot ( amin , a ) ;
    }
  /// avoid overflows (2)
  const double Q = std::pow ( 1.0 / sc1 , c2 ) ;
  const double a = c0 / ( Q + 1 ) ;
  return std::hypot ( amin , a ) ; 
}

// ============================================================================
//                                                                      The END 
// ============================================================================

