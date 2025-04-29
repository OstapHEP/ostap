// ============================================================================
// Include files 
// ============================================================================
// STD & STL 
// ============================================================================
#include <limits>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Math.h"
#include "Ostap/Fourier.h"
#include "Ostap/Clenshaw.h"
#include "Ostap/Cesaro.h"
// ============================================================================
//  Local 
// ============================================================================
#include "Exception.h"
#include "local_math.h"
#include "status_codes.h"
// ============================================================================
/** @file
 *  Implementation file for functions from the file Ostap/Fourier.h
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2010-04-19
 */
// ============================================================================
namespace
{
  // ==========================================================================
  inline void
  adjust_minmax
  ( double& xmin ,
    double& xmax )
  {
    // ========================================================================
    if ( s_equal ( -M_PI     , xmin ) ) { xmin = -M_PI     ; }
    if ( s_equal ( -1        , xmin ) ) { xmin = -1        ; }
    if ( s_equal (  0        , xmin ) ) { xmin =  0        ; }
    // ========================================================================
    if ( s_equal (  1        , xmax ) ) { xmax =  1        ; }
    if ( s_equal (      M_PI , xmax ) ) { xmax =     M_PI  ; }
    if ( s_equal (  2 * M_PI , xmax ) ) { xmax =  2 * M_PI ; }
    // ========================================================================
  }
  // ==========================================================================
} //                                              The end of anynymos namespace 
// ============================================================================
// Constructor from the degree 
// ============================================================================
Ostap::Math::FourierSum::FourierSum 
( const unsigned short N      , 
  const double         xmin   , 
  const double         xmax   ) 
  : Ostap::Math::Parameters ( 2 * N + 1 )  // ATTENTION!! 
  , m_xmin  ( std::min ( xmin , xmax ) )
  , m_xmax  ( std::max ( xmin , xmax ) )
  , m_scale ( 1 ) 
  , m_delta ( 0 ) 
{
  //
  adjust_minmax ( m_xmin , m_xmax ) ;
  //
  Ostap::Assert ( m_xmin < m_xmax                      ,
                  "Invalid xmin/xmax setting!"         ,
                  "Ostap::Math::FourierSum"            ,
                  INVALID_MINMAX , __FILE__ , __LINE__ ) ;
  //
  m_scale = 2 * M_PI / ( m_xmax - m_xmin ) ;
  m_delta = 0.5      * ( m_xmax + m_xmin ) ;
}
// ============================================================================
// protected constructor from the parameters 
// ============================================================================
Ostap::Math::FourierSum::FourierSum
( const std::vector<double>& pars  , 
  const double               xmin  , 
  const double               xmax  ) 
  : Ostap::Math::Parameters ( pars ) 
  , m_xmin  ( std::min ( xmin , xmax ) )
  , m_xmax  ( std::max ( xmin , xmax ) )
  , m_scale ( 1 ) 
  , m_delta ( 0 ) 
{
  /// ensure we always have odd number of parameters! 
  if ( 0 == npars() % 2 ) { m_pars.push_back ( 0 ) ; } ;
  //
  adjust_minmax ( m_xmin , m_xmax ) ;
  //
  Ostap::Assert ( m_xmin < m_xmax                      ,
                  "Invalid xmin/xmax setting!"         ,
                  "Ostap::Math::FourierSum"            ,
                  INVALID_MINMAX , __FILE__ , __LINE__ ) ;
  //  //
  m_scale = 2 * M_PI / ( m_xmax - m_xmin ) ;
  m_delta = 0.5      * ( m_xmax + m_xmin ) ;
}
// ============================================================================
// swap 
// ============================================================================
void Ostap::Math::FourierSum::swap ( Ostap::Math::FourierSum&  right ) 
{
  Ostap::Math::Parameters::swap ( right ) ;
  std::swap ( m_xmin  ,  right.m_xmin  ) ;
  std::swap ( m_xmax  ,  right.m_xmax  ) ;
  std::swap ( m_scale ,  right.m_scale ) ;
  std::swap ( m_delta ,  right.m_delta ) ;
  std::swap ( m_aux   ,  right.m_aux   ) ;
}
// ============================================================================
/* get the magnitude of nth-harmonic
 * \f$m_k = \sqrt( a^2_k + b^2_k) \f$
 */
// ============================================================================
double Ostap::Math::FourierSum::mag    ( const unsigned short k ) const 
{
  if      ( k > N () ) { return 0 ; }
  else if ( 0 == k   ) { return std::abs ( m_pars[0] ) ; }
  //
  return std::hypot ( m_pars[ 2 * k - 1 ] , m_pars [ 2 * k ] ) ;
}
// ============================================================================
// get the phase of nth-harmonic
// ============================================================================
double Ostap::Math::FourierSum::phase ( const unsigned short k ) const 
{
  if      ( k > N () ) { return 0 ; }
  else if ( 0 == k   ) { return 0 <= m_pars [ 0 ] ? 0. : -M_PI ; }
  //
  return std::atan2 ( m_pars [ 2 * k - 1 ] , m_pars [ 2 * k ] ) ;
} 
// ============================================================================
// calculate Fourier sum 
// ============================================================================
double Ostap::Math::FourierSum::evaluate ( const double x ) const 
{
  /// transform to "t"-representation 
  const long double tv = t ( x ) ;
  return Ostap::Math::Clenshaw::fourier_sum ( begin() , end () , tv ) ;
}
// ============================================================================
// simple  manipulations with polynoms: scale it! 
// ============================================================================
Ostap::Math::FourierSum&
Ostap::Math::FourierSum::operator*=( const double a ) 
{
  Ostap::Math::scale ( m_pars , a ) ;
  return *this ;
}
// ============================================================================
// simple  manipulations with polynoms: scale it! 
// ============================================================================
Ostap::Math::FourierSum&
Ostap::Math::FourierSum::operator/=( const double a ) 
{
  Ostap::Math::scale ( m_pars , 1/a ) ;
  return *this ;
}
// ============================================================================
// simple  manipulations with polynoms: add constant 
// ===========================================================================
Ostap::Math::FourierSum&
Ostap::Math::FourierSum::operator+=( const double a ) 
{
  m_pars [ 0 ] += 2 * a ; // NOTE FACOR 2 HERE!  
  return *this ;
}
// ============================================================================
// simple  manipulations with polynoms: subtract constant 
// ============================================================================
Ostap::Math::FourierSum&
Ostap::Math::FourierSum::operator-=( const double a ) 
{
  m_pars [ 0 ] -= 2 * a ;
  return *this ;
}
// =============================================================================
/*  sum of two Fourier series (with the same interval!) 
 *  @param other the first fourier sum
 *  @return the sum of two Fourier series 
 */
// =============================================================================
Ostap::Math::FourierSum 
Ostap::Math::FourierSum::sum
( const Ostap::Math::FourierSum& other ) const 
{
  //
  if ( this == &other ) 
  {
    FourierSum result ( *this  ) ;
    result *= 2 ;
    return result ;
  }
  //
  if      ( other.zero() ) { return *this ; } // 
  else if (       zero() ) { return other ; } // random choice 
  //
  Ostap::Assert ( s_equal ( xmin () , other.xmin() ) && 
                  s_equal ( xmax () , other.xmax() )  , 
                  "Can't sum Fourier series with different domains" ,                  
                  "Ostap::Math::FourierSum"           ,
                  INVALID_RANGE , __FILE__ , __LINE__ ) ;
  //
  //
  const unsigned short n  = std::max ( N () , other.N  () ) ;
  //
  FourierSum result ( n , xmin() , xmax() ) ;
  const unsigned np  = result.npars() ;
  for ( unsigned short i = 0 ; i < np ; ++i ) 
  { result.m_pars  [ i ] =  par ( i )  + other.par( i ) ; }
  //
  return result ;
}
// =============================================================================
/*  get "shifted" fourier sum 
 *  \f$ g(x) \equiv f ( x - a ) \f$
 *  @param a the bias aprameter 
 *  @return the shifted fourier sum 
 */
// =============================================================================
Ostap::Math::FourierSum 
Ostap::Math::FourierSum::shift ( const double a ) const
{
  if ( s_zero ( a ) ) { return *this ; }
  //
  FourierSum result ( *this ) ;
  
  const unsigned long  N = m_pars.size() ;
  for ( unsigned short k = 1 ; 2 * k < N ; ++k  ) 
  {
    //  
    const double ct  = m_pars[2*k]   ; // cosine term 
    const double st  = m_pars[2*k-1] ; // sine   term 
    //
    const double ca  = std::cos ( k * a * m_scale ) ;
    const double sa  = std::sin ( k * a * m_scale ) ;
    //
    result.m_pars [ 2*k    ] = ct * ca - st * sa ;
    result.m_pars [ 2*k -1 ] = st * ca + ct * sa ;
  }
  //
  return result ;
}
// ============================================================================
// negate it!
// ============================================================================
Ostap::Math::FourierSum 
Ostap::Math::FourierSum::operator-() const
{
  FourierSum a ( *this );
  Ostap::Math::negate ( a.m_pars ) ;
  return a ;
}
// ============================================================================
// python operators    
// ============================================================================
Ostap::Math::FourierSum 
Ostap::Math::FourierSum::__add__  ( const double  value ) const 
{ return (*this) + value ; }
Ostap::Math::FourierSum 
Ostap::Math::FourierSum::__mul__  ( const double  value ) const 
{ return (*this) * value ; }
Ostap::Math::FourierSum 
Ostap::Math::FourierSum::__sub__  ( const double  value ) const 
{ return (*this) - value ; }
Ostap::Math::FourierSum 
Ostap::Math::FourierSum::__truediv__  ( const double  value ) const 
{ return (*this) / value ; }
Ostap::Math::FourierSum 
Ostap::Math::FourierSum::__radd__ ( const double  value ) const 
{ return (*this) + value ; }
Ostap::Math::FourierSum 
Ostap::Math::FourierSum::__rmul__ ( const double  value ) const 
{ return (*this) * value ; }
Ostap::Math::FourierSum
Ostap::Math::FourierSum::__rsub__ ( const double  value ) const 
{ return value   - (*this) ; }
Ostap::Math::FourierSum 
Ostap::Math::FourierSum::__add__  ( const Ostap::Math::FourierSum& b ) const 
{ return (*this) + b ; }
Ostap::Math::FourierSum 
Ostap::Math::FourierSum::__sub__  ( const Ostap::Math::FourierSum& b ) const 
{ return (*this) - b ; }
Ostap::Math::FourierSum 
Ostap::Math::FourierSum::__neg__  () const 
{ return -(*this) ;}
// ============================================================================


// ============================================================================
// constructor from the degree 
// ============================================================================
Ostap::Math::CosineSum::CosineSum 
( const unsigned short N      , 
  const double         xmin   , 
  const double         xmax   )  
  : Ostap::Math::Parameters ( N + 1 ) 
  , m_xmin  ( std::min ( xmin , xmax ) )
  , m_xmax  ( std::max ( xmin , xmax ) )
  , m_scale ( 1 ) 
{
  //
  adjust_minmax ( m_xmin , m_xmax ) ;
  //
  Ostap::Assert ( m_xmin < m_xmax                      ,
                  "Invalid xmin/xmax setting!"         ,
                  "Ostap::Math::CosineSum"             ,
                  INVALID_MINMAX , __FILE__ , __LINE__ ) ;
  //
  m_scale = M_PI / ( m_xmax - m_xmin ) ;
}
// ============================================================================
// constructor from non-empty list of parameters 
// ============================================================================
Ostap::Math::CosineSum::CosineSum
( const std::vector<double>& pars  , 
  const double               xmin  , 
  const double               xmax  ) 
  : Ostap::Math::Parameters ( pars )
  , m_xmin  ( std::min ( xmin , xmax ) )
  , m_xmax  ( std::max ( xmin , xmax ) )
  , m_scale ( 1 ) 
{
  //
  if ( !npars() ) { m_pars.push_back ( 0 ) ; }
  //
  adjust_minmax ( m_xmin , m_xmax ) ;
  //
  Ostap::Assert ( m_xmin < m_xmax                      ,
                  "Invalid xmin/xmax setting!"         ,
                  "Ostap::Math::CosineSum"             ,
                  INVALID_MINMAX , __FILE__ , __LINE__ ) ;
  //
  m_scale = M_PI / ( m_xmax - m_xmin ) ;
}
// ============================================================================
// swap 
// ============================================================================
void Ostap::Math::CosineSum::swap ( Ostap::Math::CosineSum&  right ) 
{
  Ostap::Math::Parameters::swap ( right ) ;
  std::swap ( m_xmin  ,  right.m_xmin  ) ;
  std::swap ( m_xmax  ,  right.m_xmax  ) ;
  std::swap ( m_scale ,  right.m_scale ) ;
  std::swap ( m_aux   ,  right.m_aux   ) ;
}
// ============================================================================
// calculate Fourier sum 
// ============================================================================
double Ostap::Math::CosineSum::evaluate ( const double x ) const 
{
  /// transform to "t"-representation 
  const long double tv = t(x) ;
  return Ostap::Math::Clenshaw::cosine_sum ( begin() , end ()  , tv ) ;
}
// ============================================================================
// simple  manipulations with polynoms: scale it! 
// ============================================================================
Ostap::Math::CosineSum&
Ostap::Math::CosineSum::operator*=( const double a ) 
{
  Ostap::Math::scale ( m_pars , a ) ;
  return *this ;
}
// ============================================================================
// simple  manipulations with polynoms: scale it! 
// ============================================================================
Ostap::Math::CosineSum&
Ostap::Math::CosineSum::operator/=( const double a ) 
{
  Ostap::Math::scale ( m_pars , 1/a ) ;
  return *this ;
}
// ============================================================================
// simple  manipulations with polynoms: add constant 
// ===========================================================================
Ostap::Math::CosineSum&
Ostap::Math::CosineSum::operator+=( const double a ) 
{
  m_pars [ 0 ] += 2 * a ; // ATTENTION! Note THE FACTOR 2 HERE 
  return *this ;
}
// ============================================================================
// simple  manipulations with polynoms: subtract constant 
// ============================================================================
Ostap::Math::CosineSum&
Ostap::Math::CosineSum::operator-=( const double a )
{  
  m_pars [ 0 ] -= 2 * a ; // ATTENTION, NOTE THE FACTOR 2 here! 
  return *this ;
}
// ============================================================================
/*  sum of two Fourier series (with the same interval!) 
 *  @param other the first fourier sum
 *  @return the sum of two Fourier series 
 */
// ============================================================================
Ostap::Math::CosineSum 
Ostap::Math::CosineSum::sum
( const Ostap::Math::CosineSum& other ) const 
{
  //
  if ( this == &other ) 
  {
    CosineSum result ( *this ) ;
    result *= 2 ;
    return result ;
  }
  //
  if      ( other.zero() ) { return *this ; }
  else if (       zero() ) { return other ; }
  //
  Ostap::Assert ( s_equal ( xmin () , other.xmin() ) && 
                  s_equal ( xmax () , other.xmax() ) , 
                  "Can't sum Cosine series with different domains" ,                  
                  "Ostap::Math::CosineSum"            ,
                  INVALID_RANGE , __FILE__ , __LINE__ ) ;
  
  //
  const unsigned short idegree = std::max ( degree () , other.degree () ) ;
  //
  CosineSum result ( idegree , xmin() , xmax() ) ;
  const unsigned npars  = result.npars() ;
  for ( unsigned short i = 0 ; i < npars ; ++i ) 
  { result.m_pars[i] =  par(i)  + other.par(i) ; }
  //
  return result ;
}
// ============================================================================
// negate it!
// ============================================================================
Ostap::Math::CosineSum 
Ostap::Math::CosineSum::operator-() const
{
  CosineSum a ( *this );
  Ostap::Math::negate ( a.m_pars ) ;
  return a ;
}
// ============================================================================
// python operators    
// ============================================================================
Ostap::Math::CosineSum 
Ostap::Math::CosineSum::__add__  ( const double  value ) const 
{ return (*this) + value ; }
Ostap::Math::CosineSum 
Ostap::Math::CosineSum::__mul__  ( const double  value ) const 
{ return (*this) * value ; }
Ostap::Math::CosineSum 
Ostap::Math::CosineSum::__sub__  ( const double  value ) const 
{ return (*this) - value ; }
Ostap::Math::CosineSum 
Ostap::Math::CosineSum::__truediv__  ( const double  value ) const 
{ return (*this) / value ; }
Ostap::Math::CosineSum 
Ostap::Math::CosineSum::__radd__ ( const double  value ) const 
{ return (*this) + value ; }
Ostap::Math::CosineSum 
Ostap::Math::CosineSum::__rmul__ ( const double  value ) const 
{ return (*this) * value ; }
Ostap::Math::CosineSum
Ostap::Math::CosineSum::__rsub__ ( const double  value ) const 
{ return value   - (*this) ; }
Ostap::Math::CosineSum 
Ostap::Math::CosineSum::__add__  ( const Ostap::Math::CosineSum& b ) const 
{ return (*this) + b ; }
Ostap::Math::CosineSum 
Ostap::Math::CosineSum::__sub__  ( const Ostap::Math::CosineSum& b ) const 
{ return (*this) - b ; }
Ostap::Math::CosineSum 
Ostap::Math::CosineSum::__neg__  () const 
{ return -(*this) ;}
// ============================================================================


// ============================================================================
/*  @param degree  degree
 *  @param xmin    low  edge
 *  @param xmax    high edge
 */
// ============================================================================
Ostap::Math::SineSum::SineSum 
( const unsigned short N      ,  // degree
  const double         xmin   ,  // low edge
  const double         xmax   ) // high edge
  : Ostap::Math::Parameters ( N  ) 
  , m_xmin  ( std::min ( xmin , xmax ) )
  , m_xmax  ( std::max ( xmin , xmax ) )
  , m_scale ( 1 ) 
{
  //
  adjust_minmax ( m_xmin , m_xmax ) ;
  //
  Ostap::Assert ( m_xmin < m_xmax                      ,
                  "Invalid xmin/xmax setting!"         ,
                  "Ostap::Math::SineSum"               ,
                  INVALID_MINMAX , __FILE__ , __LINE__ ) ;
  //
  m_scale = M_PI / ( m_xmax - m_xmin ) ;
}
// ============================================================================
// constructor from non-empty list of parameters
// ============================================================================
Ostap::Math::SineSum::SineSum 
( const std::vector<double>& pars  ,
  const double               xmin  ,
  const double               xmax  ) 
  : Ostap::Math::Parameters ( pars )
  , m_xmin  ( std::min ( xmin , xmax ) )
  , m_xmax  ( std::max ( xmin , xmax ) )
  , m_scale ( 1 ) 
{
  //
  adjust_minmax ( m_xmin , m_xmax ) ;
  //
  m_scale = M_PI / ( m_xmax - m_xmin ) ;
  // ===========================================================================  
}
// =============================================================================
// calculate Fourier sum 
// =============================================================================
double Ostap::Math::SineSum::evaluate  ( const double x ) const 
{
  /// transform to "t"-representation 
  const long double tv = t(x) ;
  return Ostap::Math::Clenshaw::sine_sum ( begin() , end () , tv ) ;
}
// ============================================================================
// simple  manipulations with polynoms: scale it! 
// ============================================================================
Ostap::Math::SineSum&
Ostap::Math::SineSum::operator*=( const double a ) 
{
  Ostap::Math::scale ( m_pars , a ) ;
  return *this ;
}
// ============================================================================
// simple  manipulations with polynoms: scale it! 
// ============================================================================
Ostap::Math::SineSum&
Ostap::Math::SineSum::operator/=( const double a ) 
{
  Ostap::Math::scale ( m_pars , 1/a ) ;
  return *this ;
}
// ============================================================================
/*  sum of two Fourier series (with the same interval!) 
 *  @param other the first fourier sum
 *  @return the sum of two Fourier series 
 */
// ============================================================================
Ostap::Math::SineSum 
Ostap::Math::SineSum::sum
( const Ostap::Math::SineSum& other ) const 
{
  //
  if ( this == &other ) 
  {
    SineSum result ( *this ) ;
    result *= 2 ;
    return result ;
  }
  //
  if      ( other.zero() ) { return *this ; }
  else if (       zero() ) { return other ; }
  //
  Ostap::Assert ( s_equal ( xmin () , other.xmin() ) && 
                  s_equal ( xmax () , other.xmax() ) , 
                  "Can't sum Sine series with different domains" ,                  
                  "Ostap::Math:SineSum"            ,
                  INVALID_RANGE , __FILE__ , __LINE__ ) ;
  
  //
  const unsigned short nn = std::max ( N () , other.N () ) ;
  //
  SineSum result (  nn  , xmin() , xmax() ) ;
  const unsigned npars  = result.npars() ;
  for ( unsigned short i = 0 ; i < npars ; ++i ) 
  { result.m_pars[i] =  par(i)  + other.par(i) ; }
  //
  return result ;
}
// ============================================================================
// negate it!
// ============================================================================
Ostap::Math::SineSum 
Ostap::Math::SineSum::operator-() const
{
  SineSum a ( *this );
  Ostap::Math::negate ( a.m_pars ) ;
  return a ;
}
// ============================================================================
// python operators    
// ============================================================================
Ostap::Math::SineSum 
Ostap::Math::SineSum::__mul__  ( const double  value ) const 
{ return (*this) * value ; }
Ostap::Math::SineSum 
Ostap::Math::SineSum::__truediv__  ( const double  value ) const 
{ return (*this) / value ; }
Ostap::Math::SineSum 
Ostap::Math::SineSum::__rmul__ ( const double  value ) const 
{ return (*this) * value ; }
Ostap::Math::SineSum 
Ostap::Math::SineSum::__neg__  () const 
{ return -(*this) ;}
// ============================================================================


// ============================================================================
//  Derivatives & integrals
// ============================================================================

// ============================================================================
//  get the derivative at point x 
// ============================================================================
double Ostap::Math::FourierSum::derivative ( const double x ) const 
{
  //
  m_aux.resize ( npars () ) ;
  m_aux [ 0 ] = 0 ;
  //
  const std::size_t NP { npars () } ;
  const std::size_t K  { N ()     } ;
  for ( std::size_t k = 1 ; k <= K ; ++k ) 
  {
    m_aux [ 2 * k - 1 ] = - m_pars [ 2 * k     ] * k ;
    m_aux [ 2 * k     ] =   m_pars [ 2 * k - 1 ] * k ;
  }
  /// transform to "t"-representation 
  const long double tv = t ( x ) ;
  /// make evaluation of fourier serie
  return m_scale * Ostap::Math::Clenshaw::fourier_sum ( m_aux.begin() , m_aux.end () , tv ) ;
}
// ===============================================================================
// get derivative as object 
// ===============================================================================
Ostap::Math::FourierSum
Ostap::Math::FourierSum::the_derivative () const
{
  FourierSum result { N() , xmin () , xmax () } ;
  //
  const std::size_t K { N () } ;
  for ( std::size_t k = 1 ; k <= K  ; ++k ) 
  {
    result.m_pars [ 2 * k - 1 ] = - m_pars [ 2 * k     ] * k ;
    result.m_pars [ 2 * k     ] =   m_pars [ 2 * k - 1 ] * k ;
  }
  result *= m_scale ;
  return result ;
}
// ============================================================================
// Get the derivative 
// ============================================================================
double Ostap::Math::CosineSum::derivative ( const double x ) const 
{
  //
  const std::size_t NP { npars ()  }  ; 
  m_aux.resize ( NP - 1 ) ;
  //
  for ( std::size_t k = 0 ; k + 1 < NP ; ++k ) { m_aux [ k ] = ( k + 1 )  * m_pars [ k + 1 ] ; }  
  /// transform to "t"-representation 
  const long double tv = t ( x ) ;
  /// make evaluation of sine serie
  return - m_scale * Ostap::Math::Clenshaw::sine_sum ( m_aux.begin () , m_aux.end () , tv ) ;
}
// ===============================================================================
// get derivative as object 
// ===============================================================================
Ostap::Math::SineSum
Ostap::Math::CosineSum::the_derivative () const
{
  SineSum result { N() , xmin () , xmax () } ;
  //
  PARAMETERS& r_pars = const_cast<PARAMETERS&> ( result.pars() ) ;
  //
  const std::size_t NP { npars ()  }  ; 
  for ( std::size_t k = 0 ; k + 1 < NP ; ++k ) { r_pars [ k ] = ( k + 1 )  * m_pars [ k + 1 ] ; }  
  //
  result *= -m_scale ;
  return result ;
}
// ============================================================================
// get the derivative at point x 
// ============================================================================
double Ostap::Math::SineSum::derivative ( const double x ) const 
{
  const std::size_t NP { npars () }  ; 
  m_aux.resize ( NP + 1 ) ;
  m_aux [ 0 ] = 0     ;
  for ( std::size_t k = 0 ; k < NP ; ++k ) { m_aux [ k + 1 ] = ( k + 1 ) * m_pars [ k ] ; }
  /// transform to "t"-representation 
  const long double tv = t ( x ) ;
  /// make evaluation of sine serie
  return m_scale * Ostap::Math::Clenshaw::cosine_sum ( m_aux.begin () , m_aux.end () , tv ) ;
}
// ===============================================================================
// get derivative as object 
// ===============================================================================
Ostap::Math::CosineSum
Ostap::Math::SineSum::the_derivative () const
{
  CosineSum result { N() , xmin () , xmax () } ;
  //
  PARAMETERS& r_pars = const_cast<PARAMETERS&> ( result.pars() ) ;
  //
  const std::size_t NP { npars ()  }  ;
  r_pars [ 0 ] = 0 ;
  for ( std::size_t k = 0 ; k < NP ; ++k ) { r_pars [ k + 1 ] = ( k + 1 ) * m_pars [ k ] ; }
  //
  result *= m_scale ;
  return result ;
}
// ============================================================================
// get the integral between low and high 
// ============================================================================
double Ostap::Math::FourierSum::integral 
( const double low  ,
  const double high ) const 
{
  if ( s_equal ( low , high ) ) { return 0 ; }
  //
  m_aux.resize ( npars () ) ;
  m_aux [ 0 ] = 0 ;
  //
  const std::size_t K { N () } ;
  for ( std::size_t k = 1 ; k <= K  ; ++k  ) 
  {
    m_aux [ 2 * k - 1 ] =   m_pars [ 2 * k     ] / k ;
    m_aux [ 2 * k     ] = - m_pars [ 2 * k - 1 ] / k ;
  }
  ///
  /// transform to "t"-representation 
  const long double tl = t ( low  ) ;
  const long double th = t ( high ) ;
  /// evaluate Fourier series
  return s_equal ( tl , th ) ? 0.0 :
    ( Ostap::Math::Clenshaw::fourier_sum ( m_aux.begin() , m_aux.end() , th ) - 
      Ostap::Math::Clenshaw::fourier_sum ( m_aux.begin() , m_aux.end() , tl ) ) / m_scale    
    + 0.5 * m_pars [ 0 ] * ( high - low ) ;
}
// ===============================================================================
/*   Get integral as function object           
 *   @attention The linear term p0/2*x is not included! 
 *   It needs to be added explicitely: 
 *   @code
 *   FourierSum fs = ... ; 
 *   FourierSum ci = fs.the_integral();
 *   x = 0.12 ;
 *   double result = ci ( x ) + 0.5 * cs[0] * cs.t(x) ; 
 *   @endcode
 */
// ===============================================================================
Ostap::Math::FourierSum
Ostap::Math::FourierSum::the_integral
( const double C ) const
{
  FourierSum result { N() , xmin () , xmax () } ;
  //
  result.m_pars [ 0 ] = 0 ;
  const std::size_t K { N () } ;
  for ( std::size_t k = 1 ; k <= K ; ++k  ) 
  {
    result.m_pars [ 2 * k - 1 ] =   m_pars [ 2 * k     ] / k ;
    result.m_pars [ 2 * k     ] = - m_pars [ 2 * k - 1 ] / k ;
  }
  ///
  result              /= m_scale ;
  result.m_pars [ 0 ] -= 2 * result ( x0 () ) ;
  result.m_pars [ 0 ] += 2 * C   ;
  return result ;
}
// ============================================================================
// get the integral between low and high 
// ============================================================================
double Ostap::Math::CosineSum::integral 
( const double low  ,
  const double high ) const 
{
  if ( s_equal ( low , high ) ) { return 0 ; }
  //
  const std::size_t NP { npars ()  }  ; 
  m_aux.resize ( NP - 1 ) ;
  //
  for ( std::size_t k = 0 ; k + 1 < NP ; ++k ) { m_aux [ k ] = m_pars [ k + 1 ] / ( k + 1 ) ; }
  //
  /// transform to "t"-representation 
  const long double tl = t ( low  ) ;
  const long double th = t ( high ) ;
  /// evaluate Fourier series
  return s_equal ( tl , th ) ? 0.0 : 
    ( Ostap::Math::Clenshaw::sine_sum ( m_aux.begin() , m_aux.end () , th ) - 
      Ostap::Math::Clenshaw::sine_sum ( m_aux.begin() , m_aux.end () , tl ) ) / m_scale
    + 0.5 * m_pars [ 0 ] * ( high - low ) ;
}
// ===============================================================================
/*   Get integral as function object           
 *   @attention The term p0/2*x is not included! 
 *   It needs to be added explicitely: 
 *   @code
 *   CosineSum cs = ... ; 
 *   SineSum   ci = cs.integral();
 *   x = 0.12 ;
 *   x0 = cs.x0() ; 
 *   double result = ci ( x ) + 0.5 * cs[0] * ( x - x0 ); 
 *   @endcode
 */
// ===============================================================================
Ostap::Math::SineSum
Ostap::Math::CosineSum::the_integral () const
{
  SineSum result { N () , xmin () , xmax () } ;
  //
  PARAMETERS& r_pars = const_cast<PARAMETERS&> ( result.pars() ) ;
  //  
  const std::size_t NP { npars ()  }  ; 
  for ( std::size_t k = 0 ; k + 1 < NP ; ++k )
  { r_pars [ k ] = m_pars [ k + 1 ] / ( k + 1 ) ; }
  //
  result /= m_scale ;
  return result ;
}
// ============================================================================
// get the integral between low and high 
// ============================================================================
double Ostap::Math::SineSum::integral 
( const double low  ,
  const double high ) const 
{
  if ( s_equal ( low , high ) ) { return 0 ; }
  //
  const std::size_t NP { npars ()  }  ; 
  m_aux.resize ( NP + 1 ) ;
  m_aux [ 0 ] = 0 ;
  //
  for ( std::size_t k = 0 ; k < NP ; ++k ) { m_aux [ k + 1 ] = m_pars [ k ] / ( k + 1)  ; }
  //
  /// transform to "t"-representation 
  const long double tl = t ( low  ) ;
  const long double th = t ( high ) ;
  /// evaluate Fourier series
  return s_equal ( tl , th ) ? 0.0 :     
    - ( Ostap::Math::Clenshaw::cosine_sum ( m_aux.begin() , m_aux.end () , th ) - 
        Ostap::Math::Clenshaw::cosine_sum ( m_aux.begin() , m_aux.end () , tl ) ) / m_scale ;
}
// ===============================================================================
Ostap::Math::CosineSum
Ostap::Math::SineSum::the_integral
( const double C ) const
{
  CosineSum result { N() , xmin () , xmax () } ;
  //
  PARAMETERS& r_pars = const_cast<PARAMETERS&> ( result.pars() ) ;
  r_pars [ 0 ] = 0 ;
  //  
  const std::size_t NP { npars ()  }  ;
  for ( std::size_t k = 0 ; k < NP ; ++k )
    { r_pars [ k + 1 ] = m_pars [ k ] / ( k + 1)  ; }
  //
  result       /= -m_scale ;
  r_pars [ 0 ] -= 2 * result ( x0()  ) ;
  r_pars [ 0 ] += 2 * C ;
  return result ;
}
// ============================================================================
Ostap::Math::FourierSum
Ostap::Math::FourierSum::convolute  
( const double sigma ) const 
{
  //
  // no convolution 
  if ( sigma <= 0 || s_zero ( sigma ) ) { return *this ; }
  //
  const long double ss      =  sigma / m_scale ;
  const long double sigma2  =  ss * ss         ;
  /// create convolution object 
  FourierSum conv( m_pars , m_xmin , m_xmax ) ;
  /// fill it! 
  conv.m_pars [ 0 ] = m_pars [ 0 ]  ;
  const unsigned long  N = m_pars.size() ;
  for ( unsigned short k = 1 ; 2 * k < N ; ++k  ) 
  {
    const long double s  = std::exp ( - 0.5L * k * k * sigma2 ) ;
    //
    const long double v1 = s * m_pars[ 2 * k    ] ;
    if ( v1 && !s_zero ( v1 ) ) { conv.m_pars [ 2 * k    ] = v1 ; }
    //
    const long double v2 = s * m_pars[ 2 * k -1 ] ;
    if ( v2 && !s_zero ( v2 ) ) { conv.m_pars [ 2 * k -1 ] = v2 ; }
    //
  }
  //
  return conv ;
}
// ============================================================================
Ostap::Math::CosineSum
Ostap::Math::CosineSum::convolute 
( const double sigma ) const 
{
  //
  // no convolution 
  if ( sigma <= 0 || s_zero ( sigma ) ) { return *this ; }
  //
  const long double ss      =  sigma / m_scale ;
  const long double sigma2  =  ss*ss           ;
  // create convolution object 
  CosineSum conv( m_pars , m_xmin , m_xmax ) ;
  /// fill it! 
  conv.m_pars [ 0 ] = m_pars [ 0 ]  ;
  const unsigned long  N = m_pars.size() ;
  for ( unsigned short k = 1 ; k < N ; ++k  ) 
  {
    const long double s  = std::exp ( - 0.5L * k * k * sigma2 ) ;
    //
    const long double v1 = s * m_pars   [ k ]      ;
    if ( v1 && !s_zero ( v1 ) ) { conv.m_pars [ k ] = v1 ; }
  }
  //
  return conv ;
}

// ============================================================================
Ostap::Math::SineSum
Ostap::Math::SineSum::convolute 
( const double sigma ) const 
{
  //
  // no convolution 
  if ( sigma <= 0 || s_zero ( sigma ) ) { return *this ; }
  //
  const long double ss      =  sigma / m_scale ;
  const long double sigma2  =  ss*ss           ;
  // create convolution object 
  SineSum conv( m_pars , m_xmin , m_xmax ) ;
  /// fill it! 
  const unsigned long  N = m_pars.size() ;
  for ( unsigned short k = 0 ; k < N ; ++k  ) 
  {
    const std::size_t kk = k + 1 ;
    const long double s  = std::exp ( - 0.5L * kk * kk * sigma2 ) ;
    //
    const long double v1 = s * m_pars [ k ] ;
    if ( v1 && !s_zero ( v1 ) ) { conv.m_pars [ k ] = v1 ; }
  }
  //
  return conv ;
}

// ============================================================================
/*  Calcuate the series using Cesaro method order k 
 *  @param   k  (INPUT)  th esummation order
 *  @return the sum that will be calcualetd usnig Cesaro method        
 */
// ============================================================================
Ostap::Math::FourierSum
Ostap::Math::FourierSum::cesaro
( const unsigned short k ) const
{
  if ( !k ) { return *this ; }
  //
  FourierSum result { *this } ; 
  // calculae the coefficients 
  Ostap::Math::cesaro ( k , begin() , end () , result.m_pars.begin() ) ;
  return result ;
}
// ============================================================================
/*  Calcuate the series using Cesaro method order k 
 *  @param   k  (INPUT)  th esummation order
 *  @return the sum that will be calcualetd usnig Cesaro method        
 */
// ============================================================================
Ostap::Math::CosineSum
Ostap::Math::CosineSum::cesaro
( const unsigned short k ) const
{
  if ( !k ) { return *this ; }
  //
  CosineSum result { *this } ;
  // calculae the coefficients 
  Ostap::Math::cesaro ( k , begin() , end () , result.m_pars.begin() ) ;
  return result ; 
}
// ============================================================================
/*  Calcuate the series using Cesaro method order k 
 *  @param   k  (INPUT)  th esummation order
 *  @return the sum that will be calcualetd usnig Cesaro method        
 */
// ============================================================================
Ostap::Math::SineSum
Ostap::Math::SineSum::cesaro
( const unsigned short k ) const
{
  if ( !k ) { return *this ; }
  //
  SineSum result { *this } ; 
  // calculae the coefficients 
  Ostap::Math::cesaro ( k , begin() , end () , result.m_pars.begin() ) ;
  return result ; 
}
// // ============================================================================
// //  convolute Fourier sum with gaussian function 
// // ============================================================================
// Ostap::Math::FourierSum 
// Ostap::Math::FourierSum::deconvolve 
// ( const double sigma , 
//   const double delta ) const 
// {
//   // no convolution 
//   if ( s_zero ( sigma ) ) { return *this ; }
//   //
//   const long double ss      =  sigma / m_scale ;
//   const long double sigma2  =  ss*ss           ;
//   // create covolution object 
//   FourierSum conv( m_pars , m_xmin , m_xmax ) ;
//   /// fill it! 
//   conv.m_pars [0] = m_pars[0]  ;
//   const unsigned long  N = m_pars.size() ;
//   //
//   const bool use_delta = !s_zero ( delta ) && 0 < delta ;
//   for ( unsigned short k = 1 ; 2 * k < N ; ++k  ) 
//   {
//     //  
//     const double v_1  = m_pars[2*k] ;
//     const double v_2  = m_pars[2*k-1] ;
//     if ( s_zero ( v_1 )  && s_zero ( v_2 ) ) { continue ; }
//     //
//     long double f = my_exp ( 0.5L * k * k * sigma2 ) ;
//     //
//     if ( use_delta ) 
//     { const long double fd = f * delta ; f /= ( 1 + fd * fd ) ; }
//     //
//     const long double   v1 = f * v_1 ;
//     if ( !s_zero ( v1 ) ) { conv.m_pars [ 2 * k    ] = v1 ; }
//     else { conv.m_pars[2*k  ] = 0 ; }    
//     //
//     const long double   v2 = f * v_2 ;
//     if ( !s_zero ( v2 ) ) { conv.m_pars [ 2 * k -1 ] = v2 ; }
//     else { conv.m_pars[2*k-1] = 0 ; }    
//     //
//   }
//   //
//   return conv ;
// }
// // ============================================================================
// /* get the effective cut-off (==number of effective harmonics) 
//  * of Tikhonov's regularization 
//  * \f$ n \equiv  \sqrt{2 \ln \delta} \frac{2\pi\sigma}{L} \f$
//  * @param sigma  gaussian resolution 
//  * @param delta  regularization parameter 
//  * @return number of effective harmonic 
//  */
// // ============================================================================
// double Ostap::Math::FourierSum::regularization 
// ( const double sigma , 
//   const double delta ) const 
// {
//   if      ( 0 > delta || s_zero ( delta ) || s_zero ( sigma ) ) 
//   { return s_UL_max ; } // return
//   else if ( 1<= delta ) { return 1 ; }
//   //
//   return std::sqrt ( -2 * std::log ( delta ) ) * m_scale / std::abs ( sigma ) ;
// }

// // ============================================================================
// //  convolute Fourier sum with gaussian function 
// // ============================================================================
// Ostap::Math::CosineSum 
// Ostap::Math::CosineSum::deconvolve 
// ( const double sigma , 
//   const double delta ) const 
// {
//   // no convolution 
//   if ( s_zero ( sigma ) ) { return *this ; }
//   //
//   const long double ss      =  sigma / m_scale ;
//   const long double sigma2  =  ss*ss           ;
//   // create covolution object 
//   CosineSum conv( m_pars , m_xmin , m_xmax ) ;
//   /// fill it! 
//   conv.m_pars [0] = m_pars[0]  ;
//   const unsigned long  N = m_pars.size() ;
//   const bool use_delta = !s_zero ( delta ) && 0 < delta ;
//   for ( unsigned short k = 1 ; k < N ; ++k  ) 
//   {
//     //  
//     const double v = m_pars[k] ;
//     if ( s_zero ( v ) ) { continue ; }
//     //
//     long double f = my_exp ( 0.5L * k * k * sigma2 ) ;
//     //
//     if ( use_delta ) 
//     { const long double fd = f * delta ; f /= ( 1 + fd * fd ) ; }
//     //
//     const long double   v1 = f * v ;
//     if ( !s_zero ( v1 ) ) { conv.m_pars [ k ] = v1 ; }
//     else { conv.m_pars[k] = 0 ; }    
//     //
//   }
//   //
//   return conv ;
// }
// // ============================================================================
// /* Get the effective cut-off (==number of terms/harmonics) 
//  * of Tikhonov's regularization 
//  * \f$ n \equiv  \sqrt{2 \ln \delta} \frac{\pi\sigma}{L} \f$
//  * @param sigma  gaussian resolution 
//  * @param delta  regularization parameter 
//  * @return number of effective harmonic 
//  */
// // ============================================================================
// double Ostap::Math::CosineSum::regularization 
// ( const double sigma , 
//   const double delta ) const 
// {
//   if      ( 0 > delta || s_zero ( delta ) || s_zero ( sigma ) ) { return s_UL_max ; } 
//   else if ( 1<= delta ) { return 1 ; }
//   //
//   return std::sqrt ( -2 * std::log ( delta ) ) * m_scale / std::abs ( sigma ) ;
// }

//     public:
//       // ======================================================================
//       /** convolute with gaussian
//        *  @param sigma resoltuion parameter for gaussian
//        *  @return convolution witgh gaussian
//        */
//       FourierSum convolve 
//       ( const double sigma     ) const ;
//       /** deconvolute with optional regularization
//        *  @param sigma sigma of gaussian
//        *  @param delta parameter of Tikhonov's regularization
//        *  for delta<=0, no regularization
//        *  @return regularised deconvolution
//        */
//       FourierSum deconvolve  
//       ( const double sigma     ,
//         const double delta = 0 ) const ;
//       /**  get the effective cut-off (==number of effective harmonics)
//        *   of Tikhonov's regularization
//        *   \f$ n \equiv  \sqrt{2 \ln \delta} \frac{\pi\sigma}{L} \f$
//        *   @param sigma  gaussian resoltuion
//        *   @param delta  regularization parameter
//        *   @return number of effective harmonic
//        */
//       double     regularization
//       ( const double sigma     ,
//         const double delta     ) const ;
//       // ======================================================================

//     public:
//       // ======================================================================
//       /** convolute with gaussian
//        *  @param sigma resoltuion parameter for gaussian
//        *  @return convolution witgh gaussian
//        */
//       CosineSum   convolve
//       ( const double sigma     ) const ;
//       /** deconvolute with optional regularization
//        *  @param sigma sigma of gaussian
//        *  @param delta parameter of Tikhonov's regularization
//        *  for delta<=0, no regularization
//        *  @return regularised deconvolution
//        */
//       CosineSum deconvolve
//       ( const double sigma     ,
//         const double delta = 0 ) const ;
//       /** get the effective cut-off (==number of terms/harmonics)
//        *  of Tikhonov's regularization
//        *  \f$ n \equiv  \sqrt{2 \ln \delta} \frac{\pi\sigma}{L} \f$
//        *  @param sigma  gaussian resoltuion
//        *  @param delta  regularization parameter
//        *  @return number of effective harmonic
//        */
//       double    regularization
//       ( const double sigma     ,
//         const double delta     ) const ;
//       // ======================================================================

//       // ======================================================================
//       /** get the effective cut-off (==number of terms/harmonics)
//        *  of Tikhonov's regularization
//        *  \f$ n \equiv  \sqrt{2 \ln \delta} \frac{\pi\sigma}{L} \f$
//        *  @param sigma  gaussian resoltuion
//        *  @param delta  regularization parameter
//        *  @return number of effective harmonic
//        */
//       double    regularization
//       ( const double sigma     ,
//         const double delta     ) const ;
//       // ======================================================================



// ============================================================================
//                                                                      The END 
// ============================================================================


