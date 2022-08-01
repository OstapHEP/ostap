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
// ============================================================================
//  Local 
// ============================================================================
#include "Exception.h"
#include "local_math.h"
// ============================================================================
/** @file
 *  Implementation file for functions from the file Ostap/Fourier.h
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2010-04-19
 */
// ============================================================================
// constructor from the degree 
// ============================================================================
Ostap::Math::FourierSum::FourierSum 
( const unsigned short degree , 
  const double         xmin   , 
  const double         xmax   , 
  const bool           fejer  )
  : Ostap::Math::Parameters ( 2 * degree + 1 )
  , m_xmin  ( std::min ( xmin , xmax ) )
  , m_xmax  ( std::max ( xmin , xmax ) )
  , m_scale ( 1 ) 
  , m_delta ( 0 ) 
  , m_fejer ( fejer ) 
{
  //
  if ( s_equal ( -M_PI     , m_xmin ) ) { m_xmin = -M_PI     ; }
  if ( s_equal ( -1        , m_xmin ) ) { m_xmin = -1        ; }
  if ( s_equal (  0        , m_xmin ) ) { m_xmin =  0        ; }
  //
  if ( s_equal (  1        , m_xmax ) ) { m_xmax =  1        ; }
  if ( s_equal (      M_PI , m_xmax ) ) { m_xmax =  M_PI     ; }
  if ( s_equal (  2 * M_PI , m_xmax ) ) { m_xmax =  2 * M_PI ; }
  //
  m_scale = 2 * M_PI / ( m_xmax - m_xmin ) ;
  m_delta = 0.5      * ( m_xmax + m_xmin ) ;
}
// ============================================================================
// constructor from cosine series 
// ============================================================================
Ostap::Math::FourierSum::FourierSum 
( const Ostap::Math::CosineSum& sum )
  : Ostap::Math::Parameters ( 2 * sum.degree() + 1 )
  , m_xmin  ( 2 * sum.xmin() - sum.xmax() )
  , m_xmax  ( sum.xmax() )
  , m_scale ( 1 ) 
  , m_delta ( 0 ) 
  , m_fejer ( sum.fejer() )
{
  //
  if ( s_equal ( -M_PI     , m_xmin ) ) { m_xmin = -M_PI     ; }
  if ( s_equal ( -1        , m_xmin ) ) { m_xmin = -1        ; }
  if ( s_equal (  0        , m_xmin ) ) { m_xmin =  0        ; }
  //
  if ( s_equal (  1        , m_xmax ) ) { m_xmax =  1        ; }
  if ( s_equal (      M_PI , m_xmax ) ) { m_xmax =  M_PI     ; }
  if ( s_equal (  2 * M_PI , m_xmax ) ) { m_xmax =  2 * M_PI ; }
  //
  m_scale = 2 * M_PI / ( m_xmax - m_xmin ) ;
  m_delta = 0.5      * ( m_xmax + m_xmin ) ;
  //
  for ( unsigned short i = 0 ; i <= degree() ; ++i ) 
  { setA ( i , sum.par(i) ) ; }
  //
}
// ============================================================================
// constructor from Fourier series
// ============================================================================
Ostap::Math::FourierSum::FourierSum 
( const Ostap::Math::FourierSum& sum   , 
  const bool                     fejer )
  : Ostap::Math::Parameters ( sum.m_pars  ) 
  , m_xmin  ( sum.m_xmin  )
  , m_xmax  ( sum.m_xmax  )
  , m_scale ( sum.m_scale )
  , m_delta ( sum.m_delta ) 
  , m_fejer ( fejer   )
{}
// ============================================================================
// copy constructor from Fourier series
// ============================================================================
Ostap::Math::FourierSum::FourierSum 
( const Ostap::Math::FourierSum& sum )
  : Ostap::Math::Parameters ( sum.m_pars  ) 
  , m_xmin  ( sum.m_xmin  )
  , m_xmax  ( sum.m_xmax  )
  , m_scale ( sum.m_scale )
  , m_delta ( sum.m_delta ) 
  , m_fejer ( sum.m_fejer )
{}
// ============================================================================
// move constructor from Fourier series
// ============================================================================
Ostap::Math::FourierSum::FourierSum
(       Ostap::Math::FourierSum&& sum ) 
  : Ostap::Math::Parameters ( std::move ( sum.m_pars  ) ) 
  , m_xmin  ( std::move ( sum.m_xmin  ) )
  , m_xmax  ( std::move ( sum.m_xmax  ) )
  , m_scale ( std::move ( sum.m_scale ) )
  , m_delta ( std::move ( sum.m_delta ) )
  , m_fejer ( std::move ( sum.m_fejer ) )
{}
// ============================================================================
// protected constructor from the parameters 
// ============================================================================
Ostap::Math::FourierSum::FourierSum
( const std::vector<double>& pars  , 
  const double               xmin  , 
  const double               xmax  , 
  const double               fejer )
  : Ostap::Math::Parameters ( pars ) 
  , m_xmin  ( std::min ( xmin , xmax ) )
  , m_xmax  ( std::max ( xmin , xmax ) )
  , m_scale ( 1 ) 
  , m_delta ( 0 ) 
  , m_fejer ( fejer ) 
{
  Ostap::Assert ( 1 ==  pars.size() % 2                        ,  
                  "odd number of parameters must be supplied!" , 
                  "Ostap::Math::FourierSum"                  ) ;
  //
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
  std::swap ( m_fejer ,  right.m_fejer ) ; 
}
// ============================================================================
/* get the magnitude of nth-harmonic
 * \f$m_k = \sqrt( a^2_k + b^2_k) \f$
 */
// ============================================================================
double Ostap::Math::FourierSum::mag    ( const unsigned short k ) const 
{
  if      ( k > degree() ) { return 0 ; }
  else if ( 0 == k       ) { return std::abs ( m_pars[0] ) ; }
  //
  return std::hypot ( m_pars[ 2 * k - 1 ] , m_pars [ 2 * k ] ) ;
}
// ============================================================================
// get the phase of nth-harmonic
// ============================================================================
double Ostap::Math::FourierSum::phase ( const unsigned short k ) const 
{
  if      ( k > degree() ) { return 0 ; }
  else if ( 0 == k       ) { return 0 <= m_pars[0] ? 0. : -M_PI ; }
  //
  return std::atan2 ( m_pars[ 2 * k - 1 ] , m_pars [ 2 * k ] ) ;
} 
// ============================================================================
// calculate Fourier sum 
// ============================================================================
double Ostap::Math::FourierSum::fourier_sum ( const double x ) const 
{
  /// transform to "t"-representation 
  const long double tv = t(x) ;
  return Ostap::Math::Clenshaw::fourier_sum ( begin() , end () , tv ) ;
}
// ============================================================================
// calculate Fejer sum 
// ============================================================================
double Ostap::Math::FourierSum::fejer_sum ( const double x ) const 
{
  /// transform to "t"-representation 
  const long double tv = t(x) ;
  return Ostap::Math::Clenshaw::fejer_sum ( begin() , end ()  , tv ) ;
}
// ============================================================================
// get Fejer sum 
// ============================================================================
Ostap::Math::FourierSum 
Ostap::Math::FourierSum::fejer_sum   () const                  // get Fejer sum 
{
  // create fejer sum obejct 
  FourierSum fejer ( m_pars , m_xmin , m_xmax , false ) ;
  // fill it! 
  const unsigned long N  = m_pars.size() ;
  const double   long fd = 1.0L / ( N + 1 ) ;
  /// start scaling of harmonics 
  for ( unsigned short k = 1 ; 2 * k < N ; ++k  ) 
  {
    const long double f  = ( N + 1 - 2 * k ) * fd ;
    fejer.m_pars[ 2 * k - 1 ] *= f ;
    fejer.m_pars[ 2 * k     ] *= f ;
  }
  //
  return fejer ;
}
// ============================================================================
//  get the derivative at point x 
// ============================================================================
double Ostap::Math::FourierSum::derivative ( const double x ) const 
{
  //
  std::vector<double> deriv ( m_pars.begin() , m_pars.end() ) ; 
  deriv[0] = 0.0 ;
  //
  const unsigned long  N = m_pars.size() ;
  for ( unsigned short k = 1 ; 2 * k < N ; ++k  ) 
  {
    deriv[ 2 * k     ] =  k * m_pars[ 2 * k - 1 ] * m_scale ;
    deriv[ 2 * k - 1 ] = -k * m_pars[ 2 * k     ] * m_scale ;
    //
  }
  /// transform to "t"-representation 
  const long double tv = t(x) ;
  /// make evaluation of fourier serie
  return 
    m_fejer ? 
    Ostap::Math::Clenshaw::fejer_sum   ( deriv.begin() , deriv.end () , tv ) :
    Ostap::Math::Clenshaw::fourier_sum ( deriv.begin() , deriv.end () , tv ) ;
}
// ============================================================================
//  get the derivative 
// ============================================================================
Ostap::Math::FourierSum Ostap::Math::FourierSum::derivative () const 
{ return derivative_n ( 1 ) ; }
// ============================================================================
//  get the nth derivative 
// ============================================================================
Ostap::Math::FourierSum 
Ostap::Math::FourierSum::derivative_n ( const unsigned short n ) const 
{
  //
  if      ( 0 == n ) { return *this ; }
  // create derivate obejct 
  FourierSum deriv ( m_pars , m_xmin , m_xmax , m_fejer ) ;
  /// fill it! 
  deriv.m_pars [0] = 0.0 ;
  const unsigned long  N = m_pars.size() ;
  //
  const short ssin =   
    1 == n % 4 ?  1 : 
    2 == n % 4 ? -1 : 
    3 == n % 4 ? -1 : 1 ;
  //
  const short scos =   
    1 == n % 4 ? -1 : 
    2 == n % 4 ? -1 : 1 ;
  //
  const bool           odd =  ( 1 == n % 2 ) ;
  for ( unsigned short k = 1 ; 2 * k < N ; ++k  ) 
  {
    const double factor = Ostap::Math::POW ( 1.0L * k * m_scale , n ) ;
    if ( odd ) 
    {
      deriv.m_pars [ 2 * k     ] = m_pars [ 2 * k - 1 ] * factor * ssin ;
      deriv.m_pars [ 2 * k - 1 ] = m_pars [ 2 * k     ] * factor * scos ;
    }
    else 
    {
      deriv.m_pars [ 2 * k     ] = m_pars [ 2 * k     ] * factor * scos ;
      deriv.m_pars [ 2 * k - 1 ] = m_pars [ 2 * k - 1 ] * factor * ssin ;
    }
    //
  }
  //
  return deriv ;
}
// ============================================================================
// get integral between low and high 
// ============================================================================
double Ostap::Math::FourierSum::integral 
( const double low , const double high ) const 
{
  if ( s_equal ( low , high ) ) { return 0 ; }
  //
  // optionally shrink interval for period. but not now..
  //
  std::vector<double> integ ( m_pars.begin() , m_pars.end() ) ; 
  //
  integ[0] = 0 ;
  const unsigned long  N = m_pars.size() ;
  for ( unsigned short k = 1 ; 2 * k < N ; ++k  ) 
  {
    integ[ 2 * k     ] = -m_pars[ 2 * k - 1 ] / ( k * m_scale ) ;
    integ[ 2 * k - 1 ] =  m_pars[ 2 * k     ] / ( k * m_scale ) ;
  }
  /// transform to "t"-representation 
  const long double tl = t ( low  ) ;
  const long double th = t ( high ) ;
  /// evaluate Fourier series
  return
    m_fejer ? 
    Ostap::Math::Clenshaw::fejer_sum   ( integ.begin() , integ.end() , th ) - 
    Ostap::Math::Clenshaw::fejer_sum   ( integ.begin() , integ.end() , tl ) +
    0.5 * m_pars[0] * ( th - tl ) :
    Ostap::Math::Clenshaw::fourier_sum ( integ.begin() , integ.end() , th ) - 
    Ostap::Math::Clenshaw::fourier_sum ( integ.begin() , integ.end() , tl ) +
    0.5 * m_pars[0] * ( th - tl ) ;
}
// ============================================================================
// get integral as object 
// ============================================================================
Ostap::Math::FourierSum
Ostap::Math::FourierSum::integral ( const double c0 ) const 
{
  //
  FourierSum integ ( m_pars , m_xmin , m_xmax , m_fejer ) ;
  //
  integ.m_pars[0] = c0 ;
  const unsigned long  N   =  m_pars.size() ;
  const double         a0  =  m_pars[0]     ;
  const bool           add = !s_zero ( a0 ) ;
  for ( unsigned short k   =  1 ; 2 * k < N ; ++k  ) 
  {
    //
    const double a_cos = -m_pars[ 2*k-1 ]  / ( k * m_scale  ) ;
    const double a_sin =  m_pars[ 2*k   ]  / ( k * m_scale  ) ;
    
    integ.m_pars[ 2 * k     ] = a_cos ;
    integ.m_pars[ 2 * k - 1 ] = a_sin ;
    //
    // add a serie for f(x) = 2*a0*x 
    if ( add ) { integ.m_pars [ 2 * k - 1 ] -= ( 0 == k%2 ? 1 : -1 ) * a0 / ( k * m_scale ) ; }
  }
  // add integration constant 
  integ.setPar( 0 , 2 * c0 );
  return integ ;
}
// ============================================================================

// ============================================================================
Ostap::Math::FourierSum
Ostap::Math::FourierSum::convolve 
( const double sigma ) const 
{
  //
  // no convolution 
  if ( s_zero ( sigma ) ) { return *this ; }
  //
  const long double ss      =  sigma / m_scale ;
  const long double sigma2  =  ss*ss           ;
  // create covolution obejct 
  FourierSum conv( m_pars , m_xmin , m_xmax , m_fejer ) ;
  /// fill it! 
  conv.m_pars [0] = m_pars[0]  ;
  const unsigned long  N = m_pars.size() ;
  for ( unsigned short k = 1 ; 2 * k < N ; ++k  ) 
  {
    const long double s  = std::exp ( - 0.5L * k * k * sigma2 ) ;
    //
    const long double v1 = s * m_pars[ 2 * k    ] ;
    if ( !s_zero ( v1 ) ) { conv.m_pars [ 2 * k    ] = v1 ; }
    //
    const long double v2 = s * m_pars[ 2 * k -1 ] ;
    if ( !s_zero ( v2 ) ) { conv.m_pars [ 2 * k -1 ] = v2 ; }
    //
  }
  //
  return conv ;
}
// ============================================================================
//  convolute Fourier sum with gaussian function 
// ============================================================================
Ostap::Math::FourierSum 
Ostap::Math::FourierSum::deconvolve 
( const double sigma , 
  const double delta ) const 
{
  // no convolution 
  if ( s_zero ( sigma ) ) { return *this ; }
  //
  const long double ss      =  sigma / m_scale ;
  const long double sigma2  =  ss*ss           ;
  // create covolution object 
  FourierSum conv( m_pars , m_xmin , m_xmax , m_fejer ) ;
  /// fill it! 
  conv.m_pars [0] = m_pars[0]  ;
  const unsigned long  N = m_pars.size() ;
  //
  const bool use_delta = !s_zero ( delta ) && 0 < delta ;
  for ( unsigned short k = 1 ; 2 * k < N ; ++k  ) 
  {
    //  
    const double v_1  = m_pars[2*k] ;
    const double v_2  = m_pars[2*k-1] ;
    if ( s_zero ( v_1 )  && s_zero ( v_2 ) ) { continue ; }
    //
    long double f = my_exp ( 0.5L * k * k * sigma2 ) ;
    //
    if ( use_delta ) 
    { const long double fd = f * delta ; f /= ( 1 + fd * fd ) ; }
    //
    const long double   v1 = f * v_1 ;
    if ( !s_zero ( v1 ) ) { conv.m_pars [ 2 * k    ] = v1 ; }
    else { conv.m_pars[2*k  ] = 0 ; }    
    //
    const long double   v2 = f * v_2 ;
    if ( !s_zero ( v2 ) ) { conv.m_pars [ 2 * k -1 ] = v2 ; }
    else { conv.m_pars[2*k-1] = 0 ; }    
    //
  }
  //
  return conv ;
}
// ============================================================================
/* get the effective cut-off (==number of effective harmonics) 
 * of Tikhonov's regularization 
 * \f$ n \equiv  \sqrt{2 \ln \delta} \frac{2\pi\sigma}{L} \f$
 * @param sigma  gaussian resolution 
 * @param delta  regularization parameter 
 * @return number of effective harmonic 
 */
// ============================================================================
double Ostap::Math::FourierSum::regularization 
( const double sigma , 
  const double delta ) const 
{
  if      ( 0 > delta || s_zero ( delta ) || s_zero ( sigma ) ) 
  { return s_UL_max ; } // return
  else if ( 1<= delta ) { return 1 ; }
  //
  return std::sqrt ( -2 * std::log ( delta ) ) * m_scale / std::abs ( sigma ) ;
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
  m_pars[0] += 2.0*a ;
  return *this ;
}
// ============================================================================
// simple  manipulations with polynoms: subtract constant 
// ============================================================================
Ostap::Math::FourierSum&
Ostap::Math::FourierSum::operator-=( const double a ) 
{
  m_pars[0] -= 2.0*a ;
  return *this ;
}
// =============================================================================
/*  sum of two Fourier series (with the same interval!) 
 *  @param other the first fourier sum
 *  @return the sum of two Fourier series 
 */
// =============================================================================
Ostap::Math::FourierSum 
Ostap::Math::FourierSum::sum( const Ostap::Math::FourierSum& other ) const 
{
  //
  if ( this == &other ) 
  {
    FourierSum result ( *this  ) ;
    result  *= 2 ;
    return result ;
  }
  //
  if      ( other.zero() ) { return *this ; } // 
  else if (       zero() ) { return other ; } // random choice 
  //
  if ( !s_equal ( xmin () , other.xmin() ) || 
       !s_equal ( xmax () , other.xmax() ) ) 
  {
    Ostap::throwException ( "Can't sum Fourier series with different domains" , 
                            "Ostap::Math::FourierSum" ) ;
  }
  if ( fejer() != other.fejer () ) 
  {
    Ostap::throwException ( "Can't sum Fourier series with different 'fejer' flag" , 
                            "Ostap::Math::FourierSum"  ) ;
  }
  //
  const unsigned short idegree = std::max ( degree () , other.degree () ) ;
  //
  FourierSum result ( idegree , xmin() , xmax() , fejer() ) ;
  const unsigned npars  = result.npars() ;
  for ( unsigned short i = 0 ; i < npars ; ++i ) 
  { result.m_pars [i] =  par(i)  + other.par(i) ; }
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
( const unsigned short degree , 
  const double         xmin   , 
  const double         xmax   , 
  const bool           fejer  )
  : Ostap::Math::Parameters ( degree + 1 ) 
  , m_xmin  ( std::min ( xmin , xmax ) )
  , m_xmax  ( std::max ( xmin , xmax ) )
  , m_scale ( 1 ) 
  , m_fejer ( fejer ) 
{
  //
  if ( s_equal ( -M_PI     , m_xmin ) ) { m_xmin = -M_PI     ; }
  if ( s_equal ( -1        , m_xmin ) ) { m_xmin = -1        ; }
  if ( s_equal (  0        , m_xmin ) ) { m_xmin =  0        ; }
  //
  if ( s_equal (  1        , m_xmax ) ) { m_xmax =  1        ; }
  if ( s_equal (      M_PI , m_xmax ) ) { m_xmax =  M_PI     ; }
  if ( s_equal (  2 * M_PI , m_xmax ) ) { m_xmax =  2 * M_PI ; }
  //
  m_scale = M_PI / ( m_xmax - m_xmin ) ;
}
// ============================================================================
// constructor from Fourier sum 
// ============================================================================
Ostap::Math::CosineSum::CosineSum 
( const Ostap::Math::FourierSum& sum ) 
  : Ostap::Math::Parameters ( sum.degree() + 1  )
  , m_xmin  ( 0.5 * ( sum.xmax() + sum.xmin() ) )
  , m_xmax  (         sum.xmax()                )
  , m_scale ( 1 ) 
  , m_fejer ( sum.fejer() ) 
{
  //
  if ( s_equal ( -M_PI     , m_xmin ) ) { m_xmin = -M_PI     ; }
  if ( s_equal ( -1        , m_xmin ) ) { m_xmin = -1        ; }
  if ( s_equal (  0        , m_xmin ) ) { m_xmin =  0        ; }
  //
  if ( s_equal (  1        , m_xmax ) ) { m_xmax =  1        ; }
  if ( s_equal (      M_PI , m_xmax ) ) { m_xmax =  M_PI     ; }
  if ( s_equal (  2 * M_PI , m_xmax ) ) { m_xmax =  2 * M_PI ; }
  //
  m_scale = M_PI / ( m_xmax - m_xmin ) ;
  //
  for ( unsigned short i = 0 ; i< m_pars.size() ; ++i ) 
  { setPar ( i , sum.a(i) ) ; }
}
// ============================================================================
// constructor from Fourier series
// ============================================================================
Ostap::Math::CosineSum::CosineSum 
( const Ostap::Math::CosineSum& sum   ,
  const bool                    fejer )
  : Ostap::Math::Parameters ( sum )
  , m_xmin  ( sum.m_xmin  )
  , m_xmax  ( sum.m_xmax  )
  , m_scale ( sum.m_scale )
  , m_fejer ( fejer   )
{}
// ============================================================================
// copy constructor from Fourier series
// ============================================================================
Ostap::Math::CosineSum::CosineSum 
( const Ostap::Math::CosineSum& sum )
  : Ostap::Math::Parameters ( sum  )
  , m_xmin  ( sum.m_xmin  )
  , m_xmax  ( sum.m_xmax  )
  , m_scale ( sum.m_scale )
  , m_fejer ( sum.m_fejer )
{}
// ============================================================================
// move constructor from Fourier series
// ============================================================================
Ostap::Math::CosineSum::CosineSum
(       Ostap::Math::CosineSum&& sum ) 
  : Ostap::Math::Parameters ( std::move ( sum ) ) 
  , m_xmin  ( std::move ( sum.m_xmin  ) )
  , m_xmax  ( std::move ( sum.m_xmax  ) )
  , m_scale ( std::move ( sum.m_scale ) )
  , m_fejer ( std::move ( sum.m_fejer ) )
{}
// ============================================================================
// constructor from non-empty list of parameters 
// ============================================================================
Ostap::Math::CosineSum::CosineSum
( const std::vector<double>& pars  , 
  const double               xmin  , 
  const double               xmax  , 
  const double               fejer )
  : Ostap::Math::Parameters ( pars )
  , m_xmin  ( std::min ( xmin , xmax ) )
  , m_xmax  ( std::max ( xmin , xmax ) )
  , m_scale ( 1 ) 
  , m_fejer ( fejer ) 
{
  Ostap::Assert ( 1 <= pars.size()                       , 
                  "List ofparmaeters must be non-empty!" , 
                  "Ostap::Math::CosineSum"               ) ;
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
  std::swap ( m_fejer ,  right.m_fejer ) ; 
}
// ============================================================================
// calculate Fourier sum 
// ============================================================================
double Ostap::Math::CosineSum::fourier_sum ( const double x ) const 
{
  /// transform to "t"-representation 
  const long double tv = t(x) ;
  return Ostap::Math::Clenshaw::cosine_sum ( begin() , end ()  , tv ) ;
}
// ============================================================================
// calculate Fejer sum 
// ============================================================================
double Ostap::Math::CosineSum::fejer_sum ( const double x ) const 
{
  /// transform to "t"-representation 
  const long double tv = t(x) ;
  return Ostap::Math::Clenshaw::fejer_cosine_sum ( begin() , end ()  , tv ) ;
}
// ============================================================================
// get Fejer sum 
// ============================================================================
Ostap::Math::CosineSum 
Ostap::Math::CosineSum::fejer_sum   () const                  // get Fejer sum 
{
  // create fejer sum e obejct 
  CosineSum fejer ( m_pars , m_xmin , m_xmax , false ) ;
  // fill it! 
  const unsigned long N  = m_pars.size() ;
  const double   long fd = 1.0L / N  ;
  /// start scaling of harmonics 
  for ( unsigned short k = 0 ;  k < N ; ++k  ) 
  {
    const long double f  = ( N - k ) * fd ;
    fejer.m_pars[  k  ] *= f ;
  }
  //
  return fejer ;
}
// ============================================================================
Ostap::Math::CosineSum
Ostap::Math::CosineSum::convolve 
( const double sigma ) const 
{
  //
  // no convolution 
  if ( s_zero ( sigma ) ) { return *this ; }
  //
  const long double ss      =  sigma / m_scale ;
  const long double sigma2  =  ss*ss           ;
  // create covolution obejct 
  CosineSum conv( m_pars , m_xmin , m_xmax , m_fejer ) ;
  /// fill it! 
  conv.m_pars [0] = m_pars[0]  ;
  const unsigned long  N = m_pars.size() ;
  for ( unsigned short k = 1 ; k < N ; ++k  ) 
  {
    const long double s  = std::exp ( - 0.5L * k * k * sigma2 ) ;
    //
    const long double v1 = s * m_pars   [ k ]      ;
    if ( !s_zero ( v1 ) ) { conv.m_pars [ k ] = v1 ; }
  }
  //
  return conv ;
}
// ============================================================================
//  convolute Fourier sum with gaussian function 
// ============================================================================
Ostap::Math::CosineSum 
Ostap::Math::CosineSum::deconvolve 
( const double sigma , 
  const double delta ) const 
{
  // no convolution 
  if ( s_zero ( sigma ) ) { return *this ; }
  //
  const long double ss      =  sigma / m_scale ;
  const long double sigma2  =  ss*ss           ;
  // create covolution obejct 
  CosineSum conv( m_pars , m_xmin , m_xmax , m_fejer ) ;
  /// fill it! 
  conv.m_pars [0] = m_pars[0]  ;
  const unsigned long  N = m_pars.size() ;
  const bool use_delta = !s_zero ( delta ) && 0 < delta ;
  for ( unsigned short k = 1 ; k < N ; ++k  ) 
  {
    //  
    const double v = m_pars[k] ;
    if ( s_zero ( v ) ) { continue ; }
    //
    long double f = my_exp ( 0.5L * k * k * sigma2 ) ;
    //
    if ( use_delta ) 
    { const long double fd = f * delta ; f /= ( 1 + fd * fd ) ; }
    //
    const long double   v1 = f * v ;
    if ( !s_zero ( v1 ) ) { conv.m_pars [ k ] = v1 ; }
    else { conv.m_pars[k] = 0 ; }    
    //
  }
  //
  return conv ;
}
// ============================================================================
/* Get the effective cut-off (==number of terms/harmonics) 
 * of Tikhonov's regularization 
 * \f$ n \equiv  \sqrt{2 \ln \delta} \frac{\pi\sigma}{L} \f$
 * @param sigma  gaussian resolution 
 * @param delta  regularization parameter 
 * @return number of effective harmonic 
 */
// ============================================================================
double Ostap::Math::CosineSum::regularization 
( const double sigma , 
  const double delta ) const 
{
  if      ( 0 > delta || s_zero ( delta ) || s_zero ( sigma ) ) { return s_UL_max ; } 
  else if ( 1<= delta ) { return 1 ; }
  //
  return std::sqrt ( -2 * std::log ( delta ) ) * m_scale / std::abs ( sigma ) ;
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
  m_pars[0] += 2.0*a ;
  return *this ;
}
// ============================================================================
// simple  manipulations with polynoms: subtract constant 
// ============================================================================
Ostap::Math::CosineSum&
Ostap::Math::CosineSum::operator-=( const double a )
{  
  m_pars[0] -= 2.0*a ;
  return *this ;
}
// ============================================================================
// get the derivative at point x 
// ============================================================================
double Ostap::Math::CosineSum::derivative ( const double x ) const 
{
  //
  std::vector<double> deriv ( m_pars.size() - 1 , 0.0 ) ; 
  deriv[0] = 0.0 ;
  //
  for ( unsigned short k = 0 ; k < deriv.size()  ; ++k  )
  { 
    deriv[k]  =  -m_pars[k+1] * ( k + 1 ) * m_scale ; 
  }
  //
    
  /// transform to "t"-representation 
  const long double tv = t(x) ;
  /// make evaluation of fourier serie
  return 
    m_fejer ? 
    Ostap::Math::Clenshaw::fejer_sine_sum   ( deriv.begin() , deriv.end () , tv ) :
    Ostap::Math::Clenshaw::sine_sum         ( deriv.begin() , deriv.end () , tv ) ;
}
// ============================================================================
//  get the derivative at point x 
// ============================================================================
Ostap::Math::FourierSum Ostap::Math::CosineSum::derivative () const 
{ return derivative_n ( 1 ) ; }
// ============================================================================
//  get the nth derivative 
// ============================================================================
Ostap::Math::FourierSum 
Ostap::Math::CosineSum::derivative_n ( const unsigned short n ) const 
{
  //
  if      ( 0 == n ) { return *this ; }
  // create derivate obejct 
  FourierSum deriv ( *this ) ;
  //
  /// fill it! 
  const unsigned long  N = m_pars.size() ;
  //
  const short scos =   
    1 == n % 4 ? -1 : 
    2 == n % 4 ? -1 : 1 ;
  //
  const bool         odd =  ( 1 == n % 2 ) ;
  for ( unsigned short k = 1 ;  k < N ; ++k  ) 
  {
    const long double factor = Ostap::Math::POW ( 1.0L * k * m_scale , n ) ;
    if ( odd ) 
    {
      deriv.setPar (  2 * k - 1  , m_pars [ k ] * factor * scos ) ;
      deriv.setPar (  2 * k      , 0                            ) ;
    }
    else       
    { 
      deriv.setPar ( 2 * k       , m_pars [ k ] * factor * scos ) ;
      deriv.setPar ( 2 * k  - 1  , 0                            ) ;
    }    
    //
  }
  //
  return deriv ;
}
// ============================================================================
// get integral between low and high 
// ============================================================================
double Ostap::Math::CosineSum::integral 
( const double low , const double high ) const 
{
  if ( s_equal ( low , high ) ) { return 0 ; }
  //
  // optionally shrink interval for period. but not now..
  //
  std::vector<double> integ ( m_pars.size() - 1 , 0.0 ) ; 
  //
  for ( unsigned short k = 0 ; k < integ.size() ; ++k  ) 
  { integ[k] = m_pars[k+1] / ( k + 1 )  / m_scale ; }
  /// transform to "t"-representation 
  const long double tl = t ( low  ) ;
  const long double th = t ( high ) ;
  /// evaluate Fourier series
  return
    m_fejer ? 
    Ostap::Math::Clenshaw::fejer_sine_sum ( integ.begin() , integ.end() , th ) - 
    Ostap::Math::Clenshaw::fejer_sine_sum ( integ.begin() , integ.end() , tl ) +
    0.5 * m_pars[0] * ( th - tl ) / m_scale :
    Ostap::Math::Clenshaw::sine_sum       ( integ.begin() , integ.end() , th ) - 
    Ostap::Math::Clenshaw::sine_sum       ( integ.begin() , integ.end() , tl ) +
    0.5 * m_pars[0] * ( th - tl ) / m_scale ;
}
// ============================================================================
// get integral as object 
// ============================================================================
Ostap::Math::FourierSum
Ostap::Math::CosineSum::integral ( const double c0 ) const 
{
  //
  FourierSum integ ( *this ) ;
  //
  integ.setPar ( 0 , c0 ) ;
  const unsigned long  N   =  m_pars.size() ;
  const double         a0  =  0.5 *  m_pars[0]     ;
  const bool           add = !s_zero ( a0 ) ;
  for ( unsigned short k   =  1 ; k < N ; ++k  ) 
  {
    // integration of cosine 
    integ.setPar ( 2 * k - 1  , m_pars[k] / ( k * m_scale )  ) ; 
    //
    // integration of cconstant term 
    if ( add && 0 != k%2 ) 
    { integ.setPar  ( 2 * k , -4 * a0  / ( k * k * m_scale  * M_PI ) ) ; }    
  }  
  //
  // add integration constant 
  integ.setPar( 0 , 2*c0  + a0 * M_PI / m_scale );
  return integ ;
}
// ============================================================================
/*  sum of two Fourier series (with the same interval!) 
 *  @param other the first fourier sum
 *  @return the sum of two Fourier series 
 */
// ============================================================================
Ostap::Math::CosineSum 
Ostap::Math::CosineSum::sum ( const Ostap::Math::CosineSum& other ) const 
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
  if ( !s_equal ( xmin () , other.xmin() ) ||
       !s_equal ( xmax () , other.xmax() ) ) 
  {
    Ostap::throwException ( "Can't sum Fourier cosine series with different domains" , 
                            "Ostap::Math::CosineSum" ) ;
  }
  if ( fejer() != other.fejer () ) 
  {
    Ostap::throwException ( "Can't sum Fourier cosine series with different 'fejer' flag" , 
                            "Ostap::Math::CosineSum" ) ;
  }
  //
  const unsigned short idegree = std::max ( degree () , other.degree () ) ;
  //
  CosineSum result ( idegree , xmin() , xmax() , fejer() ) ;
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
/*  make a sum of two fourier series (with the same interval!) 
 *  @param s1 the first fourier sum
 *  @param s2 the first fourier sum 
 *  @return s1+s2 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2016-06-26
 */
// ============================================================================
Ostap::Math::FourierSum
Ostap::Math::sum 
( const Ostap::Math::FourierSum& s1 ,
  const Ostap::Math::FourierSum& s2 ) { return s1.sum ( s2 ) ; }
// ============================================================================
/*  make a sum of two fourier cosine series (with the same interval!) 
 *  @param s1 the first fourier cosine sum
 *  @param s2 the first fourier cosine sum 
 *  @return s1+s2 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2016-06-26
 */
// ============================================================================
Ostap::Math::CosineSum 
Ostap::Math::sum 
( const Ostap::Math::CosineSum& s1 , 
  const Ostap::Math::CosineSum& s2 ) { return s1.sum ( s2 ) ; }
// ============================================================================
//                                                                      The END 
// ============================================================================


