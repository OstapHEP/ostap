// ============================================================================
// Include files 
// ============================================================================
// STD&STL
// ============================================================================
#include <climits>
#include <cassert>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Math.h"
#include "Ostap/Hash.h"
#include "Ostap/Power.h"
#include "Ostap/Clenshaw.h"
#include "Ostap/Choose.h"
#include "Ostap/Bernstein.h"
#include "Ostap/MoreMath.h"
#include "Ostap/Polynomials.h"
#include "Ostap/Positive.h"
#include "Ostap/Bernstein1D.h"
// ============================================================================
// Local 
// ============================================================================
#include "Exception.h"
// ============================================================================
/** @file 
 *  Implementation file for classes from file Ostap/Polynomials.h
 *  @see LHCbMath/Polynomials.h
 *  @date 2015-02-23 
 *  @author Vanya Belyaev Ivan Belyaev@itep.ru
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  static_assert ( std::numeric_limits<double>::is_specialized           , 
                  "mumeric_limits are not specialized for doubles"      ) ;
  // ==========================================================================
  /// equality criteria for doubles
  const Ostap::Math::Equal_To<double> s_equal {} ;       // equality criteria for doubles
  /// zero for doubles  
  const Ostap::Math::Zero<double>     s_zero  {} ;       // zero for doubles
  /// zero fo vectors 
  const Ostap::Math::Zero< std::vector<double> > s_vzero {} ; // zero for vectors
  /// small element 
  const long double            s_epsilon =  2 * std::numeric_limits<double>::epsilon() ;
  /// small element 
  Ostap::Math::Small<double>    s_small ( s_epsilon ) ;
  // ==========================================================================
  /// get a factorial 
  inline long double 
  _factorial_d_ ( const unsigned short N ) 
  {
    return 
      0 == N ?  1 : 
      1 == N ?  1 : 
      2 == N ?  2 : 
      3 == N ?  6 : 
      4 == N ? 24 : N * _factorial_d_ ( N - 1 ) ;
  }
  // ==========================================================================
}
// ============================================================================
// useful local utilities 
// ============================================================================
namespace 
{
  // =========================================================================
  // Chebyshev polynomials (1st and 2nd kind) 
  // =========================================================================
  /** evaluate chebyshev polinomial (the second kind)
   *  @param N  the polynomial degree 
   *  @param x  the argument 
   *  @return value of chebyshev polynomial (2nd kind) at point x 
   */
  inline double _chebyshevU_ 
  ( const unsigned int N , 
    const double       x )
  {
    return 
      0 == N ? 1.0   :
      1 == N ? 2 * x : 
      2 == N ? 4 * x * x - 1 :
      3 == N ? 4 * x * ( 2 * x * x - 1 ) :
      s_equal ( x ,  1 ) ?                1.0 + N            :
      s_equal ( x , -1 ) ? ( 0 == N % 2 ? 1.0 + N : -1.0-N ) :
      _chebyshevU_ ( N - 2 , x ) * ( _chebyshevU_ ( 2 , x ) - 1 )  -
      _chebyshevU_ ( N - 4 , x ) ;
  } 
  // ==========================================================================
  /** evaluate the integral for chebychev polynomial (1st kind) 
   *  @param N    the polynomial degree 
   *  @param low  low edge of integration 
   *  @param high high edge of integration 
   */  
  inline double _chebyshev_int_ 
  ( const unsigned int N    ,
    const double       low  , 
    const double       high ) 
  {
    // trivial cases 
    if      ( s_equal ( low , high ) ) { return 0          ; }
    else if ( 0 == N                 ) { return high - low ; }
    else if ( 1 == N                 )
    { return 0.5 * ( high * high - low * low ) ; }
    else if ( high < low ) 
    { return -_chebyshev_int_ ( N ,  high , low )  ; }
    //
    const double ihigh = 
      Ostap::Math::chebyshev_value ( N + 1 , high ) / ( N + 1 ) -
      Ostap::Math::chebyshev_value ( N - 1 , high ) / ( N - 1 ) ;
    //
    const double ilow = 
      Ostap::Math::chebyshev_value ( N + 1 , low  ) / ( N + 1 ) -
      Ostap::Math::chebyshev_value ( N - 1 , low  ) / ( N - 1 ) ;
    //
    return 0.5 * ( ihigh - ilow ) ;
  }
  // ==========================================================================
  /** evaluate the derivative for chebychev polynomial (1st kind)
   *  @param N the polynomial degree
   *  @param x the point 
   */
  inline double _chebyshev_der_ 
  ( const unsigned int N ,
    const double       x )  
  {
    return  
      0 == N ? 0.0 : 
      1 == N ? 1.0 :
      2 == N ? 4 * x :  N * _chebyshevU_ ( N - 1 , x ) ;
  }
  // ==========================================================================
  /** evaluate the integral for chebychev polynomial (2nd kind) 
   *  @param N    the polynomial degree 
   *  @param low  low edge of integration 
   *  @param high high edge of integration 
   */  
  inline double _chebyshevU_int_
  ( const unsigned int N    ,
    const double       low  , 
    const double       high ) 
  {
    // trival cases 
    if      ( s_equal ( low , high ) ) { return 0          ; }
    else if ( 0 == N                 ) { return high - low ; }
    else if ( 1 == N                 ) { return  ( high * high - low * low ) ; }
    else if ( high < low ) 
    { return -_chebyshevU_int_ ( N ,  high , low )  ; }
    //
    const double ihigh = Ostap::Math::chebyshev_value ( N + 1 , high ) / ( N + 1 ) ;
    const double ilow  = Ostap::Math::chebyshev_value ( N + 1 , low  ) / ( N + 1 ) ;
    //
    return ihigh - ilow ;
  }
  // =========================================================================
  /** evaluate the derivative for chebychev polynomial (2nd kind)
   *  @param N the polynomial degree
   *  @param x the point 
   */
  inline double _chebyshevU_der_
  ( const unsigned int N ,
    const double       x )  
  {
    //
    if      ( 0 == N ) { return 0   ; }
    else if ( 1 == N ) { return 2   ; }
    else if ( 2 == N ) { return 8*x ; }
    //
    const unsigned long n = N + 1 ;
    if      ( s_equal ( x ,  1 ) ) { return  n * ( n * n - 1 ) / 3.0 ; }
    else if ( s_equal ( x , -1 ) ) { return  
        n * ( n * n - 1 ) / 3.0 * ( 0 == N % 2 ? 1 : -1 ) ; }
    //
    const double v1 = ( N + 1 ) * Ostap::Math::chebyshev_value  ( N + 1 , x ) ;
    const double v2 =       x   *             _chebyshevU_      ( N     , x ) ;
    const double d2 = v1 - v2 ;
    //
    return d2 / ( x * x - 1 ) ;  // ATTENTION HERE!!! it should be safe...
  }
  // ========================================================================== 
}
// ============================================================================
// evaluate the derivative of Chebyshev polynomial
// ============================================================================
double Ostap::Math::Chebyshev::derivative ( const double x ) const
{ return _chebyshev_der_ ( m_N , x ) ; }
// ============================================================================
// get integral between low and high 
// ============================================================================
double Ostap::Math::Chebyshev::integral   
( const double low  , 
  const double high ) const 
{ return _chebyshev_int_ ( m_N , low , high ) ; }
// ============================================================================
// get all roots   of the polynomial 
// ============================================================================
std::vector<double> Ostap::Math::Chebyshev::roots   () const 
{
  std::vector<double> rs ( m_N ) ;
  for ( unsigned int i = 0 ;  2 * i < m_N ; ++i ) 
  { 
    const unsigned int j = m_N - i - 1 ;
    if ( i == j ) { rs[i] = 0.0 ; }
    else 
    {
      rs[j] =  std::cos ( ( 2 * i + 1 ) * M_PIl / ( 2 * m_N ) ) ; 
      rs[i] = -rs[j] ;
    }
  }
  return rs ;
}
// ============================================================================
// get all extrema of the polynomial 
// ============================================================================
std::vector<double> Ostap::Math::Chebyshev::extrema () const 
{
  std::vector<double> e ;
  if  ( 1 <= m_N ) { e.resize ( m_N - 1 ) ; }
  for ( unsigned int i = 0 ;  2 * ( i + 1 ) < m_N ; ++i ) 
  { 
    const unsigned int j = m_N - i - 2 ;    
    if  ( i == j ) { e[i] = 0.0 ; }
    else 
    {
      e[j] = std::cos ( ( i + 1 ) * M_PIl / m_N ) ; 
      e[i] = -e[j] ;  
    }
  }
  return e ;
}
// ============================================================================
// evaluate Chebyshev polynomial of second kind 
// ============================================================================
double Ostap::Math::ChebyshevU::operator() ( const double x ) const
{ return _chebyshevU_ ( m_N , x ) ; }
// ============================================================================
// evaluate the derivative of Chebyshev polynomial of the second kind 
// ============================================================================
double Ostap::Math::ChebyshevU::derivative ( const double x ) const
{ return _chebyshevU_der_ ( m_N , x ) ; }
// ============================================================================
// get the integral between low and high 
// ============================================================================
double Ostap::Math::ChebyshevU::integral   
( const double low  , 
  const double high ) const 
{ return _chebyshevU_int_ ( m_N , low , high ) ; }
// ============================================================================
namespace 
{
  // =========================================================================
  // Legendre polynomial
  // =========================================================================
  /** evaluate the integral for legendre polynomial
   *  @param N    the polynomial degree 
   *  @param low  low edge of integration 
   *  @param high high edge of integration 
   */  
  inline long double _legendre_int_ 
  ( const unsigned int N    ,
    const long double  low  , 
    const long double  high ) 
  {
    // trivial cases 
    if      ( s_equal ( low , high ) ) { return 0          ; }
    else if ( 0 == N                 ) { return high - low ; }
    else if ( 1 == N                 )
    { return 0.5 * ( high * high - low * low ) ; }
    else if ( high < low ) 
    { return -_legendre_int_ ( N ,  high , low )  ; }
    //
    return Ostap::Math::legendre_integral ( N , low , high ) ;
  }
  // ========================================================================== 
  /** evaluate the derivative for legendre polynomial
   *  @param N the polynomial degree
   *  @param x the point 
   */
  inline long double _legendre_der_
  ( const unsigned int N ,
    const long double  x )  
  {
    //
    if      ( 0 == N ) { return 0   ; }
    else if ( 1 == N ) { return 1   ; }
    //
    if      ( s_equal ( x ,  1 ) ) { return 0.5 * N * ( N + 1 ) ; }
    else if ( s_equal ( x , -1 ) ) 
    { return 1 == N % 2 ? 0.5 * N * ( N + 1 ) : -0.5 * N * ( N + 1 ) ; }
    //
    const long double t1 = 
      x * Ostap::Math::legendre_value ( N , x ) - Ostap::Math::legendre_value ( N - 1 , x  ) ;
    //
    return N * t1 / ( x * x - 1 ) ;  // ATTENTION HERE!!! it should be safe...
  }
  // ==========================================================================
}
// ============================================================================/
// evaluate the derivative of Chebyshev polynomial
// ============================================================================
double Ostap::Math::Legendre::derivative ( const double x ) const
{ return _legendre_der_ ( m_N , x ) ; }
// ============================================================================
// get integral between low and high 
// ============================================================================
double Ostap::Math::Legendre::integral   
( const double low  , 
  const double high ) const 
{ return _legendre_int_ ( m_N , low , high ) ; }
// ============================================================================
std::map<unsigned short,std::vector<double> > s_rootmap {} ;
// ============================================================================
/** Calculate the k-th root of Legendre polynomial of order n
 *  @param k root number
 *  @param n legendre polynomial order 
 *  @return k-th root of Legendre polynomial of order n
 */
// ============================================================================
double Ostap::Math::legendre_root ( const unsigned short k , 
                                    const unsigned short n ) 
{ return Legendre ( n ).root ( k ) ; }
// ============================================================================
// get the root of the Legendre polynomial
// ============================================================================
double Ostap::Math::Legendre::calculate_root ( const unsigned short i ) const 
{
  if ( 0 == m_N ) { return -1000 ; } // convention... 
  //
  const unsigned short ii = i % m_N ;
  static const unsigned short s_maxiter = 500 ;
  /// the the first approximation 
  long double r = - std::cos ( ( 4 * ii + 3  ) * M_PIl / ( 4 * m_N + 2 ) ) ;
  /// Newton's iterations
  for ( unsigned short j = 0 ; j < s_maxiter + 1 ; ++j ) 
  {
    long double dr = evaluate ( r ) * 1.0L / derivative ( r ) ;
    if ( s_equal ( r , r - dr ) ) { break ; }                       // BREAK
    r -= dr ;
  }
  return r ;
}
// ============================================================================
/// get the root of the Legendre polynomial
double Ostap::Math::Legendre::root  ( const unsigned short i ) const 
{
  if ( 0 == m_N ) { return -1000 ; } // signal bad situation
  return  roots() [ i % m_N ] ;
}
// ============================================================================
// get the roots of the Legendre polynomial
// ============================================================================
const std::vector<double>& Ostap::Math::Legendre::roots () const 
{
  typedef std::map<unsigned int, std::vector<double> >  Roots ;
  static Roots s_rootmap {}  ;
  //
  Roots::const_iterator i = s_rootmap.find ( m_N ) ;
  if (  s_rootmap.end() != i ) { return i->second ; }   //  RETURN
  //
  // calculate and insert the entry into the map
  std::vector<double> rs ( m_N , 0 ) ;
  static const unsigned short s_maxiter = 200 ;
  for ( unsigned int i = 0 ; 2 * i  < m_N ; ++i ) 
  {
    const double ri    = calculate_root ( i ) ;
    rs [           i ] =   ri ;
    rs [ m_N - 1 - i ] = - ri ;
  }
  //
  return s_rootmap.insert ( Roots::value_type ( m_N , rs ) ).first->second ;
}
// ============================================================================
/*  update  the Legendre expansion by addition of one "event" with 
 *  the given weight
 *  @code
 *  LegendreSum sum = ... ;
 *  for ( auto x : .... ) { sum.fill ( x ) ; }
 *  @endcode
 */
// ============================================================================
bool 
Ostap::Math::LegendreSum::fill 
( const double x      , 
  const double weight ) 
{
  // no update 
  if ( x < m_xmin || x > m_xmax ) { return false ; }
  else if ( s_zero ( weight )   ) { return true  ; }
  // 
  const double tt  =  t ( x ) ;
  //
  const long double w = weight * 2.0L / ( m_xmax - m_xmin ) ;
  //
  const unsigned short N = degree() ;
  //
  m_pars[0] += w * 0.5L         ;
  if ( 0 == N ) { return true ; } //  RETURN 
  //
  m_pars[1] += w * tt * 1.5L ;
  if ( 1 == N ) { return true ; } //  RETURN 
  //
  long double p0  = 1  ;
  long double p1  = tt ;
  long double p_i = 0  ;
  //
  for ( unsigned short i = 2 ; i <= N ; ++i ) 
  {
    p_i        = ( ( 2 * i - 1 ) * tt * p1  - ( i - 1 ) * p0 ) / i ;
    m_pars[i] += w * p_i * ( i + 0.5L ) ;
    p0         = p1  ;
    p1         = p_i ;
  }
  //
  return true ;
}
// ============================================================================
/// Associated Legendre polynomials/functions
// ============================================================================
Ostap::Math::PLegendre::PLegendre 
( const unsigned int L , 
  const unsigned int M ) 
  : m_L ( L ) 
  , m_M ( L )
{}
// ============================================================================
namespace 
{
  // ==========================================================================
  //  Hermite polynomial
  // ==========================================================================
  /** evaluate the integral for Hermite polynomial
   *  @param N    the polynomial degree 
   *  @param low  low edge of integration 
   *  @param high high edge of integration 
   */  
  inline long double _hermite_int_ 
  ( const unsigned int N    ,
    const long double  low  , 
    const long double  high ) 
  {
    // trivial cases 
    if      ( s_equal ( low , high ) ) { return 0          ; }
    else if ( 0 == N                 ) { return high - low ; }
    else if ( 1 == N                 ) { return 0.5 * ( high * high - low * low ) ; }
    else if ( high < low             ) { return -_hermite_int_ ( N ,  high , low )  ; }
    //
    const long double ihigh = Ostap::Math::hermite_value ( N + 1 , high ) ;
    const long double ilow  = Ostap::Math::hermite_value ( N + 1 , low  ) ;
    //
    return ( ihigh - ilow ) / ( N + 1 ) ;
  }
  // ==========================================================================
} 
// ============================================================================
// constructor
// ============================================================================
Ostap::Math::Hermite::Hermite ( const unsigned int N ) : m_N ( N ) {} 
// ============================================================================
// get integral between low and high 
// ============================================================================
double Ostap::Math::Hermite::integral    
( const double low  , 
  const double high ) const 
{ return  _hermite_int_ ( m_N , low ,  high ) ; }
// ============================================================================


// ============================================================================
// simple  manipulations with parameters: scale it! 
// ============================================================================
// Ostap::Math::Parameters&
// Ostap::Math::Parameters::operator*=( const double a ) 
// {
//  Ostap::Math::scale ( m_pars , a ) ;
//  return *this ;
// }
// ============================================================================
// simple  manipulations with parameters scale it! 
// ============================================================================
// Ostap::Math::Parameters&
// Ostap::Math::Parameters::operator/=( const double a ) 
// {
//  Ostap::Math::scale ( m_pars , 1.0L/a ) ;
//  return *this ;
// } 
// ============================================================================
// Base class for all polynomial sums 
// ============================================================================
// constructor from polynomial degree
// ============================================================================
Ostap::Math::PolySum::PolySum ( const unsigned short degree ) 
  : Parameters ( degree + 1 ) 
{ if ( m_pars.empty() ) { m_pars.push_back ( 0 ) ; } }
// ============================================================================
// constructor from vector of parameters 
// ============================================================================
Ostap::Math::PolySum::PolySum ( const std::vector<double>& pars ) 
  : Parameters ( pars ) 
{ if ( m_pars.empty() ) { m_pars.push_back ( 0 ) ; } }
// ============================================================================
// constructor from vector of parameters 
// ============================================================================
Ostap::Math::PolySum::PolySum ( std::vector<double>&& pars ) 
  : Parameters ( std::forward<std::vector<double> >( pars ) )
{ if ( m_pars.empty() ) { m_pars.push_back ( 0 ) ; } }
// ============================================================================
// simple  manipulations with polynoms: scale it! 
// ============================================================================
// Ostap::Math::PolySum&
// Ostap::Math::PolySum::operator*=( const double a ) 
// {
//  Parameters::operator*=( a ) ;
//  return *this ;
// }
// ============================================================================
// simple  manipulations with polynoms: scale it! 
// ============================================================================
// Ostap::Math::PolySum&
// Ostap::Math::PolySum::operator/=( const double a ) 
// {
//  Parameters::operator/=( a ) ;
//  return *this ;
// }
// ============================================================================
/* Clenshaw algorithm for summation of Chebyshev polynomials 
 *  \f$ f(x) = \sum_i p_i T_i(x)\f$
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-02-10
 */
// ============================================================================
double Ostap::Math::clenshaw_chebyshev 
( const std::vector<double>& pars , 
  const double               x    ) 
{ return Ostap::Math::Clenshaw::chebyshev_sum ( pars.begin() , pars.end() , x ) ; }
// ============================================================================
/*  Clenshaw algorithm for summation of Legendre polynomials 
 *  \f$ f(x) = \sum_i p_i P_i(x) \f$
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-02-10
 */
// ============================================================================
double Ostap::Math::clenshaw_legendre
( const std::vector<double>& pars , 
  const double               x    ) 
{ return Ostap::Math::Clenshaw::legendre_sum( pars.begin() , pars.end () , x ) ; }
// ============================================================================
/*  Clenshaw algorithm for summation of Hermite series 
 *  \f$ f(x) = \sum_i p_i He_i(x) \f$
 *  @attention here we consider "probabilistic" polynomials
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-02-10
 */
// ============================================================================
double Ostap::Math::clenshaw_hermite
( const std::vector<double>& pars , 
  const double               x    ) 
{ return Ostap::Math::Clenshaw::hermite_sum( pars.begin() , pars.end () , x ) ; }
// ============================================================================
/*  Clenshaw algorithm for summation of monomial series 
 *  (aka Horner rule) 
 *  \f$ f(x) = \sum_i p_i x^i \f$
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-02-10
 */
// ============================================================================
double Ostap::Math::clenshaw_polynom
( const std::vector<double>& pars , 
  const double               x    ) 
{ return Ostap::Math::Clenshaw::monomial_sum ( pars.rbegin() , pars.rend() , x ).first ; }
// ============================================================================
/*  Clenshaw algorithm for summation of monomial series (aka Horner rule) 
 *  \f$ f(x) = \sum_i p_i x^i \f$
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-02-10
 */
// ============================================================================
double Ostap::Math::horner_a0 
( const std::vector<double>& pars , 
  const double               x    ) 
{ return Ostap::Math::Clenshaw::monomial_sum ( pars.rbegin() , pars.rend() , x ).first ; }
// ============================================================================
/*  Clenshaw algorithm for summation of monomial series (aka Horner rule) 
 *  \f$ f(x) = \sum_i p_i x^{n-i} \f$
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-02-10
 */
// ============================================================================
double Ostap::Math::horner_aN 
( const std::vector<double>& pars , 
  const double               x    ) 
{ return Ostap::Math::Clenshaw::monomial_sum ( pars.begin() , pars.end() , x ).first ; }
// ============================================================================
/*  Clenshaw algorithm for summation of cosine-series 
 *  \f$ f(x) = \frac{a_0}{2} + \sum_{i=k}^{n} a_k \cos( k x) \f$
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-02-10
 */
// ============================================================================
double Ostap::Math::clenshaw_cosine
( const std::vector<double>& pars , 
  const double               x    ) 
{ return Ostap::Math::Clenshaw::cosine_sum ( pars.begin() , pars.end() , x ) ; }
// ============================================================================
/*  Clenshaw algorithm for summation of sine-series 
 *  \f$ f(x) = \sum_{i=k}^{n} a_k \sin( k x) \f$
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-02-10
 */
// ============================================================================
double Ostap::Math::clenshaw_sine
( const std::vector<double>& pars , 
  const double               x    ) 
{ return Ostap::Math::Clenshaw::sine_sum ( pars.begin() , pars.end() , x ) ; }
// ============================================================================
/*  Clenshaw algorithm for summation of Fourier-series 
 *  \f$ f(x) = \frac{a_0}{2} + \sum_{i=k}^{n} a_{2k-1}\sin(kx)+a_{2k}\cos(kx) \f$
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-02-10
 */
// ============================================================================
double Ostap::Math::clenshaw_fourier
( const std::vector<double>& pars , 
  const double               x    ) 
{ return Ostap::Math::Clenshaw::fourier_sum ( pars.begin() , pars.end() , x ) ; }
// ============================================================================
namespace 
{
  // ==========================================================================
  /// affine tranformation of polynomial coeffficients  
  inline double _affine_ 
  ( const unsigned short j , 
    const unsigned short k , 
    const long double    a , 
    const long double    b ) 
  {
    return 
      k < j ? 0.0 : 
      Ostap::Math::choose ( k , j     ) * 
      Ostap::Math::POW    ( a , j     ) *  
      Ostap::Math::POW    ( b , k - j ) ;
  }
  // ==========================================================================
}
// ============================================================================
/* affine transformation of polynomial
 *  \f$ x ^{\prime} = \alpha x + \beta \f$
 *  @param input  (INPUT)  input polynomial coeffeicients 
 *  @param output (UPDATE) coefficinects of transformed polynomial 
 *  @param alpha  (INPUT)  alpha
 *  @param beta   (INPUT)  beta
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-03-09
 */
// ============================================================================
bool Ostap::Math::affine_transform
( const std::vector<double>& input  , 
  std::vector<double>&       result ,  
  const double               alpha  , 
  const double               beta   ) 
{
  /// invalid transform 
  if      ( s_zero  ( alpha     ) ) { return false ; } // invalid transform 
  else if ( s_equal ( alpha , 1 ) && s_zero ( beta ) ) 
  { result = input ;                  return true  ; } // trivial transformation 
  //
  result.resize ( input.size() ) ;
  std::fill ( result.begin() , result.end () , 0 ) ;
  //
  for ( unsigned int i = 0 ; i < input.size() ; ++i ) 
  { 
    for ( unsigned short k = i ; k < input.size () ; ++k ) // ATTENTNION!!! 
    { 
      const double p = input[  k ] ;
      if ( s_zero ( p ) ) { continue ; }  // SKIP nulls... 
      result[i] += _affine_ ( i , k , alpha , beta ) * p ; 
    } 
  }  
  return true ;
}
// ============================================================================




// ============================================================================
namespace 
{
  // ==========================================================================
  inline void add_root 
  ( std::vector<double>& pars , 
    const double         root ) 
  {
    //
    if ( pars.empty() ) 
    {
      pars.push_back ( -root ) ;
      pars.push_back ( 1.0   ) ;
      return ;
    }
    //
    pars.push_back ( 0 ) ;
    const std::size_t N = pars.size() ;
    for ( std::size_t k = N - 1  ;  1 <= k  ; --k )
    { pars [ k ] = pars [ k - 1 ] - pars [ k ] * root ; }
    //
    pars [ 0 ] *= -root ;
  }
  // ==========================================================================
  inline void add_roots
  ( std::vector<double>&       pars , 
    const std::complex<double> root ) 
  {
    //
    const double alpha = -2 * root.real()    ;
    const double beta  = std::norm ( root ) ;
    //
    if ( pars.empty() ) 
    {
      pars.push_back ( beta  ) ;
      pars.push_back ( alpha ) ;
      pars.push_back ( 1     ) ;
      return ;
    }
    //
    pars.push_back ( 0 ) ;
    pars.push_back ( 0 ) ;
    //
    const std::size_t N = pars.size() ;
    for ( std::size_t k = N - 1  ;  2 <= k  ; --k )
    { pars [ k ] = 
        pars [ k     ] * beta   + 
        pars [ k - 1 ] * alpha  + 
        pars [ k - 2 ] ; 
    }
    //
    pars [ 1 ]  = beta * pars[1] + alpha * pars[0] ;
    pars [ 0 ] *= beta  ;
  }
  // ==========================================================================
} 
// ============================================================================
/*  class Polynomial
 *  Trivial polynomial
 *  \f$ f(x) = \sum_i \p_i x^i\f$
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-02-22
 */
// ============================================================================
// constructor from the degree  
// ============================================================================
Ostap::Math::Polynomial::Polynomial
( const unsigned short       degree , 
  const double               xmin   , 
  const double               xmax   )  
  : Ostap::Math::PolySum ( degree ) 
  , m_xmin ( std::min ( xmin, xmax ) )
  , m_xmax ( std::max ( xmin, xmax ) )
{}
// ============================================================================
// constructor from parameters 
// ============================================================================
Ostap::Math::Polynomial::Polynomial
( const std::vector<double>& pars  , 
  const double               xmin  , 
  const double               xmax  )  
  : Ostap::Math::PolySum ( pars ) 
  , m_xmin ( std::min ( xmin, xmax ) )
  , m_xmax ( std::max ( xmin, xmax ) )
{}
// ============================================================================
// construct polynomial from different range
// ============================================================================
Ostap::Math::Polynomial::Polynomial
( const Polynomial& right ,
  const double      xmin  ,
  const double      xmax  )
  : Ostap::Math::PolySum ( right )
  , m_xmin ( std::min ( xmin, xmax ) )
  , m_xmax ( std::max ( xmin, xmax ) )
{
  if ( !s_equal ( m_xmin , right.m_xmin ) ||
       !s_equal ( m_xmax , right.m_xmax )  ) 
  {
    const double len1  =       m_xmax -       m_xmin ;
    const double len2  = right.m_xmax - right.m_xmin ;
    //
    const double scale = len1 / len2 ;
    const double shift = 2 * ( xmid () - right.xmid () ) / len2 ;
    affine_transform ( right.m_pars , 
                       m_pars       ,
                       scale        , 
                       shift        ) ;
  }
}
// ============================================================================
/* construct polynomial from its roots
 *
 *  Polinomial has a form
 *  \f$ B(x) = \prod_i (x-r_i) \prod_j (x-c_j)(x-c_j^*) \f$
 *
 *  @param xmin low  edge for polynomial
 *  @param xmax high edge for polynomial
 *  @param roots_real    the list of real  roots of the polinomial
 *  @param roots_complex the list of complex roots (only one root from cc-pair is needed)
 */
// ============================================================================
Ostap::Math::Polynomial::Polynomial
( const double xmin , 
  const double xmax , 
  const std::vector<double>&                 roots_real    ,     
  const std::vector<std::complex<double> > & roots_complex )
  : Ostap::Math::PolySum ( roots_real.size() + 2 * roots_complex.size () )
  , m_xmin ( std::min ( xmin, xmax ) )
  , m_xmax ( std::max ( xmin, xmax ) )
{
  std::vector<double> coeffs ( 1 , 1.0 )  ;
  coeffs.reserve ( degree() + 1 ) ;
  /// add real roots 
  for ( std::vector<double>::const_iterator rr = roots_real.begin() ; 
        roots_real.end()    != rr ; ++rr ) 
  { add_root ( coeffs , t ( *rr ) ) ; }
  /// add complex roots 
  for ( std::vector<std::complex<double>>::const_iterator cr = roots_complex.begin() ; 
        roots_complex.end() != cr ; ++cr ) 
  { 
    const double a = t ( cr->real () ) ;
    const double b = t ( cr->imag () ) ;
    add_roots ( coeffs , std::complex<double>( a , b ) ) ;
  }
  //
  m_pars = coeffs ;
}
// ============================================================================
/** construct polynomial from its roots
 *
 *  Polinomial has a form
 *  \f$ B(x) = \prod_i (x-r_i) \prod_j (x-c_j)(x-c_j^*) \f$
 *
 *  @param xmin low  edge for polynomial
 *  @param xmax high edge polynomial
 *  @param roots_complex the list of complex roots (only one root from cc-pair is needed)
 *  @param roots_real    the list of real  roots of the polinomial
 */
// ============================================================================
Ostap::Math::Polynomial::Polynomial
( const double xmin , 
  const double xmax , 
  const std::vector<std::complex<double> > & roots_complex ,
  const std::vector<double>&                 roots_real    ) 
  : Polynomial ( xmin , xmax , roots_real , roots_complex )
{}
// ============================================================================
// constructor from Karlin-Shapley polinomial
// ============================================================================
Ostap::Math::Polynomial::Polynomial
( const Ostap::Math::KarlinShapley& poly ) 
  : PolySum ( poly.degree() ) 
  , m_xmin  ( poly.xmin () )
  , m_xmax  ( poly.xmax () )
{
  const double alpha = poly.alpha () ;
  const double beta  = poly.beta  () ;
  if      ( s_zero ( poly.A() ) ) {}  
  else if ( 0 == degree()       ) { m_pars [0] = alpha ; }
  else if ( 1 == degree()       ) 
  {
    m_pars [ 0 ] = ( alpha + beta ) * 0.5 ;   
    m_pars [ 1 ] = ( alpha - beta ) * 0.5 ;
  }
  else
  {
    const double even =  0 == degree() % 2 ;
    const double odd  = !even ;
    //
    const unsigned int NT = degree() + 1 ;
    //
    if ( !s_zero ( alpha ) ) 
    {
      std::vector<double> aroots {} ; aroots.reserve ( degree() ) ;
      for ( unsigned k = even ? 1: 2  ; k < NT ; k+= 2 ) 
      {
        const double r = poly.x ( poly.troot ( k ) ) ;
        aroots.push_back ( r ) ;
        aroots.push_back ( r ) ;
      }
      //
      if ( odd ) { aroots.push_back ( m_xmin ) ; }
      Polynomial tmpa ( m_xmin , m_xmax , aroots ) ;
      tmpa *= alpha / std::pow ( 2.0 , degree () ) ;;
      isum ( tmpa ) ;
    }
    if ( !s_zero ( beta ) ) 
    {
      std::vector<double> broots {} ; broots.reserve ( degree() ) ;
      //
      for ( unsigned k = even ? 2 : 1  ; k + 1 < NT ; k+= 2 ) 
      {
        const double r = poly.x ( poly.troot ( k ) ) ;
        broots.push_back ( r ) ;
        broots.push_back ( r ) ;
      }
      if ( even ) { broots.push_back ( m_xmin ) ;  broots.push_back ( m_xmax ) ; }
      else        { broots.push_back ( m_xmax ) ; }
      //
      Polynomial tmpb ( m_xmin , m_xmax , broots ) ;
      tmpb *= -beta / std::pow ( 2.0 , degree () ) ;
      isum ( tmpb ) ;
    }
  }
}
// ============================================================================
// constructor from Karlin-Studden polinomial
// ============================================================================
Ostap::Math::Polynomial::Polynomial
( const Ostap::Math::KarlinStudden& poly ) 
  : PolySum ( poly.degree() ) 
  , m_xmin  ( poly.xmin ()                )
  , m_xmax  ( poly.xmin () + poly.scale() )
{
  const double alpha = poly.alpha () ;
  const double beta  = poly.beta  () ;
  if      ( s_zero ( poly.A() ) ) {}  
  else if ( 1 == degree () ) 
  {
    m_pars [ 0 ] = alpha * 0.5 + beta ;   
    m_pars [ 1 ] = alpha * 0.5  ;
  }
  else
  {
    const double even =  0 == degree() % 2 ;
    const double odd  = !even ;
    //
    const unsigned int NT = degree() + 1 ;
    //
    if ( !s_zero ( alpha ) ) 
    {
      std::vector<double> aroots {} ; aroots.reserve ( degree() ) ;
      for ( unsigned k = even ? 1: 2  ; k < NT ; k+= 2 ) 
      {
        const double r = poly.x ( poly.zroot ( k ) ) ;
        aroots.push_back ( r ) ;
        aroots.push_back ( r ) ;
      }
      //
      if ( odd ) { aroots.push_back ( m_xmin ) ; }
      Polynomial tmpa ( m_xmin , m_xmax , aroots ) ;
      tmpa *= alpha / std::pow ( 2.0 , aroots.size() ) ;;
      isum ( tmpa ) ;
    }
    if ( !s_zero ( beta ) ) 
    {
      std::vector<double> broots {} ; broots.reserve ( degree() ) ;
      //
      for ( unsigned k = even ? 2 : 1  ; k + 1 < NT ; k+= 2 ) 
      {
        const double r = poly.x ( poly.zroot ( k ) ) ;
        broots.push_back ( r ) ;
        broots.push_back ( r ) ;
      }
      if ( even ) { broots.push_back ( m_xmin ) ; }
      //
      Polynomial tmpb ( m_xmin , m_xmax , broots ) ;
      tmpb *= beta / std::pow ( 2.0 , broots.size ()  ) ;
      isum ( tmpb ) ;
    }
  }
}
// ============================================================================
//  constructor from Even Bernstein polynomial (efficient) 
// ============================================================================
Ostap::Math::Polynomial::Polynomial 
( const Ostap::Math::BernsteinEven& poly ) 
  : Polynomial ( poly.bernstein() ) 
{}
// ============================================================================
//  constructor from Positive Bernstein polynomial (efficient) 
// ============================================================================
Ostap::Math::Polynomial::Polynomial 
( const Ostap::Math::Positive& poly ) 
  : Polynomial ( poly.bernstein() ) 
{}
// ============================================================================
//  constructor from Monotonic Bernstein polynomial (efficient) 
// ============================================================================
Ostap::Math::Polynomial::Polynomial 
( const Ostap::Math::Monotonic& poly ) 
  : Polynomial ( poly.bernstein() ) 
{}
// ============================================================================
//  constructor from Convex Bernstein polynomial (efficient) 
// ============================================================================
Ostap::Math::Polynomial::Polynomial 
( const Ostap::Math::Convex& poly ) 
  : Polynomial ( poly.bernstein() ) 
{}
// ============================================================================
//  constructor from Convex-only Bernstein polynomial (efficient) 
// ============================================================================
Ostap::Math::Polynomial::Polynomial 
( const Ostap::Math::ConvexOnly& poly ) 
  : Polynomial ( poly.bernstein() ) 
{}
// ============================================================================
// get the value
// ============================================================================
double Ostap::Math::Polynomial::evaluate ( const double x ) const 
{
  // trivial cases 
  if      ( 1 == m_pars.size()       ) { return m_pars[0] ; }  // RETURN 
  else if ( zero ()                  ) { return 0         ; }
  // 1) transform argument:
  const double tx = t ( x ) ;  
  // 2) use Clenshaw's algorithm 
  return clenshaw_polynom ( m_pars , tx ) ;              // RETURN 
}
// ============================================================================
double Ostap::Math::Polynomial::integral   
( const double low  , 
  const double high ) const 
{
  //
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high < low             ) { return -integral ( high , low ) ; }
  else if ( high <= m_xmin         ) { return 0 ; }
  else if ( low  >= m_xmax         ) { return 0 ; }
  else if ( zero ()                ) { return 0 ; }
  //
  if      ( low  <  m_xmin         ) { return integral ( m_xmin , high   ) ; }
  else if ( high >  m_xmax         ) { return integral ( low    , m_xmax ) ; }
  //
  const double xl = t ( low  ) ;
  const double xh = t ( high ) ;
  //
  std::vector<double> npars ( m_pars.size() + 1 , 0 ) ;
  for ( unsigned int i = 0 ; i < m_pars.size() ; ++i ) 
  {
    if ( s_zero ( m_pars[i] ) ) { continue ; }           // CONTINUE
    npars [i+1] = m_pars[i] / ( i + 1 ) ; 
  }
  //
  const double result = 
    clenshaw_polynom ( npars , xh ) -
    clenshaw_polynom ( npars , xl ) ;     
  //
  return result  * 0.5 * ( m_xmax - m_xmin ) ;
}
// ============================================================================
// get indefinte integral 
// ============================================================================
Ostap::Math::Polynomial 
Ostap::Math::Polynomial::indefinite_integral ( const double C ) const 
{
  const double dx = 0.5 * ( m_xmax - m_xmin ) ;
  //
  Ostap::Math::Polynomial integ ( degree() + 1 , m_xmin , m_xmax ) ;
  //
  for ( unsigned int i = 0 ; i < m_pars.size() ; ++i ) 
  { 
    if ( s_zero ( m_pars[i] ) ) { continue ; }           // CONTINUE
    integ.m_pars [ i + 1 ] = m_pars [ i ] / ( i + 1 ) * dx ; 
  }
  //
  integ.m_pars[0] = C ;
  //
  return integ ;
}
// ============================================================================
// get the derivative at point "x" 
// ============================================================================
double Ostap::Math::Polynomial::derivative ( const double x     ) const 
{
  if ( x < m_xmin || x > m_xmax ) { return 0 ; }
  //
  const double tx = t ( x  ) ;
  //
  const double dx = 2 / ( m_xmax - m_xmin ) ;
  std::vector<double> npars ( m_pars.size() - 1 , 0 ) ;
  for ( unsigned int i = 0 ; i < npars.size() ; ++i ) 
  {
    const double p = m_pars[i+1] ;
    if ( s_zero ( p ) ) { continue ; }      // CONTINUE 
    npars [i] = ( i + 1 ) * p * dx ; 
  }
  //
  return clenshaw_polynom ( npars , tx ) ;
}
// ============================================================================
// get the derivative 
// ============================================================================
Ostap::Math::Polynomial 
Ostap::Math::Polynomial::derivative          () const 
{
  if ( 1 == m_pars.size() ) 
  { return Ostap::Math::Polynomial( 0 , m_xmin , m_xmax ) ; }
  //
  const double dx = 2 / ( m_xmax - m_xmin ) ;
  //
  Ostap::Math::Polynomial deriv ( degree() - 1 , m_xmin , m_xmax ) ;
  for ( unsigned int i = 0 ; i < deriv.npars () ; ++i ) 
  {    
    const double p = m_pars [ i + 1 ] ;
    if ( s_zero ( p ) ) { continue ; }       // CONTINUE
    deriv.m_pars [ i ] = ( i + 1 ) * p * dx ; 
  }
  //
  return deriv  ;
}
// ============================================================================ 
// simple  manipulations with polynoms: shift it! 
// ============================================================================
Ostap::Math::Polynomial& 
Ostap::Math::Polynomial::operator+= ( const double a ) 
{ m_pars[0] += a ; return *this ; }
// ============================================================================
// simple  manipulations with polynoms: shift it! 
// ============================================================================
Ostap::Math::Polynomial& 
Ostap::Math::Polynomial::operator-= ( const double a ) 
{ m_pars[0] -= a ; return *this ; }
// ============================================================================ 
// simple  manipulations with polynoms: scale it 
// ============================================================================
Ostap::Math::Polynomial& 
Ostap::Math::Polynomial::operator*= ( const double a ) 
{ 
  Ostap::Math::scale ( m_pars , a ) ;
  return *this ; 
}
// ============================================================================
// simple  manipulations with polynoms: shift it! 
// ============================================================================
Ostap::Math::Polynomial& 
Ostap::Math::Polynomial::operator/= ( const double a ) 
{ 
  Ostap::Math::scale ( m_pars , 1/a ) ;
  return *this ; 
}
// ============================================================================
// Add       polynomials (with the same domain!)
// ============================================================================
Ostap::Math::Polynomial&
Ostap::Math::Polynomial::isum
( const Ostap::Math::Polynomial& other ) 
{
  if ( this == &other ) { *this *= 2 ; return *this; }
  //
  Ostap::Assert ( s_equal ( xmin() , other.xmin() ) &&
                  s_equal ( xmax() , other.xmax() )               ,
                  "Cannot sum Polynomials with different domains" , 
                  "Ostap::Math::Polynomial"                       , 
                  Ostap::StatusCode ( 520 )                       )  ;
  //
  const unsigned short idegree = std::max ( degree() , other.degree() ) ;
  m_pars.resize ( idegree + 1 ) ;
  //
  for ( unsigned short i = 0 ; i < other.npars() ; ++i ) 
  { m_pars[i] += other.m_pars [ i ] ; }
  //
  return *this ;
}
// ============================================================================
// Add       polynomials (with the same domain!)
// ============================================================================
Ostap::Math::Polynomial 
Ostap::Math::Polynomial::sum
( const Ostap::Math::Polynomial& other ) const 
{
  Polynomial result(*this) ;
  result.isum ( other ) ;
  return result ;
}
// ============================================================================
// Subtract      polynomials (with the same domain!)
// ============================================================================
Ostap::Math::Polynomial&
Ostap::Math::Polynomial::isub
( const Ostap::Math::Polynomial& other ) 
{
  if ( this == &other ) { *this *= 0 ; return *this ; }
  //
  Ostap::Assert ( s_equal ( xmin() , other.xmin() ) &&
                  s_equal ( xmax() , other.xmax() )              ,
                  "Cannot subtract Polynomials with different domains" , 
                  "Ostap::Math::Polynomial"                      , 
                  Ostap::StatusCode ( 521 )                      )  ;
  //
  const unsigned short idegree = std::max ( degree() , other.degree() ) ;
  m_pars.resize ( idegree + 1 ) ;
  //
  for ( unsigned short i = 0 ; i < other.npars() ; ++i ) 
  { m_pars[i] -= other.m_pars [ i ] ; }
  //
  return *this ;
}
// ============================================================================
// Subtract      polynomials (with the same domain!)
// ============================================================================
Ostap::Math::Polynomial 
Ostap::Math::Polynomial::subtract
( const Ostap::Math::Polynomial& other ) const 
{
  Polynomial result(*this) ;
  result.isub ( other ) ;
  return result ;
}
// ============================================================================
// unary minus 
// ============================================================================
Ostap::Math::Polynomial
Ostap::Math::Polynomial::operator-() const 
{
  Polynomial a ( *this ) ;
  Ostap::Math::negate ( a.m_pars ) ;
  return a ;
}
// ============================================================================
// Python
// ============================================================================
Ostap::Math::Polynomial& 
Ostap::Math::Polynomial::__iadd__   ( const double a )  
{ (*this) += a ; return *this ; } 
// ============================================================================
Ostap::Math::Polynomial& 
Ostap::Math::Polynomial::__isub__   ( const double a ) 
{ (*this) -= a ; return *this ; } 
// ============================================================================
Ostap::Math::Polynomial& 
Ostap::Math::Polynomial::__imul__   ( const double a )
{ (*this) *= a ; return *this ; } 
// ============================================================================
Ostap::Math::Polynomial& 
Ostap::Math::Polynomial::__itruediv__   ( const double a ) 
{ (*this) /= a ; return *this ; } 
// ============================================================================
Ostap::Math::Polynomial 
Ostap::Math::Polynomial::__add__    ( const double a ) const 
{ return Polynomial(*this) += a ; } 
// ============================================================================
Ostap::Math::Polynomial 
Ostap::Math::Polynomial::__sub__    ( const double a ) const 
{ return Polynomial(*this) -= a ; } 
// ============================================================================
Ostap::Math::Polynomial 
Ostap::Math::Polynomial::__mul__   ( const double a ) const 
{ return Polynomial(*this) *= a ; } 
// ============================================================================
Ostap::Math::Polynomial 
Ostap::Math::Polynomial::__truediv__    ( const double a ) const 
{ return Polynomial(*this) /= a ; } 
// ============================================================================
Ostap::Math::Polynomial 
Ostap::Math::Polynomial::__radd__   ( const double a ) const 
{ return __add__ ( a ) ; } 
// ============================================================================
Ostap::Math::Polynomial 
Ostap::Math::Polynomial::__rsub__   ( const double a ) const 
{ return (-(*this))+a; } 
// ============================================================================
Ostap::Math::Polynomial 
Ostap::Math::Polynomial::__neg__   () const { return -(*this); } 
// ============================================================================
Ostap::Math::Polynomial 
Ostap::Math::Polynomial::__add__   
( const Ostap::Math::Polynomial& a ) const { return sum      ( a ) ; } 
// ============================================================================
Ostap::Math::Polynomial 
Ostap::Math::Polynomial::__sub__   
( const Ostap::Math::Polynomial& a ) const { return subtract ( a ) ; } 
// ============================================================================
// get unique tag 
// ============================================================================
std::size_t Ostap::Math::Polynomial::tag () const 
{
  static const std::string s_name { "Polynomial" } ;
  return Ostap::Utils::hash_combiner 
    ( s_name , Ostap::Utils::hash_range ( m_pars ) ,  
      degree () , m_xmin , m_xmax ) ;
}
// ============================================================================



// ============================================================================
// Chebyshev sum 
// ============================================================================
// constructor from parameters 
// ============================================================================
Ostap::Math::ChebyshevSum::ChebyshevSum
( const unsigned short degree , 
  const double         xmin   , 
  const double         xmax   )  
  : Ostap::Math::PolySum ( degree ) 
  , m_xmin ( std::min ( xmin, xmax ) )
  , m_xmax ( std::max ( xmin, xmax ) )
{}
// ============================================================================
// constructor from parameters 
// ============================================================================
Ostap::Math::ChebyshevSum::ChebyshevSum
( const std::vector<double>& pars  , 
  const double               xmin  , 
  const double               xmax  )  
  : Ostap::Math::PolySum ( pars ) 
  , m_xmin ( std::min ( xmin, xmax ) )
  , m_xmax ( std::max ( xmin, xmax ) )
{}
// ============================================================================
// get the value
// ============================================================================
double Ostap::Math::ChebyshevSum::evaluate ( const double x ) const 
{
  // trivial cases 
  if      ( 1 == m_pars.size()       ) { return m_pars[0] ; }  // RETURN 
  else if ( zero ()                  ) { return 0         ; }
  // 1) transform argument:
  const double tx = t ( x ) ;  
  // 2) use Clenshaw's algorithm 
  return Ostap::Math::Clenshaw::chebyshev_sum ( m_pars.begin() , m_pars.end()  , tx ) ;        
}
// ============================================================================
// get the integral between xmin and xmax
// ============================================================================
double Ostap::Math::ChebyshevSum::integral   () const 
{ return integral ( m_xmin , m_xmax ) ; }
// ============================================================================
double Ostap::Math::ChebyshevSum::integral   
( const double low  , 
  const double high ) const 
{
  //
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high < low             ) { return -integral ( high , low ) ; }
  else if ( high <= m_xmin         ) { return 0 ; }
  else if ( low  >= m_xmax         ) { return 0 ; }
  else if ( zero ()                ) { return 0 ; }
  //
  if      ( low  <  m_xmin         ) { return integral ( m_xmin , high   ) ; }
  else if ( high >  m_xmax         ) { return integral ( low    , m_xmax ) ; }
  //
  const double xl = t ( low  ) ;
  const double xh = t ( high ) ;
  //
  const double dx = 0.5 * ( m_xmax - m_xmin ) ;
  std::vector<double> npars ( m_pars.size() + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i < m_pars.size() ; ++i ) 
  {
    if ( s_zero ( m_pars[i] ) ) { continue ; }              // CONTINUE 
    //
    if      ( 0 == i ) { npars [1] +=        m_pars[0] * dx ; }
    else if ( 1 == i ) { npars [2] += 0.25 * m_pars[1] * dx ; }
    else 
    {
      npars [ i + 1] += m_pars[i] * 0.5 / ( i + 1.0 ) * dx ;  
      npars [ i - 1] -= m_pars[i] * 0.5 / ( i - 1.0 ) * dx ; 
    }
  }
  //
  const double result = 
    Ostap::Math::Clenshaw::chebyshev_sum   ( npars.begin() , npars.end()  , xh ) - 
    Ostap::Math::Clenshaw::chebyshev_sum   ( npars.begin() , npars.end()  , xl ) ;
  //
  return result ;
}
// ============================================================================
// get indefinte integral 
// ============================================================================
Ostap::Math::ChebyshevSum 
Ostap::Math::ChebyshevSum::indefinite_integral ( const double C ) const 
{
  // 
  const double dx = 0.5 * ( m_xmax - m_xmin ) ;
  //
  Ostap::Math::ChebyshevSum integ ( degree() + 1 , m_xmin , m_xmax ) ;
  //
  for ( unsigned short i = 0 ; i < m_pars.size() ; ++i ) 
  {
    if ( s_zero ( m_pars[i] ) ) { continue ; }           // CONTINUE 
    //
    if      ( 0 == i ) { integ.m_pars [ 1 ] +=        m_pars [ 0 ] * dx ; }
    else if ( 1 == i ) { integ.m_pars [ 2 ] += 0.25 * m_pars [ 1 ] * dx ; }
    else 
    {
      integ.m_pars [ i + 1 ] += m_pars [ i ] * 0.5 / ( i + 1.0 ) * dx ;  
      integ.m_pars [ i - 1 ] -= m_pars [ i ] * 0.5 / ( i - 1.0 ) * dx ; 
    }
  }
  //
  integ.m_pars [ 0 ] += C ;
  //
  return integ ;
}
// ============================================================================
// get the derivative at point "x" 
// ============================================================================
double Ostap::Math::ChebyshevSum::derivative ( const double x     ) const 
{
  if ( x < m_xmin || x > m_xmax ) { return 0 ; }
  //
  const double tx = t ( x  ) ;
  //
  // // Trvial sum to be replaced with Clenshaw 
  // double result = 0 ;
  // for ( std::vector<double>::const_iterator ip = m_pars.begin() ; m_pars.end() != ip ; ++ip ) 
  // {
  //   const double p = *ip ;
  //   if ( s_zero ( p ) ) { continue ; } // SKIP IT! 
  //   result += p * _chebyshev_der_ ( ip - m_pars.begin() , tx ) ; 
  // }
  //
  std::vector<double> npars ( m_pars.size() - 1 , 0 ) ;
  for ( unsigned short i = 1 ; i < m_pars.size() ; ++i ) 
  { 
    //
    if ( s_zero ( m_pars[i] ) ) { continue ; }           // CONTINUE 
    //
    const unsigned short id = i - 1 ;
    if ( 0 == id % 2 ) 
    {
      for ( unsigned short j = 0 ; j <= id ; j += 2 ) 
      { npars [ j ] += 2 * i * m_pars [ i ] ; } 
      npars   [ 0 ] -=     i * m_pars [ i ] ;
    }
    else 
    {
      for ( unsigned short j = 1 ; j <= id ; j += 2 ) 
      { npars [ j ] += 2 * i * m_pars [ i ] ; } 
    }
  }
  //
  const long double result = 
    Ostap::Math::Clenshaw::chebyshev_sum  ( npars.begin() , npars.end()  , tx ) ;
  //
  const long double dx = 2.0L / ( m_xmax - m_xmin ) ;
  return result * dx ;  
}
// ============================================================================
// get the derivative 
// ============================================================================
Ostap::Math::ChebyshevSum 
Ostap::Math::ChebyshevSum::derivative () const 
{
  //
  if ( 1 == m_pars.size() ) 
  { return Ostap::Math::ChebyshevSum ( 0 , m_xmin , m_xmax ) ; }
  //
  const double dx = 2 / ( m_xmax - m_xmin ) ;
  //
  // result 
  Ostap::Math::ChebyshevSum deriv  ( degree() -1 , m_xmin , m_xmax ) ;
  for ( unsigned short i = 1 ; i < m_pars.size() ; ++i ) 
  { 
    //
    if ( s_zero ( m_pars [ i ] ) ) { continue ; }           // CONTINUE 
    //
    const unsigned short id = i - 1 ;
    if ( 0 == id % 2 ) 
    {
      for ( unsigned short j = 0 ; j <= id ; j+=2 ) 
      { deriv.m_pars [ j ] += 2 * i * m_pars [ i ] * dx ; } 
      //
      deriv.m_pars [ 0 ]   -=     i * m_pars [ i ] * dx ;
    }
    else 
    {
      for ( unsigned short j = 1 ; j <= id ; j+=2 ) 
      { deriv.m_pars [ j ] += 2 * i * m_pars [ i ] * dx ; } 
    }
  }
  //
  return deriv ;
}
// ============================================================================
/*  update  the chebyshev expansion by addition of one "event" with 
 *  the given weight
 *  @code
 *  ChebyshevSum sum = ... ;
 *  for ( auto x : .... ) { sum.fill ( x ) ; }
 *  @endcode
 */
// ============================================================================
bool 
Ostap::Math::ChebyshevSum::fill 
( const double x      , 
  const double weight ) 
{
  // no update 
  if ( x <= m_xmin || x >= m_xmax ) { return false ; }
  else if ( s_zero ( weight )     ) { return true  ; }
  // 
  const long double tt =  t ( x ) ;
  if ( 1 <= std::abs ( tt )       ) { return false ; }
  //
  const long double w  = weight * 4.0L / ( m_xmax - m_xmin ) / ( std::sqrt ( 1.0L - tt * tt ) * M_PI ) ;
  if ( !std::isfinite ( w  )     ) { return false ; }
  //
  const unsigned short N = degree() ;
  //
  m_pars [ 0 ] += w * 0.5L ;
  if ( 0 == N ) { return true ; } //  RETURN 
  //
  m_pars [ 1 ] += w * tt   ;
  if ( 1 == N ) { return true ; } //  RETURN 
  //
  long double p0  = 1  ;
  long double p1  = tt ;
  long double p_i = 0  ;
  //
  for ( unsigned short i = 2 ; i <= N ; ++i ) 
  {
    p_i           = 2.0L * tt * p1 - p0 ;
    m_pars [ i ] += w * p_i ;
    p0            = p1  ;
    p1            = p_i ;
  }
  //
  return true ;
}
// ============================================================================
// simple  manipulations with polynoms: shift it! 
// ============================================================================
Ostap::Math::ChebyshevSum& 
Ostap::Math::ChebyshevSum::operator+= ( const double a ) 
{ m_pars [ 0 ] += a ; return *this ; }
// ============================================================================
// simple  manipulations with polynoms: shift it! 
// ============================================================================
Ostap::Math::ChebyshevSum& 
Ostap::Math::ChebyshevSum::operator-= ( const double a ) 
{ m_pars [ 0 ] -= a ; return *this ; }
// ============================================================================
// simple  manipulations with polynoms: scale it  
// ============================================================================
Ostap::Math::ChebyshevSum& 
Ostap::Math::ChebyshevSum::operator*= ( const double a ) 
{
  Ostap::Math::scale ( m_pars , a ) ;
  return *this ; 
}
// ============================================================================
// simple  manipulations with polynoms: scale it 
// ============================================================================
Ostap::Math::ChebyshevSum& 
Ostap::Math::ChebyshevSum::operator/= ( const double a ) 
{
  Ostap::Math::scale ( m_pars , 1/a ) ;
  return *this ; 
}
// ============================================================================
// unary minus 
// ============================================================================
Ostap::Math::ChebyshevSum
Ostap::Math::ChebyshevSum::operator-() const 
{
  ChebyshevSum a ( *this ) ;
  Ostap::Math::negate ( a.m_pars ) ;
  return a ;
}
// ============================================================================
// Python
// ============================================================================
Ostap::Math::ChebyshevSum& 
Ostap::Math::ChebyshevSum::__iadd__   ( const double a )  
{ (*this) += a ; return *this ; } 
// ============================================================================
Ostap::Math::ChebyshevSum& 
Ostap::Math::ChebyshevSum::__isub__   ( const double a ) 
{ (*this) -= a ; return *this ; } 
// ============================================================================
Ostap::Math::ChebyshevSum& 
Ostap::Math::ChebyshevSum::__imul__   ( const double a )
{ (*this) *= a ; return *this ; } 
// ============================================================================
Ostap::Math::ChebyshevSum& 
Ostap::Math::ChebyshevSum::__itruediv__   ( const double a ) 
{ (*this) /= a ; return *this ; } 
// ============================================================================
Ostap::Math::ChebyshevSum 
Ostap::Math::ChebyshevSum::__add__    ( const double a ) const 
{ return ChebyshevSum(*this) += a ; } 
// ============================================================================
Ostap::Math::ChebyshevSum 
Ostap::Math::ChebyshevSum::__sub__    ( const double a ) const 
{ return ChebyshevSum(*this) -= a ; } 
// ============================================================================
Ostap::Math::ChebyshevSum 
Ostap::Math::ChebyshevSum::__mul__    ( const double a ) const 
{ return ChebyshevSum(*this) *= a ; } 
// ============================================================================
Ostap::Math::ChebyshevSum 
Ostap::Math::ChebyshevSum::__truediv__    ( const double a ) const 
{ return ChebyshevSum(*this) /= a ; } 
// ============================================================================
Ostap::Math::ChebyshevSum 
Ostap::Math::ChebyshevSum::__radd__   ( const double a ) const 
{ return __add__ ( a ) ; } 
// ============================================================================
Ostap::Math::ChebyshevSum 
Ostap::Math::ChebyshevSum::__rsub__   ( const double a ) const 
{ return (-(*this))+a; } 
// ============================================================================
Ostap::Math::ChebyshevSum 
Ostap::Math::ChebyshevSum::__rmul__   ( const double a ) const 
{ return __mul__  ( a ) ; } 
// ============================================================================
Ostap::Math::ChebyshevSum 
Ostap::Math::ChebyshevSum::__neg__   () const { return -(*this); } 
// ============================================================================
Ostap::Math::ChebyshevSum 
Ostap::Math::ChebyshevSum::__add__   
( const Ostap::Math::ChebyshevSum& a ) const { return sum      ( a ) ; } 
// ============================================================================
Ostap::Math::ChebyshevSum 
Ostap::Math::ChebyshevSum::__sub__   
( const Ostap::Math::ChebyshevSum& a ) const { return subtract ( a ) ; } 
// ============================================================================





// ============================================================================
// Add       polynomials (with the same domain!)
// ============================================================================
Ostap::Math::ChebyshevSum&
Ostap::Math::ChebyshevSum::isum
( const Ostap::Math::ChebyshevSum& other )
{
  if ( this == &other ) { *this *= 2 ; return *this ; }
  //
  Ostap::Assert ( s_equal ( xmin() , other.xmin() ) &&
                  s_equal ( xmax() , other.xmax() )                ,
                  "Cannot sum Chebyshev with different domains" , 
                  "Ostap::Math::Chebyshev"                         , 
                  Ostap::StatusCode ( 522 )                        )  ;
  //
  const unsigned short idegree = std::max ( degree() , other.degree() ) ;
  //
  for ( unsigned short i = 0 ; i < other.npars() ; ++i ) 
  { m_pars[i] += other.m_pars [ i ] ; }
  //
  return *this ;
}
// ============================================================================
// Add       polynomials (with the same domain!)
// ============================================================================
Ostap::Math::ChebyshevSum 
Ostap::Math::ChebyshevSum::sum
( const Ostap::Math::ChebyshevSum& other ) const 
{
  ChebyshevSum result(*this) ;
  result.isum ( other ) ;
  return result ;
}
// ============================================================================
// Subtract      polynomials (with the same domain!)
// ============================================================================
Ostap::Math::ChebyshevSum&
Ostap::Math::ChebyshevSum::isub
( const Ostap::Math::ChebyshevSum& other )
{
  if ( this == &other ) { *this *= 0 ; return *this ; }
  //
  Ostap::Assert ( s_equal ( xmin() , other.xmin() ) &&
                  s_equal ( xmax() , other.xmax() )                  ,
                  "Cannot subtract Chebyshev with different domains" , 
                  "Ostap::Math::Chebyshev"                           , 
                  Ostap::StatusCode ( 523 )                          )  ;
  //
  const unsigned short idegree = std::max ( degree() , other.degree() ) ;
  //
  for ( unsigned short i = 0 ; i < other.npars() ; ++i ) 
  { m_pars[i] -= other.m_pars [ i ] ; }
  //
  return *this ;
}
// ============================================================================
// Subtract      polynomials (with the same domain!)
// ============================================================================
Ostap::Math::ChebyshevSum 
Ostap::Math::ChebyshevSum::subtract
( const Ostap::Math::ChebyshevSum& other ) const 
{
  ChebyshevSum result(*this) ;
  result.isub ( other ) ;
  return result ;
}


// ============================================================================
// get unique tag 
// ============================================================================
std::size_t Ostap::Math::ChebyshevSum::tag () const 
{
  static const std::string s_name { "ChebyshevSum" } ;
  return Ostap::Utils::hash_combiner 
    ( s_name , Ostap::Utils::hash_range ( m_pars ) ,  
      degree () , m_xmin , m_xmax ) ;
}
// ============================================================================





// ============================================================================
// LegendreSum 
// ============================================================================
// constructor from the degree 
// ============================================================================
Ostap::Math::LegendreSum::LegendreSum
( const unsigned short degree , 
  const double         xmin   , 
  const double         xmax   )  
  : Ostap::Math::PolySum ( degree ) 
  , m_xmin ( std::min ( xmin, xmax ) )
  , m_xmax ( std::max ( xmin, xmax ) )
{}
// ============================================================================
// constructor from parameters
// ============================================================================
Ostap::Math::LegendreSum::LegendreSum
( const std::vector<double>& pars , 
  const double               xmin , 
  const double               xmax )  
  : Ostap::Math::PolySum ( pars ) 
  , m_xmin ( std::min ( xmin, xmax ) )
  , m_xmax ( std::max ( xmin, xmax ) )
{}
// ============================================================================
// get the value
// ============================================================================
double Ostap::Math::LegendreSum::evaluate ( const double x ) const 
{
  if      ( 1 == m_pars.size()       ) { return m_pars[0] ; }
  else if ( zero ()                  ) { return         0 ; }
  // transform argument:
  const double tx = t ( x ) ;  
  // use Clenshaw's algorithm 
  return Ostap::Math::Clenshaw::legendre_sum ( m_pars.begin() , m_pars.end()  , tx ) ;
}
// ============================================================================
// get the integral between xmin and xmax
// ============================================================================
double Ostap::Math::LegendreSum::integral   () const 
{ return m_pars [ 0 ] * ( m_xmax - m_xmin ) ; }
// ============================================================================
// get the integral between low and high 
// ============================================================================
double Ostap::Math::LegendreSum::integral   
( const double low  , 
  const double high ) const 
{
  //
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high < low             ) { return -integral ( high , low ) ; }
  else if ( high <= m_xmin         ) { return 0 ; }
  else if ( low  >= m_xmax         ) { return 0 ; }
  else if ( zero ()                ) { return 0 ; }
  //
  if ( s_equal( low  , m_xmin ) && s_equal ( high , m_xmax ) ) { return integral () ; }
  //
  if      ( low  <  m_xmin         ) { return integral ( m_xmin , high   ) ; }
  else if ( high >  m_xmax         ) { return integral ( low    , m_xmax ) ; }
  //
  const double xl = t ( low  ) ;
  const double xh = t ( high ) ;
  //
  // // Trvial sum to be replaced with Clenshaw 
  // double result = 0 ;
  // for ( std::vector<double>::const_iterator ip = m_pars.begin() ; m_pars.end() != ip ; ++ip ) 
  // { 
  //   const double p = *ip ;
  //   if ( s_zero ( p ) ) { continue ; } // SKIP IT! 
  //   result += p * _legendre_int_ ( ip - m_pars.begin() , xl , xh  ) ;
  // }
  //
  std::vector<double> npars ( m_pars.size() + 1 , 0 ) ;
  for ( unsigned int i = 1 ; i < m_pars.size() ; ++i ) 
  { 
    if ( s_zero ( m_pars[i] ) ) { continue ; }           // CONTINUE 
    //
    npars [ i + 1 ] += m_pars[i] / ( 2*i + 1 ) ; 
    npars [ i - 1 ] -= m_pars[i] / ( 2*i + 1 ) ; 
  }
  //
  npars [ 1 ] += m_pars [ 0 ] ;
  //
  const double result = 
    Ostap::Math::Clenshaw::legendre_sum ( npars.begin() , npars.end() , xh ) -
    Ostap::Math::Clenshaw::legendre_sum ( npars.begin() , npars.end() , xl ) ;
  //
  const double dx = 0.5 *  ( m_xmax - m_xmin ) ;
  return result * dx ;
}
// ============================================================================
// get indefinite integral 
// ============================================================================
Ostap::Math::LegendreSum
Ostap::Math::LegendreSum::indefinite_integral ( const double C ) const
{
  //
  const double dx = 0.5 * ( m_xmax - m_xmin ) ;
  //
  Ostap::Math::LegendreSum integ (  degree() + 1 , m_xmin , m_xmax ) ;
  for ( unsigned int i = 1 ; i < m_pars.size() ; ++i ) 
  { 
    if ( s_zero ( m_pars [ i ] ) ) { continue ; }           // CONTINUE 
    //
    integ.m_pars [ i + 1 ] += m_pars [ i ] / ( 2*i + 1 ) * dx ; 
    integ.m_pars [ i - 1 ] -= m_pars [ i ] / ( 2*i + 1 ) * dx ; 
  }
  //
  integ.m_pars [ 1 ] += m_pars [ 0 ] * dx ;
  integ.m_pars [ 0 ] += C ; 
  //
  return integ ;
}
// ============================================================================
// get the derivative at point "x" 
// ============================================================================
double Ostap::Math::LegendreSum::derivative ( const double x     ) const 
{
  if ( x < m_xmin || x > m_xmax || 1 == m_pars.size() ) { return 0 ; }
  //
  const double tx = t ( x  ) ;
  //
  std::vector<double> npars ( m_pars.size() - 1 , 0 ) ;
  for ( unsigned short i = 1 ; i < m_pars.size() ; ++i )
  { 
    if ( s_zero ( m_pars[i] ) ) { continue ; }           // CONTINUE 
    //
    for ( int j = i-1 ; 0<=j  ; j-=2 ) 
    { npars[j] += m_pars[i] * ( 2 * j + 1 ) ; }
  }
  //
  const double dx = 2 / ( m_xmax - m_xmin ) ;
  return Ostap::Math::Clenshaw::legendre_sum ( npars.begin() , npars.end() , tx ) * dx ;  
}
// ============================================================================
// get the derivative 
// ============================================================================
Ostap::Math::LegendreSum 
Ostap::Math::LegendreSum::derivative () const 
{
  if ( 1 == m_pars.size() ) 
  { return Ostap::Math::LegendreSum ( 0 , m_xmin , m_xmax ) ; }
  //
  const double dx = 2 / ( m_xmax - m_xmin ) ;
  //
  Ostap::Math::LegendreSum deriv (  degree() - 1 , m_xmin , m_xmax ) ;
  for ( unsigned short i = 1 ; i < m_pars.size() ; ++i )
  { 
    if ( s_zero ( m_pars[i] ) ) { continue ; }           // CONTINUE 
    //
    for ( int j = i-1 ; 0<=j  ; j-=2 ) 
    { deriv.m_pars [ j ] += m_pars[i] * ( 2 * j + 1 ) * dx ; }
  }
  //
  return deriv ;
}
// ============================================================================
// simple  manipulations with polynoms: shift it! 
// ============================================================================
Ostap::Math::LegendreSum& 
Ostap::Math::LegendreSum::operator+= ( const double a ) 
{ m_pars[0] += a ; return *this ; }
// ============================================================================
// simple  manipulations with polynoms: shift it! 
// ============================================================================
Ostap::Math::LegendreSum& 
Ostap::Math::LegendreSum::operator-= ( const double a ) 
{ m_pars[0] -= a ; return *this ; }
// ============================================================================
// simple  manipulations with polynoms: scale it! 
// ============================================================================
Ostap::Math::LegendreSum& 
Ostap::Math::LegendreSum::operator*= ( const double a ) 
{ 
  Ostap::Math::scale ( m_pars , a ) ;
  return *this ; 
}
// ============================================================================
// simple  manipulations with polynoms: scale it 
// ============================================================================
Ostap::Math::LegendreSum& 
Ostap::Math::LegendreSum::operator/= ( const double a ) 
{ 
  Ostap::Math::scale ( m_pars , 1/a ) ;
  return *this ; 
}
// ============================================================================
// unary minus 
// ============================================================================
Ostap::Math::LegendreSum
Ostap::Math::LegendreSum::operator-() const 
{
  LegendreSum a ( *this ) ;
  Ostap::Math::negate ( a.m_pars ) ;
  return a ;
}
// ============================================================================
// Python
// ============================================================================
Ostap::Math::LegendreSum& 
Ostap::Math::LegendreSum::__iadd__   ( const double a )  
{ (*this) += a ; return *this ; } 
// ============================================================================
Ostap::Math::LegendreSum& 
Ostap::Math::LegendreSum::__isub__   ( const double a ) 
{ (*this) -= a ; return *this ; } 
// ============================================================================
Ostap::Math::LegendreSum& 
Ostap::Math::LegendreSum::__imul__   ( const double a )
{ (*this) *= a ; return *this ; } 
// ============================================================================
Ostap::Math::LegendreSum& 
Ostap::Math::LegendreSum::__itruediv__   ( const double a ) 
{ (*this) /= a ; return *this ; } 
// ============================================================================
Ostap::Math::LegendreSum 
Ostap::Math::LegendreSum::__add__    ( const double a ) const 
{ return LegendreSum(*this) += a ; } 
// ============================================================================
Ostap::Math::LegendreSum 
Ostap::Math::LegendreSum::__sub__    ( const double a ) const 
{ return LegendreSum(*this) -= a ; } 
// ============================================================================
Ostap::Math::LegendreSum 
Ostap::Math::LegendreSum::__mul__    ( const double a ) const 
{ return LegendreSum(*this) *= a ; } 
// ============================================================================
Ostap::Math::LegendreSum 
Ostap::Math::LegendreSum::__truediv__    ( const double a ) const 
{ return LegendreSum(*this) /= a ; } 
// ============================================================================
Ostap::Math::LegendreSum 
Ostap::Math::LegendreSum::__radd__   ( const double a ) const 
{ return __add__ ( a ) ; } 
// ============================================================================
Ostap::Math::LegendreSum 
Ostap::Math::LegendreSum::__rsub__   ( const double a ) const 
{ return (-(*this))+a; } 
// ============================================================================
Ostap::Math::LegendreSum 
Ostap::Math::LegendreSum::__rmul__   ( const double a ) const 
{ return __mul__ ( a ) ; } 
// ============================================================================
Ostap::Math::LegendreSum 
Ostap::Math::LegendreSum::__neg__   () const { return -(*this); } 
// ============================================================================
Ostap::Math::LegendreSum 
Ostap::Math::LegendreSum::__add__   
( const Ostap::Math::LegendreSum& a ) const { return sum      ( a ) ; } 
// ============================================================================
Ostap::Math::LegendreSum 
Ostap::Math::LegendreSum::__sub__   
( const Ostap::Math::LegendreSum& a ) const { return subtract ( a ) ; } 
// ============================================================================


// ============================================================================
// Add       polynomials (with the same domain!)
// ============================================================================
Ostap::Math::LegendreSum&
Ostap::Math::LegendreSum::isum
( const Ostap::Math::LegendreSum& other )
{
  // self-addition 
  if ( this == &other ) { (*this)*= 2 ; return *this ; }
  //
  Ostap::Assert ( s_equal ( xmin() , other.xmin() ) &&
                  s_equal ( xmax() , other.xmax() )            ,
                  "Cannot sum Legendre with different domains" , 
                  "Ostap::Math::LegendreSum"                   , 
                  Ostap::StatusCode ( 524 )                    )  ;
  //
  const unsigned short idegree = std::max ( degree() , other.degree() ) ;
  //
  m_pars.resize ( idegree + 1 ) ;
  //
  for ( unsigned short i = 0 ; i < other.npars() ; ++i ) 
  { m_pars[i] += other.m_pars[ i ] ; }
  //
  return *this ;
}
// ============================================================================
// Add       polynomials (with the same domain!)
// ============================================================================
Ostap::Math::LegendreSum 
Ostap::Math::LegendreSum::sum
( const Ostap::Math::LegendreSum& other ) const 
{
  LegendreSum result(*this) ;
  result.isum ( other ) ;
  return result ;
}
// ============================================================================
// Subtract      polynomials (with the same domain!)
// ============================================================================
Ostap::Math::LegendreSum&
Ostap::Math::LegendreSum::isub
( const Ostap::Math::LegendreSum& other )
{
  // self-subtraction
  if ( this == &other ) { *this *= 0.0 ; return *this ; }
  //
  Ostap::Assert ( s_equal ( xmin() , other.xmin() ) &&
                  s_equal ( xmax() , other.xmax() )                 ,
                  "Cannot subtract Legendre with different domains" , 
                  "Ostap::Math::LegendreSum"                        , 
                  Ostap::StatusCode ( 525 )                         )  ;
  //
  const unsigned short idegree = std::max ( degree() , other.degree() ) ;
  //
  m_pars.resize ( idegree + 1 ) ;
  //
  for ( unsigned short i = 0 ; i < other.npars() ; ++i ) 
  { m_pars[i] -= other.m_pars[ i ] ; }
  //
  return *this ;
}
// ============================================================================
// Subtract      polynomials (with the same domain!)
// ============================================================================
Ostap::Math::LegendreSum 
Ostap::Math::LegendreSum::subtract
( const Ostap::Math::LegendreSum& other ) const 
{
  LegendreSum result(*this) ;
  result.isub ( other ) ;
  return result ;
}
// ============================================================================
// get unique tag 
// ============================================================================
std::size_t Ostap::Math::LegendreSum::tag () const 
{
  static const std::string s_name { "LegendreSum" } ;
  return Ostap::Utils::hash_combiner 
    ( s_name , Ostap::Utils::hash_range ( m_pars ) ,  
      degree () , m_xmin , m_xmax ) ;
}
// ============================================================================




// ============================================================================
// constructor from the degree 
// ============================================================================
Ostap::Math::HermiteSum::HermiteSum 
( const unsigned short degree  ,
  const double         xmin    , 
  const double         xmax    )
  : Ostap::Math::PolySum ( degree )
  , m_xmin  ( std::min ( xmin , xmax ) )
  , m_xmax  ( std::max ( xmin , xmax ) )
  , m_scale ( 1 ) 
{
  m_scale /= ( m_xmax - m_xmin );
}
// ============================================================================
// constructor from the parameter list 
// ============================================================================
Ostap::Math::HermiteSum::HermiteSum 
( const std::vector<double>&  pars , 
  const double                xmin , 
  const double                xmax )
  : HermiteSum ( pars.begin() , 
                 pars.end  () , 
                 xmin         , 
                 xmax         )
{}
// ============================================================================
// get the value
// ============================================================================
double Ostap::Math::HermiteSum::operator() ( const double x ) const 
{ return evaluate ( x ) ; }
double Ostap::Math::HermiteSum:: evaluate  ( const double x ) const 
{
  const double tx = t ( x ) ;
  return Ostap::Math::Clenshaw::hermite_sum ( m_pars.begin() , m_pars.end() , tx ) ;
}
// ============================================================================
// get the derivative at point "x" 
// ============================================================================
double  Ostap::Math::HermiteSum::derivative ( const double x     ) const 
{
  const unsigned int d = degree() ;
  if ( 0 == d ) { return 0 ; }
  ///
  std::vector<double> deriv ( d , 0.0 ) ;
  for ( unsigned int k = 0 ; k < d ; ++k ) 
  { deriv[k] = ( k + 1 ) * m_pars[k+1] * 2 * m_scale ; }
  ///
  const double tx = t ( x ) ;
  return Ostap::Math::Clenshaw::hermite_sum ( deriv.begin() , deriv.end() , tx ) ;  
}
// ============================================================================
// get the derivative 
// ============================================================================
Ostap::Math::HermiteSum
Ostap::Math::HermiteSum::derivative () const 
{
  const unsigned int d = degree() ;
  if ( 0 == d ) { return HermiteSum ( 0 , m_xmin , m_xmax ) ; }
  ///
  HermiteSum deriv( d - 1 , m_xmin , m_xmax ) ;
  for ( unsigned int k = 0 ; k < d ; ++k ) 
  { deriv.m_pars[k] = ( k + 1 ) * m_pars[k+1] * 2 * m_scale ; }
  ///
  return deriv ;
}
// ============================================================================
// get the integral between low and high 
// ============================================================================
double Ostap::Math::HermiteSum::integral   
( const double low , const double high ) const 
{
  std::vector<double> integr ( m_pars.size() + 1 , 0.0 ) ;
  for ( unsigned int k = 0 ; k < m_pars.size() ; ++k ) 
  { integr[ k + 1 ] = m_pars[k] / ( k + 1 ) * 0.5 / m_scale ; }
  //
  const double th = t ( high ) ;
  const double tl = t ( low  ) ;
  //
  return 
    Ostap::Math::Clenshaw::hermite_sum ( integr.begin() , integr.end() , th ) -   
    Ostap::Math::Clenshaw::hermite_sum ( integr.begin() , integr.end() , tl ) ;
}
// ============================================================================
// get the integral
// ============================================================================
Ostap::Math::HermiteSum
Ostap::Math::HermiteSum::indefinite_integral ( const double c0 ) const 
{
  const unsigned int d = degree() ;
  ///
  HermiteSum integr ( d + 1 , m_xmin , m_xmax )  ;
  integr.m_pars[0] = c0 ;
  for ( unsigned int k = 0 ; k < m_pars.size() ; ++k ) 
  { integr.m_pars [ k + 1 ] = m_pars[k] / ( k + 1 ) * 0.5 / m_scale ; }
  ///
  return integr ;
}
// ============================================================================
// simple  manipulations with polynoms: shift it! 
// ============================================================================
Ostap::Math::HermiteSum& 
Ostap::Math::HermiteSum::operator+= ( const double a ) 
{ m_pars[0] += a ; return *this ; }
// ============================================================================
// simple  manipulations with polynoms: shift it! 
// ============================================================================
Ostap::Math::HermiteSum& 
Ostap::Math::HermiteSum::operator-= ( const double a ) 
{ m_pars[0] -= a ; return *this ; }
// ============================================================================
// simple  manipulations with polynoms: scal eit 
// ============================================================================
Ostap::Math::HermiteSum& 
Ostap::Math::HermiteSum::operator*= ( const double a ) 
{ 
  Ostap::Math::scale ( m_pars , a ) ;
  return *this ; 
}
// ============================================================================
// simple  manipulations with polynoms: scale it
// ============================================================================
Ostap::Math::HermiteSum& 
Ostap::Math::HermiteSum::operator/= ( const double a ) 
{ 
  Ostap::Math::scale ( m_pars , 1/a ) ;
  return *this ; 
}
// ============================================================================
// unary minus 
// ============================================================================
Ostap::Math::HermiteSum
Ostap::Math::HermiteSum::operator-() const 
{
  HermiteSum a ( *this ) ;
  Ostap::Math::negate ( a.m_pars ) ;
  return a ;
}
// ============================================================================



// ============================================================================
// Python
// ============================================================================
Ostap::Math::HermiteSum& 
Ostap::Math::HermiteSum::__iadd__   ( const double a )  
{ (*this) += a ; return *this ; } 
// ============================================================================
Ostap::Math::HermiteSum& 
Ostap::Math::HermiteSum::__isub__   ( const double a ) 
{ (*this) -= a ; return *this ; } 
// ============================================================================
Ostap::Math::HermiteSum& 
Ostap::Math::HermiteSum::__imul__   ( const double a )
{ (*this) *= a ; return *this ; } 
// ============================================================================
Ostap::Math::HermiteSum& 
Ostap::Math::HermiteSum::__itruediv__   ( const double a ) 
{ (*this) /= a ; return *this ; } 
// ============================================================================
Ostap::Math::HermiteSum 
Ostap::Math::HermiteSum::__add__    ( const double a ) const 
{ return HermiteSum(*this) += a ; } 
// ============================================================================
Ostap::Math::HermiteSum 
Ostap::Math::HermiteSum::__sub__    ( const double a ) const 
{ return HermiteSum(*this) -= a ; } 
// ============================================================================
Ostap::Math::HermiteSum 
Ostap::Math::HermiteSum::__mul__    ( const double a ) const 
{ return HermiteSum(*this) *= a ; } 
// ============================================================================
Ostap::Math::HermiteSum 
Ostap::Math::HermiteSum::__truediv__    ( const double a ) const 
{ return HermiteSum(*this) /= a ; } 
// ============================================================================
Ostap::Math::HermiteSum 
Ostap::Math::HermiteSum::__radd__   ( const double a ) const 
{ return __add__ ( a ) ; } 
// ============================================================================
Ostap::Math::HermiteSum 
Ostap::Math::HermiteSum::__rsub__   ( const double a ) const 
{ return (-(*this))+a; } 
// ============================================================================
Ostap::Math::HermiteSum 
Ostap::Math::HermiteSum::__rmul__   ( const double a ) const 
{ return __mul__ ( a ) ; } 
// ============================================================================
Ostap::Math::HermiteSum 
Ostap::Math::HermiteSum::__neg__   () const { return -(*this); } 
// ============================================================================
Ostap::Math::HermiteSum 
Ostap::Math::HermiteSum::__add__   
( const Ostap::Math::HermiteSum& a ) const { return sum      ( a ) ; } 
// ============================================================================
Ostap::Math::HermiteSum 
Ostap::Math::HermiteSum::__sub__   
( const Ostap::Math::HermiteSum& a ) const { return subtract ( a ) ; } 


// ============================================================================
// Add       polynomials (with the same domain!)
// ============================================================================
Ostap::Math::HermiteSum&
Ostap::Math::HermiteSum::isum
( const Ostap::Math::HermiteSum& other )
{
  // self-addition 
  if ( this == &other ) { (*this)*= 2 ; return *this ; }
  //
  Ostap::Assert ( s_equal ( xmin() , other.xmin() ) &&
                  s_equal ( xmax() , other.xmax() )           ,
                  "Cannot sum Hermite with different domains" , 
                  "Ostap::Math::HermiteSum"                   , 
                  Ostap::StatusCode ( 528 )                   )  ;
  //
  const unsigned short idegree = std::max ( degree() , other.degree() ) ;
  //
  m_pars.resize ( idegree + 1 ) ;
  //
  for ( unsigned short i = 0 ; i < other.npars() ; ++i ) 
  { m_pars[i] += other.m_pars[ i ] ; }
  //
  return *this ;
}
// ============================================================================
// Add       polynomials (with the same domain!)
// ============================================================================
Ostap::Math::HermiteSum 
Ostap::Math::HermiteSum::sum
( const Ostap::Math::HermiteSum& other ) const 
{
  HermiteSum result(*this) ;
  result.isum ( other ) ;
  return result ;
}
// ============================================================================
// Subtract      polynomials (with the same domain!)
// ============================================================================
Ostap::Math::HermiteSum&
Ostap::Math::HermiteSum::isub
( const Ostap::Math::HermiteSum& other )
{
  // self-subtraction
  if ( this == &other ) { *this *= 0.0 ; return *this ; }
  //
  Ostap::Assert ( s_equal ( xmin() , other.xmin() ) &&
                  s_equal ( xmax() , other.xmax() )                ,
                  "Cannot subtract Hermite with different domains" , 
                  "Ostap::Math::HermiteSum"                        , 
                  Ostap::StatusCode ( 529 )                        )  ;
  //
  const unsigned short idegree = std::max ( degree() , other.degree() ) ;
  //
  m_pars.resize ( idegree + 1 ) ;
  //
  for ( unsigned short i = 0 ; i < other.npars() ; ++i ) 
  { m_pars[i] -= other.m_pars[ i ] ; }
  //
  return *this ;
}
// ============================================================================
// Subtract      polynomials (with the same domain!)
// ============================================================================
Ostap::Math::HermiteSum 
Ostap::Math::HermiteSum::subtract
( const Ostap::Math::HermiteSum& other ) const 
{
  HermiteSum result(*this) ;
  result.isub ( other ) ;
  return result ;
}
// ============================================================================
// get unique tag 
// ============================================================================
std::size_t Ostap::Math::HermiteSum::tag () const 
{
  static const std::string s_name { "HermiteSum" } ;
  return Ostap::Utils::hash_combiner 
    ( s_name , Ostap::Utils::hash_range ( m_pars ) ,  
      degree () , m_xmin , m_xmax ) ;
}
// ============================================================================


// ============================================================================
/// legendre to bernstein transformation 
namespace 
{
  // ==========================================================================
  /** transformation matrix from bernstein to legendre basis 
   *  @see http://www.sciencedirect.com/science/article/pii/S0377042700003769 eq.21
   */
  inline 
  double b2l_mtrx 
  ( const unsigned short j , 
    const unsigned short k ,
    const unsigned short n ) 
  {
    //
    long long r = 0 ;
    for ( unsigned short i = 0 ; i <= j ; ++i ) 
    {
      0 == ( j + i ) % 2 ?
        r +=
        Ostap::Math::choose ( j              , i     ) * 
        Ostap::Math::choose ( k + i          , k     ) * 
        Ostap::Math::choose ( n - k + j - i  , n - k ) : 
        r -=
        Ostap::Math::choose ( j              , i     ) * 
        Ostap::Math::choose ( k + i          , k     ) * 
        Ostap::Math::choose ( n - k + j - i  , n - k ) ;
    }
    //
    return r * ( 2* j + 1 ) / double ( Ostap::Math::choose ( n + j , n ) ) / ( n + j + 1 ) ;
  }
  // ==========================================================================
  /** transformation matrix from bernstein basic to monomial basis 
   *  @see http://www.sciencedirect.com/science/article/pii/S0377042700003769 eq.21
   */
  inline 
  long long b2m_mtrx 
  ( const unsigned short j , 
    const unsigned short k ,
    const unsigned short n ) 
  {
    const long long r = 
      Ostap::Math::choose ( n , j ) * 
      Ostap::Math::choose ( j , k ) ;
    return 0 == ( j + k ) % 2 ? r : -r ;
  }
  // ==========================================================================
  /// affine transformation of polynomial
  inline 
  double m2m_mtrx_1
  ( const unsigned short j , 
    const unsigned short k ) 
  {
    if ( k < j ) { return 0 ; }
    const double c = Ostap::Math::choose ( k , j ) ;
    return c / Ostap::Math::POW ( 2 , k ) ;
  }
  // ==========================================================================
  /// transformation matrix from monomial to chebyshev basis 
  inline 
  double m2c_mtrx 
  ( const unsigned short j , 
    const unsigned short k )
  {
    return 
      ( 1 == ( j + k ) % 2 ) ? 0. :
      Ostap::Math::choose ( k , ( k - j ) / 2 ) * ( j == 0 ? 1. : 2. ) 
      / Ostap::Math::POW ( 2 , k ) ;
  }
  // ============================================================================
  /// transformation matrix from legendre basic to monomial basis 
  inline 
  double l2m_mtrx 
  ( const unsigned short j , 
    const unsigned short k )
  {
    return 
      1 == ( j + k ) % 2 ? 0.0 : 
      Ostap::Math::choose       ( k         , j ) * 
      Ostap::Math::choose_half  ( j + k - 1 , k ) *
      Ostap::Math::POW          ( 2         , k ) ;
  }
  // ==========================================================================
}
// ============================================================================
// constructor from bernstein polynomial
// ============================================================================
Ostap::Math::LegendreSum::LegendreSum
( const Ostap::Math::Bernstein& poly )  
  : Ostap::Math::PolySum ( poly.degree () ) 
  , m_xmin ( poly.xmin() ) 
  , m_xmax ( poly.xmax() ) 
{
  const unsigned short np = npars  () ;
  const unsigned short d  = degree () ;
  for ( unsigned short i  = 0 ; i < np ; ++i ) 
  { 
    for ( unsigned short k = 0 ; k < np ; ++k ) 
    { 
      const double p = poly.par ( k ) ;
      if ( s_zero ( p ) ) { continue ; }
      m_pars[i] += b2l_mtrx ( i , k , d ) * p ; 
    } 
  }
}
// ============================================================================
// constructor from bernstein polynomial
// ============================================================================
Ostap::Math::Polynomial::Polynomial
( const Ostap::Math::Bernstein& poly )  
  : Ostap::Math::PolySum ( poly.degree () ) 
  , m_xmin ( poly.xmin() ) 
  , m_xmax ( poly.xmax() ) 
{
  //
  const unsigned short np = npars  () ;
  const unsigned short d  = degree () ;
  //
  // 2-step tranfromation
  //
  // 1: tranfrom to 'regular' ponynomial basic 
  std::vector<double> _pars( m_pars.size() ) ;
  for ( unsigned short i = 0 ; i < np ; ++i ) 
  { 
    for ( unsigned short k = 0 ; k <= i ; ++k )  // ATTENTION!!! 
    { 
      const double p = poly.par ( k ) ;
      if ( s_zero ( p ) ) { continue ; }
      _pars[i] += b2m_mtrx ( i , k , d ) * p ; 
    } 
  }
  //
  // 2: affine tranform to use [-1,-1] range 
  for ( unsigned short i = 0 ; i < np ; ++i ) 
  { 
    for ( unsigned short k = i ; k < np ; ++k )  // ATTENTION!!!
    { 
      const double p = _pars[  k ] ;
      if ( s_zero ( p ) ) { continue ; }
      m_pars[i] += m2m_mtrx_1 ( i , k ) * p ; 
    } 
  }  
}
// ============================================================================
// constructor from legendre polynomial
// ============================================================================
Ostap::Math::Polynomial::Polynomial
( const Ostap::Math::LegendreSum& poly )  
  : Ostap::Math::PolySum ( poly.degree () ) 
  , m_xmin ( poly.xmin() ) 
  , m_xmax ( poly.xmax() ) 
{
  //
  const unsigned short np = npars  () ;
  //
  for ( unsigned short i = 0 ; i < np ; ++i ) 
  { 
    for ( unsigned short k = i ; k < np ; k+=2 )   // ATTENTION!!! 2!
    { 
      const double p = poly.par ( k ) ;
      if ( s_zero ( p ) ) { continue ; }
      m_pars[i] += l2m_mtrx ( i , k ) * p ; 
    } 
  }
}
// =============================================================================
//  constructor from Polinomial 
// =============================================================================
Ostap::Math::ChebyshevSum::ChebyshevSum
( const Ostap::Math::Polynomial& poly )  
  : Ostap::Math::PolySum ( poly.degree () ) 
  , m_xmin ( poly.xmin() ) 
  , m_xmax ( poly.xmax() ) 
{
  const unsigned short np = npars  () ;
  for ( unsigned short i = 0 ; i < np ; ++i ) 
  { 
    for ( unsigned short k = i ; k < np  ; k+=2 ) // ATTENTION !!! 2!
    { 
      const double p = poly.par ( k ) ;
      if ( s_zero ( p ) ) { continue ; }
      m_pars[i] += m2c_mtrx ( i , k ) * p ; 
    } 
  }  
}
// ============================================================================


// ============================================================================
// Delegation
// ============================================================================

// ============================================================================
// constructor from chebyshev polynomial (delegation)
// ============================================================================
Ostap::Math::Polynomial::Polynomial
( const Ostap::Math::ChebyshevSum& poly )  
  : Ostap::Math::Polynomial ( Ostap::Math::Bernstein ( poly ) ) 
{}
// ============================================================================
// constructor from bernstein polynomial (deleghation)
// ============================================================================
Ostap::Math::ChebyshevSum::ChebyshevSum
( const Ostap::Math::Bernstein& poly )  
  : Ostap::Math::ChebyshevSum ( Ostap::Math::Polynomial ( poly ) ) 
{}
// ============================================================================
// constructor from legendre polynomial (delegation)
// ============================================================================
Ostap::Math::ChebyshevSum::ChebyshevSum
( const Ostap::Math::LegendreSum& poly )  
  : Ostap::Math::ChebyshevSum ( Ostap::Math::Polynomial ( poly ) ) 
{}
// ============================================================================
// constructor from polynomial (delegation)
// ============================================================================
Ostap::Math::LegendreSum::LegendreSum
( const Ostap::Math::Polynomial& poly )  
  : Ostap::Math::LegendreSum ( Ostap::Math::Bernstein ( poly ) ) 
{}
// ============================================================================
// constructor from chebyshev polynomial (delegation)
// ============================================================================
Ostap::Math::LegendreSum::LegendreSum
( const Ostap::Math::ChebyshevSum& poly )  
  : Ostap::Math::LegendreSum ( Ostap::Math::Bernstein ( poly ) ) 
{}


// ============================================================================
// integation with an exponent 
// ============================================================================
namespace 
{
  // ===========================================================================
  // integation with an exponent 
  // ==========================================================================
  template <class POLYNOMIAL>
  inline double _integrate_
  ( const POLYNOMIAL& poly , 
    const long double tau  , 
    const long double low  , 
    const long double high ) 
  {
    const long double xlow  = std::max ( low  , (long double) poly.xmin() ) ;
    const long double xhigh = std::min ( high , (long double) poly.xmax() ) ;
    //
    // a bit esoteric way to get numerically correct resuls...
    //
    const long double eH = std::expm1 ( tau * xhigh ) ;
    const long double eL = std::expm1 ( tau * xlow  ) ;
    const long double pH = poly       (       xhigh ) ;
    const long double pL = poly       (       xlow  ) ;
    //
    const long double p1 = ( eH * pH - eL * pL ) + ( pH - pL ) ;
    //
    if ( 1 >= poly.npars  () ) { return p1 / tau ; } // RETURN 
    //
    return ( p1 - _integrate_ ( poly.derivative() , tau , xlow , xhigh ) )  / tau ;
  }
  // ==========================================================================
}
// ============================================================================
/* get the integral between low and high for a product of
 *  polynom and the exponential function with the exponent tau
 *  \f[  \int_{a}^{b} \mathcal{P} e^{\tau x } \mathrm{d}x \f] 
 *  @param poly  polynomial
 *  @param tau   slope parameter for exponential 
 *  @param low   low  integration range 
 *  @param high  high integration range 
 */
// ============================================================================
double Ostap::Math::integrate 
( const Ostap::Math::Polynomial& poly ,
  const double                   tau  ,
  const double                   low  , 
  const double                   high ) 
{
  //
  // the exponent is totally redundant...
  if      ( s_zero  ( tau )         ) { return  poly.integral ( low  , high         ) ; } 
  else if ( s_small ( tau )         ) { return  poly.integral ( low  , high         ) ; } 
  // zero range 
  else if ( s_equal ( low , high )  ) { return 0 ; }
  // ingertaion of zeros is simple 
  else if ( poly.zero ()            ) { return 0 ; }
  // invert range 
  else if ( low  >  high            ) { return -integrate ( poly , tau , high , low ) ; }
  // our the range 
  else if ( high <  poly.xmin () || 
            low  >  poly.xmax ()    ) { return  0 ; }
  // nice range to use analytic expression 
  else if ( s_equal ( low  , poly.xmin() ) && 
            s_equal ( high , poly.xmax() ) ) { return integrate ( poly , tau )  ; }
  //
  // check if special "small-tau" algorithm is needed here 
  // 
  const long double  xmin = poly.xmin () ;
  const long double  xmax = poly.xmax () ;
  const long double  _tau = ( xmax - xmin ) * tau / 2   ;
  const unsigned int N    = poly.degree() ;
  const long double  t1   = Ostap::Math::POW ( std::abs (  tau ) , N + 1 ) ;
  const long double  t2   = Ostap::Math::POW ( std::abs ( _tau ) , N + 1 ) ;
  if ( s_small ( t1 ) || s_small ( t2 ) )
  {
    const long double _fac = std::exp   ( ( xmax + xmin ) * tau / 2 ) ;
    const long double tmin = poly.t ( low ) ;
    const long double tmax = poly.t ( high ) ;    
    //
    long double result = 0 ;
    //
    const std::vector<double>& pars = poly.pars() ;
    for ( unsigned short i = 0 ; i < pars.size() ; ++i )
    { 
      const long double p = pars[i] ;
      if ( s_zero ( p ) ) { continue ; }                 // CONTINUE 
      const long double tl = Ostap::Math::POW ( tmin , i + 1 ) ;
      const long double th = Ostap::Math::POW ( tmax , i + 1 ) ;
      result += p * _factorial_d_ ( i ) * 
        ( th * Ostap::Math::gamma_star ( i + 1 , -_tau * tmax ) -
          tl * Ostap::Math::gamma_star ( i + 1 , -_tau * tmin ) ) ;
    }
    //
    result *= ( xmax - xmin ) * _fac / 2 ;
    //
    return result ;
  }
  // try generic recursive efficient, but numericaly highly unstable scheme  
  return _integrate_ ( poly , tau , low , high );
}
// ============================================================================
/*  get the integral between low and high for a product of
 *  Chebyshev polynom and the exponential function with the exponent tau
 *  \f[  \int_{a}^{b} \mathcal{T} e^{\tau x } \mathrm{d}x \f] 
 *  @param poly  chebyshev polynomial
 *  @param tau   slope parameter for exponential 
 *  @param low   low  integration range 
 *  @param high  high integration range 
 */
// ============================================================================
double Ostap::Math::integrate 
( const Ostap::Math::ChebyshevSum& poly ,
  const double                     tau  ,
  const double                     low  , 
  const double                     high ) 
{
  //
  if      ( s_zero  ( tau )         ) { return  poly.integral ( low  , high         ) ; } 
  else if ( s_small ( tau )         ) { return  poly.integral ( low  , high         ) ; } 
  else if ( s_equal ( low , high )  ) { return 0 ; }
  else if ( poly.zero ()            ) { return 0 ; }
  else if ( low  >  high            ) { return -integrate ( poly , tau , high , low ) ; }
  else if ( high <  poly.xmin () || 
            low  >  poly.xmax ()    ) { return  0 ; }
  //
  // check if special "small-tau" algorithm is needed here 
  // 
  const long double  xmin = poly.xmin () ;
  const long double  xmax = poly.xmax () ;
  const long double  _tau = ( xmax - xmin ) * tau / 2   ;
  const unsigned int N    = poly.degree() ;
  const long double  t1   = Ostap::Math::POW ( std::abs (  tau ) , N + 1 ) ;
  const long double  t2   = Ostap::Math::POW ( std::abs ( _tau ) , N + 1 ) ;
  if ( s_small ( t1 ) || s_small ( t2 ) )
  { 
    const Ostap::Math::Polynomial p ( poly ) ;
    return integrate ( p , tau , low , high ) ; 
  }
  //
  // try generic recursive, but numericaly highly unstable scheme  
  return _integrate_ ( poly , tau , low , high );
}
// ========================================================================    
/* get the integral between low and high for a product of
 *  Legendre polynom and the exponential function with the exponent tau
 *  \f[  \int_{a}^{b} \mathcal{L} e^{\tau x } \mathrm{d}x \f] 
 *  @param poly  Legendre polynomial
 *  @param tau   slope parameter for exponential 
 *  @param a     low  integration range 
 *  @param b     high integration range 
 */
// ============================================================================
double Ostap::Math::integrate 
( const Ostap::Math::LegendreSum& poly ,
  const double                    tau  ,
  const double                    low  , 
  const double                    high ) 
{
  //
  if      ( s_zero  ( tau )         ) { return  poly.integral ( low  , high         ) ; } 
  else if ( s_small ( tau )         ) { return  poly.integral ( low  , high         ) ; } 
  else if ( s_equal ( low , high )  ) { return 0 ; }
  else if ( poly.zero ()            ) { return 0 ; }
  else if ( low  >  high            ) { return -integrate ( poly , tau , high , low ) ; }
  else if ( high <  poly.xmin () || 
            low  >  poly.xmax ()    ) { return  0 ; }
  //
  //
  // check if special "small-tau" algorithm is needed here 
  // 
  const long double  xmin = poly.xmin () ;
  const long double  xmax = poly.xmax () ;
  const long double  _tau = ( xmax - xmin ) * tau / 2   ;
  const unsigned int N    = poly.degree() ;
  const long double  t1   = Ostap::Math::POW ( std::abs (  tau ) , N + 1 ) ;
  const long double  t2   = Ostap::Math::POW ( std::abs ( _tau ) , N + 1 ) ;
  if ( s_small ( t1 ) || s_small ( t2 ) )
  { 
    const Ostap::Math::Polynomial p ( poly ) ;
    return integrate ( p , tau , low , high ) ; 
  }
  //
  // try generic recursive, but numericaly highly unstable scheme  
  return _integrate_ ( poly , tau , low , high );
}
// ============================================================================


// ============================================================================
/*  get the integral between low and high for a product of
 *  polynom and the exponential function with the exponent tau
 *  \f[  \int \mathcal{P} e^{\tau x } \mathrm{d}x \f] 
 *  @param poly  polynomial
 *  @param tau   slope parameter for exponential 
 */
// ============================================================================
double Ostap::Math::integrate 
( const Ostap::Math::Polynomial& poly ,
  const double                   tau  ) 
{
  //
  if      ( s_zero  ( tau ) || s_small ( tau ) ) 
  { return  poly.integral ( poly.xmin() , poly.xmax() ) ; } 
  else if ( poly.zero()                        ) { return  0 ; }
  //
  const long double xmin = poly.xmin () ;
  const long double xmax = poly.xmax () ;
  //
  const long double _tau =              ( xmax - xmin ) * tau / 2   ;
  const long double _fac = std::exp   ( ( xmax + xmin ) * tau / 2 ) ;
  //
  if      ( s_zero  ( _tau ) || s_small ( _tau ) ) 
  { return  poly.integral ( poly.xmin() , poly.xmax() ) * _fac ; }
  //
  long double result = 0 ;
  const std::vector<double>& pars = poly.pars() ;
  for ( unsigned short i = 0 ; i < pars.size() ; ++i )
  { 
    const long double p = pars[i] ;
    if ( s_zero ( p ) ) { continue ; }                 // CONTINUE 
    //
    result += p * Ostap::Math::beta_N ( i , -_tau ) ;  // NOTE THE SIGN!!
  }
  //
  return result * ( xmax - xmin ) * _fac / 2 ;
}
// ============================================================================
/*  construct chebyshev approximation for arbitrary function 
 *  @param func the function
 *  @param N    degree of polynomial 
 *  @param x_min low edge
 *  @param x_max high edge 
 *  @return Chebyshev approximation 
 *  @see ChebyshevSum 
 *  @code 
 *  FUNC func = ...
 *  ChebyshevSum a = chebyshev_sum ( func , 10 ,  xmin , xmax ) ;
 *  @endcode 
 */
// ============================================================================
Ostap::Math::ChebyshevSum 
Ostap::Math::chebyshev_sum
( std::function<double(double)> func  , 
  const unsigned short          N     , 
  const double                  x_min , 
  const double                  x_max ) 
{ 
  // array of precomputed function values 
  std::vector<double> fv ( N ) ;
  // 
  const double xmin = std::min ( x_min , x_max ) ;
  const double xmax = std::max ( x_min , x_max ) ;
  //
  const double      xhs  = 0.5 * ( xmin + xmax ) ;
  const double      xhd  = 0.5 * ( xmax - xmin ) ;
  const long double pi_N = M_PIl / N ;
  auto _xi_ = [xhs,xhd,pi_N] ( const unsigned short k ) 
    { return std::cos ( pi_N * ( k + 0.5 ) ) * xhd + xhs ; } ;
  //
  for ( unsigned short i = 0 ; i < N ; ++i ) { fv[i] = func ( _xi_ ( i ) ) ; }
  //
  Ostap::Math::ChebyshevSum cs ( N , xmin , xmax ) ;
  for ( unsigned short i = 0 ; i < N + 1 ; ++i ) 
  {
    double c_i = 0 ;
    if ( 0 == i ) { for ( unsigned short k = 0 ; k < N ; ++k ) { c_i += fv[k] ; } }
    else 
    {
      for ( unsigned short k = 0 ; k < N ; ++k ) 
      { c_i += fv[k] * std::cos ( pi_N * i * ( k + 0.5 ) ) ; }
    }
    c_i *= 2.0 / N ;
    if ( 0 == i ) { c_i *= 0.5 ;}
    cs.setPar ( i, c_i ) ;
  }
  return cs ;
}
// ============================================================================



// ============================================================================
// The END  
// ============================================================================
