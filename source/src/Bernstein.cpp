// ============================================================================
// Include files 
// ============================================================================
// STD& STL
// ============================================================================
#include <cmath>
#include <array>
#include <climits>
#include <cassert>
#include <numeric>
// ============================================================================
// Ostap 
// ============================================================================
#include "Ostap/Math.h"
#include "Ostap/NSphere.h"
#include "Ostap/Power.h"
#include "Ostap/Choose.h"
#include "Ostap/Polynomials.h"
#include "Ostap/MoreMath.h"
#include "Ostap/Bernstein.h"
// ============================================================================
// Local
// ============================================================================
#include "Exception.h"
// ============================================================================
/** @file 
 *  Implementation file for functions, related to Bernstein's polynomnials 
 *
 *  @see http://en.wikipedia.org/wiki/Bernstein_polynomial
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2010-04-19
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  /// equality criteria for doubles
  const Ostap::Math::Equal_To<double> s_equal{} ; // equality criteria for doubles
  /// zero for doubles  
  const Ostap::Math::Zero<double>     s_zero{}  ; // zero for doubles
  /// zero fo vectors 
  const Ostap::Math::Zero< std::vector<double> > s_vzero{} ; // zero for vectors
  ///
  static_assert( std::numeric_limits<double>::is_specialized , 
                 "std::numeric_limits<double> is not specialized" ) ;
  static_assert( std::numeric_limits<long double>::is_specialized , 
                 "std::numeric_limits<long double> is not specialized" ) ;
  /// small value 
  const Ostap::Math::Small<long double> s_small
  ( 2.0L * std::numeric_limits<double>::epsilon() ) ;
  // ==========================================================================
  // De Casteljau's algorithm
  template <class ITERATOR>
  long double _casteljau_
  ( ITERATOR          first ,
    ITERATOR          last  ,
    const long double t0    ,
    const long double t1    )
  {
    // the trivial cases
    if      ( first == last    ) { return 0       ; }
    //
    const std::size_t len  = std::distance ( first , last  ) ;
    //
    if      ( 1 == len ) { return       *first                        ; }
    else if ( 2 == len ) { return t1 * (*first) + t0 * ( *(first+1) ) ; }
    //
    ITERATOR second = --last ;
    //
    // prepare recursion
    for ( ITERATOR it = first ; it != second ; ++it )
    { *it = t1 * ( *it )  + t0 * ( *( it+1 ) ) ; }
    //
    // recursion
    return _casteljau_ ( first , second , t0 , t1 ) ;
  }
  // ==========================================================================
}
// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::Bernstein::Bernstein
( const unsigned short      N    ,
  const double              xmin ,
  const double              xmax )
  : Ostap::Math::PolySum ( N ) 
  , m_xmin ( std::min ( xmin , xmax ) )
  , m_xmax ( std::max ( xmin , xmax ) )
{}
// ============================================================================
// constructor from the list of  coefficients
// ============================================================================
Ostap::Math::Bernstein::Bernstein
( const std::vector<double>& pars ,
  const double               xmin ,
  const double               xmax )
  : Ostap::Math::PolySum ( pars ) 
  , m_xmin ( std::min ( xmin , xmax ) )
  , m_xmax ( std::max ( xmin , xmax ) )
{}
// ============================================================================
// constructor from the list of  coefficients
// ============================================================================
Ostap::Math::Bernstein::Bernstein
( std::vector<double>&& ps   ,
  const double          xmin ,
  const double          xmax )
  : Ostap::Math::PolySum ( std::forward<std::vector<double> >( ps ) )
  , m_xmin ( std::min ( xmin , xmax ) )
  , m_xmax ( std::max ( xmin , xmax ) )
{}
// ============================================================================
// constructor  from Bernstein polynomial from *different* domai
// ============================================================================
namespace 
{
  //
  inline double _mjk_ ( const unsigned short    j    , 
                        const unsigned short    k    ,
                        const unsigned short    n    , 
                        Ostap::Math::Bernstein& ba   , 
                        Ostap::Math::Bernstein& bb   , 
                        const double            abar , 
                        const double            bbar ) 
  {
    if ( j > n || k > n ) { return 0 ; }
    //
    const unsigned short imin =  j + k <= n ? 0 : ( j + k - n )  ;
    const unsigned short imax =  std::min ( j , k ) ;
    //
    double m = 0 ;
    for ( unsigned short i = imin ; i <= imax ; ++i ) 
    {
      ba.setPar ( k - i , 1 ) ;      
      bb.setPar (     i , 1 ) ;      
      m += ba.evaluate ( abar ) * bb.evaluate ( bbar ) ;
      ba.setPar ( k - i , 0 ) ;      
      bb.setPar (     i , 0 ) ;      
    }
    return m ;
  }                         
}
// ============================================================================
Ostap::Math::Bernstein::Bernstein
( const Ostap::Math::Bernstein& poly , 
  const double                  xmin , 
  const double                  xmax ) 
  : Ostap::Math::PolySum ( poly       ) 
  , m_xmin ( std::min ( xmin , xmax ) )
  , m_xmax ( std::max ( xmin , xmax ) )    
{
  // recalculate domain ?
  if ( !s_equal ( this->xmin() , poly.xmin() ) ||
       !s_equal ( this->xmax() , poly.xmax() ) ) 
  {
    //
    std::vector<double> new_pars ( npars   () , 0 ) ;
    //
    const double a    = poly .xmin()  ;
    const double b    = poly .xmax()  ;
    const double abar = this->xmin()  ;
    const double bbar = this->xmax()  ;
    //
    const unsigned short N     = degree () ;
    for ( unsigned short j = 0 ; j <= N ; ++j ) 
    {
      //
      Ostap::Math::Bernstein ba ( N - j  , a , b ) ;
      Ostap::Math::Bernstein bb (     j  , a , b ) ;
      //
      for ( unsigned short k = 0 ; k <= N ; ++k ) 
      { new_pars[j] += _mjk_ ( j  , k  , N , 
                               ba , bb , abar , bbar ) * par ( k ) ; }   
    }
    //
    for ( unsigned short k = 0 ; k <= N ; ++k ) 
    { setPar ( k , new_pars[k] ) ; }
  }  
}
// ============================================================================
// construct the basic bernstein polinomial
// ============================================================================
Ostap::Math::Bernstein::Bernstein 
( const Ostap::Math::Bernstein::Basic& bb   , 
  const double                         xmin , 
  const double                         xmax ) 
  : Ostap::Math::PolySum ( bb.N()  ) 
  , m_xmin ( std::min ( xmin , xmax ) )
  , m_xmax ( std::max ( xmin , xmax ) )
    //
{
  if ( bb.k() <= bb.N() ) { m_pars[ bb.k() ] = 1 ; } 
}
// ============================================================================
/*  construct Bernstein interpolant
 *  @param x    vector of abscissas 
 *  @param y    vector of function values 
 *  @param xmin low  edge for Bernstein polynomial
 *  @param xmin high edge for Bernstein polynomial
 *  - if vector of y is longer  than vector x, extra values are ignored 
 *  - if vector of y is shorter than vector x, missing entries are assumed to be zero  
 *  It relies on Newton-Bernstein algorithm
 *  @see http://arxiv.org/abs/1510.09197
 *  @see Mark Ainsworth and Manuel A. Sanches, 
 *       "Computing of Bezier control points of Largangian interpolant 
 *       in arbitrary dimension", arXiv:1510.09197 [math.NA]
 *  @see http://adsabs.harvard.edu/abs/2015arXiv151009197A
 */
// ============================================================================
Ostap::Math::Bernstein::Bernstein 
( const std::vector<double>& x     , 
  const std::vector<double>& y     , 
  const double               xmin  ,
  const double               xmax  )
  : Ostap::Math::Bernstein ( x.begin() , x.end() , 
                             y.begin() , y.end() , 
                             xmin      , xmax    ) 
{} 
// ============================================================================
/* construct Bernstein polynomial from its roots
 *  Polinomial has a form
 *  \f$ B(x) = \prod_i (x-r_i) \prod_j (x-c_i)(x-c_i^*) \f$
 *  @param xmin low  edge for Bernstein polynomial
 *  @param xmax high edge for Bernstein polynomial
 *  @param r  the list of real  roots of the polinomial
 *  @param c  the list of complex roots (only one root from cc-pair is needed)
 */
// ============================================================================
Ostap::Math::Bernstein::Bernstein 
( const double xmin                            , 
  const double xmax                            , 
  const std::vector<double>& r                 ,
  const std::vector<std::complex<double> > & c )
  : Ostap::Math::PolySum ( r.size() + 2*c.size() ) 
  , m_xmin ( std::min ( xmin , xmax ) )
  , m_xmax ( std::max ( xmin , xmax ) )
{
  //
  Bernstein result ( std::vector<double> ( 1 , 1.0 ) , xmin , xmax ) ;
  Bernstein b1     ( std::vector<double> ( 2 , 1.0 ) , xmin , xmax ) ;
  Bernstein b2     ( std::vector<double> ( 3 , 1.0 ) , xmin , xmax ) ;
  for ( double rr : r )
  { 
    const double dmn = m_xmin - rr  ;
    const double dmx = m_xmax - rr  ;
    if      ( s_zero ( dmn ) ) 
    {
      b1.setPar ( 0 , 0 ) ;
      b1.setPar ( 1 , 1 ) ; 
    }
    else if ( s_zero ( dmx ) ) 
    {
      b1.setPar ( 0 , 1 ) ;
      b1.setPar ( 1 , 0 ) ; 
    }
    else 
    {
      b1.setPar ( 0 , dmn ) ;
      b1.setPar ( 1 , dmx ) ;
    }
    result = result * b1 ;
  }
  const double xmid = 0.5 * ( m_xmin + m_xmax );
  for ( std::complex<double> cr : c )
  {
    //  a * x * x + bx + c 
    const double a =  1               ;
    const double b = -2*cr.real()     ;
    const double c = std::norm ( cr ) ;    
    //
    const double a0 = c + m_xmin * ( b + m_xmin * a ) ;
    const double a1 = c +   xmid * ( b +   xmid * a ) ;
    const double a2 = c + m_xmax * ( b + m_xmax * a ) ;
    //    
    b2.setPar ( 0 ,     a0                      ) ;
    b2.setPar ( 1 ,  2 * a1 - 0.5 * ( a0 + a2 ) ) ;
    b2.setPar ( 2 ,     a2                      ) ;
    //
    result = result * b2 ;
  }
  //
  this->m_pars = std::move( result.m_pars ) ;
  //
  // scale it 
  const short sf = Ostap::Math::frexp2 (  norm() ).second ;
  Ostap::Math::scale_exp2 (  m_pars , -sf + 1 ) ;
} 
// ============================================================================
/*  construct Bernstein polynomial from its roots
 *
 *  Polinomial has a form
 *  \f$ B(x) = \prod_i (x-r_i) \prod_j (x-c_i)(x-c_i^*) \f$
 *
 *  @param xmin low  edge for Bernstein polynomial
 *  @param xmax high edge for Bernstein polynomial
 *  @param c  the list of complex roots (only one root from cc-pair is needed)
 *  @param r  the list of real  roots of the polinomial
 */
// ============================================================================
Ostap::Math::Bernstein::Bernstein 
( const double xmin , 
  const double xmax , 
  const std::vector<std::complex<double> > & c ,
  const std::vector<double>&                 r )
  : Ostap::Math::Bernstein::Bernstein ( xmin , xmax , r   , c ) 
{}
// ============================================================================
// copy assignement 
// ============================================================================
Ostap::Math::Bernstein&
Ostap::Math::Bernstein::operator=( const Ostap::Math::Bernstein&  right ) 
{
  if ( &right == this ) { return *this ; }
  m_xmin = right.m_xmin ;
  m_xmax = right.m_xmax ;
  Ostap::Math::PolySum::operator=( right ) ;
  return *this ;
}
// ============================================================================
// move assignement 
// ============================================================================
Ostap::Math::Bernstein&
Ostap::Math::Bernstein::operator=(       Ostap::Math::Bernstein&& right ) 
{
  if ( &right == this ) { return *this ; }
  m_xmin = right.m_xmin ;
  m_xmax = right.m_xmax ;
  Ostap::Math::PolySum::operator=( std::move ( right ) ) ;
  return *this ;
}
// ============================================================================
// assignement from the constant 
// ============================================================================
Ostap::Math::Bernstein&
Ostap::Math::Bernstein::operator=( const double right ) 
{
  std::fill ( m_pars.begin() , m_pars.end() , s_zero ( right ) ? 0.0 : right )  ;
  return *this ;
}
// ============================================================================
// all coefficients are so small that  c+p == c ? 
// ============================================================================
bool Ostap::Math::Bernstein::small ( const double c ) const 
{ 
  const static Ostap::Math::MuchSmaller<double> s_much_smaller{} ;
  return s_much_smaller ( norm() , c ) ;
}
// ============================================================================
// is it a increasing function?
// ============================================================================
bool Ostap::Math::Bernstein::increasing   () const 
{
  if ( m_pars.size() <= 1 ) { return true ; }
  for ( std::vector<double>::const_iterator it = m_pars.begin() + 1 ; 
        m_pars.end() != it ; ++it ) 
  { if (  (*(it-1)) > (*it) && !s_equal ( *(it-1) , *it ) ) { return false ; } }
  return true ;
}
// ============================================================================
// is it a decreasing function?
// ============================================================================
bool Ostap::Math::Bernstein::decreasing   () const 
{
  if ( m_pars.size() <= 1 ) { return true ; }
  for ( std::vector<double>::const_iterator it = m_pars.begin() + 1 ; 
        m_pars.end() != it ; ++it ) 
  { if (  (*(it-1)) < (*it) && !s_equal ( *(it-1) , *it ) ) { return false ; } }
  return true ;
}
// ============================================================================
// is it a constant function?
// ============================================================================
bool Ostap::Math::Bernstein::constant () const 
{
  //
  if ( m_pars.size() <= 1 ) { return true ; }
  for ( std::vector<double>::const_iterator it = m_pars.begin() + 1 ; 
        m_pars.end() != it ; ++it ) 
  { if ( !s_equal ( *(it-1) ,  *it ) ) { return false ; } }
  //
  return true ;
}
// ============================================================================
// get the integral between xmin and xmax
// ============================================================================
double Ostap::Math::Bernstein::integral () const
{
  return
    ( m_xmax - m_xmin ) *
    std::accumulate ( m_pars.begin() , m_pars.end() , 0.0 ) / npars() ;
}
// ============================================================================
/*  filter out very small terms
 *  the term is considered to be very small if
 *   - it is numerically zero
 *   - or if epsilon > 0,
 *          abs ( c(k) ) < epsilon
 *   - or if scale   > 0  , 
 *           scale + par ==  scale 
 *   - or if scale   <= 0 ,
 *           norm  + pars == norm    
 *  Since the maximum value for each term of
 *  \f$ c_k C^n_k \frac{ k^k (n-k)^{n-k}}{ n^n}\f$
 *  @param  epsilon  parameter to define "smalness" of terms
 *  @param  scale    parameter to define "smalness" of terms
 *  @return number of nullified terms
 */
// ============================================================================
unsigned short 
Ostap::Math::Bernstein::remove_noise ( const double epsilon , 
                                       const double scale   )
{
  unsigned short       num = 0           ;
  const unsigned short N   = degree()    ;
  const bool           eps = 0 < epsilon ;
  const double        n    = norm () ;
  for ( unsigned short k = 0 ; k <= N ; ++k ) 
  {
    if      (                        s_zero ( m_pars[k] )           ) { m_pars[k] = 0 ; ++num ; }
    else if ( eps && ( 0 == k ) && std::abs ( m_pars[k] ) < epsilon ) { m_pars[k] = 0 ; ++num ; }
    else if ( eps && ( N == k ) && std::abs ( m_pars[k] ) < epsilon ) { m_pars[k] = 0 ; ++num ; }
    else if ( 0 <  scale && s_equal ( scale + m_pars[k]   , scale ) ) { m_pars[k] = 0 ; ++num ; }
    else if ( 0 >= scale && s_equal ( n     + m_pars[k]   , n )     ) { m_pars[k] = 0 ; ++num ; }
  }
  return num ;
}    
// ============================================================================
/*  get indefinite integral  as function object 
 *  \f$ I(x) = \int^{x}_{x_{min}} B(t) dt + C \f$
 *  @param C the integration constant   
 */
// ============================================================================
Ostap::Math::Bernstein
Ostap::Math::Bernstein::indefinite_integral 
( const double C ) const 
{
  //
  std::vector<long double> ck ( npars() + 1 , 0.0 ) ;
  std::partial_sum   ( m_pars.begin () , m_pars.end   () ,  ck.begin() + 1 ) ;
  Ostap::Math::scale ( ck , ( m_xmax - m_xmin ) / npars() ) ;
  //
  // add the integration constant 
  if ( !s_zero ( C ) ) 
  {
    for ( std::vector<long double>::iterator ic = ck.begin() ; ck.end() != ic ; ++ic ) 
    { (*ic) += C ; }
  }
  //
  return Ostap::Math::Bernstein ( ck.begin() , ck.end  () , m_xmin , m_xmax ) ;
}
// ============================================================================
double Ostap::Math::Bernstein::integral ( const double low  ,
                                          const double high ) const 
{
  //
  if      ( s_equal ( low , high )           ) { return  0 ; }
  else if ( low  >  high                     ) { return -1*integral ( high   , low    ) ; }
  else if ( high <= xmin () || low >= xmax() ) { return  0 ; }
  else if ( s_vzero ( m_pars )               ) { return  0 ; }  
  else if ( s_equal ( low  , m_xmin ) && 
            s_equal ( high , m_xmax )        ) { return integral () ; }          
  //
  const double xlow  = std::max ( low  , m_xmin ) ;
  const double xhigh = std::min ( high , m_xmax ) ;
  if ( xlow > xhigh                          ) { return 0 ;}
  //
  if ( 1 == npars() ) { return ( xhigh - xlow ) * m_pars[0] ; }
  //
  if ( s_equal ( xlow  , m_xmin ) && 
       s_equal ( xhigh , m_xmax ) ) { return integral () ; }
  //
  // make integration: 
  //
  std::vector<long double> ck ( npars() + 1 , 0.0 ) ;
  std::partial_sum ( m_pars.begin () , m_pars.end   () ,  ck.begin() + 1 ) ;
  Ostap::Math::scale ( ck , ( m_xmax - m_xmin ) / npars() ) ;
  //
  const Ostap::Math::Bernstein b_int ( ck.begin() ,
                                       ck.end  ()  , m_xmin , m_xmax ) ;
  //
  return b_int ( xhigh ) - b_int ( xlow ) ;
}
// ============================================================================
Ostap::Math::Bernstein
Ostap::Math::Bernstein::derivative () const 
{
  //
  std::vector<long double>   ck ( npars() , 0 ) ;
  std::adjacent_difference ( m_pars.begin () , m_pars.end() , ck.begin() ) ;
  Ostap::Math::scale ( ck , ( npars() - 1 )/ ( m_xmax - m_xmin ) ) ;
  //
  return Ostap::Math::Bernstein ( ck.begin() + 1 , ck.end() ,  m_xmin  , m_xmax ) ;
}
// ============================================================================
double Ostap::Math::Bernstein::derivative ( const double x   ) const 
{
  if      ( m_pars.size() <= 1       ) { return 0 ; }
  else if ( x < m_xmin || x > m_xmax ) { return 0 ; }
  //
  std::vector<long double>   ck ( npars() , 0 ) ;
  std::adjacent_difference ( m_pars.begin () , m_pars.end() , ck.begin() ) ;
  //
  // get the t-values
  //
  const double t0 = t ( x ) ;
  const double t1 = 1 - t0  ;
  //
  return
    _casteljau_ ( ck.begin() + 1 , ck.end() , t0 , t1 ) * ( npars()-1 )  / ( m_xmax - m_xmin ) ;
}
// ============================================================================
// get the value
// ============================================================================
double Ostap::Math::Bernstein::evaluate ( const double x ) const
{
  //
  // treat the trivial cases
  //
  if      ( m_pars.empty()                             ) { return 0 ; }
  //  needed for the proper integration with an exponential 
  else if ( s_equal ( x , m_xmin )                     ) { return m_pars [0]    ; }
  else if ( s_equal ( x , m_xmax )                     ) { return m_pars.back() ; }
  else if ( 1 == npars ()                              ) { return m_pars [0]    ; }
  else if ( s_vzero ( m_pars )                         ) { return 0 ; }
  //
  // get the t-values
  //
  const long double t0 = t ( x ) ;
  const long double t1 = 1 - t0  ;
  //
  // start de casteljau algorithm:
  //
  // use fixed size: 
  if (  npars() < 16 ) 
  {
    std::array<long double,16> _pars;
    std::copy( m_pars.begin() , m_pars.end() , _pars.begin() ) ;
    return _casteljau_ ( _pars.begin() , _pars.begin() + npars() , t0 , t1 ) ;
  }
  // generic case:
  std::vector<long double> dcj ( m_pars.begin() , m_pars.end() ) ;
  return _casteljau_ ( dcj.begin() , dcj.end() , t0 , t1 ) ;
}
// ============================================================================
Ostap::Math::Bernstein&
Ostap::Math::Bernstein::operator+=( const double a ) 
{
  if   ( s_zero ( a ) ) { return *this ; }
  Ostap::Math::shift ( m_pars , a ) ;
  return *this ;
}
// ============================================================================
Ostap::Math::Bernstein&
Ostap::Math::Bernstein::operator*=( const double a ) 
{
  if      ( s_equal ( a , 1 ) ) { return *this ; }
  else if ( s_zero  ( a     ) ) { std::fill ( m_pars.begin() , m_pars.end() , 0 ) ; }
  Ostap::Math::scale ( m_pars , a ) ;
  return *this ;
}
// ============================================================================
Ostap::Math::Bernstein&
Ostap::Math::Bernstein::operator-=( const double a ) 
{
  if ( s_zero ( a ) ) { return *this ; }
  Ostap::Math::shift ( m_pars , -a ) ;
  return *this ;
}
// ============================================================================
Ostap::Math::Bernstein&
Ostap::Math::Bernstein::operator/=( const double a ) 
{
  if   ( s_equal ( a , 1 ) ) { return *this ; }
  Ostap::Math::scale ( m_pars , 1/a ) ;
  return *this ;
}
// ============================================================================
Ostap::Math::Bernstein
Ostap::Math::Bernstein::operator-() const 
{
  Bernstein b ( *this ) ;
  Ostap::Math::negate ( b.m_pars ) ;
  return b ;
}
// ============================================================================
// Sum of Bernstein polynomial and a constant 
// ============================================================================
Ostap::Math::Bernstein
Ostap::Math::Bernstein::__add__   ( const double value ) const 
{ return (*this) + value ; }
// ============================================================================
// Sum of Bernstein polynomial and a constant 
// ============================================================================
Ostap::Math::Bernstein
Ostap::Math::Bernstein::__radd__  ( const double value ) const 
{ return value + (*this) ; }
// ============================================================================
// Product of Bernstein polynomial and a constant
// ============================================================================
Ostap::Math::Bernstein
Ostap::Math::Bernstein::__mul__   ( const double value ) const 
{ return (*this) * value ; }
// ============================================================================
// Product of Bernstein polynomial and a constant
// ============================================================================
Ostap::Math::Bernstein
Ostap::Math::Bernstein::__rmul__  ( const double value ) const 
{ return value * (*this) ; }
// ============================================================================
// Subtract a constant from Benrstein polynomial
// ============================================================================
Ostap::Math::Bernstein
Ostap::Math::Bernstein::__sub__  ( const double value ) const 
{ return (*this) - value ; }
// ============================================================================
// Subtract Bernstein polynomial from a constant 
// ============================================================================
Ostap::Math::Bernstein
Ostap::Math::Bernstein::__rsub__ ( const double value ) const 
{ return value - (*this) ; }
// ============================================================================
// Divide Benrstein polynomial by a constant 
// ============================================================================
Ostap::Math::Bernstein
Ostap::Math::Bernstein:: __div__   ( const double value ) const 
{ return (*this) / value ; }
// ============================================================================
// Negate Bernstein polynomial 
// ============================================================================
Ostap::Math::Bernstein
Ostap::Math::Bernstein::__neg__ ()  const 
{ return -(*this); }
// ============================================================================
// the sum two Bernstein polynomials
// ============================================================================
Ostap::Math::Bernstein  
Ostap::Math::Bernstein::sum ( const Ostap::Math::Bernstein& other ) const 
{
  if ( this == &other ) 
  {
    Bernstein result( *this );
    result *=2 ;
    return result ;
  }
  //
  if ( !s_equal ( xmin () , other.xmin() ) || 
       !s_equal ( xmax () , other.xmax() ) ) 
  {
    const double x_min = std::min ( xmin() , other.xmin() ) ;
    const double x_max = std::max ( xmax() , other.xmax() ) ;
    Bernstein b1 ( *this , x_min , x_max ) ;
    Bernstein b2 ( other , x_min , x_max ) ;
    return b1.sum ( b2 ) ;
  }
  //
  if ( degree() < other.degree() )
  { return other.sum ( this->elevate ( other.degree() - degree       () ) ) ; }
  if ( degree() > other.degree() ) 
  { return       sum ( other.elevate (       degree() - other.degree () ) ) ; }
  //
  Bernstein result(*this) ;
  for ( unsigned short i = 0 ; i < npars() ; ++i ) 
  { result.m_pars[i] += other.par( i ) ; }
  return result ; 
}
// =============================================================================
// subtract Bernstein polynomials (with the same domain!)
// =============================================================================
Ostap::Math::Bernstein  
Ostap::Math::Bernstein::subtract ( const Ostap::Math::Bernstein& other ) const
{
  if ( this == &other ) { return Bernstein( degree() , xmin() , xmax() ) ; }
  Bernstein b ( other ) ;
  Ostap::Math::negate ( b.m_pars ) ;
  return sum ( b ) ;  
}
// ============================================================================
// Sum of Bernstein polynomials (the same domain)
// ============================================================================
Ostap::Math::Bernstein
Ostap::Math::Bernstein::__add__   
( const  Ostap::Math::Bernstein& other ) const { return sum ( other ) ; }
// ============================================================================
// Subtraction of Bernstein polynomials (the same domain)
// ============================================================================
Ostap::Math::Bernstein
Ostap::Math::Bernstein::__sub__   
( const  Ostap::Math::Bernstein& other ) const { return subtract ( other ) ; }
// ============================================================================
// Multipky twp Bernstein polynomials (the same domain)
// ============================================================================
Ostap::Math::Bernstein
Ostap::Math::Bernstein::__mul__   
( const  Ostap::Math::Bernstein& other ) const { return multiply ( other ) ; }
// ============================================================================
namespace 
{
  // ==========================================================================
  inline long double
  c_nk ( const unsigned short n , 
         const unsigned short k ) 
  {
    return 
      n < 63 ? 
      Ostap::Math::choose        ( n , k ) : 
      Ostap::Math::choose_double ( n , k ) ;  
  }
  // ==========================================================================
}
// ============================================================================
/*  elevate it: 
 *  represent as Bernstein polynomial of order N+r 
 *  @param r  INPUT increase of degree 
 *  @return new polynomial of order N+r 
 */
// ============================================================================
Ostap::Math::Bernstein
Ostap::Math::Bernstein::elevate ( const unsigned short r ) const 
{
  // no need in elevation 
  if ( 0 == r ){ return *this ; }
  //
  std::vector<long double>    _nc ( npars () + r ) ; // new coefficients
  const std::vector<double>&  _oc =  pars ()       ; // old coefficients 
  // copy it 
  const unsigned short N  = degree() ;
  //
  std::copy ( _oc.begin() , _oc.end () , _nc.begin()              ) ;
  std::fill ( _nc.begin() + _oc.size() , _nc.end  () , _oc.back() ) ;
  //
  // repeate the elevation cycles: 
  for ( unsigned short   n = N  ; n <  N + r  ; ++n ) 
  {
    // "current" degree 
    for ( unsigned short k = n ;  1<= k ; --k ) 
    {
      _nc[k]  = ( n + 1 - k ) * _nc[k] + k * _nc[k-1] ;
      _nc[k] /=   n + 1  ;
    }    
  }
  //
  return Bernstein ( _nc.begin() , _nc.end  () , xmin() , xmax() ) ;
}
// ============================================================================
/*  reduce it
 *  represent as Bernstein polynomial of order N-r 
 *  @param r  INPUT increase of degree 
 *  @return new polynomial of order N-r 
 */
// ============================================================================
Ostap::Math::Bernstein
Ostap::Math::Bernstein::reduce ( const unsigned short r ) const 
{
  // no need in reducing
  if ( 0 == r || 0 == degree() ){ return *this ; }
  //
  const unsigned short n    = degree () ;
  const unsigned short newd = r <= n ?  n - r : 0 ;
  //
  const std::vector<double>&  oc =  pars ()                ; // old coefficients 
  std::vector<long double>    nc ( oc.begin() , oc.end() ) ; // new coefficients
  //
  for ( unsigned short n = degree() ; newd < n ; --n ) 
  {
    for ( unsigned short k = 1 ; k < n ; ++k ) 
    {  nc[k] =  ( n * nc[k] - k * nc[k-1] ) / ( n - k ) ;  }
  }
  return Bernstein ( nc.begin() , nc.begin() + newd + 1 , xmin() , xmax() ) ;
}
// ============================================================================
namespace 
{
  //
  template <class ITERATOR>
  inline long double _head_ ( ITERATOR first , ITERATOR last  ) 
  {
    if ( first == last ) { return 0      ; }
    const unsigned short N = std::distance ( first , last ) - 1 ;
    //
    long double  c = 1 ;
    unsigned int i = 0 ;
    long double  h = 0 ;
    --last ;
    const unsigned int N2 = N / 2 + 1 ;
    for ( ; i < N2  ; ++first, --last, ++i ) 
    {
      if ( 0 < i ) { c *= ( N + 1 - i ) ; c /= i ; }
      if ( first == last ) 
      {
        h +=   i%2 ? c * ( *first) : -c * (*first ) ;      
        break ;
      }
      h +=    i %2 ? c * ( *first) : -c * (*first ) ;      
      h += (N-i)%2 ? c * ( *last ) : -c * (*last  ) ;
    }
    return h * ( 0 == N%2 ? -1 : 1 ) ;
  }
  inline long double _head_ ( const std::vector<double>& pars  )
  { return _head_ ( pars.begin() , pars.end() ) ; }
  inline long double _head_ ( const Ostap::Math::Bernstein& b ) 
  { return _head_ ( b.pars() ) ; }
}
// ============================================================================
/*  calculate ``nearest'' polynomial (in the sense of q-norm) of lower degree, 
 *  where q-norm is defined as:
 *  \f$ \left| f \right|_{q} = \left( \sum_i \left|c_i\right|^q\right)^{\frac{1}{q}} \f$
 *  
 *  - q_inv = 0.0 ->  \f$ max_k    \left|c_k\right|  \f$ 
 *  - q_inv = 0.5 ->  \f$ \sqrt{ \sum_k  c_k^2 }     \f$
 *  - q_inv = 1.0 ->  \f$ \sum_k \left| c_k \right|  \f$ 
 *  @see  R.M. Corless & N.Rezvani,
 *        "The nearest polynomial of lower degree",
 *        Proceedings of the 2007 international workshop on 
 *        Symbolic-numeric computation SNC'07 
 *        https://cs.uwaterloo.ca/conferences/issac2007/
 *  @see http://dl.acm.org/citation.cfm?id=1277530&CFID=799220770&CFTOKEN=25289921
 *  @see  N.Rezvani and R.M. Corless, 
 *       "The Nearest Polynomial With A Given Zero, Revisited"
 *        ACM SIGSAM Bulletin, Vol. 39, No. 3, September 2005
 *  @see http://dl.acm.org/citation.cfm?doid=1113439.1113442
 */
// ============================================================================
Ostap::Math::Bernstein
Ostap::Math::Bernstein::nearest ( const double qinv ) const 
{
  //  trivial case: 
  if ( 1 > degree() ) { return *this ; }
  //
  // get the norm
  const double iq = 0 > qinv ? 0 : 1 < qinv ? 1 : qinv ;
  // 
  const double ip = 1 - iq ;
  //
  // leading coefficients of the basis 
  //
  std::vector<long double>  u ( m_pars.size() ) ;
  u [0] = 1  ;
  const unsigned short N = m_pars.size() ;
  for ( unsigned short i = 1 ; i < N ;  ++i ) { u [i]  = ( u[i-1] * ( N  - i ) ) / i ; }
  for ( unsigned short i = 0 ; i < N ;  ++i ) 
  { if  ( 1 == ( N + 1 - i )%2 ) { u [i] *= -1 ; } }
  //
  const long double un = 1 / Ostap::Math::p_norm ( u.begin() , u.end() , ip );
  //
  long double uc = 0  ;
  for ( unsigned short i = 0 ; i < N ; ++i ) 
  {
    u[i] *= un ; // normalize it!
    uc   += u[i] * m_pars[i] ; 
  }
  //
  const signed char suc = Ostap::Math::signum ( uc ) ;
  //
  std::vector<long double>  v ( m_pars.size() ) ;
  const  bool p_inf = s_zero ( ip ) ;
  if ( !p_inf ) 
  {
    for ( unsigned short k = 0 ; k < N ; ++k ) 
    {
      const double uk = u[k] ;
      v[k] = uc * u[k] * std::pow ( std::abs ( uk ) , 1/ip - 2 ) ;
    }
  }
  else if ( 1 == N%2 ) 
  { const unsigned short k0 = (N-1)/2 ; v[ k0 ] =       uc * u[k0] ; }
  else if ( 0 == N%2 ) 
  { 
    const unsigned short k1 =  N/2    ; v[ k1 ] = 0.5 * uc * u[k1] ; 
    const unsigned short k2 =  N/2-1  ; v[ k2 ] = 0.5 * uc * u[k2] ; 
  }
  //
  std::vector<long double> nc ( m_pars.begin() , m_pars.end  () ) ;  
  for ( unsigned short i = 0 ; i < N ; ++i ) { nc[i] -= v[i] ; }
  //
  // return Bernstein ( nc.begin() , nc.end() , xmin() , xmax() ) ;
  //
  // reduce the degree:
  const unsigned short n  = degree() ;
  const unsigned short nd = 1 >= n ? 0 : n - 1 ;
  for ( unsigned short k  = 1 ; k < n ; ++k ) 
  {  nc[k] =  ( n * nc[k] - k * nc[k-1] ) / ( n - k ) ;  }
  //
  return Bernstein ( nc.begin() , nc.begin() + nd + 1 , xmin() , xmax() ) ;
}
// ============================================================================
/*  calculate q-norm of the polynomial 
 *  where q-norm is defined as:
 *  \f$ \left| f \right|_{q} = \left( \sum_i \left|c_i\right|^q\right)^{\frac{1}{q}} \f$
 *  
 *  - q_inv = 0.0 ->  \f$ max_k    \left|c_k\right|  \f$ 
 *  - q_inv = 0.5 ->  \f$ \sqrt{ \sum_k  c_k^2 }     \f$
 *  - q_inv = 1.0 ->  \f$ \sum_k \left| c_k \right|  \f$ 
 */
// ============================================================================
double Ostap::Math::Bernstein::norm   ( const double q_inv ) const 
{ return Ostap::Math::p_norm ( m_pars.begin() , m_pars.end() , q_inv ) ; }
// ============================================================================
/*  how close are two polynomials in q-norm?
 *  where q-norm is defined as:
 *  \f$ \left| f \right|_{q} = \left( \sum_i \left|c_i\right|^q\right)^{\frac{1}{q}} \f$
 *  
 *  - q_inv = 0.0 ->  \f$ max_k    \left|c_k\right|  \f$ 
 *  - q_inv = 0.5 ->  \f$ \sqrt{ \sum_k  c_k^2 }     \f$
 *  - q_inv = 1.0 ->  \f$ \sum_k \left| c_k \right|  \f$        
 */
// ============================================================================
double Ostap::Math::Bernstein::distance 
( const Ostap::Math::Bernstein& other , const double q_inv ) const 
{
  if ( &other == this ) { return 0 ; }
  //
  // 1) adjust the ranges 
  //
  if ( !s_equal ( xmin () , other.xmin() ) || !s_equal ( xmax () , other.xmax() ) ) 
  { return distance ( Bernstein ( other, xmin() , xmax() ) , q_inv ) ; }
  //
  // 2) adjust the degrees 
  //
  if      ( degree() > other.degree() ) 
  { return       distance ( other.elevate (       degree() - other.degree() ) , q_inv ) ; }
  else if ( degree() < other.degree() ) 
  { return other.distance ( this->elevate ( other.degree() -       degree() ) , q_inv ) ; }
  //
  // 3) make a real comparsion 
  //
  std::vector<long double> v ( m_pars.begin() , m_pars.end() ) ;
  const unsigned short N = degree() ;
  for ( unsigned short k = 0 ; k <= N ; ++k ) { v[k] -= other.m_pars[k] ; }
  //
  return Ostap::Math::p_norm( v.begin() , v.end() , q_inv ) ; 
}
// ============================================================================
// multiply two Bernstein polynomials
// ============================================================================
Ostap::Math::Bernstein  
Ostap::Math::Bernstein::multiply ( const Ostap::Math::Bernstein& other ) const 
{
  //
  if ( !s_equal ( xmin () , other.xmin() ) || 
       !s_equal ( xmax () , other.xmax() ) ) 
  {
    const double x_min = std::min ( xmin() , other.xmin() ) ;
    const double x_max = std::max ( xmax() , other.xmax() ) ;
    Bernstein b1 ( *this , x_min , x_max ) ;
    Bernstein b2 ( other , x_min , x_max ) ;
    return b1.multiply ( b2 ) ;
  }
  //
  if ( zero() || other.zero() ) { return Bernstein( degree() , xmin() , xmax() ) ; }
  //
  const unsigned short m =       degree() ;
  const unsigned short n = other.degree() ;
  //
  Bernstein result ( m + n , xmin() , xmax() ) ;
  //
  long double c = 1 ;
  for ( unsigned short k = 0 ; k <= m + n ; ++k ) 
  {
    if ( 0 != k ) {  c *= ( m + n - k + 1 ) ; c /= k; }
    //
    const unsigned jmax = std::min ( m , k ) ;
    const unsigned jmin = k > n ? k - n : 0 ;
    long double     cc  = 0 == jmin ? 
      c_nk ( n , k - jmin ) :
      c_nk ( m ,     jmin ) ;  
    for ( unsigned short j = jmin ; j <= jmax ; ++j ) 
    {
      if ( j != jmin ) { cc *= ( m - j + 1 )  * ( k - j + 1 ) ; cc /= j * ( n - k + j ) ; }
      result.m_pars[k] += cc * m_pars [ j]  * other.m_pars[k-j] ;
    }
    result.m_pars[k] /= c ;
  }
  return result ; 
}
// ============================================================================
// multiply two Bernstein polynomial and the basic bernstein polynomial 
// ============================================================================
Ostap::Math::Bernstein  
Ostap::Math::Bernstein::multiply 
( const Ostap::Math::Bernstein::Basic& b ) const 
{
  Bernstein   result ( multiply ( b.k() , b.N() - b.k() ) ) ;
  Ostap::Math::scale ( result.m_pars , Ostap::Math::choose ( b.N() , b.k() ) ) ;
  return result ;
}
// ============================================================================
/*  multiply Bernstein polynomial with 
 *  \f$ (x-x_{min})^i(x_{max}-x)^j \f$ 
 */
// ============================================================================
Ostap::Math::Bernstein  
Ostap::Math::Bernstein::multiply 
( const unsigned short i1 , 
  const unsigned short i2 ) const 
{
  //
  const unsigned short m = i1 + i2  ;
  const unsigned short n = degree() ;
  //
  Bernstein result ( n + m , xmin() , xmax() ) ;
  //
  const unsigned short nK = result.m_pars.size() ;
  for ( unsigned short k  = i1 ; k < nK ; ++k ) 
  {
    const unsigned short imin = k > n ? k - n : 0 ;
    const unsigned short imax = std::min ( m , k ) ;
    if ( imin <= i1 && i1 <= imax ) 
    {
      result.m_pars[k] = m_pars[k - i1] * 
        Ostap::Math::choose (     n , k - i1 ) / 
        Ostap::Math::choose ( m + n , k      ) ;      
    }
  }
  return result ;
}
// ============================================================================
namespace 
{
  // ==========================================================================
  // power function
  // ==========================================================================
  inline Ostap::Math::Bernstein _pow_ 
  ( const Ostap::Math::Bernstein& x , 
    const unsigned short          y ,  
    const Ostap::Math::Bernstein& r ) 
  { 
    return 
      0 == y ? r :
      1 == y ? ( x.degree() >= r.degree() ? x.multiply ( r ) : r.multiply ( x ) ) :
      _pow_    ( x.multiply ( x ) , y/2 , y%2 ? r * x : r ) ; 
  }
}
// ============================================================================
// power function
// ============================================================================
Ostap::Math::Bernstein
Ostap::Math::Bernstein::pow ( const unsigned short i ) const 
{
  if      ( 1 == i ) { return *this             ; }
  else if ( 2 == i ) { return multiply (*this ) ; }
  //
  Ostap::Math::Bernstein one ( 0 , xmin() , xmax() ) ;
  one.m_pars[0] = 1 ;
  return _pow_ ( *this , i , one ) ;
}
// ============================================================================
// scale all coefficients with 2**i
// ============================================================================
Ostap::Math::Bernstein  
Ostap::Math::Bernstein::ldexp ( const short i )  const 
{
  if ( 0 == i ) { return *this ; }
  Ostap::Math::Bernstein result (*this) ;
  Ostap::Math::scale_exp2 ( result.m_pars , i ) ;
  return result ;
}
// ======================================================================
namespace 
{
  // ==========================================================================
  // implmentation of polynomial division 
  // - operator Head  
  // ==========================================================================
  template <class ITERATOR>
  inline double _m_head_ 
  ( const unsigned short m     ,
    ITERATOR             first ,
    ITERATOR             last  ) 
  {
    if ( first == last ) { return 0 ; }
    //
    long double h = 0 ;
    long double c = 1 ;
    for ( unsigned short i = 0 ; first != last ; ++first , ++i ) 
    {
      if ( 0 != i ) { c *= ( m + 1 - i ) ; c /= i ; }
      h += ( 0 == i%2 ) ?  c * (*first) : -c * (*first) ;
    }
    return h ;
  }
  // ==========================================================================
  // - operator Tail
  // ==========================================================================
  template <class ITERATOR, class OUTPUT>
  inline OUTPUT _tail_ ( ITERATOR first  , 
                         ITERATOR last   , 
                         OUTPUT   output )
  {
    if ( first == last ) { return output ; } 
    const unsigned short m = std::distance ( first , last ) - 1 ;
    long double c = 1 ;
    for ( unsigned short j = 0 ; j < m ; ++j , ++output ) 
    {
      if ( 0 != j ) { c *= j ; c /= ( m - j ) ; }
      const long double t = c * _m_head_ ( m , first , first + ( j + 1 ) ) ;
      *output = 0== j%2 ? t : -t ;
    }
    return output ;
  }
  // ==========================================================================
  // - operator Match 
  // ==========================================================================
  template <class ITERATOR, class OUTPUT>
  inline OUTPUT _match_m_ ( const unsigned short m      , 
                            ITERATOR             first  , 
                            ITERATOR             last   , 
                            OUTPUT               output )
  {
    if ( first == last ) { return output ; } 
    const unsigned short n = std::distance ( first , last ) - 1 ;
    //
    long double c = 1 ;
    for ( unsigned short j = 0 ; j <= n ; ++j , ++first, ++output   )
    {
      //
      if ( 0 != j ) {  c *= ( n - j + 1 ) ; c /= ( m - j + 1 ) ; }
      *output = (*first) * c ;
    }
    //
    for ( unsigned j = n + 1 ; j <= m ; ++j, ++output ) { *output = 0 ; }
    //
    return output ;  
  }  
  // ==========================================================================
  // - operator Quot
  // ==========================================================================
  template <class ITERATOR, class OUTPUT>
  inline OUTPUT _quot_k_ ( const unsigned short k         , 
                           const unsigned short m         , 
                           ITERATOR             first     , 
                           ITERATOR             last      , 
                           OUTPUT               output    , 
                           const long double    scale = 1 )
  {
    if ( first == last ) { return output ; } 
    const unsigned short n  = std::distance ( first , last ) - 1 ;
    const unsigned short k1 = k - ( m - n ) ;
    //
    long double c = scale ;
    for ( unsigned short j = 0 ; j <= k1 ; ++j , ++first, ++output   )
    {
      //
      if ( 0 != j ) { c *= ( k1 - j + 1 ) ; c /= ( k - j + 1 )  ; }
      *output += c ;
    }
    //
    for ( unsigned j = k1 + 1 ; j <= k ; ++j , ++output ) { *output += 0;  }
    //
    return output ;  
  }  
  // ==========================================================================
  template <class ITERATORF, class ITERATORG>
  inline std::vector<long double> _divmod_ 
  ( ITERATORF f_first , 
    ITERATORF f_last  , 
    ITERATORG g_first , 
    ITERATORG g_last  ) 
  {
    //
    const unsigned short m  = std::distance ( f_first , f_last ) - 1 ;
    const unsigned short n  = std::distance ( g_first , g_last ) - 1 ;
    //
    std::vector<long double> _tail  ( m     + 1 , 0.0 ) ;
    std::vector<long double> _match ( m     + 1 , 0.0 ) ;
    std::vector<long double> _quot  ( m - n + 1 , 0.0 ) ;
    //
    typedef std::vector<long double>::iterator IT ;
    //
    ITERATORF fend         = f_last ;
    for ( unsigned short i = m ;  n <= i ; --i )
    {
      const long double h1 = _head_  ( f_first , fend ) ;
      //
      if ( !s_zero ( h1 ) ) 
      {
        // 
        IT mend = _match_m_ ( i , g_first , g_last , _match .begin() ) ;
        //
        const long double h2 = _head_ ( _match.begin() , mend ) ;
        //
        _quot_k_ ( m - n , i , g_first , g_last , _quot.begin() , h1 / h2 ) ;
        //
        for ( unsigned short j = 0 ; j < m + 1 ; ++j ) 
        { *(f_first+j) -= h1 * _match[j] / h2 ; }
        //
      }
      //
      _tail_    ( f_first , fend , _tail.begin()        ) ;
      --fend ;
      std::copy ( _tail.begin() , _tail.end() , f_first ) ;
      //
    }
    return _quot ;
    // ==========================================================================
  }
}
// ==============================================================================
// the leading power coefficient 
// ==============================================================================
double Ostap::Math::Bernstein::head   () const 
{ return _head_ ( m_pars.begin() , m_pars.end() ) ; }
// // ==============================================================================
// Ostap::Math::Bernstein 
// Ostap::Math::Bernstein::tail  () const 
// {
//   if ( 0 == degree() ) { return Bernstein ( 0 , xmin() , xmax() ) ; }
//   Bernstein result ( degree() - 1 , xmin() , xmax() ) ;
//   _tail_ ( m_pars.begin () , 
//            m_pars.end   () , 
//            result.m_pars.begin() ) ; 
//   return result ;  
// }
// // ============================================================================
// Ostap::Math::Bernstein
// Ostap::Math::Bernstein::match ( const unsigned short m ) const 
// {
//   if ( m < degree() ) { return Bernstein ( 0 , xmin() , xmax() ) ; }
//   Bernstein result ( m , xmin() , xmax() ) ;
//   _match_m_ ( m               , 
//               m_pars.begin () , 
//               m_pars.end   () , 
//               result.m_pars.begin() ) ;
//   return result ;
// }
// // ============================================================================
// Ostap::Math::Bernstein
// Ostap::Math::Bernstein::quot  ( const unsigned short k , 
//                                 const unsigned short m ) const 
// {
//   if ( m <     degree() ) { return Bernstein ( 0 , xmin() , xmax() ) ; }
//   if ( k < m - degree() ) { return Bernstein ( 0 , xmin() , xmax() ) ; }
//   Bernstein result ( k , xmin() , xmax() ) ;
//   _quot_k_ ( k               , 
//              m               , 
//              m_pars.begin () , 
//              m_pars.end   () , 
//              result.m_pars.begin() ) ; 
//   return result ;  
// }
// 
// ==============================================================================
/* polynomial division 
 * \f$  f(x) = q(z)*g(x) + r(x) \f$ 
 * @return the pair q,r 
 */
// ==============================================================================
namespace 
{
  inline 
  std::pair<Ostap::Math::Bernstein,Ostap::Math::Bernstein>
  _divmod_ ( Ostap::Math::Bernstein f , 
             Ostap::Math::Bernstein g )
  {
    using namespace Ostap::Math ;
    ///  trivial case 
    /// 1) f==0  or |f| << |g| 
    if ( f.zero () || f.small ( g.norm() ) )
    { return std::make_pair ( Bernstein ( 0 , g.xmin() , g.xmax () ) ,
                              Bernstein ( 0 , g.xmin() , g.xmax () ) ) ; }
    /// 2) g==0  or |g| << |f| 
    if ( g.zero () || g.small ( f.norm() ) ) 
    { return std::make_pair ( Bernstein ( 0 , f.xmin() , f.xmax () ) ,
                              Bernstein ( 0 , f.xmin() , f.xmax () ) ) ; }
    //
    if ( !s_equal( f.xmin() , g.xmin() ) || !s_equal( f.xmin() , g.xmin() ) ) 
    {
      const double xmin = std::min ( f.xmin () , g.xmin () ) ;
      const double xmax = std::max ( f.xmax () , g.xmax () ) ;
      return _divmod_ ( Ostap::Math::Bernstein ( f , xmin , xmax ) , 
                        Ostap::Math::Bernstein ( g , xmin , xmax ) ) ;
    }
    //
    // get the leading coefficient of "f"
    //
    double fn = f.norm()  ;
    while ( 0 < f.degree() && s_equal (  fn + _head_ ( f ) , fn ) ) 
    {
      f  = f.reduce  ( 1 ) ;
      f.remove_noise (   ) ;
      fn = f.norm    (   ) ;
    }
    //
    // get the leading coefficient of "g"
    //
    double gn = g.norm()  ;
    while ( 0 < g.degree() && s_equal (  gn + _head_ ( g ) , gn ) ) 
    {
      g  = g.reduce  ( 1 ) ;
      g.remove_noise (   ) ;
      gn = g.norm    (   ) ;
    }
    //
    const std::vector<double>& pf  = f.pars() ;
    const std::vector<double>& pg  = g.pars() ;
    const double lc_f = _head_ ( pf.begin () , pf.end () ) ;
    const double lc_g = _head_ ( pg.begin () , pg.end () ) ;
    //
    const unsigned short m = f.degree () ;
    const unsigned short n = g.degree () ;
    //
    if ( m < n  ) { return std::make_pair ( Bernstein ( 0 , f.xmin() , f.xmax () ) , f ) ; }
    //  
    std::vector<long double> _f ( pf.begin() , pf.end() ) ;
    if ( n == 0 ) 
    {
      Ostap::Math::scale ( _f , 1/g.par(0) ) ;
      return std::make_pair ( Bernstein ( _f.begin() , _f.end  () , f.xmin() , f.xmax() ) , 
                              Bernstein ( 0                       , f.xmin() , f.xmax() ) ) ;
    }
    //
    std::vector<long double> _q =
      _divmod_ ( _f.begin() , _f.end() , 
                 pg.begin() , pg.end() ) ;           
    //
    Bernstein q ( _q.begin() , _q.end  ()     , f.xmin() , f.xmax() ) ;
    Bernstein r ( _f.begin() , _f.begin() + n , f.xmin() , f.xmax() ) ;
    //
    const double qn = q.norm() ;
    const double rn = r.norm() ;
    //
    if ( r.small ( fn + qn * gn ) ) 
    { r = Bernstein( 0 , r.xmin() , r.xmax() ) ; }
    else         
    { r.remove_noise ( 0 , fn + qn * gn )      ; }
    //
    if ( s_equal ( qn * gn + fn + rn , fn + rn ) ) 
    { q = Bernstein( 0 , q.xmin() , q.xmax() ) ; }
    else   
    { q.remove_noise ( 0 , ( fn + rn ) / gn )  ; }
    //
    return std::make_pair ( q , r ) ;
  }
  /// scale the exponents for Bernstein polynomial
  inline Ostap::Math::Bernstein
  scale_exp2 ( const Ostap::Math::Bernstein& b    , 
               const int                     iexp ) 
  {
    if ( 0 ==  iexp ) { return b ; }
    std::vector<double> pars ( b.pars() ) ;
    Ostap::Math::scale_exp2  ( pars , iexp ) ;
    return Ostap::Math::Bernstein ( std::move( pars ) , b.xmin() , b.xmax() ) ;
  } 
}
// ============================================================================
/*  polynomial division 
 *  \f$  f(x) = q(z)*g(x) + r(x) \f$ 
 *  @return the pair q(x),r(x)
 */
// ============================================================================
std::pair<Ostap::Math::Bernstein,Ostap::Math::Bernstein> 
Ostap::Math::Bernstein::divmod   ( const Ostap::Math::Bernstein& g ) const 
{ return _divmod_ ( *this , g ) ; }
// ============================================================================
/*  polynomial division 
 *  \f$  f(x) = q(z)*g(x) + r(x) \f$ 
 *  @return the quotient q(x)  
 */
// ============================================================================
Ostap::Math::Bernstein
Ostap::Math::Bernstein::quotient ( const Ostap::Math::Bernstein& g ) const
{ return _divmod_ ( *this,  g ) . first ; }
// ============================================================================
/*  polynomial division 
 *  \f$  f(x) = q(z)*g(x) + r(x) \f$ 
 *  @return the reminder r(x)
 */
// ============================================================================
Ostap::Math::Bernstein
Ostap::Math::Bernstein::reminder ( const Ostap::Math::Bernstein& g ) const
{ return _divmod_ ( *this , g ) . second ; }
// ============================================================================
/* de Casteljau algorithm for summation of Bernstein polynomials 
 *  \f$ f(x) = \sum_i p_i B_ik(x) \f$
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-02-10
 */
// ============================================================================
double Ostap::Math::casteljau
( const std::vector<double>& pars , 
  const double               x    ) 
{
  std::vector<long double> _tmp ( pars.begin() , pars.end () ) ;
  //
  const long double t0 =     x  ;
  const long double t1 = 1 - t0 ;
  //
  return _casteljau_ ( _tmp.begin() , _tmp.end  () , t0 , t1 ) ;
}
// ============================================================================
namespace 
{
  // ==========================================================================
  /** transformation matrix from legendre to bernstein basis 
   *  @see http://www.sciencedirect.com/science/article/pii/S0377042700003769 eq.20 
   */
  inline 
  double l2b_mtrx ( const unsigned short j , 
                    const unsigned short k ,
                    const unsigned short n ) 
  {
    //
    const unsigned short imin = std::max ( 0 , j + k - n ) ;
    const unsigned short imax = std::min ( j ,     k     ) ;
    long long r = 0 ;
    for ( unsigned short i = imin ; i <= imax ; ++i ) 
    {
      0 == ( k + i ) % 2 ? 
        r +=
        Ostap::Math::choose ( j     , i     ) * 
        Ostap::Math::choose ( k     , i     ) * 
        Ostap::Math::choose ( n - j , k - i ) :
        r -=
        Ostap::Math::choose ( j     , i     ) * 
        Ostap::Math::choose ( k     , i     ) * 
        Ostap::Math::choose ( n - j , k - i ) ;
    }
    //
    return r / double ( Ostap::Math::choose ( n , k ) ) ;
  }
  // ==========================================================================
  /** transformation matrix from chebyshev to bernstein basis 
   *  http://www.degruyter.com/view/j/cmam.2003.3.issue-4/cmam-2003-0038/cmam-2003-0038.xml  eq. 15
   */
  inline 
  long double c2b_mtrx ( const unsigned short j , 
                         const unsigned short k ,
                         const unsigned short n ) 
  {
    const unsigned short imin = std::max ( 0 , j + k - n ) ;
    const unsigned short imax = std::min ( j ,     k     ) ;
    long long r = 0 ;
    for ( unsigned short i = imin ; i <= imax ; ++i ) 
    {
      0 == ( k - i ) % 2 ? 
        r +=
        Ostap::Math::choose ( 2 * k , 2 * i ) * 
        Ostap::Math::choose ( n - k , j - i ) :
        r -=
        Ostap::Math::choose ( 2 * k , 2 * i ) * 
        Ostap::Math::choose ( n - k , j - i ) ;
    }
    //
    return r / (  (long double)  Ostap::Math::choose ( n , j ) ) ;
  }
  // ==========================================================================
  /** transformation matrix from monomial to bernstein basis
   */
  inline 
  long double m2b_mtrx 
  ( const unsigned short j , 
    const unsigned short k ,
    const unsigned short n ) 
  {
    //
    return
      j < k ? 0.0 : 
      (double) ( Ostap::Math::choose ( j , k ) ) / 
      (double) ( Ostap::Math::choose ( n , k ) ) ; 
  }
  // ==========================================================================
  /// affine transformation of polynomial
  inline 
  long double m2m_mtrx_2
  ( const unsigned short j , 
    const unsigned short k ) 
  {
    if ( k < j ) { return 0 ; }
    const long double c = 
      Ostap::Math::choose ( k , j ) * Ostap::Math::pow ( 2 , j ) ;
    //
    return 0 == ( k - j ) % 2 ?  c : -c ;
  }
  // ==========================================================================
}
// ============================================================================
// constructor from Legendre polynomial
// ============================================================================
Ostap::Math::Bernstein::Bernstein
( const Ostap::Math::LegendreSum& poly )  
  : Ostap::Math::PolySum ( poly.degree () ) 
  , m_xmin ( poly.xmin() ) 
  , m_xmax ( poly.xmax() ) 
{
  for ( unsigned short i = 0 ; i < npars() ; ++i ) 
  { 
    for ( unsigned short k = 0 ; k < npars() ; ++k ) 
    { 
      const double p = poly.par ( k ) ;
      if ( s_zero ( p ) ) { continue ; }
      m_pars[i] += l2b_mtrx ( i , k , degree() ) * p ; 
    } 
  }
}
// ============================================================================
// constructor from Chebyshev polynomial
// ============================================================================
Ostap::Math::Bernstein::Bernstein
( const Ostap::Math::ChebyshevSum& poly )  
  : Ostap::Math::PolySum ( poly.degree () ) 
  , m_xmin ( poly.xmin() ) 
  , m_xmax ( poly.xmax() ) 
{
  //
  for ( unsigned short i = 0 ; i < npars() ; ++i ) 
  { 
    for ( unsigned short k = 0 ; k < npars() ; ++k ) 
    { 
      const double p = poly.par ( k ) ;
      if ( s_zero ( p ) ) { continue ; }
      m_pars[i] += c2b_mtrx ( i , k , degree() ) * p ; 
    } 
  }
  //
}
// ============================================================================
// constructor from simple monomial form 
// ============================================================================
Ostap::Math::Bernstein::Bernstein
( const Ostap::Math::Polynomial& poly )  
  : Ostap::Math::PolySum ( poly.degree () ) 
  , m_xmin ( poly.xmin() ) 
  , m_xmax ( poly.xmax() ) 
{
  //
  const unsigned short np = npars() ;
  // 2-step transformation
  //
  // 1: affine transform to [0,1]
  //
  std::vector<double> _pars ( np ) ;
  for ( unsigned short i = 0 ; i < np ; ++i ) 
  { 
    for ( unsigned short k = i ; k < np ; ++k )  
    { 
      const double p = poly.par ( k )  ;
      if ( s_zero ( p ) ) { continue ; }
      _pars[i] += m2m_mtrx_2 ( i , k ) * p ; 
    } 
  }  
  //
  // 2: tramsform from shifted poly basis:
  //
  for ( unsigned short i = 0 ; i < np ; ++i ) 
  { 
    for ( unsigned short k = 0 ; k <= i ;  ++k )  // ATTENTION!!
    { 
      const double p = _pars[ k ] ;
      if ( s_zero ( p ) ) { continue ; }
      m_pars[i] += m2b_mtrx ( i , k , degree() ) * p ; 
    } 
  }
  //
}
// ============================================================================

// ============================================================================
/* get the integral between 0 and 1 for a product of basic  Bernstein
 *  polynom and the exponential function with the exponent tau
 *  \f[ \int_{0}^{1} \mathcal{B}_{ik} e^{\tau x } \mathrm{d}x \f] 
 *  @param poly  bernstein polynomial
 *  @param tau   slope parameter for an exponential function
 */
// ============================================================================
double Ostap::Math::integrate 
( const Ostap::Math::Bernstein::Basic& b    ,
  const double                         tau  ) 
{
  //
  if ( b.k() > b.N()  ) { return 0                   ; }
  if ( s_zero ( tau ) ) { return 1.0 / ( b.N() + 1 ) ; }
  //
  // make use Kummer function as default scenario 
  return Ostap::Math::kummer ( b.k () + 1 , b.N () + 2 , tau ) / ( b.N() + 1 ) ;
  //
}
// ============================================================================
/* get the integral between \f$x_{min}\f$ and \f$x_{max}\f$ for a product of Bernstein
 *  polynom and the exponential function with the exponent tau
 *  \f[  \int_{x_{min}}^{x_{max}} \mathcal{B} e^{\tau x } \mathrm{d}x \f] 
 *  @param poly  bernstein polynomial
 *  @param tau   slope parameter for exponential 
 */
// ============================================================================
double Ostap::Math::integrate 
( const Ostap::Math::Bernstein& poly ,
  const double                  tau  ) 
{
  if ( s_zero ( tau ) ) { return poly.integral() ; }
  //
  const long double xmin = poly.xmin () ;
  const long double xmax = poly.xmax () ;
  //
  const long double _tau =            ( xmax - xmin ) * tau ;
  const long double _fac = std::exp   (  tau * xmin ) ;
  //
  long double result = 0 ;
  const unsigned short       N    = poly.degree () ;
  const std::vector<double>& pars = poly.pars   () ;
  //
  for ( std::vector<double>::const_iterator ip = pars.begin() ; pars.end() != ip ; ++ip ) 
  {
    if ( s_zero ( *ip ) ) { continue ; } // skip zeroes 
    const unsigned short k =  ip - pars.begin() ;
    const long double    p = *ip ;
    result +=  p * integrate ( Ostap::Math::Bernstein::Basic ( k , N ) , _tau ) ;
  }
  //
  return result * ( xmax - xmin ) * _fac ;
}
// ============================================================================
namespace 
{
  // ==========================================================================
  inline long double r_kNm  
  ( const unsigned short k , 
    const unsigned short N , 
    const unsigned short m ) 
  {
    long double r = 1.0L ;
    for ( unsigned short i = 1 ; i <= m ; ++i ) 
    {
      r *= ( k + i ) ;
      r /= ( N + i ) ;
      r /= i         ;
    }    
    return r ;
  }
  // ==========================================================================
}
// ============================================================================
/* get the integral between 0 and 1 for a product of basic  Bernstein
 *  polynom and monomial or degree m 
 *  \f[  \int_{0}^{1} \mathcal{B} \frac{x^m}{m!} \mathrm{d}x \f] 
 *  @param b     basic bernstein polynomial
 *  @param m     degree of monomial 
 */
// ============================================================================ 
double Ostap::Math::integrate_poly 
( const Ostap::Math::Bernstein::Basic& b ,
  const unsigned short                 m )
{
  //
  const unsigned short N = b.N () ;
  const unsigned short k = b.k () ;
  //
  return r_kNm ( k , N , m ) / ( N + m + 1 ) ;
}
// ============================================================================ 
/*  get the integral between xmin and xmax Bernstein
 *  polynom and monomial or degree m 
 *  \f[  \int_{x_min}^{x_max} \mathcal{B} \frac{x^m}{m!} \mathrm{d}x \f] 
 *  @param b     basic bernstein polynomial
 *  @param m     degree of monomial 
 */
// ============================================================================ 
double Ostap::Math::integrate_poly 
( const Ostap::Math::Bernstein& b ,
  const unsigned short          m ) 
{
  //
  if ( 0 == m ) { return b.integral () ; }
  //
  const std::vector<double>& pars = b.pars()   ;
  const unsigned short       N    = b.degree() ;
  std::vector<long double>   nc  ( pars.size() , 0.0L ) ;
  for ( unsigned short k = 0 ; k < nc.size() ; ++k ) 
  {
    const long double ci = pars[k] ;
    if ( s_zero ( ci ) ) { continue ; }
    nc[k] = r_kNm ( k , N , m ) * ci ;  
  }
  //
  return 
    Ostap::Math::pow ( b.xmax() - b.xmin()  , m + 1 ) * 
    std::accumulate  ( nc.begin() , nc.end() , 0.0L ) / ( N + m + 1 ) ;
}
// ============================================================================ 
namespace 
{
  // ==========================================================================
  long double _integrate_poly_ 
  ( const Ostap::Math::Bernstein& b    ,
    const unsigned short          m    , 
    const double                  low  , 
    const double                  high )
  {
    const std::vector<double>& pars = b.pars()   ;
    const unsigned short       N    = b.degree() ;
    std::vector<long double>   nc  ( pars.size() + m , 0.0L ) ;
    for ( unsigned short k = 0 ; k < pars.size() ; ++k ) 
    {
      const long double ci = pars[k] ;
      if ( s_zero ( ci ) ) { continue ; }
      nc[ k + m ] = r_kNm ( k , N , m ) * ci ;  
    }
    //
    Ostap::Math::Bernstein a ( nc.begin() , nc.end  () , b.xmin() , b.xmax() ) ;
    //
    return Ostap::Math::pow ( b.xmax() - b.xmin()  , m ) * a.integral ( low , high ) ;  
  }
  // ==========================================================================
}
// ============================================================================ 
/* get the integral between xmin and xmax Bernstein
 *  polynom and monomial or degree m 
 *  \f[  \int_{low}^{high} \mathcal{B} \frac{(x-x_min)^m}{m!} \mathrm{d}x \f] 
 *  @param b     basic bernstein polynomial
 *  @param m     degree of monomial 
 *  @param low   low  integration limit 
 *  @param high  high integtation limit 
 */
// ============================================================================ 
double Ostap::Math::integrate_poly 
( const Ostap::Math::Bernstein& b    ,
  const unsigned short          m    , 
  const double                  low  , 
  const double                  high )
{
  //
  if      ( s_equal ( low , high )      ) { return  0 ; }
  else if ( 0 == m                      ) { return  b.integral ( low ,high ) ; }
  else if ( low  > high                 ) { return -integrate_poly ( b , m , high , low ) ; }
  else if ( high < b.xmin ()            ) { return  0 ; }
  else if ( low  > b.xmax ()            ) { return  0 ; } 
  else if ( low  < b.xmin ()            ) { return  integrate_poly ( b , m , b.xmin() , high     ) ; }
  else if ( high > b.xmax ()            ) { return  integrate_poly ( b , m , low      , b.xmax() ) ; }
  else if ( s_equal ( low  , b.xmin() ) && 
            s_equal ( high , b.xmax() ) ) { return  integrate_poly ( b , m ) ; }
  //
  // make the actual integration
  return _integrate_poly_ ( b , m , low , high ) ;
}
// ============================================================================
/*  get the integral between low and high for a product of Bernstein
 *  polynom and the exponential function with the exponent tau
 *  \f[  \int_{a}^{b} \mathcal{B} e^{\tau x } \mathrm{d}x \f] 
 *  @param poly  bernstein polynomial
 *  @param tau   slope parameter for exponential 
 *  @param low   low  integration range 
 *  @param high  high integration range 
 */
// ============================================================================
double Ostap::Math::integrate 
( const Ostap::Math::Bernstein& poly ,
  const double                  tau  ,
  const double                  low  , 
  const double                  high ) 
{
  if      ( s_small ( tau )           ) { return  poly.integral ( low , high ) ; }
  else if ( s_equal ( low , high )    ) { return  0 ; }
  else if ( poly.zero ()              ) { return  0 ; }
  else if ( low  >  high              ) { return -integrate ( poly , tau , high , low ) ; }  
  else if ( high <  poly.xmin () || 
            low  >  poly.xmax ()      ) { return  0 ; }
  else if ( low  <  poly.xmin ()      ) { return  integrate ( poly , tau , poly.xmin() , high        ) ; }
  else if ( high >  poly.xmax ()      ) { return  integrate ( poly , tau , low         , poly.xmax() ) ; }
  //
  if ( s_equal ( low  , poly.xmin() ) && 
       s_equal ( high , poly.xmax() ) ) { return integrate ( poly , tau ) ; }               
  //
  // start series expansion
  // 
  long double result =  poly.integral ( low , high ) ;
  long double dd1    = 1 ;
  long double dd2    = 1 ;
  long double taum   = 1 ;
  //
  const long double xmin = poly.xmin () ;
  // const long double xmax = poly.xmax () ;
  //
  // const long double _tau =            ( xmax - xmin ) * tau ;
  const long double _fac = std::exp   (  tau * xmin ) ;
  //
  for ( unsigned int m = 1 ; m < 10000 ; ++m ) 
  {
    taum   *=  tau ;
    dd2     = _integrate_poly_ ( poly , m , low , high ) * taum  ;
    result += dd2 ;
    if ( s_small ( dd1 / result ) && s_small ( dd2 / result ) ) { break ; }
    dd1     = dd2 ;
  }
  //
  return result * _fac ; 
}
// ============================================================================


// ============================================================================
//  Deflate polynomial 
// ============================================================================
/*  deflate Bernstein polynomial at  <code>x=xmin</code>
 *  \f$ b(x)-b(x_{min})=(x-x_{min})*d(x)\f$      
 *  @param  b  berntein polynomial to be deflated 
 *  @return deflated polinomial "d"
 */ 
// ============================================================================
Ostap::Math::Bernstein 
Ostap::Math::deflate_left  ( const Ostap::Math::Bernstein& b ) 
{
  // trivial case 
  if      ( 1 >  b.degree () ) 
  { return Ostap::Math::Bernstein ( 0 , b.xmin() , b.xmax() ) ; }
  // simple  case
  const std::vector<double>& bpars = b.pars() ;
  std::vector<long double>   dpars ( bpars.size() - 1 ) ;
  //
  const double pz = bpars.front() ;
  dpars.front() =  0 ;
  const unsigned short Nd =  dpars.size() ;
  for ( unsigned short  i = 0 ; i < Nd ; ++i ) 
  {
    const long double p_i = bpars[i+1] - pz  ;
    dpars[i]              = Nd * p_i  / ( i + 1) ;
  }
  // result 
  return Ostap::Math::Bernstein ( dpars.begin () , 
                                  dpars.end   () ,
                                  b.xmin      () ,
                                  b.xmax      () ) ;
}
// ============================================================================
/*  deflate Bernstein polynomial at  <code>x=xmax</code>
 *  \f$ b(x)-b(x_{max})=(x-x_{max})*d(x)\f$      
 *  @param  b  berntein polynomial to be deflated 
 *  @return deflated polinomial "d"
 */ 
// ============================================================================
Ostap::Math::Bernstein 
Ostap::Math::deflate_right ( const Ostap::Math::Bernstein& b ) 
{
  // trivial case 
  if      ( 1 >  b.degree () ) 
  { return Ostap::Math::Bernstein ( 0 , b.xmin() , b.xmax() ) ; }
  //
  // simple  case
  const std::vector<double>& bpars = b.pars() ;
  std::vector<long double>   dpars ( bpars.size() - 1 ) ;
  //
  const double pz = bpars.back () ;
  //
  dpars.back() = 0  ;
  const unsigned short Nd = dpars.size() ;
  for ( unsigned short  i = 0 ; i < Nd ; ++i ) 
  { 
    const long double p_i = bpars[i]  - pz       ;
    dpars[i]              = Nd * p_i / ( Nd - i) ; 
  }
  // result 
  return Ostap::Math::Bernstein ( dpars.begin () , 
                                  dpars.end   () ,
                                  b.xmin      () ,
                                  b.xmax      () ) ;
}
// ============================================================================
/*  deflate Bernstein polynomial at  <code>x=x0</code>
 *  \f$ b(x)-b(x_{0})=(x-x_{0})*d(x)\f$      
 *  @param  b  berntein polynomial to be deflated 
 *  @param  x0 the delfation point 
 *  @return deflated polinomial "d"
 */ 
// ============================================================================
Ostap::Math::Bernstein 
Ostap::Math::deflate ( const Ostap::Math::Bernstein& b , 
                       const double                  x )
{
  //
  // trivial case 
  if      ( 1 >  b.degree () ) 
  { return Ostap::Math::Bernstein ( 0 , b.xmin() , b.xmax() ) ; }
  //
  if      ( s_equal ( x , b.xmin() ) ) { return deflate_left  ( b ) ; }
  else if ( s_equal ( x , b.xmax() ) ) { return deflate_right ( b ) ; } 
  // 
  const long double v   = b.evaluate ( x ) ; 
  const long double tt  =        b.t ( x ) ;
  //
  const bool   reversed = tt <= 0.5 ;
  const double tau      = reversed  ? 1 - tt : tt ;
  //
  const long double pz  = v ;
  //
  const  std::vector<double>& bpars = b.pars()        ;
  const unsigned short        Nd    = bpars.size() -1 ;
  std::vector<long double>    dpars = reversed   ? 
    std::vector<long double> ( bpars.rbegin() , bpars.rbegin() + Nd  ) :
    std::vector<long double> ( bpars. begin() , bpars. begin() + Nd  ) ;
  //
  Ostap::Math::shift ( dpars.begin() , dpars.end() , -pz ) ;
  const long double     u = ( 1 - tau ) / tau  ;
  for ( unsigned short  i = 1 ; i < Nd ; ++i ) 
  {
    const long double p_i = dpars[i] ;
    dpars[i] = ( Nd * p_i + i * u * dpars[i-1] ) / ( Nd - i ) ;
  }
  //
  if ( reversed ) { std::reverse ( dpars.begin() , dpars.end() ) ; }
  //
  return Ostap::Math::Bernstein ( dpars.begin () , 
                                  dpars.end   () , 
                                  b.xmin      () , 
                                  b.xmax      () ) ;
}
// ============================================================================
/*  get abscissas of crosssing points of the control polygon 
 *  for Bernstein polynomial
 *  @param  b bernstein polynomial
 *  @reutrn abscissas of crossing points of the control  polygon
 */
// ============================================================================
std::vector<double> 
Ostap::Math::crossing_points ( const Ostap::Math::Bernstein& b ) 
{
  // trivial case 
  if (  1 > b.degree() ) 
  {
    if ( !s_zero ( b.pars().front() ) ) { return std::vector<double> (              ) ; }
    else                                { return std::vector<double> ( 1 , b.xmin() ) ; }
  }
  //
  const double               norm  = b.norm  () ;
  const std::vector<double>& bpars = b.pars  () ;
  const unsigned  short      N     = b.npars () ;
  //
  std::vector<double> cps ; cps.reserve ( b.degree() + 1 ) ;
  //
  const double p0 = bpars[0] ;
  if ( s_zero( p0 ) || s_equal ( p0 + norm , norm ) ) { cps.push_back ( b.xmin() ) ; }
  //
  for ( unsigned short j = 1 ; j < N ; ++j ) 
  {
    const double pj = bpars[j  ] ;
    const double pi = bpars[j-1] ;
    //
    const double xj =  b.x( float(j) / ( N - 1 ) ) ;
    if ( s_zero ( pj )|| s_equal ( pj + norm , norm ) ) 
    { cps.push_back ( xj ) ; continue ; }
    //
    if ( s_zero ( pi ) || s_equal ( pi + norm , norm ) ) { continue ; }
    //
    const signed char sj = Ostap::Math::signum ( pj ) ;
    const signed char si = Ostap::Math::signum ( pi ) ;
    //
    if ( 0 > si * sj )  // there is root here! 
    {
      const double xi =  b.x( float(j-1) / ( N - 1 ) ) ;
      const double cp = ( xj * pi - xi * pj ) / ( pi - pj ) ;
      cps.push_back ( cp ) ;
    }
    //
  }
  //
  return cps ;
}
// ============================================================================
/*  get number of (strickt) sign changes in trhe sequnce of coefficients
 *  for Bernstein polynomial 
 *  if  N is number of sign changes, then the number of real roots R is 
 *  \f$ R = N - 2K\f$, where K is non-negative integer
 */
// ============================================================================
unsigned short 
Ostap::Math::sign_changes ( const Ostap::Math::Bernstein& b ) 
{
  const std::vector<double>&      bpars = b.pars();
  const Ostap::Math::Tiny<double> s_tiny { b.norm () } ;
  return Ostap::Math::sign_changes ( bpars.begin() , bpars.end() , s_tiny ) ;
}
// ============================================================================
/*  get the most left crossing  point of convex hull with  x-axis 
 *  (it is a step  towards finding the most left root, if any 
 *  if convex hull does not cross the x-axis, xmax is returned      
 */
// ============================================================================
double Ostap::Math::left_line_hull ( const Ostap::Math::Bernstein& b  ) 
{
  const double        bn = b.norm  () ;
  //
  const std::vector<double>& bpars = b.pars()      ;
  const double               p0    = bpars.front() ;
  //
  // left point is already zero 
  if ( s_zero ( p0 ) || s_equal ( p0 + bn , bn ) ) { return b.xmin() ; }
  //
  const signed  char s0 = Ostap::Math::signum ( p0 ) ;
  const bool         up = 0 > p0 ;
  //
  const unsigned short N  = b.npars () ;
  //
  // find the first element with the opposite sign
  unsigned short i = 1 ;
  for ( ;  i < N ; ++i ) 
  { 
    const double pi = bpars[i] ;
    if ( s_zero  ( pi ) || s_equal ( pi +  bn , bn ) ||  
         0 >= s0 * Ostap::Math::signum ( bpars [i] ) ) { break ; }
  }
  //
  // no  good points are found, 
  if ( i == N ) { return b.xmax() + 10 * ( b.xmax() - b.xmin() ) ; } // RETURN
  //
  double         si =  ( bpars[i] - p0 ) / i ;
  for (  unsigned short j = i + 1 ;  j < N ;  ++j ) 
  {
    const double sj = ( bpars[j] - p0 ) / j ;
    if ( ( up && sj >= si ) || ( !up && sj <= si ) )
    {
      i  = j ;
      si = sj ;  
    }
  }
  //
  const double xi = double(i) /  ( N - 1 );
  const double yi = bpars[i] ;
  //
  return b.x ( - xi * p0 / ( yi - p0 ) ) ;
}
// ============================================================================
/*  get the most right rossing  point of convex hull with  x-axis 
 *  (it is a step  towards finding the most right root, if any 
 *  if convex hull does not cross the x-axis, xmin is returned      
 */
// ============================================================================
double Ostap::Math::right_line_hull ( const Ostap::Math::Bernstein& b ) 
{  
  const double        bn = b.norm  () ;
  //
  const std::vector<double>& bpars = b.pars()      ;
  const double               p0    = bpars.back()  ; //  ATTENTION!
  //
  // right point is already zero 
  if ( s_zero ( p0 ) || s_equal ( p0 + bn , bn ) ) { return b.xmax () ; }
  //
  const signed  char s0 = Ostap::Math::signum ( p0 ) ;
  const bool         up = 0 > p0 ;
  //
  const unsigned short N  = b.npars () ;
  //
  // find the first element with the opposite sign
  unsigned short i = 0 ;
  for ( ;  i < N - 1  ; ++i ) 
  {
    const double pi = bpars[i] ;
    if ( s_zero  ( pi ) || s_equal ( pi +  bn , bn ) ||  
         0 >= s0 * Ostap::Math::signum ( bpars [i] ) ) { break ; }
  }
  //
  // no  good points are found, 
  if ( i == N - 1 ) 
  { return b.xmin() - 10 * ( b.xmax() - b.xmin() ) ; } // RETURN
  //
  double         si =  ( bpars[i] - p0 ) / ( N - i ) ;
  for (  unsigned short j = i + 1 ;  j < N ;  ++j ) 
  {
    const double sj = ( bpars[j] - p0 ) /  ( N - j  )  ;
    if ( ( up && sj >= si ) || ( !up && sj <= si ) )
    {
      i  = j ;
      si = sj ;  
    }
  }
  //
  const double xi = double(i) /  ( N - 1 );
  const double yi = bpars[i] ;
  //
  return b.x ( ( yi - xi * p0 ) / ( yi - p0 ) ) ;
}  
// ============================================================================
//  DUAL BASIC 
// ============================================================================
// destructor 
// ============================================================================
Ostap::Math::BernsteinDualBasis::~BernsteinDualBasis() {}
// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::BernsteinDualBasis::BernsteinDualBasis
( const unsigned short N ,
  const unsigned short j ) 
  : std::unary_function<double,double>()
  , m_k         ( j ) 
  , m_bernstein ( N ) 
{
  if ( j <= N ) 
  {
    const unsigned short n = N ;
    for ( unsigned short k = 0 ; k <= N ; ++k ) 
    {
      double ck = 0.0 ;
      const unsigned short imax = std::min ( j ,  k ) ;
      for ( unsigned short i = 0 ; i <= imax ; ++i ) 
      {
        long double a = 2 * i + 1 ;
        a *= c_nk (  n + i + 1 , n - j ) ;
        a *= c_nk (  n - i     , n - j ) ;
        a *= c_nk (  n + i + 1 , n - k ) ;
        a *= c_nk (  n - i     , n - k ) ;
        //
        ck += a ;
      }
      ck /= ( c_nk ( n , j ) * c_nk ( n , k ) ) ;
      if ( ( j + k )  % 2 ) { ck = -ck ; }
      m_bernstein.setPar ( k , ck ) ;
    }
  } 
}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Math::BernsteinDualBasis::BernsteinDualBasis
( const Ostap::Math::BernsteinDualBasis&  right ) 
  : std::unary_function<double,double>( right )
  , m_k         ( right.m_k         ) 
  , m_bernstein ( right.m_bernstein ) 
{}
// ============================================================================
// move constructor 
// ============================================================================
Ostap::Math::BernsteinDualBasis::BernsteinDualBasis
( Ostap::Math::BernsteinDualBasis&& right ) 
  : std::unary_function<double,double>( right )
  , m_k         (             right.m_k           ) 
  , m_bernstein ( std::move ( right.m_bernstein ) ) 
{}

// ============================================================================
// EVEN 
// ============================================================================
/*  constructor
 *  the actual degree of polynomial will be 2*N
 *  @param N  parameter that defiend the order of polynomial (2*N)
 *  @param xmin low edge 
 *  @param xmax high edge 
 */
// ============================================================================
Ostap::Math::BernsteinEven::BernsteinEven 
( const unsigned short N    , 
  const double         xmin ,
  const double         xmax ) 
  : std::unary_function<double,double> ()
  , m_N         (   N                   )   
  , m_bernstein ( 2*N + 1 , xmin , xmax )
{}
// ============================================================================
/*  constructor from list of coefficients 
 *  @param xmin low edge 
 *  @param xmax high edge 
 */
// ============================================================================
Ostap::Math::BernsteinEven::BernsteinEven 
( const std::vector<double>& pars , 
  const double               xmin ,
  const double               xmax ) 
  : std::unary_function<double,double> ()
  , m_N         (   pars.size()                   )   
  , m_bernstein ( 2*pars.size() + 1 , xmin , xmax )
{
  for ( unsigned short i = 0 ; i < pars.size() ; ++i ) { setPar ( i , pars[i] ) ; }
}
// ============================================================================
/* set k-parameter
 *  @param k index
 *  @param value new value 
 *  @return true if parameter is actually changed 
 */
// ============================================================================
bool Ostap::Math::BernsteinEven::setPar
( const unsigned short k , const double value ) 
{
  if ( npars() <= k ) { return false ; }
  const bool b1 = m_bernstein.setPar (             k , value ) ;
  const bool b2 = m_bernstein.setPar ( 2*m_N + 1 - k , value ) ;
  return b1 || b2 ;
}
// ============================================================================
// get all parameters (by value!!! COPY!!)
// ============================================================================
std::vector<double>
Ostap::Math::BernsteinEven::pars () const 
{
  return std::vector<double>( m_bernstein.pars().begin()           , 
                              m_bernstein.pars().begin() + m_N + 1 ) ;                            
}
// ============================================================================
//  Sum of Bernstein polynomial and a constant 
// ============================================================================
Ostap::Math::BernsteinEven
Ostap::Math::BernsteinEven::__add__   ( const double value ) const 
{ BernsteinEven tmp(*this) ; tmp += value ; return tmp ; }
// ============================================================================
//  Sum of Bernstein polynomial and a constant 
// ============================================================================
Ostap::Math::BernsteinEven
Ostap::Math::BernsteinEven::__radd__  ( const double value ) const 
{ BernsteinEven tmp(*this) ; tmp += value ; return tmp ; }
// ============================================================================
//  Product of Bernstein polynomial and a constant 
// ============================================================================
Ostap::Math::BernsteinEven
Ostap::Math::BernsteinEven::__mul__   ( const double value ) const 
{ BernsteinEven tmp(*this) ; tmp *= value ; return tmp ; }
// ============================================================================
//  Product of Bernstein polynomial and a constant 
// ============================================================================
Ostap::Math::BernsteinEven
Ostap::Math::BernsteinEven::__rmul__  ( const double value ) const 
{ BernsteinEven tmp(*this) ; tmp *= value ; return tmp ; }
// ============================================================================
//  Subtract a constant from Bernstein polynomial
// ============================================================================
Ostap::Math::BernsteinEven
Ostap::Math::BernsteinEven::__sub__   ( const double value ) const 
{ BernsteinEven tmp(*this) ; tmp -= value ; return tmp ; }
// ============================================================================
//  Subtract (right) a constant and Bernstein polynomial
// ============================================================================
Ostap::Math::BernsteinEven
Ostap::Math::BernsteinEven::__rsub__  ( const double value ) const 
{ BernsteinEven tmp(*this) ; tmp *= -1 ; tmp += value ; return tmp ; }
// ============================================================================
//  Division of Bernstein polynomial and constant 
// ============================================================================
Ostap::Math::BernsteinEven
Ostap::Math::BernsteinEven::__div__   ( const double value ) const 
{ BernsteinEven tmp(*this) ; tmp /= value ; return tmp ; }
// ============================================================================


    


// ============================================================================
// POSITIVE 
// ============================================================================


// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::Positive::Positive
( const unsigned short      N    ,
  const double              xmin ,
  const double              xmax )
  : std::unary_function<double,double> ()
  , m_bernstein ( N , xmin , xmax )
  , m_sphere    ( N , 3 ) 
{
  updateBernstein () ;
}
// ============================================================================
// constructor from the list of phases 
// ============================================================================
Ostap::Math::Positive::Positive
( const std::vector<double>& pars ,
  const double               xmin ,
  const double               xmax )
  : std::unary_function<double,double> ()
  , m_bernstein ( pars.size() , xmin , xmax )
  , m_sphere    ( pars , 3 ) 
{
  updateBernstein () ;
}
// ============================================================================
// constructor from the sphere with coefficients  
// ============================================================================
Ostap::Math::Positive::Positive
( const Ostap::Math::NSphere& sphere , 
  const double                xmin   , 
  const double                xmax   )
  : std::unary_function<double,double> ()
  , m_bernstein ( sphere.dim() , xmin , xmax )
  , m_sphere    ( sphere ) 
{
  updateBernstein () ;
}
// ============================================================================
// copy 
// ============================================================================
Ostap::Math::Positive::Positive
( const Ostap::Math::Positive&  right ) 
  : std::unary_function<double,double> ( right )
  , m_bernstein ( right.m_bernstein ) 
  , m_sphere    ( right.m_sphere    ) 
{}
// ============================================================================
// move 
// ============================================================================
Ostap::Math::Positive::Positive
(       Ostap::Math::Positive&& right ) 
  : std::unary_function<double,double> ( right )
  , m_bernstein ( std::move ( right.m_bernstein ) ) 
  , m_sphere    ( std::move ( right.m_sphere    ) ) 
{}
// ============================================================================
Ostap::Math::Positive::~Positive() {}
// ============================================================================
// set k-parameter
// ============================================================================
bool Ostap::Math::Positive::setPar ( const unsigned short k , const double value )
{
  //
  const bool update = m_sphere.setPhase ( k , value ) ;
  if ( !update ) { return false ; }   // no actual change 
  //
  return updateBernstein () ;
}
// =============================================================================
// update bernstein coefficients
// =============================================================================
bool Ostap::Math::Positive::updateBernstein ()
{
  ///
  bool         update = false ;
  /// degree 
  const unsigned short o = degree() ;
  //
  const double   norm    = m_bernstein.npars() / 
    ( m_bernstein.xmax() -  m_bernstein.xmin () ) ;
  //
  // few simple cases 
  //
  if       ( 0 == o ) { return m_bernstein.setPar( 0 , norm ) ; }
  else if  ( 1 == o )  
  {
    const bool updated0 = m_bernstein.setPar ( 0 , m_sphere.x2(0) * norm ) ;
    update              = updated0 || update ;
    const bool updated1 = m_bernstein.setPar ( 1 , m_sphere.x2(1) * norm ) ;
    update              = updated1 || update ;
    //
    return update ;
  }
  //
  // get the parameters of "global" parabola 
  //
  const double a0     = m_sphere.x2 ( 0 ) ;
  const double a1_    = m_sphere.x2 ( 1 ) ;
  const double a2     = m_sphere.x2 ( 2 ) ;
  //
  const double a1_min = - std::sqrt ( a0 * a2 ) ; //
  const double a1     = a1_min + a1_ ;            // positivity constraint 
  //
  // simple parabola (probably the most common case in practice) 
  //
  if ( 2 == o ) 
  {
    const double norm2  = norm / ( a0 + a1 + a2 ) ;
    //
    const bool updated0 = m_bernstein.setPar ( 0 , a0 * norm2 ) ;
    update              = updated0 || update ;
    const bool updated1 = m_bernstein.setPar ( 1 , a1 * norm2 ) ;
    update              = updated1 || update ;
    const bool updated2 = m_bernstein.setPar ( 2 , a2 * norm2 ) ;
    update              = updated2 || update ;
    //
    return update ;
  }
  //
  // generic case 
  //
  // get the coefficients from the sphere 
  // this actually represent the positive polynomial with 
  //   - f  (0)=0 
  //   - f' (0)=0  
  //   - f''(0)=0 
  std::vector<double> v ( m_sphere.nX() ) ;
  const unsigned short vs = v.size() ;
  for ( unsigned short ix = 3 ; ix < vs ; ++ix ) { v[ix] = m_sphere.x2 ( ix ) ; }
  //
  const double c0 = a0         ;
  const double c1 = 2*(a1-a0)  ;
  const double c2 = a0+a2-2*a1 ; 
  //
  for ( unsigned short k = 0 ; k < vs ; ++k ) 
  {
    double vv = c0 ;
    const double r1 =  double(k) / o ;
    if ( 0 != k ) { vv += r1             * c1             ; }
    if ( 1 <  k ) { vv += r1 * ( k - 1 ) * c2 / ( o - 1 ) ; }
    v[k] +=  vv ;
    if ( 0 != v[k] && s_zero ( v[k] ) ) {  v[k] = 0 ; }
  }
  //
  const double isum = norm / std::accumulate ( v.begin() , v.end() , 0.0 ) ;
  //
  for ( unsigned short ix = 0 ; ix < m_sphere.nX() ; ++ix ) 
  {
    const bool updated = m_bernstein.setPar ( ix , v[ix] * isum ) ;
    update = updated || update ;
  }
  //
  return update ;
}
// ============================================================================
// copy assignement 
// ============================================================================
Ostap::Math::Positive&
Ostap::Math::Positive::operator=( const Ostap::Math::Positive&  right ) 
{
  if ( &right == this ) { return *this ; }
  m_bernstein = right.m_bernstein ;
  m_sphere    = right.m_sphere    ;
  return *this ;
}
// ============================================================================
// move assignement 
// ============================================================================
Ostap::Math::Positive&
Ostap::Math::Positive::operator=(      Ostap::Math::Positive&& right ) 
{
  if ( &right == this ) { return *this ; }
  m_bernstein = std::move ( right.m_bernstein ) ;
  m_sphere    = std::move ( right.m_sphere    ) ;
  return *this ;
}
// =============================================================================
// get the integral between xmin and xmax
// =============================================================================
double Ostap::Math::Positive::integral () const { return 1 ; } 
// =============================================================================
// get the integral between low and high 
// =============================================================================
double Ostap::Math::Positive::integral
( const double low , const double high ) const 
{ 
  return 
    s_equal ( low  , xmin() ) && s_equal ( high , xmax() ) ? 1 :
    m_bernstein.integral ( low , high )  ; 
}






// ============================================================================
// POSITIVE EVEN
// ============================================================================


// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::PositiveEven::PositiveEven
( const unsigned short      N    ,
  const double              xmin ,
  const double              xmax )
  : std::unary_function<double,double> ()
  , m_even   ( N , xmin , xmax )
  , m_sphere ( N , 3 ) 
{
  updateBernstein () ;
}
// ============================================================================
// constructor from the list of phases 
// ============================================================================
Ostap::Math::PositiveEven::PositiveEven
( const std::vector<double>& pars ,
  const double               xmin ,
  const double               xmax )
  : std::unary_function<double,double> ()
  , m_even      ( pars.size() , xmin , xmax )
  , m_sphere    ( pars , 3 ) 
{
  updateBernstein () ;
}
// ============================================================================
// constructor from the sphere with coefficients  
// ============================================================================
Ostap::Math::PositiveEven::PositiveEven
( const Ostap::Math::NSphere& sphere , 
  const double                xmin   , 
  const double                xmax   )
  : std::unary_function<double,double> ()
  , m_even    ( sphere.dim() , xmin , xmax )
  , m_sphere  ( sphere ) 
{
  updateBernstein () ;
}
// ============================================================================
// set k-parameter
// ============================================================================
bool Ostap::Math::PositiveEven::setPar 
( const unsigned short k , const double value )
{
  //
  const bool update = m_sphere.setPhase ( k , value ) ;
  if ( !update ) { return false ; }   // no actual change 
  //
  return updateBernstein () ;
}
// =============================================================================
// update bernstein coefficients
// =============================================================================
bool Ostap::Math::PositiveEven::updateBernstein ()
{
  ///
  bool         update = false ;
  /// degree 
  const unsigned short o = m_even.degree() ;
  //
  const double   norm    = m_even.npars() / 
    ( m_even.xmax() -  m_even.xmin () ) ;
  //
  // few simple cases 
  //
  if       ( 0 == o ) { return m_even.setPar( 0 , norm ) ; }
  //
  // get the parameters of "global" non-negative symmetric parabola 
  //
  const double a0 = m_sphere.x2(0)      ;
  const double a1 = m_sphere.x2(1) - a0 ;
  const double a2 = a0                  ;
  //
  // "elevate to degree of bernstein"
  const unsigned short N =  m_even.bernstein().degree()  ;
  std::vector<long double> v ( N + 1 ) ;
  v[0] = a0 ;
  v[1] = a1 ;
  v[2] = a2 ;
  std::fill ( v.begin() + 3 , v.end() , a2 ) ;
  // repeate the elevation cycles: 
  for ( unsigned short   n = 2  ; n < N ; ++n ) 
  {
    // "current" degree 
    for ( unsigned short k = n ;  1<= k ; --k ) 
    {
      v[k]  = ( n + 1 - k ) * v[k] + k * v[k-1] ;
      v[k] /=   n + 1  ;
    } 
  }
  //
  // now  we have a non-negative symmetric parabola coded.
  //   - add a proper positive polynomial to it.
  const unsigned short nV = v.size() ;
  const unsigned short nX = m_sphere.nX() ;
  for ( unsigned short ix = 2 ; ix < nX ; ++ix ) 
  {
    const double x = m_sphere.x2 ( ix ) ;
    v[      ix - 2 ] += x ; 
    v[ nV - ix + 1 ] += x ; // keep symmetry 
  }
  //
  const double isum = norm / std::accumulate ( v.begin() , v.end() , 0.0L ) ;
  //
  const unsigned short nE = m_even.npars() ;
  //
  for ( unsigned short ix = 0 ; ix < nE ; ++ix ) 
  {
    const bool updated = m_even.setPar ( ix , 2 * v[ix] * isum ) ;
    update = updated || update ;
  }
  //
  return update ;
}
// =============================================================================
// get the integral between xmin and xmax
// =============================================================================
double Ostap::Math::PositiveEven::integral () const { return 1 ; } 
// =============================================================================
// get the integral between low and high 
// =============================================================================
double Ostap::Math::PositiveEven::integral
( const double low , const double high ) const 
{ 
  return 
    s_equal ( low  , xmin() ) && s_equal ( high , xmax() ) ? 1 :
    m_even.integral ( low , high )  ; 
}





// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::Monothonic::Monothonic 
( const unsigned short      N          ,
  const double              xmin       ,
  const double              xmax       , 
  const bool                increasing ) 
  : Ostap::Math::Positive ( N , xmin , xmax )  
  , m_increasing          ( increasing      )  
{
  updateBernstein () ;
}
// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::Monothonic::Monothonic 
( const std::vector<double>& pars       ,
  const double               xmin       ,
  const double               xmax       ,
  const bool                 increasing ) 
  : Ostap::Math::Positive ( pars , xmin , xmax )  
  , m_increasing          ( increasing      )  
{
  updateBernstein () ;
}
// ============================================================================
// constructor from the spline 
// ============================================================================
Ostap::Math::Monothonic::Monothonic 
( const Ostap::Math::Positive& spline   ,
  const bool                 increasing ) 
  : Ostap::Math::Positive ( spline      )  
  , m_increasing          ( increasing  )  
{
  updateBernstein () ;
}
// ============================================================================
// constructor from the spline 
// ============================================================================
Ostap::Math::Monothonic::Monothonic 
( const Ostap::Math::Monothonic& right ) 
  : Ostap::Math::Positive ( right              )  
  , m_increasing          ( right.m_increasing )
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Math::Monothonic::~Monothonic (){}
// ============================================================================
// update bernstein coefficients
// =============================================================================
bool Ostap::Math::Monothonic::updateBernstein ()
{
  //
  bool   update = false ;
  //
  // get sphere coefficients 
  std::vector<double> v ( m_sphere.nX() ) ;
  for ( unsigned short ix = 0 ; ix < m_sphere.nX() ; ++ix ) 
  { v[ix] = m_sphere.x2 ( ix ) * ( ix + 1 ) ; }
  //
  // integrate them and to get new coefficients
  if   ( m_increasing ) { std::partial_sum ( v. begin() , v. end() ,  v. begin() ) ; }
  else                  { std::partial_sum ( v.rbegin() , v.rend() ,  v.rbegin() ) ; }
  //
  const double isum = m_bernstein.npars() 
    / std::accumulate ( v.begin() , v.end() , 0.0 ) 
    / ( m_bernstein.xmax() -  m_bernstein.xmin () ) ;
  //
  for ( unsigned short ix = 0 ; ix < m_sphere.nX() ; ++ix ) 
  { 
    const bool updated = m_bernstein.setPar ( ix , v [ix] * isum ) ; 
    update = updated || update ;
  }
  //
  return update ;
}
// ============================================================================
// get the minimal value of function 
// ============================================================================
double Ostap::Math::Monothonic::fun_min () const
{
  const std::vector<double>& ps = m_bernstein.pars() ;
  return  std::min( ps.front() , ps.back() ) ;
}
// ============================================================================
// get the maximal value of function 
// ============================================================================
double Ostap::Math::Monothonic::fun_max () const
{
  const std::vector<double>& ps = m_bernstein.pars() ;
  return  std::max( ps.front() , ps.back() ) ;
}
// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::Convex::Convex
( const unsigned short      N          ,
  const double              xmin       ,
  const double              xmax       , 
  const bool                increasing ,
  const bool                convex     ) 
  : Ostap::Math::Monothonic ( N , xmin, xmax , increasing ) 
  , m_convex                ( convex )  
{
  updateBernstein () ;
}
// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::Convex::Convex
( const std::vector<double>& pars       ,
  const double               xmin       ,
  const double               xmax       ,
  const bool                 increasing ,
  const bool                 convex     ) 
  : Ostap::Math::Monothonic ( pars  , xmin, xmax , increasing ) 
  , m_convex                ( convex     )  
{
  updateBernstein () ;
}
// ============================================================================
// constructor from the 
// ============================================================================
Ostap::Math::Convex::Convex
( const Ostap::Math::Positive& poly      ,
  const bool                  increasing ,
  const bool                  convex     ) 
  : Ostap::Math::Monothonic ( poly , increasing ) 
  , m_convex                ( convex     )  
{
  updateBernstein () ;
}
// ============================================================================
// constructor from the 
// ============================================================================
Ostap::Math::Convex::Convex
( const Ostap::Math::Monothonic& poly   ,
  const bool                     convex ) 
  : Ostap::Math::Monothonic ( poly       ) 
  , m_convex                ( convex     )  
{
  updateBernstein () ;
}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Math::Convex::Convex ( const Ostap::Math::Convex& right  ) 
  : Ostap::Math::Monothonic ( right           ) 
  , m_convex                ( right.m_convex  )  
{}
// ============================================================================
Ostap::Math::Convex::~Convex (){}
// ============================================================================
// update bernstein coefficients
// =============================================================================
bool Ostap::Math::Convex::updateBernstein ()
{
  //
  bool   update = false ;
  //
  // get sphere coefficients 
  //
  std::vector<double>  v ( m_sphere.nX() ) ;
  const unsigned short vs = v.size()    ;
  //
  const std::array<double,2> a = { { m_sphere.x2(0) , m_sphere.x2(1) } };
  for ( unsigned short ix = 2 ; ix < vs ; ++ix ) 
  { v[ix] = m_sphere.x2 ( ix ) ; }
  //
  // integrate them twice and to get new coefficients
  std::partial_sum ( v.  begin() + 2 , v.  end()     ,  v.  begin() + 2 ) ; 
  std::partial_sum ( v.  begin() + 2 , v.  end()     ,  v.  begin() + 2 ) ; 
  //
  if ( !m_convex ) 
  {
    const  double last = v.back() ;
    for ( unsigned short k = 0 ; k < vs; ++k) 
    { 
      v[k] = last  - v[k] ; 
      if ( 0 != v[k] && s_zero ( v[k] ) ) {  v[k] = 0 ; }
    }
  }
  //
  if ( m_increasing != m_convex )
  { std::reverse ( v.begin() , v.end() ) ; }
  //
  // add a positive linear function 
  //
  const unsigned short d = degree() ;
  for ( unsigned short k = 0 ; k < vs ; ++k ) 
  {
    const double r1 =  double(k) / d ;
    //
    v[k] +=  
      m_increasing ?  
      a[0] +       r1   * a[1] :
      a[0] + ( 1 - r1 ) * a[1] ;
    //
    if ( 0 != v[k] && s_zero ( v[k] ) ) {  v[k] = 0 ; }
  }
  //
  const double isum = m_bernstein.npars() 
    / std::accumulate ( v.begin() , v.end() , 0.0 ) 
    / ( m_bernstein.xmax() -  m_bernstein.xmin () ) ;
  //
  for ( unsigned short ix = 0 ; ix < vs ; ++ix ) 
  { 
    const bool updated = m_bernstein.setPar ( ix , v [ix] * isum ) ; 
    update = updated || update ;
  }
  //
  return update ;  
}
// ============================================================================




// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::ConvexOnly::ConvexOnly 
( const unsigned short      N     ,
  const double              xmin   ,
  const double              xmax   , 
  const bool                convex ) 
  : Ostap::Math::Positive ( N , xmin , xmax )  
  , m_convex          ( convex )  
{
  updateBernstein () ;
}
// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::ConvexOnly::ConvexOnly 
( const std::vector<double>& pars    ,
  const double               xmin    ,
  const double               xmax    ,
  const bool                 convex  ) 
  : Ostap::Math::Positive ( pars , xmin , xmax )  
  , m_convex          ( convex      )  
{
  updateBernstein () ;
}
// ============================================================================
// constructor from the spline 
// ============================================================================
Ostap::Math::ConvexOnly::ConvexOnly 
( const Ostap::Math::Positive& poly   ,
  const bool                   convex ) 
  : Ostap::Math::Positive ( poly   )  
  , m_convex              ( convex )  
{
  updateBernstein () ;
}
// ============================================================================
// constructor from the spline 
// ============================================================================
Ostap::Math::ConvexOnly::ConvexOnly 
( const Ostap::Math::ConvexOnly& right ) 
  : Ostap::Math::Positive ( right )  
  , m_convex ( right.m_convex )  
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Math::ConvexOnly::~ConvexOnly (){}
// ============================================================================
// update bernstein coefficients
// =============================================================================
bool Ostap::Math::ConvexOnly::updateBernstein ()
{
  //
  // linear function...
  if ( 2 > degree() ) { return Ostap::Math::Positive::updateBernstein() ; }
  //
  bool   update = false ;
  //
  // get sphere coefficients 
  //
  std::vector<double>  v ( m_sphere.nX() ) ;
  const unsigned short vs = v.size()    ;
  //
  // get parameters from the sphere:
  //
  if ( !m_convex ) 
  {
    const std::array<double,2> a = { { m_sphere.x2(0) , m_sphere.x2(1) } };
    for ( unsigned short ix = 2 ; ix < vs ; ++ix ) 
    { v[ix] = m_sphere.x2 ( ix ) ; }
    //
    // integrate them twice and to get new coefficients
    std::partial_sum ( v.  begin() + 2 , v.  end()     ,  v.  begin() + 2 ) ; 
    std::partial_sum ( v.  begin() + 2 , v.  end()     ,  v.  begin() + 2 ) ; 
    //
    {
      const  double last = v.back() ;
      for ( unsigned short k = 0 ; k < vs; ++k) 
      { 
        v[k] = last  - v[k] ; 
        if ( 0 != v[k] && s_zero ( v[k] ) ) {  v[k] = 0 ; }
      }
    }
    //
    // subtract the linear component and 
    // add positive linear function
    //
    const double v1 = a[0] - v.front() ;
    const double v2 = a[1] - v.back()  ;
    const unsigned int   d = degree() ;
    for ( unsigned short k = 0 ; k < vs ; ++k ) 
    {
      const double r1 =  double(k)  / d ;
      v[k] +=  ( 1 - r1 ) * v1  + r1 * v2 ;
      if ( 0 != v[k] && s_zero ( v[k] ) ) {  v[k] = 0 ; }
    }
  }
  else 
  { 
    std::array<double,3> a = { { m_sphere.x2(0) , 
                                 m_sphere.x2(1) , 
                                 m_sphere.x2(2) } };
    for ( unsigned short ix = 3 ; ix < vs ; ++ix ) 
    { v[ix] = m_sphere.x2 ( ix ) ; }
    // integrate them twice and to get new coefficients
    std::partial_sum ( v.  begin() + 3 , v.  end()     ,  v.  begin() + 3 ) ; 
    std::partial_sum ( v.  begin() + 3 , v.  end()     ,  v.  begin() + 3 ) ; 
    //    
    const double a0 = a[0] ;
    const double a2 = a[2] ;
    const double a1_min = -1*std::sqrt ( a0 * a2 ) ;
    const double a1_max = 0.5 * ( a0 + a2 ) ;
    //
    const double a1 = a1_min + a[1] * ( a1_max - a1_min ) ;
    //
    const double c0 = a0         ;
    const double c1 = 2*(a1-a0)  ;
    const double c2 = a0+a2-2*a1 ; 
    //
    const unsigned int   d = degree() ;
    for ( unsigned short k = 0 ; k < vs ; ++k ) 
    {
      double vv = c0 ;
      const double r1 =  double(k) / d ;
      if ( 0 != k ) { vv += r1 * c1 ; }
      if ( 1 <  k ) { vv += r1 * ( k - 1 ) * c2 / ( d - 1 ) ; }
      v[k] +=  vv ;
      if ( 0 != v[k] && s_zero ( v[k] ) ) {  v[k] = 0 ; }
    }
  }
  //
  const double isum = m_bernstein.npars() 
    / std::accumulate ( v.begin() , v.end() , 0.0 ) 
    / ( m_bernstein.xmax() -  m_bernstein.xmin () ) ;
  //
  for ( unsigned short ix = 0 ; ix < vs ; ++ix ) 
  { 
    const bool updated = m_bernstein.setPar ( ix , v [ix] * isum ) ; 
    update = updated || update ;
  }
  //
  return update ;
}
// ============================================================================






// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::Bernstein2D::Bernstein2D
( const unsigned short      nX   ,
  const unsigned short      nY   ,
  const double              xmin ,
  const double              xmax ,
  const double              ymin ,
  const double              ymax )
  : std::binary_function<double,double,double> ()
//
  , m_nx   ( nX ) 
  , m_ny   ( nY )
//
  , m_pars ( ( nX + 1 ) * ( nY + 1 ) , 0.0 )
//
  , m_xmin ( std::min ( xmin , xmax ) )
  , m_xmax ( std::max ( xmin , xmax ) )
  , m_ymin ( std::min ( ymin , ymax ) )
  , m_ymax ( std::max ( ymin , ymax ) )
    //
  , m_bx   () 
  , m_by   ()
{
  //
  typedef  Ostap::Math::Bernstein::Basic BB ;
  for ( unsigned short ix = 0 ; ix <= nX ; ++ix ) 
  { m_bx.push_back ( Bernstein ( BB ( ix , nX ) , xmin , xmax ) ) ; }
  //
  for ( unsigned short iy = 0 ; iy <= nY ; ++iy ) 
  { m_by.push_back ( Bernstein ( BB ( iy , nY ) , ymin , ymax ) ) ; }
  //
}
// ============================================================================
// get the value
// ============================================================================
double Ostap::Math::Bernstein2D::operator () ( const double x ,
                                               const double y ) const
{
  /// the trivial cases
  if ( x < m_xmin || x > m_xmax ) { return 0.0        ; }
  if ( y < m_ymin || y > m_ymax ) { return 0.0        ; }
  //
  const double scalex = ( m_nx + 1 ) / ( xmax() - xmin() ) ;
  const double scaley = ( m_ny + 1 ) / ( ymax() - ymin() ) ;
  //
  if      ( 0 == npars ()       ) { return 0.0        ; }
  else if ( 1 == npars ()       ) { return m_pars [0] * scalex * scaley ; }
  ///
  std::vector<double> fy ( m_ny + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= m_ny ; ++i )
  { fy[i] = m_by[i] ( y )  ; }
  //
  std::vector<double> fx ( m_nx + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= m_nx ; ++i ) 
  { fx[i] = m_bx[i] ( x )  ; }
  //
  double       result = 0 ;
  for  ( unsigned short ix = 0 ; ix <= m_nx ; ++ix )
  { 
    for  ( unsigned short iy = 0 ; iy <= m_ny ; ++iy ) 
    { result += par ( ix , iy ) * fx[ix] * fy[iy] ; }
  }
  //
  return result * ( scalex * scaley ) ;
}
// ============================================================================
/** get the integral over 2D-region 
 *  \f[  x_min < x < x_max, y_min< y< y_max\f] 
 */
// ============================================================================
double Ostap::Math::Bernstein2D::integral() const 
{ return std::accumulate ( m_pars.begin() , m_pars.end() , 0.0 ) ; }
// ============================================================================
/* get the integral over 2D-region 
 *  \f[ \int_{x_low}^{x_high}\int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y\f] 
 *  @param xlow  low  edge in x 
 *  @param xhigh high edge in x 
 *  @param ylow  low  edge in y 
 *  @param yhigh high edge in y 
 */
// ============================================================================
double Ostap::Math::Bernstein2D::integral 
( const double xlow , const double xhigh , 
  const double ylow , const double yhigh ) const 
{
  if      ( s_equal ( xlow , xhigh ) || s_equal ( ylow  , yhigh ) ) { return 0 ; }
  //
  else if ( s_equal ( xlow  , m_xmin ) && 
            s_equal ( xhigh , m_xmax ) && 
            s_equal ( ylow  , m_ymin ) && 
            s_equal ( yhigh , m_ymax )  )  { return integral () ; }
  //
  else if ( xlow  > xhigh ) { return -1*integral ( xhigh , xlow  , ylow  , yhigh ) ; }
  else if ( ylow  > yhigh ) { return -1*integral ( xlow  , xhigh , yhigh , ylow  ) ; }
  //
  else if ( xhigh <  xmin () || xlow >  xmax() ) { return 0 ; }
  else if ( yhigh <  ymin () || ylow >  ymax() ) { return 0 ; }
  //
  const double  x_low  = std::max ( xmin() , xlow  ) ;
  const double  x_high = std::min ( xmax() , xhigh ) ;
  if ( x_low >= x_high ) { return 0 ; }
  //
  const double  y_low  = std::max ( ymin() , ylow  ) ;
  const double  y_high = std::min ( ymax() , yhigh ) ;
  if ( y_low >= y_high ) { return 0 ; }
  //
  std::vector<double> fy ( m_ny + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= m_ny ; ++i ) 
  { fy[i] = m_by[i].integral ( y_low , y_high ) ; }
  //
  std::vector<double> fx ( m_nx + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= m_nx ; ++i ) 
  { fx[i] = m_bx[i].integral ( x_low , x_high ) ; }
  //
  double result = 0 ;
  for  ( unsigned short ix = 0 ; ix <= m_nx ; ++ix ) 
  {
    for  ( unsigned short iy = 0 ; iy <= m_ny ; ++iy ) 
    { result += par ( ix , iy ) * fx[ix] * fy[iy] ; }
    //
  }
  //
  const double scalex = ( m_nx + 1 ) / ( xmax() - xmin() ) ;
  const double scaley = ( m_ny + 1 ) / ( ymax() - ymin() ) ;
  //
  return result * ( scalex * scaley ) ;
}
// ============================================================================
/*  integral over x-dimension 
 *  \f[ \int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}y\f] 
 *  @param x     variable 
 *  @param ylow  low  edge in y 
 *  @param yhigh high edge in y 
 */
// ============================================================================
double Ostap::Math::Bernstein2D::integrateX 
( const double y    , 
  const double xlow , const double xhigh ) const 
{
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; }
  else if ( xlow  > xhigh  ) { return -1*integrateX ( y , xhigh , xlow ) ; }
  else if ( xhigh <= xmin () || xlow >= xmax() ) { return 0 ; }
  else if ( y     <  ymin () || y    >  ymax() ) { return 0 ; }
  else if ( s_equal ( xlow  , m_xmin ) && 
            s_equal ( xhigh , m_xmax )         ) { return integrateX ( y ) ; }
  //
  const double  x_low  = std::max ( xmin() , xlow  ) ;
  const double  x_high = std::min ( xmax() , xhigh ) ;
  if ( x_low >= x_high ) { return 0 ; }
  //
  std::vector<double> fy ( m_ny + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= m_ny ; ++i ) 
  { fy[i] = m_by[i] ( y ) ; }
  //
  std::vector<double> fx ( m_nx + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= m_nx ; ++i ) 
  { fx[i] = m_bx[i].integral ( x_low , x_high ) ; }
  //
  double result = 0 ;
  for  ( unsigned short ix = 0 ; ix <= m_nx ; ++ix ) 
  {
    for  ( unsigned short iy = 0 ; iy <= m_ny ; ++iy ) 
    { result += par ( ix , iy ) * fx[ix] * fy[iy] ; }
    //
  }
  //
  const double scalex = ( m_nx + 1 ) / ( xmax() - xmin() ) ;
  const double scaley = ( m_ny + 1 ) / ( ymax() - ymin() ) ;
  //
  return result * ( scalex * scaley ) ;
}
// ============================================================================
/*  integral over x-dimension 
 *  \f[ \int_{x_low}^{x_high} \mathcal{B}(x,y) \mathrm{d}x\f] 
 *  @param y     variable 
 *  @param xlow  low  edge in x 
 *  @param xhigh high edge in x 
 */
// ============================================================================
double Ostap::Math::Bernstein2D::integrateY 
( const double x    , 
  const double ylow , const double yhigh ) const 
{
  if      ( s_equal ( ylow  , yhigh ) ) { return 0 ; }
  else if ( ylow  >  yhigh ) { return -1*integrateY ( x , yhigh , ylow  ) ; }
  else if ( x     <  xmin () || x    >  xmax() ) { return 0 ; }
  else if ( yhigh <= ymin () || ylow >= ymax() ) { return 0 ; }
  else if ( s_equal ( ylow  , m_ymin ) && 
            s_equal ( yhigh , m_ymax )         ) { return integrateY ( x ) ; }
  //
  const double  y_low  = std::max ( ymin() , ylow  ) ;
  const double  y_high = std::min ( ymax() , yhigh ) ;
  if ( y_low >= y_high ) { return 0 ; }
  //
  std::vector<double> fy ( m_ny + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= m_ny ; ++i ) 
  { fy[i] = m_by[i].integral ( y_low , y_high ) ; }
  //
  std::vector<double> fx ( m_nx + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= m_nx ; ++i ) 
  { fx[i] = m_bx[i] ( x ) ; }
  //
  double result = 0 ;
  for  ( unsigned short ix = 0 ; ix <= m_nx ; ++ix ) 
  {
    for  ( unsigned short iy = 0 ; iy <= m_ny ; ++iy ) 
    { result += par ( ix , iy ) * fx[ix] * fy[iy] ; }
    //
  }
  //
  const double scalex = ( m_nx + 1 ) / ( xmax() - xmin() ) ;
  const double scaley = ( m_ny + 1 ) / ( ymax() - ymin() ) ;
  //
  return result * ( scalex * scaley ) ;
}
// ============================================================================
/*  integral over x-dimension 
 *  \f[ \int_{y_{min}}^{y_{max}} \mathcal{B}(x,y) \mathrm{d}y\f] 
 *  @param x     variable 
 */
// ============================================================================
double Ostap::Math::Bernstein2D::integrateX ( const double y ) const 
{
  if ( y < ymin () || y > ymax() ) { return 0 ; }
  //
  std::vector<double> fy ( m_ny + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= m_ny ; ++i ) 
  { fy[i] = m_by[i] ( y ) ; }
  //
  const std::vector<double> fx ( m_nx + 1 , 1 ) ;
  //
  double result = 0 ;
  for  ( unsigned short ix = 0 ; ix <= m_nx ; ++ix ) 
  {
    for  ( unsigned short iy = 0 ; iy <= m_ny ; ++iy ) 
    { result += par ( ix , iy ) * fx[ix] * fy[iy] ; }
    //
  }
  //
  // const double scalex = ( m_nx + 1 ) / ( xmax() - xmin() ) ;
  const double scaley = ( m_ny + 1 ) / ( ymax() - ymin() ) ;
  //
  return result * scaley  ;
}
// ============================================================================
/*  integral over x-dimension 
 *  \f[ \int_{x_{min}}^{x_{max}} \mathcal{B}(x,y) \mathrm{d}x\f] 
 *  @param y     variable 
 */
// ============================================================================
double Ostap::Math::Bernstein2D::integrateY 
( const double x ) const 
{
  if ( x < xmin () || x > xmax() ) { return 0 ; }
  //
  const std::vector<double> fy ( m_ny + 1 , 1.0 ) ;
  //
  std::vector<double> fx ( m_nx + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= m_nx ; ++i ) 
  { fx[i] = m_bx[i] ( x ) ; }
  //
  double result = 0 ;
  for  ( unsigned short ix = 0 ; ix <= m_nx ; ++ix ) 
  {
    for  ( unsigned short iy = 0 ; iy <= m_ny ; ++iy ) 
    { result += par ( ix , iy ) * fx[ix] * fy[iy] ; }
    //
  }
  //
  const double scalex = ( m_nx + 1 ) / ( xmax() - xmin() ) ;
  // const double scaley = ( m_ny + 1 ) / ( ymax() - ymin() ) ;
  //
  return result * scalex ; // * scaley ) ;
}
// ============================================================================
// set (l,m)-parameter
// ============================================================================
bool Ostap::Math::Bernstein2D::setPar
( const unsigned short l     , 
  const unsigned short m     , 
  const double         value )
{
  if ( l > m_nx || m > m_ny )             { return false ; }
  const unsigned int k =  l * ( m_ny + 1 ) + m ;
  return setPar ( k , value ) ;
}
// ============================================================================
// set k-parameter
// ============================================================================
bool Ostap::Math::Bernstein2D::setPar
( const unsigned int   k     , 
  const double         value )
{
  if ( k >= npars() )                     { return false ; }
  if ( s_equal ( m_pars [ k ] , value ) ) { return false ; }
  m_pars [ k ] = value ;
  return true ;
}
// ============================================================================
// get (l,m)-parameter 
// ============================================================================
double  Ostap::Math::Bernstein2D::par 
( const unsigned short l ,
  const unsigned short m ) const 
{
  if ( l > m_nx || m > m_ny ) { return 0 ; }
  const unsigned int k =  l * ( m_ny + 1 ) + m ;
  return par ( k ) ;
}
// ============================================================================

  

// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::Bernstein2DSym::Bernstein2DSym
( const unsigned short      n    ,
  const double              xmin ,
  const double              xmax )
  : std::binary_function<double,double,double> ()
//
  , m_n    ( n ) 
//
  , m_pars ( ( n + 1 ) * ( n + 2 ) / 2 , 0.0 )
//
  , m_xmin ( std::min ( xmin , xmax ) )
  , m_xmax ( std::max ( xmin , xmax ) )
//
  , m_b    () 
{
  //
  typedef  Ostap::Math::Bernstein::Basic BB ;
  for ( unsigned short i = 0 ; i <= n ; ++i ) 
  { m_b.push_back ( Bernstein ( BB ( i , n ) , xmin , xmax ) ) ; }
  //
}
// ============================================================================
// get the value
// ============================================================================
double Ostap::Math::Bernstein2DSym::operator () 
  ( const double x ,
    const double y ) const
{
  /// the trivial cases
  if ( x < xmin () || x > xmax () ) { return 0.0        ; }
  if ( y < ymin () || y > ymax () ) { return 0.0        ; }
  //
  const double scale = ( m_n + 1 ) / ( xmax() - xmin() ) ;
  //
  if      ( 0 == npars ()       ) { return 0.0 ; }
  else if ( 1 == npars ()       ) { return m_pars [0] * ( scale * scale ) ; }
  ///
  std::vector<double> fy ( m_n + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= m_n ; ++i ) 
  { fy[i] = m_b[i] ( y ) ; }
  //
  std::vector<double> fx ( m_n + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= m_n ; ++i ) 
  { fx[i] = m_b[i] ( x ) ; }
  //
  double       result = 0 ;
  for  ( unsigned short ix = 0 ; ix <= m_n ; ++ix ) 
  {
    for  ( unsigned short iy = 0 ; iy <= m_n ; ++iy ) 
    { 
      result += 
        ix == iy  ? par ( ix , iy ) * fx[ix] * fy[iy] : 
        0.5       * par ( ix , iy ) * fx[ix] * fy[iy] ;
    }
  }
  //
  return result * ( scale * scale ) ;
}
// ============================================================================
/* get the integral over 2D-region 
 *  \f[ \int_{x_{low}}^{x_{high}}\int_{y_{low}}^{y_{high}} 
 *  \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y\f] 
 *  @param xlow  low  edge in x 
 *  @param xhigh high edge in x 
 *  @param ylow  low  edge in y 
 *  @param yhigh high edge in y 
 */
// ============================================================================
double Ostap::Math::Bernstein2DSym::integral 
( const double xlow , const double xhigh ,
  const double ylow , const double yhigh ) const 
{
  //
  if      ( xlow  > xhigh   ) { return -integral ( xhigh , xlow    , ylow   , yhigh  ) ; }
  else if ( ylow  > yhigh   ) { return -integral ( xlow  , xhigh   , yhigh  , ylow   ) ; }
  //
  else if ( xlow  < xmin () ) { return  integral ( xmin() , xhigh  , ylow   , yhigh  ) ; }
  else if ( xhigh > xmax () ) { return  integral ( xlow   , xmax() , ylow   , yhigh  ) ; }
  else if ( ylow  < ymin () ) { return  integral ( xlow   , xhigh  , ymin() , yhigh  ) ; }
  else if ( yhigh > ymax () ) { return  integral ( xlow   , xhigh  , ylow   , ymax() ) ; }
  //
  else if ( s_equal ( xlow , xhigh ) || s_equal ( ylow , yhigh ) ) { return 0 ; }
  //
  //
  std::vector<double> fy ( m_n + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= m_n ; ++i ) 
  { fy[i] = m_b[i].integral ( ylow , yhigh ) ; }
  //
  std::vector<double> fx ( m_n + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= m_n ; ++i ) 
  { fx[i] = m_b[i].integral ( xlow , xhigh ) ; }
  //
  double       result = 0 ;
  for  ( unsigned short ix = 0 ; ix <= m_n ; ++ix ) 
  {
    for  ( unsigned short iy = 0 ; iy <= m_n ; ++iy ) 
    {
      result += 
        ix == iy  ? par ( ix , iy ) * fx[ix] * fy[iy] : 
        0.5       * par ( ix , iy ) * fx[ix] * fy[iy] ;
    }
  }
  //
  const double scale = ( m_n + 1 ) / ( xmax() - xmin() ) ;
  return result * ( scale * scale ) ;
}
// ============================================================================
/*  integral over x-dimension 
 *  \f[ \int_{x_{low}}^{x_{high}} \mathcal{B}(x,y) \mathrm{d}x\f] 
 *  @param x     variable 
 *  @param xlow  low  edge in x 
 *  @param xhigh high edge in x 
 */
// ============================================================================
double Ostap::Math::Bernstein2DSym::integrateX 
( const double y    ,  
  const double xlow , const double xhigh ) const 
{ return integrateY ( y , xlow , xhigh ) ; }
// ============================================================================
/*  integral over y-dimension 
 *  \f[ \int_{y_{low}}^{x_{high}} \mathcal{B}(x,y) \mathrm{d}y\f] 
 *  @param y     variable 
 *  @param xlow  low  edge in x 
 *  @param xhigh high edge in x 
 */
// ============================================================================
double Ostap::Math::Bernstein2DSym::integrateY
( const double x    ,
  const double ylow , const double yhigh ) const 
{
  //
  if      ( s_equal ( ylow  , yhigh ) ) { return 0 ; }
  else if ( ylow  > yhigh ) { return -integrateY ( x , yhigh , ylow  ) ; }
  else if ( x     <  xmin () || x    >  xmax() ) { return 0 ; }
  else if ( yhigh <  ymin () || ylow >  ymax() ) { return 0 ; }
  else if ( s_equal ( ylow  , ymin () ) && 
            s_equal ( yhigh , ymax () )  ) { return integrateY ( x ) ; }
  //
  const double  y_low  = std::max ( ymin() , ylow  ) ;
  const double  y_high = std::min ( ymax() , yhigh ) ;
  if ( y_low >= y_high ) { return 0 ; }
  //
  std::vector<double> fy ( m_n + 1 , 0 ) ;
  for ( unsigned short i = 0 ; i <= m_n ; ++i ) 
  { fy[i] = m_b[i].integral ( y_low , y_high ) ; }
  //
  std::vector<double> fx ( m_n + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= m_n ; ++i ) 
  { fx[i] = m_b[i] ( x ) ; }
  //
  double       result = 0 ;
  for  ( unsigned short ix = 0 ; ix <= m_n ; ++ix ) 
  {
    for  ( unsigned short iy = 0 ; iy <= m_n ; ++iy ) 
    {
      result += 
        ix == iy  ? par ( ix , iy ) * fx[ix] * fy[iy] : 
        0.5       * par ( ix , iy ) * fx[ix] * fy[iy] ;
    } 
  }
  //
  const double scale = ( m_n + 1 ) / ( xmax() - xmin() ) ;
  return result * ( scale * scale ) ;
}
// ============================================================================
/* get the integral over 2D-region 
 *  \f[ \int_{x_{min}}^{x_{max}}\int_{y_{min}}^{y_{max}} 
 *  \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y\f] 
 */
// ============================================================================
double  Ostap::Math::Bernstein2DSym::integral   () const 
{ return std::accumulate ( m_pars.begin() , m_pars.end() , 0.0 ) ; }
// ============================================================================
/* integral over x-dimension 
 *  \f[ \int_{x_{min}}^{x_{max}} \mathcal{B}(x,y) \mathrm{d}x\f] 
 *  @param x     variable 
 */
// ============================================================================
double Ostap::Math::Bernstein2DSym::integrateX ( const double y ) const 
{ return integrateY ( y ) ; }
// ============================================================================
/*  integral over y-dimension 
 *  \f[ \int_{y_{min}}^{x_{max}} \mathcal{B}(x,y) \mathrm{d}y\f] 
 *  @param y     variable 
 */
// ============================================================================
double Ostap::Math::Bernstein2DSym::integrateY ( const double x ) const 
{
  //
  if ( x     <  xmin () || x    >  xmax() ) { return 0 ; }
  //
  std::vector<double> fx ( m_n + 1 , 0 ) ;
  for  ( unsigned short i = 0 ; i <= m_n ; ++i ) { fx[i] = m_b[i] ( x ) ; }
  //
  double       result = 0 ;
  for  ( unsigned short ix = 0 ; ix <= m_n ; ++ix ) 
  {
    for  ( unsigned short iy = 0 ; iy <= m_n ; ++iy ) 
    {
      result += 
        ix == iy  ? par ( ix , iy ) * fx[ix] : 
        0.5       * par ( ix , iy ) * fx[ix] ;
    }
  }
  //
  const double scale = ( m_n + 1 ) / ( xmax() - xmin() ) ;
  return result * scale  ;
}
// ============================================================================
// set (k)-parameter
// ============================================================================
bool Ostap::Math::Bernstein2DSym::setPar
( const unsigned int   k     , 
  const double         value )
{
  //
  if ( k >= npars() )                     { return false ; }
  if ( s_equal ( m_pars [ k ] , value ) ) { return false ; }
  m_pars [ k ] = value ;
  //
  return true ;
}
// ============================================================================
// set (l,m)-parameter
// ============================================================================
bool Ostap::Math::Bernstein2DSym::setPar
( const unsigned short l     , 
  const unsigned short m     , 
  const double         value )
{
  //
  if ( l > m_n || m > m_n )               { return false ; }
  //
  const unsigned int k = ( l < m ) ? 
    ( m * ( m + 1 ) / 2 + l ) : 
    ( l * ( l + 1 ) / 2 + m ) ;
  //
  return setPar ( k , value ) ;
}
// ============================================================================
// get (l,m)-parameter 
// ============================================================================
double Ostap::Math::Bernstein2DSym::par
( const unsigned short l ,
  const unsigned short m ) const 
{
  //
  if ( l > m_n || m > m_n )               { return 0 ; }
  //
  const unsigned int k = ( l < m ) ? 
    ( m * ( m + 1 ) / 2 + l ) : 
    ( l * ( l + 1 ) / 2 + m ) ;
  //
  return par ( k ) ;
}
// ============================================================================
// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::Positive2D::Positive2D
( const unsigned short      nX   ,
  const unsigned short      nY   ,
  const double              xmin ,
  const double              xmax ,
  const double              ymin ,
  const double              ymax )
  : std::binary_function<double,double,double> ()
//
  , m_bernstein (   nX , nY , xmin , xmax , ymin , ymax ) 
  , m_sphere    ( ( nX + 1 ) * ( nY + 1 ) - 1 )
{
  updateBernstein () ;
}
// ============================================================================
// set k-parameter
// ============================================================================
bool Ostap::Math::Positive2D::setPar 
( const unsigned int k     , 
  const double       value )
{
  //
  const bool update = m_sphere.setPhase ( k , value ) ;
  if ( !update ) { return false ; }   // no actual change 
  //
  return updateBernstein () ;
}
// =============================================================================
// update bernstein coefficients
// =============================================================================
bool Ostap::Math::Positive2D::updateBernstein ()
{
  //
  bool update = false ;
  for ( unsigned int ix = 0 ; ix < m_sphere.nX() ; ++ix ) 
  { 
    const bool updated = m_bernstein.setPar ( ix , m_sphere.x2 ( ix ) ) ;
    update = updated || update ;  
  }
  //
  return update ;
}
// ============================================================================
// get the parameter value
// ============================================================================
double Ostap::Math::Positive2D::par ( const unsigned int k ) const 
{ return m_sphere.phase ( k ) ; }
// ============================================================================
/*  get the integral over 2D-region           
 *  \f[ \int_{x_{min}}^{x_{max}}\int_{y_{min}}^{y_{max}} 
 *        \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y\f] 
 */
// ============================================================================
double  Ostap::Math::Positive2D::integral   () const { return 1 ; }
// ============================================================================
/* get the integral over 2D-region 
 *  \f[ \int_{x_low}^{x_high}\int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y\f] 
 *  @param xlow  low  edge in x 
 *  @param xhigh high edge in x 
 *  @param ylow  low  edge in y 
 *  @param yhigh high edge in y 
 */
// ============================================================================
double Ostap::Math::Positive2D::integral   
( const double xlow , const double xhigh , 
  const double ylow , const double yhigh ) const 
{ 
  return
    s_equal ( xlow  , xmin() ) && 
    s_equal ( xhigh , xmax() ) && 
    s_equal ( ylow  , ymin() ) && 
    s_equal ( yhigh , ymax() )  ?  1.0 :
    m_bernstein.integral ( xlow , xhigh , ylow , yhigh ) ; 
}
// ============================================================================




// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::Positive2DSym::Positive2DSym
( const unsigned short      N    ,
  const double              xmin ,
  const double              xmax )
  : std::binary_function<double,double,double> ()
//
  , m_bernstein (   N , xmin , xmax ) 
  , m_sphere    ( ( N + 1 ) * ( N + 2 ) / 2 - 1  )
{
  updateBernstein () ;
}
// ============================================================================
// set k-parameter
// ============================================================================
bool Ostap::Math::Positive2DSym::setPar 
( const unsigned int k     , 
  const double       value )
{
  //
  const bool update = m_sphere.setPhase ( k , value ) ;
  if ( !update ) { return false ; }   // no actual change 
  //
  return updateBernstein () ;
}
// =============================================================================
// update bernstein coefficients
// =============================================================================
bool Ostap::Math::Positive2DSym::updateBernstein ()
{
  //
  //
  bool update = false ;
  for ( unsigned int ix = 0 ; ix < m_sphere.nX() ; ++ix ) 
  { 
    const bool updated = m_bernstein.setPar ( ix , m_sphere.x2 ( ix ) ) ; 
    update = updated || update ; 
  }
  //
  return update ;
}
// ============================================================================
// get the value
// ============================================================================
double  Ostap::Math::Positive2DSym::operator () 
  ( const double x , const double y ) const 
{ return m_bernstein ( x , y ) ; }
// ============================================================================
// get the parameter value
// ============================================================================
double Ostap::Math::Positive2DSym::par ( const unsigned int  k ) const 
{ return m_sphere.phase ( k ) ; }
// ============================================================================
/*  get the integral over 2D-region 
 *  \f[ \int_{x_{min}}^{x_{max}}\int_{y_{min}}^{y_{max}} 
 *   \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y \f] 
 */
// ============================================================================
double  Ostap::Math::Positive2DSym::integral   () const { return 1 ; }
// ============================================================================
/*  get the integral over 2D-region 
 *  \f[ \int_{x_low}^{x_high}\int_{y_low}^{y_high} 
 *   \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y \f] 
 *  @param xlow  low  edge in x 
 *  @param xhigh high edge in x 
 *  @param ylow  low  edge in y 
 *  @param yhigh high edge in y 
 */
// ============================================================================
double Ostap::Math::Positive2DSym::integral
( const double xlow , const double xhigh , 
  const double ylow , const double yhigh ) const
{
  return 
    s_equal ( xlow  , xmin () ) &&
    s_equal ( xhigh , xmax () ) &&
    s_equal ( ylow  , ymin () ) &&
    s_equal ( yhigh , ymax () ) ?  1.0 :
    m_bernstein.integral ( xlow , xhigh , ylow , yhigh ) ; 
}
// ============================================================================
/* integral over x-dimension 
 *  \f[ \int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}y\f] 
 *  @param x     variable 
 *  @param ylow  low  edge in y 
 *  @param yhigh high edge in y 
 */
// ======================================================================
double  Ostap::Math::Positive2DSym::integrateX 
( const double y    , 
  const double xlow , const double xhigh ) const 
{ return m_bernstein.integrateX ( y , xlow , xhigh ) ; }
// ======================================================================
/*  integral over x-dimension 
 *  \f[ \int_{x_low}^{x_high} \mathcal{B}(x,y) \mathrm{d}x\f] 
 *  @param y     variable 
 *  @param xlow  low  edge in x 
 *  @param xhigh high edge in x 
 */
// ======================================================================
double Ostap::Math::Positive2DSym::integrateY 
( const double x    , 
  const double ylow , const double yhigh ) const 
{ return m_bernstein.integrateY ( x , ylow , yhigh ) ; }
// ======================================================================
/*  integral over x-dimension 
 *  \f[ \int_{x_{min}}^{x_{max}} \mathcal{B}(x,y) \mathrm{d}x\f] 
 *  @param x     variable 
 */
// ======================================================================
double Ostap::Math::Positive2DSym::integrateX ( const double y ) const 
{ return m_bernstein.integrateX ( y ) ; }
// ======================================================================
/*  integral over y-dimension 
 *  \f[ \int_{y_{min}}^{y_{max}} \mathcal{B}(x,y) \mathrm{d}y\f] 
 *  @param y     variable 
 */
// ======================================================================
double Ostap::Math::Positive2DSym::integrateY ( const double x ) const 
{ return m_bernstein.integrateY ( x ) ; }
// ======================================================================



// ============================================================================
// Interpolation stuff 
// ============================================================================
/*  construct interpolation polynomial (in Bernstein form)
 *  @param x       vector of abscissas 
 *  @param y       vector of function values 
 *  @param xmin low  edge for Bernstein polynomial
 *  @param xmax high edge for Bernstein polynomial       
 *  - if vector of y is longer  than vector x, extra values are ignored 
 *  - if vector of y is shorter than vector x, missing entries are assumed to be zero  
 *  It relies on Newton-Bernstein algorithm
 *  @see http://arxiv.org/abs/1510.09197
 *  @see Mark Ainsworth and Manuel A. Sanches, 
 *       "Computing of Bezier control points of Largangian interpolant 
 *       in arbitrary dimension", arXiv:1510.09197 [math.NA]
 *  @see http://adsabs.harvard.edu/abs/2015arXiv151009197A
 *  @see Ostap::Math::Bernstein 
 *  @code 
 *  std::vector<double> x = ... ; // abscissas
 *  std::vector<double> y = ... ; // functionvalues 
 *  Bernstein p = bernstein ( x , y , -1 , 1 );
 *  std::cout << " interpolant at x=0.1 is " << p(0.1) << std::endl ;
 *  std::cout << " interpolant at x=0.2 is " << p(0.2) << std::endl ;
 *  @endcode 
 */
// ============================================================================
Ostap::Math::Bernstein
Ostap::Math::Interpolation::bernstein
( const std::vector<double>& x    ,  
  const std::vector<double>& y    , 
  const double               xmin , 
  const double               xmax )
{
  return bernstein ( x.begin() , x.end() , 
                     y.begin() , y.end() , 
                     xmin      , xmax    ) ;
}
// ============================================================================
/*  construct interpolation polynomial (in Bernstein form)
 *  @param func    the function 
 *  @param x       vector of abscissas 
 *  @param xmin low  edge for Bernstein polynomial
 *  @param xmax high edge for Bernstein polynomial
 *  It relies on Newton-Bernstein algorithm
 *  @see http://arxiv.org/abs/1510.09197
 *  @see Mark Ainsworth and Manuel A. Sanches, 
 *       "Computing of Bezier control points of Largangian interpolant 
 *       in arbitrary dimension", arXiv:1510.09197 [math.NA]
 *  @see http://adsabs.harvard.edu/abs/2015arXiv151009197A
 *  @see Ostap::Math::Bernstein 
 *  @code 
 *  auto f = [] ( double t ) { return std::sin ( t ) ; }
 *  std::vector<double> x = ... ; // abscissas
 *  Bernstein p = bernstein ( f , x , -1 , 1 );
 *  std::cout << " interpolant at x=0.1 is " << p(0.1) << std::endl ;
 *  std::cout << " interpolant at x=0.2 is " << p(0.2) << std::endl ;
 *  @endcode 
 */
// ============================================================================
Ostap::Math::Bernstein
Ostap::Math::Interpolation::bernstein
( std::function<double(double)> func , 
  const std::vector<double>&    x    ,
  const double                  xmin , 
  const double                  xmax ) 
{ return bernstein ( func , x.begin() , x.end() , xmin , xmax ) ; }
// ============================================================================
/*  construct interpolation polynomial (in Bernstein form) using Gauss-Lobatto grid, 
 *  that minimises Runge's effect.
 *  @param func      the function 
 *  @param N         the interpolation  degree 
 *  @param xmin low  edge for Bernstein polynomial
 *  @param xmax high edge for Bernstein polynomial       
 *  It relies on Newton-Bernstein algorithm
 *  @see http://arxiv.org/abs/1510.09197
 *  @see Mark Ainsworth and Manuel A. Sanches, 
 *       "Computing of Bezier control points of Largangian interpolant 
 *       in arbitrary dimension", arXiv:1510.09197 [math.NA]
 *  @see http://adsabs.harvard.edu/abs/2015arXiv151009197A
 *  @see Ostap::Math::Bernstein 
 *  @code 
 *  auto f = [] ( double t ) { return std::sin ( t ) ; }
 *  Bernstein p = bernstein ( f , 5 , -1 , 1 );
 *  std::cout << " interpolant at x=0.1 is " << p(0.1) << std::endl ;
 *  std::cout << " interpolant at x=0.2 is " << p(0.2) << std::endl ;
 *  @endcode 
 */  
// ============================================================================
Ostap::Math::Bernstein
Ostap::Math::Interpolation::bernstein
( std::function<double(double)> func , 
  const unsigned short          N    , 
  const double                  xmin , 
  const double                  xmax ) 
{ return lobatto ( func , N , xmin , xmax ) ; }
// ============================================================================

// ============================================================================
// The END 
// ============================================================================
