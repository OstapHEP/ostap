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
#include "Ostap/Bernstein1D.h"
// ============================================================================
// Local
// ============================================================================
#include "Exception.h"
#include "local_math.h"
// ============================================================================
/** @file 
 *  Implementation file for functions, related to Bernstein's polynomnials 
 *
 *  @see http://en.wikipedia.org/wiki/Bernstein_polynomial
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2010-04-19
 */
// ============================================================================
// Bernstein EVEN 
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
  : m_N         (   N                   )   
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
  : m_N         (   pars.size()                   )   
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
  : m_bernstein ( N , xmin , xmax )
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
  : m_bernstein ( pars.size() , xmin , xmax )
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
  : m_bernstein ( sphere.dim() , xmin , xmax )
  , m_sphere    ( sphere ) 
{
  updateBernstein () ;
}
// ============================================================================
// copy 
// ============================================================================
Ostap::Math::Positive::Positive
( const Ostap::Math::Positive&  right ) 
  : m_bernstein ( right.m_bernstein ) 
  , m_sphere    ( right.m_sphere    ) 
{}
// ============================================================================
// move 
// ============================================================================
Ostap::Math::Positive::Positive
(       Ostap::Math::Positive&& right ) 
  : m_bernstein ( std::move ( right.m_bernstein ) ) 
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
  const double a0     =     m_sphere.x2 ( 0 ) ;
  const double a1_    = 2 * m_sphere.x2 ( 1 ) ; // NB! factor 2 is here! 
  const double a2     =     m_sphere.x2 ( 2 ) ;
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
  : m_even   ( N , xmin , xmax )
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
  : m_even      ( pars.size() , xmin , xmax )
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
  : m_even    ( sphere.dim() , xmin , xmax )
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
// The END 
// ============================================================================
