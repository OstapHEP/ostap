// ============================================================================
// Include files 
// ============================================================================
// STD & STL 
// ============================================================================
#include <numeric>
#include <vector>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Math.h"
#include "Ostap/BSpline.h"
#include "Ostap/Bernstein.h"
#include "Ostap/Lomont.h"
// ============================================================================
// Local
// ============================================================================
#include "Exception.h"
// ============================================================================
/** @file 
 *  Implementation file for class Ostap::Math::BSpline 
 *  @date 2014-11-12 
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 */  
// ============================================================================
namespace 
{
  // ==========================================================================
  const unsigned short s_ulps = Ostap::Math::mULPS_double + 200 ;
  // ==========================================================================
  /** @var s_zero 
   *  comparison with zero
   *  @code
   *  double a = ... ;
   *  if( s_zero ( a ) ) { ... }  
   *  @endcode 
   *  @see Ostap::Math::Zero 
   *  @see Ostap::Math::Lomont
   */
  const Ostap::Math::Zero<double>     s_zero{}  ;
  // ==========================================================================
  /** @var s_equal
   *  comparison of double values 
   *  @code
   *  double a = ... ;
   *  double b = ... ;
   *  if( s_equal ( a , b ) ) { ... }  
   *  @endcode 
   *  @see Ostap::Math::Equal_To 
   *  @see Ostap::Math::Lomont
   */
  const Ostap::Math::Equal_To<double> s_equal{} ;
  // ==========================================================================
  /** @var s_less 
   *  sorting criteria for splines 
   *  @see Ostap::Math::NumLess 
   *  @see Ostap::Math::Equal_To 
   */
  const Ostap::Math::NumLess<double>  s_less{}  ;
}
// ============================================================================
// de Boor-Cox 
// ============================================================================
#include <iostream>
namespace 
{
  // ==========================================================================
  // find an interval such      knot[i]<= x< knot[i+1]  
  template <class IT> 
  inline IT find_i ( IT first , IT last , const double value ) 
  {
    if ( 2 > std::distance ( first , last ) 
         || value < * first  
         || value > *(last-1) ) { return last ; }
    //
    IT l = std::upper_bound ( first , last , value  ) ; // , s_less ) ;
    //
    return 
      last  ==  l ?   l : 
      value <  *l ? --l : 
      value == *l ? ++l : l ;
  }
  // ==========================================================================
  /// get the access to the knots 
  inline double  knot ( const std::vector<double>& knots , const int index ) 
  {
    return 
      0     >         index        ? knots.front () :  
      index >=  (int) knots.size() ? knots.back  () : knots[index] ;
  }
  // ==========================================================================
  /// get the coefficient for de Boor-Cox algorithm
  inline double _sik   ( const unsigned short       i     , 
                         const unsigned short       k     ,
                         const double               x     , 
                         const std::vector<double>& knots ) 
  {
    const double ti   =  knot ( knots , i   ) ;
    const double tik  =  knot ( knots , i+k ) ;
    return ti < tik ?  ( x - ti ) / ( tik - ti )  : 0.0 ;
  }
  // ==========================================================================
  /** calculate the value of the basic spline 
   *  @param i INPUT spline number 
   *  @param k INPUT spline order 
   *  @param x INPUT the argument 
   *  @param know INPUT vector of knots (assumed to be ordered) 
   */
  inline double bspline 
  ( const          short       i     ,  // the spline number   
    const unsigned short       k     ,  // the spline order 
    const double               x     ,  // the spline argument 
    const std::vector<double>& knots )  // knots are assumed to be ordered 
  {
    // check the range 
    if ( 2 > knots.size ()                      ) { return 0 ; }
    if ( x < knots.front()  || x > knots.back() ) { return 0 ; }
    //
    if ( x < knot ( knots , i         ) ) { return 0 ; }
    if ( x > knot ( knots , i + k + 1 ) ) { return 0 ; }
    //
    //
    const double ki = knot ( knots , i ) ;
    // stopping criteria for recursion
    if ( 0 == k ) { return  ki <= x && x < knot ( knots , i + 1 )  ? 1 : 0 ; }
    //
    const double s1 = _sik ( i     , k , x , knots ) ;
    const double s2 = _sik ( i + 1 , k , x , knots ) ;
    //
    const bool null1 = s_zero  (     s1 ) ;
    const bool unit2 = s_equal ( 1 , s2 ) ;
    //
    return 
      //
      null1 && unit2 ?   0.0 : 
      //
      null1          ? ( 1 - s2 ) * bspline ( i + 1 , k - 1 , x , knots ) :
      //
      unit2          ?       s1   * bspline ( i     , k - 1 , x , knots ) :
      //
      bspline ( i      , k - 1 , x , knots ) *       s1   +
      bspline ( i + 1  , k - 1 , x , knots ) * ( 1 - s2 ) ;
  }
  // ==========================================================================
  /** calculate the value of the M-spline 
   *  @param i INPUT spline number 
   *  @param k INPUT spline order 
   *  @param x INPUT the argument 
   *  @param know INPUT vector of knots (assumed to be ordered) 
   */
  inline double mspline 
  ( const          short       i     ,  // the spline number   
    const unsigned short       k     ,  // the spline order 
    const double               x     ,  // the spline argument 
    const std::vector<double>& knots )  // knots are assumed to be ordered 
  {
    // check the range 
    if ( 2 > knots.size ()                      ) { return 0 ; }
    if ( x < knots.front()  || x > knots.back() ) { return 0 ; }
    //
    const double ki  = knot ( knots , i ) ;
    if ( x < ki  ) { return 0 ; } 
    //
    const double kio = knot ( knots , i + k + 1 ) ;
    if ( x > kio ) { return 0 ; }
    //
    // stopping criteria for recursion
    if ( 0 == k  ) { return  1.0 / ( kio - ki ) ; }
    //
    return   
      ( ( x   - ki ) * mspline ( i     , k - 1 , x , knots ) + 
        ( kio -  x ) * mspline ( i + 1 , k - 1 , x , knots ) ) 
      * ( k + 1 ) / ( k * ( kio - ki ) );
  }
  // ==========================================================================
  /** calculate the value of the I-spline 
   *  @param i INPUT spline number 
   *  @param k INPUT spline order 
   *  @param x INPUT the argument 
   *  @param know INPUT vector of knots (assumed to be ordered) 
   */
  inline double ispline 
  ( const          short       i     ,  // the spline number   
    const unsigned short       k     ,  // the spline order 
    const double               x     ,  // the spline argument 
    const std::vector<double>& knots )  // knots are assumed to be ordered 
  {
    // check the range 
    if ( 2 > knots.size ()                      ) { return 0 ; }
    if ( x < knots.front()  || x > knots.back() ) { return 0 ; }
    //
    double result = 0 ;
    for ( int j = i ; knot ( knots , j ) <= x ; ++j ) 
    { result += bspline ( j , k , x , knots ) ; }
    //
    return result ;
  }
  // ==========================================================================
  /// de-boor-cox algorithm
  inline double _deboor_
  ( const unsigned short       k      , 
    const unsigned short       order  , 
    const unsigned short       i      , 
    const double               x      ,
    const std::vector<double>& knots  , 
    const std::vector<double>& pars   ) 
  {
    //
    if ( 0 == k ) { return pars[i] ; }
    //
    const double ti  = knot ( knots , i                 ) ;
    const double tip = knot ( knots , i + 1 + order - k ) ;
    //
    if ( s_equal ( ti , tip ) ) { return 0 ; }
    //
    const double tau = ( x - ti ) / ( tip - ti ) ;
    return 
      _deboor_ ( k-1 , order , i - 1 , x , knots , pars ) * ( 1 - tau ) +
      _deboor_ ( k-1 , order , i     , x , knots , pars ) *       tau   ;  
  }
  // =====================================================================
  unsigned short _insert_ 
  ( const double         x     , 
    const unsigned short n     ,
    std::vector<double>& knots ,
    std::vector<double>& pars  , 
    const unsigned short d     ) 
  {
    const double tf = knots.front ()  ;
    const double tb = knots.back  ()  ;
    //
    if ( x < tf && !s_equal ( x, tf ) ) { return 0 ; }
    if ( x > tb && !s_equal ( x, tb ) ) { return 0 ; }
    //
    std::vector<double>::iterator il = 
      std::lower_bound ( knots.begin () , knots.end   () , x , s_less ) ;
    std::vector<double>::iterator iu = 
      std::upper_bound ( il             , knots.end   () , x , s_less ) ;
    const unsigned short nt = iu -  il ;
    //
    if ( 0 == n ) { return nt ;}
    //
    unsigned short l = ( iu - knots.begin() ) - 1 ;
    //
    //  add knot n-times 
    for  ( unsigned  short i = 0 ; i < n  ; ++i , ++l )
    {
      const unsigned short f =   l  - d  + 1  ;
      //
      pars .insert  ( pars.begin() + l , pars[l] ) ;
      for ( unsigned short j = l ; j >= f ; --j ) 
      {
        const long double tj    = knots [ j ] ;
        const long double dt    = knots [ j + d ] - tj ;
        const long double alpha = ( x - tj ) / dt ;
        pars[j] = (  1 - alpha ) * pars[j-1] + alpha*pars[j] ;
      }
      knots.insert ( knots.begin() + l + 1 , x ) ;
    }
    return  nt + n ;
  }
}
// ============================================================================
// The  Basic spline 
// ============================================================================
/** constructor from the list of knots and the order 
 *  vector of parameters will be calculated automatically 
 *  @param knots  non-empty vector of poinst/knots 
 *  @param order  the order of splines 
 *  - vector of points is not requires to be ordered 
 *  - duplicated knots will be ignored
 *  - min/max value will be used as interval boundaries 
 *  - extra knots will added at the end of interval 
 */
// ============================================================================
Ostap::Math::BSpline::BSpline 
( const std::vector<double>& knots  ,
  const unsigned short       order  ) 
  : m_knots   ( knots  ) 
  , m_pars    () 
  , m_order   ( order  ) 
  , m_inner   ( 0      )  
  , m_xmin    ( 0 ) 
  , m_xmax    ( 1 )
    //
  , m_jlast   ( 0      ) 
  , m_pars_i  ()
  , m_knots_i ()
{
  //
  std::stable_sort ( m_knots.begin() , m_knots.end() , s_less ) ; 
  std::vector<double>::iterator it = std::unique ( m_knots.begin() , m_knots.end() , s_equal );
  m_knots.erase ( it , m_knots.end() ) ;
  //
  if ( m_knots.size () < 2 ) { Ostap::throwException 
      ("Vector of knots is too short", "Ostap::Math::BSpline" , 810 ) ; } 
  //
  m_inner = m_knots.size   () - 2 ;
  // specify vector of parameters 
  m_pars = std::vector<double>( m_inner + m_order + 1 , 0.0 ) ;
  //
  m_xmin = m_knots.front () ;
  m_xmax = m_knots.back  () ;
  //
  while ( m_knots.size() + 1 < m_pars.size() + m_order + 1 ) 
  {
    m_knots.insert ( m_knots.begin () , m_xmin ) ;
    m_knots.insert ( m_knots.end   () , m_xmax ) ;    
  }
  //
  // integration cache:
  m_pars_i .resize( m_pars.size  () + 1 ) ;
  m_knots_i.resize( m_knots.size () + 2 ) ;
  std::copy ( m_knots.begin() , m_knots.end() , m_knots_i.begin() + 1 ) ;
  m_knots_i.front () = m_xmin ;
  m_knots_i.back  () = m_xmax ;
}
// ======================================================================
/*  Constructor from the list of knots and list of parameters 
 *  The spline order will be calculated automatically 
 *  @param knots  non-empty vector of poinst/knots 
 *  @param pars   non-empty vector of parameters 
 *  - vector of points is not requires to be ordered 
 *  - min/max value will be used as interval boundaries 
 *  - duplicated knots will be ignored
 *  - extra knots will added at the end of interval 
 */
// ======================================================================
Ostap::Math::BSpline::BSpline 
( const std::vector<double>& knots  ,
  const std::vector<double>& pars   ) 
  : m_knots ( knots  ) 
  , m_pars  ( pars   ) 
  , m_order ( 1      ) 
  , m_inner ( 0      )  
  , m_xmin  ( 0      ) 
  , m_xmax  ( 1      )
  , m_jlast ( 0      ) 
    //
{
  //
  std::stable_sort ( m_knots.begin() , m_knots.end() , s_less ) ; 
  std::vector<double>::iterator it = std::unique ( m_knots.begin() , m_knots.end() , s_equal );
  m_knots.erase ( it , m_knots.end() ) ;
  //
  if      ( m_knots.size()     <  2              ) 
{ Ostap::throwException
      ("Vector of knots is too short", "Ostap::Math::BSpline" , 812 ) ; } 
  else if ( m_pars .size() + 1 <  m_knots.size() ) { Ostap::throwException
      ("Vector of pars  is too short", "Ostap::Math::BSpline" , 813 ) ; }
  //
  m_inner = m_knots.size () - 2 ;
  m_order = m_pars.size  () - m_inner - 1 ;
  //
  m_xmin = m_knots.front () ;
  m_xmax = m_knots.back  () ;
  //
  while ( m_knots.size() + 1 < m_pars.size() + m_order + 1 ) 
  {
    m_knots.insert ( m_knots.begin () , m_xmin ) ;
    m_knots.insert ( m_knots.end   () , m_xmax ) ;    
  }
  //
  // integration cache:
  m_pars_i .resize ( m_pars.size  () + 1 ) ;
  m_knots_i.resize ( m_knots.size () + 2 ) ;
  std::copy ( m_knots.begin() , m_knots.end() , m_knots_i.begin() + 1 ) ;
  m_knots_i.front () = m_xmin ;
  m_knots_i.back  () = m_xmax ;
}
// ============================================================================
/* Constructor for uniform binning 
 *  @param xmin   low  edge of spline interval 
 *  @param xmax   high edge of spline interval 
 *  @param inner  number of inner points in   (xmin,xmax) interval
 *  @param order  the degree of spline 
 */
// ============================================================================
Ostap::Math::BSpline::BSpline
( const double         xmin    ,  
  const double         xmax    , 
  const unsigned short inner   , // number of inner points  
  const unsigned short order   ) 
  : m_knots () 
  , m_pars  ( order + inner + 1 , 0.0 ) 
  , m_order ( order )  
  , m_inner ( inner )  
  , m_xmin  ( std::min ( xmin  , xmax ) ) 
  , m_xmax  ( std::max ( xmin  , xmax ) ) 
  , m_jlast ( 0      ) 
{
  //
  const double dx = ( m_xmax - m_xmin ) ;
  for ( unsigned short i = 1 ; i <= inner ; ++i ) 
  { m_knots .push_back ( m_xmin + dx * double ( i )  / ( inner + 1 ) ) ; }
  //
  while ( m_knots.size() + 1 < m_pars.size() + m_order + 1 ) 
  {
    m_knots.insert ( m_knots.begin () , m_xmin ) ;
    m_knots.insert ( m_knots.end   () , m_xmax ) ;    
  }
  //
  // integration cache:
  m_pars_i .resize( m_pars.size  () + 1 ) ;
  m_knots_i.resize( m_knots.size () + 2 ) ;
  std::copy ( m_knots.begin() , m_knots.end() , m_knots_i.begin() + 1 ) ;
  m_knots_i.front () = m_xmin ;
  m_knots_i.back  () = m_xmax ;
}
// ============================================================================
/*  Constructor from the list of knots and list of parameters 
 *  The spline order will be calculated automatically 
 *  @param knots  non-empty vector of poinst/knots 
 *  @param pars   non-empty vector of parameters 
 *  - vector of points is not requires to be ordered 
 *  - min/max value will be used as interval boundaries 
 *  - duplicated knots will be ignored
 *  - extra knots will added at the end of interval 
 */
// ============================================================================
Ostap::Math::BSpline::BSpline 
( const Ostap::Math::BSpline& b   , 
  const double                xmn , 
  const double                xmx )
  : m_knots ( b.m_knots ) 
  , m_pars  ( b.m_pars  ) 
  , m_order ( b.m_order ) 
  , m_inner ( 0      )  
  , m_xmin  ( std::min ( xmn , xmx ) )  
  , m_xmax  ( std::min ( xmn , xmx ) )  
  , m_jlast ( 0      ) 
    //
{
  //
  m_xmin = std::max ( b.xmin () , std::min ( xmn ,  xmx ) ) ;
  m_xmax = std::min ( b.xmax () , std::max ( xmn ,  xmx ) ) ;
  //
  if      ( b.xmin() < m_xmin || s_equal (  b.xmin() , m_xmin ) ) { /* ok */ }
  else if ( b.xmax() > m_xmax || s_equal (  b.xmax() , m_xmax ) ) { /* ok */ }
  else { Ostap::throwException 
      ("Invalid xmin/xmax  settings", "Ostap::Math::BSpline" , 821 ) ; }
  //
  // ==========================================================================
  std::vector<double>::iterator il = 
    std::lower_bound ( m_knots.begin () , m_knots.end () , m_xmin , s_less ) ;
  std::vector<double>::iterator iu = 
    std::upper_bound ( il               , m_knots.end () , m_xmin , s_less ) ;
  //
  const unsigned short nxmin = iu - il ;
  if ( nxmin < m_order + 1 ) 
  { _insert_ ( m_xmin , 1 + m_order - nxmin , m_knots , m_pars , m_order ) ; }
  //
  il = std::lower_bound ( m_knots.begin () , m_knots.end () , m_xmax , s_less ) ;
  iu = std::upper_bound ( il               , m_knots.end () , m_xmax , s_less ) ;
  //
  const unsigned short nxmax = iu - il ;
  if ( nxmax < m_order + 1 ) 
  { _insert_ ( m_xmax , 1 + m_order - nxmax , m_knots , m_pars , m_order ) ; }
  //
  // range of knots   
  il = std::lower_bound ( m_knots.begin ()  , m_knots.end () , m_xmin , s_less ) ;
  iu = std::upper_bound ( il                , m_knots.end () , m_xmax , s_less ) ;
  //
  const unsigned  short i1 = il - m_knots.begin() ;
  const unsigned  short i2 = iu - m_knots.begin() ;
  //
  m_knots.erase ( m_knots.begin() + i2 , m_knots.end   ()      ) ;
  m_knots.erase ( m_knots.begin()      , m_knots.begin () + i1 ) ;
  //
  m_inner                  = ( i2 - i1 ) -  2 *  ( m_order + 1 ) ;
  const unsigned  short np = ( i2 - i1 ) -  m_order - 1 ;
  const unsigned  short fp =  i1 ;
  //
  m_pars = std::vector<double>( m_pars.begin() + fp      , 
                                m_pars.begin() + fp + np ) ;                              
  //
  // integration cache:
  //
  m_pars_i .resize ( m_pars.size  () + 1 ) ;
  m_knots_i.resize ( m_knots.size () + 2 ) ;
  std::copy ( m_knots.begin() , m_knots.end() , m_knots_i.begin() + 1 ) ;
  m_knots_i.front () = m_xmin ;
  m_knots_i.back  () = m_xmax ;
}
// ============================================================================(
Ostap::Math::BSpline::BSpline 
( const Ostap::Math::Bernstein& b ) 
  : Ostap::Math::BSpline ( {{ b.xmin() , b.xmax() }} , b.pars() )
{}
// ============================================================================
// move constructor 
// ============================================================================
Ostap::Math::BSpline::BSpline( Ostap::Math::BSpline&&  right ) 
  : m_knots   ( std::move ( right.m_knots   ) ) 
  , m_pars    ( std::move ( right.m_pars    ) ) 
  , m_order   ( std::move ( right.m_order   ) ) 
  , m_inner   ( std::move ( right.m_inner   ) ) 
  , m_xmin    ( std::move ( right.m_xmin    ) ) 
  , m_xmax    ( std::move ( right.m_xmax    ) ) 
  , m_jlast   ( std::move ( right.m_jlast   ) ) 
  , m_pars_i  ( std::move ( right.m_pars_i  ) )  
  , m_knots_i ( std::move ( right.m_knots_i ) )  
{}
// ============================================================================
/// assignement move operator 
// ============================================================================
Ostap::Math::BSpline& Ostap::Math::BSpline::operator=
(       Ostap::Math::BSpline&& right ) 
{
  if ( &right == this ) { return *this ; }
  //
  m_knots   = std::move ( right.m_knots   ) ;
  m_pars    = std::move ( right.m_pars    ) ;
  m_order   = std::move ( right.m_order   ) ;
  m_inner   = std::move ( right.m_inner   ) ;
  m_xmin    = std::move ( right.m_xmin    ) ;
  m_xmax    = std::move ( right.m_xmax    ) ;
  m_jlast   = std::move ( right.m_jlast   ) ;
  m_pars_i  = std::move ( right.m_pars_i  ) ;
  m_knots_i = std::move ( right.m_knots_i ) ;
  //
  return *this ;
}
// ============================================================================
/*  calculate q-norm of the spline 
 *  where q-norm is defined as:
 *  \f$ \left| f \right|_{q} = \left( \sum_i \left|c_i\right|^q\right)^{\frac{1}{q}} \f$
 *  
 *  - q_inv = 0.0 ->  \f$ max_k    \left|c_k\right|  \f$ 
 *  - q_inv = 0.5 ->  \f$ \sqrt{ \sum_k  c_k^2 }     \f$
 *  - q_inv = 1.0 ->  \f$ \sum_k \left| c_k \right|  \f$ 
 */
// ============================================================================
double Ostap::Math::BSpline::norm   ( const double q_inv ) const 
{ return Ostap::Math::p_norm ( m_pars.begin() , m_pars.end() , q_inv ) ; }
// ============================================================================
// scale all coefficients with 2**i
// ============================================================================
Ostap::Math::BSpline
Ostap::Math::BSpline::ldexp ( const short i )  const 
{
  if ( 0 == i ) { return *this ; }
  Ostap::Math::BSpline result ( *this ) ;
  Ostap::Math::scale_exp2 ( result.m_pars   , i ) ;
  Ostap::Math::scale_exp2 ( result.m_pars_i , i ) ;
  return result ;
}
// ============================================================================
// Greville's abscissas 
// ============================================================================
std::vector<double> Ostap::Math::BSpline::greville_abscissas () const 
{
  std::vector<double> ga ( npars() ) ;
  const unsigned short Na = npars()  ;
  //
  const unsigned short o  = order() ? order() : 1 ;
  for ( unsigned short i  = 0 ; i < Na ; ++i ) 
  {
    for ( unsigned short j = i ; j < i + o ;  ++j ) { ga[i] += m_knots[j+1] ; }
    ga[i] /= o ;
  }
  return ga ;
}
// ============================================================================
// Greville's abscissa
// ============================================================================
double Ostap::Math::BSpline::greville_abscissa ( const unsigned short i ) const 
{
  //
  const unsigned short o =  order() ? order() : 1 ;
  if ( i >= npars() ) { return 2 * xmin () - xmax () ; }
  long double ga = 0 ;
  for ( unsigned short j = i ; j < i + o ;  ++j ) { ga += m_knots[j+1] ; }
  ga /= o ;
  return ga ;
}
// ============================================================================
// is it a increasing function?
// ============================================================================
bool Ostap::Math::BSpline::increasing   () const 
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
bool Ostap::Math::BSpline::decreasing   () const 
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
bool Ostap::Math::BSpline::constant () const 
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
// simple  manipulations with bernstein polynoms: scale it! 
// ============================================================================
Ostap::Math::BSpline&
Ostap::Math::BSpline::operator*=( const double a ) 
{
  Ostap::Math::scale ( m_pars  , a ) ;
  return *this ;
}
// ============================================================================
// simple  manipulations with bernstein polynoms: scale it! 
// ============================================================================
Ostap::Math::BSpline&
Ostap::Math::BSpline::operator/=( const double a ) 
{
  Ostap::Math::scale ( m_pars  , 1.0/a ) ;
  return *this ;
}
// ============================================================================
// simple  manipulations with bernstein polynoms: shift it! 
// ============================================================================
Ostap::Math::BSpline&
Ostap::Math::BSpline::operator+=( const double a ) 
{
  Ostap::Math::shift ( m_pars  , a ) ;
  return *this ;
}
// ============================================================================
// simple  manipulations with bernstein polynoms: shift it! 
// ============================================================================
Ostap::Math::BSpline&
Ostap::Math::BSpline::operator-=( const double a ) 
{
  Ostap::Math::shift ( m_pars  , -a ) ;
  return *this ;
}
// ============================================================================
Ostap::Math::BSpline
Ostap::Math::BSpline::operator-() const 
{
  BSpline b ( *this ) ;
  Ostap::Math::negate ( b.m_pars  ) ;
  return b ;
}
// ============================================================================
// Sum of B-spline and a constant 
// ============================================================================
Ostap::Math::BSpline
Ostap::Math::BSpline::__add__   ( const double value ) const 
{ return (*this) + value ; }
// ============================================================================
// Sum of B-spline and a constant 
// ============================================================================
Ostap::Math::BSpline
Ostap::Math::BSpline::__radd__  ( const double value ) const 
{ return value + (*this) ; }
// ============================================================================
// Product of B-spline and a constant
// ============================================================================
Ostap::Math::BSpline
Ostap::Math::BSpline::__mul__   ( const double value ) const 
{ return (*this) * value ; }
// ============================================================================
// Product of B-spline and a constant
// ============================================================================
Ostap::Math::BSpline
Ostap::Math::BSpline::__rmul__  ( const double value ) const 
{ return value * (*this) ; }
// ============================================================================
// Subtract a constant from B-spline 
// ============================================================================
Ostap::Math::BSpline
Ostap::Math::BSpline::__sub__  ( const double value ) const 
{ return (*this) - value ; }
// ============================================================================
// Subtract B-spline from a constant 
// ============================================================================
Ostap::Math::BSpline
Ostap::Math::BSpline::__rsub__ ( const double value ) const 
{ return value - (*this) ; }
// ============================================================================
// Divide B-spline by a constant 
// ============================================================================
Ostap::Math::BSpline
Ostap::Math::BSpline::__div__   ( const double value ) const 
{ return (*this) / value ; }
// ============================================================================
// Negate Bernstein polynomial 
// ============================================================================
Ostap::Math::BSpline
Ostap::Math::BSpline::__neg__ ()  const { return -(*this); }
// ============================================================================
// set k-parameter
// ============================================================================
bool Ostap::Math::BSpline::setPar ( const unsigned short k , const double value ) 
{
  if ( k >= m_pars.size()          ) { return false ; }
  if ( s_equal ( m_pars[k] , value ) ) { return false ; }
  //
  m_pars[k ] = value ;
  return true ;
}
// ============================================================================
// get the value of the B-spline  i at point x 
// ============================================================================
double Ostap::Math::BSpline::bspline ( const short  i , 
                                       const double x )  const 
{ 
  return  
    x < m_xmin || x > m_xmax ? 0.0 :
    ::bspline ( i , m_order , x , m_knots ) ;
}
// ============================================================================
// get the value of the B-spline  (i,k) at point x
// ============================================================================
double Ostap::Math::BSpline::bspline ( const          short i , 
                                       const unsigned short k , 
                                       const double         x )  const 
{ 
  return 
    x < m_xmin || x > m_xmax ? 0.0 :
    ::bspline ( i , k , x , m_knots ) ;
}
// ============================================================================
// get the value of the M-spline  i at point x
// ============================================================================
double Ostap::Math::BSpline::mspline ( const short  i , 
                                       const double x )  const 
{ 
  return 
    x < m_xmin || x > m_xmax ? 0.0 :
    ::mspline ( i , m_order , x , m_knots ) ; 
}
// ============================================================================
// get the value of the M-spline  (i,k) at point x 
// ============================================================================
double Ostap::Math::BSpline::mspline ( const          short i , 
                                       const unsigned short k , 
                                       const double         x )  const 
{ 
  return 
    x < m_xmin || x > m_xmax ? 0.0 :
    ::mspline ( i , k , x , m_knots ) ; 
}
// ============================================================================
// get the value of the I-spline  i at point x 
// ============================================================================
double Ostap::Math::BSpline::ispline ( const short  i , 
                                       const double x )  const 
{
  return 
    x < m_xmin || m_xmax > x  ? 0.0 :
    ::ispline ( i , m_order , x , m_knots ) ; 
}
// ============================================================================
// get the value of the M-spline  (i,k) at point x 
// ============================================================================
double Ostap::Math::BSpline::ispline ( const          short i , 
                                       const unsigned short k , 
                                       const double         x )  const 
{ 
  return 
    x < m_xmin || x > m_xmax ? 0.0 :
    ::ispline ( i , k , x , m_knots ) ; 
}
// ============================================================================
// get the value
// ============================================================================
double Ostap::Math::BSpline::operator () ( const double x ) const
{
  if      ( x < m_xmin || x > m_xmax ) { return 0 ; }             // RETURN
  //
  // endpoints 
  //
  if      ( s_equal ( x , m_xmin ) ) { return m_pars.front () ; }
  else if ( s_equal ( x , m_xmax ) ) { return m_pars.back  () ; }
  //
  // const double arg =  
  //   !s_equal ( x , m_xmax ) ? x :
  //   0 <= m_xmax             ? 
  //   Ostap::Math::next_double ( m_xmax , -s_ulps ) :
  //   Ostap::Math::next_double ( m_xmax ,  s_ulps ) ;
  //
  const double arg = x ;
  //
  // find the proper "j"
  //
  // 1) try from jlast
  if ( x <  m_knots[ m_jlast ] || x >= m_knots[ m_jlast + 1 ] ) 
  {
    m_jlast = find_i ( m_knots.begin () , 
                       m_knots.end   () , arg ) - m_knots.begin() ;
  }
  //
  // straightforward calculations: 
  // double result = 0 ;
  // for ( unsigned short j = m_jlast - m_order ; j <= m_jlast ; ++j ) 
  // { result += bspline ( j , arg ) * m_pars[j] ; }
  //
  // use de Boor-Cox algorithm:
  return _deboor_ ( m_order , m_order , m_jlast , arg , m_knots , m_pars ) ;
  //
}
// ============================================================================

// ============================================================================
namespace 
{
  // ==========================================================================
  inline double   
  _spline_integral_ 
  ( const std::vector<double>& pars  , 
    const std::vector<double>& knots ,
    const unsigned short       order ) 
  {
    double             result = 0 ;
    for ( unsigned short i = 0 ; i < pars.size() ; ++i ) 
    { result +=  pars[i] * ( knots[ i + order + 1 ] - knots [ i ] ) ; }
    //
    return result / ( order + 1 ) ;
  }
  // ==========================================================================
}
// ============================================================================
// get the integral between xmin and xmax
// ============================================================================
double Ostap::Math::BSpline::integral () const
{
  //
  return _spline_integral_ ( m_pars , m_knots , m_order ) ;
  //
  // double             result = 0 ;
  // for ( unsigned short i = 0 ; i < m_pars.size() ; ++i ) 
  // { result +=  m_pars[i] * ( m_knots[ i + m_order + 1 ] - m_knots [ i ] ) ; }
  //
  // return result / ( m_order + 1 ) ;
}
// ============================================================================
// get the integral between low and high 
// ============================================================================
double Ostap::Math::BSpline::integral 
( const double low  , 
  const double high ) const 
{
  //
  if      ( high < low                      ) { return - integral ( high   , low    ) ; }
  else if ( s_equal ( low , high )          ) { return 0 ; }
  else if ( low  <  m_xmin                  ) { return   integral ( m_xmin , high   ) ; }
  else if ( high >  m_xmax                  ) { return   integral ( low    , m_xmax ) ; }
  else if ( high <= m_xmin || low >= m_xmax ) { return 0 ; }
  else if ( s_equal ( low  , m_xmin ) && s_equal ( high , m_xmax ) ) { return integral() ; }
  //
  const double xhigh =
    !s_equal ( high , m_xmax ) ? high : 
    Ostap::Math::next_double ( m_xmax , -s_ulps ) ;
  //
  // make the integration:
  m_pars_i[0] = 0 ;
  for ( unsigned int i = 0 ; i < m_pars.size() ; ++i ) 
  { m_pars_i[i+1] = m_pars_i[i] + m_pars[i] * ( m_knots[ i + m_order + 1 ] - m_knots [ i ] ) ; }
  //
  const short  jL = find_i ( m_knots_i.begin () , m_knots_i.end () ,  low  ) - m_knots_i.begin() ;  
  const short  jH = find_i ( m_knots_i.begin () , m_knots_i.end () , xhigh ) - m_knots_i.begin() ;
  //
  const double rL = _deboor_ ( m_order + 1 , m_order + 1 , jL ,  low  , m_knots_i , m_pars_i ) ;
  const double rH = _deboor_ ( m_order + 1 , m_order + 1 , jH , xhigh , m_knots_i , m_pars_i ) ;
  //
  return  ( rH - rL ) / ( m_order + 1 ) ;
}
// ============================================================================
// get the integral   as function object 
// ============================================================================
Ostap::Math::BSpline 
Ostap::Math::BSpline::indefinite_integral ( const double C ) const 
{
  // create the object 
  BSpline result ( m_xmin , m_xmax  , m_inner  , m_order + 1 ) ;
  // fill it properly 
  std::copy ( m_knots.begin () , m_knots.end  () , result.m_knots  .begin() + 1 ) ;
  std::copy ( m_knots.begin () , m_knots.end  () , result.m_knots_i.begin() + 2 ) ;
  //
  *( result.m_knots  .begin()     ) = m_knots.front () ;
  *( result.m_knots  .end  () - 1 ) = m_knots.back  () ;
  *( result.m_knots_i.begin()     ) = m_knots.front () ;
  *( result.m_knots_i.begin() + 1 ) = m_knots.front () ;
  *( result.m_knots_i.end  () - 1 ) = m_knots.back  () ;
  *( result.m_knots_i.end  () - 2 ) = m_knots.back  () ;
  //
  result.m_knots_i [0] = m_knots[0] ;
  result.m_knots_i [1] = m_knots[0] ;
  //
  result.m_pars[0] = C ;
  for ( unsigned int i = 0 ; i < m_pars.size() ; ++i ) 
  { result.m_pars[i+1] = result.m_pars[i] 
      + m_pars[i] * ( m_knots[ i + m_order + 1 ] - m_knots [ i ] ) / ( m_order + 1 ) ; }
  //
  return result ;
}
// ============================================================================
// get the derivative at point "x" 
// ============================================================================
double Ostap::Math::BSpline::derivative ( const double x   ) const 
{
  //
  if      ( x < m_xmin || x > m_xmax || 0 == m_order ) { return 0 ; }
  //
  const double arg = 
    !s_equal ( x , m_xmax ) ? x :
    Ostap::Math::next_double ( m_xmax , -s_ulps ) ;
  
  //
  // make the differentiation 
  //
  m_pars_i[0] = m_pars[0]  ;
  for ( unsigned int i = 1  ; i < m_pars.size() ; ++i ) 
  {m_pars_i[i] = ( m_pars[i] - m_pars[i-1] ) / ( m_knots [ i + m_order ] - m_knots [ i ] ) ; }
  //
  // try from jlast
  if ( arg <  m_knots[ m_jlast ] || arg >= m_knots[ m_jlast + 1 ] ) 
  { m_jlast = find_i ( m_knots.begin () , m_knots.end   () , arg ) - m_knots.begin() ; }
  //
  const double r = _deboor_ ( m_order - 1 , m_order - 1 , m_jlast , arg , m_knots , m_pars_i ) ;
  //
  return r * m_order ;
}
// ============================================================================
// get the derivative as function object 
// ============================================================================
Ostap::Math::BSpline Ostap::Math::BSpline::derivative () const 
{
  // create the object 
  BSpline result ( m_xmin , m_xmax , m_inner , 0 < m_order ? m_order -1 : 0  ) ;
  if ( 0 == m_order ) { return result ; }
  //
  std::copy ( m_knots.begin () + 1 , m_knots.end  () - 1  , result.m_knots . begin () ) ;
  //
  std::vector<double> _pars ( m_pars.size  () ) ;
  _pars[0] = m_pars[0] * m_order ;
  for ( unsigned int i = 1  ; i < m_pars.size() ; ++i ) 
  { _pars[i] = ( m_pars[i] - m_pars[i-1] ) * m_order / 
      ( m_knots [ i + m_order ] - m_knots [ i ] ) ; }
  //
  std::copy ( _pars.begin() + 1 , _pars.end() , result.m_pars.begin() ) ;
  //
  return result ;
}
// ============================================================================
/* insert new (unique) knot into the list of knots 
 * @param t new knot  to be inserted 
 * @return true if knot is indeed inserted 
 */
// ============================================================================
bool Ostap::Math::BSpline::insert ( const double t ) 
{
  //
  if ( t < xmin () || s_equal ( t  , xmin () ) ) { return false ; }
  if ( t > xmax () || s_equal ( t  , xmax () ) ) { return false ; }
  //
  std::vector<double>::iterator il = 
    std::lower_bound ( m_knots.begin () , m_knots.end   () , t , s_less ) ;
  std::vector<double>::iterator iu = 
    std::upper_bound ( il               , m_knots.end   () , t , s_less ) ;
  //
  // such not is laready in the list! 
  if ( iu != il ) { return false ; }
  //
  return _insert_ ( t , 1 , m_knots , m_pars , m_order ) ;
}
// ============================================================================
/*  calculate the value of spline defined by vector of knot and vector of 
 *  points using de-boor-cox algorithm
 *  @see https://en.wikipedia.org/wiki/De_Boor%27s_algorithm
 *  @param x     (INPUT) value of x 
 *  @param order (INPUT) the order of spline 
 *  @param knots (INPUT) the vector of knots 
 *  @param pars  (INPUT) the vector of control points 
 *  @return the valeu of b-spline at point x 
 */
// ============================================================================
double Ostap::Math::deboor
( const double               x     ,      
  const unsigned short       order , 
  const std::vector<double>& knots , 
  const std::vector<double>& pars  ) 
{
  //
  const double tf = knots.front ()  ;
  const double tb = knots.back  ()  ;
  if ( x < tf || x > tb    ||  
       s_equal ( x , tf )  ||
       s_equal ( x , tb )   ) { return 0 ; }
  //
  const unsigned short j = find_i ( knots.begin () , knots.end   () , x ) - knots.begin() ;
  //
  return _deboor_ ( order , order , j , x , knots , pars ) ;  
}
// ============================================================================
/* insert new knot at position x in the spline, defined by 
 *  knot vector knots, vector of control points pars and the order
 *  Boehm's algorithm is used 
 *  @see W.Boehm, ``Inserting new knots into B-spline curves'',
 *       Computer-Aided Design, 12, no.4, (1980) 199 
 *  @see http://dx.doi.org/10.1016/0010-4485(80)90154-2
 *  @see http://www.sciencedirect.com/science/article/pii/0010448580901542
 *  @param x     (INPUT)  position of new knot 
 *  @param knots (UPDATE) vector of knots 
 *  @param pars  (UPDATE) vector of control points 
 *  @param order (INPUT)  degree/order of spline 
 *  @param num   (INPUT)  inser knot "num"-times 
 *  @return multiplicity of inserted knot  
 */
// ============================================================================
unsigned short 
Ostap::Math::boehm ( const double         x     , 
                     std::vector<double>& knots ,
                     std::vector<double>& pars  , 
                     const unsigned short order , 
                     const unsigned short num   ) 
{ return _insert_  ( x , num , knots , pars , order ) ; }
// ============================================================================
/*  get a vector of knots from Greville abscissas 
 *  @param aabscissas (INPUT) vector of greville's abscissas 
 *  @param degree of the spline 
 *  @return vector of knots 
 */
// ============================================================================
std::vector<double> 
Ostap::Math::knots_from_abscissas ( std::vector<double> abscissas , 
                                    const  unsigned short degree  ) 
{
  // can't reconstruct for degree-0
  if       ( 0 == degree ) { return std::vector<double>() ; }              // RETURN
  // vector of abscissas is too short 
  if ( abscissas.size() < degree  + 1 ) { return std::vector<double>() ; } // RETURN
  //
  std::stable_sort ( abscissas.begin() , abscissas.end() , s_less ) ; 
  std::vector<double>::iterator it = 
    std::unique    ( abscissas.begin() , abscissas.end() , s_equal );
  abscissas.erase  ( it , abscissas.end() ) ;
  //
  // vector of abscissas is too short 
  if ( abscissas.size() < degree  + 1 ) { return std::vector<double>() ; }  // RETURN
  //
  // trivial           for degree-1
  if  ( 1 == degree ) 
  {
    abscissas.insert ( abscissas.begin() , abscissas.front () ) ;
    abscissas.insert ( abscissas.end  () , abscissas.back  () ) ;
    return abscissas             ;
  }
  //
  // get the vector of knots 
  std::vector<double> knots  ( degree + 1 , abscissas.front() ) ;
  knots.reserve ( abscissas.size() + 2 * degree + 2  ) ;
  const unsigned int N =  abscissas.size() ;
  for ( unsigned int i = 1 ; i < N ; ++i ) 
  {
    //  sum of last "degree" terms 
    const double sumt = 
      std::accumulate ( knots.crbegin() , knots.crbegin() + ( degree - 1 ) , 0.0 ) ;
    const double ti   =  abscissas[i] *  degree - sumt ;
    knots.push_back (  ti ) ;  
  }
  //
  knots.push_back ( abscissas.back() ) ;
  return knots ;
}
// ========================================================================

// ============================================================================
// ============================================================================
// POSITIVE SPLINE 
// ============================================================================
// ============================================================================
/* constructor from the list of knots and the order 
 *  vector of parameters will be calculated automatically 
 *  @param points non-empty vector of poinst/knots 
 *  @param order  the order of splines 
 *  - vector of points is not requires to be ordered 
 *  - duplicated knots will be ignored
 *  - min/max value will be used as interval boundaries 
 */
// ============================================================================
Ostap::Math::PositiveSpline::PositiveSpline 
( const std::vector<double>& points ,
  const unsigned short       order  ) 
  : m_bspline ( points , order ) 
  , m_sphere  ( 1 , 3 ) 
{
  //
  if ( 2 > m_bspline.npars() ) { Ostap::throwException
      ( "At least two spline parameters are required" , 
        "Ostap::Math::PositiveSpline"                 , 814 ) ; }
  //
  m_sphere = Ostap::Math::NSphere( m_bspline.npars() - 1 , 3 ) ;  
  updateCoefficients () ;
}
// ============================================================================
/*  Constructor from the list of knots and list of parameters 
 *  The spline order will be calculated automatically 
 *  @param points non-empty vector of poinst/knots 
 *  @param pars   non-empty vector of parameters 
 *  - vector of points is not requires to be ordered 
 *  - duplicated knots will be ignored
 *  - min/max value will be used as interval boundaries 
 */
// ============================================================================
Ostap::Math::PositiveSpline::PositiveSpline 
( const std::vector<double>& points    ,
  const std::vector<double>& pars      ) 
  : m_bspline ( points , std::vector<double>( pars.size() + 1 , 0 )  ) 
  , m_sphere  ( pars   , 3 ) 
{
  //
  if ( 2 > m_bspline.npars() ) { Ostap::throwException
      ( "At least two spline parameters are required" , 
        "Ostap::Math::PositiveSpline"                 , 814 ) ; }
  //
  m_sphere = Ostap::Math::NSphere( m_bspline.npars() - 1 , true ) ;  
  updateCoefficients () ;
}
// ============================================================================
/*  Constructor for uniform binning 
 *  @param xmin   low  edge of spline interval 
 *  @param xmax   high edge of spline interval 
 *  @param inner  number of inner points in   (xmin,xmax) interval
 *  @param order  the degree of splline 
 */
// ============================================================================
Ostap::Math::PositiveSpline::PositiveSpline 
( const double         xmin  ,  
  const double         xmax  ,  
  const unsigned short inner ,   // number of inner points 
  const unsigned short order ) 
  : m_bspline ( xmin , xmax , inner , order ) 
  , m_sphere  ( 1 , 3 ) 
{
  //
  if ( 2 > m_bspline.npars() ) { Ostap::throwException
      ( "At least two spline parameters are required" , 
        "Ostap::Math::PositiveSpline"                 , 814 ) ; }
  //
  m_sphere = Ostap::Math::NSphere( m_bspline.npars() - 1 , 3 ) ;  
  updateCoefficients() ;
}
// ============================================================================
// constructor fomr the basic spline 
// ============================================================================
Ostap::Math::PositiveSpline::PositiveSpline 
( const Ostap::Math::BSpline& spline ) 
  : m_bspline ( spline ) 
  , m_sphere  ( 1 , 3  ) 
{
  //
  if ( 2 > m_bspline.npars() ) { Ostap::throwException
      ( "At least two spline parameters are required" , 
        "Ostap::Math::PositiveSpline"                 , 814 ) ; }
  //
  m_sphere = Ostap::Math::NSphere( m_bspline.npars() - 1 , 3 ) ;  
  updateCoefficients() ;
}
// ============================================================================
Ostap::Math::PositiveSpline::~PositiveSpline(){}
// ============================================================================
// update coefficients  
// ============================================================================
bool Ostap::Math::PositiveSpline::updateCoefficients  () 
{
  //
  bool   update = false ;
  //
  // get sphere coefficients 
  std::vector<double> v( m_sphere.nX() ) ;
  for ( unsigned short ix = 0 ; ix < m_sphere.nX() ; ++ix ) 
  {  v [ix] = m_sphere.x2 ( ix ) ; }
  //
  const double isum = 1.0 / 
    _spline_integral_ ( v , m_bspline.knots() , m_bspline.order() ) ;
  //
  for ( unsigned short ix = 0 ; ix < m_sphere.nX() ; ++ix ) 
  { 
    const bool updated = m_bspline.setPar ( ix , v [ix] * isum ) ; 
    update = updated || update ;
  }
  //
  return update ;
}
// ============================================================================
// set k-parameter
// ============================================================================
bool Ostap::Math::PositiveSpline::setPar 
( const unsigned short k , const double value ) 
{
  //
  const bool update = m_sphere.setPhase ( k , value ) ;
  if ( !update ) { return false ; }   // no actual change 
  //
  return updateCoefficients () ;
}
// ============================================================================
// get the integral between xmin and xmax
// ============================================================================
double  Ostap::Math::PositiveSpline::integral   () const { return 1 ; }
// ============================================================================
// get the integral between low and high 
// ============================================================================
double  Ostap::Math::PositiveSpline::integral   
( const double low  , 
  const double high ) const 
{ 
  return 
    s_equal ( low  , xmin() ) && s_equal ( high , xmax() ) ? 1 :
    m_bspline.integral   ( low , high ) ; 
}
// ============================================================================


// ============================================================================
// MONOTONIC SPLINE 
// ============================================================================
// ============================================================================
/* constructor from the list of knots and the order 
 *  vector of parameters will be calculated automatically 
 *  @param points non-empty vector of poinst/knots 
 *  @param order  the order of splines 
 *  - vector of points is not requires to be ordered 
 *  - duplicated knots will be ignored
 *  - min/max value will be used as interval boundaries 
 */
// ============================================================================
Ostap::Math::MonotonicSpline::MonotonicSpline
( const std::vector<double>& points      ,
  const unsigned short       order       , 
  const bool                 increasing  )
  : Ostap::Math::PositiveSpline ( points , order ) 
  , m_increasing                ( increasing ) 
{
  updateCoefficients () ;
}
// ============================================================================
/*  Constructor from the list of knots and list of parameters 
 *  The spline order will be calculated automatically 
 *  @param points non-empty vector of poinst/knots 
 *  @param pars   non-empty vector of parameters 
 *  - vector of points is not requires to be ordered 
 *  - duplicated knots will be ignored
 *  - min/max value will be used as interval boundaries 
 */
// ============================================================================
Ostap::Math::MonotonicSpline::MonotonicSpline
( const std::vector<double>& points     ,
  const std::vector<double>& pars       ,
  const bool                 increasing ) 
  : Ostap::Math::PositiveSpline ( points , pars ) 
  , m_increasing                ( increasing    ) 
{
  updateCoefficients () ;
}
// ============================================================================
/*  Constructor for uniform binning 
 *  @param xmin   low  edge of spline interval 
 *  @param xmax   high edge of spline interval 
 *  @param inner  number of inner points in   (xmin,xmax) interval
 *  @param order  the degree of splline 
 */
// ============================================================================
Ostap::Math::MonotonicSpline::MonotonicSpline
( const double         xmin       ,  
  const double         xmax       , 
  const unsigned short inner      ,   // number of inner points 
  const unsigned short order      , 
  const bool           increasing ) 
  : Ostap::Math::PositiveSpline ( xmin , xmax  , inner ,order ) 
  , m_increasing                ( increasing    ) 
{
  updateCoefficients () ;
}
// ============================================================================
// constructor from the basic spline 
// ============================================================================
Ostap::Math::MonotonicSpline::MonotonicSpline
( const Ostap::Math::PositiveSpline& spline     , 
  const bool                         increasing ) 
  : Ostap::Math::PositiveSpline ( spline        ) 
  , m_increasing                ( increasing    ) 
{
  updateCoefficients () ;
}
// ============================================================================
// constructor from the basic spline 
// ============================================================================
Ostap::Math::MonotonicSpline::MonotonicSpline
( const Ostap::Math::BSpline&        spline     , 
  const bool                         increasing ) 
  : Ostap::Math::PositiveSpline ( spline        ) 
  , m_increasing                ( increasing    ) 
{
  updateCoefficients () ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::MonotonicSpline::~MonotonicSpline(){}
// ============================================================================
// update coefficients  
// ============================================================================
bool Ostap::Math::MonotonicSpline::updateCoefficients  () 
{
  //
  bool   update = false ;
  //
  // get sphere coefficients 
  std::vector<double> v ( m_sphere.nX() ) ;
  for ( unsigned short ix = 0 ; ix < m_sphere.nX() ; ++ix ) 
  { v[ix] = m_sphere.x2 ( ix ) ; }
  //
  // integrate them and to get new coefficients
  if   ( m_increasing ) { std::partial_sum ( v. begin() , v. end() ,  v. begin() ) ; }
  else                  { std::partial_sum ( v.rbegin() , v.rend() ,  v.rbegin() ) ; }
  //
  const double isum = 1.0 / 
    _spline_integral_ ( v , m_bspline.knots() , m_bspline.order() ) ;
  //
  for ( unsigned short ix = 0 ; ix < m_sphere.nX() ; ++ix ) 
  { 
    const bool updated = m_bspline.setPar ( ix , v [ix] * isum ) ; 
    update = updated || update ;
  }
  //
  return update ;
}
// ============================================================================
// ============================================================================
// COVEX ONLY  SPLINE 
// ============================================================================
// ============================================================================
/* constructor from the list of knots and the order 
 *  vector of parameters will be calculated automatically 
 *  @param points non-empty vector of poinst/knots 
 *  @param order  the order of splines 
 *  - vector of points is not requires to be ordered 
 *  - duplicated knots will be ignored
 *  - min/max value will be used as interval boundaries 
 */
// ============================================================================
Ostap::Math::ConvexOnlySpline::ConvexOnlySpline
( const std::vector<double>& points      ,
  const unsigned short       order       , 
  const bool                 convex      )
  : Ostap::Math::PositiveSpline ( points , order ) 
  , m_convex                    ( convex ) 
{
  if ( 2 > this->order() ) { Ostap::throwException
      ( "Degree of spline must be at least 2" , 
        "Ostap::Math::ConvexOnlySpline"       , 815 ) ; }
  updateCoefficients () ;
}
// ============================================================================
/*  Constructor from the list of knots and list of parameters 
 *  The spline order will be calculated automatically 
 *  @param points non-empty vector of poinst/knots 
 *  @param pars   non-empty vector of parameters 
 *  - vector of points is not requires to be ordered 
 *  - duplicated knots will be ignored
 *  - min/max value will be used as interval boundaries 
 */
// ============================================================================
Ostap::Math::ConvexOnlySpline::ConvexOnlySpline
( const std::vector<double>& points     ,
  const std::vector<double>& pars       ,
  const bool                 convex     ) 
  : Ostap::Math::PositiveSpline ( points , pars ) 
  , m_convex                    ( convex    ) 
{
  if ( 2 > order() ) { Ostap::throwException
      ( "Degree of spline must be at least 2" , 
        "Ostap::Math::ConvexOnlySpline"       , 815 ) ; }
  updateCoefficients () ;
}
// ============================================================================
/*  Constructor for uniform binning 
 *  @param xmin   low  edge of spline interval 
 *  @param xmax   high edge of spline interval 
 *  @param inner  number of inner points in   (xmin,xmax) interval
 *  @param order  the degree of splline 
 */
// ============================================================================
Ostap::Math::ConvexOnlySpline::ConvexOnlySpline
( const double         xmin      ,  
  const double         xmax      , 
  const unsigned short inner     ,   // number of inner points 
  const unsigned short order     , 
  const bool           convex    ) 
  : Ostap::Math::PositiveSpline ( xmin , xmax  , inner , order ) 
  , m_convex                    ( convex    ) 
{
  if ( this->order() < 2 ) { Ostap::throwException
      ( "Degree of spline must be at least 2" , 
        "Ostap::Math::ConvexOnlySpline"       , 815 ) ; }
  updateCoefficients () ;
}
// ============================================================================
// constructor from the basic spline 
// ============================================================================
Ostap::Math::ConvexOnlySpline::ConvexOnlySpline
( const Ostap::Math::PositiveSpline& spline     , 
  const bool                         convex     ) 
  : Ostap::Math::PositiveSpline ( spline        ) 
  , m_convex                    ( convex    ) 
{
  if ( 2 > order() ) { Ostap::throwException
      ( "Degree of spline must be at least 2" , 
        "Ostap::Math::ConvexOnlySpline"       , 815 ) ; }
  updateCoefficients () ;
}
// ============================================================================
// constructor from the basic spline 
// ============================================================================
Ostap::Math::ConvexOnlySpline::ConvexOnlySpline
( const Ostap::Math::BSpline&        spline     , 
  const bool                         convex     ) 
  : Ostap::Math::PositiveSpline ( spline        ) 
  , m_convex                    ( convex    ) 
{
  if ( 2 > order() ) { Ostap::throwException
      ( "Degree of spline must be at least 2" , 
        "Ostap::Math::ConvexOnlySpline"       , 815 ) ; }
  updateCoefficients () ;
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::ConvexOnlySpline::~ConvexOnlySpline(){}
// ============================================================================
// update coefficients  
// ============================================================================
bool Ostap::Math::ConvexOnlySpline::updateCoefficients  () 
{
  //
  if  ( order() < 2 ) { return Ostap::Math::PositiveSpline::updateCoefficients() ; }
  //
  bool   update = false ;
  //
  // get sphere coefficients 
  std::vector<double>  v ( m_sphere.nX() ) ;
  const unsigned short vs = v.size()     ;
  const unsigned short o  = order() ;  
  //
  if ( !m_convex ) 
  {
    const std::array<double,2> a = { { m_sphere.x2(0) , m_sphere.x2(1) } };
    for ( unsigned short ix = 2 ; ix < vs ; ++ix ) { v[ix] = m_sphere.x2 ( ix ) ; }
    //
    // integrate them and to get new coefficients
    std::partial_sum ( v.  begin() + 2 , v.  end() , v.  begin() + 2 ) ; 
    //
    for ( unsigned int i = 1 ; i < v.size()  ; ++i ) 
    { v[i] = v[i-1] + v[i] * ( knot_i (  i + o  + 1 ) - knot_i ( i )  ) / o ; }
    //
    const double last = v.back() ;
    for ( unsigned int i = 0 ; i < v.size() ; ++i ) { v[i] = last - v[i] ; }
    //
    const double v1 = a[0] - v.front() ;
    const double v2 = a[1] - v.back()  ;
    //
    for ( unsigned int j = 0 ; j < v.size()  ; ++j ) 
    { 
      double vj = 0 ;
      for ( unsigned short i = j + 1 ; i < j + o + 1 ; ++i ) { vj += knot_i ( i ) ; }
      v[j] += v1 + vj * ( v2 - v1 ) / o ;
    }
  }
  else 
  {
    const std::array<double,3> a = { { m_sphere.x2(0) , 
                                       m_sphere.x2(1) , 
                                       m_sphere.x2(2) } };
    for ( unsigned short ix = 3 ; ix < vs ; ++ix ) { v[ix] = m_sphere.x2 ( ix ) ; }
    // integrate them and to get new coefficients
    std::partial_sum ( v.  begin() + 3 , v.  end() , v.  begin() + 3 ) ; 
    //
    for ( unsigned int i = 3 ; i < v.size()  ; ++i ) 
    { v[i] = v[i-1] + v[i] * ( knot_i (  i + o  + 1 ) - knot_i ( i )  ) / o ; }
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
    for ( unsigned int j = 0 ; j < v.size()  ; ++j ) 
    { 
      double v1 = 0 ;
      for ( unsigned short   i = j + 1 ; i < j + o + 1 ; ++i ) { v1 += knot_i ( i ) ; }
      double v2 = 0 ;
      for ( unsigned short   i = j + 1 ; i < j + o     ; ++i ) 
      { for ( unsigned short k = i + 1 ; k < j + o + 1 ; ++k ) 
        { v2 += knot_i ( i ) * knot_i ( k ) ; } }
      v[j] += c0 + c1 * v1 / o + 2 * c2 * v2  / ( o * ( o - 1 ) ) ;
    }
  }
  //
  // normalize it! 
  //
  const double isum = 1.0 / 
    _spline_integral_ ( v , m_bspline.knots() , m_bspline.order() ) ;
  //
  for ( unsigned short ix = 0 ; ix < m_sphere.nX() ; ++ix ) 
  { 
    const bool updated = m_bspline.setPar ( ix , v [ix] * isum ) ; 
    update = updated || update ;
  }
  //
  return update ;
}
// ============================================================================
// ============================================================================
// convex SPLINE 
// ============================================================================
// ============================================================================
/* constructor from the list of knots and the order 
 *  vector of parameters will be calculated automatically 
 *  @param points non-empty vector of poinst/knots 
 *  @param order  the order of splines 
 *  - vector of points is not requires to be ordered 
 *  - duplicated knots will be ignored
 *  - min/max value will be used as interval boundaries 
 */
// ============================================================================
Ostap::Math::ConvexSpline::ConvexSpline
( const std::vector<double>& points      ,
  const unsigned short       order       , 
  const bool                 increasing  ,
  const bool                 convex      )
  : Ostap::Math::MonotonicSpline ( points , order , increasing ) 
  , m_convex                      ( convex ) 
{
  if ( 2 > this->order() ) { Ostap::throwException
      ( "Degree of spline must be at least 2" , 
        "Ostap::Math::ConvexSpline"           , 816 ) ; }
  updateCoefficients () ;
}
// ============================================================================
/*  Constructor from the list of knots and list of parameters 
 *  The spline order will be calculated automatically 
 *  @param points non-empty vector of poinst/knots 
 *  @param pars   non-empty vector of parameters 
 *  - vector of points is not requires to be ordered 
 *  - duplicated knots will be ignored
 *  - min/max value will be used as interval boundaries 
 */
// ============================================================================
Ostap::Math::ConvexSpline::ConvexSpline
( const std::vector<double>& points     ,
  const std::vector<double>& pars       ,
  const bool                 increasing ,
  const bool                 convex     ) 
  : Ostap::Math::MonotonicSpline ( points , pars , increasing ) 
  , m_convex                      ( convex ) 
{
  if ( 2 > order() ) { Ostap::throwException
      ( "Degree of spline must be at least 2" , 
        "Ostap::Math::ConvexSpline"           , 816 ) ; }
  updateCoefficients () ;
}
// ============================================================================
/*  Constructor for uniform binning 
 *  @param xmin   low  edge of spline interval 
 *  @param xmax   high edge of spline interval 
 *  @param inner  number of inner points in   (xmin,xmax) interval
 *  @param order  the degree of splline 
 */
// ============================================================================
Ostap::Math::ConvexSpline::ConvexSpline
( const double            xmin       ,  
  const double            xmax       , 
  const unsigned short    inner      ,   // number of inner points 
  const unsigned short    order      , 
  const bool              increasing ,
  const bool              convex     ) 
  : Ostap::Math::MonotonicSpline ( xmin , xmax , inner , order , increasing ) 
  , m_convex                      ( convex ) 
{
  if ( 2 > this->order() ) { Ostap::throwException
      ( "Degree of spline must be at least 2" , 
        "Ostap::Math::ConvexSpline"           , 816 ) ; }
  updateCoefficients () ;
}
// ============================================================================
// constructor from positive spline 
// ============================================================================
Ostap::Math::ConvexSpline::ConvexSpline
( const PositiveSpline&   spline      , 
  const bool              increasing  ,
  const bool              convex      ) 
  : Ostap::Math::MonotonicSpline ( spline , increasing ) 
  , m_convex                      ( convex ) 
{
  if ( 2 > order() ) { Ostap::throwException
      ( "Degree of spline must be at least 2" , 
        "Ostap::Math::ConvexSpline"           , 816 ) ; }
  updateCoefficients () ;
}
// ============================================================================
// constructor from basic spline 
// ============================================================================
Ostap::Math::ConvexSpline::ConvexSpline
( const BSpline&          spline     , 
  const bool              increasing ,
  const bool              convex     ) 
  : Ostap::Math::MonotonicSpline ( spline , increasing ) 
  , m_convex                      ( convex ) 
{
  if ( 2 > order() ) { Ostap::throwException
      ( "Degree of spline must be at least 2" , 
        "Ostap::Math::ConvexSpline"           , 816 ) ; }
  updateCoefficients () ;
}
// ============================================================================
// constructor from monotonic spline 
// ============================================================================
Ostap::Math::ConvexSpline::ConvexSpline
( const MonotonicSpline& spline , 
  const bool              convex ) 
  : Ostap::Math::MonotonicSpline ( spline ) 
  , m_convex                      ( convex ) 
{
  if ( 2 > order() ) { Ostap::throwException
      ( "Degree of spline must be at least 2" , 
        "Ostap::Math::ConvexSpline"           , 816 ) ; }
  updateCoefficients () ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Math::ConvexSpline::~ConvexSpline(){}
// ============================================================================
// update coefficients  
// ============================================================================
bool Ostap::Math::ConvexSpline::updateCoefficients  () 
{
  //
  bool   update = false ;
  //
  // get sphere coefficients  (all but 0 ) : NOTE INDICES HERE! 
  std::vector<double> v ( m_sphere.nX() - 1 ) ;
  for ( unsigned short ix = 1 ; ix < m_sphere.nX() ; ++ix ) 
  { v[ix-1] = m_sphere.x2 ( ix ) * ( ix + 1 ) ; }
  //
  // integrate them and to get new coefficients
  if   ( m_convex ) { std::partial_sum ( v. begin() , v. end() ,  v. begin() ) ; }
  else              { std::partial_sum ( v.rbegin() , v.rend() ,  v.rbegin() ) ; }
  //

  // Actual algorithm: "Safe"
  // BSpline aux1 ( m_bspline.knots() , v  )     ;
  // BSpline aux2 = aux1.indefinite_integral ( m_sphere.x2(0) ) ;  
  // v    = aux2.pars()    ;
  
  // in place: 
  //   the second integration: 
  std::vector<double> v2 ( m_sphere.nX() ) ;
  v2[0] = m_sphere.x2(0)  ;
  const unsigned short o = order() ;
  for ( unsigned int i = 0 ; i < v.size() ; ++i ) 
  { v2[i+1] = v2[i] + v[i] * ( knot_i (  i + o + 1 ) - knot_i ( i + 1 )  ) / o ; }
  //
  // revert, if needed 
  if ( !m_increasing ) { std::reverse ( v2.begin() , v2.end () ) ; }
  
  //
  const double isum = 1.0 / 
    _spline_integral_ ( v2 , m_bspline.knots() , m_bspline.order() ) ;
  //
  for ( unsigned short ix = 0 ; ix < m_sphere.nX() ; ++ix ) 
  { 
    const bool updated = m_bspline.setPar ( ix , v2 [ix] * isum ) ; 
    update = updated || update ;
  }
  //
  return update ;
}


// ============================================================================
// 2D-objects 
// ============================================================================






// ============================================================================
// Generic2D-spline
// ============================================================================
Ostap::Math::BSpline2D::BSpline2D
( const Ostap::Math::BSpline& xspline ,
  const Ostap::Math::BSpline& yspline )
  : m_xspline ( xspline        )
  , m_yspline ( yspline        )
  , m_pars    ( xspline.npars() * yspline.npars() )
{
  for ( unsigned i = 0 ; i < m_xspline.npars() ; ++i ) { m_xspline.setPar ( i , 0 ) ; }
  for ( unsigned i = 0 ; i < m_yspline.npars() ; ++i ) { m_yspline.setPar ( i , 0 ) ; }
}
// ===========================================================================
// set k-parameter
// ===========================================================================
bool Ostap::Math::BSpline2D::setParameter
( const unsigned int k , const double value ) 
{
  if ( k >= m_pars.size() )    { return false  ; }
  const double v = m_pars[k] ;
  if ( s_equal ( v , value ) ) { return false ; }
  //
  m_pars[k] = value ;
  return true ;
}
// ===========================================================================
// make the calcualtions 
// ===========================================================================
double Ostap::Math::BSpline2D::calculate
( const std::vector<double>& fx , 
  const std::vector<double>& fy ) const 
{
  const unsigned int NX = fx.size() ;
  const unsigned int NY = fy.size() ;
  //
  double result = 0.0 ;
  for ( unsigned short ix = 0 ; ix < NX ; ++ix ) 
  {
    const double vx = fx[ix] ;
    for ( unsigned short iy = 0 ; iy < NY ; ++iy ) 
    {
      const double vy = fy[iy] ;
      const double p  = par ( ix , iy ) ;
      result += p *vx * vy ;  
    }
  }
  return result * ( m_xspline.order() + 1 ) * ( m_yspline.order() + 1 ) ;
}
// ===========================================================================
// get the value
// ============================================================================
double Ostap::Math::BSpline2D::evaluate 
  ( const double x , const double y ) const
{
  //
  if ( x < xmin() || y < ymin() || x > xmax() || y > ymax() ) { return 0 ; }
  //
  const double xarg =
    !s_equal ( x , xmax ()) ? x :
    Ostap::Math::next_double ( xmax() , -s_ulps ) ;
  //
  const double yarg =
    !s_equal ( y , ymax ()) ? y :
    Ostap::Math::next_double ( ymax() , -s_ulps ) ;
  //
  const unsigned short NX = m_xspline.npars() ;
  const unsigned short NY = m_yspline.npars() ;
  //
  std::vector<double>  fx ( NX , 0 ) ;
  std::vector<double>  fy ( NY , 0 ) ;
  //
  // fill x-cache
  for ( unsigned short ix = 0 ; ix < NX ; ++ix )
  {
    m_xspline.setPar ( ix , 1.0 ) ;
    //
    double resx  = m_xspline ( xarg ) ;
    if ( 0 < resx )
    {
      const double ti  = knot ( m_xspline.knots() , ix                         ) ;
      const double tip = knot ( m_xspline.knots() , ix + m_xspline.order() + 1 ) ;
      resx /= ( tip - ti ) ;
    }
    //
    fx[ix] = resx ;
    //
    m_xspline.setPar ( ix , 0.0 ) ;
  }
  //
  // fill y-cache
  for ( unsigned short iy = 0 ; iy < NY ; ++iy )
  {
    m_yspline.setPar ( iy , 1.0 ) ;
    //
    double resy  = m_yspline ( yarg ) ;
    if ( 0 < resy )
    {
      const double ti  = knot ( m_yspline.knots() , iy                         ) ;
      const double tip = knot ( m_yspline.knots() , iy + m_yspline.order() + 1 ) ;
      resy /= ( tip - ti ) ;
    }
    //
    fy[iy] = resy ;
    //
    m_yspline.setPar ( iy , 0.0 ) ;
  }
  //
  return calculate ( fx , fy ) ;
}
// ============================================================================
/*  get the integral over 2D-region
 *  @param xlow  low  edge in x
 *  @param xhigh high edge in x
 *  @param ylow  low  edge in y
 *  @param yhigh high edge in y
 */
// ============================================================================
double Ostap::Math::BSpline2D::integral
( const double xlow  ,
  const double xhigh ,
  const double ylow  ,
  const double yhigh ) const
{
  //
  if      ( xhigh < xlow ) { return - integral ( xhigh , xlow  , ylow  , yhigh ) ; }
  else if ( yhigh < ylow ) { return - integral ( xlow  , xhigh , yhigh , ylow  ) ; }
  // boundaries
  else if ( xhigh < xmin () || yhigh < ymin () ) { return 0 ; }
  else if ( xlow  > xmax () || ylow  > ymax () ) { return 0 ; }
  else if ( s_equal ( xlow , xhigh ) || s_equal ( ylow , yhigh ) ) { return 0 ; }
  //
  else if ( s_equal ( xlow  , xmin() ) &&
            s_equal ( xhigh , xmax() ) &&
            s_equal ( ylow  , ymin() ) &&
            s_equal ( yhigh , ymax() ) ) { return integral() ; }
  // adjust
  else if ( xlow  < xmin () ) { return integral ( xmin() , xhigh   , ylow   , yhigh  ) ; }
  else if ( xhigh > xmax () ) { return integral ( xlow   , xmax () , ylow   , yhigh  ) ; }
  else if ( ylow  < ymin () ) { return integral ( xlow   , xhigh   , ymin() , yhigh  ) ; }
  else if ( yhigh > ymax () ) { return integral ( xlow   , xhigh   , ylow   , ymax() ) ; }
  //
  //
  const double xarg =
    !s_equal ( xhigh , xmax ()) ? xhigh :
    Ostap::Math::next_double ( xmax() , -s_ulps ) ;
  //
  const double yarg =
    !s_equal ( yhigh , ymax ()) ? yhigh :
    Ostap::Math::next_double ( ymax() , -s_ulps ) ;
  //
  const unsigned short NX = m_xspline.npars() ;
  const unsigned short NY = m_yspline.npars() ;
  //
  std::vector<double>  fx ( NX , 0 ) ;
  std::vector<double>  fy ( NY , 0 ) ;
  //
  // fill x-cache
  for ( unsigned short ix = 0 ; ix < NX ; ++ix )
  {
    m_xspline.setPar ( ix , 1.0 ) ;
    //
    double resx  = m_xspline.integral ( xlow , xarg ) ;
    if ( 0 < resx )
    {
      const double ti  = knot ( m_xspline.knots() , ix                         ) ;
      const double tip = knot ( m_xspline.knots() , ix + m_xspline.order() + 1 ) ;
      resx /= ( tip - ti ) ;
    }
    //
    fx[ix] = resx ;
    //
    m_xspline.setPar ( ix , 0.0 ) ;
  }
  //
  // fill y-cache
  for ( unsigned short iy = 0 ; iy < NY ; ++iy )
  {
    m_yspline.setPar ( iy , 1.0 ) ;
    //
    double resy  = m_yspline.integral ( ylow ,  yarg ) ;
    if ( 0 < resy )
    {
      const double ti  = knot ( m_yspline.knots() , iy                         ) ;
      const double tip = knot ( m_yspline.knots() , iy + m_yspline.order() + 1 ) ;
      resy /= ( tip - ti ) ;
    }
    //
    fy[iy] = resy ;
    //
    m_yspline.setPar ( iy , 0.0 ) ;
  }
  //
  return calculate ( fx , fy ) ;
}
// ============================================================================
/*  get the integral over Y  for given X
 *  @param x  (INPUT) x-value
 *  @param ylow  low  edge in y
 *  @param yhigh high edge in y
 */
// ============================================================================
double Ostap::Math::BSpline2D::integrateY
( const double x    ,
  const double ylow , const double yhigh ) const
{
  //
  if      ( x < xmin() || x > xmax()           ) { return 0 ; }
  else if ( s_equal ( ylow  , ymin() ) &&
            s_equal ( yhigh , ymax() ) ) { return integrateY ( x ) ; }
  else if ( yhigh <  ylow ) { return - integrateY ( x , yhigh , ylow ) ; }
  else if ( s_equal ( ylow , yhigh )           ) { return 0 ; }
  else if ( yhigh <  ymin() ||  ylow >  ymax() ) { return 0 ; }
  else if ( ylow  <  ymin() ) { return integrateY ( x , ymin() , yhigh  ) ; }
  else if ( yhigh >  ymax() ) { return integrateY ( x , ylow   , ymax() ) ; }
  //
  const double xarg =
    !s_equal ( x     , xmax ()) ? x :
    Ostap::Math::next_double ( xmax() , -s_ulps ) ;
  //
  const double yarg =
    !s_equal ( yhigh , ymax ()) ? yhigh :
    Ostap::Math::next_double ( ymax() , -s_ulps ) ;
  //
  const unsigned short NX = m_xspline.npars() ;
  const unsigned short NY = m_yspline.npars() ;
  //
  std::vector<double>  fx ( NX , 0 ) ;
  std::vector<double>  fy ( NY , 0 ) ;
  //
  // fill x-cache
  for ( unsigned short ix = 0 ; ix < NX ; ++ix )
  {
    m_xspline.setPar ( ix , 1.0 ) ;
    //
    double resx  = m_xspline ( xarg ) ;
    if ( 0 < resx )
    {
      const double ti  = knot ( m_xspline.knots() , ix                         ) ;
      const double tip = knot ( m_xspline.knots() , ix + m_xspline.order() + 1 ) ;
      resx /= ( tip - ti ) ;
    }
    //
    fx[ix] = resx ;
    //
    m_xspline.setPar ( ix , 0.0 ) ;
  }
  //
  // fill y-cache
  for ( unsigned short iy = 0 ; iy < NY ; ++iy )
  {
    m_yspline.setPar ( iy , 1.0 ) ;
    //
    double resy  = m_yspline.integral ( ylow ,  yarg ) ;
    if ( 0 < resy )
    {
      const double ti  = knot ( m_yspline.knots() , iy                         ) ;
      const double tip = knot ( m_yspline.knots() , iy + m_yspline.order() + 1 ) ;
      resy /= ( tip - ti ) ;
    }
    //
    fy[iy] = resy ;
    //
    m_yspline.setPar ( iy , 0.0 ) ;
  }
  //
  return calculate ( fx , fy ) ;
}
// ============================================================================
/*  get the integral over X for the given Y
 *  @param y  (INPUT) y-value
 *  @param xlow  low  eadge in x
 *  @param xhigh high edge in x
 */
// ============================================================================
double Ostap::Math::BSpline2D::integrateX
( const double y    ,
  const double xlow , const double xhigh ) const
{
  //
  if      ( y < ymin() || y > ymax()           ) { return 0 ; }
  else if ( s_equal ( xlow  , xmin() ) &&
            s_equal ( xhigh , xmax() ) ) { return integrateX ( y ) ; }
  else if ( xhigh <  xlow ) { return - integrateX ( y , xhigh , xlow ) ; }
  else if ( s_equal ( xlow , xhigh )           ) { return 0 ; }
  else if ( xhigh <= xmin() ||  xlow > xmax()  ) { return 0 ; }
  else if ( xlow  <  xmin() ) { return integrateX ( y , xmin() , xhigh  ) ; }
  else if ( xhigh >  xmax() ) { return integrateX ( y , xlow   , xmax() ) ; }
  //
  const double xarg =
    !s_equal ( xhigh , xmax ()) ? xhigh :
    Ostap::Math::next_double ( xmax() , -s_ulps ) ;
  //
  const double yarg =
    !s_equal ( y     , ymax ()) ? y     :
    Ostap::Math::next_double ( ymax() , -s_ulps ) ;
  //
  const unsigned short NX = m_xspline.npars() ;
  const unsigned short NY = m_yspline.npars() ;
  //
  std::vector<double>  fx ( NX , 0 ) ;
  std::vector<double>  fy ( NY , 0 ) ;
  //
  // fill x-cache
  for ( unsigned short ix = 0 ; ix < NX ; ++ix )
  {
    m_xspline.setPar ( ix , 1.0 ) ;
    //
    double resx  = m_xspline.integral ( xlow , xarg ) ;
    if ( 0 < resx )
    {
      const double ti  = knot ( m_xspline.knots() , ix                         ) ;
      const double tip = knot ( m_xspline.knots() , ix + m_xspline.order() + 1 ) ;
      resx /= ( tip - ti ) ;
    }
    //
    fx [ix] = resx ;
    //
    m_xspline.setPar ( ix , 0.0 ) ;
  }
  //
  // fill y-cache
  for ( unsigned short iy = 0 ; iy < NY ; ++iy )
  {
    m_yspline.setPar ( iy , 1.0 ) ;
    //
    double resy  = m_yspline( yarg ) ;
    if ( 0 < resy )
    {
      const double ti  = knot ( m_yspline.knots() , iy                         ) ;
      const double tip = knot ( m_yspline.knots() , iy + m_yspline.order() + 1 ) ;
      resy /= ( tip - ti ) ;
    }
    //
    fy [iy] = resy ;
    //
    m_yspline.setPar ( iy , 0.0 ) ;
  }
  //
  return calculate ( fx , fy ) ;
}
// ============================================================================
/*  get the integral over 2D-region
 *  \f[ x_{min}<x<x_{max}, y_{min}<y<y_{max}\f]
 *  @param xlow  low  edge in x
 *  @param xhigh high edge in x
 *  @param ylow  low  edge in y
 *  @param yhigh high edge in y
 */
// ============================================================================
double Ostap::Math::BSpline2D::integral   () const 
{ return std::accumulate ( m_pars.begin() , m_pars.end() ,  0.0 ) ; }
// ============================================================================
/*  get the integral over X  for given Y
 *  @param x  (INPUT) x-value
 */
// ============================================================================
double Ostap::Math::BSpline2D::integrateY ( const double x ) const
{
  //
  if      ( x < xmin() || x > xmax()           ) { return 0 ; }
  //
  const double xarg =
    !s_equal ( x     , xmax ()) ? x :
    Ostap::Math::next_double ( xmax() , -s_ulps ) ;
  //
  const unsigned short NX = m_xspline.npars() ;
  const unsigned short NY = m_yspline.npars() ;
  //
  std::vector<double>  fx ( NX ,  0 ) ;
  std::vector<double>  fy ( NY ,  1.0 / ( m_yspline.order() + 1 ) ) ;
  //
  // fill x-cache
  for ( unsigned short ix = 0 ; ix < NX ; ++ix )
  {
    m_xspline.setPar ( ix , 1.0 ) ;
    //
    double resx  = m_xspline ( xarg ) ;
    if ( 0 < resx )
    {
      const double ti  = knot ( m_xspline.knots() , ix                         ) ;
      const double tip = knot ( m_xspline.knots() , ix + m_xspline.order() + 1 ) ;
      resx /= ( tip - ti ) ;
    }
    //
    fx [ix] = resx ;
    //
    m_xspline.setPar ( ix , 0.0 ) ;
  }
  //
  return calculate ( fx , fy ) ;
}
// ============================================================================
/*  get the integral over X  for given Y
 *  @param y  (INPUT) y-value
 */
// ============================================================================
double Ostap::Math::BSpline2D::integrateX
( const double y    ) const
{
  //
  if      ( y < ymin() || y > ymax()           ) { return 0 ; }
  //
  const double yarg =
    !s_equal ( y     , ymax ()) ? y     :
    Ostap::Math::next_double ( ymax() , -s_ulps ) ;
  //
  const unsigned short NX = m_xspline.npars() ;
  const unsigned short NY = m_yspline.npars() ;
  //
  // fill x-cache
  std::vector<double>  fx ( NX ,  1.0 / ( m_xspline.order() +  1 ) ) ;
  std::vector<double>  fy ( NY ,  0.0 ) ;
  //
  // fill y-cache
  for ( unsigned short iy = 0 ; iy < NY ; ++iy )
  {
    m_yspline.setPar ( iy , 1.0 ) ;
    //
    double resy  = m_yspline( yarg ) ;
    if ( 0 < resy )
    {
      const double ti  = knot ( m_yspline.knots() , iy                         ) ;
      const double tip = knot ( m_yspline.knots() , iy + m_yspline.order() + 1 ) ;
      resy /= ( tip - ti ) ;
    }
    //
    fy[iy] = resy ;
    //
    m_yspline.setPar ( iy , 0.0 ) ;
  }
  //
  return calculate ( fx , fy ) ;
}
// ============================================================================
Ostap::Math::BSpline2D&
Ostap::Math::BSpline2D::operator += ( const double a ) 
{ 
  if   ( s_zero ( a ) ) { return *this ; }
  Ostap::Math::shift ( m_pars , a ) ;
  return *this ; 
}
// ============================================================================
Ostap::Math::BSpline2D&
Ostap::Math::BSpline2D::operator -= ( const double a ) 
{ 
  if   ( s_zero ( a ) ) { return *this ; }
  Ostap::Math::shift ( m_pars , -a ) ;
  return *this ; 
}
// ============================================================================
Ostap::Math::BSpline2D&
Ostap::Math::BSpline2D::operator *= ( const double a ) 
{ 
  if      ( s_equal ( a , 1 ) ) { return *this ; }
  else if ( s_zero  ( a     ) ) { std::fill ( m_pars.begin() , m_pars.end() , 0 ) ; }
  Ostap::Math::scale ( m_pars , a ) ;
  return *this ;
}
// ============================================================================
Ostap::Math::BSpline2D&
Ostap::Math::BSpline2D::operator /= ( const double a ) 
{ 
  if   ( s_equal ( a , 1 ) ) { return *this ; }
  Ostap::Math::scale ( m_pars , 1/a ) ;
  return *this ;
}
// ============================================================================
// negate it 
// ============================================================================
Ostap::Math::BSpline2D
Ostap::Math::BSpline2D::operator-() const 
{
  BSpline2D b ( *this ) ;
  Ostap::Math::negate ( b.m_pars ) ;
  return b ;
}
// ============================================================================
// Sum of BSpline polynomial and a constant
// ============================================================================
Ostap::Math::BSpline2D 
Ostap::Math::BSpline2D::__add__   ( const double value ) const 
{ return  (*this) + value ; }
// ============================================================================
// Sum of BSpline polynomial and a constant
// ============================================================================
Ostap::Math::BSpline2D 
Ostap::Math::BSpline2D::__radd__ ( const double value ) const 
{ return  (*this) + value ; }
// ============================================================================
// Subtract a constant from Benrstein polynomial
// ============================================================================
Ostap::Math::BSpline2D 
Ostap::Math::BSpline2D::__sub__   ( const double value ) const 
{ return  (*this) - value ; }
// ============================================================================
// Constant minus BSpline polynomial
// ============================================================================
Ostap::Math::BSpline2D 
Ostap::Math::BSpline2D::__rsub__  ( const double value ) const 
{ return  value - (*this) ; }
// ============================================================================
// Product of BSpline polynomial and a constant
// ============================================================================
Ostap::Math::BSpline2D 
Ostap::Math::BSpline2D::__mul__  ( const double value ) const 
{ return  (*this) * value ; }
// ============================================================================
// Product of BSpline polynomial and a constant
// ============================================================================
Ostap::Math::BSpline2D 
Ostap::Math::BSpline2D::__rmul__ ( const double value ) const 
{ return  (*this) * value ; }
// ============================================================================
// Divide Benrstein polynomial by a constant
// ============================================================================
Ostap::Math::BSpline2D 
Ostap::Math::BSpline2D::__div__  ( const double value ) const 
{ return  (*this) / value ; }
// ============================================================================
// Negate BSpline polynomial
// ============================================================================
Ostap::Math::BSpline2D 
Ostap::Math::BSpline2D::__neg__  () const { return  -(*this) ; }
// ============================================================================


// ============================================================================
// Symmetric 
// ============================================================================

// ============================================================================
// Generic2DSym-spline
// ============================================================================
Ostap::Math::BSpline2DSym::BSpline2DSym
( const Ostap::Math::BSpline& spline )
  : m_spline ( spline        )
  , m_pars   ( spline.npars() * ( spline.npars() + 1 ) / 2 )
{
  for ( unsigned i = 0 ; i < m_spline.npars() ; ++i ) { m_spline.setPar ( i , 0 ) ; }
}
// ===========================================================================
// set k-parameter
// ===========================================================================
bool Ostap::Math::BSpline2DSym::setParameter
( const unsigned int k , const double value ) 
{
  if ( k >= m_pars.size() )    { return false  ; }
  const double v = m_pars[k] ;
  if ( s_equal ( v , value ) ) { return false ; }
  //
  m_pars[k] = value ;
  return true ;
}
// ===========================================================================
// make the calcualtions 
// ===========================================================================
double Ostap::Math::BSpline2DSym::calculate
( const std::vector<double>& fx , 
  const std::vector<double>& fy ) const 
{
  //
  const unsigned int NX = fx.size() ;
  // const unsigned int NY = fy.size() ;
  //
  double result = 0.0 ;
  for  ( unsigned short ix = 0 ; ix < NX  ; ++ix )
  {
    result   += par ( ix , ix ) * fx[ix] * fy[ix] ;
    for  ( unsigned short iy = 0 ; iy < ix ; ++iy )
    { result += par ( ix , iy ) * ( fx[ix] * fy[iy] + fx[iy] * fy[ix] ) ; } 
  }
  //
  const unsigned int scale = ( m_spline.order() + 1 ) ;
  //
  return result * scale * scale ;
}
// ===========================================================================
// get the value
// ============================================================================
double Ostap::Math::BSpline2DSym::evaluate 
( const double x , const double y ) const
{
  //
  if ( x < xmin() || y < ymin() || x > xmax() || y > ymax() ) { return 0 ; }
  //
  const double xarg =
    !s_equal ( x , xmax ()) ? x :
    Ostap::Math::next_double ( xmax() , -s_ulps ) ;
  //
  const double yarg =
    !s_equal ( y , ymax ()) ? y :
    Ostap::Math::next_double ( ymax() , -s_ulps ) ;
  //
  const unsigned short NX = m_spline.npars() ;
  const unsigned short NY = m_spline.npars() ;
  //
  std::vector<double>  fx ( NX , 0 ) ;
  std::vector<double>  fy ( NY , 0 ) ;
  //
  // fill x-cache
  for ( unsigned short ix = 0 ; ix < NX ; ++ix )
  {
    m_spline.setPar ( ix , 1.0 ) ;
    //
    double resx  = m_spline ( xarg ) ;
    if ( 0 < resx )
    {
      const double ti  = knot ( m_spline.knots() , ix                        ) ;
      const double tip = knot ( m_spline.knots() , ix + m_spline.order() + 1 ) ;
      resx /= ( tip - ti ) ;
    }
    //
    fx[ix] = resx ;
    //
    m_spline.setPar ( ix , 0.0 ) ;
  }
  //
  // fill y-cache
  for ( unsigned short iy = 0 ; iy < NY ; ++iy )
  {
    m_spline.setPar ( iy , 1.0 ) ;
    //
    double resy  = m_spline ( yarg ) ;
    if ( 0 < resy )
    {
      const double ti  = knot ( m_spline.knots() , iy                        ) ;
      const double tip = knot ( m_spline.knots() , iy + m_spline.order() + 1 ) ;
      resy /= ( tip - ti ) ;
    }
    //
    fy[iy] = resy ;
    //
    m_spline.setPar ( iy , 0.0 ) ;
  }
  //
  return calculate ( fx , fy ) ;
}
// ============================================================================
/*  get the integral over 2D-region
 *  @param xlow  low  edge in x
 *  @param xhigh high edge in x
 *  @param ylow  low  edge in y
 *  @param yhigh high edge in y
 */
// ============================================================================
double Ostap::Math::BSpline2DSym::integral
( const double xlow  ,
  const double xhigh ,
  const double ylow  ,
  const double yhigh ) const
{
  //
  if      ( xhigh < xlow ) { return - integral ( xhigh , xlow  , ylow  , yhigh ) ; }
  else if ( yhigh < ylow ) { return - integral ( xlow  , xhigh , yhigh , ylow  ) ; }
  // boundaries
  else if ( xhigh < xmin () || yhigh < ymin () ) { return 0 ; }
  else if ( xlow  > xmax () || ylow  > ymax () ) { return 0 ; }
  else if ( s_equal ( xlow , xhigh ) || s_equal ( ylow , yhigh ) ) { return 0 ; }
  //
  else if ( s_equal ( xlow  , xmin() ) &&
            s_equal ( xhigh , xmax() ) &&
            s_equal ( ylow  , ymin() ) &&
            s_equal ( yhigh , ymax() ) ) { return integral() ; }
  // adjust
  else if ( xlow  < xmin () ) { return integral ( xmin() , xhigh   , ylow   , yhigh  ) ; }
  else if ( xhigh > xmax () ) { return integral ( xlow   , xmax () , ylow   , yhigh  ) ; }
  else if ( ylow  < ymin () ) { return integral ( xlow   , xhigh   , ymin() , yhigh  ) ; }
  else if ( yhigh > ymax () ) { return integral ( xlow   , xhigh   , ylow   , ymax() ) ; }
  //
  //
  const double xarg =
    !s_equal ( xhigh , xmax ()) ? xhigh :
    Ostap::Math::next_double ( xmax() , -s_ulps ) ;
  //
  const double yarg =
    !s_equal ( yhigh , ymax ()) ? yhigh :
    Ostap::Math::next_double ( ymax() , -s_ulps ) ;
  //
  const unsigned short NX = m_spline.npars() ;
  const unsigned short NY = m_spline.npars() ;
  //
  std::vector<double>  fx ( NX , 0 ) ;
  std::vector<double>  fy ( NY , 0 ) ;
  //
  // fill x-cache
  for ( unsigned short ix = 0 ; ix < NX ; ++ix )
  {
    m_spline.setPar ( ix , 1.0 ) ;
    //
    double resx  = m_spline.integral ( xlow , xarg ) ;
    if ( 0 < resx )
    {
      const double ti  = knot ( m_spline.knots() , ix                        ) ;
      const double tip = knot ( m_spline.knots() , ix + m_spline.order() + 1 ) ;
      resx /= ( tip - ti ) ;
    }
    //
    fx[ix] = resx ;
    //
    m_spline.setPar ( ix , 0.0 ) ;
  }
  //
  // fill y-cache
  for ( unsigned short iy = 0 ; iy < NY ; ++iy )
  {
    m_spline.setPar ( iy , 1.0 ) ;
    //
    double resy  = m_spline.integral ( ylow ,  yarg ) ;
    if ( 0 < resy )
    {
      const double ti  = knot ( m_spline.knots() , iy                        ) ;
      const double tip = knot ( m_spline.knots() , iy + m_spline.order() + 1 ) ;
      resy /= ( tip - ti ) ;
    }
    //
    fy[iy] = resy ;
    //
    m_spline.setPar ( iy , 0.0 ) ;
  }
  //
  return calculate ( fx , fy ) ;
}
// ============================================================================
/*  get the integral over Y  for given X
 *  @param x  (INPUT) x-value
 *  @param ylow  low  edge in y
 *  @param yhigh high edge in y
 */
// ============================================================================
double Ostap::Math::BSpline2DSym::integrateY
( const double x    ,
  const double ylow , const double yhigh ) const
{
  //
  if      ( x < xmin() || x > xmax()           ) { return 0 ; }
  else if ( s_equal ( ylow  , ymin() ) &&
            s_equal ( yhigh , ymax() ) ) { return integrateY ( x ) ; }
  else if ( yhigh <  ylow ) { return - integrateY ( x , yhigh , ylow ) ; }
  else if ( s_equal ( ylow , yhigh )           ) { return 0 ; }
  else if ( yhigh <  ymin() ||  ylow >  ymax() ) { return 0 ; }
  else if ( ylow  <  ymin() ) { return integrateY ( x , ymin() , yhigh  ) ; }
  else if ( yhigh >  ymax() ) { return integrateY ( x , ylow   , ymax() ) ; }
  //
  const double xarg =
    !s_equal ( x     , xmax ()) ? x :
    Ostap::Math::next_double ( xmax() , -s_ulps ) ;
  //
  const double yarg =
    !s_equal ( yhigh , ymax ()) ? yhigh :
    Ostap::Math::next_double ( ymax() , -s_ulps ) ;
  //
  const unsigned short NX = m_spline.npars() ;
  const unsigned short NY = m_spline.npars() ;
  //
  std::vector<double>  fx ( NX , 0 ) ;
  std::vector<double>  fy ( NY , 0 ) ;
  //
  // fill x-cache
  for ( unsigned short ix = 0 ; ix < NX ; ++ix )
  {
    m_spline.setPar ( ix , 1.0 ) ;
    //
    double resx  = m_spline ( xarg ) ;
    if ( 0 < resx )
    {
      const double ti  = knot ( m_spline.knots() , ix                        ) ;
      const double tip = knot ( m_spline.knots() , ix + m_spline.order() + 1 ) ;
      resx /= ( tip - ti ) ;
    }
    //
    fx[ix] = resx ;
    //
    m_spline.setPar ( ix , 0.0 ) ;
  }
  //
  // fill y-cache
  for ( unsigned short iy = 0 ; iy < NY ; ++iy )
  {
    m_spline.setPar ( iy , 1.0 ) ;
    //
    double resy  = m_spline.integral ( ylow ,  yarg ) ;
    if ( 0 < resy )
    {
      const double ti  = knot ( m_spline.knots() , iy                        ) ;
      const double tip = knot ( m_spline.knots() , iy + m_spline.order() + 1 ) ;
      resy /= ( tip - ti ) ;
    }
    //
    fy[iy] = resy ;
    //
    m_spline.setPar ( iy , 0.0 ) ;
  }
  //
  return calculate ( fx , fy ) ;
}

// ============================================================================
/*  get the integral over X for the given Y
 *  @param y  (INPUT) y-value
 *  @param xlow  low  eadge in x
 *  @param xhigh high edge in x
 */
// ============================================================================
double Ostap::Math::BSpline2DSym::integrateX
( const double y    ,
  const double xlow , const double xhigh ) const
{ return  integrateY ( y , xlow , xhigh ) ; }
// ============================================================================
/*  get the integral over 2D-region
 *  \f[ x_{min}<x<x_{max}, y_{min}<y<y_{max}\f]
 *  @param xlow  low  edge in x
 *  @param xhigh high edge in x
 *  @param ylow  low  edge in y
 *  @param yhigh high edge in y
 */
// ============================================================================
double Ostap::Math::BSpline2DSym::integral   () const 
{ 
  const unsigned short  NX = m_spline.npars() ;
  double result = 0 ;
  for  ( unsigned short ix = 0 ; ix < NX ; ++ix )
  {
    result   += par ( ix , ix ) ;
    for  ( unsigned short iy = 0 ; iy < ix ; ++iy )
    { result += 2 * par ( ix , iy ) ; } 
  }
  //
  return result ;
}
// ============================================================================
/*  get the integral over X  for given Y
 *  @param x  (INPUT) x-value
 */
// ============================================================================
double Ostap::Math::BSpline2DSym::integrateY ( const double x ) const
{
  //
  if      ( x < xmin() || x > xmax()           ) { return 0 ; }
  //
  const double xarg =
    !s_equal ( x     , xmax ()) ? x :
    Ostap::Math::next_double ( xmax() , -s_ulps ) ;
  //
  const unsigned short NX = m_spline.npars() ;
  const unsigned short NY = m_spline.npars() ;
  //
  std::vector<double>  fx ( NX ,  0 ) ;
  std::vector<double>  fy ( NY ,  1.0 / ( m_spline.order() + 1 ) ) ;
  //
  // fill x-cache
  for ( unsigned short ix = 0 ; ix < NX ; ++ix )
  {
    m_spline.setPar ( ix , 1.0 ) ;
    //
    double resx  = m_spline ( xarg ) ;
    if ( 0 < resx )
    {
      const double ti  = knot ( m_spline.knots() , ix                        ) ;
      const double tip = knot ( m_spline.knots() , ix + m_spline.order() + 1 ) ;
      resx /= ( tip - ti ) ;
    }
    //
    fx [ix] = resx ;
    //
    m_spline.setPar ( ix , 0.0 ) ;
  }
  //
  return calculate ( fx , fy ) ;
}
// ============================================================================
/*  get the integral over X  for given Y
 *  @param y  (INPUT) y-value
 */
// ============================================================================
double Ostap::Math::BSpline2DSym::integrateX
( const double y    ) const { return integrateY (  y ) ; }

// ============================================================================
Ostap::Math::BSpline2DSym&
Ostap::Math::BSpline2DSym::operator += ( const double a ) 
{ 
  if   ( s_zero ( a ) ) { return *this ; }
  Ostap::Math::shift ( m_pars , a ) ;
  return *this ; 
}
// ============================================================================
Ostap::Math::BSpline2DSym&
Ostap::Math::BSpline2DSym::operator -= ( const double a ) 
{ 
  if   ( s_zero ( a ) ) { return *this ; }
  Ostap::Math::shift ( m_pars , -a ) ;
  return *this ; 
}
// ============================================================================
Ostap::Math::BSpline2DSym&
Ostap::Math::BSpline2DSym::operator *= ( const double a ) 
{ 
  if      ( s_equal ( a , 1 ) ) { return *this ; }
  else if ( s_zero  ( a     ) ) { std::fill ( m_pars.begin() , m_pars.end() , 0 ) ; }
  Ostap::Math::scale ( m_pars , a ) ;
  return *this ;
}
// ============================================================================
Ostap::Math::BSpline2DSym&
Ostap::Math::BSpline2DSym::operator /= ( const double a ) 
{ 
  if   ( s_equal ( a , 1 ) ) { return *this ; }
  Ostap::Math::scale ( m_pars , 1/a ) ;
  return *this ;
}
// ============================================================================
// negate it 
// ============================================================================
Ostap::Math::BSpline2DSym
Ostap::Math::BSpline2DSym::operator-() const 
{
  BSpline2DSym b ( *this ) ;
  Ostap::Math::negate ( b.m_pars ) ;
  return b ;
}
// ============================================================================
// Sum of BSpline2D polynomial and a constant
// ============================================================================
Ostap::Math::BSpline2DSym 
Ostap::Math::BSpline2DSym::__add__   ( const double value ) const 
{ return  (*this) + value ; }
// ============================================================================
// Sum of BSpline2D polynomial and a constant
// ============================================================================
Ostap::Math::BSpline2DSym 
Ostap::Math::BSpline2DSym::__radd__ ( const double value ) const 
{ return  (*this) + value ; }
// ============================================================================
// Subtract a constant from Benrstein polynomial
// ============================================================================
Ostap::Math::BSpline2DSym 
Ostap::Math::BSpline2DSym::__sub__   ( const double value ) const 
{ return  (*this) - value ; }
// ============================================================================
// Constant minus BSpline2D polynomial
// ============================================================================
Ostap::Math::BSpline2DSym 
Ostap::Math::BSpline2DSym::__rsub__  ( const double value ) const 
{ return  value - (*this) ; }
// ============================================================================
// Product of BSpline2D polynomial and a constant
// ============================================================================
Ostap::Math::BSpline2DSym 
Ostap::Math::BSpline2DSym::__mul__  ( const double value ) const 
{ return  (*this) * value ; }
// ============================================================================
// Product of BSpline2D polynomial and a constant
// ============================================================================
Ostap::Math::BSpline2DSym 
Ostap::Math::BSpline2DSym::__rmul__ ( const double value ) const 
{ return  (*this) * value ; }
// ============================================================================
// Divide Benrstein polynomial by a constant
// ============================================================================
Ostap::Math::BSpline2DSym 
Ostap::Math::BSpline2DSym::__div__  ( const double value ) const 
{ return  (*this) / value ; }
// ============================================================================
// Negate BSpline2D polynomial
// ============================================================================
Ostap::Math::BSpline2DSym 
Ostap::Math::BSpline2DSym::__neg__  () const { return  -(*this) ; }
// ============================================================================

// ============================================================================
// Positive 2D spline 
// ============================================================================
Ostap::Math::PositiveSpline2D::PositiveSpline2D 
( const Ostap::Math::BSpline& xspline ,
  const Ostap::Math::BSpline& yspline ) 
  : m_spline  ( xspline , yspline ) 
  , m_sphere  ( xspline.npars() * yspline.npars() - 1 ) 
{
  updateSpline() ;
}
// =============================================================================
// update spline coefficients
// =============================================================================
bool Ostap::Math::PositiveSpline2D::updateSpline()  
{
  //
  bool update = false ;
  for ( unsigned int ix = 0 ; ix < m_sphere.nX() ; ++ix )
  {
    const bool updated = m_spline.setPar ( ix , m_sphere.x2 ( ix ) ) ;
    update = updated || update ;
  }
  //
  return update ;
}
// ============================================================================
// Positive symmetric 2D spline 
// ============================================================================
Ostap::Math::PositiveSpline2DSym::PositiveSpline2DSym 
( const Ostap::Math::BSpline&  spline )
  : m_spline ( spline ) 
  , m_sphere ( spline.npars() * ( spline.npars() + 1 ) / 2 - 1  )
{
  updateSpline() ;
}
// =============================================================================
// update spline coefficients
// =============================================================================
bool Ostap::Math::PositiveSpline2DSym::updateSpline()  
{
  //
  bool update = false ;
  for ( unsigned int ix = 0 ; ix < m_sphere.nX() ; ++ix )
  {
    const bool updated = m_spline.setPar ( ix , m_sphere.x2 ( ix ) ) ;
    update = updated || update ;
  }
  // 
  if ( update ) { m_spline /= m_spline.integral() ; }
  //
  return update ;
}
// ============================================================================




// ============================================================================
// Berstein polynomials 
// ============================================================================
#include "Ostap/Bernstein.h"
// ============================================================================
namespace 
{
  // =========================================================================
  /// calculate the convex hull
  template <class COMPARE> 
  inline Ostap::Math::BSpline _convex_hull_ ( const Ostap::Math::Bernstein& p , COMPARE cmp )
  {
    const std::vector<double>& bpars = p.pars() ;
    //
    std::vector<double> knots ;
    //
    if ( 2 >= bpars.size() ) // special case 
    {
      knots.push_back ( p.xmin ()        ) ;
      knots.push_back ( p.xmax ()        ) ;
      return Ostap::Math::BSpline ( knots , bpars ) ;
    }
    //
    std::vector<double> pars  ;
    const unsigned short N = bpars.size() ;
    //
    knots . push_back ( 0        ) ;
    pars  . push_back ( bpars[0] ) ;
    unsigned short icurr = 1 ;
    while (  icurr < N ) 
    {
      const double pl = pars .back() ;
      const double kl = knots.back() ;
      //
      double xi = double ( icurr ) / ( N  -  1 ) ;
      double is = ( bpars[ icurr ] - pl ) / ( xi - kl ) ;
      //
      for ( unsigned short j = icurr + 1 ; j < N ; ++j ) 
      {
        //
        const double xj = double ( j ) / ( N  -  1 ) ;
        const double js = ( bpars[j] - pl ) / ( xj - kl ) ;
        //
        if ( cmp ( js , is ) ) { is = js ; xi = xj ; icurr = j ; }
      }
      knots.push_back ( xi           ) ;
      pars .push_back ( bpars[icurr] ) ;
      ++icurr ;
    }
    //
    std::transform ( knots.begin() , knots.end() ,  
                     knots.begin() , [&p]( const double t ) { return p.x(t) ; }  );
    //
    return Ostap::Math::BSpline ( knots ,  pars ) ;
  }
  // ======================================================================
}
// ========================================================================
/*  calculate the convex hull for Bernstein Polynomial 
 *  @param p  bernstein Polynomial
 *  @return   the spline object that represents upper convex hull 
 */
// ========================================================================
Ostap::Math::BSpline
Ostap::Math::upper_convex_hull   ( const Ostap::Math::Bernstein& p ) 
{ return _convex_hull_ (  p , Ostap::Math::GreaterOrEqual<double>() ) ; }
// ========================================================================
/** calculate the convex hull for Bernstein Polynomial 
 *  @param p  bernstein Polynomial
 *  @return   the spline object that represents lower convex hull 
 */
Ostap::Math::BSpline
Ostap::Math::lower_convex_hull ( const Ostap::Math::Bernstein& p ) 
{ return _convex_hull_ (  p , Ostap::Math::LessOrEqual<double>() ) ; }
// ========================================================================
/*  get control polygon  for Bernstein polynomial
 *  @param p  bernstein Polynomial
 *  @return   the spline object that represents the control polygon
 */
// ========================================================================
Ostap::Math::BSpline
Ostap::Math::control_polygon   ( const Ostap::Math::Bernstein& p ) 
{
  const std::vector<double>& pars = p.pars() ;
  //
  if ( 1 >= pars.size() ) 
  { return Ostap::Math::BSpline( {{ p.xmin() , p.xmax() }}  , pars ) ; }
  //
  const long double dx = ( p.xmax() - p.xmin() ) / ( pars.size() - 1 ) ;
  std::vector<double>        knots ( pars.size() ) ;
  const unsigned short N =   knots.size() ;
  for ( unsigned short i = 0 ; i < N ; ++i ) 
  { knots[i] = ( i * dx  + p.xmin() ) ; } ;
  return Ostap::Math::BSpline ( knots , pars ) ;
}
// ============================================================================
/*  get control polygon  basic spline 
 *  @param p  Basis spline 
 *  @return   the spline object that represents the control polygon
 */
// ============================================================================
Ostap::Math::BSpline
Ostap::Math::control_polygon   ( const Ostap::Math::BSpline& p ) 
{ return Ostap::Math::BSpline ( p.greville_abscissas() , p.pars() ) ; }
// ============================================================================
namespace 
{
  // ==========================================================================
  inline std::vector<double> 
  _crossing_points_1_ ( const Ostap::Math::BSpline& b ) 
  {
    std::vector<double> cps ; 
    //
    const std::vector<double>&  bpars = b.pars();
    cps.reserve( bpars.size() ) ;
    //
    const double p0 =  bpars.front()  ;
    if ( s_zero ( p0 ) ) {  cps.push_back( b.xmin() ) ; }
    //
    const double norm = b.norm()  ;
    const Ostap::Math::Tiny<double> tiny { norm } ;
    //
    const unsigned short N = bpars.size() ;
    //
    for ( unsigned short j = 1  ; j < N ;  ++j ) 
    {
      //
      const double pj = bpars[j  ] ;
      const double pi = bpars[j-1] ;
      //
      const double xj = b.greville_abscissa ( j ) ;
      if ( s_zero ( pj )|| s_equal ( pj + norm , norm ) ) 
      { cps.push_back ( xj ) ; continue ;   }
      //
      if ( s_zero ( pi ) || s_equal ( pi + norm , norm ) ) { continue ; }
      //
      const signed char sj = Ostap::Math::signum ( pj ) ;
      const signed char si = Ostap::Math::signum ( pi ) ;
      //
      if ( 0 > si * sj ) // there is root here! 
      {
        const double xi =  b.greville_abscissa ( j - 1 ) ;
        const double cp = ( xj * pi - xi * pj ) / ( pi - pj ) ;
        cps.push_back ( cp ) ;
      }
    }
    //
    return cps ;
  }
  // =================================================================================
  inline std::vector<double> 
  _crossing_points_2_ ( const Ostap::Math::BSpline& b )
  {
    std::vector<double> cps ; 
    //
    const std::vector<double>&  bpars = b.pars();
    cps.reserve( bpars.size() ) ;
    //
    const double norm = b.norm()  ;
    const Ostap::Math::Tiny<double> tiny { norm } ;
    //
    const unsigned short N = bpars.size() ;
    //
    // find first non-zero 
    std::vector<double>::const_iterator         _i1 =
      std::find_if_not ( bpars. cbegin () , bpars. cend () , tiny ) ;
    // 
    if ( bpars. cend () == _i1 ) { return b.greville_abscissas() ; }
    // find last non-zero 
    std::vector<double>::const_reverse_iterator _i2 =
      std::find_if_not ( bpars.crbegin () , bpars.crend () , tiny ) ;
    if ( bpars.crend () == _i2 ) { return b.greville_abscissas() ; }
    //
    const unsigned short i1 =                _i1 - bpars. cbegin () ;
    const unsigned short i2 =  ( N - 1 ) - ( _i2 - bpars.crbegin () ) ;
    //
    if ( 0 != i1 ) { cps.push_back ( b.xmin() ) ; }
    //
    const std::vector<double>::const_iterator begin = bpars.begin() + i1     ;
    const std::vector<double>::const_iterator end   = bpars.begin() + i2 + 1 ;
    //
    std::vector<double>::const_iterator i = begin ;
    std::vector<double>::const_iterator j = begin ;
    while ( j != end ) 
    {
      j = std::find_if_not ( i + 1 , end , tiny ) ;
      if ( end == j ) { break ; }
      //
      const double pj = *j ;
      const double pi = *i ;
      //
      const signed char sj = Ostap::Math::signum ( pj ) ;
      const signed char si = Ostap::Math::signum ( pi ) ;
      //
      // there is root here! 
      if ( si * sj < 0 ) 
      {
        const double xi = b.greville_abscissa ( i - bpars.begin() ) ;
        const double xj = b.greville_abscissa ( j - bpars.begin() ) ;
        const double cp = ( xj * pi - xi * pj ) / ( pi - pj ) ;
        cps.push_back ( cp ) ;
      }
      //
      i = j;
    }
    //
    if ( N-1 != i2 ) { cps.push_back ( b.xmax() ) ; }
    //
    return cps ;
  }
}
// ============================================================================
/*  get abscissas of crossing points of the control polygon with x-axis
 *  @param  b bernstein polynomial
 *  @return abscissas of crossing points of the control  polygon
 */
// ============================================================================
std::vector<double> 
Ostap::Math::crossing_points 
( const Ostap::Math::BSpline& b   , 
  const bool                  all )  
{  
  return ( all || 1>= b.degree()  ) ? 
    _crossing_points_1_ ( b ) :
    _crossing_points_2_ ( b ) ;
}
// ============================================================================
//  Here we'll use GSL
// ============================================================================
#include <gsl/gsl_linalg.h>
// ============================================================================
/* create the interpolation spline 
 *  @param xy (INPUT)   vector of data 
 *  @param bs (UPDATE) the spline 
 *  @return status code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Math::Interpolation::bspline 
( const std::vector<double>& x  ,
  const std::vector<double>& y  ,
  Ostap::Math::BSpline&      bs ) 
{
  // 
  if  ( x.size() != y.size() ) { return 100 ; }         // RETURN 100 
  const unsigned  short N = x.size() ;
  // mismatch for number of input parameters 
  if ( N != bs.npars()       ) { return 101 ; }         // RETURN 101 
  //
  std::vector< std::pair<double,double> > xy { N } ;
  for ( unsigned short i = 0 ; i < N ;   ++i ) 
  { xy[i] = std::make_pair ( x[i] , y[i] ) ; }
  return bspline ( xy , bs ) ;
}
// ============================================================================
/* create the interpolation spline 
 *  @param xy (INPUT)   vector of data 
 *  @param bs (UPPDATE) the spline 
 *  @return status code 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::Math::Interpolation::bspline 
( std::vector< std::pair<double,double> > xy ,
  Ostap::Math::BSpline&                   bs ) 
{  
  const unsigned  short N = xy.size() ;
  // mismatch for number of input parameters 
  if ( N != bs.npars() ) { return 110 ; }             // RETURN 110 
  //
  std::sort (  xy.begin() , xy.end() ) ;
  //
  gsl_matrix      * m = gsl_matrix_alloc ( N , N );
  // 
  const double xmin = bs.xmin () ;
  const double xmax = bs.xmax () ;
  //
  for ( unsigned short i = 0 ; i < N ; ++i ) 
  { for ( unsigned short j = 0 ; j < N ; ++j ) 
    {
      double xj =  xy[j].first ;
      if      ( s_equal ( xj , xmin ) ) 
      { xj = Ostap::Math::next_double ( xmin , +s_ulps ) ; }
      else if ( s_equal ( xj , xmax ) ) 
      { xj = Ostap::Math::next_float ( xmax , -s_ulps ) ; }
      //
      const double bij = bs.bspline ( i , xj ) ;
      if  ( i == j && s_zero ( bij ) ) 
      { gsl_matrix_free  ( m ) ; { return 111 ; } }  // RETURN 111 
      gsl_matrix_set (  m , j , i , bij  ) ; 
    } 
  }
  //
  gsl_vector      *x = gsl_vector_alloc      ( N ) ;
  for (  unsigned short i = 0 ; i < N ; ++i ) 
  { gsl_vector_set ( x , i , xy[i].second ) ;  }
  //
  gsl_permutation *p = gsl_permutation_alloc ( N  );
  //
  // make LU decomposition 
  int       signum = 0 ;
  const int e1     = gsl_linalg_LU_decomp ( m , p , &signum  );
  if ( e1 )
  {
    gsl_permutation_free ( p ) ;
    gsl_matrix_free      ( m ) ;
    gsl_vector_free      ( x ) ;
    return 120 + e1 ;                       // RETURN 120 + e 
  }
  //
  const int e2 = gsl_linalg_LU_svx ( m , p , x ) ;
  if ( e2 )
  {
    gsl_permutation_free ( p ) ;
    gsl_matrix_free      ( m ) ;
    gsl_vector_free      ( x ) ;
    return 130 + e2 ;                        // RETURN 130 + e
  }
  //
  for (  unsigned short i = 0 ; i < N ; ++i ) 
  { bs.setPar ( i , gsl_vector_get ( x , i ) ) ; }
  //
  gsl_permutation_free ( p ) ;
  gsl_matrix_free      ( m ) ;
  gsl_vector_free      ( x ) ;
  //
  return Ostap::StatusCode::SUCCESS ;
}
// ============================================================================
// The END 
// ============================================================================
