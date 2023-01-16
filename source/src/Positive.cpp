// ============================================================================
// Include files 
// ============================================================================
// STD&STL
// ============================================================================
#include <algorithm>
#include <numeric>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Hash.h"
#include "Ostap/NSphere.h"
#include "Ostap/Positive.h"
#include "Ostap/Lomont.h"
// ============================================================================
// lcoal
// ============================================================================
#include "Integrator1D.h"
#include "Exception.h"
#include "local_math.h"
// ============================================================================
/** @file
 *  Implementation file the for classes defiend in the file Osta/Positive.h 
 *  @date 2023-01-14 
 *  @author Vanya Belyaev Ivan.Belyave@itep.ru
 */
// ============================================================================/
// Karlin-Shapley 
// ===========================================================================
/*  constructor from the order and interval 
 *  @param   N polynomial degree
 *  @param   xmin miniaml x 
 *  @param   xmax maximal x 
 */  
// ============================================================================
Ostap::Math::KarlinShapley::KarlinShapley
( const unsigned short N    , 
  const double         xmin , 
  const double         xmax ) 
  : m_xmin      ( std::min ( xmin , xmax ) )
  , m_xmax      ( std::max ( xmax , xmin ) )
  , m_A         ( 1.0 ) 
  , m_sphere1   ( 0 == N ? 0 : 1 )  
  , m_sphere2   ( 2 <= N ? N - 1 : 0       ) 
  , m_troots    ( 2 <= N ? N + 1 : 2 , 0.0 ) 
  , m_workspace () 
{
  //
  static const std::string s_M1 { "xmin must be smaller than xmax!"                 } ;
  static const std::string s_T  { "Ostap::Math::KarlinShapley "                     } ;
  //
  Ostap::Assert ( m_xmin < m_xmax  , s_M1 , s_T , 280 ) ;
  //
  updateRoots () ;
}
// ============================================================================
/*  constructor from all parameters 
 *  @param   pars parameters
 *  @param   xmin miniaml x 
 *  @param   xmax maximal x 
 */  
// ============================================================================
Ostap::Math::KarlinShapley::KarlinShapley
( const std::vector<double>& pars , 
  const double               xmin , 
  const double               xmax )
  : m_xmin      ( std::min ( xmin , xmax ) )
  , m_xmax      ( std::max ( xmax , xmin ) )
  , m_A         ( 1 <= pars.size() ? std::abs ( pars [0] ) : 1.0 )
  , m_sphere1   ( 1 <  pars.size() ? 1 : 0 )
  , m_sphere2   ( pars.begin() + ( pars.size() < 2 ? pars.size() : 2 ) , 
                  pars.end  () )
  , m_troots    ( pars.size() < 2 ? 2 : pars.size() , 0.0 )
  , m_workspace () 
{
  //
  static const std::string s_M1 { "xmin must be smaller than xmax!"                 } ;
  static const std::string s_T  { "Ostap::Math::KarlinShapley "                     } ;
  //
  Ostap::Assert ( m_xmin < m_xmax  , s_M1 , s_T , 282 ) ;
  //
  updateRoots () ;
}
// ============================================================================
/*  constructor from scale and phases 
 *  @param   A      overall scale 
 *  @param   phases phases of the 1st sphere 
 *  @param   phases phases of the 2nd sphere 
 *  @param   xmin minimal x 
 *  @param   xmax maximal x 
 */  
// ============================================================================
Ostap::Math::KarlinShapley::KarlinShapley
( const double               A       , 
  const std::vector<double>& phases1 , 
  const std::vector<double>& phases2 , 
  const double               xmin    , 
  const double               xmax    ) 
  : m_xmin      ( std::min ( xmin , xmax ) )
  , m_xmax      ( std::max ( xmax , xmin ) )
  , m_A         ( std::abs ( A ) )
  , m_sphere1   ( phases1 )
  , m_sphere2   ( phases2 )
  , m_troots    ( 2 + phases2.size() , 0.0 )
  , m_workspace () 
{
  //
  static const std::string s_M1 { "xmin must be smaller than xmax!"      } ;
  static const std::string s_M2 { "1st sphere shodul be atmost 1 phase!"} ;
  static const std::string s_M3 { "1st sphere cannot be empty for non-empty 2nd!"} ;
  static const std::string s_T  { "Ostap::Math::KarlinShapley "          } ;
  //
  Ostap::Assert ( m_xmin < m_xmax    , s_M1 , s_T , 282 ) ;
  Ostap::Assert ( phases1.size() <=1 , s_M2 , s_T , 283 ) ;
  Ostap::Assert ( ( !phases1.empty() ) || 
                  ( phases1.empty() && phases2.empty() ) , s_M3 , s_T , 284 );
  //
  updateRoots () ;
}
// ============================================================================
/*  constructor from the scale, phase and phases 
 *  @param   A        global  scale 
 *  @param   phi      phase of the 1st sphere 
 *  @param   phases2  phases  of the 2nd sphere 
 *  @param   xmin minimal x 
 *  @param   xmax maximal x 
 */  
// ============================================================================
Ostap::Math::KarlinShapley::KarlinShapley
( const double               A       , 
  const double               phi     , 
  const std::vector<double>& phases2 , 
  const double               xmin    , 
  const double               xmax    ) 
  : KarlinShapley ( A , std::vector<double>(1,phi) , phases2 , xmin , xmax )
{}
// ============================================================================
// set parameter alpha 
// ============================================================================
bool Ostap::Math::KarlinShapley::setA ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_A , avalue ) ) { return false ; }
  m_A = avalue ;
  return true ; 
}
// ============================================================================
// evaluate Karlin-Shapley polynomial 
// ============================================================================
double Ostap::Math::KarlinShapley::evaluate ( const double x ) const 
{
  // zero-degree polynomial 
  if ( 0 == m_sphere1.npars() ) { return m_A ; }           // RETURN  
  //
  if ( s_zero ( m_A )         ) { return 0.0 ; }           // RETURN
  //
  //  get internal variable 
  const long double tt = t ( x ) ;
  //
  long double result = 0 ;
  const long double alpha = m_sphere1 . x2 ( 0 ) ;
  const long double beta  = m_sphere1 . x2 ( 1 ) ;
  //
  // 1st degree polynomial 
  if ( 0 == m_sphere2.npars() ) 
  { return alpha * tt + beta * ( 1.0L - tt ) ; }          // RETURN 
  //
  // number of t-roots 
  const std::size_t NT = m_troots.size() ;
  //
  // even degree ? (odd #troots) 
  const bool        even = 1 == ( NT % 2 ) ;
  const bool        odd  = !even ;
  //
  // ==========================================================================
  // calculate alpha-polynomial 
  // ==========================================================================
  if ( 3 <= NT && !s_zero ( alpha ) ) 
  {
    long double ap = 1.0L ;
    ap = 1.0L;
    for ( unsigned short k = even ? 1 : 2 ; k < NT ; k += 2 ) 
    { ap *= ( tt - m_troots [ k ] ) ; }
    ap = ap * ap ;
    if ( odd ) { ap *=  ( tt - m_troots.front () ) ; }
    result += alpha * ap ;
  }
  // ==========================================================================
  // calculate beta-polynomials 
  // ===========================================================================
  if ( 3 <= NT && !s_zero ( beta ) ) 
  {
    long double bp = 1.0L ;
    for ( unsigned short k = even ? 2 : 1 ; k + 1 < NT ; k +=2 ) 
    { bp *= ( tt - m_troots [ k ] ) ; }
    bp = bp * bp ;
    if ( even ) { bp *= ( tt - m_troots.front () ) * ( m_troots.back() - tt ) ; }
    else        { bp *=                              ( m_troots.back() - tt ) ; }
    result += beta * bp ;  
  }
  //
  return m_A * result ;
}
// ============================================================================
// get the unique tag 
// ============================================================================
std::size_t Ostap::Math::KarlinShapley::tag () const 
{
  static const std::string s_name {"KarlinShapley"} ;
  return Ostap::Utils::hash_combiner
    ( s_name ,
      m_A    , m_sphere1.tag () , m_sphere2.tag () ,
      m_xmin , m_xmax ) ;
} ;
// ============================================================================
// swap two polynomials
// ============================================================================
void Ostap::Math::KarlinShapley::swap 
( Ostap::Math::KarlinShapley& right ) 
{
  std::swap         ( m_xmin      , right.m_xmin      ) ;
  std::swap         ( m_xmax      , right.m_xmax      ) ;
  std::swap         ( m_A         , right.m_A         ) ;
  Ostap::Math::swap ( m_sphere1   , right.m_sphere1   ) ;
  Ostap::Math::swap ( m_sphere2   , right.m_sphere2   ) ;
  std::swap         ( m_troots    , right.m_troots    ) ;
  Ostap::Math::swap ( m_workspace , right.m_workspace ) ;
}
// ============================================================================
// update internal Karlin-Shapley t-roots 
// ============================================================================
void Ostap::Math::KarlinShapley::updateRoots() 
{
  //
  if ( m_troots.empty() ) { return ; }
  //
  m_troots.front () = 0 ;
  m_troots.back  () = 0 ;
  //
  const unsigned short np = m_sphere2.npars () ;
  for ( unsigned short k  = 0 ; k < np ; ++k ) 
  { m_troots [ k + 1 ] = m_sphere2.x2 ( k ) ; }
  //
  std::partial_sum ( m_troots.begin ()     , m_troots.end () , m_troots.begin() ) ;
  //
  m_troots.back  () = 1 ;
  //
}
// ============================================================================
/*  get (numerical) integral 
 * \f$ f = \int^{x_{max}}_{x_{min}} P(x) dx \f$ 
 */
// ============================================================================
double Ostap::Math::KarlinShapley::integral 
( const double xmin , 
  const double xmax ) const 
{
  //
  if      ( s_equal ( xmin , xmax ) ) {return 0 ; }  
  else if (           xmax < xmin   ) { return -integral ( xmax , xmin ) ; }
  //
  static const Ostap::Math::GSL::Integrator1D<KarlinShapley> s_integrator {} ;
  static char s_message[] = "Integral(KarlinShapley)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
    ( tag () ,                        // unique tag 
      &F     ,                        // the function 
      xmin   , xmax       ,           // low & high edges
      workspace ( m_workspace ) ,     // workspace
      s_APRECISION        ,           // absolute precision
      s_RPRECISION        ,           // relative precision
      m_workspace.size()  ,           // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================



// ============================================================================
// Karlin-Studden 
// ============================================================================
/*  constructor from the order and interval 
 *  @param   N polynomial degree
 *  @param   xmin miniaml x 
 *  @param   xmax maximal x 
 */  
// ============================================================================
Ostap::Math::KarlinStudden::KarlinStudden
( const unsigned short N     , 
  const double         xmin  , 
  const double         scale ) 
  : m_xmin      ( xmin               ) 
  , m_scale     ( std::abs ( scale ) ) 
  , m_A         ( 1.0 ) 
  , m_sphere1   ( 0 == N ? 0 : 1     )  
  , m_sphere2   ( 2 <= N ? N - 1 : 0 ) 
  , m_troots    ( N , 0.0 ) 
  , m_zroots    ( N , 0.0 ) 
  , m_workspace () 
{
  //
  updateRoots () ;
}
// ============================================================================
/*  constructor from all parameters 
 *  @param   pars parameters 
 *  @param   xmin minimal x 
 *  @param   scale global scale factor (for better numerical perfoamnce)
 */  
// ============================================================================
Ostap::Math::KarlinStudden::KarlinStudden
( const std::vector<double>& pars  , 
  const double               xmin  , 
  const double               scale ) 
  : m_xmin      ( xmin               ) 
  , m_scale     ( std::abs ( scale ) ) 
  , m_A         ( 1 <= pars.size() ? std::abs ( pars [0] ) : 1.0 )
  , m_sphere1   ( 1 <  pars.size() ? 1 : 0 )
  , m_sphere2   ( pars.begin() + ( pars.size() < 2 ? pars.size() : 2 ) , 
                  pars.end  () )
  , m_troots    ( pars.size() < 2 ? 1 : pars.size() , 0.0 )
  , m_zroots    ( pars.size() < 2 ? 1 : pars.size() , 0.0 )
  , m_workspace () 
{
  updateRoots () ;
}
// ============================================================================
/** constructor from the scale and phases 
 *  @param   A        global  scale 
 *  @param   phases1  phases of the 1st sphere 
 *  @param   phases2 phases  of the 2nd sphere 
 *  @param   xmin minimal x 
 *  @param   scale global scale factor (for better numerical perfoamnce)
 */  
// ============================================================================
Ostap::Math::KarlinStudden::KarlinStudden
( const double               A         , 
  const std::vector<double>& phases1   , 
  const std::vector<double>& phases2   , 
  const double               xmin      , 
  const double               scale     ) 
  
  : m_xmin      ( xmin               ) 
  , m_scale     ( std::abs ( scale ) ) 
  , m_A         ( std::abs ( A )     )
  , m_sphere1   ( phases1 )   
  , m_sphere2   ( phases2 ) 
  , m_troots    ( phases2.size() + 1 , 0.0 ) 
  , m_zroots    ( phases2.size() + 1 , 0.0 ) 
  , m_workspace () 
{
  static const std::string s_M3 { "1st sphere cannot be empty for non-empty 2nd!"} ;
  static const std::string s_T  { "Ostap::Math::KarlinStudden"          } ;
  //
  Ostap::Assert ( ( !phases1.empty() ) || 
                  ( phases1.empty() && phases2.empty() ) , s_M3 , s_T , 285 );
  //
  updateRoots () ;
}
// ============================================================================
/** constructor from the scale, phase and phases 
 *  @param   A        global  scale 
 *  @param   phi      phase of the 1st sphere 
 *  @param   phases2  phases  of the 2nd sphere 
 *  @param   xmin minimal x 
 *  @param   scale global scale factor (for better numerical perfoamnce)
 */  
// ============================================================================
Ostap::Math::KarlinStudden::KarlinStudden
( const double               A       , 
  const double               phi     , 
  const std::vector<double>& phases2 , 
  const double               xmin    , 
  const double               scale   ) 
  : KarlinStudden ( A , {{ phi }} , phases2 , xmin , scale )
{}
// ============================================================================
// set parameter alpha 
// ============================================================================
bool Ostap::Math::KarlinStudden::setA ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_A , avalue ) ) { return false ; }
  m_A = avalue ;
  return true ; 
}
// ============================================================================
// evaluate Karlin-Studdenpolynomial 
// ============================================================================
double Ostap::Math::KarlinStudden::evaluate ( const double x ) const 
{
  // zero-degree polynomial 
  if ( 0 == m_sphere1.npars() ) { return m_A ; }             // RETURN  
  //
  if ( s_zero ( m_A )         ) { return 0.0 ; }             // RETURN
  //
  //  get internal variable 
  const long double tt = t ( x ) ;
  //
  long double result = 0 ;
  const long double alpha = m_sphere1 . x2 ( 0 ) ;
  const long double beta  = m_sphere1 . x2 ( 1 ) ;
  //
  // 1st degree polynomial 
  if ( 0 == m_sphere2.npars() ) { return alpha * tt + beta ; } // RETURN 
  //
  // number of z-roots 
  const std::size_t NT = m_zroots.size() ;
  //
  // even degree ? (odd #troots) 
  const bool        even = ( NT % 2 ) ;
  const bool        odd  = !even ;
  //
  // ==========================================================================
  // calculate alpha-polynomial 
  // ==========================================================================
  if ( !s_zero ( alpha ) ) 
  {
    long double ap = 1.0L ;
    for ( unsigned short k = even ? 1 : 2 ; k < NT ; k += 2 ) 
    { ap *= ( tt - m_zroots [ k ] ) ; }
    ap = ap * ap ;
    if ( odd ) { ap *= tt ; }
    result += alpha * ap ;
  }
  // ==========================================================================
  // calculate beta-polynomials 
  // ===========================================================================
  if ( !s_zero ( beta ) ) 
  {
    long double bp = 1.0L ;
    for ( unsigned short k = even ? 2 : 1 ; k < NT ; k +=2 ) 
    { bp *= ( tt - m_zroots [ k ] ) ; }
    bp = bp * bp ;
    if ( even ) { bp *= tt ; }
    result += beta * bp ;  
  }
  //
  return m_A * result ;
}
// ============================================================================
// get the unique tag 
// ============================================================================
std::size_t Ostap::Math::KarlinStudden::tag () const 
{
  static const std::string s_name {"KarlinStudden"} ;
  return Ostap::Utils::hash_combiner
    ( s_name ,
      m_A    , m_sphere1.tag () , m_sphere2.tag () ,
      m_xmin , m_scale ) ;
} ;
// ============================================================================
// swap two polynomials
// ============================================================================
void Ostap::Math::KarlinStudden::swap 
( Ostap::Math::KarlinStudden& right ) 
{
  std::swap         ( m_xmin      , right.m_xmin      ) ;
  std::swap         ( m_scale     , right.m_scale     ) ;
  std::swap         ( m_A         , right.m_A         ) ;
  Ostap::Math::swap ( m_sphere1   , right.m_sphere1   ) ;
  Ostap::Math::swap ( m_sphere2   , right.m_sphere2   ) ;
  std::swap         ( m_troots    , right.m_troots    ) ;
  std::swap         ( m_zroots    , right.m_zroots    ) ;
  Ostap::Math::swap ( m_workspace , right.m_workspace ) ;
}
// ============================================================================
// update internal Karlin-Studden t-roots 
// ============================================================================
void Ostap::Math::KarlinStudden::updateRoots() 
{
  if ( m_troots.empty() ) { return ; }
  //
  m_troots.front () = 0 ;
  m_zroots.front () = 0 ;
  //
  const unsigned short np = m_sphere2.npars () ;
  for ( unsigned short k  = 0 ; k < np ; ++k ) 
  { m_troots [ k + 1 ] = m_sphere2.x2 ( k ) ; }
  //
  std::reverse     ( m_troots.begin () + 1 , m_troots.end () ) ;
  std::partial_sum ( m_troots.begin ()     , m_troots.end () , m_troots.begin() ) ;
  //
  const long double r_small { Ostap::Math::next_float ( 1.0f , -2 )  } ;
  std::transform   ( m_troots.begin () + 1 , 
                     m_troots.end   () , 
                     m_zroots.begin () + 1 , 
                     [r_small]( const long double r ) -> double 
                     {
                       const long double rr = std::min ( r , r_small ) ;
                       return rr / ( 1.0L - rr ) ; 
                     } ) ;
}
// ============================================================================
/*  get (numerical) integral 
 * \f$ f = \int^{x_{max}}_{x_{min}} P(x) dx \f$ 
 */
// ============================================================================
#include <iostream> 
#include "Ostap/ToStream.h" 
double Ostap::Math::KarlinStudden::integral 
( const double xmin , 
  const double xmax ) const 
{
  //
  if      ( s_equal ( xmin , xmax ) ) {return 0 ; }  
  else if (           xmax < xmin   ) { return -integral ( xmax , xmin ) ; }
  //
  static const Ostap::Math::GSL::Integrator1D<KarlinStudden> s_integrator {} ;
  static char s_message[] = "Integral(KarlinStudden)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
    ( tag () ,                        // unique tag 
      &F     ,                        // the function 
      xmin   , xmax       ,           // low & high edges
      workspace ( m_workspace ) ,     // workspace
      s_APRECISION        ,           // absolute precision
      s_RPRECISION        ,           // relative precision
      m_workspace.size()  ,           // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  if ( ierror ) 
  {
    std::cerr << " error:" << ierror << " troots: " ;
    Ostap::Utils::toStream ( m_troots , std::cerr ) ;
    std::cerr << std::endl ;
  }
  //
  return result ;
}






// ============================================================================
//                                                                      The END 
// ============================================================================
