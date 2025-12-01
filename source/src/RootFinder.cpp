//// =============================================================================
// Include files
// =============================================================================
// STD&STL
// =============================================================================
#include <functional>
#include <array>
#include <algorithm> 
// =============================================================================
// Ostap
// =============================================================================
#include "Ostap/RootFinder.h"
// =============================================================================
// Local
// =============================================================================
#include "local_math.h"
#include "status_codes.h"
// =============================================================================
/** @file
 *  Implementation file for class Ostap::Math::RootFinder
 *  @date 2018-10-11
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 */
// =============================================================================
namespace
{
  // ===========================================================================
  // Is this function value could be taken as root ?
  inline bool is_root 
  ( const double fvalue , 
    const double froot  ) 
  { return ( 0 < froot && std::abs ( fvalue ) <= froot ) || s_zero ( fvalue ) ;  }
  // ===========================================================================
  inline bool is_root
  ( const Ostap::Math::RootFinder::Point& r      , 
    const double                          froot )
  { return is_root ( r.fx () , froot ) ; }
  // ===========================================================================
  inline double secant 
  ( const double a   , 
    const double b   , 
    const double fa  , 
    const double fb  )
  { return ( a * fb - b * fa ) / ( fb - fa ) ; }
  // ===========================================================================
  inline double secant 
  ( const Ostap::Math::RootFinder::Point& a , 
    const Ostap::Math::RootFinder::Point& b )
  { return secant ( a.x () , b.x () , a.fx () , b.fx () ) ; }
  // ===========================================================================
  inline bool bracket  
  ( const double fa , 
    const double fb )
  { return ( 0 < fa && fb < 0 ) || ( 0 < fb && fa < 0 ) ;  }
  // ============================================================================
  inline bool bracket 
  ( const Ostap::Math::RootFinder::Point& a , 
    const Ostap::Math::RootFinder::Point& b )
  { return bracket ( a.fx () , b.fx () ) ; }
  // ===========================================================================
  template <class FUNCTION, 
            class DERIVATIVE1,
            class DERIVATIVE2>
  inline bool newton_halley 
  ( const FUNCTION&                 fun    , 
    const DERIVATIVE1&              deriv1 , 
    const DERIVATIVE2&              deriv2 , 
    Ostap::Math::RootFinder::Point& r  , 
    Ostap::Math::RootFinder::Point& a  , 
    Ostap::Math::RootFinder::Point& b  , 
    std::size_t&                    ncalls )
  {
    // newton step 
    const double d1 = deriv1  ( r.x() ) ;
    ++ncalls ;
    if ( s_zero ( d1 ) ) { return false ; }
    double rn = r.fx() / d1 ;
    // halley's correction
    if ( deriv2 ) 
    {
      const double d2 = deriv2 ( r.x () ) ; 
      ++ncalls ; 
      if ( !s_zero ( d2 ) ) { rn /= ( 1.0 - 0.5 * rn * d2 / d1 ) ; }
    }
    const double x = r.x () - rn ;
    if ( a.x() <= x && x <= b.x() )
    {
      r = Ostap::Math::RootFinder::Point ( x , fun ( x ) ) ;
      ++ncalls ;
      //
      if      ( bracket ( a.fx () , r.fx () ) ) { b = r ; }
      else if ( bracket ( b.fx () , r.fx () ) ) { a = r ; }
      //
      return true ; 
    } 
    return false ;
  }
// ==========================================================================
/** Single step of STFA method: Imporved regular falsi + steffenson-like
 *  @see Xinyuan Wu, Zuhe Shen, Jianlin Xia, "An improved regula falsi method
 *  with quadratic convergence of both diameter and point for enclosing
 *  simple zeros of nonlinear equations", Applied Mathematics and Computation,
 *  144, 2 2003
 *  @see https://doi.org/10.1016/S0096-3003(02)00414-9
 *  @see https://www.sciencedirect.com/science/article/pii/S0096300302004149
 */
 template <class FUNCTION>
 inline bool SFTA
 ( const FUNCTION&fun , 
   Ostap::Math::RootFinder::Point& r      , 
   Ostap::Math::RootFinder::Point& a      , 
   Ostap::Math::RootFinder::Point& b      , 
   std::size_t&                    ncalls , 
   const double                    froot  )
 {
  //
  if ( a.x() < r.x () && r.x () < b.x() )  {}
  else 
  {
    // a bit of black magic. Make  a secant-like step, but not a secant one  
    const double aa = std::abs ( a.fx () ) ;
    const double bb = std::abs ( b.fx () ) ;
    // 
    const double x = ( aa < bb ) ?
      secant ( a.x () , b.x () , a.fx () , b.fx () * ( 0.1 * aa + 0.9 * bb ) / bb ) : 
      secant ( b.x () , a.x () , b.fx () , a.fx () * ( 0.9 * aa + 0.1 * bb ) / aa ) ; 
    //
    r = Ostap::Math::RootFinder::Point ( x , fun ( x ) ) ; 
    ++ncalls ;
    if ( is_root ( r , froot ) ) { return true ; } 
  }

  /// regular falsi step
  const double xc = secant ( a , b ) ;
  Ostap::Math::RootFinder::Point c ( xc , fun ( xc ) ) ;
  if ( is_root ( c       , froot ) )    { r = c  ; return true ; } 

  if ( s_equal ( c.fx () , r.fx () ) )  { r = c  ; return true ; }

  Ostap::Math::RootFinder::Point  abar(0,0) , bbar(0,0) ;
  if      ( bracket ( a , c ) ) { abar = a ; bbar = c ; }
  else if ( bracket ( b , c ) ) { abar = c ; bbar = b ;  }   

  const double mu  = ( b.x () - a.x () ) / ( b.fx () - a.fx () ) ;
  const double fx  = r.fx() ; 
  const double cbx = r.x () - mu * fx * fx / ( fx - c.fx () ) ;

  if ( a.x() < cbx && cbx < b.x() ) {} 
  else                              
  { 
    r = c    ;
    a = abar ; 
    b = bbar ; 
    return true ; 
  }

  Ostap::Math::RootFinder::Point  cbar( cbx , fun ( cbx ) ) ;
  ++ncalls ;
  if ( is_root ( cbar , froot ) )    { r = cbar  ; return true ; } 

  if ( abar.x() <= cbar.x() && cbar.x() <= bbar.x() )
  {
    r = cbar ;
    if      ( bracket ( abar , cbar ) ) { a = abar ; b = cbar ; }
    else if ( bracket ( bbar , cbar ) ) { a = cbar ; b = bbar ; } 
    return true ;
  }
  //
  r = c ;
  //
  if      ( bracket ( a , r ) ) { b = r ; }
  else if ( bracket ( b , r ) ) { a = r ; }
  // 
  return true ; 
 }
 // =========================================================================
 /// inverse parabolic interpolation
 inline double inverse_parabolic 
 ( const Ostap::Math::RootFinder::Point& a , 
   const Ostap::Math::RootFinder::Point& b , 
   const Ostap::Math::RootFinder::Point& c ) 
 {
   //
   if      ( s_equal ( a.x  () , b.x  () ) ) { return secant ( a , c ) ; }
   else if ( s_equal ( a.x  () , c.x  () ) ) { return secant ( a , b ) ; }
   else if ( s_equal ( b.x  () , c.x  () ) ) { return secant ( a , b ) ; }
   else if ( s_equal ( a.fx () , b.fx () ) ) { return secant ( a , c ) ; }
   else if ( s_equal ( a.fx () , c.fx () ) ) { return secant ( a , b ) ; }
   else if ( s_equal ( b.fx () , c.fx () ) ) { return secant ( a , b ) ; }
   // 
   const double x0 { a.x () } , f0 { a.fx () } ;
   const double x1 { b.x () } , f1 { b.fx () } ;
   const double x2 { c.x () } , f2 { c.fx () } ; 
   // 
   const double f01 = 1.0 / ( f0 - f1 ) ; const double f10 = - f01 ;
   const double f02 = 1.0 / ( f0 - f2 ) ; const double f20 = - f02 ;
   const double f12 = 1.0 / ( f1 - f2 ) ; const double f21 = - f12 ; 
   //
   return x0 * f1 * f2 * f01 * f02 + x1 * f0 * f2 * f10 * f12 + x2 * f0 * f1 * f20 * f21 ; 
 }
 // ==========================================================================
 /// inverse cubic  interpolation
 inline double inverse_cubic  
 ( const Ostap::Math::RootFinder::Point& a , 
   const Ostap::Math::RootFinder::Point& b , 
   const Ostap::Math::RootFinder::Point& c , 
   const Ostap::Math::RootFinder::Point& d ) 
 {
  //
  if      ( s_equal ( a.x   () , b.x   () ) ) { return inverse_parabolic ( a , c , d ) ; }
  else if ( s_equal ( a.x   () , c.x   () ) ) { return inverse_parabolic ( a , b , d ) ; }
  else if ( s_equal ( a.x   () , d.x   () ) ) { return inverse_parabolic ( a , b , c ) ; }
  else if ( s_equal ( b.x   () , c.x   () ) ) { return inverse_parabolic ( a , b , d ) ; }
  else if ( s_equal ( b.x   () , d.x   () ) ) { return inverse_parabolic ( a , b , c ) ; }
  else if ( s_equal ( c.x   () , d.x   () ) ) { return inverse_parabolic ( a , b , c ) ; }
  else if ( s_equal ( a.fx  () , b.fx  () ) ) { return inverse_parabolic ( a , c , d ) ; }
  else if ( s_equal ( a.fx  () , c.fx  () ) ) { return inverse_parabolic ( a , b , d ) ; }
  else if ( s_equal ( a.fx  () , d.fx  () ) ) { return inverse_parabolic ( a , b , c ) ; }
  else if ( s_equal ( b.fx  () , c.fx  () ) ) { return inverse_parabolic ( a , b , d ) ; }
  else if ( s_equal ( b.fx  () , d.fx  () ) ) { return inverse_parabolic ( a , b , c ) ; }
  else if ( s_equal ( c.fx  () , d.fx  () ) ) { return inverse_parabolic ( a , b , c ) ; }

  const double x0 { a.x () } , f0 { a.fx () } ;
  const double x1 { b.x () } , f1 { b.fx () } ;
  const double x2 { c.x () } , f2 { c.fx () } ; 
  const double x3 { d.x () } , f3 { d.fx () } ; 

  const double f01 = 1.0 / ( f0 - f1 ) ; const double f10 = - f01 ;
  const double f02 = 1.0 / ( f0 - f2 ) ; const double f20 = - f02 ;
  const double f03 = 1.0 / ( f0 - f3 ) ; const double f30 = - f03 ;
  const double f12 = 1.0 / ( f1 - f2 ) ; const double f21 = - f12 ;
  const double f13 = 1.0 / ( f1 - f3 ) ; const double f31 = - f13 ;
  const double f23 = 1.0 / ( f2 - f3 ) ; const double f32 = - f23 ;
  //
  return  -x0 * f1 * f2 * f3 * f01 * f02 * f03
          -x1 * f0 * f2 * f3 * f10 * f12 * f13
          -x2 * f0 * f2 * f3 * f20 * f21 * f23
          -x3 * f0 * f1 * f2 * f30 * f31 * f32  ;
 }
 // ==========================================================================
 template <class FUNCTION>
 inline bool toms748
 ( const FUNCTION&                 fun    , 
   Ostap::Math::RootFinder::Point& r      , 
   Ostap::Math::RootFinder::Point& a      , 
   Ostap::Math::RootFinder::Point& b      , 
   std::size_t&                    ncalls , 
   const double                    froot  ) 
  {
   // 
   if ( r.x() < a.x() || b.x () < r.x() )
   {
     const  double x = secant  ( a , b ) ;
     r = Ostap::Math::RootFinder::Point ( x  , fun ( x ) ) ;
     ++ncalls ; 
     if ( is_root ( r , froot ) ) { return true ; } 
   }
   // 
   const double d = inverse_parabolic ( a  , b , r ) ;
   if ( !std::isfinite ( d ) || d <= a.x() || b.x() <= d ) { return false ; }

   Ostap::Math::RootFinder::Point pd ( d , fun ( d ) ) ;
   ++ncalls ; 
   if ( is_root ( pd , froot ) ) { r = pd ; return true ; } 

   const double e = inverse_cubic ( a , b , r , pd ) ; 
   if ( !std::isfinite ( e ) || e <= a.x() || b.x() <= e ) { return false ; }

   Ostap::Math::RootFinder::Point pe ( e , fun ( e ) ) ;
   ++ncalls ; 
   if ( is_root ( pe , froot ) ) { r = pe ; return true ; } 

   std::array<Ostap::Math::RootFinder::Point,3> p { r , pe , pd } ;
   std::sort ( p.begin () , p.end () ) ;

   if      ( bracket ( a       , p [ 0 ] ) ) { b = p [ 0 ]               ; }
   else if ( bracket ( p [ 0 ] , p [ 1 ] ) ) { a = p [ 0 ] ; b = p [ 1 ] ; }
   else if ( bracket ( p [ 1 ] , p [ 2 ] ) ) { a = p [ 1 ] ; b = p [ 2 ] ; }
   else if ( bracket ( p [ 2 ] , b       ) ) { a = p [ 2 ]               ; }

   if      ( a.x () <= pe.x () && pe.x () <= b.x () ) { r = pe ; return true ; }
   else if ( a.x () <= pd.x () && pd.x () <= b.x () ) { r = pd ; return true ; }
   else if ( a.x () <= r .x () && r .x () <= b.x () ) {          return true ; }

   const double x = secant ( a , b ) ;
   r = Ostap::Math::RootFinder::Point  ( x , fun ( x ) ) ;
   ++ncalls    ;
   //
   if      ( bracket ( a , r ) ) { b = r ; }
   else if ( bracket ( b , r ) ) { a = r ; }
   // 
   return true ; 
 }
 // ===========================================================================
}
// ===========================================================================
// status code for max callas limit
// ===========================================================================
const Ostap::StatusCode Ostap::Math::RootFinder::NumCallsLimit { NUM_CALLS_LIMIT_REACHED } ; 
// ===========================================================================
/* constructor from full configuration 
 *  @param max_calls  maximal number of function calls
 *  @param froot      consider x as root if  \f$ \left| f(x) \right| < f_{root} \f$ if \f$ 0 < f_{root}\f$ 
 *  @param atolerance absolute tolerance 
 *  @param rtolerance relative tolerance  
 */
// ============================================================================
Ostap::Math::RootFinder::RootFinder
( const unsigned short max_calls  ,
  const double         froot      ,
  const double         atolerance , 
  const double         rtolerance ) 
: m_maxcalls   ( 20 < max_calls  ? max_calls  :  20    )
, m_froot      (  0 < froot      ?  froot     : -1.0   )
, m_atolerance (  0 < atolerance ? atolerance :  1.e-9 )
, m_rtolerance (  0 < rtolerance ? rtolerance :  1.e-9 )
, m_ncalls     (  0 ) 
{}
// ===========================================================================
/*  find a root in [a,b]
 *  @param fun the function
 *  @param r  (update) the root 
 *  @param a  (update) the bracket interval
 *  @param b  (update) the bracket inteval
 *  @return status code 
 */
 // ==========================================================================
Ostap::StatusCode  
Ostap::Math::RootFinder::root 
( Ostap::Math::RootFinder::function1 fun , 
  double& r   , 
  double& a   , 
  double& b   ,
  Ostap::Math::RootFinder::function1 derivative1 , 
  Ostap::Math::RootFinder::function1 derivative2 ) const 
{
  /// reset the counter  
  m_ncalls = 0 ;
  //
  Ostap::Assert ( !!fun                                  , 
                  "Invalid std::function"                ,
                  "Ostap::Math::RootFinder::root"        , 
                  INVALID_FUNCTION , __FILE__ , __LINE__ ) ; 
  //
  const double fa = fun ( a ) ;
  ++m_ncalls ; 
  if ( is_root ( fa , m_froot ) ) { r = a ; return Ostap::StatusCode::SUCCESS ; }
  // 
  const double fb = fun ( b ) ;
  ++m_ncalls ;
  if ( is_root ( fb , m_froot ) ) { r = b ; return Ostap::StatusCode::SUCCESS ; }
  //
  // valid bracket interval
  Ostap::Assert ( bracket ( fa , fb )                   ,  
                  "Invalid bracker interval"            , 
                  "Ostap::Math::RootFinder::root"       , 
                  INVALID_BRACKET , __FILE__ , __LINE__ ) ; 

  Point pa  ( a , fa ) ;
  Point pb  ( b , fb ) ; 
  if ( pb.x ()  <  pa.x() ) { Ostap::Math::swap ( pa , pb ) ; }
  //
  if ( pa.x () <= r && r <= pb.x() )  {}
  else  { r = secant ( pa , pb ) ; } 
  //
  Point pr ( r , fun ( r ) ) ;
  ++m_ncalls ;
  if ( is_root ( pr , m_froot ) ) { return Ostap::StatusCode::SUCCESS ; }
  //
  const auto cfun = std::cref ( fun ) ;
  const auto d1   = std::cref ( derivative1 ) ;
  const auto d2   = std::cref ( derivative2 ) ; 
  /// call the main method 
  Ostap::StatusCode sc = root ( cfun        , 
                                pr          , 
                                pa          , 
                                pb          , 
                                d1          , 
                                d2          ) ;
  /// copy results  
  r = pr.x () ; 
  a = pa.x () ;
  b = pb.x () ;
  // 
  return sc ;
}
// =============================================================================
/*  find a root in [a,b]
 *  @param fun the function
 *  @param r  (update) the root 
 *  @param a  (update) the bracket interval
 *  @param b  (update) the bracket innteval
 *  @return status code 
 */
 // =============================================================================
Ostap::StatusCode  
Ostap::Math::RootFinder::root
( Ostap::Math::RootFinder::function1 fun         , 
  Ostap::Math::RootFinder::Point&    r           , 
  Ostap::Math::RootFinder::Point&    a           ,  
  Ostap::Math::RootFinder::Point&    b           ,   
  Ostap::Math::RootFinder::function1 derivative1 ,
  Ostap::Math::RootFinder::function1 derivative2 ) const 
{
  // 
  if      ( is_root ( a , m_froot ) ) { r = a ; return Ostap::StatusCode::SUCCESS ; }
  else if ( is_root ( b , m_froot ) ) { r = b ; return Ostap::StatusCode::SUCCESS ; }
  // 
  if      ( b.x() < a.x () ) { Ostap::Math::swap ( a , b ) ; }
  //
  // valid bracketing interval
  Ostap::Assert ( bracket  ( a.fx () , b.fx () )       ,  
                  "Invalid bracket interval"            , 
                  "Ostap::Math::RootFinder::root_"      , 
                  INVALID_BRACKET , __FILE__ , __LINE__ ) ; 
  // ========================================================================
  // check r 
  if ( a.x () <= r.x () && r.x () <= b.x () && is_root ( r.fx () , m_froot ) ) { return Ostap::StatusCode::SUCCESS ; }
  //
  const auto cfun  = std::cref ( fun    ) ;
  const auto cder1 = std::cref ( derivative1 ) ;
  const auto cder2 = std::cref ( derivative2 ) ;

  // start iterations
  double x = ( a.x() < r.x () && r.x() < b.x() ) ? r.x() : secant ( a , b ) ; 
  //
  // iterate while the maximum number of function calls is reached 
  do  
  {
    /// (1) make a single (combined) step 
    Ostap::StatusCode sc = step ( cfun , r , a , b , cder1 , cder2 ) ;  
    if ( sc.isFailure   () ) { return sc ; }
    //
    /// (2) root is found? 
    if ( is_root ( r , m_froot ) )                                 { return Ostap::StatusCode::SUCCESS ; }

    /// (3) delta  is very small 
    const double dx = std::abs ( r.x () - x ) ;
    if ( 2 * dx <= std::max ( m_atolerance , m_rtolerance * dx ) ) { return Ostap::StatusCode::SUCCESS ; }

    /// (4) bracket is small enough  
    const double ab = std::abs ( b.x () - a.x () ) ;
    if ( 2 * ab <= std::max ( m_atolerance , m_rtolerance * ab ) ) { return Ostap::StatusCode::SUCCESS ; }

    x = r.x () ;
  }
  while ( m_ncalls < m_maxcalls ) ; 
  // 
  return NumCallsLimit ;
}
// =============================================================================
/*  find a root in [a,b]
 *  @param fun the function
 *  @param r  (update) the root 
 *  @param a  (update) the bracket interval
 *  @param b  (update) the bracket innteval
 *  @return status code 
 */
 // =============================================================================
Ostap::StatusCode  
Ostap::Math::RootFinder::step
( Ostap::Math::RootFinder::function1 fun         , 
  Ostap::Math::RootFinder::Point&    r           , 
  Ostap::Math::RootFinder::Point&    a           ,
  Ostap::Math::RootFinder::Point&    b           ,
  Ostap::Math::RootFinder::function1 derivative1 ,
  Ostap::Math::RootFinder::function1 derivative2 ) const 
{
  // 
  if      ( is_root ( a , m_froot ) ) { r = a ; return Ostap::StatusCode::SUCCESS ; }
  else if ( is_root ( b , m_froot ) ) { r = b ; return Ostap::StatusCode::SUCCESS ; }
  //
  if ( b.x () < a.x () ) { Ostap::Math::swap ( a , b ) ; }
  //
  //  root is already found? 
  if ( a.x () <= r.x () && r.x () <= b.x () && is_root ( r.fx () , m_froot ) ) 
  { return Ostap::StatusCode::SUCCESS ; }
  //
  // valid bracketing interval ? 
  Ostap::Assert ( bracket  ( a , b )                    ,  
                  "Invalid bracket interval"            , 
                  "Ostap::Math::RootFinder::step"       , 
                  INVALID_BRACKET , __FILE__ , __LINE__ ) ; 

  // initial approximation is valid, otherwise make a secant step and (re)bracket
 // if ( a.x() <= r.x() && r.x() <= b.x() ) {} 
  // else 
  {
    const double rx = secant ( a, b ) ;
    r = Point  ( rx , fun ( rx ) ) ;
    ++m_ncalls ;
    //
    if      ( bracket ( a , r ) ) { b = r ; }
    else if ( bracket ( b , r ) ) { a = r ; } 
    //   
    if ( is_root ( r , m_froot ) ) { return Ostap::StatusCode::SUCCESS ; }
  }
  //
  const double r0      = r.x() ;
  bool         updated = false ;
  const double lenght  = std::abs ( b.x () - a.x() ) ;

  // ======================================================================== 
  /// (1) make a try with Newton-Halley method
  if ( derivative1 )
  {
    updated = newton_halley ( fun , derivative1 , derivative2 , r , a , b , m_ncalls ) ;
    if ( updated && a.x () <= r.x () && r.x() <= b.x() && 
         is_root ( r , m_froot ) ) {  return Ostap::StatusCode::SUCCESS ; }
  }
  
  /// single step of SFTA method  
  if ( !updated )  
  {
    updated = SFTA ( fun , r , a , b , m_ncalls , m_froot ) ;  
    if ( updated && a.x () <= r.x () && r.x() <= b.x() && 
         is_root ( r , m_froot ) ) {  return Ostap::StatusCode::SUCCESS ; }
  }

  /// use TOM748 
  if ( !updated )
  {
    updated = toms748 ( fun , r , a , b , m_ncalls , m_froot ) ;  
    if ( updated && a.x () <= r.x () && r.x() <= b.x() && 
         is_root ( r , m_froot ) ) {  return Ostap::StatusCode::SUCCESS ; }
  }

  /// use a bullet-proof secant method 
  if ( !updated )
  {
    const double x = secant ( a , b ) ;
    r = Point ( x , fun ( x ) ) ;
    ++m_ncalls ;
    if      ( bracket ( a , r ) )  { b = r ; }
    else if ( bracket ( b , r ) )  { a = r ; }
    if ( is_root ( r , m_froot ) ) {  return Ostap::StatusCode::SUCCESS ; }
  } 

  /// use bisection as the "ultima ratio regum"
  if ( !updated 
    || r.x () < a.x () 
    || r.x () > b.x ()
    || ( lenght <= 3 * std::abs ( b.x () - a.x () ) ) )
    {
      const double x = 0.5 * ( a.x () + b.x() ) ;
      Point c ( x , fun ( x ) ) ;
      ++m_ncalls ;
      //
      if ( a.x() <= r.x () && r.x() <= b.x () )
      {
        if      ( bracket ( a , c ) && r.x () <= c.x () )  { b = c ; return Ostap::StatusCode::SUCCESS ; }
        else if ( bracket ( b , c ) && r.x () >= c.x () )  { a = c ; return Ostap::StatusCode::SUCCESS ; } 
      }
      // 
      r = c ; 
      if      ( bracket ( a , r ) )  { b = c ; }
      else if ( bracket ( b , r ) )  { a = c ; }
      //
      if ( is_root ( c , m_froot ) ) {  return Ostap::StatusCode::SUCCESS ; }
    }
  //
  return Ostap::StatusCode::SUCCESS ;
}
// =============================================================================
//                                                                      The END 
// =============================================================================
