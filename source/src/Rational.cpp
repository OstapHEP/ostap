// Include files 
// ============================================================================
// Incldue files 
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Interpolants.h"
#include "Ostap/Rational.h"
// ============================================================================
// Local 
// ===========================================================================
#include "Integrator1D.h"
#include "local_math.h"
#include "local_gsl.h"
// ===========================================================================
/** @file 
 *  Rational pole-less f
 *  A pole-free rational function at interval \f$ x_{min} \le x \le x_{max}\f$
 *  \f[ F(x) = \frac{p(x)}{q(x)} \f]
 *  Actually internally it uses 
 *  the Floater-Hormann rational barycentric interpolant 
 *  and parameters are the function valeus at Chebyshev's nodes  
 *  @
 #  @see Ostap::Math::FloaterHormann
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date   2023-09-14
 */
// ============================================================================

// ============================================================================
/*  constructor  
 *  @param n degree of numerator 
 *  @param d degree of denumerator is \f$ \max ( n - d , 0 ) \f$
 *  @param xmin  low-edge  of interval
 *  @param xmax  high-edge of interval
 */
// ============================================================================
Ostap::Math::Rational::Rational
( const unsigned short n    , 
  const unsigned short d    , 
  const double         xmin , 
  const double         xmax ) 
  : Parameters ( n ) 
  , m_abscissas ( n , xmin , xmax ) 
  , m_weights   ( n , d           )
  , m_workspace ()
{}
// ============================================================================
/*  constructor  
 *  @param pars  vector of parameters 
 *  @param d     degree defect,  
 *  @param xmin  low-edge  of interval
 *  @param xmax  high-edge of interval
 */
// ============================================================================
Ostap::Math::Rational::Rational
( const std::vector<double>& pars , 
  const unsigned short       d    , 
  const double               xmin , 
  const double               xmax ) 
  : Parameters  ( pars ) 
  , m_abscissas ( pars.size() , xmin , xmax ) 
  , m_weights   ( pars.size() , d )
  , m_workspace ()
{}
// ============================================================================
// evaluate the rational function 
// ============================================================================
double Ostap::Math::Rational::evaluate ( const double x ) const 
{
  if ( x < xmin() || xmax () < x ) { return 0 ; }
  //
  long double s1 = 0 ;
  long double s2 = 0 ;
  long double z  = 1.0L * x ;
  //
  unsigned int N = n () ;
  for ( unsigned int i = 0  ; i < N  ;  ++i )
  {
    const long double xi = m_abscissas.x ( i ) ;
    const long double yi = par ( i ) ;
    //
    if ( s_equal ( x , xi ) ) { return yi ; }  // RETURN
    //
    const long double wi = m_weights.weight ( i ) / ( z - xi ) ;
    //
    s1 += wi * yi ;
    s2 += wi ;
  }
  return s1 / s2 ;
}
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::Rational::integral () const 
{ return integral ( xmin() , xmax() ) ; }
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::Rational::integral 
(  const double xlow  , 
   const double xhigh ) const 
{
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; }
  else if ( xhigh <  xlow            ) { return - integral ( xhigh , xlow ) ; }
  else if ( xhigh <= xmin ()         ) { return 0 ; }
  else if ( xlow  >= xmax ()         ) { return 0 ; }
  //
  const double xmn = std::max ( xlow  , xmin () ) ;
  const double xmx = std::min ( xhigh , xmax () ) ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Rational> s_integrator ;
  static const char s_message[] = "Integral(Rational)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag  () , 
      &F      , 
      xmn     , xmx ,             // low & high edges
      workspace ( m_workspace ) , // workspace
      s_APRECISION              , // absolute precision
      s_RPRECISION              , // relative precision
      m_workspace.size()        , // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ; 
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Rational::tag () const 
{ 
  static const std::string s_name = "Rational" ;
  return Ostap::Utils::hash_combiner 
    ( s_name  , 
      d    () , 
      xmin () ,
      xmax () , 
      Ostap::Utils::hash_range ( pars () ) ) ;
}
// ============================================================================

  
  

// ============================================================================
Ostap::Math::RationalBernstein::RationalBernstein
( const unsigned short       p     ,
  const unsigned short       q     , 
  const double               xmin  , 
  const double               xmax  ) 
  : m_p ( p , xmin , xmax ) 
  , m_q ( q , xmin , xmax ) 
  , m_workspace () 
{}
// ============================================================================
Ostap::Math::RationalBernstein::RationalBernstein
( const std::vector<double>& p    ,
  const std::vector<double>& q    ,     
  const double               xmin , 
  const double               xmax ) 
  : m_p ( p , xmin , xmax ) 
  , m_q ( q , xmin , xmax ) 
  , m_workspace () 
{}
// ============================================================================    
Ostap::Math::RationalBernstein::RationalBernstein
( const std::vector<double>& a    ,
  const unsigned short       p    ,     
  const double               xmin , 
  const double               xmax ) 
  : m_p ( a.begin() , a.begin() + ( p < a.size() ? p : a.size() ), xmin , xmax ) 
  , m_q ( a.begin() + ( p < a.size() ? p : a.size() ) , a.end()  , xmin , xmax ) 
  , m_workspace () 
{}
// ============================================================================    
// evaluate the rational function
// ============================================================================    
double Ostap::Math::RationalBernstein::evaluate ( const double x ) const 
{ return x < xmin () ? 0.0 : x > xmax () ? 0.0 : m_p ( x ) / m_q ( x ) ; }
// ============================================================================
/// all parameters (by value!!!)
// ============================================================================    
std::vector<double> 
Ostap::Math::RationalBernstein::pars () const 
{
  const unsigned short np = npars() ;
  std::vector<double> result ( np , 0.0 ) ;
  for ( unsigned short i = 0 ; i < np ; ++i ) { result [ i ] = par ( i ) ; }
  return result ;
}
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::RationalBernstein::integral () const 
{ return integral ( xmin() , xmax() ) ; }
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::RationalBernstein::integral 
(  const double xlow  , 
   const double xhigh ) const 
{
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; }
  else if ( xhigh <  xlow            ) { return - integral ( xhigh , xlow ) ; }
  else if ( xhigh <= xmin ()         ) { return 0 ; }
  else if ( xlow  >= xmax ()         ) { return 0 ; }
  //
  const double xmn = std::max ( xlow  , xmin () ) ;
  const double xmx = std::min ( xhigh , xmax () ) ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<RationalBernstein> s_integrator ;
  static const char s_message[] = "Integral(RationalBernstein)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag  () , 
      &F      , 
      xmn     , xmx ,             // low & high edges
      workspace ( m_workspace ) , // workspace
      s_APRECISION              , // absolute precision
      s_RPRECISION              , // relative precision
      m_workspace.size()        , // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ; 
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::RationalBernstein::tag () const 
{ 
  static const std::string s_name = "RationalBernstein" ;
  return Ostap::Utils::hash_combiner 
    ( s_name     , 
      xmin    () ,
      xmax    () , 
      m_p.tag () , 
      m_q.tag () ) ;
}
// ============================================================================

// ============================================================================
Ostap::Math::RationalPositive::RationalPositive
( const unsigned short       p     ,
  const unsigned short       q     , 
  const double               xmin  , 
  const double               xmax  ) 
  : m_p ( p , xmin , xmax ) 
  , m_q ( q , xmin , xmax ) 
  , m_workspace () 
{}
// ============================================================================
Ostap::Math::RationalPositive::RationalPositive
( const std::vector<double>& p    ,
  const std::vector<double>& q    ,     
  const double               xmin , 
  const double               xmax ) 
  : m_p ( p , xmin , xmax ) 
  , m_q ( q , xmin , xmax ) 
  , m_workspace () 
{}
// ============================================================================    
Ostap::Math::RationalPositive::RationalPositive
( const std::vector<double>& a    ,
  const unsigned short       p    ,     
  const double               xmin , 
  const double               xmax ) 
  : m_p ( a.begin() , a.begin() + ( p < a.size() ? p : a.size() ), xmin , xmax ) 
  , m_q ( a.begin() + ( p < a.size() ? p : a.size() ) , a.end()  , xmin , xmax ) 
  , m_workspace () 
{}
// ============================================================================    
// evaluate the rational function
// ============================================================================    
double Ostap::Math::RationalPositive::evaluate ( const double x ) const 
{ return x < xmin () ? 0.0 : x > xmax () ? 0.0 : m_p ( x ) / m_q ( x ) ; }
// ============================================================================
/// all parameters (by value!!!)
// ============================================================================    
std::vector<double> 
Ostap::Math::RationalPositive::pars () const 
{
  const unsigned short nn = npars() ;
  std::vector<double> result ( nn , 0 ) ;
  for ( unsigned short i = 0 ; i < nn ; ++i ) { result [ i ] = par ( i ) ; }
  return result ;
}
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::RationalPositive::integral () const 
{ return integral ( xmin() , xmax() ) ; }
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::RationalPositive::integral 
(  const double xlow  , 
   const double xhigh ) const 
{
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; }
  else if ( xhigh <  xlow            ) { return - integral ( xhigh , xlow ) ; }
  else if ( xhigh <= xmin ()         ) { return 0 ; }
  else if ( xlow  >= xmax ()         ) { return 0 ; }
  //
  const double xmn = std::max ( xlow  , xmin () ) ;
  const double xmx = std::min ( xhigh , xmax () ) ;
  if ( s_equal ( xmn , xmx )         ) { return 0 ; }
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<RationalPositive> s_integrator ;
  static const char s_message[] = "Integral(RationalPositive)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag  () , 
      &F      , 
      xmn     , xmx ,             // low & high edges
      workspace ( m_workspace ) , // workspace
      s_APRECISION              , // absolute precision
      s_RPRECISION              , // relative precision
      m_workspace.size()        , // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ; 
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::RationalPositive::tag () const 
{ 
  static const std::string s_name = "RationalPositive" ;
  return Ostap::Utils::hash_combiner 
    ( s_name     , 
      xmin    () ,
      xmax    () , 
      m_p.tag () , 
      m_q.tag () ) ;  
}
// ============================================================================

// ============================================================================
//                                                                      The END 
// ============================================================================
