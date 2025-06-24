// Include files 
// ============================================================================
// GSL
// ===========================================================================
#include "gsl/gsl_linalg.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Interpolation.h"
#include "Ostap/Polynomials.h"
#include "Ostap/Rational.h"
#include "Ostap/Clenshaw.h"
#include "Ostap/LinAlg.h"
#include "Ostap/StatusCode.h"
// ============================================================================
// Local 
// ===========================================================================
#include "Integrator1D.h"
#include "local_hash.h"
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
  : Parameters  ( n ) 
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
  return ( s1 / s2 ) ;
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
      s_message                 , 
      __FILE__ , __LINE__       ,
      GSL_INTEG_GAUSS21         ) ;
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
// scale it
// ============================================================================
Ostap::Math::Rational&
Ostap::Math::Rational::scale 
( const double value ) 
{
  Ostap::Math::scale ( m_pars , value ) ;
  return *this ;
}
// ============================================================================
// add/shift  it
// ============================================================================
Ostap::Math::Rational&
Ostap::Math::Rational::add  
( const double value ) 
{
  Ostap::Math::shift ( m_pars , value ) ;
  return *this ;
}
// ============================================================================
// negate it! 
// ============================================================================
Ostap::Math::Rational
Ostap::Math::Rational::operator- () const 
{
  Rational result { *this } ;
  Ostap::Math::negate ( result.m_pars ) ;
  return result ;
}
// ============================================================================
// For python 
// ============================================================================
Ostap::Math::Rational
Ostap::Math::Rational::__add__     ( const double value ) const { return (*this)+value ; }
// ============================================================================
Ostap::Math::Rational
Ostap::Math::Rational::__sub__     ( const double value ) const { return (*this)-value ; }
// ============================================================================
Ostap::Math::Rational
Ostap::Math::Rational::__mul__     ( const double value ) const { return (*this)-value ; }
// ============================================================================
Ostap::Math::Rational
Ostap::Math::Rational::__div__     ( const double value ) const { return (*this)/value ; }
// ============================================================================
Ostap::Math::Rational
Ostap::Math::Rational::__truediv__ ( const double value ) const { return (*this)/value ; }
// ============================================================================
Ostap::Math::Rational
Ostap::Math::Rational::__radd__    ( const double value ) const { return (*this)+value ; }
// ============================================================================
Ostap::Math::Rational
Ostap::Math::Rational::__rsub__    ( const double value ) const { return value - (*this) ; }
// ============================================================================
Ostap::Math::Rational
Ostap::Math::Rational::__rmul__    ( const double value ) const { return (*this)*value ; }
// ============================================================================
Ostap::Math::Rational
Ostap::Math::Rational::__neg__     () const { return -(*this) ; }
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
{ return 
    x < xmin () ? 0.0 : 
    x > xmax () ? 0.0 : 
    ( m_p ( x ) / m_q ( x ) ) / ( m_p.xmax() - m_p.xmin () ) ; }
// ============================================================================
// all parameters (by value!!!)
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
      s_message                 , 
      __FILE__ , __LINE__       ,
      GSL_INTEG_GAUSS21         ) ;
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
// scale it
// ============================================================================
Ostap::Math::RationalBernstein&
Ostap::Math::RationalBernstein::scale 
( const double value ) { m_p *= value ; return *this ; }
// ============================================================================
// add/shift  it
// ============================================================================
Ostap::Math::RationalBernstein&
Ostap::Math::RationalBernstein::add  
( const double value ) { m_p += value ; return *this ; }
// ============================================================================
// negate it! 
// ============================================================================
Ostap::Math::RationalBernstein
Ostap::Math::RationalBernstein::operator- () const 
{
  RationalBernstein result { *this } ;
  result.m_p *= -1 ;
  return result ;
}
// ============================================================================
// Multiply with Bernstein polynomial
// ============================================================================
Ostap::Math::RationalBernstein&
Ostap::Math::RationalBernstein::operator *=  
( const Ostap::Math::Bernstein& right ) 
{
  if ( s_equal ( xmin () , right.xmin () ) && 
       s_equal ( xmax () , right.xmax () ) ) 
  { m_p = m_p * right ; }
  else 
  { m_p = m_p * Bernstein ( right , xmin () , xmax () ) ; }
  //
  return *this ;
}
// ============================================================================
// Add Bernstein polynomial
// ============================================================================
Ostap::Math::RationalBernstein&
Ostap::Math::RationalBernstein::operator +=  
( const Ostap::Math::Bernstein& right ) 
{
  if ( s_equal ( xmin () , right.xmin () ) && s_equal ( xmax () , right.xmax () ) ) 
  { 
    m_p = m_p + m_q.bernstein() * right  ;
  }
  else 
  { m_p = m_p * Bernstein ( m_q.bernstein() * right , xmin () , xmax () ) ; }
  //
  return *this ;
}
// ============================================================================
// subtract Bernstein polynomial
// ============================================================================
Ostap::Math::RationalBernstein&
Ostap::Math::RationalBernstein::operator -=  
( const Ostap::Math::Bernstein& right ) 
{ return ( *this)+= ( -right ) ; }
// ============================================================================
// For python 
// ============================================================================
Ostap::Math::RationalBernstein
Ostap::Math::RationalBernstein::__add__     ( const double value ) const { return (*this)+value ; }
// ============================================================================
Ostap::Math::RationalBernstein
Ostap::Math::RationalBernstein::__sub__     ( const double value ) const { return (*this)-value ; }
// ============================================================================
Ostap::Math::RationalBernstein
Ostap::Math::RationalBernstein::__mul__     ( const double value ) const { return (*this)-value ; }
// ============================================================================
Ostap::Math::RationalBernstein
Ostap::Math::RationalBernstein::__div__     ( const double value ) const { return (*this)/value ; }
// ============================================================================
Ostap::Math::RationalBernstein
Ostap::Math::RationalBernstein::__truediv__ ( const double value ) const { return (*this)/value ; }
// ============================================================================
Ostap::Math::RationalBernstein
Ostap::Math::RationalBernstein::__radd__    ( const double value ) const { return (*this)+value ; }
// ============================================================================
Ostap::Math::RationalBernstein
Ostap::Math::RationalBernstein::__rsub__    ( const double value ) const { return value - (*this) ; }
// ============================================================================
Ostap::Math::RationalBernstein
Ostap::Math::RationalBernstein::__rmul__    ( const double value ) const { return (*this)*value ; }
// ============================================================================
Ostap::Math::RationalBernstein
Ostap::Math::RationalBernstein::__neg__     () const { return -(*this) ; }
// ============================================================================
Ostap::Math::RationalBernstein
Ostap::Math::RationalBernstein::__add__     
( const Ostap::Math::Bernstein& right ) const { return (*this)+right ; }
// ============================================================================
Ostap::Math::RationalBernstein
Ostap::Math::RationalBernstein::__radd__     
( const Ostap::Math::Bernstein& right ) const { return (*this)+right ; }
// ============================================================================
Ostap::Math::RationalBernstein
Ostap::Math::RationalBernstein::__mul__     
( const Ostap::Math::Bernstein& right ) const { return (*this)*right ; }
// ============================================================================
Ostap::Math::RationalBernstein
Ostap::Math::RationalBernstein::__rmul__     
( const Ostap::Math::Bernstein& right ) const { return (*this)*right ; }
// ============================================================================
Ostap::Math::RationalBernstein
Ostap::Math::RationalBernstein::__sub__     
( const Ostap::Math::Bernstein& right ) const { return (*this)-right ; }
// ============================================================================
Ostap::Math::RationalBernstein
Ostap::Math::RationalBernstein::__rsub__     
( const Ostap::Math::Bernstein& right ) const { return right - (*this) ; }
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
      s_message                 , 
      __FILE__ , __LINE__       , 
      GSL_INTEG_GAUSS21         ) ;
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
// Rational Pade-like function
// ============================================================================
namespace 
{
  // =========================================================================
  /// helper function to "merge" two sequence of coefficients
  inline std::vector<double>
  pq_pars
  ( const std::vector<double>& p ,
    const std::vector<double>& q )
  {
    std::vector<double> pars{};
    if ( p.empty() )
      {
	pars.reserve   ( 1 + q.size() ) ;
	pars.push_back ( 1 ) ;               // ATTENTION! 
      }
    else
      {
	pars.reserve ( p.size() + q.size() ) ;
	pars.insert ( pars.end() , p.begin() , p.end() ) ;
      }
    pars.insert ( pars.end() , q.begin() , q.end() ) ;
    return pars ;
  }
  // =========================================================================
} // =========================================================================
// ===========================================================================
/*  simplified constructor
 *  @param pars list of "Pade"-parameters 
 *  @param n degree of P(x)
 *  @param xmin   low-x 
 *  @param xmax   high-x 
 */
// ============================================================================
Ostap::Math::Pade::Pade
( const std::vector<double>& pars ,
  const unsigned short       n    ,
  const double               xmin ,
  const double               xmax )
  : Pade ( pars , n , {} , {} , {} , {} , xmin , xmax )
{}
// ============================================================================
/*  simplified constructor
 *  @param pars list of "Pade"-parameters 
 *  @param n degree of P(x)
 *  @param zeros  list of real constituent zeroes 
 *  @param poles  list of real constituent poles         
 *  @param xmin   low-x 
 *  @param xmax   high-x 
 */
// ============================================================================
Ostap::Math::Pade::Pade
( const std::vector<double>&                pars   ,
  const unsigned short                      n      ,
  const std::vector<double>&                zeroes ,
  const std::vector<double>&                poles  ,	
  const double                              xmin   ,
  const double                              xmax   )
  : Pade ( pars , n , zeroes , poles, {} , {} , xmin , xmax )
{}  
// ============================================================================
/* full constructor
 *  @param pars list of "Pade"-parameters 
 *  @param n degree of P(x)
 *  @param zeros  list of real constituent zeroes 
 *  @param poles  list of real constituent poles         
 *  @param czeros (half) list of complex constituent zeroes 
 *  @param cpoles (half) list of complex constituent poles  
 *  @param xmin   low-x 
 *  @param xmax   high-x 
 */
// ============================================================================
Ostap::Math::Pade::Pade
( const std::vector<double>&                pars    ,
  const unsigned short                      n       ,
  const std::vector<double>&                zeroes  ,
  const std::vector<double>&                poles   ,	
  const std::vector<std::complex<double> >& czeroes ,
  const std::vector<std::complex<double> >& cpoles  ,
  const double                              xmin    ,
  const double                              xmax    ) 
  : Ostap::Math::Parameters ( std::max ( pars.size() , std::size_t ( n + 1 ) ) ) 
  , m_n       ( n )
  , m_m       ( 0 )
  , m_xmin    ( std::min ( xmin , xmax)      )
  , m_xmax    ( std::max ( xmin , xmax)      )
  , m_x0      ( 0.5 * ( xmin + xmax )        )
  , m_scale   ( 2 / std::abs ( xmax - xmin ) )
  , m_zeroes  ( zeroes  )
  , m_poles   ( poles   )
  , m_czeroes ( czeroes )
  , m_cpoles  ( cpoles  )
  , m_pnts    ( poles   ) 
{
  /// degree of Q
  m_m = npars() - m_n - 1 ;
  /// set parameters 
  setPars ( pars ) ;
  //
  auto cless = [] ( const std::complex<double>& a ,
		    const std::complex<double>& b ) -> bool
  { return std::real ( a ) < std::real ( b ) ; } ;
  //
  std::sort        ( m_poles .begin() , m_poles .end() ) ; 
  std::sort        ( m_zeroes.begin() , m_zeroes.end() ) ;
  std::stable_sort ( m_cpoles .begin() , m_cpoles .end() , cless ) ;
  std::stable_sort ( m_czeroes.begin() , m_czeroes.end() , cless ) ;     
  //
  for ( const auto z : m_cpoles )
    {
      const std::complex<double> tz = m_scale * ( z - m_x0 ) ;
      if ( std::abs ( std::imag ( tz ) ) < 0.02 )
	{ m_pnts.push_back ( std::real ( z ) ) ; }
    } 
  //
  std::sort ( m_pnts   .begin() , m_pnts   .end () ) ;
  m_pnts.erase ( std::unique( m_pnts.begin(), m_pnts.end() ), m_pnts.end() );
  //
}
// ==========================================================================
/*  simplified constructor
 *  @param ps list of P(x) -parameters, if empty interpeted as  <code>[ 1 ]</code>
 *  @param qs list of Q(x) -parameters 
 *  @param xmin   low-x 
 *  @param xmax   high-x 
 *  @attention if ps is empty it is interpreted as <code>[1]</code>
 */
// ==========================================================================
Ostap::Math::Pade::Pade
( const std::vector<double>& ps    ,
  const std::vector<double>& qs    ,
  const double               xmin  ,
  const double               xmax  )
  : Pade ( ::pq_pars ( ps , qs ) ,
	   ps.empty() ? 0 : ps.size() - 1 ,
	   xmin , xmax )
{}
// ==========================================================================
/*  simplified constructor
 *  @param ps list of P(x) -parameters, if empty interpeted as  <code>[ 1 ]</code>
 *  @param qs list of Q(x) -parameters 
 *  @param zeros  list of real constituent zeroes 
 *  @param poles  list of real constituent poles         
 *  @param xmin   low-x 
 *  @param xmax   high-x 
 *  @attention if ps is empty it is interpreted as <code>[1]</code>
 */
// ==========================================================================
Ostap::Math::Pade::Pade
( const std::vector<double>& ps     ,
  const std::vector<double>& qs     ,
  const std::vector<double>& zeroes ,
  const std::vector<double>& poles  ,	
  const double               xmin   ,
  const double               xmax   )
  : Pade ( ::pq_pars ( ps , qs ) ,
	   ps.empty() ? 0 : ps.size() - 1 ,
	   zeroes , poles , xmin , xmax )
{}
// ===========================================================================
/*  full constructor
 *  @param ps list of P(x) -parameters, if empty interpeted as  <code>[ 1 ]</code>
 *  @param qs list of Q(x) -parameters 
 *  @param zeros  list of real constituent zeroes 
 *  @param poles  list of real constituent poles         
 *  @param czeros (half) list of complex constituent zeroes 
 *  @param cpoles (half) list of complex constituent poles  
 *  @param xmin   low-x 
 *  @param xmax   high-x
 *  @attention if ps is empty it is interpreted as <code>[1]</code>
 */
// ===========================================================================
Ostap::Math::Pade::Pade
( const std::vector<double>&                ps      ,
  const std::vector<double>&                qs      ,
  const std::vector<double>&                zeroes  ,
  const std::vector<double>&                poles   ,	
  const std::vector<std::complex<double> >& czeroes ,
  const std::vector<std::complex<double> >& cpoles  ,
  const double                              xmin    ,
  const double                              xmax    ) 
  : Pade ( ::pq_pars ( ps , qs ) ,
	   ps.empty() ? 0 : ps.size() - 1 ,
	   zeroes , poles , czeroes , cpoles , xmin , xmax )
{}

// ============================================================================
// evaluate the function
// ============================================================================
double Ostap::Math::Pade::evaluate ( const double x ) const
{
  const double tx = t ( x ) ;
  //
  double result =  1 ;
  //
  // (1) get all zeroes
  if ( !m_zeroes.empty() || !m_czeroes.empty() ) { result *= Zt ( tx ) ; }
  // (2) get P&Q
  result *= Pt ( tx ) / Qt ( tx ) ;
  // (3) get all poles 
  if ( !m_poles.empty() || !m_cpoles.empty()   ) { result /= Rt ( tx ) ; }
  //
  return result ;
}
// ============================================================================
// get the integral
// ============================================================================
double Ostap::Math::Pade::integral 
(  const double xlow  , 
   const double xhigh ) const 
{
  if      ( s_equal ( xlow , xhigh ) ) { return 0 ; }
  else if ( xhigh <  xlow            ) { return - integral ( xhigh , xlow ) ; }
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Pade> s_integrator ;
  static const char s_message[] = "Integral(Pade)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  //
  // potental poles ? 
  if ( 0 < m_pnts.size () )
    {
      std::vector<double>::const_iterator il = std::lower_bound ( m_pnts.begin() , m_pnts.end(), xlow  ) ;
      std::vector<double>::const_iterator ih = std::upper_bound ( il             , m_pnts.end(), xhigh ) ;
      // points in the interval 
      if ( m_pnts.end() != il && il != ih ) 
	{
	  std::vector<double> points ( il , ih ) ;
	  std::tie ( ierror , result , error ) = s_integrator.qagp_integrate
	    ( tag  () , 
	      &F      , 
	      xlow    , xhigh ,           // low & high edges
	      points                    , // poles/special points 
	      workspace ( m_workspace ) , // workspace
	      s_APRECISION              , // absolute precision
	      s_RPRECISION              , // relative precision
	      m_workspace.size()        , // size of workspace
	      s_message                 , 
	      __FILE__ , __LINE__       ) ;
	  //
	  return result ;                                             // RETURN
	}
    }
  // regular case 
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag  () , 
      &F      , 
      xlow    , xhigh ,           // low & high edges
      workspace ( m_workspace ) , // workspace
      s_APRECISION              , // absolute precision
      s_RPRECISION              , // relative precision
      m_workspace.size()        , // size of workspace
      s_message                 , 
      __FILE__ , __LINE__       ) ;
  //
  return result ;                                                    // RETURN
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Pade::tag () const 
{ 
  static const std::string s_name = "Pade" ;
  return Ostap::Utils::hash_combiner 
    ( s_name   ,
      n     () ,
      m     () ,
      xmin  () ,
      xmax  () , 
      Ostap::Utils::hash_range ( zeroes  () ) , 
      Ostap::Utils::hash_range ( poles   () ) , 
      Ostap::Utils::hash_range ( czeroes () ) , 
      Ostap::Utils::hash_range ( cpoles  () ) , 
      Ostap::Utils::hash_range ( pars    () ) ) ;
}
// =============================================================================
// get value of P(x)
// =============================================================================
double Ostap::Math::Pade::Pt ( const double tx ) const
{
  typedef std::vector<double>::const_reverse_iterator CRI ;
  //
  const CRI b { CRI ( m_pars.begin() + ( m_n + 1 ) ) } ;
  const CRI e { CRI ( m_pars.begin()               ) } ;
  //
  /// calculate P(t)
  return Ostap::Math::Clenshaw::monomial_sum ( b , e , tx ) .first ;
}
// =============================================================================
// get value of P(x)
// =============================================================================
double Ostap::Math::Pade::Qt ( const double tx ) const
{
  typedef std::vector<double>::const_reverse_iterator CRI ;
  //
  const CRI b { CRI ( m_pars.begin() + ( m_n + 1 + m_m ) ) } ;
  const CRI e { CRI ( m_pars.begin() + ( m_n + 1       ) ) } ;
  //
  /// calculate Q(t)
  return 1.0 + tx * Ostap::Math::Clenshaw::monomial_sum ( b , e , tx ) .first ;
}
// ============================================================================
// get value of all zeroes Z(x)
// ============================================================================
double Ostap::Math::Pade::Zt ( const double tx ) const
{
  /// z -> t 
  auto t_z = [this] ( const std::complex<double>& z ) -> std::complex<double>
    { return this->m_scale * ( z - this->m_x0 ) ; } ;
  //
  double result = 1 ;
  /// (1) loop over poles 
  for ( const auto z : m_zeroes  ) { result *=           ( tx - t   ( z ) ) ; }
  /// (2) loop over pairs of complex-conjugated zeroes  
  for ( const auto z : m_czeroes ) { result *= std::norm ( tx - t_z ( z ) ) ; }
  //
  return result ;
}
// ============================================================================
// get value of all poles Rt(t)
// ============================================================================
double Ostap::Math::Pade::Rt ( const double tx ) const
{
  /// z -> t 
  auto t_z = [this] ( const std::complex<double>& z ) -> std::complex<double>
    { return this->m_scale * ( z - this->m_x0 ) ; } ;
  //
  double result = 1 ;  
  /// (1) loop over poles 
  for ( const auto z : m_poles  ) { result *=           ( tx - t   ( z ) ) ; }
  /// (2) loop over pairs of complex-conjugated poles 
  for ( const auto z : m_cpoles ) { result *= std::norm ( tx - t_z ( z ) ) ; }
  //
  return result ;
}
// ============================================================================
// swap two Pade functions 
// ============================================================================
void Ostap::Math::Pade::swap ( Ostap::Math::Pade& right )
{
  //
  Ostap::Math::Parameters::swap ( right   )  ;
  //
  std::swap ( m_n       , right.m_n       ) ;
  std::swap ( m_m       , right.m_m       ) ;
  std::swap ( m_xmin    , right.m_xmin    ) ;
  std::swap ( m_xmax    , right.m_xmax    ) ;
  std::swap ( m_x0      , right.m_x0      ) ;
  std::swap ( m_scale   , right.m_scale   ) ;
  std::swap ( m_zeroes  , right.m_zeroes  ) ;
  std::swap ( m_poles   , right.m_poles   ) ;
  std::swap ( m_czeroes , right.m_czeroes ) ;
  std::swap ( m_cpoles  , right.m_cpoles  ) ;
  std::swap ( m_pnts    , right.m_pnts    ) ;
  //
  Ostap::Math::swap ( m_workspace , right.m_workspace  ) ;
}
// ==========================================================================
/*  Interpolatory constructor 
 *  @param table interpolation table 
 *  @param n degree of P(x)
 *  @param zeros  list of real constituent zeroes 
 *  @param poles  list of real constituent poles         
 *  @param czeros (half) list of complex constituent zeroes 
 *  @param cpoles (half) list of complex constituent poles  
 */
// ==========================================================================
Ostap::Math::Pade::Pade 
( const Ostap::Math::Interpolation::Table&  table   ,
  const unsigned short                      n       ,
  const std::vector<double>&                zeroes  ,
  const std::vector<double>&                poles   ,	
  const std::vector<std::complex<double> >& czeroes ,
  const std::vector<std::complex<double> >& cpoles  )
  : Pade ( std::vector<double>( table.size() , 0.0 )              ,
	   n             , 
	   zeroes        ,
	   poles         ,
	   czeroes       ,
	   cpoles        ,
	   table.xmin () ,
	   table.xmax () )
{
  //
  Ostap::Assert ( !table.empty()               ,
		  "Empty interpolation table!" ,
		  "Ostap::Math::Pade"          ) ;
  //
  Ostap::Assert ( m_n + 1 <= table.size ()       ,
		  "Invalid size of interpolation table!" ,
		  "Ostap::Math::Pade"          ) ;
  //
  const unsigned short N = table.size() ;
  Ostap::Assert ( npars() == N  ,
		  "Mismatch interpolation table size/#pars!" ,
		  "Ostap::Math::Pade"          ) ;
  //
  // =====================================================================  
  Ostap::GSL::Matrix A { N , N , Ostap::GSL::Matrix::Zero() } ;  
  Ostap::GSL::Vector b { N     } ; // free column 
  // (2) fill the matrix A and free column b
  for ( unsigned short j = 0 ; j < N ; ++j )
    {
      const double x = table  .x ( j  ) ;
      const double y = table  .y ( j  ) ;
      const double t = this -> t ( x  ) ;
      
      const double Z = Zt ( t ) ; // all zeroes 
      const double R = Rt ( t ) ; // all poles 
      
      // =================================================================
      { // the first n + 1 columns  
        double xx = Z ;
        for ( unsigned int i = 0 ; i < m_n + 1 ; ++i ) 
          {
            A.set ( j , i , xx ) ;
            xx *= t ;	   
          }	
      }
      // =================================================================
      { // remaining columns 
        double xx = -R * y * t ;
        for ( unsigned int i = m_n + 1 ; i < N ; ++i ) 
          {
            A.set ( j , i , xx ) ;
            xx *= t ;	   
          }
      }
      // free column
      b.set ( j , R * y ) ;
    }
  // =====================================================================  
  
  // (3) solve the system Ax=b using LU decomposition with pivoting 
  
  // (3.1) make LU decomposition with pivoting 
  Ostap::GSL::Permutation  P { N } ;  
  int signum ; 
  int ierror  = gsl_linalg_LU_decomp ( A.matrix() , P.permutation() , &signum ) ;
  if ( ierror ) { gsl_error ( "Failure in LU-decomposition" , __FILE__  , __LINE__ , ierror ) ; }
  Ostap::Assert ( !ierror                        ,
                  "Failure in LU-decomposition!" ,
                  "Ostap::Math::Pade"            , 1100 + ierror ) ;
  
  // (3.2) solve the system Ax=b 
  Ostap::GSL::Vector       x { N     } ; // solution 
  ierror  = gsl_linalg_LU_solve ( A.matrix(), P.permutation() , b.vector() , x.vector() );
  if ( ierror ) { gsl_error ( "Failure in LU-solve" , __FILE__  , __LINE__ , ierror ) ; }
  Ostap::Assert ( !ierror                ,
                  "Failure in LU-solve!" ,
                  "Ostap::Math::Pade"    , 1200 + ierror ) ;
  
  // (4) Feed Pade with calculated parameters 
  for ( unsigned short k = 0 ; k < N ; ++k ) { setPar ( k , x.get ( k ) ) ; }
}

// ===========================================================================
/** constructor from the polynomial expansion 
 *  @param n degree of P(x) 
 *  @param m degree of Q(x) 
 *  @param f polynomial expansion 
 */
// ===========================================================================
Ostap::Math::Pade::Pade 
( const unsigned short       n    , 
  const unsigned short       m    , 
  const std::vector<double>& f    ,
  const double               xmin ,
  const double               xmax )
  : Pade ( std::vector<double> ( m + n + 1 , 0.0 ) ,
	   n , {} , {} , {} , {} , xmin , xmax )
{
  Ostap::Assert ( n + m + 1 <= f.size ()               ,
		  "Invalid Polynomial->Pade setting!"  , 
		  "Ostap::Math::Pade"                  ) ;
  //
  const unsigned short N = this->npars () ;
  //
  typedef std::vector<double>::const_iterator         CI  ;
  typedef std::vector<double>::const_reverse_iterator CRI ;
  //
  Ostap::GSL::Matrix A { N , N , Ostap::GSL::Matrix::Zero() } ;
  Ostap::GSL::Vector b { N     } ; // free column 
  //
  for ( unsigned short j = 0 ; j < N ; ++j )
    {
      if ( j < n + 1 ) { A.set ( j , j , 1 ) ; }
      //
      if ( 1 <= j && j < m + 1 )
	{
	  CI  i1 = f.begin()       ;
	  CRI i2 { f.begin() + j } ;
	  for ( unsigned short k = 0 ; k < j ; ++k )
	    {	     
	      A.set ( j , n + 1 + k , -1 * ( *i2++ ) ) ;
	    }
	}
      if ( m + 1 <= j )
	{
	  CRI i1 { f.begin() + ( j - m ) } ;
	  CRI i2 { f.begin() +   j       } ;	  
	  for ( unsigned short k = 0 ; k < m ; ++k )
	    {	     
	      A.set ( j , n + 1 + k , -1 * ( *i2++ ) ) ;
	    }
	}
      b.set( j , f[j] ) ;
    }
  // ===============================================================================
  // (3) solve the system Ax=b using LU decomposition with pivoting 
  
  // (3.1) make LU decomposition with pivoting 
  Ostap::GSL::Permutation  P { N } ;  
  int signum ; 
  int ierror  = gsl_linalg_LU_decomp ( A.matrix() , P.permutation() , &signum ) ;
  if ( ierror ) { gsl_error ( "Failure in LU-decomposition" , __FILE__  , __LINE__ , ierror ) ; }
  Ostap::Assert ( !ierror ,
		  "Failure in LU-decomposition!" ,
		  "Ostap::Math::Pade"            , 1100 + ierror ) ;
  
  // (3.2) solve the system Ax=b 
  Ostap::GSL::Vector       x { N     } ; // solution 
  ierror  = gsl_linalg_LU_solve ( A.matrix(), P.permutation() , b.vector() , x.vector() );
  if ( ierror ) { gsl_error ( "Failure in LU-solve" , __FILE__  , __LINE__ , ierror ) ; }
  Ostap::Assert ( !ierror                ,
		  "Failure in LU-solve!" ,
		  "Ostap::Math::Pade"    , 1200 + ierror ) ;
  
  // (4) Feed Pade with calculated parameters 
  for ( unsigned short k = 0 ; k < N ; ++k ) { setPar ( k , x.get ( k ) ) ; }
  // ==========================================================================
}
// ============================================================================
/*  constructor from polynomi expansion 
 *  @param p polynomial expansion 
 *  @param n degree of P(x) 
 */
// =============================================================================
Ostap::Math::Pade::Pade 
( const Ostap::Math::Polynomial& p ,
  const unsigned short           n ,
  const unsigned short           m )
  : Pade ( n , m , p.pars() , p.xmin() , p.xmax() )
{}
// ===============================================================================
/* constructor from the polynomial expansion 
 *  @param p polynomial expansion 
 *  @param n degree of P(x) 
 */
// ===============================================================================
Ostap::Math::Pade::Pade 
( const Ostap::Math::Polynomial& p ,
  const unsigned short           n )
  : Pade ( p , n , n + 1 <= p.npars() ? p.npars() - n - 1 : 0 )
{
  Ostap::Assert ( n  + 1 <= p.npars ()                ,
		  "Invalid Poliynomial->Pade setting!" , 
		  "Ostap::Math::Pade"                  ) ;
}



/** 
    
#include "Ostap/GSL_utils.h"
#include "GSL_helpers.h"
#include <iostream>

// Ostap::Math::Pade
void Ostap::Math::pade_ququ
( const std::vector<double>& f ,
  const unsigned short       n )
{
  Ostap::Assert ( n + 1 <= f.size() ,
		  "Insufficient array of coefficients!"
		  "Ostap::Math::pade" ) ;
  
  //
  const unsigned short N = f.size() ;
  const unsigned short m = f.size() - n - 1 ;
  //
  typedef std::vector<double>::const_iterator         CI ;
  typedef std::vector<double>::const_reverse_iterator CRI;
  //
  Ostap::GSL::Matrix A { N , N , Ostap::GSL::Matrix::Zero() } ;
  Ostap::GSL::Vector b { N     } ; // free column 
  //
  for ( unsigned short j = 0 ; j < N ; ++j )
    {
      if ( j < n + 1 ) { A.set ( j , j , 1 ) ; }
      //
      if ( 1 <= j && j < m + 1 )
	{
	  CI  i1 = f.begin()       ;
	  CRI i2 { f.begin() + j } ;
	  for ( unsigned short k = 0 ; k < j ; ++k )
	    {	     
	      A.set ( j , n + 1 + k , -1 * ( *i2++ ) ) ;
	    }
	}
      if ( m + 1 <= j )
	{
	  CRI i1 { f.begin() + ( j - m ) } ;
	  CRI i2 { f.begin() +   j       } ;	  
	  for ( unsigned short k = 0 ; k < m ; ++k )
	    {	     
	      A.set ( j , n + 1 + k , -1 * ( *i2++ ) ) ;
	    }
	}
      b.set( j , f[j] ) ;
    }
  // solve Ax=b
    
  // (3) solve the system Ax=b using LU decomposition with pivoting 
  
  // (3.1) make LU decomposition with pivoting 
  Ostap::GSL::Permutation  P { N } ;  
  int signum ; 
  int ierror  = gsl_linalg_LU_decomp ( A.matrix() , P.permutation() , &signum ) ;
  if ( ierror ) { gsl_error ( "Failure in LU-decomposition" , __FILE__  , __LINE__ , ierror ) ; }
  Ostap::Assert ( !ierror ,
		  "Failure in LU-decomposition!" ,
		  "Ostap::Math::Pade"            , 1100 + ierror ) ;
  
  // (3.2) solve the system Ax=b 
  Ostap::GSL::Vector       x { N     } ; // solution 
  ierror  = gsl_linalg_LU_solve ( A.matrix(), P.permutation() , b.vector() , x.vector() );
  if ( ierror ) { gsl_error ( "Failure in LU-solve" , __FILE__  , __LINE__ , ierror ) ; }
  Ostap::Assert ( !ierror                ,
		  "Failure in LU-solve!" ,
		  "Ostap::Math::Pade"    , 1200 + ierror ) ;
  
  // (4) Feed Pade with calculated parameters 
  // for ( unsigned short k = 0 ; k < N ; ++k ) { setPar ( k , x.get ( k ) ) ; }
  
  
  std::cout << (*A.matrix()) << std::endl ;
  std::cout << (*x.vector()) << std::endl ;
}

*/

// ============================================================================
//                                                                      The END 
// ============================================================================
