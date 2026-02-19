// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <cmath>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Hash.h"
#include "Ostap/Bernstein1D.h"
#include "Ostap/PhaseSpace.h"
#include "Ostap/MoreMath.h"
#include "Ostap/AdHocShapes.h"
// ============================================================================
//  Local 
// ============================================================================
#include "local_math.h"
#include "local_hash.h"
#include "status_codes.h"
#include "Integrator1D.h"
// ============================================================================
/** @file
 *  Implementation file for functions from the file Ostap/Models.h
 *
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2010-04-19
 */
// ============================================================================

// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::ExpoPositive::ExpoPositive
( const unsigned short      N    ,
  const double              tau  ,   
  const double              xmin ,
  const double              xmax )
  : Ostap::Math::PolyFactor1D ( N , xmin , xmax )
  , m_tau       ( tau ) 
{}
// ============================================================================
// constructor from the order
// ============================================================================
Ostap::Math::ExpoPositive::ExpoPositive
( const std::vector<double>& pars ,
  const double               tau  ,
  const double               xmin ,
  const double               xmax )
  : Ostap::Math::PolyFactor1D ( pars , xmin , xmax )
  , m_tau       ( tau ) 
{}
// ======================================================================
// constructor from polynom and exponential 
// ======================================================================
Ostap::Math::ExpoPositive::ExpoPositive
( const Ostap::Math::Positive& pol , 
  const double                 tau ) 
  : Ostap::Math::PolyFactor1D ( pol )
  , m_tau ( tau ) 
{}
// ============================================================================
// set tau-parameter
// ============================================================================
bool Ostap::Math::ExpoPositive::setTau ( const double value )
{
  if ( s_equal ( value , m_tau ) ) { return false ; }
  m_tau = value ;
  return true ;
}
// ============================================================================
// get the value
// ============================================================================
double Ostap::Math::ExpoPositive::evaluate ( const double x ) const 
{
  if ( x < xmin() || x > xmax() ) { return 0 ; }
  return my_exp ( m_tau * x ) * m_positive ( x ) ;
}
// ============================================================================
double Ostap::Math::ExpoPositive::integral () const 
{ return integral ( xmin () , xmax() )  ; }
// ============================================================================
double Ostap::Math::ExpoPositive::integral
( const double low  , 
  const double high ) const 
{ return Ostap::Math::integrate ( m_positive.bernstein() , m_tau , low , high ) ; }
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::ExpoPositive::tag () const 
{
  static const std::string s_name = "ExpoPositive" ;
  return Ostap::Utils::hash_combiner ( s_name , m_positive.tag () , m_tau ) ; 
}
// ================-===========================================================
// get the value \f$ x_{min}\$ such that  \f$ x_{min} \le p(x) \f$ 
// ================-===========================================================
double Ostap::Math::ExpoPositive::min_value () const
{
  //
  const double pmin = m_positive.min_value () ;
  const double exp1 = my_exp ( m_tau * xmin () ) ; 
  const double exp2 = my_exp ( m_tau * xmax () ) ;
  //
  return pmin * std::min ( exp1 , exp2 ) ;
}
// ================-===========================================================
// get the value \f$ x_{max}\$ such that  \f$ x_{max} \ge p(x) \f$ 
// ================-===========================================================
double Ostap::Math::ExpoPositive::max_value () const
{
  //
  const double pmax = m_positive.max_value () ;
  const double exp1 = my_exp ( m_tau * xmin () ) ; 
  const double exp2 = my_exp ( m_tau * xmax () ) ;
  //
  return pmax * std::max ( exp1 , exp2 ) ;
}

// ============================================================================
/* constructor from threshold and number of particles
 *  @param threshold_L the low-mass  threshold
 *  @param l           how many particles we consider
 *  @param N           degree of polynomial
 *  @param tau         the exponent 
 *  @param xhigh       the high edge 
 */
// ============================================================================
Ostap::Math::PhaseSpaceLeftExpoPol::PhaseSpaceLeftExpoPol
( const double         threshold_L ,   // low threshold 
  const unsigned short l           ,   // number of particles 
  const unsigned short N           ,   // degree of polynomial
  const double         tau         ,   // the exponent 
  const double         xhigh       )   // high edge 
  : PhaseSpaceLeftExpoPol ( threshold_L , l , N , tau , threshold_L , xhigh ) 
{}
// ============================================================================
/*  constructor from the phase space and polynomial degree
 *  @param ps          phase space factor
 *  @param N           degree of polynomial
 *  @param tau         the exponent 
 *  @param xhigh       the high edge 
 */
// ============================================================================
Ostap::Math::PhaseSpaceLeftExpoPol::PhaseSpaceLeftExpoPol
( const Ostap::Math::PhaseSpaceLeft& ps ,
  const unsigned short  N     ,   // degree of polynomial
  const double          tau   ,   // the exponent 
  const double          xhigh )   // high edge 
  : PhaseSpaceLeftExpoPol ( ps , N , tau , ps.threshold () , xhigh ) 
{}
// ============================================================================
/* constructor from threshold and number of particles
 *  @param threshold_L the low-mass  threshold
 *  @param l           how many particles we consider
 *  @param N           degree of polynomial
 *  @param tau         the exponent 
 *  @param xlow        the low  edge 
 *  @param xhigh       the high edge 
 */
// ============================================================================
Ostap::Math::PhaseSpaceLeftExpoPol::PhaseSpaceLeftExpoPol
( const double         threshold_L ,   // low threshold 
  const unsigned short l           ,   // number of particles 
  const unsigned short N           ,   // degree of polynomial
  const double         tau         ,   // the exponent 
  const double         xlow        ,   // low edge 
  const double         xhigh       )   // high edge 
  : PhaseSpaceLeftExpoPol ( Ostap::Math::PhaseSpaceLeft ( threshold_L , l ) ,
                            N , tau , xlow , xhigh ) 
{}
// ============================================================================
/* constructor from the phase space and polynomial degree
 *  @param ps          phase space factor
 *  @param N           degree of polynomial
 *  @param tau         the exponent 
 *  @param xlow        the low  edge 
 *  @param xhigh       the high edge 
 */
// ============================================================================
Ostap::Math::PhaseSpaceLeftExpoPol::PhaseSpaceLeftExpoPol
( const PhaseSpaceLeft& ps    ,
  const unsigned short  N     ,   // degree of polynomial
  const double          tau   ,   // the exponent 
  const double          xlow  ,   // low edge 
  const double          xhigh )  // high edge
  : Ostap::Math::PolyFactor1D ( N  ,
				std::max ( ps.threshold () , std::min ( xlow , xhigh ) ) , 
				std::max ( xlow , xhigh ) ) 
  , m_phasespace  ( ps ) 
  , m_tau         ( std::abs ( tau ) ) 
  , m_workspace   ( ) 
{
  //
  Ostap::Assert ( m_phasespace.threshold() <= m_positive.xmin () , 
                  "Invalid setting of threshold/xmin/xmax" ,
                  "Ostap::Math::PhaseSpaceLeftPol"         ,
		  INVALID_PARAMETERS , __FILE__ , __LINE__ ) ; 
  //
}
// ============================================================================
/*  constructor from the phase space and polynomial
 *  @param ps          phase space factor
 *  @param poly        polynomial
 *  @param tau         the exponent 
 */
// ============================================================================
Ostap::Math::PhaseSpaceLeftExpoPol::PhaseSpaceLeftExpoPol
( const PhaseSpaceLeft&        ps  ,   // pjase space 
  const Ostap::Math::Positive& pol ,   // polynomial 
  const double                 tau )  // the exponent 
  : Ostap::Math::PolyFactor1D ( pol ) 
  , m_phasespace  ( ps  ) 
  , m_tau         ( std::abs ( tau ) )
  , m_workspace   ( ) 
{
  //
  Ostap::Assert ( m_phasespace.threshold() < m_positive.xmax () , 
                  "Invalid setting of threshold/xmin/xmax" ,
                  "Ostap::Math::PhaseSpaceLeftPol"         , 
		  INVALID_PARAMETERS , __FILE__ , __LINE__ ) ;
  //
}
// ======================================================================
// evaluate the modulated phase space
// ============================================================================
double Ostap::Math::PhaseSpaceLeftExpoPol::evaluate ( const double x ) const 
{
  if  ( x <= xmin () || x >= xmax () ) { return 0 ; }
  const double xc = 0.5 * ( xmin() + xmax() ) ;
  return  
    m_phasespace ( x  ) / 
    m_phasespace ( xc ) * 
    m_positive   ( x  ) *  
    std::exp     ( -1 * m_tau  * ( x - xc ) ) ;
}
// ============================================================================
// set tau-parameter
// ============================================================================
bool Ostap::Math::PhaseSpaceLeftExpoPol::setTau ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( avalue , m_tau ) ) { return false ; }
  m_tau = avalue ;
  return true ;
}
// ============================================================================
// get the integral between low and high limits
// =============================================================================
double Ostap::Math::PhaseSpaceLeftExpoPol::integral 
( const double low  ,
  const double high ) const 
{
  //
  if      ( s_equal ( low , high ) ) { return 0 ; } // return 
  else if ( high < low             ) { return -1 * integral ( high , low ) ; } // RETURN
  else if ( high <= xmin ()        ) { return 0 ; }
  else if ( low  >= xmax ()        ) { return 0 ; }
  //
  const double xlow  = std::max ( low  , xmin () ) ;
  const double xhigh = std::min ( high , xmax () ) ;
  //
  // if the exponent plays important role, split the interval 
  if ( !s_zero ( m_tau ) ) 
  {
    if  ( 3 < ( xhigh - xlow ) * m_tau )  
    {
      const double xc = 0.5 * ( xhigh + xlow ) ;
      return integral ( xlow , xc ) + integral ( xc , xhigh ) ;
    }
  }
  //
  /// split near-threshold region 
  const double delta =  xmax() - threshold() ;
  const double len   =  xhigh  -  xlow  ;
  const double x1 = threshold () + 0.05 * delta ;
  if ( 0.05 * delta < len && xlow < x1 && x1 < xhigh ) 
  { return integral ( xlow , x1 ) + integral ( x1 , xhigh ) ; }
  const double x2 = threshold () + 0.15 * delta ;
  if ( 0.10 * delta < len && xlow < x2 && x2 < xhigh ) 
  { return integral ( xlow , x2 ) + integral ( x2 , xhigh ) ; }
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<PhaseSpaceLeftExpoPol> s_integrator {} ;
  static char s_message[] = "Integral(PhaseSpaceLeftExpoPol)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag () ,  
      &F     , 
      xlow   , xhigh      ,          // low & high edges
      workspace ( m_workspace ) ,    // workspace
      s_APRECISION         ,          // absolute precision
      s_RPRECISION         ,          // relative precision
      m_workspace.size()              ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::PhaseSpaceLeftExpoPol::tag () const 
{
  static const std::string s_name = "PhaseSpaceLeftExpoPol" ;
  return Ostap::Utils::hash_combiner ( s_name , m_phasespace.tag () , m_positive.tag () , m_tau ) ;
}
// ============================================================================




// ============================================================================
// constructor from polynom and parameters "alpha" and "x0"
// ============================================================================
Ostap::Math::Sigmoid::Sigmoid
( const Ostap::Math::Positive&   poly  , 
  const double                   scale ,
  const double                   x0    , 
  const double                   delta , 
  const Ostap::Math::SigmoidType st    ) 
  : Ostap::Math::PolyFactor1D ( poly  )
  , m_scale     ( scale )
  , m_x0        ( x0    )
  , m_delta     ( 0     ) 
  , m_type      ( st    )
  , m_sin2delta ( 0     )
  , m_workspace () 
{
  setDelta ( delta ) ; 
  //
  Ostap::Assert ( !s_zero ( m_scale )                     ,
		  "Parameter `scale` must be non-zero!"   ,
		  "Ostap::Math::Sigmoid"                  ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //	
  Ostap::Assert ( Ostap::Math::SigmoidType::First <= st &&
		  st <= Ostap::Math::SigmoidType::Last    ,
		  "Invalid SigmoidType!"                  ,
		  "Ostap::Math::Sigmoid"                  ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
}
// ============================================================================
// constructor from polynom and parameters "alpha" and "x0"
// ============================================================================
Ostap::Math::Sigmoid::Sigmoid 
( const unsigned short           N     , 
  const double                   xmin  , 
  const double                   xmax  , 
  const double                   scale , 
  const double                   x0    , 
  const double                   delta , 
  const Ostap::Math::SigmoidType st    )   
  : Ostap::Math::PolyFactor1D ( N , xmin , xmax )
  , m_scale     ( scale )
  , m_x0        ( x0    )
  , m_delta     ( 0     ) 
  , m_type      ( st    )
  , m_sin2delta ( 0     ) 
  , m_workspace() 
{
  //
  setDelta ( delta );
  //
  Ostap::Assert ( m_scale                                 ,
		  "Parameter `scale` nust be non-zero!"   ,
		  "Ostap::Math::Sigmoid"                  ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  Ostap::Assert ( Ostap::Math::SigmoidType::First <= st &&
		  st <= Ostap::Math::SigmoidType::Last    ,
		  "Invalid SigmoidType!"  ,
		  "Ostap::Math::Sigmoid"  ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
}
// ============================================================================
// constructor from polynom and parameters "alpha" and "x0"
// ============================================================================
Ostap::Math::Sigmoid::Sigmoid
( const std::vector<double>&     pars  ,
  const double                   xmin  , 
  const double                   xmax  , 
  const double                   scale , 
  const double                   x0    ,  
  const double                   delta , 
  const Ostap::Math::SigmoidType st    )   
  : Ostap::Math::PolyFactor1D ( pars , xmin , xmax )
  , m_scale     ( scale )
  , m_x0        ( x0    )
  , m_delta     ( 0     ) 
  , m_type      ( st    )
  , m_sin2delta ( 0     ) 
  , m_workspace () 
{
  //
  setDelta ( delta ) ;
  //
  Ostap::Assert ( !s_zero ( m_scale )                     ,
		  "Parameter `scale` must be non-zero!"   ,
		  "Ostap::Math::Sigmoid"                  ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  Ostap::Assert ( Ostap::Math::SigmoidType::First <= st &&
		  st <= Ostap::Math::SigmoidType::Last    ,
		  "Invalid SigmoidType!"  ,
		  "Ostap::Math::Sigmoid"  ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
}
// ============================================================================
// constructor from polynomial and parameters "alpha" and "x0"
// ============================================================================
Ostap::Math::Sigmoid::Sigmoid
( const std::string&             sigmoid_name , 
  const Ostap::Math::Positive&   poly         ,
  const double                   scale        ,
  const double                   x0           ,
  const double                   delta        )
  : Sigmoid ( poly  ,
	      scale ,
	      x0    ,
	      delta ,
	      Ostap::Math::sigmoid_type ( sigmoid_name ) )
{}
// ============================================================================
// constructor from polynomial and parameter "alpha"
// ============================================================================
Ostap::Math::Sigmoid::Sigmoid
( const std::string&             sigmoid_name , 
  const unsigned short           N            ,
  const double                   xmin         ,
  const double                   xmax         ,
  const double                   scale        ,
  const double                   x0           , 
  const double                   delta        )
  : Sigmoid ( N , xmin , xmax , scale , x0 , delta ,
	      Ostap::Math::sigmoid_type ( sigmoid_name ) )
{}
// ============================================================================
// constructor from polynomial and parameters "alpha" and "x0"
// ============================================================================
Ostap::Math::Sigmoid::Sigmoid
( const std::string&             sigmoid_name , 
  const std::vector<double>&     pars         ,
  const double                   xmin         ,
  const double                   xmax         ,
  const double                   scale        ,
  const double                   x0           ,
  const double                   delta        )
  : Sigmoid ( pars , xmin , xmax , scale , x0 , delta ,
	      Ostap::Math::sigmoid_type ( sigmoid_name ) )
{}
// ============================================================================
// set new value for alpha 
// ============================================================================
bool Ostap::Math::Sigmoid::setScale ( const double value )
{
  if ( s_equal ( m_scale , value ) ) { return false ; }
  //
  Ostap::Assert ( !s_zero ( m_scale )                     ,
		  "Parameter `scale` must be non-zero!"   ,
		  "Ostap::Math::Sigmoid:setScale"         ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  m_scale = value ;
  //
  return true ;
}
// ============================================================================
// set new value for x0
// ============================================================================
bool Ostap::Math::Sigmoid::setX0 ( const double value )
{
  if ( s_equal ( m_x0, value ) ) { return false ; }
  m_x0 = value ;
  //
  return true ;
}
// ============================================================================
// set new value for delta 
// ============================================================================
bool Ostap::Math::Sigmoid::setDelta ( const double value )
{
  if ( s_equal ( m_delta , value ) ) { return false ; }
  m_delta     = value ;
  m_sin2delta = std::pow ( std::sin ( m_delta ) , 2 ) ;
  m_sin2delta = std::min ( m_sin2delta , 1.0 ) ; 
  //
  return true ;
}
// ============================================================================
// get the value \f$ x_{min}\$ such that  \f$ x_{min} \le p(x) \f$ 
// ================-===========================================================
double Ostap::Math::Sigmoid::min_value () const
{
  //
  const double c2 = cos2delta () ;
  const double s2 = sin2delta () ; 
  //
  const double pmin = m_positive.min_value () ;
  const double sig1 = c2 * sigmoid ( xmin () ) + s2 ;
  const double sig2 = c2 * sigmoid ( xmax () ) + s2 ;
  //
  return pmin * std::min ( sig1 , sig2 ) ;
}
// ================-===========================================================
// get the value \f$ x_{max}\$ such that  \f$ x_{max} \ge p(x) \f$ 
// ================-===========================================================
double Ostap::Math::Sigmoid::max_value () const
{
  //
  const double c2 = cos2delta () ;
  const double s2 = sin2delta () ; 
  //
  const double pmax = m_positive.max_value () ; 
  const double sig1 = c2 * sigmoid ( xmin () ) + s2 ;
  const double sig2 = c2 * sigmoid ( xmax () ) + s2 ;
  //
  return pmax * std::max ( sig1 , sig2 ) ;
}
// ============================================================================
// get the integral between xmin and xmax 
// ============================================================================
double Ostap::Math::Sigmoid::integral   () const 
{ return integral ( m_positive.xmin () , m_positive.xmax() ) ; }
// ============================================================================
// get the integral between low and high 
// ============================================================================
double Ostap::Math::Sigmoid::integral  
( const double low  , 
  const double high ) const 
{
  //
  if      ( high < low                ) { return -integral ( high , low ) ; }
  else if ( s_equal ( low , high )    ) { return 0 ; }
  else if ( high < xmin ()            ) { return 0 ; }
  else if ( low  > xmax ()            ) { return 0 ; }
  //
  // split it, if needed 
  if ( low < m_x0 && m_x0 < high ) 
    { return integral ( low , m_x0 ) + integral ( m_x0 , high ) ; }
  // split further, if needed 
  const double a3 = m_x0 + 3 * m_scale ;
  if ( low < a3 && a3 < high ) { return integral ( low , a3 ) + integral ( a3 , high ) ; }
  // split further, if needed  
  const double a4 = m_x0 - 3 * m_scale ;
  if ( low < a4 && a4 < high ) { return integral ( low , a4 ) + integral ( a4 , high ) ; }
  //
  static const Ostap::Math::GSL::Integrator1D<Sigmoid> s_integrator {} ;
  static char s_message[] = "Integral(Sigmoid)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag  () , 
      &F      , 
      low     , high ,               // low & high edges
      workspace ( m_workspace ) ,    // workspace
      s_APRECISION         ,          // absolute precision
      s_RPRECISION         ,          // relative precision
      m_workspace.size()   ,          // size of workspace
      s_message            , 
      __FILE__ , __LINE__  ) ;
  //
  return result ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Sigmoid::tag () const 
{
  static const std::string s_name = "Sigmoid" ;
  return Ostap::Utils::hash_combiner ( s_name            , 
				       m_positive.tag () ,
				       m_scale ,
				       m_x0    ,
				       m_type  , 
				       m_delta ) ; 
}
// ============================================================================
// the name of sigmoid function 
// ============================================================================
std::string Ostap::Math::Sigmoid::sigmoid_name () const
{ return Ostap::Math::sigmoid_name ( m_type ) ; }

// ============================================================================
Ostap::Math::TwoExpos::TwoExpos
( const double alpha ,
  const double delta , 
  const double x0    ) 
  : m_alpha ( std::abs ( alpha ) ) 
  , m_delta ( std::abs ( delta ) ) 
  , m_x0    ( x0 ) 
{}
// ============================================================================
// set new value for x0
// ============================================================================
bool Ostap::Math::TwoExpos::setX0 ( const double value )
{
  if ( s_equal ( m_x0, value ) ) { return false ; }
  m_x0 = value ;
  //
  return true ;
}
// ============================================================================
// set new value for alpha
// ============================================================================
bool Ostap::Math::TwoExpos::setAlpha ( const double value )
{
  const double nv = std::abs ( value ) ;
  if ( s_equal ( m_alpha, nv ) ) { return false ; }
  m_alpha = nv ;
  //
  return true ;
}
// ============================================================================
// set new value for delta
// ============================================================================
bool Ostap::Math::TwoExpos::setDelta ( const double value )
{
  const double nv = std::abs ( value ) ;
  if ( s_equal ( m_delta, nv ) ) { return false ; }
  m_delta = nv ;
  //
  return true ;
}
// ============================================================================
// get the value 
// ============================================================================
double Ostap::Math::TwoExpos::evaluate ( const double x ) const 
{ return x < m_x0 ? 0 : derivative ( x , 0 ) ; }
// ============================================================================
// get the integral between -inf and +inf
// ============================================================================
double Ostap::Math::TwoExpos::integral   () const { return 1 ; }
// ============================================================================
// get the integral between low and high 
// ============================================================================
double Ostap::Math::TwoExpos::integral   
( const double low  , 
  const double high ) const 
{
  //
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( low  > high            ) { return -integral ( high , low  ) ; }
  else if ( high <= m_x0           ) { return 0 ; }
  else if ( low  <  m_x0           ) { return  integral ( m_x0 , high ) ; }
  //
  const double a     = m_alpha            ;
  const double b     = m_alpha + m_delta  ;
  //
  const double xlow  = low  - m_x0 ;
  const double xhigh = high - m_x0 ;  
  //
  const double norm  = 1.0 / m_alpha - 1.0 / ( m_alpha + m_delta ) ;
  return 
    ( ( std::exp ( -b * xhigh ) - std::exp ( -b * xlow ) ) / b -
      ( std::exp ( -a * xhigh ) - std::exp ( -a * xlow ) ) / a ) / norm ;
}
// ============================================================================
namespace 
{
  // ==========================================================================
  inline unsigned long long _factorial_ ( const unsigned short N ) 
  {
    return 
      0 == N ?  1 : 
      1 == N ?  1 : 
      2 == N ?  2 : 
      3 == N ?  6 : 
      4 == N ? 24 : N * _factorial_ ( N - 1 ) ;
  }
  // ==========================================================================
  /// get (un-normalized) moment 
  inline long double _moment_ 
  ( const long double    alpha , 
    const long double    delta , 
    const unsigned short N     ) 
  {
    return _factorial_ ( N ) *  
      ( 1 / Ostap::Math::POW ( alpha         , N + 1 ) - 
        1 / Ostap::Math::POW ( alpha + delta , N + 1 ) ) ;  
  }
  // ==========================================================================
}
// ============================================================================
// get normalization constant
// ============================================================================
double Ostap::Math::TwoExpos::norm () const 
{ return 1.L / _moment_ ( m_alpha , m_delta , 0 ) ; } 
// ============================================================================
// mean-value (for -inf,+inf) interval 
// ============================================================================
double Ostap::Math::TwoExpos::mean  () const 
{
  const long double n0 = _moment_ ( m_alpha , m_delta , 0 ) ;
  const long double n1 = _moment_ ( m_alpha , m_delta , 1 ) ;
  //
  return m_x0 + n1 / n0 ;
}
// ============================================================================
// mode 
// ============================================================================
double  Ostap::Math::TwoExpos::mode  () const 
{
  const long double delta = m_delta ;
  return m_x0 + std::log1p ( delta / m_alpha ) / delta ; 
}
// ============================================================================
// variance 
// ============================================================================
double Ostap::Math::TwoExpos::variance () const 
{
  const long double n0 = _moment_ ( m_alpha , m_delta , 0 ) ;
  const long double n1 = _moment_ ( m_alpha , m_delta , 1 ) ;
  const long double n2 = _moment_ ( m_alpha , m_delta , 2 ) ;
  //
  return ( n2 * n0 - n1 * n1 ) / ( n0 * n0 )  ;
}
// ============================================================================
// sigma 
// ============================================================================
double Ostap::Math::TwoExpos::sigma () const { return std::sqrt ( variance() ) ; }
// ============================================================================
// get the derivative at given value 
// ============================================================================
double Ostap::Math::TwoExpos::derivative  ( const double x    ) const 
{ return x < m_x0 ? 0 : derivative ( x , 1 ) ; }
// ============================================================================
// get the second derivative at given value
// ============================================================================
double Ostap::Math::TwoExpos::derivative2 ( const double x    ) const
{ return x < m_x0 ? 0 : derivative ( x , 2 ) ; }
// ============================================================================
// get the Nth derivative at given value
// ============================================================================
double Ostap::Math::TwoExpos::derivative
( const double   x , 
  const unsigned N ) const 
{
  if      ( x <  m_x0 ) { return            0 ; }
  //
  const long double n0 = _moment_ ( m_alpha , m_delta , 0 ) ;
  const long double dx = x - m_x0 ;
  //
  const long double a  = tau1 () ;
  const long double b  = tau2 () ;
  //
  return 
    ( Ostap::Math::POW ( a , N ) *  std::exp ( a * dx ) - 
      Ostap::Math::POW ( b , N ) *  std::exp ( b * dx ) ) / n0 ;
}
// ============================================================================
// get the value \f$ x_{max}\$ such that  \f$ x_{max} \ge p(x) \f$ 
// ============================================================================
double Ostap::Math::TwoExpos::max_value () const
{ return evaluate ( mode() ) ; }
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::TwoExpos::tag () const 
{
  static const std::string s_name = "TwoExpos" ;
  return Ostap::Utils::hash_combiner ( s_name , m_alpha , m_delta , m_x0 ) ; 
}
// ============================================================================

// ============================================================================
Ostap::Math::TwoExpoPositive::TwoExpoPositive  
( const unsigned short N     , 
  const double         alpha , 
  const double         delta , 
  const double         x0    ,
  const double         xmin  , 
  const double         xmax  ) 
  : Ostap::Math::PolyFactor1D ( N , xmin , xmax    )
  , m_2exp     ( alpha , delta , x0 )
{}
// ============================================================================
Ostap::Math::TwoExpoPositive::TwoExpoPositive  
( const std::vector<double>& pars  ,
  const double               alpha , 
  const double               delta , 
  const double               x0    ,
  const double               xmin  , 
  const double               xmax  ) 
  : Ostap::Math::PolyFactor1D ( pars  , xmin  , xmax )
  , m_2exp     ( alpha , delta , x0   )
{}
// ============================================================================
Ostap::Math::TwoExpoPositive::TwoExpoPositive  
( const Ostap::Math::Positive& poly  ,  
  const double                 alpha , 
  const double                 delta , 
  const double                 x0    ) 
  : Ostap::Math::PolyFactor1D ( poly                 )
  , m_2exp     ( alpha , delta , x0   )
{}
// ============================================================================
Ostap::Math::TwoExpoPositive::TwoExpoPositive  
( const Ostap::Math::Positive& poly   , 
  const Ostap::Math::TwoExpos& expos  ) 
  : Ostap::Math::PolyFactor1D ( poly                 )
  , m_2exp     ( expos )
{}
// ============================================================================
Ostap::Math::TwoExpoPositive::TwoExpoPositive  
( const Ostap::Math::TwoExpos& expos  , 
  const Ostap::Math::Positive& poly   )
  : Ostap::Math::PolyFactor1D ( poly                 )
  , m_2exp     ( expos )
{}
// ============================================================================
// get the value 
// ============================================================================
double Ostap::Math::TwoExpoPositive::operator() ( const double x ) const 
{
  return 
    x < x0   () ? 0 :  
    x < xmin () ? 0 : 
    x > xmax () ? 0 : m_positive ( x ) * m_2exp ( x ) ;
}
// ============================================================================
// get the integral between xmin and xmax
// ============================================================================ 
double Ostap::Math::TwoExpoPositive::integral () const
{
  const double xlow = std::max ( x0() , xmin () ) ;
  return xlow < xmax() ? integral ( xlow , xmax () ) : 0 ;
}
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::TwoExpoPositive::integral 
( const double low  , 
  const double high ) const 
{
  //
  if      ( s_equal ( low, high ) ) { return 0 ; }
  else if ( low > high            ) { return -integral ( high , low ) ; }
  //
  const long double r1 = 
    Ostap::Math::integrate ( m_positive.bernstein() , tau1 () , low , high ) ;
  const long double r2 = 
    Ostap::Math::integrate ( m_positive.bernstein() , tau2 () , low , high ) ;
  //
  return ( r1 - r2 ) / _moment_ ( alpha() , delta () , 0 ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::TwoExpoPositive::tag () const 
{ 
  static const std::string s_name = "TwoExposPositive" ;
  return Ostap::Utils::hash_combiner ( s_name , m_positive.tag () , m_2exp.tag () ) ; 
}
// ============================================================================
// get the value \f$ x_{min}\$ such that  \f$ x_{min} \le p(x) \f$ 
// ============================================================================
double Ostap::Math::TwoExpoPositive::min_value () const 
{
  if ( xmin() < x0 ()   ) { return 0 ; } 
  const double p1 = m_positive.min_value () ;
  const double p2 = std::min ( m_2exp ( xmin () ) , m_2exp ( xmax() ) ) ;
  return p1 * p2 ;
}
// ============================================================================
// get the value \f$ x_{max}\$ such that  \f$ x_{max} \ge p(x) \f$ 
// ============================================================================
double Ostap::Math::TwoExpoPositive::max_value () const
{
  const double p1 = m_positive.max_value () ;
  //
  const double emode = m_2exp.mode() ;
  const double p2    =
    xmin () <= emode && emode <= xmax () ? m_2exp.mode() :
    std::max ( m_2exp ( xmin () ) , m_2exp ( xmax () ) ) ;
  //
  return p1 * p2 ;
}


// =============================================================================
// constructor for all elements 
// =============================================================================
Ostap::Math::Argus::Argus 
( const double mu  , 
  const double  c  , 
  const double chi ) 
  : m_mu   ( mu ) 
  , m_c    ( std::abs ( c   ) )
  , m_chi  ( std::abs ( chi ) )  
  , m_norm ( -1 ) 
{ setChi ( chi ) ; }
// =============================================================================
// set mu parameter
// =============================================================================
bool Ostap::Math::Argus::setMu  ( const double value ) 
{
  if ( s_equal ( m_mu , value ) ) { return false ; }
  m_mu = value ;
  return true ;    
}
// =============================================================================
// set c parameter
// =============================================================================
bool Ostap::Math::Argus::setC   ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_c , avalue ) ) { return false ; }
  m_c = avalue ;
  return true ;    
}
// =============================================================================
// set chi parameter
// =============================================================================
bool Ostap::Math::Argus::setChi ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_chi , avalue ) && 0 < m_norm ) { return false ; }
  m_chi  = avalue ;
  m_norm = std::pow ( m_chi , 3 ) / psi ( m_chi ) * s_sqrt_1_2pi ;
  return true ;    
}
// ============================================================================
/*  helper function 
 *  \f$ \Psi ( \chi ) = \Phi(\chi )  - \chi \phi  (\chi ) - \frac{1}{2} \f$ 
 */
// ============================================================================
double Ostap::Math::Argus::psi ( const double value ) const 
{ return Ostap::Math::gauss_cdf ( value ) - value * Ostap::Math::gauss_pdf ( value ) - 0.5 ; }
// ============================================================================
// evaluate the function 
// ============================================================================
double Ostap::Math::Argus::evaluate   ( const double x ) const 
{
  if ( x + m_c <= m_mu || m_mu <= x ) { return 0 ; }
  const double dx = ( x + m_c - m_mu ) / m_c ; 
  const double dd = 1 - dx * dx ;
  return m_norm * dx * std::sqrt ( dd ) * std::exp ( -0.5 * m_chi * m_chi * dd ) / m_c ;  
}
// ============================================================================
// get CDF 
// ============================================================================
double Ostap::Math::Argus::cdf        ( const double x ) const 
{
  //
  if       ( x + m_c <= m_mu ) { return 0 ; }
  else  if ( m_mu <= x       ) { return 1 ; }
  //
  const double dx = ( x + m_c - m_mu ) / m_c ;
  const double dd = std::sqrt ( 1 - dx * dx ) ;
  //
  return 1 - psi ( m_chi * dd ) / psi ( m_chi ) ;
}
// ============================================================================
// get the integral 
// ============================================================================
double Ostap::Math::Argus::integral   () const { return 1 ; }
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::Argus::integral  
( const double low  ,
  const double high ) const 
{
  if      ( s_equal ( low , high )             ) { return 0 ; }
  else if ( high < low                         ) { return -integral ( high , low ) ; }
  else if ( high + m_c <=  m_mu                ) { return 0 ; }
  else if ( m_mu       <=  low                 ) { return 0 ; }
  else if ( low  + m_c <= m_mu && m_mu <= high ) { return 1 ; }    
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ===========================================================================
// mean of the distribution 
// ===========================================================================
double Ostap::Math::Argus::mean     () const 
{
  const double c2  = 0.25 * m_chi * m_chi ;
  return ( m_mu - m_c ) + 
    0.5 * m_c * m_chi * s_sqrt_pi_2 * std::exp ( - c2 ) * Ostap::Math::bessel_In ( 1 , c2 ) / psi ( m_chi );
}
// ===========================================================================
// mode of the distribution 
// ===========================================================================
double Ostap::Math::Argus::mode     () const 
{
  const double c2  = m_chi * m_chi ;
  return ( m_mu - m_c ) + 
    m_c * s_1_sqrt2 * std::sqrt ( ( c2 - 2 ) + std::sqrt ( c2 * c2 + 4 ) ) / m_chi ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Argus::tag () const 
{ 
  static const std::string s_name = "Argus" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_c , m_chi ) ;
}
// =============================================================================
// constructor for all elements 
// =============================================================================
Ostap::Math::GenArgus::GenArgus 
( const double mu  , 
  const double  c  , 
  const double chi ,
  const double dp  ) 
  : m_mu   ( mu ) 
  , m_c    ( std::abs ( c   ) )
  , m_chi  ( std::abs ( chi ) )  
  , m_dp   ( std::abs ( dp  ) )  
  , m_norm ( -1 ) 
{ 
  setChi ( chi ) ;
  setDp  ( dp  ) ;
}
// =============================================================================
// set mu parameter
// =============================================================================
bool Ostap::Math::GenArgus::setMu  ( const double value ) 
{
  if ( s_equal ( m_mu , value ) ) { return false ; }
  m_mu = value ;
  return true ;    
}
// =============================================================================
// set c parameter
// =============================================================================
bool Ostap::Math::GenArgus::setC   ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_c , avalue ) ) { return false ; }
  m_c = avalue ;
  return true ;    
}
// =============================================================================
// set chi parameter
// =============================================================================
bool Ostap::Math::GenArgus::setChi ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_chi , avalue ) && 0 < m_norm ) { return false ; }
  m_chi  = avalue ;
  //
  const double c2 = m_chi * m_chi ;
  const double p1 = p() + 1 ;
  m_norm = 2 * std::pow ( 0.5 * c2 , p1 ) /
  ( std::tgamma ( p() + 1 ) * ( 1 - Ostap::Math::gamma_inc_Q (   p1 , 0.5 * c2 ) ) ) ;
  //
  return true ;    
}
// =============================================================================
// set dp parameter
// =============================================================================
bool Ostap::Math::GenArgus::setDp ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_dp , avalue ) && 0 < m_norm ) { return false ; }
  m_dp   = avalue ;
  //
  const double c2 = m_chi * m_chi ;
  const double p1 = p() + 1 ;
  m_norm = 2 * std::pow ( 0.5 * c2 , p1 ) /
  ( std::tgamma ( p() + 1 ) * ( 1 - Ostap::Math::gamma_inc_Q (   p1 , 0.5 * c2 ) ) ) ;
  //
  return true ;    
}
// ============================================================================
// evaluate the function 
// ============================================================================
double Ostap::Math::GenArgus::evaluate   ( const double x ) const 
{
  if ( x + m_c <= m_mu || m_mu <= x ) { return 0 ; }
  const double dx = ( x + m_c - m_mu ) / m_c ; 
  const double dd = 1 - dx * dx ;
  return m_norm * dx * std::pow ( dd , p () ) * std::exp ( -0.5 * m_chi * m_chi * dd ) / m_c ;  
}
// ============================================================================
// get CDF 
// ============================================================================
double Ostap::Math::GenArgus::cdf        ( const double x ) const 
{
  //
  if       ( x + m_c <= m_mu ) { return 0 ; }
  else  if ( m_mu <= x       ) { return 1 ; }
  //
  const double dx = ( x + m_c - m_mu ) / m_c ;
  const double dd = 1 - dx * dx  ;
  //
  const double p1 = p() + 1 ;
  const double c2 = 0.5 * m_chi * m_chi ;
  //
  const double a1 = Ostap::Math::gamma_inc_Q ( p1 , c2 * dd ) ;
  const double a2 = Ostap::Math::gamma_inc_Q ( p1 , c2      ) ;
  //
  return ( a1 - a2 ) / ( 1 - a2 ) ;
}
// ============================================================================
// get the integral 
// ============================================================================
double Ostap::Math::GenArgus::integral   () const { return 1 ; }
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::GenArgus::integral  
( const double low  ,
  const double high ) const 
{
  if      ( s_equal ( low , high )             ) { return 0 ; }
  else if ( high < low                         ) { return -integral ( high , low ) ; }
  else if ( high + m_c <=  m_mu                ) { return 0 ; }
  else if ( m_mu       <=  low                 ) { return 0 ; }
  else if ( low  + m_c <= m_mu && m_mu <= high ) { return 1 ; }    
  //
  return cdf ( high ) - cdf ( low ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::GenArgus::tag () const 
{ 
  static const std::string s_name = "GenArgus" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mu , m_c , m_chi , m_dp ) ;
}

// ============================================================================
/*  constructor from all parameters
 *  @param a position of the left paraboilic horn
 *  @param delta half-distance fron left to right parabolic horn 
 *  @param phi   linear correction parameter ("efficiency")
 */
// ============================================================================
Ostap::Math::HORNSdini::HORNSdini
( const double a     , 
  const double delta ,
  const double phi   ) 
  : m_a        ( a   ) 
  , m_delta    ( std::abs ( delta ) )
  , m_phi      ( phi ) 
  , m_cos2_phi ( std::pow ( std::cos ( phi + s_pi_4 ) , 2 ) ) 
  , m_sin2_phi ( std::pow ( std::sin ( phi + s_pi_4 ) , 2 ) ) 
{}
// ======================================================================
// evaluate the function 
// ======================================================================
double Ostap::Math::HORNSdini::evaluate ( const double x )  const 
{
  if ( x < xmin () || xmax () < x ) { return 0 ; }
  const double z = ( x - m_a ) / m_delta - 1 ;
  return 1.5 * z * z * ( 1 + z * ( m_cos2_phi - m_sin2_phi ) ) / m_delta ;                        
}
// =============================================================================
// set a parameter
// =============================================================================
bool Ostap::Math::HORNSdini::setA   ( const double value ) 
{
  if ( s_equal ( m_a , value ) ) { return false ; }
  m_a = value ;
  return true ;    
}
// =============================================================================
// set delta parameter
// =============================================================================
bool Ostap::Math::HORNSdini::setDelta ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_delta , avalue ) ) { return false ; }
  m_delta = avalue ;
  return true ;    
}
// =============================================================================
// set phi parameter
// =============================================================================
bool Ostap::Math::HORNSdini::setPhi ( const double value ) 
{
  if ( s_equal ( m_phi , value ) ) { return false ; }
  //
  m_phi      = value ;
  m_cos2_phi = std::pow ( std::cos ( m_phi + s_pi_4 ) , 2 ) ;
  m_sin2_phi = std::pow ( std::sin ( m_phi + s_pi_4 ) , 2 ) ;
  //
  return true ;    
}
// ============================================================================
// get the integral 
// ============================================================================
double Ostap::Math::HORNSdini::integral  () const { return 1 ; }
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::HORNSdini::integral 
( const double low  ,
  const double high ) const 
{
  //
  if      ( s_equal ( low , high )           ) { return 0 ; }
  else if ( high < low                       ) { return - integral ( high , low ) ; }
  else if ( high < xmin ()                   ) { return 0 ; }
  else if ( low  > xmax ()                   ) { return 0 ; }
  else if ( low <= xmin () && high>= xmax () ) { return 1 ; }
  //
  const double xl  = std::max ( low  , xmin () ) ;
  const double xh  = std::min ( high , xmax () ) ;
  //
  const double zl  = ( xl - m_a ) / m_delta - 1  ;
  const double zh  = ( xh - m_a ) / m_delta - 1 ;
  //
  const double zh3 = std::pow ( zh , 3 ) ;
  const double zl3 = std::pow ( zl , 3 ) ;
  //
  return ( ( zh3 - zl3 ) / 3 + 
           ( m_cos2_phi - m_sin2_phi ) * ( zh3 * zh - zl3 * zl ) / 4 ) * 1.5 ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::HORNSdini::tag () const 
{ 
  static const std::string s_name = "HORNSdini" ;
  return Ostap::Utils::hash_combiner ( s_name , m_a , m_delta , m_phi ) ;
}

// ============================================================================
/*  constructor from all parameters
 *  @param a position of the left paraboilic horn
 *  @param delta half-distance fron left to right parabolic horn 
 *  @param phi   linear correction parameter ("efficiency")
 */
// ============================================================================
Ostap::Math::HILLdini::HILLdini
( const double a     , 
  const double delta ,
  const double phi   ) 
  : m_a        ( a   ) 
  , m_delta    ( std::abs ( delta ) )
  , m_phi      ( phi ) 
  , m_cos2_phi ( std::pow ( std::cos ( phi + s_pi_4 ) , 2 ) ) 
  , m_sin2_phi ( std::pow ( std::sin ( phi + s_pi_4 ) , 2 ) ) 
{}
// ======================================================================
// evaluate the function 
// ======================================================================
double Ostap::Math::HILLdini::evaluate ( const double x )  const 
{
  if ( x < xmin () || xmax () < x ) { return 0 ; }
  const double z  = ( x - m_a ) / m_delta - 1 ;
  return 0.75 * ( 1 - z * z ) * ( 1 + z * ( m_cos2_phi - m_sin2_phi ) ) / m_delta ;                        
}
// =============================================================================
// set a parameter
// =============================================================================
bool Ostap::Math::HILLdini::setA   ( const double value ) 
{
  if ( s_equal ( m_a , value ) ) { return false ; }
  m_a = value ;
  return true ;    
}
// =============================================================================
// set delta parameter
// =============================================================================
bool Ostap::Math::HILLdini::setDelta ( const double value ) 
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_delta , avalue ) ) { return false ; }
  m_delta = avalue ;
  return true ;    
}
// =============================================================================
// set phi parameter
// =============================================================================
bool Ostap::Math::HILLdini::setPhi ( const double value ) 
{
  if ( s_equal ( m_phi , value ) ) { return false ; }
  //
  m_phi      = value ;
  m_cos2_phi = std::pow ( std::cos ( m_phi + s_pi_4 ) , 2 ) ;
  m_sin2_phi = std::pow ( std::sin ( m_phi + s_pi_4 ) , 2 ) ;
  //
  return true ;    
}
// ============================================================================
// get the integral 
// ============================================================================
double Ostap::Math::HILLdini::integral  () const { return 1 ; }
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::HILLdini::integral 
( const double low  ,
  const double high ) const 
{
  //
  if      ( s_equal ( low , high )           ) { return 0 ; }
  else if ( high < low                       ) { return - integral ( high , low ) ; }
  else if ( high < xmin ()                   ) { return 0 ; }
  else if ( low  > xmax ()                   ) { return 0 ; }
  else if ( low <= xmin () && high>= xmax () ) { return 1 ; }
  //
  const double xl  = std::max ( low  , xmin () ) ;
  const double xh  = std::min ( high , xmax () ) ;
  //
  const double zh  = ( xh - m_a ) / m_delta - 1 ; 
  const double zl  = ( xl - m_a ) / m_delta - 1 ;
  //
  const double zh3 = std::pow ( zh , 3 ) ;
  const double zl3 = std::pow ( zl , 3 ) ;
  //
  const double A   = m_cos2_phi - m_sin2_phi ;
  //
  return ( (   zh       - zl       )     
           + ( zh  * zh - zl  * zl ) * A / 2 
           - ( zh3      - zl3      )     / 3
           - ( zh3 * zh - zl3 * zl ) * A / 4 ) * 0.75 ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::HILLdini::tag () const 
{ 
  static const std::string s_name = "HILLdini" ;
  return Ostap::Utils::hash_combiner ( s_name , m_a , m_delta , m_phi ) ;
}

// ============================================================================
/*  Constructor from all parameters
 *  @param right dump direction
 *  @param x0    threshold value 
 *  @param sigma sigma  
 */
// ============================================================================
Ostap::Math::CutOffGauss::CutOffGauss
( const bool   right , 
  const double x0    , 
  const double sigma ) 
  : m_right ( right ) 
  , m_x0    ( x0                 )
  , m_sigma ( std::abs ( sigma ) ) 
{}
// =========================================================================
// update sigma
// ============================================================================
bool Ostap::Math::CutOffGauss::setSigma ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_sigma , avalue ) ) { return false ; }
  m_sigma = avalue ;
  return true ;
}
// ============================================================================
// update x0
// ============================================================================
bool Ostap::Math::CutOffGauss::setX0 ( const double value )
{
  if ( s_equal ( m_x0 , value ) ) { return false ; }
  m_x0 = value ;
  return true ;
}
// ============================================================================
// the main method 
// ============================================================================
double Ostap::Math::CutOffGauss::operator() ( const double x ) const 
{
  // 
  if      (  m_right && x <= m_x0 ) { return 1 ; }
  else if ( !m_right && x >= m_x0 ) { return 1 ; }
  //
  const double dx = ( x - m_x0 ) / m_sigma ;
  return std::exp ( -0.5 * dx * dx ) ;
}
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::CutOffGauss::integral
( const double low  ,
  const double high ) const 
{
  if      ( s_equal ( low , high ) ) { return 0                        ; }
  else if ( low > high             ) { return -integral ( high , low ) ; }
  //
  if ( low < m_x0 && m_x0 < high ) 
  { return integral ( low , m_x0 ) + integral ( m_x0 , high ) ; }
  //
  if      (  m_right && high <= m_x0 ) { return high - low ; }
  else if ( !m_right && low  >= m_x0 ) { return high - low ; }
  //
  return s_sqrt_2pi * m_sigma * 
    ( Ostap::Math::gauss_cdf ( high , m_x0 , m_sigma ) -
      Ostap::Math::gauss_cdf ( low  , m_x0 , m_sigma ) ) ;  
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::CutOffGauss::tag () const 
{  
  static const std::string s_name = "CutOffGauss" ;
  return Ostap::Utils::hash_combiner ( s_name , m_right , m_x0 , m_sigma ) ; 
}
// ============================================================================

// ============================================================================
/*  Constructor from all parameters
 *  @param right dump direction
 *  @param x0    threshold value 
 *  @param nu    parameter nu 
 *  @param sigma parameter sigma  
 */
// ============================================================================
Ostap::Math::CutOffStudent::CutOffStudent 
( const bool   right , 
  const double x0    , 
  const double nu    ,
  const double sigma ) 
  : m_right  ( right ) 
  , m_x0     ( x0    ) 
  , m_nu     ( -1    ) 
  , m_sigma  ( std::abs ( sigma ) ) 
  , m_C      ( -1    )
{
  setNu ( nu ) ;
}
// =========================================================================
// update sigma
// ============================================================================
bool Ostap::Math::CutOffStudent::setSigma ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_sigma , avalue ) ) { return false ; }
  m_sigma = avalue ;
  return true ;
}
// =========================================================================
// update nu
// ============================================================================
bool Ostap::Math::CutOffStudent::setNu ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_nu , avalue ) ) { return false ; }
  m_nu = avalue ;
  //
  m_C  = std::exp ( - std::lgamma (  0.5 * ( m_nu + 1 ) )  
                    + std::lgamma (  0.5 * ( m_nu     ) ) 
                    + 0.5 * std::log ( m_nu * s_pi ) ) ;
  return true ;
}
// ============================================================================
// update x0
// ============================================================================
bool Ostap::Math::CutOffStudent::setX0 ( const double value )
{
  if ( s_equal ( m_x0 , value ) ) { return false ; }
  m_x0 = value ;
  return true ;
}
// ============================================================================
// the main method 
// ============================================================================
double Ostap::Math::CutOffStudent::operator() ( const double x ) const 
{
  // 
  if      (  m_right && x <= m_x0 ) { return 1 ; }
  else if ( !m_right && x >= m_x0 ) { return 1 ; }
  //
  const double dx = ( x - m_x0 ) / m_sigma ;
  return std::pow ( 1 + dx * dx / m_nu , -0.5 * ( m_nu + 1 ) ) ;
}
// ============================================================================
// get the integral between low and high
// ============================================================================
double Ostap::Math::CutOffStudent::integral
( const double low  ,
  const double high ) const 
{
  if      ( s_equal ( low , high ) ) { return 0                        ; }
  else if ( low > high             ) { return -integral ( high , low ) ; }
  //
  if ( low < m_x0 && m_x0 < high ) 
  { return integral ( low , m_x0 ) + integral ( m_x0 , high ) ; }
  //
  if      (  m_right && high <= m_x0 ) { return high - low ; }
  else if ( !m_right && low  >= m_x0 ) { return high - low ; }
  //
  const double xl = ( low  - m_x0 ) / m_sigma ;
  const double xh = ( high - m_x0 ) / m_sigma ;
  //
  return m_C * m_sigma * ( Ostap::Math::student_cdf ( xh , m_nu ) - 
                           Ostap::Math::student_cdf ( xl , m_nu ) )  ;
}

// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::CutOffStudent::tag () const 
{ 
  static const std::string s_name = "CutOffStudent" ;
  return Ostap::Utils::hash_combiner ( s_name , m_right , m_x0 , m_nu , m_sigma ) ;
}
// ============================================================================

// ============================================================================
//                                                                      The END 
// ============================================================================
