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
#include "Ostap/Workspace.h"
#include "Ostap/MoreMath.h"
#include "Ostap/Spectra.h"
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
// Tsallis function 
// ============================================================================
/*  constructor from all parameters 
 *  @param mass particle mass  (M>0)
 *  @param n    n-parameter    (N>1)  
 *  @param T    T-parameter    (T>0)
 */
// ============================================================================
Ostap::Math::Tsallis::Tsallis 
( const double mass  , 
  const double n     ,  
  const double T     ) 
  : m_mass ( std::abs ( mass ) ) 
  , m_n    ( std::abs ( n    ) )
  , m_T    ( std::abs ( T    ) )
  , m_workspace() 
{
  Ostap::Assert ( 0 < m_mass ,
		  "Mass  must be positive!"    , 
		  "Ostap::Math::Tsallis"       ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  Ostap::Assert ( 0 < m_T ,
		  "Temperature must be positive!" , 
		  "Ostap::Math::Tsallis"          ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
}
// ============================================================================
// set new value for mass
// ============================================================================
bool Ostap::Math::Tsallis::setMass ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_mass , avalue ) ) { return false ; }
  //
  Ostap::Assert ( avalue  ,
		  "Mass  must be positive!"    , 
		  "Ostap::Math::Tsallis"       ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  //
  m_mass = avalue ;  
  return true ;
}
// ============================================================================
// set new value for n-parameter
// ============================================================================
bool Ostap::Math::Tsallis::setN ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_n , avalue ) ) { return false ; }
  
  m_n = avalue ;  
  return true ;
}
// ============================================================================
// set new value for T-parameter
// ============================================================================
bool Ostap::Math::Tsallis::setT ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_T , avalue ) ) { return false ; }
  Ostap::Assert ( 0 < avalue  ,
		  "Temperature must be positive!" , 
		  "Ostap::Math::Tsallis"          ,
		  INVALID_PARAMETER , __FILE__ , __LINE__ ) ;
  m_T = avalue ;
  return true ;
}
// ============================================================================
//  get Tsallis PDF  
// ============================================================================
double Ostap::Math::Tsallis::evaluate ( const double x ) const 
{ return x <= 0 ? 0.0 : x * std::pow ( 1.0 + eTkin ( x ) / ( m_T * m_n ) , -m_n ) ; }
// ============================================================================
//  get Tsallis integrals  
// ============================================================================
double Ostap::Math::Tsallis::integral 
( const double low  , 
  const double high ) const
{
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high < low             ) { return - integral ( high , low ) ; }
  else if ( high <= xmin ()        ) { return 0 ; }
  //
  const double _low = std::max ( low , xmin () ) ;
  //
  // split too large intervals
  if ( 0 < m_mass ) 
  {
    // split points 
    static const std::array<int,7> s_split = {{ 1 ,  3  , 10 , 20 , 50 , 100 , 1000 }} ;
    for( const auto p : s_split )
    {
      const double middle = m_mass * p ;
      if (  _low < middle && middle < high ) 
      { return integral ( _low , middle ) + integral ( middle , high ) ; }
    }
  }
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Tsallis> s_integrator {} ;
  static char s_message[] = "Integral(Tsallis)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag  () , 
      &F      , 
      low     , high  ,              // low & high edges
      workspace ( m_workspace ) ,    // workspace
      s_APRECISION         ,          // absolute precision
      s_RPRECISION         ,          // relative precision
      m_workspace.size()   ,          // size of workspace
      s_message            , 
      __FILE__ , __LINE__  ) ;
  //
  return result ;
  //
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Tsallis::tag () const 
{ 
  static const std::string s_name = "Tsallis" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mass , m_n , m_T ) ; 
}
// ============================================================================

// ============================================================================
// QGSM function 
// ============================================================================
/*  constructor from all parameters 
 *  @param mass particle mass  (M>0)
 *  @param n    n-parameter    (N>1)  
 *  @param T    T-parameter    (T>0)
 */
// ============================================================================
Ostap::Math::QGSM::QGSM 
( const double mass , 
  const double b    ) 
  : m_mass ( std::abs ( mass ) ) 
  , m_b    ( std::abs ( b    ) )
  , m_workspace() 
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Math::QGSM::~QGSM(){}
// ============================================================================
// set new value for mass
// ============================================================================
bool Ostap::Math::QGSM::setMass ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_mass , avalue ) ) { return false ; }
  m_mass = avalue ;
  return true ;
}
// ============================================================================
// set new value for b-parameter
// ============================================================================
bool Ostap::Math::QGSM::setB ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_b , avalue ) ) { return false ; }
  m_b = avalue ;
  return true ;
}
// ============================================================================
//  get QGSM PDF  
// ============================================================================
double Ostap::Math::QGSM::pdf ( const double x ) const 
{ return x <= 0 ? 0.0 : x * std::exp ( -m_b * eTkin ( x ) ) ; }
// ============================================================================
//  get QGSM integrals  
// ============================================================================
double Ostap::Math::QGSM::integral 
( const double low  , 
  const double high ) const
{
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high < low             ) { return - integral ( high , low ) ; }
  else if ( high <= xmin()         ) { return 0 ; }
  //
  const double _low = std::max ( low , xmin() ) ;
  //
  // split too large intervals
  if ( 0 < m_mass ) 
  {
    // split points 
    static const std::array<int,5> s_split = {{ 1 ,  3  , 10 , 20 , 50 }} ;
    for( const auto p : s_split )
    {
      const double middle = m_mass * p ;
      if (  _low < middle && middle < high ) 
      { return integral ( _low , middle ) + integral ( middle , high ) ; }
    }
  }
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<QGSM> s_integrator {} ;
  static char s_message [] = "Integral(QGSM)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag  () , 
      &F      , 
      low     , high      ,          // low & high edges
      workspace ( m_workspace ) ,    // workspace
      s_APRECISION         ,          // absolute precision
      s_RPRECISION         ,          // relative precision
      m_workspace.size()              ,          // size of workspace
      s_message           , 
      __FILE__ , __LINE__ ) ;
  //
  return result ;
  //
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::QGSM::tag () const 
{ 
  static const std::string s_name = "QGSM" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mass , m_b ) ; 
}
// ============================================================================



// ============================================================================
/*  Constructor for all parameters
 *  @param mass   mas sof th eparticle 
 *  @param beta   inverse temperature
 */
// ============================================================================
Ostap::Math::Hagedorn::Hagedorn 
( const double mass , 
  const double beta ) 
  : m_mass ( std::abs ( mass ) ) 
  , m_beta ( std::abs ( beta ) ) 
{}
// ============================================================================
// set new value for mass
// ============================================================================
bool Ostap::Math::Hagedorn::setMass ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_mass , avalue ) ) { return false ; }
  m_mass = avalue ;
  return true ;
}
// ============================================================================
// set new value for inverse temporature 
// ============================================================================
bool Ostap::Math::Hagedorn::setBeta ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_beta , avalue ) ) { return false ; }
  m_beta = avalue ;
  return true ;
}
// ============================================================================
// the main method 
// ============================================================================
double Ostap::Math::Hagedorn::evaluate ( const double x ) const 
{
  if ( x <= 0 ) { return 0 ; }
  //
  const double mt  = mT ( x )    ;
  const double arg = m_beta * mt ;
  //
  return 
    // GSL_LOG_DBL_MAX < arg ? 0.0 : 
    300                < arg ? 0.0 : 
    x * mt * Ostap::Math::bessel_Kn ( 1 , arg ) / m_beta ;
}
// ============================================================================
//  get Hagedorn integral  
// ============================================================================
double Ostap::Math::Hagedorn::integral 
( const double low  , 
  const double high ) const
{
  if      ( s_equal ( low , high ) ) { return 0 ; }
  else if ( high < low             ) { return - integral ( high , low ) ; }
  else if ( high <= 0              ) { return 0 ; }
  //
  const double xmin = std::max ( 0.0 , low ) ;
  const double xmax = high ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<Hagedorn> s_integrator {} ;
  static char s_message [] = "Integral(HAgedorn)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag  ()                   , 
      &F                        , 
      xmin    , xmax            ,   // low & high edges
      workspace ( m_workspace ) ,   // workspace
      s_APRECISION              ,   // absolute precision
      s_RPRECISION              ,   // relative precision
      m_workspace.size()        ,   // size of workspace
      s_message                 , 
      __FILE__ , __LINE__       ) ;
  //
  return result ;
  //
}
// ============================================================================
// get the mean value 
// ============================================================================
double Ostap::Math::Hagedorn::mean () const 
{
  const double mt = m_mass * m_beta ;
  return 
    s_SQRTPIHALF * std::sqrt ( m_mass / m_beta ) * 
    Ostap::Math::bessel_Knu ( 5.0/2 , mt ) /
    Ostap::Math::bessel_Kn  ( 2     , mt ) ;
}
// ============================================================================
// get the tag
// ============================================================================
std::size_t Ostap::Math::Hagedorn::tag () const 
{ 
  static const std::string s_name = "Hagedorn" ;
  return Ostap::Utils::hash_combiner ( s_name , m_mass , m_beta ) ; 
}
// ============================================================================

 
// ============================================================================
//                                                                     The END 
// ============================================================================
