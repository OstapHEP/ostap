
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
#include "Ostap/PhaseSpacePol.h"
// ============================================================================
//  Local 
// ============================================================================
#include "local_math.h"
#include "local_hash.h"
#include "status_codes.h"
#include "Integrator1D.h"
// ============================================================================
/** @file
 *  Implementation file for functions from the file Ostap/PhaseSpacePol.h
 *
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2010-04-19
 */
// ============================================================================
/*  constructor from thresholds and number of particles
 *  @param threshold_L the low-mass  threshold
 *  @param threshold_H the high-mass threshold
 *  @param l           how many particles we consider
 *  @param n           total number of particles ( n>l!)
 *  @param N           degree of polynomial 
 */
// ======================================================================
Ostap::Math::PhaseSpacePol::PhaseSpacePol 
( const double         threshold1  ,
  const double         threshold2  ,
  const unsigned short l           ,
  const unsigned short n           , 
  const unsigned short N           )  // degree of polynomial
  : Ostap::Math::PolyFactor1D ( N , 
				std::min ( std::abs ( threshold1 ) , std::abs ( threshold2 ) ) , 
				std::max ( std::abs ( threshold1 ) , std::abs ( threshold2 ) ) )     
  , m_phasespace ( threshold1 , threshold2 , l , n ) 
  , m_workspace  ()
{}
// =====================================================================
/*  constructor from the phase space and polynomial degree 
 *  @param ps          phase space factor 
 *  @param N           degree of polynomial 
 */
// =====================================================================
Ostap::Math::PhaseSpacePol::PhaseSpacePol 
( const Ostap::Math::PhaseSpaceNL& ps ,
  const unsigned short             N  )  // degree of polynomial
  : Ostap::Math::PolyFactor1D ( N  , ps.lowEdge() , ps.highEdge() )     
  , m_phasespace ( ps ) 
  , m_workspace  ()
{}
// ======================================================================
/*  constructor from phase space and polynomial degree 
 *  @param ps          phase space factor 
 *  @param N           degree of polynomial 
 */
// =====================================================================
Ostap::Math::PhaseSpacePol::PhaseSpacePol 
( const Ostap::Math::PhaseSpaceNL& ps    ,
  const unsigned short             N     , 
  const double                     xlow  , 
  const double                     xhigh ) 
  : Ostap::Math::PolyFactor1D ( N  , 
				std::max ( ps. lowEdge() , std::min ( xlow , xhigh ) ) ,
				std::min ( ps.highEdge() , std::max ( xlow , xhigh ) ) )
  , m_phasespace ( ps ) 
  , m_workspace  ()
{}
// ======================================================================
// constructor from phase space and polynomial
// ======================================================================
Ostap::Math::PhaseSpacePol::PhaseSpacePol 
( const PhaseSpaceNL&          ps  ,
  const Ostap::Math::Positive& pol ) 
  : Ostap::Math::PolyFactor1D ( pol ) 
  , m_phasespace ( ps  ) 
  , m_workspace  () 
{
  Ostap::Assert ( m_phasespace.lowEdge () < m_positive.xmax      () , 
                  "Invalid setting of lowEdge/highEdge/xmin/xmax"   ,
                  "Ostap::Math::PhaseSpacePol"                      ,
		  INVALID_PARAMETERS , __FILE__ , __LINE__          ) ;
  Ostap::Assert ( m_positive.xmin      () < m_phasespace.highEdge() , 
                  "Invalid setting of lowEdge/highEdge/xmin/xmax"   ,
                  "Ostap::Math::PhaseSpacePol"                      , 
		  INVALID_PARAMETERS , __FILE__ , __LINE__          ) ;
}
// =====================================================================
// evaluate N/L-body modulated phase space
// =====================================================================
double Ostap::Math::PhaseSpacePol::evaluate ( const double x ) const 
{
  //
  if      ( x < m_phasespace . lowEdge () ) { return 0 ; }
  else if ( x > m_phasespace .highEdge () ) { return 0 ; }
  else if ( x < m_positive   .   xmin  () ) { return 0 ; }
  else if ( x > m_positive   .   xmax  () ) { return 0 ; }
  //
  return m_positive ( x ) * m_phasespace ( x ) ;
}
// =====================================================================
// get the integral
// =====================================================================
double Ostap::Math::PhaseSpacePol::integral () const 
{
  //
  if      ( m_phasespace.highEdge() <= m_positive.xmin() ) { return 0 ; }
  else if ( m_phasespace. lowEdge() >= m_positive.xmax() ) { return 0 ; }
  //
  const double mn = std::max ( m_phasespace. lowEdge() ,  m_positive.xmin () ) ;
  const double mx = std::min ( m_phasespace.highEdge() ,  m_positive.xmax () ) ;
  //
  return integral ( mn , mx ) ;
}
// =====================================================================
// get the integral between low and high limits
// =====================================================================
double  Ostap::Math::PhaseSpacePol::integral 
( const double low  ,
  const double high ) const 
{
  //
  if      ( s_equal ( low , high ) ) { return  0 ; }
  else if (           low > high   ) { return -1 * integral ( high , low ) ; }
  //
  if      ( high <= m_phasespace .  lowEdge () ) { return 0 ; }
  else if ( high <= m_positive   .     xmin () ) { return 0 ; }
  else if ( low  >= m_phasespace . highEdge () ) { return 0 ; }
  else if ( low  >= m_positive   .     xmax () ) { return 0 ; }
  //
  const double mn    = std::max ( m_phasespace. lowEdge() ,  m_positive.xmin () ) ;
  const double mx    = std::min ( m_phasespace.highEdge() ,  m_positive.xmax () ) ;
  //
  const double xlow  = std::max ( low  , mn ) ;
  const double xhigh = std::min ( high , mx ) ;
  //
  // use GSL to evaluate the integral
  //
  static const Ostap::Math::GSL::Integrator1D<PhaseSpacePol> s_integrator {} ;
  static char s_message[] = "Integral(PhaseSpacePol)" ;
  //
  const auto F = s_integrator.make_function ( this ) ;
  int    ierror   = 0   ;
  double result   = 1.0 ;
  double error    = 1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.qag_integrate
    ( tag () ,  
      &F     , 
      xlow   , xhigh       ,          // low & high edges
      workspace ( m_workspace ) ,     // workspace
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
std::size_t Ostap::Math::PhaseSpacePol::tag () const 
{
  static const std::string s_name = "PhaseSpacePol" ;
  return Ostap::Utils::hash_combiner ( s_name , m_phasespace.tag () , m_positive.tag () ) ; 
}
// ============================================================================


