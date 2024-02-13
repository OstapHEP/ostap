// ============================================================================
// Include files 
// ============================================================================
// STD & STL
// ============================================================================
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
// ============================================================================
// local
// ============================================================================
#include "Ostap/Math.h"
#include "Ostap/Hash.h"
#include "Ostap/NSphere.h"
// ============================================================================
/** @file 
 *  Implementation of class Ostap::Math::NSphere 
 *  @see Ostap::Math::NSPhere 
 *  @author Vanya BELYAEV Ivan,Belyaev@itep.ru
 *  @date 2014-01-21 
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  /// equality criteria for doubles
  const Ostap::Math::Equal_To<double> s_equal{} ; // equality criteria for doubles
  /// zero for doubles  
  const Ostap::Math::Zero<double>     s_zero{}  ; // zero for doubles
  // ==========================================================================
  inline 
  std::pair<double,double> _sincos_ ( const double phase ) 
  {
    if ( 0 == phase || s_zero ( phase ) ) { return std::make_pair ( 0 , 1 ) ; }
    //
    // return std::make_pair ( std::sin( phase ) , std::cos( phase ) ) ;
    // double sinv = 0 ;
    // double cosv = 1 ;
    // vdt::fast_sincos ( phase , sinv , cosv ) ;
    const double sinv  = std::sin ( phase ) ;
    const double cosv  = std::cos ( phase ) ;
    //
    const double acosv = std::abs ( cosv ) ;
    if ( 0 == sinv  || 1 == acosv ) 
    { return std::make_pair ( 0 , 0 < cosv ? 1 : -1 ) ; }
    //
    const double asinv = std::abs ( sinv ) ;
    if ( 0 == cosv  || 1 == asinv )
    { return std::make_pair ( 0 < sinv ? 1 : -1 , 0 ) ; }
    //
    if ( s_zero  ( asinv ) || s_equal ( acosv , 1 ) )
    { return std::make_pair ( 0 , 0 < cosv ? 1 : -1 ) ; }
    //
    if ( s_zero  ( acosv ) || s_equal ( asinv , 1 ) )
    { return std::make_pair ( 0 < sinv ? 1 : -1 , 0 ) ; }
    //
    return std::make_pair ( sinv , cosv ) ;
  }
  // ==========================================================================
  /** @var s_PIHALF 
   *  useful constant
   */
  const double s_PIHALF = 0.5 * M_PI ;
  // ==========================================================================
}
// ============================================================================
/*  Standard constructor for rotated sphere 
 *  @param nPhases  dimensionality of N-sphere 
 */
// ============================================================================
Ostap::Math::NSphere::NSphere 
( const unsigned short N ) 
  : Ostap::Math::NSphere::NSphere ( N , N + 1 ) 
{}
// ============================================================================
/*  Standard constructor
 *  @param N     dimensionality of N-sphere (number of phases)  
 *  @param bias  use the bias in phases? 
 */
// ============================================================================
Ostap::Math::NSphere::NSphere 
( const unsigned short N       ,
  const unsigned short rotated ) 
  : Ostap::Math::NSphere::NSphere ( std::vector<double> ( N , 0.0 )   , 
                                    rotated < N + 1 ? rotated : N + 1 ) 
{}
// ============================================================================
/*  Standard constructor
 *  @param nPhases  dimensionality of N-sphere 
 *  @param bias     use the rotated sphere? 
 */
// ============================================================================
Ostap::Math::NSphere::NSphere 
( const std::vector<double>& phases  )
  : Ostap::Math::NSphere::NSphere ( phases , phases.size() + 1 ) 
{}
// ============================================================================
/*  Standard constructor with deltas 
 *  @param phases  vector of phases 
 *  @param deltas  rotation deltas 
 */
// ============================================================================
Ostap::Math::NSphere::NSphere  
( const std::vector<double>& phases ,
  const std::vector<double>& deltas ) 
  : m_delta   ( phases.size () ) 
  , m_phases  ( phases         )
  , m_sin_phi ( phases.size () , 0 ) 
  , m_cos_phi ( phases.size () , 1 ) 
{
  // copy deltas 
  const unsigned int nd = deltas.size() ;
  for ( unsigned short i = 0 ; i < m_phases.size() ; ++i )
  { m_delta [ i ] =  i < nd ? deltas  [ i ] : 0.0 ; }
  //
  for ( unsigned short  i = 0 ; i < m_phases.size() ; ++i )
  {
    const double phase = m_phases [i] + m_delta[i] ;
    const std::pair<double,double> sincos = _sincos_ ( phase );
    m_sin_phi [i] = sincos.first ;
    m_cos_phi [i] = sincos.second ;
  }
}
// ============================================================================
/*  Standard constructor with deltas 
 *  @param phases  vector of phases 
 *  @param deltas  rotation deltas 
 */
// ============================================================================
Ostap::Math::NSphere::NSphere  
( const std::string&        /* fake */,
  const std::vector<double>& deltas ) 
  : m_delta   ( deltas ) 
  , m_phases  ( deltas.size () , 0.0 )
  , m_sin_phi ( deltas.size () , 0.0 ) 
  , m_cos_phi ( deltas.size () , 1.0 ) 
{
  for ( unsigned short  i = 0 ; i < m_phases.size() ; ++i )
  {
    const double phase = m_phases [i] + m_delta[i] ;
    const std::pair<double,double> sincos = _sincos_ ( phase );
    m_sin_phi [i] = sincos.first ;
    m_cos_phi [i] = sincos.second ;
  }
}
// ============================================================================
Ostap::Math::NSphere::NSphere 
( const std::vector<double>& phases , 
  const unsigned short       rotated ) 
  : m_delta   ( phases.size () , 0 ) 
  , m_phases  ( phases ) 
  , m_sin_phi ( phases.size () , 0 ) 
  , m_cos_phi ( phases.size () , 1 ) 
{ 
  // ==========================================================================
  // calculate the bias (if needed) 
  if  ( rotated )                       // ROTATE SPHERE 
  {
    const unsigned int nzero = m_phases.size() - rotated + 1 ;
    for ( unsigned short  i = 0 ; i < m_phases.size() ; ++i )
    {
      const double ni    = m_phases.size () - i ;
      //
      if ( i <  nzero ) { m_delta [i] = s_PIHALF ; }
      else 
      { m_delta [i] = std::atan2 ( std::sqrt ( ni ) , 1.0L ) ; }
      //
      const double phase = m_phases [i] + m_delta[i] ;
      //
      const std::pair<double,double> sincos = _sincos_ ( phase );
      m_sin_phi [i] = sincos.first ;
      m_cos_phi [i] = sincos.second ;
    }
  }
  //
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Math::NSphere::~NSphere() {}
// ============================================================================
// set new value for phi(i)      0 <= i < nPhi
// ============================================================================
bool Ostap::Math::NSphere::setPhase 
( const unsigned short index , 
  const double         value ) 
{
  // 
  if ( nPhi () <= index ) { return false ; } // no change in unphysical phases 
  //
  if ( s_equal ( m_phases [ index ] , value ) ) { return false ; }
  //
  const double phase  = value + m_delta [ index ] ;
  //
  const std::pair<double,double> sincos = _sincos_ ( phase );
  m_sin_phi [ index ] = sincos.first  ;
  m_cos_phi [ index ] = sincos.second ;
  m_phases  [ index ] = value         ;  // attention!! original values!! 
  //
  return true ;
}
// ============================================================================
// copy assignement 
// ============================================================================
Ostap::Math::NSphere& 
Ostap::Math::NSphere::operator=( const Ostap::Math::NSphere&  right ) 
{
  if ( &right == this ) { return *this ; }
  //
  m_delta     = right.m_delta   ;
  m_phases    = right.m_phases  ;
  m_sin_phi   = right.m_sin_phi ;
  m_cos_phi   = right.m_cos_phi ;
  //
  return *this ;
}
// ============================================================================
// move assignement 
// ============================================================================
Ostap::Math::NSphere& 
Ostap::Math::NSphere::operator=(       Ostap::Math::NSphere&& right ) 
{
  if ( &right == this ) { return *this ; }
  //
  m_delta     = std::move ( right.m_delta   ) ;
  m_phases    = std::move ( right.m_phases  ) ;
  m_sin_phi   = std::move ( right.m_sin_phi ) ;
  m_cos_phi   = std::move ( right.m_cos_phi ) ;
  //
  return *this ;
}
// ============================================================================
// swap two spheres 
// ============================================================================
void Ostap::Math::NSphere::swap ( Ostap::Math::NSphere& right )  // swap two spheres 
{
  //
  if ( &right == this ) { return ; }
  //
  std::swap ( m_delta   , right.m_delta   ) ;
  std::swap ( m_phases  , right.m_phases  ) ;
  std::swap ( m_sin_phi , right.m_sin_phi ) ;
  std::swap ( m_cos_phi , right.m_cos_phi ) ;
}
// ============================================================================
/* convert n-coordinates \f$ x_i \f$ into (n-1) phases \f$ \phi_i\f$  
 * @param  x vector  in n-dimensional space 
 * @return vector of spherical phases
 */
// ============================================================================
std::vector<double> 
Ostap::Math::NSphere::phis ( const std::vector<double>& x )
{
  //
  if      ( 0 == x.size () ) { return std::vector<double> ()          ; }
  else if ( 1 == x.size () ) { return std::vector<double> ( 1 , 0.0 ) ; }
  //
  const unsigned short nphi = x.size() - 1 ;
  //
  std::vector<long double> r2 ( x.begin() , x.end () ) ;
  std::transform ( r2.begin () , r2.end() , r2.begin () , 
                   []( long double  x ) -> long double { return x * x ; } ) ;
  std::partial_sum ( r2.rbegin() , r2.rend() , r2.rbegin() ) ;
  //
  if  ( s_zero ( r2.front () ) ) { return std::vector<double> ( nphi , 0 ) ; }
  
  std::vector<double> phis ( nphi , 0 ) ;
  
  for  ( unsigned short i = 0 ; i + 1 < nphi ; ++i  ) 
  {
    //
    const long double rr = r2 [ i ] ;
    //
    if ( rr <= 0 || s_zero ( rr ) ) { phis [ i ] = 0 ; continue ; }
    //
    const long double r =   std::sqrt ( rr ) ;
    if ( s_zero ( r )             ) { phis [ i ] = 0 ; continue ; }
    //
    const long double xi = x [ i ] ;
    const long double ax = std::abs ( xi  ) ;
    if ( r < ax || s_equal ( ax , r ) ) 
    { phis [ i ] = ( 0 < xi ) ? 0 : M_PIl ; continue ; }
    //
    phis [ i ] = std::acos ( x [ i ] / r ) ;
  }
  // the last   case 
  const unsigned short kphi = nphi - 1 ;
  const long double rr      = r2 [kphi] ;
  if   ( s_zero  ( rr ) ) { phis [ kphi ] = 0 ; }
  else 
  {
    const long double r  = std::sqrt ( rr ) ;
    const long double xk = x [ kphi ] ;
    const long double ax = std::abs ( xk  ) ;
    if ( r < ax || s_zero ( r ) ) { phis [ kphi ] = 0 ; }
    else 
    {
      phis [ kphi ] =  
        0 <= x [ nphi ] ? 
        std::acos ( x [ kphi ] / r ) : 2 * M_PI - std::acos ( x [ kphi ] / r ) ;
    }     
  }
  return phis ;
}  
// =============================================================================
// get unique tag for the given sphere 
// =============================================================================
std::size_t Ostap::Math::NSphere::tag () const 
{
  static const std::string s_name { "NSPhere" } ;
  return Ostap::Utils::hash_combiner 
    ( s_name , nPhi() , 
      Ostap::Utils::hash_range ( m_delta  ) ,  
      Ostap::Utils::hash_range ( m_phases ) ) ;
} 

// ============================================================================
//                                                                      The END 
// ============================================================================
