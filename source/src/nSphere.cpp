// ============================================================================
// Include files 
// ============================================================================
// STD & STL
// ============================================================================
#include <vector>
#include <cmath>
// ============================================================================
// local
// ============================================================================
#include "Ostap/Math.h"
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
    const double sinv = std::sin ( phase ) ;
    const double cosv = std::cos ( phase ) ;
    //
    const double acosv = std::abs ( cosv ) ;
    if      ( 0 == sinv  || 1 == acosv ) 
    { return std::make_pair ( 0  , 0 < cosv ? 1 : -1 ) ; }
    //
    const double asinv = std::abs ( sinv ) ;
    if      ( 0 == cosv  || 1 == asinv )
    { return std::make_pair ( 0 < sinv ? 1 : -1 , 0 ) ; }
    //
    if      ( s_zero  ( sinv ) || s_equal ( acosv , 1 ) )
    { return std::make_pair ( 0  , 0 < cosv ? 1 : -1 ) ; }
    //
    if      ( s_zero  ( cosv ) || s_equal ( asinv , 1 ) )
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
Ostap::Math::NSphere::NSphere 
( const std::vector<double>& phases , 
  const unsigned short       rotated ) 
  : m_rotated ( rotated < phases.size () + 1 ? rotated : phases.size() + 1 )
  , m_delta   ( phases.size () , 0 ) 
  , m_phases  ( phases ) 
  , m_sin_phi ( phases.size () , 0 ) 
  , m_cos_phi ( phases.size () , 1 ) 
{ 
  // ==============================
  // calculate the bias (if needed) 
  if  ( m_rotated )                       // ROTATE SPHERE 
  {
    const unsigned int nzero = m_phases.size() - m_rotated + 1 ;
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
// copy
// ============================================================================
Ostap::Math::NSphere::NSphere 
( const Ostap::Math::NSphere&  right ) 
  : m_rotated  ( right.m_rotated ) 
  , m_delta    ( right.m_delta   ) 
  , m_phases   ( right.m_phases  ) 
  , m_sin_phi  ( right.m_sin_phi ) 
  , m_cos_phi  ( right.m_cos_phi ) 
{}
// ============================================================================
// move
// ============================================================================
Ostap::Math::NSphere::NSphere 
(       Ostap::Math::NSphere&& right ) 
  : m_rotated  ( right.m_rotated ) 
  , m_delta    ( std::move ( right.m_delta   ) ) 
  , m_phases   ( std::move ( right.m_phases  ) ) 
  , m_sin_phi  ( std::move ( right.m_sin_phi ) )  
  , m_cos_phi  ( std::move ( right.m_cos_phi ) ) 
{}
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
  if ( nPhi() <= index ) { return false ; } // no change in unphysical phases 
  //
  if ( s_equal ( m_phases[index] , value ) ) { return false ; }
  //
  const double di     = m_rotated ? m_delta [ index ] : 0.0 ;
  const double phase  = value + di ;
  //
  const std::pair<double,double> sincos = _sincos_ ( phase );
  m_sin_phi [ index ] = sincos.first  ;
  m_cos_phi [ index ] = sincos.second ;
  m_phases  [ index ] = value ;  // attention!! original values!! 
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
  m_rotated   = right.m_rotated ;
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
  m_rotated   = right.m_rotated ;
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
  std::swap ( m_rotated , right.m_rotated ) ;
  std::swap ( m_delta   , right.m_delta   ) ;
  std::swap ( m_phases  , right.m_phases  ) ;
  std::swap ( m_sin_phi , right.m_sin_phi ) ;
  std::swap ( m_cos_phi , right.m_cos_phi ) ;
}


// ============================================================================
// The END 
// ============================================================================
