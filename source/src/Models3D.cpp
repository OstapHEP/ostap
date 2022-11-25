// Include files 
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Hash.h"
#include "Ostap/Math.h"
#include "Ostap/Vector3DTypes.h"
#include "Ostap/Models3D.h"
#include "Ostap/MoreMath.h"
// ============================================================================
//  Local 
// ============================================================================
#include "Exception.h"
#include "local_math.h"
#include "local_hash.h"
#include "Integrator1D.h"
#include "Integrator2D.h"
#include "Integrator3D.h"
// ============================================================================
/** @file 
 *  Implementation file for the classes form the file Ostap/Models3D.h 
 *  @date 2022-06-23 
 *  @author Vanya Belyaev  Ivan.Belyaev@cern.ch 
 */
// ============================================================================
namespace 
{
  // // ==========================================================================
  typedef  std::array<short,11> SPLITS ;
  /// split points for 3D Gaussian
  static const SPLITS s_SPLITS = { -15 , -9 , -5 , -3 , -1 , 0 , 1 , 3 , 5 , 9 , 15 } ;
  // ==========================================================================
}
// ============================================================================
/* constructor 
 *  @param muX    x-location 
 *  @param muY    y-location 
 *  @param muZ    z-location 
 *  @param sigmaX x-sigma
 *  @param sigmaY y-sigma
 *  @param sigmaZ z-sigma
 *  @param phi   Euler angle phi 
 *  @param theta Euler angle theta 
 *  @param psi   Euler angle psi 
 */
// ============================================================================
Ostap::Math::Gauss3D::Gauss3D 
( const double muX    , 
  const double muY    , 
  const double muZ    , 
  const double sigmaX , 
  const double sigmaY ,
  const double sigmaZ ,
  const double phi    ,
  const double theta  ,
  const double psi    ) 
  : m_muX    ( muX ) 
  , m_muY    ( muY ) 
  , m_muZ    ( muZ ) 
  , m_sigmaX ( std::abs ( sigmaX ) ) 
  , m_sigmaY ( std::abs ( sigmaY ) ) 
  , m_sigmaZ ( std::abs ( sigmaZ ) ) 
  , m_phi    ( -123456 ) 
  , m_theta  ( -234567 ) 
  , m_psi    ( -345678 ) 
{
  setEuler ( phi , theta , psi ) ;
}
// ============================================================================
// set mux-parameter
// ============================================================================
bool Ostap::Math::Gauss3D::setMuX ( const double value )
{
  if ( s_equal ( m_muX , value ) ) { return false ; }
  m_muX = value ;
  return true ;
}
// ============================================================================
// set muy-parameter
// ============================================================================
bool Ostap::Math::Gauss3D::setMuY ( const double value )
{
  if ( s_equal ( m_muY , value ) ) { return false ; }
  m_muY = value ;
  return true ;
}
// ============================================================================
// set muz-parameter
// ============================================================================
bool Ostap::Math::Gauss3D::setMuZ ( const double value )
{
  if ( s_equal ( m_muZ , value ) ) { return false ; }
  m_muZ = value ;
  return true ;
}
// ============================================================================
// set sigmax-parameter
// ============================================================================
bool Ostap::Math::Gauss3D::setSigmaX ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_sigmaX , avalue ) ) { return false ; }
  m_sigmaX = avalue ;
  return true ;
}
// ============================================================================
// set sigmay-parameter
// ============================================================================
bool Ostap::Math::Gauss3D::setSigmaY ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_sigmaY , avalue ) ) { return false ; }
  m_sigmaY = avalue ;
  return true ;
}
// ============================================================================
// set sigmaz-parameter
// ============================================================================
bool Ostap::Math::Gauss3D::setSigmaZ ( const double value )
{
  const double avalue = std::abs ( value ) ;
  if ( s_equal ( m_sigmaZ , avalue ) ) { return false ; }
  m_sigmaZ = avalue ;
  return true ;
}
// ============================================================================
// set phi-parameter
// ============================================================================
bool Ostap::Math::Gauss3D::setPhi ( const double value )
{ 
  if ( s_equal ( m_phi , value ) ) { return false ; }  
  return setEuler ( value , m_theta , m_psi ) ;
}
// ============================================================================
// set theta-parameter
// ============================================================================
bool Ostap::Math::Gauss3D::setTheta ( const double value )
{ 
  if ( s_equal ( m_theta , value ) ) { return false ; }  
  return setEuler ( m_phi , value , m_psi ) ;
}
// ============================================================================
// set psi-parameter
// ============================================================================
bool Ostap::Math::Gauss3D::setPsi ( const double value )
{ 
  if ( s_equal ( m_psi , value ) ) { return false ; }  
  return setEuler ( m_phi , m_theta , value ) ;
}
// ============================================================================
// set all Euler parameters 
// ============================================================================
bool Ostap::Math::Gauss3D::setEuler   
( const double phi    , 
  const double theta  , 
  const double psi    ) 
{
  if ( s_equal ( m_phi   , phi   ) && 
       s_equal ( m_theta , theta ) && 
       s_equal ( m_psi   , psi   ) ) { return false ; }  
  //
  const ROOT::Math::EulerAngles ea ( phi , theta , psi );
  m_phi      = ea.Phi   () ;
  m_theta    = ea.Theta () ;
  m_psi      = ea.Psi   () ;
  m_rotation = ROOT::Math::Rotation3D ( ea ) ;
  //
  return true ;
}
// ============================================================================
// set all Euler parameters 
// ============================================================================
bool Ostap::Math::Gauss3D::setEuler   
( const ROOT::Math::EulerAngles& angles ) 
{ return setEuler ( angles.Phi() , angles.Theta() , angles.Psi () ) ; }
// ============================================================================
// set all Euler parameters 
// ============================================================================
bool Ostap::Math::Gauss3D::setEuler   
( const ROOT::Math::Rotation3D& angles ) 
{ return setEuler ( ROOT::Math::EulerAngles ( angles ) ) ; }

// ============================================================================
// get the value
// ============================================================================
double Ostap::Math::Gauss3D::operator ()
  ( const double x , 
    const double y , 
    const double z ) const 
{
  // (1) rotate  
  const Ostap::Vector3D dv { m_rotation ( Ostap::Vector3D ( x - m_muX , 
                                                            y - m_muY , 
                                                            z - m_muZ ) ) } ;
  
  // (2) evaluate 3D gaussian  
  static const double s_norm = std::pow ( 2 * M_PI , 1.5 ) ;
  //
  return std::exp ( -0.5 * ( std::pow ( dv.X() / m_sigmaX , 2 )  + 
                             std::pow ( dv.Y() / m_sigmaY , 2 )  + 
                             std::pow ( dv.Z() / m_sigmaZ , 2 ) ) )
    / ( s_norm * m_sigmaX * m_sigmaY  * m_sigmaZ ) ;
}
// ============================================================================
// get the integral over the whole 3D-region
// ============================================================================
double Ostap::Math::Gauss3D::integral () const { return 1 ; }
// ============================================================================
/*  get the integral over 2D-region
 *  \f[ \int_{x_{low}}^{x_{high}}
 *      \int_{y_{low}}^{y_{high}}
 *      \int_{z_{low}}^{z_{high}}
 *        \mathcal{B}(x,y,z) \mathrm{d}x\mathrm{d}y\mathrm{d}z \f]
 *  @param xlow  low  edge in x
 *  @param xhigh high edge in x
 *  @param ylow  low  edge in y
 *  @param yhigh high edge in y
 *  @param zlow  low  edge in z
 *  @param zhigh high edge in z
 */
// ============================================================================
double Ostap::Math::Gauss3D::integral 
( const double xlow  ,
  const double xhigh ,
  const double ylow  , 
  const double yhigh ,
  const double zlow  , 
  const double zhigh ) const 
{
  if ( s_equal ( xlow , xhigh ) ||
       s_equal ( ylow , yhigh ) ||
       s_equal ( zlow , zhigh ) ) { return 0 ; }
  //
  if      ( xhigh < xlow ) { return -integral ( xhigh , xlow  , ylow  , yhigh , zlow  , zhigh ) ; }
  else if ( yhigh < ylow ) { return -integral ( xlow  , xhigh , yhigh , ylow  , zlow  , zhigh ) ; }
  else if ( zhigh < zlow ) { return -integral ( xlow  , xhigh , ylow  , yhigh , zhigh , zlow  ) ; }
  //
  if ( ( s_zero  ( m_phi   ) && 
         s_zero  ( m_theta ) && 
         s_zero  ( m_psi   ) ) || 
       ( s_equal ( m_sigmaX , m_sigmaY ) && 
         s_equal ( m_sigmaY , m_sigmaZ ) ) ) 
  {
    return 
      ( Ostap::Math::gauss_cdf ( xhigh , m_muX , m_sigmaX ) - 
        Ostap::Math::gauss_cdf ( xlow  , m_muX , m_sigmaX ) ) * 
      ( Ostap::Math::gauss_cdf ( yhigh , m_muY , m_sigmaY ) - 
        Ostap::Math::gauss_cdf ( ylow  , m_muY , m_sigmaY ) ) *
      ( Ostap::Math::gauss_cdf ( zhigh , m_muZ , m_sigmaZ ) - 
        Ostap::Math::gauss_cdf ( zlow  , m_muZ , m_sigmaZ ) ) ;  
  }
  //
  const Ostap::Vector3D v1 { m_rotation ( Ostap::Vector3D ( m_sigmaX , 0 , 0 ) ) } ;
  const Ostap::Vector3D v2 { m_rotation ( Ostap::Vector3D ( 0 , m_sigmaY , 0 ) ) } ;
  const Ostap::Vector3D v3 { m_rotation ( Ostap::Vector3D ( 0 , 0 , m_sigmaZ ) ) } ;
  //
  const double sx = std::max ( std::abs ( v1.X() ) , std::max ( std::abs ( v2.X() ) , std::abs ( v3.X() ) ) ) ;
  const double sy = std::max ( std::abs ( v1.Y() ) , std::max ( std::abs ( v2.Y() ) , std::abs ( v3.Y() ) ) ) ;
  const double sz = std::max ( std::abs ( v1.Z() ) , std::max ( std::abs ( v2.Z() ) , std::abs ( v3.Z() ) ) ) ;
  //
  if      ( xhigh <= m_muX - 50 * sx ) { return 0 ; }
  else if ( xlow  >= m_muX + 50 * sx ) { return 0 ; }
  else if ( yhigh <= m_muY - 50 * sy ) { return 0 ; }
  else if ( ylow  >= m_muY + 50 * sy ) { return 0 ; }
  else if ( zhigh <= m_muZ - 50 * sz ) { return 0 ; }
  else if ( zlow  >= m_muZ + 50 * sz ) { return 0 ; }
  //  
  // split into smaller regions 
  //
  if ( xlow < m_muX + sx * s_SPLITS.back() || xhigh > m_muX + sx * s_SPLITS.front() ) 
  {
    for ( SPLITS::const_iterator ix = s_SPLITS.begin() ; s_SPLITS.end() != ix  ; ++ix ) 
    {
      const double px = m_muX + (*ix) * sx ;
      if ( xlow < px && px < xhigh ) 
      {
      return 
        integral ( xlow , px    , ylow , yhigh , zlow , zhigh ) +
        integral ( px   , xhigh , ylow , yhigh , zlow , zhigh ) ;
      }
    } 
  }
  // ==========================================================================
  if ( ylow < m_muY + sy * s_SPLITS.back() || yhigh > m_muY + sy * s_SPLITS.front() ) 
  {
    for ( SPLITS::const_iterator iy = s_SPLITS.begin() ; s_SPLITS.end() != iy  ; ++iy ) 
    {
      const double py = m_muY + (*iy) * sy ;
      if ( ylow < py && py < yhigh ) 
      {
        return 
          integral ( xlow , xhigh , ylow , py    , zlow , zhigh ) +
          integral ( xlow , xhigh , py   , yhigh , zlow , zhigh ) ;
      }
    }
  }
  // ==========================================================================
  if ( zlow < m_muZ + sz * s_SPLITS.back() || zhigh > m_muZ + sz * s_SPLITS.front() ) 
  {
    for ( SPLITS::const_iterator iz = s_SPLITS.begin() ; s_SPLITS.end() != iz  ; ++iz ) 
    {
      const double pz = m_muZ + (*iz) * sz ;
      if ( zlow < pz && pz < zhigh ) 
      {
        return 
          integral ( xlow , xhigh , ylow , yhigh , zlow , pz    ) +
          integral ( xlow , xhigh , ylow , yhigh , pz   , zhigh ) ;
      }
    }
  }
  // ==========================================================================
  const bool in_tail = 
    ( xhigh <= m_muX + sx * s_SPLITS.front () ) || 
    ( xlow  >= m_muX + sx * s_SPLITS.back  () ) || 
    ( yhigh <= m_muY + sy * s_SPLITS.front () ) || 
    ( ylow  >= m_muY + sy * s_SPLITS.back  () ) ||
    ( zhigh <= m_muZ + sz * s_SPLITS.front () ) || 
    ( zlow  >= m_muZ + sz * s_SPLITS.back  () ) ;
  //
  // use 3D-cubature   
  //
  static const Ostap::Math::GSL::Integrator3D<Gauss3D> s_cubature{} ;
  static const char s_message[] = "Integral(Gauss3D)" ;
  const auto F = s_cubature.make_function ( this , 
                                            xlow , xhigh , 
                                            ylow , yhigh , 
                                            zlow , zhigh ) ;
  //
  int    ierror  =  0 ;
  double  result =  1 ;
  double  error  = -1 ;
  std::tie ( ierror , result , error ) = s_cubature.cubature
    ( tag () , &F , 50000 , 
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision 
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      s_message , __FILE__ , __LINE__ ) ;
  return  result ;
}
// ======================================================================
/*  integral over x&y-dimensions
 *  \f[ \int_{x_{low}}^{x_{high}}
 *      \int_{y_{low}}^{y_{high}} \mathcal{B}(x,y,z) \mathrm{d}x\mathrm{d}y\f]
 *  @param z     variable
 *  @param xlow  low  edge in x
 *  @param xhigh high edge in x
 *  @param ylow  low  edge in y
 *  @param yhigh high edge in y
 */
// ======================================================================
double Ostap::Math::Gauss3D::integrateXY 
( const double z    ,                          
  const double xlow , const double xhigh ,
  const double ylow , const double yhigh ) const 
{
  if ( s_equal ( xlow , xhigh ) || s_equal ( ylow , yhigh ) ) { return 0 ; }
  //
  if      ( xhigh < xlow ) { return -integrateXY ( z , xhigh , xlow  , ylow  , yhigh ) ; }
  else if ( yhigh < ylow ) { return -integrateXY ( z , xlow  , xhigh , yhigh , ylow  ) ; }
  //
  if ( ( s_zero  ( m_phi   ) && 
         s_zero  ( m_theta ) && 
         s_zero  ( m_psi   ) ) || 
       ( s_equal ( m_sigmaX , m_sigmaY ) && 
         s_equal ( m_sigmaY , m_sigmaZ ) ) ) 
  {
    return 
      ( Ostap::Math::gauss_cdf ( xhigh , m_muX , m_sigmaX ) - 
        Ostap::Math::gauss_cdf ( xlow  , m_muX , m_sigmaX ) ) * 
      ( Ostap::Math::gauss_cdf ( yhigh , m_muY , m_sigmaY ) - 
        Ostap::Math::gauss_cdf ( ylow  , m_muY , m_sigmaY ) ) 
      * Ostap::Math::gauss_pdf ( z     , m_muZ , m_sigmaZ ) ;                             
  }
  //
  const Ostap::Vector3D v1 { m_rotation ( Ostap::Vector3D ( m_sigmaX , 0 , 0 ) ) } ;
  const Ostap::Vector3D v2 { m_rotation ( Ostap::Vector3D ( 0 , m_sigmaY , 0 ) ) } ;
  const Ostap::Vector3D v3 { m_rotation ( Ostap::Vector3D ( 0 , 0 , m_sigmaZ ) ) } ;
  //
  const double sx = std::max ( std::abs ( v1.X() ) , std::max ( std::abs ( v2.X() ) , std::abs ( v3.X() ) ) ) ;
  const double sy = std::max ( std::abs ( v1.Y() ) , std::max ( std::abs ( v2.Y() ) , std::abs ( v3.Y() ) ) ) ;
  const double sz = std::max ( std::abs ( v1.Z() ) , std::max ( std::abs ( v2.Z() ) , std::abs ( v3.Z() ) ) ) ;
  //
  if      ( xhigh <= m_muX - 50 * sx ) { return 0 ; }
  else if ( xlow  >= m_muX + 50 * sx ) { return 0 ; }
  else if ( yhigh <= m_muY - 50 * sy ) { return 0 ; }
  else if ( ylow  >= m_muY + 50 * sy ) { return 0 ; }
  else if ( z     <= m_muZ - 50 * sz ) { return 0 ; }
  else if ( z     >= m_muZ + 50 * sz ) { return 0 ; }
  //  
  // split into smaller regions 
  //
  if ( xlow < m_muX + sx * s_SPLITS.back() || xhigh > m_muX + sx * s_SPLITS.front() ) 
  {
    for ( SPLITS::const_iterator ix = s_SPLITS.begin() ; s_SPLITS.end() != ix  ; ++ix ) 
    {
      const double xsplit = m_muX + (*ix) * sx ;
      if ( xlow < xsplit && xsplit < xhigh ) 
      {
        return 
          integrateXY ( z , xlow   , xsplit , ylow , yhigh ) +
          integrateXY ( z , xsplit , xhigh  , ylow , yhigh ) ;
      }
    } 
  }
  // ==========================================================================
  if ( ylow < m_muY + sy * s_SPLITS.back() || yhigh > m_muY + sy * s_SPLITS.front() ) 
  {
    for ( SPLITS::const_iterator iy = s_SPLITS.begin() ; s_SPLITS.end() != iy  ; ++iy ) 
    {
      const double ysplit = m_muY + (*iy) * sy ;
      if ( ylow < ysplit && ysplit < yhigh ) 
      {
        return 
          integrateXY ( z , xlow , xhigh , ylow   , ysplit ) +
          integrateXY ( z , xlow , xhigh , ysplit , yhigh  ) ;
      }
    }
  }
  // ==========================================================================
  // 
  const bool in_tail = 
    ( xhigh <= m_muX + sx * s_SPLITS.front () ) || 
    ( xlow  >= m_muX + sx * s_SPLITS.back  () ) || 
    ( yhigh <= m_muY + sy * s_SPLITS.front () ) || 
    ( ylow  >= m_muY + sy * s_SPLITS.back  () ) ||
    ( z     <= m_muZ + sz * s_SPLITS.front () ) || 
    ( z     >= m_muZ + sz * s_SPLITS.back  () ) ;
  //
  // use 2D-cubature   
  //
  typedef Ostap::Math::IntegrateXY<Gauss3D> FXY ;
  const FXY fxy ( this , z ) ;
  //
  static const Ostap::Math::GSL::Integrator2D<FXY> s_cubature{} ;
  static const char s_message[] = "IntegralXY(Gauss3D)" ;
  const auto F = s_cubature.make_function ( &fxy , 
                                            xlow , xhigh , 
                                            ylow , yhigh ) ;
  //
  int    ierror  =  0 ;
  double  result =  1 ;
  double  error  = -1 ;
  std::tie ( ierror , result , error ) = s_cubature.cubature
    ( Ostap::Utils::hash_combiner ( tag () , 'Z' , z ) , &F , 50000 , 
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      s_message , __FILE__ , __LINE__ ) ;

  return  result ;
}
// ============================================================================
/*  integral over y&z-dimensions
 *  \f[ \int_{y_{low}}^{y_{high}}
 *      \int_{z_{low}}^{z_{high}} \mathcal{B}(x,y,z) \mathrm{d}y\mathrm{d}z\f]
 *  @param x     variable
 *  @param ylow  low  edge in x
 *  @param yhigh high edge in x
 *  @param zlow  low  edge in y
 *  @param zhigh high edge in y
 */
// ============================================================================
double Ostap::Math::Gauss3D::integrateYZ 
( const double x    ,                          
  const double ylow , const double yhigh ,
  const double zlow , const double zhigh ) const 
{
  if ( s_equal ( zlow , zhigh ) ||
       s_equal ( ylow , yhigh ) ) { return 0 ; }
  //
  if      ( zhigh < zlow ) { return -integrateYZ ( x , ylow  , yhigh , zhigh , zlow  ) ; }
  else if ( yhigh < ylow ) { return -integrateYZ ( x , yhigh , ylow  , zlow  , zhigh ) ; }
  //
  if ( ( s_zero  ( m_phi   ) && 
         s_zero  ( m_theta ) && 
         s_zero  ( m_psi   ) ) || 
       ( s_equal ( m_sigmaX , m_sigmaY ) && 
         s_equal ( m_sigmaY , m_sigmaZ ) ) ) 
  {
    return 
      ( Ostap::Math::gauss_cdf ( zhigh , m_muZ , m_sigmaZ ) - 
        Ostap::Math::gauss_cdf ( zlow  , m_muZ , m_sigmaZ ) ) * 
      ( Ostap::Math::gauss_cdf ( yhigh , m_muY , m_sigmaY ) - 
        Ostap::Math::gauss_cdf ( ylow  , m_muY , m_sigmaY ) ) 
      * Ostap::Math::gauss_pdf ( x     , m_muX , m_sigmaX ) ;                             
  }
  // 
  const Ostap::Vector3D v1 { m_rotation ( Ostap::Vector3D ( m_sigmaX , 0 , 0 ) ) } ;
  const Ostap::Vector3D v2 { m_rotation ( Ostap::Vector3D ( 0 , m_sigmaY , 0 ) ) } ;
  const Ostap::Vector3D v3 { m_rotation ( Ostap::Vector3D ( 0 , 0 , m_sigmaZ ) ) } ;
  //
  const double sx = std::max ( std::abs ( v1.X() ) , std::max ( std::abs ( v2.X() ) , std::abs ( v3.X() ) ) ) ;
  const double sy = std::max ( std::abs ( v1.Y() ) , std::max ( std::abs ( v2.Y() ) , std::abs ( v3.Y() ) ) ) ;
  const double sz = std::max ( std::abs ( v1.Z() ) , std::max ( std::abs ( v2.Z() ) , std::abs ( v3.Z() ) ) ) ;
  //
  if      ( x     <= m_muX - 50 * sx ) { return 0 ; }
  else if ( x     >= m_muX + 50 * sx ) { return 0 ; }
  else if ( yhigh <= m_muY - 50 * sy ) { return 0 ; }
  else if ( ylow  >= m_muY + 50 * sy ) { return 0 ; }
  else if ( zhigh <= m_muZ - 50 * sz ) { return 0 ; }
  else if ( zlow  >= m_muZ + 50 * sz ) { return 0 ; }
  //  
  // split into smaller regions 
  //
  if ( zlow < m_muZ + sz * s_SPLITS.back() || zhigh > m_muZ + sz * s_SPLITS.front() ) 
  {
    for ( SPLITS::const_iterator iz = s_SPLITS.begin() ; s_SPLITS.end() != iz  ; ++iz ) 
    {
      const double pz = m_muZ + (*iz) * sz ;
      if ( zlow < pz && pz < zhigh ) 
      {
        return 
          integrateYZ ( x , ylow , yhigh , zlow , pz    ) +
          integrateYZ ( x , ylow , yhigh , pz   , zhigh ) ;
      }
    } 
  }
  // ==========================================================================
  if ( ylow < m_muY + sy * s_SPLITS.back() || yhigh > m_muY + sy * s_SPLITS.front() ) 
  {
    for ( SPLITS::const_iterator iy = s_SPLITS.begin() ; s_SPLITS.end() != iy  ; ++iy ) 
    {
      const double py = m_muY + (*iy) * sy ;
      if ( ylow < py && py < yhigh ) 
      {
        return 
          integrateYZ ( x , ylow , py    , zlow , zhigh ) +
          integrateYZ ( x , py   , yhigh , zlow , zhigh ) ;
      }
    }
  }
  // ==========================================================================
  const bool in_tail =
    ( x     <= m_muX + sx * s_SPLITS.front () ) || 
    ( x     >= m_muX + sx * s_SPLITS.back  () ) || 
    ( yhigh <= m_muY + sy * s_SPLITS.front () ) || 
    ( ylow  >= m_muY + sy * s_SPLITS.back  () ) ||
    ( zhigh <= m_muZ + sz * s_SPLITS.front () ) || 
    ( zlow  >= m_muZ + sz * s_SPLITS.back  () ) ;
  //
  // use 2D-cubature   
  //
  typedef Ostap::Math::IntegrateYZ<Gauss3D> FYZ ;
  const FYZ fyz ( this , x ) ;
  //
  static const Ostap::Math::GSL::Integrator2D<FYZ> s_cubature{} ;
  static const char s_message[] = "IntegralYZ(Gauss3D)" ;
  const auto F = s_cubature.make_function ( &fyz , 
                                            ylow , yhigh , 
                                            zlow , zhigh ) ;
  //
  int    ierror  =  0 ;
  double  result =  1 ;
  double  error  = -1 ;
  std::tie ( ierror , result , error ) = s_cubature.cubature
    ( Ostap::Utils::hash_combiner ( tag () , 'X' , x ) , &F , 50000 ,
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision 
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      s_message , __FILE__ , __LINE__ ) ;
  return  result ;
}
// ============================================================================
/*  integral over y&z-dimensions
 *  \f[ \int_{y_{low}}^{y_{high}}
 *      \int_{z_{low}}^{z_{high}} \mathcal{B}(x,y,z) \mathrm{d}y\mathrm{d}z\f]
 *  @param x     variable
 *  @param ylow  low  edge in x
 *  @param yhigh high edge in x
 *  @param zlow  low  edge in y
 *  @param zhigh high edge in y
 */
 // ============================================================================
double Ostap::Math::Gauss3D::integrateXZ
( const double y    ,                          
  const double xlow , const double xhigh ,
  const double zlow , const double zhigh ) const 
{
  if ( s_equal ( xlow , xhigh ) ||
       s_equal ( zlow , zhigh ) ) { return 0 ; }
  //
  if      ( xhigh < xlow ) { return -integrateXZ ( y , xhigh , xlow  , zlow  , zhigh ) ; }
  else if ( zhigh < zlow ) { return -integrateXZ ( y , xlow  , xhigh , zhigh , zlow  ) ; }
  //
  if ( ( s_zero  ( m_phi   ) && 
         s_zero  ( m_theta ) && 
         s_zero  ( m_psi   ) ) || 
       ( s_equal ( m_sigmaX , m_sigmaY ) && 
         s_equal ( m_sigmaY , m_sigmaZ ) ) ) 
  {
    return 
      ( Ostap::Math::gauss_cdf ( zhigh , m_muZ , m_sigmaZ ) - 
        Ostap::Math::gauss_cdf ( zlow  , m_muZ , m_sigmaZ ) ) * 
      ( Ostap::Math::gauss_cdf ( xhigh , m_muX , m_sigmaX ) - 
        Ostap::Math::gauss_cdf ( xlow  , m_muX , m_sigmaX ) ) 
      * Ostap::Math::gauss_pdf ( y     , m_muY , m_sigmaY ) ;                             
  }
  // 
  const Ostap::Vector3D v1 { m_rotation ( Ostap::Vector3D ( m_sigmaX , 0 , 0 ) ) } ;
  const Ostap::Vector3D v2 { m_rotation ( Ostap::Vector3D ( 0 , m_sigmaY , 0 ) ) } ;
  const Ostap::Vector3D v3 { m_rotation ( Ostap::Vector3D ( 0 , 0 , m_sigmaZ ) ) } ;
  //
  const double sx = std::max ( std::abs ( v1.X() ) , std::max ( std::abs ( v2.X() ) , std::abs ( v3.X() ) ) ) ;
  const double sy = std::max ( std::abs ( v1.Y() ) , std::max ( std::abs ( v2.Y() ) , std::abs ( v3.Y() ) ) ) ;
  const double sz = std::max ( std::abs ( v1.Z() ) , std::max ( std::abs ( v2.Z() ) , std::abs ( v3.Z() ) ) ) ;
  //
  if      ( xhigh <= m_muX - 50 * sx ) { return 0 ; }
  else if ( xlow  >= m_muX + 50 * sx ) { return 0 ; }
  else if ( y     <= m_muY - 50 * sy ) { return 0 ; }
  else if ( y     >= m_muY + 50 * sy ) { return 0 ; }
  else if ( zhigh <= m_muZ - 50 * sz ) { return 0 ; }
  else if ( zlow  >= m_muZ + 50 * sz ) { return 0 ; }
  //  
  // split into smaller regions 
  //
  if ( xlow < m_muX + sx * s_SPLITS.back() || xhigh > m_muX + sx * s_SPLITS.front() ) 
  {
    for ( SPLITS::const_iterator ix = s_SPLITS.begin() ; s_SPLITS.end() != ix  ; ++ix ) 
    {
      const double px = m_muX + (*ix) * sx ;
      if ( xlow < px && px < xhigh ) 
      {
        return 
          integrateXZ ( y , xlow , px    , zlow , zhigh ) +
          integrateXZ ( y , px   , xhigh , zlow , zhigh ) ;
      }
    }
  } 
  // ==========================================================================
  if ( zlow < m_muZ + sz * s_SPLITS.back() || zhigh > m_muZ + sz * s_SPLITS.front() ) 
  {
    for ( SPLITS::const_iterator iz = s_SPLITS.begin() ; s_SPLITS.end() != iz  ; ++iz ) 
    {
      const double pz = m_muZ + (*iz) * sz ;
      if ( zlow < pz && pz < zhigh ) 
      {
        return 
          integrateXZ ( y , xlow , xhigh , zlow , pz    ) +
          integrateXZ ( y , xlow , xhigh , pz   , zhigh ) ;
      }
    }
  }
  // ==========================================================================
  const bool in_tail =  
    ( xhigh <= m_muX + sx * s_SPLITS.front () ) || 
    ( xlow  >= m_muX + sx * s_SPLITS.back  () ) || 
    ( y     <= m_muY + sy * s_SPLITS.front () ) || 
    ( y     >= m_muY + sy * s_SPLITS.back  () ) ||
    ( zhigh <= m_muZ + sz * s_SPLITS.front () ) || 
    ( zlow  >= m_muZ + sz * s_SPLITS.back  () ) ;
  //
  // use 2D-cubature   
  //
  typedef Ostap::Math::IntegrateXZ<Gauss3D> FXZ ;
  const FXZ fxz ( this , y ) ;
  //
  static const Ostap::Math::GSL::Integrator2D<FXZ> s_cubature{} ;
  static const char s_message[] = "IntegralXZ(Gauss3D)" ;
  const auto F = s_cubature.make_function ( &fxz , 
                                            xlow , xhigh , 
                                            zlow , zhigh ) ;
  //
  int    ierror  =  0 ;
  double  result =  1 ;
  double  error  = -1 ;
  std::tie ( ierror , result , error ) = s_cubature.cubature
    ( Ostap::Utils::hash_combiner ( tag () , 'Y' , y ) , &F , 50000 ,
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision 
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      s_message , __FILE__ , __LINE__ ) ;
  return  result ;
}
// ============================================================================
/* integral over x-dimension
 *  \f[ \int_{x_{low}}^{x_{high}} \mathcal{B}(x,y,z) \mathrm{d}x\f]
 *  @param x     variable
 *  @param z     variable
 *  @param ylow  low  edge in y
 *  @param yhigh high edge in y
 */
// ============================================================================
double Ostap::Math::Gauss3D::integrateX
( const double y    ,
  const double z    ,                          
  const double xlow , const double xhigh ) const 
{
  if ( s_equal ( xlow , xhigh ) ) { return 0 ; }
  //
  if ( xhigh < xlow ) { return -integrateX ( y , z , xhigh , xlow  ) ; }
  //
  if ( ( s_zero  ( m_phi   ) && 
         s_zero  ( m_theta ) && 
         s_zero  ( m_psi   ) ) || 
       ( s_equal ( m_sigmaX , m_sigmaY ) && 
         s_equal ( m_sigmaY , m_sigmaZ ) ) ) 
  {
    return 
      ( Ostap::Math::gauss_cdf ( xhigh , m_muX , m_sigmaX ) - 
        Ostap::Math::gauss_cdf ( xlow  , m_muX , m_sigmaX ) ) 
      * Ostap::Math::gauss_pdf ( y     , m_muY , m_sigmaY )                              
      * Ostap::Math::gauss_pdf ( z     , m_muZ , m_sigmaZ ) ;                             
  }
  //
  const Ostap::Vector3D v1 { m_rotation ( Ostap::Vector3D ( m_sigmaX , 0 , 0 ) ) } ;
  const Ostap::Vector3D v2 { m_rotation ( Ostap::Vector3D ( 0 , m_sigmaY , 0 ) ) } ;
  const Ostap::Vector3D v3 { m_rotation ( Ostap::Vector3D ( 0 , 0 , m_sigmaZ ) ) } ;
  //
  const double sx = std::max ( std::abs ( v1.X() ) , std::max ( std::abs ( v2.X() ) , std::abs ( v3.X() ) ) ) ;
  const double sy = std::max ( std::abs ( v1.Y() ) , std::max ( std::abs ( v2.Y() ) , std::abs ( v3.Y() ) ) ) ;
  const double sz = std::max ( std::abs ( v1.Z() ) , std::max ( std::abs ( v2.Z() ) , std::abs ( v3.Z() ) ) ) ;
  //
  if      ( xhigh <= m_muX - 50 * sx ) { return 0 ; }
  else if ( xlow  >= m_muX + 50 * sx ) { return 0 ; }
  else if ( y     <= m_muY - 50 * sy ) { return 0 ; }
  else if ( y     >= m_muY + 50 * sy ) { return 0 ; }
  else if ( z     <= m_muZ - 50 * sz ) { return 0 ; }
  else if ( z     >= m_muZ + 50 * sz ) { return 0 ; }
  //  
  // split into smaller regions 
  //
  if ( xlow < m_muX + sx * s_SPLITS.back() || xhigh > m_muX + sx * s_SPLITS.front() ) 
  {
    for ( SPLITS::const_iterator ix = s_SPLITS.begin() ; s_SPLITS.end() != ix  ; ++ix ) 
    {
      const double px = m_muX + (*ix) * sx ;
      if ( xlow < px && px < xhigh ) 
      {
        return 
          integrateX ( y , z , xlow , px    ) +
          integrateX ( y , z , px   , xhigh ) ;
      }
    } 
  }
  // ==========================================================================
  const bool in_tail = 
    ( xhigh <= m_muX + sx * s_SPLITS.front () ) || 
    ( xlow  >= m_muX + sx * s_SPLITS.back  () ) || 
    ( y     <= m_muY + sy * s_SPLITS.front () ) || 
    ( y     >= m_muY + sy * s_SPLITS.back  () ) ||
    ( z     <= m_muZ + sz * s_SPLITS.front () ) || 
    ( z     >= m_muZ + sz * s_SPLITS.back  () ) ;
  //
  // use 1D-integration    
  //
  typedef Ostap::Math::IntegrateX3<Gauss3D> FX ;
  const FX fx ( this , y , z ) ;
  //
  static const Ostap::Math::GSL::Integrator1D<FX> s_integrator ;
  static const char message[] = "IntegrateX(Gauss3D)" ;
  const auto F = s_integrator.make_function ( &fx ) ;
  //  
  int    ierror    =  0   ;
  double result    =  1.0 ;
  double error     = -1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
    ( Ostap::Utils::hash_combiner ( tag() , 'Y' , 'Z' , y , z  ) , 
      &F                        ,   // the function
      xlow    , xhigh           ,   // low & high edges
      workspace ( m_workspace ) ,   // workspace
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision 
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      m_workspace.size()        ,   // maximum number of subintervals
      message                   ,   // message 
      __FILE__  , __LINE__      ) ; // filename & line number 
  //
  return result ;
}
// ============================================================================
/*  integral over y-dimension
 *  \f[ \int_{y_{low}}^{y_{high}} \mathcal{B}(x,y,z) \mathrm{d}y\f]
 *  @param y     variable
 *  @param z     variable
 *  @param xlow  low  edge in x
 *  @param xhigh high edge in x
 */
// ============================================================================
double Ostap::Math::Gauss3D::integrateY 
( const double x    ,
  const double z    ,
  const double ylow , const double yhigh ) const 
{
  if ( s_equal ( ylow , yhigh ) ) { return 0 ; }
  //
  if ( yhigh < ylow ) { return -integrateY ( x , z , yhigh , ylow  ) ; }
  //
  if ( ( s_zero  ( m_phi   ) && 
         s_zero  ( m_theta ) && 
         s_zero  ( m_psi   ) ) || 
       ( s_equal ( m_sigmaX , m_sigmaY ) && 
         s_equal ( m_sigmaY , m_sigmaZ ) ) ) 
  {
    return 
      ( Ostap::Math::gauss_cdf ( yhigh , m_muY , m_sigmaY ) - 
        Ostap::Math::gauss_cdf ( ylow  , m_muY , m_sigmaY ) ) 
      * Ostap::Math::gauss_pdf ( x     , m_muX , m_sigmaX )                              
      * Ostap::Math::gauss_pdf ( z     , m_muZ , m_sigmaZ ) ;                             
  }
  //
  const Ostap::Vector3D v1 { m_rotation ( Ostap::Vector3D ( m_sigmaX , 0 , 0 ) ) } ;
  const Ostap::Vector3D v2 { m_rotation ( Ostap::Vector3D ( 0 , m_sigmaY , 0 ) ) } ;
  const Ostap::Vector3D v3 { m_rotation ( Ostap::Vector3D ( 0 , 0 , m_sigmaZ ) ) } ;
  //
  const double sx = std::max ( std::abs ( v1.X() ) , std::max ( std::abs ( v2.X() ) , std::abs ( v3.X() ) ) ) ;
  const double sy = std::max ( std::abs ( v1.Y() ) , std::max ( std::abs ( v2.Y() ) , std::abs ( v3.Y() ) ) ) ;
  const double sz = std::max ( std::abs ( v1.Z() ) , std::max ( std::abs ( v2.Z() ) , std::abs ( v3.Z() ) ) ) ;
  //
  if      ( x     <= m_muX - 50 * sx ) { return 0 ; }
  else if ( x     >= m_muX + 50 * sx ) { return 0 ; }
  else if ( yhigh <= m_muY - 50 * sy ) { return 0 ; }
  else if ( ylow  >= m_muY + 50 * sy ) { return 0 ; }
  else if ( z     <= m_muZ - 50 * sz ) { return 0 ; }
  else if ( z     >= m_muZ + 50 * sz ) { return 0 ; }
  //  
  // split into smaller regions 
  //
  if ( ylow < m_muY + sy * s_SPLITS.back() || yhigh > m_muY + sy * s_SPLITS.front() ) 
  {
    for ( SPLITS::const_iterator iy = s_SPLITS.begin() ; s_SPLITS.end() != iy  ; ++iy ) 
    {
      const double py = m_muY + (*iy) * sy ;
      if ( ylow < py && py < yhigh ) 
      {
        return 
          integrateY ( x , z , ylow , py    ) +
          integrateY ( x , z , py   , yhigh ) ;
      }
    }
  } 
  // ==========================================================================
  const bool in_tail = 
    ( x     <= m_muX + sx * s_SPLITS.front () ) || 
    ( x     >= m_muX + sx * s_SPLITS.back  () ) || 
    ( yhigh <= m_muY + sy * s_SPLITS.front () ) || 
    ( ylow  >= m_muY + sy * s_SPLITS.back  () ) ||
    ( z     <= m_muZ + sz * s_SPLITS.front () ) || 
    ( z     >= m_muZ + sz * s_SPLITS.back  () ) ;
  //
  // Use 1D-integration 
  //
  typedef Ostap::Math::IntegrateY3<Gauss3D> FY ;
  const FY fy ( this , x , z ) ;
  //
  static const Ostap::Math::GSL::Integrator1D<FY> s_integrator ;
  static const char message[] = "IntegrateY(Gauss3D)" ;
  const auto F = s_integrator.make_function ( &fy ) ;
  //  
  int    ierror    =  0   ;
  double result    =  1.0 ;
  double error     = -1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
    ( Ostap::Utils::hash_combiner ( tag() , 'X' , 'Z' , x , z ) , 
      &F                        ,   // the function
      ylow    , yhigh           ,   // low & high edges
      workspace ( m_workspace ) ,   // workspace
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision 
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      m_workspace.size()        ,   // maximum number of subintervals
      message                   ,   // message 
      __FILE__  , __LINE__      ) ; // filename & line number 
  //
  return result ;
}
// ============================================================================
/*  integral over z-dimension
 *  \f[ \int_{z_{low}}^{z_{high}} \mathcal{B}(x,y,z) \mathrm{d}z\f]
 *  @param x     variable
 *  @param y     variable
 *  @param zlow  low  edge in z
 *  @param zhigh high edge in z
 */
// ============================================================================
double Ostap::Math::Gauss3D::integrateZ 
( const double x    ,
  const double y    ,
  const double zlow , const double zhigh ) const
{
  if ( s_equal ( zlow , zhigh ) ) { return 0 ; }
  //
  if ( zhigh < zlow ) { return -integrateZ ( x , y , zhigh , zlow  ) ; }
  //
  if ( ( s_zero  ( m_phi   ) && 
         s_zero  ( m_theta ) && 
         s_zero  ( m_psi   ) ) || 
       ( s_equal ( m_sigmaX , m_sigmaY ) && 
         s_equal ( m_sigmaY , m_sigmaZ ) ) ) 
  {
    return 
      ( Ostap::Math::gauss_cdf ( zhigh , m_muZ , m_sigmaZ ) - 
        Ostap::Math::gauss_cdf ( zlow  , m_muZ , m_sigmaZ ) ) 
      * Ostap::Math::gauss_pdf ( x     , m_muX , m_sigmaX )                              
      * Ostap::Math::gauss_pdf ( y     , m_muY , m_sigmaY ) ;                             
  }
  //
  const Ostap::Vector3D v1 { m_rotation ( Ostap::Vector3D ( m_sigmaX , 0 , 0 ) ) } ;
  const Ostap::Vector3D v2 { m_rotation ( Ostap::Vector3D ( 0 , m_sigmaY , 0 ) ) } ;
  const Ostap::Vector3D v3 { m_rotation ( Ostap::Vector3D ( 0 , 0 , m_sigmaZ ) ) } ;
  //
  const double sx = std::max ( std::abs ( v1.X() ) , std::max ( std::abs ( v2.X() ) , std::abs ( v3.X() ) ) ) ;
  const double sy = std::max ( std::abs ( v1.Y() ) , std::max ( std::abs ( v2.Y() ) , std::abs ( v3.Y() ) ) ) ;
  const double sz = std::max ( std::abs ( v1.Z() ) , std::max ( std::abs ( v2.Z() ) , std::abs ( v3.Z() ) ) ) ;
  //
  if      ( x     <= m_muX - 50 * sx ) { return 0 ; }
  else if ( x     >= m_muX + 50 * sx ) { return 0 ; }
  else if ( y     <= m_muY - 50 * sy ) { return 0 ; }
  else if ( y     >= m_muY + 50 * sy ) { return 0 ; }
  else if ( zhigh <= m_muZ - 50 * sz ) { return 0 ; }
  else if ( zlow  >= m_muZ + 50 * sz ) { return 0 ; }
  //  
  // split into smaller regions 
  //
  if ( zlow < m_muZ + sz * s_SPLITS.back() || zhigh > m_muZ + sz * s_SPLITS.front() ) 
  {
    for ( SPLITS::const_iterator iz = s_SPLITS.begin() ; s_SPLITS.end() != iz  ; ++iz ) 
    {
      const double pz = m_muZ + (*iz) * sz ;
      if ( zlow < pz && pz < zhigh ) 
      {
        return 
          integrateZ ( x , y , zlow , pz    ) +
          integrateZ ( x , y , pz   , zhigh ) ;
      }
    } 
  }
  // ==========================================================================
  const bool in_tail = 
    ( x     <= m_muX + sx * s_SPLITS.front () ) || 
    ( x     >= m_muX + sx * s_SPLITS.back  () ) || 
    ( y     <= m_muY + sy * s_SPLITS.front () ) || 
    ( y     >= m_muY + sy * s_SPLITS.back  () ) ||
    ( zhigh <= m_muZ + sz * s_SPLITS.front () ) || 
    ( zlow  >= m_muZ + sz * s_SPLITS.back  () ) ;
  //
  // use 1D-integration    
  //
  typedef Ostap::Math::IntegrateZ3<Gauss3D> FZ ;
  const FZ fz ( this , x , y ) ;
  //
  static const Ostap::Math::GSL::Integrator1D<FZ> s_integrator ;
  static const char message[] = "IntegrateZ(Gauss3D)" ;
  const auto F = s_integrator.make_function ( &fz ) ;
  //  
  int    ierror    =  0   ;
  double result    =  1.0 ;
  double error     = -1.0 ;
  std::tie ( ierror , result , error ) = s_integrator.gaq_integrate
    ( Ostap::Utils::hash_combiner ( tag() , 'X' , 'Y' , x , y ) , 
      &F                        ,   // the function
      zlow    , zhigh           ,   // low & high edges
      workspace ( m_workspace ) ,   // workspace
      in_tail ? s_APRECISION_TAIL : s_APRECISION , // absolute precision 
      in_tail ? s_RPRECISION_TAIL : s_RPRECISION , // relative precision
      m_workspace.size()        ,   // maximum number of subintervals
      message                   ,   // message 
      __FILE__  , __LINE__      ) ; // filename & line number 
  //
  return result ;
}
// ============================================================================
// get the tag 
// ============================================================================
std::size_t Ostap::Math::Gauss3D::tag () const 
{ return Ostap::Utils::hash_combiner ( m_muX    , 
                                       m_muY    , 
                                       m_muZ    , 
                                       m_sigmaX , 
                                       m_sigmaY , 
                                       m_sigmaZ , 
                                       m_phi    ,
                                       m_theta  ,
                                       m_psi    ) ; }

// ============================================================================
//                                                                      The END 
// ============================================================================
