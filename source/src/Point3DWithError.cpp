// ============================================================================
// Include files
// ============================================================================
// STD& STL
// ============================================================================
#include <sstream>
// ============================================================================
// local
// ============================================================================
#include "Ostap/Point3DWithError.h"
#include "Ostap/Vector3DWithError.h"
#include "Ostap/MatrixTransforms.h"
// ============================================================================
/** @file
 *  Implementation file for class Ostap::Math::PointWithError
 *  @author Vanya BELYAEV Ivane.BElyaev@nikhef.nl
 *  @date 20090603
 */
// ============================================================================
// constructor form the point and covariance matrix
// ============================================================================
Ostap::Math::Point3DWithError::Point3DWithError
( const Ostap::Math::Point3DWithError::Point3D&    point  ,
  const Ostap::Math::Point3DWithError::Covariance& matrix )
  : Ostap::XYZPoint ( point )
  , m_cov2  ( matrix )
{}
// ============================================================================
// constructor form the point and covariance matrix
// ============================================================================
Ostap::Math::Point3DWithError::Point3DWithError
( const Ostap::Math::Point3DWithError::Covariance& matrix ,
  const Ostap::Math::Point3DWithError::Point3D&    point  )
  : Ostap::XYZPoint ( point )
  , m_cov2  ( matrix )
{}
// ============================================================================
// constructor form the point and covariance matrix
// ============================================================================
Ostap::Math::Point3DWithError::Point3DWithError
( const Ostap::Math::Point3DWithError::Vector&     point  ,
  const Ostap::Math::Point3DWithError::Covariance& matrix )
  : Ostap::XYZPoint ( point[0] , point[1] , point[2] )
  , m_cov2  ( matrix )
{}
// ============================================================================
// constructor form the point and covariance matrix
// ============================================================================
Ostap::Math::Point3DWithError::Point3DWithError
( const Ostap::Math::Point3DWithError::VectorE&    point  )
  : Ostap::XYZPoint ( point.value(0) , point.value(1) , point.value(2) )
  , m_cov2  ( point.cov2() )
{}
// ============================================================================
Ostap::Math::Point3DWithError&
Ostap::Math::Point3DWithError::operator+=
( const Ostap::Math::Vector3DWithError& right )
{
  point () += right.value() ;
  m_cov2   += right.cov2 () ;
  return *this ;
}
// ============================================================================
Ostap::Math::Point3DWithError&
Ostap::Math::Point3DWithError::operator-=
( const Ostap::Math::Vector3DWithError& right )
{
  point() -= right.value() ;
  m_cov2  += right.cov2 () ;
  return *this ;
}
// ============================================================================
Ostap::Math::Point3DWithError&
Ostap::Math::Point3DWithError::operator+=
( const Ostap::Math::Point3DWithError::Vector3D& right )
{
  point() += right ;
  return *this ;
}
// ============================================================================
Ostap::Math::Point3DWithError&
Ostap::Math::Point3DWithError::operator-=
( const Ostap::Math::Point3DWithError::Vector3D& right )
{
  point() -= right;
  return *this ;
}
// ============================================================================
Ostap::Math::Point3DWithError&
Ostap::Math::Point3DWithError::operator+=
( const Ostap::Math::Point3DWithError::VectorE& right )
{
  using namespace Ostap::Math::Operators ;
  setPoint ( point() + Vector3D(right.value()[0],
                                right.value()[1],
                                right.value()[2]) ) ;
  m_cov2   += right.cov2() ;
  return *this ;
}
// ============================================================================
Ostap::Math::Point3DWithError&
Ostap::Math::Point3DWithError::operator-=
( const Ostap::Math::Point3DWithError::VectorE& right )
{
  using namespace Ostap::Math::Operators ;
  setPoint ( point() - Vector3D(right.value()[0],
                                right.value()[1],
                                right.value()[2]) ) ;
  m_cov2   += right.cov2() ;
  return *this ;
}
// ============================================================================
Ostap::Math::Point3DWithError&
Ostap::Math::Point3DWithError::operator+=
( const Ostap::Math::Point3DWithError::Vector& right )
{
  using namespace Ostap::Math::Operators ;
  setPoint ( Point3D( point() + Vector3D(right[0],right[1],right[2]) ) ) ;
  return *this ;
}
// ============================================================================
Ostap::Math::Point3DWithError&
Ostap::Math::Point3DWithError::operator-=
( const Ostap::Math::Point3DWithError::Vector& right )
{
  using namespace Ostap::Math::Operators ;
  setPoint ( Point3D( point() - Vector3D(right[0],right[1],right[2]) ) ) ;
  return *this ;
}
// ============================================================================
// scaling
// ============================================================================
Ostap::Math::Point3DWithError&
Ostap::Math::Point3DWithError::operator*= ( const double v )
{
  point () *=  v    ;
  m_cov2   *= (v*v) ;
  return *this ;
}
// ============================================================================
Ostap::Math::Point3DWithError&
Ostap::Math::Point3DWithError::operator/= ( const double v )
{
  point () /=  v    ;
  m_cov2   /= (v*v) ;
  return *this ;
}
// ============================================================================
// chi2 distance
// ============================================================================
double Ostap::Math::Point3DWithError::chi2
( const Ostap::Math::Point3DWithError& right ) const
{
  Covariance s_cov2 ( cov2() ) ;
  s_cov2 += right.cov2() ;
  // use Manuel's inverter:
  const bool ok = s_cov2.InvertChol() ;
  if  ( !ok ) { return -1 ; }                                       // RETURN
  /// calculate chi2
  return Ostap::Math::Similarity ( point() - right.point() , s_cov2 ) ;
}
// ============================================================================
// chi2 distance
// ============================================================================
double Ostap::Math::Point3DWithError::chi2
( const Ostap::XYZPoint& right ) const
{
  Covariance s_cov2 ( cov2() ) ;
  // use Manuel's inverter:
  const bool ok = s_cov2.InvertChol() ;
  if  ( !ok ) { return -1 ; }                                        // RETURN
  /// calculate chi2
  return Ostap::Math::Similarity ( point() - right , s_cov2 ) ;
}
// ============================================================================
// chi2 distance
// ============================================================================
double Ostap::Math::Point3DWithError::chi2
( const Ostap::Math::Point3DWithError::Vector& right ) const
{
  Covariance s_cov2 ( cov2() ) ;
  // use Manuel's inverter:
  const bool ok = s_cov2.InvertChol() ;
  if  ( !ok ) { return -1 ; }                                        // RETURN
  /// calculate chi2
  //
  Vector vct   ;
  Ostap::Math::geo2LA ( point() , vct ) ;
  vct -= right ;
  //
  return ROOT::Math::Similarity ( vct , s_cov2 ) ;
}
// ============================================================================
// chi2 distance
// ============================================================================
double Ostap::Math::Point3DWithError::chi2
( const Ostap::Math::Point3DWithError::VectorE& right ) const
{
  Covariance s_cov2 ( cov2() ) ;
  s_cov2 += right.cov2() ;
  // use Manuel's inverter:
  const bool ok = s_cov2.InvertChol() ;
  if  ( !ok ) { return -1 ; }                                       // RETURN
  /// calculate chi2
  //
  Vector vct ;
  Ostap::Math::geo2LA ( point() , vct ) ;
  vct -= right.value() ;
  //
  return ROOT::Math::Similarity ( vct , s_cov2 ) ;
}
// ============================================================================
// printout
// ============================================================================
namespace
{
  inline double err ( double cov )
  { return 0 <= cov ? std::sqrt ( cov ) : -std::sqrt(-cov) ; }
}
// ============================================================================
std::ostream&
Ostap::Math::Point3DWithError::fillStream ( std::ostream& s ) const // printout
{
  return s << "( "
           << X () << " +- " << err ( m_cov2(0,0) ) << " , "
           << Y () << " +- " << err ( m_cov2(1,1) ) << " , "
           << Z () << " +- " << err ( m_cov2(2,2) ) << " )";
}
// ============================================================================
// conversion to the string
// ============================================================================
std::string
Ostap::Math::Point3DWithError::toString   () const // conversion to the string
{
  std::ostringstream s ;
  fillStream ( s ) ;
  return s.str() ;
}
// ============================================================================


// ============================================================================
Ostap::Math::Point3DWithError
Ostap::Math::Point3DWithError::__add__
( const Ostap::Math::Vector3DWithError& right ) const
{
  Ostap::Math::Point3DWithError tmp ( *this ) ;
  return tmp += right ;
}
// ============================================================================
Ostap::Math::Point3DWithError
Ostap::Math::Point3DWithError::__sub__
( const Ostap::Math::Vector3DWithError& right ) const
{
  Ostap::Math::Point3DWithError tmp ( *this ) ;
  return tmp -= right ;
}
// ============================================================================
Ostap::Math::Point3DWithError
Ostap::Math::Point3DWithError::__add__ ( const Ostap::XYZVector& right ) const
{
  Ostap::Math::Point3DWithError tmp ( *this ) ;
  return tmp += right ;
}
// ============================================================================
Ostap::Math::Point3DWithError
Ostap::Math::Point3DWithError::__sub__ ( const Ostap::XYZVector& right ) const
{
  Ostap::Math::Point3DWithError tmp ( *this ) ;
  return tmp -= right ;
}


// ============================================================================
Ostap::Math::Vector3DWithError
Ostap::Math::Point3DWithError::__sub__
( const Ostap::Math::Point3DWithError& right ) const
{ return Ostap::Math::Vector3DWithError ( point() -  right.point () ,
                                          cov2 () +  right.cov2  () ) ; }
// ============================================================================
Ostap::Math::Vector3DWithError
Ostap::Math::Point3DWithError::__sub__
( const Ostap::XYZPoint& right ) const
{ return Ostap::Math::Vector3DWithError ( point() -  right , cov2() ) ; }
// ============================================================================


// ============================================================================
Ostap::Math::Vector3DWithError
Ostap::Math::Point3DWithError::__rsub__
( const Ostap::XYZPoint& right ) const
{ return Ostap::Math::Vector3DWithError ( right - point() , cov2() ) ; }
// ============================================================================


// ============================================================================
void Ostap::Math::Point3DWithError::asVector
( Ostap::Math::Point3DWithError::Vector& data ) const
{ Ostap::Math::geo2LA ( point() , data ) ; }
// ============================================================================
void Ostap::Math::Point3DWithError::asVector
( Ostap::Math::Point3DWithError::VectorE& data ) const
{
  Ostap::Math::geo2LA ( point() , data.value() ) ;
  data.setCov2 ( cov2() ) ;
}
// ============================================================================
Ostap::Math::Point3DWithError::VectorE
Ostap::Math::Point3DWithError::asVector () const
{
  Ostap::Math::Point3DWithError::VectorE data ;

  Ostap::Math::geo2LA ( point() , data.value() ) ;
  data.setCov2 ( cov2() ) ;

  return data ;
}
// ============================================================================
Ostap::Math::Point3DWithError::Vector
Ostap::Math::Point3DWithError::asVector3 () const
{
  Ostap::Math::Point3DWithError::Vector data ;
  Ostap::Math::geo2LA ( point() , data ) ;
  return data ;
}
// ============================================================================




// ============================================================================
void Ostap::Math::Point3DWithError::setValue
( const Ostap::Math::Point3DWithError::VectorE& v )
{
  Ostap::Math::la2geo ( v.value() , point() ) ;
  m_cov2 = v.cov2() ;
}
// ============================================================================
void Ostap::Math::Point3DWithError::setValue
( const Ostap::Math::Point3DWithError::Vector& v )
{
  Ostap::Math::la2geo ( v , point() ) ;
}
// ============================================================================


// ============================================================================
Ostap::Math::Point3DWithError&
Ostap::Math::Point3DWithError::__imul__ ( const double v )
{
  point() *= v ;
  m_cov2  *= (v*v) ;
  return *this ;
}
// ============================================================================
Ostap::Math::Point3DWithError&
Ostap::Math::Point3DWithError::__itruediv__ ( const double v )
{
  point() /= v ;
  m_cov2  /= (v*v) ;
  return *this ;
}
// ============================================================================

// ============================================================================
Ostap::Math::Point3DWithError
Ostap::Math::Point3DWithError::__mul__ ( const double v ) const
{
  Ostap::Math::Point3DWithError tmp (*this) ;
  return ( tmp *= v ) ;
}
// ============================================================================

// ============================================================================
Ostap::Math::Point3DWithError
Ostap::Math::Point3DWithError::__truediv__ ( const double v ) const
{
  Ostap::Math::Point3DWithError tmp (*this) ;
  return ( tmp /= v ) ;
}
// ============================================================================
Ostap::Math::Point3DWithError 
Ostap::Math::Point3DWithError::mean 
( const Ostap::Math::Point3DWithError&          right ) const 
{ return asVector ().mean ( right.asVector() ) ; }
// ============================================================================
Ostap::Math::Point3DWithError 
Ostap::Math::Point3DWithError::mean 
( const Ostap::Math::Point3DWithError::VectorE& right ) const 
{ return asVector ().mean ( right            ) ; }
// ============================================================================
/* Get symmetrized Kullback-Leibler divergency for two objects 
 *  @see https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence
 *  @see Ostap::Math::kullback_leibler 
 */
// ============================================================================
double Ostap::Math::kullback_leibler 
( const Ostap::Math::Point3DWithError& a , 
  const Ostap::Math::Point3DWithError& b ) 
{
  return Ostap::Math::kullback_leibler 
    ( a.asVector3() , a.covariance () , 
      b.asVector3() , b.covariance () ) ;
}
// ============================================================================
/* Get asymmetric Kullback-Leibler divergency for two objects 
 *  @see https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence
 *  @see Ostap::Math::asymmetric_kullback_leibler 
 */
// ============================================================================
double Ostap::Math::asymmetric_kullback_leibler 
( const Ostap::Math::Point3DWithError& a , 
  const Ostap::Math::Point3DWithError& b ) 
{
  return Ostap::Math::asymmetric_kullback_leibler 
    ( a.asVector3() , a.covariance () , 
      b.asVector3() , b.covariance () ) ;
}
// ============================================================================

// ============================================================================
//                                                                      The END
// ============================================================================

