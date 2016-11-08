// ============================================================================
// Include files 
// ============================================================================
// STD& STL 
// ============================================================================
#include <sstream>
// ============================================================================
// local
// ============================================================================
#include "Ostap/Vector3DWithError.h"
#include "Ostap/Point3DWithError.h"
#include "Ostap/MatrixTransforms.h"
#include "Ostap/Kinematics.h"
// ============================================================================
/** @file 
 *  Implementation file for class Ostap::Math::PointWithError
 *  @author Vanya BELYAEV Ivane.BElyaev@nikhef.nl
 *  @date 20090603
 */
// ============================================================================
// constructor from the point and covariance matrix
// ============================================================================
Ostap::Math::Vector3DWithError::Vector3DWithError 
( const Ostap::Math::Vector3DWithError::Vector3D&   vct    ,
  const Ostap::Math::Vector3DWithError::Covariance& matrix ) 
  : Ostap::XYZVector ( vct )  
  , m_cov2  ( matrix ) 
{}
// ============================================================================
// constructor from the point and covariance matrix
// ============================================================================
Ostap::Math::Vector3DWithError::Vector3DWithError 
( const Ostap::Math::Vector3DWithError::Covariance& matrix , 
  const Ostap::Math::Vector3DWithError::Vector3D&   vct    )
  : Ostap::XYZVector ( vct )  
  , m_cov2  ( matrix ) 
{}
// ============================================================================
// constructor from generic vector and covariance matrix
// ============================================================================
Ostap::Math::Vector3DWithError::Vector3DWithError 
( const Ostap::Math::Vector3DWithError::Vector&     vct    ,
  const Ostap::Math::Vector3DWithError::Covariance& matrix ) 
  : Ostap::XYZVector ( vct[0] , vct[1] , vct[2] )  
  , m_cov2  ( matrix ) 
{}
// ============================================================================
// constructor from generic vector and covariance matrix
// ============================================================================
Ostap::Math::Vector3DWithError::Vector3DWithError 
( const Ostap::Math::Vector3DWithError::VectorE& vct    ) 
  : Ostap::XYZVector ( vct.value(0) , vct.value(1) , vct.value(2) )  
  , m_cov2  ( vct.cov2() ) 
{}
// ============================================================================
Ostap::Math::Vector3DWithError& 
Ostap::Math::Vector3DWithError::operator+= 
( const Ostap::Math::Vector3DWithError& right ) 
{
  vector3d () += right.value() ;
  m_cov2      += right.cov2 () ;
  return *this ;
}
// ============================================================================
Ostap::Math::Vector3DWithError& 
Ostap::Math::Vector3DWithError::operator-= 
( const Ostap::Math::Vector3DWithError& right ) 
{
  vector3d () -= right.value() ;
  m_cov2      += right.cov2 () ;
  return *this ;
}
// ============================================================================
Ostap::Math::Vector3DWithError& 
Ostap::Math::Vector3DWithError::operator+= 
( const Ostap::Math::Vector3DWithError::Vector3D& right ) 
{
  vector3d () += right ;
  return *this ;
}
// ============================================================================
Ostap::Math::Vector3DWithError& 
Ostap::Math::Vector3DWithError::operator-= 
( const Ostap::Math::Vector3DWithError::Vector3D& right ) 
{
  vector3d() -= right;
  return *this ;
}
// ============================================================================
Ostap::Math::Vector3DWithError& 
Ostap::Math::Vector3DWithError::operator+= 
( const Ostap::Math::Vector3DWithError::VectorE& right ) 
{
  using namespace Ostap::Math::Operators ;
  setVector ( vector3d() + Vector3D(right.value()[0],
                                    right.value()[1],
                                    right.value()[2]) )  ;
  m_cov2 += right.cov2() ;
  return *this ;
}
// ============================================================================
Ostap::Math::Vector3DWithError& 
Ostap::Math::Vector3DWithError::operator-= 
( const Ostap::Math::Vector3DWithError::VectorE& right ) 
{
  using namespace Ostap::Math::Operators ;
  setVector ( vector3d() - Vector3D(right.value()[0],
                                    right.value()[1],
                                    right.value()[2]) )  ;
  m_cov2 += right.cov2() ;
  return *this ;
}
// ============================================================================
Ostap::Math::Vector3DWithError& 
Ostap::Math::Vector3DWithError::operator+= 
( const Ostap::Math::Vector3DWithError::Vector& right ) 
{
  using namespace Ostap::Math::Operators ;
  setVector ( vector3d() + Vector3D(right[0],right[1],right[2]) )  ;
  return *this ;
}
// ============================================================================
Ostap::Math::Vector3DWithError& 
Ostap::Math::Vector3DWithError::operator-= 
( const Ostap::Math::Vector3DWithError::Vector& right ) 
{
  using namespace Ostap::Math::Operators ;
  setVector ( vector3d() - Vector3D(right[0],right[1],right[2]) )  ;
  return *this ;
}
// ============================================================================
// scaling
// ============================================================================
Ostap::Math::Vector3DWithError& 
Ostap::Math::Vector3DWithError::operator*= ( const double v ) 
{
  vector3d () *=  v    ;
  m_cov2      *= (v*v) ;
  return *this ;
}
// ============================================================================
Ostap::Math::Vector3DWithError& 
Ostap::Math::Vector3DWithError::operator/= ( const double v ) 
{
  vector3d () /=  v    ;
  m_cov2      /= (v*v) ;
  return *this ;
}

// ============================================================================
// chi2 distance 
// ============================================================================
double Ostap::Math::Vector3DWithError::chi2 
( const Ostap::Math::Vector3DWithError& right ) const 
{
  Covariance s_cov2 ( cov2() ) ;
  s_cov2 += right.cov2() ;
  // use Manuel's inverter:
  const bool ok = s_cov2.InvertChol() ;
  if  ( !ok ) { return -1 ; }                                 // RETURN  
  /// calculate chi2 
  return Ostap::Math::Similarity ( vector3d() - right.vector3d() , s_cov2 ) ;
}

// ============================================================================
// chi2 distance 
// ============================================================================
double Ostap::Math::Vector3DWithError::chi2
( const Ostap::XYZVector& right ) const 
{
  Covariance s_cov2 ( cov2() ) ;
  // use Manuel's inverter:
  const bool ok = s_cov2.InvertChol() ;
  if  ( !ok ) { return -1 ; }                                 // RETURN  
  /// calculate chi2 
  return Ostap::Math::Similarity ( vector3d() - right , s_cov2 ) ;
}

// ============================================================================
// chi2 distance 
// ============================================================================
double Ostap::Math::Vector3DWithError::chi2 
( const Ostap::Math::Vector3DWithError::VectorE& right ) const 
{
  Covariance s_cov2 ( cov2() ) ;
  s_cov2 += right.cov2() ;
  // use Manuel's inverter:
  const bool ok = s_cov2.InvertChol() ;
  if  ( !ok ) { return -1 ; }                                 // RETURN  
  //
  /// calculate chi2 
  Vector vct ;
  Ostap::Math::geo2LA ( vector3d() , vct ) ;
  vct -= right.value() ;
  //
  return ROOT::Math::Similarity ( vct , s_cov2 ) ;
}

// ============================================================================
// chi2 distance 
// ============================================================================
double Ostap::Math::Vector3DWithError::chi2
( const Ostap::Math::Vector3DWithError::Vector& right ) const 
{
  Covariance s_cov2 ( cov2() ) ;
  // use Manuel's inverter:
  const bool ok = s_cov2.InvertChol() ;
  if  ( !ok ) { return -1 ; }                                 // RETURN  
  /// calculate chi2 
  //
  Vector vct   ;
  Ostap::Math::geo2LA ( vector3d() , vct ) ;
  vct -= right ;
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
Ostap::Math::Vector3DWithError::fillStream ( std::ostream& s ) const // printout 
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
Ostap::Math::Vector3DWithError::toString   () const // conversion to the string 
{
  std::ostringstream s ;
  fillStream ( s ) ;
  return s.str() ;
}

// ============================================================================
// unary- 
// ============================================================================
Ostap::Math::Vector3DWithError 
Ostap::Math::Vector3DWithError::operator-() const 
{ return Ostap::Math::Vector3DWithError ( -value() , cov2() ) ; }



// ============================================================================
// operators 
// ============================================================================
Ostap::Math::Vector3DWithError 
Ostap::Math::Vector3DWithError::__add__ 
( const Ostap::Math::Vector3DWithError& right ) const 
{
  Ostap::Math::Vector3DWithError tmp ( *this ) ;
  return tmp += right ;
}
// ============================================================================
Ostap::Math::Vector3DWithError 
Ostap::Math::Vector3DWithError::__add__ 
( const Ostap::XYZVector& right ) const 
{
  Ostap::Math::Vector3DWithError tmp ( *this ) ;
  return tmp += right ;
}
// ============================================================================
Ostap::Math::Point3DWithError 
Ostap::Math::Vector3DWithError::__add__ 
( const Ostap::Math::Point3DWithError& right ) const 
{
  Ostap::Math::Point3DWithError tmp ( right ) ;
  return tmp += (*this)  ;
}
// ============================================================================
Ostap::Math::Vector3DWithError 
Ostap::Math::Vector3DWithError::__sub__ 
( const Ostap::Math::Vector3DWithError& right ) const 
{
  Ostap::Math::Vector3DWithError tmp ( *this ) ;
  return tmp -= right ;
}
// ============================================================================
Ostap::Math::Vector3DWithError 
Ostap::Math::Vector3DWithError::__sub__ 
( const Ostap::XYZVector& right ) const 
{
  Ostap::Math::Vector3DWithError tmp ( *this ) ;
  return tmp -= right ;
}
// ============================================================================
Ostap::Math::Vector3DWithError 
Ostap::Math::Vector3DWithError::__rsub__ 
( const Ostap::XYZVector& right ) const 
{ return Ostap::Math::Vector3DWithError ( right - value() , cov2() ) ; }
// ============================================================================


// ============================================================================
void Ostap::Math::Vector3DWithError::setValue
( const Ostap::Math::Vector3DWithError::Vector& v ) 
{
  SetX ( v[0] ) ;
  SetY ( v[1] ) ;
  SetZ ( v[2] ) ;
}
// ============================================================================
void Ostap::Math::Vector3DWithError::setValue
( const Ostap::Math::Vector3DWithError::VectorE& v ) 
{
  setValue ( v.value () ) ;
  m_cov2 = v.cov2() ;
}
// ============================================================================

// ============================================================================
void Ostap::Math::Vector3DWithError::asVector 
( Ostap::Math::Vector3DWithError::Vector& data ) const 
{ Ostap::Math::geo2LA ( vector3d() , data ) ; }
// ============================================================================
void Ostap::Math::Vector3DWithError::asVector 
( Ostap::Math::Vector3DWithError::VectorE& data ) const 
{
  Ostap::Math::geo2LA ( vector3d() , data.value() ) ; 
  data.setCov2( cov2() ) ;
}
// ============================================================================
Ostap::Math::Vector3DWithError::VectorE 
Ostap::Math::Vector3DWithError::asVector () const 
{
  Ostap::Math::Vector3DWithError::VectorE data ;
  //
  Ostap::Math::geo2LA ( vector3d() , data.value() ) ; 
  data.setCov2( cov2() ) ;
  //
  return data ;
}
// ============================================================================

// ========================================================================
Ostap::Math::Vector3DWithError operator- 
( const Ostap::XYZPoint&               b ,
  const Ostap::Math::Point3DWithError& a ) { return a.__rsub__ ( b ) ; }
Ostap::Math::Vector3DWithError operator- 
( const Ostap::Math::Point3DWithError& a , 
  const Ostap::Math::Point3DWithError& b ) { return a.__sub__ ( b ) ; }
Ostap::Math::Vector3DWithError operator- 
( const Ostap::Math::Point3DWithError& a , 
  const Ostap::XYZPoint&               b ) { return a.__sub__ ( b ) ; }
// ========================================================================

// ============================================================================
Ostap::Math::Vector3DWithError&
Ostap::Math::Vector3DWithError::__imul__ ( const double v ) 
{
  vector3d() *= v ;
  m_cov2     *= (v*v) ;
  return *this ;
}
// ============================================================================
Ostap::Math::Vector3DWithError&
Ostap::Math::Vector3DWithError::__idiv__ ( const double v ) 
{
  vector3d() /= v ;
  m_cov2     /= (v*v) ;
  return *this ;
}
// ============================================================================

// ============================================================================
Ostap::Math::Vector3DWithError 
Ostap::Math::Vector3DWithError::__mul__ ( const double v ) const
{
  Ostap::Math::Vector3DWithError tmp (*this) ;
  return ( tmp *= v ) ;
}
// ============================================================================

// ============================================================================
Ostap::Math::Vector3DWithError 
Ostap::Math::Vector3DWithError::__div__ ( const double v ) const 
{
  Ostap::Math::Vector3DWithError tmp (*this) ;
  return ( tmp /= v ) ;
}
// ============================================================================
Ostap::Math::Vector3DWithError 
Ostap::Math::Vector3DWithError::mean 
( const Ostap::Math::Vector3DWithError&          right ) const 
{ return asVector ().mean ( right.asVector() ) ; }
// ============================================================================
Ostap::Math::Vector3DWithError 
Ostap::Math::Vector3DWithError::mean 
( const Ostap::Math::Vector3DWithError::VectorE& right ) const 
{ return asVector ().mean ( right            ) ; }
// ============================================================================


// ============================================================================
// The END 
// ============================================================================

