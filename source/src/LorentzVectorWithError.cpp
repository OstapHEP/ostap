// $Id$ 
// ============================================================================
// Include files 
// ============================================================================
// STD & STL 
// ============================================================================
#include <sstream>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/LorentzVectorWithError.h"
#include "Ostap/MatrixTransforms.h"
#include "Ostap/Math.h"
#include "Ostap/Kinematics.h"
// ============================================================================
/** @file
 *  Implementation file for class Ostap::Math::LorentVectorWithError
 *  @date 2009-06-12 
 *  @author Vanya BELYAEV  Ivan.Belyaev@nikhef.nl
 */
// ============================================================================
// constructor for lorent vector and convariance
// ============================================================================
Ostap::Math::LorentzVectorWithError::LorentzVectorWithError 
( const Ostap::LorentzVector&                            value , 
  const Ostap::Math::LorentzVectorWithError::Covariance& cov2  ) 
  : Ostap::LorentzVector ( value ) 
  , m_cov2 ( cov2 ) 
{}
// ============================================================================
// constructor for lorent vector and convariance
// ============================================================================
Ostap::Math::LorentzVectorWithError::LorentzVectorWithError 
( const Ostap::Math::LorentzVectorWithError::Covariance& cov2  , 
  const Ostap::LorentzVector&                            value )
  : Ostap::LorentzVector ( value ) 
    , m_cov2 ( cov2 ) 
{}
// ============================================================================
// constructor for lorent vector and convariance
// ============================================================================
Ostap::Math::LorentzVectorWithError::LorentzVectorWithError 
( const Ostap::Math::LorentzVectorWithError::Vector&     value , 
  const Ostap::Math::LorentzVectorWithError::Covariance& cov2  ) 
  : Ostap::LorentzVector () 
  , m_cov2 ( cov2 ) 
{
  SetPx ( value[0] ) ;
  SetPy ( value[1] ) ;
  SetPz ( value[2] ) ;
  SetE  ( value[3] ) ;
}
// ============================================================================
// constructor for lorent vector and convariance
// ============================================================================
Ostap::Math::LorentzVectorWithError::LorentzVectorWithError 
( const Ostap::Math::LorentzVectorWithError::VectorE& value )
  : Ostap::LorentzVector () 
  , m_cov2 ( value.cov2 () ) 
{
  SetPx ( value[0] ) ;
  SetPy ( value[1] ) ;
  SetPz ( value[2] ) ;
  SetE  ( value[3] ) ;
}
// ============================================================================/
void Ostap::Math::LorentzVectorWithError::setValue 
( const Ostap::Math::LorentzVectorWithError::Vector& v ) 
{
  SetPx ( v [0] ) ;
  SetPy ( v [1] ) ;
  SetPz ( v [2] ) ;
  SetE  ( v [3] ) ;
}
// ============================================================================/
void Ostap::Math::LorentzVectorWithError::setValue 
( const Ostap::Math::LorentzVectorWithError::VectorE& v ) 
{
  setValue ( v.value() ) ;
  m_cov2 = v.cov2() ;
}
// ============================================================================
Ostap::Math::LorentzVectorWithError&
Ostap::Math::LorentzVectorWithError::operator+= 
( const Ostap::Math::LorentzVectorWithError& right ) 
{
  vector4d() += right.vector4d() ;
  m_cov2     += right.m_cov2 ;
  return *this ;
}
// ============================================================================
Ostap::Math::LorentzVectorWithError&
Ostap::Math::LorentzVectorWithError::operator-= 
( const Ostap::Math::LorentzVectorWithError& right ) 
{
  vector4d() -= right.vector4d() ;
  m_cov2     += right.m_cov2 ;
  return *this ;
}
// ============================================================================
Ostap::Math::LorentzVectorWithError&
Ostap::Math::LorentzVectorWithError::operator+= 
( const Ostap::Math::LorentzVectorWithError::VectorE& right ) 
{
  Ostap::Math::add ( vector4d() , right.value () ) ;
  m_cov2     += right.cov2() ;
  return *this ;
}
// ============================================================================
Ostap::Math::LorentzVectorWithError&
Ostap::Math::LorentzVectorWithError::operator-= 
( const Ostap::Math::LorentzVectorWithError::VectorE& right ) 
{
  Ostap::Math::sub ( vector4d() , right.value () ) ;
  m_cov2     += right.cov2() ;
  return *this ;
}
// ============================================================================
Ostap::Math::LorentzVectorWithError&
Ostap::Math::LorentzVectorWithError::operator+= 
( const Ostap::Math::LorentzVectorWithError::Vector& right ) 
{
  Ostap::Math::add ( vector4d() , right ) ;
  return *this ;
}
// ============================================================================
Ostap::Math::LorentzVectorWithError&
Ostap::Math::LorentzVectorWithError::operator-= 
( const Ostap::Math::LorentzVectorWithError::Vector& right ) 
{
  Ostap::Math::sub ( vector4d() , right ) ;
  return *this ;
}
// ============================================================================
Ostap::Math::LorentzVectorWithError&
Ostap::Math::LorentzVectorWithError::operator+= 
( const Ostap::LorentzVector& right ) 
{
  vector4d() += right ;
  return *this ;
}
// ============================================================================
Ostap::Math::LorentzVectorWithError&
Ostap::Math::LorentzVectorWithError::operator-= 
( const Ostap::LorentzVector& right ) 
{
  vector4d() -= right ;
  return *this ;
}
// ============================================================================
Ostap::Math::LorentzVectorWithError&
Ostap::Math::LorentzVectorWithError::operator*= ( const double v ) 
{
  vector4d() *= v ;
  m_cov2     *= (v*v) ;
  return *this ;
}
// ============================================================================
Ostap::Math::LorentzVectorWithError&
Ostap::Math::LorentzVectorWithError::operator/= ( const double v ) 
{
  vector4d() /= v ;
  m_cov2     /= (v*v) ;
  return *this ;
}
// ============================================================================
Ostap::Math::LorentzVectorWithError
Ostap::Math::LorentzVectorWithError::__add__
( const Ostap::Math::LorentzVectorWithError& right ) const 
{
  LorentzVectorWithError tmp ( *this ) ;
  return tmp += right ;
}
// ============================================================================
Ostap::Math::LorentzVectorWithError
Ostap::Math::LorentzVectorWithError::__sub__
( const Ostap::Math::LorentzVectorWithError& right ) const 
{
  LorentzVectorWithError tmp ( *this ) ;
  return tmp -= right ;
}
// ============================================================================
Ostap::Math::LorentzVectorWithError
Ostap::Math::LorentzVectorWithError::__add__
( const Ostap::LorentzVector& right ) const 
{
  LorentzVectorWithError tmp ( *this ) ;
  return tmp += right ;
}
// ============================================================================
Ostap::Math::LorentzVectorWithError
Ostap::Math::LorentzVectorWithError::__add__
( const Ostap::Math::LorentzVectorWithError::VectorE& right ) const 
{
  LorentzVectorWithError tmp ( *this ) ;
  return tmp += right ;
}
// ============================================================================
Ostap::Math::LorentzVectorWithError
Ostap::Math::LorentzVectorWithError::__add__
( const Ostap::Math::LorentzVectorWithError::Vector& right ) const 
{
  LorentzVectorWithError tmp ( *this ) ;
  return tmp += right ;
}
// ============================================================================
Ostap::Math::LorentzVectorWithError
Ostap::Math::LorentzVectorWithError::__sub__
( const Ostap::LorentzVector& right ) const 
{
  LorentzVectorWithError tmp ( *this ) ;
  return tmp -= right ;
}
// ============================================================================
Ostap::Math::LorentzVectorWithError
Ostap::Math::LorentzVectorWithError::__sub__
( const Ostap::Math::LorentzVectorWithError::VectorE& right ) const 
{
  LorentzVectorWithError tmp ( *this ) ;
  return tmp -= right ;
}
// ============================================================================
Ostap::Math::LorentzVectorWithError
Ostap::Math::LorentzVectorWithError::__sub__
( const Ostap::Math::LorentzVectorWithError::Vector& right ) const 
{
  LorentzVectorWithError tmp ( *this ) ;
  return tmp -= right ;
}
// ============================================================================
Ostap::Math::LorentzVectorWithError
Ostap::Math::LorentzVectorWithError::__rsub__
( const Ostap::LorentzVector& right )  const 
{
  return LorentzVectorWithError( right - vector4d() , cov2() ) ;
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
Ostap::Math::LorentzVectorWithError::fillStream ( std::ostream& s ) const // printout 
{
  return s << "[" << "( " 
           << X () << " +- " << err ( m_cov2(0,0) ) << " , "
           << Y () << " +- " << err ( m_cov2(1,1) ) << " , "
           << Z () << " +- " << err ( m_cov2(2,2) ) << " ), "
           << E () << " +- " << err ( m_cov2(3,3) ) << " ]";
}  
// ============================================================================
// conversion to the string 
// ============================================================================
std::string
Ostap::Math::LorentzVectorWithError::toString   () const // conversion to the string 
{
  std::ostringstream s ;
  fillStream ( s ) ;
  return s.str() ;
}
// ============================================================================
// chi2 distance 
// ============================================================================
double Ostap::Math::LorentzVectorWithError::chi2 
( const Ostap::Math::LorentzVectorWithError& right ) const 
{
  Covariance s_cov2 ( cov2() ) ;
  s_cov2 += right.cov2() ;
  // use Manuel's inverter:
  const bool ok = s_cov2.InvertChol() ;
  if  ( !ok ) { return -1 ; }                                 // RETURN  
  /// calculate chi2 
  return Ostap::Math::Similarity ( vector4d() - right.vector4d() , s_cov2 ) ;
}
// ============================================================================
// chi2 distance 
// ============================================================================
double Ostap::Math::LorentzVectorWithError::chi2
( const Ostap::LorentzVector& right ) const 
{
  Covariance s_cov2 ( cov2() ) ;
  // use Manuel's inverter:
  const bool ok = s_cov2.InvertChol() ;
  if  ( !ok ) { return -1 ; }                                 // RETURN  
  /// calculate chi2 
  return Ostap::Math::Similarity ( vector4d() - right , s_cov2 ) ;
}

// ============================================================================
// chi2 distance 
// ============================================================================
double Ostap::Math::LorentzVectorWithError::chi2 
( const Ostap::Math::LorentzVectorWithError::VectorE& right ) const 
{
  Covariance s_cov2 ( cov2() ) ;
  s_cov2 += right.cov2() ;
  // use Manuel's inverter:
  const bool ok = s_cov2.InvertChol() ;
  if  ( !ok ) { return -1 ; }                                 // RETURN  
  /// calculate chi2
  Vector vct ;
  Ostap::Math::geo2LA ( vector4d() , vct ) ;
  vct -= right.value() ;
  //
  return ROOT::Math::Similarity ( vct , s_cov2 ) ;
}
// ============================================================================
// chi2 distance 
// ============================================================================
double Ostap::Math::LorentzVectorWithError::chi2
( const Ostap::Math::LorentzVectorWithError::Vector& right ) const 
{
  Covariance s_cov2 ( cov2() ) ;
  // use Manuel's inverter:
  const bool ok = s_cov2.InvertChol() ;
  if  ( !ok ) { return -1 ; }                                 // RETURN  
  /// calculate chi2 
  Vector vct ;
  Ostap::Math::geo2LA ( vector4d() , vct ) ;
  vct -= right ;
  //
  return ROOT::Math::Similarity ( vct , s_cov2 ) ;
}
// ============================================================================
void Ostap::Math::LorentzVectorWithError::asVector 
( Ostap::Math::LorentzVectorWithError::Vector& data ) const 
{ Ostap::Math::geo2LA ( vector4d() , data ) ; }
// ============================================================================
void Ostap::Math::LorentzVectorWithError::asVector 
( Ostap::Math::LorentzVectorWithError::VectorE& data ) const 
{
  Ostap::Math::geo2LA ( vector4d() , data.value() ) ; 
  data.setCov2( cov2() ) ;
}
// ============================================================================
Ostap::Math::LorentzVectorWithError::VectorE 
Ostap::Math::LorentzVectorWithError::asVector () const 
{
  Ostap::Math::LorentzVectorWithError::VectorE data ;
  Ostap::Math::geo2LA ( vector4d() , data.value() ) ;
  data.setCov2( cov2() ) ;
  return data ;
}
// ============================================================================
Ostap::Math::LorentzVectorWithError::Vector 
Ostap::Math::LorentzVectorWithError::asVector4 () const 
{
  Ostap::Math::LorentzVectorWithError::Vector data ;
  Ostap::Math::geo2LA ( vector4d() , data ) ;
  return data ;
}
// ============================================================================

// ============================================================================
Ostap::Math::LorentzVectorWithError&
Ostap::Math::LorentzVectorWithError::__imul__ ( const double v ) 
{
  vector4d() *= v ;
  m_cov2     *= (v*v) ;
  return *this ;
}
// ============================================================================
Ostap::Math::LorentzVectorWithError&
Ostap::Math::LorentzVectorWithError::__itruediv__ ( const double v ) 
{
  vector4d() /= v ;
  m_cov2     /= (v*v) ;
  return *this ;
}
// ============================================================================

// ============================================================================
Ostap::Math::LorentzVectorWithError 
Ostap::Math::LorentzVectorWithError::__mul__ ( const double v ) const
{
  Ostap::Math::LorentzVectorWithError tmp (*this) ;
  return ( tmp *= v ) ;
}
// ============================================================================

// ============================================================================
Ostap::Math::LorentzVectorWithError 
Ostap::Math::LorentzVectorWithError::__truediv__ ( const double v ) const
{
  Ostap::Math::LorentzVectorWithError tmp (*this) ;
  return ( tmp /= v ) ;
}
// ============================================================================

// ============================================================================
// set of helper functions 
// ============================================================================
namespace 
{
  // ==========================================================================
  /// almost zero  ? 
  const Ostap::Math::Zero<double>      s_zero { } ;
  /// almost equal ? 
  const Ostap::Math::Equal_To<double>  s_equal { } ;
  // ==========================================================================
}
// ============================================================================
/*  calculate the mass with uncertainty 
 *  @param mom 4-momentum 
 *  @param cov covariance 
 *  @return mass with uncertainty 
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::Kinematics::mass 
( const Ostap::LorentzVector& mom , 
  const Ostap::SymMatrix4x4&  cov ) 
{
  const double m = mom.M() ;
  return 
    m <= 0 || s_zero ( m ) ? 
    Ostap::Math::ValueWithError ( m ) : 
    Ostap::Math::ValueWithError ( m , Ostap::Math::sigma2mass ( mom , cov ) ) ;
}
// ============================================================================
/*  calculate the scalar momentum (p) with uncertainty 
 *  @param mom 4-momentum 
 *  @param cov covariance 
 *  @return p with uncertainty 
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::Kinematics::momentum 
( const Ostap::LorentzVector& mom , 
  const Ostap::SymMatrix4x4&  cov ) 
{
  const double p = mom.P() ;
  return 
    p <= 0 || s_zero ( p ) ?  
    Ostap::Math::ValueWithError ( p ) : 
    Ostap::Math::ValueWithError ( p , Ostap::Math::sigma2p ( mom , cov ) ) ;
}
// ============================================================================
/*  calculate the rapidity (y) with uncertainty 
 *  @param mom 4-momentum 
 *  @param cov covariance 
 *  @return y with uncertainty 
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::Kinematics::rapidity 
( const Ostap::LorentzVector& mom , 
  const Ostap::SymMatrix4x4&  cov ) 
{
  const double y  = mom.Rapidity() ;
  const double e  = std::abs ( mom.E  ()  ) ;
  const double pz = std::abs ( mom.Pz ()  ) ;
  return 
    e <= pz || s_equal ( e , pz  ) ?  
    Ostap::Math::ValueWithError ( y ) : 
    Ostap::Math::ValueWithError ( y , Ostap::Math::sigma2y ( mom , cov ) ) ;  
}
// ============================================================================
/*  calculate the pseudo-rapidity (eta) with uncertainty 
 *  @param mom 4-momentum 
 *  @param cov covariance 
 *  @return eta with uncertainty 
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::Kinematics::pseudorapidity 
( const Ostap::LorentzVector& mom , 
  const Ostap::SymMatrix4x4&  cov ) 
{
  const double eta = mom.Eta   () ;
  const double pt2 = mom.Perp2 () ;
  if ( 0 >= pt2 || s_zero ( pt2 ) ) { return Ostap::Math::ValueWithError ( eta ) ; }
  const double p   = mom.P     () ;
  if ( 0 >= p   || s_zero ( p   ) ) { return Ostap::Math::ValueWithError ( eta ) ; }
  //
  const double pz  = mom.Pz    () ;
  const double c   = - pz / ( pt2 * p ) ;
  //
  // get the vector d(Eta)/dp_i :
  ROOT::Math::SVector<double,4> dEtadP_i;
  dEtadP_i [0] =  c * mom.Px() ;
  dEtadP_i [1] =  c * mom.Py() ;
  dEtadP_i [2] =  1 / p ;
  dEtadP_i [3] =  0.0 ;
  //
  const double s2eta = ROOT::Math::Similarity ( cov , dEtadP_i ) ;
  return 
    s2eta <= 0 || s_zero ( s2eta )  ? 
    Ostap::Math::ValueWithError ( eta         ) : 
    Ostap::Math::ValueWithError ( eta , s2eta ) ;
}
// ============================================================================
/*  calculate the transverse momentum (pt) with uncertainty 
 *  @param mom 4-momentum 
 *  @param cov covariance 
 *  @return pt with uncertainty 
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::Kinematics::transverseMomentum  
( const Ostap::LorentzVector& mom , 
  const Ostap::SymMatrix4x4&  cov ) 
{
  const double pt = mom.Pt () ;
  return 
    pt <= 0 || s_zero ( pt  ) ?  
    Ostap::Math::ValueWithError ( pt ) : 
    Ostap::Math::ValueWithError ( pt , Ostap::Math::sigma2pt ( mom , cov ) ) ;  
}
// ============================================================================
/*  calculate the squared transverse mass  (mT2) with uncertainty 
 *  \f$ m_T^2 = e^2 - p_z^2  = m^2 + p_T^2 \f$ 
 *  @param mom 4-momentum 
 *  @param cov covariance 
 *  @return mT2 with uncertainty 
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::Kinematics::transverseMass2
( const Ostap::LorentzVector& mom , 
  const Ostap::SymMatrix4x4&  cov ) 
{
  const double mt2 = mom.Mt2() ;
  //
  // get the vector d(Mt2)/dp_i :
  ROOT::Math::SVector<double,4> dMt2dP_i;
  dMt2dP_i [0] =  0.0 ;
  dMt2dP_i [1] =  0.0 ;
  dMt2dP_i [2] = -2 * mom.Pz() ;
  dMt2dP_i [3] =  2 * mom.E () ;
  //
  const double s2mt2 = ROOT::Math::Similarity ( cov , dMt2dP_i ) ;
  return 
    s2mt2 <= 0 || s_zero ( s2mt2 )  ? 
    Ostap::Math::ValueWithError ( mt2         ) : 
    Ostap::Math::ValueWithError ( mt2 , s2mt2 ) ;
}
// ============================================================================
/*  calculate the transverse mass  (mT) with uncertainty 
 *  \f$ m_T  = \sqrt{e^2 - p_z^2}  = \sqrt{m^2 + p_T^2} \f$ 
 *  @param mom 4-momentum 
 *  @param cov covariance 
 *  @return mT2 with uncertainty 
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::Kinematics::transverseMass
( const Ostap::LorentzVector& mom , 
  const Ostap::SymMatrix4x4&  cov ) 
{
  const Ostap::Math::ValueWithError mt2 = transverseMass2 ( mom , cov ) ;
  return 
    !s_zero ( mt2.value() ) ?    Ostap::Math::signed_sqrt ( mt2 ) :
    Ostap::Math::ValueWithError ( Ostap::Math::signed_sqrt ( mt2.value() ) ) ;
}
// ============================================================================
/*  calculate the squared transverse energy (eT2) with uncertainty 
 *  \f$ e_T^2 = \frac{e^2p_t^2}{p^2} \f$ 
 *  @param mom 4-momentum 
 *  @param cov covariance 
 *  @return eT2 with uncertainty 
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::Kinematics::transverseEnergy2
( const Ostap::LorentzVector& mom , 
  const Ostap::SymMatrix4x4&  cov ) 
{
  const double et2 = mom.Et2() ;
  //
  const double p2  = mom.P2 () ;
  if ( 0 >= p2 || s_zero ( p2 ) ) { return Ostap::Math::ValueWithError ( et2 ) ; }
  //
  const double pt2 = mom.Perp2() ;
  const double e   = mom.E    () ;
  const double e2  = e * e       ;
  const double c1  = 2 * e2 / p2 ;
  const double c2  = pt2    / p2 ;
  //
  // get the vector d(Et2)/dp_i :
  ROOT::Math::SVector<double,4> dEt2dP_i;
  dEt2dP_i [0] =  c1 * ( 1 - c2 ) * mom.Px () ;
  dEt2dP_i [1] =  c1 * ( 1 - c2 ) * mom.Py () ;
  dEt2dP_i [2] =  c1 * ( 0 - c2 ) * mom.Pz () ;
  dEt2dP_i [3] =  2 * e * c2 ;
  //
  const double s2et2 = ROOT::Math::Similarity ( cov , dEt2dP_i ) ;
  return 
    s2et2 <= 0 || s_zero ( s2et2 )  ? 
    Ostap::Math::ValueWithError ( et2         ) : 
    Ostap::Math::ValueWithError ( et2 , s2et2 ) ;
}
// ============================================================================
/** calculate the transverse energy  (eT) with uncertainty 
 *  \f$ e_T^2 = \frac{e p_t}{p} \f$ 
 *  @param mom 4-momentum 
 *  @param cov covariance 
 *  @return eT with uncertainty 
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::Kinematics::transverseEnergy
( const Ostap::LorentzVector& mom , 
  const Ostap::SymMatrix4x4&  cov ) 
{
  const Ostap::Math::ValueWithError et2 = transverseEnergy2 ( mom , cov ) ;
  return 
    !s_zero ( et2.value() ) ?    Ostap::Math::signed_sqrt ( et2 ) :
    Ostap::Math::ValueWithError ( Ostap::Math::signed_sqrt ( et2.value() ) ) ;
}
// ============================================================================
/*  calculate the kinetic energy  (eK) with uncertainty 
 *  \f$ e_K = e - m \f$ 
 *  @param mom 4-momentum 
 *  @param cov covariance 
 *  @return eK with uncertainty 
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::Kinematics::kineticEnergy
( const Ostap::LorentzVector& mom , 
  const Ostap::SymMatrix4x4&  cov ) 
{
  const double m2 = mom.M2() ;
  if ( 0 >= m2 || s_zero ( m2 ) ) 
  {
    const double c2ee = cov( 3,3 ) ;
    return 
      0 >= c2ee || s_zero ( c2ee ) ?  
      Ostap::Math::ValueWithError ( mom.E()        ) : 
      Ostap::Math::ValueWithError ( mom.E() , c2ee ) ;
  }
  const double m   = std::sqrt  ( m2 ) ;
  const double eK  = mom.E() - m       ;
  //
  // get the vector d(Ek)/dp_i :
  ROOT::Math::SVector<double,4> dEkdP_i;
  dEkdP_i [0] =      mom.Px() / m ;
  dEkdP_i [1] =      mom.Py() / m ;
  dEkdP_i [2] =      mom.Pz() / m ;
  dEkdP_i [3] =  1 - mom.E () / m ;
  //
  const double s2ek = ROOT::Math::Similarity ( cov , dEkdP_i ) ;
  return 
    s2ek <= 0 || s_zero ( s2ek )  ? 
    Ostap::Math::ValueWithError ( eK        ) : 
    Ostap::Math::ValueWithError ( eK , s2ek ) ;
}
// ============================================================================
/* calculate the phi (asymuthal angle) with uncertainy 
 *  \f$ \tan \phi = \frac{p_y}{p_x} \f$ 
 *  @param mom 4-momentum 
 *  @param cov covariance 
 *  @return phi with uncertainty 
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::Kinematics::phi
( const Ostap::LorentzVector& mom , 
  const Ostap::SymMatrix4x4&  cov ) 
{
  const double phi = mom.Phi  () ;
  const double pt2 = mom.Perp2() ;
  //
  if ( 0 >= pt2 || s_zero ( pt2 ) ) { return Ostap::Math::ValueWithError ( phi ) ; }
  //
  // get the vector d(Phi)/dp_i :
  ROOT::Math::SVector<double,4> dPhidP_i;
  dPhidP_i [0] = - mom.Py() / pt2 ;
  dPhidP_i [1] =   mom.Px() / pt2 ;
  dPhidP_i [2] =   0.0 ;
  dPhidP_i [3] =   0.0 ;
  //
  const double s2phi = ROOT::Math::Similarity ( cov , dPhidP_i ) ;
  return 
    s2phi <= 0 || s_zero ( s2phi )  ? 
    Ostap::Math::ValueWithError ( phi         ) : 
    Ostap::Math::ValueWithError ( phi , s2phi ) ;
}
// ============================================================================
/*  calculate the theta (polar angle) with uncertainy 
 *  \f$ \tan \theta = \frac{p_T}{p_z} \f$ 
 *  @param mom 4-momentum 
 *  @param cov covariance 
 *  @return theta with uncertainty 
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::Kinematics::theta
( const Ostap::LorentzVector& mom , 
  const Ostap::SymMatrix4x4&  cov ) 
{
  const double theta = mom.Theta  () ;
  const double pt2   = mom.Perp2  () ;
  if ( 0 >= pt2 || s_zero ( pt2 ) ) 
  { return Ostap::Math::ValueWithError ( theta ) ; }
  const double p2    = mom.P2     () ;
  if ( 0 >= p2  || s_zero ( p2  ) ) 
  { return Ostap::Math::ValueWithError ( theta ) ; }
  //
  const double pt    = mom.Pt     () ;
  const double pz    = mom.Pz     () ;
  const double c     = pz / ( p2 * pt ) ;
  //
  // get the vector d(Theta)/dp_i :
  ROOT::Math::SVector<double,4> dThetadP_i;
  dThetadP_i [0] =  c * mom.Px() ;
  dThetadP_i [1] =  c * mom.Py() ;
  dThetadP_i [2] =  - pt / p2    ;
  dThetadP_i [3] =  0.0 ;
  //
  const double s2theta = ROOT::Math::Similarity ( cov , dThetadP_i ) ;
  return 
    s2theta <= 0 || s_zero ( s2theta )  ? 
    Ostap::Math::ValueWithError ( theta           ) : 
    Ostap::Math::ValueWithError ( theta , s2theta ) ;
}
// ============================================================================
/* calculate the transverse kinetic energy  (eTK) with uncertainty 
 *  \f$ e_{T,k} = m_T - m \f$ 
 *  @param mom 4-momentum 
 *  @param cov covariance 
 *  @return eTK with uncertainty 
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::Math::Kinematics::transverseKineticEnergy
( const Ostap::LorentzVector& mom , 
  const Ostap::SymMatrix4x4&  cov )
{
  //
  const double e    = mom.E  () ;
  const double px   = mom.Px () ;
  const double py   = mom.Py () ;
  const double pz   = mom.Pz () ;
  const double m    = mom.M  () ;
  //
  const double mt2  = m * m + px * px + py * py ;
  const double mt   = std::sqrt ( mt2 ) ;
  //
  const double etk  = mt - m ;
  if ( m <= 0 || s_zero ( m ) ) { return Ostap::Math::ValueWithError ( etk ) ; }
  //
  // get the vector d(etk)/dp_i :
  ROOT::Math::SVector<double,4> dEtk_dP;
  dEtk_dP [0] = px / m                  ;
  dEtk_dP [1] = py / m                  ;
  dEtk_dP [2] = pz / m * ( 1 - m / mt ) ;
  dEtk_dP [3] = -e / m * ( 1 - m / mt ) ;
  //
  const double s2etk = ROOT::Math::Similarity ( cov , dEtk_dP ) ;
  return 
    s2etk <= 0 || s_zero ( s2etk )  ? 
    Ostap::Math::ValueWithError ( etk         ) : 
    Ostap::Math::ValueWithError ( etk , s2etk ) ;
}

// ============================================================================
Ostap::Math::LorentzVectorWithError 
Ostap::Math::LorentzVectorWithError::mean 
( const Ostap::Math::LorentzVectorWithError&          right ) const 
{ return asVector ().mean ( right.asVector() ) ; }
// ============================================================================
Ostap::Math::LorentzVectorWithError 
Ostap::Math::LorentzVectorWithError::mean 
( const Ostap::Math::LorentzVectorWithError::VectorE& right ) const 
{ return asVector ().mean ( right            ) ; }
// ============================================================================

// ============================================================================
// The END 
// ============================================================================

