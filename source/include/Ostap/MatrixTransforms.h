// $Id$
// ============================================================================
#ifndef OSTAP_MATRIXTRANSFORMS_H 
#define OSTAP_MATRIXTRANSFORMS_H 1
// ============================================================================
// Include files
// ============================================================================
// ROOT
// ============================================================================
#include "Math/SMatrix.h"
#include "Math/SVector.h"
#include "Math/Point3D.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Math.h"
#include "Ostap/MatrixUtils.h"
// ============================================================================
/** @file Ostap/MatrixTransforms.h
 *  The collection of useful utils to minimize the conversion from 
 *  geometry&kinematical vectors into linear algebra vectors 
 *   In particular it includes:
 *     - conversion from geometry&kinematical vectors into Linear Algebra vectors
 *     - conversion from Linear Algebra vectors into geometry&kinematical vectors
 *     - evaluation of various "chi2"-like values, like chi2-distance between 
 *       two 3D or 4D-vectors. E.g. "vicinity" of two points, or two momenta 
 *     - conversion form "track" to 4-momenta representation (by Sean BRISBANE)
 *     - the transition matrix for the conversion form "track" to 
 *        4-momenta representation (by Sean BRISBANE)
 *
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2008-01-15
 */
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {   
    // ========================================================================
    /** fill  Linear Algebra - vector from 3D-point
     *
     *  @code
     *   
     *  const Ostap::XYZPoint xyzv = ... ;
     *  Ostap::Vector4 lav ;
     *  
     *  // fill Linear Algebra vector from 3D-point
     *  geo2LA ( xyzv , lav ) ; 
     *
     *  @endcode 
     *  
     *  @param  source (input)  3D-point
     *  @param  dest   (output) Linear Algebra vector 
     *  @return linear algebra vector 
     *  @author Vanya BELYAEV ibelyaev@physics.syr.edu
     *  @date 2006-05-24
     */
    template <class C,class T>
    inline 
    const ROOT::Math::SVector<T,3>& 
    geo2LA
    ( const ROOT::Math::PositionVector3D<C>& source , 
      ROOT::Math::SVector<T,3>&              dest   ) 
    {
      dest ( 0 ) = source.X () ;
      dest ( 1 ) = source.Y () ;
      dest ( 2 ) = source.Z () ;
      return dest ;
    } 
    // ========================================================================
    /** fill  Linear Algebra - vector from 3D-vector 
     *
     *  @code
     *   
     *  const Ostap::XYZVector xyzv = ... ;
     *  Ostap::Vector4 lav ;
     *  
     *  // fill Linear Algebra vector from 3D-Vector 
     *  geo2LA ( xyzv , lav ) ; 
     *
     *  @endcode 
     *  
     *  @param source (input)  3D-Vector 
     *  @param dest   (output) Linear Algebra vector 
     *  @return linear algebra vector 
     *  @author Vanya BELYAEV ibelyaev@physics.syr.edu
     *  @date 2006-05-24
     */
    template <class C,class T>
    inline 
    const ROOT::Math::SVector<T,3>& 
    geo2LA 
    ( const ROOT::Math::DisplacementVector3D<C>& source , 
      ROOT::Math::SVector<T,3>&                  dest   ) 
    {
      dest ( 0 ) = source.X () ;
      dest ( 1 ) = source.Y () ;
      dest ( 2 ) = source.Z () ;
      return dest ;
    }
    // ========================================================================
    /** fill  Linear Algebra - vector from 4D-vector 
     *
     *  @code
     *   
     *  const Ostap::LorenztVector lorv = ... ;
     *  Ostap::Vector4 lav ;
     *  
     *  // fill Linear Algebra vector from Lorenz Vector 
     *  geo2LA ( lorv , lav ) ; 
     *
     *  @endcode 
     *  
     *  @param source (input) Lorentz Vector 
     *  @param dest   (output) Linear Algebra vector 
     *  @return linear algebra vector 
     *  @author Vanya BELYAEV ibelyaev@physics.syr.edu
     *  @date 2006-05-24
     */
    template <class C,class T>
    inline 
    const ROOT::Math::SVector<T,4>& 
    geo2LA 
    ( const ROOT::Math::LorentzVector<C>& source , 
      ROOT::Math::SVector<T,4>&           dest  ) 
    {
      dest ( 0 ) = source.X () ;
      dest ( 1 ) = source.Y () ;
      dest ( 2 ) = source.Z () ;
      dest ( 3 ) = source.E () ;
      return dest ;
    } 
    // ========================================================================
    /** fill  Linear Algebra 3-vector from the spatial 
     *  components of 4D-(Lorentz)vector 
     *
     *  @code
     *   
     *  const Ostap::LorenztVector lorv = ... ;
     *  Ostap::Vector3 v3 ;
     *  
     *  // fill Linear Algebra vector from Lorenz Vector 
     *  geo2LA ( lorv , v3 ) ; 
     *
     *  @endcode 
     *  
     *  @param source (input) Lorentz Vector 
     *  @param dest   (output) Linear Algebra vector 
     *  @return linear algebra vector 
     *  @author Vanya BELYAEV ibelyaev@physics.syr.edu
     *  @date 2006-05-24
     */
    template <class C,class T>
    inline 
    const ROOT::Math::SVector<T,3>& 
    geo2LA 
    ( const ROOT::Math::LorentzVector<C>& source , 
      ROOT::Math::SVector<T,3>&           dest  ) 
    {
      dest ( 0 ) = source.X () ;
      dest ( 1 ) = source.Y () ;
      dest ( 2 ) = source.Z () ;
      return dest ;
    } 
    // ========================================================================
    /** fill  3D-point from Linear Algebra vector 
     *
     *  @code
     *   
     *  Ostap::XYZPoint xyzv = ... ;
     *  const Ostap::Vector4 lav ;
     *  
     *  // fill 3D-point from Linear Algebra vector 
     *  la2geo ( xyzv , lav ) ; 
     *
     *  @endcode 
     *  
     *  @param source (input)  Linear Algebra vector 
     *  @param dest   (output) 3D-point 
     *  @return linear algebra vector 
     *  @author Vanya BELYAEV ibelyaev@physics.syr.edu
     *  @date 2006-05-24
     */
    template <class C,class T>
    inline 
    const ROOT::Math::PositionVector3D<C>&
    la2geo 
    ( const ROOT::Math::SVector<T,3>&  source , 
      ROOT::Math::PositionVector3D<C>& dest   ) 
    {
      dest.SetX ( source ( 0 ) ) ;      
      dest.SetY ( source ( 1 ) ) ;
      dest.SetZ ( source ( 2 ) ) ;
      return dest ;
    } 
    // ========================================================================
    /** fill  3D-vector from Linear Algebra vector 
     *
     *  @code
     *   
     *  Ostap::XYZVector xyzv = ... ;
     *  const Ostap::Vector4 lav ;
     *  
     *  // fill 3D-Vector from Linear Algebra vector 
     *  la2geo ( xyzv , lav ) ; 
     *
     *  @endcode 
     *  
     *  @param source (input)  Linear Algebra vector 
     *  @param dest   (output) 3D-vector 
     *  @return linear algebra vector 
     *  @author Vanya BELYAEV ibelyaev@physics.syr.edu
     *  @date 2006-05-24
     */
    template <class C,class T>
    inline 
    const ROOT::Math::DisplacementVector3D<C>&
    la2geo 
    ( const ROOT::Math::SVector<T,3>&      source , 
      ROOT::Math::DisplacementVector3D<C>& dest   ) 
    {
      dest.SetX ( source ( 0 ) ) ;      
      dest.SetY ( source ( 1 ) ) ;
      dest.SetZ ( source ( 2 ) ) ;
      return dest ;
    } 
    // ========================================================================
    /** fill  Lorentz vector from Linear Algebra vector 
     *
     *  @code
     *   
     *  const Ostap::LorenztVector lorv = ... ;
     *  Ostap::Vector4 lav ;
     *  
     *  // fill Linear Algebra vector from Lorenz Vector 
     *  geo2LA ( lorv , lav ) ; 
     *
     *  @endcode 
     *  
     *  @param source (input) Lorentz Vector 
     *  @param dest   (output) Linear Algebra vector 
     *  @return linear algebra vector 
     *  @author Vanya BELYAEV ibelyaev@physics.syr.edu
     *  @date 2006-05-24
     */
    template <class C,class T>
    inline 
    const ROOT::Math::LorentzVector<C>&  
    la2geo 
    ( const ROOT::Math::SVector<T,4>& source , 
      ROOT::Math::LorentzVector<C>&   dest   )  
    {
      dest.SetPx ( source ( 0 ) ) ;
      dest.SetPy ( source ( 1 ) ) ;
      dest.SetPz ( source ( 2 ) ) ;
      dest.SetE  ( source ( 3 ) ) ;
      return dest ;
    }
    // ========================================================================
    /** construct similarity("chi2") using 3D-vector 
     *
     *  E.g. one can ask the "chi2"-distance inbetween vertices:
     *
     *  @code 
     *
     *  const LHCb::Vertex* v1 = ... ;
     *  const LHCb::Vertex* v2 = ... ;
     *  
     *  int ifail = 0 ;
     *  // evaluate the chi2 distance 
     *  const double chi2 = Ostap::Math::Similarity
     *     ( v1->position() - v2.position() , 
     *      ( v1->covMatrix() + v2->covMatrix() ).Sinverse( ifail ) ) ;
     *  if ( 0 != ifail ) { ... error here ... } ;
     *
     *  always() << " Chi2 distance between vertices is " << chi2 << endmsg ;
     *
     *  @endcode 
     *
     *  @param matrix (input) symmetric (3x3) matrix used for similarity
     *  @param delta  (input) 3D- vector 
     *  @return reult of v^T*M*v (similarity) operation
     *  @author Vanya BELYAEV ibelyaev@physics.syr.edu
     *  @date 2006-05-24
     */ 
    template <class C,class T>
    inline T 
    Similarity
    ( const ROOT::Math::DisplacementVector3D<C>&                     delta  , 
      const ROOT::Math::SMatrix<T,3,3,ROOT::Math::MatRepSym<T,3> > & matrix ) 
    {
      ROOT::Math::SVector<T,3> tmp ;
      return ROOT::Math::Similarity ( geo2LA ( delta , tmp )  , matrix ) ;
    } 
    // ========================================================================
    /** construct similarity("chi2") using 3D-vector 
     *
     *  E.g. one can ask the "chi2"-distance inbetween vertices:
     *
     *  @code 
     *
     *  const LHCb::Vertex* v1 = ... ;
     *  const LHCb::Vertex* v2 = ... ;
     *  
     *  int ifail = 0 ;
     *  // evaluate the chi2 distance 
     *  const double chi2 = Ostap::Math::Similarity
     *    ( ( v1->covMatrix() + v2->covMatrix() ).Sinverse( ifail ) ,
     *        v1->position() - v2.position() ) ;
     *  if ( 0 != ifail ) { ... error here ... } ;
     *
     *  always() << " Chi2 distance between vertices is " << chi2 << endmsg ;
     *
     *  @endcode 
     *
     *  @param matrix (input) symmetric (3x3) matrix used for similarity
     *  @param delta  (input) 3D- vector 
     *  @return result of v^T*M*v (similarity) operation
     *  @author Vanya BELYAEV ibelyaev@physics.syr.edu
     *  @date 2006-05-24
     */ 
    template <class C,class T>
    inline T 
    Similarity
    ( const ROOT::Math::SMatrix<T,3,3,ROOT::Math::MatRepSym<T,3> > & matrix , 
      const ROOT::Math::DisplacementVector3D<C>&                     delta  ) 
    { return Similarity ( delta , matrix ) ; }
    // ========================================================================
    /** construct similarity("chi2") using 4D-vector 
     *
     *  E.g. one can ask the "chi2"-distance inbetween momenta of particles
     *
     *  @code 
     *
     *  const LHCb::Particle* p1 = ... ;
     *  const LHCb::Particle* p2 = ... ;
     *  
     *  int ifail = 0 ;
     *  // evaluate the chi2 distance 
     *  const double chi2 = Ostap::Math::Similarity
     *     (  p1->momentum() - p2.momentum() ,
     *     ( p1->momCovMatrix() + p2->momCovMatrix() ).Sinverse( ifail ) ) ;
     *  if ( 0 != ifail ) { ... error here ... } ;
     *
     *  always() << " Chi2 distance in momenta is " << chi2 << endmsg ;
     *
     *  @endcode 
     *
     *  @param delta  (input) Lorentz vector 
     *  @param matrix (input) symmetric (4x4) matrix used for similarity
     *  @return result of v^T*M*v (similarity) operation
     *  @author Vanya BELYAEV ibelyaev@physics.syr.edu
     *  @date 2006-05-24
     */ 
    template <class C,class T>
    inline T 
    Similarity
    ( const ROOT::Math::LorentzVector<C>&                            delta  , 
      const ROOT::Math::SMatrix<T,4,4,ROOT::Math::MatRepSym<T,4> > & matrix ) 
    {
      ROOT::Math::SVector<T,4> tmp ;
      return ROOT::Math::Similarity ( geo2LA( delta , tmp ) , matrix ) ;
    }
    // ========================================================================
    /** construct similarity("chi2") using 4D-vector 
     *
     *  E.g. one can ask the "chi2"-distance inbetween momenta of particles
     *
     *  @code 
     *
     *  const LHCb::Particle* p1 = ... ;
     *  const LHCb::Particle* p2 = ... ;
     *  
     *  int ifail = 0 ;
     *  // evaluate the chi2 distance 
     *  const double chi2 = Ostap::Math::Similarity
     *    * ( p1->momCovMatrix() + p2->momCovMatrix() ).Sinverse( ifail ) , 
     *       p1->momentum() - p2.momentum() ) ;
     *  if ( 0 != ifail ) { ... error here ... } ;
     *
     *  always() << " Chi2 distance in momenta is " << chi2 << endmsg ;
     *
     *  @endcode 
     *
     *  @param matrix (input) symmetric (4x4) matrix used for similarity
     *  @param delta  (input) Lorentz vector 
     *  @return result of v^T*M*v (similarity) operation
     *  @author Vanya BELYAEV ibelyaev@physics.syr.edu
     *  @date 2006-05-24
     */ 
    template <class C,class T>
    inline T 
    Similarity
    ( const ROOT::Math::SMatrix<T,4,4,ROOT::Math::MatRepSym<T,4> > & matrix ,
      const ROOT::Math::LorentzVector<C>&                            delta  ) 
    { return Similarity ( delta , matrix ) ; } 
    // ========================================================================
    /** increment  Position-Vector with 3-component linear vector 
     *  
     *  @code 
     *
     *  Ostap::XYZPoint     v1 = ... ;
     *  const Ostap::Vector3       v2 = ... ;
     *  
     *  // update vector with LA vector:
     *  Ostap::Math::add ( v1 , v2 ) ;
     *
     *  @endcode 
     *  
     *  @param v1 (input/output) LorentzVector to be updated 
     *  @param v2 (input) Linear Algebra vector, to be added to LorentzVector 
     *  @return the updated LorenzVector 
     *  @author Vanya BELYAEV ibelyaev@physics.syr.edu
     *  @date 2006-05-24
     */
    template <class C, class T>
    inline 
    const ROOT::Math::PositionVector3D<C>& 
    add 
    ( ROOT::Math::PositionVector3D<C>& v1 , const ROOT::Math::SVector<T,3> & v2 ) 
    { return v1 += ROOT::Math::PositionVector3D<C>
        ( v2 ( 0 )  , v2 ( 1 )  , v2 ( 2 ) ) ; }
    // ========================================================================
    /** increment  Displacement-Vector with 3-component linear vector 
     *  
     *  @code 
     *
     *  Ostap::XYZVector v1 = ... ;
     *  const Ostap::Vector3       v2 = ... ;
     *  
     *  // update vector with LA vector:
     *  Ostap::Math::add ( v1 , v2 ) ;
     *
     *  @endcode 
     *  
     *  @param v1 (input/output) LorentzVector to be updated 
     *  @param v2 (input) Linear Algebra vector, to be added to LorentzVector 
     *  @return the updated LorenzVector 
     *  @author Vanya BELYAEV ibelyaev@physics.syr.edu
     *  @date 2006-05-24
     */
    template <class C, class T>
    inline 
    const ROOT::Math::DisplacementVector3D<C>& 
    add 
    ( ROOT::Math::DisplacementVector3D<C>& v1 , const ROOT::Math::SVector<T,3> & v2 ) 
    { return v1 += ROOT::Math::DisplacementVector3D<C>
        ( v2 ( 0 )  , v2 ( 1 )  , v2 ( 2 ) ) ; }
    // ========================================================================
    /** increment  LorentzVector with 4-component linear vector 
     *  
     *  @code 
     *
     *  Ostap::LorentzVector v1 = ... ;
     *  const Ostap::Vector4       v2 = ... ;
     *  
     *  // update Lorentz vector with LA vector:
     *  Ostap::Math::add ( v1 , v2 ) ;
     *
     *  @endcode 
     *  
     *  @param v1 (input/output) LorentzVector to be updated 
     *  @param v2 (input) Linear Algebra vector, to be added to LorentzVector 
     *  @return the updated LorenzVector 
     *  @author Vanya BELYAEV ibelyaev@physics.syr.edu
     *  @date 2006-05-24
     */
    template <class C, class T>
    inline 
    const ROOT::Math::LorentzVector<C>& 
    add 
    ( ROOT::Math::LorentzVector<C>& v1 , const ROOT::Math::SVector<T,4> & v2 ) 
    { return v1 += ROOT::Math::LorentzVector<C>
        ( v2 ( 0 )  , v2 ( 1 )  , v2 ( 2 ) , v2 ( 3 ) ) ; }
    // ========================================================================    
    /** decrement  Position-Vector with 3-component linear vector 
     *  
     *  @code 
     *
     *  Ostap::XYZPoint     v1 = ... ;
     *  const Ostap::Vector3       v2 = ... ;
     *  
     *  // update vector with LA vector:
     *  Ostap::Math::sub ( v1 , v2 ) ;
     *
     *  @endcode 
     *  
     *  @param v1 (input/output) LorentzVector to be updated 
     *  @param v2 (input) Linear Algebra vector, to be subtracted 3D-point
     *  @return the updated LorenzVector 
     *  @author Vanya BELYAEV ibelyaev@physics.syr.edu
     *  @date 2006-05-24
     */
    template <class C, class T>
    inline 
    const ROOT::Math::PositionVector3D<C>& 
    sub 
    ( ROOT::Math::PositionVector3D<C>& v1 , const ROOT::Math::SVector<T,3> & v2 ) 
    { return v1 -= ROOT::Math::PositionVector3D<C>
        ( v2 ( 0 )  , v2 ( 1 )  , v2 ( 2 ) ) ; }
    // ========================================================================
    /** decrement  Displacement-Vector with 3-component linear vector 
     *  
     *  @code 
     *
     *  Ostap::XYZVector v1 = ... ;
     *  const Ostap::Vector3       v2 = ... ;
     *  
     *  // update vector with LA vector:
     *  Ostap::Math::sub ( v1 , v2 ) ;
     *
     *  @endcode 
     *  
     *  @param v1 (input/output) LorentzVector to be updated 
     *  @param v2 (input) Linear Algebra vector, to be subtracted from Vector 
     *  @return the updated LorenzVector 
     *  @author Vanya BELYAEV ibelyaev@physics.syr.edu
     *  @date 2006-05-24
     */
    template <class C, class T>
    inline 
    const ROOT::Math::DisplacementVector3D<C>& 
    sub 
    ( ROOT::Math::DisplacementVector3D<C>& v1 , const ROOT::Math::SVector<T,3> & v2 ) 
    { return v1 -= ROOT::Math::DisplacementVector3D<C>
        ( v2 ( 0 )  , v2 ( 1 )  , v2 ( 2 ) ) ; }
    // ========================================================================
    /** decrement  LorentzVector with 4-component linear vector 
     *  
     *  @code 
     *
     *  Ostap::LorentzVector v1 = ... ;
     *  const Ostap::Vector4       v2 = ... ;
     *  
     *  // update Lorentz vector with LA vector:
     *  Ostap::Math::sub ( v1 , v2 ) ;
     *
     *  @endcode 
     *  
     *  @param v1 (input/output) LorentzVector to be updated 
     *  @param v2 (input) Linear Algebra vector, to be subtracted from  LorentzVector 
     *  @return the updated LorenzVector 
     *  @author Vanya BELYAEV ibelyaev@physics.syr.edu
     *  @date 2006-05-24
     */
    template <class C, class T>
    inline 
    const ROOT::Math::LorentzVector<C>& 
    sub 
    ( ROOT::Math::LorentzVector<C>& v1 , const ROOT::Math::SVector<T,4> & v2 ) 
    { return v1 -= ROOT::Math::LorentzVector<C>
        ( v2 ( 0 )  , v2 ( 1 )  , v2 ( 2 ) , v2 ( 3 ) ) ; }
    // ========================================================================    
    /** increment the symmetric matrix with "symmetrized" part of other matrix
     * 
     *  @code 
     * 
     *  Ostap::SymMatrix3x3 matrix = ... ;
     *  const Ostap::Matrix3x3    other  = ... ;
     * 
     *  // update "matrix" with the upper triangular part of "other"
     *  Ostap::Math::sub ( matrix , other ) ;
     *
     *  @endcode 
     * 
     *  @param matrix symmetric matrix to be updated 
     *  @param other matrix, upper traingle is used for updating
     *  @return the updated symmetric matrix 
     *  @author Vanya BELYAEV ibelyaev@physics.syr.edu
     *  @date 2006-05-24
     */
    template <class T1,class T2,unsigned int D, class R>
    inline 
    ROOT::Math::SMatrix<T1,D,D,ROOT::Math::MatRepSym<T1,D> >& 
    add
    ( ROOT::Math::SMatrix<T1,D,D,ROOT::Math::MatRepSym<T1,D> >& matrix , 
      const ROOT::Math::SMatrix<T2,D,D,R>&                      other  ) 
    {
      for   ( unsigned int i = 0 ; i < D ; ++i ) 
      { for ( unsigned int j = i ; j < D ; ++j ) { matrix(i,j) += other(i,j); } }
      return matrix ;
    }
    // ========================================================================
    /** Fill Lorentz vector from 3D displacement vector + Mass
     *     
     *  @code
     *
     *  geoLA(xyz, mass, lav)
     *
     *  @endcode
     *
     *  @param[in]  source Linear Algebra vector3D
     *  @param[in]  mass   Mass
     *  @param[out] dest   Lorentz Vector 
     *  @return Lorentz Vector
     *  @author Sean BRISBANE sean.brisbane@cern.ch
     *  @date 2007-11-27
     */
    template<class T, class C, class M>
    inline 
    const ROOT::Math::LorentzVector<C>&
    geo2LA 
    ( const ROOT::Math::SVector<T, 3>& source , 
      const M                          mass   ,
      ROOT::Math::LorentzVector<C>&    dest   ) 
    {
      // first calculate momentum in carthesian coordinates:
      const double p = 1/fabs(source(2)) ;
      const double n = sqrt( 1 + source(0)*source(0)+source(1)*source(1)) ;
      const double pz = p/n          ;
      const double px = pz*source(0) ;
      const double py = pz*source(1) ;
      dest.SetPx ( px ) ;
      dest.SetPy ( py ) ;
      dest.SetPz ( pz ) ;
      dest.SetE  ( std::sqrt ( p*p + mass*mass ) ) ;
      return dest ;
    }
    // ========================================================================
    /** Compute the jacobian for the transformation of a covariance matrix
     *  with rows representing track parameters TxTyQop and columns in xyz
     *  into a covariance matrix representing the track parameters in 
     *  PxPyPzE and columns xyz
     *
     * @code 
     *
     * ROOT::Math::SMatrix<T, 3 ,3 ,R> covTxTyQoP_xyz = ....  ;
     * ROOT::Math::SVector<C,3>& particleMomentum=...;
     * ROOT::Math::SMatrix<T, 4 , 3, R > covPxPyPzE_xyz;
     * massOfParticle = ...;
     *
     * ROOT::Math::SMatrix<T,4,3,R> Jacob;
     * JacobdP4dMom (particleMomentum, massOfParticle, Jacob) ;
     * covPxPyPzE_xyz = Jacob * covTxTyQoP_xyz;
     *
     * @endcode 
     *
     * @param mom     (input)  the txtyqop vector of track/particle momenta
     * @param mass    (input)  the particle mass
     * @param J       (output) the Jacobian for the transformation
     * @author Sean BRISBANE sean.brisbane@cern.ch
     * @date 2007-11-27
     */
    template <class T,class R, class M >
    inline void  JacobdP4dMom 
    ( const ROOT::Math::SVector<T,3>& mom  ,
      const M                         mass , 
      ROOT::Math::SMatrix<R,4,3>&     J    )
    {
      double tx = mom(0) ;
      double ty = mom(1) ;
      double qop = mom(2) ;
      double p  = 1/std::abs(qop) ;
      double n2 = 1 + tx*tx + ty*ty ;
      double n  = std::sqrt(n2) ;
      double n3 = n2*n ;
      double px = p*tx/n ;
      double py = p*ty/n ;
      double pz = p/n ;
      double E = std::hypot(p,mass) ;
      
      J(0,0) = p * (1+ty*ty)/n3 ; // dpx/dtx
      J(0,1) = p * tx * -ty/n3  ; // dpx/dty
      J(0,2) = -px/qop ;          // dpx/dqop

      J(1,0) = p * ty * -tx/n3  ; // dpy/dtx
      J(1,1) = p * (1+tx*tx)/n3 ; // dpy/dty
      J(1,2) = -py/qop ;          // dpy/dqop

      J(2,0) = pz * -tx/n2 ;      // dpz/dtx
      J(2,1) = pz * -ty/n2 ;      // dpz/dtx
      J(2,2) = -pz/qop ;          // dpz/dqop
      
      J(3,0) = 0.0          ;     // dE/dtx 
      J(3,1) = 0.0          ;     // dE/dty 
      J(3,2) = p/E * -p/qop ;     // dE/dqop
      
      return ; 
    }
    // ========================================================================
  } // end of namespace Ostap::Math
  // ==========================================================================
} // end of namespace Ostap
// ============================================================================
namespace Ostap
{
  namespace Math 
  {
    // ========================================================================
    /** @namespace Ostap::Math::Operators LHCbMath/MatrixTransforms.h
     *
     *  The helper namespace which contains inline operators for 
     *  various objects from the world of  "Geometry&Kinematics" 
     *  and the object from the world of Linear Algebra 
     *
     *  The existence of thes e operator <b>drastically simplifies</b> 
     *  the code, dealing with kinematical and/or topolofical 
     *  calculations, in particular the implementation of various 
     *  kinematical fitters
     *
     *  To make these operators vizible for the real code one 
     *  needs to use <c>using</c> directive, e.g. 
     *
     *  @code 
     *
     *   StatusCode myFunction ( ... ) 
     *   {
     *       // get access to the useful operators:
     *       using namespace Ostap::Math::Operators ; // NB! 
     *
     *        ... use the operators ... 
     *   } 
     *
     *  @endcode 
     *
     *  <table>
     *
     *  <tr>
     *     <td align="center"> <b> The First Operand  Type </b> </td> 
     *     <td align="center"> <b> The Operator            </b> </td>
     *     <td align="center"> <b> The Second Operand Type </b> </td> 
     *     <td align="center"> <b> The Result Type         </b> </td> 
     *  </tr>
     *
     *  <tr>
     *     <td>ROOT::Math::PositionVector3D<C>  </td>
     *     <td align="center"><b>+</b></td>
     *     <td>ROOT::Math::SVector&<T,3>        </td>
     *     <td>ROOT::Math::PositionVector3D<C>  </td>
     *  </tr> 
     *  <tr>
     *     <td>ROOT::Math::PositionVector3D<C>  </td>
     *     <td align="center"><b>+</b></td>
     *     <td>ROOT::Math::VecExpr<B,T,3>;      </td>
     *     <td>ROOT::Math::PositionVector3D<C>  </td>
     *  </tr> 
     *  <tr>
     *     <td>ROOT::Math::PositionVector3D<C> </td>
     *     <td align="center"><b>-</b></td>
     *     <td>ROOT::Math::SVector<T,3>        </td>
     *     <td>ROOT::Math::PositionVector3D<C> </td>
     *  </tr> 
     *  <tr>
     *     <td>ROOT::Math::PositionVector3D<C> </td>
     *     <td align="center"><b>-</b></td>
     *     <td>ROOT::Math::VecExpr<B,T,3>      </td>
     *     <td>ROOT::Math::PositionVector3D<C> </td>
     *  </tr> 
     *
     *  <tr>
     *     <td>ROOT::Math::DisplacementVector3D<C> </td>
     *     <td align="center"><b>+</b></td>
     *     <td>ROOT::Math::SVector<T,3>            </td>
     *     <td>ROOT::Math::DisplacementVector3D<C> </td>
     *  </tr> 
     *  <tr>
     *     <td>ROOT::Math::DisplacementVector3D<C> </td>
     *     <td align="center"><b>+</b></td>
     *     <td>ROOT::Math::VecExpr<B,T,3>         </td>
     *     <td>ROOT::Math::DisplacementVector3D<C></td>
     *  </tr> 
     *  <tr>
     *     <td>ROOT::Math::DisplacementVector3D<C></td>
     *     <td align="center"><b>-</b></td>
     *     <td>ROOT::Math::SVector<T,3>       </td>
     *     <td>ROOT::Math::DisplacementVector3D<C></td>
     *  </tr> 
     *  <tr>
     *     <td>ROOT::Math::DisplacementVector3D<C></td>
     *     <td align="center"><b>-</b></td>
     *     <td>ROOT::Math::VecExpr<B,T,3>         </td>
     *     <td>ROOT::Math::DisplacementVector3D<C></td>
     *  </tr> 
     *
     *  <tr>
     *     <td>ROOT::Math::LorentzVector<C> </td>
     *     <td align="center"><b>+</b></td>
     *     <td>ROOT::Math::SVector<T,4>     </td>
     *     <td>ROOT::Math::LorentzVector<C> </td>
     *  </tr> 
     *  <tr>
     *     <td>ROOT::Math::LorentzVector<C> </td>
     *     <td align="center"><b>+</b></td>
     *     <td>ROOT::Math::VecExpr<B,T,4>   </td>
     *     <td>ROOT::Math::LorentzVector<C> </td>
     *  </tr> 
     *  <tr>
     *     <td>ROOT::Math::LorentzVector<C> </td>
     *     <td align="center"><b>-</b></td>
     *     <td>ROOT::Math::SVector<T,4>     </td>
     *     <td>ROOT::Math::LorentzVector<C> </td>
     *  </tr> 
     *  <tr>
     *     <td>ROOT::Math::LorentzVector<C> </td>
     *     <td align="center"><b>-</b></td>
     *     <td>ROOT::Math::VecExpr<B,T,4>   </td>
     *     <td>ROOT::Math::LorentzVector<C> </td>
     *  </tr> 
     *
     *  <tr>
     *     <td>ROOT::Math::SVector<T,3>          </td>
     *     <td align="center"><b>+</b></td>
     *     <td>ROOT::Math::PositionVector3D<C>   </td>
     *     <td>ROOT::Math::SVector<T,3>          </td>
     *  </tr> 
     *  <tr>
     *     <td>ROOT::Math::VecExpr<B,T,3>        </td>
     *     <td align="center"><b>+</b></td>
     *     <td>ROOT::Math::PositionVector3D<C>   </td>
     *     <td>ROOT::Math::SVector<T,3>          </td>
     *  </tr> 
     *  <tr>
     *     <td>ROOT::Math::SVector<T,3>          </td>
     *     <td align="center"><b>-</b></td>
     *     <td>ROOT::Math::PositionVector3D<C>   </td>
     *     <td>ROOT::Math::SVector<T,3>          </td>
     *  </tr> 
     *  <tr>
     *     <td>ROOT::Math::VecExpr<B,T,3>        </td>
     *     <td align="center"><b>-</b></td>
     *     <td>ROOT::Math::PositionVector3D<C>   </td>
     *     <td>ROOT::Math::SVector<T,3>          </td>
     *  </tr> 
     *
     *
     *  <tr>
     *     <td>ROOT::Math::SVector<T,3>           </td>
     *     <td align="center"><b>+</b></td>
     *     <td>ROOT::Math::DisplacementVector3D<C></td>
     *     <td>ROOT::Math::SVector<T,3>           </td>
     *  </tr> 
     *  <tr>
     *     <td>ROOT::Math::VecExpr<B,T,3>         </td>
     *     <td align="center"><b>+</b></td>
     *     <td>ROOT::Math::DisplacementVector3D<C></td>
     *     <td>ROOT::Math::SVector<T,3>           </td>
     *  </tr> 
     *  <tr>
     *     <td>ROOT::Math::SVector<T,3>           </td>
     *     <td align="center"><b>-</b></td>
     *     <td>ROOT::Math::DisplacementVector3D<C></td>
     *     <td>ROOT::Math::SVector<T,3>           </td>
     *  </tr> 
     *  <tr>
     *     <td>ROOT::Math::VecExpr<B,T,3>         </td>
     *     <td align="center"><b>-</b></td>
     *     <td>ROOT::Math::DisplacementVector3D<C></td>
     *     <td>ROOT::Math::SVector<T,3>           </td>
     *  </tr> 
     *
     *  <tr>
     *     <td>ROOT::Math::SVector<T,4>       </td>
     *     <td align="center"><b>+</b></td>
     *     <td>ROOT::Math::LorentzVector<C>   </td>
     *     <td>ROOT::Math::SVector<T,4>       </td>
     *  </tr> 
     *  <tr>
     *     <td>ROOT::Math::VecExpr<B,T,4>     </td>
     *     <td align="center"><b>+</b></td>
     *     <td>ROOT::Math::LorentzVector<C>   </td>
     *     <td>ROOT::Math::SVector<T,4>       </td>
     *  </tr> 
     *  <tr>
     *     <td>ROOT::Math::SVector<T,4>       </td>
     *     <td align="center"><b>-</b></td>
     *     <td>ROOT::Math::LorentzVector<C>   </td>
     *     <td>ROOT::Math::SVector<T,4>       </td>
     *  </tr> 
     *  <tr>
     *     <td>ROOT::Math::VecExpr<B,T,4>     </td>
     *     <td align="center"><b>-</b></td>
     *     <td>ROOT::Math::LorentzVector<C>   </td>
     *     <td>ROOT::Math::SVector<T,4>       </td>
     *  </tr> 
     *
     *  <tr>
     *     <td>ROOT::Math::SMatrix<T,D,3,R>    </td>
     *     <td align="center"><b>*</b></td>
     *     <td>ROOT::Math::PositionVector3D<C> </td>
     *     <td>ROOT::Math::SVector<T,D>        </td>
     *  </tr> 
     *  <tr>
     *     <td>ROOT::Math::Expr<B,T,D,3,R>     </td>
     *     <td align="center"><b>*</b></td>
     *     <td>ROOT::Math::PositionVector3D<C> </td>
     *     <td>ROOT::Math::SVector<T,D>        </td>
     *  </tr> 
     *  <tr>
     *     <td>ROOT::Math::PositionVector3D<C> </td>
     *     <td align="center"><b>*</b></td>
     *     <td>ROOT::Math::SMatrix<T,3,D,R>    </td>
     *     <td>ROOT::Math::SVector<T,D>        </td>
     *  </tr> 
     *  <tr>
     *     <td>ROOT::Math::PositionVector3D<C> </td>
     *     <td align="center"><b>*</b></td>
     *     <td>ROOT::Math::Expr<B,T,3,D,R>     </td>
     *     <td>ROOT::Math::SVector<T,D>        </td>
     *  </tr> 
     *
     *  <tr>
     *     <td>ROOT::Math::SMatrix<T,D,3,R>        </td>
     *     <td align="center"><b>*</b></td>
     *     <td>ROOT::Math::DisplacementVector3D<C> </td>
     *     <td>ROOT::Math::SVector<T,D>            </td>
     *  </tr> 
     *  <tr>
     *     <td>ROOT::Math::Expr<B,T,D,3,R>         </td>
     *     <td align="center"><b>*</b></td>
     *     <td>ROOT::Math::DisplacementVector3D<C> </td>
     *     <td>ROOT::Math::SVector<T,D>            </td>
     *  </tr> 
     *  <tr>
     *     <td>ROOT::Math::DisplacementVector3D<C> </td>
     *     <td align="center"><b>*</b></td>
     *     <td>ROOT::Math::SMatrix<T,3,D,R>        </td>
     *     <td>ROOT::Math::SVector<T,D>            </td>
     *  </tr> 
     *  <tr>
     *     <td>ROOT::Math::DisplacementVector3D<C> </td>
     *     <td align="center"><b>*</b></td>
     *     <td>ROOT::Math::Expr<B,T,3,D,R>         </td>
     *     <td>ROOT::Math::SVector<T,D>            </td>
     *  </tr> 
     *
     *  <tr>
     *     <td>ROOT::Math::SMatrix<T,D,4,R>        </td>
     *     <td align="center"><b>*</b></td>
     *     <td>ROOT::Math::LorentzVector<C>        </td>
     *     <td>ROOT::Math::SVector<T,D>            </td>
     *  </tr> 
     *  <tr>
     *     <td>ROOT::Math::Expr<B,T,D,4,R>         </td>
     *     <td align="center"><b>*</b></td>
     *     <td>ROOT::Math::LorentzVector<C>        </td>
     *     <td>ROOT::Math::SVector<T,D>            </td>
     *  </tr> 
     *  <tr>
     *     <td>ROOT::Math::LorentzVector<C>        </td>
     *     <td align="center"><b>*</b></td>
     *     <td>ROOT::Math::SMatrix<T,4,D,R>        </td>
     *     <td>ROOT::Math::SVector<T,D>            </td>
     *  </tr> 
     *  <tr>
     *     <td>ROOT::Math::LorentzVector<C>        </td>
     *     <td align="center"><b>*</b></td>
     *     <td>ROOT::Math::Expr<B,T,4,D,R>         </td>
     *     <td>ROOT::Math::SVector<T,D>            </td>
     *  </tr> 
     *
     *   </table>
     *
     *  @author Vanya BELYAEV ibelyaev@itep.ru
     *  @date 2008-02-14
     */
    namespace Operators 
    {
      // ======================================================================
      /** addition of 3D-vector and the linear algebra vector 
       *
       *  @code
       * 
       *  const Ostap::XYZPoint& p1 = ... ;
       *  const Ostap::Vector3&  v2 = ... ;
       *
       *  Ostap::XYZPoint p = p1 + v2 ; 
       *
       *  @endcode 
       *
       *  @param  p1 the position vector (point) 
       *  @param  v2 the linear algebra vector 
       *  @return the effective position vector 
       *
       *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
       *  @date   2008-03-05
       */
      template <class C, class T>
      inline 
      ROOT::Math::PositionVector3D<C>
      operator+ 
      ( const ROOT::Math::PositionVector3D<C>& p1 ,
        const ROOT::Math::SVector<T,3>&        v2 )
      {
        ROOT::Math::PositionVector3D<C> result  ;
        result.SetXYZ(  p1 . X () + v2 ( 0 ) , 
                        p1 . Y () + v2 ( 1 ) ,
                        p1 . Z () + v2 ( 2 ) ) ;
        return result ;
      }
      // ======================================================================
      template <class C, class B, class T>
      inline 
      ROOT::Math::PositionVector3D<C>
      operator+ 
      ( const ROOT::Math::PositionVector3D<C>& p1 ,
        const ROOT::Math::VecExpr<B,T,3>&      v2 )
      {
        ROOT::Math::PositionVector3D<C> result  ;
        result.SetXYZ(  p1 . X () + v2 ( 0 ) , 
                        p1 . Y () + v2 ( 1 ) ,
                        p1 . Z () + v2 ( 2 ) ) ;
        return result ;
      }
      // ========================================================================
      /** addition of 3D-vector and the linear algebra vector 
       *
       *  @code
       * 
       *  const Ostap::XYZVector& p1 = ... ;
       *  const Ostap::Vector3&   v2 = ... ;
       *
       *  Ostap::XYZVector p = p1 + v2 ; 
       *
       *  @endcode 
       *
       *  @param  p1 the displacement vector  
       *  @param  v2 the linear algebra vector 
       *  @return the effective displacement vector 
       *
       *  @author Vanya BELYAEV ibelyaev@itep.ru
       *  @date   2008-02-14
       */
      template <class C, class T>
      inline 
      ROOT::Math::DisplacementVector3D<C>
      operator+ 
      ( const ROOT::Math::DisplacementVector3D<C>& p1 ,
        const ROOT::Math::SVector<T,3>&            v2 )
      {
        ROOT::Math::DisplacementVector3D<C> result  ;
        result.SetXYZ(  p1 . X () + v2 ( 0 ) , 
                        p1 . Y () + v2 ( 1 ) ,
                        p1 . Z () + v2 ( 2 ) ) ;
        return result ;
      }
      // ======================================================================
      template <class C, class B, class T>
      inline 
      ROOT::Math::DisplacementVector3D<C>
      operator+ 
      ( const ROOT::Math::DisplacementVector3D<C>& p1 ,
        const ROOT::Math::VecExpr<B,T,3>&          v2 )
      {
        ROOT::Math::DisplacementVector3D<C> result  ;
        result.SetXYZ(  p1 . X () + v2 ( 0 ) , 
                        p1 . Y () + v2 ( 1 ) ,
                        p1 . Z () + v2 ( 2 ) ) ;
        return result ;
      }
      // ======================================================================
      /** addition of Lorentz vector and the linear algebra vector 
       *
       *  @code
       * 
       *  const Ostap::LorentzVector& p1 = ... ;
       *  const Ostap::Vector4&       v2 = ... ;
       *
       *  Ostap::LorentzVector p = p1 + v2 ; 
       *
       *  @endcode 
       *
       *  @param  p1 Lorentz vector  
       *  @param  v2 the linear algebra vector 
       *  @return the effective Lorentz vector 
       *
       *  @author Vanya BELYAEV ibelyaev@itep.ru
       *  @date   2008-02-14
       */
      template <class C, class T>
      inline 
      ROOT::Math::LorentzVector<C>
      operator+ 
      ( const ROOT::Math::LorentzVector<C>& p1 ,
        const ROOT::Math::SVector<T,4>&     v2 )
      {
        ROOT::Math::LorentzVector<C> result  ;
        result.SetXYZT 
          (  p1 . Px () + v2 ( 0 ) , 
             p1 . Py () + v2 ( 1 ) ,
             p1 . Pz () + v2 ( 2 ) ,
             p1 . E  () + v2 ( 3 ) ) ;
        return result ;
      }
      // ======================================================================
      template <class C, class B , class T>
      inline 
      ROOT::Math::LorentzVector<C>
      operator+ 
      ( const ROOT::Math::LorentzVector<C>& p1 ,
        const ROOT::Math::VecExpr<B,T,4>&   v2 )
      {
        ROOT::Math::LorentzVector<C> result  ;
        result.SetXYZT 
          (  p1 . Px () + v2 ( 0 ) , 
             p1 . Py () + v2 ( 1 ) ,
             p1 . Pz () + v2 ( 2 ) ,
             p1 . E  () + v2 ( 3 ) ) ;
        return result ;
      }
      // ======================================================================
      /** addition of 3D-vector and the linear algebra vector 
       *
       *  @code
       * 
       *  const Ostap::Vector3&  v2 = ... ;
       *  const Ostap::XYZPoint& p1 = ... ;
       *
       *  Ostap::Vector3 p = v2 + p1 ; 
       *
       *  @endcode 
       *
       *  @param  v2 the linear algebra vector 
       *  @param  p1 the position vector (point) 
       *  @return the effective linear algebra vector 
       *
       *  @author Vanya BELYAEV ibelyaev@itep.ru
       *  @date   2008-02-14
       */
      template <class C, class T>
      inline 
      ROOT::Math::SVector<T,3>
      operator+ 
      ( const ROOT::Math::SVector<T,3>&        v2 ,
        const ROOT::Math::PositionVector3D<C>& p1 )
      {
        ROOT::Math::SVector<T,3> result ( v2 ) ;
        result ( 0 ) += p1 . X () ;
        result ( 1 ) += p1 . Y () ;
        result ( 2 ) += p1 . Z () ;
        return result ;
      }
      // ======================================================================
      template <class C, class B, class T>
      inline 
      ROOT::Math::SVector<T,3>
      operator+ 
      ( const ROOT::Math::VecExpr<B,T,3>&      v2 ,
        const ROOT::Math::PositionVector3D<C>& p1 )
      {
        ROOT::Math::SVector<T,3> result = v2 ;
        result ( 0 ) += p1 . X () ;
        result ( 1 ) += p1 . Y () ;
        result ( 2 ) += p1 . Z () ;
        return result ;
      }
      // ======================================================================
      /** addition of 3D-vector and the linear algebra vector 
       *
       *  @code
       * 
       *  const Ostap::Vector3&   v2 = ... ;
       *  const Ostap::XYZVector& p1 = ... ;
       *
       *  Ostap::Vector3 p = v2 + p1 ; 
       *
       *  @endcode 
       *
       *  @param  v2 the linear algebra vector 
       *  @param  p1 the displacement vector 
       *  @return the effective linear algebra vector 
       *
       *  @author Vanya BELYAEV ibelyaev@itep.ru
       *  @date   2008-02-14
       */
      template <class C, class T>
      inline 
      ROOT::Math::SVector<T,3>
      operator+ 
      ( const ROOT::Math::SVector<T,3>&            v2 ,
        const ROOT::Math::DisplacementVector3D<C>& p1 )
      {
        ROOT::Math::SVector<T,3> result ( v2 ) ;
        result ( 0 ) += p1 . X () ;
        result ( 1 ) += p1 . Y () ;
        result ( 2 ) += p1 . Z () ;
        return result ;
      }
      // ======================================================================
      template <class C, class B, class T>
      inline 
      ROOT::Math::SVector<T,3>
      operator+ 
      ( const ROOT::Math::VecExpr<B,T,3>&          v2 ,
        const ROOT::Math::DisplacementVector3D<C>& p1 )
      {
        ROOT::Math::SVector<T,3> result = v2 ;
        result ( 0 ) += p1 . X () ;
        result ( 1 ) += p1 . Y () ;
        result ( 2 ) += p1 . Z () ;
        return result ;
      }
      // ======================================================================
      /** addition of Lorentz vector and the linear algebra vector 
       *
       *  @code
       * 
       *  const Ostap::Vector4&        v2 = ... ;
       *  const Ostap::LorentzVector&  p1 = ... ;
       *
       *  Ostap::Vector4 p = v2 + p1 ; 
       *
       *  @endcode 
       *
       *  @param  v2 the linear algebra vector 
       *  @param  p1 Lorentz vector 
       *  @param  v2 the linear algebra vector 
       *  @return the effective linear algebra vector 
       *
       *  @author Vanya BELYAEV ibelyaev@itep.ru
       *  @date   2008-02-14
       */
      template <class C, class T>
      inline 
      ROOT::Math::SVector<T,4>
      operator+ 
      ( const ROOT::Math::SVector<T,4>&     v2 ,
        const ROOT::Math::LorentzVector<C>& p1 )
      {
        ROOT::Math::SVector<T,4> result ( v2 ) ;
        result ( 0 ) += p1 . Px () ;
        result ( 1 ) += p1 . Py () ;
        result ( 2 ) += p1 . Pz () ;
        result ( 3 ) += p1 . E  () ;
        return result ;
      }
      // ======================================================================
      template <class C, class B, class T>
      inline 
      ROOT::Math::SVector<T,4>
      operator+ 
      ( const ROOT::Math::VecExpr<B,T,4>&   v2 ,
        const ROOT::Math::LorentzVector<C>& p1 )
      {
        ROOT::Math::SVector<T,4> result = v2  ;
        result ( 0 ) += p1 . Px () ;
        result ( 1 ) += p1 . Py () ;
        result ( 2 ) += p1 . Pz () ;
        result ( 3 ) += p1 . E  () ;
        return result ;
      }
      // ======================================================================
      /** subtraction of 3D-vector and the linear algebra vector 
       *
       *  @code
       * 
       *  const Ostap::XYZPoint& p1 = ... ;
       *  const Ostap::Vector3&  v2 = ... ;
       *
       *  Ostap::XYZPoint p = p1 - v2 ; 
       *
       *  @endcode 
       *
       *  @param  p1 the position vector (point) 
       *  @param  v2 the linear algebra vector 
       *  @return the effective position vector 
       *
       *  @author Vanya BELYAEV ibelyaev@itep.ru
       *  @date   2008-02-14
       */
      template <class C, class T>
      inline 
      ROOT::Math::PositionVector3D<C>
      operator- 
      ( const ROOT::Math::PositionVector3D<C>& p1 ,
        const ROOT::Math::SVector<T,3>&        v2 )
      {
        ROOT::Math::PositionVector3D<C> result  ;
        result.SetXYZ
          (  p1 . X () - v2 ( 0 ) , 
             p1 . Y () - v2 ( 1 ) ,
             p1 . Z () - v2 ( 2 ) ) ;
        return result ;
      }
     // ======================================================================
      template <class C, class B , class T>
      inline 
      ROOT::Math::PositionVector3D<C>
      operator- 
      ( const ROOT::Math::PositionVector3D<C>& p1 ,
        const ROOT::Math::VecExpr<B,T,3>&      v2 )
      {
        ROOT::Math::PositionVector3D<C> result  ;
        result.SetXYZ
          (  p1 . X () - v2 ( 0 ) , 
             p1 . Y () - v2 ( 1 ) ,
             p1 . Z () - v2 ( 2 ) ) ;
        return result ;
      }      
      // ========================================================================
      /** subtraction of 3D-vector and the linear algebra vector 
       *
       *  @code
       * 
       *  const Ostap::XYZVector& p1 = ... ;
       *  const Ostap::Vector3&   v2 = ... ;
       *
       *  Ostap::XYZVector p = p1 - v2 ; 
       *
       *  @endcode 
       *
       *  @param  p1 the displacement vector  
       *  @param  v2 the linear algebra vector 
       *  @return the effective displacement vector 
       *
       *  @author Vanya BELYAEV ibelyaev@itep.ru
       *  @date   2008-02-14
       */
      template <class C, class T>
      inline 
      ROOT::Math::DisplacementVector3D<C>
      operator- 
      ( const ROOT::Math::DisplacementVector3D<C>& p1 ,
        const ROOT::Math::SVector<T,3>&            v2 )
      {
        ROOT::Math::DisplacementVector3D<C> result  ;
        result.SetXYZ
          (  p1 . X () - v2 ( 0 ) , 
             p1 . Y () - v2 ( 1 ) ,
             p1 . Z () - v2 ( 2 ) ) ;
        return result ;
      }
      // ======================================================================
      template <class C, class B, class T>
      inline 
      ROOT::Math::DisplacementVector3D<C>
      operator- 
      ( const ROOT::Math::DisplacementVector3D<C>& p1 ,
        const ROOT::Math::VecExpr<B,T,3>&          v2 )
      {
        ROOT::Math::DisplacementVector3D<C> result  ;
        result.SetXYZ
          (  p1 . X () - v2 ( 0 ) , 
             p1 . Y () - v2 ( 1 ) ,
             p1 . Z () - v2 ( 2 ) ) ;
        return result ;
      }      
      // ======================================================================
      /** subtraction of Lorentz vector and the linear algebra vector 
       *
       *  @code
       * 
       *  const Ostap::LorentzVector& p1 = ... ;
       *  const Ostap::Vector4&       v2 = ... ;
       *
       *  Ostap::LorentzVector p = p1 - v2 ; 
       *
       *  @endcode 
       *
       *  @param  p1 Lorentz vector  
       *  @param  v2 the linear algebra vector 
       *  @return the effective Lorentz vector 
       *
       *  @author Vanya BELYAEV ibelyaev@itep.ru
       *  @date   2008-02-14
       */
      template <class C, class T>
      inline 
      ROOT::Math::LorentzVector<C>
      operator- 
      ( const ROOT::Math::LorentzVector<C>& p1 ,
        const ROOT::Math::SVector<T,4>&     v2 )
      {
        ROOT::Math::LorentzVector<C> result  ;
        result.SetXYZT 
          (  p1 . Px () - v2 ( 0 ) , 
             p1 . Py () - v2 ( 1 ) ,
             p1 . Pz () - v2 ( 2 ) ,
             p1 . E  () - v2 ( 3 ) ) ;
        return result ;
      }
      // ======================================================================
      template <class C, class B, class T>
      inline 
      ROOT::Math::LorentzVector<C>
      operator- 
      ( const ROOT::Math::LorentzVector<C>& p1 ,
        const ROOT::Math::VecExpr<B,T,4>&   v2 )
      {
        ROOT::Math::LorentzVector<C> result  ;
        result.SetXYZT 
          (  p1 . Px () - v2 ( 0 ) , 
             p1 . Py () - v2 ( 1 ) ,
             p1 . Pz () - v2 ( 2 ) ,
             p1 . E  () - v2 ( 3 ) ) ;
        return result ;
      }
      // ======================================================================      
      /** subtract the Lorentz Vector from the Linear Algebra -vector 
       *
       *  @code
       * 
       *   const Ostap::Vector4&       vct1 = ... ;
       *   const Ostap::LorentzVector& vct2 = ... ;
       *  
       *   std::cout << " Delta is " << vct1-vct2 << std::endl ;
       *
       *  @endcode 
       *
       *  @author Vanya BEYAEV Ivan.Belyaev@nikhef.nl
       *  @date 2008-03-03
       */
      template <class C, class T> 
      inline 
      ROOT::Math::SVector<T,4>
      operator- 
      ( const ROOT::Math::SVector<T,4>&     v1 , 
        const ROOT::Math::LorentzVector<C>& v2 ) 
      {
        ROOT::Math::SVector<T,4> result = v1 ;
        result ( 0 ) -= v2 . Px () ;
        result ( 1 ) -= v2 . Py () ;
        result ( 2 ) -= v2 . Pz () ;
        result ( 3 ) -= v2 . E  () ;
        return result ;
      }
      // ======================================================================
      template <class C, class B, class T> 
      inline 
      ROOT::Math::SVector<T,4>
      operator- 
      ( const ROOT::Math::VecExpr<B,T,4>&   v1 , 
        const ROOT::Math::LorentzVector<C>& v2 ) 
      {
        ROOT::Math::SVector<T,4> result = v1 ;
        result ( 0 ) -= v2 . Px () ;
        result ( 1 ) -= v2 . Py () ;
        result ( 2 ) -= v2 . Pz () ;
        result ( 3 ) -= v2 . E  () ;
        return result ;
      }
      // ======================================================================      
      /** subtract the 3D Vector from the Linear Algebra -vector 
       *
       *  @code
       * 
       *   const Ostap::Vector3&  vct1 = ... ;
       *   const Ostap::XYZPoint& vct2 = ... ;
       *  
       *   std::cout << " Delta is " << vct1-vct2 << std::endl ;
       *
       *  @endcode 
       *
       *  @author Vanya BEYAEV Ivan.Belyaev@nikhef.nl
       *  @date 2008-03-03
       */
      template <class C, class T> 
      inline 
      ROOT::Math::SVector<T,3>
      operator- 
      ( const ROOT::Math::SVector<T,3>&        v1 , 
        const ROOT::Math::PositionVector3D<C>& v2 ) 
      {
        ROOT::Math::SVector<T,3> result = v1 ;
        result ( 0 ) -= v2 . X () ;
        result ( 1 ) -= v2 . Y () ;
        result ( 2 ) -= v2 . Z () ;
        return result ;
      }
      // ======================================================================
      template <class C, class B, class T> 
      inline 
      ROOT::Math::SVector<T,3>
      operator- 
      ( const ROOT::Math::VecExpr<B,T,3>&      v1 , 
        const ROOT::Math::PositionVector3D<C>& v2 ) 
      {
        ROOT::Math::SVector<T,3> result = v1 ;
        result ( 0 ) -= v2 . X () ;
        result ( 1 ) -= v2 . Y () ;
        result ( 2 ) -= v2 . Z () ;
        return result ;
      }
      // ======================================================================      
      /** subtract the 3D Vector from the Linear Algebra -vector 
       *    
       *  @code
       * 
       *   const Ostap::Vector3&   vct1 = ... ;
       *   const Ostap::XYZVector& vct2 = ... ;
       *  
       *   std::cout << " Delta is " << vct1-vct2 << std::endl ;
       *
       *  @endcode 
       *
       *  @author Vanya BEYAEV Ivan.Belyaev@nikhef.nl
       *  @date 2008-03-03
       */
      template <class C, class T> 
      inline 
      ROOT::Math::SVector<T,3>
      operator- 
      ( const ROOT::Math::SVector<T,3>&            v1 , 
        const ROOT::Math::DisplacementVector3D<C>& v2 ) 
      {
        ROOT::Math::SVector<T,3> result = v1 ;
        result ( 0 ) -= v2 . X () ;
        result ( 1 ) -= v2 . Y () ;
        result ( 2 ) -= v2 . Z () ;
        return result ;
      }
      // ======================================================================
      template <class C, class B, class T> 
      inline 
      ROOT::Math::SVector<T,3>
      operator- 
      ( const ROOT::Math::VecExpr<B,T,3>&          v1 , 
        const ROOT::Math::DisplacementVector3D<C>& v2 ) 
      {
        ROOT::Math::SVector<T,3> result = v1 ;
        result ( 0 ) -= v2 . X () ;
        result ( 1 ) -= v2 . Y () ;
        result ( 2 ) -= v2 . Z () ;
        return result ;
      }
      // ======================================================================      
      /** multiply the matrix and the Lorenz vector 
       *
       *  @code 
       *  
       *  const Ostap::SymMatrix4x4  mtrx  = ... ;
       *  const Ostap::LorentzVector vect  = ... ;
       *
       *  const Ostap::Vector4 resut = mrtx * vect ;
       *
       *  @endcode 
       *
       *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
       *  @date 2008-03-03
       */
      template <class T, class C, class R, unsigned int D>
      inline 
      ROOT::Math::SVector<T,D>
      operator* 
      ( const ROOT::Math::SMatrix<T,D,4,R>& mrtx , 
        const ROOT::Math::LorentzVector<C>& vect ) 
      {
        const ROOT::Math::SVector<T,4> vct 
          ( vect . Px () , vect . Py () , vect . Pz () , vect.E () ) ;
        return mrtx * vct ;
      } 
      // ======================================================================
      template <class T, class C, class B, class R, unsigned int D>
      inline 
      ROOT::Math::SVector<T,D>
      operator* 
      ( const ROOT::Math::Expr<B,T,D,4,R>&  mtrx , 
        const ROOT::Math::LorentzVector<C>& vect ) 
      {
        const ROOT::Math::SVector<T,4> vct 
          ( vect . Px () , vect . Py () , vect . Pz () , vect.E () ) ;
        return mtrx * vct ;
      }
      // ======================================================================      
      /** multiply the matrix and the Lorenz vector 
       *
       *  @code 
       *  
       *  const Ostap::Matrix4x5  mtrx  = ... ;
       *  const Ostap::LorentzVector vect  = ... ;
       *
       *  const Ostap::Vector5 result = vect * mrtx ;
       *
       *  @endcode 
       *
       *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
       *  @date 2008-03-03
       */
      template <class T, class C, class R, unsigned int D>
      inline 
      ROOT::Math::SVector<T,D>
      operator* 
      ( const ROOT::Math::LorentzVector<C>& vect , 
        const ROOT::Math::SMatrix<T,4,D,R>& mtrx )
      {
        const ROOT::Math::SVector<T,4> vct 
          ( vect . Px () , vect . Py () , vect . Pz () , vect.E () ) ;
        return vct * mtrx ;
      }
      // ======================================================================
      template <class T, class C, class B, class R, unsigned int D>
      inline 
      ROOT::Math::SVector<T,D>
      operator* 
      ( const ROOT::Math::LorentzVector<C>& vect , 
        const ROOT::Math::Expr<B,T,4,D,R>&  mtrx )
      {
        const ROOT::Math::SVector<T,4> vct 
          ( vect . Px () , vect . Py () , vect . Pz () , vect.E () ) ;
        return vct * mtrx ;
      }
      // ======================================================================      
      /** multiply the matrix and 3D-vector 
       *
       *  @code 
       *  
       *  const Ostap::SymMatrix3x3  mtrx  = ... ;
       *  const Ostap::XYZVector     vect  = ... ;
       *
       *  const Ostap::Vector3 resut = mrtx * vect ;
       *
       *  @endcode 
       *
       *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
       *  @date 2008-03-03
       */
      template <class T, class C, class R, unsigned int D>
      inline 
      ROOT::Math::SVector<T,D>
      operator* 
      ( const ROOT::Math::SMatrix<T,D,3,R>&        mtrx , 
        const ROOT::Math::DisplacementVector3D<C>& vect ) 
      {
        const ROOT::Math::SVector<T,3> vct 
          ( vect . X () , vect . Y () , vect . Z () ) ;
        return mtrx * vct ;
      } 
      // ======================================================================
      template <class T, class C, class B, class R, unsigned int D>
      inline 
      ROOT::Math::SVector<T,D>
      operator* 
      ( const ROOT::Math::Expr<B,T,D,3,R>&         mtrx , 
        const ROOT::Math::DisplacementVector3D<C>& vect ) 
      {
        const ROOT::Math::SVector<T,3> vct 
          ( vect . X () , vect . Y () , vect . Z () ) ;
        return mtrx * vct ;
      } 
      // ======================================================================      
      /** multiply the matrix and the 3D-vector 
       *
       *  @code 
       *  
       *  const Ostap::SymMatrix3x3  mtrx  = ... ;
       *  const Ostap::XYZVector vect  = ... ;
       *
       *  const Ostap::Vector3 resut = vect * mtrx ;
       *
       *  @endcode 
       *
       *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
       *  @date 2008-03-03
       */
      template <class T, class C, class R, unsigned int D>
      inline 
      ROOT::Math::SVector<T,D>
      operator* 
      ( const ROOT::Math::DisplacementVector3D<C>& vect , 
        const ROOT::Math::SMatrix<T,3,D,R>&        mtrx )
      {
        const ROOT::Math::SVector<T,3> vct 
          ( vect . X () , vect . Y () , vect . Z () ) ;
        return vct * mtrx ;
      }
      // ======================================================================      
      template <class T, class C, class B, class R, unsigned int D>
      inline 
      ROOT::Math::SVector<T,D>
      operator* 
      ( const ROOT::Math::DisplacementVector3D<C>& vect , 
        const ROOT::Math::Expr<B,T,3,D,R>&         mtrx )
      {
        const ROOT::Math::SVector<T,3> vct 
          ( vect . X () , vect . Y () , vect . Z () ) ;
        return vct * mtrx ;
      }
      // ======================================================================      
      /** multiply the matrix and 3D-vector 
       *
       *  @code 
       *  
       *  const Ostap::Matrix3x3     mtrx  = ... ;
       *  const Ostap::XYZPoint      vect  = ... ;
       *
       *  const Ostap::Vector3 result = mrtx * vect ;
       *
       *  @endcode 
       *
       *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
       *  @date 2008-03-03
       */
      template <class T, class C, class R, unsigned int D>
      inline 
      ROOT::Math::SVector<T,D>
      operator* 
      ( const ROOT::Math::SMatrix<T,D,3,R>&    mtrx , 
        const ROOT::Math::PositionVector3D<C>& vect ) 
      {
        const ROOT::Math::SVector<T,3> vct 
          ( vect . X () , vect . Y () , vect . Z () ) ;
        return mtrx * vct ;
      } 
      // ======================================================================      
      template <class T, class C, class B, class R, unsigned int D>
      inline 
      ROOT::Math::SVector<T,D>
      operator* 
      ( const ROOT::Math::Expr<B,T,D,3,R>&     mtrx , 
        const ROOT::Math::PositionVector3D<C>& vect ) 
      {
        const ROOT::Math::SVector<T,3> vct 
          ( vect . X () , vect . Y () , vect . X () ) ;
        return mtrx * vct ;
      } 
      // ======================================================================      
      /** multiply the matrix and the 3D-vector 
       *
       *  @code 
       *  
       *  const Ostap::SymMatrix3x3  mtrx  = ... ;
       *  const Ostap::XYZPoint vect  = ... ;
       *
       *  const Ostap::Vector3 resut = vect * mtrx ;
       *
       *  @endcode 
       *
       *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
       *  @date 2008-03-03
       */
      template <class T, class C, class R, unsigned int D>
      inline 
      ROOT::Math::SVector<T,D>
      operator* 
      ( const ROOT::Math::PositionVector3D<C>& vect , 
        const ROOT::Math::SMatrix<T,3,D,R>&    mtrx )
      {
        const ROOT::Math::SVector<T,3> vct 
          ( vect . X () , vect . Y () , vect . Z () ) ;
        return vct * mtrx ;
      }
      // ======================================================================      
      template <class T, class C, class B, class R, unsigned int D>
      inline 
      ROOT::Math::SVector<T,D>
      operator* 
      ( const ROOT::Math::PositionVector3D<C>& vect , 
        const ROOT::Math::Expr<B,T,3,D,R>&     mtrx )
      {
        const ROOT::Math::SVector<T,3> vct 
          ( vect . X () , vect . Y () , vect . Z () ) ;
        return vct * mtrx ;
      }
      // ======================================================================
    } //                                end of namespace Ostap::Math::Operators
    // ========================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap    
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // LHCBMATH_XXX_H
// ============================================================================
