// ============================================================================
#ifndef OSTAP_KINEMATICS_H 
#define OSTAP_KINEMATICS_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL 
// ============================================================================
#include <cmath>
// ============================================================================
// ROOT
// ============================================================================
#include "Math/SMatrix.h"
#include "Math/SVector.h"
#include "Math/Vector4D.h"
// ============================================================================
/** @file 
 *
 *  Collection of useful mathematical utilities related to the kinematics
 *
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2008-01-15
 */
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {   
    // ========================================================================
    /** helper function ("signed sqrt") 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2010-05-24
     */
    inline double _sqrt_ ( const double value ) 
    { return  0 <= value ? std::sqrt ( value ) : -std::sqrt ( - value ) ; }
    // ========================================================================
    /** evaluate the dispersion of M^2 from the particle 4-vector and 
     *  the covarinace matrix
     *
     *  @code
     *
     *   const LHCb::Particle* p = ... ;
     *   double s2m2 = sigma2mass2 ( p -> momentum() , p -> momCovMatrix() ) ; 
     *
     *  @endcode
     *  
     *  @param momentum   (in) the particle momentum
     *  @param covariance (in) 4x4 covarinnce matrix
     *  @return the estimate for dispersion of M^2
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2008-01-15
     */
    template <class C, class T>
    inline double
    sigma2mass2 
    ( const ROOT::Math::LorentzVector<C>&                           momentum   , 
      const ROOT::Math::SMatrix<T,4,4,ROOT::Math::MatRepSym<T,4> >& covariance ) 
    {
      // get the vector d(M2)/dp_i :
      ROOT::Math::SVector<T,4> dM2dp;
      dM2dp [0] = -2 * momentum.Px () ;
      dM2dp [1] = -2 * momentum.Py () ;
      dM2dp [2] = -2 * momentum.Pz () ;
      dM2dp [3] =  2 * momentum.E  () ;
      //
      return ROOT::Math::Similarity ( covariance , dM2dp ) ;
    }
    // ========================================================================
    template <class C, class T, class B , class R>
    inline double
    sigma2mass2 
    ( const ROOT::Math::LorentzVector<C>& momentum   , 
      const ROOT::Math::Expr<B,T,4,4,R> & matrix     ) 
    {
      return sigma2mass2 
        ( momentum , 
          ROOT::Math::SMatrix<T,4,4,ROOT::Math::MatRepSym<T,4> > ( matrix ) ) ;
    }
    // ========================================================================
    /** evaluate the dispersion of M from the particle 4-vector and 
     *  the covarinace matrix
     *
     *  @code
     *
     *   const LHCb::Particle* p = ... ;
     *   double s2m = sigma2mass ( p -> momentum() , p -> momCovMatrix() ) ; 
     *
     *  @endcode
     *  
     *  @attention the correct result is returned only for time-like vectors! 
     * 
     *  @param momentum   (in) the particle momentum
     *  @param covariance (in) 4x4 covarinnce matrix
     *  @return the estimate for dispersion of M
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2008-01-15
     */
    template <class C, class T>
    inline double
    sigma2mass 
    ( const ROOT::Math::LorentzVector<C>&                           momentum   , 
      const ROOT::Math::SMatrix<T,4,4,ROOT::Math::MatRepSym<T,4> >& covariance ) 
    {
      const double m2   = momentum.M2 () ;
      if ( m2   <= 0 ) { return 0 ; }                         // NB! RETURN 0
      const double s2m2 = sigma2mass2 ( momentum , covariance ) ;
      if ( s2m2 <= 0 ) { return 0 ; }                         // NB! RETURN 0  
      //
      return 0.25 * s2m2 / m2 ;
    }
    // ========================================================================
    template <class C, class T, class B, class R>
    inline double
    sigma2mass 
    ( const ROOT::Math::LorentzVector<C>& momentum   , 
      const ROOT::Math::Expr<B,T,4,4,R> & covariance ) 
    {
      const double s2m2 = sigma2mass2 ( momentum , covariance ) ;
      const double m2   = momentum.M2 () ;
      if ( 0 != m2 ) { return 0.25 * s2m2 / m2 ; }                    // RETURN
      return -1.E+24;                                                 // RETURN 
    }
    // ========================================================================
    /** evaluate sigma(M) from the particle 4-vector and 
     *  the covarinace matrix
     *
     *  @code
     *
     *   const LHCb::Particle* p = ... ;
     *   double sigma = sigmamass ( p -> momentum() , p -> momCovMatrix() ) ; 
     *
     *  @endcode
     *  
     *  @attention the correct result is returned only for time-like vectors!
     * 
     *  @param momentum   (in) the particle momentum
     *  @param covariance (in) 4x4 covarinnce matrix
     *  @return the estimate for dispersion of M
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2008-01-15
     */
    template <class C, class T>
    inline double
    sigmamass 
    ( const ROOT::Math::LorentzVector<C>&                           momentum   , 
      const ROOT::Math::SMatrix<T,4,4,ROOT::Math::MatRepSym<T,4> >& covariance ) 
    {
      const double s2m = sigma2mass ( momentum , covariance ) ;
      return _sqrt_ ( s2m ) ;
    }
    // =========================================================================
    template <class C, class T, class B, class R>
    inline double
    sigmamass 
    ( const ROOT::Math::LorentzVector<C>& momentum   , 
      const ROOT::Math::Expr<B,T,4,4,R> & covariance ) 
    {
      const double s2m = sigma2mass ( momentum , covariance ) ;
      return _sqrt_ ( s2m ) ;
    }
    // ========================================================================
    /** evaluate the chi2 of the mass 
     *
     *  @code
     *
     *   const LHCb::Particle* B = ... ;
     *
     *   const double chi2 = 
     *       chi2mass ( 5.279 * Gaudi::Units::GeV , 
     *                  B -> momentum()           , 
     *                  B -> momCovMatrix()       ) ; 
     *
     *  @endcode
     *  
     *  @param mass       (in) nominal mass
     *  @param momentum   (in) 4-momentum of the particle
     *  @param covariance (in) 4x4-covariance matrix 
     *  @return chi2 of the delta mass 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2008-01-15
     */
    template <class C, class T>
    inline double chi2mass 
    ( const double                                                  mass       , 
      const ROOT::Math::LorentzVector<C>&                           momentum   , 
      const ROOT::Math::SMatrix<T,4,4,ROOT::Math::MatRepSym<T,4> >& covariance ) 
    {
      // sigma^2(M^2):
      const double s2 = 1.0 / Ostap::Math::sigma2mass2 ( momentum , covariance ) ;
      // delta(M^2)
      const double dm2 = momentum.M2() - mass * mass ;
      //  (delta^2(M^2))/(sigma^2(M^2))
      return ( dm2 * dm2 ) * s2 ;
    }
    // ========================================================================
    template <class C, class T, class B, class R>
    inline double chi2mass 
    ( const double                        mass       , 
      const ROOT::Math::LorentzVector<C>& momentum   , 
      const ROOT::Math::Expr<B,T,4,4,R> & covariance ) 
    {
      // sigma^2(M^2):
      const double s2 = 1.0 / Ostap::Math::sigma2mass2 ( momentum , covariance ) ;
      // delta(M^2)
      const double dm2 = momentum.M2() - mass * mass ;
      //  (delta^2(M^2))/(sigma^2(M^2))
      return ( dm2 * dm2 ) * s2 ;
    }
    // ========================================================================
    /** evaluate the dispersion of p from the particle 4-vector and 
     *  the covarinace matrix
     *
     *  @code
     *
     *   const LHCb::Particle* p = ... ;
     *   double s2p = sigma2p ( p -> momentum() , p -> momCovMatrix() ) ; 
     *
     *  @endcode
     *  
     *  @param momentum   (in) the particle momentum
     *  @param covariance (in) 4x4 covarinnce matrix
     *  @return the estimate for dispersion of p
     *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
     *  @date 2010-05-25
     */
    template <class C, class T>
    inline double
    sigma2p
    ( const ROOT::Math::LorentzVector<C>&                           momentum   , 
      const ROOT::Math::SMatrix<T,4,4,ROOT::Math::MatRepSym<T,4> >& covariance ) 
    {
      // get the vector d(p)/dp_i :
      ROOT::Math::SVector<T,4> dPdP_i;
      const double P = momentum.P() ;
      dPdP_i [0] = momentum.Px () / P ;
      dPdP_i [1] = momentum.Py () / P ;
      dPdP_i [2] = momentum.Pz () / P ;
      dPdP_i [3] = 0.0 ;
      //
      return ROOT::Math::Similarity ( covariance , dPdP_i ) ;
    }
    // ========================================================================
    template <class C, class T, class B, class R>
    inline double
    sigma2p
    ( const ROOT::Math::LorentzVector<C>& momentum   , 
      const ROOT::Math::Expr<B,T,4,4,R>&  covariance ) 
    {
      // get the vector d(p)/dp_i :
      ROOT::Math::SVector<T,4> dPdP_i;
      const double P = momentum.P() ;
      dPdP_i [0] = momentum.Px () / P ;
      dPdP_i [1] = momentum.Py () / P ;
      dPdP_i [2] = momentum.Pz () / P ;
      dPdP_i [3] = 0.0 ;
      //
      return ROOT::Math::Similarity ( covariance , dPdP_i ) ;
    }
    // ========================================================================
    /** evaluate the sigma of |p| from the particle 4-vector and 
     *  the covarinace matrix
     *
     *  @code
     *
     *   const LHCb::Particle* p = ... ;
     *   double sp = sigmap ( p -> momentum() , p -> momCovMatrix() ) ; 
     *
     *  @endcode
     *  
     *  @param momentum   (in) the particle momentum
     *  @param covariance (in) 4x4 covarinnce matrix
     *  @return the estimate for sigma of p
     *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
     *  @date 2010-05-25
     */
    template <class C, class T>
    inline double
    sigmap
    ( const ROOT::Math::LorentzVector<C>&                           momentum   , 
      const ROOT::Math::SMatrix<T,4,4,ROOT::Math::MatRepSym<T,4> >& covariance ) 
    {
      const double s2p = sigma2p ( momentum , covariance ) ;
      return _sqrt_ ( s2p ) ;
    }
    // ========================================================================
    template <class C, class T, class B, class R>
    inline double
    sigmap
    ( const ROOT::Math::LorentzVector<C>& momentum   , 
      const ROOT::Math::Expr<B,T,4,4,R> & covariance ) 
    {
      const double s2p = sigma2p ( momentum , covariance ) ;
      return _sqrt_ ( s2p ) ;
    }
    // ========================================================================
    /** evaluate the dispersion of pt from the particle 4-vector and 
     *  the covariance matrix
     *
     *  @code
     *
     *   const LHCb::Particle* p = ... ;
     *   double s2pt = sigma2pt ( p -> momentum() , p -> momCovMatrix() ) ; 
     *
     *  @endcode
     *  
     *  @param momentum   (in) the particle momentum
     *  @param covariance (in) 4x4 covarinnce matrix
     *  @return the estimate for dispersion of pt
     *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
     *  @date 2010-05-25
     */
    template <class C, class T>
    inline double
    sigma2pt
    ( const ROOT::Math::LorentzVector<C>&                           momentum   , 
      const ROOT::Math::SMatrix<T,4,4,ROOT::Math::MatRepSym<T,4> >& covariance ) 
    {
      // get the vector d(pt)/dp_i :
      const double _Pt = momentum.Pt () ;
      const double _ax = momentum.Px () / _Pt ;
      const double _ay = momentum.Py () / _Pt ;
      //
      return 
        covariance ( 0 , 0 ) * _ax * _ax       + 
        covariance ( 0 , 1 ) * _ax * _ay * 2.0 +
        covariance ( 1 , 1 ) * _ay * _ay       ;  
    }
    // ========================================================================
    template <class C, class T, class B, class R>
    inline double
    sigma2pt
    ( const ROOT::Math::LorentzVector<C>& momentum   , 
      const ROOT::Math::Expr<B,T,4,4,R> & covariance ) 
    {
      // get the vector d(pt)/dp_i :
      const double _Pt = momentum.Pt () ;
      const double _ax = momentum.Px () / _Pt ;
      const double _ay = momentum.Py () / _Pt ;
      //
      return 
        covariance ( 0 , 0 ) * _ax * _ax       + 
        covariance ( 0 , 1 ) * _ay * _ay * 2.0 +
        covariance ( 1 , 1 ) * _ay * _ay       ;
    }
    // ========================================================================
    /** evaluate the sigma of pt from the particle 4-vector and 
     *  the covarinace matrix
     *
     *  @code
     *
     *   const LHCb::Particle* p = ... ;
     *   double spt = sigmapt ( p -> momentum() , p -> momCovMatrix() ) ; 
     *
     *  @endcode
     *  
     *  @param momentum   (in) the particle momentum
     *  @param covariance (in) 4x4 covarinnce matrix
     *  @return the estimate for sigma of pt
     *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
     *  @date 2010-05-25
     */
    template <class C, class T>
    inline double
    sigmapt
    ( const ROOT::Math::LorentzVector<C>&                           momentum   , 
      const ROOT::Math::SMatrix<T,4,4,ROOT::Math::MatRepSym<T,4> >& covariance ) 
    {
      const double s2pt = sigma2pt ( momentum , covariance ) ;
      return _sqrt_ ( s2pt ) ;
    }
    // ========================================================================
    template <class C, class T, class B, class R>
    inline double
    sigmapt
    ( const ROOT::Math::LorentzVector<C>& momentum   , 
      const ROOT::Math::Expr<B,T,4,4,R> & covariance ) 
    {
      const double s2pt = sigma2pt ( momentum , covariance ) ;
      return _sqrt_ ( s2pt ) ;
    }
    // ========================================================================
    /** evaluate the dispersion of rapidity from the particle 4-vector and 
     *  the covariance matrix
     *
     *  @code
     *
     *   const LHCb::Particle* p = ... ;
     *   double s2y = sigma2y ( p -> momentum() , p -> momCovMatrix() ) ; 
     *
     *  @endcode
     *  
     *  @param momentum   (in) the particle momentum
     *  @param covariance (in) 4x4 covarinnce matrix
     *  @return the estimate for dispersion of pt
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2016-02-12
     */
    template <class C, class T>
    inline double
    sigma2y
    ( const ROOT::Math::LorentzVector<C>&                           momentum   , 
      const ROOT::Math::SMatrix<T,4,4,ROOT::Math::MatRepSym<T,4> >& covariance ) 
    {
      // get the vector d(Y)/dp_i :
      ROOT::Math::SVector<T,4> dYdP_i;
      const double ePpz = 0.5 / ( momentum.E() + momentum.Pz() ) ;
      const double eMpz = 0.5 / ( momentum.E() - momentum.Pz() ) ;
      dYdP_i [0] = 0.0 ; 
      dYdP_i [1] = 0.0 ; 
      dYdP_i [2] = ePpz + eMpz ;
      dYdP_i [3] = ePpz - eMpz ;
      //
      return ROOT::Math::Similarity ( covariance , dYdP_i ) ;
    }
    // ========================================================================
    template <class C, class T, class B, class R>
    inline double
    sigma2y
    ( const ROOT::Math::LorentzVector<C>& momentum   , 
      const ROOT::Math::Expr<B,T,4,4,R> & covariance ) 
    {
      // get the vector d(Y)/dp_i :
      ROOT::Math::SVector<T,4> dYdP_i;
      const double ePpz = 0.5 / ( momentum.E() + momentum.Pz() ) ;
      const double eMpz = 0.5 / ( momentum.E() - momentum.Pz() ) ;
      dYdP_i [0] = 0.0 ; 
      dYdP_i [1] = 0.0 ; 
      dYdP_i [2] = ePpz + eMpz ;
      dYdP_i [3] = ePpz - eMpz ;
      //
      return ROOT::Math::Similarity ( covariance , dYdP_i ) ;
    }
    // ========================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_KINEMATICS_H
// ============================================================================
