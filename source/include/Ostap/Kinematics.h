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
// Ostap 
// ============================================================================
#include "Ostap/Vector3DTypes.h"
#include "Ostap/Vector4DTypes.h"
// ============================================================================
/** @file Ostap/Kinematics.h
 *  Collection of useful mathematical utilities related to the kinematics
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
     *  const LHCb::Particle* p = ... ;
     *  double s2m = sigma2mass ( p -> momentum() , p -> momCovMatrix() ) ; 
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
     *  const LHCb::Particle* p = ... ;
     *  double sigma = sigmamass ( p -> momentum() , p -> momCovMatrix() ) ; 
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
     *   const LHCb::Particle* B = ... ;
     *   const double chi2 = 
     *       chi2mass ( 5.279 * GeV , 
     *                  B -> momentum()           , 
     *                  B -> momCovMatrix()       ) ; 
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
     *   const LHCb::Particle* p = ... ;
     *   double s2p = sigma2p ( p -> momentum() , p -> momCovMatrix() ) ; 
     *  @endcode
     *  
     *  @param momentum   (in) the particle momentum
     *  @param covariance (in) 4x4 covarinnce matrix
     *  @return the estimate for dispersion of p
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
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
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
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
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
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
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
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
     *  const LHCb::Particle* p = ... ;
     *  double s2y = sigma2y ( p -> momentum() , p -> momCovMatrix() ) ; 
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
  } //                                        The end of namespace Ostyap::Math
  // ==========================================================================
  namespace Kinematics 
  {
    // ========================================================================
    /** calculate the ``triangle'' function, aka ``lambda'' or ``Kallen'' function 
     *  \f[ \lambda ( a , b, c ) = a^2 + b^2 + c^2 - 2ab - 2bc - 2ca \f]
     *  @see see https://en.wikipedia.org/wiki/K%C3%A4ll%C3%A9n_function          
     *  @param a parameter a
     *  @param b parameter b
     *  @param c parameter c
     *  @return the value of triangle function
     */
    double triangle
    ( const double a ,
      const double b ,
      const double c ) ;
    // ========================================================================
    /** calculate the ``triangle'' function, aka ``lambda'' or ``Kallen'' function 
     *  \f[ \lambda ( a , b, c ) = a^2 + b^2 + c^2 - 2ab - 2bc - 2ca \f]
     *  @see see https://en.wikipedia.org/wiki/K%C3%A4ll%C3%A9n_function          
     *  @param a parameter a
     *  @param b parameter b
     *  @param c parameter c
     *  @return the value of triangle function
     */
    inline double kallen 
    ( const double a ,
      const double b ,
      const double c ) { return triangle ( a , b , c ) ; }
    // ========================================================================
    /** universal four-particle kinematical function, aka "tetrahedron-function"
     *  @see E.Byckling, K.Kajantie, "Particle kinematics", John Wiley & Sons,
     *              London, New York, Sydney, Toronto, 1973, p.89, eq. (5.23)
     *  @see https://userweb.jlab.org/~rafopar/Book/byckling_kajantie.pdf
     *  E.g. physical range for 2->2 scattering process is defined as
     * \f$ G(s,t,m_2^2, m_a^2, m_b^2, m_1^2) \le 0 \f$
     * or the physical range  for Dalitz plot is
     * \f$ G(s_2, s_1,  m_3^2, m_1^2, s , m_2^2) \le 0 \f$
     *
     * Actually the formula in E.Byckling & K.Kajantie  has a typo.
     * 
     * See the correct formula in: 
     * @see  P. Nyborg, H.S. Song, W. Kernan, R.H. Good,
     *       Phase-Space Considerations for Four-Particle Final States"
     *       Phys.Rev. 140 (1965) B914-B920, DOI: 10.1103/PhysRev.140.B914  
     * @see https://journals.aps.org/pr/pdf/10.1103/PhysRev.140.B914
     * @see http://inspirehep.net/record/49679?ln=en
     * 
     *  The correct formula is :
     *  \f[ G(x,y,z.u,v,w) =  
     *   x^2y + xy^2 + z^2u + zu^2  + v^2w + vw^2 
     *   + xzw + zuv + yzv + yuw 
     *   - x y ( z + u + v + w ) 
     *   - z u ( x + y + v + w ) 
     *   - v w ( x + y + z + u ) \f] 
     *  - B&K gives \f$ yzw\f$ instead of \f$yzv\f$  
     */ 
    double G 
    ( const double x , 
      const double y , 
      const double z , 
      const double u ,
      const double v , 
      const double w ) ;
    // ========================================================================
    /** universal four-particle kinematical function, aka "tetrahedron-function"
     *  @see E.Byckling, K.Kajantie, "Particle kinematics", John Wiley & Sons,
     *              London, New York, Sydney, Toronto, 1973, p.89, eq. (5.23)
     *  @see https://userweb.jlab.org/~rafopar/Book/byckling_kajantie.pdf
     *  E.g. physical range for 2->2 scattering process is defined as
     * \f$ G(s,t,m_2^2, m_a^2, m_b^2, m_1^2) \le 0 \f$
     * or the physical range  for Dalitz plot is
     * \f$ G(s_2, s_1,  m_3^2, m_1^2, s , m_2^2) \le 0 \f$ 
     * 
     * Actually the formula in E.Byckling & K.Kajantie  has a typo.
     * 
     * See the correct formula in: 
     * @see  P. Nyborg, H.S. Song, W. Kernan, R.H. Good,
     *       Phase-Space Considerations for Four-Particle Final States"
     *       Phys.Rev. 140 (1965) B914-B920, DOI: 10.1103/PhysRev.140.B914  
     * @see https://journals.aps.org/pr/pdf/10.1103/PhysRev.140.B914
     * @see http://inspirehep.net/record/49679?ln=en
     */ 
    inline double tetrahedron 
    ( const double x , 
      const double y , 
      const double z , 
      const double u ,
      const double v , 
      const double w ) { return G ( x , y , z , u , v , w ) ; }
    // ========================================================================
    /** @class Gram
     *  Calculate a few simplest Gram determinants:
     *  \f[ G\left( \begin{array}{lcr} p_1, & ... & p_n \\ 
     *                                 q_1, & ... & q_n \end{array} \right) 
     *  \equiv 
     *  \left| \begin{array}{lcl}
     *    p_1q_1   & ... &   p_1q_n \\
     *    ...      & ... &   ...    \\
     *    p_nq_1   & ... &   p_nq_n  \end{array} \right| \f]          
     *
     *  @see E.Byckling, K.Kajantie, "Particle kinematics", John Wiley & Sons,
     *              London, New York, Sydney, Toronto, 1973, p.89, eq. (5.23)
     *  @see https://userweb.jlab.org/~rafopar/Book/byckling_kajantie.pdf     
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2019-07012
     */
    class Gram 
    {
      // ======================================================================
    public:
      // ======================================================================
      /** non-symmetric Gram determinant
       *  \f[ G \left(\begin{array}{ll} p_1  & p_2 \\
       *                                q_1  & q_2 \end{array}\right) 
       *  = \left| \begin{array}{ll}
       *  p_1q_1 & p_1q_2 \\
       *  p_2q_1 & p_2q_2 \end{array} \right| \f]
       */
      static double G 
      ( const Ostap::LorentzVector& p1 , 
        const Ostap::LorentzVector& p2 ,
        const Ostap::LorentzVector& q1 , 
        const Ostap::LorentzVector& q2 ) ;
      // ======================================================================
      /** non-symmetric Gram determinant
       *  \f[ G \left(\begin{array}{lll} p_1  & p_2  & p_3 \\
       *                                 q_1  & q_2  & q_3 \end{array}\right) 
       *  = \left| \begin{array}{lll}
       *  p_1q_1 & p_1q_2 & p_1q_3 \\
       *  p_2q_1 & p_2q_2 & p_2q_3 \\
       *  p_3q_1 & p_3q_2 & p_3q_3 \end{array} \right| \f]
       */
      static double G 
      ( const Ostap::LorentzVector& p1 , 
        const Ostap::LorentzVector& p2 ,
        const Ostap::LorentzVector& p3 ,
        const Ostap::LorentzVector& q1 ,
        const Ostap::LorentzVector& q2 ,
        const Ostap::LorentzVector& q3 ) ;
      // ======================================================================
      /** non-symmetric Gram determinant
       *  \f[ G \left(\begin{array}{llll} p_1  & p_2  & p_3 & p_4 \\
       *                                  q_1  & q_2  & q_3 & q_4 \end{array}\right) 
       *  = \left| \begin{array}{llll}
       *  p_1q_1 & p_1q_2 & p_1q_3 & p_1q_4\\
       *  p_2q_1 & p_2q_2 & p_2q_3 & p_2q_4\\
       *  p_3q_1 & p_3q_2 & p_3q_3 & p_3q_4\\ 
       *  p_4q_1 & p_4q_2 & p_4q_3 & p_4q_4\end{array} \right| \f]
       */
      static double G 
      ( const Ostap::LorentzVector& p1 , 
        const Ostap::LorentzVector& p2 ,
        const Ostap::LorentzVector& p3 ,
        const Ostap::LorentzVector& p4 ,
        const Ostap::LorentzVector& q1 ,
        const Ostap::LorentzVector& q2 ,
        const Ostap::LorentzVector& q3 ,
        const Ostap::LorentzVector& q4 ) ;
      // ======================================================================
    public: 
      // ======================================================================
      /** symmetric Gram determinant
       *  @see E.Byckling, K.Kajantie, "Particle kinematics", John Wiley & Sons,
       *              London, New York, Sydney, Toronto, 1973, p.89, eq. (5.23)
       *  @see https://userweb.jlab.org/~rafopar/Book/byckling_kajantie.pdf     
       */  
      static double Delta  
      ( const Ostap::LorentzVector& p1 ) ;
      // ======================================================================
      /** symmetric Gram determinant
       *  \f[ Delta( p_1, p_2) \equiv 
       *   G \left( \begin{array}{ll}p_1 & p_2 \\ 
       *                             p_1 & p_2 \end{array} \right) \f] 
       *  @see E.Byckling, K.Kajantie, "Particle kinematics", John Wiley & Sons,
       *              London, New York, Sydney, Toronto, 1973, p.89, eq. (5.23)
       *  @see https://userweb.jlab.org/~rafopar/Book/byckling_kajantie.pdf     
       */  
      static double Delta  
      ( const Ostap::LorentzVector& p1 , 
        const Ostap::LorentzVector& p2 ) ;
      // ======================================================================
      /** symmetric Gram determinant
       *  \f[ Delta( p_1, p_2, p_3 ) \equiv 
       *   G \left( \begin{array}{lll}p_1 & p_2 & p_3 \\ 
       *                              p_1 & p_2 & p_3 \end{array} \right) \f] 
       *  @see E.Byckling, K.Kajantie, "Particle kinematics", John Wiley & Sons,
       *              London, New York, Sydney, Toronto, 1973, p.89, eq. (5.23)
       *  @see https://userweb.jlab.org/~rafopar/Book/byckling_kajantie.pdf     
       */
      static double Delta
      ( const Ostap::LorentzVector& p1 , 
        const Ostap::LorentzVector& p2 ,
        const Ostap::LorentzVector& p3 ) ;
      // ======================================================================
      /** symmetric Gram determinant
       *  \f[ Delta( p_1, p_2, p_3, p_4 ) \equiv 
       *   G \left( \begin{array}{llll}p_1 & p_2 &p_3 &p_4 \\ 
       *                               p_1 & p_2 &p_3 &p_4 \end{array} \right) \f] 
       *  @see E.Byckling, K.Kajantie, "Particle kinematics", John Wiley & Sons,
       *              London, New York, Sydney, Toronto, 1973, p.89, eq. (5.23)
       *  @see https://userweb.jlab.org/~rafopar/Book/byckling_kajantie.pdf     
       */
      static double Delta
      ( const Ostap::LorentzVector& p1 , 
        const Ostap::LorentzVector& p2 ,
        const Ostap::LorentzVector& p3 ,
        const Ostap::LorentzVector& p4 ) ;
      // ======================================================================
    } ;    
    // ========================================================================
    /** boost Lorentz vector into  rest-frame of another Lorentz vector 
     *  @param what   the vextro to be bosted 
     *  @param frame  the 4-vector of the frame 
     *  @return boosted vector 
     */
      Ostap::LorentzVector boost 
      ( const Ostap::LorentzVector& what  ,
        const Ostap::LorentzVector& frame ) ;  
    // ========================================================================
    /** simple function which evaluates the magnitude of 3-momentum
     *  of particle "v" in the rest system of particle "M"
     *
     *  \f[ \left|\vec{p}\right| =
     *     \sqrt{  \frac{\left(v\cdot M\right)^2}{M^2} -v^2} \f]
     *
     *  Note that this is clear Lorentz invarinat expresssion.
     *
     *  @attention particle M must be time-like particle: M^2 > 0 !
     *  For invalid configurations large negative number is returned 
     *
     *  @param v the vector to be checked
     *  @param M the defintion of "rest"-system
     *  @return the magnitude of 3D-momentum of v in rest-frame of M
     *  @date 2008-07-27
     */
    double restMomentum 
    ( const Ostap::LorentzVector& v ,
      const Ostap::LorentzVector& M ) ;
    // ========================================================================
    /** simple function which evaluates the energy
     *  of particle "v" in the rest system of particle "M"
     *
     *  \f[ e = \frac{ v \cdot M }{\sqrt{ M^2 } } \f]
     *
     *  Note that this is clear Lorentz invarinat expresssion.
     *
     *  @attention particle M must be time-like particle: M^2 > 0 !
     *
     *  @param v the vector to be checked
     *  @param M the defintion of "rest"-system
     *  @return the energy of v in rest-frame of M
     *  @date 2008-07-27
     */
    double restEnergy 
    ( const Ostap::LorentzVector& v ,
      const Ostap::LorentzVector& M ) ;
    // =======================================================================
    /** simple function for evaluation of the euclidiam norm
     *  for LorentzVectors (E**2+Px**2+Py**2+Pz**2)
     *  @param vct the vector
     *  @return euclidian norm squared
     *  @date 2006-01-17
     */
    double euclidianNorm2 ( const Ostap::LorentzVector& vct ) ;
    // ========================================================================
    /** simple function which evaluates the transverse
     *  momentum with respect a certain 3D-direction:
     *
     * \f[ r_T = \left| \vec{r} \right| = \left| \vec{v} -
     *    \vec{d}\frac{\left(\vec{v}\vec{d}\right)}
     *  { \left| \vec{d} \right|^2 } \right| \f]
     *
     *  @param mom the momentum
     *  @param dir the direction
     *  @return the transverse moementum with respect to the direction
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-01-17
     */
    double transverseMomentumDir
    ( const Ostap::Vector3D& mom, 
      const Ostap::Vector3D& dir );
    // ========================================================================
    /** simple function which evaluates the transverse
     *  momentum with respect a certain 3D-direction:
     *
     * \f[ r_T = \left| \vec{r} \right| =  = \left| \vec{v} -
     * \vec{d}\frac{\left(\vec{v}\vec{d}\right)}
     *  { \left| \vec{d} \right|^2 } \right| \f]
     *
     *  @param mom the momentum
     *  @param dir the direction
     *  @return the transverse moementum with respect to the direction
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-01-17
     */
    inline double transverseMomentumDir
    ( const Ostap::LorentzVector& mom, 
      const Ostap::Vector3D&      dir ) 
    { return transverseMomentumDir( mom.Vect(), dir ); }
    // ========================================================================
    /** This routine returns the cosine angle theta
     *  The decay angle calculated  is that between
     *  the flight direction of the daughter meson, "D",
     *  in the rest frame of "Q" (the parent of "D"),
     *  with respect to "Q"'s flight direction in "P"'s
     *  (the parent of "Q") rest frame
     *
     *  \f[
     *  \cos \theta = \frac
     *  { \left(P \cdot D\right)Q^2 -
     *    \left(P \cdot Q\right)\left(D \cdot Q \right) }
     *  {\sqrt{ \left[ \left( P \cdot Q \right)^2 - Q^2 P^2 \right]
     *          \left[ \left( D \cdot Q \right)^2 - Q^2 D^2 \right] } }
     *  \f]
     *
     *  Note that the expression has the symmetry: \f$ P \leftrightarrow D \f$
     *
     *  Essentially it is a rewritten <c>EvtDecayAngle(P,Q,D)</c>
     *  routine from EvtGen package
     *
     *  @param D 4-momentum of the daughter particle
     *  @param Q 4-momentum of mother particle
     *  @param P "rest frame system"
     *  @return cosine of decay angle
     *
     *  @see Ostap::LorentzVector
     *
     *  - The function <c>Ostap::Kinematics::cos_theta</c> provides  
     *    the alternative way to calculate the same quantity as 
     *    <code>Ostap::Kinematics::cos_theta ( P , -D , Q )</code>
     *
     *  - The function <c>Ostap::Kinematics::cosThetaRest</c> provides  
     *    the alternative way to calculate the same quantity as 
     *    <code>Ostap::Kinematics::cosThetaRest ( P , -D , Q )</code>
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2004-12-03
     */
    double decayAngle 
    ( const Ostap::LorentzVector& P , 
      const Ostap::LorentzVector& Q ,
      const Ostap::LorentzVector& D ) ;    
    // ========================================================================
    /** This routine returns the cosine angle theta
     *  The decay angle calculated  is that between
     *  the flight direction of the daughter meson, "D",
     *  in the rest frame of "M" (the parent of "D"),
     *  with respect to the boost direction from
     *  "M"'s rest frame
     *
     *  @param D 4-momentum of the daughter particle
     *  @param M 4-momentum of mother particle
     *  @return cosine of decay angle
     *
     *  Clearly it is a variant of 3-argument with the
     *  P-argument to be of type (0,0,0,E)
     *  (=="laborator frame")
     *
     *  @see Ostap::LorentzVector
     *  @see Ostap::Kinematics::decayAngle 
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2004-12-03
     */
    double decayAngle
    ( const Ostap::LorentzVector& D , 
      const Ostap::LorentzVector& M ) ;
    // ========================================================================
    /** simple function to evaluate the cosine angle between
     *  two directions (v1 and v2) in the rest system of M
     *
     * \f[
     * \cos\theta =
     * \frac{\vec{p}_1\vec{p}_2}{\left|\vec{p}_1\right|
     * \left|\vec{p}_2\right|} =
     * \frac{1}{\left|\vec{p}_1\right|\left|\vec{p}_2\right|}
     * \left( E_1E_2 -\frac{1}{2}
     * \left(\left(v_1+v_2\right)^2-v_1^2-v_2^2 \right) \right),
     * \f]
     *
     *  where
     *  \f[
     *  E_1 E_2 = \frac{ \left ( v_1 \cdot M\right) \left (v_2 \cdot M \right ) }{M^2}
     *  \f]
     *  and
     *  \f[
     * \left|\vec{p}_1\right|\left|\vec{p}_2\right| =
     * \sqrt{
     * \left( \frac{\left(v_1\cdot M\right)^2}{M^2}-v_1^2 \right)
     *      \left( \frac{\left(v_2\cdot M\right)^2}{M^2}-v_2^2 \right) }
     * \f]
     *
     *  Note that the expressions are clear Lorentz invariant
     *
     *  @attention the particle M must be time-like particle: M^2 > 0 !
     *  @param v1 the first vector
     *  @param v2 the last vector
     *  @param M  the defintion of rest-system
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2008-07-27
     *  - The function <c>Ostap::Kinematics::cos_theta</c> provides  
     *    the alternative way to calcualet the same quantity 
     *  @see Ostap::Kinematics::cos_theta     
     */
    double cosThetaRest
    ( const Ostap::LorentzVector& v1 , 
      const Ostap::LorentzVector& v2 ,
      const Ostap::LorentzVector& M  ) ;
    // ========================================================================
    /** Cosine of the angle between p1 and p2 in the rest frame of M
     *  \f$ \cos \theta = - \frac 
     *  { G   \left( \begin{array}{ll} M, &p_1 \\ M,& p_2 \end{array}\right) }
     *  { \left[ \Delta_2(M,p_1)\Delta_2(M,p_2)\right]^{1/2} }\f$  
     *  @param v1 the first vector
     *  @param v2 the last vector
     *  @param M  the defintion of rest-system
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2019-07-17
     *  - The function <c>Ostap::Kinematics::cosThetaRest</c> provides  
     *    the alternative way to calcualet the same quantity 
     *  @see Ostap::Kinematics::cosThetaRest
     */
    double cos_theta 
    ( const Ostap::LorentzVector& p1 ,
      const Ostap::LorentzVector& p2 ,
      const Ostap::LorentzVector& M  ) ;  
    // ========================================================================
    /** Sine squared  of the angle between p1 and p2 in the rest frame of M
     *  \f$ \sin^2 \theta = \frac { \Delta ( M ) \Delta ( M , p_1 , p_2 ) }
     *  { \Delta ( M, p_1  ) \Delta ( M , p_2 ) } \f$ 
     *  @param v1 the first vector
     *  @param v2 the last vector
     *  @param M  the defintion of rest-system
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2019-07-17
     */
    double sin2_theta 
    ( const Ostap::LorentzVector& p1 ,
      const Ostap::LorentzVector& p2 ,
      const Ostap::LorentzVector& M  ) ;  
    // ========================================================================
    /** evaluate the angle \f$\chi\f$
     *  beween two decay planes,
     *  formed by particles v1&v2 and h1&h2 correspondingly.
     *  The angle is evaluated in the rest frame
     *  of "mother" particles (defined as v1+v2+h1+h2)
     *
     *  @param d1 the 1st daughter
     *  @param d2 the 2nd daughter
     *  @param h1 the 3rd daughter
     *  @param h2 the 4th daughter
     *  @return angle chi
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2008-07-27
     */
    double decayAngleChi 
    ( const Ostap::LorentzVector& d1 , 
      const Ostap::LorentzVector& d2 ,
      const Ostap::LorentzVector& h1 , 
      const Ostap::LorentzVector& h2 );
    // ========================================================================
    /** evaluate \f$\cos \chi\f$, where \f$\chi\f$ if the angle
     *  beween two decay planes, formed by particles d1&d2
     *  and h1&h2 correspondingly.
     *
     *  The angle is evaluated in the rest frame
     *  of "mother" particles (defined as d1+d2+h1+h2)
     *
     *  The angle is evaluated using the explicit
     *  Lorentz-invariant expression:
     *  \f[
     *  \cos \chi =
     *   - \frac{ L_D^{\mu} L_H^{\mu} }
     *     { \sqrt{ \left[ -L_D^2 \right]\left[ -L_H^2 \right] }},
     *   =
     *   - \frac{
     *     \epsilon_{ijkl}d_1^{j}d_2^{k}\left(h_1+h_2\right)^l
     *     \epsilon_{imnp}h_1^{m}h_2^{n}\left(d_1+d_2\right)^p }
     *     { \sqrt{ \left[ -L_D^2 \right]\left[ -L_H^2 \right] }},
     *  \f],
     *  where "4-normales" are defined as:
     *  \f[
     *   L_D^{\mu} = \epsilon_{\mu\nu\lambda\kappa}
     *                d_1^{\nu}d_2^{\lambda}\left(h_1+h_2\right)^{\kappa}
     *  \f]
     *   and
     *  \f[
     *   L_H^{\mu} = \epsilon_{\mu\lambda\delta\rho}
     *                h_1^{\lambda}h_2^{\delta}\left(d_1+d_2\right)^{\rho}
     *   \f]
     *
     *  @param d1 the 1st daughter
     *  @param d2 the 2nd daughter
     *  @param h1 the 3rd daughter
     *  @param h2 the 4th daughter
     *  @return cos(chi)
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2008-07-27
     */
    double cosDecayAngleChi
    ( const Ostap::LorentzVector& d1 , 
      const Ostap::LorentzVector& d2 ,
      const Ostap::LorentzVector& h1 , 
      const Ostap::LorentzVector& h2 ) ;
    // ========================================================================    
    /** evaluate \f$\sin\chi\f$, where \f$\chi\f$ is the angle
     *  beween two decay planes,
     *  formed by particles v1&v2 and h1&h2 correspondingly.
     *  The angle is evaluated in the rest frame
     *  of "mother" particles (defined as v1+v2+h1+h2)
     *
     *  The angle is  calculated using the explicit
     *   Lorentz-invariant expression:
     *  \f[ \sin \chi =
     *   \frac  {
     *   \epsilon_{\mu\nu\lambda\delta}
     *   L_D^{\mu}L_H^{\nu}H^{\lambda}M^{\delta} }
     *   { \sqrt{
     *   \left[ -L_D^2 \right]\left[ -L_H^2 \right]
     *   \left[ \left( H\ cdot M\right)^2-H^2M^2 \right]
     *   }} = \frac {
     *   \epsilon_{\mu\nu\lambda\delta}
     *   d_1^{\mu}d_2^{\nu}h_1^{\lambda}h_2^{\delta}
     *   \left( \left( D \cdot H \right)^2 - D^2H^2 \right) }
     *   { \sqrt{
     *   \left[ -L_D^2 \right]\left[ -L_H^2    \right]
     *   \left[ \left(H\cdot M\right)^2-H^2M^2 \right]
     *   }} \f],
     *  where "4-normales" are defined as:
     *  \f$ L_D^{\mu} \equiv \epsilon_{\mu\nu\lambda\kappa}
     *                d_1^{\nu}d_2^{\lambda}\left(h_1+h_2\right)^{\kappa} \f$,
     *  \f$ L_H^{\mu} \equiv \epsilon_{\mu\lambda\delta\rho}
     *  h_1^{\lambda}h_2^{\delta}\left(d_1+d_2\right)^{\rho} \f$
     *  and   \f$ D = d_1 + d_2 \f$,
     *        \f$ H = h_1 + h_2 \f$,
     *        \f$ M = D + H = d_1 + d_2 + h_1+h_2 \f$.
     *
     *  @param d1 the 1st daughter
     *  @param d2 the 2nd daughter
     *  @param h1 the 3rd daughter
     *  @param h2 the 4th daughter
     *  @return angle chi
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2008-07-27
     */
    double sinDecayAngleChi 
    ( const Ostap::LorentzVector& d1 , 
      const Ostap::LorentzVector& d2 ,
      const Ostap::LorentzVector& h1 , 
      const Ostap::LorentzVector& h2 ) ;
    // ========================================================================
    /** evaluate the Armenteros-Podolanski variable \f$\mathbf{\alpha}\f$,
     *  defined as:
     *  \f[
     *  \mathbf{\alpha} = \frac
     *  { \mathrm{p}^{\mathrm{L},1} - \mathrm{p}^{\mathrm{L},1} }
     *  { \mathrm{p}^{\mathrm{L},1} + \mathrm{p}^{\mathrm{L},1} },
     *  \f]
     *  where
     *   \f$ \mathrm{p}^{\mathrm{L},1}\f$ and
     *   \f$ \mathrm{p}^{\mathrm{L},2}\f$ are longitudinal momentum
     *   components for the first and the second daughter particles
     *   with respect to the total momentum direction.
     *
     *  Clearly this expression could be rewritten in an equivalent
     *  form which however much more easier for calculation:
     *  \f[
     *  \mathbf{\alpha} = \frac
     *  { \vec{\mathbf{p}}_1^2 - \vec{\mathbf{p}}_2^2 }
     *  { \left( \vec{\mathbf{p}}_1 + \vec{\mathbf{p}}_2 \right)^2 }
     *  \f]
     *
     *  @attention instead of
     *     2D \f$\left(\mathrm{p_T},\mathbf{\alpha}\right)\f$ diagram,
     *     in the case of twobody decays it is much better to
     *     use 2D diagram \f$\left(\cos \theta, \mathrm{m} \right)\f$
     *     diagram, where \f$\cos\theta\f$-is the decay
     *     angle, and \f$\mathrm{m}\f$ is an
     *     invariant evalauted for some (fixed) mass prescription,
     *     e.g. \f$\pi^+\pi^-\f$.
     *
     *  @param d1  three momentum of the first  daughter particle
     *  @param d2  three momentum of the second daughter particle
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2008-09-21
     */
    double armenterosPodolanskiX 
    ( const Ostap::Vector3D& d1 , 
      const Ostap::Vector3D& d2 );
    // ========================================================================
    /** evaluate the Armenteros-Podolanski variable \f$\mathbf{\alpha}\f$,
     *  defined as:
     *  \f[
     *  \mathbf{\alpha} = \frac
     *  { \mathrm{p}^{\mathrm{L},1} - \mathrm{p}^{\mathrm{L},1} }
     *  { \mathrm{p}^{\mathrm{L},1} + \mathrm{p}^{\mathrm{L},1} },
     *  \f]
     *  where
     *   \f$ \mathrm{p}^{\mathrm{L},1}\f$ and
     *   \f$ \mathrm{p}^{\mathrm{L},2}\f$ are longitudinal momentum
     *   components for the first and the second daughter particles
     *   with respect to the total momentum direction.
     *
     *  Clearly this expression could be rewritten in an equivalent
     *  form which however much more easier for calculation:
     *  \f[
     *  \mathbf{\alpha} = \frac
     *  { \vec{\mathbf{p}}_1^2 - \vec{\mathbf{p}}_2^2 }
     *  { \left( \vec{\mathbf{p}}_1 + \vec{\mathbf{p}}_2 \right)^2 }
     *  \f]
     *
     *  @attention instead of
     *     2D \f$\left(\mathrm{p_T},\mathbf{\alpha}\right)\f$ diagram,
     *     in the case of twobody decays at it is much better to
     *     use 2D diagram \f$\left(\cos \theta, \mathrm{m} \right)\f$
     *     diagram, where \f$\cos\theta\f$-is the decay
     *     angle, and \f$\mathrm{m}\f$ is an
     *     invariant evalauted for some (fixed) mass prescription,
     *     e.g. \f$\pi^+\pi^-\f$.
     *
     *  @param d1  four momentum of the first  daughter particle
     *  @param d2  four momentum of the second daughter particle
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2008-09-21
     */
    inline double armenterosPodolanskiX
    ( const Ostap::LorentzVector& d1 , 
      const Ostap::LorentzVector& d2 ) 
    { return armenterosPodolanskiX ( d1.Vect() , d2.Vect() ) ; }
    // ========================================================================
    /** trivial function to get the component of "a", transverse to "b"
     *  @param a (INPUT)  three vector
     *  @param b (INPUT)  reference direction
     *  @return component of "a", transverse to "b"
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-02-04
     */
    Ostap::Vector3D transverse
    ( const Ostap::ThreeVector& a, 
      const Ostap::ThreeVector& b );
    // ========================================================================
    /** trivial function to get the component of "a" parallel to "b"
     *  @param a (INPUT)  three vector
     *  @param b (INPUT)  reference direction
     *  @return component of "a" parallel to "b"
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-02-04
     */
    Ostap::Vector3D parallel
    ( const Ostap::Vector3D& a  , 
      const Ostap::Vector3D& b  ) ;
    // ========================================================================
    /** momentum of the first partle form two-body decay
     *  \f$ m\rightarrow m_1 m_2 \f$ in th eret frasme of \f$ m \f$.
     *  \f[ q ( m , m_1 , m_2 )  \equiv 
     *  \frac{\lambda^{1/2}\left( m^2, m_1^2, m_2^2\right)}{2m} \f]
     */
    double q 
    ( const double m  , 
      const double m1 ,
      const double m2 ) ;
    // ========================================================================
    /** two-body phase space: 
     *  \f[ \Phi_2( m ) = \dfrac{1}{8\pi} 
     *   \dfrac{ \sqrt{\lambda \left( m^2 , m_1^2, m_2^2 \right) }}{m^2} \f]
     * Note that there exist also an alternative normalization:
     *  \f[ \Phi_2^{\prime}( m ) = \dfrac{\pi}{2} 
     *   \dfrac{ \sqrt{\lambda \left( m^2 , m_1^2, m_2^2 \right) }}{m^2} \f],
     *  this one is used e.g. in  
     *  E.Byckling, K.Kajantie, "Particle kinematics", John Wiley & Sons,
     *              London, New York, Sydney, Toronto, 1973, Eq. (V.1.9)
     *  @see Ostap::Kinematics::phasespace2_bk 
     */
    double phasespace2
    ( const double x  , 
      const double m1 , 
      const double m2 ) ; 
    // ========================================================================
    /** two-body phase space: 
     *  \f[ \Phi_2^{\prime}( m ) = \dfrac{\pi}{2} 
     *   \dfrac{ \sqrt{\lambda \left( m^2 , m_1^2, m_2^2 \right) }}{m^2} \f],
     *  This one is used e.g. in  
     *  E.Byckling, K.Kajantie, "Particle kinematics", John Wiley & Sons,
     *              London, New York, Sydney, Toronto, 1973, Eq. (V.1.9)
     * Note that there exist also an alternative normalization:
     *  \f[ \Phi_2( m ) = \dfrac{1}{8\pi} 
     *   \dfrac{ \sqrt{\lambda \left( m^2 , m_1^2, m_2^2 \right) }}{m^2} \f]
     *  @see Ostap::Kinematics::phasespace2_bk 
     */
    double phasespace2_bk 
    ( const double x  , 
      const double m1 , 
      const double m2 ) ; 
    // ========================================================================
    /** three-body phase space, analytic symmetric expression via 
     *  elliptic  integrals 
     *  @see https://indico.cern.ch/event/368497/contributions/1786992/attachments/1134067/1621999/davydychev.PDF
     *  @see http://cds.cern.ch/record/583358/files/0209233.pdf
     *  @see https://www.researchgate.net/publication/2054534_Three-body_phase_space_symmetrical_treatments
     */
    double phasespace3 ( const double x  , 
                         const double m1 , 
                         const double m2 , 
                         const double m3 ) ; 
    // ========================================================================
  } //                                      end of namespace Ostap::Kinenmatics 
  // ==========================================================================
  namespace Math 
  {
    using Ostap::Kinematics::q        ;
    using Ostap::Kinematics::kallen   ;
    using Ostap::Kinematics::triangle ;    
  }
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
namespace ROOT 
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    /** scalar product of two 4-vectors 
     */
    template <class TYPE>
    inline typename TYPE::Scalar 
    operator* ( const LorentzVector<TYPE>& a , 
                const LorentzVector<TYPE>& b ) { return a.Dot ( b ) ; }
    // ========================================================================
  } //                                          The end of namespace ROOT::Math
  // ==========================================================================
} //                                                  The end of namespace ROOT
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_KINEMATICS_H
// ============================================================================
