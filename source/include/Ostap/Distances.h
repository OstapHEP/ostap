// ==========================================================================
#ifndef OSTAP_MATRIX_DISTANCES_H 
#define OSTAP_MATRIX_DISTANCES_H 1
// ==========================================================================
// Include files
// ==========================================================================
// STD& STL
// ==========================================================================
#include <cmath>
// ==========================================================================
// Ostap::
// ==========================================================================
#include "Ostap/MatrixUtils.h"
#include "Ostap/EigenSystem.h"
// ==========================================================================
/** @file Ostap/Distances.h
 *  Collection of functions defining "distances" between two datasets 
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date   2020-11-21
 */
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    /** Get the (asymmetrical)  Kullback-Leibler divergency for two objects 
     *  @see https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence
     *  @param v1 the first  data vector 
     *  @param c1 the first  covariance matrix 
     *  @param v2 the second data vector 
     *  @param c2 the second covariance matrix 
     *  @return (asymmetric) Kullback-Leibler divergency, or -999 
     */
    template <unsigned int N, typename SCALAR>
    inline double 
    kullback_leibler
    ( const ROOT::Math::SVector<SCALAR,N>&                                    v1 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c1 , 
      const ROOT::Math::SVector<SCALAR,N>&                                    v2 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c2 )
    {
      /// the actual type of covariance matrix
      typedef typename ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> > COV ;
      ///
      static const double s_bad = -999 ;
      ///
      // Specialization for N = 1 to avoid Cling warnings and overhead
      if constexpr ( N == 1 ) 
      {
        const double var1 = c1 ( 0 , 0 ) ;
        const double var2 = c2 ( 0 , 0 ) ;
        if ( var1 <= 0 || var2 <= 0 ) { return s_bad ; }

        const double diff     = v2 [ 0 ] - v1 [ 0 ];
        const double inv_var2 = 1.0 / var2;

        // For N=1: trace(inv(C2) * C1) becomes var1 / var2, and N becomes 1
        return 0.5 * ( (var1 * inv_var2) - 1.0 + (diff * diff * inv_var2) + std::log ( var2 / var1 ) ) ;
      }
      //
      SCALAR det1 = 1 ;
      if ( !c1.Det2 ( det1 ) || det1 <= 0 ) { return s_bad ; }
      SCALAR det2 = 1 ;
      if ( !c2.Det2 ( det2 ) || det2 <= 0 ) { return s_bad ; }
      //
      /// try to invert matrices 
      COV g2 { c2 } ;
      if  ( !g2.InvertChol () ) { return s_bad ; }
      //
      return 0.5 * ( Ostap::Math::trace     ( g2 * c1      ) - N  + 
                     ROOT::Math::Similarity ( g2 , v2 - v1 ) + 
                     std::log ( det2 / det1 )              ) ;
    }
    // ========================================================================
    /** Get the symmetrized Kullback-Leibler divergency,
     *  aka Jeffrey's divergence, for two objects 
     *  @see https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence
     *  \f[ f(v_1, C_1, v_2, C_2) = 
     *  (v_1-v_2)^{T} \left  ( C_1^{-1} + C_2^{-1} \right) (v_1 - v_2)  
     *   + Sp \left ( C_1 - C_2 \right ) 
     *    \times  \left ( C_2^{-1} - C_1^{-1} \right ) \f]
     *  @param v1 the first  data vector 
     *  @param c1 the first  covariance matrix 
     *  @param v2 the second data vector 
     *  @param c2 the second covariance matrix 
     *  @return Symmetrised Kullback-Leibler/Jeffrey divergency, or -999 
     */
    template <unsigned int N, typename SCALAR>
    inline double jeffrey 
    ( const ROOT::Math::SVector<SCALAR,N>&                                    v1 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c1 , 
      const ROOT::Math::SVector<SCALAR,N>&                                    v2 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c2 )
    {
      /// the actual type of covariance matrix
      typedef typename ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> > COV ;
      ///
      static const double s_bad = -999 ;
      /// try to invert matrices 
      COV g1 { c1 } ;
      if  ( !g1.InvertChol () ) { return s_bad ; }
      COV g2 { c2 } ;
      if  ( !g2.InvertChol () ) { return s_bad ; }
      ///
      return 0.5 * ( ROOT::Math::Similarity ( g1 + g2   , v2 - v1   ) + 
                     Ostap::Math::trace     ( g2 * c1 ) +
                     Ostap::Math::trace     ( g1 * c2 ) ) - N ;
    }
    // ========================================================================
    /** Get Jensen-Shannon divergency
     */
    template <unsigned int N, typename SCALAR>
    inline double jensen_shannon
    ( const ROOT::Math::SVector<SCALAR,N>&                                    v1 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c1 , 
      const double                                                            n1 ,  
      const ROOT::Math::SVector<SCALAR,N>&                                    v2 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c2 , 
      const double                                                            n2 ) 
    {
      //
      static const double s_bad = -999 ;
      //
      if ( n1 <= 1.0 || n2 <= 1.0 ) { return s_bad ; }
      //
      const double w1 = ( n1 - 1.0 ) / ( n1 + n2 - 2.0 ) ;
      const double w2 = ( n2 - 1.0 ) / ( n1 + n2 - 2.0 ) ;
      //
      /// the actual type of covariance matrix
      typedef typename ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> > COV ;
      /// the actual type of vector 
      typedef typename ROOT::Math::SVector<SCALAR,N>                                    VCT ;
      //
      const double w12 = w1 * w2 ;
      //
      const VCT v { w1 * v1 + w2 * v2 } ;
      COV       c { w1 * c1 + w2 * c2 } ;
      const VCT d { v2 - v1 } ;
      for ( unsigned int i = 0 ; i < N ; ++i )
      {
        const double di    = d [ i ]   ;
        const double diw   = di  * w12 ;
        c ( i , i )       += diw * di  ;
        for ( unsigned int j = i + 1 ; j < N ; ++j ) { c ( i , j ) +=  diw * d [ j ]  ; }
      }
      ///
      const double kl1 = kullback_leibler ( v1 , c1 , v , c ) ;
      const double kl2 = kullback_leibler ( v2 , c2 , v , c ) ;
      //
      return ( s_bad == kl1 || s_bad == kl2 ) ? s_bad : kl1 + kl2 ;
    }
    // ========================================================================
    /** Get Jensen-Shannon divergency
     */
    template <unsigned int N, typename SCALAR>
    inline double jensen_shannon
    ( const ROOT::Math::SVector<SCALAR,N>&                                    v1 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c1 , 
      const ROOT::Math::SVector<SCALAR,N>&                                    v2 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c2 )
    { return jensen_shannon ( v1 , v1 , 2.0 , v2 , c2 , 2.0 ) ; }
    // ========================================================================
    /** get Mahalanobis' distance
     *  \f[ D_M = \sqrt { (v_2-v_1)^T S^{-1} (v_2-v1) } \f]  , whre
     *  \f$ S = \frac{ (n_1-1)\Sigma_1 + ( n_2 -1 ) Sigma_2} {n_1+n_2-2} \f$ 
     *
     *  @param v1 (INPUT) the first data vector
     *  @param c1 (INPUT) the covariance matrix for the first data vector
     *  @param n1 (INPUT) sum of weights fot the first data vector 
     *  @param v2 (INPUT) the second data vector
     *  @param c2 (INPUT) the covariance matrix for the second data vector
     *  @param n2 (INPUT) sum of weights fot the second data vector 
     *
     *  @author Vanya BELYUAEV Ivan.Belyaev@itep.ru
     *  @date 2023-03-07
     */
    template <unsigned int N, typename SCALAR>
    inline double mahalanobis
    ( const ROOT::Math::SVector<SCALAR,N>&                                    v1 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c1 , 
      const double                                                            n1 ,  
      const ROOT::Math::SVector<SCALAR,N>&                                    v2 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c2 , 
      const double                                                            n2 ) 
    {
      //
      static const double s_bad = -999 ;
      //
      if ( n1 <= 1.0 || n2 <= 1.0 ) { return s_bad ; } 
      /// the actual type of covariance matrix
      typedef typename ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> > COV ;
      //
      const double w1 = ( n1 - 1.0 ) / ( n1 + n2 - 2.0 ) ;
      const double w2 = ( n2 - 1.0 ) / ( n1 + n2 - 2.0 ) ;
      //
      // Specialization for N = 1 to avoid Cling warnings and overhead
      if constexpr ( N == 1 ) 
      {
        const double si_val = w1 * c1 ( 0 , 0 ) + w2 * c2 ( 0 , 0 );
        if ( si_val <= 0 ) { return s_bad ; }

        const double diff   = v1 [ 0 ] - v2 [ 0 ];
        const double result = ( diff * diff ) / si_val;
        return 0 <= result ? std::sqrt ( result ) : s_bad ;
      }
      //      
      // pooled covariace matrix 
      COV si { w1 * c1 + w2 * c2 } ;
      if ( !si.InvertChol () ) { return s_bad ; }      
      //
      const double result = ROOT::Math::Similarity ( si , v1 - v2 ) ;
      return 0 <= result ? std::sqrt ( result ) : s_bad ;
    }
    // ========================================================================
    /** get Mahalanobis' distance
     *  \f[ D_M = \sqrt { (v_2-v_1)^T S^{-1} (v_2-v1) } \f]  , whre
     *  \f$ S = \frac{ (n_1-1)\Sigma_1 + ( n_2 -1 ) Sigma_2} {n_1+n_2-2} \f$ 
     *
     *  @param v1 (INPUT) the first data vector
     *  @param c1 (INPUT) the covariance matrix for the first data vector
     *  @param v2 (INPUT) the second data vector
     *  @param c2 (INPUT) the covariance matrix for the second data vector
     *
     *  @author Vanya BELYUAEV Ivan.Belyaev@itep.ru
     *  @date 2023-03-07
     */
    template <unsigned int N, typename SCALAR>
    inline double mahalanobis
    ( const ROOT::Math::SVector<SCALAR,N>&                                    v1 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c1 , 
      const ROOT::Math::SVector<SCALAR,N>&                                    v2 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c2 )
    { return mahalanobis ( v1 , c1 , 2.0 , v2 , c2 , 2.0 ) ; }
    // ========================================================================    
    /** get Hotelling's t-squared statistics 
     *  @see https://en.wikipedia.org/wiki/Hotelling%27s_T-squared_distribution#Two-sample_statistic
     *  \f[ t^2 = \frac{n_1 n_2}{n_1+n_2} \left(v_1-v_2\right)^T \Sigma^{-1} \left( v1-v2) \sim
     *   T^2 ( p , n_1 + n_2 -2 \f]
     *
     *  @param v1 (INPUT) the first data vector
     *  @param c1 (INPUT) the covariance matrix for the first data vector
     *  @param n1 (INPUT) sum of weights fot the first data vector 
     *  @param v2 (INPUT) the second data vector
     *  @param c2 (INPUT) the covariance matrix for the second data vector
     *  @param n2 (INPUT) sum of weights fot the second data vector 
     *
     *  @author Vanya BELYUAEV Ivan.Belyaev@itep.ru
     *  @date 2023-03-07
     */
    template <unsigned int N, typename SCALAR>
    inline double hotelling 
    ( const ROOT::Math::SVector<SCALAR,N>&                                    v1 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c1 , 
      const double                                                            n1 ,  
      const ROOT::Math::SVector<SCALAR,N>&                                    v2 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c2 , 
      const double                                                            n2 ) 
    {
      //
      static const double s_bad = -999 ;
      //
      if ( n1 <= 1.0 || n2 <= 1.0 ) { return s_bad ; } 
      /// the actual type of covariance matrix
      typedef typename ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> > COV ;
      //
      const double w1 = ( n1 - 1.0 ) / ( n1 + n2 - 2.0 ) ;
      const double w2 = ( n2 - 1.0 ) / ( n1 + n2 - 2.0 ) ;
      //
      // Specialization for N = 1 to avoid Cling warnings and overhead
      if constexpr ( N == 1 ) 
      {
        const double si_val = w1 * c1 ( 0 , 0 ) + w2 * c2 ( 0 , 0 ) ;
        if ( si_val <= 0 ) { return s_bad ; }

        const double diff = v1 [ 0 ] - v2 [ 0 ];
        const double quad = ( diff * diff ) / si_val;
        return ( n1 * n2 / ( n1 + n2 ) ) * quad;
      }
      //
      // pooled covariace matrix 
      COV si { w1 * c1 + w2 * c2 } ;
      if ( !si.InvertChol () ) { return s_bad ; }      
      //
      return n1 * n2 / ( n1 + n2 ) * ROOT::Math::Similarity ( si , v1 - v2 ) ;
    }
    // ========================================================================
    /** get Bhattacharyya's distance
     *  \f[ D_B = \frac{1}{8} (v_2-v_1)^T\left(\frac{\Sigma_1+\Sigma_2}{2}\right)^{-1}(v_2-v_1)
     *   + \frac{1}{2} \log \left( \frac{ det (\frac{\Sigma_1+\Sigma_2}{2}) }
     *   { \sqrt{ det \sigma_1 det \Sigma_2 }} \right) \f] 
     *
     *  @param v1 (INPUT) the first data vector
     *  @param c1 (INPUT) the covariance matrix for the first data vector
     *  @param v2 (INPUT) the second data vector
     *  @param c2 (INPUT) the covariance matrix for the second data vector
     *
     *  @author Vanya BELYUAEV Ivan.Belyaev@itep.ru
     *  @date 2023-03-07
     */
    template <unsigned int N, typename SCALAR>
    inline double bhattacharyya 
    ( const ROOT::Math::SVector<SCALAR,N>&                                    v1 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c1 ,
      const double                                                            n1 , 
      const ROOT::Math::SVector<SCALAR,N>&                                    v2 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c2 ,
      const double                                                            n2 )
    {
      //
      static const double s_bad = -999 ;
      //
      if ( n1 <= 1.0 || n2 <= 1.0 ) { return s_bad ; } 
      //
      const double w1 = ( n1 - 1.0 ) / ( n1 + n2 - 2.0 ) ;
      const double w2 = ( n2 - 1.0 ) / ( n1 + n2 - 2.0 ) ;
      //
      if constexpr ( N == 1 ) 
      {
        const double diff   = v2[0] - v1[0];
        //
        const double c100   = c1 ( 0 , 0 ) ;
        const double c200   = c2 ( 0 , 0 ) ;
        //
        const double si_val = w1 * c100 + w2 * c200 ;
        
        if ( si_val <= 0 || c100 <= 0 || c200 <= 0 ) { return s_bad ; }
        //
        const double quad = ( diff * diff ) / si_val;
        
        // using pooled variance depending on definition, keeping consistent with matrix logic       
        return 0.5 * w1 * w2 * quad + 0.5 * std::log ( si_val / std::sqrt ( c100 * c200 ) ) ;
      }
      /// the actual type of covariance matrix
      typedef typename ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> > COV ;
      // 
      // get "average"  covariance matrix 
      COV si { w1 * c1 + w2 * c2 } ;
      //
      // (1) calculate the determinant
      SCALAR det  = 1 ;
      if ( !si.Det2 ( det ) || 0 >= det ) { return s_bad ; }
      // (2) invert it! 
      if ( !si.InvertChol ()            ) { return s_bad ; }      
      //
      SCALAR det1 = 1 ;
      if ( !c1.Det2 ( det1 ) || det1 < 0 ) { return s_bad ; }
      SCALAR det2 = 1 ;
      if ( !c2.Det2 ( det2 ) || det2 < 0 ) { return s_bad ; }
      //
      return
        0.5 * w1 * w2 * ROOT::Math::Similarity ( si , v2 - v1      ) + 
        0.5           * std::log ( det / std::sqrt ( det1 * det2 ) ) ;
    }
    // ========================================================================
    /** get Bhattacharyya's distance
     *  \f[ D_B = \frac{1}{8} (v_2-v_1)^T\left(\frac{\Sigma_1+\Sigma_2}{2}\right)^{-1}(v_2-v_1)
     *   + \frac{1}{2} \log \left( \frac{ det (\frac{\Sigma_1+\Sigma_2}{2}) }
     *   { \sqrt{ det \sigma_1 det \Sigma_2 }} \right) \f] 
     *
     *  @param v1 (INPUT) the first data vector
     *  @param c1 (INPUT) the covariance matrix for the first data vector
     *  @param v2 (INPUT) the second data vector
     *  @param c2 (INPUT) the covariance matrix for the second data vector
     *
     *  @author Vanya BELYUAEV Ivan.Belyaev@itep.ru
     *  @date 2023-03-07
     */
    template <unsigned int N, typename SCALAR>
    inline double bhattacharyya 
    ( const ROOT::Math::SVector<SCALAR,N>&                                    v1 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c1 ,
      const ROOT::Math::SVector<SCALAR,N>&                                    v2 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c2 )
      { return bhattachariyya  ( v1 , c1 , 2.0 , v2 , c2 , 2.0 ) ; }    
    // ========================================================================     
    /** get Wasserstain' distance
     *  \f[ W^2 = (v_2-v_1)^T \Sigma ^{-1} (v_2-v1) } 
     *          + tr \left( \Sigma_1 + \Sigma_2
     *          - 2  \left( \Sigma_2^{1/2}\Sigma_1 \Sigma_2^{1/2}\right)^{1/2} \right) \f]
     *
     *  @param v1 (INPUT) the first data vector
     *  @param c1 (INPUT) the covariance matrix for the first data vector
     *  @param n1 (INPUT) sum of weights fot the first data vector 
     *  @param v2 (INPUT) the second data vector
     *  @param c2 (INPUT) the covariance matrix for the second data vector
     *  @param n2 (INPUT) sum of weights fot the second data vector 
     *
     *  @author Vanya BELYUAEV Ivan.Belyaev@itep.ru
     *  @date 2023-03-07
     */
    template <unsigned int N, typename SCALAR>
    inline double wasserstein 
    ( const ROOT::Math::SVector<SCALAR,N>&                                    v1 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c1 , 
      const double                                                            n1 ,  
      const ROOT::Math::SVector<SCALAR,N>&                                    v2 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c2 , 
      const double                                                            n2 ) 
    {
      //
      static const double s_bad = -999 ;
      //
      if ( n1 <= 1.0 || n2 <= 1.0 ) { return s_bad ; } 
      /// the actual type of covariance matrix
      typedef typename ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >   COV ;
      /// the actual type of cholesky matrix 
      typedef typename ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepStd<SCALAR,N,N> > CHOL ;
      /// Vector of eigenvalues 
      typedef typename ROOT::Math::SVector<SCALAR,N> EVCT ;      
      //
      const double w1 = ( n1 - 1.0 ) / ( n1 + n2 - 2.0 ) ;
      const double w2 = ( n2 - 1.0 ) / ( n1 + n2 - 2.0 ) ;
      //
      if constexpr ( N == 1 ) 
      {
        const double diff = v1 [ 0 ] - v2 [ 0 ];
        //
        const double c100   = c1 ( 0 , 0 ) ;
        const double c200   = c2 ( 0 , 0 ) ;
        //
        const double si      = w1 * c100 + w2 * c200 ;
        if ( si <= 0 ) { return s_bad ; }
        //
        const double result1 = ( diff * diff ) / si ;
        const double cross   = std::sqrt ( c100 * c200 ) ;
        //
        const double wsq = result1 - 2.0 * cross + c100 + c200 ;
        return 0 <= wsq ? wsq : 0.0 ;
      }
      //      
      // pooled covariace matrix 
      COV si { w1 * c1 + w2 * c2 } ;
      if ( !si.InvertChol () ) { return s_bad ; }      
      //
      /// the first term 
      const double result1 = ROOT::Math::Similarity ( si , v1 - v2 ) ;
      ///
      CHOL L {} ;
      if ( !cholesky ( c1 , L ) ) { return s_bad ; }
      /// helper matrix M 
      const COV M { ROOT::Math::Similarity ( L , c2 ) } ;
      ///
      EVCT V {} ;
      const Ostap::Math::GSL::EigenSystem eigen {} ;
      if ( !eigen.eigenValues ( M , V ).isSuccess() ) { return s_bad ; }
      //
      double result2 = 0 ;
      for ( unsigned int i = 0 ; i < N ; ++i )
      {
        const double ev = V [ i ] ;
        if ( 0 < ev ) { result2 += std::sqrt ( ev ) ; }
      }
      //
      const double wsq = result1 - 2 * result2 +
        Ostap::Math::trace ( c1 ) +
        Ostap::Math::trace ( c2 ) ;
      //
      return 0 <= wsq ? wsq : 0.0 ;
    }
    // ========================================================================
    /** get Wasserstain' distance
     *  \f[ W^2 = (v_2-v_1)^T \Sigma ^{-1} (v_2-v1) } 
     *          + tr \left( \Sigma_1 + \Sigma_2
     *          - 2  \left( \Sigma_2^{1/2}\Sigma_1 \Sigma_2^{1/2}\right)^{1/2} \right) \f]
     *
     *  @param v1 (INPUT) the first data vector
     *  @param c1 (INPUT) the covariance matrix for the first data vector
     *  @param n1 (INPUT) sum of weights fot the first data vector 
     *  @param v2 (INPUT) the second data vector
     *  @param c2 (INPUT) the covariance matrix for the second data vector
     *  @param n2 (INPUT) sum of weights fot the second data vector 
     *
     *  @author Vanya BELYUAEV Ivan.Belyaev@itep.ru
     *  @date 2023-03-07
     */
    template <unsigned int N, typename SCALAR>
    inline double wasserstein 
    ( const ROOT::Math::SVector<SCALAR,N>&                                    v1 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c1 , 
      const ROOT::Math::SVector<SCALAR,N>&                                    v2 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c2 )
      { return wasserstein ( v1 , c1 , 2.0 , v2 , c2 , 2.0 ) ; }       
    // ========================================================================

    // ========================================================================
    /** get Hellinger's distance (squared) between two multivariate Gaussians
     *  @see https://en.wikipedia.org/wiki/Hellinger_distance
     *  \f[ H^2 = 1 - \frac{(\det C_1)^{1/4} (\det C_2)^{1/4}}{(\det \Sigma_{avg})^{1/4}} 
     *      \exp \left( - \frac{1}{8} (v_2 - v_1)^T \Sigma_{avg}^{-1} (v_2 - v_1) \right) \f]
     *
     *  @param v1 (INPUT) the first data vector
     *  @param c1 (INPUT) the covariance matrix for the first data vector
     *  @param n1 (INPUT) sum of weights for the first data vector 
     *  @param v2 (INPUT) the second data vector
     *  @param c2 (INPUT) the covariance matrix for the second data vector
     *  @param n2 (INPUT) sum of weights for the second data vector 
     *
     *  @author Vanya BELYUAEV Ivan.Belyaev@itep.ru
     *  @date 2026-06-06
     */
    template <unsigned int N, typename SCALAR>
    inline double hellinger 
    ( const ROOT::Math::SVector<SCALAR,N>&                                    v1 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c1 , 
      const double                                                            n1 ,  
      const ROOT::Math::SVector<SCALAR,N>&                                    v2 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c2 , 
      const double                                                            n2 ) 
    {
      //
      static const double s_bad = -999 ;
      //
      if ( n1 <= 1.0 || n2 <= 1.0 ) { return s_bad ; } 
      //
      const double w1 = ( n1 - 1.0 ) / ( n1 + n2 - 2.0 ) ;
      const double w2 = ( n2 - 1.0 ) / ( n1 + n2 - 2.0 ) ;
      //

      if constexpr ( N == 1 ) 
      {
        const double diff = v2[0] - v1[0];
        
        const double c100 = c1 ( 0 , 0 ) ;
        const double c200 = c2 ( 0 , 0 ) ;
        
        const double si_val = w1 * c100 + w2 * c200 ;
        if ( si_val <= 0 || c100 < 0 || c200 < 0 ) { return s_bad ; }

        const double quad      = w1 * w2 * ( diff * diff ) / si_val;
        const double exp_term  = std::exp ( -0.125 * quad );
        
        const double det_ratio = std::sqrt ( std::sqrt ( c100 * c200 / si_val ) ) ;        
        const double h_sq      = 1.0 - det_ratio * exp_term ;

        return ( h_sq >= 0.0 ) ? h_sq : 0.0 ;
      }
      
      /// the actual type of covariance matrix
      typedef typename ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> > COV ;
      
      // pooled covariance matrix 
      COV si { w1 * c1 + w2 * c2 } ;
      //
      // (1) calculate the determinant of pooled matrix
      SCALAR det  = 1 ;
      if ( !si.Det2 ( det ) || 0 >= det ) { return s_bad ; }
      // (2) invert it for quadratic form
      if ( !si.InvertChol ()            ) { return s_bad ; }      
      //
      SCALAR det1 = 1 ;
      if ( !c1.Det2 ( det1 ) || det1 < 0 ) { return s_bad ; }
      SCALAR det2 = 1 ;
      if ( !c2.Det2 ( det2 ) || det2 < 0 ) { return s_bad ; }
      //
      // Calculate exponential term via quadratic form and weights
      const double quad      = w1 * w2 * ROOT::Math::Similarity ( si , v2 - v1 ) ;
      const double exp_term  = std::exp ( -0.125 * quad ) ;
      
      // Calculate determinant ratio: (det1 * det2)^(1/4) / det^(1/4)
      const double det_ratio = std::sqrt ( std::sqrt ( det1 * det2 / det ) ) ;
      
      const double h_sq      = 1.0 - det_ratio * exp_term ;
      
      // Hellinger squared distance is bound in [0, 1]
      return ( h_sq >= 0.0 ) ? h_sq : 0.0 ;
    }
    // ========================================================================
    /** get Hellinger's distance (squared) between two multivariate Gaussians
     *  @param v1 (INPUT) the first data vector
     *  @param c1 (INPUT) the covariance matrix for the first data vector
     *  @param v2 (INPUT) the second data vector
     *  @param c2 (INPUT) the covariance matrix for the second data vector
     */
    template <unsigned int N, typename SCALAR>
    inline double hellinger 
    ( const ROOT::Math::SVector<SCALAR,N>&                                    v1 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c1 , 
      const ROOT::Math::SVector<SCALAR,N>&                                    v2 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c2 )
    { return hellinger ( v1 , c1 , 2.0 , v2 , c2 , 2.0 ) ; }
    // ========================================================================    
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
#endif // OSTAP_MATRIX_DISTANCES_H
// ============================================================================
//                                                                      The END 
// ============================================================================



