// ============================================================================
#ifndef OSTAP_MATRIXUTILS2_H
#define OSTAP_MATRIXUTILS2_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL 
// ============================================================================
#include <algorithm>
#include <numeric>
#include <functional>
#include <utility>
#include <cmath>
#include <array>
// ============================================================================
// ROOT
// ============================================================================
#include "Math/SMatrix.h"
#include "Math/SVector.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Constants.h"
#include "Ostap/Epsilon.h"
#include "Ostap/Math.h"
#include "Ostap/Norms.h"
#include "Ostap/MatrixAsBuffer.h"
#include "Ostap/LinAlg.h"
#include "Ostap/EigenSystem.h"
#include "Ostap/StatusCode.h"
#include "Ostap/MatrixUtils.h"
// ============================================================================
/** @file Ostap/MatrixUtils2.h
 *  The collection of functions for manipulation with matrices and vectors.
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 */
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================

    // ========================================================================
    // spectral matrix norms
    // ========================================================================

    // ========================================================================
    /** @brief Compute spectral matrix norm (maximum singular value) 
     *  @param m (INPUT) Input general matrix
     *  @return Spectral norm value
     */
    template <typename T,
              unsigned int D1,
              unsigned int D2,
              typename R>            
    inline T norm_spectral
    ( const ROOT::Math::SMatrix<T,D1,D2,R>& m )
    {
      //
      // 1. Wrap SMatrix data into a GSL matrix via Ostap::Utils::Buffer
      Ostap::Math::GSL::Matrix A ( D1 , D2 ,  Ostap::Utils::buffer ( m  ) ) ;
      //
      // 2. Allocate structures for SVD[cite: 1]
      const  std::size_t K = std::min ( D1 , D2 ) ;
      Ostap::Math::GSL::Vector S ( K      ) ; // singular values 
      Ostap::Math::GSL::Matrix U ( D1 , K ) ; 
      Ostap::Math::GSL::Matrix V ( D2 , K ) ;
      //
      // 3. Compute SVD using the GSL module
      const Ostap::StatusCode sc = Ostap::Math::GSL::SVD( A , S , U , V ) ;
      if ( sc.isFailure () ) { return static_cast<T> ( Ostap::v_INVALID_NORM ) ; }
      //      
      // 4. Spectral norm is the maximum singular value
      return static_cast<T> ( norm_max ( S ) ) ;
    }
    
    // ========================================================================
    /** @brief Compute spectral matrix norm (maximum singular value) for symmetric matrices 
     *  @param m (INPUT) Input general matrix
     *  @return Spectral norm value
     */
    template <typename T,
              unsigned int D>
    inline T norm_spectral
    ( const ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >& m  )
    {
      // allocate eigensystem 
      Ostap::Math::GSL::EigenSystem eigen_system {}  ;      
      ROOT::Math::SVector<T, D>     values;
      //
      // compute eigenvalues 
      const Ostap::StatusCode sc = eigen_system.eigenValues ( m , values , false ) ;
      if ( sc.isFailure () ) { return static_cast<T> ( Ostap::v_INVALID_NORM ) ; }
      //
      //  Spectral norm is the maximum singular value
      return static_cast<T> ( norm_max ( values ) ) ;
    }

    // ========================================================================
    /** @brief Compute nuclear norm (sum of singular value) 
     *  @param m (INPUT) Input general matrix
     *  @return Nuclear norm value
     */
    template <typename T,
              unsigned int D1,
              unsigned int D2,
              typename R>            
    inline T norm_nuclear 
    ( const ROOT::Math::SMatrix<T,D1,D2,R>& m )
    {
      //
      // 1. Wrap SMatrix data into a GSL matrix via Ostap::Utils::Buffer
      Ostap::Math::GSL::Matrix A ( D1 , D2 ,  Ostap::Utils::buffer ( m ) ) ;
      //
      // 2. Allocate structures for SVD[cite: 1]
      const  std::size_t K = std::min ( D1 , D2 ) ;
      Ostap::Math::GSL::Vector S ( K      ) ; // singular values 
      Ostap::Math::GSL::Matrix U ( D1 , K ) ; 
      Ostap::Math::GSL::Matrix V ( D2 , K ) ;
      //
      // 3. Compute SVD using the GSL module
      const Ostap::StatusCode sc = Ostap::Math::GSL::SVD( A , S , U , V ) ;
      if ( sc.isFailure () ) { return static_cast<T> ( Ostap::v_INVALID_NORM ) ; }
      //      
      // 4. Nuclear norm is a sum of singular values 
      return norm_L1 ( S );
    }

    // ========================================================================
    /** @brief Compute nuclear norm (sum of singular value) foe symmetric matrix 
     *  @param m (INPUT) Input general matrix
     *  @return Nuclear norm value
     */
    template <typename T,
              unsigned int D>
    inline T norm_nuclear
    ( const ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >& m )
    {
      // allocate eigensystem 
      Ostap::Math::GSL::EigenSystem eigen_system {}  ;      
      ROOT::Math::SVector<T, D>     values;
      //
      // compute eigenvalues 
      const Ostap::StatusCode sc = eigen_system.eigenValues ( m , values , false ) ;
      if ( sc.isFailure () ) { return static_cast<T> ( Ostap::v_INVALID_NORM ) ; }
      //
      //  Nuclear norm is a sum of singular value = sum of eigenvalues moduli 
      return static_cast<T> ( norm_L1( values ) ) ;
    }

    // ========================================================================
    /** @brief Compute Schatten matrix norm \f$ \left(\Sum \left|\sigma_i\right|^p\right)^{1/p}\f$
     *  @param m (INPUT) Input general matrix
     *  @param p (INPUT) power parameter 
     *  @return Spectral norm value
     */
    template <typename T,
              unsigned int D1,
              unsigned int D2,
              typename R>            
    inline T norm_schatten 
    ( const ROOT::Math::SMatrix<T,D1,D2,R>& m , 
      const double                          p = 2 )
    {
      //
      // 1. Wrap SMatrix data into a GSL matrix via Ostap::Utils::Buffer
      Ostap::Math::GSL::Matrix A ( D1 , D2 ,  Ostap::Utils::buffer ( m ) ) ;
      //
      // 2. Allocate structures for SVD[cite: 1]
      const  std::size_t K = std::min ( D1 , D2 ) ;
      Ostap::Math::GSL::Vector S ( K      ) ; // singular values 
      Ostap::Math::GSL::Matrix U ( D1 , K ) ; 
      Ostap::Math::GSL::Matrix V ( D2 , K ) ;
      //
      // 3. Compute SVD using the GSL module
      const Ostap::StatusCode sc = Ostap::Math::GSL::SVD ( A , S , U , V ) ;
      if ( sc.isFailure () ) { return static_cast<T> ( Ostap::v_INVALID_NORM ) ; }
      //      
      // 4. Schatten' norm is Lp norm fro vector of singular values 
      return norm_Lp ( S , p );
    }

    // ========================================================================
    /** @brief Compute Schatten matrix norm \f$ \left(\Sum \left|\sigma_i\right|^p\right)^{1/p}\f$
     *  @param m (INPUT) Input general matrix
     *  @param p (INPUT) power parameter 
     *  @return Schatten norm value
     */
    template <typename T,
              unsigned int D>
    inline T norm_schatten 
    ( const ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >& m    , 
      const double                                                  p = 2 )
    {
      // allocate eigensystem 
      Ostap::Math::GSL::EigenSystem eigen_system {}  ;      
      ROOT::Math::SVector<T, D>     values;
      //
      // compute eigenvalues 
      const Ostap::StatusCode sc = eigen_system.eigenValues ( m , values , false ) ;
      if ( sc.isFailure () ) { return static_cast<T> ( Ostap::v_INVALID_NORM ) ; }
      //
      //  Schatten' norm is a Lp norm of vector of eigenvalues 
      return static_cast<T> ( norm_Lp ( values , p ) ) ;
    }

    // ========================================================================
    /** Get the rank of general matrix (SMatrix) via GSL rank function.
     *  @param  m   (INPUT) input matrix
     *  @param  eps (INPUT) tolerance for rank determination (negative for default)
     *  @return matrix rank
     */
    template <typename     T ,
              unsigned int D1,
              unsigned int D2,
              typename     R>
    inline std::size_t rank
    ( const ROOT::Math::SMatrix<T, D1, D2, R>& m    ,
      const double                             eps  = epsilon_v<T> )
    {      
      Ostap::Math::GSL::Matrix A ( D1, D2, Ostap::Utils::buffer ( m ) ) ;
      return Ostap::Math::rank ( A , static_cast<double> ( eps ) ) ;
    }
    
    // ========================================================================
    /** Get the rank of symmetrical matrix (SMatrix) via EigenSystem.
     *  @param  m   (INPUT) input symmetrical matrix
     *  @param  eps (INPUT) tolerance for rank determination (negative for default)
     *  @return matrix rank
     */
    template <typename T, unsigned int D>
    inline std::size_t rank
    ( const ROOT::Math::SMatrix<T, D, D, ROOT::Math::MatRepSym<T, D>>& m   ,
      const T                                                          eps = epsilon_v<T> )
    {
      //
      Ostap::Math::GSL::EigenSystem eigen_system {} ;
      ROOT::Math::SVector<T, D>     values          ;
      //
      const Ostap::StatusCode sc = eigen_system.eigenValues ( m , values , false ) ;
      if ( sc.isFailure() ) { return 0 ; }
      //
      return norm_L0 ( values ,  eps ) ;
    }

    // ========================================================================
    /** Compute Moore-Penrose Pseudoinverse using SVD: A^+ = V * Sigma^+ * U^T
     *  @param a     (INPUT)  Input matrix A (m x n)
     *  @param a_pinv(OUTPUT) Pseudoinverse matrix A^+ (n x m)
     *  @param tol   (INPUT)  Tolerance for zeroing small singular values (< 0 for default)
     *  @return status code
     */
    template <typename T,
              unsigned int D1,
              unsigned int D2,
              typename     R>
    inline Ostap::StatusCode PINV 
    ( const ROOT::Math::SMatrix<T, D1, D2, R>& a      ,
      ROOT::Math::SMatrix<T, D2, D1, R>&       a_pinv ,      
      const T                                  eps    = epsilon_v<T> )
    {
      
      // Wrap SMatrix data into a GSL matrix
      Ostap::Math::GSL::Matrix A ( D1, D2 , Ostap::Utils::buffer ( a ) );
      
      // Allocate structures for SVD
      std::size_t K = std::min ( D1, D2 );
      Ostap::Math::GSL::Vector S ( K );
      Ostap::Math::GSL::Matrix U ( D1, K );
      Ostap::Math::GSL::Matrix V ( D2, K );
      
      // Compute SVD: A = U * S * V^T
      const Ostap::StatusCode sc1 = Ostap::Math::GSL::SVD ( A , S, U, V  ) ;
      if ( sc1.isFailure() ) { return sc1 ; }

      //
      const double tol = ( 0 < eps ? eps : epsilon_v<T> ) * S ( 0 ) * std::max ( D1 , D2 ) ;
      //
      // Compute pseudo-inverse: A^+ = V * S^+ * U^T
      for ( std::size_t k = 0 ; k < K ; ++ k ) 
      {
        const double si = S ( k ) ;
        const double v  = ( si && tol < std::abs ( si ) ) ? 1.0/si : 0.0 ;
        S.set ( k , v ) ;
      }
      //
      // final result 
      Ostap::Math::GSL::Matrix r ( D2, D1 ) ;       
      const Ostap::StatusCode sc2 = Ostap::Math::GSL::MDM ( V , false , S , U , true , r ) ;
      if ( sc2.isFailure() ) { return sc2 ; }
      //
      // Convert the result to S-matrix            
      for ( std::size_t i = 0 ; i < D2 ; ++i )
      { for ( std::size_t j = 0; j < D1; ++j )
        { a_pinv ( i, j ) = r.get ( i , j ) ;} }
      //
      return Ostap::StatusCode::SUCCESS ;
    }

    // ========================================================================
    /** Compute Moore-Penrose Pseudoinverse using SVD: A^+ = V * Sigma^+ * U^T
     *  @param a     (INPUT)  Input matrix A (m x n)
     *  @param a_pinv(OUTPUT) Pseudoinverse matrix A^+ (n x m)
     *  @param tol   (INPUT)  Tolerance for zeroing small singular values (< 0 for default)
     *  @return status code
     */
    template <typename T,
              unsigned int D>
    inline Ostap::StatusCode PINV 
    ( const ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T, D>>& a      ,
      ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T, D>>&       a_pinv ,      
      const T                                                       eps    = epsilon_v<T> )
    {
      //
      Ostap::Math::GSL::EigenSystem                            eigen_system ( D ) ;
      
      // Wrap SMatrix data into a GSL matrix
      Ostap::Math::GSL::Matrix A ( D , D , Ostap::Utils::buffer ( a ) ) ;

      // eigenvalues 
      Ostap::Math::GSL::Vector S ( D ) ;
      
      // eigenvectors 
      Ostap::Math::GSL::Matrix V ( D , D ) ; 
      //
      const Ostap::StatusCode sc = eigen_system.eigenVectors ( A , S , V , false ) ; 
      if ( sc.isFailure() ) { return sc ; } 

      const double max_val  = norm_max ( D ) ;
      const double tol = ( 0 < eps ? eps : epsilon_v<T> ) * max_val * D ;
      
      // Compute pseudo-inverse of eigenvalues: S^+
      for ( std::size_t k = 0; k < D; ++k )  
      {
        const double si = S.get( k );
        const double s  = ( si && tol < std::abs ( si ) ) ? 1.0 / si : 0.0 ;
        S.set ( k , s ) ;
      }
      
      // final result: A^+ = V * S^+ * V^T      
      Ostap::Math::GSL::Matrix r ( D, D );      
      const Ostap::StatusCode sc2 = Ostap::Math::GSL::MDM ( V , false, S, V, true, r );
      if ( sc2.isFailure() ) { return sc2; }
      
      // Convert the result to S-matrix            
      for ( std::size_t i = 0; i < D; ++i )
      { for ( std::size_t j = 0; j <= i ; ++j )
        { a_pinv ( i, j ) = r.get ( i, j ) ; } }
      //
      return Ostap::StatusCode::SUCCESS ;
    }
    
    // ========================================================================
    /** @brief Compute Variance Inflation Factors (VIF) for a covariance matrix.
     *
     *  Calculates the Variance Inflation Factor (VIF) vector \f$ \vec{v} \f$ 
     *  for a symmetric $D \times D$ covariance matrix \f$ \Sigma \f$:
     *  \f[
     *      v_i = \Sigma_{ii} \cdot (\Sigma^+)_{ii}
     *  \f]
     *  where \f$ \Sigma_{ii} \f$ is the variance of variable $i$, and 
     *  \f$ (\Sigma^+)_{ii} \f$ is the corresponding diagonal element of the 
     *  Moore-Penrose pseudoinverse matrix \f$ \Sigma^+ \f$.
     *
     *  @par Connection to Global Correlation Coefficient:
     *  In classical linear regression, the VIF of variable $i$ measures how much 
     *  the variance of the estimated regression coefficient is inflated due to 
     *  multicollinearity. It is strictly related to the **Global Correlation 
     *  Coefficient** \f$ R_i \f$ (the coefficient of determination when regressing 
     *  variable $i$ against all other $D-1$ variables):
     *  \f[
     *      v_i = \frac{1}{1 - R_i^2} \quad \Longleftrightarrow \quad R_i = \sqrt{1 - \frac{1}{v_i}}
     *  \f]
     *  - \f$ R_i = 0 \implies v_i = 1 \f$: Variable $i$ is orthogonal (uncorrelated) to all others.
     *  - \f$ R_i \to 1 \implies v_i \to \infty \f$: Variable $i$ is a linear combination of other variables.
     *
     *  @par Numerical Robustness & Fallback Architecture:
     *  1. **Fast Path**: Attempts in-place Cholesky decomposition (\f$ \Sigma = L L^T \f$). 
     *     This is $O(D^3)$ with minimal overhead for well-behaved, Positive Definite matrices.
     *  2. **Fallback Path**: If Cholesky fails (due to zero or negative eigenvalues 
     *     arising from exact linear dependencies or negative \f$sPlot\f$ event weights), 
     *     it gracefully falls back to spectral pseudoinversion (`PINV`). Zero and negative 
     *     eigenvalues (\f$\lambda_k \le \text{eps} \cdot \lambda_{\max}\f$) are zeroed out, 
     *     preventing division-by-zero crashes while preserving the subspace mapping.
     *  3. **Non-positive Variances**: Variables with \f$ \Sigma_{ii} \le 0 \f$ (constants or 
     *     severe \f$sPlot\f$ noise) are explicitly assigned `std::numeric_limits<T>::infinity()`, 
     *     marking them as primary targets for elimination.
     *
     *  @par Interpretation Thresholds:
     *  - \f$ v_i \approx 1 \f$: No collinearity.
     *  - \f$ 1 < v_i < 5 \f$: Moderate, acceptable correlation.
     *  - \f$ v_i > 10 \f$: High collinearity (\f$ R_i > 0.95 \f$), feature removal recommended.
     *  - \f$ v_i > 10^4 \text{ or } \infty \f$: Critical geometric degeneracy or constant feature.
     *
     *  @tparam T Data type (`double`, `float`).
     *  @tparam D Dimension of the covariance matrix.
     *  @param[in]  cov Input $D \times D$ symmetric covariance matrix \f$ \Sigma \f$.
     *  @param[out] vif Output $D$-dimensional vector containing VIF values.
     *  @param[in]  eps Numerical tolerance ratio for truncating small eigenvalues in `PINV`.
     *  @return `Ostap::StatusCode::SUCCESS` if computation completed successfully.
     */
    template <typename T, unsigned int D>
    inline Ostap::StatusCode VIF 
    ( const ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >& cov ,
      ROOT::Math::SVector<T,D>&                                     vif ,
      const T                                                       eps = epsilon_v<T> ) 
    {
      // (symmetric) matrix type 
      typedef typename ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >  MTRX ;
      //
      //  Inverse/pseudoinverse matrix 
      MTRX sp { cov };
      // (1) try Cholesky' inversion 
      if  ( !sp.InvertChol() )
      {
        // (2) if Cholesky fails - use Moore-Penrouse'
        const Ostap::StatusCode sc = PINV ( cov , sp , eps ) ;
        if ( sc.isFailure() ) { return sc ; } 
      }
      //
      // (3) finally calculate VIFs
      for  ( std::size_t i = 0 ; i < D ; ++i )
      {
        const T cii = cov ( i , i ) ;
        // (4) mark the pathological components with infinities 
        vif [ i ] = cii <= 0 ? std::numeric_limits<T>::infinity () : cii * sp ( i , i )  ;
      }
      return Ostap::StatusCode::SUCCESS ;
    }


    // ========================================================================
    // helper functions to allow proper operations in PyROOT
    // - need to bypass expressions  (no easy way to use them in PyROOT)
    // ========================================================================

    // ========================================================================
    // Vector operations 
    // ========================================================================
    template <class VECTOR>
    struct VctrOps ;
    
    // ========================================================================
    template <class T, unsigned int D>
    struct VctrOps < ROOT::Math::SVector<T,D> >
    {
      // a + b 
      static
      ROOT::Math::SVector<T,D> 
      add  ( const ROOT::Math::SVector<T,D>& a ,
             const ROOT::Math::SVector<T,D>& b ) { return a + b ; }
      // a + c 
      static
      ROOT::Math::SVector<T,D>
      add  ( const ROOT::Math::SVector<T,D>& a ,
             const double                    c ) { return a + c ; }
      // a - b 
      static
      ROOT::Math::SVector<T,D> 
      sub  ( const ROOT::Math::SVector<T,D>& a ,
             const ROOT::Math::SVector<T,D>& b ) { return a - b ; }
      // a - c 
      static
      ROOT::Math::SVector<T,D> 
      sub  ( const ROOT::Math::SVector<T,D>& a ,
             const double                    b ) { return a + b ; }
      // c - a  
      static
      ROOT::Math::SVector<T,D> 
      rsub ( const ROOT::Math::SVector<T,D>& a ,
             const double                    c ) { return c - a  ; }
    };
    
    // ========================================================================
    // matrix operations 
    // ========================================================================
    template <class MATRIX>
    struct MtrxOps ;

    /// generic matrices 
    template <class T, unsigned int D1, unsigned int D2>
    struct MtrxOps< ROOT::Math::SMatrix<T,D1,D2,ROOT::Math::MatRepStd<T,D1,D2> > > 
    {
      typedef ROOT::Math::SMatrix<T,D1,D2,ROOT::Math::MatRepStd<T,D1,D2> > MATRIX ;
      //
      // ======================================================================
      // m + m 
      static MATRIX add  ( const MATRIX& a ,
                           const MATRIX& b ) { return a + b ; }
      // m + c 
      static MATRIX add  ( const MATRIX& a ,
                           const double  c ) { return a + c ; }
      // ======================================================================
      // m - m 
      static MATRIX sub  ( const MATRIX& a ,
                           const MATRIX& b ) { return a - b ; }
      // m - c
      static MATRIX sub  ( const MATRIX& a ,
                           const double  c ) { return a - c ; }
      // ======================================================================
      // c - m
      static MATRIX rsub ( const MATRIX& a ,
                           const double  c ) { return c - a ; }
      // ======================================================================
    };

    /// square matrices 
    template <class T, unsigned int D>
    struct MtrxOps<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepStd<T,D,D> > > 
    {
      typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepStd<T,D,D> >      MATRIX ;
      typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >     SYMMATRIX ;
      //
      // ======================================================================
      // m + m
      static MATRIX add  ( const    MATRIX& a ,
                           const    MATRIX& b ) { return a + b ; }
      // m + s 
      static MATRIX add  ( const    MATRIX& a ,
                           const SYMMATRIX& b ) { return a + b ; }
      // m + c
      static MATRIX add  ( const    MATRIX& a ,
                           const    double  b ) { return a + b ; }
      // ======================================================================
      // m - m 
      static MATRIX sub  ( const    MATRIX& a ,
                           const    MATRIX& b ) { return a - b ; }
      // m - s 
      static MATRIX sub  ( const    MATRIX& a ,
                           const SYMMATRIX& b ) { return a - b ; }
      // m - c 
      static MATRIX sub  ( const    MATRIX& a ,
                           const    double  c ) { return a - c ; }
      // ======================================================================
      // s - m 
      static MATRIX rsub ( const    MATRIX& a ,
                           const SYMMATRIX& b ) { return b - a ; }
      // c - m  
      static MATRIX rsub ( const    MATRIX& a ,
                           const    double  c ) { return c - a ; }
      
      // ======================================================================
    };

    /// symmetric matrices
    template <class T, unsigned int D>
    struct MtrxOps<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > > 
    {
      typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >      MATRIX ;
      typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepStd<T,D,D> > GENMATRIX ;
      //
      // ======================================================================
      // s + s 
      static MATRIX    add  ( const    MATRIX& a ,
                              const    MATRIX& b ) { return a + b ; }
      // s + m      
      static GENMATRIX add  ( const    MATRIX& a ,
                              const GENMATRIX& b ) { return a + b ; }
      // s + c
      static MATRIX    add  ( const    MATRIX& a ,
                              const     double c ) { return a + c ; }      
      // ======================================================================
      // s - s 
      static MATRIX    sub  ( const    MATRIX& a ,
                              const    MATRIX& b ) { return a - b ; }
      // s - m      
      static GENMATRIX sub  ( const    MATRIX& a ,
                              const GENMATRIX& b ) { return a - b ; }
      // s - c
      static MATRIX    sub  ( const    MATRIX& a ,
                              const     double c ) { return a - c ; }
      // ======================================================================      
      // m - s      
      static GENMATRIX rsub ( const    MATRIX& a ,
                              const GENMATRIX& b ) { return b - a ; }
      // c - s 
      static MATRIX    rsub ( const    MATRIX& a ,
                              const     double c ) { return c - a ; }
      // ======================================================================      
    };

    
    // ========================================================================
    /// multiplication
    template <class OBJ1, class OBJ2>
    struct MultiplyOp ;
    // ========================================================================
    // vector * vector  
    template <class T, unsigned int D>
    struct MultiplyOp < ROOT::Math::SVector<T,D>, ROOT::Math::SVector<T,D> >
    {
      // dot:
      static
      double 
      dot ( const ROOT::Math::SVector<T,D> & a , 
            const ROOT::Math::SVector<T,D> & b ) 
      {
        double result = 0 ;
        for ( unsigned short i = 0 ; i < D ; ++i ) { result +=  double( a[i] ) * b[i] ; }
        return result ;
      }
      // cross: 
      static
      ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepStd<T,D,D> >
      cross ( const ROOT::Math::SVector<T,D> & a , 
              const ROOT::Math::SVector<T,D> & b ) 
      {
        ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepStd<T,D,D> > result ;
        for ( unsigned short i = 0 ; i < D ; ++i ) 
        { for ( unsigned short j = 0 ; j < D ; ++j ) 
          { result(i,j) = a[i] * b[j] ; } }
        return result ;
      }
      // multiply
      static 
      double 
      multiply ( const ROOT::Math::SVector<T,D> & a , 
                 const ROOT::Math::SVector<T,D> & b ) { return dot ( a , b ) ; }
    } ;
    // ========================================================================    
    /// cross/tensor : vector * vector  
    template <class T, unsigned int D1, unsigned int D2>
    struct MultiplyOp < ROOT::Math::SVector<T,D1>, ROOT::Math::SVector<T,D2> >
    {
      //
      static
      ROOT::Math::SMatrix<T,D1,D2,ROOT::Math::MatRepStd<T,D1,D2> >
      cross ( const ROOT::Math::SVector<T,D1> & a , 
              const ROOT::Math::SVector<T,D2> & b ) 
      {
        ROOT::Math::SMatrix<T,D1,D2,ROOT::Math::MatRepStd<T,D1,D2> > result ;
        for ( unsigned short i = 0 ; i < D1 ; ++i ) 
        { for ( unsigned short j = 0 ; j < D2 ; ++j ) 
          { result(i,j) = a[i] * b[j] ; } }
        return result ;
      }
    } ;
    // ========================================================================
    // vector * matrix 
    template <class T, unsigned int D, unsigned D2, class  R>
    struct MultiplyOp<ROOT::Math::SVector<T,D>,ROOT::Math::SMatrix<T,D,D2,R> >
    {
      static 
      ROOT::Math::SVector<T,D2> 
      multiply ( const ROOT::Math::SVector<T,D>      & a , 
                 const ROOT::Math::SMatrix<T,D,D2,R> & b ) { return a * b ; }
    } ;
    // =======================================================================
    // matrix * matrix  
    template <class T, unsigned int D1, unsigned D2, unsigned D3, class R1, class R2>
    struct MultiplyOp<ROOT::Math::SMatrix<T,D1,D2,R1>,ROOT::Math::SMatrix<T,D2,D3,R2> > 
    {
      static 
      ROOT::Math::SMatrix<T,D1,D3,ROOT::Math::MatRepStd<T,D1,D3> >  
      multiply  ( const ROOT::Math::SMatrix<T,D1,D2,R1> & a , 
                  const ROOT::Math::SMatrix<T,D2,D3,R2> & b ) { return a * b ; }
    } ;
    // =======================================================================
    // matrix * vector  
    template <class T, unsigned int D, unsigned D2, class  R>
    struct MultiplyOp <ROOT::Math::SMatrix<T,D,D2,R>,ROOT::Math::SVector<T,D2> >
    {
      static 
      ROOT::Math::SVector<T,D> 
      multiply  ( const ROOT::Math::SMatrix<T,D,D2,R> & a , 
                  const ROOT::Math::SVector<T,D2>     & b ) { return a * b ; }
    } ;
    // ========================================================================

    // ========================================================================
    template <class OBJ1, class OBJ2> 
    struct EqualityOp ;
    // ========================================================================
    // vector == vector 
    template <class T1,class T2, unsigned int D>
    struct EqualityOp< ROOT::Math::SVector<T1,D>,ROOT::Math::SVector<T2,D> >
    {
      static 
      bool 
      equal ( const ROOT::Math::SVector<T1,D>& v1 , 
              const ROOT::Math::SVector<T2,D>& v2 ) 
      {
        static const Equal_To<ROOT::Math::SVector<T1,D> > m_cmp{} ;
        return m_cmp ( v1 , v2 ) ;
      }
    } ;
    // ========================================================================
    // matrix == matrix 
    template <class T1,class T2, unsigned int D1,unsigned int D2, class R1, class R2>
    struct EqualityOp< ROOT::Math::SMatrix<T1,D1,D2,R1>,ROOT::Math::SMatrix<T2,D1,D2,R2> >
    {
      static 
      bool 
      equal ( const ROOT::Math::SMatrix<T1,D1,D2,R1>& v1 , 
              const ROOT::Math::SMatrix<T2,D1,D2,R2>& v2 ) 
      {
        static const Equal_To<ROOT::Math::SMatrix<T1,D1,D2,R1> > m_cmp{} ;
        return m_cmp ( v1 , v2 ) ;
      }
    } ;
    // ========================================================================
    
    // ========================================================================
    namespace  Ops
    {
      // ======================================================================
      // CHECKERS
      // ======================================================================
      
      template <class M1, class M2>
      struct CanAdd   { static bool operation ( const M1& /* m1 */ , const M2& /* m2 */ ) { return false ; } } ;
      
      template <class M1, class M2>
      struct CanMul   { static bool operation ( const M1& /* m1 */ , const M2& /* m2 */ ) { return false ; } } ;
      template <class M1, class M2>
      struct CanIMul  { static bool operation ( const M1& /* m1 */ , const M2& /* m2 */ ) { return false ; } } ;
      
      template <class M1, class M2>
      struct CanDiv   { static bool operation ( const M1& /* m1 */ , const M2& /* m2 */ ) { return false ; } } ;
      template <class M1, class M2>
      struct CanIDiv  { static bool operation ( const M1& /* m1 */ , const M2& /* m2 */ ) { return false ; } } ;
      
      template <class M1, class M2>
      struct CanDot   { static bool operation ( const M1& /* m1 */ , const M2& /* m2 */ ) { return false ; } } ;
      template <class M1, class M2>
      struct CanCross { static bool operation ( const M1& /* m1 */ , const M2& /* m2 */ ) { return false ; } } ;
      template <class M1, class M2>
      struct CanSim   { static bool operation ( const M1& /* m1 */ , const M2& /* m2 */ ) { return false ; } } ;
      template <class M1, class M2>
      struct CanSimT  { static bool operation ( const M1& /* m1 */ , const M2& /* m2 */ ) { return false ; } } ;

      template <class M1, class M2>
      struct CanEq    { static bool operation ( const M1& /* m1 */ , const M2& /* m2 */ ) { return false ; } } ;
      
      template <class M1>
      struct CanPow   { static bool operation ( const M1& /* m1 */ , const double /* p */ ) { return false ; } } ;

      template <class M1>
      struct CanSym   { static bool operation ( const M1& /* m1 */ ) { return false ; } } ;

      template <class M1>
      struct CanASym  { static bool operation ( const M1& /* m1 */ ) { return false ; } } ;

      template <class M1>
      struct CanInvert { static bool operation ( const M1& /* m1 */ ) { return false ; } } ;

      template <class M1, class M2>
      struct CanRMul
      { static bool operation ( const M1& m1 , const M2& m2 ) { return CanMul<M2,M1>::operation ( m2 , m1 ) ; } } ;
      
      // ======================================================================
      // partial specializations with scalar/double 
      // ======================================================================
      template <class M1>
      struct CanMul<M1,double>  
      { static bool operation ( const M1&    /* m1 */ , const double /* m2 */ ) { return true ; } } ;
      template <class M1>
      struct CanMul<double,M1> 
      { static bool operation ( const double /* m2 */ , const M1&    /* m1 */ ) { return true ; } } ;
      template <class M1>
      struct CanRMul<M1,double> 
      { static bool operation ( const M1&    /* m1 */ , const double /* m2 */ ) { return true ; } } ;
      template <class M1>
      struct CanIMul<M1,double> 
      { static bool operation ( const M1&    /* m1 */ , const double /* m2 */ ) { return true ; } } ;      
      template <class M1>
      struct CanDiv<M1,double> 
      { static bool operation ( const M1&    /* m1 */ , const double /* m2 */ ) { return true ; } } ;
      template <class M1>
      struct CanIDiv<M1,double> 
      { static bool operation ( const M1&    /* m1 */ , const double /* m2 */ ) { return true ; } } ;
      
      template <class T, unsigned int D, class R1>
      struct CanInvert< ROOT::Math::SMatrix<T,D,D,R1> > 
      {
        typedef ROOT::Math::SMatrix<T,D,D,R1> M1 ;
        typedef bool                          R  ;
        //
        static R operation ( const M1& /* m1 */ ) { return true ; } 
      } ;
            
      // ======================================================================
      // new cases  with "almost" scalar
      // ======================================================================      
      template <class M1, class T, class R1>
      struct CanMul<M1, ROOT::Math::SMatrix<T,1,1,R1> > 
      {
        typedef ROOT::Math::SMatrix<T,1,1,R1> M2 ;
        static bool operation ( const M1& /* m1 */ , const M2& /* m2 */ ) { return true ; } 
      } ;
      // ======================================================================
      template <class M1, class T, class R1>
      struct CanRMul<M1, ROOT::Math::SMatrix<T,1,1,R1> > 
      {
        typedef ROOT::Math::SMatrix<T,1,1,R1> M2 ;
        static bool operation ( const M1& /* m1 */ , const M2& /* m2 */ ) { return true ; } 
      } ;
      // ======================================================================
      template <class M1, class T, class R1>
      struct CanIMul<M1, ROOT::Math::SMatrix<T,1,1,R1> > 
      {
        typedef ROOT::Math::SMatrix<T,1,1,R1> M2 ;
        static bool operation ( const M1& /* m1 */ , const M2& /* m2 */ ) { return true ; } 
      } ;
      // ======================================================================
      template <class M1, class T, class R1>
      struct CanMul<ROOT::Math::SMatrix<T,1,1,R1> , M1 > 
      {
        typedef ROOT::Math::SMatrix<T,1,1,R1> M2 ;
        static bool operation ( const M2& /* m1 */ , 
                                const M1& /* m2 */ ) { return true ; } 
      } ;      
      // ======================================================================
      template <class M1, class T, class R1>
      struct CanRMul<ROOT::Math::SMatrix<T,1,1,R1> , M1 > 
      {
        typedef ROOT::Math::SMatrix<T,1,1,R1> M2 ;
        static bool operation ( const M2& /* m1 */ , const M1& /* m2 */ ) { return true ; } 
      } ;
      // ======================================================================

      
      // ======================================================================
      // Can be added 
      // ======================================================================
      
      template <class T,unsigned int D1,unsigned int D2,class R1,
                class R2>
      struct CanAdd<ROOT::Math::SMatrix<T,D1,D2,R1> ,
                    ROOT::Math::SMatrix<T,D1,D2,R2> >
      {
        static bool operation
        ( const ROOT::Math::SMatrix<T,D1,D2,R1>& /* m1 */ , 
          const ROOT::Math::SMatrix<T,D1,D2,R2>& /* m2 */ ) { return true ; }  
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct CanAdd<ROOT::Math::SVector<T,D> ,
                    ROOT::Math::SVector<T,D> >
      {
        static bool operation
        ( const ROOT::Math::SVector<T,D>& /* m1 */ , 
          const ROOT::Math::SVector<T,D>& /* m2 */ ) { return true ; }  
      } ;
      // ======================================================================
      template <class T,unsigned int D,class R1>
      struct CanAdd<ROOT::Math::SMatrix<T,D,D,R1> , double>
      {
        static bool operation
        ( const ROOT::Math::SMatrix<T,D,D,R1>& /* m1 */ , 
          const double                         /* m2 */ ) { return true ; }  
      } ;
      
      // ======================================================================
      // Can be multiplied
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2,class R1,
                unsigned int D3,class R2>
      struct CanMul<ROOT::Math::SMatrix<T,D1,D2,R1> ,
                    ROOT::Math::SMatrix<T,D2,D3,R2> >
      {
        static bool operation
        ( const ROOT::Math::SMatrix<T,D1,D2,R1>& /* m1 */ , 
          const ROOT::Math::SMatrix<T,D2,D3,R2>& /* m2 */ ) { return true ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2,class R1>
      struct CanMul<ROOT::Math::SMatrix<T,D1,D2,R1> ,
                    ROOT::Math::SVector<T,D2> >
      {
        static bool operation
        ( const ROOT::Math::SMatrix<T,D1,D2,R1>& /* m1 */ , 
          const ROOT::Math::SVector<T,D2>&       /* m2 */ ) { return true ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2,class R1>
      struct CanMul<ROOT::Math::SVector<T,D1>       , 
                    ROOT::Math::SMatrix<T,D1,D2,R1> >
      {
        static bool operation
        ( const ROOT::Math::SVector<T,D1>&       /* m2 */ ,
          const ROOT::Math::SMatrix<T,D1,D2,R1>& /* m1 */ )  { return true ; } 
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct CanMul<ROOT::Math::SVector<T,D>       , 
                    ROOT::Math::SVector<T,D> >
      {
        static bool operation
        ( const ROOT::Math::SVector<T,D>& /* m2 */ ,
          const ROOT::Math::SVector<T,D>& /* m1 */ )  { return true ; } 
      } ;
      // ========================================================================
      
      // ========================================================================
      template <class T,unsigned int D1, unsigned int D2,
                class R2>
      struct CanIMul<ROOT::Math::SMatrix<T,D1,D2>    ,
                     ROOT::Math::SMatrix<T,D2,D2,R2> >
      {
        static bool operation
        ( const ROOT::Math::SMatrix<T,D1,D2>&    /* m1 */ , 
          const ROOT::Math::SMatrix<T,D2,D2,R2>& /* m2 */ ) { return true ; }
      } ;
      // ======================================================================
      
      // ======================================================================
      template <class T,unsigned int D>
      struct CanDot<ROOT::Math::SVector<T,D> , 
                    ROOT::Math::SVector<T,D> >
      {
        static bool operation 
        ( const ROOT::Math::SVector<T,D>& /* m2 */ ,
          const ROOT::Math::SVector<T,D>& /* m1 */ )  { return true ; } 
      } ;
      // ======================================================================
      
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2>
      struct CanCross<ROOT::Math::SVector<T,D1> , 
                      ROOT::Math::SVector<T,D2> >
      {
        static bool operation
        ( const ROOT::Math::SVector<T,D1>& /* m2 */ ,
          const ROOT::Math::SVector<T,D2>& /* m1 */ )  { return true ; } 
      } ;
      // ======================================================================
      
      // ======================================================================
      template <class T,unsigned int D,unsigned int D2, class R2>
      struct CanSim<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > ,
                    ROOT::Math::SMatrix<T,D2,D,R2> >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef ROOT::Math::SMatrix<T,D2,D,R2>                         M2 ;
        // check
        static bool operation
        ( const M1&  /* m1 */ ,
          const M2&  /* m2 */ ) { return true ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct CanSim<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D>> ,
                    ROOT::Math::SVector<T,D> >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef ROOT::Math::SVector<T,D>                               M2 ;
        // check 
        static bool operation 
          ( const M1& /* m1 */ ,
            const M2& /* m2 */ ) { return true ; }
      } ;
      // =====================================================================
      template <class T,unsigned int D> 
      struct CanSim<ROOT::Math::SVector<T,D> ,
                    ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > >
      {
        typedef ROOT::Math::SVector<T,D>                               M1 ;
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M2 ;
        // check 
        static bool operation 
        ( const M1& /* m1 */ ,
          const M2& /* m2 */ ) { return true ; }
      } ;
      // ======================================================================
      
      // ======================================================================
      template <class T,unsigned int D,unsigned int D2, class R2> 
      struct CanSimT<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > ,
                     ROOT::Math::SMatrix<T,D,D2,R2> >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >    M1 ;
        typedef ROOT::Math::SMatrix<T,D,D2,R2>                            M2 ;
        // check 
        static bool operation ( const M1& /* m1 */ , const M2& /* m2 */ ) { return true ; }
      } ;
      // ======================================================================
      
      template <class T, unsigned int D, class R1>
      struct CanPow<ROOT::Math::SMatrix<T,D,D,R1> >
      { 
        static bool operation 
        ( const ROOT::Math::SMatrix<T,D,D,R1>& /* m1 */ ,
          const unsigned short                 /* p  */ ) { return true ; } 
      } ; 

      template <class T, unsigned int D, class R1>
      struct CanSym<ROOT::Math::SMatrix<T,D,D,R1> >
      { static bool operation ( const ROOT::Math::SMatrix<T,D,D,R1>& /* m1 */ ) { return true ; } } ;
      
      template <class T, unsigned int D, class R1>
      struct CanASym<ROOT::Math::SMatrix<T,D,D,R1> >
      { static bool operation ( const ROOT::Math::SMatrix<T,D,D,R1>& /* m1 */ ) { return true ; } } ;

      // ======================================================================
      // Operations
      // ======================================================================
      template <class T1,class T2>
      struct  Add  ;
      template <class T1,class T2>
      struct IAdd ;
      template <class T1,class T2>
      struct  Sub  ;
      template <class T1,class T2>
      struct ISub ;
      template <class T1,class T2>
      struct  Mul  ;
      template <class T1,class T2>
      struct IMul ;
      template <class T1,class T2>
      struct  Div  ;
      template <class T1,class T2>
      struct IDiv ;
      
      template <class T1,class T2>
      struct RAdd ;
      template <class T1,class T2>
      struct RSub ;
      template <class T1,class T2>
      struct RDiv ;
      
      template <class T1,class T2>
      struct Dot   ;
      template <class T1,class T2>
      struct Cross ;
      template <class T1,class T2>
      struct Sim   ;
      template <class T1,class T2>
      struct SimT  ;


      template <class T1,class T2>
      struct Eq ;

      template <class T1>
      struct Pow ;

      template <class T1>
      struct Sym  ;
      
      template <class T1>
      struct ASym ;

      template <class T1>
      struct Invert ;

      // ======================================================================
      // "Right" operations
      // ======================================================================
      template <class M1, class M2>
      struct RAdd
      {
        typedef Add<M2,M1>    A ;
        typedef typename A::R R ;
        static R operation ( const M1& m1 , const  M2& m2 ) { return A::operation ( m2 , m1 ) ; }
      } ;
      // ======================================================================
      template <class M1, class M2>
      struct RSub
      {
        typedef Sub<M2,M1>    S ;
        typedef typename S::R R ;
        static R operation ( const M1& m1 , const  M2& m2 ) { return S::operation ( m2 , m1 ) ; }
      } ;
      // ======================================================================
      template <class M1, class M2>
      struct RMul
      {
        typedef Mul<M2,M1>    M ;
        typedef typename M::R R ;
        static R operation ( const M1& m1 , const  M2& m2 ) { return M::operation ( m2 , m1 ) ; }
      } ;
      
      // ======================================================================
      // scaling
      // ======================================================================
      template <class M1>
      struct IMul<M1,double>
      { static void operation ( M1& m1 , const double m2 ) 
        {  m1 *= m2 ; } } ;
      // ======================================================================
      template <class M1>
      struct IDiv<M1,double>
      { static void operation ( M1& m1 , const double m2 )
        { IMul<M1,double>::operation ( m1 , 1 / m2 ) ; } } ;
      // ======================================================================
      template <class M1>
      struct Mul<M1,double>
      {
        typedef M1 R ;
        static R operation ( const M1& m1 , const double m2 ) { return m1 * m2  ; }
      } ;
      // ======================================================================
      template <class M1>
      struct RMul<M1,double>
      {
        typedef M1 R ;
        static R operation ( const M1& m1 , const double m2 ) 
        { return Mul<M1,double>::operation ( m1 , m2 ) ; }
      } ;
      // ======================================================================
      template <class M1>
      struct Div<M1,double>
      {
        typedef M1 R ;
        static R operation ( const M1& m1 , const double m2 )
        { return Mul<M1,double>::operation ( m1 , 1 / m2 ) ; }
      } ;
      // ======================================================================

      // ======================================================================
      // Trivia 
      // ======================================================================
      template <class M>
      struct IAdd<M,M>     
      { static void operation ( M& m1 , const M& m2 ) { m1 += m2 ; } } ;
      // ======================================================================
      template <class M>
      struct ISub<M,M>
      { static void operation ( M& m1 , const M& m2 ) { m1 -= m2 ; } } ;
      

      // ======================================================================
      /// "almost constants"
      // ======================================================================
      template <class M1, class T, class R1>
      struct Mul<M1, ROOT::Math::SMatrix<T,1,1,R1> > 
      {
        typedef ROOT::Math::SMatrix<T,1,1,R1> M2 ;
        typedef Mul<M1,double>                O  ;
        typedef typename O::R                 R  ;
        //
        static R operation ( const M1& m1 , const M2& m2 )   
        { return O::operation (  m1 , m2 ( 0 , 0 ) ) ; }
      } ;
      // ======================================================================
      template <class M1, class T, class R1>
      struct RMul<M1, ROOT::Math::SMatrix<T,1,1,R1> > 
      {
        typedef ROOT::Math::SMatrix<T,1,1,R1> M2 ;
        typedef Mul<M1,double>                O  ;
        typedef typename O::R                 R  ;
        //
        static R operation ( const M1& m1 , const M2& m2 )
        { return O::operation ( m1 , m2 ( 0 , 0 ) ) ; }
      } ;      
      // ======================================================================
      template <class M1, class T, class R1>
      struct IMul<M1, ROOT::Math::SMatrix<T,1,1,R1> > 
      {
        typedef ROOT::Math::SMatrix<T,1,1,R1> M2 ;
        typedef IMul<M1,double>               O  ;
        static void operation ( M1& m1 , const M2& m2 ) { O::operation ( m1 , m2 ( 0 , 0 ) ) ; }
      } ;
      // ======================================================================
      template <class M1, class T, class R1>
      struct Mul<ROOT::Math::SMatrix<T,1,1,R1>, M1 > 
      {
        typedef ROOT::Math::SMatrix<T,1,1,R1> M2 ;
        typedef Mul<M1,double>                O  ;
        typedef typename O::R                 R  ;
        //
        static R operation ( const M2& m2 , const M1& m1 ) 
        { return O::operation ( m1 , m2 ( 0 , 0 ) ) ; } 
      } ;
      // ======================================================================
      template <class M1, class T, class R1>
      struct RMul<ROOT::Math::SMatrix<T,1,1,R1>, M1> 
      {
        typedef ROOT::Math::SMatrix<T,1,1,R1> M2 ;
        typedef RMul<M1,double>               O  ;
        typedef typename O::R                 R  ;
        //
        static R operation ( const M2& m2 , const M1& m1 ) 
        { return O::operation ( m1 , m2 ( 0 , 0 ) ) ; } 
      } ;
      // ======================================================================

      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2>
      struct IAdd<ROOT::Math::SMatrix<T,D1,D2> ,
                  ROOT::Math::SMatrix<T,D1,D2> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2> M ;
        static void operation ( M& m1 , const M& m2 ) { m1 += m2 ; }
      }; 
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2>
      struct ISub<ROOT::Math::SMatrix<T,D1,D2> ,
                  ROOT::Math::SMatrix<T,D1,D2> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2> M ;
        static void operation ( M& m1 , const M& m2 ) { m1 -= m2 ; }
      }; 
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2,class R>
      struct IAdd<ROOT::Math::SMatrix<T,D1,D2,R> ,
                  ROOT::Math::SMatrix<T,D1,D2,R> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R> M ;
        static void operation ( M& m1 , const M& m2 ) { m1 += m2 ; }
      }; 
      // =====================================================================
      template <class T,unsigned int D1,unsigned int D2,class R>
      struct ISub<ROOT::Math::SMatrix<T,D1,D2,R> ,
                  ROOT::Math::SMatrix<T,D1,D2,R> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R> M ;
        static void operation ( M& m1 , const M& m2 ) { m1 -= m2 ; }
      };
      // ======================================================================
      template <class T,unsigned int D>
      struct IAdd<ROOT::Math::SVector<T,D> ,
                  ROOT::Math::SVector<T,D> >
      {
        typedef ROOT::Math::SVector<T,D> M ;
        static void operation ( M& m1 , const M& m2 ) { m1 += m2 ; }
      }; 
      // ======================================================================
      template <class T,unsigned int D>
      struct ISub<ROOT::Math::SVector<T,D> ,
                  ROOT::Math::SVector<T,D> >
      {
        typedef ROOT::Math::SVector<T,D> M ;
        static void operation ( M& m1 , const M& m2 ) { m1 -= m2 ; }
      };

      // ======================================================================
      // ADD
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2,class R1,
                class R2>
      struct Add<ROOT::Math::SMatrix<T,D1,D2,R1> ,
                 ROOT::Math::SMatrix<T,D1,D2,R2> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M1 ;
        typedef ROOT::Math::SMatrix<T,D1,D2,R2> M2 ;
        typedef ROOT::Math::SMatrix<T,D1,D2>    R  ;
        // addition
        static R operation ( const M1& m1 , const M2& m2 ) { return m1 + m2 ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct Add<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > ,
                 ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M2 ;
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > R  ;
        // addition
        static R operation ( const M1 & m1 , const M2 & m2 ) { return m1 + m2 ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct Add<ROOT::Math::SVector<T,D> ,
                 ROOT::Math::SVector<T,D> >
      {
        typedef ROOT::Math::SVector<T,D> M1 ;
        typedef ROOT::Math::SVector<T,D> M2 ;
        typedef ROOT::Math::SVector<T,D> R  ;
        // addition
        static R operation ( const M1 & m1 , const M2 & m2 ) { return m1 + m2 ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D,class R1>
      struct Add<ROOT::Math::SMatrix<T,D,D,R1>, double> 
      {
        typedef ROOT::Math::SMatrix<T,D,D,R1> M1 ;
        typedef double                        M2 ;
        typedef ROOT::Math::SMatrix<T,D,D,R1> R  ;
        // addition
        static R operation ( const M1& m1 , const double m2 ) 
        { 
          R result  { m1 } ;
          for  ( unsigned int i = 0 ; i < D ; ++i ) { result ( i , i )  += m2 ; }
          return result ; 
        }
      } ;
      // ======================================================================
      template <class T,unsigned int D,class R1>
      struct RAdd<ROOT::Math::SMatrix<T,D,D,R1>, double> 
        : public Add<ROOT::Math::SMatrix<T,D,D,R1>, double> {} ;
      

      // ======================================================================
      // IADD
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2, class R2>
      struct IAdd<ROOT::Math::SMatrix<T,D1,D2>    ,
                  ROOT::Math::SMatrix<T,D1,D2,R2> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2>    M1 ;
        typedef ROOT::Math::SMatrix<T,D1,D2,R2> M2 ;
        // addition
        static void operation ( M1 & m1 , const M2 & m2 ) { m1 += m2 ; }
      } ;
      // ======================================================================      
      template <class T,unsigned int D>
      struct IAdd<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > ,
                 ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > >
      {
       typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
       typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M2 ;
       // addition
        static void operation ( M1 & m1 , const M2 & m2 ) { m1 += m2 ; }
      } ;
      // ======================================================================      
      template <class T,unsigned int D,class R1>
      struct IAdd<ROOT::Math::SMatrix<T,D,D,R1> , double> 
      {
        typedef ROOT::Math::SMatrix<T,D,D,R1> M1 ;
        typedef double                        M2 ;
        // addition
        static void operation ( M1& m1 , const double m2 ) 
        { for  ( unsigned int i = 0 ; i < D ; ++i ) { m1 ( i , i )  += m2 ; } }
      } ;

      
      // ======================================================================
      // SUB
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2,class R1,
                class R2>
      struct Sub<ROOT::Math::SMatrix<T,D1,D2,R1> ,
                 ROOT::Math::SMatrix<T,D1,D2,R2> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M1 ;
        typedef ROOT::Math::SMatrix<T,D1,D2,R2> M2 ;
        typedef ROOT::Math::SMatrix<T,D1,D2>    R  ;
        // subtraction 
        static R operation ( const M1& m1 , const M2& m2 ) { return m1 - m2 ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct Sub<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > ,
                 ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M2 ;
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > R  ;
        // subtraction 
        static R operation ( const M1 & m1 , const M2 & m2 ) { return m1 - m2 ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct Sub<ROOT::Math::SVector<T,D> ,
                 ROOT::Math::SVector<T,D> >
      {
        typedef ROOT::Math::SVector<T,D> M1 ;
        typedef ROOT::Math::SVector<T,D> M2 ;
        typedef ROOT::Math::SVector<T,D> R  ;
        // subtraction
        static R operation ( const M1 & m1 , const M2 & m2 ) { return m1 - m2 ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D,class R1>
      struct Sub<ROOT::Math::SMatrix<T,D,D,R1> , double> 
      {
        typedef ROOT::Math::SMatrix<T,D,D,R1> M1 ;
        typedef double                        M2 ;
        typedef ROOT::Math::SMatrix<T,D,D,R1> R  ;
        // addition
        static R operation ( const M1& m1 , const double m2 ) 
        { 
          R result  { m1 } ;
          for  ( unsigned int i = 0 ; i < D ; ++i ) { result ( i , i ) -= m2 ; }
          return result ; 
        }
      } ;
      // ======================================================================
      template <class T,unsigned int D,class R1>
      struct RSub<ROOT::Math::SMatrix<T,D,D,R1> , double> 
      {
        typedef ROOT::Math::SMatrix<T,D,D,R1> M1 ;
        typedef double                        M2 ;
        typedef ROOT::Math::SMatrix<T,D,D,R1> R  ;
        // addition
        static R operation ( const M1& m1 , const double m2 ) 
        { 
          R result  { m1 } ; result *= -1 ; // ATTENTION! 
          for  ( unsigned int i = 0 ; i < D ; ++i ) { result ( i , i ) += m2 ; }
          return result ; 
        }
      } ;
      
      // ======================================================================
      // ISUB
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2, class R2>
      struct ISub<ROOT::Math::SMatrix<T,D1,D2>    ,
                  ROOT::Math::SMatrix<T,D1,D2,R2> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2>    M1 ;
        typedef ROOT::Math::SMatrix<T,D1,D2,R2> M2 ;
        // subtraction
        static void operation ( M1 & m1 , const M2 & m2 ) { m1 -= m2 ; }
      } ;
      // ======================================================================

      template <class T,unsigned int D>
      struct ISub<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > ,
                  ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M2 ;
        // addition
        static void operation ( M1 & m1 , const M2 & m2 ) { m1 -= m2 ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D,class R1>
      struct ISub<ROOT::Math::SMatrix<T,D,D,R1> , double> 
      {
        typedef ROOT::Math::SMatrix<T,D,D,R1> M1 ;
        typedef double                        M2 ;
        // addition
        static void operation ( M1& m1 , const double m2 ) 
        { for  ( unsigned int i = 0 ; i < D ; ++i ) { m1 ( i , i ) -= m2 ; } }
      } ;

      
      // ======================================================================
      // MUL 
      // ======================================================================      
      template <class T,unsigned int D1,unsigned int D2,class R1,
                unsigned int D3,class R2>
      struct Mul<ROOT::Math::SMatrix<T,D1,D2,R1> ,
                 ROOT::Math::SMatrix<T,D2,D3,R2> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M1 ;
        typedef ROOT::Math::SMatrix<T,D2,D3,R2> M2 ;
        typedef ROOT::Math::SMatrix<T,D1,D3>     R  ;
        // multiplication
        static R operation ( const M1 & m1 , const M2 & m2 ) { return m1 * m2 ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2,class R1>
      struct Mul<ROOT::Math::SMatrix<T,D1,D2,R1> ,
                 ROOT::Math::SVector<T,D2>       >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M1 ;
        typedef ROOT::Math::SVector<T,D2>       M2 ;
        typedef ROOT::Math::SVector<T,D1>       R  ;
        // multiplication
        static R operation ( const M1 & m1 , const M2 & m2 ) { return m1 * m2 ; }
      } ;
      // ===================================================================
      template <class T,unsigned int D1,unsigned int D2,class R1>
      struct Mul<ROOT::Math::SVector<T,D1>       ,
                 ROOT::Math::SMatrix<T,D1,D2,R1> >
      {
        typedef ROOT::Math::SVector<T,D1>       M1 ;
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M2 ;
        typedef ROOT::Math::SVector<T,D2>       R  ;
        // multiplication
        static R operation ( const M1& m1 , const M2 & m2 ) { return m1 * m2 ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2,class R1>
      struct Mul<ROOT::Math::SMatrix<T,D1,D2,R1>,double>
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M1 ;
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> R  ;
        // multiplication
        static R operation ( const M1& m1 , const double m2 ) { return m1 * m2  ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2,class R1>
      struct Mul<double,ROOT::Math::SMatrix<T,D1,D2,R1>>
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M1 ;
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> R  ;
        // multiplication
        static R operation ( const double m2 , const M1& m1) { return m1 * m2  ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct Mul<ROOT::Math::SVector<T,D> ,
                 ROOT::Math::SVector<T,D> >
      {
        typedef ROOT::Math::SVector<T,D>       M1 ;
        typedef ROOT::Math::SVector<T,D>       M2 ;
        typedef double                         R  ;
        // multiplication
        static double operation ( const M1 & m1 , const M2 & m2 )
        { return std::inner_product ( m1.begin() , m1.end() , m2.begin() , 0.0 ); }
      } ;
      // ===============================================================
      template <class T,unsigned int D>
      struct Mul<ROOT::Math::SVector<T,D>,double>
      {
        typedef ROOT::Math::SVector<T,D> M1 ;
        typedef ROOT::Math::SVector<T,D> R  ;
        // multiplication
        static R operation ( const M1& m1 , const double m2 ) { return m1 * m2  ; }
      } ;
      // ===============================================================      
      template <class T,unsigned int D>
      struct Mul<double,ROOT::Math::SVector<T,D>>
      {
        typedef ROOT::Math::SVector<T,D> M1 ;
        typedef ROOT::Math::SVector<T,D> R  ;
        // multiplication
        static R operation ( const double m2,  const M1& m1 ) { return m1 * m2  ; }
      } ;
      // ======================================================================

      
      // ======================================================================
      // IMUL 
      // ======================================================================      
      template <class T,unsigned int D1,unsigned int D2,
                class R2>
      struct IMul<ROOT::Math::SMatrix<T,D1,D2>    ,
                  ROOT::Math::SMatrix<T,D2,D2,R2> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2>    M1 ;
        typedef ROOT::Math::SMatrix<T,D2,D2,R2> M2 ;
        // in-place multiplication
        static void operation ( M1& m1 , const M2&m2 ) { m1 *= m2 ; }
      } ;
      // ======================================================================


      // ======================================================================
      // DOT
      // ======================================================================
      template <class T1,unsigned int D, class T2>
      struct Dot<ROOT::Math::SVector<T1,D> ,
                 ROOT::Math::SVector<T2,D> >
      {
        typedef ROOT::Math::SVector<T1,D> M1 ;
        typedef ROOT::Math::SVector<T2,D> M2 ;
        typedef double                    R  ;
        // multiplication
        static double operation ( const M1 & m1 , const M2 & m2 )
        { return std::inner_product ( m1.begin() , m1.end() , m2.begin() , 0.0 ) ; }
      } ;            
      // ======================================================================

      // ======================================================================
      // CROSS
      // ======================================================================
      template <class T,unsigned int D1,unsigned D2>
      struct Cross<ROOT::Math::SVector<T,D1> ,
                   ROOT::Math::SVector<T,D2> >
      {
        typedef ROOT::Math::SVector<T,D1>       M1 ;
        typedef ROOT::Math::SVector<T,D2>       M2 ;
        typedef ROOT::Math::SMatrix<T,D1,D2>    R  ;
        // multiplication
        static R operation ( const M1 & m1 , const M2 & m2 )
        {
          R r ;
          for ( unsigned int i = 0 ; i < D1 ; ++i ) 
          { for ( unsigned int j = 0 ; j < D2 ; ++j ) 
            { r(i,j) = m1[i] * m2[j] ; } }
          return r;
        }
      } ;      
      // ======================================================================

      // ======================================================================
      // SIM 
      // ======================================================================
      template <class T,unsigned int D,
                unsigned int D2, class R2> 
      struct Sim<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > ,
                 ROOT::Math::SMatrix<T,D2,D,R2> >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >    M1 ;
        typedef ROOT::Math::SMatrix<T,D2,D,R2>                            M2 ;
        typedef ROOT::Math::SMatrix<T,D2,D2,ROOT::Math::MatRepSym<T,D2> > R  ;        
        // similarity
        static R operation ( const M1& A , const M2& U  )
        { return ROOT::Math::Similarity ( U , A ) ; }
      } ;
      // ======================================================================      
      template <class T,unsigned int D> 
      struct Sim<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > ,
                 ROOT::Math::SVector<T,D> >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef ROOT::Math::SVector<T,D>                               M2 ;
        typedef double                                                 R  ;
        // similarity 
        static double operation ( const M1& A , const M2& V )
        { return ROOT::Math::Similarity ( A , V ) ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D> 
      struct Sim<ROOT::Math::SVector<T,D> ,
                 ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > >
      {
        typedef ROOT::Math::SVector<T,D>                               M1 ;
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M2 ;
        typedef double                                                 R  ;
        // check 
        static double operation ( const M1& V , const M2& A  )
        { return ROOT::Math::Similarity ( A , V ) ; }
      } ;
      // ======================================================================
      
      // ======================================================================
      // SIMT 
      // ======================================================================
      template <class T,unsigned int D,
                unsigned int D2, class R2> 
      struct SimT<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > ,
                  ROOT::Math::SMatrix<T,D,D2,R2> >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >    M1 ;
        typedef ROOT::Math::SMatrix<T,D,D2,R2>                            M2 ;
        typedef ROOT::Math::SMatrix<T,D2,D2,ROOT::Math::MatRepSym<T,D2> > R  ;        
        // similarity
        static R operation ( const M1& A , const M2& U  )
        { return ROOT::Math::SimilarityT ( U , A ) ; }
      } ;
      // ======================================================================

      
      // ======================================================================
      template <class T, unsigned int D, class R1>
      struct Invert<ROOT::Math::SMatrix<T,D,D,R1> > 
      {
        typedef ROOT::Math::SMatrix<T,D,D,R1> M1 ;
        typedef ROOT::Math::SMatrix<T,D,D,R1> R  ;
        //
        static R operation ( const M1& m1 , int& flag  )
        { return m1.Inverse ( flag ) ; }
      } ;      
      // ======================================================================
      
      // ======================================================================
      /// can be compared ?
      // ======================================================================
      template <class T1, unsigned D, class T2> 
      struct CanEq < ROOT::Math::SVector<T1,D> , 
                     ROOT::Math::SVector<T2,D> >
      {
        typedef ROOT::Math::SVector<T1,D> M1 ;
        typedef ROOT::Math::SVector<T2,D> M2 ;
        // equality 
        static bool operation 
        ( const M1& /* m1 */ , 
          const M2& /* m2 */ ) { return true ; }
      } ;
      // ======================================================================
      template <class T, unsigned D> 
      struct CanEq < ROOT::Math::SVector<T,D> , 
                     ROOT::Math::SVector<T,D> >
      {
        typedef ROOT::Math::SVector<T,D> M1 ;
        typedef ROOT::Math::SVector<T,D> M2 ;
        // equality 
        static bool operation 
        ( const M1& /* m1 */ , 
          const M2& /* m2 */ ) { return true ; }
      } ;
      // ======================================================================
      template <class T, unsigned D1, unsigned D2, class R1> 
      struct CanEq < ROOT::Math::SMatrix<T,D1,D2,R1> , 
                     ROOT::Math::SMatrix<T,D1,D2,R1> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M1 ;
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M2 ;
        // equality 
        static bool operation 
        ( const M1& /* m1 */ , 
          const M2& /* m2 */ ) { return true ; }
      } ;
      // ======================================================================
      template <class T, unsigned D1, unsigned D2, class R1, class R2> 
      struct CanEq < ROOT::Math::SMatrix<T,D1,D2,R1> , 
                     ROOT::Math::SMatrix<T,D1,D2,R2> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M1 ;
        typedef ROOT::Math::SMatrix<T,D1,D2,R2> M2 ;
        // equality 
        static bool operation 
        ( const M1& /* m1 */ , 
          const M2& /* m2 */ ) { return true ; }
      } ;
      // ======================================================================
      template <class T1, unsigned D1, unsigned D2, class R1, class T2, class R2> 
      struct CanEq < ROOT::Math::SMatrix<T1,D1,D2,R1> , 
                     ROOT::Math::SMatrix<T2,D1,D2,R2> >
      {
        typedef ROOT::Math::SMatrix<T1,D1,D2,R1> M1 ;
        typedef ROOT::Math::SMatrix<T2,D1,D2,R2> M2 ;
        // equality 
        static bool operation 
        ( const M1& /* m1 */ , 
          const M2& /* m2 */ ) { return true ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D,class R1>
      struct CanEq<ROOT::Math::SMatrix<T,D,D,R1> , double >
      {
        typedef ROOT::Math::SMatrix<T,D,D,R1> M1 ;
        typedef double                        M2 ;
        typedef bool                          R  ;
        // 
        static R operation
        ( const M1& /* m1 */ ,
          const M2  /* m2 */ ) { return true ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct CanEq<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > , double >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >  M1 ;
        typedef double                                                  M2 ;
        typedef bool                                                    R  ;
        // 
        static R operation
        ( const M1& /* m1 */ ,
          const M2  /* m2 */ ) { return true ; }
      } ;

      // ======================================================================
      // Equality
      // ======================================================================
      template <class T1,unsigned int D1,unsigned int D2,class R1,
                class T2, class R2>
      struct Eq<ROOT::Math::SMatrix<T1,D1,D2,R1> ,
                ROOT::Math::SMatrix<T2,D1,D2,R2> >
      {
        typedef ROOT::Math::SMatrix<T1,D1,D2,R1> M1 ;
        typedef ROOT::Math::SMatrix<T2,D1,D2,R2> M2 ;
        typedef bool                      R  ;
        // addition
        static R operation ( const M1& m1 , const M2& m2 )
        {
          static const Ostap::Math::Equal_To<M1> s_cmp ;
          return s_cmp ( m1 ,  m2 ) ;
        }
      } ;
      // ======================================================================
      template <class T1,unsigned int D,
                class T2>
      struct Eq<ROOT::Math::SVector<T1,D> ,
                ROOT::Math::SVector<T2,D> >
      {
        typedef ROOT::Math::SVector<T1,D> M1 ;
        typedef ROOT::Math::SVector<T2,D> M2 ;
        typedef bool                      R  ;
        // addition
        static R operation ( const M1& m1 , const M2& m2 )
        {
          static const Ostap::Math::Equal_To<M1> s_cmp ;
          return s_cmp ( m1 ,  m2 ) ;
        }
      } ;
      // ======================================================================

      
      template <class T,unsigned int D,class R1>
      struct Eq<ROOT::Math::SMatrix<T,D,D,R1> , double >
      {
        typedef ROOT::Math::SMatrix<T,D,D,R1> M1 ;
        typedef double                        M2 ;
        typedef bool                          R  ;
        // 
        static R operation ( const M1& m1 , const M2 m2 )
        {
          static const Ostap::Math::Equal_To<T>  s_cmp  ;
          static const Ostap::Math::Zero<T>      s_zero ;
          for ( unsigned int i = 0 ; i < D ; ++i ) 
          {
            if ( !s_cmp  ( m1 ( i , i ) , m2 )        ) { return false ; }            
            for ( unsigned int j = 0 ; j < D ; ++j ) 
            { 
              if ( i != j && !s_zero ( m1 ( i , j ) ) ) { return false ; }
            }
          }
          return true ;
        }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct Eq<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > , double >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef double                                                 M2 ;
        typedef bool                                                   R  ;
        //
        static R operation ( const M1& m1 , const M2 m2 )
        {
          static const Ostap::Math::Equal_To<T>  s_cmp  ;
          static const Ostap::Math::Zero<T>      s_zero ;
          for ( unsigned int i = 0 ; i < D ; ++i ) 
          {
            if ( !s_cmp  ( m1 ( i , i ) , m2 )   ) { return false ; }            
            for ( unsigned int j = i +1  ; j < D ; ++j ) 
            {
              if ( !s_zero ( m1 ( i , j )      ) ) { return false ; }
            }
          }
          return true ;
        }
      } ;

      // =======================================================================
      // EXTRA
      // =======================================================================
      template <class T>
      struct Eigen ;
          
      template <class T,unsigned int D> 
      struct Eigen<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef ROOT::Math::SVector<T,D>                               M2 ;
        typedef ROOT::Math::SMatrix<T,D,D>                             M3 ;
        // get eigen values 
        static Ostap::StatusCode operation 
        ( const M1&  m             , 
          M2&        values        ,
          const bool sorted = true )
        {
          Ostap::Math::GSL::EigenSystem eigen {} ;
          return eigen.eigenValues  ( m , values , sorted ) ;
        }
        // get eigen values and eigenvectors 
        static Ostap::StatusCode operation 
        ( const M1&  m                 , 
          M2&        values            , 
          M3&        vectors           ,
          const bool sorted     = true ,
          const bool ascending  = true )
        {
          Ostap::Math::GSL::EigenSystem eigen {} ;
          return eigen.eigenVectors ( m , values , vectors , sorted , ascending ) ;
        }
      } ;
      
      // ======================================================================
      template <class T, unsigned int D, class R1>
      struct Pow<ROOT::Math::SMatrix<T,D,D,R1> >
      {
        //
        typedef ROOT::Math::SMatrix<T,D,D,R1> M ;
        typedef ROOT::Math::SMatrix<T,D,D>    R ;
        //
        static R operation ( const M& m , const unsigned short n )
        {
          //
          if      ( 0 == n ) { return M ( ROOT::Math::SMatrixIdentity () ) ; }
          else if ( 1 == n ) { return         m ; }
          else if ( 2 == n ) { return     m * m ; }
          else if ( 3 == n ) { return m * m * m ; }
          //
          R r = operation ( m , n / 2 )  ;
          if ( 0 == n / 2 ) { return r * r ; }
          //
          return r * r * m ;
        }
      } ;
      // ======================================================================
      template <class T, class R1>
      struct Pow<ROOT::Math::SMatrix<T,1,1,R1> >
      {
        //
        typedef ROOT::Math::SMatrix<T,1,1,R1> M ;
        typedef double   R ;
        //
        static R operation ( const M& m , const int    n )
        { return 0 == n ? 1.0 : std::pow ( m ( 0 , 0 )  , n ) ; }
      } ;
      // ======================================================================
      template <class T, unsigned int D>
      struct Sym<ROOT::Math::SMatrix<T,D,D> >
      {
        //
        typedef ROOT::Math::SMatrix<T,D,D>                             M ;
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > R ;
        //
        static R operation ( const M& m )
        {
          //
          R r ;
          for ( unsigned int i = 0 ; i < D ; ++i )
          { r ( i , i ) = m ( i , i ) ;
            for ( unsigned int j = i + 1  ; j < D ; ++j )
            { r ( i , j ) = 0.5 * ( m ( i , j ) + m ( j , i ) ) ; } }
          //
          return r ;
        }
      } ;
      // ======================================================================
      template <class T, unsigned int D>
      struct Sym<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > >
      {
        //
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M ;
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > R ;
        //
        static const R& operation ( const M& m ) { return m ; }
      } ;
      // ======================================================================
      template <class T, unsigned int D>
      struct ASym<ROOT::Math::SMatrix<T,D,D> >
      {
        //
        typedef ROOT::Math::SMatrix<T,D,D> M ;
        typedef ROOT::Math::SMatrix<T,D,D> R ;
        //
        static R operation ( const M& m )
        {
          //
          R r ;
          for ( unsigned int i = 0 ; i < D ; ++i )
          {
            r ( i , i ) = 0 ;
            for ( unsigned int j = i + 1 ; j < D ; ++j )
            {
              const T v = 0.5 * ( m ( i , j ) - m ( j , i ) ) ;
              r ( i , j ) =  v ;
              r ( j , i ) = -v ;          
            } 
          }
          //
          return r ;
        }
      } ;
      // ======================================================================
      template <class T, unsigned int D>
      struct ASym<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > >
      {
        //
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M ;
        typedef ROOT::Math::SMatrix<T,D,D> R ;
        //
        static R operation ( const M& /* m */ ) { return R () ; }  
      } ;      
      // ======================================================================
    } //                                  The end of namespace Ostap::Math::Ops
    // ========================================================================
  } //                                        The end of namespace  Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // LHCBMATH_MATRIXUTILS2_H
// ============================================================================
