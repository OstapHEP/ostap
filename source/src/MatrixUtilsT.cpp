// ============================================================================
// Include files 
// ============================================================================
// ROOT
// ============================================================================
#include "TVectorT.h"
#include "TMatrixT.h"
#include "TMatrixTSym.h"
#include "TDecompSVD.h"
#include "TMatrixDEigen.h"
#include "TMatrixDSymEigen.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatusCode.h"
#include "Ostap/Math.h"
#include "Ostap/MatrixUtilsT.h"
// ============================================================================
// local
// ============================================================================
#include "format.h"
#include "status_codes.h"
#include "local_math.h"
// ============================================================================
/** @file 
 *  Non-inline functions from Ostap/MatrixUtilsT.h 
 *  @date 2020-09-22 
 *  @author Vanya BELYAEV IvanBelyaev@iter.ru
 */
// ============================================================================
// get the rank of symmetrix  matrix
// ===========================================================================
std::size_t Ostap::Math::rank
( const TMatrixTSym<float>&     matrix , 
  const float                   eps    )
{
  if ( !matrix.IsValid () || matrix.GetNcols() < 1 ) { return 0 ; }
  /// convert to double 
  TMatrixTSym<double> m { matrix } ;
  return rank ( m , static_cast<double> ( eps ) ) ;
}
// ===========================================================================
// get the rank of symmetrix  matrix
// ===========================================================================
std::size_t Ostap::Math::rank
( const TMatrixTSym<double>&    matrix , 
  const double                  eps    )
{
  if ( !matrix.IsValid () || matrix.GetNcols() < 1 ) { return 0 ; }
  // Eigenvalues 
  TMatrixDSymEigen eigen ( matrix ) ;
  return norm_L0 ( eigen.GetEigenValues() , eps ) ;
}
// ===========================================================================
// get the rank of the matrix
// ===========================================================================
std::size_t Ostap::Math::rank
( const TMatrixT<float>&     matrix , 
  const float                eps    )
{
  if ( !matrix.IsValid () || matrix.GetNcols() < 1 || matrix.GetNrows () < 1 ) { return 0 ; }  
  /// convert to double 
  TMatrixT<double> m { matrix } ;
  return rank ( m , static_cast<double> ( eps ) ) ;
}
// ===========================================================================
// get the rank of the matrix
// ===========================================================================
std::size_t Ostap::Math::rank
( const TMatrixT<double>&    matrix , 
  const double               eps    )
{
  if ( !matrix.IsValid () || matrix.GetNcols() < 1 || matrix.GetNrows () < 1 ) { return 0 ; }  
  /// SVD decomposition
  TDecompSVD svd ( matrix );
  if ( !svd.Decompose () ) { return 0 ; }
  //
  return norm_L0 ( svd.GetSig() ) ; 
}
// ============================================================================

// ============================================================================
/*  @brief Compute spectral norm (maximum singular value) 
 *  @param m (INPUT) Input general matrix
 *  @return Spectral norm value
 */
// ============================================================================
double Ostap::Math::norm_spectral
( const TMatrixT<float>&     matrix )
{
  if ( !matrix.IsValid () || matrix.GetNcols() < 1 || matrix.GetNrows () < 1 ) { return INVALID_NORM_v ; } 
  /// convert to double 
  TMatrixT<double> m { matrix } ;
  return norm_spectral ( m ) ;
}   
// ============================================================================
/*  @brief Compute spectral norm (maximum singular value) 
 *  @param m (INPUT) Input general matrix
 *  @return Spectral norm value
 */
// ============================================================================
double Ostap::Math::norm_spectral
( const TMatrixT<double>&     matrix )
{
  if ( !matrix.IsValid () || matrix.GetNcols() < 1 || matrix.GetNrows () < 1 ) { return INVALID_NORM_v ; } 
  /// SVD decomposition
  TDecompSVD svd ( matrix );
  if ( !svd.Decompose () ) { return INVALID_NORM_v ; }

  //
  return norm_Linf ( svd.GetSig() ) ;    
}
// ============================================================================
/*  @brief Compute spectral norm (maximum singular value) 
 *  @param m (INPUT) Input general matrix
 *  @return Spectral norm value
 */
// ============================================================================
double Ostap::Math::norm_spectral
( const TMatrixTSym<float>&     matrix )
{
  if ( !matrix.IsValid () || matrix.GetNcols() < 1 || matrix.GetNrows () < 1 ) { return INVALID_NORM_v ; } 
  /// convert to double 
  TMatrixTSym<double> m { matrix } ;
  return norm_spectral ( m ) ;  
}
// ============================================================================
/*  @brief Compute spectral norm (maximum singular value) 
 *  @param m (INPUT) Input general matrix
 *  @return Spectral norm value
 */
// ============================================================================
double Ostap::Math::norm_spectral
( const TMatrixTSym<double>&     matrix )
{
  if ( !matrix.IsValid () || matrix.GetNcols() < 1 || matrix.GetNrows () < 1 ) { return INVALID_NORM_v ; } 
  // 
  TMatrixDSymEigen eigen ( matrix ) ;
  return norm_Linf ( eigen.GetEigenValues () ) ;
}


// ============================================================================
/*  @brief Compute nuclear norm (sum of singular value) 
 *  @param m (INPUT) Input general matrix
 *  @return Nuclear norm value
 */
// ============================================================================
double Ostap::Math::norm_nuclear 
( const TMatrixT<float>&     matrix )
{
  if ( !matrix.IsValid () || matrix.GetNcols() < 1 || matrix.GetNrows () < 1 ) { return INVALID_NORM_v ; } 
  /// convert to double 
  TMatrixT<double> m { matrix } ;
  return norm_nuclear ( m ) ;
}   
// ============================================================================
/*  @brief Compute nuclear norm (sum of singular value) 
 *  @param m (INPUT) Input general matrix
 *  @return Nuclear norm value
 */
// ============================================================================
double Ostap::Math::norm_nuclear 
( const TMatrixTSym<float>&     matrix )
{
  if ( !matrix.IsValid () || matrix.GetNcols() < 1 || matrix.GetNrows () < 1 ) { return INVALID_NORM_v ; } 
  /// convert to double 
  TMatrixTSym<double> m { matrix } ;
  return norm_nuclear ( m ) ;
}   
// ============================================================================
/*  @brief Compute nuclear norm (sum of singular value) 
 *  @param m (INPUT) Input general matrix
 *  @return Nuclear norm value
 */
// ============================================================================
double Ostap::Math::norm_nuclear 
( const TMatrixT<double>&     matrix )
{
  if ( !matrix.IsValid () || matrix.GetNcols() < 1 || matrix.GetNrows () < 1 ) { return INVALID_NORM_v ; } 
  /// SVD decomposition
  TDecompSVD svd ( matrix );
  if ( !svd.Decompose () ) { return INVALID_NORM_v ; }
  //
  return sum ( svd.GetSig() ) ;    
}   
// ============================================================================
/*  @brief Compute nuclear norm (sum of singular value) 
 *  @param m (INPUT) Input general matrix
 *  @return Nuclear norm value
 */
// ============================================================================
double Ostap::Math::norm_nuclear 
( const TMatrixTSym<double>&     matrix )
{
  if ( !matrix.IsValid () || matrix.GetNcols() < 1 || matrix.GetNrows () < 1 ) { return INVALID_NORM_v ; } 
  //
  TMatrixDSymEigen eigen ( matrix ) ;
  return sum1 ( eigen.GetEigenValues () ) ; 
}   

// ============================================================================
/*  @brief Compute Schatten' norm \f$ \left(\Sum \left|\sigma_i\right|^p\right)^{1/p}\f$
 *  @param m (INPUT) Input general matrix
 *  @return Schatten's norm value
 */
// ============================================================================
double Ostap::Math::norm_schatten 
( const TMatrixT<float>&    matrix , 
  const double              p      )
{
  if ( !matrix.IsValid () || matrix.GetNcols() < 1 || matrix.GetNrows () < 1 ) { return INVALID_NORM_v ; } 
  /// special cases:
  static const Ostap::Math::Equal_To<double> s_equal {} ;
  static const Ostap::Math::Zero    <double> s_zero  {} ;
  //  
  if      ( std::isinf ( p )            ) { return norm_spectral ( matrix ) ; }
  else if ( 1 == p || s_equal ( 1 , p ) ) { return norm_nuclear  ( matrix ) ; }
  else if ( 2 == p || s_equal ( 2 , p ) ) { return norm_L2       ( matrix ) ; }
  else if ( 0 >= p || s_zero  (     p ) ) { return rank          ( matrix ) ; }
  //
  /// convert to double 
  TMatrixT<double> m { matrix } ;
  return norm_schatten ( m , p ) ;
}
// ============================================================================
/*  @brief Compute Schatten' norm \f$ \left(\Sum \left|\sigma_i\right|^p\right)^{1/p}\f$
 *  @param m (INPUT) Input general matrix
 *  @return Schatten's norm value
 */
// ============================================================================
double Ostap::Math::norm_schatten 
( const TMatrixTSym<float>& matrix , 
  const double              p      )
{
  if ( !matrix.IsValid () || matrix.GetNcols() < 1 || matrix.GetNrows () != matrix.GetNcols () ) { return INVALID_NORM_v ; } 
  /// convert to double 
  TMatrixTSym<double> m { matrix } ;
  return norm_schatten ( m , p ) ;
}
// ============================================================================
/*  @brief Compute Schatten' norm \f$ \left(\Sum \left|\sigma_i\right|^p\right)^{1/p}\f$
 *  @param m (INPUT) Input general matrix
 *  @return Schatten's norm value
 */
// ============================================================================
double Ostap::Math::norm_schatten 
( const TMatrixT<double>&   matrix , 
  const double              p      )
{
  if ( !matrix.IsValid () || matrix.GetNcols() < 1 || matrix.GetNrows () < 1 ) { return INVALID_NORM_v ; } 
  ///
  static const Ostap::Math::Equal_To<double> s_equal {} ;
  static const Ostap::Math::Zero    <double> s_zero  {} ;
  //  
  if      ( std::isinf ( p )            ) { return norm_spectral ( matrix ) ; }
  else if ( 1 == p || s_equal ( 1 , p ) ) { return norm_nuclear  ( matrix ) ; }
  else if ( 2 == p || s_equal ( 2 , p ) ) { return norm_L2       ( matrix ) ; }
  else if ( 0 >= p || s_zero  (     p ) ) { return rank          ( matrix ) ; }
  //  
  /// SVD decomposition
  TDecompSVD svd ( matrix );
  if ( !svd.Decompose () ) { return INVALID_NORM_v ; }  
  //
  return std::pow ( sum_pow ( svd.GetSig() , p ) , 1 / p ) ;
}

// ============================================================================
/* Compute Moore-Penrose Pseudoinverse using SVD: A^+ = V * Sigma^+ * U^T
 *  @param a     (INPUT)  Input matrix A (m x n)
 *  @param a_pinv(OUTPUT) Pseudoinverse matrix A^+ (n x m)
 *  @param tol   (INPUT)  Tolerance for zeroing small singular values ( 0>= for default)
 *  @return status code
 */
// ============================================================================
Ostap::StatusCode Ostap::Math::PINV
( const TMatrixT<float>& a      ,
  TMatrixT<float>&       a_pinv ,
  const float            tol    )
{
  if ( !a.IsValid () || a.GetNcols() < 1 || a.GetNrows () < 1 ) { return INVALID_TMATRIX ; }
  ///
  const Int_t M = a.GetNrows () ;
  const Int_t N = a.GetNcols () ;
  //
  /// convert to double 
  TMatrixT<double> b  { a      } ;
  TMatrixT<double> bi { N , M  } ;
  const Ostap::StatusCode sc = PINV ( b , bi , static_cast<double> ( tol ) ) ;
  if ( sc.isFailure() ) { return sc ; }
  //
  a_pinv.ResizeTo ( N , M ) ; 
  a_pinv = bi ;
  return sc ;   
}
// ============================================================================
/* Compute Moore-Penrose Pseudoinverse using SVD: A^+ = V * Sigma^+ * U^T
 *  @param a     (INPUT)  Input matrix A (m x n)
 *  @param a_pinv(OUTPUT) Pseudoinverse matrix A^+ (n x m)
 *  @param tol   (INPUT)  Tolerance for zeroing small singular values ( 0>= for default)
 *  @return status code
 */
// ============================================================================
Ostap::StatusCode Ostap::Math::PINV
( const TMatrixT<double>& a      ,
  TMatrixT<double>&       a_pinv ,
  const double            tol    )
{
  //
  if ( !a.IsValid () || a.GetNcols() < 1 || a.GetNrows () < 1 ) { return INVALID_TMATRIX ; }
  //
  const std::size_t M  = a.GetNrows () ;
  const std::size_t N  = a.GetNcols () ;
  const std::size_t K  = std::min ( M , N ) ; 
  
  // SVD decomposition: A = U * Sigma * V^T
  TDecompSVD svd ( a ) ;
  if ( !svd.Decompose () ) { return INVALID_SVD_DECOMPOSITION ; }
  //  
  const auto& U = svd.GetU   () ;
  const auto& V = svd.GetV   () ;
  const auto& S = svd.GetSig () ;
  
  // Determine default tolerance if not specified (tol <= 0)
  const double max_sig = std::abs ( norm_max ( S ) ) ;
  const double atol    = ( 0 < tol ? tol : std::numeric_limits<double>::epsilon() ) * max_sig * std::max ( M , N ) ;
  
  // Construct Sigma^+ (n x m matrix)
  TMatrixD sigma_pinv ( N , M ) ;
  sigma_pinv.Zero();
  
  for ( std::size_t i = 0 ; i < K ; ++i )
  {
    const double si = S [ i ] ;
    if ( si && ( atol < std::abs ( si ) ) ) { sigma_pinv ( i , i ) = 1.0 / si ; }
  }
  
  // A^+ = V * Sigma^+ * U^T
  // V is (n x n), sigma_pinv is (n x m), U^T is (m x m)
  TMatrixD u_transpose ( TMatrixD::kTransposed , U );
  
  // Resize output matrix to (n x m)
  a_pinv.ResizeTo ( N , M ) ;
  a_pinv = V * sigma_pinv * u_transpose;
  //
  return Ostap::StatusCode::SUCCESS;
}
// ============================================================================
/* Compute Moore-Penrose Pseudoinverse using SVD: A^+ = V * Sigma^+ * U^T
 *  @param a     (INPUT)  Input matrix A (m x n)
 *  @param a_pinv(OUTPUT) Pseudoinverse matrix A^+ (n x m)
 *  @param tol   (INPUT)  Tolerance for zeroing small singular values ( 0>= for default)
 *  @return status code
 */
// ============================================================================
Ostap::StatusCode Ostap::Math::PINV
( const TMatrixTSym<float>& a      ,
  TMatrixTSym<float>&       a_pinv ,
  const float               tol    )  
{
  //
  if ( !a.IsValid () || a.GetNcols() < 1 || a.GetNrows () < 1 ) { return INVALID_TMATRIX ; }
  //
  const Int_t N  = a.GetNrows () ;
  //
  TMatrixTSym<double> b  { a } ; 
  TMatrixTSym<double> bi { N , N } ;
  
  const Ostap::StatusCode sc = PINV ( b , bi , static_cast<double> ( tol ) ) ;
  if ( sc.isFailure() ) { return sc ; }
  //
  a_pinv.ResizeTo ( N , N ) ; 
  a_pinv = bi ;
  return sc ;     
}

// ============================================================================
/* Compute Moore-Penrose Pseudoinverse using SVD: A^+ = V * Sigma^+ * U^T
 *  @param a     (INPUT)  Input matrix A (m x n)
 *  @param a_pinv(OUTPUT) Pseudoinverse matrix A^+ (n x m)
 *  @param tol   (INPUT)  Tolerance for zeroing small singular values ( 0>= for default)
 *  @return status code
 */
// ============================================================================
Ostap::StatusCode Ostap::Math::PINV
( const TMatrixTSym<double>& a      ,
  TMatrixTSym<double>&       a_pinv ,
  const double               tol    )  
{
  //
  if ( !a.IsValid () || a.GetNcols() < 1 || a.GetNrows () < 1 ) { return INVALID_TMATRIX ; }
  //
  const std::size_t N  = a.GetNrows () ;
  //
  // Symmetric eigenvalue decomposition: A = V * Lambda * V^T
  TMatrixDSymEigen eigen ( a ) ;
  //
  const auto& V        = eigen.GetEigenVectors () ;
  const auto& lambda   = eigen.GetEigenValues  () ;
  
  const double max_sig = std::abs ( norm_max ( lambda ) ) ;
  const double atol    = ( 0 < tol ? tol : std::numeric_limits<double>::epsilon() ) * max_sig * N ;
  
  // Construct Lambda^+ (diagonal matrix of inverted eigenvalues)
  TMatrixTSym<double> lambda_pinv ( N , N  ) ;  
  lambda_pinv.Zero();
  for ( std::size_t i = 0 ; i < N; ++i )
  {
    const double li = lambda [ i ] ;
    if ( li && ( atol < std::abs ( li ) ) ) { lambda_pinv(i, i) = 1.0 / li ; }
  }
  //  
  a_pinv = lambda_pinv.Similarity ( V ) ;
  //
  return Ostap::StatusCode::SUCCESS; 
}

// =============================================================================
namespace
{
  // ===========================================================================
  /** @brief Compute Variance Inflation Factors (VIF) for a TMatrixTSym covariance matrix.
   *
   *  Calculates the Variance Inflation Factor (VIF) vector \f$ \vec{v} \f$ 
   *  for a ROOT ROOT::TMatrixTSym<T> covariance matrix \f$ \Sigma \f$:
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
   *  1. **Fast Path**: Attempts in-place fast Cholesky/LU inversion on a local copy. 
   *     Optimal for well-behaved Positive Definite matrices.
   *  2. **Fallback Path**: If fast inversion fails (due to zero or negative eigenvalues 
   *     arising from exact linear dependencies or negative \f$sPlot\f$ event weights), 
   *     it gracefully falls back to spectral pseudoinversion (`PINV`). Zero and negative 
   *     eigenvalues (\f$\lambda_k \le \text{eps} \cdot \lambda_{\max}\f$) are zeroed out.
   *  3. **Non-positive Variances**: Variables with \f$ \Sigma_{ii} \le 0 \f$ (constants or 
   *     severe \f$sPlot\f$ noise) are explicitly assigned `std::numeric_limits<T>::infinity()`, 
   *     marking them as primary targets for elimination.
   *
   *  @tparam T Data type (`double`, `float`).
   *  @param[in]  cov Input symmetric covariance matrix \f$ \Sigma \f$ (TMatrixTSym<T>).
   *  @param[out] vif Output vector containing VIF values (TVectorT<T>).
   *  @param[in]  eps Numerical tolerance ratio for truncating small eigenvalues in `PINV`.
   *  @return `Ostap::StatusCode::SUCCESS` if computation completed successfully.
   */
   template <typename T>
   inline Ostap::StatusCode _VIF_ 
   ( const TMatrixTSym<T>& cov  ,
     TVectorT<T>&          vif  ,
     const T               eps = std::numeric_limits<T>::epsilon() ) 
    {
      //
      if ( !cov.IsValid () || cov.GetNcols() < 1 || cov.GetNrows () != cov.GetNcols () ) 
      { return INVALID_TMATRIX ; } 
      // 
      const Int_t nrows = cov.GetNrows() ;

      // Resize output vector if necessary
      if ( vif.GetNrows() != nrows ) { vif.ResizeTo( nrows ) ; }

      // Working copy for inversion
      TMatrixTSym<T> sp { cov } ;
      
      // Fast path inversion (Cholesky / Fast inversion)
      sp.InvertFast() ;
      
      // Check if InvertFast failed (sp becomes invalid or non-invertible)
      if ( !sp.IsValid() )
      {
        // Reset copy to original matrix before pseudoinversion
        sp = cov ;
        Ostap::StatusCode sc = Ostap::Math::PINV ( cov , sp , eps ) ; 
        if ( sc.isFailure () ) { return sc ; }
      }

      // Compute VIF array
      for ( Int_t i = 0 ; i < nrows ; ++i )
      {
        const T cii = cov ( i , i ) ;
        vif [ i ]   = cii <= 0 ? std::numeric_limits<T>::infinity () : cii * sp ( i , i ) ;
      }
      //
      return Ostap::StatusCode::SUCCESS ;
    }
  // ==========================================================================
}
// ============================================================================

Ostap::StatusCode Ostap::Math::VIF 
( const TMatrixTSym<float>& cov ,
  TVectorT<float>&          vif ,
  const float               eps )
{ return _VIF_ ( cov , vif , eps ) ; } 
Ostap::StatusCode VIF 
( const TMatrixTSym<double>& cov ,
  TVectorT<double>&          vif ,
  const double               eps ) 
{ return _VIF_ ( cov , vif , eps ) ; } 

// ============================================================================
//                                                                      The END 
// ============================================================================
