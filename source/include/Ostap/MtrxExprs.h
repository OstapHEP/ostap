// ============================================================================
#ifndef OSTAP_MATRIX_EXPR_H
#define OSTAP_MATRIX_EXPR_H 1
// ============================================================================
// Include files
// ============================================================================
// ROOT 
// ============================================================================
#include "Math/SMatrix.h"
#include  "Math/SVector.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    /** @class DiagMatrixExpr
     *  @brief Lightweight expression adapter for representing a diagonal matrix from an SVector.
     * 
     *  This class wraps a reference to a ROOT::Math::SVector and exposes a matrix-like
     *  interface with $O(1)$ lazy evaluation. It avoids memory allocation and dynamic 
     *  copying of vector data.
     * 
     *  @tparam T      Numeric data type (e.g., double, float).
     *  @tparam D_dim  Dimension of the square diagonal matrix (\f$D\_dim \times D\_dim\f$).
     *
     *  @author Vanya BELYAEV
     *  @date   2026
     */
    template <typename T, unsigned int D_dim>
    class DiagMatrixExpr
    {
    public:
      /// Data type contained in the matrix expression
      using value_type = T;

      /** @brief Construct a diagonal matrix expression wrapping a constant vector reference.
       *  @param[in] vec The vector representing the diagonal elements.
       */
      explicit DiagMatrixExpr ( const ROOT::Math::SVector<T, D_dim>& vec ) 
        : m_vec ( vec ) {}

      /** @brief Get a constant reference to the underlying diagonal vector.
       *  @return Constant reference to ROOT::Math::SVector.
       */
      inline const ROOT::Math::SVector<T, D_dim>& vector() const { return m_vec; }

      /** @brief Element access operator simulating matrix indexing \f$M(i, j)\f$.
       *  @param[in] i Row index.
       *  @param[in] j Column index.
       *  @return Diagonal element \f$v_i\f$ if \f$i == j\f$, otherwise \f$0\f$.
       */
      inline T operator() ( unsigned int i, unsigned int j ) const 
      {
        return ( i == j ) ? m_vec[i] : T(0) ;
      }

      /// Number of rows in the diagonal matrix
      static constexpr unsigned int kRows = D_dim;
      /// Number of columns in the diagonal matrix
      static constexpr unsigned int kCols = D_dim;

    private:
      /// Constant reference to the underlying vector
      const ROOT::Math::SVector<T, D_dim>& m_vec ;
    };

    // ========================================================================
    /** @brief Factory helper function to instantiate a DiagMatrixExpr with automatic template argument deduction.
     *  @tparam T      Numeric data type.
     *  @tparam D_dim  Dimension of the vector.
     *  @param[in] v   The SVector to adapt as a diagonal matrix.
     *  @return A DiagMatrixExpr wrapping vector @p v.
     */
    template <typename T, unsigned int D_dim>
    inline DiagMatrixExpr<T, D_dim> 
    Diag ( const ROOT::Math::SVector<T, D_dim>& v ) 
    {
      return DiagMatrixExpr<T, D_dim>( v );
    }

    // ========================================================================
    /** @class DiagMulExpr
     *  @brief Lazy expression evaluation object for matrix-diagonal multiplication (\f$M \cdot D\f$ or \f$D \cdot M\f$).
     * 
     *  Implements element-wise computation directly without intermediate matrix copies, 
     *  reducing multiplication complexity from \f$O(N^3)\f$ to \f$O(N^2)\f$.
     * 
     *  @tparam A         Expression representation type inside ROOT::Math::Expr.
     *  @tparam T         Numeric data type.
     *  @tparam D1        Number of rows in the input matrix.
     *  @tparam D2        Number of columns in the input matrix / dimension of diagonal matrix.
     *  @tparam R         Matrix storage representation type (e.g., MatRepStd).
     *  @tparam RightDiag Flag indicating if diagonal is on the right (@c true for \f$M \cdot D\f$, @c false for \f$D \cdot M\f$).
     */
    template <typename A, typename T, unsigned int D1, unsigned int D2, typename R, bool RightDiag = true>
    class DiagMulExpr
    {
    public:
      /// Data type contained in the matrix expression
      using value_type = T;
      /// Alias for the specific ROOT::Math::Expr being wrapped
      using MatrixExpr = ROOT::Math::Expr<A, T, D1, D2, R>;

      /** @brief Constructor initializing references to the matrix expression and diagonal adapter.
       *  @param[in] mat  Reference to the ROOT matrix expression.
       *  @param[in] diag Reference to the diagonal matrix adapter expression.
       */
      DiagMulExpr ( const MatrixExpr&            mat  , 
                    const DiagMatrixExpr<T, D2>& diag )
        : m_mat ( mat ) , m_diag ( diag ) {}

      /** @brief Lazy evaluation operator for element \f$(i, j)\f$.
       *  @param[in] i Row index.
       *  @param[in] j Column index.
       *  @return Calculated element value \f$M_{ij} \cdot d_j\f$ or \f$d_i \cdot M_{ij}\f$.
       */
      inline T operator() ( unsigned int i, unsigned int j ) const 
      {
        return RightDiag ? ( m_mat(i, j) * m_diag.vector()[j] ) 
                         : ( m_diag.vector()[i] * m_mat(i, j) ) ;
      }

      /// Number of rows of the resulting expression
      static constexpr unsigned int kRows = D1;
      /// Number of columns of the resulting expression
      static constexpr unsigned int kCols = D2;

    private:
      /// Constant reference to the input matrix expression
      const MatrixExpr&            m_mat  ;
      /// Constant reference to the diagonal matrix adapter
      const DiagMatrixExpr<T, D2>& m_diag ;
    };

    // ========================================================================
    /** @brief Lazy multiplication operator: \f$M \cdot D\f$ (Matrix Expression * Diagonal).
     *  @tparam A  Expression template type.
     *  @tparam T  Data scalar type.
     *  @tparam D1 Number of rows in matrix M.
     *  @tparam D2 Number of columns in matrix M / dimension of diagonal D.
     *  @tparam R  Storage representation.
     *  @param[in] M Input ROOT::Math::Expr matrix.
     *  @param[in] D Input diagonal matrix expression.
     *  @return A wrapped ROOT::Math::Expr encapsulating lazy multiplication.
     */
    template <typename A, typename T, unsigned int D1, unsigned int D2, typename R>
    inline ROOT::Math::Expr<
        DiagMulExpr<A, T, D1, D2, R, true>,
        T, D1, D2,
        ROOT::Math::MatRepStd<T, D1, D2>
    >
    operator* ( const ROOT::Math::Expr<A, T, D1, D2, R>& M ,
                const DiagMatrixExpr<T, D2>&             D )
    {
      using RhsExpr = DiagMulExpr<A, T, D1, D2, R, true>;
      using RepType = ROOT::Math::MatRepStd<T, D1, D2>;
      
      return ROOT::Math::Expr<RhsExpr, T, D1, D2, RepType>( RhsExpr( M , D ) );
    }

    // ========================================================================
    /** @brief Lazy multiplication operator: \f$D \cdot M\f$ (Diagonal * Matrix Expression).
     *  @tparam A  Expression template type.
     *  @tparam T  Data scalar type.
     *  @tparam D1 Dimension of diagonal D / number of rows in matrix M.
     *  @tparam D2 Number of columns in matrix M.
     *  @tparam R  Storage representation.
     *  @param[in] D Input diagonal matrix expression.
     *  @param[in] M Input ROOT::Math::Expr matrix.
     *  @return A wrapped ROOT::Math::Expr encapsulating lazy multiplication.
     */
    template <typename A, typename T, unsigned int D1, unsigned int D2, typename R>
    inline ROOT::Math::Expr<
        DiagMulExpr<A, T, D1, D2, R, false>,
        T, D1, D2,
        ROOT::Math::MatRepStd<T, D1, D2>
    >
    operator* ( const DiagMatrixExpr<T, D1>&             D ,
                const ROOT::Math::Expr<A, T, D1, D2, R>& M )
    {
      using RhsExpr = DiagMulExpr<A, T, D1, D2, R, false>;
      using RepType = ROOT::Math::MatRepStd<T, D1, D2>;
      
      return ROOT::Math::Expr<RhsExpr, T, D1, D2, RepType>( RhsExpr( M , D ) );
    }

    // ========================================================================
    /** @brief Overloaded operator* for direct ROOT::Math::SMatrix (Right multiplication: \f$M \cdot D\f$).
     *  Converts SMatrix to ROOT::Math::Expr via @c .asExpr() automatically.
     *  @param[in] M Input dense SMatrix.
     *  @param[in] D Input diagonal matrix expression.
     */
    template <typename T, unsigned int D1, unsigned int D2, typename R>
    inline auto operator* ( const ROOT::Math::SMatrix<T, D1, D2, R>& M ,
                            const DiagMatrixExpr<T, D2>&             D )
    {
      return M.asExpr() * D;
    }

    // ========================================================================
    /** @brief Overloaded operator* for direct ROOT::Math::SMatrix (Left multiplication: \f$D \cdot M\f$).
     *  Converts SMatrix to ROOT::Math::Expr via @c .asExpr() automatically.
     *  @param[in] D Input diagonal matrix expression.
     *  @param[in] M Input dense SMatrix.
     */
    template <typename T, unsigned int D1, unsigned int D2, typename R>
    inline auto operator* ( const DiagMatrixExpr<T, D1>&             D ,
                            const ROOT::Math::SMatrix<T, D1, D2, R>& M )
    {
      return D * M.asExpr();
    }

    // ========================================================================
    /** @brief Fast element-wise vector scaling operator: \f$D \cdot v\f$ (Hadamard product).
     *  Performs linear time \f$O(N)\f$ diagonal-vector multiplication.
     *  @tparam T     Data scalar type.
     *  @tparam D_dim Dimension of vector and matrix.
     *  @param[in] D  Input diagonal matrix adapter.
     *  @param[in] v  Input vector.
     *  @return Resulting vector containing element-wise product \f$d_i \cdot v_i\f$.
     */
    template <typename T, unsigned int D_dim>
    inline ROOT::Math::SVector<T, D_dim> 
    operator* ( const DiagMatrixExpr<T, D_dim>&       D ,
                const ROOT::Math::SVector<T, D_dim>& v ) 
    {
      ROOT::Math::SVector<T, D_dim> result ;
      const auto& d = D.vector() ;
      for ( unsigned int i = 0 ; i < D_dim ; ++i ) 
      {
        result[i] = d[i] * v[i] ;
      }
      return result ;
    }

    // ========================================================================
    /** @brief Specialized similarity transformation \f$R = M \cdot D \cdot M^T\f$ for ROOT::Math::Expr.
     * 
     *  Computes the transformation directly into a symmetric matrix representation (@c MatRepSym),
     *  reducing the operation count to \f$O(D_1^2 \cdot D_2)\f$ while skipping off-diagonal zero terms.
     * 
     *  @tparam A  Expression representation type.
     *  @tparam T  Numeric data type.
     *  @tparam D1 Number of rows in matrix M (dimension of output symmetric matrix).
     *  @tparam D2 Number of columns in matrix M / dimension of diagonal matrix D.
     *  @tparam R  Storage representation type.
     *  @param[in] M Input matrix expression.
     *  @param[in] D Diagonal matrix expression.
     *  @return Symmetric matrix of size \f$D_1 \times D_1\f$.
     */
    template <typename A, typename T, unsigned int D1, unsigned int D2, typename R>
    inline auto Similarity ( const ROOT::Math::Expr<A, T, D1, D2, R>& M ,
                             const DiagMatrixExpr<T, D2>&              D ) 
      -> ROOT::Math::SMatrix<T, D1, D1, ROOT::Math::MatRepSym<T, D1>>
    {
      ROOT::Math::SMatrix<T, D1, D1, ROOT::Math::MatRepSym<T, D1>> result ;
      const auto& d = D.vector() ;

      for ( unsigned int i = 0 ; i < D1 ; ++i ) 
      {
        for ( unsigned int j = i ; j < D1 ; ++j ) 
        {
          T sum = 0 ;
          for ( unsigned int k = 0 ; k < D2 ; ++k ) 
          {
            sum += M(i, k) * d[k] * M(j, k) ;
          }
          result(i, j) = sum ;
        }
      }
      return result ;
    }

    // ========================================================================
    /** @brief Overload of Similarity transformation for direct ROOT::Math::SMatrix inputs.
     *  @param[in] M Input dense SMatrix.
     *  @param[in] D Input diagonal matrix adapter.
     */
    template <typename T, unsigned int D1, unsigned int D2, typename R>
    inline auto Similarity ( const ROOT::Math::SMatrix<T, D1, D2, R>& M ,
                             const DiagMatrixExpr<T, D2>&             D )
    {
      return Similarity( M.asExpr(), D );
    }

    // ========================================================================
    /** @brief Symmetric syntax overload for Similarity transformation: @c Similarity(Diag(d), M).
     *  Delegates call to @c Similarity(M, Diag(d)).
     *  @param[in] D Input diagonal matrix adapter.
     *  @param[in] M Input matrix (SMatrix or Expr).
     */
    template <typename MatrixType, typename T, unsigned int D2>
    inline auto Similarity ( const DiagMatrixExpr<T, D2>& D ,
                             const MatrixType&            M ) 
      -> decltype(Similarity(M, D))
    {
      return Similarity ( M , D ) ;
    }

    // =======================================================================
  } //                                         The end of namespace Otap::Math
  // =========================================================================
} //                                                The end of namespace Ostap
// ===========================================================================
#endif // OSTAP_MATRIX_EXPR_H
// ===========================================================================
//                                                                     The END
// =========================================================================== 

