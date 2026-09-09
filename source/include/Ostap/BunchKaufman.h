#ifndef OSTAP_BUNCHKAUFMAN_H 
#define OSTAP_BUNCHKAUFMAN_H 1
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
#include "TVectorT.h"
#include "TMatrixT.h"
#include "TMatrixTSym.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatusCode.h"
// ============================================================================
/** @file Ostap/BunchKaufman.h
 *  Bunhc-Kaufman decomposition of symmetric matrices
 *  @see TDecompBK
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 */
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    /** @class TPermutation
     *  Trival permutation vector
     */
    /*
    class TPermutation
    {
    public  :
      // ======================================================================
      /// constructor  from vector of indices
      TPermutation ( const std::vector<std::size_t>& indices ) ;
      // ======================================================================
      TPermutation () = delete ;
      // ======================================================================
    public :
      // ======================================================================
      /// get the size
      inline std::size_t                     size    () const { return m_indices.size () ; }
      //// get indices
      inline const std::vector<std::size_t>& indices () const { return m_indices         ; }
      // ======================================================================
    private :
      // ======================================================================
      std::vector<std::size_t> m_indices ;
      // ======================================================================
    };
    */
    
    // ========================================================================
    /** @brief Applies permutation vector p to symmetric matrix M: A = P * M * P^T
     *  @param [in]  M Input symmetric matrix N x N
     *  @param [in]  p Permutation vector of size N
     *  @param [out] A Permuted symmetric matrix A
     *  @return status code 
     */
    Ostap::StatusCode ApplyPermutation
      ( const TMatrixTSym<double>& M ,
        const TVectorT<double>&    P ,
        TMatrixTSym<double>&       A ) ;
    
    // ========================================================================
    /** @brief Applies permutation vector p to symmetric matrix M: A = P * M * P^T
     *  @param [in]  M Input symmetric matrix N x N
     *  @param [in]  p Permutation vector of size N
     *  @param [out] A Permuted symmetric matrix A
     */
    Ostap::StatusCode ApplyPermutation
      ( const TMatrixTSym<float>& M ,
        const TVectorT<float>&    P ,
        TMatrixTSym<float>&       A ) ;
    
    // ========================================================================
    /** @brief Convert a permutation vector \f$ p \f$  into an explicit permutation matrix \f$ P \f$ 
     *
     *  Constructs an \f$ N \times N \f$ matrix \f$ P \f$,
     *  where \f$ P( i , p ( i ) ) = 1.0 \f$ and all other 
     *  entries are zero 
     *
     *  @param[in]  p Input permutation vector of size \f$ N \f$ .
     *  @param[out] P Output  \f$ N \times N \f$  orthogonal permutation matrix.
     *  @return status code.
     */
    Ostap::StatusCode PermutationMatrix
    ( const TVectorT<double>& p ,
      TMatrixT<double>&       P ) ; 
    
    // ========================================================================
    /** @brief Convert a permutation vector \f$ p \f$  into an explicit permutation matrix \f$ P \f$ 
     *
     *  Constructs an \f$ N \times N \f$ matrix \f$ P \f$,
     *  where \f$ P( i , p ( i ) ) = 1.0 \f$ and all other 
     *  entries are zero 
     *
     *  @param[in]  p Input permutation vector of size \f$ N \f$ .
     *  @param[out] P Output  \f$ N \times N \f$  orthogonal permutation matrix.
     *  @return status code.
     */
     Ostap::StatusCode PermutationMatrix
     ( const TVectorT<float>& p ,
       TMatrixT<float>&       P ) ; 
    
    // ========================================================================    
    /** Bunch-Kaufman decompositon of symmetric matrices 
     * @param[in]  A Input symmetric matrix to decompose.
     * @param[out] U triangular factor matrix U.
     * @param[out] D block-diagonal symmetric matrix D containing 1x1 and 2x2 blocks.
     * @return status code 
     */
    Ostap::StatusCode BunchKaufman
    ( const TMatrixTSym<double>& A ,
      TMatrixT<double>&          U ,
      TMatrixTSym<double>&       D ) ;
  
    // ========================================================================    
    /** Bunch-Kaufman decompositon of symmetric matrices 
     * @param[in]  A Input symmetric matrix to decompose.
     * @param[out] U triangular factor matrix U.
     * @param[out] D block-diagonal symmetric matrix D containing 1x1 and 2x2 blocks.
     * @return status code 
     */
    Ostap::StatusCode BunchKaufman
    ( const TMatrixTSym<float>& A ,
      TMatrixT<float>&          U ,
      TMatrixTSym<float>&       D ) ;
    
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_BUNCHKAUFMAN_H
// ============================================================================
