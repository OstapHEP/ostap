// ============================================================================
// Include files 
// ============================================================================
// STD&STL 
// ============================================================================
#include <vector>
#include <utility>
// ============================================================================
// ROOT
// ============================================================================
#include "TVectorT.h"
#include "TMatrixT.h"
#include "TMatrixTSym.h"
#include "TDecompBK.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatusCode.h"
#include "Ostap/Math.h"
#include "Ostap/BunchKaufman.h"
// ============================================================================
// local
// ============================================================================
#include "status_codes.h"
// ============================================================================
/** @file
 *  Implementation of Bunch-Kaufman decomposition
 *  @date 2020-09-22 
 *  @author Vanya BELYAEV IvanBelyaev@iter.ru
 */
// ============================================================================
namespace
{
  // =========================================================================
  /// check if vector p corresponds to a valid permutation 
  template <typename T>
  inline bool valid_permutation
  ( const TVectorT<T> & p )
  {
    if ( !p.IsValid() || p.GetNrows() < 1 ) { return false ; }
    //
    const std::size_t N = p.GetNrows() ;
    std::vector<bool> seen ( N , false ) ;
    //
    for ( Int_t i = 0 ; i < N ; ++i )
    {      
      const T v = p ( i ) ;
      if ( !Ostap::Math::isuint ( v )   ) { return false ; }
      const std::size_t index = Ostap::Math::round ( v  ) ;
      if ( N <= index || seen [ index ] ) { return false ; }
      seen [ index ] = true ;
    }
    //
    return true ;
  }
  // =========================================================================
  /** @brief Applies permutation vector p to symmetric matrix M: A = P * M * P^T
   *  @param [in]  M Input symmetric matrix N x N
   *  @param [in]  p Permutation vector of size N
   *  @param [out] A Permuted symmetric matrix A
   */
  // ============================================================================
  template <typename T> 
  inline Ostap::StatusCode _AP_ 
  ( const TMatrixTSym<T>& M ,
    const TVectorT<T>&    P ,
    TMatrixTSym<T>&       A )
  {
    //
    /// Validate arguments
    if ( !M.IsValid() || M.GetNrows() < 1 || M.GetNrows() != M.GetNcols () ) { return INVALID_TMATRIX     ; }
    if ( !P.IsValid() || P.GetNrows() < 1 || P.GetNrows() != M.GetNrows () ) { return INVALID_TVECTOR     ; }
    if ( !valid_permutation ( P )                                          ) { return INVALID_PERMUTATION ; }
    //
    const Int_t N = M.GetNrows() ;
    //
    // Handle in-place operation (aliasing protection when &A == &M)
    if ( &A == &M ) 
    {
      TMatrixTSym<T> tmp ( N ) ;
      Ostap::StatusCode sc = _AP_ ( M , P , tmp ) ;
      if ( sc.isSuccess() ) { A = tmp ; }
      return sc ;
    }
    else {  A.ResizeTo ( N , N ) ; } 
    //
    for ( Int_t i = 0 ; i < N ; ++i )
    {
      const Int_t pi = static_cast<Int_t> ( P ( i ) ) ;       
      // Diagonal element
      A ( i , i ) = M ( pi , pi ) ;      
      // Off-diagonal elements
      for ( Int_t j = i + 1 ; j < N ; ++j )
      {
        const Int_t pj  = static_cast<Int_t> ( P ( j ) ) ;            
        const T     val = M ( pi , pj ) ; 
        A ( i , j ) = val ;
        A ( j , i ) = val ; 
      }
    }
    //
    return Ostap::StatusCode::SUCCESS ;
  }
  // =========================================================================
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
  template <typename T>
  inline Ostap::StatusCode _PM_ 
  ( const TVectorT<T>& p ,
    TMatrixT<T>&       P )
  {
    if ( !p.IsValid() || p.GetNrows() < 1 ) { return INVALID_TVECTOR     ; }
    if ( !valid_permutation ( p )         ) { return INVALID_PERMUTATION ; }
    // 
    const Int_t n = p.GetNrows() ;
    //
    P.ResizeTo ( n , n ) ;
    P.Zero     (       ) ;
    //
    for ( Int_t i = 0 ; i < n ; ++i )
    {
      const T     p_i = p ( i ) ;
      /// it must be "integer": check for unsigned int   
      if ( !Ostap::Math::isuint ( p_i )     ) { return INVALID_PERMUTATION_INDEX ; }
      const Int_t idx = static_cast<Int_t>( p ( i ) ) ;
      if ( idx < 0 || idx >= n              ) { return INVALID_PERMUTATION_INDEX ; }
      //
      P ( i , idx ) = 1.0 ;
    }
    //
    return Ostap::StatusCode::SUCCESS ;
  }
  // ==========================================================================
  /**
   * @brief Helper class to access the protected fIpiv array of TDecompBK.
   */
  class TDecompBKSpy : public TDecompBK 
  {
  public:
    // ========================================================================
    using TDecompBK::TDecompBK;    
    /// Get the internal pivoting and block index array.
    const Int_t* GetIpiv() const { return fIpiv; }
    // ========================================================================
  };
  // ==========================================================================
}
// ============================================================================
/** @brief Applies permutation vector p to symmetric matrix M: A = P * M * P^T
 *  @param [in]  M Input symmetric matrix N x N
 *  @param [in]  p Permutation vector of size N
 *  @param [out] A Permuted symmetric matrix A
 *  @return status code 
 */
// ============================================================================
Ostap::StatusCode Ostap::Math::ApplyPermutation
( const TMatrixTSym<double>& M ,
  const TVectorT<double>&    P ,
  TMatrixTSym<double>&       A )
{ return ::_AP_ ( M , P , A ) ; }

// ============================================================================
/** @brief Applies permutation vector p to symmetric matrix M: A = P * M * P^T
 *  @param [in]  M Input symmetric matrix N x N
 *  @param [in]  p Permutation vector of size N
 *  @param [out] A Permuted symmetric matrix A
 */
// ============================================================================
Ostap::StatusCode Ostap::Math::ApplyPermutation
( const TMatrixTSym<float>& M ,
  const TVectorT<float>&    P ,
  TMatrixTSym<float>&       A ) 
{ return ::_AP_ ( M , P , A ) ; }

// ============================================================================
/*  @brief Convert a permutation vector <code>p</code> into an explicit permutation matrix <code>P</code>.
 *
 *  Constructs an \f$ N \times N \f$ matrix \f$ P \f$,
 *  where \f$ P( i , p ( i ) ) = 1.0 \f$ and all other 
 *  entries are zero 
 *
 *  @param[in]  p Input permutation vector of size \f$ N \f$ .
 *  @param[out] P Output  \f$ N \times N \f$  orthogonal permutation matrix.
 *  @return status code.
 */
// ============================================================================
Ostap::StatusCode Ostap::Math::PermutationMatrix
( const TVectorT<double>& p ,
  TMatrixT<double>&       P )
{ return ::_PM_ ( p , P ) ; }
// ============================================================================
/*  @brief Convert a permutation vector \f$ p \f$  into an explicit permutation matrix \f$ P \f$ 
 *
 *  Constructs an \f$ N \times N \f$ matrix \f$ P \f$,
 *  where \f$ P( i , p ( i ) ) = 1.0 \f$ and all other 
 *  entries are zero 
 *
 *  @param[in]  p Input permutation vector of size \f$ N \f$ .
 *  @param[out] P Output  \f$ N \times N \f$  orthogonal permutation matrix.
 *  @return status code.
 */
// ============================================================================
Ostap::StatusCode Ostap::Math::PermutationMatrix
( const TVectorT<float>& p ,
  TMatrixT<float>&       P )
{ return ::_PM_ ( p , P ) ; }
// ============================================================================

// ============================================================================
// Bunch-Kaufman decomposition 
// ============================================================================

// ========================================================================
/* Bunch-Kaufman decomposition of symmetric matrices A = U D U^T
 * @param[in]  A Input symmetric matrix to decompose.
 * @param[out] U matrix U.
 * @param[out] D block-diagonal symmetric matrix D containing 1x1 and 2x2 blocks.
 * @return status code 
 */
// ========================================================================
Ostap::StatusCode Ostap::Math::BunchKaufman
( const TMatrixTSym<float>& A ,
  TMatrixT<float>&          U ,
  TMatrixTSym<float>&       D )
{
  if ( !A.IsValid () || A.GetNrows() < 1 || A.GetNrows () != A.GetNcols () ) { return INVALID_TMATRIX ; } 
  //
  const Int_t N = A.GetNrows() ;
  //
  const TMatrixTSym<double> a { A } ;
  TMatrixT<double>          u { N , N } ;
  TMatrixTSym<double>       d { N } ;
  //
  const Ostap::StatusCode sc = BunchKaufman ( a , u , d ) ;
  if ( sc.isFailure() ) { return sc ; } ;
  //
  U = u ;
  D = d ;
  //
  return Ostap::StatusCode::SUCCESS ;  
}
// ========================================================================
/** @brief Bunch-Kaufman decomposition with implicit permutation \f$ A = U D U^T \f$.
 *
 *  Decomposes a real symmetric matrix \f$ A \f$ into \f$ A = U D U^T \f$.
 *  The matrix \f$ U \f$ explicitly accumulates unit upper triangular multipliers
 *  and row permutations in exact LAPACK sequence order.
 *
 *  @param[in]  A Input real symmetric matrix.
 *  @param[out] U Output factor matrix \f$ U \f$ such that \f$ A = U D U^T \f$.
 *  @param[out] D Output symmetric block-diagonal matrix \f$ D \f$.
 *  @return Ostap::StatusCode Status code (SUCCESS if factorization succeeded).
 */
// ========================================================================
Ostap::StatusCode Ostap::Math::BunchKaufman
( const TMatrixTSym<double>& A ,
  TMatrixT<double>&          U ,
  TMatrixTSym<double>&       D ) 
{
  if ( !A.IsValid() || A.GetNrows() < 1 || A.GetNrows() != A.GetNcols() ) { return INVALID_TMATRIX ; } 
  
  ::TDecompBKSpy bk ( A ) ;
  if ( !bk.Decompose() ) { return INVALID_BK_DECOMPOSITION ; } 
  
  const TMatrixD& rawU = bk.GetU() ;
  const Int_t*    ipiv = bk.GetIpiv() ;
  if ( !ipiv ) { return INVALID_BK_DECOMPOSITION ; }
  
  const Int_t n = A.GetNrows() ;
  
  /// Helper structure to represent LAPACK pivot steps
  struct BKStep {
    Int_t k;     /// Pivot column/row index
    Int_t type;  /// Pivot block type (1 for 1x1, 2 for 2x2)
    Int_t kp;    /// Permutation target index
  };
  
  std::vector<BKStep> steps ;
  steps.reserve( n ) ;
  
  D.ResizeTo ( n , n ) ;
  D.Zero     () ;
  
  /// Parse LAPACK ipiv array and extract block diagonal D
  Int_t k = n - 1 ;
  while ( k >= 0 )
  {
    if ( ipiv[k] > 0 ) // 1x1 pivot block
    {
      const Int_t kp = ipiv[k] - 1 ;
      steps.push_back( { k , 1 , kp } ) ;
      D ( k , k ) = rawU ( k , k ) ;
      k -= 1 ;
    }
    else // 2x2 pivot block
    {
      const Int_t kp = -ipiv[k] - 1 ;
      steps.push_back( { k , 2 , kp } ) ;
      D ( k - 1 , k - 1 ) = rawU ( k - 1 , k - 1 ) ;
      D ( k - 1 , k     ) = rawU ( k - 1 , k     ) ;
      D ( k     , k - 1 ) = rawU ( k - 1 , k     ) ;
      D ( k     , k     ) = rawU ( k     , k     ) ;
      k -= 2 ;
    }
  }
  
  U.ResizeTo   ( n , n ) ;
  U.UnitMatrix () ;
  
  /// Construct U = P_{n-1} U_{n-1} ... P_0 U_0 by applying steps in forward order
  for ( auto it = steps.rbegin() ; it != steps.rend() ; ++it )
  {
    const BKStep& s = *it ;
    if ( s.type == 1 )
    {
      const Int_t sk  = s.k ;
      const Int_t skp = s.kp ;
      
      for ( Int_t i = 0 ; i < sk ; ++i )
      {
        const double mult = rawU ( i , sk ) ;
        for ( Int_t c = 0 ; c < n ; ++c )
        {
          U ( i , c ) += mult * U ( sk , c ) ;
        }
      }
      if ( skp != sk )
      {
        for ( Int_t c = 0 ; c < n ; ++c )
        {
          std::swap ( U ( sk , c ) , U ( skp , c ) ) ;
        }
      }
    }
    else // 2x2 pivot block
    {
      const Int_t sk  = s.k ;
      const Int_t skp = s.kp ;
      
      for ( Int_t i = 0 ; i < sk - 1 ; ++i )
      {
        const double mult_k1 = rawU ( i , sk - 1 ) ;
        const double mult_k  = rawU ( i , sk     ) ;
        for ( Int_t c = 0 ; c < n ; ++c )
        {
          U ( i , c ) += mult_k1 * U ( sk - 1 , c ) + mult_k * U ( sk , c ) ;
        }
      }
      if ( skp != sk - 1 )
      {
        for ( Int_t c = 0 ; c < n ; ++c )
        {
          std::swap ( U ( sk - 1 , c ) , U ( skp , c ) ) ;
        }
      }
    }
  }
  //
  return Ostap::StatusCode::SUCCESS ;
}
// ========================================================================


#include "Ostap/LinAlg.h"
#include "Ostap/LinAlgUtils.h"


// ======================================================================
/* Bunch-Kaufman decomposition of symmetric matrices 
 * @param[in]  A Input symmetric matrix to decompose.
 * @param[out] U triangular factor matrix U.
 * @param[out] D block-diagonal symmetric matrix D containing 1x1 and 2x2 blocks.
 * @return status code
 * @attention  here we use TDecompBK from ROOT 
 * @see  TDecompBK 
 */
// ======================================================================
Ostap::StatusCode Ostap::Math::GSL::BK
( const Ostap::Math::GSL::Matrix&  A ,  
  Ostap::Math::GSL::Matrix&        U , 
  Ostap::Math::GSL::Matrix&        D )
{
  //
  if ( A.nRows() != A.nCols() ) { return INVALID_GMATRIX ; }
  //
  const std::size_t N = A.nRows() ;
  const Int_t       n = static_cast<Int_t> ( N ) ;
  //
  TMatrixTSym<double>  a { n     } ;
  TMatrixT<double>     u { n , n } ;
  TMatrixTSym<double>  d { n     } ;
  //
  // copy GSL matrix into symmetric T-matrix 
  for ( Int_t i = 0 ; i < n ; ++i )
  {
    const double aii = A ( i , i ) ;
    a ( i , i ) = aii ;
    // read only low-triangular part 
    for ( Int_t j = 0 ; j < i ; ++j )
    {
      const double aij = A ( i , j ) ;
      a ( i , j ) = aij ;
      a ( j , i ) = aij ;      
    }
  }
  //
  const Ostap::StatusCode sc = Ostap::Math::BunchKaufman ( a , u , d ) ;
  if ( sc.isFailure() ) { return sc ; }
  //
  U = Ostap::Math::GSL::matrix ( u ) ;
  D = Ostap::Math::GSL::matrix ( d ) ;
  //
  return Ostap::StatusCode::SUCCESS ; 
}
// ======================================================================

// ============================================================================
//                                                                     The END 
// ============================================================================
