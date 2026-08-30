// ============================================================================
#ifndef OSTAP_MATRIXFUN_H
#define OSTAP_MATRIXFUN_H 1
// ============================================================================
// Include files
// ============================================================================
// STD& STL 
// ============================================================================
#include <algorithm>
#include <cmath>
// ============================================================================
// ROOT
// ============================================================================
#include "Math/SMatrix.h"
#include "Math/SVector.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Power.h"
#include "Ostap/EigenSystem.h"
#include "Ostap/StatusCode.h"
// ============================================================================
/** @file Ostap/MatrixFun.h
 *  Useful matrix functions
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 */
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    /** function from vector: apply arbitrary function/functor to the vector
     *  @code
     *  ROOT::Math::SVector<...> input  = ... ;
     *  ROOT::Math::SVector<...> result = apply ( input , []( double v ) -> double { return std::sin ( v ) ; } ) ;     
     *  @endcode
     */
    template <class T, unsigned int D, typename F>
    inline ROOT::Math::SVector<T,D>
    apply 
    ( const ROOT::Math::SVector<T,D>& vct , F&& func ) 
    {
      ROOT::Math::SVector<T,D> result {} ;
      std::transform ( vct.begin() , vct.end() , result.begin() , std::forward<F>( func ) ) ;
      return result ;
    }
    // ========================================================================

    // ========================================================================
    /** Generic function of symmetric matrix
     *  @code
     *  MATRIX input  = ... ;
     *  MATRIX result {} ;
     *  auto func = [] ( double v ) -> double { return std::sin ( v ) ; } ;
     *  if ( apply ( input , result , func ) ) {  OK ; } 
     *  @endcode 
     */
    template <class T, unsigned int D, typename R, typename F>
    inline 
    bool apply 
    ( const ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >& mtrx   ,
      ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >&       result , F&& func)
    {
      // 
      if constexpr ( 1 == D )
      {
        const T w { mtrx ( 0 , 0 ) } ;
        result ( 0 , 0 ) = static_cast<T> ( func ( w ) ) ;
        return std::isfinite ( result ( 0 , 0 ) ) ;
      }
      //
      /// self-type 
      typedef typename ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > MTRX ;
      /// vector of eigenvalues
      typedef typename ROOT::Math::SVector<T,D>                             VCT  ;
      /// matrix of eigenvectors
      typedef typename ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > VCTR ;
      //
      VCT  V {} ;
      VCTR M {} ;
      /// get the eigen-values and matrix of eigen-vectors 
      const Ostap::Math::GSL::EigenSystem eigen {} ;
      if ( !eigen.eigenVectors ( mtrx , V , M ).isSuccess() ) { return false ; }
      //
      for ( unsigned int i = 0 ; i < D ; ++i )
      {
        const T vv = V [ i ] ;
        const T rr = static_cast<T> ( func ( vv ) ) ;
        if ( !std::isfinite ( rr ) ) { return false ; }
        V [ i ] = rr ; 
      }
      //
      result = ROOT::Math::SMatrixIdentity();
      for ( unsigned int i = 0 ; i < D ; ++i ) { result ( i , i ) = V ( i ) ; }
      //
      result = ROOT::Math::Similarity ( result , M ) ;
      //
      return true ;
    }     
    // ========================================================================
    
    // ========================================================================
    /// pow-function for square matrices 
    // ========================================================================

    // ========================================================================
    /** pow-function for square matrices with unsigned int compile-time power 
     *  @code
     *  MATRIX mtrx   = ...
     *  MATRIX result = pow<5> ( mtrx ) ;
     *  @endcode
     */
    template <class T, unsigned int D, typename R, unsigned int N>
    inline ROOT::Math::SMatrix<T,D,D,R>
    pow
    ( const ROOT::Math::SMatrix<T,D,D,R>& mtrx ) 
    {
      /// trivial cases
      if      constexpr ( 0 == N ) { return ROOT::Math::SMatrix<T,D,D,R> ( ROOT::Math::SMatrixIdentity () ) ; }
      else if constexpr ( 1 == N ) { return mtrx        ; }
      else if constexpr ( 2 == N ) { return mtrx * mtrx ; }
      else if constexpr ( 3 == N ) { return mtrx * mtrx * mtrx ; }
      //      
      if constexpr ( 1 == D )
      {
        const T w { mtrx ( 0 , 0 ) } ;
        ROOT::Math::SMatrix<T,D,D,R> result {} ;
        result ( 0 , 0 ) = static_cast<T> ( std::pow ( w , N ) ) ;
        return result ;
      }
      //
      return Ostap::Math::POW  ( mtrx , N ) ; 
    }
    // ========================================================================

    // ========================================================================
    /** pow-function for square matrices with unsigned int power (runtime)
     *  @code
     *  MATRIX  mtrx  = ...
     *  MATRIX result = pow ( mtrx , 5 ) ;
     *  @endcode
     */
    template <class T, unsigned int D, typename R>
    inline ROOT::Math::SMatrix<T,D,D,R>
    pow
    ( const ROOT::Math::SMatrix<T,D,D,R>& mtrx ,
      const unsigned int                  n    ) 
    {
      /// trivial cases
      if      ( 0 == n ) { return ROOT::Math::SMatrix<T,D,D,R> ( ROOT::Math::SMatrixIdentity () ) ; }
      else if ( 1 == n ) { return        mtrx ; }
      else if ( 2 == n ) { return mtrx * mtrx ; }
      else if ( 3 == n ) { return mtrx * mtrx * mtrx ; } 
      ///
      if constexpr ( 1 == D )
      {
        const T w { mtrx ( 0 , 0 ) } ;
        ROOT::Math::SMatrix<T,D,D,R> result {} ;
        result ( 0 , 0 ) = static_cast<T> ( std::pow ( w , n ) ) ; 
        return result ;
      }
      //
      return Ostap::Math::POW ( mtrx , n ) ; 
    }
    // ========================================================================

    // ========================================================================
    /** pow-function for symmetric square matrices and arbitrary integer argument 
     *  @code
     *  MATRIX  mtrx   = ...
     *  MTATRIX result {} ; 
     *  if ( pow ( mtrx , +5 , result ) ) {  OK ; } 
     *  if ( pow ( mtrx , -5 , result ) ) {  OK ; } 
     *  @endcode
     */
    template <class T, unsigned int D, typename R >
    inline 
    bool pow
    ( const ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >& mtrx   ,
      const int                                                     n      ,
      ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >&       result )
    {
      /// result type 
      typedef typename ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > MTRX ;
      /// trivial cases
      if      ( 0 == n ) { result = MTRX ( ROOT::Math::SMatrixIdentity () ) ; return true ; }
      else if ( 1 == n ) { result =               mtrx                      ; return true ; } 
      else if ( 2 == n ) { result =        mtrx * mtrx                      ; return true ; } 
      else if ( 3 == n ) { result = mtrx * mtrx * mtrx                      ; return true ; } 
      /// another trivial case 
      if constexpr ( 1 == D )
      {
        const T w { mtrx ( 0 , 0 ) } ;
        if ( w == 0 && n < 0 ) { return false ; }
        result ( 0 , 0 ) = static_cast<T> ( std::pow ( w , n ) ) ;
        return std::isfinite ( result ( 0 , 0 ) ) ;
      }
      //
      if ( n < 0 )
      {
        /// here matrix needs to be inverted
        int  ifail { 0 };
        /// 1) Try fast cholesky inversion 
        MTRX m { mtrx.InverseChol ( ifail ) } ; 
        /// 2) Fall to regular inversion 
        if ( 0 != ifail ) { m = mtrx.Inverse ( ifail ) ; }
        /// 3) matrix cannot be inverted 
        if ( 0 != ifail ) { return false ; }
        //
        return pow ( m , std::abs ( n ) , result ) ; 
      }
      //
      result = Ostap::Math::POW  ( mtrx , n ) ;
      return true ; 
    }     
    // ========================================================================

    // ========================================================================
    /** pow-function for symmetric square matrices with real power 
     *  @code
     *  MATRIX mtrx   = ...
     *  MATRIX result {} ; 
     *  if ( pow ( mtrx , +4.1 , result ) ) {  OK ; } 
     *  if ( pow ( mtrx , -4.1 , result ) ) {  OK ; } 
     *  @endcode
     */
    template <class T, unsigned int D, typename R>
    inline 
    bool pow 
    ( const ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >& mtrx   ,
      const double                                                  p      , 
      ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >&       result ) 
    {
      //
      if ( Ostap::Math::isshort ( p ) )
      {
        const int n = static_cast<short> ( Ostap::Math::round ( p ) ) ;
        return pow ( mtrx , n , result ) ; 
      }
      //
      auto func = [p]( T x ) -> T { return std::pow ( x , p ) ; } ;
      return apply ( mtrx , result , func ) ;
    }
    // ========================================================================

    // ========================================================================
    /** square root for symmetric positively semidefinite matrices 
     *  @code
     *  MATRIX mtrx   = ...
     *  MATRIX result {} ;  
     *  if ( sqrt ( mtrx , result ) ) {  OK ; } 
     *  @endcode
     */
    template <class T, unsigned int D, typename R>
    inline 
    bool sqrt
    ( const ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >& mtrx   ,
      ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >&       result ) 
    {
      //
      auto func = []( T x ) -> T { return std::sqrt ( x ) ; } ;
      return apply ( mtrx , result , func ) ;
    }
    // ========================================================================
    
    // ========================================================================
    /** matrix exponent 
     *  @code
     *  MATRIX mtrx   = ...
     *  MATRIX result {} ;  
     *  if ( exp ( mtrx , result ) ) {  OK ; } 
     *  @endcode
     */
    template <class T, unsigned int D, typename R>
    inline 
    bool exp 
    ( const ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >& mtrx   ,
      ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >&       result ) 
    {
      //
      auto func = []( T x ) -> T { return std::exp ( x ) ; } ;
      return apply ( mtrx , result , func ) ;
    }
    // ========================================================================
    
    // ========================================================================
    /** matrix logarithm (for symmetric positive deinite matrix)
     *  @code
     *  MATRIX mtrx   = ...
     *  MATRIX result {} ;  
     *  if ( log ( mtrx , result ) ) {  OK ; } 
     *  @endcode
     */
    template <class T, unsigned int D, typename R>
    inline 
    bool log 
    ( const ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >& mtrx   ,
      ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >&       result ) 
    {
      //
      auto func = []( T x ) -> T { return std::log ( x ) ; } ;
      return apply ( mtrx , result , func ) ;
    }
    // ========================================================================
    
    // ========================================================================    
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
#endif // OSTAP_MATRIXFUN_H
// ============================================================================
//                                                                      The END
// ============================================================================

