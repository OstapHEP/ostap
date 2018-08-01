// ============================================================================
#ifndef OSTAP_CHI2SOLUTION_H 
#define OSTAP_CHI2SOLUTION_H 1
// ============================================================================
// Include files
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatusCode.h"
#include "Ostap/SVectorWithError.h"
// ============================================================================
// ROOT 
// ============================================================================
#include "Math/SMatrix.h"
#include "Math/SVector.h"
// ============================================================================
/** @file Ostap/Chi2Solution.h
 *  Generic solution for N-dimensional chi2-problem with R-constraints 
 */ 
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    /** @class Chi2Solution Ostap/Chi2Solution.h
     *  Generic solution for N-dimensional chi2-problem with R-constraints 
     *
     *  All Formulae and notation from Paul Avery :
     *  - Applied Fitting theory I: General Least Squares Theory
     *  - CBX 92-72, October 18, 1991 
     *
     *  @author Vanya Belyaev
     *  @date   2012-04-27
     */
    template <unsigned int N, unsigned int R, class T = double>
    class Chi2Solution 
    {
    public:
      // ======================================================================
      typedef ROOT::Math::SVector<T,N>                               DATA    ;
      typedef ROOT::Math::SMatrix<T,N,N,ROOT::Math::MatRepSym<T,N> > COV2    ;
      typedef ROOT::Math::SMatrix<T,R,N>                             CMTRX1  ;
      typedef ROOT::Math::SVector<T,R>                               COFF    ;
      typedef Ostap::Math::SVectorWithError<N,T>                     VECT    ;
      // ======================================================================
      // backup solution 
      typedef std::vector<DATA>                                      CMTRX2  ;
      // ======================================================================
    public:
      // ======================================================================
      /** make N-dimensional chi2-solution with R-constraints
       *  @param data (UPDATE) input approximation for the data vector 
       *  @param cov2 (UPDATE) the covariance matrix for input data 
       *  @param D    (UPDATE) the matrix   of constraints
       *  @param d    (UPDATE) the offsets for constraints
       *  @param chi2 (UPDATE) the chi2 
       *  @return status code 
       */
      static Ostap::StatusCode solve 
      ( DATA&         data , 
        COV2&         cov2 , 
        const CMTRX1& D    , 
        const COFF &  d    , 
        double&       chi2 ) 
      {
        //
        typedef ROOT::Math::SMatrix<T,R,R,ROOT::Math::MatRepSym<T,R> > MRxR ;
        typedef ROOT::Math::SMatrix<T,N,R>                             MNxR ;
        //
        MRxR vD = ROOT::Math::Similarity ( D , cov2 ) ;
        if ( !vD.Invert() ) 
        {
          chi2 = 1.e+100 ;
          return Ostap::StatusCode::FAILURE ;                             // RETURN
        } 
        //
        const ROOT::Math::SVector<T,R> alpha  =  D * data + d ;
        const ROOT::Math::SVector<T,R> lambda = vD * alpha    ;
        //
        const MNxR vTimesDt = cov2 * ROOT::Math::Transpose ( D ) ;
        //
        // make the solution 
        data -= vTimesDt * lambda ;
        cov2 -= ROOT::Math::Similarity ( vTimesDt , vD ) ;
        //
        // get chi2 
        chi2  = ROOT::Math::Dot        ( alpha  , lambda ) ;
        // chi2  = ROOT::Math::Similarity ( , vD ) ;
        //
        return StatusCode::SUCCESS ;                            // RETURN 
      }
      // ==================================================================
      /** make N-dimensional chi2-solution with R-constraints
       *  @param data (UPDATE) input approximation for the data vector 
       *  @param cov2 (UPDATE) the covariance matrix for input data 
       *  @param D2   (UPDATE) the matrix   of constraints
       *  @param d    (UPDATE) the offsets for constraints
       *  @param chi2 (UPDATE) the chi2 
       *  @return status code 
       */
      static Ostap::StatusCode solve 
      ( DATA&         data , 
        COV2&         cov2 , 
        const CMTRX2& D2   , 
        const COFF  & d    , 
        double&       chi2 ) 
      {
        //
        CMTRX1 D ;
        for ( unsigned int i = 0 ; i < R ; ++i )
        {
          for ( unsigned int j = 0 ; j < N ; ++j ) 
          {
            if ( i < D2.size() ) { D( i , j ) = D2[i][j] ; }
          }
        }
        //
        return solve ( data , cov2 , D , d, chi2 ) ; // REUTRN 
      }
      // ======================================================================
      /** make N-dimensional chi2-solution with R-constraints
       *  @param data (UPDATE) for the data&covarinace  
       *  @param D    (UPDATE) the matrix   of constraints
       *  @param d    (UPDATE) the offsets for constraints
       *  @param chi2 (UPDATE) the chi2 
       *  @return status code 
       */
      static Ostap::StatusCode solve
      ( VECT&         data , 
        const CMTRX1& D    , 
        const COFF&   d    , 
        double&       chi2 ) 
      {
        return solve ( data.value () , data.cov2  () , D , d , chi2 ) ;
      }
      // ======================================================================
      /** make N-dimensional chi2-solution with R-constraints
       *  @param data (UPDATE) for the data&covarinace  
       *  @param D    (UPDATE) the matrix   of constraints
       *  @param d    (UPDATE) the offsets for constraints
       *  @param chi2 (UPDATE) the chi2 
       *  @return status code 
       */
      static Ostap::StatusCode solve
      ( VECT&         data , 
        const CMTRX2& D    , 
        const COFF&   d    , 
        double&       chi2 ) 
      {
        return solve ( data.value () , data.cov2  () , D , d , chi2 ) ;
      }
      // ======================================================================
    }; 
    // ========================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_CHI2SOLUTION_H
// ============================================================================
