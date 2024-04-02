// ============================================================================
#ifndef OSTAP_COMBINE_H 
#define OSTAP_COMBINE_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL 
// ============================================================================
#include <array>
// ============================================================================
// ROOT
// ============================================================================
#include  "Math/SVector.h"
#include  "Math/SMatrix.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/SymmetricMatrixTypes.h"
#include "Ostap/MatrixUtils.h"
#include "Ostap/ValueWithError.h"
#include "Ostap/SVectorWithError.h"
#include "Ostap/StatusCode.h"
// ============================================================================
/** @file Ostap/Combine.h 
 *  Helper utility to combine   correlated measurements:
 *   "BLUE" : Best Linear Unbiased Estimator 
 * 
 *  @see P.Avery "Combining measurements with correlated errors", CBX 95 55
 *  @see http://www.phys.ufl.edu/~avery/fitting/error_correl.ps.gz
 *  @see http://www.researchgate.net.publication/2345482_Combining_Measurements_with_Correlated_Errors
 * 
 *  @see Louis Lyons, Duncan Gibaut, Peter Clifford, 
 *       "How to combine correlated estimates of a single physical quantity",
 *       Nuclear Instruments and Methods in Physics Research Section A:
 *        Accelerators, Spectrometers, Detectors and Associated Equipment
 *        Volume 270, Issue 1, 1 July 1988, Pages 110-117
 *  @see https://doi.org/10.1016/0168-9002(88)90018-6
 */  
// ============================================================================
namespace  Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    /** @class Combine
     *  Helper utility to combine   correlated measurements 
     *  "BLUE" : Best Linear Unbiased Estimator 
     *
     *  @see P.Avery "Combining measurements with correlated errors", CBX 95 55
     *  @see http://www.phys.ufl.edu/~avery/fitting/error_correl.ps.gz
     *  @see http://www.researchgate.net.publication/2345482_Combining_Measurements_with_Correlated_Errors
     *  @see Louis Lyons, Duncan Gibaut, Peter Clifford, 
     *       "How to combine correlated estimates of a single physical quantity",
     *       Nuclear Instruments and Methods in Physics Research Section A:
     *        Accelerators, Spectrometers, Detectors and Associated Equipment
     *        Volume 270, Issue 1, 1 July 1988, Pages 110-117
     *  @see https://doi.org/10.1016/0168-9002(88)90018-6
     *
     *  @author  Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-09-28
     */
    template <unsigned int D,
              class T=double,
              typename std::enable_if<(D>1),int>::type = 1 >
    class Combine
    {
    public:
      // ============================================================================
      typedef  ROOT::Math::SVector<T,D>                                 Data          ;
      typedef  ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >   Covariance    ;
      typedef  Ostap::Math::SVectorWithError<D,T>                       DataWithError ;
      // ============================================================================
    public:
      // =======================================================================
      // constructor from the vector of data and cov matrix 
      Combine
      ( const Data&       data , 
        const Covariance& cov2 ) 
        : m_data ( data ) 
        , m_cov2 ( cov2 )
        , m_vxi  ()  
        , m_w    () 
      {
        if ( Ostap::Math::inverse ( m_cov2 , m_vxi ) ) 
        {
          Ostap::throwException ( "Covariance matrix is not innvertile!" ,
                                  "Ostap::Math::Combine<>" , 730 ) ;
        }
        const Data& vone = this->units() ;
        m_w = ( m_vxi * vone ) / ROOT::Math::Similarity( m_vxi , vone ) ; 
      }
      // ======================================================================
      // constructor from the data
      Combine ( const DataWithError& data ) 
        : Combine ( data.value() , data.cov2() )
      {}
      // ======================================================================
      Combine
      ( const std::array<double,D>& data , 
        const Covariance&           cov2 ) 
        : Combine ( Data ( data.begin() , data.end() ), cov2 )
      {}
      // ======================================================================
      // constructor from the vector of data and cov matrices
      template <typename ... ARGS> 
      Combine
      ( const Data&       data , 
        const Covariance& cov1 , 
        const Covariance& cov2 ,
        ARGS...           args ) 
        : Combine ( data , cov1 + cov2 , args... )
      {}
      // ======================================================================
      // constructor from the data
      template <typename ... ARGS> 
      Combine
      ( const DataWithError& data ,
        ARGS...              args ) 
        : Combine ( data.value() , data.cov2() , args... )
      {}
      // ======================================================================
      // constructor from the data
      template <typename ... ARGS> 
      Combine
      ( const std::array<double,D>& data , 
        const Covariance&           cov2 ,
        ARGS...                     args )
        : Combine ( Data ( data.begin() , data.end() ), cov2 , args... )
      {}
      // ======================================================================
    public:
      // ======================================================================
      /// the main method:  get a combined value using the calculated weights
      Ostap::Math::ValueWithError result () const 
      { 
        const double r  = ROOT::Math::Dot        ( m_data , m_w ) ;
        const double e2 = ROOT::Math::Similarity ( m_cov2 , m_w ) ;
        return Ostap::Math::ValueWithError ( r , e2 ) ;
      }
      /// get the calculated weights
      const Data&       weights () const { return m_w    ; }
      /// get the data
      const Data&       data    () const { return m_data ; }
      /// get the covarinace 
      const Covariance& cov2    () const { return m_data ; }
      /// get the chi2
      double            chi2    () const
      {
        const double r  = ROOT::Math::Dot        ( m_data , m_w ) ;
        Data delta = m_data - r  ;
        return ROOT::Math::Similarity ( m_vxi , delta ) ;
      }
      // ======================================================================
    private:
      // ======================================================================
      /// get vector of units 
      const Data& units () const 
      {
        static Data s_units ;
        if ( 1 != s_units[0] ) { Ostap::Math::setToScalar ( s_units , T(1) ) ; }
        return s_units ;
      }
      // ======================================================================
    private:
      // ======================================================================
      /// input data vector 
      Data        m_data  {} ;            // input data vector 
      /// the overall covariance matrix 
      Covariance  m_cov2  {} ;            // the overall covariance matrix
      // ======================================================================
    private:
      // ======================================================================
      /// inverse  covariance matrix 
      Covariance  m_vxi   {} ;            // inverse covariance matrix
      /// weights 
      Data        m_w     {} ; // weights
      // ======================================================================
    } ;  
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
namespace Ostap
{
  // =========================================================================
  namespace  Math
  {
    // ========================================================================
    /** combine two measurements <code>x</code> and <code>y</code>
     *  with covarinace matrix <code>cov</code>
     *  @param x   (INPUT) the first  measurement 
     *  @param y   (INPUT) the second measurement 
     *  @param cov (INPUT) covariance matrix 
     *  @return combined result
     *  @author  Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-09-28
     */
    Ostap::Math::ValueWithError 
    combine  ( const double               x   , 
               const double               y   , 
               const Ostap::SymMatrix2x2& cov ) ;
    // ========================================================================
    /** combine two measurements <code>x1</code> and <code>x2</code>
     *  using correlation coefficient <code>rho</code>:  \f$-1\le\rho\le1\f$
     *  @param x1  (INPUT) the first  measurement 
     *  @param x2  (INPUT) the second measurement 
     *  @param rho (INPUT) correlation coefficient 
     *  @return combined result
     *  @author  Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-09-28
     */
    Ostap::Math::ValueWithError 
    combine 
    ( const Ostap::Math::ValueWithError& x1  ,
      const Ostap::Math::ValueWithError& x2  , 
      const double                       rho ) ;
    // ========================================================================
    /** combine two measurements <code>x1</code> and <code>x2</code>
     *  using their "statistical" uncertainties (assumed to be uncorrelated) 
     *  and a covariance matrix of "systematic" uncertainties
     *  @param x1   (INPUT) the first  measurement 
     *  @param x2   (INPUT) the second measurement 
     *  @param syst (INPUT) covariance matrix of systematic uncertainties  
     *  @return combined result
     *  @author  Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-09-28
     */
    Ostap::Math::ValueWithError 
    combine 
    ( const Ostap::Math::ValueWithError& x1   ,
      const Ostap::Math::ValueWithError& x2   , 
      const Ostap::SymMatrix2x2&         syst ) ;
    // ========================================================================
    /** combine three measurements:
     *  - <code>x1</code>, 
     *  - <code>x2</code> and 
     *  - <code>x3</code>
     *  using their "statistical" uncertainties (assumed to be uncorrelated) 
     *  and a covariance matrix of "systematic" uncertainties
     *  @param x1   (INPUT) the first  measurement 
     *  @param x2   (INPUT) the second measurement 
     *  @param x3   (INPUT) the third  measurement 
     *  @param syst (INPUT) covariance matrix of systematic uncertainties  
     *  @return combined result
     *  @author  Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-09-28
     */
    Ostap::Math::ValueWithError 
    combine 
    ( const Ostap::Math::ValueWithError& x1   ,
      const Ostap::Math::ValueWithError& x2   , 
      const Ostap::Math::ValueWithError& x3   , 
      const Ostap::SymMatrix3x3&         syst ) ;
    // ========================================================================
    /** combine four measurements:
     *  - <code>x1</code>, 
     *  - <code>x2</code>,
     *  - <code>x3</code> and 
     *  - <code>x4</code>
     *  using their "statistical" uncertainties (assumed to be uncorrelated) 
     *  and a covariance matrix of "systematic" uncertainties
     *  @param x1   (INPUT) the first  measurement 
     *  @param x2   (INPUT) the second measurement 
     *  @param x3   (INPUT) the third  measurement 
     *  @param x4   (INPUT) the fourth measurement 
     *  @param syst (INPUT) covariance matrix of systematic uncertainties  
     *  @return combined result
     *  @author  Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2015-09-28
     */
    Ostap::Math::ValueWithError 
    combine 
    ( const Ostap::Math::ValueWithError& x1   ,
      const Ostap::Math::ValueWithError& x2   , 
      const Ostap::Math::ValueWithError& x3   , 
      const Ostap::Math::ValueWithError& x4   , 
      const Ostap::SymMatrix4x4&         syst ) ;
    // ========================================================================
  }
  // ==========================================================================
}
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // LHCBMATH_COMBINE_H
// ============================================================================
