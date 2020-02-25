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
    template <unsigned int D, class T=double>
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
      Combine ( const Data&       data , 
                const Covariance& cov2 ) 
        : m_data ( data ) 
        , m_cov2 ( cov2 ) 
        , m_w    () 
      { m_w = this->getWeights() ; }
      // ======================================================================
      // constructor from the vector of data and cov matrices 
      Combine ( const Data&       data , 
                const Covariance& cov1 , 
                const Covariance& cov2 )
        : m_data ( data ) 
        , m_cov2 ( cov1 + cov2  ) 
        , m_w    () 
      { m_w = this->getWeights() ; }
      // ======================================================================
      // constructor from the vector of data and cov matrices 
      Combine ( const Data&       data , 
                const Covariance& cov1 , 
                const Covariance& cov2 ,
                const Covariance& cov3 )
        : m_data ( data ) 
        , m_cov2 ( cov1 + cov2 + cov3 ) 
        , m_w    () 
      { m_w = this->getWeights() ; }
      // ======================================================================
      // constructor from the vector of data and cov matrices 
      Combine ( const Data&       data , 
                const Covariance& cov1 , 
                const Covariance& cov2 ,
                const Covariance& cov3 ,
                const Covariance& cov4 )
        : m_data ( data ) 
        , m_cov2 ( cov1 + cov2 + cov3 + cov4 ) 
        , m_w    () 
      { m_w = this->getWeights() ; }
      // ======================================================================
      // constructor from the data
      Combine ( const DataWithError& data ) 
        : m_data ( data.value() ) 
        , m_cov2 ( data.cov2 () ) 
        , m_w    () 
      { m_w = this->getWeights() ; }
      // ======================================================================
      // constructor from the vector of data and cov matrix 
      Combine ( const std::array<double,D>& data , 
                const Covariance&           cov2 ) 
        : m_data ( data.begin() , data.end ()  ) 
        , m_cov2 ( cov2 ) 
        , m_w    () 
      { m_w = this->getWeights() ; }
      // ======================================================================
    public:
      // ======================================================================
      /// calculate weights
      Data getWeights() const 
      {
        // the inverse covariance matrix 
        Covariance vxi ;
        if ( Ostap::Math::inverse ( m_cov2 , vxi ) ) 
        { 
          /* try to regularize the matrix */
          const double tr    = std::abs ( Ostap::Math::trace           ( m_cov2 ) ) / D     ;
          const double mnd   = std::abs ( Ostap::Math::min_diagonal    ( m_cov2 ) )         ;
          const double mxd   = std::abs ( Ostap::Math::max_diagonal    ( m_cov2 ) ) / D     ;
          const double mnad  = std::abs ( Ostap::Math::minabs_diagonal ( m_cov2 ) )         ;
          const double mxad  = std::abs ( Ostap::Math::maxabs_diagonal ( m_cov2 ) ) / D     ;
          const double mn    = std::abs ( Ostap::Math::min_element     ( m_cov2 ) ) / D / D ;
          const double mx    = std::abs ( Ostap::Math::max_element     ( m_cov2 ) ) / D / D ;
          const double mna   = std::abs ( Ostap::Math::minabs_element  ( m_cov2 ) ) / D / D ;
          const double mxa   = std::abs ( Ostap::Math::maxabs_element  ( m_cov2 ) ) / D / D ;
          //
          double small = 1.e-4 * mxa ;
          //
          if ( 0 < tr   )  {  small = std::min ( small , 1.e-3 * tr   ) ; }
          if ( 0 < mnd  )  {  small = std::min ( small , 1.e-3 * mnd  ) ; }
          if ( 0 < mxd  )  {  small = std::min ( small , 1.e-3 * mxd  ) ; }
          if ( 0 < mnad )  {  small = std::min ( small , 1.e-3 * mnad ) ; }
          if ( 0 < mxad )  {  small = std::min ( small , 1.e-3 * mxad ) ; }
          if ( 0 < mn   )  {  small = std::min ( small , 1.e-3 * mn   ) ; }
          if ( 0 < mx   )  {  small = std::min ( small , 1.e-3 * mx   ) ; }
          if ( 0 < mna  )  {  small = std::min ( small , 1.e-3 * mna  ) ; }
          if ( 0 < mxa  )  {  small = std::min ( small , 1.e-3 * mxa  ) ; }
          //
          Covariance fixed ( m_cov2 ) ;
          for ( unsigned int i = 0 ; i < D ; ++i ) { fixed ( i, i ) += small ; }
          Ostap::Math::inverse ( fixed , vxi ) ;
        }
        //
        const Data& vone = this->units() ;
        return ( vxi * vone ) / ROOT::Math::Similarity( vxi , vone ) ;
      }
      /// the main method:  get a combined value using the calculated weights
      Ostap::Math::ValueWithError result () const 
      { 
        const double r  = ROOT::Math::Dot        ( m_data , m_w ) ;
        const double e2 = ROOT::Math::Similarity ( m_cov2 , m_w ) ;
        return Ostap::Math::ValueWithError ( r , e2 ) ;
      }
      /// get the calculated weights
      const Data& weights() const { return m_w ; }
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
      Data        m_data  ;            // input data vector 
      /// the overall covariance matrix 
      Covariance  m_cov2  ;            // the overall covariance matrix
      // ======================================================================
    private:
      // ======================================================================
      /// weights 
      Data m_w            ; // weights 
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
