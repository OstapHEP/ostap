// ============================================================================
#ifndef OSTAP_EPDF_H 
#define OSTAP_EPDF_H 1
// ============================================================================
// Include files
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Power.h"
#include "Ostap/ECDF.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    // Kernels
    // https://en.wikipedia.org/wiki/Kernel_(statistics)
    // ========================================================================
    /// uniform Kernel 
    inline double k_uniform      ( const double u )
    { return 0.5 * ( std::abs ( u ) <= 1 ) ; }
    /// Triangular kernel
    inline double k_triangular   ( const double u )
    { return std::abs ( u ) <= 1 ? ( 1 - std::abs ( u ) ) : 0.0 ; }
    /// Epanechnikov/parabolic  kernel 
    inline double k_epanechnikov ( const double u )
    { return std::abs ( u ) <= 1 ? ( 1 - u * u )          : 0.0 ; }
    /// Epanechnikov Parabolic kernel
    inline double k_parabolic    ( const double u )
    { return k_epanechnikov ( u ) ; }
    /// Quartic/biweight kernel  
    inline double k_quartic      ( const double u )
    { return std::abs ( u ) <= 1 ? 15 * Ostap::Math::POW ( 1.0 - u * u ,  2  ) / 16 : 0.0 ; } 
    /// Quartic/biweight kernel
    inline double k_biweight     ( const double u )
    { return k_quartic ( u ) ; }
    /// Triweight kernel  
    inline double k_triweight    ( const double u )
    { return std::abs ( u ) <= 1 ? 35 * Ostap::Math::POW ( 1.0 - u * u ,  3  ) / 32 : 0.0 ; } 
    /// Tricube kernel
    inline double k_tricube      ( const double u )
    { return std::abs ( u ) <= 1 ? 70 * Ostap::Math::POW ( 1.0 - std::abs ( u * u * u ) , 3 ) / 81 : 0.0 ; } 
    /// Gaussian kernel 
    double        k_gaussian     ( const double u ) ; 
    /// Cosine kernel 
    double        k_cosine       ( const double u ) ; 
    /// Logistic Kernel
    double        k_logistic     ( const double u ) ;
    /// sigmoid kernel 
    double        k_sigmoid      ( const double u ) ; 
    // ====================================================================
    /** @class DensityEstimator
     *  Helper class for non-parametetric density estimators
     *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
     */ 
    class DensityEstimator
    {
      //====================================================================
    public: 
      //====================================================================
      enum Kernel
	{
	  Uniform      , 
	  Rectangular  = Uniform      , 
	  Boxcar       = Uniform      , 
	  Triangular   ,
	  Epanechnikov , 
	  Parabolic    = Epanechnikov , 
	  Quartic      , 
	  Biweight     = Quartic      , 
	  Triweight    , 
	  Tricube      ,
	  Gaussian     , 
	  Cosine       , 
	  Logistic     , 
	  Sigmoid      , 
	  Last         = Sigmoid    
	};
      // ======================================================================
    public:
      // ======================================================================
      /// get the kernel estimate
      static double kernel ( const double u , const Kernel k ) ;
      // ======================================================================
      /// get the "optimal" value for smoothing parameter 
      static double hopt   ( const  ECDF& data ) ;
      /// get the "optimal" value for smoothing parameter 
      static double hopt   ( const WECDF& data ) ;
      // ======================================================================
    };
    // ========================================================================
    
    // =======================================================================
    /** @class EPDF
     *  Helper utility to eatimate the PDF for emprical data using 
     *  Kernel estimators
     */
    class EPDF 
    {
    public:
      // =====================================================================
      /// create the emppirical PDF from empirical CDF 
      EPDF
      ( const ECDF&                                 ecdf  ,
	const Ostap::Math::DensityEstimator::Kernel k     ,
	const double                                h = 0 ) ;
      // =====================================================================
    public:
      // =====================================================================
      /// get the PDF 
      inline double operator () ( const double x ) const { return evaluate ( x ) ; }
      /// get the PDF 
      inline double pdf         ( const double x ) const { return evaluate ( x ) ; }
      /// get the PDF 
      double evaluate ( const double x ) const ;
      /// get the CDF 
      inline double cdf         ( const double x ) const { return m_cdf ( x ) ; }
      // =====================================================================
    public:
      // =====================================================================
      /// get ECDF 
      const  ECDF&   cdf () const { return m_cdf  ; }
      // =====================================================================
      /// get the kernel 
      inline Ostap::Math::DensityEstimator::Kernel kernel () const { return m_k ; }
      // =====================================================================
      /// get the smoothing parameter
      inline double h    () const { return m_h ; }
      // =====================================================================
    public:
      // =====================================================================
      /// update smoothing parameters  (
      bool setH      ( const double h ) ;
      /// set hew kernel
      bool setKernel ( const Ostap::Math::DensityEstimator::Kernel k ) ;
      // =====================================================================
    private:
      // =====================================================================
      ECDF                                  m_cdf ;
      Ostap::Math::DensityEstimator::Kernel m_k    { Ostap::Math::DensityEstimator::Epanechnikov } ;
      double                                m_h    { 0       } ; 
      // =====================================================================
    } ;
    // ========================================================================
    /** @class EPDF
     *  Helper utility to eatimate the PDF for emprical data using 
     *  Kernel estimators
     */
    class WEPDF 
    {
    public:
      // =====================================================================
      /// create the emppirical PDF from empirical CDF 
      WEPDF
      ( const WECDF&                                ecdf  ,
	const Ostap::Math::DensityEstimator::Kernel k     ,
	const double                                h = 0 ) ;
      // =====================================================================
    public:
      // =====================================================================
      /// get the PDF 
      inline double operator () ( const double x ) const { return evaluate ( x ) ; }
      /// get the PDF 
      inline double pdf         ( const double x ) const { return evaluate ( x ) ; }
      /// get the PDF 
      double evaluate ( const double x ) const ;
      /// get the CDF 
      inline double cdf         ( const double x ) const { return m_cdf ( x ) ; }
      // =====================================================================
    public:
      // =====================================================================
      /// get ECDF 
      const WECDF& cdf () const { return m_cdf  ; }
      // =====================================================================
      /// get the kernel 
      inline Ostap::Math::DensityEstimator::Kernel kernel () const { return m_k ; }
      // =====================================================================
      /// get the smoothing parameter
      inline double h  () const { return m_h ; }
      // =====================================================================
    public:
      // =====================================================================
      /// update smoothing parameters  (
      bool setH      ( const double h ) ;
      /// set hew kernel
      bool setKernel ( const Ostap::Math::DensityEstimator::Kernel k ) ;
      // =====================================================================
    private: 
      // =====================================================================
      WECDF                                 m_cdf ;
      Ostap::Math::DensityEstimator::Kernel m_k    { Ostap::Math::DensityEstimator::Epanechnikov } ;
      double                                m_h    { 0       } ; 
      // =====================================================================
    } ;    
    // =======================================================================    
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_EPDF_H
// ============================================================================
