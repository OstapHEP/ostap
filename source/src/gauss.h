// ===========================================================================
#ifndef OSTAP_GAUSS_H 
#define OSTAP_GAUSS_H 1
// ===========================================================================
/** @file gauss.h
 *  (local) collection of useful integrals 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 */
// ===========================================================================
namespace Ostap 
{
  // =========================================================================
  namespace Math
  {
    // =======================================================================
    namespace details 
    {
      // =====================================================================
      /** get the gaussian integral numerically
       *  \f[ f = \int_a^b \mathrm{e}^ {-\alpha x^2 + \beta x } \mathrm{d}x \f]
       *  @attention note the sign for alpha-term!
       *  @param alpha the alpha parameter
       *  @param beta  the beta  parameter
       *  @param a     the low  integration limit
       *  @param b     the high integration limit
       *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
       *  @date 2010-05-23
       */
      double gaussian_int_num
      ( const double alpha ,
        const double beta  ,
        const double a     ,
        const double b     ) ;
      // ======================================================================
      /** get the gaussian integral:
       *  \f[ f = \int_a^{\inf} \exp { -\alpha x^2 + \beta x } \mathrm{d}x \f]
       *  @attention note the sign for alpha-term!
       *  @param alpha the alpha parameter
       *  @param beta  the beta  parameter
       *  @param a     the low integration limit
       *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
       *  @date 2010-05-23
       */
      double gaussian_int_R ( const double alpha ,
                              const double beta  ,
                              const double a     ) ; 
      // ======================================================================
      /** get the gaussian integral:
       *  \f[ f = \int_{-\inf}^{b} \exp { -\alpha x^2 + \beta x } \mathrm{d}x \f]
       *  @attention note the sign for alpha-term!
       *  @param alpha the alpha parameter
       *  @param beta  the beta  parameter
       *  @param b     the high integration limit
       *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
       *  @date 2010-05-23
       */
      inline double gaussian_int_L ( const double alpha ,
                                     const double beta  ,
                                     const double b     ) 
      { return gaussian_int_R ( alpha , -beta , -b ) ; }        
      // ======================================================================
      /** get the gaussian integral:
       *  \f[ f = \int_a^b \exp { -\alpha x^2 + \beta x } \mathrm{d}x \f]
       *  @attention note the sign for alpha-term!
       *  @param alpha the alpha parameter
       *  @param beta  the beta  parameter
       *  @param a     the low  integration limit
       *  @param b     the high integration limit
       *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
       *  @date 2010-05-23
       */
      double gaussian_int
      ( const double alpha ,
        const double beta  ,
        const double a     ,
        const double b     ) ;
      // ======================================================================
      /** get the exponential integral 
       *  \f[ f = \int_a^b \exp { \beta x } \mathrm{d}x \f]
       *  @param beta  the beta  parameter
       *  @param a     the low  integration limit
       *  @param b     the high integration limit
       */
      double exponent_int
      ( const double beta , 
        const double a    , 
        const double b    ) ;
      // ======================================================================
    } //                              The end of namespace Ostap::Math::details 
    // ========================================================================
  } //                                         The end of namespace Ostap::Math 
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // GAUSS_H
// ============================================================================
