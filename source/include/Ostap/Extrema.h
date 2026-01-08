// ============================================================================
#ifndef OSTAP_EXTREMA_H
#define OSTAP_EXTREMA_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <functional>
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================    
    /** Search for the minimum of function at [a,b] using 
     *  golden-section-search rule  
     *  @aparam f the function 
     *  @param low    lower interval edge 
     *  @param high   upper interval edge
     *  @param guess  initial guess 
     *  @param abseps the required absolute precision   
     *  @param releps the required relative precision   
     *  @author Vanya BELYAEV  Ivan.Belyaev@cern.ch
     */
    double  minimum_golden_section
    ( std::function<double(double)> f      ,
      const double                  low    ,
      const double                  high   ,
      const double                  guess  ,
      const double                  abseps ,
      const double                  releps ) ;
    // ========================================================================
    /** Search for the maximum of function at [a,b] using 
     *  golden-section-search rule  
     *  @aparam f the function 
     *  @param low    lower interval edge 
     *  @param high   upper interval edge
     *  @param guess  initial guess 
     *  @param abseps the required absolute precision   
     *  @param releps the required relative precision   
     *  @author Vanya BELYAEV  Ivan.Belyaev@cern.ch
     */
    double maximum_golden_section
    ( std::function<double(double)> f      ,
      const double                  low    ,
      const double                  high   ,
      const double                  guess  ,
      const double                  abseps ,
      const double                  releps ) ;
    // ========================================================================
    /** Search for the minimum of function at [a,b] using 
     *  Brent' rule  
     *  @aparam f the function 
     *  @param low    lower interval edge 
     *  @param high   upper interval edge
     *  @param guess  initial guess 
     *  @param abseps the required absolute precision   
     *  @param releps the required relative precision   
     *  @author Vanya BELYAEV  Ivan.Belyaev@cern.ch
     */
    double  minimum_brent 
    ( std::function<double(double)> f      ,
      const double                  low    ,
      const double                  high   ,
      const double                  guess  ,
      const double                  abseps ,
      const double                  releps ) ;
    // ========================================================================
    /** Search for the maximum of function at [a,b] using 
     *  Brent' rule  
     *  @aparam f the function 
     *  @param low    lower interval edge 
     *  @param high   upper interval edge
     *  @param guess  initial guess 
     *  @param abseps the required absolute precision   
     *  @param releps the required relative precision   
     *  @author Vanya BELYAEV  Ivan.Belyaev@cern.ch
     */
    double maximum_brent 
    ( std::function<double(double)> f      ,
      const double                  low    ,
      const double                  high   ,
      const double                  guess  ,
      const double                  abseps ,
      const double                  releps ) ;
    // ========================================================================
    /** Search for the minimum of function at [a,b] using 
     *  variant of Brent’s algorithm which uses 
     *  the safe-guarded step-length algorithm of Gill and Murray.
     *  @aparam f the function 
     *  @param low    lower interval edge 
     *  @param high   upper interval edge
     *  @param guess  initial guess 
     *  @param abseps the required absolute precision   
     *  @param releps the required relative precision   
     *  @author Vanya BELYAEV  Ivan.Belyaev@cern.ch
     */
    double  minimum_quad_golden
    ( std::function<double(double)> f      ,
      const double                  low    ,
      const double                  high   ,
      const double                  guess  ,
      const double                  abseps ,
      const double                  releps ) ;
    // ========================================================================
    /** Search for the maximum of function at [a,b] using 
     *  variant of Brent’s algorithm which uses 
     *  the safe-guarded step-length algorithm of Gill and Murray.
     *  @aparam f the function 
     *  @param low    lower interval edge 
     *  @param high   upper interval edge
     *  @param guess  initial guess 
     *  @param abseps the required absolute precision   
     *  @param releps the required relative precision   
     *  @author Vanya BELYAEV  Ivan.Belyaev@cern.ch
     */
    double maximum_quad_golden
    ( std::function<double(double)> f      ,
      const double                  low    ,
      const double                  high   ,
      const double                  guess  ,
      const double                  abseps ,
      const double                  releps ) ;
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================   
} //                                                 The end of namespace Ostap
// ============================================================================
#endif // OSTAP_EXTREMA_H
// ============================================================================ 
