// ============================================================================
#ifndef LOCAL_GSL_H 
#define LOCAL_GSL_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&   STL
// ============================================================================
#include <limits>
// ============================================================================
// GSL
// ============================================================================
#include "gsl/gsl_errno.h"
#include "gsl/gsl_integration.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Workspace.h"
// ============================================================================
// Local
// ============================================================================
#include "GSL_sentry.h"
// ============================================================================
namespace 
{
  // ==========================================================================  
  // GSL
  // ==========================================================================
  typedef Ostap::Math::GSL::GSL_Error_Handler Sentry ;
  // ==========================================================================
  /// get GSL-workspace
  inline gsl_integration_workspace* workspace
  ( const Ostap::Math::WorkSpace& ws )
  {
    void* _ws =  ws.workspace() ;
    return (gsl_integration_workspace*) _ws ;
  }
  // ==========================================================================
  // get size of GSL-workspace 
  // ==========================================================================
  /** @var s_SIZE
   *  the workspace size parameter for GSL-integration
   *  @see https://www.gnu.org/software/gsl/doc/html/integration.html
   *  quote: "The maximum number of subintervals is given by limit,
   *  which may not exceed the allocated size of the workspace."  
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const std::size_t s_SIZE  = 2000 ;
  // ==========================================================================
  /** @var s_PRECISION
   *  the default precision for various calculations,
   *  in particular GSL integration
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_PRECISION  = 1.e-8 ;
  // ==========================================================================
  /** @var s_APRECISION
   *  the default absolute precision for various calculations,
   *  in particular GSL integration
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_APRECISION  = s_PRECISION ;
  // ==========================================================================
  /** @var s_RPRECISION
   *  the default relative precision for various calculations,
   *  in particular GSL integration
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_RPRECISION  = s_PRECISION ;
  // ==========================================================================
  /** @var s_PRECISION_TAIL
   *  the low-relative precision for "tails" in GSL integration
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_PRECISION_TAIL = 1.e-5 ;
  // ===========================================================================
  /** @var s_PRECISION_QAWC
   *  the default QAWC precision for various calculations,
   *  in particular GSL integration
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_PRECISION_QAWC  = 1.e-7 ;
  // ==========================================================================
  /** @var s_APRECISION_QAWC
   *  the default absolute QAWC precision for various calculations,
   *  in particular GSL integration
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_APRECISION_QAWC  = s_PRECISION_QAWC ;
  // ==========================================================================
  /** @var s_EPRECISION_QAWC
   *  the default relative QAWC precision for various calculations,
   *  in particular GSL integration
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_RPRECISION_QAWC  = s_PRECISION_QAWC ;
  // ==========================================================================
}
// ============================================================================
//                                                                     The END 
// ============================================================================
#endif // LOCAL_GSL_H
// ============================================================================
