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
  const std::size_t s_SIZE  = 5000 ;
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
  const double s_APRECISION  = 1.e-9 ;
  // ==========================================================================
  /** @var s_RPRECISION
   *  the default relative precision for various calculations,
   *  in particular GSL integration
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_RPRECISION  = 1.e-8 ;
  // ==========================================================================

  // ==========================================================================
  /** @var s_PRECISION_TAIL
   *  the low-relative precision for "tails" in GSL integration
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_PRECISION_TAIL = 1.e-7 ;
  // ===========================================================================

  // ==========================================================================
  /** @var s_APRECISION_TAIL
   *  the default absolute precision for various calculations,
   *  in particular GSL integration
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_APRECISION_TAIL  = s_APRECISION ;
  // ==========================================================================

  // ==========================================================================
  /** @var s_RPRECISION_TAIL
   *  the default relative precision for various calculations,
   *  in particular GSL integration
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_RPRECISION_TAIL  = s_PRECISION_TAIL ;
  // ==========================================================================

  // ===========================================================================
  /** @var s_APRECISION_GAQ
   *  the default absolute precision for GAQ calculations,
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_APRECISION_GAQ  = s_APRECISION ;
  // ==========================================================================
  /** @var s_RPRECISION_QAG
   *  the default relative precision for GAQ calculations,
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_RPRECISION_GAQ  = s_RPRECISION ;
  // ==========================================================================

  // ==========================================================================
  /** @var s_APRECISION_GAQI
   *  the default absolute precision for GAQI calculations,
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_APRECISION_GAQI  = s_APRECISION ;
  // ==========================================================================
  /** @var s_RPRECISION_QAGI
   *  the default relative precision for GAQI calculations,
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_RPRECISION_GAQI  = 1.e-7 ;
  // ==========================================================================

  // ==========================================================================
  /** @var s_APRECISION_GAQIU
   *  the default absolute precision for GAQIU calculations,
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_APRECISION_GAQIU  = s_APRECISION_GAQI ;
  // ==========================================================================
  /** @var s_RPRECISION_QAGIU
   *  the default relative precision for GAQIU calculations,
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_RPRECISION_GAQIU  = s_RPRECISION_GAQI ;
  // ==========================================================================
  
  // ==========================================================================
  /** @var s_APRECISION_GAQIL
   *  the default absolute precision for GAQIL calculations,
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_APRECISION_GAQIL  = s_APRECISION_GAQIU ;
  // ==========================================================================
  /** @var s_RPRECISION_QAGIL
   *  the default relative precision for GAQIL calculations,
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_RPRECISION_GAQIL  = s_RPRECISION_GAQIU ;
  // ==========================================================================

  // ==========================================================================
  /** @var s_APRECISION_GAQP
   *  the default absolute precision for GAQP calculations,
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_APRECISION_GAQP  = s_APRECISION ;
  // ==========================================================================
  /** @var s_RPRECISION_QAGP
   *  the default relative precision for GAQP calculations,
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_RPRECISION_GAQP  = 1.e-7 ;
  // ==========================================================================

  // ==========================================================================
  /** @var s_APRECISION_QAWC
   *  the default absolute QAWC precision for various calculations,
   *  in particular GSL integration
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_APRECISION_QAWC  = s_APRECISION ;
  // ==========================================================================
  /** @var s_EPRECISION_QAWC
   *  the default relative QAWC precision for various calculations,
   *  in particular GSL integration
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_RPRECISION_QAWC  = 2.e-7 ;
  // ==========================================================================
}
// ============================================================================
//                                                                     The END 
// ============================================================================
#endif // LOCAL_GSL_H
// ============================================================================
