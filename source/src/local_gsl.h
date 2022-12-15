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
  /// get GSL-workspace for CQUAD integration 
  inline gsl_integration_cquad_workspace* workspace_cquad 
  ( const Ostap::Math::WorkSpace& ws )
  {
    void* _ws =  ws.workspace_cquad () ;
    return (gsl_integration_cquad_workspace*) _ws ;
  }
  // ==========================================================================
  /// get GSL-workspace for Romberg integration 
  inline gsl_integration_romberg_workspace* workspace_romberg 
  ( const Ostap::Math::WorkSpace& ws )
  {
    void* _ws =  ws.workspace_romberg () ;
    return (gsl_integration_romberg_workspace*) _ws ;
  }
  // ==========================================================================
  /// get GSL-workspace for CQUAD integration 
  inline gsl_integration_cquad_workspace* cquad_workspace
  ( const Ostap::Math::WorkSpace& ws )
  { return workspace_cquad ( ws ) ; }
  // ==========================================================================
  /// get GSL-workspace for Romberg integration 
  inline gsl_integration_romberg_workspace* romberg_workspace
  ( const Ostap::Math::WorkSpace& ws )
  { return workspace_romberg ( ws ) ; }
  // ==========================================================================
  // get the default size of the main GSL-workspace 
  // ==========================================================================
  /** @var s_SIZE
   *  the default workspace size parameter for GSL-integration
   *  @see https://www.gnu.org/software/gsl/doc/html/integration.html
   *  quote: "The maximum number of subintervals is given by limit,
   *  which may not exceed the allocated size of the workspace."  
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const std::size_t s_SIZE = 10000 ;
  // ==========================================================================
  /** @var s_SIZE_CQUAD
   *  the default workspace size parameter for CQUAD double adaptive GSL-integration
   *  @see https://www.gnu.org/software/gsl/doc/html/integration.html
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2022-12-13
   */
  const std::size_t s_SIZE_CQUAD = 5000 ;
  // ==========================================================================
  /** @var s_SIZE_ROMBERG
   *  the default workspace size parameter for Romberg GSL-integration
   *  Number of divisions is \f$ 2^n + 1 \f$ 
   *  @see https://www.gnu.org/software/gsl/doc/html/integration.html
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2022-12-13
   */
  const std::size_t s_SIZE_ROMBERG  = 26 ;
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
  /** @var s_APRECISION_QAG
   *  the default absolute precision for QAG calculations,
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_APRECISION_QAG  = s_APRECISION ;
  // ==========================================================================
  /** @var s_RPRECISION_QAG
   *  the default relative precision for QAG calculations,
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_RPRECISION_QAG  = s_RPRECISION ;
  // ==========================================================================

  // ==========================================================================
  /** @var s_APRECISION_QAGI
   *  the default absolute precision for QAGI calculations,
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_APRECISION_QAGI  = s_APRECISION ;
  // ==========================================================================
  /** @var s_RPRECISION_QAGI
   *  the default relative precision for QAGI calculations,
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_RPRECISION_QAGI  = 1.e-7 ;
  // ==========================================================================

  // ==========================================================================
  /** @var s_APRECISION_QAGIU
   *  the default absolute precision for QAGIU calculations,
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_APRECISION_QAGIU  = s_APRECISION_QAGI ;
  // ==========================================================================
  /** @var s_RPRECISION_QAGIU
   *  the default relative precision for QAGIU calculations,
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_RPRECISION_QAGIU  = s_RPRECISION_QAGI ;
  // ==========================================================================
  
  // ==========================================================================
  /** @var s_APRECISION_QAGIL
   *  the default absolute precision for GAGIL calculations,
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_APRECISION_QAGIL  = s_APRECISION_QAGIU ;
  // ==========================================================================
  /** @var s_RPRECISION_QAGIL
   *  the default relative precision for GAGIL calculations,
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_RPRECISION_QAGIL  = s_RPRECISION_QAGIU ;
  // ==========================================================================

  // ==========================================================================
  /** @var s_APRECISION_QAGP
   *  the default absolute precision for QAGP calculations,
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_APRECISION_QAGP  = s_APRECISION ;
  // ==========================================================================
  /** @var s_RPRECISION_QAGP
   *  the default relative precision for QAGP calculations,
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_RPRECISION_QAGP  = 1.e-7 ;
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
  /** @var s_RPRECISION_QAWC
   *  the default relative QAWC precision for various calculations,
   *  in particular GSL integration
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_RPRECISION_QAWC  = 2.e-7 ;
  // ==========================================================================

  // ==========================================================================
  /** @var s_APRECISION_CQUAD
   *  the default absolute precision for CQUAD double adaptive integrator 
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_APRECISION_CQUAD  = s_APRECISION ;
  // ==========================================================================
  /** @var s_RPRECISION_CQUAD
   *  the default relative precision for CQUAD double adaptive integrator 
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_RPRECISION_CQUAD  = s_RPRECISION;
  // ==========================================================================

  // ==========================================================================
  /** @var s_APRECISION_ROMBERG
   *  the default absolute precision for Romberg integrator 
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_APRECISION_ROMBERG  = s_APRECISION ;
  // ==========================================================================
  /** @var s_RPRECISION_ROMBERG
   *  the default relative precision for Romberg integrator 
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_RPRECISION_ROMBERG  = s_RPRECISION;
  // ==========================================================================

  // ==========================================================================
  /** @var s_APRECISION_CUBE2D
   *  the default absolute precision for 2D-cubatures 
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_APRECISION_CUBE2D = s_APRECISION ;
  // ==========================================================================
  /** @var s_RPRECISION_CUBE2D
   *  the default relative precision for 2D-cubatures 
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_RPRECISION_CUBE2D = s_RPRECISION ;
  // ==========================================================================

  // ==========================================================================
  /** @var s_APRECISION_CUBE3D
   *  the default absolute precision for 3D-cubatures 
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_APRECISION_CUBE3D = s_APRECISION ;
  // ==========================================================================
  /** @var s_RPRECISION_CUBE3D
   *  the default relative precision for 3D-cubatures 
   *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
   *  @date 2010-05-23
   */
  const double s_RPRECISION_CUBE3D  = s_RPRECISION ;
  // ==========================================================================
}
// ============================================================================
//                                                                     The END 
// ============================================================================
#endif // LOCAL_GSL_H
// ============================================================================
