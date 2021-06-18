// ============================================================================
#ifndef OSTAP_PLOTON_H 
#define OSTAP_PLOTON_H 1
// ============================================================================
// Include files
// ============================================================================
// ROOT 
// ============================================================================
#include "RooCmdArg.h"
// ============================================================================
// Forward declarations 
// ============================================================================
class RooPlot    ;
class RooAbsReal ;
class RooAbsData ;
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Utils 
  {
    // ========================================================================
    RooPlot* plotOn	(	const RooAbsReal* value                      , 
                      RooPlot*  	      frame                      ,
                      const RooCmdArg&  arg1  = RooCmdArg::none () ,
                      const RooCmdArg&  arg2  = RooCmdArg::none () ,
                      const RooCmdArg&  arg3  = RooCmdArg::none () , 
                      const RooCmdArg&  arg4  = RooCmdArg::none () ,
                      const RooCmdArg&  arg5  = RooCmdArg::none () ,
                      const RooCmdArg&  arg6  = RooCmdArg::none () ,
                      const RooCmdArg&  arg7  = RooCmdArg::none () ,
                      const RooCmdArg&  arg8  = RooCmdArg::none () ) ;
    // ========================================================================
    RooPlot* plotOn	(	const RooAbsData* data                       , 
                      RooPlot* 	        frame                      ,
                      const RooCmdArg&  arg1  = RooCmdArg::none () ,
                      const RooCmdArg&  arg2  = RooCmdArg::none () ,
                      const RooCmdArg&  arg3  = RooCmdArg::none () ,
                      const RooCmdArg&  arg4  = RooCmdArg::none () ,
                      const RooCmdArg&  arg5  = RooCmdArg::none () ,
                      const RooCmdArg&  arg6  = RooCmdArg::none () ,
                      const RooCmdArg&  arg7  = RooCmdArg::none () ,
                      const RooCmdArg&  arg8  = RooCmdArg::none () ) ;
    // ========================================================================
  } 
  // ==========================================================================
}
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_PLOTON_
// ============================================================================
