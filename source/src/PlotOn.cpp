// ============================================================================
// Include files 
// ============================================================================
// ROOT 
// ============================================================================
#include "RooCmdArg.h"
#include "RooPlot.h"
#include "RooAbsData.h"
#include "RooAbsReal.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/PlotOn.h"
// ============================================================================
/** @file
 *  Implementation file for functions from file Ostap/PlotOn.h
 *  @date 2021-06-18 
 *  @author Vanya Belyaev  Ivan.Belyaev@itep.ru
 */
// ============================================================================
RooPlot* Ostap::Utils::plotOn	(	const RooAbsReal* value , 
                                RooPlot*  	      frame ,
                                const RooCmdArg&  arg1  ,
                                const RooCmdArg&  arg2  ,
                                const RooCmdArg&  arg3  , 
                                const RooCmdArg&  arg4  ,
                                const RooCmdArg&  arg5  ,
                                const RooCmdArg&  arg6  ,
                                const RooCmdArg&  arg7  ,
                                const RooCmdArg&  arg8  )
{
  return  
    ( nullptr == value ) ? frame :
    ( nullptr == frame ) ? frame :
    value -> plotOn ( frame , 
                      arg1 , arg2 , arg3 , arg4 ,
                      arg5 , arg6 , arg7 , arg8 ) ;                
}
// =====================================================================
RooPlot* Ostap::Utils::plotOn	(	const RooAbsData* data  , 
                                RooPlot*  	      frame ,
                                const RooCmdArg&  arg1  ,
                                const RooCmdArg&  arg2  ,
                                const RooCmdArg&  arg3  , 
                                const RooCmdArg&  arg4  ,
                                const RooCmdArg&  arg5  ,
                                const RooCmdArg&  arg6  ,
                                const RooCmdArg&  arg7  ,
                                const RooCmdArg&  arg8  ) 
{
  return  
    ( nullptr == data  ) ? frame :
    ( nullptr == frame ) ? frame :
    data -> plotOn ( frame , 
                     arg1 , arg2 , arg3 , arg4 ,
                     arg5 , arg6 , arg7 , arg8 ) ;                
}
// ============================================================================
//                                                                      The END 
// ============================================================================
