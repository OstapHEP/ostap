// $Id$
// ==================================================================
// Include files 
// ==================================================================
// ROOT 
// ==================================================================
#include "RooAbsData.h"
#include "RooCmdArg.h"
// ==================================================================
// Ostap
// ==================================================================
#include "Ostap/HistoMake.h"
// ==================================================================
/** @file
 *  Implementation file for class Analysis::RooMakeHistos
 *  @see Analysis::RooMakeHistos 
 *  @see RooAbsData::createHistos 
 *  @date   2011-07-16 
 *  @author Vanya Belyaev
 */
// ==================================================================
// the only one method 
// ==================================================================
TH1* Ostap::HistoMake::create_histo 
( const RooAbsData&       dataset , 
  const std::string&      name    , 
  const RooAbsRealLValue& xvar    ,
  const RooCmdArg&        arg1    , 
  const RooCmdArg&        arg2    , 
  const RooCmdArg&        arg3    , 
  const RooCmdArg&        arg4    , 
  const RooCmdArg&        arg5    , 
  const RooCmdArg&        arg6    , 
  const RooCmdArg&        arg7    , 
  const RooCmdArg&        arg8    ) 
{
  return dataset.createHistogram ( name.c_str() , 
                                   xvar         , 
                                   arg1         , 
                                   arg2         , 
                                   arg3         , 
                                   arg4         , 
                                   arg5         , 
                                   arg6         , 
                                   arg7         , 
                                   arg8         ) ;
}
// ==================================================================
// The END 
// ==================================================================

