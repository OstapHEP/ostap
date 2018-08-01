// ===============================================================================
#ifndef OSTAP_HISTOMAKE_H 
#define OSTAP_HISTOMAKE_H 1
// ===============================================================================
// Include files
// ===============================================================================
// ROOT
// ===============================================================================
#include "RooCmdArg.h"
// ===============================================================================
// forward declaration 
// ===============================================================================
class RooAbsData       ;
class TH1              ; 
class RooAbsRealLValue ;
// ===============================================================================
namespace Ostap
{
  // ==============================================================================
  /** @class HistoMake  Ostap/HistoMake.h
   *  Helper class to "fix" the problem with "masked" 
   *  RooAbsData::createHistogram method
   *
   *  @author Vanya Belyaev
   *  @date   2011-07-16
   */
  class HistoMake  
  {
  public:
    // ================
    /// the only one method 
    static TH1* create_histo 
    ( const RooAbsData&       dataset                   , 
      const std::string&      name                      , 
      const RooAbsRealLValue& xvar                      ,
      const RooCmdArg&        arg1 = RooCmdArg::none () , 
      const RooCmdArg&        arg2 = RooCmdArg::none () , 
      const RooCmdArg&        arg3 = RooCmdArg::none () , 
      const RooCmdArg&        arg4 = RooCmdArg::none () , 
      const RooCmdArg&        arg5 = RooCmdArg::none () , 
      const RooCmdArg&        arg6 = RooCmdArg::none () , 
      const RooCmdArg&        arg7 = RooCmdArg::none () , 
      const RooCmdArg&        arg8 = RooCmdArg::none () ) ;
    // ========================================================================
  } ;
  // ==========================================================================
} //                                                  End of namespace Analysis 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_HISTOMAKE_H
// ============================================================================
