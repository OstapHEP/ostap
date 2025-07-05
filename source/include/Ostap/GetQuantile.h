// ============================================================================
#ifndef OSTAP_GETQUANTILE_H 
#define OSTAP_GETQUANTILE_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <limits>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatusCode.h"
#include "Ostap/StatVar.h"
// ============================================================================
// forward declarations 
// ============================================================================
namespace Ostap { namespace Math { class Quantile  ; } } 
namespace Ostap { namespace Math { class Quantiles ; } } 
// ============================================================================
/** @file Ostap/GetQuantile.h
 *  Collecton of useful metghdo to get quantiles 
 * @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 * @data 2025-06-19
 */
namespace Ostap
{  
  // ===========================================================================
  /* @class GetQuantile 
   * helper class to get quantiles form TTree/RooAbsData  
   * @see TTree 
   * @see RooAbsData 
   */
  class GetQuantile : public StatVar
  {
    // ========================================================================
  public :
    // ========================================================================
    /// construtctor with the progress flag
    GetQuantile ( const Ostap::Utils::ProgressConf& progress = false) ;
    // ========================================================================    
  public: // 
    // ========================================================================
    /** get (approximate) quantile using P2 algorithm 
     *  @param data       (INPUT)  input data source 
     *  @param quantile   (UPDATE) quantile  
     *  @param expression (INOUT)  expression 
     *  @param selection  (INPUT)  slection (treatd as boolean!)
     *  @param first      (INPUT)  the first event to process (inclusive) 
     *  @param last       (INPUT)  the last event to process (exclusive) 
     *  @param xmin       (INPUT)  low  limit for expressoon 
     *  @param xmax       (INPUT)  high limit for expressoon 
     *  @return status code 
     */
    Ostap::StatusCode quantile
      ( TTree*                       data                           , 
        Ostap::Math::Quantile&       quantile                       ,   
        const std::string&           expression                     ,  
        const std::string&           selection  = ""                ,
        const Ostap::EventIndex      first      = Ostap::FirstEvent ,
        const Ostap::EventIndex      last       = Ostap::LastEvent  ,
        const Ostap::DataType        xmin       = Ostap::MinValue   ,
        const Ostap::DataType        xmax       = Ostap::MaxValue   ) const ;
    // =======================================================================
    /** get (approximate) quantile using P2 algorithm  
     *  @param data       (INPUT)   input data source 
     *  @param quantile   (UPDATE)  quantile  
     *  @param expression (INOUT)   expression 
     *  @param selection  (INPUT)   selection (treatd as boolean!)
     *  @param cut_range  (INPUT)   cut-range 
     *  @param first      (INPUT)  the first event to process (inclusive) 
     *  @param last       (INPUT)  the last event to process (exclusive) 
     *  @param xmin       (INPUT)   low  limit for expressoon 
     *  @param xmax       (INPUT)  high limit for expressoon 
     *  @return status code 
     *  @attention data must be non-weighted!!
     */
    Ostap::StatusCode quantile
      ( const RooAbsData*            data                           , 
        Ostap::Math::Quantile&       quantile                       ,   
        const std::string&           expression                     ,  
        const std::string&           selection  = ""                ,
        const std::string&           cut_range  = ""                ,
        const Ostap::EventIndex      first      = Ostap::FirstEvent ,
        const Ostap::EventIndex      last       = Ostap::LastEvent  ,
        const Ostap::DataType        xmin       = Ostap::MinValue   ,
        const Ostap::DataType        xmax       = Ostap::MaxValue   ) const ;
    // =======================================================================
  public: // 
    // ========================================================================
    /** get (approximate) quantile using P2 algorithm 
     *  @param data       (INPUT)  input data source 
     *  @param qunaitles  (UPDATE) quantile  
     *  @param expression (INOUT)  expression 
     *  @param selection  (INPUT)  slection (treatd as boolean!)
     *  @param first      (INPUT)  the first event to process (inclusive) 
     *  @param last       (INPUT)  the last event to process (exclusive) 
     *  @param xmin       (INPUT)  low  limit for expressoon 
     *  @param xmax       (INPUT)  high limit for expressoon 
     *  @return status code 
     */
    Ostap::StatusCode quantiles
      ( TTree*                       data                           , 
        Ostap::Math::Quantiles&      quantiles                      ,   
        const std::string&           expression                     ,  
        const std::string&           selection  = ""                ,
        const Ostap::EventIndex      first      = Ostap::FirstEvent ,
        const Ostap::EventIndex      last       = Ostap::LastEvent  ,
        const Ostap::DataType        xmin       = Ostap::MinValue   ,
        const Ostap::DataType        xmax       = Ostap::MaxValue   ) const ;
    // =======================================================================
    /** get (approximate) quantile using P2 algorithm  
     *  @param data       (INPUT)   input data source 
     *  @param quantiles  (UPDATE)  quantiles
     *  @param expression (INOUT)   expression 
     *  @param selection  (INPUT)   selection (treatd as boolean!)
     *  @param cut_range  (INPUT)   cut-range 
     *  @param first      (INPUT)  the first event to process (inclusive) 
     *  @param last       (INPUT)  the last event to process (exclusive) 
     *  @param xmin       (INPUT)   low  limit for expressoon 
     *  @param xmax       (INPUT)  high limit for expressoon 
     *  @return status code 
     *  @attention data must be non-weighted!!
     */
    Ostap::StatusCode quantiles
      ( const RooAbsData*            data                           , 
        Ostap::Math::Quantiles&      quantiles                      ,   
        const std::string&           expression                     ,  
        const std::string&           selection  = ""                ,
        const std::string&           cut_range  = ""                ,
        const Ostap::EventIndex      first      = Ostap::FirstEvent ,
        const Ostap::EventIndex      last       = Ostap::LastEvent  ,
        const Ostap::DataType        xmin       = Ostap::MinValue   ,
        const Ostap::DataType        xmax       = Ostap::MaxValue   ) const ;
    // =======================================================================
  } ;                                   // The end of class Ostap::GetQuantile
  // ==========================================================================
} //                                                 The enf of namesoace Ostap 
// ============================================================================
#endif // OSTAP_GETQUANTILE_H
// ============================================================================
//                                                                      The END 
// ============================================================================
