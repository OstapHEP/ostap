// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <cstring>
// ============================================================================
// ROOT
// ============================================================================
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
// ============================================================================
// Ostap
// ===========================================================================
#include "Ostap/GetQuantile.h"
#include "Ostap/Quantile.h"
// ============================================================================
// Local 
// ============================================================================
#include "local_utils.h"
#include "status_codes.h"
// ============================================================================
/** @file
 *  Implementation file for class Ostap::GetQuantile
 *  @see Ostap::HistoProject
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date   2015-10-08
 */
// ============================================================================
// construtctor with the progress flag
// ============================================================================
Ostap::GetQuantile::GetQuantile
( const Ostap::Utils::ProgressConf& progress )
  : Ostap::StatVar ( progress )
{}
// ============================================================================
// single quantile: 
// ============================================================================
/*  get (approximate) quantile using P2 algorithm 
 *  @param data       (INPUT)  input data source 
 *  @param quantile   (UPDATE) quantile  
 *  @param expression (INOUT)  expression 
 *  @param selection  (INPUT)  slection (treated as boolean!)
 *  @param first      (INPUT)  the first event to process (inclusive) 
 *  @param last       (INPUT)  the last event to process (exclusive) 
 *  @param xmin       (INPUT)  low  limit for expressoon 
 *  @param xmax       (INPUT)  high limit for expressoon 
 *  @return status code 
 */
// ============================================================================
Ostap::StatusCode Ostap::GetQuantile::quantile
( TTree*                       data       , 
  Ostap::Math::Quantile&       quantile   ,   
  const std::string&           expression ,  
  const std::string&           selection  ,
  const Ostap::EventIndex      first      ,
  const Ostap::EventIndex      last       ,
  const Ostap::DataType        xmin       ,
  const Ostap::DataType        xmax       ) const
{
  quantile.reset() ;  
  return get_stat ( data       ,
                    quantile   ,
                    expression ,
                    selection  ,
                    first      ,
                    last       ,
                    xmin       ,
                    xmax       ) ;
}
// ============================================================================
/*  get (approximate) quantile using P2 algorithm  
 *  @param data       (INPUT)   input data source 
 *  @param qunaitle   (UPDATE)  quantile  
 *  @param expression (INOUT)   expression 
 *  @param selection  (INPUT)   selection (treatd as boolean!)
 *  @param cut_range  (INPUT)   cut-range 
 *  @param first       (INPUT)  the first event to process (inclusive) 
 *  @param last        (INPUT)  the last event to process (exclusive) 
 *  @param xmin       (INPUT)   low  limit for expressoon 
 *  @param xmax        (INPUT)  high limit for expressoon 
 *  @return status code 
 *  @attention data must be non-weighted!!
 */
// ============================================================================
Ostap::StatusCode Ostap::GetQuantile::quantile
( const RooAbsData*            data       , 
  Ostap::Math::Quantile&       quantile   ,   
  const std::string&           expression ,  
  const std::string&           selection  ,
  const std::string&           cut_range  ,
  const Ostap::EventIndex      first      ,
  const Ostap::EventIndex      last       ,
  const Ostap::DataType        xmin       ,
  const Ostap::DataType        xmax       ) const
{
  quantile.reset() ;  
  return get_stat ( data       ,
                    quantile   ,
                    expression ,
                    selection  ,
                    cut_range  ,
                    first      ,
                    last       ,
                    xmin       ,
                    xmax       ) ;
}
// ============================================================================
// several quantiles: 
// ============================================================================
/*  get (approximate) quantile using P2 algorithm 
 *  @param data       (INPUT)  input data source 
 *  @param quantile   (UPDATE) quantile  
 *  @param expression (INOUT)  expression 
 *  @param selection  (INPUT)  slection (treated as boolean!)
 *  @param first      (INPUT)  the first event to process (inclusive) 
 *  @param last       (INPUT)  the last event to process (exclusive) 
 *  @param xmin       (INPUT)  low  limit for expressoon 
 *  @param xmax       (INPUT)  high limit for expressoon 
 *  @return status code 
 */
// ============================================================================
Ostap::StatusCode Ostap::GetQuantile::quantiles
( TTree*                       data       , 
  Ostap::Math::Quantiles&      quantiles  ,   
  const std::string&           expression ,  
  const std::string&           selection  ,
  const Ostap::EventIndex      first      ,
  const Ostap::EventIndex      last       ,
  const Ostap::DataType        xmin       ,
  const Ostap::DataType        xmax       ) const
{
  quantiles.reset() ;  
  return get_stat ( data       ,
                    quantiles  ,
                    expression ,
                    selection  ,
                    first      ,
                    last       ,
                    xmin       ,
                    xmax       ) ;
}
// ============================================================================
/*  get (approximate) quantile using P2 algorithm  
 *  @param data       (INPUT)   input data source 
 *  @param qunaitle   (UPDATE)  quantile  
 *  @param expression (INOUT)   expression 
 *  @param selection  (INPUT)   selection (treatd as boolean!)
 *  @param cut_range  (INPUT)   cut-range 
 *  @param first       (INPUT)  the first event to process (inclusive) 
 *  @param last        (INPUT)  the last event to process (exclusive) 
 *  @param xmin       (INPUT)   low  limit for expressoon 
 *  @param xmax        (INPUT)  high limit for expressoon 
 *  @return status code 
 *  @attention data must be non-weighted!!
 */
// ============================================================================
Ostap::StatusCode Ostap::GetQuantile::quantiles
( const RooAbsData*            data       , 
  Ostap::Math::Quantiles&      quantiles  ,   
  const std::string&           expression ,  
  const std::string&           selection  ,
  const std::string&           cut_range  ,
  const Ostap::EventIndex      first      ,
  const Ostap::EventIndex      last       ,
  const Ostap::DataType        xmin       ,
  const Ostap::DataType        xmax       ) const
{
  quantiles.reset() ;  
  return get_stat ( data       ,
                    quantiles  ,
                    expression ,
                    selection  ,
                    cut_range  ,
                    first      ,
                    last       ,
                    xmin       ,
                    xmax       ) ;
}


// ============================================================================
//                                                                      The END 
// ============================================================================
