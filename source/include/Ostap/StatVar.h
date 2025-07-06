// ============================================================================
#ifndef OSTAP_STATVAR_H
#define OSTAP_STATVAR_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <limits>
#include <vector>
// ============================================================================
// Forward declarations 
// =============================================================================
class TTree      ; // ROOT 
class TChain     ; // ROOT 
class RooAbsData ; // RooFit
// =============================================================================
// Ostap
// ============================================================================
#include "Ostap/Types.h"
#include "Ostap/StatEntity.h"
#include "Ostap/WStatEntity.h"
#include "Ostap/StatusCode.h"
#include "Ostap/ProgressConf.h"
// ============================================================================
// forward declarations 
// ============================================================================
namespace Ostap { namespace Math { class  Statistic   ; } }
namespace Ostap { namespace Math { class WStatistic   ; } }
namespace Ostap { namespace Math { class  Statistic2  ; } }
namespace Ostap { namespace Math { class WStatistic2  ; } }
namespace Ostap { namespace Math { class  Statistic3  ; } }
namespace Ostap { namespace Math { class WStatistic3  ; } }
namespace Ostap { namespace Math { class  Statistic4  ; } }
namespace Ostap { namespace Math { class WStatistic4  ; } }
namespace Ostap { namespace Math { class  Covariance  ; } }
namespace Ostap { namespace Math { class WCovariance  ; } }
namespace Ostap { namespace Math { class  Covariances ; } }
namespace Ostap { namespace Math { class WCovariances ; } }
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  /** @class StatVar Ostap/StatVar.h
   *  Helper class to get statistical 
   *  infomation  about the variable/expression 
   *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
   *  @date   2013-10-13
   */
  class StatVar 
  {
    // ========================================================================
  public: // helper typedef 
    // ========================================================================
    typedef std::vector<Ostap::StatEntity>   StatVector ;
    typedef std::vector<Ostap::WStatEntity> WStatVector ;
    typedef Ostap::Dict<Ostap::StatEntity>   StatMap    ;
    typedef Ostap::Dict<Ostap::WStatEntity> WStatMap    ;    
    // ========================================================================
  public :
    // ========================================================================
    /// construtctor with the progress flag
    StatVar ( const Ostap::Utils::ProgressConf& progress = false) ;
    // ========================================================================    
  public:  // Statistics&WStatistics 
    // ========================================================================
    /** Fill/update statistical counter 
     *  @param data       (input) data 
     *  @param stat       (UDATE) statistical counter 
     *  @param expression (INPUT) the variable 
     *  @param selection  (INPUT) selection/cut (treated as boolean!)
     *  @param first      (INPUT) the first event to process (inclusive) 
     *  @param last       (INPUT) the last event to process (exclusive) 
     *  @param xmin       (INPUT)  low  limit for expressoon 
     *  @param xmax       (INPUT)  high limit for expressoon 
     *  @return status code 
     *  @see Ostap::Math::Statistics
     *  @attention selection/cut is treated as boolean!
     */
    Ostap::StatusCode get_stat
    ( TTree*                    data                           ,
      Ostap::Math::Statistic&   stat                           ,
      const std::string&        expression                     , 
      const std::string&        selection  = ""                ,
      const Ostap::EventIndex   first      = Ostap::FirstEvent ,
      const Ostap::EventIndex   last       = Ostap::LastEvent  ,
      const Ostap::DataType     xmin       = Ostap::MinValue   ,
      const Ostap::DataType     xmax       = Ostap::MaxValue   ) const ;
    // ========================================================================
    /** Fill/update statistical counter 
     *  @param data       (input) data 
     *  @param stat       (UDATE) statistical counter 
     *  @param expression (INPUT) the variable 
     *  @param selection  (INPUT) selection/cut (treated as weight!)
     *  @param first      (INPUT) the first event to process (inclusibe) 
     *  @param last       (INPUT) the last event to process (exclusive) 
     *  @param xmin       (INPUT)  low  limit for expressoon 
     *  @param xmax       (INPUT)  high limit for expressoon 
     *  @return status code 
     *  @see Ostap::Math::WStatistics
     *  @attention selection/cut is treated as weight!
     */
    Ostap::StatusCode get_stat
    ( TTree*                    data                           ,
      Ostap::Math::WStatistic&  stat                           ,
      const std::string&        expression                     , 
      const std::string&        selection  = ""                ,
      const Ostap::EventIndex   first      = Ostap::FirstEvent ,
      const Ostap::EventIndex   last       = Ostap::LastEvent  ,
      const Ostap::DataType     xmin       = Ostap::MinValue   ,
      const Ostap::DataType     xmax       = Ostap::MaxValue   ) const ;
    // ========================================================================
    /** Fill/update statistical counter 
     *  @param data       (input) data 
     *  @param stat       (UDATE) statistical counter 
     *  @param expression (INPUT) the variable 
     *  @param selection  (INPUT) selection/cut (treated as weight!)
     *  @param cut_range  (INPUT) if non empty: use evene only fomthis cut-range
     *  @param first      (INPUT) the first event to process (inclusibe) 
     *  @param last       (INPUT) the last event to process (exclusive) 
     *  @param xmin       (INPUT)  low  limit for expressoon 
     *  @param xmax       (INPUT)  high limit for expressoon 
     *  @return status code 
     *  @see Ostap::Math::WStatistics
     *  @attention selection/cut is treated as boolean 
     *  @attention datatse must be nonn-weighted 
     */
    Ostap::StatusCode get_stat
    ( const RooAbsData*         data                           ,
      Ostap::Math::Statistic&   stat                           ,
      const std::string&        expression                     , 
      const std::string&        selection  = ""                ,
      const std::string&        cut_range  = ""                ,
      const Ostap::EventIndex   first      = Ostap::FirstEvent ,
      const Ostap::EventIndex   last       = Ostap::LastEvent  , 
      const Ostap::DataType     xmin       = Ostap::MinValue   ,
      const Ostap::DataType     xmax       = Ostap::MaxValue   ) const ;
   // ========================================================================    
    /** Fill/update statistical counter 
     *  @param data       (input) data 
     *  @param stat       (UDATE) statistical counter 
     *  @param expression (INPUT) the variable 
     *  @param selection  (INPUT) selection/cut (treated as weight!)
     *  @param cut_range  (INPUT) if non empty: use evene only fomthis cut-range
     *  @param first      (INPUT) the first event to process (inclusibe) 
     *  @param last       (INPUT) the last event to process (exclusive) 
     *  @param xmin       (INPUT)  low  limit for expressoon 
     *  @param xmax       (INPUT)  high limit for expressoon 
     *  @return status code 
     *  @see Ostap::Math::WStatistics
     *  @attention selection/cut is treated as weight!
     */
    Ostap::StatusCode get_stat
    ( const RooAbsData*         data                           ,
      Ostap::Math::WStatistic&  stat                           ,
      const std::string&        expression                     , 
      const std::string&        selection  = ""                ,
      const std::string&        cut_range  = ""                ,
      const Ostap::EventIndex   first      = Ostap::FirstEvent ,
      const Ostap::EventIndex   last       = Ostap::LastEvent  , 
      const Ostap::DataType     xmin       = Ostap::MinValue   ,
      const Ostap::DataType     xmax       = Ostap::MaxValue   ) const ;
   // ========================================================================    
  public:  // Statistics2 & WStatistics2 
    // ========================================================================
    /** Fill/update statistical counter 
     *  @param data       (input) data 
     *  @param stat       (UDATE) statistical counter 
     *  @param expr1      (INPUT) the 1st variable 
     *  @param expr1      (INPUT) the 2nd variable 
     *  @param selection  (INPUT) selection/cut (treated as boolean!)
     *  @param first      (INPUT) the first event to process (inclusibe) 
     *  @param last       (INPUT) the last event to process (exclusive) 
     *  @param xmin       (INPUT)  low  limit for expressoon 
     *  @param xmax       (INPUT)  high limit for expressoon 
     *  @param ymin       (INPUT)  low  limit for expressoon 
     *  @param ymax       (INPUT)  high limit for expressoon 
     *  @return status code 
     *  @see Ostap::Math::Statistics2
     *  @attention selection/cut is treated as boolean!
     */
    Ostap::StatusCode get_stat
    ( TTree*                    data                           ,
      Ostap::Math::Statistic2&  stat                           ,
      const std::string&        expr1                          , 
      const std::string&        expr2                          ,       
      const std::string&        selection  = ""                ,
      const Ostap::EventIndex   first      = Ostap::FirstEvent ,
      const Ostap::EventIndex   last       = Ostap::LastEvent  ,      
      const Ostap::DataType     xmin       = Ostap::MinValue   ,
      const Ostap::DataType     xmax       = Ostap::MaxValue   , 
      const Ostap::DataType     ymin       = Ostap::MinValue   ,
      const Ostap::DataType     ymax       = Ostap::MaxValue   ) const ;    
    // ========================================================================
    /** Fill/update statistical counter 
     *  @param data       (input) data 
     *  @param stat       (UDATE) statistical counter 
     *  @param expr1      (INPUT) the 1st variable 
     *  @param expr1      (INPUT) the 2nd variable 
     *  @param selection  (INPUT) selection/cut (treated as weight!)
     *  @param first      (INPUT) the first event to process (inclusibe) 
     *  @param last       (INPUT) the last event to process (exclusive) 
     *  @param xmin       (INPUT)  low  limit for expressoon 
     *  @param xmax       (INPUT)  high limit for expressoon 
     *  @param ymin       (INPUT)  low  limit for expressoon 
     *  @param ymax       (INPUT)  high limit for expressoon 
     *  @return status code 
     *  @see Ostap::Math::WStatistics2
     *  @attention selection/cut is treated as weight!
     */
    Ostap::StatusCode get_stat
    ( TTree*                    data                           ,
      Ostap::Math::WStatistic2& stat                           ,
      const std::string&        expr1                          , 
      const std::string&        expr2                          ,       
      const std::string&        selection  = ""                ,
      const Ostap::EventIndex   first      = Ostap::FirstEvent ,
      const Ostap::EventIndex   last       = Ostap::LastEvent  , 
      const Ostap::DataType     xmin       = Ostap::MinValue   ,
      const Ostap::DataType     xmax       = Ostap::MaxValue   , 
      const Ostap::DataType     ymin       = Ostap::MinValue   ,
      const Ostap::DataType     ymax       = Ostap::MaxValue   ) const ;
    // ========================================================================
    /** Fill/update statistical counter 
     *  @param data       (input) data 
     *  @param stat       (UDATE) statistical counter 
     *  @param expr1       (INPUT) the 1st variable 
     *  @param expr1      (INPUT) the 2nd variable 
     *  @param selection  (INPUT) selection/cut (treated as boolean!)
     *  @param cut_range  (INPUT) if non empty: use evene only fomthis cut-range
     *  @param first      (INPUT) the first event to process (inclusibe) 
     *  @param last       (INPUT) the last event to process (exclusive) 
     *  @param xmin       (INPUT)  low  limit for expressoon 
     *  @param xmax       (INPUT)  high limit for expressoon 
     *  @param ymin       (INPUT)  low  limit for expressoon 
     *  @param ymax       (INPUT)  high limit for expressoon 
     *  @return status code 
     *  @see Ostap::Math::WStatistics2
     *  @attention selection/cut is treated as boolean 
     *  @attention dat amust be non-weighted 
     */
    Ostap::StatusCode get_stat
    ( const RooAbsData*         data                           ,
      Ostap::Math::Statistic2&  stat                           ,
      const std::string&        expr1                          , 
      const std::string&        expr2                          ,       
      const std::string&        selection  = ""                ,
      const std::string&        cut_range  = ""                ,
      const Ostap::EventIndex   first      = Ostap::FirstEvent ,
      const Ostap::EventIndex   last       = Ostap::LastEvent  ,       
      const Ostap::DataType     xmin       = Ostap::MinValue   ,
      const Ostap::DataType     xmax       = Ostap::MaxValue   , 
      const Ostap::DataType     ymin       = Ostap::MinValue   ,
      const Ostap::DataType     ymax       = Ostap::MaxValue   ) const ;
    // ========================================================================    
    /** Fill/update statistical counter 
     *  @param data       (input) data 
     *  @param stat       (UDATE) statistical counter 
     *  @param expr1       (INPUT) the 1st variable 
     *  @param expr1      (INPUT) the 2nd variable 
     *  @param selection  (INPUT) selection/cut (treated as weight!)
     *  @param cut_range  (INPUT) if non empty: use evene only fomthis cut-range
     *  @param first      (INPUT) the first event to process (inclusibe) 
     *  @param last       (INPUT) the last event to process (exclusive) 
     *  @param xmin       (INPUT)  low  limit for expressoon 
     *  @param xmax       (INPUT)  high limit for expressoon 
     *  @param ymin       (INPUT)  low  limit for expressoon 
     *  @param ymax       (INPUT)  high limit for expressoon 
     *  @return status code 
     *  @see Ostap::Math::WStatistics2
     *  @attention selection/cut is treated as weight!
     */
    Ostap::StatusCode get_stat
    ( const RooAbsData*         data                           ,
      Ostap::Math::WStatistic2& stat                           ,
      const std::string&        expr1                          , 
      const std::string&        expr2                          ,       
      const std::string&        selection  = ""                ,
      const std::string&        cut_range  = ""                ,
      const Ostap::EventIndex   first      = Ostap::FirstEvent ,
      const Ostap::EventIndex   last       = Ostap::LastEvent  ,       
      const Ostap::DataType     xmin       = Ostap::MinValue   ,
      const Ostap::DataType     xmax       = Ostap::MaxValue   , 
      const Ostap::DataType     ymin       = Ostap::MinValue   ,
      const Ostap::DataType     ymax       = Ostap::MaxValue   ) const ;
    // ========================================================================    
  public:  // Statistics3 & WStatistics3
    // ========================================================================
    /** Fill/update statistical counter 
     *  @param data       (input) data 
     *  @param stat       (UDATE) statistical counter 
     *  @param expr1      (INPUT) the 1st variable 
     *  @param expr1      (INPUT) the 2nd variable 
     *  @param expr3      (INPUT) the 3rd variable 
     *  @param selection  (INPUT) selection/cut (treated as boolean!)
     *  @param first      (INPUT) the first event to process (inclusibe) 
     *  @param last       (INPUT) the last event to process (exclusive) 
     *  @param xmin       (INPUT)  low  limit for expressoon 
     *  @param xmax       (INPUT)  high limit for expressoon 
     *  @param ymin       (INPUT)  low  limit for expressoon 
     *  @param ymax       (INPUT)  high limit for expressoon 
     *  @param zmin       (INPUT)  low  limit for expressoon 
     *  @param zmax       (INPUT)  high limit for expressoon 
     *  @return status code 
     *  @see Ostap::Math::Statistics3
     *  @attention selection/cut is treated as boolean!
     */
    Ostap::StatusCode get_stat
    ( TTree*                    data                           ,
      Ostap::Math::Statistic3&  stat                           ,
      const std::string&        expr1                          , 
      const std::string&        expr2                          ,       
      const std::string&        expr3                          , 
      const std::string&        selection  = ""                ,
      const Ostap::EventIndex   first      = Ostap::FirstEvent ,
      const Ostap::EventIndex   last       = Ostap::LastEvent  , 
      const Ostap::DataType     xmin       = Ostap::MinValue   ,
      const Ostap::DataType     xmax       = Ostap::MaxValue   , 
      const Ostap::DataType     ymin       = Ostap::MinValue   ,
      const Ostap::DataType     ymax       = Ostap::MaxValue   , 
      const Ostap::DataType     zmin       = Ostap::MinValue   ,
      const Ostap::DataType     zmax       = Ostap::MaxValue   ) const ;
    // ========================================================================
    /** Fill/update statistical counter 
     *  @param data       (input) data 
     *  @param stat       (UDATE) statistical counter 
     *  @param expr1      (INPUT) the 1st variable 
     *  @param expr1      (INPUT) the 2nd variable 
     *  @param expr3      (INPUT) the 3rd variable  
     *  @param selection  (INPUT) selection/cut (treated as weight!)
     *  @param first      (INPUT) the first event to process (inclusibe) 
     *  @param last       (INPUT) the last event to process (exclusive) 
     *  @param xmin       (INPUT)  low  limit for expressoon 
     *  @param xmax       (INPUT)  high limit for expressoon 
     *  @param ymin       (INPUT)  low  limit for expressoon 
     *  @param ymax       (INPUT)  high limit for expressoon 
     *  @param zmin       (INPUT)  low  limit for expressoon 
     *  @param zmax       (INPUT)  high limit for expressoon 
     *  @return status code 
     *  @see Ostap::Math::WStatistics3
     *  @attention selection/cut is treated as weight!
     */
    Ostap::StatusCode get_stat
    ( TTree*                    data                           ,
      Ostap::Math::WStatistic3& stat                           ,
      const std::string&        expr1                          , 
      const std::string&        expr2                          ,
      const std::string&        expr3                          ,       
      const std::string&        selection  = ""                ,
      const Ostap::EventIndex   first      = Ostap::FirstEvent ,
      const Ostap::EventIndex   last       = Ostap::LastEvent  , 
      const Ostap::DataType     xmin       = Ostap::MinValue   ,
      const Ostap::DataType     xmax       = Ostap::MaxValue   , 
      const Ostap::DataType     ymin       = Ostap::MinValue   ,
      const Ostap::DataType     ymax       = Ostap::MaxValue   , 
      const Ostap::DataType     zmin       = Ostap::MinValue   ,
      const Ostap::DataType     zmax       = Ostap::MaxValue   ) const ;
    // ========================================================================
    /** Fill/update statistical counter 
     *  @param data       (input) data 
     *  @param stat       (UDATE) statistical counter 
     *  @param expr1      (INPUT) the 1st variable 
     *  @param expr1      (INPUT) the 2nd variable 
     *  @param expr3      (INPUT) the 3rd variable 
     *  @param selection  (INPUT) selection/cut (treated as boolean !)
     *  @param cut_range  (INPUT) if non empty: use evene only fomthis cut-range
     *  @param first      (INPUT) the first event to process (inclusibe) 
     *  @param last       (INPUT) the last event to process (exclusive) 
     *  @param xmin       (INPUT)  low  limit for expressoon 
     *  @param xmax       (INPUT)  high limit for expressoon 
     *  @param ymin       (INPUT)  low  limit for expressoon 
     *  @param ymax       (INPUT)  high limit for expressoon 
     *  @param zmin       (INPUT)  low  limit for expressoon 
     *  @param zmax       (INPUT)  high limit for expressoon 
     *  @return status code 
     *  @see Ostap::Math::WStatistics3
     *  @attention selection/cut is treated as weight!
     *  @attention data MUST be non-weighted 
     */
    Ostap::StatusCode get_stat
    ( const RooAbsData*         data                           ,
      Ostap::Math::Statistic3&  stat                           ,
      const std::string&        expr1                          , 
      const std::string&        expr2                          , 
      const std::string&        expr3                          ,      
      const std::string&        selection  = ""                ,
      const std::string&        cut_range  = ""                ,
      const Ostap::EventIndex   first      = Ostap::FirstEvent ,
      const Ostap::EventIndex   last       = Ostap::LastEvent  , 
      const Ostap::DataType     xmin       = Ostap::MinValue   ,
      const Ostap::DataType     xmax       = Ostap::MaxValue   , 
      const Ostap::DataType     ymin       = Ostap::MinValue   ,
      const Ostap::DataType     ymax       = Ostap::MaxValue   , 
      const Ostap::DataType     zmin       = Ostap::MinValue   ,
      const Ostap::DataType     zmax       = Ostap::MaxValue   ) const ;
    // ======================================================================    
    /** Fill/update statistical counter 
     *  @param data       (input) data 
     *  @param stat       (UDATE) statistical counter 
     *  @param expr1      (INPUT) the 1st variable 
     *  @param expr1      (INPUT) the 2nd variable 
     *  @param expr3      (INPUT) the 3rd variable 
     *  @param selection  (INPUT) selection/cut (treated as weight!)
     *  @param cut_range  (INPUT) if non empty: use evene only fomthis cut-range
     *  @param first      (INPUT) the first event to process (inclusibe) 
     *  @param last       (INPUT) the last event to process (exclusive) 
     *  @param xmin       (INPUT)  low  limit for expressoon 
     *  @param xmax       (INPUT)  high limit for expressoon 
     *  @param ymin       (INPUT)  low  limit for expressoon 
     *  @param ymax       (INPUT)  high limit for expressoon 
     *  @param zmin       (INPUT)  low  limit for expressoon 
     *  @param zmax       (INPUT)  high limit for expressoon 
     *  @return status code 
     *  @see Ostap::Math::WStatistics3
     *  @attention selection/cut is treated as weight!
     */
    Ostap::StatusCode get_stat
    ( const RooAbsData*         data                           ,
      Ostap::Math::WStatistic3& stat                           ,
      const std::string&        expr1                          , 
      const std::string&        expr2                          , 
      const std::string&        expr3                          ,      
      const std::string&        selection  = ""                ,
      const std::string&        cut_range  = ""                ,
      const Ostap::EventIndex   first      = Ostap::FirstEvent ,
      const Ostap::EventIndex   last       = Ostap::LastEvent  , 
      const Ostap::DataType     xmin       = Ostap::MinValue   ,
      const Ostap::DataType     xmax       = Ostap::MaxValue   , 
      const Ostap::DataType     ymin       = Ostap::MinValue   ,
      const Ostap::DataType     ymax       = Ostap::MaxValue   , 
      const Ostap::DataType     zmin       = Ostap::MinValue   ,
      const Ostap::DataType     zmax       = Ostap::MaxValue   ) const ;
    // ======================================================================    
  public:  // Statistics4 & WStatistics4
    // ======================================================================
    /** Fill/update statistical counter 
     *  @param data       (input) data 
     *  @param stat       (UDATE) statistical counter 
     *  @param expr1      (INPUT) the 1st variable 
     *  @param expr1      (INPUT) the 2nd variable 
     *  @param expr3      (INPUT) the 3rd variable 
     *  @param expr4      (INPUT) the 4th variable   
     *  @param selection  (INPUT) selection/cut (treated as boolean!)
     *  @param first      (INPUT) the first event to process (inclusibe) 
     *  @param last       (INPUT) the last event to process (exclusive) 
     *  @param xmin       (INPUT)  low  limit for expressoon 
     *  @param xmax       (INPUT)  high limit for expressoon 
     *  @param ymin       (INPUT)  low  limit for expressoon 
     *  @param ymax       (INPUT)  high limit for expressoon 
     *  @param zmin       (INPUT)  low  limit for expressoon 
     *  @param zmax       (INPUT)  high limit for expressoon 
     *  @param tmin       (INPUT)  low  limit for expressoon 
     *  @param tmax       (INPUT)  high limit for expressoon 
     *  @return status code 
     *  @see Ostap::Math::Statistics4
     *  @attention selection/cut is treated as boolean!
     */
    Ostap::StatusCode get_stat
    ( TTree*                    data                           ,
      Ostap::Math::Statistic4&  stat                           ,
      const std::string&        expr1                          , 
      const std::string&        expr2                          ,       
      const std::string&        expr3                          , 
      const std::string&        expr4                          , 
      const std::string&        selection  = ""                ,
      const Ostap::EventIndex   first      = Ostap::FirstEvent ,
      const Ostap::EventIndex   last       = Ostap::LastEvent  , 
      const Ostap::DataType     xmin       = Ostap::MinValue   ,
      const Ostap::DataType     xmax       = Ostap::MaxValue   , 
      const Ostap::DataType     ymin       = Ostap::MinValue   ,
      const Ostap::DataType     ymax       = Ostap::MaxValue   , 
      const Ostap::DataType     zmin       = Ostap::MinValue   ,
      const Ostap::DataType     zmax       = Ostap::MaxValue   , 
      const Ostap::DataType     tmin       = Ostap::MinValue   ,
      const Ostap::DataType     tmax       = Ostap::MaxValue   ) const ;
    // ========================================================================
    /** Fill/update statistical counter 
     *  @param data       (input) data 
     *  @param stat       (UDATE) statistical counter 
     *  @param expr1      (INPUT) the 1st variable 
     *  @param expr1      (INPUT) the 2nd variable 
     *  @param expr3      (INPUT) the 3rd variable  
     *  @param expr4      (INPUT) the 4th variable 
     *  @param selection  (INPUT) selection/cut (treated as weight!)
     *  @param first      (INPUT) the first event to process (inclusibe) 
     *  @param last       (INPUT) the last event to process (exclusive) 
     *  @param xmin       (INPUT)  low  limit for expressoon 
     *  @param xmax       (INPUT)  high limit for expressoon 
     *  @param ymin       (INPUT)  low  limit for expressoon 
     *  @param ymax       (INPUT)  high limit for expressoon 
     *  @param zmin       (INPUT)  low  limit for expressoon 
     *  @param zmax       (INPUT)  high limit for expressoon 
     *  @param tmin       (INPUT)  low  limit for expressoon 
     *  @param tmax       (INPUT)  high limit for expressoon 
     *  @return status code 
     *  @see Ostap::Math::WStatistics4
     *  @attention selection/cut is treated as weight!
     */
    Ostap::StatusCode get_stat
    ( TTree*                    data                           ,
      Ostap::Math::WStatistic4& stat                           ,
      const std::string&        expr1                          , 
      const std::string&        expr2                          ,
      const std::string&        expr3                          ,
      const std::string&        expr4                          ,        
      const std::string&        selection  = ""                ,
      const Ostap::EventIndex   first      = Ostap::FirstEvent ,
      const Ostap::EventIndex   last       = Ostap::LastEvent  , 
      const Ostap::DataType     xmin       = Ostap::MinValue   ,
      const Ostap::DataType     xmax       = Ostap::MaxValue   , 
      const Ostap::DataType     ymin       = Ostap::MinValue   ,
      const Ostap::DataType     ymax       = Ostap::MaxValue   , 
      const Ostap::DataType     zmin       = Ostap::MinValue   ,
      const Ostap::DataType     zmax       = Ostap::MaxValue   , 
      const Ostap::DataType     tmin       = Ostap::MinValue   ,
      const Ostap::DataType     tmax       = Ostap::MaxValue   ) const ;
    // ========================================================================
    /** Fill/update statistical counter 
     *  @param data       (input) data 
     *  @param stat       (UDATE) statistical counter 
     *  @param expr1      (INPUT) the 1st variable 
     *  @param expr1      (INPUT) the 2nd variable 
     *  @param expr3      (INPUT) the 3rd variable 
     *  @param expr4      (INPUT) the 4th variable 
     *  @param selection  (INPUT) selection/cut (treated as boolean)
     *  @param cut_range  (INPUT) if non empty: use evene only fomthis cut-range
     *  @param first      (INPUT) the first event to process (inclusibe) 
     *  @param last       (INPUT) the last event to process (exclusive) 
     *  @param xmin       (INPUT)  low  limit for expressoon 
     *  @param xmax       (INPUT)  high limit for expressoon 
     *  @param ymin       (INPUT)  low  limit for expressoon 
     *  @param ymax       (INPUT)  high limit for expressoon 
     *  @param zmin       (INPUT)  low  limit for expressoon 
     *  @param zmax       (INPUT)  high limit for expressoon 
     *  @param tmin       (INPUT)  low  limit for expressoon 
     *  @param tmax       (INPUT)  high limit for expressoon 
     *  @return status code 
     *  @see Ostap::Math::WStatistics2
     *  @attention selection/cut is treated as boolean 
     *  @attention data must be non-weighted 
     */
    Ostap::StatusCode get_stat
    ( const RooAbsData*         data                           ,
      Ostap::Math::Statistic4&  stat                           ,
      const std::string&        expr1                          , 
      const std::string&        expr2                          , 
      const std::string&        expr3                          ,
      const std::string&        expr4                          ,       
      const std::string&        selection  = ""                ,
      const std::string&        cut_range  = ""                ,
      const Ostap::EventIndex   first      = Ostap::FirstEvent ,
      const Ostap::EventIndex   last       = Ostap::LastEvent  , 
      const Ostap::DataType     xmin       = Ostap::MinValue   ,
      const Ostap::DataType     xmax       = Ostap::MaxValue   , 
      const Ostap::DataType     ymin       = Ostap::MinValue   ,
      const Ostap::DataType     ymax       = Ostap::MaxValue   , 
      const Ostap::DataType     zmin       = Ostap::MinValue   ,
      const Ostap::DataType     zmax       = Ostap::MaxValue   , 
      const Ostap::DataType     tmin       = Ostap::MinValue   ,
      const Ostap::DataType     tmax       = Ostap::MaxValue   ) const ;
    // ========================================================================
    /** Fill/update statistical counter 
     *  @param data       (input) data 
     *  @param stat       (UDATE) statistical counter 
     *  @param expr1      (INPUT) the 1st variable 
     *  @param expr1      (INPUT) the 2nd variable 
     *  @param expr3      (INPUT) the 3rd variable 
     *  @param expr4      (INPUT) the 4th variable 
     *  @param selection  (INPUT) selection/cut (treated as weight!)
     *  @param cut_range  (INPUT) if non empty: use evene only fomthis cut-range
     *  @param first      (INPUT) the first event to process (inclusibe) 
     *  @param last       (INPUT) the last event to process (exclusive) 
     *  @param xmin       (INPUT)  low  limit for expressoon 
     *  @param xmax       (INPUT)  high limit for expressoon 
     *  @param ymin       (INPUT)  low  limit for expressoon 
     *  @param ymax       (INPUT)  high limit for expressoon 
     *  @param zmin       (INPUT)  low  limit for expressoon 
     *  @param zmax       (INPUT)  high limit for expressoon 
     *  @param tmin       (INPUT)  low  limit for expressoon 
     *  @param tmax       (INPUT)  high limit for expressoon 
     *  @return status code 
     *  @see Ostap::Math::WStatistics2
     *  @attention selection/cut is treated as weight!
     */
    Ostap::StatusCode get_stat
    ( const RooAbsData*         data                           ,
      Ostap::Math::WStatistic4& stat                           ,
      const std::string&        expr1                          , 
      const std::string&        expr2                          , 
      const std::string&        expr3                          ,
      const std::string&        expr4                          ,       
      const std::string&        selection  = ""                ,
      const std::string&        cut_range  = ""                ,
      const Ostap::EventIndex   first      = Ostap::FirstEvent ,
      const Ostap::EventIndex   last       = Ostap::LastEvent  , 
      const Ostap::DataType     xmin       = Ostap::MinValue   ,
      const Ostap::DataType     xmax       = Ostap::MaxValue   , 
      const Ostap::DataType     ymin       = Ostap::MinValue   ,
      const Ostap::DataType     ymax       = Ostap::MaxValue   , 
      const Ostap::DataType     zmin       = Ostap::MinValue   ,
      const Ostap::DataType     zmax       = Ostap::MaxValue   , 
      const Ostap::DataType     tmin       = Ostap::MinValue   ,
      const Ostap::DataType     tmax       = Ostap::MaxValue   ) const ;
    // ========================================================================
  public: /// statistc for single variable
    // =========================================================================
    /** Is there at leats one good entry in dataset ?
     *  @param data       (INPUT) input data 
     *  @param selecttion (INPUT) selection criteria 
     *  @param first      (INPUT) the first event to process (inclusibe) 
     *  @param last       (INPUT) the last event to process (exclusive) 
     */
    bool hasEntry 
    ( TTree*                    data                           , 
      const std::string         selection  = ""                , 
      const Ostap::EventIndex   first      = Ostap::FirstEvent ,
      const Ostap::EventIndex   last       = Ostap::LastEvent  ) const ;
    // ========================================================================
    /** Is there at leats one good entry in dataset ?
     *  @param data       (INPUT) input data 
     *  @param selecttion (INPUT) selection criteria 
     *  @param first      (INPUT) the first event to process (inclusibe) 
     *  @param last       (INPUT) the last event to process (exclusive) 
     */
    bool hasEntry 
    ( const RooAbsData*         data                           , 
      const std::string&        selection  = ""                ,
      const std::string&        cut_range  = ""                ,   
      const Ostap::EventIndex   first      = Ostap::FirstEvent ,
      const Ostap::EventIndex   last       = Ostap::LastEvent  ) const ;    
    // ========================================================================
  public: /// statistc for single variable
    // =========================================================================
    /** get the nmnber of good entries 
     *  @param data (INPUT) input dat 
     *  @param selecttion (INPUT) selection criteria 
     *  @param first      (INPUT) the first event to process (inclusibe) 
     *  @param last       (INPUT) the last event to process (exclusive) 
     */
    Ostap::EventIndex size
    ( TTree*                  data           ,
      const std::string&      selection = "" ,
      const Ostap::EventIndex first     = Ostap::FirstEvent ,
      const Ostap::EventIndex last      = Ostap::LastEvent  ) const ;
    // ======================================================================
    /** get the nmnber of good entries 
     *  @param data (INPUT) input dat 
     *  @param selecttion (INPUT) selection criteria 
     *  @param first      (INPUT) the first event to process (inclusibe) 
     *  @param last       (INPUT) the last event to process (exclusive) 
     */
    Ostap::EventIndex size
    ( const RooAbsData*       data           ,
      const std::string&      selection = "" ,
      const std::string&      cut_range = "" , 
      const Ostap::EventIndex first     = Ostap::FirstEvent ,
      const Ostap::EventIndex last      = Ostap::LastEvent  ) const ;    
    // ========================================================================
  public: /// statistc for single variable
    // =========================================================================
    /** build statistic for the <code>expression</code>
     *  @param tree       (INPUT) the tree 
     *  @param expression (INPUT) the expression
     *  @param first      (INPUT) the first entry 
     *  @param last       (INPUT) the last entry      *
     *  @param xmin       (INPUT)  low  limit for expressoon 
     *  @param xmax       (INPUT)  high limit for expressoon 
     *
     *  @code
     *  tree = ... 
     *  stat = tree.statVar( 'S_sw' ) 
     *  @endcode 
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2013-10-13
     */
    Ostap::StatEntity statVar
    ( TTree*                  tree                           , 
      const std::string&      expression                     , 
      const Ostap::EventIndex first      = Ostap::FirstEvent ,
      const Ostap::EventIndex last       = Ostap::LastEvent  , 
      const Ostap::DataType   xmin       = Ostap::MinValue   ,
      const Ostap::DataType   xmax       = Ostap::MaxValue   ) const ;
    // =================================================================
    /** build statistic for the <code>expression</code>
     *  @param tree       (INPUT) the tree 
     *  @param expression (INPUT) the expression
     *  @param selection  (INPUT) selection/ (as boolean) 
     *  @param first      (INPUT) the first entry 
     *  @param last       (INPUT) the last entry      *
     *  @param xmin       (INPUT)  low  limit for expressoon 
     *  @param xmax       (INPUT)  high limit for expressoon 
     *
     *  @code
     *  tree = ... 
     *  stat = tree.statVar( 'S_sw' ) 
     *  @endcode 
     *  @attetion sleection is treated as boolean! 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2013-10-13
     */
    Ostap::StatEntity statVar_cut
    ( TTree*                  tree                           , 
      const std::string&      expression                     , 
      const std::string&      selection  = ""                , 
      const Ostap::EventIndex first      = Ostap::FirstEvent ,
      const Ostap::EventIndex last       = Ostap::LastEvent  , 
      const Ostap::DataType   xmin       = Ostap::MinValue   ,
      const Ostap::DataType   xmax       = Ostap::MaxValue   ) const ;      
    // =================================================================
    /** build statistic for the <code>expression</code>
     *  @param tree       (INPUT) the tree 
     *  @param expression (INPUT) the expression
     *  @param selection  (INPUT) selection/weight 
     *  @param first      (INPUT) the first entry 
     *  @param last       (INPUT) the last entry      *
     *  @param xmin       (INPUT)  low  limit for expressoon 
     *  @param xmax       (INPUT)  high limit for expressoon 
     *
     *  @code
     *  tree = ... 
     *  stat = tree.statVar( 'S_sw' ) 
     *  @endcode 
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2013-10-13
     */
    Ostap::WStatEntity statVar
    ( TTree*                  tree                           , 
      const std::string&      expression                     , 
      const std::string&      selection                      , 
      const Ostap::EventIndex first      = Ostap::FirstEvent ,
      const Ostap::EventIndex last       = Ostap::LastEvent  , 
      const Ostap::DataType   xmin       = Ostap::MinValue   ,
      const Ostap::DataType   xmax       = Ostap::MaxValue   ) const ;
    // ===============================================================
    /** build statistic for the <code>expression</code>
     *  @param tree       (INPUT) the tree 
     *  @param expression (INPUT) the expression
     *  @param selection  (INPUT) seleciton/weight 
     *  @param cut_range  (INPUT) cut-range 
     *  @param first      (INPUT) the first entry 
     *  @param last       (INPUT) the last entry      *
     *  @param xmin       (INPUT)  low  limit for expressoon 
     *  @param xmax       (INPUT)  high limit for expressoon 
     *
     *  @code
     *  data = ... 
     *  data = data.statVar( 'S_sw' ) 
     *  @endcode 
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2013-10-13
     */
    Ostap::WStatEntity statVar 
    ( const RooAbsData*       data                           , 
      const std::string&      expression                     , 
      const std::string&      selection                      ,
      const std::string&      cut_range  = ""                ,  
      const Ostap::EventIndex first      = Ostap::FirstEvent ,
      const Ostap::EventIndex last       = Ostap::LastEvent  , 
      const Ostap::DataType   xmin       = Ostap::MinValue   ,
      const Ostap::DataType   xmax       = Ostap::MaxValue   ) const ;
    // =========================================================================
  public: /// get infomation about several variables simultabneously 
    // =========================================================================
    /** build statistic for the <code>expressions</code>
     *  @param data        (INPUT)  input data 
     *  @param result      (UPDATE) the output statistics for specified expressions 
     *  @param expressions (INPUT)  the list of  expressions
     *  @param selection   (INPUT)  selection 
     *  @param first       (INPUT)  the first entry to process 
     *  @param last        (INPUT)  the last entry to process (not including!)
     *  @return status code 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @attentoon selection is treated as boolean 
     *  @date   2018-11-04
     */
    Ostap::StatusCode statVars
    ( TTree*                       data                      ,  
      Ostap::StatVar::StatVector&  result                    ,  
      const Ostap::Strings&        expressions               ,
      const std::string&           selection   = ""            , 
      const Ostap::EventIndex      first       = Ostap::FirstEvent ,
      const Ostap::EventIndex      last        = Ostap::LastEvent  ) const ;
    // =========================================================================
    /**  build statistic for the <code>expressions</code>
     *  @param data        (INPUT)  iiput data 
     *  @param result      (UPDATE) the output statistics for specified expressions 
     *  @param expressions (INPUT)  the list of  expressions
     *  @param selection   (INPUT)  selection 
     *  @param first       (INPUT)  the first entry to process 
     *  @param last        (INPUT)  the last entry to process (not including!)
     *  @return statsu code 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @attentoon selection is treated as weight
     *  @date   2018-11-04
     */
    Ostap::StatusCode statVars
    ( TTree*                       data                            ,  
      Ostap::StatVar::WStatVector& result                          ,  
      const Ostap::Strings&        expressions                     ,
      const std::string&           selection   = ""                , 
      const Ostap::EventIndex      first       = Ostap::FirstEvent ,
      const Ostap::EventIndex      last        = Ostap::LastEvent  ) const ;
    // ========================================================================
    /**  build statistic for the <code>expressions</code>
     *  @param data        (INPUT)  iiput data 
     *  @param result      (UPDATE) the output statistics for specified expressions 
     *  @param expressions (INPUT)  the list of  expressions
     *  @param selection   (INPUT)  selection 
     *  @param cut_range   (INPUT)  the cut-range 
     *  @param first       (INPUT)  the first entry to process 
     *  @param last        (INPUT)  the last entry to process (not including!)
     *  @return statsu code 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2018-11-04
     */
    Ostap::StatusCode statVars
    ( const RooAbsData*            data        , 
      Ostap::StatVar::WStatVector& result      , 
      const Ostap::Strings&        expressions ,
      const std::string&           cuts        = "" ,
      const std::string&           cut_range   = "" ,
      const Ostap::EventIndex      first       = Ostap::FirstEvent ,
      const Ostap::EventIndex      last        = Ostap::LastEvent  ) const ;
    // ========================================================================
  public : // covariances
    // =========================================================================
    /** calculate the covariance of two expressions 
     *  @param tree  (INPUT)  the inpout tree 
     *  @param exp1  (INPUT)  the first  expression
     *  @param exp2  (INPUT)  the second expression
     *  @param selection 
     *  @param first (INPUT)  the first  event to process (inclusive) 
     *  @param last  (INPUT)  the last   event to process (exc;usive) 
     *  @param xmin  (INPUT)  low  limit for expressoon 
     *  @param xmax  (INPUT)  high limit for expressoon 
     *  @param ymin  (INPUT)  low  limit for expressoon 
     *  @param ymax  (INPUT)  high limit for expressoon 
     *  @return Covariance 
     *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
     *  @date   2024-07-22-
     */
    Ostap::StatusCode statCov
    ( TTree*                   tree                          , 
      Ostap::Math::Covariance& stat                          , 
      const std::string&       exp1                          , 
      const std::string&       exp2                          , 
      const std::string&       selection = ""                , 
      const Ostap::EventIndex  first     = Ostap::FirstEvent ,
      const Ostap::EventIndex  last      = Ostap::LastEvent  , 
      const Ostap::DataType    xmin      = Ostap::MinValue   ,
      const Ostap::DataType    xmax      = Ostap::MaxValue   , 
      const Ostap::DataType    ymin      = Ostap::MinValue   ,
      const Ostap::DataType    ymax      = Ostap::MaxValue   ) const ;
    // ======================================================================    
    /** calculate the covariance of two expressions 
     *  @param tree  (INPUT)  the inpout tree 
     *  @param exp1  (INPUT)  the first  expression
     *  @param exp2  (INPUT)  the second expression
     *  @param selection 
     *  @param first (INPUT)  the first  event to process (inclusive) 
     *  @param last  (INPUT)  the last   event to process (exc;usive) 
     *  @param xmin  (INPUT)  low  limit for expressoon 
     *  @param xmax  (INPUT)  high limit for expressoon 
     *  @param ymin  (INPUT)  low  limit for expressoon 
     *  @param ymax  (INPUT)  high limit for expressoon 
     *  @return Covariance 
     *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
     *  @date   2024-07-22-
     */
    Ostap::StatusCode statCov 
    ( TTree*                    tree       , 
      Ostap::Math::WCovariance& stat       , 
      const std::string&        exp1       , 
      const std::string&        exp2       ,
      const std::string&        selection  , 
      const Ostap::EventIndex   first      = Ostap::FirstEvent ,
      const Ostap::EventIndex   last       = Ostap::LastEvent  , 
      const Ostap::DataType     xmin       = Ostap::MinValue   ,
      const Ostap::DataType     xmax       = Ostap::MaxValue   , 
      const Ostap::DataType     ymin       = Ostap::MinValue   ,
      const Ostap::DataType     ymax       = Ostap::MaxValue   ) const ;
    // ========================================================================
    /** calculate the covariance of two expressions 
     *  @param tree  (INPUT)  the inpout tree 
     *  @param exp1  (INPUT)  the first  expression
     *  @param exp2  (INPUT)  the second expression
     *  @param selection 
     *  @param first (INPUT)  the first  event to process (inclusive) 
     *  @param last  (INPUT)  the last   event to process (exc;usive) 
     *  @param xmin  (INPUT)  low  limit for expressoon 
     *  @param xmax  (INPUT)  high limit for expressoon 
     *  @param ymin  (INPUT)  low  limit for expressoon 
     *  @param ymax  (INPUT)  high limit for expressoon 
     *  @return Covariance 
     *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
     *  @date   2024-07-22-
     */
    Ostap::StatusCode statCov 
    ( const RooAbsData*         data                           , 
      Ostap::Math::WCovariance& stat                           , 
      const std::string&        exp1                           , 
      const std::string&        exp2                           ,
      const std::string&        selection  = ""                ,
      const std::string&        cut_range  = ""                ,      
      const Ostap::EventIndex   first      = Ostap::FirstEvent ,
      const Ostap::EventIndex   last       = Ostap::LastEvent  , 
      const Ostap::DataType     xmin       = Ostap::MinValue   ,
      const Ostap::DataType     xmax       = Ostap::MaxValue   , 
      const Ostap::DataType     ymin       = Ostap::MinValue   ,
      const Ostap::DataType     ymax       = Ostap::MaxValue   ) const ;
    // =========================================================================
  public: // Covarinaces for many variables 
    // =========================================================================
    /** calculate the covariance of several expressions 
     *  @param data        (INPUT)  the inpout tree 
     *  @param stats       (UPDATE) the statistics  
     *  @param expressions (INPUT)  expressions 
     *  @param selection   (INPUT)  the selection criteria 
     *  @return status code 
     *  @attention selection is treate as boolean! 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2023-02-28
     */
    Ostap::StatusCode statCov
    ( TTree*                          data       , 
      Ostap::Math::Covariances&       stats      , 
      const Ostap::Strings&           exressions , 
      const std::string&              sekection  = ""                ,
      const Ostap::EventIndex         first      = Ostap::FirstEvent ,
      const Ostap::EventIndex         last       = Ostap::LastEvent  ) const ;
    // ========================================================================
    /** calculate the covariance of several expressions 
     *  @param data        (INPUT)  the inpout tree 
     *  @param stats       (UPDATE) the statistics 
     *  @param expressions (INPUT)  expressions 
     *  @param selection   (INPUT)  the selection criteria 
     *  @return status code 
     *  @attention selection is treate as weight!
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2023-02-28
     */
    Ostap::StatusCode statCov
    ( TTree*                          data        , 
      Ostap::Math::WCovariances&      stats       , 
      const Ostap::Strings&           expressions , 
      const std::string&              sekection   = ""                ,
      const Ostap::EventIndex         first       = Ostap::FirstEvent ,
      const Ostap::EventIndex         last        = Ostap::LastEvent  ) const ;          
    // =========================================================================
    /** calculate the covariance of several expressions 
     *  @param data        (INPUT)  the inpout tree 
     *  @param stats       (UPDATE) the statistics 
     *  @param cov2        (UPDATE) the covariance matrix 
     *  @param expressions (INPUT)  expressions 
     *  @param selection   (INPUT)  the selection criteria 
     *  @return status code 
     *  @attention selection is treate as weight!
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2023-02-28
     */
    Ostap::StatusCode statCov
    ( const RooAbsData*               data        , 
      Ostap::Math::WCovariances&      stats       , 
      const Ostap::Strings&           expressions , 
      const std::string&              sekection   = ""                ,
      const std::string&              cut_range   = ""                , 
      const Ostap::EventIndex         first       = Ostap::FirstEvent ,
      const Ostap::EventIndex         last        = Ostap::LastEvent  ) const ;          
    // ========================================================================
  public:     
    // ========================================================================
    /// the table column type 
    typedef std::vector<double> Column ;
    /// The table type
    typedef Ostap::Dict<Column> Table  ;
    // ========================================================================    
    /** Get data table from RooAbsData 
     *  @param data       (input)  data 
     *  @param table      (UPDATE) tale 
     *  @param selection  (INPUT)  selection/cut (treated as weight!)
     *  @param cut_range  (INPUT) if non empty: use events only from this cut-range
     *  @param first      (INPUT) the first event to process (inclusibe) 
     *  @param last       (INPUT) the last event to process (exclusive) 
     *  @return status code 
     */
    Ostap::StatusCode get_table
    ( const RooAbsData*         data                  ,
      Ostap::StatVar::Table&    table                 ,
      const std::string&        selection  = ""       ,
      const std::string&        cut_range  = ""       ,
      const Ostap::EventIndex   first      = Ostap::FirstEvent ,
      const Ostap::EventIndex   last       = Ostap::LastEvent  ) const ; 
    // ========================================================================
  public: 
    // ========================================================================
    /// congfiguration of the progress bar 
    inline const Ostap::Utils::ProgressConf& progress () const
    { return m_progress ; }
    // ========================================================================
  private :
    // ========================================================================
    /// congfiguration of the progress bar 
    Ostap::Utils::ProgressConf m_progress { false } ; 
    // ========================================================================
  } ;  
  // ==========================================================================


  
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_STATVAR_H
// ============================================================================
