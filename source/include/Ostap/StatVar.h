// ============================================================================
#ifndef OSTAP_STATVAR_H
#define OSTAP_STATVAR_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <limits>
// ============================================================================
// ROOT
// ============================================================================
#include "TMatrixDfwd.h"
// ============================================================================
// Forward declarations 
// =============================================================================
class TTree      ; // ROOT 
class TChain     ; // ROOT 
class RooAbsData ; // RooFit
// =============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatEntity.h"
#include "Ostap/WStatEntity.h"
#include "Ostap/ValueWithError.h"
#include "Ostap/Covariance.h"
#include "Ostap/StatusCode.h"
#include "Ostap/Types.h"
#include "Ostap/ProgressConf.h"
// ============================================================================
// forward declarations 
// ============================================================================
template <class SCALAR> class TMatrixTSym ; // ROOT 
// ============================================================================
namespace Ostap { namespace Math { class  Statistic   ; } }
namespace Ostap { namespace Math { class WStatistic   ; } }
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
  public: // helper structures 
    // ========================================================================
    /** @struct Interval 
     *  the actual type for interval
     */
    struct Interval
    {
      Interval
      ( const double l = 0 ,
	const double h = 0 ) ; 
      // ======================================================================
      /// low edge of the interval 
      double low  { 0 } ; // low edge of the interval 
      /// high edge of the interval 
      double high { 0 } ; // high edge of the interval 
      // ======================================================================
    };
    // ========================================================================
    /** @struct Quantile 
     *  The actual type for quantile and statitsics 
     */    
    struct Quantile 
    {
      Quantile
      ( const double      q = 0 ,
	const std::size_t n = 0 ) ;
      // ======================================================================
      /// quantile value  
      double      quantile { 0 } ; // quantile value  
      /// number of events used for estimation
      std::size_t nevents  { 0 } ; // number of events used for estimation
      // ======================================================================
    };
    // ========================================================================
    /** @struct Quantiles 
     *  The actual type for quantiles and statitsics 
     */    
    struct Quantiles 
    {
      // =======================================================================
      Quantiles
      ( const std::vector<double>& q = std::vector<double>() , 
	const std::size_t          n = 0  );
      // ======================================================================
      /// quantile values  
      std::vector<double> quantiles {} ; // quantile values  
      /// number of events used for estimation
      std::size_t         nevents  { 0 } ; // number of events used for estimation
      // ======================================================================
    };
    // ========================================================================
    /** @struct QInterval
     *  The actual type for interval with statitsics 
     */    
    struct QInterval    
    {
      // =====================================================================
      QInterval
      ( const Interval&     i = Interval () ,
	const std::size_t   n = 0           ) ;
      // ======================================================================
      /// the interval
      Interval      interval {}    ; // the interval
      /// number of events used for estimation
      std::size_t   nevents  { 0 } ; // number of events used for estimation
      // ======================================================================
    };    
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
     *  @param first      (INPIUT) the first event to process (inclusibe) 
     *  @param last       (INPIUT) the last event to process (exclusive) 
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
      const Ostap::EventIndex   last       = Ostap::LastEvent  ) const ;
    // ========================================================================
    /** Fill/update statistical counter 
     *  @param data       (input) data 
     *  @param stat       (UDATE) statistical counter 
     *  @param expression (INPUT) the variable 
     *  @param selection  (INPUT) selection/cut (treated as weight!)
     *  @param first      (INPIUT) the first event to process (inclusibe) 
     *  @param last       (INPIUT) the last event to process (exclusive) 
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
      const Ostap::EventIndex   last       = Ostap::LastEvent  ) const ;
    // ========================================================================
    /** Fill/update statistical counter 
     *  @param data       (input) data 
     *  @param stat       (UDATE) statistical counter 
     *  @param expression (INPUT) the variable 
     *  @param selection  (INPUT) selection/cut (treated as weight!)
     *  @param cut_range  (INPUT) if non empty: use evene only fomthis cut-range
     *  @param first      (INPIUT) the first event to process (inclusibe) 
     *  @param last       (INPIUT) the last event to process (exclusive) 
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
      const Ostap::EventIndex   last       = Ostap::LastEvent  ) const ;
   // ========================================================================    
  public:  // Statistics2 & WStatistics2 
    // ========================================================================
    /** Fill/update statistical counter 
     *  @param data       (input) data 
     *  @param stat       (UDATE) statistical counter 
     *  @param expr1      (INPUT) the 1st variable 
     *  @param expr1      (INPUT) the 2nd variable 
     *  @param selection  (INPUT) selection/cut (treated as boolean!)
     *  @param first      (INPIUT) the first event to process (inclusibe) 
     *  @param last       (INPIUT) the last event to process (exclusive) 
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
      const Ostap::EventIndex   last       = Ostap::LastEvent  ) const ;
    // ========================================================================
    /** Fill/update statistical counter 
     *  @param data       (input) data 
     *  @param stat       (UDATE) statistical counter 
     *  @param expr1      (INPUT) the 1st variable 
     *  @param expr1      (INPUT) the 2nd variable 
     *  @param selection  (INPUT) selection/cut (treated as weight!)
     *  @param first      (INPIUT) the first event to process (inclusibe) 
     *  @param last       (INPIUT) the last event to process (exclusive) 
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
      const Ostap::EventIndex   last       = Ostap::LastEvent  ) const ;
    // ========================================================================
    /** Fill/update statistical counter 
     *  @param data       (input) data 
     *  @param stat       (UDATE) statistical counter 
     *  @param expr1       (INPUT) the 1st variable 
     *  @param expr1      (INPUT) the 2nd variable 
     *  @param selection  (INPUT) selection/cut (treated as weight!)
     *  @param cut_range  (INPUT) if non empty: use evene only fomthis cut-range
     *  @param first      (INPIUT) the first event to process (inclusibe) 
     *  @param last       (INPIUT) the last event to process (exclusive) 
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
      const Ostap::EventIndex   last       = Ostap::LastEvent  ) const ;
    // ========================================================================





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
     *  @param first      (INPIUT) the first event to process (inclusibe) 
     *  @param last       (INPIUT) the last event to process (exclusive) 
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
      const Ostap::EventIndex   last       = Ostap::LastEvent  ) const ;
    // ========================================================================
    /** Fill/update statistical counter 
     *  @param data       (input) data 
     *  @param stat       (UDATE) statistical counter 
     *  @param expr1      (INPUT) the 1st variable 
     *  @param expr1      (INPUT) the 2nd variable 
     *  @param expr3      (INPUT) the 3rd variable  
     *  @param selection  (INPUT) selection/cut (treated as weight!)
     *  @param first      (INPIUT) the first event to process (inclusibe) 
     *  @param last       (INPIUT) the last event to process (exclusive) 
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
      const Ostap::EventIndex   last       = Ostap::LastEvent  ) const ;
    // ========================================================================
    /** Fill/update statistical counter 
     *  @param data       (input) data 
     *  @param stat       (UDATE) statistical counter 
     *  @param expr1      (INPUT) the 1st variable 
     *  @param expr1      (INPUT) the 2nd variable 
     *  @param expr3      (INPUT) the 3rd variable 
     *  @param selection  (INPUT) selection/cut (treated as weight!)
     *  @param cut_range  (INPUT) if non empty: use evene only fomthis cut-range
     *  @param first      (INPIUT) the first event to process (inclusibe) 
     *  @param last       (INPIUT) the last event to process (exclusive) 
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
      const Ostap::EventIndex   last       = Ostap::LastEvent  ) const ;
    // ======================================================================    
  public:  // Statistics3 & WStatistics3
    // ======================================================================
    /** Fill/update statistical counter 
     *  @param data       (input) data 
     *  @param stat       (UDATE) statistical counter 
     *  @param expr1      (INPUT) the 1st variable 
     *  @param expr1      (INPUT) the 2nd variable 
     *  @param expr3      (INPUT) the 3rd variable 
     *  @param expr4      (INPUT) the 4th variable   
     *  @param selection  (INPUT) selection/cut (treated as boolean!)
     *  @param first      (INPIUT) the first event to process (inclusibe) 
     *  @param last       (INPIUT) the last event to process (exclusive) 
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
      const Ostap::EventIndex   last       = Ostap::LastEvent  ) const ;
    // ========================================================================
    /** Fill/update statistical counter 
     *  @param data       (input) data 
     *  @param stat       (UDATE) statistical counter 
     *  @param expr1      (INPUT) the 1st variable 
     *  @param expr1      (INPUT) the 2nd variable 
     *  @param expr3      (INPUT) the 3rd variable  
     *  @param expr4      (INPUT) the 4th variable 
     *  @param selection  (INPUT) selection/cut (treated as weight!)
     *  @param first      (INPIUT) the first event to process (inclusibe) 
     *  @param last       (INPIUT) the last event to process (exclusive) 
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
      const Ostap::EventIndex   last       = Ostap::LastEvent  ) const ;
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
     *  @param first      (INPIUT) the first event to process (inclusibe) 
     *  @param last       (INPIUT) the last event to process (exclusive) 
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
      const Ostap::EventIndex   last       = Ostap::LastEvent  ) const ;
    // ========================================================================
  public: /// get infomation about several variables simultabneously 
    // =========================================================================
    /**  build statistic for the <code>expressions</code>
     *  @param data        (INPUT)  input data 
     *  @param result      (UPDATE) the output statistics for specified expressions 
     *  @param expressions (INPUT)  the list of  expressions
     *  @param selection   (INPUT)  selection 
     *  @param first       (INPUT)  the first entry to process 
     *  @param last        (INPUT)  the last entry to process (not including!)
     *  @return statsu code 
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
  public : //  
    // ========================================================================
    /** build statistic for the <code>expression</code>
     *  @param tree       (INPUT) the tree 
     *  @param expression (INPUT) the expression
     *  @param first      (INPUT) the first entry 
     *  @param last       (INPUT) the last entry      *
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
      const Ostap::EventIndex last       = Ostap::LastEvent  ) const ;
    // =================================================================
    /** build statistic for the <code>expression</code>
     *  @param tree       (INPUT) the tree 
     *  @param expression (INPUT) the expression
     *  @param selection  (INPUT) selection/ (as boolean) 
     *  @param first      (INPUT) the first entry 
     *  @param last       (INPUT) the last entry      *
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
      const Ostap::EventIndex last       = Ostap::LastEvent  ) const ;
    // =================================================================
    /** build statistic for the <code>expression</code>
     *  @param tree       (INPUT) the tree 
     *  @param expression (INPUT) the expression
     *  @param selection  (INPUT) selection/weight 
     *  @param first      (INPUT) the first entry 
     *  @param last       (INPUT) the last entry      *
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
      const Ostap::EventIndex last       = Ostap::LastEvent  ) const ;
    // ===============================================================
    /** build statistic for the <code>expression</code>
     *  @param tree       (INPUT) the tree 
     *  @param expression (INPUT) the expression
     *  @param selection  (INPUT) seleciton/weight 
     *  @param cut_range  (INPUT) cut-range 
     *  @param first      (INPUT) the first entry 
     *  @param last       (INPUT) the last entry      *
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
      const Ostap::EventIndex last       = Ostap::LastEvent  ) const ;
    // =========================================================================
  public : // covariances
    // =========================================================================
    /** calculate the covariance of two expressions 
     *  @param tree  (INPUT)  the inpout tree 
     *  @param exp1  (INPUT)  the first  expression
     *  @param exp2  (INPUT)  the second expression
     *  @param first (INPUT)  the first  event to process (inclusive) 
     *  @param last  (INPUT)  the last   event to process (exc;usive) 
     *  @return Covariance 
     *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
     *  @date   2024-07-22-
     */
    Ostap::Math::Covariance statCov 
    ( TTree*                  tree                      , 
      const std::string&      exp1                      , 
      const std::string&      exp2                      , 
      const Ostap::EventIndex first = Ostap::FirstEvent ,
      const Ostap::EventIndex last  = Ostap::LastEvent  ) const ;
    // ======================================================================
    /** calculate the covariance of two expressions 
     *  @param tree  (INPUT)  the inpout tree 
     *  @param exp1  (INPUT)  the first  expression
     *  @param exp2  (INPUT)  the second expression
     *  @param selection 
     *  @param first (INPUT)  the first  event to process (inclusive) 
     *  @param last  (INPUT)  the last   event to process (exc;usive) 
     *  @return Covariance 
     *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
     *  @date   2024-07-22-
     */
    Ostap::Math::Covariance statCov_cut 
    ( TTree*                  tree                          , 
      const std::string&      exp1                          , 
      const std::string&      exp2                          , 
      const std::string&      selection = ""                , 
      const Ostap::EventIndex first     = Ostap::FirstEvent ,
      const Ostap::EventIndex last      = Ostap::LastEvent  ) const ;
    // ======================================================================    
    /** calculate the covariance of two expressions 
     *  @param tree  (INPUT)  the inpout tree 
     *  @param exp1  (INPUT)  the first  expression
     *  @param exp2  (INPUT)  the second expression
     *  @param selection 
     *  @param first (INPUT)  the first  event to process (inclusive) 
     *  @param last  (INPUT)  the last   event to process (exc;usive) 
     *  @return Covariance 
     *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
     *  @date   2024-07-22-
     */
    Ostap::Math::WCovariance statCov 
    ( TTree*                  tree       , 
      const std::string&      exp1       , 
      const std::string&      exp2       ,
      const std::string&      selection  , 
      const Ostap::EventIndex first      = Ostap::FirstEvent ,
      const Ostap::EventIndex last       = Ostap::LastEvent  ) const ;
    // ========================================================================
    /** calculate the covariance of two expressions 
     *  @param tree  (INPUT)  the inpout tree 
     *  @param exp1  (INPUT)  the first  expression
     *  @param exp2  (INPUT)  the second expression
     *  @param selection 
     *  @param first (INPUT)  the first  event to process (inclusive) 
     *  @param last  (INPUT)  the last   event to process (exc;usive) 
     *  @return Covariance 
     *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
     *  @date   2024-07-22-
     */
    Ostap::Math::WCovariance statCov 
    ( const RooAbsData*       data                           , 
      const std::string&      exp1                           , 
      const std::string&      exp2                           ,
      const std::string&      selection  = ""                ,
      const std::string&      cut_range  = ""                ,      
      const Ostap::EventIndex first      = Ostap::FirstEvent ,
      const Ostap::EventIndex last       = Ostap::LastEvent  ) const ;    
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
  private:
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
