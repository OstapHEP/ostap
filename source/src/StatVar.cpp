// ============================================================================
// Include files
// ============================================================================
//   STD&STL
// ============================================================================
#include <cmath>
#include <numeric>
#include <algorithm>
#include <set>
#include <random>
#include <tuple>
// ============================================================================
// ROOT
// ============================================================================
#include "TTree.h"
#include "TMatrixTSym.h"
#include "RooDataSet.h"
#include "RooAbsReal.h"
#include "RooFormulaVar.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Names.h"
#include "Ostap/Combiner.h"
#include "Ostap/Formula.h"
#include "Ostap/Notifier.h"
#include "Ostap/StatVar.h"
#include "Ostap/FormulaVar.h"
#include "Ostap/Covariances.h"
#include "Ostap/ProgressBar.h"
// ============================================================================
// Local
// ============================================================================
#include "Formulae.h"
#include "Exception.h"
#include "status_codes.h"
// ============================================================================
/** @file
 *  Implementation file for class Ostap::StatVar
 *  @date 2013-10-13
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 */
// ============================================================================
namespace
{
  // ==========================================================================
  /// make FormulaVar 
  std::unique_ptr<Ostap::FormulaVar>
  make_formula 
  ( const std::string& expression           , 
    const RooAbsData*  data                 , 
    const bool         allow_empty = false  , 
    const bool         allow_null  = false  ) 
  { 
    if ( allow_empty && Ostap::trivial ( expression ) ) { return nullptr ; }  // RETURN!
    //
    const RooArgSet*  aset = data -> get() ;
    if ( allow_null && nullptr == aset        ) { return nullptr ; }
    //
    Ostap::Assert ( nullptr != aset                ,  
                    "Invalid varset"               , 
                    "Ostap::StatVar::make_formula" ,
		    INVALID_ARGSET , __FILE__ , __LINE__ ) ;
    //
    RooArgList        alst { *aset } ;
    auto result = std::make_unique<Ostap::FormulaVar> ( expression , alst , false ) ;
    if ( allow_null && ( !result || !result->ok() ) ) { return nullptr ; } 
    //
    Ostap::Assert ( result && result->ok()                   , 
                    "Invalid formula:\"" + expression + "\"" , 
                    "Ostap::StatVar::make_formula"           ,
		    INVALID_FORMULA , __FILE__ , __LINE__    ) ;
    //
    return result ;
  }
  // ==========================================================================
  /// make FormulaVar 
  std::unique_ptr<Ostap::Formula>
  make_formula 
  ( const std::string& expression           , 
    const TTree*       data                 , 
    const bool         allow_empty = false  , 
    const bool         allow_null  = false  ) 
  { 
    if ( allow_empty && Ostap::trivial ( expression ) ) { return nullptr ; }  // RETURN!
    // 
    auto result = std::make_unique<Ostap::Formula> ( expression , data  ) ;
    if ( allow_null && ( !result || !result->ok() ) ) { return nullptr ; } 
    //
    Ostap::Assert ( result && result->ok()                   , 
                    "Invalid formula:\"" + expression + "\"" , 
                    "Ostap::StatVar::make_formula"           ,
		    INVALID_FORMULA , __FILE__ , __LINE__    ) ;
    //
    return result ;
  }
  // ==========================================================================
} //                                                 end of anonymous namespace
//  ============================================================================
// constructor with ProgressBar configuration 
// =============================================================================
Ostap::StatVar::StatVar
( const Ostap::Utils::ProgressConf& progress)
  : m_progress ( progress )
{}
// ============================================================================


// =============================================================================
// The basic methods 
// =============================================================================
/* Fill/update statistical counter 
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
// =============================================================================
Ostap::StatusCode Ostap::StatVar::get_stat
( TTree*                    data       ,
  Ostap::Math::Statistic&   stat       ,
  const std::string&        expression , 
  const std::string&        selection  ,
  const Ostap::EventIndex   first      ,
  const Ostap::EventIndex   last       ) const 
{
  /// reset the statistic
  stat.reset() ; // reset the statistic
  //
  if ( nullptr == data    ) { return INVALID_DATA    ; }
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  //
  const Ostap::EventIndex  num_entries = data -> GetEntries () ;
  const Ostap::EventIndex  the_last    = std::min ( last , num_entries ) ;
  if ( the_last <= first ) { return Ostap::StatusCode::SUCCESS ; }
  //
  Ostap::Formula formula ( expression , data ) ;
  if ( !formula.ok ()     ) { return INVALID_FORMULA ; }
  //
  const std::unique_ptr<Ostap::Formula> cuts { make_formula ( selection , data , true ) } ;
  const bool                with_cuts   = cuts && cuts->ok () ;
  //
  Ostap::Utils::ProgressBar bar    ( the_last - first , m_progress ) ;
  Ostap::Utils::Notifier    notify ( data , &formula , cuts.get () ) ;
  //
  std::vector<double>       results {} ;
  for ( Ostap::EventIndex   entry = first ; entry < the_last  ; ++entry , ++bar )
    {
      //
      const long ievent = data -> GetEntryNumber ( entry ) ;
      if ( 0 > ievent )                      { return INVALID_ENTRY ; }  // RETURN 
      if ( 0 > data -> LoadTree ( ievent ) ) { return INVALID_EVENT ; }  // RETURN 
      //
      const double  cut = with_cuts ? cuts->evaluate() : 1.0 ;
      if ( !cut ) { continue ; }                                         // ATTENTION! 
      //
      formula.evaluate ( results ) ;
      for  ( const double value : results ) { stat.update ( value ) ; }
    }
  //
  return Ostap::StatusCode::SUCCESS ;
}
// =========================================================================
/*  Fill/update statistical counter 
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
// =============================================================================
Ostap::StatusCode Ostap::StatVar::get_stat
( TTree*                    data       ,
  Ostap::Math::WStatistic&  stat       ,
  const std::string&        expression , 
  const std::string&        selection  ,
  const Ostap::EventIndex   first      ,
  const Ostap::EventIndex   last       ) const
{
  /// reset the statistic
  stat.reset() ; // reset the statistic
  //
  if ( nullptr == data    ) { return INVALID_DATA    ; }
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  //
  const Ostap::EventIndex   num_entries = data -> GetEntries()  ;
  const Ostap::EventIndex   the_last    = std::min ( last , num_entries ) ;
  if ( the_last <= first   ) { return Ostap::StatusCode::SUCCESS ; }
  //
  Ostap::Formula formula ( expression , data ) ;
  if ( !formula.ok()      ) { return INVALID_FORMULA ; }
  //
  const std::unique_ptr<Ostap::Formula> cuts { make_formula ( expression , data , true ) } ; 
  const bool                with_cuts   = cuts && cuts->ok () ;
  //
  Ostap::Utils::Notifier    notify ( data , &formula  , cuts.get() ) ;
  Ostap::Utils::ProgressBar bar    ( the_last - first , m_progress ) ;
  //
  std::vector<double>       results {} ;
  for ( Ostap::EventIndex   entry = first ; entry < the_last ; ++entry , ++bar )
    {
      //
      const long ievent = data -> GetEntryNumber ( entry ) ;
      if ( 0 > ievent                      ) { return INVALID_ENTRY ; } // RETURN 
      if ( 0 > data -> LoadTree ( ievent ) ) { return INVALID_EVENT ; } // RETURN 
      //
      const double weight = with_cuts ? cuts->evaluate() : 1.0 ;
      if ( !weight ) { continue ; }                                    // ATTENTION! 
      //
      formula.evaluate ( results ) ;
      for  ( const double value : results ) { stat.update ( value , weight ) ; }
    }
  //
  return Ostap::StatusCode::SUCCESS ;
}
// =============================================================================
/*  Fill/update statistical counter 
 *  @param data       (input) data 
 *  @param stat       (UDATE) statistical counter 
 *  @param expression (INPUT) the variable 
 *  @param selection  (INPUT) selection/cut (treated as weight!)
 *  @param cut_range  (INPUT) if non empty: use evene only fomthis cut-range
 *  @param first      (INPIUT) the first event to process (inclusibe) 
 *  @param last       (INPIUT) the last event to process (exclusive) 
 *  @return status code 
 *  @see Ostap::Math::WStatistic
 *  @attention selection/cut is treated as weight!
 */
// =============================================================================
Ostap::StatusCode Ostap::StatVar::get_stat
( const RooAbsData*         data       ,
  Ostap::Math::WStatistic&  stat       ,
  const std::string&        expression , 
  const std::string&        selection  ,
  const std::string&        cut_range  ,
  const Ostap::EventIndex   first      ,
  const Ostap::EventIndex   last       ) const
{
  /// reset the statistic
  stat.reset() ; // reset the statistic
  // 
  if ( nullptr == data    ) { return INVALID_DATA    ; }
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  //
  const Ostap::EventIndex num_entries = data -> numEntries() ;
  const Ostap::EventIndex the_last    = std::min ( num_entries , last ) ;
  if ( the_last <= first  ) { return Ostap::StatusCode::SUCCESS ; }
  //
  const bool  weighted = data -> isWeighted() ;
  const char* cutrange = cut_range.empty() ? nullptr : cut_range.c_str() ;
  const std::unique_ptr<Ostap::FormulaVar> expr { make_formula ( expression , data       ) } ;
  if ( !expr || !expr->ok()            ) { return INVALID_FORMULA ; }
  //
  const std::unique_ptr<Ostap::FormulaVar> cuts { make_formula ( selection  , data , true ) } ;
  const bool                with_cuts   = cuts && cuts->ok () ;
  //
  Ostap::Utils::ProgressBar bar { the_last - first , m_progress } ; 
  for ( Ostap::EventIndex entry = first ; entry < the_last ; ++entry , ++bar )
    {
      const RooArgSet* vars = data -> get ( entry ) ;
      if ( nullptr == vars )                             { break    ; } // BREAK 
      //
      if ( cutrange && !vars->allInRange ( cutrange ) )  { continue ; } // CONTINUE    
      // apply weight:
      const double wd = weighted  ? data -> weight () : 1.0 ;
      if ( !wd ) { continue ; }                                   // CONTINUE    
      // apply cuts:
      const double wc = with_cuts ? cuts -> getVal () : 1.0 ;
      if ( !wc ) { continue ; }                                   // CONTINUE  
      // Total: cuts & weight:
      const double weight  = wd *  wc ;
      if ( !weight  ) { continue ; }                             // CONTINUE        
      //
      const double value = expr -> getVal () ;
      //
      stat.update ( value , weight ) ;
    }
  //
  return Ostap::StatusCode::SUCCESS ; 
}
// =============================================================================

// =============================================================================
// Get information about several varibales 
// =============================================================================
/* Fill/update statistical counter 
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
// =============================================================================
Ostap::StatusCode Ostap::StatVar::get_stat
( TTree*                    data       ,
  Ostap::Math::Statistic2&  stat       ,
  const std::string&        expr1      ,
  const std::string&        expr2      ,  
  const std::string&        selection  , 
  const Ostap::EventIndex   first      ,
  const Ostap::EventIndex   last       ) const 
{
  /// reset the statistic
  stat.reset() ; // reset the statistic
  // 
  if ( nullptr == data    ) { return INVALID_DATA    ; }
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  //
  const Ostap::EventIndex  num_entries = data -> GetEntries () ;
  const Ostap::EventIndex  the_last    = std::min ( last , num_entries ) ;
  if ( the_last <= first ) { return Ostap::StatusCode::SUCCESS ; }
  //
 /// formulae for exressons
  const Ostap::Formulae                 formulae { data , { expr1 , expr2 } } ;
  const std::unique_ptr<Ostap::Formula> cuts     { make_formula ( selection , data , true ) } ;
  const bool                with_cuts = cuts && cuts->ok() ;
  //
  Ostap::Utils::Notifier    notify  ( formulae.begin() , formulae.end() , cuts.get() , data ) ;
  //
  Ostap::Utils::ProgressBar bar    ( the_last - first , m_progress ) ;
  //
  typedef std::vector<double>  RESULT  ;
  typedef std::array<RESULT,2> RESULTS ; 
  RESULTS results ;
  ///  
  for ( Ostap::EventIndex   entry = first ; entry < the_last  ; ++entry , ++bar )
    {
      //
      const long ievent = data -> GetEntryNumber ( entry ) ;
      if ( 0 > ievent )                      { return INVALID_ENTRY ; }  // RETURN 
      if ( 0 > data -> LoadTree ( ievent ) ) { return INVALID_EVENT ; }  // RETURN 
      //
      const double  cut = with_cuts ? cuts->evaluate() : 1.0 ;
      if ( !cut ) { continue ; }                                         // ATTENTION! 
      //
      formulae.evaluate ( 0 , results [ 0 ] ) ;
      formulae.evaluate ( 1 , results [ 1 ] ) ;
      //
      for ( const double x : results [ 0 ] )
      { for ( const double y : results [ 1 ] )
          { stat.update ( x , y ) ; }}
    }
  //
  return Ostap::StatusCode::SUCCESS ;
}
// =========================================================================
/*  Fill/update statistical counter 
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
// =============================================================================
Ostap::StatusCode Ostap::StatVar::get_stat
( TTree*                    data       ,
  Ostap::Math::WStatistic2& stat       ,
  const std::string&        expr1      ,
  const std::string&        expr2      ,  
  const std::string&        selection  ,
  const Ostap::EventIndex   first      ,
  const Ostap::EventIndex   last       ) const
{
  /// reset the statistic
  stat.reset() ; // reset the statistic
  // 
  if ( nullptr == data    ) { return INVALID_DATA    ; }
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  //
  const Ostap::EventIndex   num_entries = data -> GetEntries()  ;
  const Ostap::EventIndex   the_last    = std::min ( last , num_entries ) ;
  if ( the_last <= first   ) { return Ostap::StatusCode::SUCCESS ; }
  //
  /// formulae for exressons
  const Ostap::Formulae                 formulae { data , { expr1 , expr2 } } ;
  const std::unique_ptr<Ostap::Formula> cuts     { make_formula ( selection , data , true ) } ;
  const bool                with_cuts = cuts && cuts->ok() ;
  //
  Ostap::Utils::Notifier    notify ( formulae.begin() , formulae.end()  , cuts.get() , data ) ;
  Ostap::Utils::ProgressBar bar    ( the_last - first , m_progress ) ;
  //
  typedef std::vector<double>  RESULT  ;
  typedef std::array<RESULT,2> RESULTS ; 
  RESULTS results ;
  //
  for ( Ostap::EventIndex   entry = first ; entry < the_last ; ++entry , ++bar )
    {
      //
      const long ievent = data -> GetEntryNumber ( entry ) ;
      if ( 0 > ievent                      ) { return INVALID_ENTRY ; } // RETURN 
      if ( 0 > data -> LoadTree ( ievent ) ) { return INVALID_EVENT ; } // RETURN 
      //
      const double weight = with_cuts ? cuts->evaluate() : 1.0 ;
      if ( !weight ) { continue ; }                                    // ATTENTION! 
      //
      formulae.evaluate ( 0 , results [ 0 ] ) ;
      formulae.evaluate ( 1 , results [ 1 ] ) ;
      //
      for ( const double x : results [ 0 ] )
      { for ( const double y : results [ 1 ] )
          { stat.update ( x , y , weight ) ; }}
    }
  //
  return Ostap::StatusCode::SUCCESS ;
}
// =============================================================================
/*  Fill/update statistical counter 
 *  @param data       (input) data 
 *  @param stat       (UDATE) statistical counter 
 *  @param expression (INPUT) the variable 
 *  @param selection  (INPUT) selection/cut (treated as weight!)
 *  @param cut_range  (INPUT) if non empty: use evene only fomthis cut-range
 *  @param first      (INPIUT) the first event to process (inclusibe) 
 *  @param last       (INPIUT) the last event to process (exclusive) 
 *  @return status code 
 *  @see Ostap::Math::WStatistic
 *  @attention selection/cut is treated as weight!
 */
// =============================================================================
Ostap::StatusCode Ostap::StatVar::get_stat
( const RooAbsData*         data       ,
  Ostap::Math::WStatistic2& stat       ,
  const std::string&        expr1      ,
  const std::string&        expr2      ,  
  const std::string&        selection  , 
  const std::string&        cut_range  ,
  const Ostap::EventIndex   first      ,
  const Ostap::EventIndex   last       ) const
{
  /// reset the statistic
  stat.reset() ; // reset the statistic
  //
  if ( nullptr == data    ) { return INVALID_DATA    ; }
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  //
  const Ostap::EventIndex num_entries = data -> numEntries() ;
  const Ostap::EventIndex the_last    = std::min ( num_entries , last ) ;
  if ( the_last <= first  ) { return Ostap::StatusCode::SUCCESS ; }
  //
  const bool  weighted = data -> isWeighted() ;
  const char* cutrange = cut_range.empty() ? nullptr : cut_range.c_str() ;
  //
  // formulae for exressons
  const Ostap::FormulaVars                 formulae { data , { expr1 , expr2 } } ;
  const std::unique_ptr<Ostap::FormulaVar> cuts { make_formula ( selection  , data , true ) } ;
  const bool                with_cuts   = cuts && cuts->ok () ;
  //
  Ostap::Utils::ProgressBar bar { the_last - first , m_progress } ; 
  for ( Ostap::EventIndex entry = first ; entry < the_last ; ++entry , ++bar )
    {
      const RooArgSet* vars = data -> get ( entry ) ;
      if ( nullptr == vars )                             { break    ; } // BREAK 
      //
      if ( cutrange && !vars->allInRange ( cutrange ) )  { continue ; } // CONTINUE    
      // apply weight:
      const double wd = weighted  ? data -> weight () : 1.0 ;
      if ( !wd ) { continue ; }                                   // CONTINUE    
      // apply cuts:
      const double wc = with_cuts ? cuts -> getVal () : 1.0 ;
      if ( !wc ) { continue ; }                                   // CONTINUE  
      // Total: cuts & weight:
      const double weight  = wd *  wc ;
      if ( !weight  ) { continue ; }                             // CONTINUE        
      //
      const double x = formulae.evaluate ( 0 ) ; 
      const double y = formulae.evaluate ( 1 ) ; 
      //
      stat.update ( x , y  , weight ) ;
    }
  //
  return Ostap::StatusCode::SUCCESS ; 
}
// =============================================================================

// ============================================================================
// Get information about several varibales 
// =============================================================================
/* Fill/update statistical counter 
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
// =============================================================================
Ostap::StatusCode Ostap::StatVar::get_stat
( TTree*                    data       ,
  Ostap::Math::Statistic3&  stat       ,
  const std::string&        expr1      ,
  const std::string&        expr2      ,
  const std::string&        expr3      ,   
  const std::string&        selection  , 
  const Ostap::EventIndex   first      ,
  const Ostap::EventIndex   last       ) const 
{
  /// reset the statistic
  stat.reset() ; // reset the statistic
  // 
  if ( nullptr == data    ) { return INVALID_DATA    ; }
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  //
  const Ostap::EventIndex  num_entries = data -> GetEntries () ;
  const Ostap::EventIndex  the_last    = std::min ( last , num_entries ) ;
  if ( the_last <= first ) { return Ostap::StatusCode::SUCCESS ; }
  //
 /// formulae for exressons
  const Ostap::Formulae                 formulae { data , { expr1 , expr2 , expr3 } } ;
  const std::unique_ptr<Ostap::Formula> cuts     { make_formula ( selection , data , true ) } ;
  const bool                with_cuts = cuts && cuts->ok() ;
  //
  Ostap::Utils::Notifier    notify  ( formulae.begin() , formulae.end() , cuts.get() , data ) ;
  //
  Ostap::Utils::ProgressBar bar    ( the_last - first , m_progress ) ;
  //
  typedef std::vector<double>  RESULT  ;
  typedef std::array<RESULT,3> RESULTS ; 
  RESULTS results ;
  ///  
  for ( Ostap::EventIndex   entry = first ; entry < the_last  ; ++entry , ++bar )
    {
      //
      const long ievent = data -> GetEntryNumber ( entry ) ;
      if ( 0 > ievent )                      { return INVALID_ENTRY ; }  // RETURN 
      if ( 0 > data -> LoadTree ( ievent ) ) { return INVALID_EVENT ; }  // RETURN 
      //
      const double  cut = with_cuts ? cuts->evaluate() : 1.0 ;
      if ( !cut ) { continue ; }                                         // ATTENTION! 
      //
      formulae.evaluate ( 0 , results [ 0 ] ) ;
      formulae.evaluate ( 1 , results [ 1 ] ) ;
      formulae.evaluate ( 2 , results [ 2 ] ) ;
      //
      for ( const double x : results [ 0 ] )
      { for ( const double y : results [ 1 ] )
        { for ( const double z : results [ 2 ] )
          { stat.update ( x , y , z ) ; } } }
    }
  //
  return Ostap::StatusCode::SUCCESS ;
}
// =========================================================================
/*  Fill/update statistical counter 
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
// =============================================================================
Ostap::StatusCode Ostap::StatVar::get_stat
( TTree*                    data       ,
  Ostap::Math::WStatistic3& stat       ,
  const std::string&        expr1      ,
  const std::string&        expr2      ,  
  const std::string&        expr3      ,
  const std::string&        selection  ,
  const Ostap::EventIndex   first      ,
  const Ostap::EventIndex   last       ) const
{
  /// reset the statistic
  stat.reset() ; // reset the statistic
  // 
  if ( nullptr == data    ) { return INVALID_DATA    ; }
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  //
  const Ostap::EventIndex   num_entries = data -> GetEntries()  ;
  const Ostap::EventIndex   the_last    = std::min ( last , num_entries ) ;
  if ( the_last <= first   ) { return Ostap::StatusCode::SUCCESS ; }
  //
  /// formulae for exressons
  const Ostap::Formulae                 formulae { data , { expr1 , expr2 , expr3 } } ;
  const std::unique_ptr<Ostap::Formula> cuts     { make_formula ( selection , data , true ) } ;
  const bool                with_cuts = cuts && cuts->ok() ;
  //
  Ostap::Utils::Notifier    notify ( formulae.begin() , formulae.end()  , cuts.get() , data ) ;
  Ostap::Utils::ProgressBar bar    ( the_last - first , m_progress ) ;
  //
  typedef std::vector<double>  RESULT  ;
  typedef std::array<RESULT,3> RESULTS ; 
  RESULTS results ;
  //
  for ( Ostap::EventIndex   entry = first ; entry < the_last ; ++entry , ++bar )
    {
      //
      const long ievent = data -> GetEntryNumber ( entry ) ;
      if ( 0 > ievent                      ) { return INVALID_ENTRY ; } // RETURN 
      if ( 0 > data -> LoadTree ( ievent ) ) { return INVALID_EVENT ; } // RETURN 
      //
      const double weight = with_cuts ? cuts->evaluate() : 1.0 ;
      if ( !weight ) { continue ; }                                    // ATTENTION! 
      //
      formulae.evaluate ( 0 , results [ 0 ] ) ;
      formulae.evaluate ( 1 , results [ 1 ] ) ;
      formulae.evaluate ( 2 , results [ 2 ] ) ;
      //
      for ( const double x : results [ 0 ] )
      { for ( const double y : results [ 1 ] )
        { for ( const double z : results [ 2 ] )
          { stat.update ( x , y , z , weight ) ; } } }
    }
  //
  return Ostap::StatusCode::SUCCESS ;
}
// =============================================================================
/*  Fill/update statistical counter 
 *  @param data       (input) data 
 *  @param stat       (UDATE) statistical counter 
 *  @param expression (INPUT) the variable 
 *  @param selection  (INPUT) selection/cut (treated as weight!)
 *  @param cut_range  (INPUT) if non empty: use evene only fomthis cut-range
 *  @param first      (INPIUT) the first event to process (inclusibe) 
 *  @param last       (INPIUT) the last event to process (exclusive) 
 *  @return status code 
 *  @see Ostap::Math::WStatistic
 *  @attention selection/cut is treated as weight!
 */
// =============================================================================
Ostap::StatusCode Ostap::StatVar::get_stat
( const RooAbsData*         data       ,
  Ostap::Math::WStatistic3& stat       ,
  const std::string&        expr1      ,
  const std::string&        expr2      ,  
  const std::string&        expr3      , 
  const std::string&        selection  , 
  const std::string&        cut_range  ,
  const Ostap::EventIndex   first      ,
  const Ostap::EventIndex   last       ) const
{
  /// reset the statistic
  stat.reset() ; // reset the statistic
  // 
  if ( nullptr == data    ) { return INVALID_DATA    ; }
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  //
  const Ostap::EventIndex num_entries = data -> numEntries() ;
  const Ostap::EventIndex the_last    = std::min ( num_entries , last ) ;
  if ( the_last <= first  ) { return Ostap::StatusCode::SUCCESS ; }
  //
  const bool  weighted = data -> isWeighted() ;
  const char* cutrange = cut_range.empty() ? nullptr : cut_range.c_str() ;
  //
  // formulae for exressons
  const Ostap::FormulaVars                 formulae { data , { expr1 , expr2 , expr3} } ;
  const std::unique_ptr<Ostap::FormulaVar> cuts { make_formula ( selection  , data , true ) } ;
  const bool                with_cuts   = cuts && cuts->ok () ;
  //
  Ostap::Utils::ProgressBar bar { the_last - first , m_progress } ; 
  for ( Ostap::EventIndex entry = first ; entry < the_last ; ++entry , ++bar )
    {
      const RooArgSet* vars = data -> get ( entry ) ;
      if ( nullptr == vars )                             { break    ; } // BREAK 
      //
      if ( cutrange && !vars->allInRange ( cutrange ) )  { continue ; } // CONTINUE    
      // apply weight:
      const double wd = weighted  ? data -> weight () : 1.0 ;
      if ( !wd ) { continue ; }                                   // CONTINUE    
      // apply cuts:
      const double wc = with_cuts ? cuts -> getVal () : 1.0 ;
      if ( !wc ) { continue ; }                                   // CONTINUE  
      // Total: cuts & weight:
      const double weight  = wd *  wc ;
      if ( !weight  ) { continue ; }                             // CONTINUE        
      //
      const double x = formulae.evaluate ( 0 ) ; 
      const double y = formulae.evaluate ( 1 ) ; 
      const double z = formulae.evaluate ( 2 ) ; 
      //
      stat.update ( x , y , z , weight ) ;
    }
  //
  return Ostap::StatusCode::SUCCESS ; 
}
// =============================================================================

// ============================================================================
// Get information about several varibales 
// =============================================================================
/* Fill/update statistical counter 
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
// =============================================================================
Ostap::StatusCode Ostap::StatVar::get_stat
( TTree*                    data       ,
  Ostap::Math::Statistic4&  stat       ,
  const std::string&        expr1      ,
  const std::string&        expr2      ,
  const std::string&        expr3      ,   
  const std::string&        expr4      , 
  const std::string&        selection  , 
  const Ostap::EventIndex   first      ,
  const Ostap::EventIndex   last       ) const 
{
  /// reset the statistic
  stat.reset() ; // reset the statistic
  // 
  if ( nullptr == data    ) { return INVALID_DATA    ; }
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  //
  const Ostap::EventIndex  num_entries = data -> GetEntries () ;
  const Ostap::EventIndex  the_last    = std::min ( last , num_entries ) ;
  if ( the_last <= first ) { return Ostap::StatusCode::SUCCESS ; }
  //
 /// formulae for exressons
  const Ostap::Formulae                 formulae { data , { expr1 , expr2 , expr3 , expr4 } } ; 
  const std::unique_ptr<Ostap::Formula> cuts     { make_formula ( selection , data , true ) } ;
  const bool                with_cuts = cuts && cuts->ok() ;
  //
  Ostap::Utils::Notifier    notify  ( formulae.begin() , formulae.end() , cuts.get() , data ) ;
  //
  Ostap::Utils::ProgressBar bar    ( the_last - first , m_progress ) ;
  //
  typedef std::vector<double>  RESULT  ;
  typedef std::array<RESULT,3> RESULTS ; 
  RESULTS results ;
  ///  
  for ( Ostap::EventIndex   entry = first ; entry < the_last  ; ++entry , ++bar )
    {
      //
      const long ievent = data -> GetEntryNumber ( entry ) ;
      if ( 0 > ievent )                      { return INVALID_ENTRY ; }  // RETURN 
      if ( 0 > data -> LoadTree ( ievent ) ) { return INVALID_EVENT ; }  // RETURN 
      //
      const double  cut = with_cuts ? cuts->evaluate() : 1.0 ;
      if ( !cut ) { continue ; }                                         // ATTENTION! 
      //
      formulae.evaluate ( 0 , results [ 0 ] ) ;
      formulae.evaluate ( 1 , results [ 1 ] ) ;
      formulae.evaluate ( 2 , results [ 2 ] ) ;
      formulae.evaluate ( 3 , results [ 3 ] ) ;
      //
      for ( const double x : results [ 0 ] )
      { for ( const double y : results [ 1 ] )
        { for ( const double z : results [ 2 ] )
          { for  ( const double t :  results [ 3 ] ) 
            { stat.update ( x , y , z , t ) ; } } } } 
    }
  //
  return Ostap::StatusCode::SUCCESS ;
}
// =========================================================================
/*  Fill/update statistical counter 
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
// =============================================================================
Ostap::StatusCode Ostap::StatVar::get_stat
( TTree*                    data       ,
  Ostap::Math::WStatistic4& stat       ,
  const std::string&        expr1      ,
  const std::string&        expr2      ,  
  const std::string&        expr3      ,
  const std::string&        expr4      , 
  const std::string&        selection  ,
  const Ostap::EventIndex   first      ,
  const Ostap::EventIndex   last       ) const
{
  /// reset the statistic
  stat.reset() ; // reset the statistic
  // 
  if ( nullptr == data    ) { return INVALID_DATA    ; }
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  //
  const Ostap::EventIndex   num_entries = data -> GetEntries()  ;
  const Ostap::EventIndex   the_last    = std::min ( last , num_entries ) ;
  if ( the_last <= first   ) { return Ostap::StatusCode::SUCCESS ; }
  //
  /// formulae for exressons
  const Ostap::Formulae                 formulae { data , { expr1 , expr2 , expr3 , expr4 } } ; 
  const std::unique_ptr<Ostap::Formula> cuts     { make_formula ( selection , data , true ) } ;
  const bool                with_cuts = cuts && cuts->ok() ;
  //
  Ostap::Utils::Notifier    notify ( formulae.begin() , formulae.end()  , cuts.get() , data ) ;
  Ostap::Utils::ProgressBar bar    ( the_last - first , m_progress ) ;
  //
  typedef std::vector<double>  RESULT  ;
  typedef std::array<RESULT,4> RESULTS ; 
  RESULTS results ;
  //
  for ( Ostap::EventIndex   entry = first ; entry < the_last ; ++entry , ++bar )
    {
      //
      const long ievent = data -> GetEntryNumber ( entry ) ;
      if ( 0 > ievent                      ) { return INVALID_ENTRY ; } // RETURN 
      if ( 0 > data -> LoadTree ( ievent ) ) { return INVALID_EVENT ; } // RETURN 
      //
      const double weight = with_cuts ? cuts->evaluate() : 1.0 ;
      if ( !weight ) { continue ; }                                    // ATTENTION! 
      //
      formulae.evaluate ( 0 , results [ 0 ] ) ;
      formulae.evaluate ( 1 , results [ 1 ] ) ;
      formulae.evaluate ( 2 , results [ 2 ] ) ;
      formulae.evaluate ( 3 , results [ 3 ] ) ;
      //
      for ( const double x : results [ 0 ] )
      { for ( const double y : results [ 1 ] )
        { for ( const double z : results [ 2 ] )
          { for ( const double t : results [ 3 ] )
            { stat.update ( x , y , z , t , weight ) ; } } } }
    }
  //
  return Ostap::StatusCode::SUCCESS ;
}
// =============================================================================
/*  Fill/update statistical counter 
 *  @param data       (input) data 
 *  @param stat       (UDATE) statistical counter 
 *  @param expression (INPUT) the variable 
 *  @param selection  (INPUT) selection/cut (treated as weight!)
 *  @param cut_range  (INPUT) if non empty: use evene only fomthis cut-range
 *  @param first      (INPIUT) the first event to process (inclusibe) 
 *  @param last       (INPIUT) the last event to process (exclusive) 
 *  @return status code 
 *  @see Ostap::Math::WStatistic
 *  @attention selection/cut is treated as weight!
 */
// =============================================================================
Ostap::StatusCode Ostap::StatVar::get_stat
( const RooAbsData*         data       ,
  Ostap::Math::WStatistic4& stat       ,
  const std::string&        expr1      ,
  const std::string&        expr2      ,  
  const std::string&        expr3      ,
  const std::string&        expr4      ,  
  const std::string&        selection  , 
  const std::string&        cut_range  ,
  const Ostap::EventIndex   first      ,
  const Ostap::EventIndex   last       ) const
{
  /// reset the statistic
  stat.reset() ; // reset the statistic
  // 
  if ( nullptr == data    ) { return INVALID_DATA    ; }
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  //
  const Ostap::EventIndex num_entries = data -> numEntries() ;
  const Ostap::EventIndex the_last    = std::min ( num_entries , last ) ;
  if ( the_last <= first  ) { return Ostap::StatusCode::SUCCESS ; }
  //
  const bool  weighted = data -> isWeighted() ;
  const char* cutrange = cut_range.empty() ? nullptr : cut_range.c_str() ;
  //
  // formulae for exressons
  const Ostap::FormulaVars                 formulae { data , { expr1 , expr2 , expr3 , expr4 } } ;
  const std::unique_ptr<Ostap::FormulaVar> cuts { make_formula ( selection  , data , true ) } ;
  const bool                with_cuts   = cuts && cuts->ok () ;
  //
  Ostap::Utils::ProgressBar bar { the_last - first , m_progress } ; 
  for ( Ostap::EventIndex entry = first ; entry < the_last ; ++entry , ++bar )
    {
      const RooArgSet* vars = data -> get ( entry ) ;
      if ( nullptr == vars )                             { break    ; } // BREAK 
      //
      if ( cutrange && !vars->allInRange ( cutrange ) )  { continue ; } // CONTINUE    
      // apply weight:
      const double wd = weighted  ? data -> weight () : 1.0 ;
      if ( !wd ) { continue ; }                                   // CONTINUE    
      // apply cuts:
      const double wc = with_cuts ? cuts -> getVal () : 1.0 ;
      if ( !wc ) { continue ; }                                   // CONTINUE  
      // Total: cuts & weight:
      const double weight  = wd *  wc ;
      if ( !weight  ) { continue ; }                             // CONTINUE        
      //
      const double x = formulae.evaluate ( 0 ) ; 
      const double y = formulae.evaluate ( 1 ) ; 
      const double z = formulae.evaluate ( 2 ) ; 
      const double t = formulae.evaluate ( 3 ) ; 
      //
      stat.update ( x , y , z , t ,  weight ) ;
    }
  //
  return Ostap::StatusCode::SUCCESS ; 
}
// =============================================================================


// =============================================================================
/*  build statistic for the <code>expressions</code>
 *  @param tree        (INPUT)  the tree 
 *  @param result      (UPDATE) the output statistics for specified expressions 
 *  @param expressions (INPUT)  the list of  expressions
 *  @param cuts        (INPUT)  the selection criteria 
 *  @param first       (INPUT)  the first entry to process 
 *  @param last        (INPUT)  the last entry to process (not including!)
 *  @return number of processed entries 
 *  @attentoon   selection is treated as boolean!
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2018-11-04
 */
// ============================================================================
Ostap::StatusCode Ostap::StatVar::statVars
( TTree*                       data        ,  
  Ostap::StatVar::StatVector&  result      ,  
  const Ostap::Strings&        expressions ,
  const std::string&           selection   ,  
  const Ostap::EventIndex      first       ,
  const Ostap::EventIndex      last        ) const
{
  //
  const std::size_t N = expressions.size() ;
  result.resize ( N ) ; 
  for ( auto& r : result ) { r.reset () ; }
  //
  if ( nullptr == data     ) { return INVALID_TREE ; }
  if ( expressions.empty() ) { return Ostap::StatusCode::SUCCESS ; }
  //
  const Ostap::EventIndex num_entries =  data->GetEntries() ;
  const Ostap::EventIndex the_last    = std::min ( first , num_entries ) ;
  if ( the_last <= first   ) { return Ostap::StatusCode::SUCCESS ; }
  //
  /// formulae for exressons
  const Ostap::Formulae                 formulae { data , expressions } ;  
  const std::unique_ptr<Ostap::Formula> cuts     { make_formula ( selection , data , true ) } ;
  const bool                with_cuts = cuts && cuts->ok() ;
  //
  Ostap::Utils::Notifier    notify  ( formulae.begin() , formulae.end() , cuts.get() , data ) ;
  //
  std::vector<double>       results {} ;
  Ostap::Utils::ProgressBar bar     ( the_last - first , m_progress ) ;
  for ( Ostap::EventIndex   entry = first ; entry < the_last  ; ++entry , ++bar )
    {
      //
      const long ievent = data -> GetEntryNumber ( entry ) ;
      if ( 0 > ievent                    ) { return INVALID_ENTRY ; } // RETURN
      if ( 0 > data->LoadTree ( ievent ) ) { return INVALID_EVENT ; } // RETURN
      //
      const double cut = with_cuts ? cuts ->evaluate() : 1.0 ;
      if ( !cut ) { continue  ; }
      //
      for ( std::size_t i = 0 ; i < N ; ++i ) 
	{
	  formulae.evaluate ( i , results ) ;
	  for ( const double value : results ) { result[i].add ( value ) ; }
	}
    }
  //
  return Ostap::StatusCode::SUCCESS ;
}
// =============================================================================
/*  build statistic for the <code>expressions</code>
 *  @param tree        (INPUT)  the tree 
 *  @param result      (UPDATE) the output statistics for specified expressions 
 *  @param expressions (INPUT)  the list of  expressions
 *  @param cuts        (INPUT)  the selection criteria 
 *  @param first       (INPUT)  the first entry to process 
 *  @param last        (INPUT)  the last entry to process (not including!)
 *  @return number of processed entries 
 *  @attentoon   selection is treated as weight 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2018-11-04
 */
// ============================================================================
Ostap::StatusCode Ostap::StatVar::statVars
( TTree*                       data        ,  
  Ostap::StatVar::WStatVector& result      ,  
  const Ostap::Strings&        expressions ,
  const std::string&           selection   ,  
  const Ostap::EventIndex      first       ,
  const Ostap::EventIndex      last        ) const
{
  //
  const std::size_t N = expressions.size() ;
  result.resize ( N ) ; 
  for ( auto& r : result ) { r.reset () ; }
  //
  if ( nullptr == data     ) { return INVALID_TREE ; }
  if ( expressions.empty() ) { return Ostap::StatusCode::SUCCESS ; }
  //
  const Ostap::EventIndex num_entries =  data->GetEntries() ;
  const Ostap::EventIndex the_last    = std::min ( first , num_entries ) ;
  if ( the_last <= first   ) { return Ostap::StatusCode::SUCCESS ; }
  //
  /// formulae for exressons
  const Ostap::Formulae                 formulae { data , expressions } ;  
  const std::unique_ptr<Ostap::Formula> cuts     { make_formula ( selection , data , true ) } ;  
  // 
  Ostap::Utils::Notifier    notify  ( formulae.begin() , formulae.end() , cuts.get() , data ) ;
  std::vector<double>       results {} ;
  Ostap::Utils::ProgressBar bar     ( the_last - first , m_progress ) ;
  for ( Ostap::EventIndex   entry = first ; entry < the_last  ; ++entry , ++bar )
    {
      //
      const long ievent = data -> GetEntryNumber ( entry ) ;
      if ( 0 > ievent                    ) { return INVALID_ENTRY ; } // RETURN
      if ( 0 > data->LoadTree ( ievent ) ) { return INVALID_EVENT ; } // RETURN
      //
      const double weight = cuts ? cuts ->evaluate() : 1.0 ;
      if ( !weight ) { continue  ; }
      //
      for ( std::size_t i = 0 ; i < N ; ++i ) 
	{
	  formulae.evaluate ( i , results ) ;
	  for ( const double value : results ) { result[i].add ( value , weight ) ; }
	}
    }
  //
  return Ostap::StatusCode::SUCCESS ;
}
// ============================================================================
/** build statistic for the <code>expressions</code>
 *  @param data        (INPUT)  input data 
 *  @param result      (UPDATE) the output statistics for specified expressions 
 *  @param expressions (INPUT)  the list of  expressions
 *  @param cuts        (INPUT)  the selection 
 *  @param cut_range   (INPUT)  cut range  
 *  @param first       (INPUT)  the first entry to process 
 *  @param last        (INPUT)  the last entry to process (not including!)
 *  @return number of processed entries 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2021-06-04
 */
// ============================================================================
Ostap::StatusCode 
Ostap::StatVar::statVars
( const RooAbsData*               data        , 
  Ostap::StatVar::WStatVector&    result      , 
  const Ostap::Strings&           expressions ,
  const std::string&              selection   ,
  const std::string&              cut_range   ,
  const Ostap::EventIndex         first       ,
  const Ostap::EventIndex         last        ) const
{
  // 
  const std::size_t N = expressions.size() ;
  result.resize ( N ) ; 
  for ( auto& r : result ) { r.reset () ; }
  //
  if ( nullptr == data     ) { return INVALID_DATA ; }
  if ( expressions.empty() ) { return Ostap::StatusCode::SUCCESS ; }
  //
  const Ostap::EventIndex num_entries =  data->numEntries() ;
  const Ostap::EventIndex the_last    = std::min ( first , num_entries ) ;
  if ( the_last <= first   ) { return Ostap::StatusCode::SUCCESS ; }
  //
  /// formulae for exressons
  const Ostap::FormulaVars                 formulae { data , expressions } ;    
  const std::unique_ptr<Ostap::FormulaVar> cuts     { make_formula ( selection , data , true ) } ;
  //
  const char* cutrange = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  const bool  weighted = data->isWeighted() ;
  //
  // start the loop
  Ostap::Utils::ProgressBar bar ( the_last - first , m_progress ) ;
  for ( Ostap::EventIndex   entry = first ; entry < the_last ; ++entry , ++bar )
    {
      //
      const RooArgSet* vars = data->get( entry ) ;
      if ( nullptr == vars  )                           { break    ; } // RETURN
      if ( cutrange && !vars->allInRange ( cutrange ) ) { continue ; } // CONTINUE    
      //
      // apply weight:
      const double wd = weighted  ? data->weight()        : 1.0 ;
      if ( !wd ) { continue ; }                                   // CONTINUE    
      // apply cuts:
      const double wc = cuts ? cuts -> getVal() : 1.0 ;
      if ( !wc ) { continue ; }                                   // CONTINUE  
      // cuts & weight:
      const double weight  = wd *  wc ;
      if ( !weight ) { continue ; }                                   // CONTINUE        
      //
      for ( std::size_t i = 0 ; i < N ; ++i ) 
	{
	  const double value = formulae.evaluate ( i ) ;
	  result [ i ].add ( value , weight ) ;
	}
      //
    }
  //
  return Ostap::StatusCode::SUCCESS ; 
}
// ==========================================================================


// =========================================================================
// Collect information for single variable 
// =========================================================================

// =========================================================================
/* build statistic for the <code>expression</code>
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
// =========================================================================
Ostap::StatEntity Ostap::StatVar::statVar
( TTree*                  data        , 
  const std::string&      expression  , 
  const Ostap::EventIndex first       ,
  const Ostap::EventIndex last        ) const
{
  /// 
  Ostap::StatVar::StatVector vct         { 1 , Ostap::StatEntity() } ;
  const Ostap::Strings       expressions { 1 , expression          } ; 
  const Ostap::StatusCode    sc = statVars ( data , vct , expressions , "" , first , last ) ;
  Ostap::Assert ( sc.isSuccess ()  ,
		  "Error from Ostap::StatVar::statVars" ,
		  "Ostap::StatVar::statVar" , sc , __FILE__ , __LINE__ ) ;
  Ostap::Assert ( 1 == vct.size()  ,
		  "Ivvalid size of counters array " ,
		  "Ostap::StatVar::statVar" ,
		  INVALID_COUNTERS , __FILE__ , __LINE__ ) ;
  return vct.front () ;		  
}
// =========================================================================
/*  build statistic for the <code>expression</code>
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
// =========================================================================
Ostap::StatEntity Ostap::StatVar::statVar_cut
( TTree*                  data       , 
  const std::string&      expression , 
  const std::string&      selection  , 
  const Ostap::EventIndex first      ,
  const Ostap::EventIndex last       ) const
{
  /// 
  Ostap::StatVar::StatVector vct         { 1 , Ostap::StatEntity() } ;
  const Ostap::Strings       expressions { 1 , expression          } ; 
  const Ostap::StatusCode    sc = statVars ( data , vct , expressions , selection , first , last ) ;
  Ostap::Assert ( sc.isSuccess ()  ,
		  "Error from Ostap::StatVar::statVars" ,
		  "Ostap::StatVar::statVar_cut" , sc , __FILE__ , __LINE__ ) ;
  Ostap::Assert ( 1 == vct.size()  ,
		  "Ivvalid size of counters array " ,
		  "Ostap::StatVar::statVar_cut" ,
		  INVALID_COUNTERS , __FILE__ , __LINE__ ) ;
  return vct.front () ;		  
}
// =========================================================================
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
// =========================================================================
Ostap::WStatEntity Ostap::StatVar::statVar
( TTree*                  data       , 
  const std::string&      expression , 
  const std::string&      selection  , 
  const Ostap::EventIndex first      ,
  const Ostap::EventIndex last       ) const
{
  /// 
  Ostap::StatVar::WStatVector vct         { 1 , Ostap::WStatEntity() } ;
  const Ostap::Strings        expressions { 1 , expression           } ; 
  const Ostap::StatusCode     sc = statVars ( data , vct , expressions , selection , first , last ) ;
  Ostap::Assert ( sc.isSuccess ()  ,
		  "Error from Ostap::StatVar::statVars" ,
		  "Ostap::StatVar::statVar" , sc , __FILE__ , __LINE__ ) ;
  Ostap::Assert ( 1 == vct.size()  ,
		  "Ivvalid size of counters array " ,
		  "Ostap::StatVar::statVar" ,
		  INVALID_COUNTERS , __FILE__ , __LINE__ ) ;
  return vct.front () ;		  
}
// =========================================================================
/* build statistic for the <code>expression</code>
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
// =========================================================================
Ostap::WStatEntity Ostap::StatVar::statVar 
( const RooAbsData*       data       , 
  const std::string&      expression , 
  const std::string&      selection  ,
  const std::string&      cut_range  ,  
  const Ostap::EventIndex first      ,
  const Ostap::EventIndex last       ) const
{
  /// 
  Ostap::StatVar::WStatVector vct         { 1 , Ostap::WStatEntity() } ;
  const Ostap::Strings        expressions { 1 , expression           } ; 
  const Ostap::StatusCode     sc = statVars ( data , vct , expressions , selection , cut_range , first , last ) ;
  Ostap::Assert ( sc.isSuccess ()  ,
		  "Error from Ostap::StatVar::statVars" ,
		  "Ostap::StatVar::statVar" , sc , __FILE__ , __LINE__ ) ;
  Ostap::Assert ( 1 == vct.size()  ,
		  "Ivvalid size of counters array " ,
		  "Ostap::StatVar::statVar" ,
		  INVALID_COUNTERS , __FILE__ , __LINE__ ) ;
  return vct.front () ;		  
}


// ============================================================================
// Covariances for two variables 
// ============================================================================

// ============================================================================
/*  calculate the covariance of two expressions
 *  @param tree  (INPUT)  the input tree
 *  @param exp1  (INPUT)  the first  expression
 *  @param exp2  (INPUT)  the second expression
 *  @return covariance
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date   2024-07-22
 */
// ============================================================================
Ostap::Math::Covariance
Ostap::StatVar::statCov
( TTree*                  data  ,
  const std::string&      exp1  ,
  const std::string&      exp2  ,
  const Ostap::EventIndex first ,
  const Ostap::EventIndex last  ) const 
{
  //
  Ostap::Assert ( nullptr != data                    ,  
		  "Invalid TTree"                    , 
		  "Ostap::StatVar::statCov"          ,
		  INVALID_TREE , __FILE__ , __LINE__ ) ;
  //
  // prepare the result 
  Ostap::Math::Covariance result {} ;
  //
  const Ostap::EventIndex num_entries = data -> GetEntries() ;
  const Ostap::EventIndex the_last    = std::min ( last , num_entries ) ;
  if ( the_last <= first ) { return result ; } 
  //
  // formulas for expressins
  const Ostap::Formulae formulae { data , { exp1 , exp2 } } ;
  //
  Ostap::Utils::Notifier notify ( formulae.begin()  , formulae.end()  , data ) ;
  //
  std::vector<double> results1 {} ;
  std::vector<double> results2 {} ;
  //
  Ostap::Utils::ProgressBar bar ( the_last - first , m_progress ) ;
  for ( Ostap::EventIndex entry = first ; entry < the_last ; ++entry , ++bar )
    {
      //
      const long ievent = data->GetEntryNumber ( entry ) ;
      if ( 0 > ievent                    ) { break ; } // BREAK
      if ( 0 > data->LoadTree ( ievent ) ) { break ; } // BREAK
      //
      formulae.evaluate ( 0 , results1 ) ;
      formulae.evaluate ( 1 , results2 ) ;
      //
      for ( const long double value1 : results1 ) 
	{  for ( const long double value2 : results2 ) 
	    { result.add ( value1 , value2 ) ; } }
      //
    }
  //
  return result ;
}
// ============================================================================
/*  calculate the covariance of two expressions 
 *  @param tree  (INPUT)  the inpout tree 
 *  @param exp1  (INPUT)  the first  expression
 *  @param exp2  (INPUT)  the second expression
 *  @param cuts  (INPUT)  cuts/weigths expression 
 *  @return Covariance 
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date   2024-07-22
 */
// ============================================================================
Ostap::Math::Covariance
Ostap::StatVar::statCov_cut 
( TTree*                  data      ,
  const std::string&      exp1      ,
  const std::string&      exp2      ,
  const std::string&      selection ,
  const Ostap::EventIndex first     ,
  const Ostap::EventIndex last      ) const 
{
  //
  Ostap::Assert ( nullptr != data                    ,  
		  "Invalid TTree"                    , 
		  "Ostap::StatVar::statCov"          ,
		  INVALID_TREE , __FILE__ , __LINE__ ) ;
  //
  // prepare the result 
  Ostap::Math::Covariance result {} ;
  //
  const Ostap::EventIndex num_entries = data -> GetEntries() ;
  const Ostap::EventIndex the_last    = std::min ( last , num_entries ) ;
  if ( the_last <= first ) { return result ; } 
  //
  // formulas for expressins
  const Ostap::Formulae formulae             { data , { exp1 , exp2 } } ;
  const std::unique_ptr<Ostap::Formula> cuts { make_formula ( selection , data , true ) } ;
  //
  Ostap::Utils::Notifier notify ( formulae.begin()  , formulae.end()  , cuts.get() , data ) ;
  //
  const bool with_cuts = cuts && cuts ->ok () ;

  std::vector<double> results1 {} ;
  std::vector<double> results2 {} ;
  //
  Ostap::Utils::ProgressBar bar ( the_last - first , m_progress ) ;
  for ( Ostap::EventIndex entry = first ; entry < the_last ; ++entry , ++bar )
    {
      //
      const long ievent = data->GetEntryNumber ( entry ) ;
      if ( 0 > ievent                    ) { break ; } // BREAK
      if ( 0 > data->LoadTree ( ievent ) ) { break ; } // BREAK
      //
      const double cut = with_cuts ? cuts->evaluate() : 1.0 ;
      if ( !cut ) { continue ;}                      // ATTENTION
      //
      formulae.evaluate ( 0 , results1 ) ;
      formulae.evaluate ( 1 , results2 ) ;      
      //
      for ( const long double value1 : results1 ) 
	{  for ( const long double value2 : results2 ) 
	    { result.add ( value1 , value2 ) ; } }
      //
    }
  //
  return result ;
}

// ============================================================================
/*  calculate the covariance of two expressions 
 *  @param tree  (INPUT)  the inpout tree 
 *  @param exp1  (INPUT)  the first  expression
 *  @param exp2  (INPUT)  the second expression
 *  @param cuts  (INPUT)  cuts/weigths expression 
 *  @return Covariance 
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date   2024-07-22
 */
// ============================================================================
Ostap::Math::WCovariance
Ostap::StatVar::statCov
( TTree*                  data      ,
  const std::string&      exp1      ,
  const std::string&      exp2      ,
  const std::string&      selection ,
  const Ostap::EventIndex first     ,
  const Ostap::EventIndex last      ) const
{
  //
  Ostap::Assert ( nullptr != data                    ,  
		  "Invalid TTree"                    , 
		  "Ostap::StatVar::statCov"          ,
		  INVALID_TREE , __FILE__ , __LINE__ ) ;
  //
  // prepare the result 
  Ostap::Math::WCovariance result {} ;
  //
  const Ostap::EventIndex num_entries = data -> GetEntries() ;
  const Ostap::EventIndex the_last    = std::min ( last , num_entries ) ;
  if ( the_last <= first ) { return result ; } 
  //
  // formulas for expressins
  const Ostap::Formulae formulae             { data , { exp1 , exp2 } } ;
  const std::unique_ptr<Ostap::Formula> cuts { make_formula ( selection , data , true ) } ;
  //
  Ostap::Utils::Notifier notify ( formulae.begin() , formulae.end() , cuts.get() , data ) ;
  //
  const bool with_cuts = cuts && cuts ->ok () ;
  //
  std::vector<double> results1 {} ;
  std::vector<double> results2 {} ;
  //
  Ostap::Utils::ProgressBar bar ( the_last - first , m_progress ) ;
  for ( Ostap::EventIndex entry = first ; entry < the_last ; ++entry , ++bar )
    {
      //
      const long ievent = data->GetEntryNumber ( entry ) ;
      if ( 0 > ievent                    ) { break ; } // BREAK
      if ( 0 > data->LoadTree ( ievent ) ) { break ; } // BREAK
      //
      const double weight = with_cuts ? cuts->evaluate() : 1.0 ;
      if ( !weight ) { continue ;}                      // ATTENTION
      //
      formulae.evaluate ( 0 , results1 ) ;
      formulae.evaluate ( 1 , results2 ) ;      
      //
      for ( const long double value1 : results1 ) 
	{  for ( const long double value2 : results2 ) 
	    { result.add ( value1 , value2 , weight ) ; } }
      //
    }
  //
  return result ;  
}
// ============================================================================
/*  calculate the covariance of two expressions 
 *  @param tree  (INPUT)  the input  tree 
 *  @param exp1  (INPUT)  the first  expresiion
 *  @param exp2  (INPUT)  the second expresiion
 *  @param cuts  (INPUT)  selection/weight expression 
 *  @return number of processed events 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2014-03-27
 */
// ============================================================================
Ostap::Math::WCovariance
Ostap::StatVar::statCov
( const RooAbsData*       data      , 
  const std::string&      exp1      , 
  const std::string&      exp2      , 
  const std::string&      selection , 
  const std::string&      cut_range , 
  const Ostap::EventIndex first     ,
  const Ostap::EventIndex last      ) const 
{
  // ===========================================================================
  Ostap::Assert ( nullptr != data                    ,  
		  "Invalid RotAbsData"               , 
		  "Ostap::StatVar::statCov"          ,
		  INVALID_DATA , __FILE__ , __LINE__ ) ;
  //
  // prepare the result 
  Ostap::Math::WCovariance result {} ;
  //
  const Ostap::EventIndex num_entries = data->numEntries() ;
  const Ostap::EventIndex the_last    = std::min ( last , num_entries ) ;
  if ( the_last <= first ) { return result ; }         // RETURN
  //
  const Ostap::FormulaVars                 formulae { data , { exp1 , exp2 } } ;
  const std::unique_ptr<Ostap::FormulaVar> cuts     { make_formula ( selection , data , true ) } ;
  //
  const char* cutrange  = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  const bool  weighted  = data->isWeighted () ;
  const bool  with_cuts = cuts && cuts->ok () ;
  //
  // start the loop
  Ostap::Utils::ProgressBar bar ( the_last - first , m_progress ) ;
  for ( unsigned long entry = first ; entry < the_last ; ++entry , ++bar )
    {
      //
      const RooArgSet* vars = data->get( entry ) ;
      if ( nullptr == vars                            ) { break    ; } // RETURN
      if ( cutrange && !vars->allInRange ( cutrange ) ) { continue ; } // CONTINUE    
      //
      // apply weight:
      const double wd = weighted  ? data->weight()        : 1.0 ;
      if ( !wd ) { continue ; }                                   // CONTINUE    
      // apply cuts:
      const double wc = with_cuts ? cuts  -> getVal() : 1.0 ;
      if ( !wc ) { continue ; }                                   // CONTINUE  
      // cuts & weight:
      const double weight = wd *  wc ;
      if ( !weight ) { continue ; }                                   // CONTINUE        
      //
      const double value1 = formulae.evaluate ( 0 ) ;
      const double value2 = formulae.evaluate ( 1 ) ;
      //
      result.add ( value1 , value2 , weight ) ;
    }
  //
  return result ;
}


// ============================================================================
// Covarinaes for several variables 
// ============================================================================

// ============================================================================
/*  calculate the covariance of several expressions 
 *  @param data        (INPUT)  the inpout tree 
 *  @param stats       (UPDATE) the statistics 
 *  @param cov2        (UPDATE) the covariance matrix 
 *  @param expressions (INPUT)  expressions 
 *  @param selection   (INPUT)  the selection criteria 
 *  @return status code 
 *  @attention selection is treate as boolen! 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2023-02-28
 */
// ============================================================================
Ostap::StatusCode Ostap::StatVar::statCov
( TTree*                          data        ,   
  Ostap::Math::Covariances&       stats       ,  
  const Ostap::Strings&           expressions , 
  const std::string&              selection   ,
  const Ostap::EventIndex         first       ,
  const Ostap::EventIndex         last        ) const
{
  //
  const std::size_t N = expressions.size() ;
  stats = Ostap::Math::Covariances (  N  ) ; 
  //
  if ( nullptr == data     ) { return INVALID_TREE ; }
  if ( expressions.empty() ) { return Ostap::StatusCode::SUCCESS ; }
  //
  const Ostap::EventIndex num_entries =  data->GetEntries() ;
  const Ostap::EventIndex the_last    = std::min ( first , num_entries ) ;
  if ( the_last <= first   ) { return Ostap::StatusCode::SUCCESS ; }
  //
  /// formulae for exressons
  const Ostap::Formulae                 formulae { data , expressions } ;  
  const std::unique_ptr<Ostap::Formula> cuts     { make_formula ( selection , data , true ) } ;    
  // 
  Ostap::Utils::Notifier    notify  ( formulae.begin() , formulae.end() , cuts.get() , data ) ;
  // 
  const bool                with_cuts = cuts && cuts->ok() ;
  //
  typedef std::vector<double> RESULT  ;
  typedef std::vector<RESULT> RESULTS ;
  RESULTS results {} ;
  RESULT  input   ( N , 0.0 ) ;
  //
  Ostap::Utils::ProgressBar bar {the_last - first , m_progress } ;
  for ( Ostap::EventIndex    entry = first ; entry < the_last  ; ++entry , ++bar )
    {
      //
      const long ievent = data -> GetEntryNumber ( entry ) ;
      if ( 0 > ievent                    ) { return INVALID_ENTRY ; } // RETURN
      if ( 0 > data->LoadTree ( ievent ) ) { return INVALID_EVENT ; } // RETURN
      //
      const double cut = with_cuts ? cuts ->evaluate() : 1.0 ;
      if ( !cut ) { continue  ; }
      //
      for ( std::size_t i = 0 ; i < N ; ++i )  { formulae.evaluate ( i , results [ i ] ) ; }
      //
      Ostap::Combiner_<RESULT> combiner { results.begin() , results.end() };
      while ( combiner.valid() ) 
      {
        //get the curren combination 
        combiner.current ( input.begin() ) ;
        // update counter 
        stats.add ( input ) ; 
        // go to next combination
        ++combiner ; 
      } 
    }
  //
  return Ostap::StatusCode::SUCCESS ; 
}
// ============================================================================
/*  calculate the covariance of several expressions 
 *  @param data        (INPUT)  the inpout tree 
 *  @param stats       (UPDATE) the statistics 
 *  @param cov2        (UPDATE) the covariance matrix 
 *  @param expressions (INPUT)  expressions 
 *  @param selection   (INPUT)  the selection criteria 
 *  @return status code 
 *  @attention selection is treate as weight 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2023-02-28
 */
// ============================================================================
Ostap::StatusCode 
Ostap::StatVar::statCov
( TTree*                          data        ,   
  Ostap::Math::WCovariances&      stats       ,  
  const Ostap::Strings&           expressions , 
  const std::string&              selection   ,
  const Ostap::EventIndex         first       ,
  const Ostap::EventIndex         last        ) const
{
  const std::size_t N = expressions.size() ;
  stats = Ostap::Math::WCovariances (  N  ) ; 
  //
  if ( nullptr == data     ) { return INVALID_TREE ; }
  if ( expressions.empty() ) { return Ostap::StatusCode::SUCCESS ; }
  //
  const Ostap::EventIndex num_entries =  data->GetEntries() ;
  const Ostap::EventIndex the_last    = std::min ( first , num_entries ) ;
  if ( the_last <= first   ) { return Ostap::StatusCode::SUCCESS ; }
  //
  /// formulae for exressons
  const Ostap::Formulae                 formulae { data , expressions } ;  
  const std::unique_ptr<Ostap::Formula> cuts     { make_formula ( selection , data , true ) } ;    
  // 
  Ostap::Utils::Notifier    notify  ( formulae.begin() , formulae.end() , cuts.get() , data ) ;
  // 
  const bool                with_cuts = cuts && cuts->ok() ;
  //
  typedef std::vector<double> RESULT  ;
  typedef std::vector<RESULT> RESULTS ;
  RESULTS results {} ;
  RESULT  input   ( N , 0.0 ) ;
  //
  Ostap::Utils::ProgressBar bar {the_last - first , m_progress } ;
  for ( Ostap::EventIndex    entry = first ; entry < the_last  ; ++entry , ++bar )
    {
      //
      const long ievent = data -> GetEntryNumber ( entry ) ;
      if ( 0 > ievent                    ) { return INVALID_ENTRY ; } // RETURN
      if ( 0 > data->LoadTree ( ievent ) ) { return INVALID_EVENT ; } // RETURN
      //
      const double weight  = with_cuts ? cuts ->evaluate() : 1.0 ;
      if ( !weight ) { continue  ; }
      //
      for ( std::size_t i = 0 ; i < N ; ++i )  { formulae.evaluate ( i , results [ i ] ) ; }
      //
      Ostap::Combiner_<RESULT> combiner { results.begin() , results.end() };
      while ( combiner.valid() ) 
      {
        //get the curren combination 
        combiner.current ( input.begin() ) ;
        // update counter 
        stats.add ( input , weight ) ; 
        // go to next combination
        ++combiner ; 
      } 
    }
  //
  return Ostap::StatusCode::SUCCESS ; 
}
// ============================================================================
/*  calculate the covariance of several expressions 
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
// ============================================================================
Ostap::StatusCode 
Ostap::StatVar::statCov
( const RooAbsData*               data        , 
  Ostap::Math::WCovariances&      stats       ,   
  const Ostap::Strings&           expressions , 
  const std::string&              selection   ,
  const std::string&              cut_range   , 
  const Ostap::EventIndex         first       , 
  const Ostap::EventIndex         last        ) const
{
  // 
  const std::size_t N = expressions.size() ;
  stats = Ostap::Math::WCovariances (  N  ) ; 
  //
  if ( nullptr == data     ) { return INVALID_DATA ; }
  if ( expressions.empty() ) { return Ostap::StatusCode::SUCCESS ; }
  //
  const Ostap::EventIndex num_entries =  data -> numEntries () ;
  const Ostap::EventIndex the_last    = std::min ( first , num_entries ) ;
  if ( the_last <= first   ) { return Ostap::StatusCode::SUCCESS ; }
  //
  /// formulae for exressons
  const Ostap::FormulaVars                 formulae { data , expressions } ;  
  const std::unique_ptr<Ostap::FormulaVar> cuts     { make_formula ( selection , data , true ) } ;  
  // 
  const bool  with_cuts = cuts && cuts->ok() ;
  const char* cutrange  = cut_range.empty() ? nullptr : cut_range.c_str() ;
  const bool  weighted  = data->isWeighted() ;
  //
  typedef std::vector<double> RESULTS ;
  RESULTS result ( N , 0.0 ) ;
  //
  Ostap::Utils::ProgressBar bar {the_last - first , m_progress } ;
  for ( Ostap::EventIndex   entry = first ; entry < the_last  ; ++entry , ++bar )
    {
      //
      const RooArgSet* vars = data->get( entry ) ;
      if ( nullptr == vars  )                           { break    ; } // BREAK
      /// 
      if ( cutrange && !vars->allInRange ( cutrange ) ) { continue ; } // CONTINUE
      /// apply weight
      const double wd = weighted ? data->weight() : 1.0 ;
      if ( !wd )                                        { continue ; } // CONTINUE  
      //
      const double wc = with_cuts ? cuts->getVal() : 1.0 ;
      if ( !wc )                                        { continue ; } // CONTINUE  
      //
      const double weight = wd * wc ;
      if ( !weight )                                    { continue ; } // CONTINUE  
      //
      // evaluate all expressions 
      for ( std::size_t i = 0 ; i < N ; ++i ) { result[i] = formulae.evaluate ( i ) ; }
      //
      stats.add ( result , weight ) ;
    } 
  //
  return Ostap::StatusCode::SUCCESS ; 
}
// ============================================================================

// ============================================================================
//                                                                      The END
// ============================================================================
