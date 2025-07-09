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
#include "Ostap/StatusCode.h"
#include "Ostap/Combiner.h"
#include "Ostap/Formula.h"
#include "Ostap/Notifier.h"
#include "Ostap/StatVar.h"
#include "Ostap/FormulaVar.h"
#include "Ostap/Covariance.h"
#include "Ostap/Covariances.h"
#include "Ostap/ProgressBar.h"
#include "Ostap/GetWeight.h"
// ============================================================================
// Local
// ============================================================================
#include "Formulae.h"
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
  inline bool in_range
  ( const Ostap::DataType value                  ,
    const Ostap::DataType xmin = Ostap::MinValue , 
    const Ostap::DataType xmax = Ostap::MaxValue )
  { return std::isfinite ( value )
      && ( xmin  <= value )
      && ( value <= xmax  ) ; }
  // ==========================================================================
} //                                             The end of anonymous namespace 
// ============================================================================
// constructor with ProgressBar configuration
// =============================================================================
Ostap::StatVar::StatVar
( const Ostap::Utils::ProgressConf& progress)
  : m_progress ( progress )
{}
// =============================================================================
// The basic methods 
// =============================================================================
/*  Fill/update statistical counter 
 *  @param data       (input) data 
 *  @param stat       (UDATE) statistical counter 
 *  @param expression (INPUT) the variable 
 *  @param selection  (INPUT) selection/cut (treated as boolean!)
 *  @param first      (INPIUT) the first event to process (inclusibe) 
 *  @param last       (INPIUT) the last event to process (exclusive) 
 *  @param xmin       (INPUT)  low  limit for expressoon 
 *  @param xmax       (INPUT)  high limit for expressoon 
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
  const Ostap::EventIndex   last       , 
  const Ostap::DataType     xmin       ,
  const Ostap::DataType     xmax       ) const 
{
  /// reset the statistic
  stat.reset() ; // reset the statistic
  //
  // check the ranges
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  if ( xmax <  xmin       ) { return Ostap::StatusCode::SUCCESS ; }
  //
  if ( nullptr == data    ) { return INVALID_DATA    ; }
  //
  const Ostap::EventIndex  num_entries = data -> GetEntries () ;
  const Ostap::EventIndex  the_last    = std::min ( last , num_entries ) ;
  if ( the_last <= first ) { return Ostap::StatusCode::SUCCESS ; }
  //
  Ostap::Formula formula ( expression , data ) ;
  if ( !formula.ok ()     ) { return INVALID_FORMULA ; }
  //
  const std::unique_ptr<Ostap::Formula> cuts { Ostap::makeFormula ( selection , data , true ) } ;
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
      const bool cut = with_cuts ? cuts->evaluate() : 1.0  ;
      if ( !cut ) { continue ; }                                         // ATTENTION! 
      //
      formula.evaluate ( results ) ;
      for ( const double value : results )
        { if ( in_range ( value , xmin , xmax ) ) { stat.update ( value ) ; } }
    }
  //
  return Ostap::StatusCode::SUCCESS ;
} ;
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
  const Ostap::EventIndex   last       , 
  const Ostap::DataType     xmin       ,
  const Ostap::DataType     xmax       ) const 
{
  /// reset the statistic
  stat.reset() ; // reset the statistic
  //
  // check the ranges
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  if ( xmax <  xmin       ) { return Ostap::StatusCode::SUCCESS ; }
  //
  if ( nullptr == data    ) { return INVALID_DATA    ; }
  //
  const Ostap::EventIndex   num_entries = data -> GetEntries()  ;
  const Ostap::EventIndex   the_last    = std::min ( last , num_entries ) ;
  if ( the_last <= first   ) { return Ostap::StatusCode::SUCCESS ; }
  //
  Ostap::Formula formula ( expression , data ) ;
  if ( !formula.ok()      ) { return INVALID_FORMULA ; }
  //
  const std::unique_ptr<Ostap::Formula> cuts { Ostap::makeFormula ( selection , data , true ) } ; 
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
      if ( !weight || !std::isfinite ( weight ) ) { continue ; }        // ATTENTION! 
      //
      formula.evaluate ( results ) ;
      for  ( const double value : results )
        { if ( in_range ( value , xmin , xmax ) ) { stat.update ( value , weight ) ; } }
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
 *  @attention selection is treated as boolean 
 *  @attention for weighted datasets, the weigth is treated as boolean! 
 */
// =============================================================================
Ostap::StatusCode Ostap::StatVar::get_stat
( const RooAbsData*         data       ,
  Ostap::Math::Statistic&   stat       ,
  const std::string&        expression , 
  const std::string&        selection  ,
  const std::string&        cut_range  ,
  const Ostap::EventIndex   first      ,
  const Ostap::EventIndex   last       , 
  const Ostap::DataType     xmin       ,
  const Ostap::DataType     xmax       ) const 
{
  /// reset the statistic
  stat.reset() ; // reset the statistic
  // 
  // check the ranges
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  if ( xmax <  xmin       ) { return Ostap::StatusCode::SUCCESS ; }
  //
  if ( nullptr == data    ) { return INVALID_DATA        ; }
  //
  const Ostap::EventIndex num_entries = data -> numEntries() ;
  const Ostap::EventIndex the_last    = std::min ( num_entries , last ) ;
  if ( the_last <= first  ) { return Ostap::StatusCode::SUCCESS ; }
  //
  const bool  weighted = data -> isWeighted() ;
  const char* cutrange = cut_range.empty() ? nullptr : cut_range.c_str() ;
  const std::unique_ptr<Ostap::FormulaVar> expr { Ostap::makeFormula ( expression , data       ) } ;
  if ( !expr || !expr->ok()            ) { return INVALID_FORMULA ; }
  //
  const std::unique_ptr<Ostap::FormulaVar> cuts { Ostap::makeFormula ( selection  , data , true ) } ;
  const bool                with_cuts   = cuts && cuts->ok () ;
  //
  Ostap::Utils::ProgressBar bar { the_last - first , m_progress } ; 
  for ( Ostap::EventIndex entry = first ; entry < the_last ; ++entry , ++bar )
    {
      const RooArgSet* vars = data -> get ( entry ) ;
      if ( nullptr == vars )                             { break    ; } // BREAK 
      //
      if ( cutrange && !vars->allInRange ( cutrange ) )  { continue ; } // CONTINUE
      //
      /// data weight?
      const bool keep = weighted ? data-> weight() : 1.0 ;
      if ( !keep ) { continue ; }                                       // CONTINUE
      //
      /// apply cuts:
      const bool cut  = with_cuts ? cuts -> getVal () : 1.0 ;
      if ( !cut ) { continue ; }                                        // CONTINUE  
      //
      const double value = expr -> getVal () ;
      if ( !in_range ( value , xmin , xmax ) ) { continue ; }          // CONTINUE 
      //
      stat.update ( value ) ;
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
  const Ostap::EventIndex   last       , 
  const Ostap::DataType     xmin       ,
  const Ostap::DataType     xmax       ) const 
{
  /// reset the statistic
  stat.reset() ; // reset the statistic
  // 
  // check the ranges
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  if ( xmax <  xmin       ) { return Ostap::StatusCode::SUCCESS ; }
  //
  if ( nullptr == data    ) { return INVALID_DATA    ; }
  //
  const Ostap::EventIndex num_entries = data -> numEntries() ;
  const Ostap::EventIndex the_last    = std::min ( num_entries , last ) ;
  if ( the_last <= first  ) { return Ostap::StatusCode::SUCCESS ; }
  //
  const bool  weighted = data -> isWeighted() ;
  const char* cutrange = cut_range.empty() ? nullptr : cut_range.c_str() ;
  const std::unique_ptr<Ostap::FormulaVar> expr { Ostap::makeFormula ( expression , data       ) } ;
  if ( !expr || !expr->ok()            ) { return INVALID_FORMULA ; }
  //
  const std::unique_ptr<Ostap::FormulaVar> cuts { Ostap::makeFormula ( selection  , data , true ) } ;
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
      if ( !wc ) { continue ; }                                    // CONTINUE  
      // Total: cuts & weight:
      const double weight  = wd *  wc ;
      if ( !weight  || !std::isfinite ( weight ) ) { continue ; }  // CONTINUE        
      //
      const double value = expr -> getVal () ;
      if ( !in_range ( value , xmin , xmax ) ) { continue ; }      // CONTINUE 
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
  const Ostap::EventIndex   last       , 
  const Ostap::DataType     xmin       ,
  const Ostap::DataType     xmax       , 
  const Ostap::DataType     ymin       , 
  const Ostap::DataType     ymax       ) const 
{
  /// reset the statistic
  stat.reset() ; // reset the statistic
  //
  // check the ranges
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  if ( xmax <  xmin       ) { return Ostap::StatusCode::SUCCESS ; }
  if ( ymax <  ymin       ) { return Ostap::StatusCode::SUCCESS ; }
  //
  if ( nullptr == data    ) { return INVALID_DATA    ; }
  //
  const Ostap::EventIndex  num_entries = data -> GetEntries () ;
  const Ostap::EventIndex  the_last    = std::min ( last , num_entries ) ;
  if ( the_last <= first ) { return Ostap::StatusCode::SUCCESS ; }
  //
 /// formulae for exressons
  const Ostap::Formulae                 formulae { data , { expr1 , expr2 } } ;
  const std::unique_ptr<Ostap::Formula> cuts     { Ostap::makeFormula ( selection , data , true ) } ;
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
      const bool cut = with_cuts ? cuts->evaluate() : 1.0 ;
      if ( !cut ) { continue ; }                                         // ATTENTION! 
      //
      formulae.evaluate ( 0 , results [ 0 ] ) ;
      formulae.evaluate ( 1 , results [ 1 ] ) ;
      //
      for ( const double x : results [ 0 ] )
	{ if ( !in_range ( x , xmin , xmax ) ) { continue ; } 
	  for ( const double y : results [ 1 ] )
	    { if ( !in_range ( y , ymin , ymax ) ) { continue ; } 
	      stat.update ( x , y ) ; } }       
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
  const Ostap::EventIndex   last       , 
  const Ostap::DataType     xmin       ,
  const Ostap::DataType     xmax       , 
  const Ostap::DataType     ymin       , 
  const Ostap::DataType     ymax       ) const 
{
  /// reset the statistic
  stat.reset() ; // reset the statistic
  //
  // check the ranges
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  if ( xmax <  xmin       ) { return Ostap::StatusCode::SUCCESS ; }
  if ( ymax <  ymin       ) { return Ostap::StatusCode::SUCCESS ; }
  //
  if ( nullptr == data    ) { return INVALID_DATA    ; }
  //
  const Ostap::EventIndex   num_entries = data -> GetEntries()  ;
  const Ostap::EventIndex   the_last    = std::min ( last , num_entries ) ;
  if ( the_last <= first   ) { return Ostap::StatusCode::SUCCESS ; }
  //
  /// formulae for exressons
  const Ostap::Formulae                 formulae { data , { expr1 , expr2 } } ;
  const std::unique_ptr<Ostap::Formula> cuts     { Ostap::makeFormula ( selection , data , true ) } ;
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
      if ( !weight || !std::isfinite ( weight ) ) { continue ; }         // ATTENTION! 
      //
      formulae.evaluate ( 0 , results [ 0 ] ) ;
      formulae.evaluate ( 1 , results [ 1 ] ) ;
      //
      for ( const double x : results [ 0 ] )
	{ if ( !in_range ( x , xmin , xmax ) ) { continue ; } 
	  for ( const double y : results [ 1 ] )
	    { if ( !in_range ( y , ymin , ymax ) ) { continue ; } 
	      stat.update ( x , y , weight ) ; } } 
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
 *  @attention selection is treated as boolean 
 *  @attention for weighted datasets, the weigth is treated as boolean! 
 
 */
// =============================================================================
Ostap::StatusCode Ostap::StatVar::get_stat
( const RooAbsData*         data       ,
  Ostap::Math::Statistic2&  stat       ,
  const std::string&        expr1      ,
  const std::string&        expr2      ,  
  const std::string&        selection  , 
  const std::string&        cut_range  ,
  const Ostap::EventIndex   first      ,
  const Ostap::EventIndex   last       , 
  const Ostap::DataType     xmin       ,
  const Ostap::DataType     xmax       , 
  const Ostap::DataType     ymin       , 
  const Ostap::DataType     ymax       ) const
{
  /// reset the statistic
  stat.reset() ; // reset the statistic
  //
  // check the ranges
  if ( last <= first       ) { return Ostap::StatusCode::SUCCESS ; }
  if ( xmax <  xmin        ) { return Ostap::StatusCode::SUCCESS ; }
  if ( ymax <  ymin        ) { return Ostap::StatusCode::SUCCESS ; }
  //
  if ( nullptr == data     ) { return INVALID_DATA        ; }
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
  const std::unique_ptr<Ostap::FormulaVar> cuts { Ostap::makeFormula ( selection  , data , true ) } ;
  const bool with_cuts   = cuts && cuts->ok () ;
  //
  Ostap::Utils::ProgressBar bar { the_last - first , m_progress } ; 
  for ( Ostap::EventIndex entry = first ; entry < the_last ; ++entry , ++bar )
    {
      const RooArgSet* vars = data -> get ( entry ) ;
      if ( nullptr == vars )                             { break    ; } // BREAK 
      if ( cutrange && !vars->allInRange ( cutrange ) )  { continue ; } // CONTINUE
      //
      /// data weight?
      const bool keep = weighted ? data-> weight() : 1.0 ;
      if ( !keep ) { continue ; }                                       // CONTINUE
      //
      /// apply cuts:
      const bool cut = with_cuts ? cuts -> getVal () : 1.0 ;
      if ( !cut ) { continue ; }                                   // CONTINUE
      //
      const double x = formulae.evaluate ( 0 ) ;
      if ( !in_range ( x , xmin , xmax ) ) { continue ; }  // CONTINUE 
      //
      const double y = formulae.evaluate ( 1 ) ;  
      if ( !in_range ( y , ymin , ymax ) ) { continue ; }  // CONTINUE 
      //
      stat.update ( x , y ) ;
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
  const Ostap::EventIndex   last       , 
  const Ostap::DataType     xmin       ,
  const Ostap::DataType     xmax       , 
  const Ostap::DataType     ymin       , 
  const Ostap::DataType     ymax       ) const
{
  /// reset the statistic
  stat.reset() ; // reset the statistic
  //
  // check the ranges
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  if ( xmax <  xmin       ) { return Ostap::StatusCode::SUCCESS ; }
  if ( ymax <  ymin       ) { return Ostap::StatusCode::SUCCESS ; }
  //
  if ( nullptr == data    ) { return INVALID_DATA    ; }
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
  const std::unique_ptr<Ostap::FormulaVar> cuts { Ostap::makeFormula ( selection  , data , true ) } ;
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
      if ( !weight || !std::isfinite ( weight ) ) { continue ; }  // CONTINUE        
      //
      const double x = formulae.evaluate ( 0 ) ;
      if ( !in_range ( x , xmin , xmax ) ) { continue ; }  // CONTINUE 
      //
      const double y = formulae.evaluate ( 1 ) ;  
      if ( !in_range ( y , ymin , ymax ) ) { continue ; }  // CONTINUE 
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
  const Ostap::EventIndex   last       , 
  const Ostap::DataType     xmin       ,
  const Ostap::DataType     xmax       , 
  const Ostap::DataType     ymin       , 
  const Ostap::DataType     ymax       , 
  const Ostap::DataType     zmin       , 
  const Ostap::DataType     zmax       ) const
{
  /// reset the statistic
  stat.reset() ; // reset the statistic
  //
  // check the ranges
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  if ( xmax <  xmin       ) { return Ostap::StatusCode::SUCCESS ; }
  if ( ymax <  ymin       ) { return Ostap::StatusCode::SUCCESS ; }
  if ( zmax <  zmin       ) { return Ostap::StatusCode::SUCCESS ; }
  //
  if ( nullptr == data    ) { return INVALID_DATA    ; }
  //
  const Ostap::EventIndex  num_entries = data -> GetEntries () ;
  const Ostap::EventIndex  the_last    = std::min ( last , num_entries ) ;
  if ( the_last <= first ) { return Ostap::StatusCode::SUCCESS ; }
  //
 /// formulae for exressons
  const Ostap::Formulae                 formulae { data , { expr1 , expr2 , expr3 } } ;
  const std::unique_ptr<Ostap::Formula> cuts     { Ostap::makeFormula ( selection , data , true ) } ;
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
      const bool cut = with_cuts ? cuts->evaluate() : 1.0 ;
      if ( !cut ) { continue ; }                                         // ATTENTION! 
      //
      formulae.evaluate ( 0 , results [ 0 ] ) ;
      formulae.evaluate ( 1 , results [ 1 ] ) ;
      formulae.evaluate ( 2 , results [ 2 ] ) ;
      //
      for ( const double x : results [ 0 ] )
        { if ( !in_range ( x , xmin , xmax ) ) { continue ; } // CONTINUE 
          for ( const double y : results [ 1 ] )
            { if ( !in_range ( y , ymin , ymax ) ) { continue ; } // CONTINUE 
              for ( const double z : results [ 2 ] )
                { if ( !in_range ( z , zmin , zmax ) ) { continue ; } // CONTINUE 
                  stat.update ( x , y , z ) ; } } }
      //
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
  const Ostap::EventIndex   last       , 
  const Ostap::DataType     xmin       ,
  const Ostap::DataType     xmax       , 
  const Ostap::DataType     ymin       , 
  const Ostap::DataType     ymax       , 
  const Ostap::DataType     zmin       , 
  const Ostap::DataType     zmax       ) const
{
  /// reset the statistic
  stat.reset() ; // reset the statistic
  //
  // check the ranges
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  if ( xmax <  xmin       ) { return Ostap::StatusCode::SUCCESS ; }
  if ( ymax <  ymin       ) { return Ostap::StatusCode::SUCCESS ; }
  if ( zmax <  zmin       ) { return Ostap::StatusCode::SUCCESS ; }
  //
  if ( nullptr == data    ) { return INVALID_DATA    ; }
  //
  const Ostap::EventIndex   num_entries = data -> GetEntries()  ;
  const Ostap::EventIndex   the_last    = std::min ( last , num_entries ) ;
  if ( the_last <= first   ) { return Ostap::StatusCode::SUCCESS ; }
  //
  /// formulae for exressons
  const Ostap::Formulae                 formulae { data , { expr1 , expr2 , expr3 } } ;
  const std::unique_ptr<Ostap::Formula> cuts     { Ostap::makeFormula ( selection , data , true ) } ;
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
      if ( !weight || !std::isfinite ( weight ) ) { continue ; }       // ATTENTION! 
      //
      formulae.evaluate ( 0 , results [ 0 ] ) ;
      formulae.evaluate ( 1 , results [ 1 ] ) ;
      formulae.evaluate ( 2 , results [ 2 ] ) ;
      //
      for ( const double x : results [ 0 ] )
	{ if ( !in_range ( x , xmin , xmax ) ) { continue ; } // CONTINUE 
	  for ( const double y : results [ 1 ] )
	    { if ( !in_range ( y , ymin , ymax ) ) { continue ; } // CONTINUE 
	      for ( const double z : results [ 2 ] )
		{ if ( !in_range ( z , zmin , zmax ) ) { continue ; } // CONTINUE 
		  stat.update ( x , y , z , weight ) ; } } }
      //
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
 *  @attention selection is treated as boolean 
 *  @attention for weighted datasets, the weigth is treated as boolean! 
 */
// =============================================================================
Ostap::StatusCode Ostap::StatVar::get_stat
( const RooAbsData*         data       ,
  Ostap::Math::Statistic3&  stat       ,
  const std::string&        expr1      ,
  const std::string&        expr2      ,  
  const std::string&        expr3      , 
  const std::string&        selection  , 
  const std::string&        cut_range  ,
  const Ostap::EventIndex   first      ,
  const Ostap::EventIndex   last       , 
  const Ostap::DataType     xmin       ,
  const Ostap::DataType     xmax       , 
  const Ostap::DataType     ymin       , 
  const Ostap::DataType     ymax       , 
  const Ostap::DataType     zmin       , 
  const Ostap::DataType     zmax       ) const
{
  /// reset the statistic
  stat.reset() ; // reset the statistic
  //
  // check the ranges
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  if ( xmax <  xmin       ) { return Ostap::StatusCode::SUCCESS ; }
  if ( ymax <  ymin       ) { return Ostap::StatusCode::SUCCESS ; }
  if ( zmax <  zmin       ) { return Ostap::StatusCode::SUCCESS ; }
  //
  if ( nullptr == data    ) { return INVALID_DATA        ; }
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
  const std::unique_ptr<Ostap::FormulaVar> cuts { Ostap::makeFormula ( selection  , data , true ) } ;
  const bool                with_cuts   = cuts && cuts->ok () ;
  //
  Ostap::Utils::ProgressBar bar { the_last - first , m_progress } ; 
  for ( Ostap::EventIndex entry = first ; entry < the_last ; ++entry , ++bar )
    {
      const RooArgSet* vars = data -> get ( entry ) ;
      if ( nullptr == vars )                             { break    ; } // BREAK 
      //
      if ( cutrange && !vars->allInRange ( cutrange ) )  { continue ; } // CONTINUE
      //
      /// data weight?
      const bool keep = weighted ? data-> weight() : 1.0 ;
      if ( !keep ) { continue ; }                                       // CONTINUE
      //
      // apply cuts:
      const bool cut  = with_cuts ? cuts -> getVal () : 1.0 ;
      if ( !cut ) { continue ; }                                   // CONTINUE  
      //
      const double x = formulae.evaluate ( 0 ) ;
      if ( !in_range ( x , xmin , xmax ) ) { continue ; }  // CONTINUE 
      //
      const double y = formulae.evaluate ( 1 ) ;
      if ( !in_range ( y , ymin , ymax ) ) { continue ; }  // CONTINUE 
      //
      const double z = formulae.evaluate ( 2 ) ;
      if ( !in_range ( z , zmin , zmax ) ) { continue ; }  // CONTINUE 
      //
      stat.update ( x , y , z ) ;
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
  const Ostap::EventIndex   last       , 
  const Ostap::DataType     xmin       ,
  const Ostap::DataType     xmax       , 
  const Ostap::DataType     ymin       , 
  const Ostap::DataType     ymax       , 
  const Ostap::DataType     zmin       , 
  const Ostap::DataType     zmax       ) const
{
  /// reset the statistic
  stat.reset() ; // reset the statistic
  //
  // check the ranges
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  if ( xmax <  xmin       ) { return Ostap::StatusCode::SUCCESS ; }
  if ( ymax <  ymin       ) { return Ostap::StatusCode::SUCCESS ; }
  if ( zmax <  zmin       ) { return Ostap::StatusCode::SUCCESS ; }
  //
  if ( nullptr == data    ) { return INVALID_DATA    ; }
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
  const std::unique_ptr<Ostap::FormulaVar> cuts { Ostap::makeFormula ( selection  , data , true ) } ;
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
      if ( !weight || !std::isfinite ( weight )  ) { continue ; }   // CONTINUE        
      //
      const double x = formulae.evaluate ( 0 ) ;
      if ( !in_range ( x , xmin , xmax ) ) { continue ; }  // CONTINUE 
      //
      const double y = formulae.evaluate ( 1 ) ;
      if ( !in_range ( y , ymin , ymax ) ) { continue ; }  // CONTINUE 
      //
      const double z = formulae.evaluate ( 2 ) ;
      if ( !in_range ( z , zmin , zmax ) ) { continue ; }  // CONTINUE 
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
  const Ostap::EventIndex   last       , 
  const Ostap::DataType     xmin       ,
  const Ostap::DataType     xmax       , 
  const Ostap::DataType     ymin       , 
  const Ostap::DataType     ymax       , 
  const Ostap::DataType     zmin       , 
  const Ostap::DataType     zmax       , 
  const Ostap::DataType     tmin       , 
  const Ostap::DataType     tmax       ) const
{
  /// reset the statistic
  stat.reset() ; // reset the statistic
  //
  // check the ranges
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  if ( xmax <  xmin       ) { return Ostap::StatusCode::SUCCESS ; }
  if ( ymax <  ymin       ) { return Ostap::StatusCode::SUCCESS ; }
  if ( zmax <  zmin       ) { return Ostap::StatusCode::SUCCESS ; }
  if ( tmax <  tmin       ) { return Ostap::StatusCode::SUCCESS ; }
  //
  if ( nullptr == data    ) { return INVALID_DATA    ; }
  //
  const Ostap::EventIndex  num_entries = data -> GetEntries () ;
  const Ostap::EventIndex  the_last    = std::min ( last , num_entries ) ;
  if ( the_last <= first ) { return Ostap::StatusCode::SUCCESS ; }
  //
 /// formulae for exressons
  const Ostap::Formulae                 formulae { data , { expr1 , expr2 , expr3 , expr4 } } ; 
  const std::unique_ptr<Ostap::Formula> cuts     { Ostap::makeFormula ( selection , data , true ) } ;
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
      const bool cut = with_cuts ? cuts->evaluate() : 1.0 ;
      if ( !cut ) { continue ; }                                         // ATTENTION! 
      //
      formulae.evaluate ( 0 , results [ 0 ] ) ;
      formulae.evaluate ( 1 , results [ 1 ] ) ;
      formulae.evaluate ( 2 , results [ 2 ] ) ;
      formulae.evaluate ( 3 , results [ 3 ] ) ;
      //
      for ( const double x : results [ 0 ] )
	{ if ( !in_range ( x , xmin , xmax ) ) { continue ; } // CONTINUE 
	  for ( const double y : results [ 1 ] )
	    { if ( !in_range ( y , ymin , ymax ) ) { continue ; } // CONTINUE 
	      for ( const double z : results [ 2 ] )
		{ if ( !in_range ( z , zmin , zmax ) ) { continue ; } // CONTINUE
		  for ( const double t : results [ 3 ] )
		    { if ( !in_range ( t , tmin , tmax ) ) { continue ; } // CONTINUE
		      stat.update ( x , y , z , t ) ; } } } } 
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
  const Ostap::EventIndex   last       , 
  const Ostap::DataType     xmin       ,
  const Ostap::DataType     xmax       , 
  const Ostap::DataType     ymin       , 
  const Ostap::DataType     ymax       , 
  const Ostap::DataType     zmin       , 
  const Ostap::DataType     zmax       , 
  const Ostap::DataType     tmin       , 
  const Ostap::DataType     tmax       ) const
{
  /// reset the statistic
  stat.reset() ; // reset the statistic
  //
  // check the ranges
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  if ( xmax <  xmin       ) { return Ostap::StatusCode::SUCCESS ; }
  if ( ymax <  ymin       ) { return Ostap::StatusCode::SUCCESS ; }
  if ( zmax <  zmin       ) { return Ostap::StatusCode::SUCCESS ; }
  if ( tmax <  tmin       ) { return Ostap::StatusCode::SUCCESS ; }
  //
  if ( nullptr == data    ) { return INVALID_DATA    ; }
  //
  const Ostap::EventIndex   num_entries = data -> GetEntries()  ;
  const Ostap::EventIndex   the_last    = std::min ( last , num_entries ) ;
  if ( the_last <= first   ) { return Ostap::StatusCode::SUCCESS ; }
  //
  /// formulae for exressons
  const Ostap::Formulae                 formulae { data , { expr1 , expr2 , expr3 , expr4 } } ; 
  const std::unique_ptr<Ostap::Formula> cuts     { Ostap::makeFormula ( selection , data , true ) } ;
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
      if ( !weight || !std::isfinite ( weight ) ) { continue ; }       // ATTENTION! 
      //
      formulae.evaluate ( 0 , results [ 0 ] ) ;
      formulae.evaluate ( 1 , results [ 1 ] ) ;
      formulae.evaluate ( 2 , results [ 2 ] ) ;
      formulae.evaluate ( 3 , results [ 3 ] ) ;
      //
      for ( const double x : results [ 0 ] )
	{ if ( !in_range ( x , xmin , xmax ) ) { continue ; } // CONTINUE 
	  for ( const double y : results [ 1 ] )
	    { if ( !in_range ( y , ymin , ymax ) ) { continue ; } // CONTINUE 
	      for ( const double z : results [ 2 ] )
		{ if ( !in_range ( z , zmin , zmax ) ) { continue ; } // CONTINUE
		  for ( const double t : results [ 3 ] )
		    { if ( !in_range ( t , tmin , tmax ) ) { continue ; } // CONTINUE
		      stat.update ( x , y , z , t , weight ) ; } } } } 
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
 *  @attention selection is treated as boolean 
 *  @attention for weighted datasets, the weigth is treated as boolean!  
 */
// =============================================================================
Ostap::StatusCode Ostap::StatVar::get_stat
( const RooAbsData*         data       ,
  Ostap::Math::Statistic4&  stat       ,
  const std::string&        expr1      ,
  const std::string&        expr2      ,  
  const std::string&        expr3      ,
  const std::string&        expr4      ,  
  const std::string&        selection  , 
  const std::string&        cut_range  ,
  const Ostap::EventIndex   first      ,
  const Ostap::EventIndex   last       , 
  const Ostap::DataType     xmin       ,
  const Ostap::DataType     xmax       , 
  const Ostap::DataType     ymin       , 
  const Ostap::DataType     ymax       , 
  const Ostap::DataType     zmin       , 
  const Ostap::DataType     zmax       , 
  const Ostap::DataType     tmin       , 
  const Ostap::DataType     tmax       ) const
{
  /// reset the statistic
  stat.reset() ; // reset the statistic
  //
  // check the ranges
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  if ( xmax <  xmin       ) { return Ostap::StatusCode::SUCCESS ; }
  if ( ymax <  ymin       ) { return Ostap::StatusCode::SUCCESS ; }
  if ( zmax <  zmin       ) { return Ostap::StatusCode::SUCCESS ; }
  if ( tmax <  tmin       ) { return Ostap::StatusCode::SUCCESS ; }
  //
  if ( nullptr == data    ) { return INVALID_DATA        ; }
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
  const std::unique_ptr<Ostap::FormulaVar> cuts { Ostap::makeFormula ( selection  , data , true ) } ;
  const bool                with_cuts   = cuts && cuts->ok () ;
  //
  Ostap::Utils::ProgressBar bar { the_last - first , m_progress } ; 
  for ( Ostap::EventIndex entry = first ; entry < the_last ; ++entry , ++bar )
    {
      const RooArgSet* vars = data -> get ( entry ) ;
      if ( nullptr == vars )                             { break    ; } // BREAK 
      //
      if ( cutrange && !vars->allInRange ( cutrange ) )  { continue ; } // CONTINUE
      //
      /// data weight?
      const bool keep = weighted ? data-> weight() : 1.0 ;
      if ( !keep ) { continue ; }                                       // CONTINUE
      //      
      // apply cuts:
      const bool  cut = with_cuts ? cuts -> getVal () : 1.0 ;
      if ( !cut ) { continue ; }                                   // CONTINUE  
      //
      const double x = formulae.evaluate ( 0 ) ;
      if ( !in_range ( x , xmin , xmax ) ) { continue ; }  // CONTINUE 
      //
      const double y = formulae.evaluate ( 1 ) ;
      if ( !in_range ( y , ymin , ymax ) ) { continue ; }  // CONTINUE 
      //
      const double z = formulae.evaluate ( 2 ) ;
      if ( !in_range ( z , zmin , zmax ) ) { continue ; }  // CONTINUE
      //
      const double t = formulae.evaluate ( 3 ) ;
      if ( !in_range ( t , tmin , tmax ) ) { continue ; }  // CONTINUE 
      //
      stat.update ( x , y , z , t ) ;
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
  const Ostap::EventIndex   last       , 
  const Ostap::DataType     xmin       ,
  const Ostap::DataType     xmax       , 
  const Ostap::DataType     ymin       , 
  const Ostap::DataType     ymax       , 
  const Ostap::DataType     zmin       , 
  const Ostap::DataType     zmax       , 
  const Ostap::DataType     tmin       , 
  const Ostap::DataType     tmax       ) const
{
  /// reset the statistic
  stat.reset() ; // reset the statistic
  //
  // check the ranges
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  if ( xmax <  xmin       ) { return Ostap::StatusCode::SUCCESS ; }
  if ( ymax <  ymin       ) { return Ostap::StatusCode::SUCCESS ; }
  if ( zmax <  zmin       ) { return Ostap::StatusCode::SUCCESS ; }
  if ( tmax <  tmin       ) { return Ostap::StatusCode::SUCCESS ; }
  //
  if ( nullptr == data    ) { return INVALID_DATA    ; }
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
  const std::unique_ptr<Ostap::FormulaVar> cuts { Ostap::makeFormula ( selection  , data , true ) } ;
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
      if ( !weight || !std::isfinite ( weight )  ) { continue ; } // CONTINUE        
      //
      const double x = formulae.evaluate ( 0 ) ;
      if ( !in_range ( x , xmin , xmax ) ) { continue ; }  // CONTINUE 
      //
      const double y = formulae.evaluate ( 1 ) ;
      if ( !in_range ( y , ymin , ymax ) ) { continue ; }  // CONTINUE 
      //
      const double z = formulae.evaluate ( 2 ) ;
      if ( !in_range ( z , zmin , zmax ) ) { continue ; }  // CONTINUE
      //
      const double t = formulae.evaluate ( 3 ) ;
      if ( !in_range ( t , tmin , tmax ) ) { continue ; }  // CONTINUE 
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
  const Ostap::EventIndex the_last    = std::min ( last , num_entries ) ;
  if ( the_last <= first   ) { return Ostap::StatusCode::SUCCESS ; }
  //
  /// formulae for exressons
  const Ostap::Formulae                 formulae { data , expressions } ;  
  const std::unique_ptr<Ostap::Formula> cuts     { Ostap::makeFormula ( selection , data , true ) } ;
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
          for ( const double value : results )
            { if ( std::isfinite ( value ) ) { result[i].add ( value ) ; } } 
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
  const Ostap::EventIndex the_last    = std::min ( last , num_entries ) ;
  if ( the_last <= first   ) { return Ostap::StatusCode::SUCCESS ; }
  //
  /// formulae for exressons
  const Ostap::Formulae                 formulae { data , expressions } ;  
  const std::unique_ptr<Ostap::Formula> cuts     { Ostap::makeFormula ( selection , data , true ) } ;  
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
      if ( !weight || !std::isfinite ( weight ) ) { continue  ; }     // CONTNUE 
      //
      for ( std::size_t i = 0 ; i < N ; ++i ) 
	{
	  formulae.evaluate ( i , results ) ;
	  for ( const double value : results )
	    { if ( std::isfinite ( value ) ) {  result[i].add ( value , weight ) ; } }
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
  const Ostap::EventIndex the_last    = std::min ( last , num_entries ) ;
  if ( the_last <= first   ) { return Ostap::StatusCode::SUCCESS ; }
  //
  /// formulae for exressons
  const Ostap::FormulaVars                 formulae { data , expressions } ;    
  const std::unique_ptr<Ostap::FormulaVar> cuts     { Ostap::makeFormula ( selection , data , true ) } ;
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
      if ( !weight || !std::isfinite ( weight ) ) { continue ; } // CONTINUE        
      //
      for ( std::size_t i = 0 ; i < N ; ++i ) 
        {
          const double value = formulae.evaluate ( i ) ;
          if ( !std::isfinite ( value ) ) { continue ; } // CONTINUE 
          result [ i ].add ( value , weight ) ;
        }
      //
    }
  //
  return Ostap::StatusCode::SUCCESS ; 
}
// ==========================================================================

// =========================================================================
/// Is there at leats one good entry in dataset ?
// =========================================================================

// =========================================================================
/*  Is there at leats one good entry in dataset ?
 *  @param data       (INPUT) input data 
 *  @param selecttion (INPUT) selection criteria 
 *  @param first      (INPUT) the first event to process (inclusibe) 
 *  @param last       (INPUT) the last event to process (exclusive) 
 */
// =========================================================================
bool Ostap::StatVar::hasEntry 
( TTree*                    data                           , 
  const std::string         selection  , 
  const Ostap::EventIndex   first      ,
  const Ostap::EventIndex   last       ) const 
{
  if ( !data ) { return false ; } // NB!
  //
  const Ostap::EventIndex num_entries =  data->GetEntries() ;
  const Ostap::EventIndex the_last    = std::min ( last , num_entries ) ;
  if ( the_last <= first  ) { return false ; }
  //
  /// no cuts ? 
  if ( Ostap::trivial ( selection ) ) { return first < the_last ; }
  //
  /// formulae for exression
  const std::unique_ptr<Ostap::Formula> cuts { Ostap::makeFormula ( selection , data , true ) } ;
  const bool with_cuts = cuts && cuts->ok() ;
  //
  if  ( !with_cuts ) { return first < the_last ; }
  //
  Ostap::Utils::Notifier    notify ( data , cuts.get() ) ;
  Ostap::Utils::ProgressBar bar    ( the_last - first , m_progress ) ;
  for ( Ostap::EventIndex   entry = first ; entry < the_last  ; ++entry , ++bar )
    {
      //
      const long ievent = data -> GetEntryNumber ( entry ) ;
      if ( 0 > ievent                    ) { return false ; } // RETURN
      if ( 0 > data->LoadTree ( ievent ) ) { return false ; } // RETURN
      //
      const bool cut = with_cuts ? cuts ->evaluate() : 1.0 ;
      if ( !cut ) { continue  ; }                        
      //
      return true ;                                            // RETURN 
    }
  return false ;                                               // RETURN
}
// ====================================================================================
/*  Is there at leats one good entry in dataset ?
 *  @param data       (INPUT) input data 
 *  @param selecttion (INPUT) selection criteria 
 *  @param first      (INPUT) the first event to process (inclusibe) 
 *  @param last       (INPUT) the last event to process (exclusive) 
 */
// ====================================================================================
bool Ostap::StatVar::hasEntry 
( const RooAbsData*         data       , 
  const std::string&        selection  ,
  const std::string&        cut_range  , 
  const Ostap::EventIndex   first      ,
  const Ostap::EventIndex   last       ) const
{
  if ( !data ) { return false ; } // NB!                                 // RETURN 
  //
  const Ostap::EventIndex num_entries =  data->numEntries() ;
  const Ostap::EventIndex the_last    = std::min ( last , num_entries ) ;
  if ( the_last <= first   ) { return false ; ; }                         // RETURN 
  //
  /// formulae for exressons
  const std::unique_ptr<Ostap::FormulaVar> cuts     { Ostap::makeFormula ( selection , data , true ) } ;
  //
  const char* cutrange = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  const bool  weighted = data->isWeighted() ;
  //
  // No need to run the loop 
  if ( !cuts && !cutrange && !weighted ) { return first < the_last ; }   // RETURN
  //
  // start the loop
  Ostap::Utils::ProgressBar bar ( the_last - first , m_progress ) ;
  for ( Ostap::EventIndex   entry = first ; entry < the_last ; ++entry , ++bar )
    {
      //
      const RooArgSet* vars = data->get( entry ) ;
      if ( nullptr == vars  )                           { return false ; } // RETURN
      if ( cutrange && !vars->allInRange ( cutrange ) ) { continue     ; } // CONTINUE    
      //
      // apply weight:
      const double wd = weighted  ? data->weight()        : 1.0 ;
      if ( !wd ) { continue ; }                                   // CONTINUE    
      // apply cuts:
      const double wc = cuts ? cuts -> getVal() : 1.0 ;
      if ( !wc ) { continue ; }                                   // CONTINUE  
      // cuts & weight:
      const double weight  = wd *  wc ;
      if ( !weight || !std::isfinite ( weight ) ) { continue ; }  // CONTINUE        
      //
      return true ;                                               // RETURN 
    }
  return false ;
}
// =========================================================================


// =========================================================================
/*  Number of good entries 
 *  @param data       (INPUT) input data 
 *  @param selecttion (INPUT) selection criteria 
 *  @param first      (INPUT) the first event to process (inclusibe) 
 *  @param last       (INPUT) the last event to process (exclusive) 
 */
// =========================================================================
Ostap::EventIndex Ostap::StatVar::size 
( TTree*                    data                           , 
  const std::string&        selection  , 
  const Ostap::EventIndex   first      ,
  const Ostap::EventIndex   last       ) const 
{
  if ( !data ) { return 0  ; } // NB!
  //
  const Ostap::EventIndex num_entries =  data->GetEntries() ;
  const Ostap::EventIndex the_last    = std::min ( last , num_entries ) ;
  if ( the_last <= first  ) { return 0 ; }
  //
  /// no cuts ? 
  if ( Ostap::trivial ( selection ) ) { return the_last - first  ; }
  //
  /// formulae for exression
  const std::unique_ptr<Ostap::Formula> cuts { Ostap::makeFormula ( selection , data , true ) } ;
  const bool with_cuts = cuts && cuts->ok() ;
  //
  if  ( !with_cuts ) { return the_last - first ; }
  //
  Ostap::Utils::Notifier    notify ( data , cuts.get() ) ;
  Ostap::EventIndex         result { 0 } ; 
  Ostap::Utils::ProgressBar bar    ( the_last - first , m_progress ) ;
  for ( Ostap::EventIndex   entry = first ; entry < the_last  ; ++entry , ++bar )
    {
      //
      const long ievent = data -> GetEntryNumber ( entry ) ;
      if ( 0 > ievent                    ) { return result ; } // RETURN
      if ( 0 > data->LoadTree ( ievent ) ) { return result ; } // RETURN
      //
      const bool cut = with_cuts ? cuts ->evaluate() : 1.0 ;
      if ( !cut ) { continue  ; }                        
      //
      ++result ;
    }
  return result ;                                               // RETURN
}
// ====================================================================================
/*  Is there at leats one good entry in dataset ?
 *  @param data       (INPUT) input data 
 *  @param selecttion (INPUT) selection criteria 
 *  @param first      (INPUT) the first event to process (inclusibe) 
 *  @param last       (INPUT) the last event to process (exclusive) 
 */
// ====================================================================================
Ostap::EventIndex Ostap::StatVar::size 
( const RooAbsData*         data       , 
  const std::string&        selection  ,
  const std::string&        cut_range  , 
  const Ostap::EventIndex   first      ,
  const Ostap::EventIndex   last       ) const
{
  if ( !data ) { return 0  ; } // NB!                                 // RETURN 
  //
  const Ostap::EventIndex num_entries =  data->numEntries() ;
  const Ostap::EventIndex the_last    = std::min ( last , num_entries ) ;
  if ( the_last <= first   ) { return 0  ; ; }                         // RETURN 
  //
  /// formulae for exressons
  const std::unique_ptr<Ostap::FormulaVar> cuts     { Ostap::makeFormula ( selection , data , true ) } ;
  //
  const char* cutrange = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  const bool  weighted = data->isWeighted() ;
  //
  // No need to run the loop 
  if ( !cuts && !cutrange && !weighted ) { return the_last - first ; }   // RETURN
  //
  Ostap::EventIndex result { 0 } ;
  // start the loop  
  Ostap::Utils::ProgressBar bar ( the_last - first , m_progress ) ;
  for ( Ostap::EventIndex   entry = first ; entry < the_last ; ++entry , ++bar )
    {
      //
      const RooArgSet* vars = data->get( entry ) ;
      if ( nullptr == vars  )                           { return  result ; } // RETURN
      if ( cutrange && !vars->allInRange ( cutrange ) ) { continue       ; } // CONTINUE    
      //
      // apply weight:
      const double wd = weighted  ? data->weight()        : 1.0 ;
      if ( !wd ) { continue ; }                                   // CONTINUE    
      // apply cuts:
      const double wc = cuts ? cuts -> getVal() : 1.0 ;
      if ( !wc ) { continue ; }                                   // CONTINUE  
      // cuts & weight:
      const double weight  = wd *  wc ;
      if ( !weight || !std::isfinite ( weight ) ) { continue ; }  // CONTINUE        
      //
      ++result ;
    }
  return result;
}
// =========================================================================

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
( TTree*                  data       , 
  const std::string&      expression , 
  const Ostap::EventIndex first      ,
  const Ostap::EventIndex last       , 
  const Ostap::DataType   xmin       ,
  const Ostap::DataType   xmax       ) const 
{ return statVar_cut ( data , expression , std::string () , first , last , xmin , xmax ) ; }
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
  const Ostap::EventIndex last       , 
  const Ostap::DataType   xmin       ,
  const Ostap::DataType   xmax       ) const 
{
  ///
  Ostap::StatEntity result {} ;
  const Ostap::StatusCode sc = get_stat ( data       ,
                                          result     , 
                                          expression ,
                                          selection  , 
                                          first      ,
                                          last       ,
                                          xmin       ,
                                          xmax       ) ;
  // check status code 
  Ostap::Assert ( sc.isSuccess ()  ,
                  "Error from Ostap::StatVar::get_stat " ,
                  "Ostap::StatVar::statVar_cut" , sc , __FILE__ , __LINE__ ) ;
  //
  return result ;
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
  const Ostap::EventIndex last       , 
  const Ostap::DataType   xmin       ,
  const Ostap::DataType   xmax       ) const 
{
  ///
  Ostap::WStatEntity result {} ;
  const Ostap::StatusCode sc = get_stat ( data       ,
                                          result     , 
                                          expression ,
                                          selection  , 
                                          first      ,
                                          last       ,
                                          xmin       ,
                                          xmax       ) ;
  // check status code 
  Ostap::Assert ( sc.isSuccess ()  ,
                  "Error from Ostap::StatVar::get_stat " ,
                  "Ostap::StatVar::statVar" , sc , __FILE__ , __LINE__ ) ;
  //
  return result ;
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
  const Ostap::EventIndex last       ,
  const Ostap::DataType   xmin       ,
  const Ostap::DataType   xmax       ) const 
{
  ///
  Ostap::WStatEntity result {} ;
  const Ostap::StatusCode sc = get_stat ( data       ,
					  result     , 
					  expression ,
					  selection  ,
					  cut_range  , 
					  first      ,
					  last       ,
					  xmin       ,
					  xmax       ) ;
  // check status code 
  Ostap::Assert ( sc.isSuccess ()  ,
		  "Error from Ostap::StatVar::get_stat " ,
		  "Ostap::StatVar::statVar" , sc , __FILE__ , __LINE__ ) ;
  //
  return result ;
}

// ============================================================================
// Covariances for two variables 
// ============================================================================

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
Ostap::StatusCode
Ostap::StatVar::statCov
( TTree*                   data      ,
  Ostap::Math::Covariance& stat      , 
  const std::string&       exp1      ,
  const std::string&       exp2      ,
  const std::string&       selection ,
  const Ostap::EventIndex  first     ,
  const Ostap::EventIndex  last      , 
  const Ostap::DataType    xmin      ,
  const Ostap::DataType    xmax      , 
  const Ostap::DataType    ymin      ,
  const Ostap::DataType    ymax      ) const 
{
  return get_stat ( data  ,
		    stat  , 
		    exp1  , 
		    exp2  , 
		    selection ,
		    first ,
		    last  ,
		    xmin  ,
		    xmax  ,
		    ymin  ,
		    ymax  ) ;
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
Ostap::StatusCode
Ostap::StatVar::statCov
( TTree*                    data      ,
  Ostap::Math::WCovariance& stat      , 
  const std::string&        exp1      ,
  const std::string&        exp2      ,
  const std::string&        selection ,
  const Ostap::EventIndex   first     ,
  const Ostap::EventIndex   last      , 
  const Ostap::DataType     xmin      ,
  const Ostap::DataType     xmax      , 
  const Ostap::DataType     ymin      ,
  const Ostap::DataType     ymax      ) const 
{
  return get_stat ( data  ,
		    stat  , 
		    exp1  , 
		    exp2  , 
		    selection ,
		    first ,
		    last  ,
		    xmin  ,
		    xmax  ,
		    ymin  ,
		    ymax  ) ;
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
Ostap::StatusCode 
Ostap::StatVar::statCov
( const RooAbsData*         data      , 
  Ostap::Math::WCovariance& stat      , 
  const std::string&        exp1      , 
  const std::string&        exp2      , 
  const std::string&        selection , 
  const std::string&        cut_range , 
  const Ostap::EventIndex   first     ,
  const Ostap::EventIndex   last      , 
  const Ostap::DataType     xmin      ,
  const Ostap::DataType     xmax      , 
  const Ostap::DataType     ymin      ,
  const Ostap::DataType     ymax      ) const
{
  return get_stat ( data      ,
                    stat      , 
                    exp1      , 
                    exp2      , 
                    selection ,
                    cut_range , 
                    first     ,
                    last      ,
                    xmin      ,
                    xmax      ,
                    ymin      ,
                    ymax      ) ;
}
// ============================================================================

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
  if ( N < 2               ) { return INVALID_SIZE ; } 
  if ( nullptr == data     ) { return INVALID_TREE ; }
  //
  const Ostap::EventIndex num_entries =  data->GetEntries() ;
  const Ostap::EventIndex the_last    = std::min ( last , num_entries ) ;
  if ( the_last <= first   ) { return Ostap::StatusCode::SUCCESS ; }
  //
  /// formulae for exressons
  const Ostap::Formulae                 formulae { data , expressions } ;  
  const std::unique_ptr<Ostap::Formula> cuts     { Ostap::makeFormula ( selection , data , true ) } ;    
  //
  Ostap::Utils::Notifier    notify  ( formulae.begin() , formulae.end() , cuts.get() , data ) ;
  // 
  const bool                with_cuts = cuts && cuts->ok() ;
  //
  typedef std::vector<double> RESULT  ;
  typedef std::vector<RESULT> RESULTS ;
  RESULTS results { N } ;
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
      const bool cut = with_cuts ? cuts ->evaluate() : 1.0 ;
      if ( !cut ) { continue  ; }
      //
      results.resize ( N ) ;
      for ( std::size_t i = 0 ; i < N ; ++i )  { formulae.evaluate ( i , results [ i ] ) ; }
      //
      Ostap::Combiner_<RESULT> combiner { results.begin() , results.end() };
      //
      while ( combiner.valid() ) 
        {
          // get the current combination 
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
Ostap::StatusCode  Ostap::StatVar::statCov
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
  if ( N < 2               ) { return INVALID_SIZE ; } 
  if ( nullptr == data     ) { return INVALID_TREE ; }
  //
  const Ostap::EventIndex num_entries =  data->GetEntries() ;
  const Ostap::EventIndex the_last    = std::min ( last , num_entries ) ;
  if ( the_last <= first   ) { return Ostap::StatusCode::SUCCESS ; }
  //
  /// formulae for exressons
  const Ostap::Formulae                 formulae { data , expressions } ;  
  const std::unique_ptr<Ostap::Formula> cuts     { Ostap::makeFormula ( selection , data , true ) } ;    
  // 
  Ostap::Utils::Notifier    notify  ( formulae.begin() , formulae.end() , cuts.get() , data ) ;
  // 
  const bool                with_cuts = cuts && cuts->ok() ;
  //
  typedef std::vector<double> RESULT  ;
  typedef std::vector<RESULT> RESULTS ;
  RESULTS results { N } ;
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
      if ( !weight || !std::isfinite ( weight ) ) { continue  ; }
      //
      results.resize ( N ) ;
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
Ostap::StatusCode Ostap::StatVar::statCov
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
  if ( N < 2               ) { return INVALID_SIZE ; } 
  if ( nullptr == data     ) { return INVALID_DATA ; }
  //
  const Ostap::EventIndex num_entries =  data -> numEntries () ;
  const Ostap::EventIndex the_last    = std::min ( last , num_entries ) ;
  if ( the_last <= first   ) { return Ostap::StatusCode::SUCCESS ; }
  //
  /// formulae for exressons
  const Ostap::FormulaVars                 formulae { data , expressions } ;  
  const std::unique_ptr<Ostap::FormulaVar> cuts     { Ostap::makeFormula ( selection , data , true ) } ;  
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
      if ( !weight || !std::isfinite ( weight ) )       { continue ; } // CONTINUE  
      //
      // evaluate all expressions
      result.resize ( N ) ;
      for ( std::size_t i = 0 ; i < N ; ++i )
        { result [ i ] = formulae.evaluate ( i ) ; }
      //
      stats.add ( result , weight ) ;
    } 
  //
  return Ostap::StatusCode::SUCCESS ; 
}
// ============================================================================

// ============================================================================
// get tanle fomr RooAbsData 
// ============================================================================
/* Get data table from RooAbsData 
 *  @param data       (input)  data 
 *  @param table      (UPDATE) table 
 *  @param selection  (INPUT)  selection/cut (treated as weight!)
 *  @param cut_range  (INPUT)  if non empty: use events only from this cut-range
 *  @param first      (INPUT)  the first event to process (inclusibe) 
 *  @param last       (INPUT)  the last event to process (exclusive) 
 *  @return status code  
 */
// ============================================================================
Ostap::StatusCode Ostap::StatVar::get_table
( const RooAbsData*         data       ,
  Ostap::StatVar::Table&    table      ,
  const std::string&        selection  ,
  const std::string&        cut_range  ,
  const Ostap::EventIndex   first      ,
  const Ostap::EventIndex   last       ) const
{
  // clear data dable ;
  for ( Table::iterator i = table.begin () ; table.end() != i ; ++i ) { i->second.clear() ; }
  //
  if  ( nullptr == data     ) { table.clear() ; return INVALID_DATA ; }
  //
  const Ostap::EventIndex num_entries =  data -> numEntries () ;
  const Ostap::EventIndex the_last    = std::min ( last , num_entries ) ;
  if ( the_last <= first   ) { table.clear() ; return Ostap::StatusCode::SUCCESS ; }
  //
  std::vector<std::string> expressions{} ; expressions.reserve ( table.size() + 1 ) ;
  for ( Table::const_iterator i = table.begin () ; table.end() != i ; ++i )
    { expressions.push_back ( i->first ) ; }
  //
  /// formulae for expressons
  const Ostap::FormulaVars                 formulae { data , expressions } ;  
  const std::unique_ptr<Ostap::FormulaVar> cuts     { Ostap::makeFormula ( selection , data , true ) } ;  
  ///
  const bool  with_cuts = cuts && cuts->ok() ;
  const char* cutrange  = cut_range.empty() ? nullptr : cut_range.c_str() ;
  const bool  weighted  = data->isWeighted() ;
  //
  /// the name of weight variable
  const std::string wname   { weighted ? Ostap::Utils::getWeight ( data ) :  "" } ;  
  const bool        wsep =  ( weighted || with_cuts ) ;
  ///
  const std::size_t N = formulae.size () ;
  /// 
  typedef std::vector<Column>  TABLE ;
  TABLE results { wsep ? N + 1 : N } ;
  //
  /// reserved the space 
  for ( auto& r : results ) { r.reserve ( ( last - first ) / 3 ) ; } 
  ///
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
      if ( !weight || !std::isfinite ( weight )  ) { continue ; }   // CONTINUE        
      //
      /// fill the table 
      for ( std::size_t i = 0 ; i < N ; ++i )
        { results [ i ].push_back ( formulae.evaluate ( i ) ) ; }
      /// add the combined weight 
      if ( wsep ) { results.back().push_back ( weight  ) ; }  
      //     
    }
  //
  // move data to the output table
  std::size_t i = 0 ;
  for ( Table::iterator c = table.begin ()  ; table.end() != c ; ++c , ++i )
    { std::swap ( c->second , results [ i ] ) ; }
  //
  /// add the weight/cuts column 
  if ( wsep ) { table [ wname ] = results.back () ;  }
  //
  return Ostap::StatusCode::SUCCESS ;
}
// ============================================================================
//                                                                      The END
// ============================================================================
