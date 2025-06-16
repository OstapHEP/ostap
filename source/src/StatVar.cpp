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
#include "RVersion.h"
#include "TTree.h"
#include "TMatrixTSym.h"
#include "RooDataSet.h"
#include "RooAbsReal.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Formula.h"
#include "Ostap/Iterator.h"
#include "Ostap/Notifier.h"
#include "Ostap/MatrixUtils.h"
#include "Ostap/StatVar.h"
#include "Ostap/FormulaVar.h"
#include "Ostap/P2Quantile.h"
#include "Ostap/Moments.h"
#include "Ostap/GetWeight.h"
#include "Ostap/Moments.h"
#include "Ostap/ECDF.h"
#include "Ostap/DataFrameUtils.h"
#include "Ostap/ProgressBar.h"
// ============================================================================
// Local
// ============================================================================
#include "OstapDataFrame.h"
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
    const RooAbsData&  data                 , 
    const bool         allow_empty = false  , 
    const bool         allow_null  = false  ) 
  { 
    if ( allow_empty && expression.empty() ) { return nullptr ; }  // RETURN!
    //
    const RooArgSet*  aset = data.get() ;
    //
    if ( allow_null && nullptr == aset     ) { return nullptr ; }
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
} //                                                 end of anonymous namespace
//  ============================================================================
// constructor with ProgressBar configuration 
// =============================================================================
Ostap::StatVar::StatVar
( const Ostap::Utils::ProgressConf& progress)
  : m_progress ( progress )
{}
// ============================================================================
// Helper structures 
// ============================================================================
Ostap::StatVar::Interval::Interval
(  const double l ,
   const double h )
  : low  ( std::min ( l , h ) )
  , high ( std::max ( l , h ) )
{}
// ============================================================================
Ostap::StatVar::Quantile::Quantile
( const double      q ,
  const std::size_t n )
  : quantile ( q ) 
  , nevents  ( n ) 
{}
// ============================================================================
Ostap::StatVar::Quantiles::Quantiles
( const std::vector<double>& q ,  
  const std::size_t          n )
  : quantiles ( q ) 
  , nevents   ( n ) 
{}
// ============================================================================
Ostap::StatVar::QInterval::QInterval
( const Ostap::StatVar::Interval& i ,
  const std::size_t               n )
  : interval ( i ) 
  , nevents  ( n ) 
{}
// =============================================================================


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
  // 
  if ( nullptr == data    ) { return INVALID_DATA    ; }
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  //
  Ostap::Formula formula ( expression , data ) ;
  if ( !formula.ok ()     ) { return INVALID_FORMULA ; }
  //
  std::unique_ptr<Ostap::Formula> cuts { nullptr } ;
  if  ( !trivial ( selection ) ) 
    { 
      cuts = std::make_unique<Ostap::Formula> ( selection , data ) ; 
      if ( !cuts || !cuts->ok() ) { return INVALID_SELECTION ; }
    }
  //
  const bool               with_cuts   = cuts && cuts->ok () ;
  const Ostap::EventIndex  num_entries = data -> GetEntries () ;
  const Ostap::EventIndex  the_last    = std::min ( last , num_entries ) ;
  if ( the_last <= first ) { return Ostap::StatusCode::SUCCESS ; }
  //
  Ostap::Utils::ProgressBar bar    ( the_last - first , m_progress ) ;
  Ostap::Utils::Notifier    notify ( data , &formula , cuts.get () ) ;
  //
  std::vector<double>       results {} ;
  for ( Ostap::EventIndex   entry = first ; entry < the_last  ; ++entry , +bar )
    {
      //
      const long ievent = data -> GetEntryNumber ( entry ) ;
      if ( 0 > ievent ) { return INVALID_ENTRY ; }  // RETURN 
      //
      const ievent      = data -> LoadTree ( ievent ) ;
      if ( 0 > ievent ) { return INVALID_EVENT ; }  // RETURN 
      //
      const double w = with_cuts ? cuts->evaluate() : 1.0 ;
      //
      if ( !w ) { continue ; }                      // ATTENTION! 
      //
      formula.evaluate ( results ) ;
      for  ( const double r : results ) { moment.update ( r ) ; }
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
  // 
  if ( nullptr == data    ) { return INVALID_DATA    ; }
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  //
  Ostap::Formula formula ( expression , data ) ;
  if ( !formula.ok()      ) { return INVALID_FORMULA ; }
  //
  std::unique_ptr<Ostap::Formula> cuts { nullptr } ;
  if ( !trivial ( selection ) ) 
    { 
      cuts = std::make_unique<Ostap::Formula> ( selection , data ) ; 
      if ( !cuts || !cuts->ok() ) { return INVALID_SELECTION ; }
    }
  //
  const bool                with_cuts = cuts && cuts->ok () ;
  const Ostap::EventIndex   the_last  = std::min ( last , (Ostap::EventIndex) tree->GetEntries() ) ;
  if ( the_last <= first   ) { return Ostap::StatusCode::SUCCESS ; }
  //
  Ostap::Utils::Notifier    notify ( data , &formula  , cuts ) ;
  Ostap::Utils::ProgressBar bar    ( the_last - first , m_progress ) ;
  //
  std::vector<double>       results {} ;
  for ( Ostap::EventIndex   entry = first ; entry < the_last ; ++entry , ++bar )
    {
      //
      const long ievent = data -> GetEntryNumber ( entry ) ;
      if ( 0 > ievent ) { return INVALID_ENTRY ; }          // RETURN 
      //
      const ievent      = data -> LoadTree ( ievent ) ;
      if ( 0 > ievent ) { return INVALID_EVENT ; }          // RETURN 
      //
      const double w = with_cuts ? cuts->evaluate() : 1.0 ;
      //
      if ( !w ) { continue ; }                             // ATTENTION! 
      //
      formula.evaluate ( results ) ;
      for  ( const double r : results ) { moment.update ( r , w ) ; }
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
( const RooAbdData*         data       ,
  Ostap::Math::WStatistic&  stat       ,
  const std::string&        expression , 
  const std::string&        selection  ,
  const std::string&        cut_range  ,
  const Ostap::EventIndex   first      ,
  const Ostap::EventIndex   last       ) const
{
  // 
  if ( nullptr == data    ) { return INVALID_DATA    ; }
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  //
  const Ostap::EventIndex num_entries = data.numEntries() ;
  const Ostap::EventIndex the_last    = std::min ( num_entries , last ) ;
  if ( the_last <= first  ) { return Ostap::StatusCode::SUCCESS ; }
  //
  const bool  weighted = data -> isWeighted() ;
  const char* cutrange = cut_range.empty() ? nullptr : cut_range.c_str() ;
  const std::unique_ptr<Ostap::FormulaVar> expr { make_formula ( expr      , *data        ) } ;
  //
  const std::unique_ptr<Ostap::FormulaVar> cuts =
    trivial ( selection ) ? nullptr : make_formula ( selection , *data , true ) ;
  const bool with_cuts = nullptr != cuts ;
  //
  if ( !expr || !expr=>ok()            ) { return INVALID_FORMULA ; }
  //
  Ostap::Utils::ProgressBar bar { the_last - first , m_progres } ; 
  for ( unsigned long entry = first ; entry < the_last ; ++entry , ++bar )
    {
      const RooArgSet* vars = data -> get ( entry ) ;
      if ( nullptr == vars )                              { break    ; } // BREAK 
      //
      if ( cut_range && !vars->allInRange ( cut_range ) ) { continue ; } // CONTINUE    
      // apply cuts:
      const double wc = with_cuts ? cuts -> getVal () : 1.0 ;
      if ( !wc ) { continue ; }                                   // CONTINUE  
      // apply weight:
      const double wd = weighted  ? data -> weight() : 1.0 ;
      if ( !wd ) { continue ; }                                   // CONTINUE    
      // cuts & weight:
      const double w  = wd *  wc ;
      if ( !w  ) { continue ; }                                   // CONTINUE        
      //
      const double v = expr -> getVal () ;
      //
      stat.update ( value , w ) ;
    }
  //
  return Ostap::StatusCode::SUCCESS ; 
}
// =============================================================================




// =============================================================================
// Get information about several varibales 
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
  const std::string&           selecton    ,  
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
  const Ostap::EventIndex numentries =  data->GetEntries() ;
  const Ostap::EventIndex the_last = std::min ( first , num_entries ) ;
  if ( the_last <= first   ) { return Ostap::StatusCode::SUCCESS ; }
  //
  std::unique_ptr<Ostap::Formula> cuts { nullptr } ;
  if  ( !trivial ( selection ) ) 
    { 
      cuts = std::make_unique<Ostap::Formula> ( selection , data ) ; 
      if ( !cuts || !cuts->ok() ) { return INVALID_FORMULA ; }
    }
  //
  typedef std::unique_ptr<Ostap::Formula> UOF ;
  std::vector<UOF> formulas ; formulas.reserve ( N ) ;
  for ( const auto& e : expressions  ) 
    {
      auto p = std::make_unique<Ostap::Formula>( e , data ) ;
      if ( !p || !p->ok() ) { return 0 ; }
      formulas.push_back ( std::move ( p ) ) ;
    }
  //
  Ostap::Assert ( N == formulas.size()                 , 
                  "Inconsistent size of structures"    , 
                  "Ostap::StatVar::statVars"           ,
		  INVALID_FORMULAE , __FILE__ , __LINE ) ;
  //
  // 
  Ostap::Utils::Notifier notify ( formulas.begin() , formulas.end() , cuts , tree ) ;
  //
  std::vector<double>      results {} ;
  Ostap::Utils::ProressBar bar ( the_last - first , m_progress ) ;
  for ( Ostap::EventIndex  entry = first ; entry < the_last  ; ++entry , ++bar )
    {
      //
      const long ievent = tree->GetEntryNumber ( entry ) ;
      if ( 0 > ievent ) { return INVALID_ENTRY ; }                // RETURN
      //
      const long ievent = tree->LoadTree ( ievent ) ;
      const if ( 0 > ievent ) { return INVALID_EVENT ; }          // RETURN
      //
      const double w = cuts ? cuts ->evaluate() ; 1.0 ;
      if ( !w ) { continue  ; }
      //
      for ( std::size_t i = 0 ; i < N ; ++i ) 
	{
	  formulas [ i ]->evaluate ( results ) ;
	  for ( const double r : results ) { result[i].add ( r ) ; }
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
unsigned long Ostap::StatVar::statVars
( TTree*                       data        ,  
  Ostap::StatVar::WStatVector& result      ,  
  const Ostap::Strings&        expressions ,
  const std::string&           selecton    ,  
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
  const Ostap::EventIndex numentries =  data->GetEntries() ;
  const Ostap::EventIndex the_last = std::min ( first , num_entries ) ;
  if ( the_last <= first   ) { return Ostap::StatusCode::SUCCESS ; }
  //
  std::unique_ptr<Ostap::Formula> cuts { nullptr } ;
  if  ( !trivial ( selection ) ) 
    { 
      cuts = std::make_unique<Ostap::Formula> ( selection , data ) ; 
      if ( !cuts || !cuts->ok() ) { return INVALID_FORMULA ; }
    }
  //
  typedef std::unique_ptr<Ostap::Formula> UOF ;
  std::vector<UOF> formulas ; formulas.reserve ( N ) ;
  for ( const auto& e : expressions  ) 
    {
      auto p = std::make_unique<Ostap::Formula>( e , data ) ;
      if ( !p || !p->ok() ) { return 0 ; }
      formulas.push_back ( std::move ( p ) ) ;
    }
  //
  Ostap::Assert ( N == formulas.size()                 , 
                  "Inconsistent size of structures"    , 
                  "Ostap::StatVar::statVars"           ,
		  INVALID_FORMULAE , __FILE__ , __LINE ) ;
  //
  Ostap::Utils::Notifier   notify ( formulas.begin() , formulas.end() , cuts , tree ) ;
  //
  std::vector<double>      results {} ;
  Ostap::Utils::ProressBar bar ( the_last - first , m_progress ) ;
  for ( Ostap::EventIndex  entry = first ; entry < the_last  ; ++entry , ++bar )
    {
      //
      const long ievent = tree->GetEntryNumber ( entry ) ;
      if ( 0 > ievent ) { return INVALID_ENTRY ; }                // RETURN
      //
      const long ievent = tree->LoadTree ( ievent ) ;
      const if ( 0 > ievent ) { return INVALID_EVENT ; }          // RETURN
      //
      const double w = cuts ? cuts ->evaluate() ; 1.0 ;
      if ( !w ) { continue  ; }
      //
      for ( unsigned int i = 0 ; i < N ; ++i ) 
	{
	  formulas [ i ]->evaluate ( results ) ;
	  for ( const double r : results ) { result[i].add ( r , w ) ; }
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
  const unsigned long             first       ,
  const unsigned long             last        ) 
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
  const std::unique_ptr<Ostap::FormulaVar> cuts { make_formula ( selection , *data , true ) } ;
  //
  typedef std::unique_ptr<Ostap::FormulaVar> UOF ;
  std::vector<UOF> formulas ; formulas.reserve ( N ) ;
  //
  for ( const auto& e : expressions  ) 
    {
      auto p = make_formula ( e , *data , false ) ;
      if ( !p ) { return 0 ; }
      formulas.push_back ( std::move ( p ) ) ;  
    }
  //
  const char* cutrange = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  const bool  weighted = data->isWeighted() ;
  //
  // start the loop
  Ostap::Utils::ProgressBar bar ( the_last - first , m_progress ) ;
  for ( unsigned long entry = first ; entry < the_last ; ++entry , ++bar )
    {
      //
      const RooArgSet* vars = data->get( entry ) ;
      if ( nullptr == vars  )                           { break    ; } // RETURN
      if ( cutrange && !vars->allInRange ( cutrange ) ) { continue ; } // CONTINUE    
      //
      // apply cuts:
      const double wc = cuts ? cuts -> getVal() : 1.0 ;
      if ( !wc ) { continue ; }                                   // CONTINUE  
      // apply weight:
      const double wd = weighted  ? data->weight()        : 1.0 ;
      if ( !wd ) { continue ; }                                   // CONTINUE    
      // cuts & weight:
      const double w  = wd *  wc ;
      if ( !w  ) { continue ; }                                   // CONTINUE        
      //
      for ( unsigned int i = 0 ; i < N ; ++i ) 
	{
	  const double v = formulas[i]->getVal () ;
	  result[i].add ( v , w ) ;
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
  Ostap::StatVector      vct         { 1 , Ostap::StatEntity() } ;
  const Ostap::Strings   expressions { 1 , expression          } ; 
  const Ostap::SatusCode sc = statVars ( data , vct , expressions , "" , first , last ) ;
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
  Ostap::StatVector      vct         { 1 , Ostap::StatEntity() } ;
  const Ostap::Strings   expressions { 1 , expression          } ; 
  const Ostap::SatusCode sc = statVars ( tree , vct , expressions , selection , first , last ) ;
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
( TTree*                  tree       , 
  const std::string&      expression , 
  const std::string&      selection  , 
  const Ostap::EventIndex first      ,
  const Ostap::EventIndex last       ) const
{
  /// 
  Ostap::WStatVector     vct         { 1 , Ostap::WStatEntity() } ;
  const Ostap::Strings   expressions { 1 , expression           } ; 
  const Ostap::SatusCode sc = statVars ( tree , vct , expressions , selection , first , last ) ;
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
  const syd::string&      cut_range  ,  
  const Ostap::EventIndex first      ,
  const Ostap::EventIndex last       ) const
{
  /// 
  Ostap::WStatVector     vct         { 1 , Ostap::WStatEntity() } ;
  const Ostap::Strings   expressions { 1 , expression           } ; 
  const Ostap::SatusCode sc = statVars ( tree , vct , expression , selection , cut_range , first , last ) ;
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
( TTree*                  tree  ,
  const std::string&      exp1  ,
  const std::string&      exp2  ,
  const Ostap::EventIndex first ,
  const Ostap::EventIndex last  ) const 
{
  //
  Ostap::Assert ( nullptr != tree                    ,  
		  "Invalid TTree"                    , 
		  "Ostap::StatVar::statCov"          ,
		  INVALID_TREE , __FILE__ , __LINE__ ) ;
  //
  // prepare the result 
  Ostap::Math::Covariance result {} ;
  //
  const Ostap::EventInex num_entries = tree->GetEntries() ;
  const Ostap::EventIndex the_last   = std::min ( last , num_entries ) ;
  if ( the_last <= first ) { return result ; } 
  //
  Ostap::Formula formula1 ( exp1 , tree ) ;
  Ostap::Assert ( formula1.ok()                                         ,
		  std::string ( "Invalid first  expression: " ) + exp1  ,  
		  "Ostap::StatVar::statCov"                             ,
		  INVALID_FORMULA , __FILE__ , __LINE__ ) ;
  Ostap::Formula formula2 ( exp2 , tree ) ;
  Ostap::Assert ( formula2.ok()                                        ,
		  std::string ( "Invalid second expression: " ) + exp2 , 
		  "Ostap::StatVar::statCov"                            , 
		  INVALID_FORMULA , __FILE__ , __LINE__ ) ;
  //
  Ostap::Utils::Notifier notify ( tree , &formula1 , &formula2 ) ;
  //
  std::vector<double> results1 {} ;
  std::vector<double> results2 {} ;
  //
  Ostap::Utils::ProgressBar bar ( the_last - first , m_progress ) ;
  for ( Ostap::EventIndex entry = first ; entry < the_last ; ++entry , ++bar )
    {
      //
      long ievent = tree->GetEntryNumber ( entry ) ;
      if ( 0 > ievent ) { break ; }                        // BREAK
      //
      ievent      = tree->LoadTree ( ievent ) ;
      if ( 0 > ievent ) { break ; }                        // BREAK
      //
      formula1.evaluate ( results1 ) ;
      formula2.evaluate ( results2 ) ;
      //
      for ( const long double v1 : results1 ) 
	{ 
	  for ( const long double v2 : results2 ) 
	    { 
	      result.add ( v1 , v2 ) ;
	    }
	}
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
( TTree*                  tree      ,
  const std::string&      exp1      ,
  const std::string&      exp2      ,
  const std::string&      selection ,
  const Ostap::EventIndex first     ,
  const Ostap::EventIndex last      ) const 
{
  //
  Ostap::Assert ( nullptr != tree                    ,  
		  "Invalid TTree"                    , 
		  "Ostap::StatVar::statCov"          ,
		  INVALID_TREE , __FILE__ , __LINE__ ) ;
  //
  // prepare the result 
  Ostap::Math::Covariance result {} ;
  //
  const Ostap::EventInex num_entries = tree->GetEntries() ;
  const Ostap::EventIndex the_last   = std::min ( last , num_entries ) ;
  if ( the_last <= first ) { return result ; } 
  //
  Ostap::Formula formula1 ( exp1 , tree ) ; 
  Ostap::Assert ( formula1.ok()                                        ,
		  std::string ( "Invalid first expression: " ) + exp1  , 
		  "Ostap::StatVar::statCov"                            , 
		  INVALID_FORMULA , __FILE__ , __LINE__ ) ;
  //
  Ostap::Formula formula2 ( exp2 , tree ) ;
  Ostap::Assert ( formula2.ok()                                        ,
		  std::string ( "Invalid second expression: " ) + exp2 , 
		  "Ostap::StatVar::statCov"                            , 
		  INVALID_FORMULA , __FILE__ , __LINE__ ) ;
  //
  std::unique_ptr<Ostap::Formula> cuts { nullptr } ;
  if  ( !trivial ( selection ) ) 
    { 
      cuts = std::make_unique<Ostap::Formula> ( selection , data ) ; 
      Ostap::Assert ( cuts && cuts->ok ()                              ,
		      std::string ( "Invalid selction " ) + selection  ,  
		      "Ostap::StatVar::statCov"                        ,
		      INVALID_SELECTION , __FILE__ , __LINE__ ) ;
    }
  //
  Ostap::Utils::Notifier notify ( tree , &formula1 , &formula2 , cuts ) ;
  //
  std::vector<double> results1 {} ;
  std::vector<double> results2 {} ;
  //
  Ostap::Utils::ProgressBar bar ( the_last - first , m_progress ) ;
  for ( Ostap::EventIndex entry = first ; entry < the_last ; ++entry , ++bar )
    {
      //
      long ievent = tree->GetEntryNumber ( entry ) ;
      if ( 0 > ievent ) { break ; }                        // BREAK
      //
      ievent      = tree->LoadTree ( ievent ) ;
      if ( 0 > ievent ) { break ; }                        // BREAK
      //
      const double w = cuts ? cuts->evaluate() : 1.0 ;
      if ( !w ) { continue ; }
      //
      formula1.evaluate ( results1 ) ;
      formula2.evaluate ( results2 ) ;
      //
      for ( const long double v1 : results1 ) 
	{ 
	  for ( const long double v2 : results2 ) 
	    { 
	      result.add ( v1 , v2 ) ;
	    }
	}
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
( TTree*                  tree      ,
  const std::string&      exp1      ,
  const std::string&      exp2      ,
  const std::string&      selection ,
  const Ostap::EventIndex first     ,
  const Ostap::EventIndex last      ) const 
{
  //
  Ostap::Assert ( nullptr != tree                    ,  
		  "Invalid TTree"                    , 
		  "Ostap::StatVar::statCov"          ,
		  INVALID_TREE , __FILE__ , __LINE__ ) ;
  //
  // prepare the result 
  Ostap::Math::WCovariance result {} ;
  //
  const Ostap::EventInex num_entries = tree->GetEntries() ;
  const Ostap::EventIndex the_last   = std::min ( last , num_entries ) ;
  if ( the_last <= first ) { return result ; } 
  //
  Ostap::Formula formula1 ( exp1 , tree ) ; 
  Ostap::Assert ( formula1.ok()                                        ,
		  std::string ( "Invalid first expression: " ) + exp1 , 
		  "Ostap::StatVar::statCov"                            , 
		  INVALID_FORMULA , __FILE__ , __LINE__ ) ;
  Ostap::Formula formula2 ( exp2 , tree ) ;
  Ostap::Assert ( formula2.ok()                                        ,
		  std::string ( "Invalid second expression: " ) + exp2 , 
		  "Ostap::StatVar::statCov"                            , 
		  INVALID_FORMULA , __FILE__ , __LINE__ ) ;
  //
  std::unique_ptr<Ostap::Formula> cuts { nullptr } ;
  if  ( !trivial ( selection ) ) 
    { 
      cuts = std::make_unique<Ostap::Formula> ( selection , data ) ; 
      Ostap::Assert ( cuts && cuts->ok ()                              ,
		      std::string ( "Invalid selction " ) + selection  ,  
		      "Ostap::StatVar::statCov"                        ,
		      INVALID_SELECTION , __FILE__ , __LINE__ ) ;
    }
  //
  Ostap::Utils::Notifier notify ( tree , &formula1 , &formula2 , cuts ) ;
  //
  std::vector<double> results1 {} ;
  std::vector<double> results2 {} ;
  //
  Ostap::Utils::ProgressBar bar ( the_last - first , m_progress ) ;
  for ( Ostap::EventIndex entry = first ; entry < the_last ; ++entry , ++bar )
    {
      //
      long ievent = tree->GetEntryNumber ( entry ) ;
      if ( 0 > ievent ) { break ; }                        // BREAK
      //
      ievent      = tree->LoadTree ( ievent ) ;
      if ( 0 > ievent ) { break ; }                        // BREAK
      //
      const double w = cuts ? cuts->evaluate() : 1.0 ;
      if ( !w ) { continue ; }
      //
      formula1.evaluate ( results1 ) ;
      formula2.evaluate ( results2 ) ;
      //
      for ( const long double v1 : results1 ) 
	{ for ( const long double v2 : results2 ) 
	    { result.add ( v1 , v2 , w ) ; } }       // ADD WITH WEIGHT 
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
  const std::unique_ptr<Ostap::FormulaVar> formula1  { make_formula ( exp1 , *data        ) } ;
  const std::unique_ptr<Ostap::FormulaVar> formula2  { make_formula ( exp2 , *data        ) } ;
  const std::unique_ptr<Ostap::FormulaVar> selection { make_formula ( cuts , *data , true ) } ;
  //
  const char* cutrange = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  const bool  weighted = data->isWeighted() ;
  //
  //
  // start the loop
  Ostap::Utils::ProgressBar bar ( the_last - first , m_progress ) ;
  for ( unsigned long entry = first ; entry < the_last ; ++entry , ++bar )
    {
      //
      const RooArgSet* vars = data->get( entry ) ;
      if ( nullptr == vars  )                           { break    ; } // RETURN
      if ( cutrange && !vars->allInRange ( cutrange ) ) { continue ; } // CONTINUE    
      //
      // apply cuts:
      const double wc = selection ? selection -> getVal() : 1.0 ;
      if ( !wc ) { continue ; }                                   // CONTINUE  
      // apply weight:
      const double wd = weighted  ? data->weight()        : 1.0 ;
      if ( !wd ) { continue ; }                                   // CONTINUE    
      // cuts & weight:
      const double w  = wd *  wc ;
      if ( !w  ) { continue ; }                                   // CONTINUE        
      //
      const double v1 = formula1->getVal () ;
      const double v2 = formula2->getVal () ;
      //
      result.add ( v1 , v2  , w ) ;
    }
  //
  return result ;
}

















// ============================================================================
/*  calculate the covariance of several expressions 
 *  @param tree  (INPUT)  the inpout tree 
 *  @param vars  (INPUT)  expressions 
 *  @param cuts  (INPUT)  the selection criteria 
 *  @param stats (UPDATE) the statistics 
 *  @param cov2  (UPDATE) the covariance matrix 
 *  @return number of processed events 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2023-02-28
 */
// ============================================================================
unsigned long
Ostap::StatVar::statCov
( TTree*                          tree  , 
  const std::vector<std::string>& vars  , 
  const std::string&              cuts  ,
  Ostap::StatVar::WStatVector&    stats ,  
  TMatrixTSym<double>&            cov2  , 
  const unsigned long             first ,
  const unsigned long             last  ) 
{
  //
  cov2 *= 0.0 ;
  //
  if ( 0 == tree || last <= first ) { stats.clear() ; return 0 ; }//
  std::vector< std::unique_ptr<Ostap::Formula> > formulas ;
  std::vector< std::vector<double> >             results  ;
  std::vector< TObject*>                         objects  ;
  //
  formulas.reserve ( vars.size() ) ;
  results .reserve ( vars.size() ) ;
  //
  for ( std::vector<std::string>::const_iterator ie = vars.begin() ; vars.end() != ie ; ++ie ) 
  {
    std::unique_ptr<Ostap::Formula> expr { new Ostap::Formula ( *ie , tree ) } ;
    if ( !expr || !expr->ok() ) { stats.clear() ; return 0 ; }
    formulas.push_back ( std::move ( expr )    ) ;
    results.push_back  ( std::vector<double>() ) ;
    objects.push_back  ( formulas.back().get() ) ;
  }
  //
  const unsigned int N = formulas.size() ;
  if ( N < 1 ) { stats.clear() ; return 0 ; }
  //
  std::unique_ptr<Ostap::Formula> selection ;
  if ( !cuts.empty() ) 
  {
    selection.reset ( new Ostap::Formula ( cuts , tree ) ) ;
    if ( !selection || !selection->ok() ) { stats.clear() ; return 0 ; }
    objects.push_back ( selection.get() ) ;
  }
  //
  Ostap::Utils::Notifier notify ( objects.begin() , objects.end() , tree ) ;
  //
  const unsigned long nEntries =
    std::min ( last , (unsigned long) tree->GetEntries() ) ;
  //
  cov2 = TMatrixTSym<double>( N ) ; 
  stats.resize( N ) ;
  const bool with_cuts = selection && selection->ok () ;
  for ( WStatVector::iterator s = stats.begin() ; stats.end() != s ; ++s ) { s->reset() ; }
  //
  for ( unsigned long entry = first ; entry < nEntries ; ++entry )
  {
    //
    long ievent = tree->GetEntryNumber ( entry ) ;
    if ( 0 > ievent ) { break ; }                              // RETURN
    //
    ievent      = tree->LoadTree ( ievent ) ;
    if ( 0 > ievent ) { break ; }                              // RETURN
    //
    const double w = with_cuts ? selection->evaluate() : 1.0 ;
    //
    if ( !w ) { continue ; }                                   // ATTENTION
    //
    for ( unsigned int i = 0 ; i < N ; ++i ) 
    { formulas[i]->evaluate( results[i] ) ; }
    //
    for ( unsigned int i = 0 ; i < N ; ++i ) 
    {
      for ( const double ri : results[i] ) 
      {
        stats[i].add ( ri , w ) ;
        for ( unsigned int j = i  ; j < N ; ++j ) 
        {
          for ( const double rj : results[j] ) 
          { 
            const double val = w * ri * rj ;
            cov2 ( i , j ) += val ;
          }  
        }
      }
    }
  }
  //
  if ( 0 == stats[0].nEntries() ) { return 0 ; }
  //
  cov2 *= 1.0 / stats[0].weights().sum()  ;
  //
  for  (unsigned int i = 0 ; i < N ; ++i ) 
  {
    const double vi_mean = stats[i].mean() ;
    for  (unsigned int j = i ; j < N ; ++j ) 
      {
        const double vj_mean = stats[j].mean() ;
        const double val     = vi_mean * vj_mean ;
        cov2 ( i , j ) -= val ;
      }
  }
  //
  /// strange lines.... due to ROOT 
  for ( unsigned int i = 0 ; i < N ; ++i ) 
    { for ( unsigned int j = 0 ; j < i ; ++j ) 
        { if ( !cov2 ( i , j ) ) { cov2 ( i , j ) = cov2 ( j , i ) ; } } }
  //
  return stats[0].nEntries() ;
}
// ============================================================================


// ============================================================================
/*  calculate the covariance of several expressions 
 *  @param tree      (INPUT)  the inpout tree 
 *  @param vars      (INPUT)  expressions 
 *  @param cuts      (INPUT)  the selection criteria 
 *  @param stats     (UPDATE) the statistics 
 *  @param cov2      (UPDATE) the covariance matrix 
 *  @param cut_range (INPUT)  range  
 *  @return number of processed events 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2023-02-28
 */
// ============================================================================
unsigned long
Ostap::StatVar::statCov
( const RooAbsData*               data      , 
  const std::vector<std::string>& vars      ,  
  const std::string&              cuts      ,
  Ostap::StatVar::WStatVector&    stats     ,  
  TMatrixTSym<double>&            cov2      , 
  const std::string&              cut_range , 
  const unsigned long             first     ,
  const unsigned long             last      ) 
{
  //
  cov2 *= 0.0 ;
  //
  if ( 0 == data || last <= first ) { stats.clear() ; return 0 ; }
  //
  const bool  weighted = data->isWeighted() ;
  //
  std::vector<std::unique_ptr<Ostap::FormulaVar> > formulas ;
  //
  for ( std::vector<std::string>::const_iterator ie = vars.begin() ; vars.end() != ie ; ++ie ) 
  { formulas.push_back ( make_formula ( *ie , *data ) ) ; }
  //
  const std::unique_ptr<Ostap::FormulaVar> selection { make_formula ( cuts , *data ,  true ) } ;
  
  const unsigned int N = formulas.size() ;
  if ( 1 > N ) { stats.clear() ; return 0 ; }
  //
  cov2 = TMatrixTSym<double>( N ) ; 
  stats.resize( N ) ;
  //
  const char* cutrange = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  std::vector<double> results ( N , 0.0 ) ;
  //
  const bool with_cuts  = !(!selection) ;
  const unsigned long nEntries = std::min ( last , (unsigned long) data->numEntries() ) ;
  for ( unsigned long entry = first ; entry < nEntries ; ++entry )
  {
    //
    const RooArgSet* vars = data->get( entry ) ;
    if ( nullptr == vars  )                           { break    ; } // BREAK
    if ( cutrange && !vars->allInRange ( cutrange ) ) { continue ; } // CONTINUE    
    //
    // apply weight:
    const long double w = weighted  ? data->weight()        : 1.0L ;
    if ( !w      ) { continue ; }                                   // CONTINUE    
    //
    const double weight = with_cuts ? w * selection->getVal() : w ;
    if ( !weight ) { continue ; }                                   // CONTINUE    
    
    //
    for ( unsigned int i = 0 ; i < N ; ++i ) 
    { results [i] = formulas[i]->getVal() ; }
    
    for ( unsigned int i = 0 ; i < N ; ++i ) 
      {
        const double ri = results [i] ;
        stats[i].add ( ri , weight ) ;
        for ( unsigned int j = i ; j < N ; ++j ) 
          { cov2 ( i , j ) += weight * ri * results [ j ] ; }
      }
    //
  }
  //
  if ( 0 == stats[0].nEntries() ) { return 0 ; }
  //
  cov2 *= 1.0 / stats[0].weights().sum()  ;
  //
  for  (unsigned int i = 0 ; i < N ; ++i ) 
    {
      const double vi_mean = stats[i].mean() ;
      for  (unsigned int j = i ; j < N ; ++j ) 
        {
          const double vj_mean = stats[j].mean() ;
          cov2 ( i , j ) -= vi_mean * vj_mean ;
        }
    }
  //
  /// strange lines.... due to ROOT 
  for ( unsigned int i = 0 ; i < N ; ++i ) 
    { for ( unsigned int j = 0 ; j < i ; ++j ) 
        { if ( !cov2 ( i , j ) ) { cov2 ( i , j ) = cov2 ( j , i ) ; } } }
  //
  return stats[0].nEntries() ;
}
// ============================================================================

// ============================================================================
/*  Get the empirical cumulative distribution function 
 *  @param data  (INPUT) data 
 *  @param ecdf  (UDATE) cumulative distribtion function 
 *  @param expression (INPUT) the variable 
 *  @param first  (INPUT) the first event to process (inclusive)
 *  @param last   (INPUT) the last  event to process (non-inclusive)
 *  @returs status code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::StatVar::ECDF
( TTree*                  data       ,
  Ostap::Math::ECDF&      ecdf       ,
  const std::string&      expression ,
  const std::string&      selection  ,
  const Ostap::EventIndex first      ,
  const Ostap::EventIndex last       ) ;
{
  // reset ECDF 
  ecdf = Ostap::Math::ECDF () ;
  const Ostap::StatusCode sc = get_stat ( data , ecdf , expression , seelction , first , last ) ;
  if       ( sc.isFailure() ) { return sc ; }
  else if  ( !ecdf.ok()     ) { return INVALID_ECDF ; } 
  return Ostap::StatusCode::SUCCESS ; 
}
// ==========================================================================
/* Get the empirical cumulative distribtion function 
 *  @param data       (INPUT) data 
 *  @param ecdf       (UDATE) cumulative distribtion function 
 *  @param expression (INPUT) the variable 
 *  @param selection  (INOUT) selectgion/weight 
 *  @param first      (INPUT) the first event to process (inclusive)
 *  @param last       (INPUT) the last  event to process (non-inclusive)
 *  @returs status code 
 */
// ==========================================================================
Ostap::StatusCode
Ostap::StatVar::ECDF
( TTree*                  data       ,
  Ostap::Math::WECDF&     ecdf       ,
  const std::string&      expression , 
  const std::string&      selection  , 
  const Ostap::EventIndex first      ,
  const Ostap::EventIndex last       ) ;
{
  // reset ECDF 
  ecdf = Ostap::Math::WECDF () ;
  const Ostap::StatusCode sc = get_stat  ( data , ecdf , expression , selection , first , last ) ;
  if       ( sc.isFailure() ) { return sc ; }
  else if  ( !ecdf.ok()     ) { return INVALID_WECDF ; } 
  return Ostap::StatusCode::SUCCESS ; 
}
// ==========================================================================
/* Get the empirical cumulative distribtion function 
 *  @param data       (INPUT) data 
 *  @param ecdf       (UDATE) cumulative distribtion function 
 *  @param expression (INPUT) the variable 
 *  @param selection  (INOUT) selectgion/weight 
 *  @param first      (INPUT) the first event to process (inclusive)
 *  @param last       (INPUT) the last  event to process (non-inclusive)
 *  @returs status code 
 */
// ==========================================================================
Ostap::StatusCode
Ostap::StatVar::ECDF
( const RooAbsData*       data       ,
  Ostap::Math::WECDF&     ecdf       ,
  const std::string&      expression , 
  const std::string&      selection  , 
  const std::string&      cut_range  , 
  const Ostap::EventIndex first      ,
  const Ostap::EventIndex last       ) ;
{
  // reset ECDF 
  ecdf = Ostap::Math::WECDF () ;
  const Ostap::StatusCode sc = get_stat ( data , ecdf , expression , selection , cut_range , first , last ) ;
  if       ( sc.isFailure() ) { return sc ; }
  else if  ( !ecdf.ok()     ) { return INVALID_WECDF ; } 
  return Ostap::StatusCode::SUCCESS ; 
}



// ============================================================================
//                                                                      The END
// ============================================================================
