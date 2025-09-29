// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <cstring>
#include <memory>
// ============================================================================
// ROOT
// ============================================================================
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"
// ============================================================================
// ROOT/RooFit
// ============================================================================
#include "RooAbsData.h"
// ============================================================================
// Ostap
// ===========================================================================
#include "Ostap/StatVar.h"
#include "Ostap/GetWeight.h"
#include "Ostap/Project.h"
#include "Ostap/Exception.h"
#include "Ostap/ECDF.h"
#include "Ostap/Polynomials.h"
#include "Ostap/Bernstein.h"
#include "Ostap/Bernstein2D.h"
#include "Ostap/Bernstein3D.h"
#include "Ostap/Parameterization.h"
#include "Ostap/ProgressBar.h"
#include "Ostap/FormulaVar.h"
// ============================================================================
// Local 
// ============================================================================
#include "local_utils.h"
#include "status_codes.h"
#include "hstats.h"
// ============================================================================
/** @file
 *  Implementation file for class Ostap::HistoProject
 *  @see Ostap::HistoProject
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2015-10-08
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
// construtctor with the progress flag
// ============================================================================
Ostap::Project::Project					
( const Ostap::Utils::ProgressConf& progress )
  : Ostap::StatVar ( progress )
{}
// ============================================================================

// ============================================================================
// 1D-histos 
// ============================================================================
/*  Project data in the 1D-ihstogram
 *  @param data       (INPUT)  input data 
 *  @param histo      (UPDATE) the !D-histogram 
 *  @param expression (INPUT)  expression 
 *  @param selection  (INPUT)  selection 
 *  @param first      (INOPUT) the first event to process (inclusive) 
 *  @param last       (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 *  @see   TH1::Fill
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project1
( TTree*                  data        ,
  TH1*                    histo       ,
  const std::string&      expression  ,
  const std::string&      selection   ,
  const Ostap::EventIndex first       ,
  const Ostap::EventIndex last        ) const 
{
  if ( !histo || 1 != histo->GetDimension () ) { return INVALID_TH1 ; }
  //
  Ostap::Utils::H1 h1 ( histo )    ;
  //
  const TAxis* xaxis = histo->GetXaxis () ;
  if ( !xaxis ) { return INVALID_XAXIS ; }
  //
  return get_stat ( data                ,
		    h1                  ,
		    expression          ,
		    selection           ,
		    first               ,
		    last                ,
		    xaxis -> GetXmin () ,
		    xaxis -> GetXmax () ) ;
}
// ============================================================================
/*  Project data in the 1D-ihstogram
 *  @param data       (INPUT)  input data 
 *  @param histo      (UPDATE) the 1D-histogram 
 *  @param expression (INPUT)  expression 
 *  @param selection  (INPUT)  selection 
 *  @param first      (INOPUT) the first event to process (inclusive) 
 *  @param last       (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 *  @see   TH1::Fill
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project1
( const RooAbsData*       data        ,
  TH1*                    histo       ,
  const std::string&      expression  ,
  const std::string&      selection   ,
  const std::string&      cut_range   ,  
  const Ostap::EventIndex first       ,
  const Ostap::EventIndex last        ) const 
{
  if ( !data )                                 { return INVALID_DATA ; }
  if ( !histo || 1 != histo->GetDimension () ) { return INVALID_TH1  ; }
  //
  histo->Reset () ;
  if ( !histo->GetSumw2() ) { histo->Sumw2() ; }
  //
  // check the ranges
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  //
  const Ostap::EventIndex num_entries = data -> numEntries() ;
  const Ostap::EventIndex the_last    = std::min ( num_entries , last ) ;
  if ( the_last <= first  ) { return Ostap::StatusCode::SUCCESS ; }  
  //
  const TAxis* xaxis = histo->GetXaxis () ;
  if ( !xaxis ) { return INVALID_XAXIS ; }
  //
  const double xmin = xaxis->GetXmin () ;
  const double xmax = xaxis->GetXmax () ;
  //
  // special processing where the weighth has errors
  //
  if ( data->isWeighted() && Ostap::Utils::storeError ( data ) )
    {
      //
      const char* cutrange = cut_range.empty() ? nullptr : cut_range.c_str() ;
      const std::unique_ptr<Ostap::FormulaVar> expr { Ostap::makeFormula ( expression , data       ) } ;
      if ( !expr || !expr->ok()            ) { return INVALID_FORMULA ; }
      //
      const std::unique_ptr<Ostap::FormulaVar> cuts { Ostap::makeFormula ( selection  , data , true ) } ;
      const bool                with_cuts   = cuts && cuts->ok () ;
      //
      Ostap::Utils::ProgressBar bar { the_last - first , progress () } ; 
      for ( Ostap::EventIndex entry = first ; entry < the_last ; ++entry , ++bar )	
	{
	  const RooArgSet* vars = data -> get ( entry ) ;
	  if ( nullptr == vars )                             { break    ; } // BREAK 
	  //
	  if ( cutrange && !vars->allInRange ( cutrange ) )  { continue ; } // CONTINUE
	  //
	  // data weight?
	  const double wd  = data-> weight() ;
	  if ( !wd ) { continue ; }                                       // CONTINUE
	  //
	  // apply cuts:
	  const double wc = with_cuts ? cuts -> getVal () : 1.0 ;
	  if ( !wc ) { continue ; }                                        // CONTINUE  
	  //
	  // total weight
	  const double wt = wd * wc ;
	  if ( !wt ) { continue ; } 
	  //
	  const double value = expr -> getVal () ;
	  if ( !in_range ( value , xmin , xmax ) ) { continue ; }          // CONTINUE
	  //
	  histo->Fill ( value , wt ) ;
	  //
	  // correct the errors 
	  const double we  = data -> weightError() * wc ;
	  if ( we )
	    {
	      const int    bin = histo->FindBin     ( value ) ;
	      const double he  = histo->GetBinError ( bin   ) ;
	      const double e2  = he * he - wt * wt + we * we ;
	      histo->SetBinError ( bin , std::sqrt ( std::abs ( e2 ) ) ) ;
	    }
	}
      return Ostap::StatusCode::SUCCESS ;
    }
  //
  // regular processing 
  //  
  Ostap::Utils::H1 h1 ( histo )    ;
  //
  return get_stat ( data       ,
		    h1         ,
		    expression ,
		    selection  ,
		    cut_range  , 
		    first      ,
		    last       ,
		    xmin       ,
		    xmax       ) ;
}
// ============================================================================

// ============================================================================
// 2D-histos 
// ============================================================================
/*  Project data in the 2D-ihstoram
 *  @param data        (INPUT)  input data 
 *  @param histo       (UPDATE) the !D-histogram 
 *  @param expression1 (INPUT)  expression1 
 *  @param expression2 (INPUT)  expression2 
 *  @param selection   (INPUT)  selection/weight  
 *  @param first       (INOPUT) the first event to process (inclusive) 
 *  @param last        (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 *  @see   TH1::Fill
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project2
( TTree*                  data        ,
  TH2*                    histo       ,
  const std::string&      expression1 ,
  const std::string&      expression2 ,
  const std::string&      selection   ,
  const Ostap::EventIndex first       ,
  const Ostap::EventIndex last        ) const
{
  if ( !histo || 2 != histo->GetDimension () ) { return INVALID_TH2 ; }
  //
  Ostap::Utils::H2 h2 ( histo )    ;
  //
  const TAxis* xaxis = histo->GetXaxis () ;
  if ( !xaxis ) { return INVALID_XAXIS ; }
  //
  const TAxis* yaxis = histo->GetYaxis () ;
  if ( !yaxis ) { return INVALID_YAXIS ; }
  //
  return get_stat ( data                ,
		    h2                  ,
		    expression1         ,
		    expression2         ,
		    selection           ,
		    first               ,
		    last                ,
		    xaxis -> GetXmin () ,
		    xaxis -> GetXmax () ,
		    yaxis -> GetXmin () ,
		    yaxis -> GetXmax () ) ;
}
// ============================================================================
/*  Project data in the 2D-ihstoram
 *  @param data        (INPUT)  input data 
 *  @param histo       (UPDATE) the !D-histogram 
 *  @param expression1 (INPUT)  expression1 
 *  @param expression2 (INPUT)  expression2 
 *  @param selection   (INPUT)  selection/weight  
 *  @param first       (INOPUT) the first event to process (inclusive) 
 *  @param last        (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 *  @see   TH1::Fill
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project2
( const RooAbsData*       data        ,
  TH2*                    histo       ,
  const std::string&      expression1 ,
  const std::string&      expression2 ,
  const std::string&      selection   ,
  const std::string&      cut_range   ,
  const Ostap::EventIndex first       ,
  const Ostap::EventIndex last        ) const
{
  if ( !data )                                 { return INVALID_DATA ; }
  if ( !histo || 2 != histo->GetDimension () ) { return INVALID_TH2 ; }
  //
  histo->Reset () ;
  if ( !histo->GetSumw2() ) { histo->Sumw2() ; }
  //
  // check the ranges
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }  
  //
  const Ostap::EventIndex num_entries = data -> numEntries() ;
  const Ostap::EventIndex the_last    = std::min ( num_entries , last ) ;
  if ( the_last <= first  ) { return Ostap::StatusCode::SUCCESS ; }  
  //
  const TAxis* xaxis = histo->GetXaxis () ;
  if ( !xaxis ) { return INVALID_XAXIS ; }
  //
  const double xmin = xaxis->GetXmin () ;
  const double xmax = xaxis->GetXmax () ;
  //
  const TAxis* yaxis = histo->GetYaxis () ;
  if ( !yaxis ) { return INVALID_YAXIS ; }
  //
  const double ymin = yaxis->GetXmin () ;
  const double ymax = yaxis->GetXmax () ;
  //   
  // special processing where the weighth has errors
  //
  if ( data->isWeighted() && Ostap::Utils::storeError ( data ) )
    {
      //
      const char* cutrange = cut_range.empty() ? nullptr : cut_range.c_str() ;
      //
      const std::unique_ptr<Ostap::FormulaVar> xexpr { Ostap::makeFormula ( expression1 , data       ) } ;
      if ( !xexpr || !xexpr->ok()            ) { return INVALID_FORMULA ; }
      //
      const std::unique_ptr<Ostap::FormulaVar> yexpr { Ostap::makeFormula ( expression2 , data       ) } ;
      if ( !yexpr || !yexpr->ok()            ) { return INVALID_FORMULA ; }
      //
      const std::unique_ptr<Ostap::FormulaVar> cuts { Ostap::makeFormula ( selection  , data , true ) } ;
      const bool                with_cuts   = cuts && cuts->ok () ;
      //
      Ostap::Utils::ProgressBar bar { the_last - first , progress () } ; 
      for ( Ostap::EventIndex entry = first ; entry < the_last ; ++entry , ++bar )	
	{
	  const RooArgSet* vars = data -> get ( entry ) ;
	  if ( nullptr == vars )                             { break    ; } // BREAK 
	  //
	  if ( cutrange && !vars->allInRange ( cutrange ) )  { continue ; } // CONTINUE
	  //
	  // data weight?
	  const double wd  = data-> weight() ;
	  if ( !wd ) { continue ; }                                       // CONTINUE
	  //
	  // apply cuts:
	  const double wc = with_cuts ? cuts -> getVal () : 1.0 ;
	  if ( !wc ) { continue ; }                                        // CONTINUE  
	  //
	  // total weight
	  const double wt = wd * wc ;
	  if ( !wt ) { continue ; } 
	  //
	  const double xvalue = xexpr -> getVal () ;
	  if ( !in_range ( xvalue , xmin , xmax ) ) { continue ; }          // CONTINUE
	  //
	  const double yvalue = yexpr -> getVal () ;
	  if ( !in_range ( yvalue , ymin , ymax ) ) { continue ; }          // CONTINUE
	  //
	  histo->Fill ( xvalue , yvalue , wt ) ;
	  //
	  // correct the errors 
	  const double we  = data -> weightError() * wc ;
	  if ( we )
	    {
	      const int    xbin = xaxis->FindBin    ( xvalue ) ;
	      const int    ybin = yaxis->FindBin    ( yvalue ) ;
	      const double he  = histo->GetBinError ( xbin , ybin  ) ;
	      const double e2  = he * he - wt * wt + we * we ;
	      histo->SetBinError ( xbin , ybin , std::sqrt ( std::abs ( e2 ) ) ) ;
	    }
	}
      return Ostap::StatusCode::SUCCESS ;
    }
  // 
  // regular processing 
  // 
  Ostap::Utils::H2 h2 ( histo ) ;
  //
  return get_stat ( data        ,
		    h2          ,
		    expression1 ,
		    expression2 ,
		    selection   ,
		    cut_range   , 
		    first       ,
		    last        ,
		    xmin        ,
		    xmax        ,
		    ymin        ,
		    ymax        ) ; 
}
// ============================================================================

// ============================================================================
// 1D-profile 
// ============================================================================
/*  Project data in the 1D-profile 
 *  @param data        (INPUT)  input data 
 *  @param histo       (UPDATE) the 1D-profile 
 *  @param expression1 (INPUT)  expression1 
 *  @param expression2 (INPUT)  expression2 
 *  @param selection   (INPUT)  selection/weight  
 *  @param first       (INOPUT) the first event to process (inclusive) 
 *  @param last        (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 *  @see   TH1::Fill
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project2
( TTree*                  data        ,
  TProfile*               histo       ,
  const std::string&      expression1 ,
  const std::string&      expression2 ,
  const std::string&      selection   ,
  const Ostap::EventIndex first       ,
  const Ostap::EventIndex last        ) const
{
  if ( !histo || 1 != histo->GetDimension () ) { return INVALID_TPROFILE ; }
  //
  Ostap::Utils::P1 p1 ( histo )    ;
  //
  const TAxis* xaxis = histo->GetXaxis () ;
  if ( !xaxis ) { return INVALID_XAXIS ; }
  //
  return get_stat ( data                ,
		    p1                  ,
		    expression1         ,
		    expression2         ,
		    selection           ,
		    first               ,
		    last                ,
		    xaxis -> GetXmin () ,
		    xaxis -> GetXmax () ) ;
}
// ============================================================================
/*  Project data in the 1D-profile 
 *  @param data        (INPUT)  input data 
 *  @param histo       (UPDATE) the 1D-profile 
 *  @param expression1 (INPUT)  expression1 
 *  @param expression2 (INPUT)  expression2 
 *  @param selection   (INPUT)  selection/weight  
 *  @param first       (INOPUT) the first event to process (inclusive) 
 *  @param last        (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 *  @see   TH1::Fill
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project2
( const RooAbsData*       data        ,
  TProfile*               histo       ,
  const std::string&      expression1 ,
  const std::string&      expression2 ,
  const std::string&      selection   ,
  const std::string&      cut_range   ,
  const Ostap::EventIndex first       ,
  const Ostap::EventIndex last        ) const
{
  if ( !histo || 1 != histo->GetDimension () ) { return INVALID_TPROFILE ; }
  //
  Ostap::Utils::P1 p1 ( histo )    ;
  //
  const TAxis* xaxis = histo->GetXaxis () ;
  if ( !xaxis ) { return INVALID_XAXIS ; }
  //
  return get_stat ( data                ,
		    p1                  ,
		    expression1         ,
		    expression2         ,
		    selection           ,
		    cut_range           , 
		    first               ,
		    last                ,
		    xaxis -> GetXmin () ,
		    xaxis -> GetXmax () ) ;
}
// ============================================================================

// ============================================================================
// 3D-histos 
// ============================================================================
/* Project data in the 3D-ihstoram
 *  @param data        (INPUT)  input data 
 *  @param histo       (UPDATE) the !D-histogram 
 *  @param expression1 (INPUT)  expression1 
 *  @param expression2 (INPUT)  expression2 
 *  @param expression2 (INPUT)  expression3 
 *  @param selection   (INPUT)  selection/weight  
 *  @param first       (INOPUT) the first event to process (inclusive) 
 *  @param last        (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 *  @see   TH1::Fill
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project3
( TTree*                  data        ,
  TH3*                    histo       ,
  const std::string&      expression1 ,
  const std::string&      expression2 ,
  const std::string&      expression3 ,
  const std::string&      selection   ,
  const Ostap::EventIndex first       ,
  const Ostap::EventIndex last        ) const
{
  if ( !histo || 3 != histo->GetDimension () ) { return INVALID_TH3 ; }
  //
  Ostap::Utils::H3 h3 ( histo )    ;
  //
  const TAxis* xaxis = histo->GetXaxis () ;
  if ( !xaxis ) { return INVALID_XAXIS ; }
  //
  const TAxis* yaxis = histo->GetYaxis () ;
  if ( !yaxis ) { return INVALID_YAXIS ; }
  //
  const TAxis* zaxis = histo->GetZaxis () ;
  if ( !zaxis ) { return INVALID_ZAXIS ; }
  //
  return get_stat ( data                ,
		    h3                  ,
		    expression1         ,
		    expression2         ,
		    expression3         ,
		    selection           ,
		    first               ,
		    last                ,
		    xaxis -> GetXmin () ,
		    xaxis -> GetXmax () ,
		    yaxis -> GetXmin () ,
		    yaxis -> GetXmax () ,
		    zaxis -> GetXmin () ,
		    zaxis -> GetXmax () ) ;
} 
// ============================================================================
/* Project data in the 2D-ihstoram
 *  @param data        (INPUT)  input data 
 *  @param histo       (UPDATE) the !D-histogram 
 *  @param expression1 (INPUT)  expression1 
 *  @param expression2 (INPUT)  expression2 
 *  @param expression2 (INPUT)  expression3 
 *  @param selection   (INPUT)  selection/weight  
 *  @param first       (INOPUT) the first event to process (inclusive) 
 *  @param last        (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 *  @see   TH1::Fill
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project3
( const RooAbsData*       data        ,
  TH3*                    histo       ,
  const std::string&      expression1 ,
  const std::string&      expression2 ,
  const std::string&      expression3 ,
  const std::string&      selection   ,
  const std::string&      cut_range   , 
  const Ostap::EventIndex first       ,
  const Ostap::EventIndex last        ) const
{
  
  if ( !data )                                 { return INVALID_DATA ; }
  if ( !histo || 2 != histo->GetDimension () ) { return INVALID_TH2 ; }
  //
  histo->Reset () ;
  if ( !histo->GetSumw2() ) { histo->Sumw2() ; }
  //
  // check the ranges
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }  
  //
  const Ostap::EventIndex num_entries = data -> numEntries() ;
  const Ostap::EventIndex the_last    = std::min ( num_entries , last ) ;
  if ( the_last <= first  ) { return Ostap::StatusCode::SUCCESS ; }  
  //
  const TAxis* xaxis = histo->GetXaxis () ;
  if ( !xaxis ) { return INVALID_XAXIS ; }
  //
  const double xmin = xaxis->GetXmin () ;
  const double xmax = xaxis->GetXmax () ;
  //
  const TAxis* yaxis = histo->GetYaxis () ;
  if ( !yaxis ) { return INVALID_YAXIS ; }
  //
  const double ymin = yaxis->GetXmin () ;
  const double ymax = yaxis->GetXmax () ;
  //   
  const TAxis* zaxis = histo->GetZaxis () ;
  if ( !zaxis ) { return INVALID_ZAXIS ; }
  //
  const double zmin = zaxis->GetXmin () ;
  const double zmax = zaxis->GetXmax () ;
  //   
  // special processing where the weighth has errors
  //
  if ( data->isWeighted() && Ostap::Utils::storeError ( data ) )
    {
      //
      const char* cutrange = cut_range.empty() ? nullptr : cut_range.c_str() ;
      //
      const std::unique_ptr<Ostap::FormulaVar> xexpr { Ostap::makeFormula ( expression1 , data       ) } ;
      if ( !xexpr || !xexpr->ok()            ) { return INVALID_FORMULA ; }
      //
      const std::unique_ptr<Ostap::FormulaVar> yexpr { Ostap::makeFormula ( expression2 , data       ) } ;
      if ( !yexpr || !yexpr->ok()            ) { return INVALID_FORMULA ; }
      //
      const std::unique_ptr<Ostap::FormulaVar> zexpr { Ostap::makeFormula ( expression3 , data       ) } ;
      if ( !zexpr || !zexpr->ok()            ) { return INVALID_FORMULA ; }
      //
      const std::unique_ptr<Ostap::FormulaVar> cuts { Ostap::makeFormula ( selection  , data , true ) } ;
      const bool                with_cuts   = cuts && cuts->ok () ;
      //
      Ostap::Utils::ProgressBar bar { the_last - first , progress () } ; 
      for ( Ostap::EventIndex entry = first ; entry < the_last ; ++entry , ++bar )	
	{
	  const RooArgSet* vars = data -> get ( entry ) ;
	  if ( nullptr == vars )                             { break    ; } // BREAK 
	  //
	  if ( cutrange && !vars->allInRange ( cutrange ) )  { continue ; } // CONTINUE
	  //
	  // data weight?
	  const double wd  = data-> weight() ;
	  if ( !wd ) { continue ; }                                       // CONTINUE
	  //
	  // apply cuts:
	  const double wc = with_cuts ? cuts -> getVal () : 1.0 ;
	  if ( !wc ) { continue ; }                                        // CONTINUE  
	  //
	  // total weight
	  const double wt = wd * wc ;
	  if ( !wt ) { continue ; } 
	  //
	  const double xvalue = xexpr -> getVal () ;
	  if ( !in_range ( xvalue , xmin , xmax ) ) { continue ; }          // CONTINUE
	  //
	  const double yvalue = yexpr -> getVal () ;
	  if ( !in_range ( yvalue , ymin , ymax ) ) { continue ; }          // CONTINUE
	  //
	  const double zvalue = zexpr -> getVal () ;
	  if ( !in_range ( zvalue , zmin , zmax ) ) { continue ; }          // CONTINUE
	  //
	  histo->Fill ( xvalue , yvalue , zvalue , wt ) ;
	  //
	  // correct the errors 
	  const double we  = data -> weightError() * wc ;
	  if ( we )
	    {
	      const int    xbin = xaxis->FindBin    ( xvalue ) ;
	      const int    ybin = yaxis->FindBin    ( yvalue ) ;
	      const int    zbin = yaxis->FindBin    ( zvalue ) ;
	      const double he  = histo->GetBinError ( xbin , ybin , zbin ) ;
	      const double e2  = he * he - wt * wt + we * we ;
	      histo->SetBinError ( xbin , ybin , zbin , std::sqrt ( std::abs ( e2 ) ) ) ;
	    }
	}
      return Ostap::StatusCode::SUCCESS ;
    }
  // 
  Ostap::Utils::H3 h3 ( histo ) ;
  //
  return get_stat ( data        ,
		    h3          ,
		    expression1 ,
		    expression2 ,
		    expression3 ,
		    selection   ,
		    cut_range   , 
		    first       ,		    
		    last        ,
		    xmin        ,
		    xmax        , 
		    ymin        ,
		    ymax        ,
		    zmin        ,
		    zmax        ) ;
}
// ============================================================================

// ============================================================================
// 2D-profiles 
// ============================================================================
/* Project data in the 3D-ihstoram
 *  @param data        (INPUT)  input data 
 *  @param histo       (UPDATE) the !D-histogram 
 *  @param expression1 (INPUT)  expression1 
 *  @param expression2 (INPUT)  expression2 
 *  @param expression2 (INPUT)  expression3 
 *  @param selection   (INPUT)  selection/weight  
 *  @param first       (INOPUT) the first event to process (inclusive) 
 *  @param last        (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 *  @see   TH1::Fill
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project3
( TTree*                  data        ,
  TProfile2D*             histo       ,
  const std::string&      expression1 ,
  const std::string&      expression2 ,
  const std::string&      expression3 ,
  const std::string&      selection   ,
  const Ostap::EventIndex first       ,
  const Ostap::EventIndex last        ) const
{
  if ( !histo || 2 != histo->GetDimension () ) { return INVALID_TPROFILE2D ; }
  //
  Ostap::Utils::P2 p2 ( histo )    ;
  //
  const TAxis* xaxis = histo->GetXaxis () ;
  if ( !xaxis ) { return INVALID_XAXIS ; }
  //
  const TAxis* yaxis = histo->GetYaxis () ;
  if ( !yaxis ) { return INVALID_YAXIS ; }
  //
  return get_stat ( data                ,
		    p2                  ,
		    expression1         ,
		    expression2         ,
		    expression3         ,
		    selection           ,
		    first               ,
		    last                ,
		    xaxis -> GetXmin () ,
		    xaxis -> GetXmax () ,
		    yaxis -> GetXmin () ,
		    yaxis -> GetXmax () ) ;
} 
// ============================================================================
/* Project data in the 2D-ihstoram
 *  @param data        (INPUT)  input data 
 *  @param histo       (UPDATE) the !D-histogram 
 *  @param expression1 (INPUT)  expression1 
 *  @param expression2 (INPUT)  expression2 
 *  @param expression2 (INPUT)  expression3 
 *  @param selection   (INPUT)  selection/weight  
 *  @param first       (INOPUT) the first event to process (inclusive) 
 *  @param last        (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 *  @see   TH1::Fill
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project3
( const RooAbsData*       data        ,
  TProfile2D*             histo       ,
  const std::string&      expression1 ,
  const std::string&      expression2 ,
  const std::string&      expression3 ,
  const std::string&      selection   ,
  const std::string&      cut_range   , 
  const Ostap::EventIndex first       ,
  const Ostap::EventIndex last        ) const
{
  if ( !histo || 2 != histo->GetDimension () ) { return INVALID_TPROFILE2D ; }
  //
  Ostap::Utils::P2 p2 ( histo )    ;
  //
  const TAxis* xaxis = histo->GetXaxis () ;
  if ( !xaxis ) { return INVALID_XAXIS ; }
  //
  const TAxis* yaxis = histo->GetYaxis () ;
  if ( !yaxis ) { return INVALID_YAXIS ; }
  //
  return get_stat ( data                ,
		    p2                  ,
		    expression1         ,
		    expression2         ,
		    expression3         ,
		    selection           ,
		    cut_range           , 
		    first               ,
		    last                ,
		    xaxis -> GetXmin () ,
		    xaxis -> GetXmax () ,
		    yaxis -> GetXmin () ,
		    yaxis -> GetXmax () ) ;
}
// ============================================================================

// ============================================================================
// ECDF & WECDF 
// ============================================================================
/*  Project data in ECDS/WECDF 
 *  @param data       (INPUT)  input data 
 *  @param poly       (UPDATE) the polynoimal 
 *  @param expression (INPUT)  expression 
 *  @param selection  (INPUT)  selection (as boolean)
 *  @param first      (INOPUT) the first event to process (inclusive) 
 *  @param last       (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project1
( TTree*                    data       ,
  Ostap::Math::ECDF&        ecdf       ,
  const std::string&        expression ,
  const std::string&        selection  ,
  const Ostap::EventIndex   first      ,
  const Ostap::EventIndex   last       , 
  const Ostap::DataType     xmin       , 
  const Ostap::DataType     xmax       ) const
{
  return get_stat ( data                ,
		    ecdf                ,
		    expression          ,
		    selection           ,
		    first               ,
		    last                ,
		    xmin                ,
		    xmax                ) ;
}
// ============================================================================
/*  Project data in ECDS/WECDF 
 *  @param data       (INPUT)  input data 
 *  @param poly       (UPDATE) the polynoimal 
 *  @param expression (INPUT)  expression 
 *  @param selection  (INPUT)  selection (as weighjt )
 *  @param first      (INOPUT) the first event to process (inclusive) 
 *  @param last       (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project1
( TTree*                    data       ,
  Ostap::Math::WECDF&       ecdf       ,
  const std::string&        expression ,
  const std::string&        selection  ,
  const Ostap::EventIndex   first      ,
  const Ostap::EventIndex   last       , 
  const Ostap::DataType     xmin       , 
  const Ostap::DataType     xmax       ) const
{
  return get_stat ( data                ,
		    ecdf                ,
		    expression          ,
		    selection           ,
		    first               ,
		    last                ,
		    xmin                ,
		    xmax                ) ;
}
// ============================================================================
/*  Project data in ECDS/WECDF 
 *  @param data       (INPUT)  input data 
 *  @param poly       (UPDATE) the polynoimal 
 *  @param expression (INPUT)  expression 
 *  @param selection  (INPUT)  selection (as weight) 
 *  @param first      (INOPUT) the first event to process (inclusive) 
 *  @param last       (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project1
( const RooAbsData*         data       ,
  Ostap::Math::WECDF&       ecdf       ,
  const std::string&        expression ,
  const std::string&        selection  ,
  const std::string&        cut_range  , 
  const Ostap::EventIndex   first      , 
  const Ostap::EventIndex   last       , 
  const Ostap::DataType     xmin       ,
  const Ostap::DataType     xmax       ) const 
{
  return get_stat ( data                ,
		    ecdf                ,
		    expression          ,
		    selection           ,
		    cut_range           , 
		    first               ,
		    last                ,
		    xmin                ,
		    xmax                ) ;
}
// ============================================================================

// ============================================================================
// 1D-chebyshev polynoials: on-flight parameterization
// ============================================================================

// ============================================================================
/*  Project data in 1D-polyom=ominal: on-flight parameterisation 
 *  @param data       (INPUT)  input data 
 *  @param poly       (UPDATE) the polynoimal 
 *  @param expression (INPUT)  expression 
 *  @param selection  (INPUT)  selection (as boolean)
 *  @param first      (INOPUT) the first event to process (inclusive) 
 *  @param last       (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project1
( TTree*                     data       ,
  Ostap::Math::ChebyshevSum& poly       ,
  const std::string&         expression ,
  const std::string&         selection  ,
  const Ostap::EventIndex    first      ,
  const Ostap::EventIndex    last       ) const
{
  return get_stat ( data                ,
		    poly                ,
		    expression          ,
		    selection           ,
		    first               ,
		    last                ,
		    poly.xmin ()        ,
		    poly.xmax ()        ) ;
}
// =========================================================================
/*  Project data in 1D-polyom=ominal: on-flight parameterisation 
 *  @param data       (INPUT)  input data 
 *  @param poly       (UPDATE) the polynoimal 
 *  @param expression (INPUT)  expression 
 *  @param selection  (INPUT)  selection (as boolean)
 *  @param first      (INOPUT) the first event to process (inclusive) 
 *  @param last       (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project1
( const RooAbsData*          data       ,
  Ostap::Math::ChebyshevSum& poly       ,
  const std::string&         expression ,
  const std::string&         selection  ,
  const std::string&         cut_range  , 
  const Ostap::EventIndex    first      ,
  const Ostap::EventIndex    last       ) const
{
  return get_stat ( data                ,
		    poly                ,
		    expression          ,
		    selection           ,
		    cut_range           , 
		    first               ,
		    last                ,
		    poly.xmin ()        ,
		    poly.xmax ()        ) ;
}
// =========================================================================
  

// ============================================================================
// 1D-Legendre polynoials: on-flight parameterization
// ============================================================================

// ============================================================================
/*  Project data in 1D-polyom=ominal: on-flight parameterisation 
 *  @param data       (INPUT)  input data 
 *  @param poly       (UPDATE) the polynoimal 
 *  @param expression (INPUT)  expression 
 *  @param selection  (INPUT)  selection (as boolean)
 *  @param first      (INOPUT) the first event to process (inclusive) 
 *  @param last       (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project1
( TTree*                     data       ,
  Ostap::Math::LegendreSum&  poly       ,
  const std::string&         expression ,
  const std::string&         selection  ,
  const Ostap::EventIndex    first      ,
  const Ostap::EventIndex    last       ) const
{
  return get_stat ( data                ,
		    poly                ,
		    expression          ,
		    selection           ,
		    first               ,
		    last                ,
		    poly.xmin ()        ,
		    poly.xmax ()        ) ;
}
// =========================================================================
/*  Project data in 1D-polyom=ominal: on-flight parameterisation 
 *  @param data       (INPUT)  input data 
 *  @param poly       (UPDATE) the polynoimal 
 *  @param expression (INPUT)  expression 
 *  @param selection  (INPUT)  selection (as boolean)
 *  @param first      (INOPUT) the first event to process (inclusive) 
 *  @param last       (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project1
( const RooAbsData*          data       ,
  Ostap::Math::LegendreSum&  poly       ,
  const std::string&         expression ,
  const std::string&         selection  ,
  const std::string&         cut_range  , 
  const Ostap::EventIndex    first      ,
  const Ostap::EventIndex    last       ) const
{
  return get_stat ( data                ,
		    poly                ,
		    expression          ,
		    selection           ,
		    cut_range           , 
		    first               ,
		    last                ,
		    poly.xmin ()        ,
		    poly.xmax ()        ) ;
}
// ============================================================================

// ============================================================================
// 1D-Legendre polynoials: on-flight parameterization
// ============================================================================

// ============================================================================
/*  Project data in 1D-polyom=ominal: on-flight parameterisation 
 *  @param data       (INPUT)  input data 
 *  @param poly       (UPDATE) the polynoimal 
 *  @param expression (INPUT)  expression 
 *  @param selection  (INPUT)  selection (as boolean)
 *  @param first      (INOPUT) the first event to process (inclusive) 
 *  @param last       (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project1
( TTree*                     data       ,
  Ostap::Math::Bernstein&    poly       ,
  const std::string&         expression ,
  const std::string&         selection  ,
  const Ostap::EventIndex    first      ,
  const Ostap::EventIndex    last       ) const
{
  return get_stat ( data                ,
		    poly                ,
		    expression          ,
		    selection           ,
		    first               ,
		    last                ,
		    poly.xmin ()        ,
		    poly.xmax ()        ) ;
}
// =========================================================================
/*  Project data in 1D-polyom=ominal: on-flight parameterisation 
 *  @param data       (INPUT)  input data 
 *  @param poly       (UPDATE) the polynoimal 
 *  @param expression (INPUT)  expression 
 *  @param selection  (INPUT)  selection (as boolean)
 *  @param first      (INOPUT) the first event to process (inclusive) 
 *  @param last       (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project1
( const RooAbsData*          data       ,
  Ostap::Math::Bernstein&    poly       ,
  const std::string&         expression ,
  const std::string&         selection  ,
  const std::string&         cut_range  , 
  const Ostap::EventIndex    first      ,
  const Ostap::EventIndex    last       ) const
{
  return get_stat ( data                ,
		    poly                ,
		    expression          ,
		    selection           ,
		    cut_range           , 
		    first               ,
		    last                ,
		    poly.xmin ()        ,
		    poly.xmax ()        ) ;
}
// =========================================================================


// =========================================================================
// 2D-Bernstein polynomial: on-flight parameterisation
// =========================================================================
/*  Project data in 2D-polyoominal: on-flight parameterisation 
 *  @param data        (INPUT)  input data 
 *  @param poly        (UPDATE) the polynoimal 
 *  @param expression1 (INPUT)  expression1 
 *  @param expression2 (INPUT)  expression2 
 *  @param selection   (INPUT)  selection (as boolean)
 *  @param first       (INOPUT) the first event to process (inclusive) 
 *  @param last        (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// =========================================================================
Ostap::StatusCode
Ostap::Project::project2
( TTree*                     data        ,
  Ostap::Math::Bernstein2D&  poly        ,
  const std::string&         expression1 ,
  const std::string&         expression2 ,
  const std::string&         selection   ,
  const Ostap::EventIndex    first       ,
  const Ostap::EventIndex    last        ) const
{
  return get_stat ( data                ,
		    poly                ,
		    expression1         ,
		    expression2         ,
		    selection           ,
		    first               ,
		    last                ,
		    poly.xmin ()        ,
		    poly.xmax ()        , 
		    poly.ymin ()        ,
		    poly.ymax ()        ) ;
}
// =========================================================================
/*  Project data in 2D-polyoominal: on-flight parameterisation 
 *  @param data        (INPUT)  input data 
 *  @param poly        (UPDATE) the polynoimal 
 *  @param expression1 (INPUT)  expression1 
 *  @param expression2 (INPUT)  expression2 
 *  @param selection   (INPUT)  selection (as boolean)
 *  @param first       (INOPUT) the first event to process (inclusive) 
 *  @param last        (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// =========================================================================
Ostap::StatusCode
Ostap::Project::project2
( const RooAbsData*          data        ,
  Ostap::Math::Bernstein2D&  poly        ,
  const std::string&         expression1 ,
  const std::string&         expression2 ,
  const std::string&         selection   ,
  const std::string&         cut_range   , 
  const Ostap::EventIndex    first       ,
  const Ostap::EventIndex    last        ) const
{
  return get_stat ( data                ,
		    poly                ,
		    expression1         ,
		    expression2         ,
		    selection           ,
		    cut_range           , 
		    first               ,
		    last                ,
		    poly.xmin ()        ,
		    poly.xmax ()        , 
		    poly.ymin ()        ,
		    poly.ymax ()        ) ;
}
// =========================================================================


// =========================================================================
// 2D-Legendre polynomial: on-flight parameterisation
// =========================================================================
/*  Project data in 2D-polyoominal: on-flight parameterisation 
 *  @param data        (INPUT)  input data 
 *  @param poly        (UPDATE) the polynoimal 
 *  @param expression1 (INPUT)  expression1 
 *  @param expression2 (INPUT)  expression2 
 *  @param selection   (INPUT)  selection (as boolean)
 *  @param first       (INOPUT) the first event to process (inclusive) 
 *  @param last        (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// =========================================================================
Ostap::StatusCode
Ostap::Project::project2
( TTree*                     data        ,
  Ostap::Math::LegendreSum2& poly        ,
  const std::string&         expression1 ,
  const std::string&         expression2 ,
  const std::string&         selection   ,
  const Ostap::EventIndex    first       ,
  const Ostap::EventIndex    last        ) const
{
  return get_stat ( data                ,
		    poly                ,
		    expression1         ,
		    expression2         ,
		    selection           ,
		    first               ,
		    last                ,
		    poly.xmin ()        ,
		    poly.xmax ()        , 
		    poly.ymin ()        ,
		    poly.ymax ()        ) ;
}
// =========================================================================
/*  Project data in 2D-polyoominal: on-flight parameterisation 
 *  @param data        (INPUT)  input data 
 *  @param poly        (UPDATE) the polynoimal 
 *  @param expression1 (INPUT)  expression1 
 *  @param expression2 (INPUT)  expression2 
 *  @param selection   (INPUT)  selection (as boolean)
 *  @param first       (INOPUT) the first event to process (inclusive) 
 *  @param last        (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// =========================================================================
Ostap::StatusCode
Ostap::Project::project2
( const RooAbsData*          data        ,
  Ostap::Math::LegendreSum2& poly        ,
  const std::string&         expression1 ,
  const std::string&         expression2 ,
  const std::string&         selection   ,
  const std::string&         cut_range   , 
  const Ostap::EventIndex    first       ,
  const Ostap::EventIndex    last        ) const
{
  return get_stat ( data                ,
		    poly                ,
		    expression1         ,
		    expression2         ,
		    selection           ,
		    cut_range           , 
		    first               ,
		    last                ,
		    poly.xmin ()        ,
		    poly.xmax ()        , 
		    poly.ymin ()        ,
		    poly.ymax ()        ) ;
}
// ============================================================================


// ============================================================================
// 3D-Bernstein polynomials 
// ============================================================================
/** Project data in 2D-polyoominal: on-flight parameterisation 
 *  @param data        (INPUT)  input data 
 *  @param poly        (UPDATE) the polynoimal 
 *  @param expression1 (INPUT)  expression1 
 *  @param expression2 (INPUT)  expression2 
 *  @param expression3 (INPUT)  expression2 
 *  @param selection   (INPUT)  selection (as boolean)
 *  @param first       (INOPUT) the first event to process (inclusive) 
 *  @param last        (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project3
( TTree*                     data        ,
  Ostap::Math::Bernstein3D&  poly        ,
  const std::string&         expression1 ,
  const std::string&         expression2 ,
  const std::string&         expression3 ,
  const std::string&         selection   ,
  const Ostap::EventIndex    first       ,
  const Ostap::EventIndex    last        ) const
{
  return get_stat ( data                ,
		    poly                ,
		    expression1         ,
		    expression2         ,
		    expression3         ,
		    selection           ,
		    first               ,
		    last                ,
		    poly.xmin ()        ,
		    poly.xmax ()        , 
		    poly.ymin ()        ,
		    poly.ymax ()        , 
		    poly.zmin ()        ,
		    poly.zmax ()        ) ;
}
// =========================================================================
/*  Project data in 2D-polyoominal: on-flight parameterisation 
 *  @param data        (INPUT)  input data 
 *  @param poly        (UPDATE) the polynoimal 
 *  @param expression1 (INPUT)  expression1 
 *  @param expression2 (INPUT)  expression2 
 *  @param selection   (INPUT)  selection (as boolean)
 *  @param first       (INOPUT) the first event to process (inclusive) 
 *  @param last        (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// =========================================================================
Ostap::StatusCode
Ostap::Project::project3
( const RooAbsData*          data        ,
  Ostap::Math::Bernstein3D&  poly        ,
  const std::string&         expression1 ,
  const std::string&         expression2 ,
  const std::string&         expression3 ,
  const std::string&         selection   ,
  const std::string&         cut_range   , 
  const Ostap::EventIndex    first       ,
  const Ostap::EventIndex    last        ) const
{
  return get_stat ( data                ,
		    poly                ,
		    expression1         ,
		    expression2         ,
		    expression3         ,
		    selection           ,
		    cut_range           , 
		    first               ,
		    last                ,
		    poly.xmin ()        ,
		    poly.xmax ()        , 
		    poly.ymin ()        ,
		    poly.ymax ()        ,
		    poly.zmin ()        ,
		    poly.zmax ()        ) ;
}
// ============================================================================

// ============================================================================
// 3D-Legendre polynomials 
// ============================================================================
/** Project data in 2D-polyoominal: on-flight parameterisation 
 *  @param data        (INPUT)  input data 
 *  @param poly        (UPDATE) the polynoimal 
 *  @param expression1 (INPUT)  expression1 
 *  @param expression2 (INPUT)  expression2 
 *  @param expression3 (INPUT)  expression2 
 *  @param selection   (INPUT)  selection (as boolean)
 *  @param first       (INOPUT) the first event to process (inclusive) 
 *  @param last        (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project3
( TTree*                     data        ,
  Ostap::Math::LegendreSum3& poly        ,
  const std::string&         expression1 ,
  const std::string&         expression2 ,
  const std::string&         expression3 ,
  const std::string&         selection   ,
  const Ostap::EventIndex    first       ,
  const Ostap::EventIndex    last        ) const
{
  return get_stat ( data                ,
		    poly                ,
		    expression1         ,
		    expression2         ,
		    expression3         ,
		    selection           ,
		    first               ,
		    last                ,
		    poly.xmin ()        ,
		    poly.xmax ()        , 
		    poly.ymin ()        ,
		    poly.ymax ()        , 
		    poly.zmin ()        ,
		    poly.zmax ()        ) ;
}
// =========================================================================
/*  Project data in 2D-polyoominal: on-flight parameterisation 
 *  @param data        (INPUT)  input data 
 *  @param poly        (UPDATE) the polynoimal 
 *  @param expression1 (INPUT)  expression1 
 *  @param expression2 (INPUT)  expression2 
 *  @param selection   (INPUT)  selection (as boolean)
 *  @param first       (INOPUT) the first event to process (inclusive) 
 *  @param last        (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// =========================================================================
Ostap::StatusCode
Ostap::Project::project3
( const RooAbsData*          data        ,
  Ostap::Math::LegendreSum3& poly        ,
  const std::string&         expression1 ,
  const std::string&         expression2 ,
  const std::string&         expression3 ,
  const std::string&         selection   ,
  const std::string&         cut_range   , 
  const Ostap::EventIndex    first       ,
  const Ostap::EventIndex    last        ) const
{
  return get_stat ( data                ,
		    poly                ,
		    expression1         ,
		    expression2         ,
		    expression3         ,
		    selection           ,
		    cut_range           , 
		    first               ,
		    last                ,
		    poly.xmin ()        ,
		    poly.xmax ()        , 
		    poly.ymin ()        ,
		    poly.ymax ()        ,
		    poly.zmin ()        ,
		    poly.zmax ()        ) ;
}
// ============================================================================


// ============================================================================
// 4D-Legendre polynomials 
// ============================================================================
/** Project data in 2D-polyoominal: on-flight parameterisation 
 *  @param data        (INPUT)  input data 
 *  @param poly        (UPDATE) the polynoimal 
 *  @param expression1 (INPUT)  expression1 
 *  @param expression2 (INPUT)  expression2 
 *  @param expression3 (INPUT)  expression2 
 *  @param selection   (INPUT)  selection (as boolean)
 *  @param first       (INOPUT) the first event to process (inclusive) 
 *  @param last        (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Project::project4
( TTree*                     data        ,
  Ostap::Math::LegendreSum4& poly        ,
  const std::string&         expression1 ,
  const std::string&         expression2 ,
  const std::string&         expression3 ,
  const std::string&         expression4 ,
  const std::string&         selection   ,
  const Ostap::EventIndex    first       ,
  const Ostap::EventIndex    last        ) const
{
  return get_stat ( data                ,
		    poly                ,
		    expression1         ,
		    expression2         ,
		    expression3         ,
		    expression4         ,
		    selection           ,
		    first               ,
		    last                ,
		    poly.xmin ()        ,
		    poly.xmax ()        , 
		    poly.ymin ()        ,
		    poly.ymax ()        , 
		    poly.zmin ()        ,
		    poly.zmax ()        , 
		    poly.umin ()        ,
		    poly.umax ()        ) ;
}
// =========================================================================
/*  Project data in 2D-polyoominal: on-flight parameterisation 
 *  @param data        (INPUT)  input data 
 *  @param poly        (UPDATE) the polynoimal 
 *  @param expression1 (INPUT)  expression1 
 *  @param expression2 (INPUT)  expression2 
 *  @param selection   (INPUT)  selection (as boolean)
 *  @param first       (INOPUT) the first event to process (inclusive) 
 *  @param last        (INOPUT) the last  event to process (exclusive)
 *  @return statis code 
 */
// =========================================================================
Ostap::StatusCode
Ostap::Project::project4
( const RooAbsData*          data        ,
  Ostap::Math::LegendreSum4& poly        ,
  const std::string&         expression1 ,
  const std::string&         expression2 ,
  const std::string&         expression3 ,
  const std::string&         expression4 ,
  const std::string&         selection   ,
  const std::string&         cut_range   , 
  const Ostap::EventIndex    first       ,
  const Ostap::EventIndex    last        ) const
{
  return get_stat ( data                ,
		    poly                ,
		    expression1         ,
		    expression2         ,
		    expression3         ,
		    expression4         ,
		    selection           ,
		    cut_range           , 
		    first               ,
		    last                ,
		    poly.xmin ()        ,
		    poly.xmax ()        , 
		    poly.ymin ()        ,
		    poly.ymax ()        ,
		    poly.zmin ()        ,
		    poly.zmax ()        , 
		    poly.umin ()        ,
		    poly.umax ()        ) ;
}
// ============================================================================




// ============================================================================
//                                                                      The END 
// ============================================================================
