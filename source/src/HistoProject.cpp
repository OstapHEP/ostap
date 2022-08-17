// ============================================================================
// Include files
// ============================================================================
// ROOT 
// ============================================================================
#include "RooDataSet.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
// ============================================================================
// Local: 
// ============================================================================
#include "Ostap/StatVar.h"
#include "Ostap/Formula.h"
#include "Ostap/FormulaVar.h"
#include "Ostap/HistoProject.h"
#include "Ostap/Iterator.h"
// ============================================================================
#include "OstapDataFrame.h"
#include "local_math.h"
#include "local_utils.h"
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
  static_assert (std::numeric_limits<unsigned long>::is_specialized   , 
                 "Numeric_limist<unsigned long> are not specialized!" ) ;
  // ========================================================================== 
  /// get variable by name from RooArgSet
  RooAbsReal* get_var ( const RooArgSet&   aset , 
                        const std::string& name ) 
  {
    RooAbsArg* arg = aset.find ( name.c_str() ) ;
    return 0 != arg ? dynamic_cast<RooAbsReal*> ( arg ) : nullptr ;
  }
  // ==========================================================================
}
// ============================================================================
/** make a projection of RooDataSet into the histogram 
 *  @param data  (INPUT)  input data 
 *  @param histo (UPDATE) histogram 
 *  @param expression (INPUT) expression
 *  @param selection  (INPUT) selection criteria/weight 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project
( const RooAbsData*   data       , 
  TH1*                histo      ,
  const RooAbsReal&   expression ,
  const RooAbsReal*   selection  ,
  const unsigned long first      ,
  const unsigned long last       ) 
{
  //
  if ( 0 == histo ) { return Ostap::StatusCode ( 301 ) ; }
  else { histo->Reset() ; } // reset the historgam 
  if ( 0 == data  ) { return Ostap::StatusCode ( 300 ) ; }
  //
  const unsigned long nEntries = 
    std::min ( last , (unsigned long) data->numEntries() ) ;
  if ( nEntries <= first  ) { return Ostap::StatusCode::RECOVERABLE ; }
  //
  const bool weighted = data->isWeighted() ;
  //
  const double xmin = histo->GetXaxis()->GetXmin () ;
  const double xmax = histo->GetXaxis()->GetXmax () ;
  //
  for ( unsigned long entry = first ; entry < nEntries ; ++entry )   
  {
    //
    if ( 0 == data->get( entry)  ) { break ; }                    // BREAK
    //
    // data weight 
    const double dw = weighted  ? data      -> weight () : 1.0 ;
    if ( !dw ) { continue ; }                         // SKIP    
    //
    // selection weight 
    const double sw = selection ? selection -> getVal () : 1.0 ;
    if ( !sw ) { continue ; }                         // SKIP    
    //
    // calculate the total weight 
    const double w = sw * dw ;
    if ( !w  ) { continue ; }                          // SKIP 
    //
    // calculate the x-value (only for non-zero weights)
    const double xvalue = expression.getVal()  ;
    //
    // check the range 
    if ( xmax <= xvalue || xvalue < xmin ) { continue ; }
    //
    // fill the histogram  (only for non-zero weights and in-range entries!)
    histo -> Fill ( xvalue , w ) ; 
    //
    if  ( weighted )
    {
      const double dwe = data -> weightError ( RooAbsData::SumW2 ) ;
      const double we  = ( dwe ? dwe : dw ) * sw ;
      if  ( !s_equal ( we , w ) )
      {
        const int    bin    = histo -> FindBin     ( xvalue ) ;
        const double binerr = histo -> GetBinError ( bin    ) ;
        const double err2   = binerr * binerr - w * w + we * we  ;
        histo -> SetBinError ( bin , std::sqrt ( err2 ) ) ;  
      }
    }
  }
  //
  return StatusCode::SUCCESS ;  
}
// ============================================================================
/*  make a projection of RooDataSet into the histogram 
 *  @param data  (INPUT)  input data 
 *  @param histo (UPDATE) histogram 
 *  @param xexpression (INPUT) expression for x-axis 
 *  @param yexpression (INPUT) expression for y-axis 
 *  @param selection  (INPUT) selection criteria/weight 
 */
// ============================================================================
Ostap::StatusCode
Ostap::HistoProject::project2
( const RooAbsData*   data        , 
  TH2*                histo       ,
  const RooAbsReal&   xexpression ,
  const RooAbsReal&   yexpression ,
  const RooAbsReal*   selection   ,
  const unsigned long first       ,
  const unsigned long last        ) 
{
  //
  if ( 0 == histo ) { return Ostap::StatusCode ( 301 ) ; }
  else { histo->Reset() ; } // reset the historgam 
  if ( 0 == data  ) { return Ostap::StatusCode ( 300 ) ; }
  //
  const unsigned long nEntries = 
    std::min ( last , (unsigned long) data->numEntries() ) ;
  if ( nEntries <= first  ) { return Ostap::StatusCode::RECOVERABLE ; }
  //
  const bool weighted = data->isWeighted() ;
  //
  const double xmin = histo -> GetXaxis () -> GetXmin () ;
  const double xmax = histo -> GetXaxis () -> GetXmax () ;
  const double ymin = histo -> GetYaxis () -> GetXmin () ;
  const double ymax = histo -> GetYaxis () -> GetXmax () ;
  //
  for ( unsigned long entry = first ; entry < nEntries ; ++entry )   
  {
    //
    if ( 0 == data->get( entry)  ) { break ; }                    // BREAK
    //
    // data weight 
    const double dw = weighted  ? data      -> weight () : 1.0 ;
    if ( !dw ) { continue ; }   // SKIP
    //
    // selection weight 
    const double sw = selection ? selection -> getVal () : 1.0 ;
    if ( !sw ) { continue ; }   // SKIP    
    //
    // calculate the total weight 
    const double w = sw * dw ;
    if ( !w  ) { continue ; }   // SKIP 
    //
    // calculate the x-value (only for non-zero weights)
    const double xvalue = xexpression.getVal()  ;
    // check the range 
    if ( xmax <= xvalue || xvalue < xmin ) { continue ; }
    // calculate the y-value (only for non-zero weights)
    const double yvalue = yexpression.getVal()  ;
    // check the range 
    if ( ymax <= yvalue || yvalue < ymin ) { continue ; }
    //
    // fill the histogram  (only for non-zero weights)
    histo->Fill ( xvalue , yvalue , w ) ; 
    //
    if  ( weighted )
    {
      const double dwe = data -> weightError ( RooAbsData::SumW2 ) ;
      const double we  = ( dwe ? dwe : dw ) * sw ;
      if  ( !s_equal ( we , w ) ) 
      {
        const int    bin    = histo -> FindBin     ( xvalue , yvalue ) ;
        const double binerr = histo -> GetBinError ( bin    ) ;
        const double err2   = binerr * binerr - w * w + we * we  ;
        histo -> SetBinError ( bin , std::sqrt ( err2 ) ) ;  
      }
    }
  }
  //
  return StatusCode::SUCCESS ;  
}
// ============================================================================
/*  make a projection of RooDataSet into the histogram 
 *  @param data  (INPUT)  input data 
 *  @param histo (UPDATE) histogram 
 *  @param xexpression (INPUT) expression for x-axis 
 *  @param yexpression (INPUT) expression for y-axis 
 *  @param zexpression (INPUT) expression for z-axis 
 *  @param selection  (INPUT) selection criteria/weight 
 */
// ============================================================================
Ostap::StatusCode
Ostap::HistoProject::project3
( const RooAbsData*   data        , 
  TH3*                histo       ,
  const RooAbsReal&   xexpression ,
  const RooAbsReal&   yexpression ,
  const RooAbsReal&   zexpression ,
  const RooAbsReal*   selection   ,
  const unsigned long first       ,
  const unsigned long last        ) 
{
  //
  if ( 0 == histo ) { return Ostap::StatusCode ( 301 ) ; }
  else { histo->Reset() ; } // reset the historgam 
  if ( 0 == data  ) { return Ostap::StatusCode ( 300 ) ; }
  //
  const unsigned long nEntries = 
    std::min ( last , (unsigned long) data->numEntries() ) ;
  if ( nEntries <= first  ) { return Ostap::StatusCode::RECOVERABLE ; }
  //
  const bool weighted = data->isWeighted() ;
  //
  const double xmin = histo -> GetXaxis () -> GetXmin () ;
  const double xmax = histo -> GetXaxis () -> GetXmax () ;
  const double ymin = histo -> GetYaxis () -> GetXmin () ;
  const double ymax = histo -> GetYaxis () -> GetXmax () ;
  const double zmin = histo -> GetZaxis () -> GetXmin () ;
  const double zmax = histo -> GetZaxis () -> GetXmax () ;
  //
  for ( unsigned long entry = first ; entry < nEntries ; ++entry )   
  {
    //
    if ( 0 == data->get( entry)  ) { break ; }                    // BREAK
    //
    // data weight 
    const double dw = weighted  ? data      -> weight () : 1.0 ;
    if ( !dw ) { continue ; }    // SKIP 
    //
    // selection weight 
    const double sw = selection ? selection -> getVal () : 1.0 ;
    if ( !sw ) { continue ; }    // SKIP 
    //
    // calculate the total weight 
    const double w = sw * dw ;
    if ( !w  ) { continue ; }    // SKIP 
    //
    // calculate the x-value (only for non-zero weights)
    const double xvalue = xexpression.getVal()  ;
    if ( xmax <= xvalue || xvalue < xmin ) { continue ; }
    // calculate the y-value (only for non-zero weights)
    const double yvalue = yexpression.getVal()  ;
    if ( ymax <= yvalue || yvalue < ymin ) { continue ; }
    // calculate the z-value (only for non-zero weights)
    const double zvalue = zexpression.getVal()  ;
    if ( zmax <= zvalue || zvalue < zmin ) { continue ; }
    //
    // fill the histogram  (only for non-zero weights)
    histo->Fill ( xvalue , yvalue , zvalue , w ) ; 
    //
    if  ( weighted )
    {
      const double dwe = data -> weightError ( RooAbsData::SumW2 ) ;
      const double we  = ( dwe ? dwe : dw ) * sw ;
      if  ( !s_equal ( we , w ) )
      {
        const int    bin    = histo -> FindBin     ( xvalue , yvalue , zvalue ) ;
        const double binerr = histo -> GetBinError ( bin    ) ;
        const double err2   = binerr * binerr - w * w + we * we  ;
        histo -> SetBinError ( bin , std::sqrt ( err2 ) ) ;  
      }
    }
  }
  //
  return StatusCode::SUCCESS ;  
}
// ============================================================================
/* make a projection of RooDataSet into the histogram 
 * @param data  (INPUT)  input data 
 * @param histo (UPDATE) histogram 
 * @param expression (INPUT) expression
 * @param selection  (INPUT) selection criteria/weight 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project
( const RooAbsData*   data       , 
  TH1*                histo      ,
  const std::string&  expression ,
  const std::string&  selection  ,
  const unsigned long first      ,                                          
  const unsigned long last       ) 
{
  //
  if ( 0 == histo ) { return Ostap::StatusCode ( 301 ) ; }
  else { histo->Reset() ; } // reset the historgam 
  if ( 0 == data  ) { return Ostap::StatusCode ( 300 ) ; }
  //
  const unsigned long nEntries = 
    std::min ( last , (unsigned long) data->numEntries() ) ;
  if ( nEntries <= first  ) { return Ostap::StatusCode::RECOVERABLE ; }
  //
  RooArgList        alst ;
  const RooArgSet*  aset = data->get() ;
  if ( 0 == aset       ) { return  0 ; }                          // RETURN
  Ostap::Utils::Iterator iter ( *aset );
  RooAbsArg*   coef = 0 ;
  while ( ( coef = (RooAbsArg*) iter.next() ) ){ alst.add ( *coef ); }
  //
  // convert expressions into FormulaVar 
  const bool        with_cuts   = !selection.empty() ;
  const RooAbsReal* cut_var     = 0 ;
  if ( with_cuts ) { cut_var =get_var ( *aset , selection ) ; }
  std::unique_ptr<Ostap::FormulaVar> cuts ;
  if ( with_cuts && 0 == cut_var ) 
  {
    cuts.reset ( new Ostap::FormulaVar( selection , alst , false ) ) ;
    if ( !cuts->ok () ) { return Ostap::StatusCode(302) ; } //        // RETURN 
  }
  //
  const RooAbsReal* x_var = get_var ( *aset , expression ) ;
  std::unique_ptr<Ostap::FormulaVar> xwhat ;
  if ( 0 == x_var ) 
  {
    xwhat.reset( new Ostap::FormulaVar( expression , alst , false ) ) ;
    if ( !xwhat->ok()   ) { return Ostap::StatusCode(303)  ; }             // RETURN
  }
  //
  return project ( data                                   , 
                   histo                                  , 
                   0 !=   x_var ?   *x_var : *xwhat       , 
                   0 != cut_var ?  cut_var :   cuts.get() , first , last ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into 2D-histogram 
 *  @param data  (INPUT)  input data 
 *  @param histo (UPDATE) histogram 
 *  @param xexpression (INPUT) expression for x-axis 
 *  @param yexpression (INPUT) expression for y-axis 
 *  @param selection  (INPUT) selection criteria/weight 
 */
// ============================================================================
Ostap::StatusCode
Ostap::HistoProject::project2
( const RooAbsData*   data        , 
  TH2*                histo       ,
  const std::string&  xexpression ,
  const std::string&  yexpression ,
  const std::string&  selection   ,
  const unsigned long first       ,
  const unsigned long last        ) 
{
  //
  if ( 0 == histo ) { return Ostap::StatusCode ( 301 ) ; }
  else { histo->Reset() ; } // reset the historgam 
  if ( 0 == data  ) { return Ostap::StatusCode ( 300 ) ; }
  //
  const unsigned long nEntries = 
    std::min ( last , (unsigned long) data->numEntries() ) ;
  if ( nEntries <= first  ) { return Ostap::StatusCode::SUCCESS ; }
  //
  RooArgList        alst ;
  const RooArgSet*  aset = data->get() ;
  if ( 0 == aset       ) { return  0 ; }  // RETURN
  //
  Ostap::Utils::Iterator iter ( *aset );
  RooAbsArg*   coef = 0 ;
  while ( ( coef = (RooAbsArg*) iter.next() ) ){ alst.add ( *coef ); }
  //
  // convert expressions into FormulaVar 
  const bool        with_cuts   = !selection.empty() ;
  const RooAbsReal* cut_var     = 0 ;
  if ( with_cuts ) { cut_var =get_var ( *aset , selection ) ; }
  std::unique_ptr<Ostap::FormulaVar> cuts ;
  if ( with_cuts && 0 == cut_var ) 
  {
    cuts.reset ( new Ostap::FormulaVar ( selection , alst , false ) ) ;
    if ( !cuts->ok () ) { return Ostap::StatusCode(302) ; } //        // RETURN 
  }
  //
  const RooAbsReal* x_var = get_var ( *aset , xexpression ) ;
  std::unique_ptr<Ostap::FormulaVar> xwhat ;
  if ( 0 == x_var ) 
  {
    xwhat.reset( new Ostap::FormulaVar( xexpression , alst , false ) ) ;
    if ( !xwhat->ok()   ) { return Ostap::StatusCode(303)  ; }             // RETURN
  }
  //
  const RooAbsReal* y_var = get_var ( *aset , yexpression ) ;
  std::unique_ptr<Ostap::FormulaVar> ywhat ;
  if ( 0 == y_var ) 
  { 
    ywhat.reset( new Ostap::FormulaVar( yexpression , alst , false ) ) ;
    if ( !ywhat->ok()   ) { return Ostap::StatusCode(304)  ; }             // RETURN
  }
  //
  return project2 ( data                                   , 
                    histo                                  , 
                    0 !=   x_var ?   *x_var : *xwhat       , 
                    0 !=   y_var ?   *y_var : *ywhat       , 
                    0 != cut_var ?  cut_var :   cuts.get() , first , last ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into the histogram 
 *  @param data  (INPUT)  input data 
 *  @param histo (UPDATE) histogram 
 *  @param xexpression (INPUT) expression for x-axis 
 *  @param yexpression (INPUT) expression for y-axis 
 *  @param zexpression (INPUT) expression for z-axis 
 *  @param selection  (INPUT) selection criteria/weight 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project3
( const RooAbsData*   data        , 
  TH3*                histo       ,
  const std::string&  xexpression ,
  const std::string&  yexpression ,
  const std::string&  zexpression ,
  const std::string&  selection   ,
  const unsigned long first       ,
  const unsigned long last        ) 
{
  //
  if ( 0 == histo ) { return Ostap::StatusCode ( 301 ) ; }
  else { histo->Reset() ; } // reset the historgam 
  if ( 0 == data  ) { return Ostap::StatusCode ( 300 ) ; }
  //
  const unsigned long nEntries = 
    std::min ( last , (unsigned long) data->numEntries() ) ;
  if ( nEntries <= first  ) { return Ostap::StatusCode::SUCCESS ; }
  //
  RooArgList        alst ;
  const RooArgSet*  aset = data->get() ;
  if ( 0 == aset       ) { return  0 ; }  // RETURN
  //
  Ostap::Utils::Iterator iter ( *aset );
  RooAbsArg*   coef = 0 ;
  while ( ( coef = (RooAbsArg*) iter.next() ) ){ alst.add ( *coef ); }
  //
  // convert expressions into FormulaVar 
  const bool        with_cuts   = !selection.empty() ;
  const RooAbsReal* cut_var     = 0 ;
  if ( with_cuts ) { cut_var =get_var ( *aset , selection ) ; }
  std::unique_ptr<Ostap::FormulaVar> cuts ;
  if ( with_cuts && 0 == cut_var ) 
  {
    cuts.reset ( new Ostap::FormulaVar( selection , alst , false ) ) ;
    if ( !cuts->ok () ) { return Ostap::StatusCode(302) ; } //        // RETURN 
  }
  //
  const RooAbsReal* x_var = get_var ( *aset , xexpression ) ;
  std::unique_ptr<Ostap::FormulaVar> xwhat ;
  if ( 0 == x_var ) 
  {
    xwhat.reset( new Ostap::FormulaVar( xexpression , alst , false ) ) ;
    if ( !xwhat->ok()   ) { return Ostap::StatusCode(303)  ; }             // RETURN
  }
  //
  const RooAbsReal* y_var = get_var ( *aset , yexpression ) ;
  std::unique_ptr<Ostap::FormulaVar> ywhat ;
  if ( 0 == y_var ) 
  {
    ywhat.reset( new Ostap::FormulaVar( yexpression , alst , false ) ) ;
    if ( !ywhat->ok()   ) { return Ostap::StatusCode(304)  ; }             // RETURN
  }
  //
  const RooAbsReal* z_var = get_var ( *aset , zexpression ) ;
  std::unique_ptr<Ostap::FormulaVar> zwhat ;
  if ( 0 == z_var ) 
  {
    zwhat.reset( new Ostap::FormulaVar( zexpression , alst , false ) ) ;
    if ( !zwhat->ok()   ) { return Ostap::StatusCode(305)  ; }             // RETURN
  }
  //
  return project3 ( data                                   , 
                    histo                                  , 
                    0 !=   x_var ?   *x_var : *xwhat       , 
                    0 !=   y_var ?   *y_var : *ywhat       , 
                    0 !=   z_var ?   *z_var : *zwhat       , 
                    0 != cut_var ?  cut_var :   cuts.get() , first , last ) ;
}
// ============================================================================
/*  make a projection of DataFrame into the histogram 
 *  @param data  (INPUT)  input data 
 *  @param histo (UPDATE) histogram 
 *  @param expression (INPUT) expression
 *  @param selection  (INPUT) selection criteria/weight 
 */
// ============================================================================
Ostap::StatusCode Ostap::HistoProject::project
( Ostap::FrameNode    data       , 
  TH1*                histo      ,
  const std::string&  expression ,
  const std::string&  selection  ) 
{
  //
  if ( nullptr  == histo ) { return Ostap::StatusCode ( 301 ) ; }
  else if ( dynamic_cast<TH2*>( histo ) ) { return Ostap::StatusCode ( 301 ) ; }  
  else { histo->Reset() ; } // reset the histogram 
  //
  TH1D model {} ; histo->Copy ( model ) ;
  //
  const bool no_cuts = trivial ( selection ) ;
  //
  const std::string xvar   = Ostap::tmp_name ( "vx_" , expression ) ;
  const std::string weight = Ostap::tmp_name ( "w_"  , selection  ) ;
  //
  auto h = data
    .Define  ( xvar   ,                   "1.0*(" + expression + ")" )
    .Define  ( weight , no_cuts ? "1.0" : "1.0*(" + selection  + ")" ) 
    .Histo1D ( model  , xvar , weight ) ;
  //
  h->Copy ( *histo ) ;
  //
  return Ostap::StatusCode::SUCCESS ;
}
// ========================================================================
/*  make a projection of RooDataSet into the histogram 
 *  @param data  (INPUT)  input data 
 *  @param histo (UPDATE) histogram 
 *  @param xexpression (INPUT) expression for x-axis 
 *  @param yexpression (INPUT) expression for y-axis 
 *  @param selection  (INPUT) selection criteria/weight 
 */
// ========================================================================
Ostap::StatusCode Ostap::HistoProject::project2
( Ostap::FrameNode    data        , 
  TH2*                histo       ,
  const std::string&  xexpression ,
  const std::string&  yexpression ,
  const std::string&  selection   )
{
  //
  if      ( nullptr == histo            ) { return Ostap::StatusCode ( 301 ) ; }
  else if ( dynamic_cast<TH3*>( histo ) ) { return Ostap::StatusCode ( 301 ) ; }  
  else { histo->Reset() ; } // reset the historgam 
  //
  const bool no_cuts = trivial ( selection ) ;
  //
  const std::string xvar   = Ostap::tmp_name ( "vx_" , xexpression ) ;
  const std::string yvar   = Ostap::tmp_name ( "vy_" , yexpression ) ;
  const std::string weight = Ostap::tmp_name ( "w_"  , selection   ) ;
  //
  TH2D model {} ; histo->Copy ( model ) ;
  //
  auto h = data
    .Define  ( xvar   ,                   "1.0*(" + xexpression + ")" )
    .Define  ( yvar   ,                   "1.0*(" + yexpression + ")" )
    .Define  ( weight , no_cuts ? "1.0" : "1.0*(" + selection   + ")" ) 
    .Histo2D ( model  , xvar , yvar , weight ) ;
  //
  h->Copy ( *histo ) ;
  //
  return Ostap::StatusCode::SUCCESS ;
}
// ========================================================================
/*  make a projection of RooDataSet into the histogram 
 *  @param data        (INPUT)  input data 
 *  @param histo       (UPDATE) histogram 
 *  @param xexpression (INPUT)  expression for x-axis 
 *  @param yexpression (INPUT)  expression for y-axis 
 *  @param zexpression (INPUT)  expression for z-axis 
 *  @param selection   (INPUT)  selection criteria/weight 
 */
// ========================================================================
Ostap::StatusCode Ostap::HistoProject::project3
( Ostap::FrameNode    data        , 
  TH3*                histo       ,
  const std::string&  xexpression ,
  const std::string&  yexpression ,
  const std::string&  zexpression ,
  const std::string&  selection   )
{
  //
  if ( nullptr  == histo ) { return Ostap::StatusCode ( 301 ) ; }
  else { histo->Reset() ; } // reset the historgam 
  //
  const bool no_cuts = trivial ( selection  ) ; 
  //
  const std::string xvar   = Ostap::tmp_name ( "vx_" , xexpression ) ;
  const std::string yvar   = Ostap::tmp_name ( "vy_" , yexpression ) ;
  const std::string zvar   = Ostap::tmp_name ( "vz_" , yexpression ) ;
  const std::string weight = Ostap::tmp_name ( "w_"  , selection   ) ;
  //
  //
  TH3D model {} ; histo->Copy ( model ) ;
  //
  auto h = data
    .Define  ( xvar   ,                   "1.0*(" + xexpression + ")" )
    .Define  ( yvar   ,                   "1.0*(" + yexpression + ")" )
    .Define  ( zvar   ,                   "1.0*(" + zexpression + ")" )
    .Define  ( weight , no_cuts ? "1.0" : "1.0*(" + selection   + ")" ) 
    .Histo3D ( model  , xvar , yvar , zvar ,  weight ) ;
  //
  h->Copy ( *histo ) ;
  //
  return Ostap::StatusCode::SUCCESS ;
}

// ============================================================================
//                                                                      The END 
// ============================================================================
