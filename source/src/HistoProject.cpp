// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <cstring>
// ============================================================================
// ROOT
// ============================================================================
#include "RooDataSet.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
// ============================================================================
// Local: 
// ============================================================================
#include "Ostap/StatVar.h"
#include "Ostap/Formula.h"
#include "Ostap/FormulaVar.h"
#include "Ostap/HistoProject.h"
#include "Ostap/Iterator.h"
#include "Ostap/ProgressBar.h"
#include "Ostap/Notifier.h"
#include "Ostap/Exception.h"
#include "Ostap/Params.h"
#include "Ostap/Parameterization.h"
#include "Ostap/Polynomials.h"
#include "Ostap/Bernstein.h"
#include "Ostap/Bernstein2D.h"
#include "Ostap/Bernstein3D.h"
// ============================================================================
#include "OstapDataFrame.h"
#include "Exception.h"
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
  RooAbsReal* get_var
  ( const RooArgSet&   aset , 
    const std::string& name ) 
  {
    RooAbsArg* arg = aset.find ( name.c_str() ) ;
    return 0 != arg ? dynamic_cast<RooAbsReal*> ( arg ) : nullptr ;
  }
  // ==========================================================================
  // RooAbsData 
  // ==========================================================================
  /// project data into object 
  template <class WHAT>
  Ostap::StatusCode _project_ 
  ( const RooAbsData*                 data       , 
    const Ostap::Utils::ProgressConf& progress   ,
    WHAT&                             what       , 
    const RooAbsReal&                 expression ,
    const RooAbsReal*                 selection  ,
    const char*                       range      , 
    const double                      xmin       ,
    const double                      xmax       ,
    const unsigned long               first      ,
    const unsigned long               last       ) 
  {
    if ( 0 == data  ) { return Ostap::StatusCode ( 300 ) ; }
    //
    const unsigned long nEntries = 
      std::min ( last , (unsigned long) data->numEntries() ) ;
    if ( nEntries <= first  ) { return Ostap::StatusCode::RECOVERABLE ; }
    //
    const bool weighted  = data->isWeighted() ;
    const bool has_range = range && 1 <= std::strlen ( range ) ;
    //
    Ostap::Utils::ProgressBar bar ( nEntries - first , progress ) ;
    for ( unsigned long entry = first ; entry < nEntries ; ++entry , ++bar )   
    {
      //
      const RooArgSet* vars = data->get( entry ) ;
      if ( nullptr == vars )                          { break    ; } // BREAK 
      //
      if ( has_range && !vars->allInRange ( range ) ) { continue ; } // CONTINUE    
      //
      // data weight 
      const double dw = weighted  ? data      -> weight () : 1.0 ;
      if ( !dw ) { continue ; }                                      // SKIP    
      //
      // selection weight 
      const double sw = selection ? selection -> getVal () : 1.0 ;
      if ( !sw ) { continue ; }                                      // SKIP    
      //
      // calculate the total weight 
      const double w = sw * dw ;
      if ( !w  ) { continue ; }                                      // SKIP 
      //
      // calculate the x-value (only for non-zero weights)
      const double xvalue = expression.getVal()  ;
      //
      // check the range 
      if ( xmax <= xvalue || xvalue < xmin ) { continue ; }
      //
      // fill the histogram  (only for non-zero weights and in-range entries!)
      what.fill ( xvalue , w ) ; 
      //
    }
    //
    return Ostap::StatusCode::SUCCESS ;
  } ;
  // ==========================================================================
  /// project data into object 
  template <class WHAT>
  Ostap::StatusCode _project_ 
  ( const RooAbsData*                 data       , 
    const Ostap::Utils::ProgressConf& progress   ,
    WHAT&                             what       , 
    const std::string&                expression ,
    const std::string&                selection  ,
    const char*                       range      , 
    const double                      xmin       ,
    const double                      xmax       ,
    const unsigned long               first      ,
    const unsigned long               last       ) 
  {
    if ( 0 == data  ) { return Ostap::StatusCode ( 300 ) ; }
    //
    RooArgList        alst ;
    const RooArgSet*  aset = data->get() ;
    if ( 0 == aset       ) { return  0 ; }                            // RETURN
    Ostap::Utils::Iterator iter ( *aset );
    RooAbsArg*   coef = 0 ;
    while ( ( coef = (RooAbsArg*) iter.next() ) ){ alst.add ( *coef ); }
    //
    // convert expressions into FormulaVar 
    const bool        with_cuts = !selection.empty() ;
    const RooAbsReal* cut_var   = 0 ;
    if ( with_cuts ) { cut_var  = get_var ( *aset , selection ) ; }
    std::unique_ptr<Ostap::FormulaVar> cuts ;
    if ( with_cuts && 0 == cut_var ) 
    {
      cuts.reset ( new Ostap::FormulaVar ( selection , alst , false ) ) ;
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
    return _project_ ( data                                  , 
                       progress                              ,
                       what                                  , 
                       0 !=   x_var ?   *x_var : *xwhat      , 
                       0 != cut_var ?  cut_var :  cuts.get() , 
                       range , xmin , xmax , first , last    ) ;
  } ;
  // ==========================================================================
  /// project data into object 
  template <class WHAT>
  Ostap::StatusCode _project2_ 
  ( const RooAbsData*                 data        , 
    const Ostap::Utils::ProgressConf& progress    ,
    WHAT&                             what        , 
    const RooAbsReal&                 xexpression ,
    const RooAbsReal&                 yexpression ,
    const RooAbsReal*                 selection   ,
    const char*                       range       , 
    const double                      xmin        ,
    const double                      xmax        ,
    const double                      ymin        ,
    const double                      ymax        ,
    const unsigned long               first       ,
    const unsigned long               last        ) 
  {
    if ( 0 == data  ) { return Ostap::StatusCode ( 300 ) ; }
    //
    const unsigned long nEntries = 
      std::min ( last , (unsigned long) data->numEntries() ) ;
    if ( nEntries <= first  ) { return Ostap::StatusCode::RECOVERABLE ; }
    //
    const bool weighted  = data->isWeighted() ;
    const bool has_range = range && 1 <= std::strlen ( range ) ;
    //
    Ostap::Utils::ProgressBar bar ( nEntries - first , progress ) ;
    for ( unsigned long entry = first ; entry < nEntries ; ++entry , ++bar )   
    {
      //
      const RooArgSet* vars = data->get( entry ) ;
      if ( nullptr == vars )                          { break    ; } // BREAK 
      //
      if ( has_range && !vars->allInRange ( range ) ) { continue ; } // CONTINUE    
      //
      // data weight 
      const double dw = weighted  ? data      -> weight () : 1.0 ;
      if ( !dw ) { continue ; }                                      // SKIP    
      //
      // selection weight 
      const double sw = selection ? selection -> getVal () : 1.0 ;
      if ( !sw ) { continue ; }                                      // SKIP    
      //
      // calculate the total weight 
      const double w = sw * dw ;
      if ( !w  ) { continue ; }                                      // SKIP 
      //
      // calculate the x-value (only for non-zero weights)
      const double xvalue = xexpression.getVal()  ;
      // check the x-range 
      if ( xmax <= xvalue || xvalue < xmin ) { continue ; }
      // calculate the y-value (only for non-zero weights)
      const double yvalue = yexpression.getVal()  ;
      // check the y-range 
      if ( ymax <= yvalue || yvalue < ymin ) { continue ; }
      //
      // fill the histogram  (only for non-zero weights and in-range entries!)
      what.fill ( xvalue , yvalue , w ) ; 
      //
    }
    //
    return Ostap::StatusCode::SUCCESS ;
  } ;
  // ==========================================================================
  /// project data into object 
  template <class WHAT>
  Ostap::StatusCode _project2_ 
  ( const RooAbsData*                 data        , 
    const Ostap::Utils::ProgressConf& progress    ,
    WHAT&                             what        , 
    const std::string&                xexpression ,
    const std::string&                yexpression ,
    const std::string&                selection   ,
    const char*                       range       , 
    const double                      xmin        ,
    const double                      xmax        ,
    const double                      ymin        ,
    const double                      ymax        ,
    const unsigned long               first       ,
    const unsigned long               last        ) 
  {
    if ( 0 == data  ) { return Ostap::StatusCode ( 300 ) ; }
    //
    RooArgList        alst ;
    const RooArgSet*  aset = data->get() ;
    if ( 0 == aset       ) { return  0 ; }                            // RETURN
    Ostap::Utils::Iterator iter ( *aset );
    RooAbsArg*   coef = 0 ;
    while ( ( coef = (RooAbsArg*) iter.next() ) ){ alst.add ( *coef ); }
    //
    // convert expressions into FormulaVar 
    const bool        with_cuts = !selection.empty() ;
    const RooAbsReal* cut_var   = 0 ;
    if ( with_cuts ) { cut_var  = get_var ( *aset , selection ) ; }
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
    return _project2_ ( data                                  , 
                        progress                              ,
                        what                                  , 
                        0 !=   x_var ?   *x_var : *xwhat      , 
                        0 !=   y_var ?   *y_var : *ywhat      , 
                        0 != cut_var ?  cut_var :  cuts.get() , 
                        range ,
                        xmin  , xmax , 
                        ymin  , ymax , 
                        first , last ) ;
  } ;
  // ==========================================================================
  /// project data into object 
  template <class WHAT>
  Ostap::StatusCode _project3_ 
  ( const RooAbsData*                 data        , 
    const Ostap::Utils::ProgressConf& progress    ,
    WHAT&                             what        , 
    const RooAbsReal&                 xexpression ,
    const RooAbsReal&                 yexpression ,
    const RooAbsReal&                 zexpression ,
    const RooAbsReal*                 selection   ,
    const char*                       range       , 
    const double                      xmin        ,
    const double                      xmax        ,
    const double                      ymin        ,
    const double                      ymax        ,
    const double                      zmin        ,
    const double                      zmax        ,
    const unsigned long               first       ,
    const unsigned long               last        ) 
  {
    if ( 0 == data  ) { return Ostap::StatusCode ( 300 ) ; }
    //
    const unsigned long nEntries = 
      std::min ( last , (unsigned long) data->numEntries() ) ;
    if ( nEntries <= first  ) { return Ostap::StatusCode::RECOVERABLE ; }
    //
    const bool weighted  = data->isWeighted() ;
    const bool has_range = range && 1 <= std::strlen ( range ) ;
    //
    Ostap::Utils::ProgressBar bar ( nEntries - first , progress ) ;
    for ( unsigned long entry = first ; entry < nEntries ; ++entry , ++bar )   
    {
      //
      const RooArgSet* vars = data->get( entry ) ;
      if ( nullptr == vars )                          { break    ; } // BREAK 
      //
      if ( has_range && !vars->allInRange ( range ) ) { continue ; } // CONTINUE    
      //
      // data weight 
      const double dw = weighted  ? data      -> weight () : 1.0 ;
      if ( !dw ) { continue ; }                                      // SKIP    
      //
      // selection weight 
      const double sw = selection ? selection -> getVal () : 1.0 ;
      if ( !sw ) { continue ; }                                      // SKIP    
      //
      // calculate the total weight 
      const double w = sw * dw ;
      if ( !w  ) { continue ; }                                      // SKIP 
      //
      // calculate the x-value (only for non-zero weights)
      const double xvalue = xexpression.getVal()  ;
      // check the x-range 
      if ( xmax <= xvalue || xvalue < xmin ) { continue ; }
      // calculate the y-value (only for non-zero weights)
      const double yvalue = yexpression.getVal()  ;
      // check the y-range 
      if ( ymax <= yvalue || yvalue < ymin ) { continue ; }
      // calculate the z-value (only for non-zero weights)
      const double zvalue = zexpression.getVal()  ;
      // check the y-range 
      if ( zmax <= zvalue || zvalue < zmin ) { continue ; }
      //
      // fill the histogram  (only for non-zero weights and in-range entries!)
      what.fill ( xvalue , yvalue , zvalue , w ) ; 
      //
    }
    //
    return Ostap::StatusCode::SUCCESS ;
  } ;
  // ==========================================================================
  /// project data into object 
  template <class WHAT>
  Ostap::StatusCode _project3_ 
  ( const RooAbsData*                 data        , 
    const Ostap::Utils::ProgressConf& progress    ,
    WHAT&                             what        , 
    const std::string&                xexpression ,
    const std::string&                yexpression ,
    const std::string&                zexpression ,
    const std::string&                selection   ,
    const char*                       range       , 
    const double                      xmin        ,
    const double                      xmax        ,
    const double                      ymin        ,
    const double                      ymax        ,
    const double                      zmin        ,
    const double                      zmax        ,
    const unsigned long               first       ,
    const unsigned long               last        ) 
  {
    if ( 0 == data  ) { return Ostap::StatusCode ( 300 ) ; }
    //
    RooArgList        alst ;
    const RooArgSet*  aset = data->get() ;
    if ( 0 == aset       ) { return  0 ; }                            // RETURN
    Ostap::Utils::Iterator iter ( *aset );
    RooAbsArg*   coef = 0 ;
    while ( ( coef = (RooAbsArg*) iter.next() ) ){ alst.add ( *coef ); }
    //
    // convert expressions into FormulaVar 
    const bool        with_cuts = !selection.empty() ;
    const RooAbsReal* cut_var   = 0 ;
    if ( with_cuts ) { cut_var  = get_var ( *aset , selection ) ; }
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
    const RooAbsReal* z_var = get_var ( *aset , zexpression ) ;
    std::unique_ptr<Ostap::FormulaVar> zwhat ;
    if ( 0 == z_var ) 
    {
      zwhat.reset( new Ostap::FormulaVar( zexpression , alst , false ) ) ;
      if ( !zwhat->ok()   ) { return Ostap::StatusCode(305)  ; }             // RETURN
    }
    //
    return _project3_ ( data                                  , 
                        progress                              ,
                        what                                  , 
                        0 !=   x_var ?   *x_var : *xwhat      , 
                        0 !=   y_var ?   *y_var : *ywhat      , 
                        0 !=   z_var ?   *z_var : *zwhat      , 
                        0 != cut_var ?  cut_var :  cuts.get() , 
                        range ,
                        xmin  , xmax , 
                        ymin  , ymax , 
                        zmin  , zmax , 
                        first , last ) ;
  } ;
  
  // ==========================================================================
  /// project data into object 
  template <class WHAT>
  Ostap::StatusCode _project4_ 
  ( const RooAbsData*                 data        , 
    const Ostap::Utils::ProgressConf& progress    ,
    WHAT&                             what        , 
    const RooAbsReal&                 xexpression ,
    const RooAbsReal&                 yexpression ,
    const RooAbsReal&                 zexpression ,
    const RooAbsReal&                 uexpression ,
    const RooAbsReal*                 selection   ,
    const char*                       range       , 
    const double                      xmin        ,
    const double                      xmax        ,
    const double                      ymin        ,
    const double                      ymax        ,
    const double                      zmin        ,
    const double                      zmax        ,
    const double                      umin        ,
    const double                      umax        ,
    const unsigned long               first       ,
    const unsigned long               last        ) 
  {
    if ( 0 == data  ) { return Ostap::StatusCode ( 300 ) ; }
    //
    const unsigned long nEntries = 
      std::min ( last , (unsigned long) data->numEntries() ) ;
    if ( nEntries <= first  ) { return Ostap::StatusCode::RECOVERABLE ; }
    //
    const bool weighted  = data->isWeighted() ;
    const bool has_range = range && 1 <= std::strlen ( range ) ;
    //
    Ostap::Utils::ProgressBar bar ( nEntries - first , progress ) ;
    for ( unsigned long entry = first ; entry < nEntries ; ++entry , ++bar )   
    {
      //
      const RooArgSet* vars = data->get( entry ) ;
      if ( nullptr == vars )                          { break    ; } // BREAK 
      //
      if ( has_range && !vars->allInRange ( range ) ) { continue ; } // CONTINUE    
      //
      // data weight 
      const double dw = weighted  ? data      -> weight () : 1.0 ;
      if ( !dw ) { continue ; }                                      // SKIP    
      //
      // selection weight 
      const double sw = selection ? selection -> getVal () : 1.0 ;
      if ( !sw ) { continue ; }                                      // SKIP    
      //
      // calculate the total weight 
      const double w = sw * dw ;
      if ( !w  ) { continue ; }                                      // SKIP 
      //
      // calculate the x-value (only for non-zero weights)
      const double xvalue = xexpression.getVal()  ;
      // check the x-range 
      if ( xmax <= xvalue || xvalue < xmin ) { continue ; }
      // calculate the y-value (only for non-zero weights)
      const double yvalue = yexpression.getVal()  ;
      // check the y-range 
      if ( ymax <= yvalue || yvalue < ymin ) { continue ; }
      // calculate the z-value (only for non-zero weights)
      const double zvalue = zexpression.getVal()  ;
      // check the y-range 
      if ( zmax <= zvalue || zvalue < zmin ) { continue ; }
      // calculate the u-value (only for non-zero weights)
      const double uvalue = uexpression.getVal()  ;
      // check the y-range 
      if ( umax <= uvalue || uvalue < umin ) { continue ; }
      //
      // fill the histogram  (only for non-zero weights and in-range entries!)
      what.fill ( xvalue , yvalue , zvalue , uvalue , w ) ; 
      //
    }
    //
    return Ostap::StatusCode::SUCCESS ;
  } ;
  // ==========================================================================
  /// project data into object 
  template <class WHAT>
  Ostap::StatusCode _project4_ 
  ( const RooAbsData*                 data        , 
    const Ostap::Utils::ProgressConf& progress    ,
    WHAT&                             what        , 
    const std::string&                xexpression ,
    const std::string&                yexpression ,
    const std::string&                zexpression ,
    const std::string&                uexpression ,
    const std::string&                selection   ,
    const char*                       range       , 
    const double                      xmin        ,
    const double                      xmax        ,
    const double                      ymin        ,
    const double                      ymax        ,
    const double                      zmin        ,
    const double                      zmax        ,
    const double                      umin        ,
    const double                      umax        ,
    const unsigned long               first       ,
    const unsigned long               last        ) 
  {
    if ( 0 == data  ) { return Ostap::StatusCode ( 300 ) ; }
    //
    RooArgList        alst ;
    const RooArgSet*  aset = data->get() ;
    if ( 0 == aset       ) { return  0 ; }                            // RETURN
    Ostap::Utils::Iterator iter ( *aset );
    RooAbsArg*   coef = 0 ;
    while ( ( coef = (RooAbsArg*) iter.next() ) ){ alst.add ( *coef ); }
    //
    // convert expressions into FormulaVar 
    const bool        with_cuts = !selection.empty() ;
    const RooAbsReal* cut_var   = 0 ;
    if ( with_cuts ) { cut_var  = get_var ( *aset , selection ) ; }
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
    const RooAbsReal* z_var = get_var ( *aset , zexpression ) ;
    std::unique_ptr<Ostap::FormulaVar> zwhat ;
    if ( 0 == z_var ) 
    {
      zwhat.reset( new Ostap::FormulaVar( zexpression , alst , false ) ) ;
      if ( !zwhat->ok()   ) { return Ostap::StatusCode(305)  ; }             // RETURN
    }
    //
    const RooAbsReal* u_var = get_var ( *aset , uexpression ) ;
    std::unique_ptr<Ostap::FormulaVar> uwhat ;
    if ( 0 == u_var ) 
    {
      uwhat.reset( new Ostap::FormulaVar( uexpression , alst , false ) ) ;
      if ( !uwhat->ok()   ) { return Ostap::StatusCode(306)  ; }             // RETURN
    }
    //
    return _project4_ ( data                                  , 
                        progress                              ,
                        what                                  , 
                        0 !=   x_var ?   *x_var : *xwhat      , 
                        0 !=   y_var ?   *y_var : *ywhat      , 
                        0 !=   z_var ?   *z_var : *zwhat      , 
                        0 !=   u_var ?   *u_var : *uwhat      , 
                        0 != cut_var ?  cut_var :  cuts.get() , 
                        range ,
                        xmin  , xmax , 
                        ymin  , ymax , 
                        zmin  , zmax , 
                        umin  , umax , 
                        first , last ) ;
  } ;
  // ==========================================================================
  // TTree 
  // ==========================================================================
  template <class OBJECT> 
  Ostap::StatusCode 
  _project_ 
  ( TTree*                            data       , 
    const Ostap::Utils::ProgressConf& progress   ,
    OBJECT&                           obj        ,
    const std::string&                expression ,
    const std::string&                selection  ,
    const double                      xmin       , 
    const double                      xmax       ,
    const unsigned long               first      ,
    const unsigned long               last       )
  {
    if ( 0 == data          ) { return Ostap::StatusCode(300 ) ; }
    //
    const unsigned long nEntries = std::min ( last , (unsigned long) data->GetEntries() ) ;
    if ( nEntries <= first  ) { return Ostap::StatusCode::RECOVERABLE ; } 
    //
    std::unique_ptr<Ostap::Formula> cut { nullptr } ;
    if  ( !selection.empty() ) 
    { 
      cut = std::make_unique<Ostap::Formula> ( selection , data ) ; 
      if ( !cut || !cut->ok () ) { return Ostap::StatusCode ( 302 ) ; }
    }
    //
    Ostap::Formula xvar ( expression , data ) ;
    if ( !xvar || !xvar.ok() ) { return Ostap::StatusCode ( 303 ) ; }
    //
    Ostap::Utils::Notifier notify ( data , &xvar , cut.get() ) ;
    std::vector<double> results {} ;
    /// make an explicit loop 
    Ostap::Utils::ProgressBar bar (  nEntries - first , progress ) ;
    for ( unsigned long entry = first ; entry < nEntries ; ++entry , ++bar )
    {
      long ievent = data->GetEntryNumber ( entry ) ;
      if ( ievent < 0 ) { return Ostap::StatusCode( 400 ) ; }
      //
      ievent      = data->LoadTree ( ievent ) ;      
      if ( ievent < 0 ) { return Ostap::StatusCode( 401 ) ; }
      //
      const double weight = cut ? cut -> evaluate() : 1.0 ;
      if ( !weight    ) { continue ; }
      //
      xvar.evaluate ( results ) ;
      for ( double x : results ) 
      { 
        // check the x-range  & fill 
        if ( xmin <= x && x <= xmax ) { obj.fill ( x , weight  ) ; }
      }
    }
    return Ostap::StatusCode::SUCCESS ;
  }
  // ==========================================================================
  template <class OBJECT> 
  Ostap::StatusCode 
  _project2_ 
  ( TTree*                            data        , 
    const Ostap::Utils::ProgressConf& progress    ,
    OBJECT&                           obj         ,
    const std::string&                xexpression ,
    const std::string&                yexpression ,
    const std::string&                selection   ,
    const double                      xmin        , 
    const double                      xmax        ,
    const double                      ymin        , 
    const double                      ymax        ,
    const unsigned long               first       ,
    const unsigned long               last        )
  {
    if ( 0 == data          ) { return Ostap::StatusCode(300 ) ; }
    //
    const unsigned long nEntries = std::min ( last , (unsigned long) data->GetEntries() ) ;
    if ( nEntries <= first  ) { return Ostap::StatusCode::RECOVERABLE ; } 
    //
    std::unique_ptr<Ostap::Formula> cut { nullptr } ;
    if  ( !selection.empty() ) 
    { 
      cut = std::make_unique<Ostap::Formula> ( selection , data ) ; 
      if ( !cut || !cut->ok () ) { return Ostap::StatusCode ( 302 ) ; }
    }
    //
    Ostap::Formula xvar ( xexpression , data ) ;
    if ( !xvar || !xvar.ok() ) { return Ostap::StatusCode ( 303 ) ; }
    //
    Ostap::Formula yvar ( yexpression , data ) ;
    if ( !yvar || !yvar.ok() ) { return Ostap::StatusCode ( 304 ) ; }
    //
    Ostap::Utils::Notifier notify ( data , &xvar , &yvar, cut.get() ) ;
    //
    std::vector<double> xresults {} ;
    std::vector<double> yresults {} ;
    /// make an explicit loop 
    Ostap::Utils::ProgressBar bar (  nEntries - first , progress ) ;
    for ( unsigned long entry = first ; entry < nEntries ; ++entry , ++bar )
    {
      long ievent = data->GetEntryNumber ( entry ) ;
      if ( ievent < 0 ) { return Ostap::StatusCode( 400 ) ; }
      //
      ievent      = data->LoadTree ( ievent ) ;      
      if ( ievent < 0 ) { return Ostap::StatusCode( 401 ) ; }
      //
      const double weight = cut ? cut -> evaluate() : 1.0 ;
      if ( !weight    ) { continue ; }
      //
      xvar.evaluate ( xresults ) ;
      yvar.evaluate ( yresults ) ;      
      for ( double y : yresults ) 
      {
        if ( ymax < y || y < ymin ) { continue ; }
        for ( double x : xresults ) 
        { 
          if ( xmin <= x && x <= xmax ) { obj.fill ( x , y , weight  ) ; }
        }
      }
    }
    return Ostap::StatusCode::SUCCESS ;
  }
  // ==========================================================================
  template <class OBJECT> 
  Ostap::StatusCode 
  _project3_ 
  ( TTree*                            data        , 
    const Ostap::Utils::ProgressConf& progress    ,
    OBJECT&                           obj         ,
    const std::string&                xexpression ,
    const std::string&                yexpression ,
    const std::string&                zexpression ,
    const std::string&                selection   ,
    const double                      xmin        , 
    const double                      xmax        ,
    const double                      ymin        , 
    const double                      ymax        ,
    const double                      zmin        , 
    const double                      zmax        ,
    const unsigned long               first       ,
    const unsigned long               last        )
  {
    if ( 0 == data          ) { return Ostap::StatusCode(300 ) ; }
    //
    const unsigned long nEntries = std::min ( last , (unsigned long) data->GetEntries() ) ;
    if ( nEntries <= first  ) { return Ostap::StatusCode::RECOVERABLE ; } 
    //
    std::unique_ptr<Ostap::Formula> cut { nullptr } ;
    if  ( !selection.empty() ) 
    { 
      cut = std::make_unique<Ostap::Formula> ( selection , data ) ; 
      if ( !cut || !cut->ok () ) { return Ostap::StatusCode ( 302 ) ; }
    }
    //
    Ostap::Formula xvar ( xexpression , data ) ;
    if ( !xvar || !xvar.ok() ) { return Ostap::StatusCode ( 303 ) ; }
    //
    Ostap::Formula yvar ( yexpression , data ) ;
    if ( !yvar || !yvar.ok() ) { return Ostap::StatusCode ( 304 ) ; }
    //
    Ostap::Formula zvar ( zexpression , data ) ;
    if ( !zvar || !zvar.ok() ) { return Ostap::StatusCode ( 305 ) ; }
    //
    Ostap::Utils::Notifier notify ( data , &xvar , &yvar, &zvar , cut.get() ) ;
    //
    std::vector<double> xresults {} ;
    std::vector<double> yresults {} ;
    std::vector<double> zresults {} ;
    /// make an explicit loop 
    Ostap::Utils::ProgressBar bar (  nEntries - first , progress ) ;
    for ( unsigned long entry = first ; entry < nEntries ; ++entry , ++bar )
    {
      long ievent = data->GetEntryNumber ( entry ) ;
      if ( ievent < 0 ) { return Ostap::StatusCode( 400 ) ; }
      //
      ievent      = data->LoadTree ( ievent ) ;      
      if ( ievent < 0 ) { return Ostap::StatusCode( 401 ) ; }
      //
      const double weight = cut ? cut -> evaluate() : 1.0 ;
      if ( !weight    ) { continue ; }
      //
      xvar.evaluate ( xresults ) ;
      yvar.evaluate ( yresults ) ;
      zvar.evaluate ( zresults ) ;
      //
      for ( double z : zresults ) 
      {
        if ( zmax < z || z < zmin ) { continue ; }
        for ( double y : yresults ) 
        {
          if ( ymax < y || y < ymin ) { continue ; }
          for ( double x : xresults ) 
          { 
            if ( xmin <= x && x <= xmax ) { obj.fill ( x , y , z , weight  ) ; }
          }
        }
      }
    }
    return Ostap::StatusCode::SUCCESS ;
  }
  // ==========================================================================
  template <class OBJECT> 
  Ostap::StatusCode 
  _project4_ 
  ( TTree*                            data        , 
    const Ostap::Utils::ProgressConf& progress    ,
    OBJECT&                           obj         ,
    const std::string&                xexpression ,
    const std::string&                yexpression ,
    const std::string&                zexpression ,
    const std::string&                uexpression ,
    const std::string&                selection   ,
    const double                      xmin        , 
    const double                      xmax        ,
    const double                      ymin        , 
    const double                      ymax        ,
    const double                      zmin        , 
    const double                      zmax        ,
    const double                      umin        , 
    const double                      umax        ,
    const unsigned long               first       ,
    const unsigned long               last        )
  {
    if ( 0 == data          ) { return Ostap::StatusCode(300 ) ; }
    //
    const unsigned long nEntries = std::min ( last , (unsigned long) data->GetEntries() ) ;
    if ( nEntries <= first  ) { return Ostap::StatusCode::RECOVERABLE ; } 
    //
    std::unique_ptr<Ostap::Formula> cut { nullptr } ;
    if  ( !selection.empty() ) 
    { 
      cut = std::make_unique<Ostap::Formula> ( selection , data ) ; 
      if ( !cut || !cut->ok () ) { return Ostap::StatusCode ( 302 ) ; }
    }
    //
    Ostap::Formula xvar ( xexpression , data ) ;
    if ( !xvar || !xvar.ok() ) { return Ostap::StatusCode ( 303 ) ; }
    //
    Ostap::Formula yvar ( yexpression , data ) ;
    if ( !yvar || !yvar.ok() ) { return Ostap::StatusCode ( 304 ) ; }
    //
    Ostap::Formula zvar ( zexpression , data ) ;
    if ( !zvar || !zvar.ok() ) { return Ostap::StatusCode ( 305 ) ; }
    //
    Ostap::Formula uvar ( uexpression , data ) ;
    if ( !uvar || !uvar.ok() ) { return Ostap::StatusCode ( 306 ) ; }
    //
    Ostap::Utils::Notifier notify ( data , &xvar , &yvar, &zvar , &uvar, cut.get() ) ;
    //
    std::vector<double> xresults {} ;
    std::vector<double> yresults {} ;
    std::vector<double> zresults {} ;
    std::vector<double> uresults {} ;
    /// make an explicit loop 
    Ostap::Utils::ProgressBar bar (  nEntries - first , progress ) ;
    for ( unsigned long entry = first ; entry < nEntries ; ++entry , ++bar )
    {
      long ievent = data->GetEntryNumber ( entry ) ;
      if ( ievent < 0 ) { return Ostap::StatusCode( 400 ) ; }
      //
      ievent      = data->LoadTree ( ievent ) ;      
      if ( ievent < 0 ) { return Ostap::StatusCode( 401 ) ; }
      //
      const double weight = cut ? cut -> evaluate() : 1.0 ;
      if ( !weight    ) { continue ; }
      //
      xvar.evaluate ( xresults ) ;
      yvar.evaluate ( yresults ) ;
      zvar.evaluate ( zresults ) ;
      uvar.evaluate ( uresults ) ;
      //
      for ( double u : uresults ) 
      {
        if ( umax < u || u < umin ) { continue ; }
        for ( double z : zresults ) 
        {
          if ( zmax < z || z < zmin ) { continue ; }
          for ( double y : yresults ) 
          {
            if ( ymax < y || y < ymin ) { continue ; }
            for ( double x : xresults ) 
            { 
              if ( xmin <= x && x <= xmax ) { obj.fill ( x , y , z , u , weight  ) ; }
            }
          }
        }
      }
    }
    return Ostap::StatusCode::SUCCESS ;
  }
  // ==========================================================================
} // the end of anonymous namespace
// ============================================================================
/** make a projection of RooDataSet into the histogram 
 *  @param data       (INPUT)  input data
 *  @param progress   (INPUT)  configuration of progres bar 
 *  @param histo      (UPDATE) histogram 
 *  @param expression (INPUT) expression
 *  @param selection  (INPUT) selection criteria/weight 
 *  @param first      (INPUT) the first event to process 
 *  @param last       (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project
( const RooAbsData*                 data       , 
  const Ostap::Utils::ProgressConf& progress   ,
  TH1*                              histo      ,
  const RooAbsReal&                 expression ,
  const RooAbsReal*                 selection  ,
  const unsigned long               first      ,
  const unsigned long               last       ) 
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
  Ostap::Utils::ProgressBar bar ( nEntries - first , progress ) ;
  for ( unsigned long entry = first ; entry < nEntries ; ++entry , ++bar )   
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
/*  make a projection of RooDataSet into the Legendre object 
 *  @param data       (INPUT)  input data 
 *  @param object     (UPDATE) the object 
 *  @param expression (INPUT) expression
 *  @param selection  (INPUT) selection criteria/weight 
 *  @param first      (INPUT) the first event to process 
 *  @param last       (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project
( const RooAbsData*         data       , 
  Ostap::Math::LegendreSum& object     ,
  const RooAbsReal&         expression ,
  const RooAbsReal*         selection  ,
  const unsigned long       first      ,
  const unsigned long       last       ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project ( data       , 
                   progress   , 
                   object     , 
                   expression , 
                   selection  ,
                   first      ,
                   last       ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into the Legendre object 
 *  @param data       (INPUT)  input data 
 *  @param progress   (INPUT)  configuration of progres bar 
 *  @param object     (UPDATE) the object 
 *  @param expression (INPUT) expression
 *  @param selection  (INPUT) selection criteria/weight 
 *  @param first      (INPUT) the first event to process 
 *  @param last       (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project
( const RooAbsData*                 data       , 
  const Ostap::Utils::ProgressConf& progress   ,
  Ostap::Math::LegendreSum&         object     ,
  const RooAbsReal&                 expression ,
  const RooAbsReal*                 selection  ,
  const unsigned long               first      ,
  const unsigned long               last       ) 
{
  /// reset the object 
  object *= 0 ;
  return _project_ ( data           , 
                     progress       , 
                     object         ,   
                     expression     , 
                     selection      , 
                     nullptr        , 
                     object.xmin () ,
                     object.xmax () , 
                     first          ,
                     last           ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into ChebyshevSum object 
 *  @param data       (INPUT)  input data 
 *  @param object     (UPDATE) the object 
 *  @param expression (INPUT) expression
 *  @param selection  (INPUT) selection criteria/weight 
 *  @param first      (INPUT) the first event to process 
 *  @param last       (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project
( const RooAbsData*          data       , 
  Ostap::Math::ChebyshevSum& object     ,
  const RooAbsReal&          expression ,
  const RooAbsReal*          selection  ,
  const unsigned long        first      ,
  const unsigned long        last       ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project ( data       , 
                   progress   , 
                   object     , 
                   expression , 
                   selection  ,
                   first      ,
                   last       ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into ChebyshevSum  object 
 *  @param data       (INPUT)  input data 
 *  @param progress   (INPUT)  configuration of progres bar 
 *  @param object     (UPDATE) the object 
 *  @param expression (INPUT) expression
 *  @param selection  (INPUT) selection criteria/weight 
 *  @param first      (INPUT) the first event to process 
 *  @param last       (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project
( const RooAbsData*                 data       , 
  const Ostap::Utils::ProgressConf& progress   ,
  Ostap::Math::ChebyshevSum&        object     ,
  const RooAbsReal&                 expression ,
  const RooAbsReal*                 selection  ,
  const unsigned long               first      ,
  const unsigned long               last       ) 
{
  /// reset the object 
  object *= 0 ;
  return _project_ ( data           , 
                     progress       , 
                     object         ,   
                     expression     , 
                     selection      , 
                     nullptr        , 
                     object.xmin () ,
                     object.xmax () , 
                     first          ,
                     last           ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into Bernstein object 
 *  @param data       (INPUT)  input data 
 *  @param object     (UPDATE) the object 
 *  @param expression (INPUT) expression
 *  @param selection  (INPUT) selection criteria/weight 
 *  @param first      (INPUT) the first event to process 
 *  @param last       (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project
( const RooAbsData*          data       , 
  Ostap::Math::Bernstein&    object     ,
  const RooAbsReal&          expression ,
  const RooAbsReal*          selection  ,
  const unsigned long        first      ,
  const unsigned long        last       ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project ( data       , 
                   progress   , 
                   object     , 
                   expression , 
                   selection  ,
                   first      ,
                   last       ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into Bernstein object 
 *  @param data       (INPUT)  input data 
 *  @param progress   (INPUT)  configuration of progres bar 
 *  @param object     (UPDATE) the object 
 *  @param expression (INPUT) expression
 *  @param selection  (INPUT) selection criteria/weight 
 *  @param first      (INPUT) the first event to process 
 *  @param last       (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project
( const RooAbsData*                 data       , 
  const Ostap::Utils::ProgressConf& progress   ,
  Ostap::Math::Bernstein&           object     ,
  const RooAbsReal&                 expression ,
  const RooAbsReal*                 selection  ,
  const unsigned long               first      ,
  const unsigned long               last       ) 
{
  /// reset the object 
  object *= 0 ;
  return _project_ ( data           , 
                     progress       , 
                     object         ,   
                     expression     , 
                     selection      , 
                     nullptr        , 
                     object.xmin () ,
                     object.xmax () , 
                     first          ,
                     last           ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into the LegendreSum2 object 
 *  @param data       (INPUT)  input data 
 *  @param object     (UPDATE) the object 
 *  @param xexpression (INPUT) expression
 *  @param yexpression (INPUT) expression
 *  @param selection  (INPUT) selection criteria/weight 
 *  @param first      (INPUT) the first event to process 
 *  @param last       (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project2
( const RooAbsData*          data        , 
  Ostap::Math::LegendreSum2& object      ,
  const RooAbsReal&          xexpression ,
  const RooAbsReal&          yexpression ,
  const RooAbsReal*          selection   ,
  const unsigned long        first       ,
  const unsigned long        last        ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project2 ( data        , 
                    progress    , 
                    object      , 
                    xexpression , 
                    yexpression , 
                    selection   ,
                    first       ,
                    last        ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into the Legendre object 
 *  @param data       (INPUT)  input data 
 *  @param progress   (INPUT)  configuration of progres bar 
 *  @param object     (UPDATE) the object 
 *  @param expression (INPUT) expression
 *  @param selection  (INPUT) selection criteria/weight 
 *  @param first      (INPUT) the first event to process 
 *  @param last       (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project2
( const RooAbsData*                 data        , 
  const Ostap::Utils::ProgressConf& progress    ,
  Ostap::Math::LegendreSum2&        object      ,
  const RooAbsReal&                 xexpression ,
  const RooAbsReal&                 yexpression ,
  const RooAbsReal*                 selection   ,
  const unsigned long               first       ,
  const unsigned long               last        ) 
{
  /// reset the object 
  object *= 0 ;
  return _project2_ ( data           , 
                      progress       , 
                      object         ,   
                      xexpression    , 
                      yexpression    , 
                      selection      , 
                      nullptr        , 
                      object.xmin () ,
                      object.xmax () , 
                      object.ymin () ,
                      object.ymax () , 
                      first          ,
                      last           ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into the Bernstein2D  object 
 *  @param data       (INPUT)  input data 
 *  @param object     (UPDATE) the object 
 *  @param xexpression (INPUT) expression
 *  @param yexpression (INPUT) expression
 *  @param selection  (INPUT) selection criteria/weight 
 *  @param first      (INPUT) the first event to process 
 *  @param last       (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project2
( const RooAbsData*          data        , 
  Ostap::Math::Bernstein2D&  object      ,
  const RooAbsReal&          xexpression ,
  const RooAbsReal&          yexpression ,
  const RooAbsReal*          selection   ,
  const unsigned long        first       ,
  const unsigned long        last        ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project2 ( data        , 
                    progress    , 
                    object      , 
                    xexpression , 
                    yexpression , 
                    selection   ,
                    first       ,
                    last        ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into the Bernstein2D object 
 *  @param data       (INPUT)  input data 
 *  @param progress   (INPUT)  configuration of progres bar 
 *  @param object     (UPDATE) the object 
 *  @param expression (INPUT) expression
 *  @param selection  (INPUT) selection criteria/weight 
 *  @param first      (INPUT) the first event to process 
 *  @param last       (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project2
( const RooAbsData*                 data        , 
  const Ostap::Utils::ProgressConf& progress    ,
  Ostap::Math::Bernstein2D&  object      ,
  const RooAbsReal&                 xexpression ,
  const RooAbsReal&                 yexpression ,
  const RooAbsReal*                 selection   ,
  const unsigned long               first       ,
  const unsigned long               last        ) 
{
  /// reset the object 
  object *= 0 ;
  return _project2_ ( data           , 
                      progress       , 
                      object         ,   
                      xexpression    , 
                      yexpression    , 
                      selection      , 
                      nullptr        , 
                      object.xmin () ,
                      object.xmax () , 
                      object.ymin () ,
                      object.ymax () , 
                      first          ,
                      last           ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into the LegendreSum3 object 
 *  @param data       (INPUT)  input data 
 *  @param object     (UPDATE) the object 
 *  @param xexpression (INPUT) expression
 *  @param yexpression (INPUT) expression
 *  @param zexpression (INPUT) zxpression
 *  @param selection  (INPUT) selection criteria/weight 
 *  @param first      (INPUT) the first event to process 
 *  @param last       (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project3
( const RooAbsData*          data        , 
  Ostap::Math::LegendreSum3& object      ,
  const RooAbsReal&          xexpression ,
  const RooAbsReal&          yexpression ,
  const RooAbsReal&          zexpression ,
  const RooAbsReal*          selection   ,
  const unsigned long        first       ,
  const unsigned long        last        ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project3 ( data        , 
                    progress    , 
                    object      , 
                    xexpression , 
                    yexpression , 
                    zexpression , 
                    selection   ,
                    first       ,
                    last        ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into the LegendreSum3 object 
 *  @param data       (INPUT)  input data 
 *  @param progress   (INPUT)  configuration of progres bar 
 *  @param object     (UPDATE) the object 
 *  @param xexpression (INPUT) expression
 *  @param yexpression (INPUT) expression
 *  @param zexpression (INPUT) zxpression
 *  @param selection  (INPUT) selection criteria/weight 
 *  @param first      (INPUT) the first event to process 
 *  @param last       (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project3
( const RooAbsData*                 data        , 
  const Ostap::Utils::ProgressConf& progress    ,
  Ostap::Math::LegendreSum3&        object      ,
  const RooAbsReal&                 xexpression ,
  const RooAbsReal&                 yexpression ,
  const RooAbsReal&                 zexpression ,
  const RooAbsReal*                 selection   ,
  const unsigned long               first       ,
  const unsigned long               last        ) 
{
  /// reset the object 
  object *= 0 ;
  return _project3_ ( data           , 
                      progress       , 
                      object         ,   
                      xexpression    , 
                      yexpression    , 
                      zexpression    , 
                      selection      , 
                      nullptr        , 
                      object.xmin () ,
                      object.xmax () , 
                      object.ymin () ,
                      object.ymax () , 
                      object.zmin () ,
                      object.zmax () , 
                      first          ,
                      last           ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into the Bernstein3D  object 
 *  @param data       (INPUT)  input data 
 *  @param object     (UPDATE) the object 
 *  @param xexpression (INPUT) x-expression
 *  @param yexpression (INPUT) y-expression
 *  @param zexpression (INPUT) z-expression
 *  @param selection  (INPUT) selection criteria/weight 
 *  @param first      (INPUT) the first event to process 
 *  @param last       (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project3
( const RooAbsData*          data        , 
  Ostap::Math::Bernstein3D&  object      ,
  const RooAbsReal&          xexpression ,
  const RooAbsReal&          yexpression ,
  const RooAbsReal&          zexpression ,
  const RooAbsReal*          selection   ,
  const unsigned long        first       ,
  const unsigned long        last        ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project3 ( data        , 
                    progress    , 
                    object      , 
                    xexpression , 
                    yexpression , 
                    zexpression , 
                    selection   ,
                    first       ,
                    last        ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into the Bernstein3D object 
 *  @param data       (INPUT)  input data 
 *  @param progress   (INPUT)  configuration of progres bar 
 *  @param object     (UPDATE) the object 
 *  @param xexpression (INPUT) x-expression
 *  @param yexpression (INPUT) y-expression
 *  @param zexpression (INPUT) z-expression
 *  @param selection  (INPUT) selection criteria/weight 
 *  @param first      (INPUT) the first event to process 
 *  @param last       (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project3
( const RooAbsData*                 data        , 
  const Ostap::Utils::ProgressConf& progress    ,
  Ostap::Math::Bernstein3D&  object      ,
  const RooAbsReal&                 xexpression ,
  const RooAbsReal&                 yexpression ,
  const RooAbsReal&                 zexpression ,
  const RooAbsReal*                 selection   ,
  const unsigned long               first       ,
  const unsigned long               last        ) 
{
  /// reset the object 
  object *= 0 ;
  return _project3_ ( data           , 
                      progress       , 
                      object         ,   
                      xexpression    , 
                      yexpression    , 
                      zexpression    , 
                      selection      , 
                      nullptr        , 
                      object.xmin () ,
                      object.xmax () , 
                      object.ymin () ,
                      object.ymax () , 
                      object.zmin () ,
                      object.zmax () , 
                      first          ,
                      last           ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into the LegendreSum4 object 
 *  @param data       (INPUT)  input data 
 *  @param object     (UPDATE) the object 
 *  @param xexpression (INPUT) xexpression
 *  @param yexpression (INPUT) yexpression
 *  @param zexpression (INPUT) zxpression
 *  @param uexpression (INPUT) uxpression
 *  @param selection  (INPUT) selection criteria/weight 
 *  @param first      (INPUT) the first event to process 
 *  @param last       (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project4
( const RooAbsData*          data        , 
  Ostap::Math::LegendreSum4& object      ,
  const RooAbsReal&          xexpression ,
  const RooAbsReal&          yexpression ,
  const RooAbsReal&          zexpression ,
  const RooAbsReal&          uexpression ,
  const RooAbsReal*          selection   ,
  const unsigned long        first       ,
  const unsigned long        last        ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project4 ( data        , 
                    progress    , 
                    object      , 
                    xexpression , 
                    yexpression , 
                    zexpression , 
                    uexpression , 
                    selection   ,
                    first       ,
                    last        ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into the LegendreSum3 object 
 *  @param data       (INPUT)  input data 
 *  @param progress   (INPUT)  configuration of progres bar 
 *  @param object     (UPDATE) the object 
 *  @param xexpression (INPUT) expression
 *  @param yexpression (INPUT) expression
 *  @param zexpression (INPUT) zxpression
 *  @param selection  (INPUT) selection criteria/weight 
 *  @param first      (INPUT) the first event to process 
 *  @param last       (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project4
( const RooAbsData*                 data        , 
  const Ostap::Utils::ProgressConf& progress    ,
  Ostap::Math::LegendreSum4&        object      ,
  const RooAbsReal&                 xexpression ,
  const RooAbsReal&                 yexpression ,
  const RooAbsReal&                 zexpression ,
  const RooAbsReal&                 uexpression ,
  const RooAbsReal*                 selection   ,
  const unsigned long               first       ,
  const unsigned long               last        ) 
{
  /// reset the object 
  object *= 0 ;
  return _project4_ ( data           , 
                      progress       , 
                      object         ,   
                      xexpression    , 
                      yexpression    , 
                      zexpression    , 
                      uexpression    , 
                      selection      , 
                      nullptr        , 
                      object.xmin () ,
                      object.xmax () , 
                      object.ymin () ,
                      object.ymax () , 
                      object.zmin () ,
                      object.zmax () , 
                      object.umin () ,
                      object.umax () , 
                      first          ,
                      last           ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into LegendreSum object 
 *  @param data       (INPUT)  input data 
 *  @param object     (UPDATE) the object 
 *  @param expression (INPUT) expression
 *  @param selection  (INPUT) selection criteria/weight 
 *  @param first      (INPUT) the first event to process 
 *  @param last       (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project
( const RooAbsData*         data       , 
  Ostap::Math::LegendreSum& object     ,
  const std::string&        expression , 
  const std::string&        selection  , 
  const unsigned long       first      ,
  const unsigned long       last       ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project ( data        , 
                   progress    , 
                   object      , 
                   expression  , 
                   selection   ,
                   first       ,
                   last        ) ;
}
// ========================================================================
/*  make a projection of RooDataSet into LegendreSum object 
 *  @param data       (INPUT)  input data 
 *  @param progress   (INPUT)  configuration of progres bar 
 *  @param object     (UPDATE) the object 
 *  @param expression (INPUT) expression
 *  @param selection  (INPUT) selection criteria/weight 
 *  @param first      (INPUT) the first event to process 
 *  @param last       (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project
( const RooAbsData*                 data       , 
  const Ostap::Utils::ProgressConf& progress   ,
  Ostap::Math::LegendreSum&         object     ,
  const std::string&                expression , 
  const std::string&                selection  , 
  const unsigned long               first      ,
  const unsigned long               last       ) 
{
  /// reset the object 
  object *= 0 ;
  return _project_ ( data           , 
                     progress       , 
                     object         ,   
                     expression     , 
                     selection      , 
                     nullptr        , 
                     object.xmin () ,
                     object.xmax () , 
                     first          ,
                     last           ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into LegendreSum object 
 *  @param data       (INPUT)  input data 
 *  @param object     (UPDATE) the object 
 *  @param expression (INPUT) expression
 *  @param selection  (INPUT) selection criteria/weight 
 *  @param first      (INPUT) the first event to process 
 *  @param last       (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project
( const RooAbsData*          data       , 
  Ostap::Math::ChebyshevSum& object     ,
  const std::string&         expression , 
  const std::string&         selection  , 
  const unsigned long        first      ,
  const unsigned long        last       ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project ( data        , 
                   progress    , 
                   object      , 
                   expression  , 
                   selection   ,
                   first       ,
                   last        ) ;
}
// ========================================================================
/*  make a projection of RooDataSet into LegendreSum object 
 *  @param data       (INPUT)  input data 
 *  @param progress   (INPUT)  configuration of progres bar 
 *  @param object     (UPDATE) the object 
 *  @param expression (INPUT) expression
 *  @param selection  (INPUT) selection criteria/weight 
 *  @param first      (INPUT) the first event to process 
 *  @param last       (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project
( const RooAbsData*                 data       , 
  const Ostap::Utils::ProgressConf& progress   ,
  Ostap::Math::ChebyshevSum&        object     ,
  const std::string&                expression , 
  const std::string&                selection  , 
  const unsigned long               first      ,
  const unsigned long               last       ) 
{
  /// reset the object 
  object *= 0 ;
  return _project_ ( data           , 
                     progress       , 
                     object         ,   
                     expression     , 
                     selection      , 
                     nullptr        , 
                     object.xmin () ,
                     object.xmax () , 
                     first          ,
                     last           ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into LegendreSum object 
 *  @param data       (INPUT)  input data 
 *  @param object     (UPDATE) the object 
 *  @param expression (INPUT) expression
 *  @param selection  (INPUT) selection criteria/weight 
 *  @param first      (INPUT) the first event to process 
 *  @param last       (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project
( const RooAbsData*         data       , 
  Ostap::Math::Bernstein&   object     ,
  const std::string&        expression , 
  const std::string&        selection  , 
  const unsigned long       first      ,
  const unsigned long       last       ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project ( data        , 
                   progress    , 
                   object      , 
                   expression  , 
                   selection   ,
                   first       ,
                   last        ) ;
}
// ========================================================================
/*  make a projection of RooDataSet into LegendreSum object 
 *  @param data       (INPUT)  input data 
 *  @param progress   (INPUT)  configuration of progres bar 
 *  @param object     (UPDATE) the object 
 *  @param expression (INPUT) expression
 *  @param selection  (INPUT) selection criteria/weight 
 *  @param first      (INPUT) the first event to process 
 *  @param last       (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project
( const RooAbsData*                 data       , 
  const Ostap::Utils::ProgressConf& progress   ,
  Ostap::Math::Bernstein&           object     ,
  const std::string&                expression , 
  const std::string&                selection  , 
  const unsigned long               first      ,
  const unsigned long               last       ) 
{
  /// reset the object 
  object *= 0 ;
  return _project_ ( data           , 
                     progress       , 
                     object         ,   
                     expression     , 
                     selection      , 
                     nullptr        , 
                     object.xmin () ,
                     object.xmax () , 
                     first          ,
                     last           ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into LegendreSum2 object 
 *  @param data        (INPUT)  input data 
 *  @param object      (UPDATE) the object 
 *  @param xexpression (INPUT) x-expression
 *  @param yexpression (INPUT) y-expression
 *  @param selection   (INPUT) selection criteria/weight 
 *  @param first       (INPUT) the first event to process 
 *  @param last        (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project2
( const RooAbsData*          data        , 
  Ostap::Math::LegendreSum2& object      ,
  const std::string&         xexpression , 
  const std::string&         yexpression , 
  const std::string&         selection  , 
  const unsigned long        first      ,
  const unsigned long        last       ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project2 ( data        , 
                    progress    , 
                    object      , 
                    xexpression ,
                    yexpression , 
                    selection   ,
                    first       ,
                    last        ) ;
}
// ========================================================================
/*  make a projection of RooDataSet into LegendreSum2 object 
 *  @param data        (INPUT)  input data 
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param object      (UPDATE) the object 
 *  @param xexpression (INPUT) x-expression
 *  @param yexpression (INPUT) y-expression
 *  @param selection   (INPUT) selection criteria/weight 
 *  @param first       (INPUT) the first event to process 
 *  @param last        (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project2
( const RooAbsData*                 data        , 
  const Ostap::Utils::ProgressConf& progress    ,
  Ostap::Math::LegendreSum2&        object      ,
  const std::string&                xexpression , 
  const std::string&                yexpression , 
  const std::string&                selection   , 
  const unsigned long               first       ,
  const unsigned long               last        ) 
{
  /// reset the object 
  object *= 0 ;
  return _project2_ ( data           , 
                      progress       , 
                      object         ,   
                      xexpression    , 
                      yexpression    , 
                      selection      , 
                      nullptr        , 
                      object.xmin () ,
                      object.xmax () , 
                      object.ymin () ,
                      object.ymax () , 
                      first          ,
                      last           ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into Bernstein2D object 
 *  @param data        (INPUT)  input data 
 *  @param object      (UPDATE) the object 
 *  @param xexpression (INPUT) x-expression
 *  @param yexpression (INPUT) y-expression
 *  @param selection   (INPUT) selection criteria/weight 
 *  @param first       (INPUT) the first event to process 
 *  @param last        (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project2
( const RooAbsData*          data        , 
  Ostap::Math::Bernstein2D&  object      ,
  const std::string&         xexpression , 
  const std::string&         yexpression , 
  const std::string&         selection  , 
  const unsigned long        first      ,
  const unsigned long        last       ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project2 ( data        , 
                    progress    , 
                    object      , 
                    xexpression ,
                    yexpression , 
                    selection   ,
                    first       ,
                    last        ) ;
}
// ========================================================================
/*  make a projection of RooDataSet into Bernstein2D object 
 *  @param data        (INPUT)  input data 
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param object      (UPDATE) the object 
 *  @param xexpression (INPUT) x-expression
 *  @param yexpression (INPUT) y-expression
 *  @param selection   (INPUT) selection criteria/weight 
 *  @param first       (INPUT) the first event to process 
 *  @param last        (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project2
( const RooAbsData*                 data        , 
  const Ostap::Utils::ProgressConf& progress    ,
  Ostap::Math::Bernstein2D&         object      ,
  const std::string&                xexpression , 
  const std::string&                yexpression , 
  const std::string&                selection   , 
  const unsigned long               first       ,
  const unsigned long               last        ) 
{
  /// reset the object 
  object *= 0 ;
  return _project2_ ( data           , 
                      progress       , 
                      object         ,   
                      xexpression    , 
                      yexpression    , 
                      selection      , 
                      nullptr        , 
                      object.xmin () ,
                      object.xmax () , 
                      object.ymin () ,
                      object.ymax () , 
                      first          ,
                      last           ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into LegendreSum3 object 
 *  @param data        (INPUT)  input data 
 *  @param object      (UPDATE) the object 
 *  @param xexpression (INPUT) x-expression
 *  @param yexpression (INPUT) y-expression
 *  @param zexpression (INPUT) z-expression
 *  @param selection   (INPUT) selection criteria/weight 
 *  @param first       (INPUT) the first event to process 
 *  @param last        (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project3
( const RooAbsData*          data        , 
  Ostap::Math::LegendreSum3& object      ,
  const std::string&         xexpression , 
  const std::string&         yexpression , 
  const std::string&         zexpression , 
  const std::string&         selection  , 
  const unsigned long        first      ,
  const unsigned long        last       ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project3 ( data        , 
                    progress    , 
                    object      , 
                    xexpression ,
                    yexpression , 
                    zexpression , 
                    selection   ,
                    first       ,
                    last        ) ;
}
// ========================================================================
/*  make a projection of RooDataSet into LegendreSum3 object 
 *  @param data        (INPUT)  input data 
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param object      (UPDATE) the object 
 *  @param xexpression (INPUT) x-expression
 *  @param yexpression (INPUT) y-expression
 *  @param zexpression (INPUT) z-expression
 *  @param selection   (INPUT) selection criteria/weight 
 *  @param first       (INPUT) the first event to process 
 *  @param last        (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project3
( const RooAbsData*                 data        , 
  const Ostap::Utils::ProgressConf& progress    ,
  Ostap::Math::LegendreSum3&        object      ,
  const std::string&                xexpression , 
  const std::string&                yexpression , 
  const std::string&                zexpression , 
  const std::string&                selection   , 
  const unsigned long               first       ,
  const unsigned long               last        ) 
{
  /// reset the object 
  object *= 0 ;
  return _project3_ ( data           , 
                      progress       , 
                      object         ,   
                      xexpression    , 
                      yexpression    , 
                      zexpression    , 
                      selection      , 
                      nullptr        , 
                      object.xmin () ,
                      object.xmax () , 
                      object.ymin () ,
                      object.ymax () , 
                      object.zmin () ,
                      object.zmax () , 
                      first          ,
                      last           ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into LegendreSum3 object 
 *  @param data        (INPUT)  input data 
 *  @param object      (UPDATE) the object 
 *  @param xexpression (INPUT) x-expression
 *  @param yexpression (INPUT) y-expression
 *  @param zexpression (INPUT) z-expression
 *  @param selection   (INPUT) selection criteria/weight 
 *  @param first       (INPUT) the first event to process 
 *  @param last        (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project3
( const RooAbsData*          data        , 
  Ostap::Math::Bernstein3D&  object      ,
  const std::string&         xexpression , 
  const std::string&         yexpression , 
  const std::string&         zexpression , 
  const std::string&         selection  , 
  const unsigned long        first      ,
  const unsigned long        last       ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project3 ( data        , 
                    progress    , 
                    object      , 
                    xexpression ,
                    yexpression , 
                    zexpression , 
                    selection   ,
                    first       ,
                    last        ) ;
}
// ========================================================================
/*  make a projection of RooDataSet into LegendreSum3 object 
 *  @param data        (INPUT)  input data 
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param object      (UPDATE) the object 
 *  @param xexpression (INPUT) x-expression
 *  @param yexpression (INPUT) y-expression
 *  @param zexpression (INPUT) z-expression
 *  @param selection   (INPUT) selection criteria/weight 
 *  @param first       (INPUT) the first event to process 
 *  @param last        (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project3
( const RooAbsData*                 data        , 
  const Ostap::Utils::ProgressConf& progress    ,
  Ostap::Math::Bernstein3D&         object      ,
  const std::string&                xexpression , 
  const std::string&                yexpression , 
  const std::string&                zexpression , 
  const std::string&                selection   , 
  const unsigned long               first       ,
  const unsigned long               last        ) 
{
  /// reset the object 
  object *= 0 ;
  return _project3_ ( data           , 
                      progress       , 
                      object         ,   
                      xexpression    , 
                      yexpression    , 
                      zexpression    , 
                      selection      , 
                      nullptr        , 
                      object.xmin () ,
                      object.xmax () , 
                      object.ymin () ,
                      object.ymax () , 
                      object.zmin () ,
                      object.zmax () , 
                      first          ,
                      last           ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into LegendreSum4 object 
 *  @param data        (INPUT)  input data 
 *  @param object      (UPDATE) the object 
 *  @param xexpression (INPUT) x-expression
 *  @param yexpression (INPUT) y-expression
 *  @param zexpression (INPUT) z-expression
 *  @param uexpression (INPUT) u-expression
 *  @param selection   (INPUT) selection criteria/weight 
 *  @param first       (INPUT) the first event to process 
 *  @param last        (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project4
( const RooAbsData*          data        , 
  Ostap::Math::LegendreSum4& object      ,
  const std::string&         xexpression , 
  const std::string&         yexpression , 
  const std::string&         zexpression , 
  const std::string&         uexpression , 
  const std::string&         selection   , 
  const unsigned long        first       ,
  const unsigned long        last        ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project4 ( data        , 
                    progress    , 
                    object      , 
                    xexpression ,
                    yexpression , 
                    zexpression , 
                    uexpression , 
                    selection   ,
                    first       ,
                    last        ) ;
}
// ========================================================================
/*  make a projection of RooDataSet into LegendreSum4 object 
 *  @param data        (INPUT)  input data 
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param object      (UPDATE) the object 
 *  @param xexpression (INPUT) x-expression
 *  @param yexpression (INPUT) y-expression
 *  @param zexpression (INPUT) z-expression
 *  @param uexpression (INPUT) u-expression
 *  @param selection   (INPUT) selection criteria/weight 
 *  @param first       (INPUT) the first event to process 
 *  @param last        (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project4
( const RooAbsData*                 data        , 
  const Ostap::Utils::ProgressConf& progress    ,
  Ostap::Math::LegendreSum4&        object      ,
  const std::string&                xexpression , 
  const std::string&                yexpression , 
  const std::string&                zexpression , 
  const std::string&                uexpression , 
  const std::string&                selection   , 
  const unsigned long               first       ,
  const unsigned long               last        ) 
{
  /// reset the object 
  object *= 0 ;
  return _project4_ ( data           , 
                      progress       , 
                      object         ,   
                      xexpression    , 
                      yexpression    , 
                      zexpression    , 
                      uexpression    , 
                      selection      , 
                      nullptr        , 
                      object.xmin () ,
                      object.xmax () , 
                      object.ymin () ,
                      object.ymax () , 
                      object.zmin () ,
                      object.zmax () , 
                      object.umin () ,
                      object.umax () , 
                      first          ,
                      last           ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into the histogram 
 *  @param data       (INPUT)  input data 
 *  @param histo      (UPDATE) histogram 
 *  @param expression (INPUT) expression
 *  @param selection  (INPUT) selection criteria/weight 
 *  @param first      (INPUT) the first event to process 
 *  @param last       (INPUT) the last event to process 
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
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project ( data       , 
                   progress   , 
                   histo      , 
                   expression , 
                   selection  , 
                   first      ,
                   last       ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into the histogram 
 *  @param data        (INPUT)  input data 
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param histo       (UPDATE) histogram 
 *  @param xexpression (INPUT) expression for x-axis 
 *  @param yexpression (INPUT) expression for y-axis 
 *  @param selection   (INPUT) selection criteria/weight 
 *  @param first       (INPUT) the first event to process 
 *  @param last        (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode
Ostap::HistoProject::project2
( const RooAbsData*                 data        , 
  const Ostap::Utils::ProgressConf& progress    ,
  TH2*                              histo       ,
  const RooAbsReal&                 xexpression ,
  const RooAbsReal&                 yexpression ,
  const RooAbsReal*                 selection   ,
  const unsigned long               first       ,
  const unsigned long               last        ) 
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
  Ostap::Utils::ProgressBar bar ( nEntries - first , progress ) ;
  for ( unsigned long entry = first ; entry < nEntries ; ++entry , ++bar )   
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
 *  @param data        (INPUT)  input data 
 *  @param histo       (UPDATE) histogram 
 *  @param xexpression (INPUT) expression for x-axis 
 *  @param yexpression (INPUT) expression for y-axis 
 *  @param selection   (INPUT) selection criteria/weight 
 *  @param first       (INPUT) the first event to process 
 *  @param last        (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode
Ostap::HistoProject::project2
( const RooAbsData*                 data        , 
  TH2*                              histo       ,
  const RooAbsReal&                 xexpression ,
  const RooAbsReal&                 yexpression ,
  const RooAbsReal*                 selection   ,
  const unsigned long               first       ,
  const unsigned long               last        ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project2 ( data        , 
                    progress    , 
                    histo       , 
                    xexpression , 
                    yexpression , 
                    selection   , 
                    first       ,
                    last        ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into the histogram 
 *  @param data        (INPUT)  input data 
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param histo       (UPDATE) histogram 
 *  @param xexpression (INPUT) expression for x-axis 
 *  @param yexpression (INPUT) expression for y-axis 
 *  @param zexpression (INPUT) expression for z-axis 
 *  @param selection   (INPUT) selection criteria/weight 
 *  @param first       (INPUT) the first event to process 
 *  @param last        (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode
Ostap::HistoProject::project3
( const RooAbsData*                 data        , 
  const Ostap::Utils::ProgressConf& progress    ,
  TH3*                              histo       ,
  const RooAbsReal&                 xexpression ,
  const RooAbsReal&                 yexpression ,
  const RooAbsReal&                 zexpression ,
  const RooAbsReal*                 selection   ,
  const unsigned long               first       ,
  const unsigned long               last        ) 
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
  Ostap::Utils::ProgressBar bar ( nEntries - first , progress ) ;
  for ( unsigned long entry = first ; entry < nEntries ; ++entry , ++bar )   
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
/*  make a projection of RooDataSet into the histogram 
 *  @param data        (INPUT)  input data 
 *  @param histo       (UPDATE) histogram 
 *  @param xexpression (INPUT) expression for x-axis 
 *  @param yexpression (INPUT) expression for y-axis 
 *  @param zexpression (INPUT) expression for z-axis 
 *  @param selection   (INPUT) selection criteria/weight 
 *  @param first       (INPUT) the first event to process 
 *  @param last        (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode
Ostap::HistoProject::project3
( const RooAbsData*                 data        , 
  TH3*                              histo       ,
  const RooAbsReal&                 xexpression ,
  const RooAbsReal&                 yexpression ,
  const RooAbsReal&                 zexpression ,
  const RooAbsReal*                 selection   ,
  const unsigned long               first       ,
  const unsigned long               last        ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project3 ( data        , 
                    progress    , 
                    histo       , 
                    xexpression , 
                    yexpression , 
                    zexpression , 
                    selection   , 
                    first       ,
                    last        ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into the histogram 
 *  @param data        (INPUT)  input data 
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param histo       (UPDATE) histogram 
 *  @param expression  (INPUT)  expression
 *  @param selection   (INPUT)  selection criteria/weight 
 *  @param first       (INPUT)  the first event to process 
 *  @param last        (INPUT)  the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project
( const RooAbsData*                 data       , 
  const Ostap::Utils::ProgressConf& progress   ,
  TH1*                              histo      ,
  const std::string&                expression ,
  const std::string&                selection  ,
  const unsigned long               first      ,                                          
  const unsigned long               last       ) 
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
  if ( 0 == aset       ) { return  0 ; }                            // RETURN
  Ostap::Utils::Iterator iter ( *aset );
  RooAbsArg*   coef = 0 ;
  while ( ( coef = (RooAbsArg*) iter.next() ) ){ alst.add ( *coef ); }
  //
  // convert expressions into FormulaVar 
  const bool        with_cuts = !selection.empty() ;
  const RooAbsReal* cut_var   = 0 ;
  if ( with_cuts ) { cut_var  = get_var ( *aset , selection ) ; }
  std::unique_ptr<Ostap::FormulaVar> cuts ;
  if ( with_cuts && 0 == cut_var ) 
  {
    cuts.reset ( new Ostap::FormulaVar ( selection , alst , false ) ) ;
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
  return project ( data                                  , 
                   progress                              ,
                   histo                                 , 
                   0 !=   x_var ?   *x_var : *xwhat      , 
                   0 != cut_var ?  cut_var :  cuts.get() , first , last ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into the histogram 
 *  @param data        (INPUT)  input data 
 *  @param histo       (UPDATE) histogram 
 *  @param expression  (INPUT)  expression
 *  @param selection   (INPUT)  selection criteria/weight 
 *  @param first       (INPUT)  the first event to process 
 *  @param last        (INPUT)  the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project
( const RooAbsData*                 data       , 
  TH1*                              histo      ,
  const std::string&                expression ,
  const std::string&                selection  ,
  const unsigned long               first      ,                                          
  const unsigned long               last       ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project ( data        , 
                   progress    , 
                   histo       , 
                   expression  , 
                   selection   , 
                   first       ,
                   last        ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into the histogram 
 *  @param data        (INPUT)  input data 
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param histo       (UPDATE) histogram 
 *  @param xexpression (INPUT) expression for x-axis 
 *  @param yexpression (INPUT) expression for y-axis 
 *  @param selection   (INPUT) selection criteria/weight 
 *  @param first       (INPUT) the first event to process 
 *  @param last        (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode
Ostap::HistoProject::project2
( const RooAbsData*                 data        , 
  const Ostap::Utils::ProgressConf& progress    ,
  TH2*                              histo       ,
  const std::string&                xexpression ,
  const std::string&                yexpression ,
  const std::string&                selection   ,
  const unsigned long               first       ,
  const unsigned long               last        ) 
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
                    progress                               ,
                    histo                                  , 
                    0 !=   x_var ?   *x_var : *xwhat       , 
                    0 !=   y_var ?   *y_var : *ywhat       , 
                    0 != cut_var ?  cut_var :   cuts.get() , first , last ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into the histogram 
 *  @param data        (INPUT)  input data 
 *  @param histo       (UPDATE) histogram 
 *  @param xexpression (INPUT) expression for x-axis 
 *  @param yexpression (INPUT) expression for y-axis 
 *  @param selection   (INPUT) selection criteria/weight 
 *  @param first       (INPUT) the first event to process 
 *  @param last        (INPUT) the last event to process 
 */
// ============================================================================
Ostap::StatusCode
Ostap::HistoProject::project2
( const RooAbsData*                 data        , 
  TH2*                              histo       ,
  const std::string&                xexpression ,
  const std::string&                yexpression ,
  const std::string&                selection   ,
  const unsigned long               first       ,
  const unsigned long               last        ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project2 ( data        , 
                    progress    , 
                    histo       , 
                    xexpression , 
                    yexpression , 
                    selection   , 
                    first       ,
                    last        ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into the histogram 
 *  @param data        (INPUT)  input data 
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param histo       (UPDATE) histogram 
 *  @param xexpression (INPUT)  expression for x-axis 
 *  @param yexpression (INPUT)  expression for y-axis 
 *  @param zexpression (INPUT)  expression for z-axis 
 *  @param selection   (INPUT)  selection criteria/weight 
 *  @param first       (INPUT)  the first event to process 
 *  @param last        (INPUT)  the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project3
( const RooAbsData*                 data        , 
  const Ostap::Utils::ProgressConf& progress    ,
  TH3*                              histo       ,
  const std::string&                xexpression ,
  const std::string&                yexpression ,
  const std::string&                zexpression ,
  const std::string&                selection   ,
  const unsigned long               first       ,
  const unsigned long               last        ) 
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
                    progress                               , 
                    histo                                  , 
                    0 !=   x_var ?   *x_var : *xwhat       , 
                    0 !=   y_var ?   *y_var : *ywhat       , 
                    0 !=   z_var ?   *z_var : *zwhat       , 
                    0 != cut_var ?  cut_var :   cuts.get() , first , last ) ;
}
// ============================================================================
/*  make a projection of RooDataSet into the histogram 
 *  @param data        (INPUT)  input data 
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param histo       (UPDATE) histogram 
 *  @param xexpression (INPUT)  expression for x-axis 
 *  @param yexpression (INPUT)  expression for y-axis 
 *  @param zexpression (INPUT)  expression for z-axis 
 *  @param selection   (INPUT)  selection criteria/weight 
 *  @param first       (INPUT)  the first event to process 
 *  @param last        (INPUT)  the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project3
( const RooAbsData*                 data        , 
  TH3*                              histo       ,
  const std::string&                xexpression ,
  const std::string&                yexpression ,
  const std::string&                zexpression ,
  const std::string&                selection   ,
  const unsigned long               first       ,
  const unsigned long               last        ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project3 ( data        , 
                    progress    , 
                    histo       , 
                    xexpression , 
                    yexpression , 
                    zexpression , 
                    selection   , 
                    first       ,
                    last        ) ;
}
// ============================================================================
/*  make a projection of DataFrame into the histogram 
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param data        (INPUT)  input data 
 *  @param histo       (UPDATE) histogram 
 *  @param expression  (INPUT)  expression
 *  @param selection   (INPUT)  selection criteria/weight 
 */
// ============================================================================
Ostap::StatusCode Ostap::HistoProject::project
( Ostap::FrameNode                     data        , 
  const Ostap::Utils::ProgressConf& /* progress */ ,
  TH1*                                 histo       ,
  const std::string&                   expression  ,
  const std::string&                   selection   ) 
{ return project ( data , histo , expression , selection ) ; }
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
#if   ROOT_VERSION(6,12,0) <= ROOT_VERSION_CODE
  //
  TH1D model {} ; histo->Copy ( model ) ;
  //
#else
  //
  TH1F model {} ; histo->Copy ( model ) ;
  //  
#endif
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
 *  @param data        (INPUT)  input data 
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param histo       (UPDATE) histogram 
 *  @param xexpression (INPUT)  expression for x-axis 
 *  @param yexpression (INPUT)  expression for y-axis 
 *  @param selection   (INPUT)  selection criteria/weight 
 */
// ========================================================================
Ostap::StatusCode Ostap::HistoProject::project2
( Ostap::FrameNode                     data        , 
  const Ostap::Utils::ProgressConf& /* progress */ ,
  TH2*                                 histo       ,
  const std::string&                   xexpression ,
  const std::string&                   yexpression ,
  const std::string&                   selection   )
{ return project2 ( data , histo , xexpression , yexpression , selection ) ; }
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
#if   ROOT_VERSION(6,12,0) <= ROOT_VERSION_CODE
  //
  TH2D model {} ; histo->Copy ( model ) ;
  //
#else
  //
  TH2F model {} ; histo->Copy ( model ) ;
  //  
#endif
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
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param histo       (UPDATE) histogram 
 *  @param xexpression (INPUT)  expression for x-axis 
 *  @param yexpression (INPUT)  expression for y-axis 
 *  @param zexpression (INPUT)  expression for z-axis 
 *  @param selection   (INPUT)  selection criteria/weight 
 */
// ========================================================================
Ostap::StatusCode Ostap::HistoProject::project3
( Ostap::FrameNode                     data        , 
  const Ostap::Utils::ProgressConf& /* progress */ ,
  TH3*                                 histo       ,
  const std::string&                   xexpression ,
  const std::string&                   yexpression ,
  const std::string&                   zexpression ,
  const std::string&                   selection   )
{ return project3 ( data , histo , xexpression , yexpression , zexpression , selection ) ; }
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
#if   ROOT_VERSION(6,12,0) <= ROOT_VERSION_CODE
  //
  TH3D model {} ; histo->Copy ( model ) ;
  //
#else
  //
  TH3F model {} ; histo->Copy ( model ) ;
  //  
#endif
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


// ============================================================================
//  TTree
// ============================================================================

// ============================================================================
/*  make a projection of TTree into the histogram 
 *  @param data        (INPUT)  input data 
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param histo       (UPDATE) histogram 
 *  @param expression  (INPUT)  expression
 *  @param selection   (INPUT)  selection criteria/weight 
 *  @param first       (INPUT)  the first event to process 
 *  @param last        (INPUT)  the last event to process 
 */
// ============================================================================
Ostap::StatusCode Ostap::HistoProject::project
( TTree*                            data       , 
  const Ostap::Utils::ProgressConf& progress   ,
  TH1*                              histo      ,
  const std::string&                expression ,
  const std::string&                selection  ,
  const unsigned long               first      ,
  const unsigned long               last       ) 
{
  if ( 0 == histo     ) { return Ostap::StatusCode ( 301 ) ; }
  else { histo->Reset() ; } // reset the historgam 
  if ( 0 == data      ) { return Ostap::StatusCode ( 300 ) ; }
  //
  const unsigned long nEntries = std::min ( last , (unsigned long) data->GetEntries() ) ;
  if ( nEntries <= first  ) { return Ostap::StatusCode::RECOVERABLE ; }
  //
  Ostap::Formula xvar ( expression , data ) ;
  if ( !xvar || !xvar.ok() ) { return Ostap::StatusCode ( 303 ) ; }  // RETURN 
  //
  std::unique_ptr<Ostap::Formula> cut { nullptr } ;
  if  ( !selection.empty() ) 
  { 
    cut = std::make_unique<Ostap::Formula>( selection , data ) ; 
    if ( !cut || !cut->ok () ) { return Ostap::StatusCode ( 303 ) ; }  // RETURN
  }
  //
  Ostap::Utils::Notifier notify ( data , &xvar , cut.get() ) ;
  std::vector<double> results {} ;
  /// make an explicit  loop 
  Ostap::Utils::ProgressBar bar (  nEntries - first , progress ) ;
  for ( unsigned long entry = first ; entry < nEntries ; ++entry , ++bar )
  {
    long ievent = data->GetEntryNumber ( entry ) ;
    if ( 0 > ievent ) { return Ostap::StatusCode ( 400 ) ; } // RETURN 
    //
    ievent      = data->LoadTree ( ievent ) ;
    if ( 0 > ievent ) { return Ostap::StatusCode ( 401 ) ; } // RETURN 
    //
    const double weight = cut ? cut -> evaluate() : 1.0 ;
    if ( !weight    ) { continue ;}
    //
    xvar.evaluate ( results ) ;
    for ( double x : results ) { histo -> Fill ( x , weight ) ; }
  }
  return StatusCode::SUCCESS ;
}
// ========================================================================
/** make a projection of DataFrame into the histogram 
 *  @param data        (INPUT)  input data 
 *  @param histo       (UPDATE) histogram 
 *  @param expression  (INPUT)  expression
 *  @param selection   (INPUT)  selection criteria/weight 
 */
// ============================================================================
Ostap::StatusCode Ostap::HistoProject::project
( TTree*                            data       , 
  TH1*                              histo      ,
  const std::string&                expression ,
  const std::string&                selection  ,
  const unsigned long               first      ,
  const unsigned long               last       ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project ( data        , 
                   progress    , 
                   histo       ,
                   expression , 
                   selection   , 
                   first       ,
                   last        ) ;
}
// ============================================================================
/*  make a projection of TTree into the histogram 
 *  @param data        (INPUT)  input data 
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param histo       (UPDATE) histogram 
 *  @param xexpression (INPUT)  x-expression
 *  @param yexpression (INPUT)  y-expression
 *  @param selection   (INPUT)  selection criteria/weight 
 *  @param first       (INPUT)  the first event to process 
 *  @param last        (INPUT)  the last event to process 
 */
// ============================================================================
Ostap::StatusCode Ostap::HistoProject::project2
( TTree*                            data        , 
  const Ostap::Utils::ProgressConf& progress    ,
  TH2*                              histo       ,
  const std::string&                xexpression ,
  const std::string&                yexpression ,
  const std::string&                selection   ,
  const unsigned long               first       ,
  const unsigned long               last        ) 
{
  if ( 0 == histo     ) { return Ostap::StatusCode ( 301 ) ; }
  else { histo->Reset() ; } // reset the historgam 
  if ( 0 == data      ) { return Ostap::StatusCode ( 300 ) ; }
  //
  const unsigned long nEntries = std::min ( last , (unsigned long) data->GetEntries() ) ;
  if ( nEntries <= first  ) { return Ostap::StatusCode::RECOVERABLE ; }
  //
  Ostap::Formula xvar ( xexpression , data ) ;
  if ( !xvar || !xvar.ok() ) { return Ostap::StatusCode ( 303 ) ; }  // RETURN 
  //
  Ostap::Formula yvar ( yexpression , data ) ;
  if ( !yvar || !yvar.ok() ) { return Ostap::StatusCode ( 304 ) ; }  // RETURN 
  //
  std::unique_ptr<Ostap::Formula> cut { nullptr } ;
  if  ( !selection.empty() ) 
  { 
    cut = std::make_unique<Ostap::Formula>( selection , data ) ; 
    if ( !cut || !cut->ok () ) { return Ostap::StatusCode ( 302 ) ; }  // RETURN
  }
  // 
  Ostap::Utils::Notifier notify ( data , &xvar , &yvar , cut.get() ) ;
  std::vector<double> xresults {} ;
  std::vector<double> yresults {} ;
  /// make an explicit  loop 
  Ostap::Utils::ProgressBar bar (  nEntries - first , progress ) ;
  for ( unsigned long entry = first ; entry < nEntries ; ++entry , ++bar )
  {
    long ievent = data->GetEntryNumber ( entry ) ;
    if ( 0 > ievent ) { return Ostap::StatusCode ( 400 ) ; } // RETURN 
    //
    ievent      = data->LoadTree ( ievent ) ;
    if ( 0 > ievent ) { return Ostap::StatusCode ( 401 ) ; } // RETURN 
    //
    const double weight = cut ? cut -> evaluate() : 1.0 ;
    if ( !weight    ) { continue ;}
    //
    xvar.evaluate ( xresults ) ;
    yvar.evaluate ( yresults ) ;
    for ( double x : xresults )
    { for ( double y : yresults )
      { histo -> Fill ( x , y , weight ) ; } }
  }
  return StatusCode::SUCCESS ;
}
// ========================================================================
/*  make a projection of TTree into the histogram 
 *  @param data        (INPUT)  input data 
 *  @param histo       (UPDATE) histogram 
 *  @param xexpression (INPUT)  x-expression
 *  @param yexpression (INPUT)  y-expression
 *  @param selection   (INPUT)  selection criteria/weight 
 *  @param first       (INPUT)  the first event to process 
 *  @param last        (INPUT)  the last event to process 
 */
 // ============================================================================
 Ostap::StatusCode Ostap::HistoProject::project2
( TTree*                            data       , 
  TH2*                              histo       ,
  const std::string&                xexpression ,
  const std::string&                yexpression ,
  const std::string&                selection   ,
  const unsigned long               first       ,
  const unsigned long               last        ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project2 ( data        , 
                    progress    , 
                    histo       ,
                    xexpression ,
                    yexpression , 
                    selection   , 
                    first       ,
                    last        ) ;
}
// ============================================================================
/*  make a projection of TTree into the histogram 
 *  @param data        (INPUT)  input data 
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param histo       (UPDATE) histogram 
 *  @param xexpression (INPUT)  x-expression
 *  @param yexpression (INPUT)  y-expression
 *  @param zexpression (INPUT)  z-expression
 *  @param selection   (INPUT)  selection criteria/weight 
 *  @param first       (INPUT)  the first event to process 
 *  @param last        (INPUT)  the last event to process 
 */
// ============================================================================
Ostap::StatusCode Ostap::HistoProject::project3
( TTree*                            data        , 
  const Ostap::Utils::ProgressConf& progress    ,
  TH3*                              histo       ,
  const std::string&                xexpression ,
  const std::string&                yexpression ,
  const std::string&                zexpression ,
  const std::string&                selection   ,
  const unsigned long               first       ,
  const unsigned long               last        ) 
{
  if ( 0 == histo     ) { return Ostap::StatusCode ( 301 ) ; }
  else { histo->Reset() ; } // reset the historgam 
  if ( 0 == data      ) { return Ostap::StatusCode ( 300 ) ; }
  //
  const unsigned long nEntries = std::min ( last , (unsigned long) data->GetEntries() ) ;
  if ( nEntries <= first  ) { return Ostap::StatusCode::RECOVERABLE ; }
  //
  Ostap::Formula xvar ( xexpression , data ) ;
  if ( !xvar || !xvar.ok() ) { return Ostap::StatusCode ( 303 ) ; }  // RETURN 
  //
  Ostap::Formula yvar ( yexpression , data ) ;
  if ( !yvar || !yvar.ok() ) { return Ostap::StatusCode ( 304 ) ; }  // RETURN 
  //
  Ostap::Formula zvar ( zexpression , data ) ;
  if ( !zvar || !zvar.ok() ) { return Ostap::StatusCode ( 305 ) ; }  // RETURN 
  //
  std::unique_ptr<Ostap::Formula> cut { nullptr } ;
  if  ( !selection.empty() ) 
  { 
    cut = std::make_unique<Ostap::Formula>( selection , data) ; 
    if ( !cut || !cut->ok () ) { return Ostap::StatusCode ( 302 ) ; }  // RETURN
  }
  //
  Ostap::Utils::Notifier notify ( data , &xvar , &yvar , &zvar , cut.get() ) ;
  std::vector<double> xresults {} ;
  std::vector<double> yresults {} ;
  std::vector<double> zresults {} ;
  /// make an explicit  loop 
  Ostap::Utils::ProgressBar bar (  nEntries - first , progress ) ;
  for ( unsigned long entry = first ; entry < nEntries ; ++entry , ++bar )
  {
    long ievent = data->GetEntryNumber ( entry ) ;
    if ( 0 > ievent ) { return Ostap::StatusCode ( 400 ) ; } // RETURN 
    //
    ievent      = data->LoadTree ( ievent ) ;
    if ( 0 > ievent ) { return Ostap::StatusCode ( 401 ) ; } // RETURN 
    //
    const double weight = cut ? cut -> evaluate() : 1.0 ;
    if ( !weight    ) { continue ;}
    //
    xvar.evaluate ( xresults ) ;
    yvar.evaluate ( yresults ) ;
    yvar.evaluate ( zresults ) ;
    for ( double x : xresults )
    { for ( double y : yresults )
      { for ( double z : zresults )
        { histo -> Fill ( x , y , z , weight ) ; } } }
  }
  return StatusCode::SUCCESS ;
}
// ========================================================================
/*  make a projection of TTree into the histogram 
 *  @param data        (INPUT)  input data 
 *  @param histo       (UPDATE) histogram 
 *  @param xexpression (INPUT)  x-expression
 *  @param yexpression (INPUT)  y-expression
 *  @param zexpression (INPUT)  z-expression
 *  @param selection   (INPUT)  selection criteria/weight 
 *  @param first       (INPUT)  the first event to process 
 *  @param last        (INPUT)  the last event to process 
 */
 // ============================================================================
 Ostap::StatusCode Ostap::HistoProject::project3
( TTree*                            data       , 
  TH3*                              histo       ,
  const std::string&                xexpression ,
  const std::string&                yexpression ,
  const std::string&                zexpression ,
  const std::string&                selection   ,
  const unsigned long               first       ,
  const unsigned long               last        ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project3 ( data        , 
                    progress    , 
                    histo       ,
                    xexpression ,
                    yexpression , 
                    zexpression , 
                    selection   , 
                    first       ,
                    last        ) ;
}
// ============================================================================






// ============================================================================
// TTree -> non-histograms  (1D) 
// ============================================================================
/** make a projection of TTree into the object 
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param data        (INPUT)  input data 
 *  @param sum         (UPDATE) sum  
 *  @param expression  (INPUT)  expression
 *  @param selection   (INPUT)  selection criteria/weight 
 *  @param first       (INPUT)  the first event to process 
 *  @param last        (INPUT)  the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project
( TTree*                     data       , 
  Ostap::Math::LegendreSum&  sum        ,
  const std::string&         expression ,
  const std::string&         selection  ,
  const unsigned long        first      ,
  const unsigned long        last       ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project ( data       ,
                   progress   , 
                   sum        , 
                   expression , 
                   selection  , 
                   first      , 
                   last       ) ;
}
// ============================================================================
/** make a projection of TTree into the object 
 *  @param data        (INPUT)  input data 
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param sum         (UPDATE) sum 
 *  @param expression  (INPUT)  expression
 *  @param selection   (INPUT)  selection criteria/weight 
 *  @param first       (INPUT)  the first event to process 
 *  @param last        (INPUT)  the last event to process 
 */

// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project
( TTree*                            data           , 
  const Ostap::Utils::ProgressConf& progress       ,
  Ostap::Math::LegendreSum&         sum            ,
  const std::string&                expression     ,
  const std::string&                selection      ,
  const unsigned long               first          ,
  const unsigned long               last           ) 
{
  /// reset the sum 
  sum *= 0 ;
  return _project_ ( data        , 
                     progress    , 
                     sum         , 
                     expression  , 
                     selection   ,
                     sum.xmin () , 
                     sum.xmax () , 
                     first       , 
                     last        ) ;
}
// ============================================================================
/** make a projection of TTree into the object 
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param data        (INPUT)  input data 
 *  @param sum         (UPDATE) sum  
 *  @param expression  (INPUT)  expression
 *  @param selection   (INPUT)  selection criteria/weight 
 *  @param first       (INPUT)  the first event to process 
 *  @param last        (INPUT)  the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project
( TTree*                      data       , 
  Ostap::Math::ChebyshevSum&  sum        ,
  const std::string&          expression ,
  const std::string&          selection  ,
  const unsigned long         first      ,
  const unsigned long         last       ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project ( data       ,
                   progress   , 
                   sum        , 
                   expression , 
                   selection  , 
                   first      , 
                   last       ) ;
}
// ============================================================================
/** make a projection of TTree into the object 
 *  @param data        (INPUT)  input data 
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param sum         (UPDATE) sum 
 *  @param expression  (INPUT)  expression
 *  @param selection   (INPUT)  selection criteria/weight 
 *  @param first       (INPUT)  the first event to process 
 *  @param last        (INPUT)  the last event to process 
 */

// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project
( TTree*                            data           , 
  const Ostap::Utils::ProgressConf& progress       ,
  Ostap::Math::ChebyshevSum&        sum            ,
  const std::string&                expression     ,
  const std::string&                selection      ,
  const unsigned long               first          ,
  const unsigned long               last           ) 
{
  /// reset the sum 
  sum *= 0 ;
  return _project_ ( data        , 
                     progress    , 
                     sum         , 
                     expression  , 
                     selection   ,
                     sum.xmin () , 
                     sum.xmax () , 
                     first       , 
                     last        ) ;
}// ============================================================================
/** make a projection of TTree into the object 
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param data        (INPUT)  input data 
 *  @param sum         (UPDATE) sum  
 *  @param expression  (INPUT)  expression
 *  @param selection   (INPUT)  selection criteria/weight 
 *  @param first       (INPUT)  the first event to process 
 *  @param last        (INPUT)  the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project
( TTree*                      data       , 
  Ostap::Math::Bernstein&     sum        ,
  const std::string&          expression ,
  const std::string&          selection  ,
  const unsigned long         first      ,
  const unsigned long         last       ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project ( data       ,
                   progress   , 
                   sum        , 
                   expression , 
                   selection  , 
                   first      , 
                   last       ) ;
}
// ============================================================================
/** make a projection of TTree into the object 
 *  @param data        (INPUT)  input data 
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param sum         (UPDATE) sum 
 *  @param expression  (INPUT)  expression
 *  @param selection   (INPUT)  selection criteria/weight 
 *  @param first       (INPUT)  the first event to process 
 *  @param last        (INPUT)  the last event to process 
 */

// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project
( TTree*                            data           , 
  const Ostap::Utils::ProgressConf& progress       ,
  Ostap::Math::Bernstein&           sum            ,
  const std::string&                expression     ,
  const std::string&                selection      ,
  const unsigned long               first          ,
  const unsigned long               last           ) 
{
  /// reset the sum 
  sum *= 0 ;
  return _project_ ( data        , 
                     progress    , 
                     sum         , 
                     expression  , 
                     selection   ,
                     sum.xmin () , 
                     sum.xmax () , 
                     first       , 
                     last        ) ;
}
// ============================================================================
// TTeee -> 2D non-histograms 
// ============================================================================
/** make a projection of TTree into the object
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param data        (INPUT)  input data 
 *  @param sum         (UPDATE) sum 
 *  @param xexpression (INPUT)  x-expression
 *  @param yexpression (INPUT)  y-expression
 *  @param selection   (INPUT)  selection criteria/weight 
 *  @param first       (INPUT)  the first event to process 
 *  @param last        (INPUT)  the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project2
( TTree*                     data         , 
  Ostap::Math::LegendreSum2& sum          ,
  const std::string&         xexpression  ,
  const std::string&         yexpression  ,
  const std::string&         selection    ,
  const unsigned long        first        ,
  const unsigned long        last         ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project2 ( data        ,
                    progress    , 
                    sum         , 
                    xexpression , 
                    yexpression , 
                    selection   , 
                    first       , 
                    last        ) ;
}
// ========================================================================
/** make a projection of TTree into the object 
 *  @param data        (INPUT)  input data 
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param sum         (UPDATE) sum 
 *  @param xexpression (INPUT)  x-expression
 *  @param yexpression (INPUT)  y-expression
 *  @param selection   (INPUT)  selection criteria/weight 
 *  @param first       (INPUT)  the first event to process 
 *  @param last        (INPUT)  the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project2
( TTree*                             data        , 
  const Ostap::Utils::ProgressConf&  progress    ,
  Ostap::Math::LegendreSum2&         sum         ,
  const std::string&                 xexpression ,
  const std::string&                 yexpression ,
  const std::string&                 selection   ,
  const unsigned long                first       ,
  const unsigned long                last        ) 
{
  /// reset the sum 
  sum *= 0 ;
  return _project2_ ( data        , 
                      progress    , 
                      sum         , 
                      xexpression , 
                      yexpression , 
                      selection   ,
                      sum.xmin () , 
                      sum.xmax () , 
                      sum.ymin () , 
                      sum.ymax () , 
                      first       , 
                      last        ) ;
}
// ============================================================================
/** make a projection of TTree into the object
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param data        (INPUT)  input data 
 *  @param sum         (UPDATE) sum 
 *  @param xexpression (INPUT)  x-expression
 *  @param yexpression (INPUT)  y-expression
 *  @param selection   (INPUT)  selection criteria/weight 
 *  @param first       (INPUT)  the first event to process 
 *  @param last        (INPUT)  the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project2
( TTree*                     data         , 
  Ostap::Math::Bernstein2D&  sum          ,
  const std::string&         xexpression  ,
  const std::string&         yexpression  ,
  const std::string&         selection    ,
  const unsigned long        first        ,
  const unsigned long        last         ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project2 ( data        ,
                    progress    , 
                    sum         , 
                    xexpression , 
                    yexpression , 
                    selection   , 
                    first       , 
                    last        ) ;
}
// ========================================================================
/** make a projection of TTree into the object 
 *  @param data        (INPUT)  input data 
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param sum         (UPDATE) sum 
 *  @param xexpression (INPUT)  x-expression
 *  @param yexpression (INPUT)  y-expression
 *  @param selection   (INPUT)  selection criteria/weight 
 *  @param first       (INPUT)  the first event to process 
 *  @param last        (INPUT)  the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project2
( TTree*                               data        , 
  const Ostap::Utils::ProgressConf&    progress    ,
  Ostap::Math::Bernstein2D&            sum         ,
  const std::string&                   xexpression ,
  const std::string&                   yexpression ,
  const std::string&                   selection   ,
  const unsigned long                  first       ,
  const unsigned long                  last        ) 
{
  /// reset the sum 
  sum *= 0 ;
  return _project2_ ( data        , 
                      progress    , 
                      sum         , 
                      xexpression , 
                      yexpression , 
                      selection   ,
                      sum.xmin () , 
                      sum.xmax () , 
                      sum.ymin () , 
                      sum.ymax () , 
                      first       , 
                      last        ) ;
}
// ============================================================================
// TTeee -> 3D non-histograms 
// ============================================================================
/** make a projection of TTree into the object
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param data        (INPUT)  input data 
 *  @param sum         (UPDATE) sum 
 *  @param xexpression (INPUT)  x-expression
 *  @param yexpression (INPUT)  y-expression
 *  @param zexpression (INPUT)  y-expression
 *  @param selection   (INPUT)  selection criteria/weight 
 *  @param first       (INPUT)  the first event to process 
 *  @param last        (INPUT)  the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project3
( TTree*                     data         , 
  Ostap::Math::LegendreSum3& sum          ,
  const std::string&         xexpression  ,
  const std::string&         yexpression  ,
  const std::string&         zexpression  ,
  const std::string&         selection    ,
  const unsigned long        first        ,
  const unsigned long        last         ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project3 ( data        ,
                    progress    , 
                    sum         , 
                    xexpression , 
                    yexpression , 
                    zexpression , 
                    selection   , 
                    first       , 
                    last        ) ;
}
// ========================================================================
/** make a projection of TTree into the object 
 *  @param data        (INPUT)  input data 
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param sum         (UPDATE) sum 
 *  @param xexpression (INPUT)  x-expression
 *  @param yexpression (INPUT)  y-expression
 *  @param yexpression (INPUT)  y-expression
 *  @param selection   (INPUT)  selection criteria/weight 
 *  @param first       (INPUT)  the first event to process 
 *  @param last        (INPUT)  the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project3
( TTree*                               data        , 
  const Ostap::Utils::ProgressConf&    progress    ,
  Ostap::Math::LegendreSum3&           sum         ,
  const std::string&                   xexpression ,
  const std::string&                   yexpression ,
  const std::string&                   zexpression ,
  const std::string&                   selection   ,
  const unsigned long                  first       ,
  const unsigned long                  last        ) 
{
  /// reset the sum 
  sum *= 0 ;
  return _project3_ ( data        , 
                      progress    , 
                      sum         , 
                      xexpression , 
                      yexpression , 
                      zexpression , 
                      selection   ,
                      sum.xmin () , 
                      sum.xmax () , 
                      sum.ymin () , 
                      sum.ymax () , 
                      sum.zmin () , 
                      sum.zmax () , 
                      first       , 
                      last        ) ;
}
// ============================================================================
/*  make a projection of TTree into the object
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param data        (INPUT)  input data 
 *  @param sum         (UPDATE) sum 
 *  @param xexpression (INPUT)  x-expression
 *  @param yexpression (INPUT)  y-expression
 *  @param zexpression (INPUT)  z-expression
 *  @param selection   (INPUT)  selection criteria/weight 
 *  @param first       (INPUT)  the first event to process 
 *  @param last        (INPUT)  the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project3
( TTree*                     data         , 
  Ostap::Math::Bernstein3D&  sum          ,
  const std::string&         xexpression  ,
  const std::string&         yexpression  ,
  const std::string&         zexpression  ,
  const std::string&         selection    ,
  const unsigned long        first        ,
  const unsigned long        last         ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project3 ( data        ,
                    progress    , 
                    sum         , 
                    xexpression , 
                    yexpression , 
                    zexpression , 
                    selection   , 
                    first       , 
                    last        ) ;
}


// ========================================================================
/*  make a projection of TTree into the object 
 *  @param data        (INPUT)  input data 
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param sum         (UPDATE) sum 
 *  @param xexpression (INPUT)  x-expression
 *  @param yexpression (INPUT)  y-expression
 *  @param zexpression (INPUT)  z-expression
 *  @param selection   (INPUT)  selection criteria/weight 
 *  @param first       (INPUT)  the first event to process 
 *  @param last        (INPUT)  the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project3
( TTree*                               data        , 
  const Ostap::Utils::ProgressConf&    progress    ,
  Ostap::Math::Bernstein3D&            sum         ,
  const std::string&                   xexpression ,
  const std::string&                   yexpression ,
  const std::string&                   zexpression ,
  const std::string&                   selection   ,
  const unsigned long                  first       ,
  const unsigned long                  last        ) 
{
  /// reset the sum 
  sum *= 0 ;
  return _project3_ ( data        , 
                      progress    , 
                      sum         , 
                      xexpression , 
                      yexpression , 
                      zexpression , 
                      selection   ,
                      sum.xmin () , 
                      sum.xmax () , 
                      sum.ymin () , 
                      sum.ymax () , 
                      sum.zmin () , 
                      sum.zmax () , 
                      first       , 
                      last        ) ;
}
// ============================================================================
// TTeee -> 4D non-histograms 
// ============================================================================
/** make a projection of TTree into the object
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param data        (INPUT)  input data 
 *  @param sum         (UPDATE) sum 
 *  @param xexpression (INPUT)  x-expression
 *  @param yexpression (INPUT)  y-expression
 *  @param zexpression (INPUT)  y-expression
 *  @param uexpression (INPUT)  u-expression
 *  @param selection   (INPUT)  selection criteria/weight 
 *  @param first       (INPUT)  the first event to process 
 *  @param last        (INPUT)  the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project4
( TTree*                     data         , 
  Ostap::Math::LegendreSum4& sum          ,
  const std::string&         xexpression  ,
  const std::string&         yexpression  ,
  const std::string&         zexpression  ,
  const std::string&         uexpression  ,
  const std::string&         selection    ,
  const unsigned long        first        ,
  const unsigned long        last         ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return project4 ( data        ,
                    progress    , 
                    sum         , 
                    xexpression , 
                    yexpression , 
                    zexpression , 
                    uexpression , 
                    selection   , 
                    first       , 
                    last        ) ;
}
// ========================================================================
/*  make a projection of TTree into the object 
 *  @param data        (INPUT)  input data 
 *  @param progress    (INPUT)  configuration of progres bar 
 *  @param sum         (UPDATE) sum 
 *  @param xexpression (INPUT)  x-expression
 *  @param yexpression (INPUT)  y-expression
 *  @param yexpression (INPUT)  y-expression
 *  @param selection   (INPUT)  selection criteria/weight 
 *  @param first       (INPUT)  the first event to process 
 *  @param last        (INPUT)  the last event to process 
 */
// ============================================================================
Ostap::StatusCode 
Ostap::HistoProject::project4
( TTree*                               data        , 
  const Ostap::Utils::ProgressConf&    progress    ,
  Ostap::Math::LegendreSum4&           sum         ,
  const std::string&                   xexpression ,
  const std::string&                   yexpression ,
  const std::string&                   zexpression ,
  const std::string&                   uexpression ,
  const std::string&                   selection   ,
  const unsigned long                  first       ,
  const unsigned long                  last        ) 
{
  /// reset the sum 
  sum *= 0 ;
  return _project4_ ( data        , 
                      progress    , 
                      sum         , 
                      xexpression , 
                      yexpression , 
                      zexpression , 
                      uexpression , 
                      selection   ,
                      sum.xmin () , 
                      sum.xmax () , 
                      sum.ymin () , 
                      sum.ymax () , 
                      sum.zmin () , 
                      sum.zmax () , 
                      sum.umin () , 
                      sum.umax () , 
                      first       , 
                      last        ) ;
}

// ============================================================================
//                                                                      The END 
// ============================================================================
