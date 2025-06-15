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
  static_assert ( std::numeric_limits<unsigned long>::is_specialized   ,
                  "Numeric_limist<unsigned long> are not specialized!" ) ;
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
                    "Ostap::StatVar::make_formula" ) ;
    //
    RooArgList        alst { *aset } ;
    //
    auto result = std::make_unique<Ostap::FormulaVar> ( expression , alst , false ) ;
    //
    if ( allow_null && ( !result || !result->ok() ) ) { return nullptr ; } 
    //
    Ostap::Assert ( result && result->ok()                   , 
                    "Invalid formula:\"" + expression + "\"" , 
                    "Ostap::StatVar::make_formula"           ) ;
    return result ;
  }
  // ==========================================================================
  /** get the number of equivalent entries 
   *  \f$ n_{eff} \equiv = \frac{ (\sum w)^2}{ \sum w^2} \f$
   */
  double _neff_  
  ( TTree&               tree  , 
    Ostap::Formula*      cuts  , 
    const unsigned long  first ,
    const unsigned long  last  ) 
  {
    //
    const unsigned long nEntries = std::min ( last , (unsigned long) tree.GetEntries() ) ;
    if ( last <= first  ) { return 0 ; }                // RETURN 
    // 
    if ( !cuts ) { return nEntries - first ; }          // RETURN    
    //
    Ostap::Utils::Notifier notify ( &tree , cuts ) ;
    //
    long double sumw  = 0     ;
    long double sumw2 = 0     ;
    bool        empty = false ;
    // 
    for ( unsigned long entry = first ; entry < nEntries ; ++entry ) 
    {      
      long ievent = tree.GetEntryNumber ( entry ) ;
      if ( 0 > ievent ) { break ; }                        // BREAK
      //

      ievent      = tree.LoadTree ( ievent ) ;
      if ( 0 > ievent ) { break ; }                        // BREAK
      //
      const long double w = cuts->evaluate() ;
      //
      if ( !w  ) { continue ; }                           // CONTINUE       
      //
      sumw  += w     ;
      sumw2 += w * w ;
      empty  = false ;
    }
    //
    return empty ? 0 : sumw *  sumw / sumw2 ; // RETURN 
  }
  // ==========================================================================
  /** calculate the moment of order "order" relative to the center "center"
   *  @param  tree   (INPUT) input tree 
   *  @param  expr   (INPUT) expression  (must  be valid TFormula!)
   *  @param  cuts   (INPUT) cuts 
   *  @param  order  (INPUT) the order 
   *  @param  center (INPUT) the center 
   *  @param  first  (INPUT) the first  event to process 
   *  @param  last   (INPUT) the last event to  process
   *  @return the moment 
   */
  double _moment1_ ( TTree&               tree   , 
                     Ostap::Formula&   var    ,
                     Ostap::Formula*   cuts   , 
                     const unsigned short order  , 
                     const double         center , 
                     const unsigned long  first  , 
                     const unsigned long  last   )
  {
    // 
    if ( 0 == order     ) { return  1 ; }                // RETURN
    // the loop 
    const unsigned long nEntries = std::min ( last , (unsigned long) tree.GetEntries() ) ;
    if ( last <= first  ) { return  0 ; }                // RETURN ???    
    //
    Ostap::Utils::Notifier notify ( &tree , &var , cuts ) ;
    const bool with_cuts = nullptr != cuts ? true : false ;
    //
    long double         mom     = 0      ;
    long double         sumw    = 0      ;
    bool                empty   = true   ;
    const long double   v0      = center ;
    std::vector<double> results {} ;
    for ( unsigned long entry = first ; entry < nEntries ; ++entry ) 
    {      
      long ievent = tree.GetEntryNumber ( entry ) ;
      if ( 0 > ievent ) { break ; }                        // BREAK
      //
      ievent      = tree.LoadTree ( ievent ) ;
      if ( 0 > ievent ) { break ; }                        // BREAK
      //
      const long double w = with_cuts ? cuts->evaluate() : 1.0L ;
      //
      if  ( !w ) { continue ; }                            // ATTENTION!
      //
      var.evaluate ( results ) ;
      for ( const long double r : results ) 
      {
        const long double dx =  r -  v0 ;
        //
        mom   += w  * std::pow ( dx , order ) ;
        sumw  += w ;
        empty  = false ;
      } 
    }
    //
    return empty ? 0 : mom / sumw  ;
    // ========================================================================
  }
  // ==========================================================================
  double _moment_ ( const RooAbsData&    data      , 
                    const RooAbsReal&    expr      , 
                    const RooAbsReal*    cuts      , 
                    const unsigned short order     , 
                    const double         center    , 
                    const unsigned long  first     , 
                    const unsigned long  last      , 
                    const char*          cut_range )
  {
    //
    if (  0 == order       ) { return  1 ; }    // RETURN
    //
    const bool  weighted = data.isWeighted () ;
    //
    long double mom   = 0    ;
    long double sumw  = 0    ;
    bool        empty = true ;
    //  
    for ( unsigned long entry = first ; entry < last ; ++entry )
    {
      const RooArgSet* vars = data.get( entry ) ;
      if ( nullptr == vars )                              { break    ; } // BREAK 
      //
      if ( cut_range && !vars->allInRange ( cut_range ) ) { continue ; } // CONTINUE    
      // apply cuts:
      const long double wc = nullptr != cuts ? cuts -> getVal() : 1.0L ;
      if ( !wc ) { continue ; }                                   // CONTINUE  
      // apply weight:
      const long double wd = weighted  ? data.weight()   : 1.0L ;
      if ( !wd ) { continue ; }                                   // CONTINUE    
      // cuts & weight:
      const long double w  = wd *  wc ;
      if ( !w  ) { continue ; }                                   // CONTINUE        
      //
      const double dx = expr.getVal() - center ;
      //
      mom   += w * std::pow ( dx , order ) ;
      sumw  += w     ;
      empty  = false ;
      //
    }
    //
    return empty ? 0.0 : mom / sumw ;
  }  
  // ==========================================================================
  /** calculate the moment of order "order"
   *  @param  tree  (INPUT) input tree 
   *  @param  order (INPUT) the order 
   *  @param  expr  (INPUT) expression  (must  be valid TFormula!)
   *  @param  cuts  (INPUT) cuts 
   *  @param  first (INPUT) the first  event to process 
   *  @param  last  (INPUT) the last event to  process
   *  @return the moment 
   */ 
  Ostap::Math::ValueWithError
  _moment2_
  ( TTree&               tree  ,  
    const unsigned short order ,
    Ostap::Formula&   var   ,
    Ostap::Formula*   cuts  , 
    const unsigned long  first , 
    const unsigned long  last  )
  {
    // 
    if ( 0 == order     ) { return  1 ; }                // RETURN
    // the loop 
    const unsigned long nEntries = std::min ( last , (unsigned long) tree.GetEntries() ) ;
    if ( last <= first  )
    { return Ostap::Math::ValueWithError ( -1 , -1 ) ;  } // RETURN ???    
    //
    Ostap::Utils::Notifier notify ( &tree , &var , cuts ) ;
    const bool with_cuts = nullptr != cuts ? true : false ;
    //
    long double         mom   = 0    ;
    long double         sumw  = 0    ; // sum of weights 
    // for uncertainties 
    long double         sumw2 = 0    ; // sum of weights^2
    long double         c2    = 0    ;
    double              empty = true ;
    std::vector<double> results {}   ;
    for ( unsigned long entry = first ; entry < nEntries ; ++entry ) 
    {      
      long ievent = tree.GetEntryNumber ( entry ) ;
      if ( 0 > ievent ) { break ; }                        // BREAK
      //
      ievent      = tree.LoadTree ( ievent ) ;
      if ( 0 > ievent ) { break ; }                        // BREAK
      //
      const long double w = with_cuts ? cuts->evaluate() : 1.0L ;
      //
      if  ( !w ) { continue ; }                            // ATTENTION!
      //
      var.evaluate (  results ) ;
      for ( const long double r : results ) 
      {
        const long double x =  r ;
        //
        mom   += w  * std::pow ( x , order ) ;
        sumw  += w     ;
        // for uncertainty:
        sumw2 += w * w ;
        c2    += w * std::pow ( x , 2 * order ) ;
        // 
        empty  =  false ;
      } 
    }
    //
    if  ( empty ) { return 0 ; }   //    RETURN
    //
    const long double v = mom / sumw ;
    //
    c2 /= sumw  ; // the moment of "2*order"
    c2 -= v * v ; // m(2*order) - m(order)**2
    //
    const long double n = sumw * sumw / sumw2 ;
    c2 /= n     ; // 
    //
    return Ostap::Math::ValueWithError ( v , c2 ) ;            // RETURN
  }
  // ==========================================================================
  /** calculate the central moment of order "order"
   *  @param  tree  (INPUT) input tree 
   *  @param  order (INPUT) the order 
   *  @param  expr  (INPUT) expression  (must  be valid TFormula!)
   *  @param  cuts  (INPUT) cuts 
   *  @param  first (INPUT) the first  event to process 
   *  @param  last  (INPUT) the last event to  process
   *  @return the moment 
   */ 
  Ostap::Math::ValueWithError
  _moment3_
  ( TTree&               tree  ,  
    const unsigned short order ,
    Ostap::Formula&      var   ,
    Ostap::Formula*      cuts  , 
    const unsigned long  first , 
    const unsigned long  last  )
  {
    // 
    if ( 0 == order     ) { return  1 ; }                // RETURN
    // the loop 
    const unsigned long nEntries = std::min ( last , (unsigned long) tree.GetEntries() ) ;
    if ( last <= first  )
    { return Ostap::Math::ValueWithError ( -1 , -1 ) ;  } // RETURN ???    
    //
    // calculate mean value
    const long double mean = _moment1_  ( tree , var , cuts , 1 , 0 , first , last ) ;
    //
    Ostap::Utils::Notifier notify ( &tree , &var , cuts ) ;
    const bool with_cuts = nullptr != cuts ? true : false ;
    //
    long double         mom   = 0    ;
    long double         sumw  = 0    ; // sum of weights 
    // for uncertainty:
    long double         sumw2 = 0    ; // sum of weights^2
    long double         m2o   = 0    ; // moment of 2*order 
    long double         mm1   = 0    ; // moment of   order-1
    long double         mp1   = 0    ; // moment of   order+1
    long double         m2    = 0    ; // moment of 2
    bool                empty = true ;
    std::vector<double> results ;
    for ( unsigned long entry = first ; entry < nEntries ; ++entry ) 
    {      
      long ievent = tree.GetEntryNumber ( entry ) ;
      if ( 0 > ievent ) { break ; }                        // BREAK
      //
      ievent      = tree.LoadTree ( ievent ) ;
      if ( 0 > ievent ) { break ; }                        // BREAK
      //
      const long double w = with_cuts ? cuts->evaluate() : 1.0L ;
      //
      if  ( !w ) { continue ; }                            // ATTENTION!
      //
      var.evaluate ( results ) ;
      for ( const long double r : results ) 
      {
        const long double dx =  r -  mean ;
        //
        mom   += w  * std::pow ( dx , order ) ;
        sumw  += w     ;
        // for uncertainty:
        sumw2 += w * w ;
        m2o   += w * std::pow ( dx , 2 * order     ) ; 
        mm1   += w * std::pow ( dx ,     order - 1 ) ; 
        mp1   += w * std::pow ( dx ,     order + 1 ) ; 
        m2    += w * std::pow ( dx , 2             ) ;
        //
        empty  = false ;
      }
    }
    //
    if  ( empty ) { return 0 ; } //  RETURN
    //
    // number of effective entries:
    const long double n = sumw * sumw / sumw2 ;
    long double v = mom / sumw ;
    /// correct O(1/n) bias  for   3rd and 4th order moments :
    if      ( 3 == order ) { v *=  n * n / ( ( n - 1  ) * ( n - 2 ) ) ; }
    else if ( 4 == order ) 
    {
      const double n0 =  ( n - 1 ) * ( n - 2 ) * ( n - 3 ) ;
      const double n1 =  n * ( n * n - 2 * n +  3 ) / n0   ;
      const double n2 =  3 * n * ( 2 * n - 3 )      / n0   ;
      v = n1 * v - n2 * m2 * m2 / ( sumw * sumw ) ;
    }
    //
    m2o /= sumw ;
    mm1 /= sumw ;
    mp1 /= sumw ;
    m2  /= sumw ;
    //
    long double c2 = m2o  ;
    c2 -= 2 * order * mm1 * mp1 ;
    c2 -= v * v ;
    c2 += order * order * m2 * mm1 * mm1 ;
    c2 /= n  ;
    //
    return Ostap::Math::ValueWithError ( v , c2 ) ;            // RETURN
  }
  // ==========================================================================
  /*  calculate the skewness of the  distribution
   *  @param  tree  (INPUT) input tree 
   *  @param  expr  (INPUT) expression  (must  be valid TFormula!)
   *  @param  cuts  (INPUT) cuts 
   *  @param  first (INPUT) the first  event to process 
   *  @param  last  (INPUT) the last event to  process
   *  @return the skewness 
   */
  Ostap::Math::ValueWithError
  _skewness_
  ( TTree&               tree  ,  
    Ostap::Formula&   var   ,
    Ostap::Formula*   cuts  , 
    const unsigned long  first , 
    const unsigned long  last  )
  {
    // the loop 
    const unsigned long nEntries = std::min ( last , (unsigned long) tree.GetEntries() ) ;
    if ( last <= first  ) { return 0 ;  } // RETURN ???    
    //
    // calculate mean value
    const long double mean = _moment1_  ( tree , var , cuts , 1 , 0 , first , last ) ;
    //
    Ostap::Utils::Notifier notify ( &tree , &var , cuts ) ;
    const bool with_cuts = nullptr != cuts ? true : false ;
    //
    long double          mom   = 0    ;
    long double          sumw  = 0    ; // sum of weights 
    // for uncertainty:
    long double          sumw2 = 0    ; // sum of weights^2
    long double          m2    = 0    ; // moment of 2
    bool                 empty = true ;
    std::vector<double>  results {} ;
    for ( unsigned long entry = first ; entry < nEntries ; ++entry ) 
    {      
      long ievent = tree.GetEntryNumber ( entry ) ;
      if ( 0 > ievent ) { break ; }                        // BREAK
      //
      ievent      = tree.LoadTree ( ievent ) ;
      if ( 0 > ievent ) { break ; }                        // BREAK
      //
      const long double w = with_cuts ? cuts->evaluate() : 1.0L ;
      //
      if  ( !w ) { continue ; }                            // ATTENTION!
      //
      var.evaluate ( results ) ;
      for ( const long double r : results ) 
      {
        const long double dx =  r -  mean ;
        //
        mom   += w * std::pow ( dx , 3 ) ;
        sumw  += w     ;
        // for uncertainty:
        sumw2 += w * w ;
        m2    += w * std::pow ( dx , 2 ) ;
        //
        empty  = false ;
      }
    }
    //
    if  (  empty ) { return 0 ; }
    //
    // number of effective entries:
    const long double n = sumw * sumw / sumw2 ;
    //
    long double v = mom / sumw ;
    /// correct O(1/n) bias  for 3rd moment 
    v *=  n * n / ( ( n - 1  ) * ( n - 2 ) ) ; 
    //
    m2 /= sumw  ;
    v  /= std::pow ( m2 , 1.5 ) ;
    //
    long double c2 = 6  ;
    c2 *= ( n - 2 )  ;
    c2 /= ( n + 1 ) * ( n + 3 ) ;    
    //
    return Ostap::Math::ValueWithError ( v , c2 ) ;            // RETURN
  }
  // ==============================================================================
  /*  calculate the (excess) kurtosis of the  distribution
   *  @param  tree  (INPUT) input tree 
   *  @param  expr  (INPUT) expression  (must  be valid TFormula!)
   *  @param  cuts  (INPUT) cuts 
   *  @param  first (INPUT) the first  event to process 
   *  @param  last  (INPUT) the last event to  process
   *  @return the (excess) kurtosis
   */
  // ============================================================================
  Ostap::Math::ValueWithError
  _kurtosis_
  ( TTree&               tree  ,  
    Ostap::Formula&   var   ,
    Ostap::Formula*   cuts  , 
    const unsigned long  first , 
    const unsigned long  last  )
  {
    // the loop 
    const unsigned long nEntries = std::min ( last , (unsigned long) tree.GetEntries() ) ;
    if ( last <= first  ) { return 0 ;  } // RETURN ???    
    //
    // calculate mean value
    const long double mean = _moment1_  ( tree , var , cuts , 1 , 0 , first , last ) ;
    //
    Ostap::Utils::Notifier notify ( &tree , &var , cuts ) ;
    const bool with_cuts = nullptr != cuts ? true : false ;
    //
    long double         mom   = 0    ;
    long double         sumw  = 0    ; // sum of weights 
    // for uncertainty:
    long double         sumw2 = 0    ; // sum of weights^2
    long double         m2    = 0    ; // moment of 2
    bool                empty = true ;
    std::vector<double> results {} ;
    for ( unsigned long entry = first ; entry < nEntries ; ++entry ) 
    {      
      long ievent = tree.GetEntryNumber ( entry ) ;
      if ( 0 > ievent ) { break ; }                        // BREAK
      //
      ievent      = tree.LoadTree ( ievent ) ;
      if ( 0 > ievent ) { break ; }                        // BREAK
      //
      const long double w = with_cuts ? cuts->evaluate() : 1.0L ;
      //
      if  ( !w ) { continue ; }                            // ATTENTION!
      //
      var.evaluate ( results ) ;
      for ( const long double r : results ) 
      { 
        const long double dx =  r -  mean ;
        //
        mom   += w * std::pow ( dx , 4 ) ;
        sumw  += w     ;
        // for uncertainty:
        sumw2 += w * w ;
        m2    += w * std::pow ( dx , 2 ) ;
        //
        empty  = false ;
      } 
    }
    //
    if ( empty ) { return 0 ; } // RETURN
    //
    // number of effective entries:
    const long double n = sumw * sumw / sumw2 ;
    //
    long double v = mom / sumw ;
    m2 /= sumw  ; // second order moment
    /// correct for O(1/n) bias:
    const double n0 =  ( n - 1 ) * ( n - 2 ) * ( n - 3 ) ;
    const double n1 =  n * ( n * n - 2 * n +  3 ) / n0   ;
    const double n2 =  3 * n * ( 2 * n - 3 )      / n0   ;
    v  = n1 * v - n2 * m2 * m2 ;
    /// normalize  it: 
    v /= std::pow ( m2 , 2 ) ;
    ///
    long double c2 = 24 * n ;
    c2 *= ( n - 2 ) * ( n - 3 ) ;
    c2 /= ( n + 1 ) * ( n + 1 ) ;
    c2 /= ( n + 3 ) * ( n + 5 ) ;
    //
    return Ostap::Math::ValueWithError ( v , c2 ) ;            // RETURN
  }
  // ==========================================================================
  /*   get quantile of the distribution  
   *   @param tree  (INPUT) the input tree 
   *   @param q     (INPUT) quantile value   0 < q < 1  
   *   @param expr  (INPUT) the expression 
   *   @param cuts  (INPUT) selection cuts 
   *   @param  first (INPUT) the first  event to process 
   *   @param  last  (INPUT) the last event to  process
   *   @return the quantile value 
   */
  Ostap::StatVar::Quantiles
  _quantiles_
  ( TTree&                  tree      ,
    const std::set<double>& quantiles , //  0<q<1 
    Ostap::Formula&         var       ,
    Ostap::Formula*         cuts      , 
    const unsigned long     first     ,
    const unsigned long     last      ) 
  {
    // the loop 
    const unsigned long the_last = std::min ( last , (unsigned long) tree.GetEntries() ) ;
    //
    Ostap::Utils::Notifier notify ( &tree , &var , cuts ) ;
    const bool with_cuts = nullptr != cuts ? true : false ;
    //
    unsigned long num = 0 ;
    for ( unsigned long entry = first ; entry < the_last ; ++entry ) 
    {      
      long ievent = tree.GetEntryNumber ( entry ) ;
      if ( 0 > ievent ) { break ; }                        // BREAK
      //
      ievent      = tree.LoadTree ( ievent ) ;
      if ( 0 > ievent ) { break ; }                        // BREAK
      //
      const long double w = with_cuts ? cuts->evaluate() : 1.0L ;
      //
      if ( !w  ) { continue ; }                           // CONTINUE       
      //
      ++num ;
    }
    //
    if ( 0 == num ) { return Ostap::StatVar::Quantiles ( std::vector<double>() , num ) ; }
    //
    typedef std::vector<double> VALUES ;
    VALUES values{} ; values.reserve ( num ) ;
    //
    std::vector<double> results {} ;
    for ( unsigned long entry = first ; entry < the_last ; ++entry ) 
    {      
      long ievent = tree.GetEntryNumber ( entry ) ;
      if ( 0 > ievent ) { break ; }                        // BREAK
      //
      ievent      = tree.LoadTree ( ievent ) ;
      if ( 0 > ievent ) { break ; }                        // BREAK
      //
      const long double w = with_cuts ? cuts->evaluate() : 1.0L ;
      //
      if ( !w  ) { continue ; }                           // CONTINUE       
      //
      var.evaluate  ( results ) ;
      values.insert ( values.end() , results.begin() , results.end() ) ;
      //
    }
    //
    std::vector<double> result ; result.reserve ( quantiles.size() ) ;
    //
    VALUES::iterator start = values.begin() ;
    for ( const double q : quantiles )
    {
      const unsigned long current = values.size() * q ;
      std::nth_element  ( start , values.begin () + current , values.end () ) ;
      start = values.begin() + current ;
      result.push_back ( *start ) ;
    }
    //
    return Ostap::StatVar::Quantiles ( result , values.size () ) ; 
  }
  // ==========================================================================
  /*   get quantile of the distribution  
   *   @param tree  (INPUT) the input tree 
   *   @param q     (INPUT) quantile value   0 < q < 1  
   *   @param expr  (INPUT) the expression 
   *   @param cuts  (INPUT) selection cuts 
   *   @param  first (INPUT) the first  event to process 
   *   @param  last  (INPUT) the last event to  process
   *   @return the quantile value 
   */
  Ostap::StatVar::Quantiles
  _p2quantiles_
  ( TTree&                  tree      ,
    const std::set<double>& quantiles , //  0<q<1 
    Ostap::Formula&         var       ,
    Ostap::Formula*         cuts      , 
    const unsigned long     first     ,
    const unsigned long     last      ) 
  {
    // the loop 
    const unsigned long the_last = std::min ( last , (unsigned long) tree.GetEntries() ) ;
    //
    Ostap::Utils::Notifier notify ( &tree , &var , cuts ) ;
    const bool with_cuts = nullptr != cuts ? true : false ;
    //
    std::vector<Ostap::Math::GSL::P2Quantile> qs ( quantiles.begin() , quantiles.end() ) ;
    //
    unsigned long num = 0  ;
    std::vector<double> results {} ;
    for ( unsigned long entry = first ; entry < the_last ; ++entry ) 
    {      
      long ievent = tree.GetEntryNumber ( entry ) ;
      if ( 0 > ievent ) { break ; }                        // BREAK
      //
      ievent      = tree.LoadTree ( ievent ) ;
      if ( 0 > ievent ) { break ; }                        // BREAK
      //
      const long double w = with_cuts ? cuts->evaluate() : 1.0L ;
      //
      if ( !w  ) { continue ; }                           // CONTINUE       
      //
      var.evaluate  ( results ) ;
      for ( auto& q : qs ) { q.add  ( results.begin() , results.end () ) ; }
      num += results.size();
    }
    //
    return Ostap::StatVar::Quantiles ( std::vector<double>( qs.begin (), qs.end () )  , num ) ;
  }
  // ==========================================================================
  /*   get quantile of the distribution  
   *   @param tree  (INPUT) the input tree 
   *   @param q     (INPUT) quantile value   0 < q < 1  
   *   @param expr  (INPUT) the expression 
   *   @param cuts  (INPUT) selection cuts 
   *   @param  first (INPUT) the first  event to process 
   *   @param  last  (INPUT) the last event to  process
   *   @return the quantile value 
   */
  Ostap::StatVar::Quantiles 
  _quantiles_
  ( const RooAbsData&       data      ,
    const std::set<double>& quantiles , //  0<q<1 
    const RooAbsReal&       var       ,
    const RooAbsReal*       cuts      , 
    const unsigned long     first     ,
    const unsigned long     last      , 
    const char*             cut_range ) 
  {
    // the loop 
    const unsigned long the_last = std::min ( last , (unsigned long) data.numEntries() ) ;
    //
    const bool  weighted = data.isWeighted () ;
    //
    unsigned long num = 0 ;
    for ( unsigned long entry = first ; entry < the_last ; ++entry )
    {
      const RooArgSet* vars = data.get( entry ) ;
      if ( nullptr == vars )                              { break    ; } // BREAK 
      //
      if ( cut_range && !vars->allInRange ( cut_range ) ) { continue ; } // CONTINUE    
      // apply cuts:
      const long double wc = nullptr != cuts ? cuts -> getVal() : 1.0L ;
      if ( !wc ) { continue ; }                                          // CONTINUE  
      // apply weight:
      const long double wd = weighted  ? data.weight()   : 1.0L ;
      if ( !wd ) { continue ; }                                          // CONTINUE    
      // cuts & weight:
      const long double w  = wd *  wc ; 
      if ( !w  ) { continue ; }                                          // CONTINUE        
      //
      ++num ;
    }
    //
    typedef std::vector<double> VALUES ;
    VALUES values{} ; values.reserve ( num ) ;
    //
    for ( unsigned long entry = first ; entry < the_last ; ++entry )
    {
      const RooArgSet* vars = data.get( entry ) ;
      if ( nullptr == vars )                              { break    ; } // BREAK 
      //
      if ( cut_range && !vars->allInRange ( cut_range ) ) { continue ; } // CONTINUE    
      // apply cuts:
      const long double wc = nullptr != cuts ? cuts -> getVal() : 1.0L ;
      if ( !wc ) { continue ; }                                          // CONTINUE  
      // apply weight:
      const long double wd = weighted  ? data.weight()   : 1.0L ;
      if ( !wd ) { continue ; }                                          // CONTINUE    
      // cuts & weight:
      const long double w  = wd *  wc ; 
      if ( !w  ) { continue ; }                                          // CONTINUE        
      //
      values.push_back ( var.getVal() ) ;
    }
    //
    std::vector<double> result ; result.reserve ( quantiles.size() ) ;
    //
    VALUES::iterator start = values.begin() ;
    for ( const double q : quantiles )
    {
      const unsigned long current = values.size() * q ;
      std::nth_element  ( start , values.begin () + current , values.end () ) ;
      start = values.begin() + current ;
      result.push_back ( *start ) ;
    }
    //
    return Ostap::StatVar::Quantiles ( result , values.size () ) ; 
  }
  // ==========================================================================
  /*   get (Appeoximate) quantile of the distribution using  P^2 algorithm  
   *   @param tree  (INPUT) the input tree 
   *   @param q     (INPUT) quantile value   0 < q < 1  
   *   @param expr  (INPUT) the expression 
   *   @param cuts  (INPUT) selection cuts 
   *   @param  first (INPUT) the first  event to process 
   *   @param  last  (INPUT) the last event to  process
   *   @return the quantile value 
   */
  Ostap::StatVar::Quantiles 
  _p2quantiles_
  ( const RooAbsData&       data      ,
    const std::set<double>& quantiles , //  0<q<1 
    const RooAbsReal&       var       ,
    const RooAbsReal*       cuts      , 
    const unsigned long     first     ,
    const unsigned long     last      , 
    const char*             cut_range ) 
  {
    // the loop 
    const unsigned long the_last = std::min ( last , (unsigned long) data.numEntries() ) ;
    //
    const bool  weighted = data.isWeighted () ;
    //
    std::vector<Ostap::Math::GSL::P2Quantile> qs ( quantiles.begin() , quantiles.end() ) ;
    //
    unsigned long num = 0 ;
    for ( unsigned long entry = first ; entry < the_last ; ++entry )
    {
      const RooArgSet* vars = data.get( entry ) ;
      if ( nullptr == vars )                              { break    ; } // BREAK 
      //
      if ( cut_range && !vars->allInRange ( cut_range ) ) { continue ; } // CONTINUE    
      // apply cuts:
      const long double wc = nullptr != cuts ? cuts -> getVal() : 1.0L ;
      if ( !wc ) { continue ; }                                          // CONTINUE  
      // apply weight:
      const long double wd = weighted  ? data.weight()   : 1.0L ;
      if ( !wd ) { continue ; }                                          // CONTINUE    
      // cuts & weight:
      const long double w  = wd *  wc ; 
      if ( !w  ) { continue ; }                                          // CONTINUE        
      //
      for ( auto& q : qs ) { q.add ( var.getVal() ) ; }
      ++num ;
    }
    //
    return Ostap::StatVar::Quantiles ( std::vector<double>( qs.begin (), qs.end () )  , num ) ;
  }
  // ==========================================================================
  /** calculate the moment of order "order" relative to the center "center"
   *  @param  tree   (INPUT) input tree 
   *  @param  expr   (INPUT) expression  (must  be valid TFormula!)
   *  @param  cuts   (INPUT) cuts 
   *  @param  order  (INPUT) the order 
   *  @param  center (INPUT) the center 
   *  @param  first  (INPUT) the first  event to process 
   *  @param  last   (INPUT) the last event to  process
   *  @return the moment 
   */
  void _moment_
  ( TTree&                  tree    ,
    Ostap::Math::Statistic& counter ,
    Ostap::Formula&         var     ,
    Ostap::Formula*         cuts    , 
    const unsigned long     first   , 
    const unsigned long     last    )
  {
    // 
    // the loop
    const unsigned long nEntries = std::min ( last , (unsigned long) tree.GetEntries() ) ;
    if ( last <= first  ) { return ; } // RETURN ???    
    //
    Ostap::Utils::Notifier notify ( &tree , &var , cuts ) ;
    const bool with_cuts = nullptr != cuts ? true : false ;
    //
    std::vector<double> results {} ;
    for ( unsigned long entry = first ; entry < nEntries ; ++entry ) 
    {      
      long ievent = tree.GetEntryNumber ( entry ) ;
      if ( 0 > ievent ) { break ; }                        // BREAK
      //
      ievent      = tree.LoadTree ( ievent ) ;
      if ( 0 > ievent ) { break ; }                        // BREAK
      //
      const long double w = with_cuts ? cuts->evaluate() : 1.0L ;
      //
      if  ( !w ) { continue ; }                            // ATTENTION!
      //
      var.evaluate ( results ) ;
      for ( const long double r : results ) { counter.update ( r ) ; } 
    }
    //
  }
  // ==========================================================================
  /** calculate the moment of order "order" relative to the center "center"
   *  @param  tree   (INPUT) input tree 
   *  @param  expr   (INPUT) expression  (must  be valid TFormula!)
   *  @param  cuts   (INPUT) cuts 
   *  @param  order  (INPUT) the order 
   *  @param  center (INPUT) the center 
   *  @param  first  (INPUT) the first  event to process 
   *  @param  last   (INPUT) the last event to  process
   *  @return the moment 
   */
  void _moment_
  ( TTree&                   tree    ,
    Ostap::Math::WStatistic& counter ,
    Ostap::Formula&          var     ,
    Ostap::Formula*          weight  ,
    Ostap::Formula*          cuts    , 
    const unsigned long      first   , 
    const unsigned long      last    )
  {
    // 
    // the loop
    const unsigned long nEntries = std::min ( last , (unsigned long) tree.GetEntries() ) ;
    if ( last <= first  ) { return ; } // RETURN ???    
    //
    Ostap::Utils::Notifier notify ( &tree , &var , weight , cuts ) ;
    const bool with_cuts   = nullptr != cuts   ? true : false ;
    const bool with_weight = nullptr != weight ? true : false ;
    //
    std::vector<double> results {} ;
    for ( unsigned long entry = first ; entry < nEntries ; ++entry ) 
    {      
      long ievent = tree.GetEntryNumber ( entry ) ;
      if ( 0 > ievent ) { break ; }                        // BREAK
      //
      ievent      = tree.LoadTree ( ievent ) ;
      if ( 0 > ievent ) { break ; }                        // BREAK
      //
      const double c = with_cuts   ? cuts  ->evaluate() : 1.0 ;
      if  ( !c ) { continue ; }                            // ATTENTION!
      //
      const double w = with_weight ? weight->evaluate() : 1.0 ;
      var.evaluate ( results ) ;
      for ( const long double r : results ) { counter.update ( r , w  ) ; } 
    }
    //
  }
  // ==========================================================================
} //                                                 end of anonymous namespace
// ============================================================================
/*  check if there is at least one entry that satisfies criteria 
 *  @param tree       (INPUT) the tree 
 *  @param cuts       (INPUT) criteria 
 *  @param first      (INPUT) the first entry 
 *  @param last       (INPUT) the last entry
 *  @return true if there exist at leats one entry 
 */
// ============================================================================
bool Ostap::StatVar::hasEntry 
( TTree*              tree  , 
  const std::string&  cuts  ,   
  const unsigned long first ,
  const unsigned long last  )
{
  // check arguments 
  if ( nullptr == tree || last <= first || tree->GetEntries() < first ) { return false ; }
  //
  Ostap::Formula formula ( cuts , tree ) ;
  if ( !formula.GetNdim() )         { return false  ; }  // RETURN
  //
  Ostap::Utils::Notifier notify ( tree , &formula ) ;
  //
  const unsigned long nEntries =
    std::min ( last , (unsigned long) tree->GetEntries() ) ;
  //
  std::vector<double>  results {} ;
  for ( unsigned long entry = first ; entry < nEntries ; ++entry )
  {
    long ievent = tree->GetEntryNumber ( entry ) ;
    if ( 0 > ievent ) { return false  ; }                // RETURN
    //
    ievent      = tree->LoadTree ( ievent ) ;
    if ( 0 > ievent ) { return false  ; }                // RETURN
    //
    formula.evaluate ( results ) ;
    for  ( const double r : results ) 
    { if ( r ) { return true ; } }
    //
  }
  return false ;
}
// ============================================================================
/** check if there is at least one entry that satisfies criteria 
 *  @param data       (INPUT) data 
 *  @param cuts       (INPUT) criteria 
 *  @param cut_range  (INPUT) cut range 
 *  @param first      (INPUT) the first entry 
 *  @param last       (INPUT) the last entry
 *  @return true if there exist at leats one entry 
 */
// ============================================================================
bool Ostap::StatVar::hasEntry 
( const RooAbsData*   data      , 
  const std::string&  cuts      , 
  const std::string&  cut_range , 
  const unsigned long first     ,
  const unsigned long last      ) 
{
  if ( nullptr == data || last <= first || data->numEntries() < first ) { return false ; }
  //
  const std::unique_ptr<Ostap::FormulaVar> selection { make_formula ( cuts , *data , true ) } ;
  //
  const unsigned long the_last  = std::min ( last , (unsigned long) data->numEntries() ) ;
  //
  const char* cutrange = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  //
  // start the loop
  for ( unsigned long entry = first ; entry < the_last ; ++entry )
  {
    //
    const RooArgSet* vars = data->get( entry ) ;
    if ( nullptr == vars  )                           { break    ; } // RETURN
    if ( cutrange && !vars->allInRange ( cutrange ) ) { continue ; } // CONTINUE    
    //
    // apply cuts:
    const long double wc = selection -> getVal() ;
    if ( wc ) { return true ; }
  }
  //
  return false ;
}
// ============================================================================
/** check if there is at least one entry that satisfies criteria 
 *  @param data       (INPUT) data 
 *  @param cuts       (INPUT) criteria 
 *  @param first      (INPUT) the first entry 
 *  @param last       (INPUT) the last entry
 *  @return true if there exist at leats one entry 
 */
// ============================================================================
bool Ostap::StatVar::hasEntry 
( const RooAbsData*   data      , 
  const std::string&  cuts      , 
  const unsigned long first     ,
  const unsigned long last      ) 
{ return hasEntry ( data , cuts , std::string() , first , last ) ; }
// ============================================================================
/*  build statistic for the <code>expression</code>
 *  @param tree (INPUT) the tree
 *  @param expression (INPUT) the expression
 *
 *  @code
 *  tree = ...
 *  stat = tree.statVar( 'S_sw' )
 *  @endcode
 *
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2013-10-13
 */
// ============================================================================
Ostap::StatEntity
Ostap::StatVar::statVar
( TTree*              tree       ,
  const std::string&  expression ,
  const unsigned long first      ,
  const unsigned long last       )
{
  Ostap::StatEntity result {} ;
  if ( 0 == tree || last <= first ) { return result ; }  // RETURN
  Ostap::Formula formula ( expression , tree ) ;
  if ( !formula.GetNdim() )         { return result ; }  // RETURN
  //
  Ostap::Utils::Notifier notify ( tree , &formula ) ;
  //
  const unsigned long nEntries =
    std::min ( last , (unsigned long) tree->GetEntries() ) ;
  //
  std::vector<double>  results {} ;
  for ( unsigned long entry = first ; entry < nEntries ; ++entry )
  {
    //
    long ievent = tree->GetEntryNumber ( entry ) ;
    if ( 0 > ievent ) { return result ; }                // RETURN
    //
    ievent      = tree->LoadTree ( ievent ) ;
    if ( 0 > ievent ) { return result ; }                // RETURN
    //
    formula.evaluate ( results ) ;
    for  ( const double r : results ) { result += r ; }
  }
  //
  return result ;
}
// ============================================================================
/*  build statistic for the <code>expression</code>
 *  @param tree       (INPUT) the tree
 *  @param expression (INPUT) the expression
 *  @param cuts       (INPUT) the selection criteria
 *
 *  @code
 *  tree = ...
 *  stat = tree.statVar( 'S_sw' ,'pt>1000')
 *  @endcode
 *
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2013-10-13
 */
// ============================================================================
Ostap::WStatEntity
Ostap::StatVar::statVar
( TTree*              tree       ,
  const std::string&  expression ,
  const std::string&  cuts       ,
  const unsigned long first      ,
  const unsigned long last       )
{
  //
  if ( cuts.empty() ) { return statVar( tree , expression , first , last ) ; }
  //
  Ostap::WStatEntity result ;
  if ( 0 == tree || last <= first ) { return result ; }  // RETURN
  Ostap::Formula selection ( cuts      , tree ) ;
  if ( !selection.ok () ) { return result ; }            // RETURN
  Ostap::Formula formula   ( expression , tree ) ;
  if ( !formula  .ok () ) { return result ; }            // RETURN
  //
  Ostap::Utils::Notifier notify ( tree , &selection,  &formula ) ;
  //
  const unsigned long nEntries =
    std::min ( last , (unsigned long) tree->GetEntries() ) ;
  //
  std::vector<double> results {} ;
  for ( unsigned long entry = first ; entry < nEntries ; ++entry )
  {
    //
    long ievent = tree->GetEntryNumber ( entry ) ;
    if ( 0 > ievent ) { return result ; }                // RETURN
    //
    ievent      = tree->LoadTree ( ievent ) ;
    if ( 0 > ievent ) { return result ; }                // RETURN
    //
    const double w = selection.evaluate() ;
    //
    if  ( !w ) { continue ; }                            // ATTENTION!
    //
    formula.evaluate ( results ) ;
    for  ( const double r : results ) { result.add (  r , w ) ; }
    //
  }
  //
  return result ;
}
// ============================================================================
/*  build statistic for the <code>expressions</code>
 *  @param tree        (INPUT)  the tree 
 *  @param result      (UPDATE) the output statistics for specified expressions 
 *  @param expressions (INPUT)  the list of  expressions
 *  @param first       (INPUT)  the first entry to process 
 *  @param last        (INPUT)  the last entry to process (not including!)
 *  @return number of processed entries 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2018-11-04
 */
// ============================================================================
unsigned long Ostap::StatVar::statVars
( TTree*                                  tree        ,  
  Ostap::StatVar::WStatVector& result      ,  
  const Ostap::Strings&            expressions ,
  const unsigned long                     first       ,
  const unsigned long                     last        ) 
{
  //
  const unsigned int N = expressions.size() ;
  //
  result.resize ( N ) ;
  for ( auto& r : result ) { r.reset () ; }
  //
  if ( 0 == tree || last <= first ) { return 0 ; }  // RETURN
  if ( expressions.empty()        ) { return 0 ; }  // RETURN  
  //
  typedef std::unique_ptr<Ostap::Formula> UOF ;
  std::vector<UOF> formulas ; formulas.reserve ( N ) ;
  //
  for ( const auto& e : expressions  ) 
  {
    auto p = std::make_unique<Ostap::Formula>( e , tree ) ;
    if ( !p || !p->ok() ) { return 0 ; }
    formulas.push_back ( std::move ( p ) ) ;  
  }
  //
  Ostap::Assert ( N == formulas.size()           , 
                  "Inconsistent size of structures" , 
                  "Ostap::StatVar::statVars"        ) ;
  //
  Ostap::Utils::Notifier notify ( formulas.begin() , formulas.end() , tree ) ;
  //
  const unsigned long nEntries =
    std::min ( last , (unsigned long) tree->GetEntries() ) ;
  //
  std::vector<double>  results {} ;
  for ( unsigned long entry = first ; entry < nEntries ; ++entry )
  {
    //
    long ievent = tree->GetEntryNumber ( entry ) ;
    if ( 0 > ievent ) { return entry - first ; }                // RETURN
    //
    ievent      = tree->LoadTree ( ievent ) ;
    if ( 0 > ievent ) { return entry - first  ; }                // RETURN
    //
    for ( unsigned int i = 0 ; i < N ; ++i ) 
    {
      formulas[i]->evaluate ( results ) ;
      for ( const double r : results ) { result[i] += r ; }
    }
  }
  //
  return results.empty() ? 0 : result[0].nEntries() ;
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
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2018-11-04
 */
// ============================================================================
unsigned long Ostap::StatVar::statVars
( TTree*                                  tree        ,  
  Ostap::StatVar::WStatVector&            result      ,  
  const Ostap::Strings&                   expressions ,
  const std::string&                      cuts        ,
  const unsigned long                     first       ,
  const unsigned long                     last        ) 
{
  //
  if ( cuts.empty() ) { return statVars ( tree , result , expressions , first , last ) ; }
  //
  const unsigned int N = expressions.size() ;
  //
  result.resize ( N ) ; 
  for ( auto& r : result ) { r.reset () ; }
  //
  if ( 0 == tree || last <= first ) { return 0 ; }  // RETURN
  if ( expressions.empty()        ) { return 0 ; }  // RETURN  
  //
  Ostap::Formula selection ( cuts , tree ) ;
  if ( !selection .ok ()          ) { return 0 ; }  // RETURN
  //
  typedef std::unique_ptr<Ostap::Formula> UOF ;
  std::vector<UOF> formulas ; formulas.reserve ( N ) ;
  for ( const auto& e : expressions  ) 
  {
    auto p = std::make_unique<Ostap::Formula>( e , tree ) ;
    if ( !p || !p->ok() ) { return 0 ; }
    formulas.push_back ( std::move ( p ) ) ;
  }
  //
  Ostap::Assert ( N == formulas.size()              , 
                  "Inconsistent size of structures" , 
                  "Ostap::StatVar::statVars"        ) ;
  //
  Ostap::Utils::Notifier notify ( formulas.begin() , formulas.end() , &selection , tree ) ;
  //
  const unsigned long nEntries =
    std::min ( last , (unsigned long) tree->GetEntries() ) ;
  //
  std::vector<double>  results {} ;
  for ( unsigned long entry = first ; entry < nEntries ; ++entry )
  {
    //
    long ievent = tree->GetEntryNumber ( entry ) ;
    if ( 0 > ievent ) { return entry - first ; }                // RETURN
    //
    ievent      = tree->LoadTree ( ievent ) ;
    if ( 0 > ievent ) { return entry - first ; }                // RETURN
    //
    const double w = selection.evaluate() ;
    if ( !w ) { continue  ; }
    //
    for ( unsigned int i = 0 ; i < N ; ++i ) 
    {
      formulas[i]->evaluate ( results ) ;
      for ( const double r : results ) { result[i].add ( r , w ) ; }
    }
  }
  //
  return results.empty() ? 0 : result[0].nEntries() ;
}
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
( TTree*                     tree    ,
  const std::string&         exp1    ,
  const std::string&         exp2    ,
  const unsigned long        first   ,
  const unsigned long        last    )
{
  //
  Ostap::Assert ( nullptr != tree           ,  
		  "Invalid TTree"           , 
		  "Ostap::StatVar::statCov" ) ;
  // 
  // prepare the result 
  Ostap::Math::Covariance result {} ;
  ///
  if ( last <= first ) { return result ; }  // RETURN
  
  Ostap::Formula formula1 ( exp1 , tree ) ;
  Ostap::Assert ( formula1.ok()                                        ,
		  std::string ( "Invalid first  expression: " ) + exp1  ,  
		  "Ostap::StatVar::statCov"                            ) ;
  Ostap::Formula formula2 ( exp2 , tree ) ;
  Ostap::Assert ( formula2.ok()                                        ,
		  std::string ( "Invalid second expression: " ) + exp2  , 
		  "Ostap::StatVar::statCov"                            ) ;
  //
  Ostap::Utils::Notifier notify ( tree , &formula1 , &formula2 ) ;
  //
  const unsigned long nEntries = std::min ( last , (unsigned long) tree->GetEntries() ) ;
  //
  std::vector<double> results1 {} ;
  std::vector<double> results2 {} ;
  //
  for ( unsigned long entry = first ; entry < nEntries ; ++entry )
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
Ostap::Math::WCovariance
Ostap::StatVar::statCov
( TTree*                     tree    ,
  const std::string&         exp1    ,
  const std::string&         exp2    ,
  const std::string&         cuts    ,
  const unsigned long        first   ,
  const unsigned long        last    )
{
  //
  Ostap::Assert ( nullptr != tree           ,  
		  "Invalid TTree"           , 
		  "Ostap::StatVar::statCov" ) ;
  // 
  // prepare the result 
  Ostap::Math::WCovariance result {} ;
  ///
  if ( last <= first ) { return result ; }
  //
  Ostap::Formula formula1 ( exp1 , tree ) ;
  Ostap::Assert ( formula1.ok()                        ,
		  std::string ( "Invalid first  expression: " ) + exp1  ,  
		  "Ostap::StatVar::statCov"            ) ;
  //
  Ostap::Formula formula2 ( exp2 , tree ) ;
  Ostap::Assert ( formula2.ok()                        ,
		  std::string ( "Invalid second expression: " ) + exp2  , 
		  "Ostap::StatVar::statCov"            ) ;
  //
  const std::unique_ptr<Ostap::Formula> selection { !cuts.empty() ? new Ostap::Formula ( cuts , tree ) : nullptr } ;
  Ostap::Assert ( !selection || selection->ok ()       ,
		  std::string ( "Invalid selection/weight: " ) + cuts , 
		  "Ostap::StatVar::statCov"            ) ;
  ///
  Ostap::Utils::Notifier notify ( tree , &formula1 , selection.get() ) ;
  //
  const bool with_cuts = selection && selection->ok () ;
  //
  const unsigned long nEntries = std::min ( last , (unsigned long) tree->GetEntries() ) ;
  //
  std::vector<double> results1 {} ;
  std::vector<double> results2 {} ;
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
      formula1.evaluate ( results1 ) ;
      formula2.evaluate ( results2 ) ;
      //
      for ( const long double v1 : results1 ) 
	{ 
	  for ( const long double v2 : results2 ) 
	    {
	      result.add ( v1 , v2 , w ) ;
	    }
	}
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
/*  calculate the covariance of two expressions 
 *  @param tree  (INPUT)  the inpout tree 
 *  @param vars  (INPUT)  expressions 
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
  Ostap::StatVar::WStatVector&    stats ,  
  TMatrixTSym<double>&            cov2  , 
  const unsigned long             first ,
  const unsigned long             last  ) 
{
  static const std::string cuts{} ;
  return statCov ( tree , vars ,  cuts , stats , cov2 , first , last  ) ;
}
// ============================================================================
Ostap::WStatEntity
Ostap::StatVar::statVar
( const RooAbsData*   data        ,
  const std::string&  expression  ,
  const std::string&  cuts        ,
  const std::string&  cut_range   , 
  const unsigned long first       ,
  const unsigned long last        )
{
    
  // ==========================================================================
  Ostap::Assert ( nullptr != data           ,  
		  "Invalid RotAbsData"      , 
		  "Ostap::StatVar::statVar" ) ;
  // ==========================================================================
  Ostap::WStatEntity result{} ;
  if ( 0 == data || last <= first ) { return result ; }         // RETURN
  //
  const std::unique_ptr<Ostap::FormulaVar> formula   { make_formula ( expression , *data        ) } ;
  const std::unique_ptr<Ostap::FormulaVar> selection { make_formula ( cuts       , *data , true ) } ;
  //
  const bool  weighted = data->isWeighted() ;
  //
  const unsigned long the_last  = std::min ( last , (unsigned long) data->numEntries() ) ;
  //
  const char* cutrange = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  //
  // start the loop
  for ( unsigned long entry = first ; entry < the_last ; ++entry )
  {
    //
    const RooArgSet* vars = data->get( entry ) ;
    if ( nullptr == vars  )                           { break    ; } // RETURN
    if ( cutrange && !vars->allInRange ( cutrange ) ) { continue ; } // CONTINUE
    //
    // apply cuts:
    const long double wc = selection ? selection -> getVal() : 1.0L ;
    if ( !wc ) { continue ; }                                   // CONTINUE  
    // apply weight:
    const long double wd = weighted  ? data->weight()        : 1.0L ;
    if ( !wd ) { continue ; }                                   // CONTINUE    
    // cuts & weight:
    const long double w  = wd *  wc ;
    if ( !w  ) { continue ; }                                   // CONTINUE        
    //
    const double v = formula->getVal () ;
    //      
    result.add ( v , w ) ;
  }
  //
  return result ;
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
unsigned long
Ostap::StatVar::statVars
( const RooAbsData*               data        , 
  Ostap::StatVar::WStatVector&    result      , 
  const Ostap::Strings&           expressions ,
  const std::string&              cuts        ,
  const std::string&              cut_range   , 
  const unsigned long             first       ,
  const unsigned long             last        ) 
{
  // 
  const unsigned int N = expressions.size() ;
  //
  result.resize ( N ) ; 
  for ( auto& r : result ) { r.reset () ; }
  //
  if ( expressions.empty()              ) { return 0 ; }
  if ( nullptr == data || last <= first ) { return 0 ; }
  if ( data->numEntries() <= first      ) { return 0 ; }
  //
  const std::unique_ptr<Ostap::FormulaVar> selection { make_formula ( cuts       , *data , true ) } ;
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
  const unsigned long the_last  = std::min ( last , (unsigned long) data->numEntries() ) ;
  //
  // start the loop
  for ( unsigned long entry = first ; entry < the_last ; ++entry )
  {
    //
    const RooArgSet* vars = data->get( entry ) ;
    if ( nullptr == vars  )                           { break    ; } // RETURN
    if ( cutrange && !vars->allInRange ( cutrange ) ) { continue ; } // CONTINUE    
    //
    // apply cuts:
    const long double wc = selection ? selection -> getVal() : 1.0L ;
    if ( !wc ) { continue ; }                                   // CONTINUE  
    // apply weight:
    const long double wd = weighted  ? data->weight()        : 1.0L ;
    if ( !wd ) { continue ; }                                   // CONTINUE    
    // cuts & weight:
    const long double w  = wd *  wc ;
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
  return result.empty() ? 0 : result[0].nEntries() ;
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
( const RooAbsData*    data      , 
  const std::string&   exp1      , 
  const std::string&   exp2      , 
  const std::string&   cuts      , 
  const std::string&   cut_range , 
  const unsigned long  first     ,
  const unsigned long  last      )
{
  // ===========================================================================
  Ostap::Assert ( nullptr != data           ,  
		  "Invalid RotAbsData"      , 
		  "Ostap::StatVar::statCov" ) ;
  //
  // prepare the result 
  Ostap::Math::WCovariance result {} ;
  //
  if ( last <= first ) { return result ; }         // RETURN
  //
  const std::unique_ptr<Ostap::FormulaVar> formula1  { make_formula ( exp1 , *data        ) } ;
  const std::unique_ptr<Ostap::FormulaVar> formula2  { make_formula ( exp2 , *data        ) } ;
  const std::unique_ptr<Ostap::FormulaVar> selection { make_formula ( cuts , *data , true ) } ;
  //
  const char* cutrange = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  const bool  weighted = data->isWeighted() ;
  //
  const unsigned long the_last  = std::min ( last , (unsigned long) data->numEntries() ) ;
  //
  // start the loop
  for ( unsigned long entry = first ; entry < the_last ; ++entry )
    {
      //
      const RooArgSet* vars = data->get( entry ) ;
      if ( nullptr == vars  )                           { break    ; } // RETURN
      if ( cutrange && !vars->allInRange ( cutrange ) ) { continue ; } // CONTINUE    
      //
      // apply cuts:
      const long double wc = selection ? selection -> getVal() : 1.0L ;
      if ( !wc ) { continue ; }                                   // CONTINUE  
      // apply weight:
      const long double wd = weighted  ? data->weight()        : 1.0L ;
      if ( !wd ) { continue ; }                                   // CONTINUE    
      // cuts & weight:
      const long double w  = wd *  wc ;
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
/*  calculate the covariance of several expressions 
 *  @param tree      (INPUT)  the inpout tree 
 *  @param vars      (INPUT)  expressions 
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
  Ostap::StatVar::WStatVector&    stats     ,  
  TMatrixTSym<double>&            cov2      ,  
  const std::string&              cut_range , 
  const unsigned long             first     ,
  const unsigned long             last      )
{
  static const std::string _cuts{} ;
  return statCov ( data , vars , _cuts , stats , cov2 , cut_range , first , last ) ;  
}
// ============================================================================
/*  get the number of equivalent entries 
 *  \f$ n_{eff} \equiv = \frac{ (\sum w)^2}{ \sum w^2} \f$
 *  @param tree  (INPUT) the tree 
 *  @param cuts  (INPUT) selection  criteria 
 *  @param  first  (INPUT) the first  event to process 
 *  @param  last   (INPUT) the last event to  process
 *  @return number of equivalent entries 
 */
// ========================================================================
double Ostap::StatVar::nEff 
( TTree&               tree  , 
  const std::string&   cuts  , 
  const unsigned long  first ,
  const unsigned long  last  ) 
{
  //
  std::unique_ptr<Ostap::Formula> cut { nullptr } ;
  if  ( !cuts.empty() ) 
  { 
    cut = std::make_unique<Ostap::Formula>( cuts , &tree ) ; 
    Ostap::Assert ( cut && cut->ok()               , 
                    "Invalid cut:\"" + cuts + '\"' ,
                    "Ostap::StatVar::nEff"         ) ;
  }
  //
  return  _neff_ ( tree , cut.get() , first , last ) ;
}
// ============================================================================
/*  calculate the moment of order "order" relative to the center "center"
 *  @param  tree   (INPUT) input tree 
 *  @param  expr   (INPUT) expression  (must  be valid TFormula!)
 *  @param  order  (INPUT) the order 
 *  @param  center (INPUT) the center 
 *  @param  cuts   (INPUT) cuts 
 *  @param  first  (INPUT) the first  event to process 
 *  @param  last   (INPUT) the last event to  process
 *  @return the moment 
 */ 
// ============================================================================
double Ostap::StatVar::get_moment 
( TTree&               tree   ,  
  const unsigned short order  , 
  const std::string&   expr   , 
  const double         center ,
  const std::string&   cuts   , 
  const unsigned long  first  ,
  const unsigned long  last   ) 
{
  //
  if ( 0 == order ){ return 1 ; } // RETURN 
  //
  Ostap::Formula var ( expr , &tree ) ;
  Ostap::Assert ( var.ok()                            , 
                  "Invalid expression:'" + expr + "'" ,
                  "Ostap::StatVar::moment"            ) ;
  //
  std::unique_ptr<Ostap::Formula> cut { nullptr } ;
  if  ( !cuts.empty() ) 
  { 
    cut = std::make_unique<Ostap::Formula>( cuts , &tree ) ; 
    Ostap::Assert ( cut && cut->ok()               , 
                    "Invalid cut:\"" + cuts + "\"" ,
                    "Ostap::StatVar::moment"       ) ;
  }
  //
  return _moment1_ ( tree , var , cut.get() , order , center , first , last ) ;
}
// ============================================================================
/*  calculate the moment of order "order"
 *  @param  tree  (INPUT) input tree 
 *  @param  order (INPUT) the order 
 *  @param  expr  (INPUT) expression  (must  be valid TFormula!)
 *  @param  cuts  (INPUT) cuts 
 *  @param  error (INPUT) calculate the uncertainty?
 *  @param  first (INPUT) the first  event to process 
 *  @param  last  (INPUT) the last event to  process
 *  @return the moment 
 */ 
// ============================================================================
Ostap::Math::ValueWithError
Ostap::StatVar::moment
( TTree&               tree  ,  
  const unsigned short order ,
  const std::string&   expr  , 
  const std::string&   cuts  , 
  const unsigned long  first ,
  const unsigned long  last  ) 
{
  //
  if ( 0 == order ){ return 1 ; } // RETURN 
  //
  Ostap::Formula var ( expr , &tree ) ;
  Ostap::Assert ( var.ok()                              ,
                  "Invalid expression:\"" + expr + "\"" ,
                  "Ostap::StatVar::moment"              ) ;  
  //
  std::unique_ptr<Ostap::Formula> cut { nullptr } ;
  if  ( !cuts.empty() ) 
  { 
    cut = std::make_unique<Ostap::Formula>( cuts , &tree ) ; 
    Ostap::Assert ( cut && cut->ok()               ,   
                    "Invalid cut:\"" + cuts + "\"" ,
                    "Ostap::StatVar::moment"       ) ;
  }
  //
  return _moment2_ ( tree , order , var , cut.get() , first , last ) ;
}
// ============================================================================
/* calculate the central moment of order "order"
 *  @param  tree  (INPUT) input tree 
 *  @param  order (INPUT) the order 
 *  @param  expr  (INPUT) expression  (must  be valid TFormula!)
 *  @param  cuts  (INPUT) cuts 
 *  @param  error (INPUT) calculate the uncertainty?
 *  @param  first (INPUT) the first  event to process 
 *  @param  last  (INPUT) the last event to  process
 *  @return the moment 
 */ 
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::StatVar::central_moment 
( TTree&               tree  ,  
  const unsigned short order , 
  const std::string&   expr  , 
  const std::string&   cuts  , 
  const unsigned long  first ,
  const unsigned long  last  ) 
{
  //
  if      ( 0 == order ){ return 1 ; } // RETURN 
  else if ( 1 == order ){ return 0 ; } // RETURN 
  //
  Ostap::Formula var ( expr , &tree ) ;
  Ostap::Assert ( var.ok()                              , 
                  "Invalid expression:\"" + expr + "\"" ,
                  "Ostap::StatVar::central_moment"      ) ;
  //
  std::unique_ptr<Ostap::Formula> cut { nullptr } ;
  if  ( !cuts.empty() ) 
  { 
    cut = std::make_unique<Ostap::Formula>( cuts , &tree ) ; 
    Ostap::Assert ( cut && cut->ok()                 , 
                    "Invalid cut:\"" + cuts + "\""   ,
                    "Ostap::StatVar::central_moment" ) ;
  }
  //
  return _moment3_ ( tree , order , var , cut.get() , first , last ) ;
}
// ============================================================================
/*  calculate the skewness of the  distribution
 *  @param  tree  (INPUT) input tree 
 *  @param  expr  (INPUT) expression  (must  be valid TFormula!)
 *  @param  cuts  (INPUT) cuts 
 *  @param  first (INPUT) the first  event to process 
 *  @param  last  (INPUT) the last event to  process
 *  @return the skewness 
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::StatVar::skewness 
( TTree&               tree  ,  
  const std::string&   expr  , 
  const std::string&   cuts  , 
  const unsigned long  first ,
  const unsigned long  last  ) 
{
  //
  Ostap::Formula var ( expr , &tree ) ;
  Ostap::Assert ( var.ok()                              , 
                  "Invalid expression:\"" + expr + "\"" ,
                  "Ostap::StatVar::skewness"            ) ;
  //
  std::unique_ptr<Ostap::Formula> cut { nullptr } ;
  if  ( !cuts.empty() ) 
  { 
    cut = std::make_unique<Ostap::Formula>( cuts , &tree ) ; 
    Ostap::Assert ( cut && cut->ok()               , 
                    "Invalid cut:\"" + cuts + "\"" ,
                    "Ostap::StatVar::skewness"     ) ;
  }
  //
  return _skewness_ ( tree , var , cut.get() , first , last ) ;
}
// ============================================================================
/*  calculate the (excess) kurtosis of the  distribution
 *  @param  tree  (INPUT) input tree 
 *  @param  expr  (INPUT) expression  (must  be valid TFormula!)
 *  @param  cuts  (INPUT) cuts 
 *  @param  first (INPUT) the first  event to process 
 *  @param  last  (INPUT) the last event to  process
 *  @return the (excess) kurtosis
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::StatVar::kurtosis
( TTree&               tree  ,  
  const std::string&   expr  , 
  const std::string&   cuts  , 
  const unsigned long  first ,
  const unsigned long  last  ) 
{
  //
  Ostap::Formula var ( expr , &tree ) ;
  Ostap::Assert ( var.ok()                             , 
                  "Invalid expression:\"" + expr + "\"" , 
                  "Ostap::StatVar::kurtosis"            ) ;  
  //
  std::unique_ptr<Ostap::Formula> cut { nullptr } ;
  if  ( !cuts.empty() ) 
  { 
    cut = std::make_unique<Ostap::Formula>( "", cuts , &tree ) ; 
    Ostap::Assert ( cut && cut->ok()               ,  
                    "Invalid cut:\"" + cuts + "\"" ,
                    "Ostap::StatVar::kurtosis"     ) ;
  }
  //
  return _kurtosis_ ( tree , var , cut.get() , first , last ) ;
} 
// ============================================================================
/*   get quantile of the distribution  
 *   @param tree  (INPUT) the input tree 
 *   @param q     (INPUT) quantile value   0 < q < 1  
 *   @param expr  (INPUT) the expression 
 *   @param cuts  (INPUT) selection cuts 
 *   @param  first (INPUT) the first  event to process 
 *   @param  last  (INPUT) the last event to  process
 *   @return the quantile value 
 */
// ============================================================================
Ostap::StatVar::Quantile
Ostap::StatVar::quantile
( TTree&              tree  ,
  const double        q     , //  0<q<1 
  const std::string&  expr  , 
  const std::string&  cuts  , 
  const unsigned long first ,
  const unsigned long last  ) 
{
  Ostap::Assert ( 0 < q && q < 1             , 
                  "Invalid quantile"         ,
                  "Ostap::StatVar::quantile" ) ;
  //
  Ostap::Formula var ( expr , &tree ) ;
  Ostap::Assert ( var.ok()                              , 
                  "Invalid expression:\"" + expr + "\"" ,
                  "Ostap::StatVar::quantile"            ) ;
  //
  std::unique_ptr<Ostap::Formula> cut { nullptr } ;
  if  ( !cuts.empty() ) 
  { 
    cut = std::make_unique<Ostap::Formula>( "", cuts , &tree ) ; 
    Ostap::Assert ( cut && cut->ok()               ,
                    "Invalid cut:\"" + cuts + "\"" ,
                    "Ostap::StatVar::quantile"     ) ;
  }
  //
  std::set<double> qset {} ;
  qset.insert ( q )  ;
  //
  auto result = _quantiles_ ( tree , qset , var , cut.get() , first , last ) ; 
  //
  Ostap::Assert ( 1 == result.quantiles.size()         , 
                  "Invalid quantiles size"   ,
                  "Ostap::StatVar::interval" ) ;
  //
  return Quantile ( result.quantiles[0] , result.nevents ) ;
}
// ============================================================================
/*   get approximate  quantile of the distribution  using p^2 algortihm
 *   @param tree  (INPUT) the input tree 
 *   @param q     (INPUT) quantile value   0 < q < 1  
 *   @param expr  (INPUT) the expression 
 *   @param cuts  (INPUT) selection cuts 
 *   @param  first (INPUT) the first  event to process 
 *   @param  last  (INPUT) the last event to  process
 *   @return the quantile value 
 */
// ============================================================================
Ostap::StatVar::Quantile
Ostap::StatVar::p2quantile
( TTree&              tree  ,
  const double        q     , //  0<q<1 
  const std::string&  expr  , 
  const std::string&  cuts  , 
  const unsigned long first ,
  const unsigned long last  ) 
{
  Ostap::Assert ( 0 < q && q < 1             , 
                  "Invalid quantile"         ,
                  "Ostap::StatVar::quantile" ) ;
  //
  Ostap::Formula var ( expr , &tree ) ;
  Ostap::Assert ( var.ok()                              , 
                  "Invalid expression:\"" + expr + "\"" ,
                  "Ostap::StatVar::quantile"            ) ;
  //
  std::unique_ptr<Ostap::Formula> cut { nullptr } ;
  if  ( !cuts.empty() ) 
  { 
    cut = std::make_unique<Ostap::Formula>( "", cuts , &tree ) ; 
    Ostap::Assert ( cut && cut->ok()               ,
                    "Invalid cut:\"" + cuts + "\"" ,
                    "Ostap::StatVar::quantile"     ) ;
  }
  //
  std::set<double> qset {} ;
  qset.insert ( q )  ;
  //
  auto result = _p2quantiles_ ( tree , qset , var , cut.get() , first , last ) ; 
  //
  Ostap::Assert ( 1 == result.quantiles.size()         , 
                  "Invalid quantiles size"   ,
                  "Ostap::StatVar::interval" ) ;
  //
  return Quantile ( result.quantiles[0] , result.nevents ) ;
}
// ============================================================================


// ============================================================================
/*   get quantiles of the distribution  
 *   @param tree  (INPUT) the input tree 
 *   @param q     (INPUT) quantile value   0 < q < 1  
 *   @param expr  (INPUT) the expression 
 *   @param cuts  (INPUT) selection cuts 
 *   @param  first (INPUT) the first  event to process 
 *   @param  last  (INPUT) the last event to  process
 *   @return the quantile value 
 */
// ============================================================================
Ostap::StatVar::Quantiles
Ostap::StatVar::quantiles
( TTree&                     tree      ,
  const std::vector<double>& quantiles , 
  const std::string&         expr      , 
  const std::string&         cuts      , 
  const unsigned long        first     ,
  const unsigned long        last      ) 
{
  //
  std::set<double> qs ;
  for ( double v : quantiles ) { qs.insert ( v ) ; }
  Ostap::Assert ( !qs.empty ()                ,
                  "Invalid quantiles"         ,
                  "Ostap::StatVar::quantiles" ) ;
  Ostap::Assert ( 0 < *qs. begin ()           , 
                  "Invalid quantile"          ,
                  "Ostap::StatVar::quantiles" ) ;  
  Ostap::Assert ( 1 > *qs.rbegin ()           , 
                  "Invalid quantile"          ,
                  "Ostap::StatVar::quantiles" ) ;
  //
  Ostap::Formula var ( expr , &tree ) ;
  Ostap::Assert ( var.ok()                              ,
                  "Invalid expression:\"" + expr + "\"" ,
                  "Ostap::StatVar::quantile"            ) ;
  //
  std::unique_ptr<Ostap::Formula> cut { nullptr } ;
  if  ( !cuts.empty() ) 
  { 
    cut = std::make_unique<Ostap::Formula>( cuts , &tree ) ; 
    Ostap::Assert ( cut && cut->ok()               , 
                    "Invalid cut:\"" + cuts + "\"" ,
                    "Ostap::StatVar::quantile"     ) ;
  }
  //
  return _quantiles_ ( tree , qs , var  ,  cut.get() , first , last ) ; 
}
// ============================================================================
/*   get (approximate) quantiles of the distribution using P^2 algorithm   
 *   @param tree  (INPUT) the input tree 
 *   @param q     (INPUT) quantile value   0 < q < 1  
 *   @param expr  (INPUT) the expression 
 *   @param cuts  (INPUT) selection cuts 
 *   @param  first (INPUT) the first  event to process 
 *   @param  last  (INPUT) the last event to  process
 *   @return the quantile value 
 */
// ============================================================================
Ostap::StatVar::Quantiles
Ostap::StatVar::p2quantiles
( TTree&                     tree      ,
  const std::vector<double>& quantiles , 
  const std::string&         expr      , 
  const std::string&         cuts      , 
  const unsigned long        first     ,
  const unsigned long        last      ) 
{
  //
  std::set<double> qs ;
  for ( double v : quantiles ) { qs.insert ( v ) ; }
  Ostap::Assert ( !qs.empty ()                ,
                  "Invalid quantiles"         ,
                  "Ostap::StatVar::quantiles" ) ;
  Ostap::Assert ( 0 < *qs. begin ()           , 
                  "Invalid quantile"          ,
                  "Ostap::StatVar::quantiles" ) ;  
  Ostap::Assert ( 1 > *qs.rbegin ()           , 
                  "Invalid quantile"          ,
                  "Ostap::StatVar::quantiles" ) ;
  //
  Ostap::Formula var ( expr , &tree ) ;
  Ostap::Assert ( var.ok()                              ,
                  "Invalid expression:\"" + expr + "\"" ,
                  "Ostap::StatVar::quantile"            ) ;
  //
  std::unique_ptr<Ostap::Formula> cut { nullptr } ;
  if  ( !cuts.empty() ) 
  { 
    cut = std::make_unique<Ostap::Formula>( cuts , &tree ) ; 
    Ostap::Assert ( cut && cut->ok()               , 
                    "Invalid cut:\"" + cuts + "\"" ,
                    "Ostap::StatVar::quantile"     ) ;
  }
  //
  return _p2quantiles_ ( tree , qs , var  ,  cut.get() , first , last ) ; 
}
// ============================================================================
/*  get the interval of the distribution  
 *   @param tree  (INPUT) the input tree 
 *   @param q1    (INPUT) quantile value   0 < q1 < 1  
 *   @param q2    (INPUT) quantile value   0 < q2 < 1  
 *   @param expr  (INPUT) the expression 
 *   @param cuts  (INPUT) selection cuts 
 *   @param  first (INPUT) the first  event to process 
 *   @param  last  (INPUT) the last event to  process
 *   @return the quantile value 
 *   @code
 *   Tree& tree = ... ;
 *   /// get 90% interval:
 *   Interval ab = interval ( tree , 0.05 , 0.95 , 'mass' , 'pt>3' ) ;
 *   @code 
 */
// ============================================================================
Ostap::StatVar::QInterval 
Ostap::StatVar::interval 
( TTree&              tree  ,
  const double        q1    , //  0<q1<1 
  const double        q2    , //  0<q2<1 
  const std::string&  expr  , 
  const std::string&  cuts  , 
  const unsigned long first ,
  const unsigned long last  ) 
{
  Ostap::Assert ( 0 < q1 && q1 < 1 , 
                  "Invalid quantile1"        ,
                  "Ostap::StatVar::interval" ) ;
  Ostap::Assert ( 0 < q2 && q2 < 1 , 
                  "Invalid quantile2"        ,
                  "Ostap::StatVar::interval" ) ;
  //
  Ostap::Formula var ( expr , &tree ) ;
  Ostap::Assert ( var.ok()                              , 
                  "Invalid expression:\"" + expr + "\"" ,
                  "Ostap::StatVar::interval"            ) ;
  //
  std::unique_ptr<Ostap::Formula> cut { nullptr } ;
  if  ( !cuts.empty() ) 
  { 
    cut = std::make_unique<Ostap::Formula>( "", cuts , &tree ) ; 
    Ostap::Assert ( cut && cut->ok()               ,
                    "Invalid cut:\"" + cuts + "\"" ,
                    "Ostap::StatVar::interval"     ) ;
  }
  //
  auto result = 
    _quantiles_ ( tree , std::set<double>{{ q1 , q2 }} , 
                  var  ,  cut.get() , first , last ) ; 
  Ostap::Assert ( 2 == result.quantiles.size()         ,
                  "Invalid interval"         ,
                  "Ostap::StatVar::interval" ) ;
  //
  return QInterval ( Interval ( result.quantiles[0] , result.quantiles[1] ) , result.nevents ) ;
}
// ============================================================================
/*  get the approximate  interval of the distribution  usnig P^2 algorithm
 *   @param tree  (INPUT) the input tree 
 *   @param q1    (INPUT) quantile value   0 < q1 < 1  
 *   @param q2    (INPUT) quantile value   0 < q2 < 1  
 *   @param expr  (INPUT) the expression 
 *   @param cuts  (INPUT) selection cuts 
 *   @param  first (INPUT) the first  event to process 
 *   @param  last  (INPUT) the last event to  process
 *   @return the quantile value 
 *   @code
 *   Tree& tree = ... ;
 *   /// get 90% interval:
 *   Interval ab = p2interval ( tree , 0.05 , 0.95 , 'mass' , 'pt>3' ) ;
 *   @code 
 */
// ============================================================================
Ostap::StatVar::QInterval 
Ostap::StatVar::p2interval 
( TTree&              tree  ,
  const double        q1    , //  0<q1<1 
  const double        q2    , //  0<q2<1 
  const std::string&  expr  , 
  const std::string&  cuts  , 
  const unsigned long first ,
  const unsigned long last  ) 
{
  Ostap::Assert ( 0 < q1 && q1 < 1 , 
                  "Invalid quantile1"        ,
                  "Ostap::StatVar::interval" ) ;
  Ostap::Assert ( 0 < q2 && q2 < 1 , 
                  "Invalid quantile2"        ,
                  "Ostap::StatVar::interval" ) ;
  //
  Ostap::Formula var ( expr , &tree ) ;
  Ostap::Assert ( var.ok()                              , 
                  "Invalid expression:\"" + expr + "\"" ,
                  "Ostap::StatVar::interval"            ) ;
  //
  std::unique_ptr<Ostap::Formula> cut { nullptr } ;
  if  ( !cuts.empty() ) 
  { 
    cut = std::make_unique<Ostap::Formula>( "", cuts , &tree ) ; 
    Ostap::Assert ( cut && cut->ok()               ,
                    "Invalid cut:\"" + cuts + "\"" ,
                    "Ostap::StatVar::interval"     ) ;
  }
  //
  auto result = 
    _p2quantiles_ ( tree , std::set<double>{{ q1 , q2 }} , 
                    var  ,  cut.get() , first , last ) ; 
  Ostap::Assert ( 2 == result.quantiles.size()         ,
                  "Invalid interval"         ,
                  "Ostap::StatVar::interval" ) ;
  //
  return QInterval ( Interval ( result.quantiles[0] , result.quantiles[1] ) , result.nevents ) ;
}
// ============================================================================
/** get the number of equivalent entries 
 *  \f$ n_{eff} \equiv = \frac{ (\sum w)^2}{ \sum w^2} \f$
 *  @param tree  (INPUT) the tree 
 *  @param cuts  (INPUT) selection  criteria 
 *  @param  first  (INPUT) the first  event to process 
 *  @param  last   (INPUT) the last event to  process
 *  @return number of equivalent entries 
 */
// ============================================================================
double Ostap::StatVar::nEff 
( const RooAbsData&    data      , 
  const std::string&   cuts      , 
  const std::string&   cut_range , 
  const unsigned long  first     ,
  const unsigned long  last      )
{
  //
  const unsigned long num_entries = data.numEntries() ;
  const unsigned long the_last    = std::min ( num_entries , last ) ;
  if ( the_last <= first ) { return 0 ; }                                     //  RETURN
  //
  const bool  weighted   = data.isWeighted () ;
  const bool  with_cuts  = !cuts    .empty () ;
  //
  const char* cutrange = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  // the simplest case :
  if ( !with_cuts && !cutrange && !weighted ) { return the_last - first ; } // RETURN
  //
  const std::unique_ptr<Ostap::FormulaVar> cut { make_formula ( cuts , data , true ) } ;
  //
  long double sumw  = 0    ;
  long double sumw2 = 0    ;
  bool        empty = true ; //  empty  dataset (after cuts&selection) ? 
  for ( unsigned long entry = first ; entry < the_last ; ++entry )
  {
    const RooArgSet* vars = data.get( entry ) ;
    if ( nullptr == vars )                            { break    ; } // BREAK 
    //
    if ( cutrange && !vars->allInRange ( cutrange ) ) { continue ; } // CONTINUE    
    // apply cuts:
    const long double wc = cut       ? cut -> getVal() : 1.0L ;
    if ( !wc ) { continue ; }                                   // CONTINUE  
    // apply weight:
    const long double wd = weighted  ? data.weight()   : 1.0L ;
    if ( !wd ) { continue ; }                                   // CONTINUE    
    // cuts & weight:
    const long double w  = wd *  wc ;
    if ( !w  ) { continue ; }                                   // CONTINUE        
    //
    sumw  += w     ;
    sumw2 += w * w ;
    empty  = false ;
    //
  }
  //
  return empty ? 0.0 : sumw * sumw / sumw2 ;
}
// ============================================================================
/* calculate the moment of order "order" relative to the center "center"
 *  @param  data      (INPUT) input data
 *  @param  expr      (INPUT) expression  (must  be valid TFormula!)
 *  @param  order     (INPUT) the order 
 *  @param  center    (INPUT) the center 
 *  @param  cuts      (INPUT) cuts 
 *  @param  first     (INPUT) the first  event to process 
 *  @param  last      (INPUT) the last event to  process
 *  @param  cut_range (INPUT) cut range 
 *  @return the moment 
 */
// ============================================================================
double Ostap::StatVar::get_moment 
( const RooAbsData&    data      ,  
  const unsigned short order     , 
  const std::string&   expr      , 
  const double         center    ,
  const std::string&   cuts      , 
  const std::string&   cut_range , 
  const unsigned long  first     ,
  const unsigned long  last      )
{
  //
  if (  0 == order       ) { return  1 ; }    // RETURN
  //
  const unsigned long num_entries = data.numEntries() ;
  const unsigned long the_last    = std::min ( num_entries , last ) ;
  if ( the_last <= first ) { return  0 ; }    // RETURN
  //
  const char* cutrange = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  const std::unique_ptr<Ostap::FormulaVar> expression { make_formula ( expr , data        ) } ;
  const std::unique_ptr<Ostap::FormulaVar> cut        { make_formula ( cuts , data , true ) } ;
  //  
  return _moment_ ( data      , *expression ,
                    cut.get() ,  
                    order     , center   , 
                    first     , the_last , cutrange ) ;
}
// ============================================================================
/*  calculate the moment of order "order"
 *  @param  data  (INPUT) input data
 *  @param  order (INPUT) the order 
 *  @param  expr  (INPUT) expression  (must  be valid TFormula!)
 *  @param  cuts  (INPUT) cuts 
 *  @param  first (INPUT) the first  event to process 
 *  @param  last  (INPUT) the last event to  process
 *  @return the moment 
 */ 
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::StatVar::moment
( const  RooAbsData&   data      ,  
  const unsigned short order     ,
  const std::string&   expr      , 
  const std::string&   cuts      , 
  const std::string&   cut_range , 
  const unsigned long  first     ,
  const unsigned long  last      )
{
  if ( 0 == order )        { return  1 ; }    // RETURN
  //
  const unsigned long num_entries = data.numEntries() ;
  const unsigned long the_last    = std::min ( num_entries , last ) ;
  if ( the_last <= first ) { return  0 ; }    // RETURN
  //
  const bool  weighted = data.isWeighted () ;
  const bool  with_cuts  = !cuts    .empty () ;
  //
  const std::unique_ptr<Ostap::FormulaVar> expression { make_formula ( expr , data        ) } ;
  const std::unique_ptr<Ostap::FormulaVar> cut        { make_formula ( cuts , data , true ) } ;
  //  
  long double mom   = 0 ;
  long double sumw  = 0 ; // sum of weights 
  // for uncertainties 
  long double sumw2 = 0 ; // sum of weights^2
  long double c2    = 0 ;
  //
  const char* cutrange = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  bool        empty = true ;
  for ( unsigned long entry = first ; entry < the_last ; ++entry )
  {
    const RooArgSet* vars = data.get( entry ) ;
    if ( nullptr == vars )                            { break    ; } // BREAK 
    //
    if ( cutrange && !vars->allInRange ( cutrange ) ) { continue ; } // CONTINUE    
    // apply cuts:
    const long double wc = with_cuts ? cut -> getVal() : 1.0L ;
    if ( !wc ) { continue ; }                                   // CONTINUE  
    // apply weight:
    const long double wd = weighted  ? data.weight()   : 1.0L ;
    if ( !wd ) { continue ; }                                   // CONTINUE    
    // cuts & weight:
    const long double w  = wd *  wc ;
    if ( !w  ) { continue ; }                                   // CONTINUE        
    //
    const double x = expression->getVal() ;
    //
    mom   += w  * std::pow ( x , order ) ;
    sumw  += w     ;
    // for uncertainty:
    sumw2 += w * w ;
    c2    += w * std::pow  ( x , 2 * order ) ;
    // 
    empty  =  false ;
    //
  }
  //
  if ( empty ) { return 0 ; }  //  RETURN 
  //
  const long double v = mom / sumw ;
  //
  c2 /= sumw  ; // the moment of "2*order"
  c2 -= v * v ; // m(2*order) - m(order)**2
  //
  const long double n = sumw * sumw / sumw2 ;
  c2 /= n     ; // 
  //
  return Ostap::Math::ValueWithError ( v , c2 ) ;            // RETURN
}
// ============================================================================
/*  calculate the central moment of order "order"
 *  @param  data  (INPUT) input data
 *  @param  order (INPUT) the order 
 *  @param  expr  (INPUT) expression  (must  be valid TFormula!)
 *  @param  cuts  (INPUT) cuts 
 *  @param  first (INPUT) the first  event to process 
 *  @param  last  (INPUT) the last event to  process
 *  @param  cut_range (INPUT) cut range 
 *  @return the moment 
 */ 
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::StatVar::central_moment 
( const RooAbsData&    data      ,  
  const unsigned short order     ,
  const std::string&   expr      , 
  const std::string&   cuts      , 
  const std::string&   cut_range , 
  const unsigned long  first     ,
  const unsigned long  last      )
{
  if      ( 0 == order )        { return  1 ; }    // RETURN
  else if ( 1 == order )        { return  0 ; }    // RETURN
  //
  const unsigned long num_entries = data.numEntries() ;
  const unsigned long the_last    = std::min ( num_entries , last ) ;
  if ( the_last <= first ) { return  0 ; }    // RETURN
  //
  const bool  weighted = data.isWeighted () ;
  const bool  with_cuts  = !cuts    .empty () ;
  //
  const std::unique_ptr<Ostap::FormulaVar> expression { make_formula ( expr , data        ) } ;
  const std::unique_ptr<Ostap::FormulaVar> cut        { make_formula ( cuts , data , true ) } ;
  //
  const char* cutrange = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  const long double mu =
    _moment_ ( data ,  *expression , cut.get() , 1  , 0 , first ,  the_last , cutrange );
  //
  
  long double mom   = 0    ;
  long double sumw  = 0    ; // sum of weights 
  // for uncertainty:
  long double sumw2 = 0    ; // sum of weights^2
  long double m2o   = 0    ; // moment of 2*order 
  long double mm1   = 0    ; // moment of   order-1
  long double mp1   = 0    ; // moment of   order+1
  long double m2    = 0    ; // moment of 2
  bool        empty = true ;  
  for ( unsigned long entry = first ; entry < the_last ; ++entry )
  {
    const RooArgSet* vars = data.get( entry ) ;
    if ( nullptr == vars )                            { break    ; } // BREAK 
    //
    if ( cutrange && !vars->allInRange ( cutrange ) ) { continue ; } // CONTINUE    
    // apply cuts:
    const long double wc = with_cuts ? cut -> getVal() : 1.0L ;
    if ( !wc ) { continue ; }                                   // CONTINUE  
    // apply weight:
    const long double wd = weighted  ? data.weight()   : 1.0L ;
    if ( !wd ) { continue ; }                                   // CONTINUE    
    // cuts & weight:
    const long double w  = wd *  wc ;
    if ( !w  ) { continue ; }                                   // CONTINUE        
    //
    const double dx = expression->getVal() - mu ;
    //
    mom   += w  * std::pow ( dx , order ) ;
    sumw  += w     ;
    // for uncertainty:
    sumw2 += w * w ;
    m2o   += w * std::pow ( dx , 2 * order     ) ; 
    mm1   += w * std::pow ( dx ,     order - 1 ) ; 
    mp1   += w * std::pow ( dx ,     order + 1 ) ; 
    m2    += w * std::pow ( dx , 2             ) ;
    //
    empty  = false ;
  }
  //
  if ( empty ) { return 0 ; }  //  RETURN 
  //
    // number of effective entries:
    const long double n = sumw * sumw / sumw2 ;
    long double v = mom / sumw ;
    /// correct O(1/n) bias  for   3rd and 4th order moments :
    if      ( 3 == order ) { v *=  n * n / ( ( n - 1  ) * ( n - 2 ) ) ; }
    else if ( 4 == order ) 
    {
      const double n0 =  ( n - 1 ) * ( n - 2 ) * ( n - 3 ) ;
      const double n1 =  n * ( n * n - 2 * n +  3 ) / n0   ;
      const double n2 =  3 * n * ( 2 * n - 3 )      / n0   ;
      v = n1 * v - n2 * m2 * m2 / ( sumw * sumw ) ;
    }
    //
    m2o /= sumw ;
    mm1 /= sumw ;
    mp1 /= sumw ;
    m2  /= sumw ;
    //
    long double c2 = m2o  ;
    c2 -= 2 * order * mm1 * mp1 ;
    c2 -= v * v ;
    c2 += order * order * m2 * mm1 * mm1 ;
    c2 /= n  ;
    //
    return Ostap::Math::ValueWithError ( v , c2 ) ;            // RETURN
}
// ============================================================================
/*  calculate the skewness of the  distribution
 *  @param  tree      (INPUT) input tree 
 *  @param  expr      (INPUT) expression  (must  be valid TFormula!)
 *  @param  cuts      (INPUT) cuts 
 *  @param  first     (INPUT) the first  event to process 
 *  @param  last      (INPUT) the last event to  process
 *  @param  cut_range (INPUT) cut range 
 *  @return the skewness 
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::StatVar::skewness 
( const RooAbsData&    data      ,  
  const std::string&   expr      , 
  const std::string&   cuts      , 
  const std::string&   cut_range , 
  const unsigned long  first     ,
  const unsigned long  last      )
{ 
  //
  const unsigned long num_entries = data.numEntries() ;
  const unsigned long the_last    = std::min ( num_entries , last ) ;
  if ( the_last <= first ) { return  0 ; }    // RETURN
  //
  const bool  weighted = data.isWeighted () ;
  const bool  with_cuts  = !cuts    .empty () ;
  //
  const std::unique_ptr<Ostap::FormulaVar> expression { make_formula ( expr , data        ) } ;
  const std::unique_ptr<Ostap::FormulaVar> cut        { make_formula ( cuts , data , true ) } ;
  //
  const char* cutrange = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  const long double mu =
    _moment_ ( data ,  *expression , cut.get() , 1  , 0 , first ,  the_last , cutrange );
  //
  long double mom   = 0 ;
  long double sumw  = 0 ; // sum of weights 
  // for uncertainty:
  long double sumw2 = 0 ; // sum of weights^2
  long double m2    = 0 ; // moment of 2
  //
  bool        empty = true ;
  //
  for ( unsigned long entry = first ; entry < the_last ; ++entry )
  {
    const RooArgSet* vars = data.get( entry ) ;
    if ( nullptr == vars )                            { break    ; } // BREAK 
    //
    if ( cutrange && !vars->allInRange ( cutrange ) ) { continue ; } // CONTINUE    
    // apply cuts:
    const long double wc = with_cuts ? cut -> getVal() : 1.0L ;
    if ( !wc ) { continue ; }                                   // CONTINUE  
    // apply weight:
    const long double wd = weighted  ? data.weight()   : 1.0L ;
    if ( !wd ) { continue ; }                                   // CONTINUE    
    // cuts & weight:
    const long double w  = wd *  wc ;
    if ( !w  ) { continue ; }                                   // CONTINUE        
    //
    const double dx = expression->getVal() - mu ;
    //
    mom   += w * std::pow ( dx , 3 ) ;
    sumw  += w     ;
    // for uncertainty:
    sumw2 += w * w ;
    m2    += w * std::pow ( dx , 2 ) ;
    //
    empty  = false ;
  }
  //
  if ( empty ) {  return 0 ; }
  // number of effective entries:
  const long double n = sumw * sumw / sumw2 ;
  //
  long double v = mom / sumw ;
  /// correct O(1/n) bias  for 3rd moment 
  v *=  n * n / ( ( n - 1  ) * ( n - 2 ) ) ; 
  //
  m2 /= sumw  ;
  v  /= std::pow ( m2 , 1.5 ) ;
  //
  long double c2 = 6  ;
  c2 *= ( n - 2 )  ;
  c2 /= ( n + 1 ) * ( n + 3 ) ;    
  //
  return Ostap::Math::ValueWithError ( v , c2 ) ;            // RETURN
}
// ============================================================================
/*  calculate the (excess) kurtosis of the  distribution
 *  @param  data  (INPUT) input data
 *  @param  expr  (INPUT) expression  (must  be valid TFormula!)
 *  @param  cuts  (INPUT) cuts 
 *  @param  first (INPUT) the first  event to process 
 *  @param  last  (INPUT) the last event to  process
 *  @return the (excess) kurtosis
 */
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::StatVar::kurtosis
( const RooAbsData&    data      ,  
  const std::string&   expr      , 
  const std::string&   cuts      , 
  const std::string&   cut_range , 
  const unsigned long  first     ,
  const unsigned long  last      ) 
{
  //
  const unsigned long num_entries = data.numEntries() ;
  const unsigned long the_last    = std::min ( num_entries , last ) ;
  if ( the_last <= first ) { return  0 ; }    // RETURN
  //
  const bool  weighted  = data.isWeighted () ;
  const bool  with_cuts = !cuts    .empty () ;
  //
  const std::unique_ptr<Ostap::FormulaVar> expression { make_formula ( expr , data        ) } ;
  const std::unique_ptr<Ostap::FormulaVar> cut        { make_formula ( cuts , data , true ) } ;
  //
  const char* cutrange = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  const long double mu =
    _moment_ ( data ,  *expression , cut.get() , 1  , 0 , first ,  the_last , cutrange );
  //
  long double mom   = 0 ;
  long double sumw  = 0 ; // sum of weights 
  // for uncertainty:
  long double sumw2 = 0 ; // sum of weights^2
  long double m2    = 0 ; // moment of 2
  //
  bool        empty = true ;
  //
  for ( unsigned long entry = first ; entry < the_last ; ++entry )
  {
    const RooArgSet* vars = data.get( entry ) ;
    if ( nullptr == vars )                            { break    ; } // BREAK 
    //
    if ( cutrange && !vars->allInRange ( cutrange ) ) { continue ; } // CONTINUE    
    // apply cuts:
    const long double wc = with_cuts ? cut -> getVal() : 1.0L ;
    if ( !wc ) { continue ; }                                   // CONTINUE  
    // apply weight:
    const long double wd = weighted  ? data.weight()   : 1.0L ;
    if ( !wd ) { continue ; }                                   // CONTINUE    
    // cuts & weight:
    const long double w  = wd *  wc ;
    if ( !w  ) { continue ; }                                   // CONTINUE        
    //
    const double dx = expression->getVal() - mu ;
    //
    mom   += w * std::pow ( dx , 4 ) ;
    sumw  += w     ;
    // for uncertainty:
    sumw2 += w * w ;
    m2    += w * std::pow ( dx , 2 ) ;
    //
    empty  = false ;
  }
  //
  if ( empty ) { return 0 ; } // RETURN
  //
  // number of effective entries:
  const long double n = sumw * sumw / sumw2 ;
  //
  long double v = mom / sumw ;
  m2 /= sumw  ; // second order moment
  /// correct for O(1/n) bias:
  const double n0 =  ( n - 1 ) * ( n - 2 ) * ( n - 3 ) ;
  const double n1 =  n * ( n * n - 2 * n +  3 ) / n0   ;
  const double n2 =  3 * n * ( 2 * n - 3 )      / n0   ;
  v  = n1 * v - n2 * m2 * m2 ;
  /// normalize  it: 
  v /= std::pow ( m2 , 2 ) ;
  ///
  long double c2 = 24 * n ;
  c2 *= ( n - 2 ) * ( n - 3 ) ;
  c2 /= ( n + 1 ) * ( n + 1 ) ;
  c2 /= ( n + 3 ) * ( n + 5 ) ;
  //
  return Ostap::Math::ValueWithError ( v , c2 ) ;            // RETURN
}
// ============================================================================
/*   get quantile of the distribution  
 *   @param data   (INPUT) the input data
 *   @param q      (INPUT) quantile value   0 < q < 1  
 *   @param expr   (INPUT) the expression 
 *   @param cuts   (INPUT) selection cuts 
 *   @param  first (INPUT) the first  event to process 
 *   @param  last  (INPUT) the last event to  process
 *   @param  cut_range (INPUT) cut range 
 *   @return the quantile value 
 */
// ============================================================================
Ostap::StatVar::Quantile
Ostap::StatVar::quantile
( const RooAbsData&   data      ,
  const double        q         , //  0<q<1 
  const std::string&  expr      , 
  const std::string&  cuts      , 
  const std::string&  cut_range , 
  const unsigned long first     ,
  const unsigned long last      )
{
  Ostap::Assert ( 0 < q && q < 1             , 
                  "Invalid quantile"         ,
                  "Ostap::StatVar::quantile" ) ;
  //
  const unsigned long num_entries = data.numEntries() ;
  const unsigned long the_last    = std::min ( num_entries , last ) ;
  if ( the_last <= first ) { return  0 ; }    // RETURN
  //
  const std::unique_ptr<Ostap::FormulaVar> expression { make_formula ( expr , data        ) } ;
  const std::unique_ptr<Ostap::FormulaVar> cut        { make_formula ( cuts , data , true ) } ;
  //  
  std::set<double> qset {} ;
  qset.insert ( q )  ;
  //
  const char* cutrange = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  auto result = _quantiles_ 
    (  data,  qset , 
       *expression ,  cut.get() , first , the_last , cutrange ) ;
  //
  Ostap::Assert ( 1 == result.quantiles.size()         ,
                  "Invalid quantile size"    ,
                  "Ostap::StatVar::quantile" ) ;
  //
  return Quantile ( result.quantiles[0] , result.nevents ) ;
}
// ============================================================================
/*   get (approximate) quantile of the distribution using P^2 algorithm
 *   @param data   (INPUT) the input data
 *   @param q      (INPUT) quantile value   0 < q < 1  
 *   @param expr   (INPUT) the expression 
 *   @param cuts   (INPUT) selection cuts 
 *   @param  first (INPUT) the first  event to process 
 *   @param  last  (INPUT) the last event to  process
 *   @param  cut_range (INPUT) cut range 
 *   @return the quantile value 
 */
// ============================================================================
Ostap::StatVar::Quantile
Ostap::StatVar::p2quantile
( const RooAbsData&   data      ,
  const double        q         , //  0<q<1 
  const std::string&  expr      , 
  const std::string&  cuts      , 
  const std::string&  cut_range , 
  const unsigned long first     ,
  const unsigned long last      )
{
  Ostap::Assert ( 0 < q && q < 1             , 
                  "Invalid quantile"         ,
                  "Ostap::StatVar::quantile" ) ;
  //
  const unsigned long num_entries = data.numEntries() ;
  const unsigned long the_last    = std::min ( num_entries , last ) ;
  if ( the_last <= first ) { return  0 ; }    // RETURN
  //
  const std::unique_ptr<Ostap::FormulaVar> expression { make_formula ( expr , data        ) } ;
  const std::unique_ptr<Ostap::FormulaVar> cut        { make_formula ( cuts , data , true ) } ;
  //  
  const char* cutrange = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  std::set<double> qset {} ;
  qset.insert ( q )  ;
  //
  auto result = _p2quantiles_ 
    (  data,  qset, 
       *expression ,  cut.get() , first , the_last , cutrange ) ;
  //
  Ostap::Assert ( 1 == result.quantiles.size()         ,
                  "Invalid quantile size"    ,
                  "Ostap::StatVar::quantile" ) ;
  //
  return Quantile ( result.quantiles[0] , result.nevents ) ;
}
// ============================================================================
/*   get the interval of the distribution  
 *   @param data  (INPUT) the input data
 *   @param q1    (INPUT) quantile value   0 < q1 < 1  
 *   @param q2    (INPUT) quantile value   0 < q2 < 1  
 *   @param expr  (INPUT) the expression 
 *   @param cuts  (INPUT) selection cuts 
 *   @param  first (INPUT) the first  event to process 
 *   @param  last  (INPUT) the last event to  process
 *   @return the quantile value 
 *   @code
 *   const RooAbsData& data = ... ;
 *   /// get 90% interval:
 *   Interval ab = interval ( data , 0.05 , 0.95 , 'mass' , 'pt>3' ) ;
 *   @code 
 */
// ============================================================================
Ostap::StatVar::QInterval 
Ostap::StatVar::interval 
( const RooAbsData&   data      ,
  const double        q1        , //  0<q1<1 
  const double        q2        , //  0<q2<1 
  const std::string&  expr      , 
  const std::string&  cuts      , 
  const std::string&  cut_range , 
  const unsigned long first     ,
  const unsigned long last      )
{
  Ostap::Assert ( 0 < q1 && q1 < 1           , 
                  "Invalid quantile1"        ,
                  "Ostap::StatVar::quantile" ) ;
  Ostap::Assert ( 0 < q2 && q2 < 1           , 
                  "Invalid quantile2"        ,
                  "Ostap::StatVar::quantile" ) ;
  //
  const unsigned long num_entries = data.numEntries() ;
  const unsigned long the_last    = std::min ( num_entries , last ) ;
  if ( the_last <= first ) { return QInterval() ; }    // RETURN
  //
  const std::unique_ptr<Ostap::FormulaVar> expression { make_formula ( expr , data        ) } ;
  const std::unique_ptr<Ostap::FormulaVar> cut        { make_formula ( cuts , data , true ) } ;
  //  
  const char* cutrange = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  auto result = _quantiles_ 
    (  data,  std::set<double>{{ q1 , q2 }} , 
       *expression ,  cut.get() , first , the_last , cutrange ) ;
  Ostap::Assert ( 2 ==  result.quantiles.size()        ,
                  "Invalid quantile size"    ,
                  "Ostap::StatVar::quantile" ) ;
  //
  return QInterval ( Interval ( result.quantiles[0] , result.quantiles [1] ) , result.nevents ) ;
}
// ============================================================================
/*   get the approximate  interval of the distribution  usnig  P^2 algorithm
 *   @param data  (INPUT) the input data
 *   @param q1    (INPUT) quantile value   0 < q1 < 1  
 *   @param q2    (INPUT) quantile value   0 < q2 < 1  
 *   @param expr  (INPUT) the expression 
 *   @param cuts  (INPUT) selection cuts 
 *   @param  first (INPUT) the first  event to process 
 *   @param  last  (INPUT) the last event to  process
 *   @return the quantile value 
 *   @code
 *   const RooAbsData& data = ... ;
 *   /// get 90% interval:
 *   Interval ab = p2interval ( data , 0.05 , 0.95 , 'mass' , 'pt>3' ) ;
 *   @code 
 */
// ============================================================================
Ostap::StatVar::QInterval 
Ostap::StatVar::p2interval 
( const RooAbsData&   data      ,
  const double        q1        , //  0<q1<1 
  const double        q2        , //  0<q2<1 
  const std::string&  expr      , 
  const std::string&  cuts      , 
  const std::string&  cut_range ,   
  const unsigned long first     ,
  const unsigned long last      )
{
  Ostap::Assert ( 0 < q1 && q1 < 1           , 
                  "Invalid quantile1"        ,
                  "Ostap::StatVar::quantile" ) ;
  Ostap::Assert ( 0 < q2 && q2 < 1           , 
                  "Invalid quantile2"        ,
                  "Ostap::StatVar::quantile" ) ;
  //
  const unsigned long num_entries = data.numEntries() ;
  const unsigned long the_last    = std::min ( num_entries , last ) ;
  if ( the_last <= first ) { return  QInterval()  ; }    // RETURN
  //
  const std::unique_ptr<Ostap::FormulaVar> expression { make_formula ( expr , data        ) } ;
  const std::unique_ptr<Ostap::FormulaVar> cut        { make_formula ( cuts , data , true ) } ;
  //  
  const char* cutrange = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  auto result = _p2quantiles_ 
    (  data,  std::set<double>{{ q1 , q2 }} , 
       *expression ,  cut.get() , first , the_last , cutrange ) ;
  Ostap::Assert ( 2 ==  result.quantiles.size()        ,
                  "Invalid quantile size"    ,
                  "Ostap::StatVar::quantile" ) ;
  //
  return QInterval ( Interval ( result.quantiles[0] , result.quantiles [1] ) , result.nevents ) ;
}
// ============================================================================
/*   get quantiles of the distribution  
 *   @param data  (INPUT) the input data
 *   @param q     (INPUT) quantile value   0 < q < 1  
 *   @param expr  (INPUT) the expression 
 *   @param cuts  (INPUT) selection cuts 
 *   @param  first (INPUT) the first  event to process 
 *   @param  last  (INPUT) the last event to  process
 *   @return the quantile value 
 */
// ============================================================================
Ostap::StatVar::Quantiles
Ostap::StatVar::quantiles
( const RooAbsData&          data      ,
  const std::vector<double>& quantiles , 
  const std::string&         expr      , 
  const std::string&         cuts      , 
  const std::string&         cut_range ,   
  const unsigned long        first     ,
  const unsigned long        last      )
{
  std::set<double> qs ;
  for ( double v : quantiles ) { qs.insert( v ) ; }
  Ostap::Assert (  !qs.empty()                 , 
                   "Invalid quantiles"         ,
                   "Ostap::StatVar::quantiles" ) ;
  Ostap::Assert ( 0 < *qs. begin ()            ,  
                  "Invalid quantile"           ,
                  "Ostap::StatVar::quantiles"  ) ;
  Ostap::Assert ( 1 > *qs.rbegin ()            ,  
                  "Invalid quantile"           ,
                  "Ostap::StatVar::quantiles"  ) ;
  //
  const unsigned long num_entries = data.numEntries() ;
  const unsigned long the_last    = std::min ( num_entries , last ) ;
  if ( the_last <= first ) { return  std::vector<double>() ; }    // RETURN
  //
  const char* cutrange = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  const std::unique_ptr<Ostap::FormulaVar> expression { make_formula ( expr , data        ) } ;
  const std::unique_ptr<Ostap::FormulaVar> cut        { make_formula ( cuts , data , true ) } ;
  //  
  return _quantiles_ ( data  , qs  , 
                       *expression , cut.get() , 
                       first , the_last , cutrange ) ;
}
// ============================================================================
/*   get (approximate) quantiles of the distribution using P^2 algorithm  
 *   @param data  (INPUT) the input data
 *   @param q     (INPUT) quantile value   0 < q < 1  
 *   @param expr  (INPUT) the expression 
 *   @param cuts  (INPUT) selection cuts 
 *   @param  first (INPUT) the first  event to process 
 *   @param  last  (INPUT) the last event to  process
 *   @return the quantile value 
 */
// ============================================================================
Ostap::StatVar::Quantiles
Ostap::StatVar::p2quantiles
( const RooAbsData&          data      ,
  const std::vector<double>& quantiles , 
  const std::string&         expr      , 
  const std::string&         cuts      , 
  const std::string&         cut_range ,   
  const unsigned long        first     ,
  const unsigned long        last      )
{
  std::set<double> qs ;
  for ( double v : quantiles ) { qs.insert( v ) ; }
  Ostap::Assert (  !qs.empty()                 , 
                   "Invalid quantiles"         ,
                   "Ostap::StatVar::quantiles" ) ;
  Ostap::Assert ( 0 < *qs. begin ()            ,  
                  "Invalid quantile"           ,
                  "Ostap::StatVar::quantiles"  ) ;
  Ostap::Assert ( 1 > *qs.rbegin ()            ,  
                  "Invalid quantile"           ,
                  "Ostap::StatVar::quantiles"  ) ;
  //
  const unsigned long num_entries = data.numEntries() ;
  const unsigned long the_last    = std::min ( num_entries , last ) ;
  if ( the_last <= first )
  { return  Ostap::StatVar::Quantiles  ( std::vector<double>() , 0 ) ; }    // RETURN
  //
  const char* cutrange = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  const std::unique_ptr<Ostap::FormulaVar> expression { make_formula ( expr , data        ) } ;
  const std::unique_ptr<Ostap::FormulaVar> cut        { make_formula ( cuts , data , true ) } ;
  //  
  return _p2quantiles_ ( data  , qs  , 
                         *expression , cut.get() , 
                         first , the_last , cutrange ) ;
}
// ============================================================================
// Actions with frames 
// ============================================================================
/*  get the number of equivalent entries 
 *  \f$ n_{eff} \equiv = \frac{ (\sum w)^2}{ \sum w^2} \f$
 *  @param tree      (INPUT) the frame 
 *  @param cuts      (INPUT) selection criteria   (used as "weight")
 *  @return number of equivalent entries 
 */
// ============================================================================
double Ostap::StatVar::nEff 
( Ostap::FrameNode     frame ,
  const std::string&   cuts  )
{
  //
  const bool no_cuts = trivial ( cuts ) ;
  if ( no_cuts ) { return  *frame.Count(); }                     // RETURN
  //
  /// define temporary columns 
  const std::string weight  { Ostap::tmp_name ( "w_"  , cuts ) } ;
  const std::string weight2 { Ostap::tmp_name ( "w2_" , cuts ) } ;
  const std::string bcut    { Ostap::tmp_name ( "b_"  , cuts ) } ;
  /// decorate the frame
  auto t = frame
    .Define ( bcut     , "(bool)   ( " + cuts + " ) ;" ) 
    .Filter ( bcut     )
    .Define ( weight   , "(double) ( " + cuts + " ) ;" )
    .Define ( weight2  , []( double v ) { return v * v ; } , { weight } ) ;
  //
  const double zero = 0.0 ;
  auto  sumw_  = t.Reduce ( std::plus<double>() ,  weight   , zero ) ;
  auto  sumw2_ = t.Reduce ( std::plus<double>() ,  weight2  , zero ) ;
  //
  const double sumw  = *sumw_  ;
  const double sumw2 = *sumw2_ ;
  //
  return !sumw2 ? 0.0 : sumw * sumw / sumw2 ;
}
// ============================================================================
/*  build statistic for the <code>expression</code>
 *  @param data       (INPUT) the data 
 *  @param expression (INPUT) the expression
 *  @param cuts       (INPUT) the selection
 *
 *  @code
 *  data = ... 
 *  stat = data.statVar( 'S_sw' , 'pt>10') 
 *  @endcode 
 *
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2018-06-18
 */
// ============================================================================
Ostap::WStatEntity
Ostap::StatVar::statVar 
( Ostap::FrameNode    frame      , 
  const std::string&  expression , 
  const std::string&  cuts       ) 
{
  //
  const bool no_cuts = trivial ( cuts ) ; 
  //
  /// define the temporary columns 
  const std::string var    { Ostap::tmp_name ( "v_" , expression ) } ;
  const std::string weight { Ostap::tmp_name ( "w_" , cuts       ) } ;
  const std::string bcut   { Ostap::tmp_name ( "b_" , cuts       ) } ;
  /// define actions 
  auto t = frame
    .Define ( bcut   , no_cuts ? "true" : "(bool)   ( " + cuts + " ) ;" ) 
    .Filter ( bcut   ) 
    .Define ( var    ,  "1.0*(" + expression + ")"   )
    .Define ( weight , no_cuts ? "1.0"  : "1.0*(" + cuts + ")" ) ;
  //
  const unsigned int nSlots = std::max ( 1u , Ostap::Utils::mt_pool_size() ) ;
  //
  WStatVector _stat ( nSlots ? nSlots : 1 ) ;
  //
  auto fun = [&_stat,nSlots] ( unsigned int slot , double v , double w ) 
    { _stat [ slot % nSlots ].add ( v , w ) ; } ;
  t.ForeachSlot ( fun ,  { var , weight } ) ; 
  //
  Ostap::WStatEntity stat ; for ( const auto& s : _stat ) { stat += s ; }
  return stat ;
}
// ============================================================================
/*  calculate the covariance of two expressions 
 *  @param frame (INPUT)  data frame 
 *  @param exp1  (INPUT)  the first  expresiion
 *  @param exp2  (INPUT)  the second expresiion
 *  @return covariance 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2024-07-22
 */
// ============================================================================
Ostap::Math::Covariance
Ostap::StatVar::statCov
( Ostap::FrameNode     frame , 
  const std::string&   exp1  , 
  const std::string&   exp2  )  
{
  /// define the temporary columns 
  const std::string var1   = Ostap::tmp_name ( "v_" , exp1 ) ;
  const std::string var2   = Ostap::tmp_name ( "v_" , exp2 ) ;
  auto t = frame
    .Define ( var1   ,                   "1.0*(" + exp1 + ")" ) 
    .Define ( var2   ,                   "1.0*(" + exp2 + ")" ) ;
  ///
  const unsigned int nSlots = std::max ( 1u , Ostap::Utils::mt_pool_size() ) ;
  //
  std::vector<Ostap::Math::Covariance>           _covs ( nSlots ? nSlots : 1 ) ;
  //
  auto fun = [&_covs,nSlots] ( unsigned int slot , double v1 , double v2 )
  { _covs [ slot % nSlots ].add ( v1 , v2 ) ; } ;
  t.ForeachSlot ( fun , { var1 , var2 } ); 
  //
  Ostap::Math::Covariance result ;
  for ( const auto& s : _covs ) { result += s ; }
  //
  return result ;
}
// ============================================================================
/*  calculate the covariance of two expressions 
 *  @param frame (INPUT)  data frame 
 *  @param exp1  (INPUT)  the first  expresiion
 *  @param exp2  (INPUT)  the second expresiion
 *  @param cuts  (INPUT)  the selection/weight  
 *  @return covariance 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2024-07-22
 */
// ============================================================================
Ostap::Math::WCovariance
Ostap::StatVar::statCov
( Ostap::FrameNode     frame , 
  const std::string&   exp1  , 
  const std::string&   exp2  , 
  const std::string&   cuts  ) 
{
  const bool no_cuts = trivial ( cuts ) ; 
  /// define the temporary columns 
  const std::string var1   = Ostap::tmp_name ( "v_" , exp1 ) ;
  const std::string var2   = Ostap::tmp_name ( "v_" , exp2 ) ;
  const std::string bcut   = Ostap::tmp_name ( "b_" , cuts ) ;
  const std::string weight = Ostap::tmp_name ( "w_" , cuts ) ;
  auto t = frame
    .Define ( bcut   , no_cuts ? "true" : "(bool) ( " + cuts + " ) ;" ) 
    .Filter ( bcut   )    
    .Define ( var1   ,                   "1.0*(" + exp1 + ")" ) 
    .Define ( var2   ,                   "1.0*(" + exp2 + ")" ) 
    .Define ( weight , no_cuts ? "1.0" : "1.0*(" + cuts + ")" ) ;
  ///
  const unsigned int nSlots = std::max ( 1u , Ostap::Utils::mt_pool_size() ) ;
  //
  std::vector<Ostap::Math::WCovariance>           _covs ( nSlots ? nSlots : 1 ) ;
  //
  auto fun = [&_covs,nSlots] ( unsigned int slot , double v1 , double v2 , double w )
  { if ( w )  { _covs [ slot % nSlots ].add ( v1 , v2 , w ) ;  } } ;
  t.ForeachSlot ( fun , { var1 , var2 , weight } ); 
  //
  Ostap::Math::WCovariance result ;
  for ( const auto& s : _covs ) { result += s ; }
  //
  return result ;
}
// ============================================================================
/*  calculate the moment of order "order" relative to the center "center"
 *  @param  frame  (INPUT) input frame 
 *  @param  expr   (INPUT) expression 
 *  @param  order  (INPUT) the order 
 *  @param  center (INPUT) the center 
 *  @param  cuts   (INPUT) cuts 
 *  @return the moment 
 */
// ============================================================================
double Ostap::StatVar::get_moment 
( Ostap::FrameNode     frame  ,  
  const unsigned short order  , 
  const std::string&   expr   , 
  const double         center ,    
  const std::string&   cuts   ) 
{
  //
  if ( 0 == order ) { return 1 ; } // RETURN 
  //
  const bool no_cuts =  trivial ( cuts ) ;
  ///
  /// define the temporary columns 
  const std::string var    = Ostap::tmp_name ( "v_" , expr ) ;
  const std::string bcut   = Ostap::tmp_name ( "b_" , cuts ) ;
  const std::string weight = Ostap::tmp_name ( "w_" , cuts ) ;
  const std::string mom    = Ostap::tmp_name ( "m_" , expr ) ;
  //
  auto t = frame
    .Define ( bcut   , no_cuts ? "true" : "(bool) ( " + cuts + " ) ;" ) 
    .Filter ( bcut   )    
    .Define ( var    ,                   "1.0*(" + expr + ")" ) 
    .Define ( weight , no_cuts ? "1.0" : "1.0*(" + cuts + ")" ) 
    .Define ( mom    , [order,center]( double v , double w )->double {
        const long double lv = v ;
        const long double lw = w ;
        return lw ? lw * std::pow ( lv -  center , order ) : 0.0L ;
      } , { var , weight } ) ;
  //
  const double zero = 0 ;
  auto _sumv = t.Reduce( std::plus<double>() , mom    ) ;
  auto _sumw = t.Reduce( std::plus<double>() , weight ) ;
  //
  const double sumv = *_sumv ; 
  const double sumw = *_sumw ;
  //
  return !sumw ? 0 : sumv / sumw ;
}
// ============================================================================
/*  calculate the moment of order "order"
 *  @param  tree  (INPUT) input tree 
 *  @param  order (INPUT) the order 
 *  @param  expr  (INPUT) expression 
 *  @param  cuts  (INPUT) cuts 
 *  @return the moment 
 */ 
// ============================================================================
Ostap::Math::ValueWithError Ostap::StatVar::moment
( Ostap::FrameNode     frame  ,  
  const unsigned short order  ,
  const std::string&   expr   , 
  const std::string&   cuts   ) 
{
  //
  if ( 0 == order ){ return 1 ; } // RETURN 
  //
  const bool no_cuts =  trivial ( cuts ) ; 
  ///
  /// define the temporary columns 
  const std::string var     = Ostap::tmp_name ( "v_"  , expr ) ;
  const std::string bcut    = Ostap::tmp_name ( "b_"  , cuts ) ;
  const std::string weight  = Ostap::tmp_name ( "w_"  , cuts ) ;
  const std::string weight2 = Ostap::tmp_name ( "w2_" , cuts ) ;
  const std::string vmom    = Ostap::tmp_name ( "m_"  , expr ) ;
  const std::string vmom2   = Ostap::tmp_name ( "m2_" , expr ) ;
  //
  // decorate the frame 
  auto t = frame
    .Define ( bcut    , no_cuts ? "true" : "(bool) ( " + cuts + " ) ;" ) 
    .Filter ( bcut    )    
    .Define ( var     ,                   "1.0*(" + expr + ")" ) 
    .Define ( weight  , no_cuts ? "1.0" : "1.0*(" + cuts + ")" ) 
    .Define ( weight2 , [] ( double w ) { return w * w ; } , { weight } ) 
    // 
    .Define ( vmom    , [order]( double v , double w )->double {
        const long double lv = v ;
        const long double lw = w ;
        return lw ? lw * std::pow ( lv , order ) : 0.0 ;
      } , { var , weight } )   
    //
    .Define ( vmom2   , [order]( double v , double w )->double {
        const long double lv = v ;
        const long double lw = w ;
        return lw ? lw * std::pow ( lv , 2 * order ) : 0.0 ;
      } , { var , weight } ) ;
  //
  auto _sum   = t.Reduce ( std::plus<double> () ,  vmom    ) ;
  auto _sum2  = t.Reduce ( std::plus<double> () ,  vmom2   ) ;
  auto _sumw  = t.Reduce ( std::plus<double> () ,  weight  ) ;
  auto _sumw2 = t.Reduce ( std::plus<double> () ,  weight2 ) ;
  //
  const long double sumw  = *_sumw  ;
  const long double sumw2 = *_sumw2 ;
  const long double sum   = *_sum   ;
  long double       sum2  = *_sum2  ;
  //
  if ( !sumw ) { return 0 ; }  // RETURN 
  //
  const long double v = sum / sumw ; // the moment of "order"
  sum2 /= sumw  ;                    // the moment of "2*order"
  sum2 -= v * v ;                    // m(2*order) - m(order)**2
  //
  const long double n = sumw * sumw / sumw2 ;
  sum2 /= n ; 
  //
  return Ostap::Math::ValueWithError ( v , sum2 ) ;
}
// ============================================================================
/*  calculate the central moment of order "order"
 *  @param  tree  (INPUT) input tree 
 *  @param  order (INPUT) the order 
 *  @param  expr  (INPUT) expression
 *  @param  cuts  (INPUT) cuts 
 *  @return the moment 
 */ 
// ============================================================================
Ostap::Math::ValueWithError 
Ostap::StatVar::central_moment 
( Ostap::FrameNode     frame ,  
  const unsigned short order ,
  const std::string&   expr  , 
  const std::string&   cuts  ) 
{
  //
  if      ( 0 == order ){ return 1 ; } // RETURN 
  else if ( 1 == order ){ return 0 ; } // RETURN 
  //
  const bool no_cuts = trivial (  cuts ) ; 
  /// get the mean-value  (1-loop) 
  const long double mu = get_moment ( frame , 1 , expr , 0.0 , cuts ) ;
  //
  /// define the temporary columns 
  const std::string var     = Ostap::tmp_name ( "v_"   , expr ) ;
  const std::string bcut    = Ostap::tmp_name ( "b_"   , cuts ) ;
  const std::string weight  = Ostap::tmp_name ( "w_"   , cuts ) ;
  const std::string weight2 = Ostap::tmp_name ( "w2_"  , cuts ) ;
  const std::string vmom    = Ostap::tmp_name ( "m_"   , expr ) ;
  const std::string vmom2   = Ostap::tmp_name ( "m2_"  , expr ) ;
  const std::string vmp1    = Ostap::tmp_name ( "mp1_" , expr ) ;
  const std::string vmm1    = Ostap::tmp_name ( "mm1_" , expr ) ;
  const std::string vm2     = Ostap::tmp_name ( "mo2_" , expr ) ;
  //
  // decorate the frame 
  auto t = frame
    .Define ( bcut    , no_cuts ? "true" : "(bool) ( " + cuts + " ) ;" ) 
    .Filter ( bcut    )    
    .Define ( var     ,                   "1.0*(" + expr + ")" ) 
    .Define ( weight  , no_cuts ? "1.0" : "1.0*(" + cuts + ")" ) 
    .Define ( weight2 , [] ( double w ) { return w * w ; } , { weight } ) 
    // 
    .Define ( vmom    , [order,mu]( double v , double w )->double 
              {
                const long double lv = v ;
                const long double lw = w ;
                return lw ? lw * std::pow ( lv - mu , order ) : 0.0 ;
              } , { var , weight } )   
    //
    .Define ( vmom2   , [order,mu]( double v , double w )->double 
              {
                const long double lv = v ;
                const long double lw = w ;
                return lw ? lw * std::pow ( lv - mu , 2 * order ) : 0.0 ;
              } , { var , weight } ) 
    //
    .Define ( vmp1    , [order,mu]( double v , double w )->double 
              {
                const long double lv = v ;
                const long double lw = w ;
                return lw ? lw * std::pow ( lv - mu , order + 1 ) : 0.0 ;
              } , { var , weight } )   
    //
    .Define ( vmm1   , [order,mu]( double v , double w )->double 
              {
                const long double lv = v ;
                const long double lw = w ;
                return lw ? lw * std::pow ( lv - mu , order - 1 ) : 0.0 ;
              } , { var , weight } ) 
    //
    .Define ( vm2   , [order,mu]( double v , double w )->double 
              {
                const long double lv = v ;
                const long double lw = w ;
                return lw ? lw * std::pow ( lv - mu , 2 ) : 0.0 ;
              } , { var , weight } ) ;
  //
  auto _mom   = t.Reduce ( std::plus<double>() , vmom     ) ;
  auto _mom2  = t.Reduce ( std::plus<double>() , vmom2    ) ;
  auto _mp1   = t.Reduce ( std::plus<double>() , vmp1     ) ;
  auto _mm1   = t.Reduce ( std::plus<double>() , vmm1     ) ;
  auto _m2    = t.Reduce ( std::plus<double>() , vm2      ) ;
  auto _sumw  = t.Reduce ( std::plus<double>() , weight   ) ;
  auto _sumw2 = t.Reduce ( std::plus<double>() , weight2  ) ;
  //
  const long double sumw  = *_sumw  ;
  //
  if ( !sumw ) { return  0 ; }  // RETURN 
  //
  const long double mom   = *_mom   ;
  const long double sumw2 = *_sumw2 ;
  //
  long double m2o         = *_mom2  ;
  long double mm1         = *_mm1   ; 
  long double mp1         = *_mp1   ;
  long double m2          = *_m2    ;
  //
  // number of effective entries:
  const long double n = sumw * sumw / sumw2 ;
  //  
  long double v = mom / sumw ;
  /// correct O(1/n) bias  for  3rd and 4th order moments :
  if      ( 3 == order ) { v *= n * n / ( ( n - 1 ) * ( n - 2 ) ) ; }
  else if ( 4 == order ) 
  {
    const double n0 =  ( n - 1 ) * ( n - 2 ) * ( n - 3 ) ;
    const double n1 =  n * ( n * n - 2 * n + 3 ) / n0 ;
    const double n2 =  3 * n * ( 2 * n - 3 )     / n0 ;
    v = n1 * v - n2 * m2 * m2 / ( sumw * sumw ) ;
  }
  //
  m2o /= sumw ;
  mm1 /= sumw ;
  mp1 /= sumw ;
  m2  /= sumw ;
  //
  long double c2 = m2o  ;
  c2 -= 2 * order * mm1 * mp1 ;
  c2 -= v * v ;
  c2 += order * order * m2 * mm1 * mm1 ;
  c2 /= n  ;
  //
  return Ostap::Math::ValueWithError ( v , c2 ) ;            // RETURN
}  
// ============================================================================
/*  calculate the skewness of the  distribution
 *  @param  frame (INPUT) frame
 *  @param  expr  (INPUT) expression
 *  @param  cuts  (INPUT) cuts 
 *  @return the skewness 
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::StatVar::skewness 
( Ostap::FrameNode     frame ,  
  const std::string&   expr  , 
  const std::string&   cuts  ) 
{
  //
  const bool no_cuts =  trivial (  cuts ) ;
  //
  /// get the mean-value  (1-loop) 
  const long double mu   = get_moment ( frame , 1 , expr , 0.0 , cuts ) ;
  //  
  /// define the temporary columns 
  const std::string var     = Ostap::tmp_name ( "v_"   , expr ) ;
  const std::string bcut    = Ostap::tmp_name ( "b_"   , cuts ) ;
  const std::string weight  = Ostap::tmp_name ( "w_"   , cuts ) ;
  const std::string weight2 = Ostap::tmp_name ( "w2_"  , cuts ) ;
  const std::string vmom3   = Ostap::tmp_name ( "m3_"  , expr ) ;
  const std::string vmom2   = Ostap::tmp_name ( "m2_"  , expr ) ;
  //
  // decorate the frame 
  auto t = frame
    .Define ( bcut    , no_cuts ? "true" : "(bool) ( " + cuts + " ) ;" ) 
    .Filter ( bcut    )    
    .Define ( var     ,                   "1.0*(" + expr + ")" ) 
    .Define ( weight  , no_cuts ? "1.0" : "1.0*(" + cuts + ")" ) 
    .Define ( weight2 , [] ( double w ) { return w * w ; } , { weight } ) 
    // 
    .Define ( vmom3   , [mu]( double v , double w )->double {
        const long double lv = v ;
        const long double lw = w ;
        return lw ? lw * std::pow ( lv - mu , 3 ) : 0.0 ;
      } , { var , weight } ) 
    //
    .Define ( vmom2   , [mu]( double v , double w )->double {
        const long double lv = v ;
        const long double lw = w ;
        return lw ? lw * std::pow ( lv - mu , 2 ) : 0.0 ;
      } , { var , weight } ) ;
  //
  auto _mom3  = t.Reduce ( std::plus<double>() , vmom3   ) ;
  auto _mom2  = t.Reduce ( std::plus<double>() , vmom2   ) ;
  auto _sumw  = t.Reduce ( std::plus<double>() , weight  ) ;
  auto _sumw2 = t.Reduce ( std::plus<double>() , weight2 ) ;
  //
  const long double sumw  = *_sumw ;
  
  //
  if ( !sumw ) { return  0 ; }  // RETURN 
  //
  const long double mom3  = *_mom3  ; 
  const long double sumw2 = *_sumw2 ;
  long double       mom2  = *_mom2  ;
  //
  // number of effective entries:
  const long double n = sumw * sumw / sumw2 ;
  //
  long double v = mom3 / sumw ;
  /// correct O(1/n) bias  for 3rd moment 
  v *=  n * n / ( ( n - 1  ) * ( n - 2 ) ) ; 
  //
  mom2 /= sumw  ;
  v  /= std::pow ( mom2 , 1.5 ) ;
  //
  long double c2 = 6  ;
  c2 *= ( n - 2 )  ;
  c2 /= ( n + 1 ) * ( n + 3 ) ;    
  //
  return Ostap::Math::ValueWithError ( v , c2 ) ;            // RETURN 
}
// ============================================================================
/*  calculate the (excess) kurtosis of the  distribution
 *  @param  frame (INPUT) input frame 
 *  @param  expr  (INPUT) expression  (must  be valid TFormula!)
 *  @param  cuts  (INPUT) cuts 
 *  @param  first (INPUT) the first  event to process 
 *  @param  last  (INPUT) the last event to  process
 *  @return the (excess) kurtosis
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::StatVar::kurtosis
( Ostap::FrameNode     frame ,  
  const std::string&   expr  , 
  const std::string&   cuts  )
{
  //
  const bool no_cuts = trivial (  cuts ) ; 
  //
  /// get the mean-value  (1-loop) 
  const long double mu   = get_moment ( frame , 1 , expr , 0.0 , cuts ) ;
  //  
  const std::string var     = Ostap::tmp_name ( "v_"   , expr ) ;
  const std::string bcut    = Ostap::tmp_name ( "b_"   , cuts ) ;
  const std::string weight  = Ostap::tmp_name ( "w_"   , cuts ) ;
  const std::string weight2 = Ostap::tmp_name ( "w2_"  , cuts ) ;
  const std::string vmom4   = Ostap::tmp_name ( "m4_"  , expr ) ;
  const std::string vmom2   = Ostap::tmp_name ( "m2_"  , expr ) ;
  //
  // decorate the frame 
  auto t = frame
    .Define ( bcut    , no_cuts ? "true" : "(bool) ( " + cuts + " ) ;" ) 
    .Filter ( bcut    )    
    .Define ( var     ,                   "1.0*(" + expr + ")" ) 
    .Define ( weight  , no_cuts ? "1.0" : "1.0*(" + cuts + ")" ) 
    .Define ( weight2 , [] ( double w ) { return w * w ; } , { weight } ) 
    // 
    .Define ( vmom4   , [mu]( double v , double w )->double {
        const long double lv = v ;
        const long double lw = w ;
        return lw ? lw * std::pow ( lv - mu , 4 ) : 0.0 ;
      } , { var , weight } ) 
    //
    .Define ( vmom2   , [mu]( double v , double w )->double {
        const long double lv = v ;
        const long double lw = w ;
        return lw ? lw * std::pow ( lv - mu , 2 ) : 0.0 ;
      } , { var , weight } ) ;
  //
  //
  auto _mom4  = t.Reduce ( std::plus<double>() , vmom4   ) ;
  auto _mom2  = t.Reduce ( std::plus<double>() , vmom2   ) ;
  auto _sumw  = t.Reduce ( std::plus<double>() , weight  ) ;
  auto _sumw2 = t.Reduce ( std::plus<double>() , weight2 ) ;
  //
  //
  const long double sumw  = *_sumw ; 
  //
  if ( !sumw ) { return  0 ; }  // RETURN 
  //
  const long double mom4  = *_mom4  ; 
  const long double sumw2 = *_sumw2 ;
  long double       mom2  = *_mom2  ; 
  //
  // number of effective entries:
  const long double n = sumw * sumw / sumw2 ;
  //
  long double v = mom4 / sumw ;
  mom2 /= sumw  ; // second order moment
  /// correct for O(1/n) bias:
  const double n0 =  ( n - 1 ) * ( n - 2 ) * ( n - 3 ) ;
  const double n1 =  n * ( n * n - 2 * n +  3 ) / n0   ;
  const double n2 =  3 * n * ( 2 * n - 3 )      / n0   ;
  v  = n1 * v - n2 * mom2 * mom2 ;
  /// normalize  it: 
  v /= std::pow ( mom2 , 2 ) ;
  ///
  long double c2 = 24 * n ;
  c2 *= ( n - 2 ) * ( n - 3 ) ;
  c2 /= ( n + 1 ) * ( n + 1 ) ;
  c2 /= ( n + 3 ) * ( n + 5 ) ;
  //
  return Ostap::Math::ValueWithError ( v , c2 ) ;            // RETURN 
}
// ============================================================================
namespace 
{
  // ==========================================================================
  /// get quantiles
  Ostap::StatVar::Quantiles _quantiles_ 
  ( Ostap::FrameNode        frame , 
    const std::set<double>& qs    , 
    const std::string&      expr  , 
    const std::string&      cuts  ) 
  {
    const bool no_cuts = trivial (  cuts ) ; 
    //
    const std::string var   = Ostap::tmp_name ( "v_"   , expr ) ;
    const std::string bcut  = Ostap::tmp_name ( "b_"   , cuts ) ;
    //
    // decorate the frame 
    auto t = frame
      .Define       ( bcut    , no_cuts ? "true" : "(bool) ( " + cuts + " ) ;" ) 
      .Filter       ( bcut    )    
      .Define       ( var     , "1.0*(" + expr + ")" ) 
      .Take<double> ( var     ) ;
    // 
    typedef  std::vector<double> DATA ;
    DATA values { *t } ;
    //
    std::vector<double> result ; result.reserve ( qs.size() ) ;
    //
    DATA::iterator start = values.begin() ;
    for ( const double q : qs )
    {
      const unsigned long current = values.size() * q ;
      std::nth_element  ( start , values.begin () + current , values.end () ) ;
      start = values.begin() + current ;
      result.push_back ( *start ) ;
    }
    return  Ostap::StatVar::Quantiles ( result , values.size() ) ;
  }
  // ==========================================================================
  /// get approximate  quantiles using P^2 algorithm  
  Ostap::StatVar::Quantiles
  _p2quantiles_ 
  ( Ostap::FrameNode        frame     , 
    const std::set<double>& quantiles , 
    const std::string&      expr      , 
    const std::string&      cuts      ) 
  {
    const bool no_cuts = trivial (  cuts ) ; 
    //
    const std::string var   = Ostap::tmp_name ( "v_"   , expr ) ;
    const std::string bcut  = Ostap::tmp_name ( "b_"   , cuts ) ;
    //
    std::vector<Ostap::Math::GSL::P2Quantile> qs ( quantiles.begin() , quantiles.end() ) ;
    // decorate the frame 
    auto t = frame
      .Define   ( bcut    , no_cuts ? "true" : "(bool) ( " + cuts + " ) ;" ) 
      .Filter   ( bcut    )    
      .Define   ( var     , "1.0*(" + expr + ")" ) ;
    auto l = t.Count() ;
    t.Foreach  (  [&qs] ( double v ) { for  ( auto&q : qs ) { q.add ( v ) ; } } , { var } ) ;
    // 
    return  Ostap::StatVar::Quantiles ( std::vector<double>( qs.begin() , qs.end() ) , *l ) ;
  }
  // ==========================================================================
}
// ============================================================================
/**  get quantile of the distribution  
 *   @param frame  (INPUT) the input data
 *   @param q      (INPUT) quantile value   0 < q < 1  
 *   @param expr   (INPUT) the expression 
 *   @param cuts   (INPUT) selection cuts 
 *   @return the quantile value 
 */
// ============================================================================
Ostap::StatVar::Quantile
Ostap::StatVar::quantile
( Ostap::FrameNode    frame  ,
  const double        q      , //  0<q<1 
  const std::string&  expr   , 
  const std::string&  cuts   ) 
{
  Ostap::Assert ( 0 < q && q < 1             , 
                  "Invalid quantile"         ,
                  "Ostap::StatVar::quantile" ) ;
  //
  std::set<double> qset {} ;
  qset.insert ( q )  ;
  //
  auto result = _quantiles_ ( frame , qset , expr , cuts ) ;
  Ostap::Assert ( 1 == result.quantiles.size()         , 
                  "Invalid quantiles size"   ,
                  "Ostap::StatVar::interval" ) ;
  return Ostap::StatVar::Quantile ( result.quantiles[0] , result.nevents ) ;
}
// ============================================================================
/**  get approximate quantile of the distribution using  P^2 algorithm  
 *   @param frame  (INPUT) the input data
 *   @param q      (INPUT) quantile value   0 < q < 1  
 *   @param expr   (INPUT) the expression 
 *   @param cuts   (INPUT) selection cuts 
 *   @return the quantile value 
 */
// ============================================================================
Ostap::StatVar::Quantile
Ostap::StatVar::p2quantile
( Ostap::FrameNode    frame  ,
  const double        q      , //  0<q<1 
  const std::string&  expr   , 
  const std::string&  cuts   ) 
{
  Ostap::Assert ( 0 < q && q < 1             , 
                  "Invalid quantile"         ,
                  "Ostap::StatVar::quantile" ) ;

  std::set<double> qset {} ;
  qset.insert ( q )  ;
  //
  auto result = _p2quantiles_ ( frame , qset , expr , cuts ) ;
  Ostap::Assert ( 1 == result.quantiles.size()         , 
                  "Invalid quantiles size"   ,
                  "Ostap::StatVar::interval" ) ;
  return Ostap::StatVar::Quantile ( result.quantiles[0] , result. nevents ) ;
}
// ========================================================================    
/*  get quantiles of the distribution  
 *   @param frame (INPUT) the input frame
 *   @param q     (INPUT) quantile value   0 < q < 1  
 *   @param expr  (INPUT) the expression 
 *   @param cuts  (INPUT) selection cuts 
 *   @return the quantile value 
 */
// ========================================================================    
Ostap::StatVar::Quantiles
Ostap::StatVar::quantiles
( Ostap::FrameNode           frame     ,
  const std::vector<double>& quantiles , 
  const std::string&         expr      , 
  const std::string&         cuts      ) 
{
  Ostap::Assert ( 1 <= quantiles.size()          , 
                  "Invalid vector of quantiles"  ,
                  "Ostap::StatVar::quantile"     ) ;
  std::set<double> qs ;
  for ( double v : quantiles ) { qs.insert ( v ) ; }
  Ostap::Assert ( !qs.empty ()                ,
                  "Invalid quantiles"         ,
                  "Ostap::StatVar::quantiles" ) ;
  Ostap::Assert ( 0 < *qs. begin ()           , 
                  "Invalid quantile"          ,
                  "Ostap::StatVar::quantiles" ) ;  
  Ostap::Assert ( 1 > *qs.rbegin ()           , 
                  "Invalid quantile"          ,
                  "Ostap::StatVar::quantiles" ) ;
  //
  return _quantiles_ ( frame , qs , expr , cuts ) ;
}
// ========================================================================    
/*   get approximate  quantiles of the distribution  using P^2 algorithm 
 *   @param frame (INPUT) the input frame
 *   @param q     (INPUT) quantile value   0 < q < 1  
 *   @param expr  (INPUT) the expression 
 *   @param cuts  (INPUT) selection cuts 
 *   @return the quantile value 
 */
// ========================================================================    
Ostap::StatVar::Quantiles
Ostap::StatVar::p2quantiles
( Ostap::FrameNode           frame     ,
  const std::vector<double>& quantiles , 
  const std::string&         expr      , 
  const std::string&         cuts      ) 
{
  Ostap::Assert ( 1 <= quantiles.size()          , 
                  "Invalid vector of quantiles"  ,
                  "Ostap::StatVar::quantile"     ) ;
  std::set<double> qs ;
  for ( double v : quantiles ) { qs.insert ( v ) ; }
  Ostap::Assert ( !qs.empty ()                ,
                  "Invalid quantiles"         ,
                  "Ostap::StatVar::quantiles" ) ;
  Ostap::Assert ( 0 < *qs. begin ()           , 
                  "Invalid quantile"          ,
                  "Ostap::StatVar::quantiles" ) ;  
  Ostap::Assert ( 1 > *qs.rbegin ()           , 
                  "Invalid quantile"          ,
                  "Ostap::StatVar::quantiles" ) ;
  //
  return _p2quantiles_ ( frame , qs , expr , cuts ) ;
}
// ============================================================================
/* Get the interval of the distribution  
 * @param tree  (INPUT) the input tree 
 * @param q1    (INPUT) quantile value   0 < q1 < 1  
 * @param q2    (INPUT) quantile value   0 < q2 < 1  
 * @param expr  (INPUT) the expression 
 * @param cuts  (INPUT) selection cuts 
 * @return the quantile value 
 * @code
 * FRAME& frame = ... ;
 * /// get 90% interval:
 * Interval ab = interval ( frame , 0.05 , 0.95 , 'mass' , 'pt>3' ) ;
 * @code 
 */
// ============================================================================
Ostap::StatVar::QInterval 
Ostap::StatVar::interval 
( Ostap::FrameNode    frame ,
  const double        q1    , //  0<q1<1 
  const double        q2    , //  0<q2<1 
  const std::string&  expr  , 
  const std::string&  cuts  )
{
  Ostap::Assert ( 0 < q1 && q1 < 1 , 
                  "Invalid quantile1"        ,
                  "Ostap::StatVar::interval" ) ;
  Ostap::Assert ( 0 < q2 && q2 < 1 , 
                  "Invalid quantile2"        ,
                  "Ostap::StatVar::interval" ) ;
  //
  std::set<double> qset {} ;
  qset.insert ( q1 )  ;
  qset.insert ( q2 )  ;
  //
  auto result = _quantiles_( frame , qset ,  expr , cuts ) ; 
  Ostap::Assert ( 2 == result.quantiles.size()         ,
                  "Invalid interval"         ,
                  "Ostap::StatVar::interval" ) ;
  //
  return QInterval ( Interval ( result.quantiles[0] , result.quantiles[1] ) , result.nevents ) ;
}
// ============================================================================
/* Get the approximate  interval of the distribution  using P^2 algorithm
 * @param tree  (INPUT) the input tree 
 * @param q1    (INPUT) quantile value   0 < q1 < 1  
 * @param q2    (INPUT) quantile value   0 < q2 < 1  
 * @param expr  (INPUT) the expression 
 * @param cuts  (INPUT) selection cuts 
 * @return the quantile value 
 * @code
 * FRAME& frame = ... ;
 * /// get 90% interval:
 * Interval ab = interval ( frame , 0.05 , 0.95 , 'mass' , 'pt>3' ) ;
 * @code 
 */
// =======================x=====================================================
Ostap::StatVar::QInterval 
Ostap::StatVar::p2interval 
( Ostap::FrameNode    frame ,
  const double        q1    , //  0<q1<1 
  const double        q2    , //  0<q2<1 
  const std::string&  expr  , 
  const std::string&  cuts  )
{
  Ostap::Assert ( 0 < q1 && q1 < 1 , 
                  "Invalid quantile1"        ,
                  "Ostap::StatVar::interval" ) ;
  Ostap::Assert ( 0 < q2 && q2 < 1 , 
                  "Invalid quantile2"        ,
                  "Ostap::StatVar::interval" ) ;
  //
  std::set<double> qset {} ;
  qset.insert ( q1 )  ;
  qset.insert ( q2 )  ;
  //
  auto result = _p2quantiles_ ( frame , qset ,  expr , cuts ) ; 
  Ostap::Assert ( 2 == result.quantiles.size()         ,
                  "Invalid interval"         ,
                  "Ostap::StatVar::interval" ) ;
  //
  return QInterval ( Interval ( result.quantiles[0] , result.quantiles[1] ) , result.nevents ) ;
}
// ============================================================================


// ============================================================================
/*  get the moment as Ostap::Math::Moment_<N>
 *  @see Ostap::Math::Moment_
 *  @see Ostap::Math::Moment
 */
// ============================================================================
Ostap::StatusCode 
Ostap::StatVar::the_moment
( TTree*                  tree       , 
  Ostap::Math::Statistic& moment     , 
  const std::string&      expression , 
  const unsigned long     first      ,
  const unsigned long     last       ) 
{
  if ( nullptr == tree    ) { return INVALID_DATA     ; }
  //
  Ostap::Formula formula ( expression , tree ) ;
  if ( !formula.ok()      ) { return INVALID_FORMULA ; }
  //
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  //  
  Ostap::Utils::Notifier notify ( tree , &formula ) ;
  //
  const unsigned long nEntries =
    std::min ( last , (unsigned long) tree->GetEntries() ) ;
  //
  std::vector<double>  results {} ;
  for ( unsigned long entry = first ; entry < nEntries ; ++entry )
  {
    //
    long ievent = tree->GetEntryNumber ( entry ) ;
    if ( 0 > ievent ) { return INVALID_ENTRY ; }
    //
    ievent      = tree->LoadTree ( ievent ) ;
    if ( 0 > ievent ) { return INVALID_EVENT  ; } 
    //
    formula.evaluate ( results ) ;
    for  ( const double r : results ) { moment.update ( r ) ; }
  }
  //
  return Ostap::StatusCode::SUCCESS ;
}
// ========================================================================    
/*  get the moment as Ostap::Math::WMoment_<N>
 *  @see Ostap::Math::WMoment_
 *  @see Ostap::Math::WMoment
 */
// ============================================================================
Ostap::StatusCode 
Ostap::StatVar::the_moment
( TTree*                   tree       , 
  Ostap::Math::WStatistic& moment     , 
  const std::string&       expression , 
  const std::string&       selection  , 
  const unsigned long      first      ,
  const unsigned long      last       ) 
{
  if ( nullptr == tree    ) { return Ostap::StatusCode ( INVALID_DATA       ) ; }
  //
  Ostap::Formula formula ( expression , tree ) ;
  if ( !formula.ok()      ) { return INVALID_FORMULA ; }
  //
  std::unique_ptr<Ostap::Formula> cuts { nullptr } ;
  if  ( !selection.empty() ) 
  { 
    cuts = std::make_unique<Ostap::Formula>( selection , tree ) ; 
    if ( !cuts || !cuts->ok() ) { return INVALID_FORMULA ; }
  }
  //
  const bool with_cuts = !(!cuts) ;
  //
  if ( last <= first      ) { return Ostap::StatusCode::SUCCESS ; }
  //  
  Ostap::Utils::Notifier notify ( tree , &formula ) ;
  //
  const unsigned long nEntries =
    std::min ( last , (unsigned long) tree->GetEntries() ) ;
  //
  std::vector<double>  results {} ;
  for ( unsigned long entry = first ; entry < nEntries ; ++entry )
  {
    //
    long ievent = tree->GetEntryNumber ( entry ) ;
    if ( 0 > ievent ) { return INVALID_ENTRY ; }
    //
    ievent      = tree->LoadTree ( ievent ) ;
    if ( 0 > ievent ) { return INVALID_EVENT ; } 
    //
    const long double w = with_cuts ? cuts->evaluate() : 1.0L ;
    //
    if ( !w ) { continue ; } // ATTENTION! 
    //
    formula.evaluate ( results ) ;
    for  ( const double r : results ) { moment.update ( r , w ) ; }
  }
  //
  return Ostap::StatusCode::SUCCESS ;
}
// ========================================================================
/*  get the moment as Ostap::Math::WMoment_<N>
 *  @see Ostap::Math::WMoment_
 *  @see Ostap::Math::WMoment
 */
// ========================================================================
Ostap::StatusCode 
Ostap::StatVar::the_moment
( const RooAbsData*        data       , 
  Ostap::Math::WStatistic& moment     , 
  const std::string&       expression , 
  const std::string&       selection  , 
  const std::string&       cut_range  , 
  const unsigned long      first      ,
  const unsigned long      last       ) 
{
  if ( nullptr == data ) { return INVALID_DATA ; }
  //
  const std::unique_ptr<Ostap::FormulaVar> expr { ::make_formula ( expression , *data , false , true ) } ;
  const std::unique_ptr<Ostap::FormulaVar> cuts { ::make_formula ( selection  , *data , true  , true ) } ;
  //
  if ( !expr || !expr->ok() )                          { return INVALID_FORMULA ; }
  if ( selection.empty() && ( !cuts || !cuts->ok() ) ) { return INVALID_FORMULA  ; }
  //
  if ( last <= first ) { return Ostap::StatusCode::SUCCESS ; }
  //
  const char* cutrange = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  const bool  weighted = data->isWeighted() ;
  //
  const unsigned long the_last = std::min ( last , (unsigned long) data->numEntries() ) ;
  //
  // start the loop
  for ( unsigned long entry = first ; entry < the_last ; ++entry )
  {
    //
    const RooArgSet* vars = data->get( entry ) ;
    if ( nullptr == vars  ) { return Ostap::StatusCode ( INVALID_ENTRY ) ; }
    if ( cutrange && !vars->allInRange ( cutrange ) ) { continue ; } // CONTINUE    
    //
    // apply cuts:
    const long double wc = cuts      ? cuts->getVal() : 1.0L ;
    if ( !wc ) { continue ; }                                   // CONTINUE  
    // apply weight:
    const long double wd = weighted  ? data->weight()        : 1.0L ;
    if ( !wd ) { continue ; }                                   // CONTINUE    
    // cuts & weight:
    const long double w  = wd *  wc ;
    if ( !w  ) { continue ; }                                   // CONTINUE        
    //
    const double v = expr->getVal () ;
    //
    moment.update ( v , w ) ;
  }
  //
  return Ostap::StatusCode::SUCCESS ;
  // ==========================================================================
}
// ========================================================================
/*  get the moment as Ostap::Math::WMoment_<N>
 *  @see Ostap::Math::WMoment_
 *  @see Ostap::Math::WMoment
 */
// ========================================================================
Ostap::StatusCode 
Ostap::StatVar::the_moment
( const RooAbsData*        data       , 
  Ostap::Math::WStatistic& moment     , 
  const std::string&       expression , 
  const std::string&       selection  , 
  const unsigned long      first      ,
  const unsigned long      last       ) 
{
  return the_moment ( data       , 
                      moment     , 
                      expression , 
                      selection  ,
                      ""         , 
                      first      , 
                      last       ) ;
}
// ========================================================================
/*  get the moment as Ostap::Math::WMoment_<N>
 *  @see Ostap::Math::WMoment_
 *  @see Ostap::Math::WMoment
 */
// ========================================================================
Ostap::StatusCode 
Ostap::StatVar::the_moment
( const RooAbsData*        data       , 
  Ostap::Math::WStatistic& moment     , 
  const std::string&       expression , 
  const unsigned long      first      ,
  const unsigned long      last       ) 
{
  return the_moment ( data       , 
                      moment     , 
                      expression ,
                      ""         , 
                      ""         , 
                      first      , 
                      last       ) ;
}
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
( TTree*              data       ,
  Ostap::Math::ECDF&  ecdf       ,
  const std::string&  expression ,
  const unsigned long first      ,
  const unsigned long last       )
{
  // reset ECDF 
  ecdf = Ostap::Math::ECDF () ;
  const Ostap::StatusCode sc = the_moment ( data , ecdf , expression , first , last ) ;
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
( TTree*              data       ,
  Ostap::Math::WECDF& ecdf       ,
  const std::string&  expression , 
  const std::string&  selection  , 
  const unsigned long first      ,
  const unsigned long last       )
{
  // reset ECDF 
  ecdf = Ostap::Math::WECDF () ;
  const Ostap::StatusCode sc = the_moment ( data , ecdf , expression , selection , first , last ) ;
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
( const RooAbsData*   data       ,
  Ostap::Math::WECDF& ecdf       ,
  const std::string&  expression , 
  const std::string&  selection  , 
  const std::string&  cut_range  , 
  const unsigned long first      ,
  const unsigned long last       )
{
  // reset ECDF 
  ecdf = Ostap::Math::WECDF () ;
  const Ostap::StatusCode sc = the_moment ( data , ecdf , expression , selection , cut_range , first , last ) ;
  if       ( sc.isFailure() ) { return sc ; }
  else if  ( !ecdf.ok()     ) { return INVALID_WECDF ; } 
  return Ostap::StatusCode::SUCCESS ; 
}




// ============================================================================
/** Get the empirical cumulative distribtion function 
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
( Ostap::FrameNode    data        ,
  Ostap::Math::ECDF&  ecdf        ,
  const std::string&  expression  )
{
  // (1) reset ECDF 
  ecdf = Ostap::Math::ECDF () ;
  //
  /// define the temporary columns 
  const std::string var    { Ostap::tmp_name ( "v_" , expression ) } ;
  auto t = data.Define ( var    ,  "1.0*(" + expression + ")" ) ;
  //
  const unsigned int nSlots = std::max ( 1u , Ostap::Utils::mt_pool_size() ) ;
  std::vector<Ostap::Math::ECDF> _stat ( nSlots ? nSlots : 1 ) ;
  //
  auto fun = [&_stat,nSlots] ( unsigned int slot , double v ) 
  { _stat [ slot % nSlots ].add ( v ) ; } ;
  //
  t.ForeachSlot ( fun ,  { var } ) ;
  //
  // merge results 
  for ( unsigned int i = 1 ; i < nSlots ; ++i )
    { _stat [ 0 ] += _stat [ i ] ; }
  //
  if ( !_stat [ 0 ].ok () ) {  return INVALID_ECDF ; }
  ecdf = _stat [ 0 ] ;
  //
  return Ostap::StatusCode::SUCCESS ; 
}
// ============================================================================
/** Get the empirical cumulative distribtion function 
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
( Ostap::FrameNode    data        ,
  Ostap::Math::WECDF& ecdf        ,
  const std::string&  expression  ,
  const std::string&  selection   )
{
  //
  const bool no_cuts = trivial ( selection ) ; 
  /// define the temporary columns 
  const std::string var    { Ostap::tmp_name ( "v_" , expression ) } ;
  const std::string weight { Ostap::tmp_name ( "w_" , selection  ) } ;
  const std::string bcut   { Ostap::tmp_name ( "b_" , selection  ) } ;
  //
  auto t = data
  .Define ( bcut   , no_cuts ? "true" : "(bool)   ( " + selection + " ) ;" ) 
  .Filter ( bcut   ) 
  .Define ( var    ,  "1.0*(" + expression + ")"   )
  .Define ( weight , no_cuts ? "1.0"  : "1.0*(" + selection + ")" ) ;
  //
  const unsigned int nSlots = std::max ( 1u , Ostap::Utils::mt_pool_size() ) ;
  std::vector<Ostap::Math::WECDF> _stat ( nSlots ? nSlots : 1 ) ;
  //  
  auto fun = [&_stat,nSlots] ( unsigned int slot , double v , double w ) 
  { _stat [ slot % nSlots ].add ( v , w ) ; } ;
  //
  t.ForeachSlot ( fun ,  { var , weight } ) ; 
  //
  // merge results 
  for ( unsigned int i = 1 ; i < nSlots ; ++i )
    { _stat [ 0 ] += _stat [ i ] ; }
  //
  if ( !_stat [ 0 ].ok () ) {  return INVALID_WECDF ; }
  ecdf = _stat [ 0 ] ;
  //
  return Ostap::StatusCode::SUCCESS ;   
}
// ============================================================================
//                                                                      The END
// ============================================================================
