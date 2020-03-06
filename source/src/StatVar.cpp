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
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Formula.h"
#include "Ostap/Iterator.h"
#include "Ostap/Notifier.h"
#include "Ostap/MatrixUtils.h"
#include "Ostap/StatVar.h"
// ============================================================================
// ROOT
// ============================================================================
#include "TTree.h"
#include "TCut.h"
#include "RooDataSet.h"
// ============================================================================
// Local
// ============================================================================
#include "OstapDataFrame.h"
#include "Exception.h"
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
  /// make RooFormulaVar 
  std::unique_ptr<RooFormulaVar>
  make_formula ( const std::string& expression           , 
                 const RooAbsData&  data                 , 
                 const bool         allow_empty = false  ) 
  { 
    if ( allow_empty && expression.empty() ) { return nullptr ; }  // RETURN!
    //
    RooArgList        alst ;
    const RooArgSet*  aset = data.get() ;
    Ostap::Assert ( nullptr != aset                ,  
                    "Invalid varset"               , 
                    "Ostap::StatVar::make_formula" ) ;
    //
    Ostap::Utils::Iterator iter ( *aset );
    //
    RooAbsArg* coef = 0 ;
    while ( ( coef = (RooAbsArg*) iter.next() ) ) { alst.add ( *coef ); }
    //
    const auto fname = Ostap::tmp_name ( "formula_" , expression ) ;
    auto result = std::make_unique<RooFormulaVar> ( fname.c_str () , expression.c_str() , alst , false ) ;
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
  std::vector<double> 
  _quantiles_
  ( TTree&                  tree      ,
    const std::set<double>& quantiles , //  0<q<1 
    Ostap::Formula&      var       ,
    Ostap::Formula*      cuts      , 
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
    if ( 0 == num ) { return std::vector<double>() ; }
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
    return result ;
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
  std::vector<double> 
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
    return result ;
  }
  // ==========================================================================
} //                                                 end of anonymous namespace
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
Ostap::StatVar::Statistic
Ostap::StatVar::statVar
( TTree*              tree       ,
  const std::string&  expression ,
  const unsigned long first      ,
  const unsigned long last       )
{
  Statistic result ;
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
Ostap::StatVar::Statistic
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
  Ostap::StatVar::Statistic result ;
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
Ostap::StatVar::Statistic
Ostap::StatVar::statVar
( TTree*              tree       ,
  const std::string&  expression ,
  const TCut&         cuts       ,
  const unsigned long first      ,
  const unsigned long last       )
{
  const std::string _cuts = cuts.GetTitle() ;
  return statVar ( tree , expression , _cuts , first , last ) ;
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
  std::vector<Ostap::StatVar::Statistic>& result      ,  
  const std::vector<std::string>&         expressions ,
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
  std::vector<Ostap::StatVar::Statistic>& result      ,  
  const std::vector<std::string>&         expressions ,
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
/** build statistic for the <code>expressions</code>
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
  std::vector<Ostap::StatVar::Statistic>& result      ,  
  const std::vector<std::string>&         expressions ,
  const TCut&                             cuts        ,
  const unsigned long                     first       ,
  const unsigned long                     last        ) 
{
  const std::string _cuts = cuts.GetTitle() ;
  return statVars ( tree , result , expressions , _cuts , first , last ) ;
}
// ============================================================================
/*  calculate the covariance of two expressions
 *  @param tree  (INPUT)  the input tree
 *  @param exp1  (INPUT)  the first  expresiion
 *  @param exp2  (INPUT)  the second expresiion
 *  @param stat1 (UPDATE) the statistic for the first  expression
 *  @param stat2 (UPDATE) the statistic for the second expression
 *  @param cov2  (UPDATE) the covariance matrix
 *  @return number of processed events
 *
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2014-03-27
 */
// ============================================================================
unsigned long
Ostap::StatVar::statCov
( TTree*               tree    ,
  const std::string&   exp1    ,
  const std::string&   exp2    ,
  Ostap::StatVar::Statistic& stat1 ,
  Ostap::StatVar::Statistic& stat2 ,
  Ostap::SymMatrix2x2& cov2    ,
  const unsigned long  first   ,
  const unsigned long  last    )
{
  //
  stat1.reset () ;
  stat2.reset () ;
  Ostap::Math::setToScalar ( cov2 , 0.0 ) ;
  //
  if ( 0 == tree || last <= first ) { return 0 ; }         // RETURN
  Ostap::Formula formula1 ( exp1 , tree ) ;
  if ( !formula1 .ok () ) { return 0 ; }                   // RETURN
  Ostap::Formula formula2 ( exp2 , tree ) ;
  if ( !formula2 .ok () ) { return 0 ; }                   // RETURN
  //
  Ostap::Utils::Notifier notify ( tree , &formula1 , &formula2 ) ;
  //
  const unsigned long nEntries =
    std::min ( last , (unsigned long) tree->GetEntries() ) ;
  //
  std::vector<double> results1 {} ;
  std::vector<double> results2 {} ;
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
        //
        stat1 += v1 ;
        stat2 += v2 ;
        //
        cov2 ( 0 , 0 ) += v1*v1 ;
        cov2 ( 0 , 1 ) += v1*v2 ;
        cov2 ( 1 , 1 ) += v2*v2 ;
      }
    }
    //
  }
  //
  if ( 0 == stat1.nEntries() ) { return 0 ; }          // RETURN
  //
  cov2 /= stat1.nEntries () ;
  //
  const double v1_mean = stat1.mean() ;
  const double v2_mean = stat2.mean() ;
  //
  cov2 ( 0 , 0 ) -= v1_mean * v1_mean ;
  cov2 ( 0 , 1 ) -= v1_mean * v2_mean ;
  cov2 ( 1 , 1 ) -= v2_mean * v2_mean ;
  //
  return stat1.nEntries () ;
}
// ============================================================================
/*  calculate the covariance of two expressions
 *  @param tree  (INPUT)  the inpout tree
 *  @param exp1  (INPUT)  the first  expresiion
 *  @param exp2  (INPUT)  the second expresiion
 *  @param cuts  (INPUT)  the selection criteria
 *  @param stat1 (UPDATE) the statistic for the first  expression
 *  @param stat2 (UPDATE) the statistic for the second expression
 *  @param cov2  (UPDATE) the covariance matrix
 *  @return number of processed events
 *
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2014-03-27
 */
// ============================================================================
unsigned long
Ostap::StatVar::statCov
( TTree*               tree    ,
  const std::string&   exp1    ,
  const std::string&   exp2    ,
  const std::string&   cuts    ,
  Ostap::StatVar::Statistic& stat1 ,
  Ostap::StatVar::Statistic& stat2 ,
  Ostap::SymMatrix2x2& cov2    ,
  const unsigned long  first   ,
  const unsigned long  last    )
{
  if ( cuts.empty() ) 
  { return statCov ( tree , exp1 , exp2 , stat1 , stat2 , cov2 , first , last ) ; }
  //
  stat1.reset () ;
  stat2.reset () ;
  Ostap::Math::setToScalar ( cov2 , 0.0 ) ;
  //
  if ( 0 == tree || last <= first ) { return 0 ; }              // RETURN
  Ostap::Formula formula1 ( exp1 , tree ) ;
  if ( !formula1 .ok () ) { return 0 ; }                        // RETURN
  Ostap::Formula formula2 ( exp2 , tree ) ;
  if ( !formula2 .ok () ) { return 0 ; }                        // RETURN
  Ostap::Formula selection ( cuts      , tree ) ;
  if ( !selection.ok () ) { return 0 ; }                        // RETURN
  //
  Ostap::Utils::Notifier notify ( tree , &formula1 , &formula2 ) ;
  //
  const unsigned long nEntries =
    std::min ( last , (unsigned long) tree->GetEntries() ) ;
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
    const double w = selection.evaluate() ;
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
        //
        stat1.add ( v1 , w ) ;
        stat2.add ( v2 , w ) ;
        //
        cov2 ( 0 , 0 ) += w*v1*v1 ;
        cov2 ( 0 , 1 ) += w*v1*v2 ;
        cov2 ( 1 , 1 ) += w*v2*v2 ;
      }
    }
  }
  //
  if ( 0 == stat1.nEntries() || 0 == stat1.nEff () ) { return 0 ; }
  //
  cov2 /= stat1.weights().sum()  ;
  //
  const double v1_mean = stat1.mean() ;
  const double v2_mean = stat2.mean() ;
  //
  cov2 ( 0 , 0 ) -= v1_mean * v1_mean ;
  cov2 ( 0 , 1 ) -= v1_mean * v2_mean ;
  cov2 ( 1 , 1 ) -= v2_mean * v2_mean ;
  //
  return stat1.nEntries() ;
}
// ============================================================================
/*  calculate the covariance of two expressions
 *  @param tree  (INPUT)  the inpout tree
 *  @param exp1  (INPUT)  the first  expresiion
 *  @param exp2  (INPUT)  the second expresiion
 *  @param cuts  (INPUT)  the selection criteria
 *  @param stat1 (UPDATE) the statistic for the first  expression
 *  @param stat2 (UPDATE) the statistic for the second expression
 *  @param cov2  (UPDATE) the covariance matrix
 *  @return number of processed events
 *
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2014-03-27
 */
// ============================================================================
unsigned long
Ostap::StatVar::statCov
( TTree*               tree    ,
  const std::string&   exp1    ,
  const std::string&   exp2    ,
  const TCut&          cuts    ,
  Ostap::StatVar::Statistic& stat1 ,
  Ostap::StatVar::Statistic& stat2 ,
  Ostap::SymMatrix2x2& cov2    ,
  const unsigned long  first   ,
  const unsigned long  last    )
{
  const std::string _cuts = cuts.GetTitle() ;
  //
  return statCov ( tree  ,
                   exp1  , exp2    , _cuts ,
                   stat1 , stat2   , cov2  ,
                   first , last    ) ;
}
// ============================================================================

// ============================================================================
Ostap::StatVar::Statistic
Ostap::StatVar::statVar
( const RooAbsData*   data        ,
  const std::string&  expression  ,
  const TCut&         cuts        ,
  const std::string&  cut_range   ,
  const unsigned long first       ,
  const unsigned long last        )
{
  const std::string _cuts = cuts.GetTitle() ;
  return statVar ( data , expression , _cuts , cut_range , first , last ) ;
}
// ============================================================================
Ostap::StatVar::Statistic
Ostap::StatVar::statVar
( const RooAbsData*   data        ,
  const std::string&  expression  ,
  const TCut&         cuts        ,
  const unsigned long first       ,
  const unsigned long last        )
{
  const std::string _cuts = cuts.GetTitle() ;
  return statVar ( data , expression , _cuts , "" , first , last ) ;
}
// ============================================================================
Ostap::StatVar::Statistic
Ostap::StatVar::statVar
( const RooAbsData*   data        ,
  const std::string&  expression  ,
  const std::string&  cuts        ,
  const std::string&  cut_range   ,
  const unsigned long first       ,
  const unsigned long last        )
{
  Statistic result ;
  if ( 0 == data || last <= first ) { return result ; }         // RETURN
  //
  const std::unique_ptr<RooFormulaVar> formula   { make_formula ( expression , *data        ) } ;
  const std::unique_ptr<RooFormulaVar> selection { make_formula ( cuts       , *data , true ) } ;
  //
  const bool  weighted = data->isWeighted() ;
  const char* cutrange = cut_range.empty() ?  nullptr : cut_range.c_str() ;
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
    const double v = formula->getVal () ;
    //
    result.add ( v , w ) ;
  }
  //
  return result ;
}
// ============================================================================
/*  calculate the covariance of two expressions
 *  @param tree  (INPUT)  the input tree
 *  @param exp1  (INPUT)  the first  expresiion
 *  @param exp2  (INPUT)  the second expresiion
 *  @param stat1 (UPDATE) the statistic for the first  expression
 *  @param stat2 (UPDATE) the statistic for the second expression
 *  @param cov2  (UPDATE) the covariance matrix
 *  @return number of processed events
 *
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2014-03-27
 */
// ============================================================================
unsigned long
Ostap::StatVar::statCov
( const RooAbsData*             data      ,
  const std::string&            exp1      ,
  const std::string&            exp2      ,
  Ostap::StatVar::Statistic& stat1     ,
  Ostap::StatVar::Statistic& stat2     ,
  Ostap::SymMatrix2x2&          cov2      ,
  const std::string&            cut_range ,
  const unsigned long           first     ,
  const unsigned long           last      )
{
  //
  stat1.reset () ;
  stat2.reset () ;
  Ostap::Math::setToScalar ( cov2 , 0.0 ) ;
  //
  if ( 0 == data || last <= first ) { return 0 ; }               // RETURN
  //
  const std::unique_ptr<RooFormulaVar> formula1 { make_formula ( exp1 , *data ) } ;
  const std::unique_ptr<RooFormulaVar> formula2 { make_formula ( exp2 , *data ) } ;
  //
  const bool  weighted = data->isWeighted() ;
  const char* cutrange = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  //
  const unsigned long nEntries = std::min ( last , (unsigned long) data->numEntries() ) ;
  //
  for ( unsigned long entry = first ; entry < nEntries ; ++entry )
  {
    //
    const RooArgSet* vars = data->get( entry ) ;
    if ( nullptr == vars  )                             { break    ; } // BREAK
    if ( cutrange && !vars->allInRange ( cutrange ) ) { continue ; } // CONTINUE    
    //
    // apply weight:
    const long double w = weighted  ? data->weight()        : 1.0L ;
    if ( !w ) { continue ; }                                   // CONTINUE    
    //
    const double v1 = formula1->getVal() ;
    const double v2 = formula2->getVal() ;
    //
    stat1.add ( v1 , w ) ;
    stat2.add ( v2 , w ) ;
    //
    cov2 ( 0 , 0 ) += w * v1 * v1 ;
    cov2 ( 0 , 1 ) += w * v1 * v2 ;
    cov2 ( 1 , 1 ) += w * v2 * v2 ;
    //
  }
  //
  if ( 0 == stat1.nEntries() || 0 == stat1.nEff () ) { return 0 ; }
  //
  cov2 /= stat1.weights().sum()  ;
  //
  const double v1_mean = stat1.mean() ;
  const double v2_mean = stat2.mean() ;
  //
  cov2 ( 0 , 0 ) -= v1_mean * v1_mean ;
  cov2 ( 0 , 1 ) -= v1_mean * v2_mean ;
  cov2 ( 1 , 1 ) -= v2_mean * v2_mean ;
  //
  return stat1.nEntries() ;
}
// ============================================================================
/*  calculate the covariance of two expressions
 *  @param tree  (INPUT)  the input tree
 *  @param exp1  (INPUT)  the first  expresiion
 *  @param exp2  (INPUT)  the second expresiion
 *  @param cuts  (INPUT)  selection
 *  @param stat1 (UPDATE) the statistic for the first  expression
 *  @param stat2 (UPDATE) the statistic for the second expression
 *  @param cov2  (UPDATE) the covariance matrix
 *  @return number of processed events
 *
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2014-03-27
 */
// ============================================================================
unsigned long
Ostap::StatVar::statCov
( const RooAbsData*             data      ,
  const std::string&            exp1      ,
  const std::string&            exp2      ,
  const std::string&            cuts      ,
  Ostap::StatVar::Statistic& stat1     ,
  Ostap::StatVar::Statistic& stat2     ,
  Ostap::SymMatrix2x2&          cov2      ,
  const std::string&            cut_range ,
  const unsigned long           first     ,
  const unsigned long           last      )
{
  //
  stat1.reset () ;
  stat2.reset () ;
  Ostap::Math::setToScalar ( cov2 , 0.0 ) ;
  //
  if ( 0 == data || last <= first ) { return 0 ; }               // RETURN
  //
  const std::unique_ptr<RooFormulaVar> formula1  { make_formula ( exp1 , *data         ) } ;
  const std::unique_ptr<RooFormulaVar> formula2  { make_formula ( exp2 , *data         ) } ;
  const std::unique_ptr<RooFormulaVar> selection { make_formula ( cuts , *data ,  true ) } ;
  //
  const bool weighted = data->isWeighted() ;
  const char* cutrange = cut_range.empty() ?  nullptr : cut_range.c_str() ;
 //
  const unsigned long nEntries = std::min ( last , (unsigned long) data->numEntries() ) ;
  //
  for ( unsigned long entry = first ; entry < nEntries ; ++entry )
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
    const double v1 = formula1->getVal() ;
    const double v2 = formula2->getVal() ;
    //
    stat1.add ( v1 , w ) ;
    stat2.add ( v2 , w ) ;
    //
    cov2 ( 0 , 0 ) += w * v1 * v1 ;
    cov2 ( 0 , 1 ) += w * v1 * v2 ;
    cov2 ( 1 , 1 ) += w * v2 * v2 ;
    //
  }
  //
  if ( 0 == stat1.nEntries() || 0 == stat1.nEff () ) { return 0 ; }
  //
  cov2 /= stat1.weights().sum()  ;
  //
  const double v1_mean = stat1.mean() ;
  const double v2_mean = stat2.mean() ;
  //
  cov2 ( 0 , 0 ) -= v1_mean * v1_mean ;
  cov2 ( 0 , 1 ) -= v1_mean * v2_mean ;
  cov2 ( 1 , 1 ) -= v2_mean * v2_mean ;
  //
  return stat1.nEntries() ;
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
double Ostap::StatVar::quantile
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
  const std::vector<double> result = _quantiles_ 
    ( tree , std::set<double> {{ q }} , var , cut.get() , first , last ) ; 
  //
  Ostap::Assert ( 1 == result.size()         , 
                  "Invalid quantiles size"   ,
                  "Ostap::StatVar::interval" ) ;
  //
  return result[0] ;
}
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
std::vector<double> 
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
Ostap::StatVar::Interval 
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
  const std::vector<double> result = 
    _quantiles_ ( tree , std::set<double>{{ q1 , q2 }} , 
                  var  ,  cut.get() , first , last ) ; 
  Ostap::Assert ( 2 == result.size()         ,
                  "Invalid interval"         ,
                  "Ostap::StatVar::interval" ) ;
  //
  return std::make_pair( result[0] , result[1] ) ;
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
  const char* cutrange   = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  //
  // the simplest case :
  if ( !with_cuts && cut_range.empty() && !weighted ) { return the_last - first ; } // RETURN
  //
  const std::unique_ptr<RooFormulaVar> cut { make_formula ( cuts , data , true ) } ;
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
  const char* cutrange   = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  //
  const std::unique_ptr<RooFormulaVar> expression { make_formula ( expr , data        ) } ;
  const std::unique_ptr<RooFormulaVar> cut        { make_formula ( cuts , data , true ) } ;
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
  const char* cutrange   = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  //
  const std::unique_ptr<RooFormulaVar> expression { make_formula ( expr , data        ) } ;
  const std::unique_ptr<RooFormulaVar> cut        { make_formula ( cuts , data , true ) } ;
  //  
  long double mom   = 0 ;
  long double sumw  = 0 ; // sum of weights 
  // for uncertainties 
  long double sumw2 = 0 ; // sum of weights^2
  long double c2    = 0 ;
  //
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
  const char* cutrange   = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  //
  const std::unique_ptr<RooFormulaVar> expression { make_formula ( expr , data        ) } ;
  const std::unique_ptr<RooFormulaVar> cut        { make_formula ( cuts , data , true ) } ;
  //  
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
  const char* cutrange   = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  //
  const std::unique_ptr<RooFormulaVar> expression { make_formula ( expr , data        ) } ;
  const std::unique_ptr<RooFormulaVar> cut        { make_formula ( cuts , data , true ) } ;
  //
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
  const char* cutrange  = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  //
  const std::unique_ptr<RooFormulaVar> expression { make_formula ( expr , data        ) } ;
  const std::unique_ptr<RooFormulaVar> cut        { make_formula ( cuts , data , true ) } ;
  //
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
double Ostap::StatVar::quantile
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
  const char* cutrange  = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  //
  const std::unique_ptr<RooFormulaVar> expression { make_formula ( expr , data        ) } ;
  const std::unique_ptr<RooFormulaVar> cut        { make_formula ( cuts , data , true ) } ;
  //  
  const std::vector<double> result = _quantiles_ 
    (  data,  std::set<double>{{ q }} , 
       *expression ,  cut.get() , first , the_last , cutrange ) ;
  //
  Ostap::Assert ( 1 == result.size()         ,
                  "Invalid quantile size"    ,
                  "Ostap::StatVar::quantile" ) ;
  //
  return result[0] ;
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
Ostap::StatVar::Interval 
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
  if ( the_last <= first ) { return  std::make_pair(0.0,0.0) ; }    // RETURN
  //
  const char* cutrange  = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  //
  const std::unique_ptr<RooFormulaVar> expression { make_formula ( expr , data        ) } ;
  const std::unique_ptr<RooFormulaVar> cut        { make_formula ( cuts , data , true ) } ;
  //  
  const std::vector<double> result = _quantiles_ 
    (  data,  std::set<double>{{ q1 , q2 }} , 
       *expression ,  cut.get() , first , the_last , cutrange ) ;
  Ostap::Assert ( 2 ==  result.size()        ,
                  "Invalid quantile size"    ,
                  "Ostap::StatVar::quantile" ) ;
  //
  return std::make_pair( result[0] , result[1] ) ;
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
std::vector<double> 
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
  const char* cutrange  = cut_range.empty() ?  nullptr : cut_range.c_str() ;
  //
  const std::unique_ptr<RooFormulaVar> expression { make_formula ( expr , data        ) } ;
  const std::unique_ptr<RooFormulaVar> cut        { make_formula ( cuts , data , true ) } ;
  //  
  return _quantiles_ ( data  , qs  , 
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
( Ostap::DataFrame     frame ,
  const std::string&   cuts  )
{
  //
  const bool no_cuts = trivial ( cuts ) ;
  if ( no_cuts ) { return  *frame.Count(); }                     // RETURN
  //
  /// define temporary columns 
  const std::string weight  = Ostap::tmp_name ( "w_"  , cuts ) ;
  const std::string weight2 = Ostap::tmp_name ( "w2_" , cuts ) ;
  const std::string bcut    = Ostap::tmp_name ( "b_"  , cuts ) ;
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
Ostap::StatVar::Statistic 
Ostap::StatVar::statVar 
( Ostap::DataFrame    frame      , 
  const std::string&  expression , 
  const std::string&  cuts       ) 
{
  //
  const bool no_cuts = trivial ( cuts ) ; 
  //
  /// define the temporary columns 
  const std::string var    = Ostap::tmp_name ( "v_" , expression ) ;
  const std::string weight = Ostap::tmp_name ( "w_" , cuts       ) ;
  const std::string bcut   = Ostap::tmp_name ( "b_" , cuts       ) ;
  /// define actions 
  auto t = frame
    .Define ( bcut   , no_cuts ? "true" : "(bool)   ( " + cuts + " ) ;" ) 
    .Filter ( bcut   ) 
    .Define ( var    ,  "1.0*(" + expression + ")"   )
    .Define ( weight , no_cuts ? "1.0"  : "1.0*(" + cuts + ")" ) ;
  //
  const unsigned int nSlots = ROOT::GetImplicitMTPoolSize();
  std::vector<Statistic> _stat ( nSlots ? nSlots : 1 ) ;
  //
  auto fun = [&_stat] ( unsigned int slot , double v , double w ) 
    { _stat[slot].add ( v , w ) ; } ;
  t.ForeachSlot ( fun ,  { var , weight } ) ; 
  //
  Statistic stat ; for ( const auto& s : _stat ) { stat += s ; }
  return stat ;
}
// ============================================================================
/*  calculate the covariance of two expressions 
 *  @param frame (INPUT)  data frame 
 *  @param exp1  (INPUT)  the first  expresiion
 *  @param exp2  (INPUT)  the second expresiion
 *  @param stat1 (UPDATE) the statistic for the first  expression
 *  @param stat2 (UPDATE) the statistic for the second expression
 *  @param cov2  (UPDATE) the covariance matrix 
 *  @return number of processed events 
 *  
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2018-06-18
 */
// ============================================================================
unsigned long Ostap::StatVar::statCov
( Ostap::DataFrame     frame  , 
  const std::string&   exp1   , 
  const std::string&   exp2   , 
  Statistic&           stat1  ,  
  Statistic&           stat2  ,  
  Ostap::SymMatrix2x2& cov2   ) 
{ return statCov ( frame , exp1 , exp2 , "" , stat1 , stat2  ,  cov2 ) ; }
// ============================================================================
/*  calculate the covariance of two expressions 
 *  @param frame (INPUT)  data frame 
 *  @param exp1  (INPUT)  the first  expresiion
 *  @param exp2  (INPUT)  the second expresiion
 *  @param cuts  (INPUT)  the selection 
 *  @param stat1 (UPDATE) the statistic for the first  expression
 *  @param stat2 (UPDATE) the statistic for the second expression
 *  @param cov2  (UPDATE) the covariance matrix 
 *  @return number of processed events 
 *  
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2018-06-18
 */
// ============================================================================
unsigned long Ostap::StatVar::statCov
( Ostap::DataFrame     frame , 
  const std::string&   exp1  , 
  const std::string&   exp2  , 
  const std::string&   cuts  , 
  Statistic&           stat1 ,  
  Statistic&           stat2 ,  
  Ostap::SymMatrix2x2& cov2  ) 
{
  //
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
  const unsigned int nSlots = ROOT::GetImplicitMTPoolSize();
  std::vector<Statistic>           _sta1 ( nSlots ? nSlots : 1 ) ;
  std::vector<Statistic>           _sta2 ( nSlots ? nSlots : 1 ) ;
  std::vector<Ostap::SymMatrix2x2> _cov2 ( nSlots ? nSlots : 1 ) ;
  //
  auto fun = [&_sta1,&_sta2,&_cov2] 
    ( unsigned int slot , double v1 , double v2 , double w )  { 
    if ( w ) 
    {
      _sta1[slot].add ( v1 , w ) ; 
      _sta2[slot].add ( v2 , w ) ; 
      _cov2[slot]( 0 , 0 ) += w * v1 * v1 ;
      _cov2[slot]( 0 , 1 ) += w * v1 * v2 ;
      _cov2[slot]( 1 , 1 ) += w * v2 * v2 ;        
    }
  } ;
  t.ForeachSlot ( fun , { var1 , var2 , weight } ); 
  // 
  for ( const auto& s : _sta1 ) { stat1 += s ; }
  for ( const auto& s : _sta2 ) { stat2 += s ; }
  for ( const auto& s : _cov2 ) { cov2  += s ; }
  //
  if  ( 0 == stat1.nEntries() ) { return 0 ; }
  //
  cov2 /= stat1.nEntries() ;
  //
  const double v1_mean = stat1.mean() ;
  const double v2_mean = stat2.mean() ;
  //
  cov2 ( 0 , 0 ) -= v1_mean * v1_mean ;
  cov2 ( 0 , 1 ) -= v1_mean * v2_mean ;
  cov2 ( 1 , 1 ) -= v2_mean * v2_mean ;
  //
  return stat1.nEntries () ;  
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
( Ostap::DataFrame     frame  ,  
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
( Ostap::DataFrame     frame  ,  
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
( DataFrame            frame ,  
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
( DataFrame            frame ,  
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
( DataFrame            frame ,  
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
  std::vector<double> _quantiles_ 
  ( Ostap::DataFrame        frame , 
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
    values.clear() ;
    
    return  result ;
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
double Ostap::StatVar::quantile
( Ostap::DataFrame    frame  ,
  const double        q      , //  0<q<1 
  const std::string&  expr   , 
  const std::string&  cuts   ) 
{
  Ostap::Assert ( 0 < q && q < 1             , 
                  "Invalid quantile"         ,
                  "Ostap::StatVar::quantile" ) ;
  auto result = _quantiles_ ( frame , std::set<double>{{ q }} , expr , cuts ) ;
  Ostap::Assert ( 1 == result.size()         , 
                  "Invalid quantiles size"   ,
                  "Ostap::StatVar::interval" ) ;
  return result [0] ;
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
std::vector<double> Ostap::StatVar::quantiles
( Ostap::DataFrame           frame     ,
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
Ostap::StatVar::Interval 
Ostap::StatVar::interval 
( Ostap::DataFrame    frame ,
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
  const std::vector<double> result = 
    _quantiles_ ( frame , std::set<double>{{ q1 , q2 }} ,  expr , cuts ) ; 
  Ostap::Assert ( 2 == result.size()         ,
                  "Invalid interval"         ,
                  "Ostap::StatVar::interval" ) ;
  //
  return std::make_pair( result[0] , result[1] ) ;
}
// ============================================================================
// The END
// ============================================================================
