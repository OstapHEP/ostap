// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <vector>
#include <memory>
// ============================================================================
// ROOT 
// ============================================================================
#include "TTree.h"
#include "TCut.h"
#include "RooDataSet.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Formula.h"
#include "Ostap/MatrixUtils.h"
#include "Ostap/StatVar.h"
#include "Ostap/Iterator.h"
// ============================================================================
/** @file 
 *  Implementation file for class Analysis::StatVar
 *  @date 2013-10-13 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  static_assert (std::numeric_limits<unsigned long>::is_specialized   , 
                 "Numeric_limits<unsigned long> are not specialized!" ) ;
  // ==========================================================================
  /** @class Notifier
   *  Local helper class to keep the proper notifications for TTree
   *  @date 2013-10-13 
   *  @author Vanya BELYAEV Ivan.Brlyaev@itep.ru
   */
  class Notifier : public TObject 
  {
    // ======================================================================== 
 public:
    // ========================================================================
    Notifier 
    ( TTree*   tree     , 
      TObject* obj0     , 
      TObject* obj1 = 0 , 
      TObject* obj2 = 0 , 
      TObject* obj3 = 0 , 
      TObject* obj4 = 0 ,
      TObject* obj5 = 0 , 
      TObject* obj6 = 0 , 
      TObject* obj7 = 0 , 
      TObject* obj8 = 0 , 
      TObject* obj9 = 0 ) ;
    // templated constructor 
    template <class ITERATOR>
    Notifier ( ITERATOR  begin ,
               ITERATOR  end   , 
               TTree*    tree  ) ;
    /// virtual destructor 
    virtual        ~Notifier () ; // virtual destructor 
    /// the main method 
    virtual Bool_t  Notify   () ;
    // ========================================================================
  private:
    // ========================================================================
    Notifier () ;
    Notifier ( const Notifier & ) ;
    // ========================================================================
  private:
    // ========================================================================
    TTree*   m_tree ;
    TObject* m_old  ; // old notifier  
    // list of fobject to be notified 
    std::vector<TObject*> m_objects ;
    // ========================================================================
  } ;  
  // ==========================================================================
  Notifier::Notifier 
  ( TTree*   tree , 
    TObject* obj0 , 
    TObject* obj1 , 
    TObject* obj2 , 
    TObject* obj3 , 
    TObject* obj4 ,
    TObject* obj5 , 
    TObject* obj6 , 
    TObject* obj7 , 
    TObject* obj8 , 
    TObject* obj9 ) 
    : TObject   () 
    , m_tree    ( tree    ) 
    , m_old     ( nullptr )
    , m_objects () 
  {
    //
    if ( nullptr != m_tree ) { m_old = m_tree->GetNotify ()  ; }
    //
    if ( nullptr != m_old  ) { m_objects.push_back ( m_old ) ; } // prepend it! 
    //
    if ( nullptr != obj0   ) { m_objects.push_back ( obj0  ) ; }
    if ( nullptr != obj1   ) { m_objects.push_back ( obj1  ) ; }
    if ( nullptr != obj2   ) { m_objects.push_back ( obj2  ) ; }
    if ( nullptr != obj3   ) { m_objects.push_back ( obj3  ) ; }
    if ( nullptr != obj4   ) { m_objects.push_back ( obj4  ) ; }
    if ( nullptr != obj5   ) { m_objects.push_back ( obj5  ) ; }
    if ( nullptr != obj6   ) { m_objects.push_back ( obj6  ) ; }
    if ( nullptr != obj7   ) { m_objects.push_back ( obj7  ) ; }
    if ( nullptr != obj8   ) { m_objects.push_back ( obj8  ) ; }
    if ( nullptr != obj9   ) { m_objects.push_back ( obj9  ) ; }
    //
    if ( nullptr != m_tree ) { m_tree->SetNotify   ( this  ) ; }
    //
  }    
  // ==========================================================================
  template <class ITERATOR>
  Notifier::Notifier ( ITERATOR  begin ,
                       ITERATOR  end   , 
                       TTree*    tree  ) 
    : TObject   () 
    , m_tree    ( tree    ) 
    , m_old     ( nullptr )
    , m_objects () 
  {
    if ( nullptr != m_tree ) { m_old = m_tree->GetNotify ()  ; }
    if ( nullptr != m_old  ) { m_objects.push_back ( m_old ) ; } // prepend it! 
    for ( ; begin != end ; ++begin ) 
    {
      TObject* o = *begin ;
      if ( nullptr != o ) { m_objects.push_back ( o ) ; }
    }
    if ( nullptr != m_tree ) { m_tree->SetNotify   ( this  ) ; }
  }
  // ==========================================================================
  Notifier::~Notifier () { if ( nullptr != m_tree ) { m_tree->SetNotify ( m_old ) ; } }
  // ==========================================================================
  Bool_t  Notifier::Notify   () 
  {
    for ( TObject* o : m_objects ) { if ( nullptr != o ) { o->Notify() ; } }    
    return kTRUE ;
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
( TTree*             tree       , 
  const std::string& expression , 
  const unsigned long first     ,
  const unsigned long last      )
{
  Statistic result ;
  if ( 0 == tree || last <= first ) { return result ; }  // RETURN 
  Ostap::Formula formula ( "" , expression , tree ) ;
  if ( !formula.GetNdim() )         { return result ; }  // RETURN 
  //
  Notifier notify ( tree , &formula ) ;
  //
  const unsigned long nEntries = 
    std::min ( last , (unsigned long) tree->GetEntries() ) ;
  //
  for ( unsigned long entry = first ; entry < nEntries ; ++entry )   
  {
    //
    long ievent = tree->GetEntryNumber ( entry ) ;
    if ( 0 > ievent ) { return result ; }                // RETURN
    //
    ievent      = tree->LoadTree ( ievent ) ;
    if ( 0 > ievent ) { return result ; }                // RETURN 
    //
    result += formula.evaluate() ;
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
  Statistic result ;
  if ( 0 == tree || last <= first ) { return result ; }  // RETURN 
  Ostap::Formula selection ( "" , cuts      , tree ) ;
  if ( !selection.ok () ) { return result ; }            // RETURN 
  Ostap::Formula formula   ( "" , expression , tree ) ;
  if ( !formula  .ok () ) { return result ; }            // RETURN 
  //
  Notifier notify ( tree , &selection,  &formula ) ;
  //
  const unsigned long nEntries = 
    std::min ( last , (unsigned long) tree->GetEntries() ) ;
  //
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
    if   (  !w ) { continue ; }                          // ATTENTION
    //
    const double v = formula.evaluate() ;
    //
    result.add ( v , w ) ;
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
  //
  const std::string _cuts = cuts.GetTitle() ;
  //
  return statVar ( tree , expression , _cuts , first , last ) ;
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
  Ostap::Formula formula1 ( "" , exp1 , tree ) ;
  if ( !formula1 .ok () ) { return 0 ; }                   // RETURN 
  Ostap::Formula formula2 ( "" , exp2 , tree ) ;
  if ( !formula2 .ok () ) { return 0 ; }                   // RETURN 
  //
  Notifier notify ( tree , &formula1 , &formula2 ) ;
  //
  const unsigned long nEntries = 
    std::min ( last , (unsigned long) tree->GetEntries() ) ;
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
    const double v1 = formula1.evaluate() ;
    const double v2 = formula2.evaluate() ;
    //
    stat1 += v1 ;
    stat2 += v2 ;
    //
    cov2 ( 0 , 0 ) += v1*v1 ;
    cov2 ( 0 , 1 ) += v1*v2 ;
    cov2 ( 1 , 1 ) += v2*v2 ;
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
  //
  stat1.reset () ;
  stat2.reset () ;
  Ostap::Math::setToScalar ( cov2 , 0.0 ) ;
  //
  if ( 0 == tree || last <= first ) { return 0 ; }              // RETURN 
  Ostap::Formula formula1 ( "" , exp1 , tree ) ;
  if ( !formula1 .ok () ) { return 0 ; }                        // RETURN 
  Ostap::Formula formula2 ( "" , exp2 , tree ) ;
  if ( !formula2 .ok () ) { return 0 ; }                        // RETURN 
  Ostap::Formula selection ( "" , cuts      , tree ) ;
  if ( !selection.ok () ) { return 0 ; }                        // RETURN 
  //
  Notifier notify ( tree , &formula1 , &formula2 ) ;
  //
  const unsigned long nEntries = 
    std::min ( last , (unsigned long) tree->GetEntries() ) ;
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
    const double w = selection.evaluate() ;
    //
    if  ( !w ) { continue ; }                                  //  RETURN
    //
    const double v1 = formula1.evaluate() ;
    const double v2 = formula2.evaluate() ;
    //
    stat1.add ( v1 , w ) ;
    stat2.add ( v2 , w ) ;
    //
    cov2 ( 0 , 0 ) += w*v1*v1 ;
    cov2 ( 0 , 1 ) += w*v1*v2 ;
    cov2 ( 1 , 1 ) += w*v2*v2 ;
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
Ostap::StatVar::Statistic
Ostap::StatVar::statVar 
( const RooAbsData*   data       , 
  const std::string&  expression , 
  const unsigned long first      ,
  const unsigned long last       )
{
  Statistic result ;
  if ( 0 == data || last <= first ) { return result ; }             // RETURN 
  //
  RooArgList        alst ;
  const RooArgSet*  aset = data->get() ;
  if ( 0 == aset    )               { return result ; }            // RETURN 
  Ostap::Utils::Iterator iter ( *aset );
  //
  RooAbsArg*   coef = 0 ;
  while ( ( coef = (RooAbsArg*) iter.next() ) ) 
  { alst.add ( *coef ); }
  //
  RooFormulaVar formula ( "" ,  expression.c_str() , alst ) ;
  if ( !formula.ok()   ) { return result ; }                        // RETURN
  //
  const bool weighted = data->isWeighted() ;
  //
  const unsigned long nEntries = 
    std::min ( last , (unsigned long) data->numEntries() ) ;
  // start the loop 
  for ( unsigned long entry = first ; nEntries > entry ; ++entry ) 
  {
    //
    if ( 0 == data->get( entry)  ) { return result ; }             // RETURN 
    //
    if ( !weighted ) { result += formula.getVal() ; }
    else 
    {
      const double w = data->weight() ;
      if ( w ) { result.add ( formula.getVal () , w ) ; }    //  ATTENTION
    } 
  }
  //
  return result ;
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
  return statVar ( data , expression , _cuts , first , last ) ;
}
// ============================================================================
Ostap::StatVar::Statistic
Ostap::StatVar::statVar 
( const RooAbsData*   data        , 
  const std::string&  expression  , 
  const std::string&  cuts        , 
  const unsigned long first       ,
  const unsigned long last        ) 
{
  Statistic result ;
  if ( 0 == data || last <= first ) { return result ; }         // RETURN 
  //
  RooArgList        alst ;
  const RooArgSet*  aset = data->get() ;
  if ( 0 == aset    ) { return result ; }                       // RETURN
  Ostap::Utils::Iterator iter ( *aset );
  //
  RooAbsArg*   coef = 0 ;
  while ( ( coef = (RooAbsArg*) iter.next() ) ) 
  { alst.add ( *coef ); }
  //
  RooFormulaVar formula   ( "" ,  expression.c_str() , alst ) ;
  if ( !formula.ok()   ) { return result ; }                     // RETURN
  //
  RooFormulaVar selection ( "" ,  cuts      .c_str() , alst ) ;
  if ( !selection.ok() ) { return result ; }                     // RETURN 
  //
  const bool weighted = data->isWeighted() ;
  //
  const unsigned long nEntries = 
    std::min ( last , (unsigned long) data->numEntries() ) ;
  // start the loop 
  for ( unsigned long entry = first ; nEntries > entry ; ++entry ) 
  {
    //
    if ( 0 == data->get( entry)  ) { return result ; }            // RETURN
    //
    const double w = 
      weighted ? 
      selection.getVal () * data->weight() :
      selection.getVal ()                  ;
    //
    if  ( !w ) { continue ; }                                  // ATTENTION
    //
    const double v = formula.getVal () ;
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
( const RooAbsData*    data    , 
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
  if ( 0 == data || last <= first ) { return 0 ; }               // RETURN 
  //
  RooArgList        alst ;
  const RooArgSet*  aset = data->get() ;
  if ( 0 == aset       ) { return  0 ; }                          // RETURN
  Ostap::Utils::Iterator iter ( *aset );
  //
  RooAbsArg*   coef = 0 ;
  while ( ( coef = (RooAbsArg*) iter.next() ) ) 
  { alst.add ( *coef ); }
  //
  RooFormulaVar formula1 ( "" ,  exp1.c_str() , alst ) ;
  if ( !formula1.ok()   ) { return 0 ; }                          // RETURN
  RooFormulaVar formula2 ( "" ,  exp2.c_str() , alst ) ;
  if ( !formula2.ok()   ) { return 0 ; }                          // RETURN
  //
  const bool weighted = data->isWeighted() ;
  //
  const unsigned long nEntries = 
    std::min ( last , (unsigned long) data->numEntries() ) ;
  //
  for ( unsigned long entry = first ; entry < nEntries ; ++entry )   
  {
    //
    if ( 0 == data->get( entry)  ) { break ; }                    // BREAK
    //    
    const double w  = weighted ? data->weight() : 1.0 ;
    //
    if ( !w ) { continue ; }                                     //  ATTENTION
    //
    const double v1 = formula1.getVal() ;
    const double v2 = formula2.getVal() ;
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
( const RooAbsData*    data    , 
  const std::string&   exp1    , 
  const std::string&   exp2    , 
  const std::string&   cuts    , 
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
  if ( 0 == data || last <= first ) { return 0 ; }               // RETURN 
  //
  RooArgList        alst ;
  const RooArgSet*  aset = data->get() ;
  if ( 0 == aset       ) { return  0 ; }                          // RETURN
  Ostap::Utils::Iterator iter ( *aset );
  //
  RooAbsArg*   coef = 0 ;
  while ( ( coef = (RooAbsArg*) iter.next() ) ) 
  { alst.add ( *coef ); }
  //
  RooFormulaVar formula1  ( "" ,  exp1.c_str() , alst ) ;
  if ( !formula1.ok()   ) { return 0 ; }                          // RETURN
  RooFormulaVar formula2  ( "" ,  exp2.c_str() , alst ) ;
  if ( !formula2.ok()   ) { return 0 ; }                          // RETURN
  RooFormulaVar selection ( "" ,  cuts.c_str() , alst ) ;
  if ( !selection.ok()  ) { return 0 ; }                          // RETURN
  //
  const bool weighted = data->isWeighted() ;
  //
  const unsigned long nEntries = 
    std::min ( last , (unsigned long) data->numEntries() ) ;
  //
  for ( unsigned long entry = first ; entry < nEntries ; ++entry )   
  {
    //
    if ( 0 == data->get( entry)  ) { break ; }                    // BREAK
    //    
    const double w = 
      weighted ? 
      selection.getVal () * data->weight() :
      selection.getVal ()                  ;
    //
    if ( !w ) { continue ; }                                    // ATTENTION
    //
    const double v1 = formula1.getVal() ;
    const double v2 = formula2.getVal() ;
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
/*  calculate the covariances for generic case 
 *  @param tree  (INPUT)  the input T-tree 
 *  @param vars  (INPUT)  the list of variables 
 *  @param cuts  (INPUT)  the selection criteria/weights 
 *  @param stats (UPDATE) list of statistics 
 *  @param covs  (UPDATE) elements of "covariance matrix" 
 *  @return number of processed events 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2017-02-17
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  inline unsigned short index 
  ( const unsigned short i , 
    const unsigned short j ) { return i * ( i + 1 ) / 2 + j ; } 
  // ==========================================================================
}  
// ============================================================================
unsigned long Ostap::StatVar::_statCov
( TTree*                                  tree  ,  
  const std::vector<std::string>&         vars  ,
  const std::string&                      cuts  ,
  std::vector<Ostap::StatVar::Statistic>& stats , 
  std::vector<double>&                    covs  , 
  const unsigned long                     first ,
  const unsigned long                     last  )
{
  //
  //
  if ( 0 == tree || last <= first ) { return 0 ; }              // RETURN 
  //
  std::vector<std::unique_ptr<Ostap::Formula> > _vars ;
  for ( const std::string& v : vars ) 
  {
    // auto f = std::make_unique<Ostap::Formula>( "" , v , tree ) ;
    auto f = std::unique_ptr<Ostap::Formula>( new Ostap::Formula( "" , v , tree ) ) ;
    if ( !f || !f->ok() ) { return 0 ; }                        // RETURN
    _vars.push_back ( std::move ( f ) ) ;
  }
  //
  const unsigned short l = _vars.size() ;
  //
  stats.resize( l         ) ;
  covs.resize ( l*(l+1)/2 ) ;
  if ( 0 == l          ) { return 0 ; }  // RETURN
  //
  std::vector<double> vals  ( l , 0.0 ) ;
  //
  Ostap::Formula selection ( "" , cuts      , tree ) ;
  if ( !selection.ok () ) { return 0 ; }                        // RETURN
  //
  std::vector<TObject*> tobjs ;
  for ( auto& _v : _vars ) { tobjs.push_back ( _v.get() ) ; }
  tobjs.push_back( &selection ) ;
  //
  Notifier notifier ( tobjs.begin() , tobjs.end() , tree ) ;
  //
  const unsigned long nEntries = 
    std::min ( last , (unsigned long) tree->GetEntries() ) ;
  //
  for ( unsigned long entry = first ; entry < nEntries ; ++entry )   
  {
    long ievent = tree->GetEntryNumber ( entry ) ;
    if ( 0 > ievent ) { break ; }                              // RETURN
    //
    ievent      = tree->LoadTree ( ievent ) ;
    if ( 0 > ievent ) { break ; }                              // RETURN
    
    //
    const double w = selection.evaluate() ;
    //
    if  ( !w ) { continue ; }                                  // ATTENTION!
    //
    // fill satistics and covarinaces:
    for ( unsigned short i = 0 ; i < l ; ++i ) 
    {
      const double vi = _vars[i]->evaluate() ;
      vals [i] = vi ;
      stats[i].add ( vi , w ) ;  
      for ( unsigned short j = 0 ; j <= i ; ++j ) 
      {
        const double vj = vals[j] ;
        covs [ index ( i , j ) ] += w*vi*vj ;
      } 
    } 
  } //                                    the end of loop over entries in Ttree
  //
  if ( stats.empty() || 0 == stats[0].nEntries() || 0 == stats[0].nEff () ) { return 0 ; }
  const double sumw = stats[0].weights().sum() ; 
  //
  for ( unsigned short i = 0 ; i < l ; ++i ) 
  {
    const double vi_mean = stats[i].mean() ;
    for ( unsigned short j = 0 ; j <= i ; ++j ) 
    {
      const double          vj_mean = stats[j].mean() ;
      const unsigned short  ij      = index ( i , j ) ; 
      covs [ ij ] /= sumw ;
      covs [ ij ] -= vi_mean * vj_mean ;
    }
  }
  //
  return stats[0].nEntries() ;
}
// ============================================================================
/** calculate the covariances for generic case 
 *  @param tree  (INPUT)  the input T-tree 
 *  @param vars  (INPUT)  the list of variables 
 *  @param stats (UPDATE) list of statistics 
 *  @param covs  (UPDATE) elements of "covariance matrix" 
 *  @return number of processed events 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2017-02-17
 */
// ============================================================================
unsigned long Ostap::StatVar::_statCov
( TTree*                                  tree  ,  
  const std::vector<std::string>&         vars  ,
  std::vector<Ostap::StatVar::Statistic>& stats , 
  std::vector<double>&                    covs  , 
  const unsigned long                     first ,
  const unsigned long                     last  )
{
  //
  if ( 0 == tree || last <= first ) { return 0 ; }              // RETURN 
  //
  std::vector<std::unique_ptr<Ostap::Formula> > _vars ;
  for ( const std::string& v : vars ) 
  {
    // auto f = std::make_unique<Ostap::Formula>( "" , v , tree ) ;
    auto f = std::unique_ptr<Ostap::Formula> ( new Ostap::Formula( "" , v , tree ) ) ;
    if ( !f || !f->ok() ) { return 0 ; }                        // RETURN
    _vars.push_back ( std::move ( f ) ) ;
  }
  //
  const unsigned short l = _vars.size() ;
  //
  stats.resize( l         ) ;
  covs.resize ( l*(l+1)/2 ) ;
  if ( 0 == l          ) { return 0 ; }  // RETURN
  //
  std::vector<double> vals ( l , 0.0 ) ;
  //
  std::vector<TObject*> tobjs ;
  for ( auto& _v : _vars ) { tobjs.push_back ( _v.get() ) ; }
  //
  Notifier notifier ( tobjs.begin() , tobjs.end() , tree ) ;
  //
  const unsigned long nEntries = 
    std::min ( last , (unsigned long) tree->GetEntries() ) ;
  //
  for ( unsigned long entry = first ; entry < nEntries ; ++entry )   
  {
    long ievent = tree->GetEntryNumber ( entry ) ;
    if ( 0 > ievent ) { break ; }                              // RETURN
    //
    ievent      = tree->LoadTree ( ievent ) ;
    if ( 0 > ievent ) { break ; }                              // RETURN
    //
    const double w = 1 ;
    // fill satistics and covarinaces:
    for ( unsigned short i = 0 ; i < l ; ++i ) 
    {
      const double vi = _vars[i]->evaluate() ;
      vals [i] = vi ;
      stats[i].add ( vi , w ) ;  
      for ( unsigned short j = 0 ; j <= i ; ++j ) 
      {
        const double vj = vals[j] ;
        covs [ index ( i , j ) ] +=  w * vi * vj ;
      } 
    } 
  } //                                    the end of loop over entries in Ttree
  //
  if ( stats.empty() || 0 == stats[0].nEntries() || 0 == stats[0].nEff () ) { return 0 ; }
  const double sumw = stats[0].weights().sum() ; 
  //
  for ( unsigned short i = 0 ; i < l ; ++i ) 
  {
    const double vi_mean = stats[i].mean() ;
    for ( unsigned short j = 0 ; j <= i ; ++j ) 
    {
      const double          vj_mean = stats[j].mean() ;
      const unsigned short  ij      = index ( i , j ) ; 
      covs [ ij ] /= sumw ;
      covs [ ij ] -= vi_mean * vj_mean ;
    }
  }
  //
  return stats[0].nEntries() ;
}
// ============================================================================
/* calculate the covariances for generic case 
 *  @param tree  (INPUT)  the input tree 
 *  @param vars  (INPUT)  the list of variables 
 *  @param cuts  (INPUT)  the selection criteria/weights 
 *  @param stats (UPDATE) list of statistics 
 *  @param covs  (UPDATE) elements of "covariance matrix" 
 *  @return number of processed events 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2017-02-17
 */
// ============================================================================
unsigned long Ostap::StatVar::_statCov
( TTree*               tree   , 
  const std::vector<std::string>&         vars  , 
  const TCut&                             cuts  ,
  std::vector<Ostap::StatVar::Statistic>& stats ,  
  std::vector<double>&                    covs  ,
  const unsigned long                     first ,
  const unsigned long                     last  ) 
{
  const std::string _cuts = cuts.GetTitle() ;
  return _statCov ( tree , vars , _cuts , stats , covs , first , last ) ;
}
// ========================================================================
/*  calculate the covariances for generic case 
 *  @param data  (INPUT)  the inpout dataset 
 *  @param vars  (INPUT)  the list of variables 
 *  @param cuts  (INPUT)  the selection criteria/weights 
 *  @param stats (UPDATE) list of statistics 
 *  @param covs  (UPDATE) elements of "covariance matrix" 
 *  @return number of processed events 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2017-02-17
 */
// ========================================================================    
unsigned long Ostap::StatVar::_statCov 
( const RooAbsData*               data  , 
  const std::vector<std::string>& vars  , 
  const std::string&              cuts  ,
  std::vector<Statistic>&         stats ,  
  std::vector<double>&            covs  ,
  const unsigned long             first ,
  const unsigned long             last  ) 
{
  //
  if ( 0 == data || last <= first ) { return 0 ; }               // RETURN 
  //
  RooArgList        alst ;
  const RooArgSet*  aset = data->get() ;
  if ( 0 == aset       ) { return  0 ; }                          // RETURN
  Ostap::Utils::Iterator iter ( *aset );
  //
  RooAbsArg*   coef = 0 ;
  while ( ( coef = (RooAbsArg*) iter.next() ) ) 
  { alst.add ( *coef ); }
  // variables?
  std::vector<std::unique_ptr<RooFormulaVar> > _vars ;
  for ( const std::string& v : vars ) 
  {
    // auto f = std::make_unique<RooFormulaVar>( "" , v.c_str() , alst ) ;
    auto f = std::unique_ptr<RooFormulaVar> ( new RooFormulaVar( "" , v.c_str() , alst ) ) ;
    if ( !f || !f->ok() ) { return 0 ; }                        // RETURN
    _vars.push_back ( std::move ( f ) ) ;
  }
  // cuts ? 
  RooFormulaVar selection ( "" ,  cuts.c_str() , alst ) ;
  if ( !selection.ok()  ) { return 0 ; }                          // RETURN
  //
  const unsigned short l = _vars.size() ;
  //
  stats.resize( l         ) ;
  covs.resize ( l*(l+1)/2 ) ;
  if ( 0 == l          ) { return 0 ; }  // RETURN
  //
  std::vector<double> vals ( l , 0.0 ) ;
  //
  const bool weighted = data->isWeighted() ;
  //
  const unsigned long nEntries = 
    std::min ( last , (unsigned long) data->numEntries() ) ;
  //
  for ( unsigned long entry = first ; entry < nEntries ; ++entry )   
  {
    //
    if ( 0 == data->get( entry)  ) { break ; }                    // BREAK
    //    
    const double w = 
      weighted ? 
      selection.getVal () * data->weight() :
      selection.getVal ()                  ;
    //
    if   ( !w ) {  continue ; }                                   // ATTENTION!
    //
    // fill satisticst and covarinaces:
    for ( unsigned short i = 0 ; i < l ; ++i ) 
    {
      const double vi = _vars[i]->getVal() ;
      vals [i] = vi ;
      stats[i].add ( vi , w ) ;  
      for ( unsigned short j = 0 ; j <= i ; ++j ) 
      {
        const double vj = vals[j] ;
        covs [ index ( i , j ) ] += w*vi*vj ;
      } 
    }
  }
  //
  if ( stats.empty() || 0 == stats[0].nEntries() || 0 == stats[0].nEff () ) { return 0 ; }
  const double sumw = stats[0].weights().sum() ; 
  //
  for ( unsigned short i = 0 ; i < l ; ++i ) 
  {
    const double vi_mean = stats[i].mean() ;
    for ( unsigned short j = 0 ; j <= i ; ++j ) 
    {
      const double          vj_mean = stats[j].mean() ;
      const unsigned short  ij      = index ( i , j ) ; 
      covs [ ij ] /= sumw ;
      covs [ ij ] -= vi_mean * vj_mean ;
    }
  }
  //
  return stats[0].nEntries() ;
}
// ========================================================================
/*  calculate the covariances for generic case 
 *  @param data  (INPUT)  the inpout dataset 
 *  @param vars  (INPUT)  the list of variables 
 *  @param stats (UPDATE) list of statistics 
 *  @param covs  (UPDATE) elements of "covariance matrix" 
 *  @return number of processed events 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2017-02-17
 */
// ========================================================================    
unsigned long Ostap::StatVar::_statCov 
( const RooAbsData*               data  , 
  const std::vector<std::string>& vars  , 
  std::vector<Statistic>&         stats ,  
  std::vector<double>&            covs  ,
  const unsigned long             first ,
  const unsigned long             last  ) 
{
  //
  if ( 0 == data || last <= first ) { return 0 ; }               // RETURN 
  //
  RooArgList        alst ;
  const RooArgSet*  aset = data->get() ;
  if ( 0 == aset       ) { return  0 ; }                          // RETURN
  Ostap::Utils::Iterator iter ( *aset );
  //
  RooAbsArg*   coef = 0 ;
  while ( ( coef = (RooAbsArg*) iter.next() ) ) 
  { alst.add ( *coef ); }
  // variables?
  std::vector<std::unique_ptr<RooFormulaVar> > _vars ;
  for ( const std::string& v : vars ) 
  {
    // auto f = std::make_unique<RooFormulaVar>( "" , v.c_str() , alst ) ;
    auto f = std::unique_ptr<RooFormulaVar> ( new RooFormulaVar( "" , v.c_str() , alst ) ) ;
    if ( !f || !f->ok() ) { return 0 ; }                        // RETURN
    _vars.push_back ( std::move ( f ) ) ;
  }
  //
  const unsigned short l = _vars.size() ;
  //
  stats.resize( l         ) ;
  covs.resize ( l*(l+1)/2 ) ;
  if ( 0 == l          ) { return 0 ; }  // RETURN
  //
  std::vector<double> vals ( l , 0.0 ) ;
  //
  const bool weighted = data->isWeighted() ;
  //
  const unsigned long nEntries = 
    std::min ( last , (unsigned long) data->numEntries() ) ;
  //
  for ( unsigned long entry = first ; entry < nEntries ; ++entry )   
  {
    //
    if ( 0 == data->get( entry)  ) { break ; }                    // BREAK
    //    
    const double w = weighted ? data->weight() : 1.0 ;
    //
    if  ( !w ) { continue ; }                                     // ATTENTION
    //
    // fill satisticst and covarinaces:
    for ( unsigned short i = 0 ; i < l ; ++i ) 
    {
      const double vi = _vars[i]->getVal() ;
      vals [i] = vi ;
      stats[i].add ( vi , w ) ;  
      for ( unsigned short j = 0 ; j <= i ; ++j ) 
      {
        const double vj = vals[j] ;
        covs [ index ( i , j ) ] += w*vi*vj ;
      } 
    }
  }
  //
  if ( stats.empty() || 0 == stats[0].nEntries() || 0 == stats[0].nEff () ) { return 0 ; }
  const double sumw = stats[0].weights().sum() ; 
  //
  for ( unsigned short i = 0 ; i < l ; ++i ) 
  {
    const double vi_mean = stats[i].mean() ;
    for ( unsigned short j = 0 ; j <= i ; ++j ) 
    {
      const double          vj_mean = stats[j].mean() ;
      const unsigned short  ij      = index ( i , j ) ; 
      covs [ ij ] /= sumw ;
      covs [ ij ] -= vi_mean * vj_mean ;
    }
  }
  //
  return stats[0].nEntries() ;
}
// ========================================================================    
 

// ============================================================================
// The END 
// ============================================================================
