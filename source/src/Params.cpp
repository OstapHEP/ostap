// ============================================================================
// Include files 
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Params.h"
#include "Ostap/Formula.h"
#include "Ostap/Notifier.h"
#include "Ostap/Parameterization.h"
#include "Ostap/Polynomials.h"
// ============================================================================
// ROOT
// ============================================================================
#include "TTree.h"
// ============================================================================
// Local 
// ============================================================================
#include "Exception.h"
// ============================================================================
/** @file 
 *  Implementation file for header Ostap/Params.h
 *  @date 2019-07-03 
 *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
 */
// ============================================================================
/*  fill Legendre sum with data from the Tree 
 *  @see Ostap::Math::LegendreSum 
 *  @see Ostap::Math::LegendreSum::fill
 *  @param tree       (INPUT)  the input tree 
 *  @param sum        (UPDATE) the parameterization object 
 *  @param expression (INPUT)  expresion to be parameterized
 *  @param first      (INPUT)  the first event in Tree 
 *  @param last       (INPUT)  the last  event in Tree 
 *  @return number of events sused in parameterization 
 *  @code
 *  Tree*  tree = ...
 *  LegendreSum s ( 5 , -1 , 1 ) ;
 *  DataParam::parameterize ( tree , s , "x" ) ;
 *  @endcode
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2019-07-3
 */ 
// ============================================================================
unsigned long Ostap::DataParam::parameterize 
( TTree*                    tree       , 
  Ostap::Math::LegendreSum& sum        , 
  const std::string&        expression , 
  const unsigned long       first      ,
  const unsigned long       last       ) 
{
  // invalid tree or emptry range 
  if ( nullptr == tree || last <= first ) { return 0 ; }
  //
  Ostap::Formula var ( expression , tree ) ;
  Ostap::Assert ( var.ok()                                    , 
                  "Invalid expression:\"" + expression + "\"" ,
                  "Ostap::DataParams::parameterize"           ) ;
  //
  Ostap::Utils::Notifier notify ( tree , &var ) ;
  //
  const unsigned long nEntries = std::min ( last , (unsigned long) tree->GetEntries() ) ;
  unsigned long filled = 0 ;
  for ( unsigned long entry = first ; entry < nEntries ; ++entry ) 
  { 
    long ievent = tree->GetEntryNumber ( entry ) ;
    if ( 0 > ievent ) { break ; }                        // BREAK
    //
    ievent      = tree->LoadTree ( ievent ) ;
    if ( 0 > ievent ) { break ; }                        // BREAK
    //
    // get the variable   
    const double v = var.evaluate() ;
    //
    // fill the sum 
    filled += ( sum.fill ( v , 1 ) ? 1 : 0 ) ;
  }
  //
  return filled ;
}
// ==================================================================================
/*  fill Legendre sum with data from the Tree 
 *  @see Ostap::Math::LegendreSum 
 *  @see Ostap::Math::LegendreSum::fill
 *  @param tree       (INPUT)  the input tree 
 *  @param sum        (UPDATE) the parameterization object 
 *  @param expression (INPUT)  expression to be parameterized
 *  @param seelction  (INPUT)  selection/weight to be used 
 *  @param first      (INPUT)  the first event in Tree 
 *  @param last       (INPUT)  the last  event in Tree 
 *  @return  sum of weigths  used in parameterization
 *  @code
 *  Tree*  tree = ...
 *  LegendreSum s ( 5 , -1 , 1 ) ;
 *  DataParam::parameterize ( tree , s , "x"  , "y>10" ) ;
 *  @endcode
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2019-07-3
 */ 
// ==================================================================================
double Ostap::DataParam::parameterize 
( TTree*                    tree       , 
  Ostap::Math::LegendreSum& sum        , 
  const std::string&        expression , 
  const std::string&        selection  , 
  const unsigned long       first      ,
  const unsigned long       last       ) 
{
  // invalid tree or emptry range 
  if ( nullptr == tree || last <= first ) { return 0 ; }
  //
  if ( 0 == selection.size() ) 
  { return parameterize ( tree ,  sum , expression , first , last ) ; }
  //
  //
  Ostap::Formula var ( expression , tree ) ;
  Ostap::Assert ( var.ok()                                    , 
                  "Invalid expression:\"" + expression + "\"" ,
                  "Ostap::DataParams::parameterize"           ) ;
  //
  Ostap::Formula weight ( selection , tree ) ;
  Ostap::Assert ( weight.ok()                                 , 
                  "Invalid selection:\"" + selection   + "\"" ,
                  "Ostap::DataParams::parameterize"           ) ;
  //
  Ostap::Utils::Notifier notify ( tree , &var ,  &weight ) ;
  //
  const unsigned long nEntries = std::min ( last , (unsigned long) tree->GetEntries() ) ;
  //
  long double wsum = 0 ;
  for ( unsigned long entry = first ; entry < nEntries ; ++entry ) 
  { 
    long ievent = tree->GetEntryNumber ( entry ) ;
    if ( 0 > ievent ) { break ; }                        // BREAK
    //
    ievent      = tree->LoadTree ( ievent ) ;
    if ( 0 > ievent ) { break ; }                        // BREAK
    //
    // get the weight  
    const double w = weight.evaluate () ;
    if ( !w        ) { continue ; }
    //
    // get the variable   
    const double v =    var.evaluate() ;
    //
    // fill the sum 
    wsum += ( sum.fill ( v , w ) ? w : 0 ) ;
  }
  //
  return wsum ;
}
// ============================================================================
/*  fill Legendre sum with data from the Tree 
 *  @see Ostap::Math::LegendreSum2 
 *  @see Ostap::Math::LegendreSum2::fill
 *  @param tree        (INPUT)  the input tree 
 *  @param sum         (UPDATE) the parameterization object 
 *  @param xexpression (INPUT)  x-expression to be parameterized
 *  @param yexpression (INPUT)  y-expression to be parameterized
 *  @param first       (INPUT)  the first event in Tree 
 *  @param last        (INPUT)  the last  event in Tree 
 *  @return number of events used in parameterization 
 *  @code
 *  Tree*  tree = ...
 *  LegendreSum2 s ( 5 , 3 , -1 , 1 , -2 , 2 ) ;
 *  DataParam::parameterize ( tree , s , "x" , "y/z") ;
 *  @endcode
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2019-07-3
 */ 
// ============================================================================
unsigned long Ostap::DataParam::parameterize 
( TTree*                     tree        , 
  Ostap::Math::LegendreSum2& sum         , 
  const std::string&         xexpression , 
  const std::string&         yexpression , 
  const unsigned long        first       ,
  const unsigned long        last        ) 
{
  // invalid tree or emptry range 
  if ( nullptr == tree || last <= first ) { return 0 ; }
  //
  Ostap::Formula xvar ( xexpression , tree ) ;
  Ostap::Assert ( xvar.ok()                                     , 
                  "Invalid x-expression:\"" + xexpression + "\"" ,
                  "Ostap::DataParams::parameterize  "           ) ;
  Ostap::Formula yvar ( yexpression , tree ) ;
  Ostap::Assert ( yvar.ok()                                     , 
                  "Invalid y-expression:\"" + yexpression + "\"" ,
                  "Ostap::DataParams::parameterize  "           ) ;
  //
  Ostap::Utils::Notifier notify ( tree , &xvar, &yvar ) ;
  //
  const unsigned long nEntries = std::min ( last , (unsigned long) tree->GetEntries() ) ;
  unsigned long filled = 0 ;
  for ( unsigned long entry = first ; entry < nEntries ; ++entry ) 
  { 
    long ievent = tree->GetEntryNumber ( entry ) ;
    if ( 0 > ievent ) { break ; }                        // BREAK
    //
    ievent      = tree->LoadTree ( ievent ) ;
    if ( 0 > ievent ) { break ; }                        // BREAK
    //
    // get the variable   
    const double x = xvar.evaluate() ;
    const double y = yvar.evaluate() ;
    //
    // fill the sum 
    filled += ( sum.fill ( x , y  , 1 ) ? 1 : 0 ) ;
  }
  //
  return filled ;
}
// ============================================================================
/*  fill Legendre sum with data from the Tree 
 *  @see Ostap::Math::LegendreSum2 
 *  @see Ostap::Math::LegendreSum2::fill
 *  @param tree        (INPUT)  the input tree 
 *  @param sum         (UPDATE) the parameterization object 
 *  @param xexpression (INPUT)  expression to be parameterized
 *  @param yexpression (INPUT)  expression to be parameterized
 *  @param selection   (INPUT)  selection/weight to be used 
 *  @param first       (INPUT)  the first event in Tree 
 *  @param last        (INPUT)  the last  event in Tree 
 *  @return  sum of weigths  used in parameterization
 *  @code
 *  Tree*  tree = ...
 *  LegendreSum2 s ( 5 , 2 ,  -1 , 1 , -4  , -5 ) ;
 *  DataParam::parameterize ( tree , s , "x"  , "z" , "y>10" ) ;
 *  @endcode
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2019-07-3
 */ 
// ============================================================================
double Ostap::DataParam::parameterize 
( TTree*                     tree        ,  
  Ostap::Math::LegendreSum2& sum         , 
  const std::string&         xexpression , 
  const std::string&         yexpression , 
  const std::string&         selection   , 
  const unsigned long        first       ,
  const unsigned long        last        ) 
{
  // invalid tree or emptry range 
  if ( nullptr == tree || last <= first ) { return 0 ; }
  //
  if ( 0 == selection.size() ) 
  { return parameterize ( tree , sum , xexpression , yexpression , first , last ) ; }
  //
  Ostap::Formula xvar ( xexpression , tree ) ;
  Ostap::Assert ( xvar.ok()                                     , 
                  "Invalid x-expression:\"" + xexpression + "\"" ,
                  "Ostap::DataParams::parameterize  "           ) ;
  Ostap::Formula yvar ( yexpression , tree ) ;
  Ostap::Assert ( yvar.ok()                                     , 
                  "Invalid y-expression:\"" + yexpression + "\"" ,
                  "Ostap::DataParams::parameterize  "           ) ;

  Ostap::Formula weight  ( selection , tree ) ;
  Ostap::Assert ( weight.ok()                                 , 
                  "Invalid selection:\"" + selection   + "\"" ,
                  "Ostap::DataParams::parameterize"           ) ;
  //
  Ostap::Utils::Notifier notify ( tree , &xvar, &yvar , &weight ) ;
  //
  const unsigned long nEntries = std::min ( last , (unsigned long) tree->GetEntries() ) ;
  //
  long double sumw = 0 ;
  for ( unsigned long entry = first ; entry < nEntries ; ++entry ) 
  { 
    long ievent = tree->GetEntryNumber ( entry ) ;
    if ( 0 > ievent ) { break ; }                        // BREAK
    //
    ievent      = tree->LoadTree ( ievent ) ;
    if ( 0 > ievent ) { break ; }                        // BREAK
    //
    // get the weight  
    const double w = weight.evaluate () ;
    if ( !w        ) { continue ; }                      // CONTINUE 
    //
    // get the variable   
    const double x = xvar.evaluate() ;
    const double y = yvar.evaluate() ;
    //
    // fill the sum 
    sumw += ( sum.fill ( x , y  , w ) ? w : 0 ) ;
  }
  //
  return sumw ;
}
// ============================================================================
/*  fill Legendre sum with data from the Tree 
 *  @see Ostap::Math::LegendreSum3 
 *  @see Ostap::Math::LegendreSum3::fill
 *  @param tree        (INPUT)  the input tree 
 *  @param sum         (UPDATE) the parameterization object 
 *  @param xexpression (INPUT)  x-expression to be parameterized
 *  @param yexpression (INPUT)  y-expression to be parameterized
 *  @param zexpression (INPUT)  z-expression to be parameterized
 *  @param first       (INPUT)  the first event in Tree 
 *  @param last        (INPUT)  the last  event in Tree 
 *  @return number of events used in parameterization 
 *  @code
 *  Tree*  tree = ...
 *  LegendreSum3 s ( 5 , 3 , 2 , -1 , 1 , -2 , 2 , 0 , 4 ) ;
 *  DataParam::parameterize ( tree , s , "x" , "y" , "y/z") ;
 *  @endcode
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2019-07-3
 */ 
// ============================================================================
unsigned long Ostap::DataParam::parameterize 
( TTree*                     tree        , 
  Ostap::Math::LegendreSum3& sum         , 
  const std::string&         xexpression , 
  const std::string&         yexpression , 
  const std::string&         zexpression , 
  const unsigned long        first       ,
  const unsigned long        last        ) 
{
  // invalid tree or emptry range 
  if ( nullptr == tree || last <= first ) { return 0 ; }
  //
  Ostap::Formula xvar ( xexpression , tree ) ;
  Ostap::Assert ( xvar.ok()                                      , 
                  "Invalid x-expression:\"" + xexpression + "\"" ,
                  "Ostap::DataParams::parameterize  "            ) ;
  Ostap::Formula yvar ( yexpression , tree ) ;
  Ostap::Assert ( yvar.ok()                                      , 
                  "Invalid y-expression:\"" + yexpression + "\"" ,
                  "Ostap::DataParams::parameterize  "            ) ;
  Ostap::Formula zvar ( zexpression , tree ) ;
  Ostap::Assert ( zvar.ok()                                      , 
                  "Invalid z-expression:\"" + zexpression + "\"" ,
                  "Ostap::DataParams::parameterize  "            ) ;
  //
  Ostap::Utils::Notifier notify ( tree , &xvar, &yvar , &zvar ) ;
  //
  const unsigned long nEntries = std::min ( last , (unsigned long) tree->GetEntries() ) ;
  unsigned long filled = 0 ;
  for ( unsigned long entry = first ; entry < nEntries ; ++entry ) 
  { 
    long ievent = tree->GetEntryNumber ( entry ) ;
    if ( 0 > ievent ) { break ; }                        // BREAK
    //
    ievent      = tree->LoadTree ( ievent ) ;
    if ( 0 > ievent ) { break ; }                        // BREAK
    //
    // get the variable   
    const double x = xvar.evaluate() ;
    const double y = yvar.evaluate() ;
    const double z = zvar.evaluate() ;
    //
    // fill the sum 
    filled += ( sum.fill ( x , y  , z , 1 ) ? 1 : 0 ) ;
  }
  //
  return filled ;
}
// ============================================================================
/*  fill Legendre sum with data from the Tree 
 *  @see Ostap::Math::LegendreSum3 
 *  @see Ostap::Math::LegendreSum3::fill
 *  @param tree        (INPUT)  the input tree 
 *  @param sum         (UPDATE) the parameterization object 
 *  @param xexpression (INPUT)  x-expression to be parameterized
 *  @param yexpression (INPUT)  y-expression to be parameterized
 *  @param zexpression (INPUT)  z-expression to be parameterized
 *  @param selection   (INPUT)  selection/weight to be used 
 *  @param first       (INPUT)  the first event in Tree 
 *  @param last        (INPUT)  the last  event in Tree 
 *  @return  sum of weigths  used in parameterization
 *  @code
 *  Tree*  tree = ...
 *  LegendreSum3 s ( 5 , 2 , 4 ,  -1 , 1 , -4  , -5 , 0 , 3 ) ;
 *  DataParam::parameterize ( tree , s , "x"  , "y" , "z" , "t>10" ) ;
 *  @endcode
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2019-07-3
 */ 
// ============================================================================
double Ostap::DataParam::parameterize 
( TTree*                     tree        ,  
  Ostap::Math::LegendreSum3& sum         , 
  const std::string&         xexpression , 
  const std::string&         yexpression , 
  const std::string&         zexpression , 
  const std::string&         selection   , 
  const unsigned long        first       ,
  const unsigned long        last        ) 
{
  // invalid tree or emptry range 
  if ( nullptr == tree || last <= first ) { return 0 ; }
  //
  if ( 0 == selection.size() ) 
  { return parameterize ( tree , sum , 
                          xexpression , yexpression , zexpression , first , last ) ; }
  //
  Ostap::Formula xvar ( xexpression , tree ) ;
  Ostap::Assert ( xvar.ok()                                     , 
                  "Invalid x-expression:\"" + xexpression + "\"" ,
                  "Ostap::DataParams::parameterize  "           ) ;
  Ostap::Formula yvar ( yexpression , tree ) ;
  Ostap::Assert ( yvar.ok()                                     , 
                  "Invalid y-expression:\"" + yexpression + "\"" ,
                  "Ostap::DataParams::parameterize  "           ) ;
  Ostap::Formula zvar ( zexpression , tree ) ;
  Ostap::Assert ( zvar.ok()                                     , 
                  "Invalid z-expression:\"" + zexpression + "\"" ,
                  "Ostap::DataParams::parameterize  "           ) ;

  Ostap::Formula weight  ( selection , tree ) ;
  Ostap::Assert ( weight.ok()                                 , 
                  "Invalid selection:\"" + selection   + "\"" ,
                  "Ostap::DataParams::parameterize"           ) ;
  //
  Ostap::Utils::Notifier notify ( tree , &xvar, &yvar , &zvar , &weight ) ;
  //
  const unsigned long nEntries = std::min ( last , (unsigned long) tree->GetEntries() ) ;
  //
  long double sumw = 0 ;
  for ( unsigned long entry = first ; entry < nEntries ; ++entry ) 
  { 
    long ievent = tree->GetEntryNumber ( entry ) ;
    if ( 0 > ievent ) { break ; }                        // BREAK
    //
    ievent      = tree->LoadTree ( ievent ) ;
    if ( 0 > ievent ) { break ; }                        // BREAK
    //
    // get the weight  
    const double w = weight.evaluate () ;
    if ( !w        ) { continue ; }                      // CONTINUE 
    //
    // get the variable   
    const double x = xvar.evaluate() ;
    const double y = yvar.evaluate() ;
    const double z = zvar.evaluate() ;
    //
    // fill the sum 
    sumw += ( sum.fill ( x , y , z , w ) ? w : 0 ) ;
  }
  //
  return sumw ;
}
// ============================================================================
/*  fill Legendre sum with data from the Tree 
 *  @see Ostap::Math::LegendreSum4 
 *  @see Ostap::Math::LegendreSum4::fill
 *  @param tree        (INPUT)  the input tree 
 *  @param sum         (UPDATE) the parameterization object 
 *  @param xexpression (INPUT)  x-expression to be parameterized
 *  @param yexpression (INPUT)  y-expression to be parameterized
 *  @param zexpression (INPUT)  z-expression to be parameterized
 *  @param uexpression (INPUT)  u-expression to be parameterized
 *  @param first       (INPUT)  the first event in Tree 
 *  @param last        (INPUT)  the last  event in Tree 
 *  @return number of events used in parameterization 
 *  @code
 *  Tree*  tree = ...
 *  LegendreSum4 s ( 5 , 3 , 2 , 2 , -1 , 1 , -2 , 2 , 0 , 4 , 0 , 1 ) ;
 *  DataParam::parameterize ( tree , s , "x" , "y" , "z" , "t" ) ;
 *  @endcode
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2019-07-3
 */ 
// ============================================================================
unsigned long Ostap::DataParam::parameterize 
( TTree*                     tree        , 
  Ostap::Math::LegendreSum4& sum         , 
  const std::string&         xexpression , 
  const std::string&         yexpression , 
  const std::string&         zexpression , 
  const std::string&         uexpression , 
  const unsigned long        first       ,
  const unsigned long        last        ) 
{
  // invalid tree or emptry range 
  if ( nullptr == tree || last <= first ) { return 0 ; }
  //
  Ostap::Formula xvar ( xexpression , tree ) ;
  Ostap::Assert ( xvar.ok()                                      , 
                  "Invalid x-expression:\"" + xexpression + "\"" ,
                  "Ostap::DataParams::parameterize  "            ) ;
  Ostap::Formula yvar ( yexpression , tree ) ;
  Ostap::Assert ( yvar.ok()                                      , 
                  "Invalid y-expression:\"" + yexpression + "\"" ,
                  "Ostap::DataParams::parameterize  "            ) ;
  Ostap::Formula zvar ( zexpression , tree ) ;
  Ostap::Assert ( zvar.ok()                                      , 
                  "Invalid z-expression:\"" + zexpression + "\"" ,
                  "Ostap::DataParams::parameterize  "            ) ;
  Ostap::Formula uvar ( uexpression , tree ) ;
  Ostap::Assert ( uvar.ok()                                      , 
                  "Invalid u-expression:\"" + uexpression + "\"" ,
                  "Ostap::DataParams::parameterize  "            ) ;
  //
  Ostap::Utils::Notifier notify ( tree , &xvar, &yvar , &zvar , &uvar ) ;
  //
  const unsigned long nEntries = std::min ( last , (unsigned long) tree->GetEntries() ) ;
  unsigned long filled = 0 ;
  for ( unsigned long entry = first ; entry < nEntries ; ++entry ) 
  { 
    long ievent = tree->GetEntryNumber ( entry ) ;
    if ( 0 > ievent ) { break ; }                        // BREAK
    //
    ievent      = tree->LoadTree ( ievent ) ;
    if ( 0 > ievent ) { break ; }                        // BREAK
    //
    // get the variable   
    const double x = xvar.evaluate() ;
    const double y = yvar.evaluate() ;
    const double z = zvar.evaluate() ;
    const double u = uvar.evaluate() ;
    //
    // fill the sum 
    filled += ( sum.fill ( x , y  , z , u , 1 ) ? 1 : 0 ) ;
  }
  //
  return filled ;
}
// ============================================================================
/*  fill Legendre sum with data from the Tree 
 *  @see Ostap::Math::LegendreSum4 
 *  @see Ostap::Math::LegendreSum4::fill
 *  @param tree        (INPUT)  the input tree 
 *  @param sum         (UPDATE) the parameterization object 
 *  @param xexpression (INPUT)  x-expression to be parameterized
 *  @param yexpression (INPUT)  y-expression to be parameterized
 *  @param zexpression (INPUT)  z-expression to be parameterized
 *  @param uexpression (INPUT)  u-expression to be parameterized
 *  @param selection   (INPUT)  selection/weight to be used 
 *  @param first       (INPUT)  the first event in Tree 
 *  @param last        (INPUT)  the last  event in Tree 
 *  @return  sum of weigths  used in parameterization
 *  @code
 *  Tree*  tree = ...
 *  LegendreSum4 s ( 5 , 2 , 4 ,  -1 , 1 , -4  , -5 , 0 , 3 ) ;
 *  DataParam::parameterize ( tree , s , "x"  , "y" , "z" , "t" , "q>10" ) ;
 *  @endcode
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2019-07-3
 */ 
// ============================================================================
double Ostap::DataParam::parameterize 
( TTree*                     tree        ,  
  Ostap::Math::LegendreSum4& sum         , 
  const std::string&         xexpression , 
  const std::string&         yexpression , 
  const std::string&         zexpression , 
  const std::string&         uexpression , 
  const std::string&         selection   , 
  const unsigned long        first       ,
  const unsigned long        last        ) 
{
  // invalid tree or emptry range 
  if ( nullptr == tree || last <= first ) { return 0 ; }
  //
  if ( 0 == selection.size() ) 
  { return parameterize ( tree , sum , 
                          xexpression , yexpression , 
                          zexpression , uexpression , first , last ) ; }
  //
  Ostap::Formula xvar ( xexpression , tree ) ;
  Ostap::Assert ( xvar.ok()                                     , 
                  "Invalid x-expression:\"" + xexpression + "\"" ,
                  "Ostap::DataParams::parameterize  "           ) ;
  Ostap::Formula yvar ( yexpression , tree ) ;
  Ostap::Assert ( yvar.ok()                                     , 
                  "Invalid y-expression:\"" + yexpression + "\"" ,
                  "Ostap::DataParams::parameterize  "           ) ;
  Ostap::Formula zvar ( zexpression , tree ) ;
  Ostap::Assert ( zvar.ok()                                     , 
                  "Invalid z-expression:\"" + zexpression + "\"" ,
                  "Ostap::DataParams::parameterize  "           ) ;
  Ostap::Formula uvar ( uexpression , tree ) ;
  Ostap::Assert ( uvar.ok()                                     , 
                  "Invalid u-expression:\"" + uexpression + "\"" ,
                  "Ostap::DataParams::parameterize  "           ) ;

  Ostap::Formula weight ( selection , tree ) ;
  Ostap::Assert ( weight.ok()                                 , 
                  "Invalid selection:\"" + selection   + "\"" ,
                  "Ostap::DataParams::parameterize"           ) ;
  //
  Ostap::Utils::Notifier notify ( tree , &xvar, &yvar , &zvar , &uvar , &weight ) ;
  //
  const unsigned long nEntries = std::min ( last , (unsigned long) tree->GetEntries() ) ;
  //
  long double sumw = 0 ;
  for ( unsigned long entry = first ; entry < nEntries ; ++entry ) 
  { 
    long ievent = tree->GetEntryNumber ( entry ) ;
    if ( 0 > ievent ) { break ; }                        // BREAK
    //
    ievent      = tree->LoadTree ( ievent ) ;
    if ( 0 > ievent ) { break ; }                        // BREAK
    //
    // get the weight  
    const double w = weight.evaluate () ;
    if ( !w        ) { continue ; }                      // CONTINUE 
    //
    // get the variable   
    const double x = xvar.evaluate() ;
    const double y = yvar.evaluate() ;
    const double z = zvar.evaluate() ;
    const double u = uvar.evaluate() ;
    //
    // fill the sum 
    sumw += ( sum.fill ( x , y , z , u , w ) ? w : 0 ) ;
  }
  //
  return sumw ;
}
// ============================================================================
/*  fill chebyshev sum with data from the Tree 
 *  @see Ostap::Math::ChebyshevSum 
 *  @see Ostap::Math::ChebyshevSum::fill
 *  @param tree       (INPUT)  the input tree 
 *  @param sum        (UPDATE) the parameterization object 
 *  @param expression (INPUT)  expresion to be parameterized
 *  @param first      (INPUT)  the first event in Tree 
 *  @param last       (INPUT)  the last  event in Tree 
 *  @return number of events sused in parameterization 
 *  @code
 *  Tree*  tree = ...
 *  ChebyshevSum s ( 5 , -1 , 1 ) ;
 *  DataParam::parameterize ( tree , s , "x" ) ;
 *  @endcode
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2019-07-3
 */ 
// ============================================================================
unsigned long Ostap::DataParam::parameterize 
( TTree*                     tree       , 
  Ostap::Math::ChebyshevSum& sum        , 
  const std::string&         expression , 
  const unsigned long        first      ,
  const unsigned long        last       ) 
{
  // invalid tree or emptry range 
  if ( nullptr == tree || last <= first ) { return 0 ; }
  //
  Ostap::Formula var ( expression , tree ) ;
  Ostap::Assert ( var.ok()                                    , 
                  "Invalid expression:\"" + expression + "\"" ,
                  "Ostap::DataParams::parameterize"           ) ;
  //
  Ostap::Utils::Notifier notify ( tree , &var ) ;
  //
  const unsigned long nEntries = std::min ( last , (unsigned long) tree->GetEntries() ) ;
  unsigned long filled = 0 ;
  for ( unsigned long entry = first ; entry < nEntries ; ++entry ) 
  { 
    long ievent = tree->GetEntryNumber ( entry ) ;
    if ( 0 > ievent ) { break ; }                        // BREAK
    //
    ievent      = tree->LoadTree ( ievent ) ;
    if ( 0 > ievent ) { break ; }                        // BREAK
    //
    // get the variable   
    const double v = var.evaluate() ;
    //
    // fill the sum 
    filled += ( sum.fill ( v , 1 ) ? 1 : 0 ) ;
  }
  //
  return filled ;
}
// ==================================================================================
/*  fill chebyshev sum with data from the Tree 
 *  @see Ostap::Math::ChebyshevSum 
 *  @see Ostap::Math::ChebyshevSum::fill
 *  @param tree       (INPUT)  the input tree 
 *  @param sum        (UPDATE) the parameterization object 
 *  @param expression (INPUT)  expression to be parameterized
 *  @param seelction  (INPUT)  selection/weight to be used 
 *  @param first      (INPUT)  the first event in Tree 
 *  @param last       (INPUT)  the last  event in Tree 
 *  @return  sum of weigths  used in parameterization
 *  @code
 *  Tree*  tree = ...
 *  ChebyshevSum s ( 5 , -1 , 1 ) ;
 *  DataParam::parameterize ( tree , s , "x"  , "y>10" ) ;
 *  @endcode
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2019-07-3
 */ 
// ==================================================================================
double Ostap::DataParam::parameterize 
( TTree*                     tree       , 
  Ostap::Math::ChebyshevSum& sum        , 
  const std::string&         expression , 
  const std::string&         selection  , 
  const unsigned long        first      ,
  const unsigned long        last       ) 
{
  // invalid tree or emptry range 
  if ( nullptr == tree || last <= first ) { return 0 ; }
  //
  if ( 0 == selection.size() ) 
  { return parameterize ( tree ,  sum , expression , first , last ) ; }
  //
  //
  Ostap::Formula var ( expression , tree ) ;
  Ostap::Assert ( var.ok()                                    , 
                  "Invalid expression:\"" + expression + "\"" ,
                  "Ostap::DataParams::parameterize"           ) ;
  //
  Ostap::Formula weight ( selection , tree ) ;
  Ostap::Assert ( weight.ok()                                 , 
                  "Invalid selection:\"" + selection   + "\"" ,
                  "Ostap::DataParams::parameterize"           ) ;
  //
  Ostap::Utils::Notifier notify ( tree , &var ,  &weight ) ;
  //
  const unsigned long nEntries = std::min ( last , (unsigned long) tree->GetEntries() ) ;
  //
  long double wsum = 0 ;
  for ( unsigned long entry = first ; entry < nEntries ; ++entry ) 
  { 
    long ievent = tree->GetEntryNumber ( entry ) ;
    if ( 0 > ievent ) { break ; }                        // BREAK
    //
    ievent      = tree->LoadTree ( ievent ) ;
    if ( 0 > ievent ) { break ; }                        // BREAK
    //
    // get the weight  
    const double w = weight.evaluate () ;
    if ( !w        ) { continue ; }
    //
    // get the variable   
    const double v =    var.evaluate() ;
    //
    // fill the sum 
    wsum += ( sum.fill ( v , w ) ? w : 0 ) ;
  }
  //
  return wsum ;
}


// ============================================================================
//                                                                      The END 
// ============================================================================
