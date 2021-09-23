// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <set>
#include <string>
#include <memory>
#include <regex>
#include <sstream>
// ============================================================================
// ROOT/RooFit 
// ============================================================================
#include "RVersion.h"
#include "RooFormulaVar.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/FormulaVar.h"
#include "Ostap/Iterator.h"
// ============================================================================
// local
// ============================================================================
#include "Exception.h"
#include "local_utils.h"
// ============================================================================
/** @file 
 *  Implementation file for class Ostap::FormulaVar
 *  @date 2020-03-06
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 */
// ============================================================================
ClassImp(Ostap::FormulaVar) ;
// ============================================================================
/** make formula (skip unnesessary dependents)
 *  @param name  formula name 
 *  @param title formula title 
 *  @param expression formula expression
 *  @param dependent  formula dependents 
 *  @return the formula  
 */
// ============================================================================
std::unique_ptr<Ostap::FormulaVar>  
Ostap::makeFormula 
( const std::string& name       , 
  const std::string& title      , 
  const std::string& expression , 
  const RooArgList & dependents ) 
{
  const std::string vname { Ostap::tmp_name ( "formula1_" , expression ) } ;
  //
#if ROOT_VERSION_CODE < ROOT_VERSION(6,20,0)
  //
  std::unique_ptr<RooFormula> ptr { new RooFormula ( vname      .c_str () ,
                                                     expression.c_str  () , 
                                                     dependents         ) };
  //
#else 
  //
  std::unique_ptr<RooFormula> ptr { new RooFormula ( vname     .c_str () ,
                                                     expression.c_str () , 
                                                     dependents          , 
                                                     false            ) };
  //
#endif
  //
  if ( !ptr || !ptr->ok() ) { return nullptr ; }
  //
  std::unique_ptr<Ostap::FormulaVar> result 
  { new Ostap::FormulaVar ( name       , 
                            title      , 
                            expression , 
                            usedVariables ( *ptr , dependents ) , 
                            true       ) } ;
  //
  if ( !result || !result->ok() ) { return nullptr ; }
  //
  return result ;
}
// ============================================================================
/** make formula (skip unnesessary dependents)
 *  @param name  formula name 
 *  @param expression formula expression
 *  @param dependent  formula dependents 
 *  @return the formula  
 */
// ============================================================================
std::unique_ptr<Ostap::FormulaVar>  
Ostap::makeFormula 
( const std::string& name       , 
  const std::string& expression , 
  const RooArgList & dependents ) 
{ return makeFormula ( name , expression , expression , dependents ) ; }
// ============================================================================
/** make formula (skip unnesessary dependents)
 *  @param expression formula expression
 *  @param dependent  formula dependents 
 *  @return the formula  
 */
// ============================================================================
std::unique_ptr<Ostap::FormulaVar>  
Ostap::makeFormula 
( const std::string& expression , 
  const RooArgList & dependents ) 
{ return makeFormula ( Ostap::tmp_name (  "formula_" , expression ) ,
                        expression , expression  , dependents ) ; }
// ============================================================================
/* get the list of used variables for the given formula
 *  @code
 *  @endcode 
 *  @param formula the actual formula 
 *  @param variables list of variables 
 *  @return list of actually used variables (subset of <code>varibales</code>
 *  @see RooFormula 
 */
// ============================================================================
RooArgList 
Ostap::usedVariables
( const std::string& formula    , 
  const RooArgList&  variables  ) 
{
  //
  const std::string vname { Ostap::tmp_name ( "formula2_" , formula ) } ;
  //
#if ROOT_VERSION_CODE < ROOT_VERSION(6,20,0)
  //
  std::unique_ptr<RooFormula> ptr { new RooFormula ( vname  .c_str () ,
                                                     formula.c_str () , 
                                                     variables        ) };
  //
#else 
  //
  std::unique_ptr<RooFormula> ptr { new RooFormula ( vname  .c_str () ,
                                                     formula.c_str () , 
                                                     variables        , 
                                                     false            ) };
  //
#endif
  //
  if ( !ptr || !ptr->ok() ) { return RooArgList() ; }
  //
  return usedVariables ( *ptr , variables ) ;
}
// ============================================================================
/*  get the list of used variables for the given formula
 *  @code
 *  @endcode 
 *  @param formula the actual formula 
 *  @param variables list of variables 
 *  @return list of actually used variables (subset of <code>varibales</code>
 *  @see RooFormula 
 */
// ============================================================================
RooArgList 
Ostap::usedVariables
( const RooFormula&  formula    , 
  const RooArgList&  variables  ) 
{
  RooArgList used {};
  //
#if ROOT_VERSION_CODE < ROOT_VERSION(6,22,0)
  if ( !const_cast<RooFormula&>(formula).ok() ) { return used; }   // RETURN 
#else 
  if ( !formula.ok()                          ) { return used; }   // RETURN 
#endif
  //
  const RooArgSet actual { formula.actualDependents() } ;
  //
#if ROOT_VERSION_CODE < ROOT_VERSION(6,18,0)
  //   
  Ostap::Utils::Iterator tmp ( variables ) ;
  RooAbsArg* c = 0 ;
  while ( c = (RooAbsArg*) tmp.next() )
  { if ( c && actual.contains ( *c ) ) { used.add ( *c ) ; } }
  //
#else
  //
  for ( const auto* arg : variables ) 
  { if ( arg && actual.contains ( *arg ) ) { used.add ( *arg ) ; } }
  //
#endif 
  //
  return used ;
}
// ============================================================================
/*  get the list of used variables for the given formula
 *  @code
 *  @endcode 
 *  @param formula the actual formula 
 *  @param variables list of variables 
 *  @return list of actually used variables (subset of <code>varibales</code>
 *  @see RooFormulaVar
 *  @see RooFormulaVar::formula  
 */
// ============================================================================
RooArgList 
Ostap::usedVariables
( const RooFormulaVar& formula    , 
  const RooArgList&    variables  )
{
  //
#if ROOT_VERSION_CODE < ROOT_VERSION(6,20,0)
  //
  std::ostringstream os {} ;
  formula.printMetaArgs( os ) ;
  std::string expression {os.str() } ;
  //
  std::size_t pos = expression.find ( "formula=\"" ) ;
  Ostap::Assert ( 0 == pos ,
                  "Invalid formula expression/1: " + expression ,
                  "Ostap::usedVariables" );
  //
  expression = expression.substr ( 9 ) ;
  pos = expression.find('"');
  Ostap::Assert ( 0 <= pos && pos != std::string::npos , 
                  "Invalid formula expression/2: " + expression ,
                  "Ostap::usedVariables" );
  expression = expression.substr ( 0 , pos ) ;
  // 
  std::unique_ptr<RooFormula> ptr { new RooFormula ( formula.GetName  () ,
                                                     expression.c_str () , 
                                                     variables           ) } ;  
  if ( !ptr || !ptr->ok() ) { return RooArgList() ; }
  return usedVariables ( *ptr , variables ) ;
  //
#elif ROOT_VERSION_CODE < ROOT_VERSION(6,22,0)
  //
  std::ostringstream os {} ;
  formula.printMetaArgs( os ) ;
  std::string expression {os.str() } ;
  //
  std::size_t pos = expression.find ( "formula=\"" ) ;
  Ostap::Assert ( 0 == pos ,
                  "Invalid formula expression/1: " + expression ,
                  "Ostap::usedVariables" );
  //
  expression = expression.substr ( 9 ) ;
  pos = expression.find('"');
  Ostap::Assert ( 0 <= pos && pos != std::string::npos , 
                  "Invalid formula expression/2: " + expression ,
                  "Ostap::usedVariables" );
  expression = expression.substr ( 0 , pos ) ;
  // 
  std::unique_ptr<RooFormula> ptr { new RooFormula ( formula.GetName  () ,
                                                     expression.c_str () , 
                                                     variables           , 
                                                     false               ) } ;  
  if ( !ptr || !ptr->ok() ) { return RooArgList() ; }
  return usedVariables ( *ptr , variables ) ;
  //
#else 
  //
  return usedVariables ( formula.formula() , variables ) ; 
  //
#endif
}
// ============================================================================
/*  full constructor 
 *  @param name       formula name 
 *  @param title      formula title 
 *  @param expression formula expression 
 *  @param dependent  formula dependents 
 *  @param check      check dependents?
 */
// ============================================================================
Ostap::FormulaVar::FormulaVar
( const std::string& name          , 
  const std::string& title         , 
  const std::string& expression    , 
  const RooArgList & dependents    ,
  const bool 	      check          ) 
  : RooFormulaVar ( name       . c_str () , 
                    title      . c_str () , 
                    expression . c_str () , 
                    dependents            
#if    ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
                    , check
#endif 
                    ) 
{}
// ============================================================================
/*  full constructor 
 *  @param name       formula name 
 *  @param expression formula expression 
 *  @param dependent  formula dependents 
 *  @param check      check dependents?
 */
// ============================================================================
Ostap::FormulaVar::FormulaVar
( const std::string& name       , 
  const std::string& expression , 
  const RooArgList & dependents ,
  const bool 	      check       ) 
  : RooFormulaVar ( name       . c_str () , 
                    expression . c_str () , 
                    dependents            
#if    ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
                    , check
#endif 
                    ) 
{}
// ============================================================================
/*  full constructor 
 *  @param expression formula expression 
 *  @param dependent  formula dependents 
 *  @param check      check dependents?
 */
// ============================================================================
Ostap::FormulaVar::FormulaVar
( const std::string& expression , 
  const RooArgList & dependents ,
  const bool 	      check       ) 
  : RooFormulaVar ( Ostap::tmp_name ( "formula_" , expression ) . c_str () , 
                    expression                                  . c_str () , 
                    dependents            
#if    ROOT_VERSION_CODE >= ROOT_VERSION(6,20,0)
                    , check
#endif 
                    ) 
{}
// ============================================================================
// copy constuctor 
// ============================================================================
Ostap::FormulaVar::FormulaVar
( const Ostap::FormulaVar& right , 
  const char*              name  ) 
  : RooFormulaVar ( right , name ) 
{}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::FormulaVar::FormulaVar
( const RooFormulaVar&     right , 
  const char*              name  ) 
 : RooFormulaVar ( right , name  ) 
{}
// ============================================================================
//  destructor 
// ============================================================================
Ostap::FormulaVar::~FormulaVar(){}
// ============================================================================

// ============================================================================
//                                                                      The END 
// ============================================================================
