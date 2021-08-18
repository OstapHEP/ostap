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
  const std::string vname { Ostap::tmp_name ( "formula_" , formula ) } ;
  std::unique_ptr<RooFormula> ptr { new RooFormula ( vname  .c_str () ,
                                                     formula.c_str () , 
                                                     variables        , 
                                                     false            ) } ;
  
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
  if ( !formula.ok() ) { return used; }   // RETURN 
  //
  const RooArgSet actual { formula.actualDependents() } ;
  //
#if ROOT_VERSION_CODE < ROOT_VERSION(6,18,0)
  //   
  Ostap::Utils::Iterator tmp ( from ) ;
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
{ return usedVariables ( formula.formula() , variables ) ; }
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
