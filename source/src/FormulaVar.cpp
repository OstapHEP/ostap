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
#include <exception>
// ============================================================================
// ROOT/RooFit 
// ============================================================================
#include "RVersion.h"
#include "RooFormulaVar.h"
#include "TError.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/FormulaVar.h"
#include "Ostap/Iterator.h"
#include "Ostap/Mute.h"
// ============================================================================
// local
// ============================================================================
#include "Exception.h"
#include "local_utils.h"
#include "local_roofit.h"
// ============================================================================
/** @file 
 *  Implementation file for class Ostap::FormulaVar
 *  @date 2020-03-06
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 */
// ============================================================================
ClassImp(Ostap::FormulaVar) ;
// ============================================================================
namespace
{
  // ==========================================================================
  class ESentry
  {
  public :
    ESentry  () ;
    ~ESentry () ;
  private:
    int m_level ;    
  };
  //===========================================================================    
  ESentry:: ESentry()
    : m_level ( gErrorIgnoreLevel )
  {
    gErrorIgnoreLevel = kFatal + 1 ;
  }  
  ESentry::~ESentry()
  {
    gErrorIgnoreLevel = m_level ;    
  }
  // ==========================================================================
#if ROOT_VERSION_CODE < ROOT_VERSION(6,29,0)
  // ==========================================================================
  /** get the list of used variables for the given formula
   *  @code
   *  @endcode 
   *  @param formula the actual formula 
   *  @param variables list of variables 
   *  @return list of actually used variables (subset of <code>variables</code>
   *  @see RooFormula 
   *  @see RooFormula::actualDependents 
   */
  RooArgList 
  usedVariables
  ( const RooFormula&       formula    , 
    const RooArgList&       variables  ) 
  {
    RooArgList used {};
    //
    if ( !formula.ok()                          ) { return used; }   // RETURN 
    //
    const RooArgSet actual { formula.actualDependents() } ;
    //
    for ( const auto* arg : variables ) 
    { if ( arg && actual.contains ( *arg ) ) { used.add ( *arg ) ; } }
    //
    return used ;
  }
  //
#endif
  // ==========================================================================
}
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
( const std::string&      name        ,
  const std::string&      title       ,
  const std::string&      expression  , 
  const RooArgList&       dependents  ) 
{
  //
  const std::string vname { Ostap::tmp_name ( "test_formula1_" , expression ) } ;
  //
  // ========================================================================
#if ROOT_VERSION_CODE < ROOT_VERSION(6,29,0)
  //
  std::unique_ptr<RooFormula> ptr ;
  //
  try
  {
    ESentry sentry {} ;
    ptr.reset ( new RooFormula ( vname     .c_str () ,                                                        
                                 expression.c_str () , 
                                 dependents          , 
                                 false             ) ) ;
  }
  catch ( std::invalid_argument& /* e1 */ ){ return nullptr ;  }
  catch ( std::runtime_error&    /* e2 */ ){ return nullptr ;  }
  //
  if ( !ptr || !ptr->ok() ) { return nullptr ; }
  const RooArgList used { ::usedVariables ( *ptr , dependents ) } ;
  //
#else
  //
  std::unique_ptr<RooFormulaVar> ptr ;
  //
  try
  {
    ESentry sentry {} ;
    Ostap::Utils::Mute m1  ( true  ) ;
    Ostap::Utils::Mute m2  ( false ) ;    
    ptr.reset ( new RooFormulaVar ( vname     .c_str () ,                                                        
                                    expression.c_str () , 
                                    expression.c_str () , 
                                    dependents          , 
                                    false             ) ) ;
  }
  catch ( std::invalid_argument& /* e1 */ ){ return nullptr ;  }
  catch ( std::runtime_error&    /* e2 */ ){ return nullptr ;  }
  //    
  if ( !ptr || !ptr->ok() ) { return nullptr ; }
  const RooArgList used { usedVariables ( *ptr , dependents ) } ;
  // ========================================================================
#endif
  // ==========================================================================
  //
  if ( !ptr || !ptr->ok() ) { return nullptr ; }
  //
  std::unique_ptr<Ostap::FormulaVar> result
                        ( new Ostap::FormulaVar ( name       , 
                                                  title      , 
                                                  expression , 
                                                  used       , 
                                                  true       ) ) ;
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
( const std::string&   name       , 
  const std::string&   expression , 
  const RooArgList&    dependents ) 
{ return makeFormula ( name       ,
                       expression , 
                       expression , 
                       dependents ) ; }
// ============================================================================
/** make formula (skip unnesessary dependents)
 *  @param expression formula expression
 *  @param dependent  formula dependents 
 *  @return the formula  
 */
// ============================================================================
std::unique_ptr<Ostap::FormulaVar>  
Ostap::makeFormula 
( const std::string&   expression ,
  const RooArgList&    dependents ) 
{ return makeFormula ( Ostap::tmp_name (  "formula_" , expression ) ,
                       expression , 
                       expression , 
                       dependents ) ; }
// ============================================================================
/** valid formula expression ?
 *  @param expression formula expression
 *  @param dependent  formula dependents 
 *  @return true for valid formula
 */
// ============================================================================
bool 
Ostap::validFormula 
( const std::string& expression , 
  const RooArgList&  dependents ) 
{
  const std::unique_ptr<Ostap::FormulaVar> ptr
  { makeFormula ( expression , dependents ) } ;
  return ptr && ptr->ok () ;
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
( const std::string& formula    ,
  const RooArgList&  variables  ) 
{
  //
  const std::string vname { Ostap::tmp_name ( "formula2_" , formula ) } ;
  //
#if ROOT_VERSION_CODE < ROOT_VERSION(6,29,0)
  //
  std::unique_ptr<RooFormula> ptr ;
  try
  {
    ESentry sentry {} ;
    ptr.reset ( new RooFormula ( vname  .c_str () ,
                                 formula.c_str () , 
                                 variables        , 
                                 false            ) ) ;
  }
  catch ( std::invalid_argument& /* e1 */ ){ return RooArgList () ; }
  catch ( std::runtime_error&    /* e2 */ ){ return RooArgList () ; }
  //
  if ( !ptr || !ptr->ok() ) { return RooArgList() ; }
  return ::usedVariables ( *ptr , variables ) ;
  //
#else 
  //
  ESentry sentry {} ;
  std::unique_ptr<RooFormulaVar> ptr ;
  try
  {
    ESentry sentry {} ;
    ptr.reset ( new RooFormulaVar ( vname.c_str () ,
                                    formula.c_str () , 
                                    variables        , 
                                    false            ) ) ;
  }
  catch ( std::invalid_argument& /* e1 */ ){ return RooArgList () ; }
  catch ( std::runtime_error&    /* e2 */ ){ return RooArgList () ; }
  //
  if ( !ptr || !ptr->ok() ) { return RooArgList() ; }
  return usedVariables ( *ptr , variables ) ;
  //
#endif
  //
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
( const RooFormulaVar&    formula    , 
  const RooArgList&       variables  )
{
  //
#if ROOT_VERSION_CODE < ROOT_VERSION(6,29,0)
  //
  return ::usedVariables ( formula.formula() , variables ) ; 
  //
#else 
  // ==========================================================================
  /// FIX ME!!!! @see https://github.com/root-project/root/commit/a470a3d85e8b5c93917d3e84c39e9d5f0066da97
  return formula.dependents() ; // FIX ME!!!
  // ==========================================================================
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
#if  ROOT_VERSION_CODE < ROOT_VERSION(6,29,0)
  const std::string&    title      , 
#else        
  const std::string& /* title */   , 
#endif 
  const std::string& expression    , 
  const RooArgList & dependents    ,
  const bool            check         
  )
  : RooFormulaVar ( name       . c_str () 
#if  ROOT_VERSION_CODE < ROOT_VERSION(6,29,0)
  , title      . c_str () 
#else        
  , expression . c_str ()  
#endif 
  , expression . c_str () 
  , dependents            
  , check
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
  const bool 	       check      )
  : FormulaVar ( name ,  expression , expression ,  dependents , check )
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
  const bool         check      ) 
  : FormulaVar ( Ostap::tmp_name ( "formula_" , expression ) , 
                 expression , 
                 expression ,
                 dependents ,
                 check      ) 
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
// get true formula expression 
// ============================================================================
std::string Ostap::FormulaVar::expression () const 
{
#if ROOT_VERSION_CODE < ROOT_VERSION(6,29,0)
  return formula().GetTitle() ;
#else 
  // ==========================================================================
  /// FIX ME!!!! @see https://github.com/root-project/root/commit/a470a3d85e8b5c93917d3e84c39e9d5f0066da97
  return GetTitle ()  ;
  // ==========================================================================
#endif 
} 
// =========================================================================

// ============================================================================
//                                                                      The END 
// ============================================================================
