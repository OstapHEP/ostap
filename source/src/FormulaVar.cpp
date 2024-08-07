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
    Ostap::Utils::Iterator tmp ( variables ) ; // only for ROOT < 6.18
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
#if ROOT_VERSION_CODE < ROOT_VERSION(6,20,0)
  //
  std::unique_ptr<RooFormula> ptr ;
  try
  {
    ESentry sentry {} ;
    ptr.reset ( new RooFormula ( vname     .c_str () ,
                                 expression.c_str () , 
                                 dependents          ) ) ;
  }
  catch ( std::invalid_argument& /* e1 */ ){ return nullptr ; }
  catch ( std::runtime_error&    /* e2 */ ){ return nullptr ; }
  //
  if ( !ptr || !ptr->ok() ) { return nullptr ; }
  const RooArgList used { ::usedVariables ( *ptr , dependents ) } ;
  //
#elif ROOT_VERSION_CODE < ROOT_VERSION(6,29,0)
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
#if ROOT_VERSION_CODE < ROOT_VERSION(6,20,0)
  //
  std::unique_ptr<RooFormula> ptr ;
  try
  {
    ESentry sentry {} ;
    ptr.reset (  new RooFormula ( vname  .c_str () ,
                                  formula.c_str () ,
                                  variables        ) ) ;
  }
  catch ( std::invalid_argument& /* e1 */ ){ return RooArgList () ; }
  catch ( std::runtime_error&    /* e2 */ ){ return RooArgList () ; }
  //
  if ( !ptr || !ptr->ok() ) { return RooArgList() ; }
  return ::usedVariables ( *ptr , variables ) ;
  //
#elif ROOT_VERSION_CODE < ROOT_VERSION(6,29,0)
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
  return ::usedVariables ( *ptr , variables ) ;
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
  return ::usedVariables ( *ptr , variables ) ;
  //
#elif ROOT_VERSION_CODE < ROOT_VERSION(6,29,0)
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
#if    ROOT_VERSION_CODE < ROOT_VERSION(6,22,0)
  const std::string& /* title */   , 
#elif  ROOT_VERSION_CODE < ROOT_VERSION(6,29,0)
  const std::string&    title      , 
#else        
  const std::string& /* title */   , 
#endif 
  const std::string& expression    , 
  const RooArgList & dependents    ,
#if ROOT_VERSION_CODE < ROOT_VERSION(6,20,0)
  const bool         /* check */   
#else 
  const bool            check         
#endif 
  )
  : RooFormulaVar ( name       . c_str () , 
#if    ROOT_VERSION_CODE < ROOT_VERSION(6,22,0)
  expression . c_str () , 
#elif  ROOT_VERSION_CODE < ROOT_VERSION(6,29,0)
  title      . c_str () , 
#else        
  expression . c_str () , 
#endif 
  expression . c_str () , 
  dependents            
#if ROOT_VERSION(6,20,0) <= ROOT_VERSION_CODE
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
#if   ROOT_VERSION_CODE < ROOT_VERSION(6,22,0)
  return           GetTitle() ;
#elif ROOT_VERSION_CODE < ROOT_VERSION(6,29,0)
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
