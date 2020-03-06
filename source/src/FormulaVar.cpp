// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <string>
// ============================================================================
// ROOT/RooFit 
// ============================================================================
#include "RVersion.h"
#include "RooFormulaVar.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/FormulaVar.h"
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
