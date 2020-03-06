// ============================================================================
#ifndef OSTAP_FORMULAVAR_H 
#define OSTAP_FORMULAVAR_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <string>
// ============================================================================
// ROOT/RooFit 
// ============================================================================
#include "RooFormulaVar.h"
// ============================================================================
namespace Ostap 
{
  // ==========================================================================
  /** @class FormulaVar Ostap/FormulaVar.h
   *  Tiny extension of class FormulaVar 
   *  @author Vanya Belyaev
   *  @date   2020-03-06
   */
  class FormulaVar : public RooFormulaVar
  { 
  public:
    // ========================================================================
    ClassDef ( Ostap::FormulaVar , 1 ) ;
    // ========================================================================    
  public:
    // ========================================================================
    /** full constructor 
     *  @param name       formula name 
     *  @param title      formula title 
     *  @param expression formula expression 
     *  @param dependent  formula dependents 
     *  @param check      check dependents?
     */
    FormulaVar ( const std::string& name         , 
                 const std::string& title        , 
                 const std::string& expression   , 
                 const RooArgList & dependents   ,
                 const bool 	      check = true ) ;
    // ========================================================================
    /** full constructor 
     *  @param name       formula name 
     *  @param expression formula expression 
     *  @param dependent  formula dependents 
     *  @param check      check dependents?
     */
    FormulaVar ( const std::string& name         , 
                 const std::string& expression   , 
                 const RooArgList & dependents   ,
                 const bool 	      check = true ) ;
    // ========================================================================
    /** full constructor 
     *  @param expression formula expression 
     *  @param dependent  formula dependents 
     *  @param check      check dependents?
     */    
    FormulaVar ( const std::string& name         , 
                 const RooArgList & dependents   ,
                 const bool 	      check = true ) ;
    // ========================================================================
    /// copy constructor 
    FormulaVar ( const FormulaVar&    right , 
                 const char*          name = nullptr ) ;
    // =========================================================================
    /// copy constructor 
    FormulaVar ( const RooFormulaVar& right , 
                 const char*          name = nullptr ) ;
    // =========================================================================
    ///  destructor 
    virtual ~FormulaVar() ;
    // =========================================================================
  };
  // ===========================================================================
} //                                                  The end of namespace Ostap 
// =============================================================================
//                                                                       The END 
// =============================================================================
#endif // OSTAP_FORMULAVAR_H
// =============================================================================
