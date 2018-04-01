// ============================================================================
#ifndef OSTAP_FUNCS_H 
#define OSTAP_FUNCS_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <string>
#include <memory>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/IFuncs.h"
// ============================================================================
class RooFormulaVar ; // formm RooFit 
// ============================================================================
namespace Ostap 
{
  // ==========================================================================
  /// from Ostap 
  class  Formula ;
  // ==========================================================================
  namespace Functions 
  {
    // ========================================================================
    /** @class FuncFormula 
     *  simple implementation of TTRee-function based on Ostap::Formula
     */
    class FuncFormula : public Ostap::IFuncTree 
    {
    public :
      // ======================================================================
      /** constructor from the formula expression 
       *  @param expression the  formula expression 
       *  @param tree       the tree 
       *  @param name       the name for the formula 
       */
      FuncFormula ( const std::string& expression            , 
                    const TTree*       tree       =  nullptr ,
                    const std::string& name       = ""       ) ;
      // ======================================================================
    public:
      // ======================================================================
      ///  evaluate the formula for  TTree
      double evaluate ( const TTree* tree = nullptr ) const override ;
      // ======================================================================
    private:
      // ======================================================================
      mutable const TTree*                    m_tree    ;
      mutable std::unique_ptr<Ostap::Formula> m_formula ;
      // ======================================================================
      /// the  expression itself 
      std::string m_expression ; // the  expression itself 
      /// the name  
      std::string m_name       ; // the name  
      // ======================================================================
    } ;
    // ========================================================================
    /** @class FuncRooFormula
     *  simple implementation of 'RooAbsData'-function based on RooFormulaVar
     */
    class FuncRooFormula : public Ostap::IFuncData
    {
    public :
      // ======================================================================
      /** constructor from the formula expression 
       *  @param expression the formula expression 
       *  @param data       the data
       *  @param name       the name for the formula 
       */
      FuncRooFormula ( const std::string& expression            , 
                       const RooAbsData*  data       =  nullptr ,
                       const std::string& name       = ""       ) ;
      // ======================================================================
    public:
      // ======================================================================
      ///  evaluate the formula for  data
      double evaluate ( const RooAbsData* data = nullptr ) const override ;
      // ======================================================================
    private:
      // ======================================================================
      mutable const RooAbsData*              m_data    ;
      mutable std::unique_ptr<RooFormulaVar> m_formula ;
      // ======================================================================
      /// the  expression itself 
      std::string m_expression ; // the  expression itself 
      /// the name  
      std::string m_name       ; // the name  
      // ======================================================================
    } ;
    // ========================================================================
  } //                                   The END of  namespace Ostap::Functions
  // ==========================================================================
} //                                                The END of  namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_FUNCS_H
// ============================================================================

