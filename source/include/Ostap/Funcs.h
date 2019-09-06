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
#include "Ostap/Formula.h"
// ============================================================================
// ROOT
// ============================================================================
#include  "TObject.h"
#include  "RooFormulaVar.h"
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
    /** @class FuncFormula Ostap/Funcs.h
     *  simple implementation of TTree-function based on Ostap::Formula
     */
    class FuncFormula : public Ostap::IFuncTree, public TObject
    {
    public :
      // ======================================================================
      ClassDef(Ostap::Functions::FuncFormula,1) ;
      // ======================================================================
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
      /// copy constructor 
      FuncFormula ( const  FuncFormula& right )  ;
      // ======================================================================
      /// default constructor, needed for serialization 
      FuncFormula () = default ;
      /// destructor 
      virtual ~FuncFormula() ;
      // ======================================================================
    public:
      // ======================================================================
      FuncFormula* Clone ( const char* newname = "" ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      ///  evaluate the formula for  TTree
      double operator() ( const TTree* tree ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      Bool_t Notify   () override ; 
      // ======================================================================
   private:
      // ======================================================================
      /// make formula 
      bool make_formula() const ;
      // ======================================================================
    private:
      // ======================================================================
      mutable const TTree*                    m_tree    { nullptr } ; //! 
      mutable std::unique_ptr<Ostap::Formula> m_formula { nullptr } ; //!
      // ======================================================================
      /// the  expression itself 
      std::string m_expression {} ; // the  expression itself 
      /// the name  
      std::string m_name       {} ; // the name  
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
      /// copy contructor
      FuncRooFormula ( const FuncRooFormula&  right ) ;
      // ======================================================================
      /// default constructor, needed for serialization 
      FuncRooFormula () = default ;
      // ======================================================================
      /// destructor 
      virtual ~FuncRooFormula() ;
      // ======================================================================
    public:
      // ======================================================================
      ///  evaluate the formula for  data
      double operator () ( const RooAbsData* data ) const override ;
      // ======================================================================
    private:
      // ======================================================================
      /// make formula 
      bool make_formula() const ;
      // ======================================================================
    private:
      // ======================================================================
      mutable const RooAbsData*              m_data    { nullptr } ;
      mutable std::unique_ptr<RooFormulaVar> m_formula { nullptr } ;
      // ======================================================================
      /// the  expression itself 
      std::string m_expression {} ; // the  expression itself 
      /// the name  
      std::string m_name       {} ; // the name  
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

