// ============================================================================
#ifndef OSTAP_FORMULA_H 
#define OSTAP_FORMULA_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL 
// ============================================================================
#include <memory>
#include <string>
// ============================================================================
// ROOT 
// ============================================================================
#include "TTreeFormula.h"
// ============================================================================
class TCut ; // ROOT 
// ============================================================================
/** file Ostap/Formula.h
 *  Simple extention of class TTreeFormula 
 *  for easier usage in python 
 */
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  /** @class Formula Ostap/Formula.h
   *  Simple extension of class TTreeFormula for a bit easier usage in python 
   *  @see TTreeFormula
   *  @author Vanya Belyaev
   *  @date   2013-05-06
   */
  class Formula : public TTreeFormula 
  {
    // ========================================================================
  public:
    // ========================================================================
    ClassDefOverride(Ostap::Formula, 3) ;
    // ========================================================================
  public:
    // ========================================================================
    /// constructor from name, expression and the tree 
    Formula 
    ( const std::string& name       , 
      const std::string& expression ,
      const TTree*       tree       ) ;
    /// constructor from name, expression and the tree 
    Formula 
    ( const std::string& name       , 
      const TCut&        expression ,
      const TTree*       tree       ) ;
    /// constructor from name, expression and the tree 
    Formula 
    ( const std::string& expression ,
      const TTree*       tree       ) ;
    /// constructor from name, expression and the tree 
    Formula 
    ( const TCut&        expression ,
      const TTree*       tree       ) ;
    /// default constructor, needed for serialisationn 
    Formula () ;
    /// virtual destructor 
    virtual ~Formula () ;
    // ========================================================================
  public:
    // ========================================================================
    /// validity check 
    bool operator!() const { return !ok() ; }
    // ========================================================================
  public:
    // ========================================================================
    /// evaluate the formula 
    double evaluate () ;       // evaluate the formula 
    /// evaluate the specified instance of the formula 
    double evaluate ( const unsigned short i ) ; // evaluate the formula 
    /// evaluate all instances of the formula 
    Int_t  evaluate ( std::vector<double>& results ) ;
    // is formula OK?
    bool   ok       () const { return this->GetNdim() ; } // is formula OK ? 
    // ========================================================================    
  };
  // ==========================================================================
  /** make Formula
   *  @param expression  (input) formula expresson  
   *  @param daat        (INPUT) input data 
   *  @param allow_empty (INPUT) return nullptr for "trivial" formula 
   *  @para, allow_null  (INPUT) return nullptr instead of exceptios 
   *  @see Ostap::trivial 
   */
  std::unique_ptr<Ostap::Formula>
  makeFormula 
  ( const std::string& expression           , 
    const TTree*       data                 , 
    const bool         allow_empty = false  , 
    const bool         allow_null  = false  ) ;
  // ==========================================================================
} //                                                     End of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_FORMULA_H
// ============================================================================
