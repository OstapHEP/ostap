// ============================================================================
#ifndef OSTYP_FORMULAE_H 
#define OSTAP_FORMULAE_H 1
// ============================================================================
// Include files 
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Formula.h"
#include "Ostap/FormulaVar.h"
// ============================================================================
namespace Ostap 
{
  // ==========================================================================
  /** @class Formular 
   *  Helper class to keep several formulae togather 
   */
  class Formulae
  {
  public:
    // ========================================================================
    typedef std::unique_ptr<Ostap::Formula> FormulaT ;
    typedef std::vector<FormulaT>           FormulaV ;
    typedef FormulaV::const_iterator        iterator ;
    // ========================================================================
  public:
    // ========================================================================
    /// creates several formulae at ne go 
    Formulae
    ( const TTree*                    tree        ,
      const std::vector<std::string>& expressions ) ;
    // ========================================================================
  public:
    // ========================================================================
    inline void evaluate 
    ( const std::size_t    index   ,
      std::vector<double>& results ) const
    { m_formulae [ index ] -> evaluate ( results ) ; }
    // ========================================================================
  public:
    // ========================================================================
    iterator begin () const { return m_formulae.begin () ; }
    iterator end   () const { return m_formulae.end   () ; }
    // ========================================================================
  private: 
    // ========================================================================
    FormulaV m_formulae {} ;
    // ========================================================================
  } ;
  // ==========================================================================
  /** @class FormulaVars  
   *  Helper class to keep several formulae togather 
   */
  class FormulaVars 
  {
  public:
    // ========================================================================
    typedef std::unique_ptr<Ostap::FormulaVar> FormulaT ;
    typedef std::vector<FormulaT>              FormulaV ;
    typedef FormulaV::const_iterator           iterator ;
    // ========================================================================
  public:
    // ========================================================================
    /// create several formulae at once     
    FormulaVars
    ( const RooArgList&               vars        ,
      const std::vector<std::string>& expressions ) ;
    // ========================================================================
    /// create several formulae at once     
    FormulaVars
    ( const RooAbsCollection*         vars        ,
      const std::vector<std::string>& expressions ) ;
    // ========================================================================
    /// create several formulae at once     
    FormulaVars
    ( const RooAbsData*               vars        ,
      const std::vector<std::string>& expressions ) ;
    // ========================================================================
  public:
    // ========================================================================
    iterator begin () const { return m_formulae.begin () ; }
    iterator end   () const { return m_formulae.end   () ; }
    // ========================================================================
  public:
    // ========================================================================
    inline double evaluate 
    ( const std::size_t index ) const
    { return m_formulae [ index ] -> getVal() ; }
    // ========================================================================
  private :
    // ========================================================================
    void make_vars
    ( const RooArgList&               vars        ,
      const std::vector<std::string>& expressions ) ;
    // ========================================================================
  private: 
    // ========================================================================
    FormulaV m_formulae {} ;
    // ========================================================================
  } ;
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
#endif // OSTAP_FORMULAE_H 1
// ============================================================================
//                                                                     The END 
// ============================================================================
