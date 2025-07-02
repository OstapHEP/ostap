// ============================================================================
// Include files 
// ============================================================================
// ROOT/RooFit 
// ============================================================================
#include "RooAbsData.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Formula.h"
#include "Ostap/FormulaVar.h"
// ============================================================================
// local 
// ============================================================================
#include "Formulae.h"
#include "status_codes.h"
// ============================================================================
// creates several formulae in one go 
// ============================================================================
Ostap::Formulae::Formulae
( const TTree*                    tree        ,
  const std::vector<std::string>& expressions )
  : m_formulae () 
{
  Ostap::Assert ( nullptr != tree                    ,
		  "Invalid TTree!"                   ,
		  "Ostap::Formulae"                  ,
		  INVALID_TREE , __FILE__ , __LINE__ ) ;
  //
  m_formulae.reserve ( expressions.size() ) ;
  for ( const auto& expr : expressions )
    {
      auto f = std::make_unique<Ostap::Formula> ( expr , tree ) ;
      Ostap::Assert ( f && f->ok ()                         ,
		      "Invalid expression:" + expr          ,
		      "Ostap::Formulae"                     ,
		      INVALID_FORMULA , __FILE__ , __LINE__ ) ;
      //
      m_formulae.push_back ( std::move ( f ) ) ;
    }
}
// ============================================================================
// create several formulae in one go     
// ============================================================================
Ostap::FormulaVars::FormulaVars
( const RooArgList&               vars        , 
  const std::vector<std::string>& expressions )
  : m_formulae()
{
  make_vars ( vars, expressions ) ;
}
// ============================================================================
// create several formulae in one go  
// ============================================================================
Ostap::FormulaVars::FormulaVars
( const RooAbsCollection*         vset        , 
  const std::vector<std::string>& expressions )
  : m_formulae()
{
  Ostap::Assert ( nullptr != vset                       ,
                  "Invalid list of dependents"          ,
                  "Ostap::FormulaVars"                  ,
                  INVALID_ARGSET, __FILE__ , __LINE__ ) ;
  //
  const RooArgList vlst { *vset } ;
  make_vars ( vlst , expressions ) ;
}
// ============================================================================
// create several formulae in one go 
// ============================================================================
Ostap::FormulaVars::FormulaVars
( const RooAbsData*               data        , 
  const std::vector<std::string>& expressions )
  : m_formulae()
{
  Ostap::Assert ( nullptr != data                    ,
                  "Invalid list of dependents"       ,
                  "Ostap::FormulaVars"               ,
                  INVALID_DATA , __FILE__ , __LINE__ ) ;
  const RooArgSet* vset = data->get () ;
  Ostap::Assert ( nullptr != vset                       ,
                  "Invalid list of dependents"          ,
                  "Ostap::FormulaVars"                  ,
                  INVALID_ARGSET, __FILE__ , __LINE__ ) ;  
  //
  const RooArgList vlst { *vset } ; 
  make_vars ( vlst , expressions ) ;
}
// ============================================================================
// create several formulae in one go     
// ============================================================================
void Ostap::FormulaVars::make_vars
( const RooArgList&               vars        ,
  const std::vector<std::string>& expressions )
{
  m_formulae.clear() ;
  m_formulae.reserve ( expressions.size() ) ;
  for ( const auto& expr : expressions )
    {
      auto f = makeFormula ( expr , vars ) ;
      Ostap::Assert ( f && f->ok ()                         ,
                      "Invalid expression:" + expr          ,
                      "Ostap::FormulaVars"                  ,
                      INVALID_FORMULA , __FILE__ , __LINE__ ) ;
      //
      m_formulae.push_back ( std::move ( f ) ) ;
    }
}
// ============================================================================
//                                                                     The END 
// ============================================================================
