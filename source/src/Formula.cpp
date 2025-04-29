// ============================================================================
// Include files 
// ============================================================================
// STD&STL
// ============================================================================
#include <string>
// ============================================================================
// ROOT 
// ============================================================================
#include "RVersion.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TCut.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Formula.h"
// ============================================================================
// Local
// ============================================================================
#include "Exception.h"
#include "local_utils.h"
// ============================================================================
/** Implementation file for class Ostap::Formula
 *  @see Ostap::Formula
 *  @see TTreeFormula
 *  @date 2013-05-06 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  std::string formula_name ( const std::string& prefix     ,
                             const std::string& expression ,
                             const TTree*       tree       ) 
  { return Ostap::tmp_name ( prefix , expression , tree , true ) ; }
  // ==========================================================================
} //                                             The end of anonymous namespace 
// ============================================================================
#if ROOT_VERSION_CODE < ROOT_VERSION(6,36,0)
// ============================================================================
ClassImp(Ostap::Formula)
// ============================================================================
#endif 
// ============================================================================
// constructor from name, expression and the tree 
// ============================================================================
Ostap::Formula::Formula
( const std::string& name       , 
  const std::string& expression ,
  const TTree*       tree       ) 
: TTreeFormula ( name.c_str()                ,
                 expression.c_str()          ,
                 const_cast<TTree*> ( tree ) ) 
{}
// ============================================================================
Ostap::Formula::Formula
( const std::string& name       , 
  const TCut&        expression ,
  const TTree*       tree       ) 
  : TTreeFormula ( name.c_str()                ,
                   expression                  , 
                   const_cast<TTree*> ( tree ) ) 
{}
// ============================================================================
Ostap::Formula::Formula
( const std::string& expression ,
  const TTree*       tree       ) 
  : Formula ( formula_name ( "formula_" , expression , tree ) , expression , tree )
{}
// ============================================================================
Ostap::Formula::Formula
( const TCut&        expression ,
  const TTree*       tree       ) 
  : Formula ( formula_name ( "formula_" , expression.GetName() , tree ) , expression , tree ) 
{}
// ============================================================================
Ostap::Formula::Formula() 
  : Formula ( std::string ( "1" ) , nullptr )
{}
// ============================================================================

// ============================================================================
// destructor 
// ============================================================================
Ostap::Formula::~Formula()
{
  TTree* tree = GetTree() ;
  if ( nullptr != tree && this == tree->GetNotify() ) { tree -> SetNotify ( 0 ) ; }
}
// ============================================================================
// evaluate the formula 
// ============================================================================
double Ostap::Formula::evaluate () // evaluate the formula 
{ 
  const Int_t d = GetNdata() ; 
  Ostap::Assert ( 1 == d , 
                  "evaluate: scalar call for GetNdata()!=1 function" , 
                  "Ostap::Formula"           , 
                  Ostap::StatusCode::FAILURE ) ;
  return EvalInstance () ; 
}
// ============================================================================
// evaluate the given instance of the formula 
// ============================================================================
double Ostap::Formula::evaluate ( const unsigned short i ) // evaluate the formula 
{ 
  const Int_t d = GetNdata() ; 
  Ostap::Assert ( i  < d ,
                  "evaluate: invalid instance counter" , 
                  "Ostap::Formula"           , 
                  Ostap::StatusCode::FAILURE ) ;
  return EvalInstance ( i ) ; 
}
// ============================================================================
// evaluate all the instances of the formula 
// ============================================================================
Int_t Ostap::Formula::evaluate ( std::vector<double>& results ) 
{ 
  const Int_t d = GetNdata() ; 
  results.resize ( d ) ;
  for ( Int_t i = 0 ; i < d ; ++i ) { results [ i ] = EvalInstance ( i ) ; }
  return d ;  
}
// ============================================================================
//                                                                      The END 
// ============================================================================
