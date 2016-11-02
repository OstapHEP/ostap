// $Id$
// ============================================================================
// Include files 
// ============================================================================
// ROOT 
// ============================================================================
#include "TTree.h"
#include "TCut.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Formula.h"
// ============================================================================
/** Implementation file for class Ostap::Formula
 *  @see Ostap::Formula
 *  @see TTreeFormula
 *  @date 2013-05-06 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 */
// ============================================================================
// constructor from name, expression and the tree 
// ============================================================================
Ostap::Formula::Formula
( const std::string& name       , 
  const std::string& expression ,
  TTree*             tree       ) 
  : TTreeFormula ( name.c_str() , expression.c_str() , tree )
{}
// ============================================================================
Ostap::Formula::Formula
( const std::string& name       , 
  const TCut&        expression ,
  TTree*             tree       ) 
  : TTreeFormula ( name.c_str() , expression , tree )
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Formula::~Formula()
{
  TTree* tree = GetTree() ;
  if ( 0 != tree && this == tree->GetNotify() ) { tree -> SetNotify ( 0 ) ; }
}
// ============================================================================
// evaluate the formula 
// ============================================================================
double Ostap::Formula::evaluate () // evaluate the formula 
{ 
  GetNdata() ; 
  return EvalInstance () ; 
}
// ============================================================================
// The END 
// ============================================================================
