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
// Local
// ============================================================================
#include "Exception.h"
// ============================================================================
/** Implementation file for class Ostap::Formula
 *  @see Ostap::Formula
 *  @see TTreeFormula
 *  @date 2013-05-06 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 */
// ============================================================================
ClassImp(Ostap::Formula)
// ============================================================================
// constructor from name, expression and the tree 
// ============================================================================
Ostap::Formula::Formula
( const std::string& name       , 
  const std::string& expression ,
  TTree*             tree       ) 
: TTreeFormula ( name.c_str() , expression.c_str() , tree )
{
  
}
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
// The END 
// ============================================================================
