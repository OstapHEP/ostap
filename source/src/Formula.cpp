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
#include "Ostap/Names.h"
#include "Ostap/StatusCode.h"
#include "Ostap/Formula.h"
// ============================================================================
// Local
// ============================================================================
#include "status_codes.h"
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
  std::string formula_name
  ( const std::string& prefix     ,
    const std::string& expression ,
    const TTree*       tree       ) 
  { return Ostap::tmp_name ( prefix , Ostap::strip ( expression ) , tree , true ) ; }
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
: TTreeFormula ( name                        . c_str() ,
		 Ostap::strip ( expression ) . c_str() ,
		 const_cast<TTree*> ( tree )           ) 
{}
// ============================================================================
Ostap::Formula::Formula
( const std::string& name   , 
  const TCut&        cut   ,
  const TTree*       tree  ) 
  : TTreeFormula ( name                             . c_str() ,
		   Ostap::strip ( cut.GetTitle () ) . c_str() ,
                   const_cast<TTree*> ( tree ) ) 
{}
// ============================================================================
Ostap::Formula::Formula
( const std::string& expression ,
  const TTree*       tree       ) 
  : Formula ( formula_name ( "formula_" , expression , tree ) ,
	      Ostap::strip ( expression ) , tree )
{}
// ============================================================================
Ostap::Formula::Formula
( const TCut&   cut  ,
  const TTree*  tree ) 
  : Formula ( formula_name ( "formula_" , cut.GetTitle() , tree ) , cut , tree ) 
{}
// ============================================================================
Ostap::Formula::Formula() 
  : Formula ( std::string ( "1" ) , nullptr )
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
                  "evaluate: scalar call for vector [ GetNdata()!=1 ]  function" , 
                  "Ostap::Formula"     , 
		  INVALID_FORMULA_CALL , __FILE__ , __LINE__ ) ;
  return EvalInstance () ; 
}
// ============================================================================
// evaluate the given instance of the formula 
// ============================================================================
double Ostap::Formula::evaluate ( const unsigned short i ) // evaluate the formula 
{ 
  const Int_t d = GetNdata() ; 
  Ostap::Assert ( i < d ,
                  "evaluate: invalid instance counter" , 
                  "Ostap::Formula"           ,
		  INVALID_FORMULA_CALL , __FILE__ , __LINE__ ) ;
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
/*  make Formula
 *  @param expression  (input) formula expresson  
 *  @param daat        (INPUT) input data 
 *  @param allow_empty (INPUT) return nullptr for "trivial" formula 
 *  @para, allow_null  (INPUT) return nullptr instead of exceptios 
 *  @see Ostap::trivial 
 */
// ============================================================================
std::unique_ptr<Ostap::Formula>
Ostap::makeFormula 
( const std::string& expression  , 
  const TTree*       data        , 
  const bool         allow_empty , 
  const bool         allow_null  )
{
  if ( allow_empty && Ostap::trivial ( expression ) ) { return nullptr ; }  // RETURN!
  // 
  auto result = std::make_unique<Ostap::Formula> ( expression , data  ) ;
  if ( allow_null && ( !result || !result->ok() ) ) { return nullptr ; } 
  //
  Ostap::Assert ( result && result->ok()                   , 
		  "Invalid formula:\'" + expression + "\'" , 
		  "Ostap::Formula::makeFormula"            ,
		  INVALID_FORMULA , __FILE__ , __LINE__    ) ;
  //
  return result ;
}
// ============================================================================
//                                                                      The END 
// ============================================================================
