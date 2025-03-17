// ============================================================================
// Include files 
// ============================================================================
// ROOT
// ============================================================================
#include "TBranch.h"
#include "TTree.h"
// ============================================================================
// ROOT/RooFit 
// ============================================================================
#include "RooAbsReal.h"
#include "RooAbsData.h"
#include "RooAddPdf.h"
#include "RooRealVar.h"
#include "RooAbsRealLValue.h"
#include "RooAbsCategoryLValue.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/RooFun.h"
#include "Ostap/FitResult.h"
#include "Ostap/MoreRooFit.h"
#include "Ostap/ProgressBar.h"
#include "Ostap/AddBranch.h"
#include "Ostap/TreeGetter.h"
#include "Ostap/SPLOT.h"
#include "Ostap/RooFitUtils.h"
// ============================================================================
// local
// ============================================================================
#include "Exception.h"
#include "local_roofit.h"
#include "status_codes.h"
// ============================================================================
/** @file
 *  implementaton of various small additions to RooFit 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
 *  @date 2019-11-21
 */
// ============================================================================
/*  @param addpdf input extended RooAddPdf 
 *  @param observabels list of observables 
 *  @param normalzation normalizartion set 
 */
// ============================================================================
Ostap::Utils::SPLOT::SPLOT
( const RooAddPdf&        addpdf        ,
  const RooAbsCollection& observables   ,
  const RooFitResult&     fitresult     ,   
  const RooAbsCollection* normalization )
  : COWs ( addpdf , observables , normalization )
  , m_coefs  {}
  , m_result { std::make_unique<Ostap::Utils::FitResults>( fitresult ) }
{
  Ostap::Assert ( pdf().canBeExtended()              ,
                  "PDF must be extended!"            ,
                  "Ostap::Utils::SPLOT"              ,
                  INVALID_PDF , __FILE__ , __LINE__  ) ;
  /// get original fractions 
  bool recursive ;
  m_coefs = std::make_unique<RooArgList> ( Ostap::MoreRooFit::fractions ( pdf () , recursive ) ) ;
  Ostap::Assert ( ::size ( *m_cmps ) == ::size ( *m_coefs )  ,
		  "Mismatch in components/coefficients size" ,
		  "Ostap::Utils::COWs"                       , 
		  INVALID_PDF , __FILE__ , __LINE__          ) ;
  Ostap::Assert ( !recursive                                ,
		  "Fractions cannot be recursive"           ,
                 "Ostap::Utils::COWs"                      , 
		  INVALID_PDF , __FILE__ , __LINE__         ) ;
  // ==========================================================================
  // check validity of coefficiencts & get the notal integral 
  // ==========================================================================
  double total = 0 ;
  { // ========================================================================
    for ( auto* c : *m_coefs )
      {
        // ====================================================================
        Ostap::Assert ( nullptr != c                         ,
                        "Invalid/nullptr coefficient"        , 
                        "Ostap::Utils::SPLOT"                ,
                        INVALID_ABSARG , __FILE__ , __LINE__ ) ;
        const RooAbsRealLValue*  rv = dynamic_cast<RooAbsRealLValue*> ( c ) ;
        Ostap::Assert ( nullptr != rv , 
                        "Illegal coefficient: " + Ostap::Utils::toString ( *c ) , 
                        "Ostap::Utils::SPLOT"                                  ,
                        INVALID_ABSARG , __FILE__ , __LINE__                   ) ;
	total += rv->getVal() ;
        // ====================================================================
      } //                                                  The enf of the loop
    // ========================================================================
  } //                                                      The end of if-block
  // ==========================================================================
  // Check that all floaitng parameters are coefficients  
  // ==========================================================================
  auto fpars = std::make_unique<RooArgList> ( m_result->floatParsFinal() ) ;
  // ==========================================================================
  { // ========================================================================    
    for ( auto* p : *fpars )
      {
        // ====================================================================
        Ostap::Assert ( nullptr != p                                ,
                        "Invalid/nullptr parameter in RooFitResult" , 
                        "Ostap::Utils::SPLOT"                       ,
                        INVALID_ABSARG , __FILE__ , __LINE__        ) ;
        Ostap::Assert ( m_coefs->contains ( *p ) ,
                        "Parameter `" + std::string ( p->GetName() ) + "' is not coefficient!" ,
                        "Ostap::Utils::SPLOT"                ,
                        INVALID_ABSARG , __FILE__ , __LINE__ ) ;
        // ====================================================================
      } // The end of the loop
    // ========================================================================
  } //                                                      The enf of if-block 
  // ==========================================================================
  // load parameters from fit result
  set_parameters ( m_result->constPars      () ) ;
  set_parameters ( m_result->floatParsFinal () ) ;
  // ==========================================================================      
  // get the matrix (and extend it if needed) 
  // ==========================================================================
  const TMatrixDSym& cm = m_result->covarianceMatrix () ; 
  // extended covariance matrix 
  const std::size_t N { size() } ;
  TMatrixDSym cov { static_cast<Int_t> ( N ) } ;
  if ( N == cm.GetNrows() ) { cov = cm ; }
  else
    {
      // fill the matrix with the content & zeroes 
      const RooArgList& floatParsFinal = m_result->floatParsFinal() ;
      const RooArgList& coefficients   = this->coefficients() ;
      for ( const RooAbsArg* ci : coefficients  )
	{
	  const int i   = floatParsFinal.index ( ci->GetName () ) ;
	  const int row = coefficients  .index ( ci ) ; 
	  for ( const RooAbsArg* cj : coefficients)
	    {
	      const int j   = floatParsFinal.index ( cj->GetName()  ) ;
	      const int col = coefficients  .index ( cj ) ;
	      //
	      const double cij = ( i < 0 || j < 0 ) ? 0.0 : cm ( i , j ) ;
	      //
	      cov ( row , col ) = cij  ;
	      if  ( row != col ) { cov ( col , row ) = cij ; }
	    }
        }
    }
  //
  Ostap::Assert ( cov.IsValid () &&
		  ( N == cov.GetNcols () ) &&
		  ( N == cov.GetNrows () )               ,
		  "Invalid covariance matrix"            ,
		  "Ostap::Utils::SPLOT"                  ,                      
		  INVALID_ARGSET , __FILE__ , __LINE__ ) ;
  // scale it! 
  cov *= 1.0L / total ; 
  // copy matrux to COWs :
  m_A.ResizeTo ( cov ) ;  
  m_A = cov ;
} // ==========================================================================
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Utils::SPLOT::SPLOT
( const Ostap::Utils::SPLOT& right )
  : COWs ( right ) 
  , m_coefs  {}
  , m_result { right.m_result ? right.m_result->Clone() : nullptr }
{
  /// get original fractios 
  bool recursive ;
  m_coefs = std::make_unique<RooArgList> ( Ostap::MoreRooFit::fractions ( pdf () , recursive ) ) ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Utils::SPLOT::~SPLOT() {}
// ============================================================================
// clone 
// ============================================================================
Ostap::Utils::SPLOT* 
Ostap::Utils::SPLOT::clone () const 
{ return new Ostap::Utils::SPLOT  ( *this ) ; } 

// ==========================================================================
/*  Add sPlot information to the tree 
 *  @param tree  input tree 
 *  @param splot sPlot object 
 *  @param prefix prefix for the branch names 
 *  @param suffix suffix for the branch names 
 *  @return StatusCode
 */
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_branch
( TTree*                            tree     ,
  const Ostap::Utils::SPLOT&        splot    ,
  const std::string&                prefix   ,
  const std::string&                suffix   , 
  const Ostap::Trees::DCT&          mapping  ,
  const Ostap::Utils::ProgressConf& progress )
{
  //
  if ( !tree ) { return INVALID_TREE ; }
  //
  const std::size_t  N  { splot.size() } ; 
  // ==========================================================================
  std::vector<std::string> names ; names.reserve ( N ) ; 
  // ==========================================================================
  for ( const RooAbsArg* c : splot.coefficients () )
    {
      Ostap::Assert ( nullptr != c                         ,
                      "Invalid coefficient"                ,
                      "Ostap::Trees::add_branch"           ,
                      INVALID_ABSARG , __FILE__ , __LINE__ ) ;
      const RooAbsArg*       aa = splot.fitresult().floatParsFinal().find ( c->GetName() ) ;
      if ( nullptr == aa ) { aa = splot.fitresult().constPars()     .find ( c->GetName() ) ; }
      Ostap::Assert ( nullptr != aa  ,
		      "Coefficient is not found:" + Ostap::Utils::toString ( *c ) ,
		      "Ostap::Trees::add_branch"           ,                      
		      INVALID_ABSARG , __FILE__ , __LINE__ ) ;
      names.push_back (  prefix + aa->GetName() + suffix ) ;
    }
  // =========================================================================
  return add_branch ( tree     , 
                      splot    , 
                      names    ,   
                      mapping  ,  
                      progress ) ;  
}
// ============================================================================
//                                                                      The END 
// ============================================================================
