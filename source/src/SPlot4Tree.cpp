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
#include "Ostap/ToStream.h"
#include "Ostap/RooFun.h"
#include "Ostap/FitResult.h"
#include "Ostap/MoreRooFit.h"
#include "Ostap/ProgressBar.h"
#include "Ostap/AddBranch.h"
#include "Ostap/TreeGetter.h"
#include "Ostap/SPlot4Tree.h"
#include "Ostap/RooFitUtils.h"
#include "Ostap/ToStream.h"
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
// get the pdf 
// ============================================================================
const RooAddPdf&
Ostap::MoreRooFit::SPlot4Tree::pdf () const
{ return static_cast<const RooAddPdf&> ( *m_fun ) ; }
// ============================================================================
// size of object: number of ocmponnent
// ============================================================================
std::size_t
Ostap::MoreRooFit::SPlot4Tree::size      () const
{ return  m_cmps->getSize () ; }
// ============================================================================
/*  @param addpdf input extended RooAddPdf 
 *  @param observabels list of observables 
 *  @param normalzation normalizartion set 
 */
// ============================================================================
Ostap::MoreRooFit::SPlot4Tree::SPlot4Tree
( const RooAddPdf&        addpdf        ,
  const RooAbsCollection& observables   ,
  const RooFitResult&     fitresult     ,   
  const RooAbsCollection* normalization )
  : RooFun ( addpdf , observables , normalization )
  , m_cmps        {} 
  , m_coefs       {} 
  , m_result      { std::make_unique<Ostap::Utils::FitResults>( fitresult ) }
{
  // 
  Ostap::Assert ( pdf().canBeExtended()                   ,
                  "PDF canot be extended"                 ,
                  "Ostap::MoreRooFit::Splot2Tree"         ,
                  INVALID_PDF , __FILE__ , __LINE__       ) ;
  /// get all components 
  m_cmps = std::make_unique<RooArgList> ( pdf().pdfList () ) ;
  /// get original fractios 
  bool recursive ;
  m_coefs = std::make_unique<RooArgList> ( fractions ( pdf () , recursive ) ) ;
  Ostap::Assert ( ::size ( *m_cmps ) == ::size ( *m_coefs )  ,
                  "Mismatch in cmponent/coefficients size" ,
                  "Ostap::MoreRoofit:Splot2Tree"           , 
                  INVALID_PDF , __FILE__ , __LINE__        ) ;
  Ostap::Assert ( !recursive                               ,
                  "Fractions cannot be recursive"          ,
                  "Ostap::MoreRoofit:Splot2Tree"           , 
                  INVALID_PDF , __FILE__ , __LINE__        ) ;
  // ==========================================================================
  // check validity of coefficiencts
  // ==========================================================================
  { // ========================================================================
    for ( auto* c : *m_coefs )
      {
        // ====================================================================
        Ostap::Assert ( nullptr != c                                ,
                        "Invalid/nullptr coefficient"               , 
                        "Ostap::MoreRoofit::SPlot4Tree"             ,
                        INVALID_ABSARG , __FILE__ , __LINE__        ) ;
        //
        const RooAbsReal*  rv = dynamic_cast<RooAbsReal*> ( c ) ;
        Ostap::Assert ( nullptr != rv , 
                        "Illegal coefficient " + Ostap::Utils::toString ( *c ) , 
                        "Ostap::MoreRoofit::SPlot4Tree"                        ,
                        INVALID_ABSARG , __FILE__ , __LINE__                   ) ;
        // ====================================================================
      } //                                                  The enf of the loop
    // ========================================================================
  } //                                                      The end of if-block
  // ==========================================================================
  // Check that all floaitng parameters are correficients
  // ==========================================================================
  auto fparams = std::make_unique<RooArgList> ( m_result->floatParsFinal() ) ;
  // ==========================================================================
  { // ========================================================================    
    for ( auto* p : *fparams )
      {
        // ====================================================================
        Ostap::Assert ( nullptr != p                                ,
                        "Invalid/nullptr parameter in RooFitResult" , 
                        "Ostap::MoreRoofit::SPlot4Tree"             ,
                        INVALID_ABSARG , __FILE__ , __LINE__        ) ;
        Ostap::Assert ( m_coefs->contains ( *p ) ,
                        "Parameteter '" + std::string ( p->GetName() ) + "' is not coefficient!" ,
                        "Ostap::MoreRoofit::SPlot4Tree"      ,
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
} // ==========================================================================
// ============================================================================
// copy constructor
// ============================================================================
Ostap::MoreRooFit::SPlot4Tree::SPlot4Tree
( const Ostap::MoreRooFit::SPlot4Tree& right )
  : RooFun   ( right ) 
  , m_cmps   { std::make_unique<RooArgList>() } 
  , m_coefs  { std::make_unique<RooArgList>() } 
  , m_result { right.m_result ? right.m_result->Clone() : nullptr }
{
  ::copy ( *right.m_cmps  , *m_cmps  ) ;
  ::copy ( *right.m_coefs , *m_coefs ) ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::MoreRooFit::SPlot4Tree::~SPlot4Tree() {}
// ============================================================================
/*  Add sPlot information to the tree 
 *  @param tree  input tree 
 *  @param splot sPlot object 
 *  @param prefix prefix for the branch names 
 *  @param suffix suffix for the branch names 
 *  @return StatusCode
 */
// ========================================================================
Ostap::StatusCode
Ostap::Trees::add_branch
( TTree*                               tree     ,
  const Ostap::MoreRooFit::SPlot4Tree& splot    ,
  const std::string&                   prefix   ,
  const std::string&                   suffix   , 
  const Ostap::Trees::DCT&             mapping  ,
  const Ostap::Utils::ProgressConf&    progress )
{
  //
  if ( !tree )                       { return INVALID_TREE ; }
  //
  // keep a local copy
  const Ostap::MoreRooFit::SPlot4Tree the_splot { splot                   } ;
  const Ostap::Trees::RooGetter       getter    ( mapping                 ,
                                                  the_splot.observables() ,
                                                  tree                    ) ;
  //
  const std::size_t N = the_splot.size()   ;
  //                  name/0      c/1   d/2    b/3 
  typedef std::tuple<std::string,double,double, double> ITEM  ;
  typedef std::vector<ITEM>                             ITEMS ; 
  ITEMS items {} ; items.reserve ( N ) ;
  // ==========================================================================
  { // ========================================================================
    // ========================================================================
    for ( const RooAbsArg* c : the_splot.coefficients () )
      {
        // =====================================================================
        Ostap::Assert ( nullptr != c                         ,
                        "Invalid coefficient"                ,
                        "Ostap::Trees::add_branch"           ,
                        INVALID_ABSARG , __FILE__ , __LINE__ ) ;
        const RooAbsArg*       aa = the_splot.fitresult().floatParsFinal().find ( c->GetName() ) ;
        if ( nullptr == aa ) { aa = the_splot.fitresult().constPars()     .find ( c->GetName() ) ; }
        Ostap::Assert ( nullptr != aa  ,
                        "Coeffcient is not found:" + Ostap::Utils::toString ( *c ) ,
                        "Ostap::Trees::add_branch"           ,                      
                        INVALID_ABSARG , __FILE__ , __LINE__ ) ;
        const RooAbsReal*      rv = dynamic_cast<const RooAbsReal*>     ( c ) ;
        const RooAbsCategory*  cv = nullptr ;
        if ( nullptr == rv ) { cv = dynamic_cast<const RooAbsCategory*> ( c ) ; }
        Ostap::Assert ( ( nullptr != rv ) || ( nullptr != cv ) ,
                        "Invalid Coeffcient:" + Ostap::Utils::toString ( *c ) , 
                        "Ostap::Trees::add_branch"           ,                      
                        INVALID_ABSARG , __FILE__ , __LINE__ ) ;
        items.push_back ( ITEM () ) ;
        if      ( rv ) { items.back() = std::make_tuple ( c->GetName() , rv->getVal (     ) , 0.0 , 0.0 ) ; } 
        else if ( cv ) { items.back() = std::make_tuple ( c->GetName() , ::getValue ( *cv ) , 0.0 , 0.0 ) ; } 
        // ====================================================================
      } // ====================================================================
    // ========================================================================
  } // ========================================================================
  // ==========================================================================
  Ostap::Assert ( items.size () == N                   , 
                  "Invalid coefficients"               ,
                  "Ostap::Trees::add_branch"           ,                      
                  INVALID_ARGSET , __FILE__ , __LINE__ ) ;
  //
  const TMatrixDSym cov { the_splot.fitresult().conditionalCovarianceMatrix ( the_splot.coefficients() ) } ;
  Ostap::Assert (  ( N == cov.GetNcols () ) &&
                   ( N == cov.GetNrows () )            ,
                  "Invalid covarinace matrix"          ,
                  "Ostap::Trees::add_branch"           ,                      
                  INVALID_ARGSET , __FILE__ , __LINE__ ) ;
  //
  // create branches 
  typedef std::vector<TBranch*> BRANCHES ; 
  BRANCHES branches {} ; branches.reserve ( N ) ;
  { // ========================================================================
    for ( auto& item  : items ) 
      {
        const std::string  cname = std::get<0> ( item  ) ;
        const std::string  bname { prefix + cname + suffix } ;
        const std::string  bspec { bname  + "/D" } ;
        TBranch* branch = tree->Branch ( bname.c_str() , &std::get<3> ( item ) , bspec.c_str () ) ; // get<3>!
        Ostap::Assert ( branch                         ,
                        "Cannot create branch " + bname             ,
                        "Ostap::Trees::add_branch"                 , 
                        CANNOT_CREATE_BRANCH , __FILE__ , __LINE__ ) ;
        //
        branches.push_back ( branch ) ;
      } // ====================================================================
    // ========================================================================
  } // ========================================================================
  // check branches 
  Ostap::Assert ( N == branches.size()                       ,
                  "Missing branch"                           , 
                  "Ostap::Trees::add_branch"                 , 
                  CANNOT_CREATE_BRANCH , __FILE__ , __LINE__ ) ;
  //
  // normalization 
  const RooArgSet* normset = the_splot.normalization() ;
  //
  // loop over the TTree entries
  const Long64_t nentries = tree->GetEntries(); 
  Ostap::Utils::ProgressBar bar ( nentries , progress  ) ;
  for ( Long64_t entry = 0 ; entry < nentries ; ++entry , ++bar )
    {
      if ( tree->GetEntry ( entry ) < 0 ) { break ; };
      //
      /// assign the observables 
      getter.assign ( the_splot.observables () , tree ) ;
      /// (1) get thc omponent values and the total 
      double total = 0 ; // sum_i ci*d_i
      for  ( unsigned int i = 0 ; i < N ; ++i )
        {
          const RooAbsReal* cmp = static_cast<RooAbsReal*> ( the_splot.components().at ( i ) ) ;
          Ostap::Assert ( nullptr != cmp                       ,
                          "Invalid fit ocmponent"              , 
                          "Ostap::Trees::add_branch"           , 
                          INVALID_ABSARG , __FILE__ , __LINE__ ) ;
          // get the component value 
          const double dval = normset ? cmp->getVal ( normset ) : cmp->getVal() ;
          total += dval * std::get<1> ( items [ i ] ) ;                  // get<1>
          std::get<2> ( items [ i ] ) = dval ;                           // get<2> 
        }
      /// (2) calculate s-weights   
      for ( unsigned int i = 0 ; i < N ; ++i )
        {
          double sum = 0 ;
          for ( unsigned int j = 0 ; j < N ; ++j )
            { sum += cov ( i , j ) * std::get<2> ( items [ j ] ) ; }    // get<2> 
          std::get<3> ( items [ i ] ) = sum / total ;                    // get<3> 
        }
      /// (3) commit branches 
      for ( auto* branch : branches ) { branch->Fill() ; }
      // ======================================================================
    } //                                   The end of th loop over Tree-entries
  // ==========================================================================
  //
  return Ostap::StatusCode::SUCCESS ;
  // ===========================================================================
}

// ============================================================================
//                                                                      The END 
// ============================================================================
