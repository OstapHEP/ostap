// ============================================================================
// Include files 
// ============================================================================
// ROOT
// ============================================================================
#include "TTree.h"
#include "TBranch.h"
// ============================================================================
// ROOT/RooFit 
// ============================================================================
#include "RooAbsReal.h"
#include "RooArgSet.h"
#include "RooAbsData.h"
#include "RooAbsRealLValue.h"
#include "RooAbsCategoryLValue.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatusCode.h"
#include "Ostap/ProgressBar.h"
#include "Ostap/AddBranch.h"
#include "Ostap/Notifier.h"
#include "Ostap/TreeGetter.h"
#include "Ostap/RooFun.h"
#include "Ostap/RooFitUtils.h"
// ============================================================================
// local
// ============================================================================
#include "local_math.h"
#include "local_roofit.h"
#include "status_codes.h"
// ============================================================================
/** @file
 *  implementaton of various small additions to RooFit 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
 *  @date 2019-11-21
 */
// ============================================================================
// destructor
// ============================================================================
Ostap::Utils::RooFun::~RooFun() {}
// ============================================================================
// evaluate it!
// ============================================================================
double
Ostap::Utils::RooFun::evaluate () const
{ return m_normset ? m_fun->getVal ( *m_normset ) : m_fun->getVal() ; }
// ============================================================================
/// assign observables 
// ============================================================================
void
Ostap::Utils::RooFun::set_observables
( const RooAbsCollection& obs  ) const
{ m_observables -> assignValueOnly ( obs  ) ; }
// ============================================================================
void
Ostap::Utils::RooFun::set_parameters
( const RooAbsCollection& pars ) const
{ m_parameters  -> assignValueOnly ( pars ) ; }
// ============================================================================
/*  @param fun the function 
 *  @param observabels list of observables 
 *  @param normalzation normalization set 
 */
// ============================================================================
Ostap::Utils::RooFun::RooFun 
( const RooAbsReal&       fun           ,
  const RooAbsCollection& observables   ,
  const RooAbsCollection* normalization )
  : m_fun         { static_cast<RooAbsReal*> ( fun.clone () ) }
  , m_observables () 
  , m_parameters  () 
  , m_normset     () 
{
  this->Init ( observables , normalization ) ;
}
// ============================================================================
/*  @param fun the function 
 *  @param data input data 
 *  @param normalzation normalization set 
 */
// ============================================================================
Ostap::Utils::RooFun::RooFun 
( const RooAbsReal&       fun           ,
  const RooAbsData&       data          ,
  const RooAbsCollection* normalization )
  : m_fun         { static_cast<RooAbsReal*> ( fun.clone () ) }
  , m_observables () 
  , m_parameters  () 
  , m_normset     () 
{
  // ==========================================================================
  const std::unique_ptr<RooArgSet> observables { m_fun->getObservables ( data ) } ; 
  this->Init ( *observables , normalization ) ;
  // ==========================================================================
}
// ============================================================================
// perform initialization 	  
// ============================================================================
void Ostap::Utils::RooFun::Init
( const RooAbsCollection& observables   ,
  const RooAbsCollection* normalization )
{  
  /// the get observables
  // ==========================================================================
  RooArgSet obsset { observables } ;
  // ==========================================================================
  // get obserevables 
  m_observables = std::unique_ptr<RooArgSet> ( m_fun->getObservables ( obsset ) ) ;
  // ==========================================================================      
  Ostap::Assert ( ::size (    observables ) == ::size ( obsset ) &&
                  ::size ( *m_observables ) == ::size ( obsset )  ,
                  "Invalid input observables"                     ,
                  "Ostap::Utils::Roofun"                          ,
                  INVALID_OBSERVABLES , __FILE__ , __LINE__       ) ;                   
  // ==========================================================================
  { // ========================================================================
    for ( auto* o : *m_observables )
      {
        // ====================================================================
        Ostap::Assert ( nullptr != o                             ,
                        "Invalid/nullptr observable"             , 
                        "Ostap::MoreRoofit::RooFun"              ,
                        INVALID_OBSERVABLE , __FILE__ , __LINE__ ) ;
	//
	const RooAbsRealLValue*     rv = dynamic_cast<RooAbsRealLValue*> ( o ) ;
	const RooAbsCategoryLValue* cv = nullptr ;
	if (  nullptr == rv ) { cv = dynamic_cast<RooAbsCategoryLValue*> ( o ) ; }
	Ostap::Assert ( ( nullptr != rv ) || ( nullptr != cv ) , 
			            "Illegal observable " + Ostap::Utils::toString ( *o ) , 
			            "Ostap::Utils::RooFun"                           ,
			            INVALID_OBSERVABLE , __FILE__ , __LINE__              ) ;
	      // ====================================================================
      } //                                                  The end of the loop
    // ========================================================================
  } //                                                       The en of if-block
  // ==========================================================================
  // get parameters
  m_parameters = std::unique_ptr<RooArgSet> ( m_fun->getParameters ( *m_observables ) ) ;
  // ==========================================================================
  // normalization
  // ==========================================================================
  if ( normalization )
    {
      // ======================================================================
      m_normset = std::make_unique<RooArgSet> ( *normalization ) ;
      // ======================================================================
    }
  // ==========================================================================
}
// ============================================================================
// Copy constuctor
// ============================================================================
Ostap::Utils::RooFun::RooFun
  ( const Ostap::Utils::RooFun& right )
  : RooFun ( *right.m_fun , *right.m_observables , right.m_normset.get() )
{}
// ============================================================================
// clone
// ============================================================================
Ostap::Utils::RooFun* 
Ostap::Utils::RooFun::clone () const 
{ return new RooFun ( *this ) ; }


// ============================================================================
// RooFun -> TTree 
// ============================================================================
/*  add new branch to the tree from RooFit function
 *  @param tree input tree 
 *  @param bname branch name 
 *  @param fun   the function 
 *  @param mapping : observables -> brnaches 
 */
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_branch
( TTree*                            tree      ,
  const std::string&                bname     , 
  const Ostap::Utils::RooFun&       fun       ,
  const DCT&                        mapping   ,
  const Ostap::Utils::ProgressConf& progress  ) 
{
  //
  if ( !tree )                       { return INVALID_TREE ; }
  //
  // Pair of helper objects  
  std::unique_ptr<Ostap::Utils::RooFun> the_fun { fun.clone() } ;
  Ostap::Trees::RooGetter               getter  ( mapping , the_fun->observables() , tree ) ;
  //
  // create the branch 
  Double_t bvalue = 0  ;
  TBranch* branch = tree->Branch ( bname.c_str() , &bvalue , ( bname + "/D" ).c_str() );
  if ( !branch ) { return CANNOT_CREATE_BRANCH ; }
  //
  /// create the notifies 
  Ostap::Utils::Notifier notifier { tree , &getter } ; 
  // due to  some strange reasons we need to invoke the Notifier explicitely.
  notifier.Notify() ;
  //
  // loop over the TTree entries
  const Long64_t nentries = tree->GetEntries(); 
  Ostap::Utils::ProgressBar bar ( nentries , progress  ) ;
  for ( Long64_t i = 0 ; i < nentries ; ++i , ++bar )
    {
      if ( tree->GetEntry ( i ) < 0 ) { break ; };
      //
      /// assign the observables 
      getter.assign ( the_fun->observables () , tree ) ;
      /// assign the branch  
      bvalue = the_fun->evaluate () ;
      //
      branch->Fill () ;
    }
  //
  return Ostap::StatusCode::SUCCESS ;  
}
// ============================================================================
/** add new branch to the tree from RooFir function
 *  @param tree input tree 
 *  @param bname branch name 
 *  @param fun   the function 
 *  @param observables   function observables 
 *  @oaram normalization normalization set 
 *  @param mapping : observables -> brnaches 
 */
// ========================================================================
Ostap::StatusCode
Ostap::Trees::add_branch
( TTree*                            tree          ,
  const std::string&                bname         , 
  const RooAbsReal&                 fun           ,
  const RooAbsCollection&           observables   ,
  const RooAbsCollection*           normalization , 
  const Ostap::Trees::DCT&          mapping       , 
  const Ostap::Utils::ProgressConf& progress      ) 
{
  //
  if ( !tree )                       { return INVALID_TREE ; }
  //
  const std::unique_ptr<Ostap::Utils::RooFun> 
  the_fun { std::make_unique<Ostap::Utils::RooFun>( fun , observables , normalization ) } ;
  return add_branch ( tree     ,
                      bname    , 
                      *the_fun ,
                      mapping  ,
                      progress ) ;
}

// ============================================================================
//                                                                      The END 
// ============================================================================
  
