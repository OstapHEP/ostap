// ============================================================================
// Include files 
// ============================================================================
// ROOT
// ============================================================================
#include "TDecompChol.h"
#include "TTree.h"
#include "TBranch.h"
// ============================================================================
// ROOT/RooFit 
// ============================================================================
#include "RooAbsReal.h"
#include "RooAddPdf.h"
#include "RooRealVar.h"
#include "RooAbsRealLValue.h"
#include "RooAbsCategoryLValue.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/AddBranch.h"
#include "Ostap/COWs.h"
#include "Ostap/MoreRooFit.h"
#include "Ostap/TreeGetter.h"
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
// get the pdf 
// ============================================================================
const RooAddPdf&
Ostap::MoreRooFit::COWs::pdf () const
{ return static_cast<const RooAddPdf&> ( *m_fun ) ; }
// ============================================================================
// size of object: number of ocmponnent
// ============================================================================
std::size_t
Ostap::MoreRooFit::COWs::size      () const
{ return  m_cmps->getSize () ; }
// ============================================================================
/* @param addpdf input extended RooAddPdf 
 *  @param datax   input data
 *  @param normalzation normalisation set 
 */
// ============================================================================
Ostap::MoreRooFit::COWs::COWs
( const RooAddPdf&        addpdf        ,
  const RooAbsData&       data          ,
  const RooAbsCollection* normalization )
  : RooFun ( addpdf , data , normalization )
  , m_cmps        {} 
  , m_coefs       {} 
{
  /// get all components 
  m_cmps = std::make_unique<RooArgList> ( pdf().pdfList () ) ;
  /// get original fractios 
  bool recursive ;
  m_coefs = std::make_unique<RooArgList> ( fractions ( pdf () , recursive ) ) ;
  Ostap::Assert ( ::size ( *m_cmps ) == ::size ( *m_coefs ) ,
                  "Mismatch in component/coefficients size" ,
                  "Ostap::MoreRoofit:COWs"                  , 
                  INVALID_PDF , __FILE__ , __LINE__         ) ;
  Ostap::Assert ( !recursive                                ,
                  "Fractions cannot be recursive"           ,
                  "Ostap::MoreRoofit::COWs"                 , 
                  INVALID_PDF , __FILE__ , __LINE__         ) ;
  ///
  // ==========================================================================
  // check validity of coefficiencts
  // ==========================================================================
  { // ========================================================================
    for ( auto* c : *m_coefs )
      {
        // ====================================================================
        Ostap::Assert ( nullptr != c                         ,
                        "Invalid/nullptr coefficient"        , 
                        "Ostap::MoreRooFit::COWs"            ,
                        INVALID_ABSARG , __FILE__ , __LINE__ ) ;
        //
        const RooAbsReal*  rv = dynamic_cast<RooAbsReal*> ( c ) ;
        Ostap::Assert ( nullptr != rv , 
                        "Illegal coefficient " + Ostap::Utils::toString ( *c ) , 
                        "Ostap::MoreRooFit::COWs"                              ,
                        INVALID_ABSARG , __FILE__ , __LINE__                   ) ;
        // ====================================================================
      } //                                                  The enf of the loop
    // ========================================================================
  } //                                                      The end of if-block
  // ==========================================================================
  
  /// get obs 
  std::unique_ptr<RooArgSet> obs  { pdf().getObservables ( data ) } ;
  const RooArgSet*           norm { this->normalization () } ;
  const RooAbsPdf*           fun  { &pdf() } ;
  const std::size_t          N    { size() } ;
  //
  
  // value of N-components 
  std::vector<double>  cmp_val ( N , 0.0 ) ;
  // the matrix W
  Int_t NN = N ;
  TMatrixDSym W { NN } ;
  
  // loop over entries
  const bool        weighted = data.isWeighted () ;
  const std::size_t nEntries = data.numEntries () ;
  for ( std::size_t entry = 0 ; entry < nEntries ; ++entry )
    {
      //
      const RooArgSet* item = data.get  ( entry ) ;
      Ostap::Assert ( nullptr != item                       ,
		      "Invalid/null  entry in datatset!"    ,
		      "Ostap::MoreRooFit::COWs"             ,
		      INVALID_ABSDATA , __FILE__ , __LINE__ ) ;
      //
      const double weight = weighted ? data.weight() : 1.0 ;
      if ( !weight ) { continue ; } 
      //
      ::assign ( *obs , *item ) ;
      // evalute pdf 
      const double pdf_val = norm ? fun->getVal ( norm ) : fun->getVal () ;
      //
      if ( !pdf_val ) { continue ; } // skip it...  Is it correct? 
      //      
      // evaluate all components
      std::size_t i = 0 ;
      for ( const RooAbsArg* c : *m_cmps )
	{
	  const RooAbsReal* r = static_cast<const RooAbsReal*> ( c ) ;
	  cmp_val [ i ] = norm ? r->getVal ( norm ) : r->getVal () ;	  
	  ++i ;
	}
      //
      const double factor = weight / ( pdf_val * pdf_val ) ;
      for ( std::size_t k = 0 ; k < N ; ++k )
	{
	  const double kval = cmp_val [ k ] ;
	  if ( !kval ) { continue ; }                    // CONTINUE 
	  W ( k , k ) += factor * kval * kval ;
	  for ( std::size_t l = k + 1  ; l < N ; ++l  )
	    {
	      const double wkl = factor * kval * cmp_val [ l ] ;
	      W ( k , l ) += 0.5 * wkl ;
	      W ( l , k ) += 0.5 * wkl ;							      
	    }
	}
    }
  //
  W.Print ( "vvv" ) ;
  // invert the matrix
  double det = 0 ;
  W.Invert ( &det ) ;
  Ostap::Assert ( W.IsValid () ,
		  "Matrix W cannot be inverted!" ,
		  "Ostap::MoreRooFot::COWs" ,
		  INVALID_ABSDATA , __FILE__ , __LINE__ ) ;
  /// finally store matrix A
  m_A.ResizeTo ( W ) ;
  m_A = W ;
  m_A.Print ( "vvv" ) ;    
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::MoreRooFit::COWs::COWs
( const Ostap::MoreRooFit::COWs& right )
  : RooFun ( right ) 
  , m_cmps   { std::make_unique<RooArgList>() } 
  , m_coefs  { std::make_unique<RooArgList>() } 
  , m_A ( right.m_A ) 
{
  ::copy ( *right.m_cmps  , *m_cmps  ) ;
  ::copy ( *right.m_coefs , *m_coefs ) ;
}
// ============================================================================
/* "Recovery constructor
 *  @param addpdf input extended RooAddPdf 
 *  @param observables  observables set  
 *  @param normalzation normalisation set 
 *  @param cows         the symmetric matrix of weights 
 */
// ============================================================================
Ostap::MoreRooFit::COWs::COWs
( const RooAddPdf&        addpdf        ,
  const RooAbsCollection& observables   , 
  const TMatrixDSym&      cows          , 
  const RooAbsCollection* normalization )
  : COWs ( addpdf , observables , normalization , cows )
{}
// ============================================================================
/* "Recovery constructor
 *  @param addpdf input extended RooAddPdf 
 *  @param observables  observables set  
 *  @param normalzation normalisation set 
 *  @param cows         the symmetric matrix of weights 
 */
// ============================================================================
Ostap::MoreRooFit::COWs::COWs
( const RooAddPdf&        addpdf        ,
  const RooAbsCollection& observables   , 
  const RooAbsCollection* normalization ,
  const TMatrixDSym&      cows          )
  : RooFun ( addpdf , observables , normalization )
  , m_cmps        {} 
  , m_coefs       {}
  , m_A           { cows } 
{
  /// get all components 
  m_cmps = std::make_unique<RooArgList> ( pdf().pdfList () ) ;
  /// get original fractios 
  bool recursive ;
  m_coefs = std::make_unique<RooArgList> ( fractions ( pdf () , recursive ) ) ;
  Ostap::Assert ( ::size ( *m_cmps ) == ::size ( *m_coefs ) ,
                  "Mismatch in component/coefficients size" ,
                  "Ostap::MoreRoofit:COWs"                  , 
                  INVALID_PDF , __FILE__ , __LINE__         ) ;
  Ostap::Assert ( !recursive                                ,
                  "Fractions cannot be recursive"           ,
                  "Ostap::MoreRoofit::COWs"                 , 
                  INVALID_PDF , __FILE__ , __LINE__         ) ;
  ///
  // ==========================================================================
  // check validity of coefficiencts
  // ==========================================================================
  { // ========================================================================
    for ( auto* c : *m_coefs )
      {
        // ====================================================================
        Ostap::Assert ( nullptr != c                         ,
                        "Invalid/nullptr coefficient"        , 
                        "Ostap::MoreRooFit::COWs"            ,
                        INVALID_ABSARG , __FILE__ , __LINE__ ) ;
        //
        const RooAbsReal*  rv = dynamic_cast<RooAbsReal*> ( c ) ;
        Ostap::Assert ( nullptr != rv , 
                        "Illegal coefficient " + Ostap::Utils::toString ( *c ) , 
                        "Ostap::MoreRooFit::COWs"                              ,
                        INVALID_ABSARG , __FILE__ , __LINE__                   ) ;
        // ====================================================================
      } //                                                  The enf of the loop
    // ========================================================================
  } //                                                      The end of if-block
  // ==========================================================================
  // check the matrix
  Ostap::Assert ( m_A.IsValid()                          ,
		  "Input matrix is invalid"              ,
		  "Ostap::MoreRooFit::COWs"              ,
		  INVALID_TMATRIX  , __FILE__ , __LINE__ ) ;
  // check the matrix
  Ostap::Assert ( m_A.GetNrows() == size () &&
		  m_A.GetNcols() == size ()              , 
		  "Wrong structure of input matrix"      ,


		  "Ostap::MoreRooFit::COWs"              ,
		  INVALID_TMATRIX  , __FILE__ , __LINE__ ) ;
  // check that it is positive definite
  TDecompChol check { m_A } ;
  Ostap::Assert ( check.GetU().IsValid() ,
		  "Input matrix is not posiitve definite" ,
		  "Ostap::MoreRooFit::COWs"               ,
		  INVALID_TMATRIX  , __FILE__ , __LINE__  ) ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::MoreRooFit::COWs::~COWs() {} ;




// ============================================================================
/*  Add sPlot/COWs  information to the tree 
 *  @param tree  input tree 
 *  @param cows  COWs object 
 *  @param prefix prefix for the branch names 
 *  @param suffix suffix for the branch names 
 *  @return StatusCode
 */
// ============================================================================
Ostap::StatusCode
Ostap::Trees::add_branch
( TTree*                            tree     ,
  const Ostap::MoreRooFit::COWs&    cows     ,
  const std::string&                prefix   ,
  const std::string&                suffix   , 
  const Ostap::Trees::DCT&          mapping  , 
  const Ostap::Utils::ProgressConf& progress )
{ 
  if ( !tree )                       { return INVALID_TREE ; }
  //
  // keep a local copy
  const Ostap::MoreRooFit::COWs the_cows { cows                  } ;
  const Ostap::Trees::RooGetter getter   ( mapping                ,
					   the_cows.observables() ,
					   tree                   ) ;
  //
  const std::size_t N = the_cows.size()   ;
  
  //                  name/0      c/1    d/2    b/3 
  typedef std::tuple<std::string,double,double, double> ITEM  ;
  typedef std::vector<ITEM>                             ITEMS ; 
  ITEMS items {} ; items.reserve ( N ) ;
  // ==========================================================================
  { // ========================================================================
    // ========================================================================
    for ( const RooAbsArg* c : the_cows.coefficients () )
      {
        // =====================================================================
        Ostap::Assert ( nullptr != c                         ,
                        "Invalid coefficient"                ,
                        "Ostap::Trees::add_branch"           ,
                        INVALID_ABSARG , __FILE__ , __LINE__ ) ;
	//
        const RooAbsReal*      rv = dynamic_cast<const RooAbsReal*>     ( c ) ;
        const RooAbsCategory*  cv = nullptr ;
        if ( nullptr == rv ) { cv = dynamic_cast<const RooAbsCategory*> ( c ) ; }
        Ostap::Assert ( ( nullptr != rv ) || ( nullptr != cv ) ,
                        "Invalid coeffcient:" + Ostap::Utils::toString ( *c ) , 
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
                  "Invalid coefficients!"              ,
                  "Ostap::Trees::add_branch"           ,                      
                  INVALID_ARGSET , __FILE__ , __LINE__ ) ;
  // ========================================================================== 
  // create branches
  // ==========================================================================
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

  
  return 1 ; 
}


// ============================================================================
//                                                                      The END 
// ============================================================================
