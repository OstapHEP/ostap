// ============================================================================
// Include files 
// ============================================================================
// ROOT
// ============================================================================
#include "TDecompChol.h"
#include "TTree.h"
#include "TBranch.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"
// ============================================================================
// ROOT/RooFit 
// ============================================================================
#include "RooAbsReal.h"
#include "RooAddPdf.h"
#include "RooRealVar.h"
#include "RooAbsRealLValue.h"
#include "RooAbsCategoryLValue.h"
#include "RooAbsData.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/ProgressConf.h" 
#include "Ostap/ProgressBar.h"
#include "Ostap/AddBranch.h" 
#include "Ostap/COWs.h"
#include "Ostap/MoreRooFit.h"
#include "Ostap/TreeGetter.h"
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
Ostap::Utils::COWs::pdf () const
{ return static_cast<const RooAddPdf&> ( *m_fun ) ; }
// ============================================================================
/// size of object: number of componnents
// ============================================================================
std::size_t
Ostap::Utils::COWs::size () const { return m_cmps->getSize () ; }
// ============================================================================
/* @param addpdf input extended RooAddPdf 
 *  @param datax   input data 
 *  @param progress   configraton of progreess bar 
 */
// ============================================================================
Ostap::Utils::COWs::COWs
( const RooAddPdf&                  addpdf        ,
  const RooAbsData&                 data          , 
  const Ostap::Utils::ProgressConf& progress      ) 
  : COWs ( addpdf , data , nullptr , progress ) 
{}
// =============================================================================
/* @param addpdf input extended RooAddPdf 
 *  @param datax   input data
 *  @param normalzation normalisation set 
 *  @param progress   configraton of progreess bar 
 */
// ============================================================================
Ostap::Utils::COWs::COWs
( const RooAddPdf&                  addpdf        ,
  const RooAbsData&                 data          ,
  const RooAbsCollection*           normalization , 
  const Ostap::Utils::ProgressConf& progress      ) 
  : RooFun ( addpdf , data , normalization )
  , m_cmps        {} 
{
  /// get all components 
  m_cmps = std::make_unique<RooArgList> ( pdf().pdfList () ) ;
  /// check the components
  for ( const RooAbsArg* cmp : *m_cmps )
    {
      Ostap::Assert ( nullptr != cmp ,
                      "Invaild component!" , 
                      "Ostap::Utils::COWs" ,
                      INVALID_ABSARG  , __FILE__ , __LINE__ ) ;
      const RooAbsPdf*  a = dynamic_cast<const RooAbsPdf*> ( cmp ) ; 
      Ostap::Assert ( nullptr != a                         , 
                      "Invalid component!"                 , 
                      "Ostap::Utils::COWs"                 , 
                      INVALID_ABSPDF , __FILE__ , __LINE__ ) ;
    }
  //
  const std::size_t          N    { size() } ;
  /// get obs 
  std::unique_ptr<RooArgSet> obs     { pdf().getObservables ( data ) } ;
  const RooArgSet*           normset { this->normalization () ? this->normalization() : obs.get() } ;
  const RooAbsPdf*           pdffun  { &pdf() } ;
  //
  // value of N-components 
  std::vector<double>  cmp_val ( N , 0.0 ) ;
  // the matrix W
  TMatrixDSym W { static_cast<Int_t> ( N ) };
  // loop over the entries
  const bool        weighted = data.isWeighted () ;
  const std::size_t nEntries = data.numEntries () ;
  //
  long double sumw = 0 ;
  Ostap::Utils::ProgressBar bar ( nEntries , progress ) ;
  for ( std::size_t entry = 0 ; entry < nEntries ; ++entry, ++bar )
    {
      //
      const RooArgSet* item = data.get  ( entry ) ;
      Ostap::Assert ( nullptr != item                       ,
		      "Invalid/null  entry in datatset!"    ,
		      "Ostap::Utils::COWs"                  ,
		      INVALID_ABSDATA , __FILE__ , __LINE__ ) ;
      //
      const double weight = weighted ? data.weight() : 1.0 ;
      if ( !weight ) { continue ; }
      //
      sumw += weight ;
      //
      ::assign ( *obs , *item ) ;
      //
      // evalute the total pdf 
      const double pdf_val = normset ? pdffun->getVal ( normset ) : pdffun->getVal () ;
      //
      if ( !pdf_val ) { continue ; } // skip it...  Is it correct? 
      //      
      // evaluate all individual components
      for ( std::size_t i = 0 ; i < N ; ++i ) 
	{
	  const RooAbsReal* r  = static_cast<const RooAbsReal*> ( components().at ( i ) ) ;
	  cmp_val [ i ] = normset ? r->getVal ( normset ) : r->getVal () ;	  
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
	      const double lval = cmp_val [ l ] ;
	      const double wkl  = factor * kval * lval ;
	      // strange "feature" of TMatrixTSym ....
	      W ( l , k ) += 1.0 * wkl ;
	      W ( k , l ) += 1.0 * wkl ;
	    }
	} 
    }
  //
  // W.Print ( "vvv" ) ;
  // ensure it is symmetrical 
  for ( std::size_t k = 0 ; k < N ; ++k )
    { for ( std::size_t l = 0 ; l < k ; ++ l )
	{ W ( l , k ) = 1.0 * W ( k , l ) ;  } }
  //
  // W.Print ( "vvv" ) ;
  //
  if ( weighted  ) { W *= ( 1.0L / sumw     ) ; }
  else             { W *= ( 1.0L / nEntries ) ; }
  //
  // W.Print ( "vvv" ) ;
  // invert the matrix
  double det = 0 ;
  W.Invert ( &det ) ;
  Ostap::Assert ( W.IsValid ()                          ,
		  "The matrix `W' cannot be inverted!"  ,
		  "Ostap::Utils::COWs"                  ,
		  INVALID_ABSDATA , __FILE__ , __LINE__ ) ;
  /// finally store matrix A
  m_A.ResizeTo ( W ) ;
  m_A = W ;
  // m_A.Print ( "vvv" ) ;    
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Utils::COWs::COWs
( const Ostap::Utils::COWs& right )
  : RooFun   ( right ) 
  , m_cmps   {} 
  , m_A ( right.m_A ) 
{
   /// get all components 
   m_cmps = std::make_unique<RooArgList> ( pdf().pdfList () ) ;
}
// ============================================================================
/* "Recovery constructor
 *  @param addpdf input extended RooAddPdf 
 *  @param observables  observables set  
 *  @param normalzation normalisation set 
 *  @param cows         the symmetric matrix of weights 
 */
// ============================================================================
Ostap::Utils::COWs::COWs
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
Ostap::Utils::COWs::COWs
( const RooAddPdf&        addpdf        ,
  const RooAbsCollection& observables   , 
  const RooAbsCollection* normalization ,
  const TMatrixDSym&      cows          )
  : RooFun ( addpdf , observables , normalization )
  , m_cmps        {} 
  , m_A           { cows } 
{
  /// get all components 
  m_cmps = std::make_unique<RooArgList> ( pdf().pdfList () ) ;
  // ==========================================================================
  // check the components
  for ( const RooAbsArg* cmp : *m_cmps )
    {
      Ostap::Assert ( nullptr != cmp ,
		      "Invaild component!" , 
		      "Ostap::Utils::COWs" ,
		      INVALID_ABSARG  , __FILE__ , __LINE__ ) ;
      const RooAbsPdf*  a = dynamic_cast<const RooAbsPdf*> ( cmp ) ; 
      Ostap::Assert ( nullptr != a                         , 
                      "Invalid component!"                 , 
                      "Ostap::Utils::COWs"                 , 
                      INVALID_ABSPDF , __FILE__ , __LINE__ ) ;
    }
  // ==========================================================================
  // check the matrix
  Ostap::Assert ( m_A.IsValid()                          ,
		  "Input matrix is invalid"              ,
		  "Ostap::Utils::COWs"                   ,
		  INVALID_TMATRIX  , __FILE__ , __LINE__ ) ;
  // check the matrix
  Ostap::Assert ( m_A.GetNrows() == size () &&
		  m_A.GetNcols() == size ()              , 
		  "Wrong structure of input matrix"      ,
		  "Ostap::Utils::COWs"                   ,
		  INVALID_TMATRIX  , __FILE__ , __LINE__ ) ;
}
// ============================================================================
/* protected constructor fro SPLOT 
 *  @param addpdf input extended RooAddPdf 
 *  @param observables  observables set  
 *  @param normalzation normalisation set 
 *  @param cows         the symmetric matrix of weights 
 */
// ============================================================================
Ostap::Utils::COWs::COWs
( const RooAddPdf&        addpdf        ,
	const RooAbsCollection& observables   , 
	const RooAbsCollection* normalization ) 
  : RooFun ( addpdf , observables , normalization )
  , m_cmps        {} 
{
  /// get all components 
  m_cmps = std::make_unique<RooArgList> ( pdf().pdfList () ) ;
  /// check the components
  for ( const RooAbsArg* cmp : *m_cmps )
    {
      Ostap::Assert ( nullptr != cmp ,
                      "Invaild component!" , 
                      "Ostap::Utils::COWs" ,
                      INVALID_ABSARG  , __FILE__ , __LINE__ ) ;
      const RooAbsPdf*  a = dynamic_cast<const RooAbsPdf*> ( cmp ) ; 
      Ostap::Assert ( nullptr != a                         , 
                      "Invalid component!"                 , 
                      "Ostap::Utils::COWs"                 , 
                      INVALID_ABSPDF , __FILE__ , __LINE__ ) ;
    }
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Utils::COWs::~COWs() {} ;
// ============================================================================
// cl0one 
// ============================================================================
Ostap::Utils::COWs* 
Ostap::Utils::COWs::clone() const 
{ return new Ostap::Utils::COWs ( *this ) ; }

// =============================================================================
/* Add sPlot/COWs  information to the tree 
 *  @param tree  input tree 
 *  @param cows  COWs object 
 *  @param names names for branches
 *  @param mapping mapping of branch names to RooFit varibabls
 *  @param progress  progress bar 
 *  @return StatusCode
 */
// ==============================================================================
Ostap::StatusCode
Ostap::Trees::add_branch
( TTree*                            tree     ,     
  const Ostap::Utils::COWs&         cows     ,
  const std::vector<std::string>&   names    ,  
  const Ostap::Trees::DCT&          mapping  ,
  const Ostap::Utils::ProgressConf& progress )
{
  if ( nullptr == tree ) { return INVALID_TREE  ; } 
  /// use local version 
  const std::unique_ptr<Ostap::Utils::COWs> the_cows { cows.clone() } ;
  /// the size of the problem
  const std::size_t N { the_cows->size() } ; 
  Ostap::Assert ( N == names.size()          , 
                  "Invalid vector of names"  , 
                  "Ostap::Trees::add_branch" ,
                  INVALID_NAME , __FILE__  , __LINE__ ) ; 
  //                   bname/0   bvalue/1 
  typedef std::tuple<std::string,double> ITEM  ;
  typedef std::vector<ITEM>                             ITEMS ; 
  ITEMS items {} ; items.reserve ( N ) ;
  // ==========================================================================
  for ( const auto& name : names ) { items.emplace_back ( name , 0.0 ) ; }
  // create branches 
  typedef std::vector<TBranch*> BRANCHES ; 
  BRANCHES branches {} ; branches.reserve ( N ) ;
  for ( auto& item : items ) 
    {
      const std::string& bname = std::get<0> ( item  ) ;
      const std::string  bspec { bname  + "/D" } ;
      TBranch* branch = tree->Branch ( bname.c_str() , &std::get<1> ( item ) , bspec.c_str () ) ; // get<1>
      Ostap::Assert ( nullptr != branch                          ,
		      "Cannot create branch '" + bname  + "'"    ,
		      "Ostap::Trees::add_branch"                 , 
		      CANNOT_CREATE_BRANCH , __FILE__ , __LINE__ ) ;
      branches.push_back ( branch ) ;
    } //
  // =========================================================================
  /// observables 
  std::unique_ptr<RooArgSet>    obsset  { std::make_unique<RooArgSet> ( the_cows->observables() ) } ;
  const Ostap::Trees::RooGetter getter  ( mapping , *obsset , tree ) ;
  const RooArgSet*              normset { the_cows->normalization ()   ? the_cows->normalization () : obsset.get() } ;
  const RooAddPdf*              pdffun  { &the_cows->pdf() } ;  
  // 
  // loop over the TTree entries
  const Long64_t nentries = tree->GetEntries(); 
  Ostap::Utils::ProgressBar bar ( nentries , progress  ) ;
  TVectorD cmpvals { static_cast<Int_t> ( N ) } ; 
  for ( Long64_t entry = 0 ; entry < nentries ; ++entry , ++bar )
    {
      if ( tree->GetEntry ( entry ) < 0 ) { break ; };
      /// assign the observables 
      getter.assign ( *obsset, tree ) ;
      // evalaute all individual components:
      for ( std::size_t i = 0 ; i < N ; ++i )
	{
	  const RooAbsReal* cmp = static_cast<const RooAbsReal*> ( the_cows->components().at ( i ) ) ;
	  // get the value of componet 
	  cmpvals [ i ] = normset ? cmp->getVal ( normset ) : cmp->getVal() ;
	}
      /// get the total PDF
      const double total = normset ? pdffun->getVal ( normset ) : pdffun->getVal() ;
      TVectorD sweights { the_cows->A() * cmpvals } ; sweights *= 1.0L / total ; 
      for ( std::size_t i = 0 ; i < N ; ++i ) { std::get<1> ( items [ i ] ) = sweights [ i ] ; }
      /// (3) commit branches 
      for ( auto* branch : branches ) { branch->Fill() ; }
      // ====================================================================
    } //                                 The end of th loop over Tree-entries
  // ========================================================================
  return Ostap::StatusCode::SUCCESS ;
}

// ============================================================================
//                                                                      The END 
// ============================================================================
