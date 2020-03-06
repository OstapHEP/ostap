// ============================================================================
// Include files 
// ============================================================================
// STD&STL
// ============================================================================
#include <cmath>
#include <climits>
#include <tuple>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Math.h"
#include "Ostap/Tmva.h"
#include "Ostap/Iterator.h"
#include "Ostap/Formula.h"
#include "Ostap/FormulaVar.h"
#include "Ostap/Notifier.h"
// ============================================================================
// TMVA
// ============================================================================
#include "TMVA/Reader.h"
// ============================================================================
// ROOT
// ============================================================================
#include "TTree.h"
#include "TBranch.h"
// ============================================================================
// RooFit 
// ============================================================================
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooDataSet.h"
// ============================================================================
namespace
{
  // ===========================================================================
  static_assert ( std::numeric_limits<float>::is_specialized , 
                  "std::numeric_limits<float> is not specialized" ) ;
  // ==========================================================================
  constexpr double s_max =   std::numeric_limits<float>::max() ;  
  constexpr double s_min = - std::numeric_limits<float>::max() ;
  // ==========================================================================
  static_assert ( s_max > 0 , "std::numeric_limits<float>::max is too small" );
  static_assert ( s_min < 0 , "std::numeric_limits<float>::max is too small" );
  // ==========================================================================
  /** @typedef VARIABLE 
   *  helper structure to keep "variable":  name, accessor and placeholder
   */
  typedef std::tuple<std::string,RooAbsReal*,float> VARIABLE ;
  // ==========================================================================
  /** @typedef VARIABLES 
   *  list of all input variables 
   */
  typedef std::vector<VARIABLE>  VARIABLES  ;
  // ==========================================================================
  /// actual type for the reader 
  typedef TMVA::Reader           TMVAReader ;
  // ==========================================================================
  class READER 
  {
    // ========================================================================
  public: 
    // ========================================================================
    READER ( RooDataSet& data                     ,
             const Ostap::TMVA::MAP& inputs       , 
             const Ostap::TMVA::MAP& weight_files )
      : m_data         ( &data        ) 
      , m_inputs       ( inputs       )
      , m_weight_files ( weight_files )
    {}
    // prepare it  for usage 
    Ostap::StatusCode build ( const std::string& options = "" ) 
    {
      //
      RooArgList       varlst ;
      const RooArgSet* varset = m_data->get() ;
      if ( nullptr == varset ) { return Ostap::TMVA::InvalidDataSet ; }
      Ostap::Utils::Iterator iter ( *varset );
      // RooAbsArg*   coef = 0 ;
      while ( RooAbsArg* coef = iter.static_next<RooAbsArg>() ) { varlst.add ( *coef ); }
      //
      // 1)  create variables 
      for ( const auto& i : m_inputs ) 
      {
        const std::string& name   = i.first  ;
        const std::string formula { i.second.empty() ? name : i.second } ;
        RooAbsReal* var = nullptr ;
        /// primitive variable
        if ( std::string::npos == formula.find_first_of( "+-*/&|%()[] ") ) 
        { var = (RooAbsReal*) varset->find ( name.c_str ()  ) ; }
        else 
        {
          auto _v = std::make_unique<Ostap::FormulaVar>( formula , varlst , false ) ;
          if ( !_v->ok() ) { return Ostap::TMVA:: InvalidFormula ;}      
          var    = _v.get() ;
          _vars.push_back ( std::move ( _v ) ) ;
        }
        //
        if ( nullptr == var ) { return Ostap::TMVA::InvalidVariable ; }
        //
        m_variables.push_back ( std::make_tuple ( name , var , 0.0f ) ) ;  
      }
      //
      // 2) create the actual reader 
      m_reader = std::make_unique<TMVAReader>( options ) ;
      //
      //
      // 3) connect the reader with names&placeholders
      for ( auto& v : m_variables ) 
      { 
        m_reader->AddVariable ( std::get<0> ( v ) , &std::get<2>( v ) ) ;
      }
      //
      // 4) book   TMVA methods 
      for ( const auto& p : m_weight_files ) 
      { 
        auto m = m_reader->BookMVA   ( p.first , p.second ) ; 
        if  ( nullptr == m ) { return Ostap::TMVA::InvalidBookTMVA ; }
        m_methods.push_back ( p.first ) ;
      }
      //
      return Ostap::StatusCode::SUCCESS ;
    }
    // 
  public:
    // ========================================================================
    const std::vector<std::string> methods      () const { return m_methods      ; }
    TMVAReader*                    reader       () const { return m_reader.get() ; }
    const Ostap::TMVA::MAP&        inputs       () const { return m_inputs       ; }
    const Ostap::TMVA::MAP&        weight_files () const { return m_weight_files ; }
    VARIABLES&                     variables    ()       { return m_variables    ; }
    // ========================================================================
  private:
    // ========================================================================     
    Ostap::TMVA::MAP         m_inputs                   ;
    Ostap::TMVA::MAP         m_weight_files             ;
    std::vector<std::string> m_methods      {}          ;
    const RooAbsData*        m_data         { nullptr } ;
    // ========================================================================    
  private: // cache 
    // ========================================================================
    std::vector<std::unique_ptr<Ostap::FormulaVar> > _vars {} ;
    // ========================================================================
  private:  
    // ========================================================================    
    VARIABLES                   m_variables  {}           ;
    std::unique_ptr<TMVAReader> m_reader     { nullptr }  ;
    // ========================================================================
  } ;  
  // ==========================================================================
  Ostap::StatusCode _add_response_ 
  ( RooDataSet&                     data      , 
    READER&                         reader    ,
    const std::string&              prefix    , 
    const std::string&              suffix    , 
    const double                    aux       )
  {
    //
    const unsigned long long nEntries = data.numEntries() ;
    if  ( 0 == nEntries || reader.methods().empty() ) { return Ostap::StatusCode::SUCCESS ; }
    //
    RooArgSet tmva_vars;
    std::map<std::string , std::unique_ptr<RooRealVar> > varmap ;
    for ( const auto& m : reader.methods() )
    {
      const std::string vname = prefix + m + suffix ;
      const std::string vdesc = "Response of TMVA/" + m + " method" ;
      auto v = std::make_unique<RooRealVar> ( vname.c_str() , vdesc.c_str() , 0 , s_min  , s_max ) ;
      varmap [ m ] = std::move ( v ) ;
      tmva_vars.add( *varmap[m] ) ;
    }
    //
    auto tmva_ds = std::make_unique<RooDataSet>( "",  "" , tmva_vars ) ;
    //
    for ( unsigned long long entry = 0 ; entry < nEntries ; ++entry ) 
    {
      if ( 0 == data.get( entry ) ) { return Ostap::TMVA::InvalidEntry ; }
      //
      for ( auto& e : reader.variables() ) 
      { std::get<2>(e) = std::get<1> ( e )->getVal() ; }
      // 
      // call TMVA here ... 
      for ( auto& e : varmap ) 
      {
        const double r = 
          reader.reader()->EvaluateMVA ( e.first.c_str () , aux ) ; // EVALUATE TMVA! 
        e.second->setVal ( r )  ;                                         // ATTENTION HERE! 
      }
      //
      tmva_ds->add ( tmva_vars ) ;
    }
    //
    if ( 0 < tmva_ds->numEntries() ) { data.merge ( tmva_ds.get () ) ; }
    //
    return Ostap::StatusCode::SUCCESS ;
  }
  // ==========================================================================
  // Chopping 
  // ==========================================================================
  typedef std::vector<READER>   READERS ;
  // ==========================================================================
  Ostap::StatusCode _add_chopping_response_ 
  ( RooDataSet&        data     , 
    RooAbsReal&        chopping ,
    RooCategory&       category ,
    READERS&           readers  ,
    const std::string& prefix   , 
    const std::string& suffix   ,
    const double       aux      )
  {
    //
    const unsigned long long nEntries = data.numEntries() ;
    if  ( 0 == nEntries ) { return Ostap::StatusCode::SUCCESS ; }
    //
    RooArgSet tmva_vars;
    std::map<std::string , std::unique_ptr<RooRealVar> > varmap ;
    //
    for ( const auto& m : readers[0].methods() )
    {
      const std::string vname = prefix + m + suffix ;
      const std::string vdesc = "Response of TMVA/" + m + " method" ;
      auto v = std::make_unique<RooRealVar> ( vname.c_str() , vdesc.c_str() , 0 , s_min  , s_max ) ;
      varmap [ m ] = std::move ( v ) ;
      tmva_vars.add( *varmap[m] ) ;
    }
    //
    tmva_vars.add ( category ) ;
    const unsigned int N = readers.size() ;
    //
    auto tmva_ds = std::make_unique<RooDataSet>( "",  "" , tmva_vars ) ;
    //
    for ( unsigned long long entry = 0 ; entry < nEntries ; ++entry ) 
    {
      if ( nullptr == data.get( entry ) ) { return Ostap::TMVA::InvalidEntry ; }
      //
      const double chopval  = chopping.getVal() ;
      if ( !Ostap::Math::islong ( chopval ) ) { return Ostap::TMVA::InvalidChoppingCategory ; }
      const long     choplong = std::lround ( chopval ) ;
      const unsigned index    = choplong % N ;
      //
      category.setIndex ( index ) ;
      //
      READER& reader = readers[index] ;
      for ( auto& e : reader.variables() ) 
      { std::get<2>(e) = std::get<1> ( e )->getVal() ; }
      // 
      // call TMVA here ... 
      for ( auto& e : varmap ) 
      {
        const double r = 
          reader.reader()->EvaluateMVA ( e.first.c_str () , aux ) ; // EVALUATE TMVA! 
        e.second->setVal ( r )  ;                                           // ATTENTION HERE! 
      }
      //
      tmva_ds->add ( tmva_vars ) ;
    }
    //
    if ( 0 < tmva_ds->numEntries() ) { data.merge ( tmva_ds.get () ) ; }
    //
    return Ostap::StatusCode::SUCCESS ;
  }
  // ==========================================================================
  /** @typedef VARIABLE2 
   *  helper structure to keep "variable":  name, accessor and placeholder
   */
  typedef std::tuple<std::string,Ostap::Formula*,float> VARIABLE2;
  // ==========================================================================
  /** @typedef VARIABLES2 
   *  list of all input variables 
   */
  typedef std::vector<VARIABLE2> VARIABLES2 ;
  // ==========================================================================
  class READER2
  {
    // ========================================================================
  public: 
    // ========================================================================
    READER2 ( TTree*                  data         ,
              const Ostap::TMVA::MAP& inputs       , 
              const Ostap::TMVA::MAP& weight_files )
      : m_data         ( data         ) 
      , m_inputs       ( inputs       )
      , m_weight_files ( weight_files )
    {}
    // prepare it  for usage 
    Ostap::StatusCode build ( const std::string& options = "" ) 
    {
      //
      // 1)  create variables 
      for ( const auto& i : m_inputs ) 
      {
        const std::string& name = i.first  ;
        const std::string  formula { i.second.empty() ? name : i.second } ;
        //
        std::string fname = name + "_formula" ;
        auto _v = std::make_unique<Ostap::Formula>( fname , formula , m_data ) ;
        if ( !_v->ok() ) { return Ostap::TMVA:: InvalidFormula ; }      
        Ostap::Formula* var = _v.get() ;
        _vars.push_back ( std::move ( _v ) ) ;        
        //
        if ( nullptr == var ) { return Ostap::TMVA::InvalidVariable ; }
        //
        m_variables.push_back ( std::make_tuple ( name , var , 0.0f ) ) ;  
      }
      //
      // 2) create the actual reader 
      m_reader = std::make_unique<TMVAReader>( options ) ;
      //
      //
      // 3) connect the reader with names&placeholders
      for ( auto& v : m_variables ) 
      { 
        m_reader->AddVariable ( std::get<0> ( v ) , &std::get<2>( v ) ) ;
      }
      //
      // 4) book   TMVA methods 
      for ( const auto& p : m_weight_files ) 
      { 
        auto m = m_reader->BookMVA   ( p.first , p.second ) ; 
        if  ( nullptr == m ) { return Ostap::TMVA::InvalidBookTMVA ; }
        m_methods.push_back ( p.first ) ;
      }
      //
      return Ostap::StatusCode::SUCCESS ;
    }
    // 
  public:
    // ========================================================================
    const std::vector<std::string> methods      () const { return m_methods      ; }
    TMVAReader*                    reader       () const { return m_reader.get() ; }
    const Ostap::TMVA::MAP&        inputs       () const { return m_inputs       ; }
    const Ostap::TMVA::MAP&        weight_files () const { return m_weight_files ; }
    VARIABLES2&                    variables    ()       { return m_variables    ; }
    // ========================================================================
  private:
    // ========================================================================     
    Ostap::TMVA::MAP         m_inputs                   ;
    Ostap::TMVA::MAP         m_weight_files             ;
    std::vector<std::string> m_methods      {}          ;
    TTree*                   m_data         { nullptr } ;
    // ========================================================================    
  private:
    // ========================================================================
    std::vector<std::unique_ptr<Ostap::Formula> > _vars {} ;
    // ========================================================================
  private:  
    // ========================================================================    
    VARIABLES2                  m_variables  {}           ;
    std::unique_ptr<TMVAReader> m_reader     { nullptr }  ;
    // ========================================================================
  } ;  
  // ===========================================================================
  Ostap::StatusCode _add_response_ 
  ( TTree*                          tree      , 
    READER2&                        reader    ,
    const std::string&              prefix    , 
    const std::string&              suffix    , 
    const double                    aux       )
  {
    //
    const Long64_t nEntries = tree->GetEntries() ;
    if  ( 0 == nEntries || reader.methods().empty() ) { return Ostap::StatusCode::SUCCESS ; }
    //
    typedef std::tuple<TBranch*,std::string,double> Branch   ;
    typedef std::vector<Branch>                     Branches ;
    //
    Branches branches { reader.methods().size() } ;
    unsigned short index = 0 ;
    for(  auto& branch : branches ) 
    {
      const std::string method = reader.methods()[ index ] ;
      const std::string bname  = prefix + method + suffix       ;
      //
      std::get<1>( branch ) = method ;
      std::get<2>( branch ) = 0.0    ;
      std::get<0>( branch ) = tree->Branch ( bname.c_str() , &std::get<2> ( branch) , ( bname+"/D").c_str() ) ;
      //
      if ( !std::get<0>( branch ) ) { return Ostap::TMVA::InvalidBranch ; } 
      //
      ++index ;
    }
    //
    Ostap::Utils::Notifier  notifier { tree } ;
    for ( auto& e : reader.variables () ) { notifier.add ( std::get<1> ( e ) ) ; }
    //
    for ( Long64_t entry = 0 ; entry < nEntries ; ++entry ) 
    {
      if ( tree->GetEntry ( entry ) < 0 ) { break ; }
      //
      // prepare TMVA input 
      for ( auto& e : reader.variables() ) 
      { std::get<2>(e) = std::get<1> ( e )->evaluate () ; }
      //
      // evalaute TMVA and fill braches
      for ( auto& branch : branches ) 
      {
        // call TMVA here ... 
        const double result    = reader.reader()->EvaluateMVA ( std::get<1>( branch ). c_str () , aux ) ;
        std::get<2> ( branch ) = result ;
        std::get<0> ( branch ) -> Fill () ;  
      }
      // 
    }
    //
    return Ostap::StatusCode::SUCCESS ;
  } ;
  // ==========================================================================
  // Chopping 
  // ==========================================================================
  typedef std::vector<READER2>   READERS2 ;
  // ==========================================================================
  Ostap::StatusCode _add_chopping_response_ 
  ( TTree*             tree     , 
    Ostap::Formula&    chopping ,
    const std::string& category ,
    READERS2&          readers  ,
    const std::string& prefix   , 
    const std::string& suffix   ,
    const double       aux      )
  {
    //
    const Long64_t nEntries = tree->GetEntries() ;
    if  ( 0 == nEntries ) { return Ostap::StatusCode::SUCCESS ; }
    //    
    typedef std::tuple<TBranch*,std::string,double> Branch   ;
    typedef std::vector<Branch>                     Branches ;
    //
    // Variables
    //
    Branches branches {readers[0].methods().size() } ;
    unsigned short index = 0 ;
    for(  auto& branch : branches ) 
    {
      const std::string method = readers[0].methods()[ index ] ;
      const std::string bname  = prefix + method + suffix       ;
      //
      std::get<1>( branch ) = method ;
      std::get<2>( branch ) = 0.0    ;
      std::get<0>( branch ) = tree->Branch ( bname.c_str() , &std::get<2> ( branch) , ( bname+"/D").c_str() ) ;
      //
      if ( !std::get<0>( branch ) ) { return Ostap::TMVA::InvalidBranch ; } 
      //
      ++index ;
    }
    //
    Ostap::Utils::Notifier  notifier{ tree , &chopping } ;
    //
    for ( auto&  reader :readers ) 
    { for ( auto& e : readers[0].variables () ) { notifier.add ( std::get<1> ( e ) ) ; } }
    //
    // category in Tree:
    //
    UInt_t i_category = 0 ;
    TBranch* bcat = tree->Branch( category.c_str() , &i_category , ( category + "/i" ).c_str() ) ;
    //
    if ( !bcat ) { return Ostap::TMVA::InvalidBranch ; } 
    //
    const unsigned int N = readers.size() ;
    //
    for ( Long64_t entry = 0 ; entry < nEntries ; ++entry ) 
    {
      if ( tree->GetEntry ( entry ) < 0 ) { break ; }
      //
      const double  chopval = chopping.evaluate() ;
      if ( !Ostap::Math::islong ( chopval ) ) { return Ostap::TMVA::InvalidChoppingCategory ; }
      const long         choplong = std::lround ( chopval ) ;
      const unsigned int index    = choplong % N ;
      i_category =  index  ;
      bcat       -> Fill() ;
      //
      // prepare TMVA input 
      READER2& reader = readers[index] ;
      for ( auto& e : reader.variables() ) 
      { std::get<2>(e) = std::get<1> ( e )->evaluate () ; }
      //
      // evaluate TMVA and fill braches
      for ( auto& branch : branches ) 
      {
        // call TMVA here ... 
        const double result    = reader.reader()->EvaluateMVA ( std::get<1>( branch ). c_str () , aux ) ;
        std::get<2> ( branch ) = result ;
        std::get<0> ( branch ) -> Fill () ;  
      }
    }
    //
    return Ostap::StatusCode::SUCCESS ;
  }
  // ==========================================================================
}
// ============================================================================
/*  Add TMVA response to dataset 
 *  The  function add variables  "prefix+methods+suffix" that 
 *  are the responses of TMVA. TMNVA  configurtaion for methods 
 *  is read from the trained (xml) weight-files 
 *  @param data         (UPDATE) dataset 
 *  @param inputs       (INPUT) map  { varname : formula     }  
 *  @param weight_files (INPUT) map  { method  : weight_file }  
 *  @param prefix       (INPUT) the prefix for added varibales 
 *  @param suffix       (INPUT) the suffix for added varibales 
 */ 
// ============================================================================
Ostap::StatusCode Ostap::TMVA::addResponse
( RooDataSet& data                     ,
  const Ostap::TMVA::MAP& inputs       , 
  const Ostap::TMVA::MAP& weight_files , 
  const std::string&       options      , 
  const std::string&      prefix       , 
  const std::string&      suffix       ,
  const double            aux          )
{
  // create the helper structure  
  READER reader  ( data , inputs , weight_files ) ;
  Ostap::StatusCode sc =  reader.build ( options ) ;
  if ( sc.isFailure() ) { return sc ; }
  //
  return _add_response_ ( data    ,
                          reader  , 
                          prefix  , 
                          suffix  , 
                          aux     ) ;
}
// ============================================================================
/*  Add TMVA response to TTree
 *  The  function add branches <code>prefix+method+suffix</code> that 
 *  are the responses of TMVA. 
 *  - TMVA  configuration for all methods 
 *  is read from the trained (xml) weight-files 
 *  @param tree         (UPDATE) input TTree
 *  @param inputs       (INPUT) map  { varname : formula     }  
 *  @param weight_files (INPUT) map  { method  : weight_file }  
 *  @param prefix       (INPUT) the prefix for added variables 
 *  @param suffix       (INPUT) the suffix for added variables 
 *  @param aux          (INPUT) obligatory for the cuts method
 *                              where it represents the efficiency cutoff
 */ 
// ============================================================================
Ostap::StatusCode Ostap::TMVA::addResponse
( TTree*                  tree          ,
  const Ostap::TMVA::MAP& inputs        , 
  const Ostap::TMVA::MAP& weight_files  ,
  const std::string&      options       ,
  const std::string&      prefix        , 
  const std::string&      suffix        , 
  const double            aux           ) 
{
  if ( nullptr == tree ) { return InvalidTree ; }
  // create the helper structure  
  READER2 reader  ( tree , inputs , weight_files ) ;
  Ostap::StatusCode sc =  reader.build ( options ) ;
  if ( sc.isFailure() ) { return sc ; }
  //
  return _add_response_ ( tree    ,
                          reader  , 
                          prefix  , 
                          suffix  , 
                          aux     ) ; 
}
// ========================================================================
// Chopping 
// ============================================================================
Ostap::StatusCode Ostap::TMVA::addChoppingResponse 
( RooDataSet&              data         ,
  RooAbsReal&              chopping     , // category function 
  RooCategory&             category     , // category variable  
  const unsigned short     N            , // number of categories 
  const Ostap::TMVA::MAP&  inputs       , // mapping of input variables 
  const Ostap::TMVA::MAPS& weight_files ,
  const std::string&       options      , 
  const std::string&       prefix       , 
  const std::string&       suffix       ,
  const double             aux          )
{
  // ==========================================================================
  if  ( 0 == N || N != weight_files.size() ) { return InvalidChoppingWeightFiles ; }
  // ==========================================================================
  //
  READERS readers ;
  readers.reserve ( N ) ;
  //
  for ( const auto& wfs : weight_files )
  { readers.emplace_back ( data , inputs ,  wfs ) ;}
  //
  // initialize the  readers:
  bool first = true ;
  for ( auto& r : readers ) 
  {
    Ostap::StatusCode sc =  r.build( first ? options : "" ) ;
    if   ( sc.isFailure () ) { return sc ; }  
    first = false ;
  }
  //
  return _add_chopping_response_ ( data      ,
                                   chopping  , 
                                   category  ,
                                   readers   , 
                                   prefix    , 
                                   suffix    ,
                                   aux       )  ;
} 
// ============================================================================
/*  Add TMVA/Chopping response to TTree
 *  The  function add branches <code>prefix+method+suffix</code> that 
 *  are the responses of TMVA. 
 *  - TMVA  configuration for all methods 
 *  is read from the trained (xml) weight-files 
 *  @param tree         (UPDATE) input tree 
 *  @param chopping     (INPUT) chopping variable/expression
 *  @param chopping     (INPUT) chopping category 
 *  @param N            (INPUT) number of categories 
 *  @param inputs       (INPUT) map  { varname : formula     }  
 *  @param weight_files (INPUT) map  { method  : weight_file }  
 *  @param prefix       (INPUT) the prefix for added variables 
 *  @param suffix       (INPUT) the suffix for added variables 
 *  @param aux          (INPUT) obligatory for the cuts method
 *                              where it represents the efficiency cutoff
 */ 
// ============================================================================
Ostap::StatusCode Ostap::TMVA::addChoppingResponse 
( TTree*                   tree          ,
  const std::string&       chopping      , // category function 
  const std::string&       category_name , // category variable 
  const unsigned short     N             , // number of categories 
  const Ostap::TMVA::MAP&  inputs        , // mapping of input variables 
  const Ostap::TMVA::MAPS& weight_files  ,
  const std::string&       options       ,
  const std::string&       prefix        , 
  const std::string&       suffix        ,
  const double             aux           ) 
{
  if ( nullptr == tree ) { return InvalidTree ; }
  // ==========================================================================
  if  ( 0 == N || N != weight_files.size() ) { return InvalidChoppingWeightFiles ; }
  // ==========================================================================
  //
  auto chop_var = std::make_unique<Ostap::Formula> ( chopping , chopping  , tree ) ;
  if ( !chop_var || !chop_var->ok() ) { return Ostap::TMVA:: InvalidFormula ; }      
  //
  READERS2 readers ;
  readers.reserve ( N ) ;
  //
  for ( const auto& wfs : weight_files )
  { readers.emplace_back ( tree , inputs ,  wfs ) ;}
  //
  // initialize the  readers:
  bool first = true ;
  for ( auto& r : readers ) 
  {
    Ostap::StatusCode sc =  r.build( first ? options : "" ) ;
    if   ( sc.isFailure () ) { return sc ; }  
    first = false ;
  }
  //
  return _add_chopping_response_ ( tree          ,
                                   *chop_var     , 
                                   category_name ,
                                   readers       , 
                                   prefix        , 
                                   suffix        ,
                                   aux           );
}



// ============================================================================
//                                                                      The END 
// ============================================================================



