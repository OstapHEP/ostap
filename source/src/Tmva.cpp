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
// ============================================================================
// TMVA
// ============================================================================
#include "TMVA/Reader.h"
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
    Ostap::StatusCode build () 
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
          std::string fname = name + "_formula" ;
          auto _v = std::make_unique<RooFormulaVar>( fname.c_str(), formula.c_str() , varlst ) ;
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
      m_reader = std::make_unique<TMVAReader>() ;
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
    std::vector<std::unique_ptr<RooFormulaVar> > _vars {} ;
    // ========================================================================
  private:  
    // ========================================================================    
    VARIABLES                   m_variables  {}           ;
    std::unique_ptr<TMVAReader> m_reader     { nullptr }  ;
    // ========================================================================
  } ;  
  // ==========================================================================
  typedef std::vector<READER>   READERS ;
  // ==========================================================================
  inline Ostap::StatusCode _add_response_ 
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
  inline Ostap::StatusCode _add_chopping_response_ 
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
      if ( 0 == data.get( entry ) ) { return Ostap::TMVA::InvalidEntry ; }
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
  //
}
// ============================================================================
/*  Add TMVA response to dataset 
 *  The  function add variables  "prefix+methos+suffix" that 
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
  const std::string&      prefix       , 
  const std::string&      suffix       ,
  const double            aux          )
{
  // create the helper structure  
  READER reader  ( data , inputs , weight_files ) ;
  Ostap::StatusCode sc =  reader.build() ;
  if ( sc.isFailure() ) { return sc ; }
  //
  return _add_response_ ( data    ,
                          reader  , 
                          prefix  , 
                          suffix  , 
                          aux     ) ;
}
// ========================================================================
/*  Add TMVA response to dataset 
 *  The  function add variables  "prefix+methos+suffix" that 
 *  are the responses of TMVA. TMNVA  configurtaion for methods 
 *  is read from the trained (xml) weight-files 
 *  @param data         (UPDATE) dataset 
 *  @param inputs       (INPUT) [ (varnname,formula)   , ...  ] 
 *  @param weight_files (INPUT) [ (method,weight_file) , ...  ]
 *  @param prefix       (INPUT) the prefix for added varibales 
 *  @param suffix       (INPUT) the suffix for added varibales 
 */ 
// ========================================================================
Ostap::StatusCode Ostap::TMVA::addResponse
( RooDataSet& data                       ,
  const Ostap::TMVA::PAIRS& inputs       , 
  const Ostap::TMVA::PAIRS& weight_files , 
  const std::string&        prefix       , 
  const std::string&        suffix       ,
  const double              aux          )
{
  MAP _i ;
  for ( const auto& p : inputs     ) { _i[p.first] = p.second ; }
  if  ( _i.size() != inputs.size() ) { return InvalidInputVariables ; }
  //
  return addResponse (  data , _i , weight_files , prefix , suffix , aux ) ;
}
// =======================================================================-----
/*  Add TMVA response to dataset 
 *  The  function add variables  "prefix+methos+suffix" that 
 *  are the responses of TMVA. TMNVA  configurtaion for methods 
 *  is read from the trained (xml) weight-files 
 *  @param data         (UPDATE) dataset 
 *  @param inputs       (INPUT) [ (varnname,formula)   , ...  ] 
 *  @param weight_files (INPUT) map  { method  : weight_file }  
 *  @param prefix       (INPUT) the prefix for added varibales 
 *  @param suffix       (INPUT) the suffix for added varibales 
 */
// =======================================================================-----
Ostap::StatusCode Ostap::TMVA::addResponse
( RooDataSet& data                       ,
  const Ostap::TMVA::PAIRS& inputs       , 
  const Ostap::TMVA::MAP&   weight_files , 
  const std::string&        prefix       , 
  const std::string&        suffix       ,
  const double              aux          )
{
  MAP _i;
  for ( const  auto& p : inputs    ) { _i[p.first] = p.second ; }
  if  ( _i.size() != inputs.size() ) { return InvalidInputVariables ; }
  //
  return addResponse (  data , _i , weight_files , prefix , suffix , aux ) ; 
}
// =======================================================================-----
/*  Add TMVA response to dataset 
 *  The  function add variables  "prefix+methos+suffix" that 
 *  are the responses of TMVA. TMNVA  configurtaion for methods 
 *  is read from the trained (xml) weight-files 
 *  @param data         (UPDATE) dataset 
 *  @param inputs       (INPUT) map  { varname : formula     }  
 *  @param weight_files (INPUT) [ (method,weight_file) , ...  ]
 *  @param prefix       (INPUT) the prefix for added varibales 
 *  @param suffix       (INPUT) the suffix for added varibales 
 */
// =======================================================================-----
Ostap::StatusCode Ostap::TMVA::addResponse
( RooDataSet& data                       ,
  const Ostap::TMVA::MAP&   inputs       , 
  const Ostap::TMVA::PAIRS& weight_files , 
  const std::string&        prefix       , 
  const std::string&        suffix       ,
  const double              aux          )
{
  MAP _w ;
  for ( const auto& p : weight_files    ) { _w[p.first] = p.second ; }
  if ( _w.size() != weight_files.size() ) { return InvalidWeightFiles  ; }
  //
  return addResponse (  data , inputs , _w , prefix , suffix , aux ) ;
}
// ============================================================================
// Chopping 
// ============================================================================
Ostap::StatusCode Ostap::TMVA::addChoppingResponse 
( RooDataSet&              data         ,
  RooAbsReal&              chopping     , // category function 
  RooCategory&             category     , // category variable  
  const unsigned short     N            , // number of categories 
  const Ostap::TMVA::MAP&  inputs       , // mapping of input variables 
  const Ostap::TMVA::MAPS& weight_files ,
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
  { readers.emplace_back ( data , inputs ,  wfs ) ; }
  //
  // initialize the  readers:
  for ( auto& r : readers ) 
  {
    Ostap::StatusCode sc =  r.build() ;
    if   ( sc.isFailure () ) { return sc ; }  
  }
  //
  return _add_chopping_response_ ( data      ,
                                   chopping  , 
                                   category  ,
                                   readers   , 
                                   prefix    , 
                                   suffix    ,
                                   aux       )  ;
  //
  return  Ostap::StatusCode::SUCCESS ;
} 
// ============================================================================
//                                                                      The END 
// ============================================================================



