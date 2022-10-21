// ============================================================================
// Include files 
// ============================================================================
// ROOT/RooFit
// ============================================================================
#include "RVersion.h"
#include "RooAbsPdf.h"
#include "RooConstVar.h"
#include "RooAddPdf.h"
#include "RooRecursiveFraction.h"
#include "RooGaussian.h"
#include "RooFFTConvPdf.h"
#include "RooArgList.h"
#include "RooEfficiency.h"
#include "RooPolyVar.h"
#include "RooPolynomial.h"
#include "RooMultiVarGaussian.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/MoreVars.h"
#include "Ostap/Iterator.h"
// ============================================================================
// local
// ============================================================================
#include "Exception.h"
#include "local_roofit.h"
// ============================================================================
/** @file
 *  Implementation file for objects from the file Ostap/MoreVars.h
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
 *  @date 2020-03-06     
 */
// ============================================================================
namespace 
{
  // ===========================================================================
  const std::string s_EMPTYPARS   = "Vector of coefficients is empty!" ;
  const std::string s_INVALIDPARS = "Invalid parameters!"              ;
  const std::string s_INVALIDPAR  = "Invalid parameter!"               ;
  const std::string s_v1          = "Ostap::MoreRooFit::Bernstein"     ;
  const std::string s_v2          = "Ostap::MoreRooFit::Monotonic"     ;
  const std::string s_v3          = "Ostap::MoreRooFit::Convex"        ;
  const std::string s_v4          = "Ostap::MoreRooFit::ConvexOnly"    ;
  const std::string s_v5          = "Ostap::MoreRooFit::BSpline"       ;
  // ===========================================================================
  class FakeRecursiveFraction : public RooRecursiveFraction 
  {
  public:
    // ========================================================================
    /// constructor from RooRecuriveFraction object 
    FakeRecursiveFraction 
    ( const RooRecursiveFraction& right , const char* newname = 0 )
      : RooRecursiveFraction ( right , newname ) 
    {}
    /// virtual desttructor 
    virtual ~FakeRecursiveFraction () {}
    // ========================================================================
  public: 
    // ========================================================================
    /** get the original fractions form the <code>RooRecursiveFraction</code>
     *  @see RooRecursiveFraction
     */
    RooArgList fractions () const 
    {
      RooArgList result {} ;
      ::copy_real ( _list , result ) ;
      return result ;
    }
    // ======================================================================
  } ;
  // ========================================================================
  /** get the original fractions form the <code>RooRecursiveFraction</code>
   *  @see RooRecursiveFraction
   *  @attention the list has inverse order and the last efraction if
   *  prepended with the constant 1 
   */
  RooArgList fractions
  ( const RooRecursiveFraction& rf ) 
  {
    std::unique_ptr<::FakeRecursiveFraction> fake { new ::FakeRecursiveFraction ( rf ) } ;
    return fake->fractions () ; 
  }
  // ==========================================================================
  class FakeAddPdf : public RooAddPdf 
  {
  public: 
    // ========================================================================
    /// constructor from RooAddPdf object 
    FakeAddPdf ( const RooAddPdf& right  , const char* newname = 0 )
      : RooAddPdf ( right , newname ) 
    {}
    /// virtual desttructor 
    virtual ~FakeAddPdf() {}
    // ========================================================================
  public: 
    // ========================================================================
    /// were recursive fractions used for creation ?
    bool recursive () const { return _recursive ; }
    /// get the original fractions 
    RooArgList fractions ( bool& resursive )  const 
    {
      resursive = _recursive ;
      if ( !_recursive || ::size( _coefList ) <= 1 ) { return _coefList ; }
      //
      RooArgList result {} ;
      //
      // get the last coefficienct 
      //
      const int index_last  = ::size ( _coefList ) - 1 ;
      const RooAbsArg* last = _coefList.at ( index_last )  ;
      Ostap::Assert  ( nullptr != last                , 
                       "Invalid size of _coefList!"   , 
                       "Ostap::MoreRooFit::fractions" ) ;
      //
      const RooRecursiveFraction* rf = 
        dynamic_cast<const RooRecursiveFraction*>( last ) ;
      Ostap::Assert  ( nullptr != rf                , 
                       "Last fraction is not RooRecursiveFraction!"   , 
                       "Ostap::MoreRooFit::fractions" ) ;
      //
      RooArgList tmp { ::fractions ( *rf ) } ;
      Ostap::Assert  ( nullptr != rf                , 
                       "Invalid size of fractions!"  , 
                       "Ostap::MoreRooFit::fractions" ) ;
      //
      for  ( int i = index_last ; 1<= i ; --i ) 
      { result.add ( *tmp.at ( i ) ) ; }
      //
      return result ;
    }
    // ========================================================================
  } ;
  // ========================================================================
  class FakeGaussian : public RooGaussian 
  {
  public: 
    // ======================================================================
    FakeGaussian ( const RooGaussian& pdf , const char* newname = 0) 
      : RooGaussian ( pdf , newname ) 
    {}
    /// virtual desttructor 
    virtual ~FakeGaussian() {}
    // ========================================================================
  public: 
    // ========================================================================
    const RooAbsReal& get_x     () const { return x     .arg() ; }
    const RooAbsReal& get_mean  () const { return mean  .arg() ; }
    const RooAbsReal& get_sigma () const { return sigma .arg() ; }
    // ========================================================================
  } ;  
  // ==========================================================================
  class FakeFFTConvPdf : public RooFFTConvPdf 
  {
  public: 
    // ========================================================================
    FakeFFTConvPdf ( const RooFFTConvPdf& pdf , const char* newname = 0 ) 
      : RooFFTConvPdf ( pdf , newname ) 
    {}
    // ========================================================================
    virtual ~FakeFFTConvPdf() {}
    // ========================================================================
  public:
    // ========================================================================
    RooArgList 
    getpars
    ( double&              shift1 ,
      double&              shift2 ) const 
    {
      shift1 = _shift1 ;
      shift2 = _shift2 ;      
      RooArgList result {} ;
      if ( _xprime.absArg() ) { result.add ( _xprime.arg() ) ; }
      result.add ( _x.arg      () ) ;
      result.add ( _pdf1.arg   () ) ;
      result.add ( _pdf2.arg   () ) ;
      return result ;
    }
  } ;  
  // ==========================================================================
  class FakeEfficiency : public RooEfficiency 
  {
  public: 
    // ========================================================================
    FakeEfficiency ( const RooEfficiency& pdf , const char* newname = 0 ) 
      : RooEfficiency( pdf , newname ) 
    {}
    // ========================================================================
    virtual ~FakeEfficiency () {}
    // ========================================================================
  public:
    // ========================================================================
    const RooAbsReal&     get_eff () const { return _effFunc   .arg () ; }
    const RooAbsCategory& get_cat () const { return _cat       .arg () ; }
    std::string           get_acc () const { return _sigCatName.Data() ; }
    // ========================================================================
  } ;  
  // ==========================================================================
  class FakePolyVar : public RooPolyVar
  {
  public: 
    // ========================================================================
    FakePolyVar ( const RooPolyVar& var , const char* newname = 0 ) 
      : RooPolyVar ( var , newname ) 
    {}
    // ========================================================================
    virtual ~FakePolyVar () {}
    // ========================================================================
  public:
    // ========================================================================
    const RooArgList& coefficients () const { return _coefList ; }
    // ========================================================================
  } ;    
  // ==========================================================================
  class FakePolynomial : public RooPolynomial
  {
  public: 
    // ========================================================================
    FakePolynomial ( const RooPolynomial& var , const char* newname = 0 ) 
      : RooPolynomial ( var , newname ) 
    {}
    // ========================================================================
    virtual ~FakePolynomial () {}
    // ========================================================================
  public:
    // ========================================================================
    const RooArgList& coefficients () const { return _coefList ; }
    // ========================================================================
  } ;    
  // ==========================================================================
  class FakeMultiVarGaussian : public RooMultiVarGaussian
  {
  public: 
    // ========================================================================
    FakeMultiVarGaussian 
    ( const RooMultiVarGaussian& var , const char* newname = 0 ) 
      : RooMultiVarGaussian ( var , newname ) 
    {}
    // ========================================================================
    virtual ~FakeMultiVarGaussian () {}
    // ========================================================================
  public:
    // ========================================================================
    const RooArgList& observables () const { return _x ; }
    TVectorD          mu_vec      () const 
    { syncMuVec() ; return _muVec ; }
    // ========================================================================
  } ;    
  // ==========================================================================
} //                                             The end of anonymous namespace 
// ============================================================================
ClassImp ( Ostap::MoreRooFit::Bernstein  ) ;
ClassImp ( Ostap::MoreRooFit::Monotonic  ) ;
ClassImp ( Ostap::MoreRooFit::Convex     ) ;
ClassImp ( Ostap::MoreRooFit::ConvexOnly ) ;
ClassImp ( Ostap::MoreRooFit::BSpline    ) ;
// ============================================================================
// constructor from the variable, range and list of coefficients
// ============================================================================
Ostap::MoreRooFit::Bernstein::Bernstein 
( const std::string& name  ,
  const std::string& title ,
  RooAbsReal&        xvar  ,
  const RooArgList&  pars  ,
  const double       xmin  , 
  const double       xmax  )
  : RooAbsReal ( name.c_str  () , title.c_str () )
  , m_xvar      ( "x"    , "Dependent"  , this , xvar ) 
  , m_pars      ( "pars" , "Parameters" , this ) 
  , m_bernstein ( ::size ( pars ) - 1 , xmin , xmax ) 
{
  //
  Ostap::Assert ( 1 <= ::size ( pars ) , s_EMPTYPARS  , s_v1 , 510 ) ;
  //
  ::copy_real   ( pars , m_pars , s_INVALIDPAR , s_v1 ) ;
  //
  Ostap::Assert ( m_bernstein.npars() == ::size ( m_pars ), s_INVALIDPARS , s_v1 , 512 ) ;
  //
}
// =============================================================================
// copy constructor 
// =============================================================================
Ostap::MoreRooFit::Bernstein::Bernstein
( const Ostap::MoreRooFit::Bernstein& right , 
  const char*                         name  )
  : RooAbsReal  ( right , name ) 
  , m_xvar      ( "x"    , this , right.m_xvar ) 
  , m_pars      ( "pars" , this , right.m_pars ) 
  , m_bernstein ( right.m_bernstein ) 
{}
// =============================================================================
// default constructor 
// ============================================================================
Ostap::MoreRooFit::Bernstein::Bernstein(){}
// =============================================================================
// destructor  
// ============================================================================
Ostap::MoreRooFit::Bernstein::~Bernstein(){}
// ============================================================================
// clone it!
// ============================================================================
Ostap::MoreRooFit::Bernstein*
Ostap::MoreRooFit::Bernstein::clone ( const char* newname ) const
{ return new Bernstein( *this , newname ) ; }
// ============================================================================
void Ostap::MoreRooFit::Bernstein::setPars() const 
{ ::set_pars ( m_pars , m_bernstein ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::Bernstein::evaluate  () const 
{
  setPars () ;
  const double x = m_xvar ;
  return m_bernstein ( x ) ;
}
// ============================================================================
Int_t    Ostap::MoreRooFit::Bernstein::getAnalyticalIntegral
( RooArgSet&  allVars   , 
  RooArgSet&  analVars  , 
  const char* /* rangeName */ ) const 
{
  return matchArgs  ( allVars , analVars , m_xvar ) ? 1 : 0 ;
}
// ===========================================================================
Double_t Ostap::MoreRooFit::Bernstein::analyticalIntegral 
( Int_t code            , 
  const char* rangeName ) const 
{
  setPars() ;
  const double xmin = m_xvar.min ( rangeName ) ;
  const double xmax = m_xvar.max ( rangeName ) ;
  return m_bernstein.integral ( xmin , xmax ) ;
}
// ============================================================================


// ============================================================================
// constructor from the variable, range and list of coefficients
// ============================================================================
Ostap::MoreRooFit::Monotonic::Monotonic
( const std::string& name       ,
  const std::string& title      ,
  RooAbsReal&        xvar       ,
  const RooArgList&  pars       ,
  const bool         increasing ,
  const double       xmin       , 
  const double       xmax       ,
  RooAbsReal&        a          ,
  RooAbsReal&        b          )
  : RooAbsReal ( name.c_str  () , title.c_str () )
  , m_xvar      ( "x"    , "Dependent"  , this , xvar ) 
  , m_a         ( "a"    , "shift/bias" , this , a    ) 
  , m_b         ( "b"    , "scale"      , this , b    ) 
  , m_pars      ( "pars" , "Parameters" , this ) 
  , m_monotonic ( ::size ( pars ) , xmin , xmax , increasing ) 
{
  //
  Ostap::Assert ( 0 <= ::size ( pars ) , s_EMPTYPARS  , s_v2 , 510 ) ;
  //
  ::copy_real   ( pars , m_pars , s_INVALIDPAR , s_v2 ) ;
  //
  Ostap::Assert ( m_monotonic.npars() == ::size ( m_pars ), s_INVALIDPARS , s_v2 , 512 ) ;
  //
}
// ============================================================================
// constructor from the variable, range and list of coefficients
// ============================================================================
Ostap::MoreRooFit::Monotonic::Monotonic
( const std::string& name       ,
  const std::string& title      ,
  RooAbsReal&        xvar       ,
  const RooArgList&  pars       ,
  const bool         increasing ,
  const double       xmin       , 
  const double       xmax       ,
  const double       a          ,
  const double       b          ) 
  : Monotonic ( name , title , xvar , pars , increasing , xmin , xmax ,
                RooFit::RooConst ( a ) , 
                RooFit::RooConst ( b ) )
{}
// =============================================================================
// copy constructor 
// =============================================================================
Ostap::MoreRooFit::Monotonic::Monotonic
( const Ostap::MoreRooFit::Monotonic& right , 
  const char*                         name  )
  : RooAbsReal  ( right , name ) 
  , m_xvar      ( "x"    , this , right.m_xvar ) 
  , m_a         ( "a"    , this , right.m_a    ) 
  , m_b         ( "b"    , this , right.m_b    ) 
  , m_pars      ( "pars" , this , right.m_pars ) 
  , m_monotonic ( right.m_monotonic ) 
{}
// =============================================================================
// default constructor 
// ============================================================================
Ostap::MoreRooFit::Monotonic::Monotonic (){}
// =============================================================================
// destructor  
// ============================================================================
Ostap::MoreRooFit::Monotonic::~Monotonic(){}
// ============================================================================
// clone it!
// ============================================================================
Ostap::MoreRooFit::Monotonic*
Ostap::MoreRooFit::Monotonic::clone ( const char* newname ) const
{ return new Monotonic ( *this , newname ) ; }
// ============================================================================
void Ostap::MoreRooFit::Monotonic::setPars() const 
{ ::set_pars (  m_pars , m_monotonic ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::Monotonic::evaluate  () const 
{
  setPars () ;
  //
  const double x = m_xvar ;
  if ( x < m_monotonic.xmin () || x > m_monotonic.xmax() ) { return 0 ; }
  const double a = m_a    ;
  const double b = m_b    ;
  //
  return a + b * m_monotonic ( x ) ;
}
// ============================================================================
Int_t    Ostap::MoreRooFit::Monotonic::getAnalyticalIntegral
( RooArgSet&  allVars   , 
  RooArgSet&  analVars  , 
  const char* /* rangeName */ ) const 
{
  return matchArgs  ( allVars , analVars , m_xvar ) ? 1 : 0 ;
}
// ===========================================================================
Double_t Ostap::MoreRooFit::Monotonic::analyticalIntegral 
( Int_t code            , 
  const char* rangeName ) const 
{
  setPars() ;
  //
  const double xmin = std::max ( m_xvar.min ( rangeName ) , m_monotonic.xmin () ) ;
  const double xmax = std::min ( m_xvar.max ( rangeName ) , m_monotonic.xmax () ) ;
  //
  const double a = m_a    ;
  const double b = m_b    ;
  //
  return a * ( xmax - xmin ) + b * m_monotonic.integral ( xmin , xmax ) ;
}
// ============================================================================



// ============================================================================
// constructor from the variable, range and list of coefficients
// ============================================================================
Ostap::MoreRooFit::Convex::Convex
( const std::string& name       ,
  const std::string& title      ,
  RooAbsReal&        xvar       ,
  const RooArgList&  pars       ,
  const bool         increasing ,
  const bool         convex     ,
  const double       xmin       , 
  const double       xmax       ,
  RooAbsReal&        a          ,
  RooAbsReal&        b          )
  : RooAbsReal ( name.c_str  () , title.c_str () )
  , m_xvar      ( "x"    , "Dependent"  , this , xvar ) 
  , m_a         ( "a"    , "shift/bias" , this , a    ) 
  , m_b         ( "b"    , "scale"      , this , b    ) 
  , m_pars      ( "pars" , "Parameters" , this ) 
  , m_convex    ( ::size ( pars ) , xmin , xmax , increasing , convex ) 
{
  //
  Ostap::Assert ( 0 <= ::size ( pars ) , s_EMPTYPARS  , s_v3 , 510 ) ;
  //
  ::copy_real   ( pars , m_pars , s_INVALIDPAR , s_v3 ) ;
  //
  Ostap::Assert ( m_convex.npars() == ::size ( m_pars ), s_INVALIDPARS , s_v3 , 512 ) ;
  //
}
// ============================================================================
// constructor from the variable, range and list of coefficients
// ============================================================================
Ostap::MoreRooFit::Convex::Convex
( const std::string& name       ,
  const std::string& title      ,
  RooAbsReal&        xvar       ,
  const RooArgList&  pars       ,
  const bool         increasing ,
  const bool         convex     ,
  const double       xmin       , 
  const double       xmax       ,
  const double       a          , 
  const double       b          )
  : Convex ( name , title , xvar , pars , increasing , convex , xmin , xmax , 
             RooFit::RooConst ( a ) , 
             RooFit::RooConst ( b ) )
{}
// =============================================================================
// copy constructor 
// =============================================================================
Ostap::MoreRooFit::Convex::Convex
( const Ostap::MoreRooFit::Convex& right , 
  const char*                         name  )
  : RooAbsReal  ( right , name ) 
  , m_xvar      ( "x"    , this , right.m_xvar ) 
  , m_a         ( "a"    , this , right.m_a    ) 
  , m_b         ( "b"    , this , right.m_b    ) 
  , m_pars      ( "pars" , this , right.m_pars ) 
  , m_convex    ( right.m_convex ) 
{}
// =============================================================================
// default constructor 
// ============================================================================
Ostap::MoreRooFit::Convex::Convex (){}
// =============================================================================
// destructor  
// ============================================================================
Ostap::MoreRooFit::Convex::~Convex (){}
// ============================================================================
// clone it!
// ============================================================================
Ostap::MoreRooFit::Convex*
Ostap::MoreRooFit::Convex::clone ( const char* newname ) const
{ return new Convex ( *this , newname ) ; }
// ============================================================================
void Ostap::MoreRooFit::Convex::setPars() const 
{ ::set_pars (  m_pars , m_convex ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::Convex::evaluate  () const 
{
  setPars () ;
  //
  const double x = m_xvar ;
  if ( x < m_convex.xmin () || x > m_convex.xmax() ) { return 0 ; }
  const double a = m_a    ;
  const double b = m_b    ;
  //
  return a + b * m_convex ( x ) ;
}
// ============================================================================
Int_t    Ostap::MoreRooFit::Convex::getAnalyticalIntegral
( RooArgSet&  allVars   , 
  RooArgSet&  analVars  , 
  const char* /* rangeName */ ) const 
{
  return matchArgs  ( allVars , analVars , m_xvar ) ? 1 : 0 ;
}
// ===========================================================================
Double_t Ostap::MoreRooFit::Convex::analyticalIntegral 
( Int_t code            , 
  const char* rangeName ) const 
{
  setPars() ;
  //
  const double xmin = std::max ( m_xvar.min ( rangeName ) , m_convex.xmin () ) ;
  const double xmax = std::min ( m_xvar.max ( rangeName ) , m_convex.xmax () ) ;
  //
  const double a = m_a    ;
  const double b = m_b    ;
  //
  return a * ( xmax - xmin ) + b * m_convex.integral ( xmin , xmax ) ;
}
// ============================================================================



// ============================================================================
// constructor from the variable, range and list of coefficients
// ============================================================================
Ostap::MoreRooFit::ConvexOnly::ConvexOnly
( const std::string& name       ,
  const std::string& title      ,
  RooAbsReal&        xvar       ,
  const RooArgList&  pars       ,
  const bool         convex     ,
  const double       xmin       , 
  const double       xmax       ,
  RooAbsReal&        a          ,
  RooAbsReal&        b          )
  : RooAbsReal ( name.c_str  () , title.c_str () )
  , m_xvar      ( "x"    , "Dependent"  , this , xvar ) 
  , m_a         ( "a"    , "shift/bias" , this , a    ) 
  , m_b         ( "b"    , "scale"      , this , b    ) 
  , m_pars      ( "pars" , "Parameters" , this ) 
  , m_convex    ( ::size ( pars ) , xmin , xmax , convex ) 
{
  //
  Ostap::Assert ( 0 <= :: size ( pars ) , s_EMPTYPARS  , s_v4 , 510 ) ;
  //
  ::copy_real   ( pars , m_pars , s_INVALIDPAR , s_v4 ) ;
  //
  Ostap::Assert ( m_convex.npars() == ::size ( m_pars ), s_INVALIDPARS , s_v4 , 512 ) ;
  //
}
// ============================================================================
// constructor from the variable, range and list of coefficients
// ============================================================================
Ostap::MoreRooFit::ConvexOnly::ConvexOnly
( const std::string& name       ,
  const std::string& title      ,
  RooAbsReal&        xvar       ,
  const RooArgList&  pars       ,
  const bool         convex     ,
  const double       xmin       , 
  const double       xmax       ,
  const double       a          ,
  const double       b          )
  : ConvexOnly ( name , title , xvar , pars , convex , xmin , xmax ,
                 RooFit::RooConst ( a ) , 
                 RooFit::RooConst ( b ) )
{}
// =============================================================================
// copy constructor 
// =============================================================================
Ostap::MoreRooFit::ConvexOnly::ConvexOnly
( const Ostap::MoreRooFit::ConvexOnly& right , 
  const char*                         name  )
  : RooAbsReal  ( right , name ) 
  , m_xvar      ( "x"    , this , right.m_xvar ) 
  , m_a         ( "a"    , this , right.m_a    ) 
  , m_b         ( "b"    , this , right.m_b    ) 
  , m_pars      ( "pars" , this , right.m_pars ) 
  , m_convex    ( right.m_convex ) 
{}
// =============================================================================
// default constructor 
// ============================================================================
Ostap::MoreRooFit::ConvexOnly::ConvexOnly (){}
// =============================================================================
// destructor  
// ============================================================================
Ostap::MoreRooFit::ConvexOnly::~ConvexOnly (){}
// ============================================================================
// clone it!
// ============================================================================
Ostap::MoreRooFit::ConvexOnly*
Ostap::MoreRooFit::ConvexOnly::clone ( const char* newname ) const
{ return new ConvexOnly ( *this , newname ) ; }
// ============================================================================
void Ostap::MoreRooFit::ConvexOnly::setPars() const 
{ ::set_pars (  m_pars , m_convex ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::ConvexOnly::evaluate  () const 
{
  setPars () ;
  //
  const double x = m_xvar ;
  if ( x < m_convex.xmin () || x > m_convex.xmax() ) { return 0 ; }
  const double a = m_a    ;
  const double b = m_b    ;
  //
  return a + b * m_convex ( x ) ;
}
// ============================================================================
Int_t    Ostap::MoreRooFit::ConvexOnly::getAnalyticalIntegral
( RooArgSet&  allVars   , 
  RooArgSet&  analVars  , 
  const char* /* rangeName */ ) const 
{
  return matchArgs  ( allVars , analVars , m_xvar ) ? 1 : 0 ;
}
// ===========================================================================
Double_t Ostap::MoreRooFit::ConvexOnly::analyticalIntegral 
( Int_t code            , 
  const char* rangeName ) const 
{
  setPars() ;
  //
  const double xmin = std::max ( m_xvar.min ( rangeName ) , m_convex.xmin () ) ;
  const double xmax = std::min ( m_xvar.max ( rangeName ) , m_convex.xmax () ) ;
  //
  const double a = m_a    ;
  const double b = m_b    ;
  //
  return a * ( xmax - xmin ) + b * m_convex.integral ( xmin , xmax ) ;
}












// ============================================================================
// constructor from the variable, range and list of coefficients
// ============================================================================
Ostap::MoreRooFit::BSpline::BSpline
( const std::string&         name  ,
  const std::string&         title ,
  RooAbsReal&                xvar  ,
  const std::vector<double>& knots ,
  const RooArgList&          pars  )
  : RooAbsReal  ( name.c_str  () , title.c_str () )
  , m_xvar      ( "!x"    , "observable" , this , xvar ) 
  , m_pars      ( "!pars" , "parameters" , this )
  , m_bspline   ( knots , std::vector<double>( ::size( pars ) , 0 ) )
{
  //
  Ostap::Assert ( 1 <= ::size ( pars ) , s_EMPTYPARS  , s_v5 , 510 ) ;
  ::copy_real   ( pars , m_pars , s_INVALIDPAR , s_v1 ) ;
  Ostap::Assert ( m_bspline.npars () == ::size ( m_pars ), s_INVALIDPARS , s_v5 , 512 ) ;
  //
}
// =============================================================================
// copy constructor 
// =============================================================================
Ostap::MoreRooFit::BSpline::BSpline
( const Ostap::MoreRooFit::BSpline& right , 
  const char*                         name  )
  : RooAbsReal  ( right , name ) 
  , m_xvar      ( "!x"    , this , right.m_xvar ) 
  , m_pars      ( "!pars" , this , right.m_pars ) 
  , m_bspline   ( right.m_bspline ) 
{}
// =============================================================================
// default constructor 
// ============================================================================
Ostap::MoreRooFit::BSpline::BSpline(){}
// =============================================================================
// destructor  
// ============================================================================
Ostap::MoreRooFit::BSpline::~BSpline(){}
// ============================================================================
// clone it!
// ============================================================================
Ostap::MoreRooFit::BSpline*
Ostap::MoreRooFit::BSpline::clone ( const char* newname ) const
{ return new BSpline( *this , newname ) ; }
// ============================================================================
void Ostap::MoreRooFit::BSpline::setPars() const 
{ ::set_pars ( m_pars , m_bspline ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::BSpline::evaluate  () const 
{
  setPars () ;
  const double x = m_xvar ;
  return m_bspline ( x ) ;
}
// ============================================================================
Int_t    Ostap::MoreRooFit::BSpline::getAnalyticalIntegral
( RooArgSet&  allVars   , 
  RooArgSet&  analVars  , 
  const char* /* rangeName */ ) const 
{
  return matchArgs  ( allVars , analVars , m_xvar ) ? 1 : 0 ;
}
// ===========================================================================
Double_t Ostap::MoreRooFit::BSpline::analyticalIntegral 
( Int_t code            , 
  const char* rangeName ) const 
{
  setPars() ;
  const double xmin = m_xvar.min ( rangeName ) ;
  const double xmax = m_xvar.max ( rangeName ) ;
  return m_bspline.integral ( xmin , xmax ) ;
}
// ============================================================================




// ============================================================================
/*  Helper method to check if recursive fractions were 
 *  used for creation of RooAddPdf object
 *  @see RooAddPdf 
 */
// ============================================================================
bool Ostap::MoreRooFit::recursive ( const RooAddPdf& pdf ) 
{
  std::unique_ptr<::FakeAddPdf> fake { new ::FakeAddPdf ( pdf ) } ;
  return fake->recursive() ; 
}
// ============================================================================
/*  get the original fractions from the <code>RooAddPdf</code>
 *  @see RooAddPdf
 */
// ============================================================================
RooArgList Ostap::MoreRooFit::fractions
( const RooAddPdf& pdf       , 
  bool&            recursive )
{
  std::unique_ptr<::FakeAddPdf> fake { new ::FakeAddPdf ( pdf ) } ;
  return fake->fractions ( recursive ) ; 
}
// ============================================================================
/*  get the original fractions from the <code>RooAddPdf</code>
 *  @see RooAddPdf
 */
// ============================================================================
RooArgList Ostap::MoreRooFit::fractions
( const RooAddPdf& pdf ) 
{
  bool recursive ;
  return fractions ( pdf , recursive ) ;
}
// ============================================================================
/* get x-observable
 *  @see RooGauissian
 */
// ============================================================================
const RooAbsReal& Ostap::MoreRooFit::getX ( const RooGaussian& pdf ) 
{
#if ROOT_VERSION(6,26,0)<=ROOT_VERSION_CODE
  return pdf.getX() ;
#else 
  std::unique_ptr<::FakeGaussian> fake { new ::FakeGaussian ( pdf ) } ;
  return fake->get_x() ;
#endif 
}
// ============================================================================
/*  get mean value 
 *  @see RooGauissian
 */
// ============================================================================
const RooAbsReal& Ostap::MoreRooFit::getMean ( const RooGaussian& pdf ) 
{
#if ROOT_VERSION(6,26,0)<=ROOT_VERSION_CODE
  return pdf.getMean() ;
#else 
  std::unique_ptr<::FakeGaussian> fake { new ::FakeGaussian ( pdf ) } ;
  return fake->get_mean() ;
#endif 
}
// ============================================================================
/*  get sigma 
 *  @see RooGauissian
 */
// ============================================================================
const RooAbsReal& Ostap::MoreRooFit::getSigma ( const RooGaussian& pdf ) 
{
#if ROOT_VERSION(6,26,0)<=ROOT_VERSION_CODE
  return pdf.getSigma() ;
#else 
  std::unique_ptr<::FakeGaussian> fake { new ::FakeGaussian ( pdf ) } ;
  return fake->get_sigma() ;
#endif 
}
// ============================================================================
/*  get parameters from RooFFTConvPdf 
 *  @see RooFFTConvPdf 
 */
// ============================================================================
RooArgList 
Ostap::MoreRooFit::fft_pars 
( const RooFFTConvPdf& pdf    , 
  double&              shift1 ,
  double&              shift2 ) 
{
  std::unique_ptr<::FakeFFTConvPdf> fake { new ::FakeFFTConvPdf ( pdf , "QUQU") } ;
  return fake->getpars ( shift1 , shift2 ) ; 
}
// ============================================================================
/*  get the efficiency function from the RooEfficiency object
 *  @see RooEfficiency
 */
// ============================================================================
const RooAbsReal&    
Ostap::MoreRooFit::get_eff ( const RooEfficiency& pdf ) 
{
  std::unique_ptr<::FakeEfficiency> fake { new ::FakeEfficiency ( pdf ) } ;
  return fake->get_eff() ;
}
// ============================================================================
/*  get the category from the RooEfficiency object
 *  @see RooEfficiency
 */
// ============================================================================
const RooAbsCategory& 
Ostap::MoreRooFit::get_cat ( const RooEfficiency& pdf ) 
{
  std::unique_ptr<::FakeEfficiency> fake { new ::FakeEfficiency ( pdf ) } ;
  return fake->get_cat() ;
}
// ============================================================================
/** get the name of the 'accept' category from RooEfficiency object
 *  @see RooEfficiency
 */
// ============================================================================
std::string          
Ostap::MoreRooFit::get_acc ( const RooEfficiency& pdf ) 
{
  std::unique_ptr<::FakeEfficiency> fake { new ::FakeEfficiency ( pdf ) } ;
  return fake->get_acc() ;  
}
// ============================================================================
/*  get the coefficiencts from the <code>RooPolyVar</code>
 *  @see RooPolyVar 
 */
// ============================================================================
RooArgList Ostap::MoreRooFit::coefficients 
( const RooPolyVar& var      )
{
  std::unique_ptr<::FakePolyVar> fake { new ::FakePolyVar( var ) } ;
  return fake->coefficients () ;  
}
// ============================================================================
/*  get the coefficiencts from the <code>RooPolynomial</code>
 *  @see RooPolynomial 
 */
// ============================================================================
RooArgList Ostap::MoreRooFit::coefficients 
( const RooPolynomial& var      )
{
  std::unique_ptr<::FakePolynomial> fake { new ::FakePolynomial( var ) } ;
  return fake->coefficients () ;  
}
// ============================================================================
/*  get the observables from <code>RooMultiVarGaussian</code>
 *  @see RooMultiVarGaussian
 */
// ============================================================================
RooArgList Ostap::MoreRooFit::observables 
( const RooMultiVarGaussian& pdf ) 
{
  std::unique_ptr<::FakeMultiVarGaussian> fake { new ::FakeMultiVarGaussian ( pdf ) } ;
  return fake->observables () ;  
}
// ============================================================================
/*  get vectro of mu-values from <code>RooMultiVarGaussian</code>
 *  @see RooMultiVarGaussian
 */
// ============================================================================
TVectorD Ostap::MoreRooFit::mu_vec
( const RooMultiVarGaussian& pdf ) 
{
  std::unique_ptr<::FakeMultiVarGaussian> fake { new ::FakeMultiVarGaussian ( pdf ) } ;
  return fake->mu_vec () ;  
}
// ========================================================================




// ============================================================================
//                                                                      The END 
// ============================================================================



