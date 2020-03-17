// ============================================================================
// Include files 
// ============================================================================
// ROOT/RooFit 
// ============================================================================
#include "RooAbsReal.h"
#include "RooArgList.h"
#include "RooAddition.h"
#include "RooGlobalFunc.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/MoreMath.h"
#include "Ostap/MoreRooFit.h"
#include "Ostap/Power.h"
// ============================================================================
// local
// ============================================================================
#include "Exception.h"
#include "local_math.h"
// ============================================================================
/** @file
 *  implementaton of various small additions to RooFit 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru 
 *  @date 2019-11-21
 */
// ============================================================================
ClassImp(Ostap::MoreRooFit::Addition      )
ClassImp(Ostap::MoreRooFit::Subtraction   )
ClassImp(Ostap::MoreRooFit::Division      )
ClassImp(Ostap::MoreRooFit::Fraction      )
ClassImp(Ostap::MoreRooFit::Asymmetry     )
ClassImp(Ostap::MoreRooFit::Power         )
ClassImp(Ostap::MoreRooFit::Exp           )
ClassImp(Ostap::MoreRooFit::Log           )
ClassImp(Ostap::MoreRooFit::Erf           )
ClassImp(Ostap::MoreRooFit::Sin           )
ClassImp(Ostap::MoreRooFit::Cos           )
ClassImp(Ostap::MoreRooFit::Tan           )
ClassImp(Ostap::MoreRooFit::Tanh          )
ClassImp(Ostap::MoreRooFit::Atan2         )
ClassImp(Ostap::MoreRooFit::Gamma         )
ClassImp(Ostap::MoreRooFit::LGamma        )
ClassImp(Ostap::MoreRooFit::IGamma        )
ClassImp(Ostap::MoreRooFit::Id            )
// ============================================================================
namespace 
{
  //
  const std::string s_division_by_zero      { "division by zero"       } ;
  const std::string s_negative_log          { "log of negativenumber"  } ;
  //
  inline std::string  name_  ( const std::string&  name , 
                               const std::string&  oper , 
                               const TNamed&       a    ,
                               const TNamed&       b    )
  { 
    return 
      name.empty () ? ( oper + "_" + a.GetName() + "_" + b.GetName() ) : name ; 
  }  
  // ==========================================================================
  inline std::string  title_ ( const std::string& title , 
                               const std::string& oper  , 
                               const TNamed&      a    ,
                               const TNamed&      b    )
  { return 
      title.empty () ? std::string ( "(" ) + a.GetName() + oper + b.GetName() + ")" : title  ; }
  // ==========================================================================
  inline std::string  title1_ ( const std::string& title , 
                                const std::string& oper  , 
                                const TNamed&      a    ,
                                const TNamed&      b    )
  { return 
      title.empty () ? oper + "(" + a.GetName() + "," + b.GetName() + ")" : title  ; }
  // ==========================================================================
  inline std::string  name2_ ( const std::string&  name , 
                               const std::string&  oper , 
                               const TNamed&       b    )
  { 
    return 
      name.empty () ? ( oper + "_" + b.GetName() ) : name ; 
  }  
  // ==========================================================================
  inline std::string  title2_ ( const std::string& title , 
                               const std::string&  oper  , 
                               const TNamed&        b    )
  { return 
      title.empty () ? oper + "(" + b.GetName() + ")" : title  ; }
  // ==========================================================================
}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Addition::Addition
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     )
  : RooAddition ( name_  ( name  , "sum" , a , b ).c_str() , 
                  title_ ( title , "+"   , a , b ).c_str() , 
                  RooArgList  ( a , b ) )
{}
// ============================================================================
//  copy constructor 
// ============================================================================
Ostap::MoreRooFit::Addition::Addition
( const Ostap::MoreRooFit::Addition& right   , 
  const char*                        newname ) 
  : RooAddition ( right , newname )
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::MoreRooFit::Addition::~Addition(){}
// ============================================================================
Ostap::MoreRooFit::Addition* 
Ostap::MoreRooFit::Addition::clone ( const char* newname ) const 
{ return new Addition ( *this , newname ) ; }
// ============================================================================



// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Product::Product
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     )
  : RooProduct 
    ( name_  ( name  , "mult" , a , b ).c_str() , 
      title_ ( title , "*"    , a , b ).c_str() , RooArgList ( a , b ) ) 
{}
// ============================================================================
//  copy constructor 
// ============================================================================
Ostap::MoreRooFit::Product::Product
( const Ostap::MoreRooFit::Product& right   , 
  const char*                        newname ) 
  : RooProduct ( right , newname )
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::MoreRooFit::Product::~Product(){}
// ============================================================================
Ostap::MoreRooFit::Product* 
Ostap::MoreRooFit::Product::clone ( const char* newname ) const 
{ return new Product ( *this , newname ) ; }
// ============================================================================



// ============================================================================
// constructor with two variables
// ============================================================================
Ostap::MoreRooFit::Subtraction::Subtraction
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     )
  : Addition 
    ( name_  ( name  , "subtract" , a , b )  ,
      title_ ( title , "-"        , a , b )  , a , b )
{}
// ============================================================================
Ostap::MoreRooFit::Subtraction::Subtraction
( const Ostap::MoreRooFit::Subtraction& right   , 
  const char*                           newname ) 
  : Addition ( right , newname )
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::MoreRooFit::Subtraction::~Subtraction(){}
// ============================================================================
Ostap::MoreRooFit::Subtraction* 
Ostap::MoreRooFit::Subtraction::clone ( const char* newname ) const 
{ return new Subtraction ( *this , newname ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::Subtraction::analyticalIntegral 
( Int_t code            , 
  const char* rangeName ) const 
{
  // note: rangeName implicit encoded in code: see _cacheMgr.setObj in getPartIntList...
  CacheElem *cache = (CacheElem*) _cacheMgr.getObjByIndex(code-1);
  if (cache==0) {
    // cache got sterilized, trigger repopulation of this slot, then try again...
    std::unique_ptr<RooArgSet> vars( getParameters(RooArgSet()) );  
    std::unique_ptr<RooArgSet> iset(  _cacheMgr.nameSet2ByIndex(code-1)->select(*vars) );
    RooArgSet dummy;
    Int_t code2 = getAnalyticalIntegral(*iset,dummy,rangeName);
    assert(code==code2); // must have revived the right (sterilized) slot...
    return analyticalIntegral(code2,rangeName);
  }
  assert(cache!=0);
  
#if ROOT_VERSION_CODE <= ROOT_VERSION(6,18,0)
  {
    // loop over cache, and sum...
    std::unique_ptr<TIterator> iter( cache->_I.createIterator() );
    RooAbsReal *I;
    double result(0);
    bool   first = true  ;
    while ( ( I=(RooAbsReal*)iter->Next() ) != 0 ) 
    { 
      if ( first )  { result += I->getVal() ; first = false ; }
      else          { result -= I->getVal() ; }
      
    }
    return result;
  }
#else
  {
    // loop over cache, and sum...
    double result = 0   ;
    bool   first  = true ;
    for (auto I : cache->_I) 
    {
      //
      const double tmp =  static_cast<const RooAbsReal*>(I)->getVal();
      if ( first ) { result += tmp ; first = false ; }
      else         { result -= tmp ;                 }
      //
    }
    return result;
  }
#endif 
}
// ============================================================================
Double_t Ostap::MoreRooFit::Subtraction::evaluate() const
{
#if ROOT_VERSION_CODE <= ROOT_VERSION(6,18,0)
  {
    Double_t sum(0);
    bool first = true ;
    const RooArgSet* nset = _set.nset() ;
    RooFIter setIter = _set.fwdIterator() ;
    RooAbsReal* comp ;
    while((comp=(RooAbsReal*)setIter.next())) 
    {
      Double_t tmp = comp->getVal(nset) ;
      if ( first ) { sum += tmp ; first = false ; }
      else         { sum -= tmp ; }
    }
    return sum ;
  }
#else 
  {
    Double_t result = 0 ;
    const RooArgSet* nset = _set.nset() ;
    //
    bool first = true ;
    for ( const auto arg : _set) 
    {
      const auto comp = static_cast<RooAbsReal*>(arg);
      const Double_t tmp = comp->getVal(nset);
      if ( first ) { result += tmp ;  first = false ; }
      else         { result -= tmp ;                  }
    }
  return result ;
  }
#endif 
}
// ============================================================================


// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Division::Division 
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : RooAbsReal 
    ( name_  ( name  , "divide" , a , b ).c_str() ,
      title_ ( title , "/"      , a , b ).c_str() )
  , m_A ( "!A" , "A" , this , a ) 
  , m_B ( "!B" , "B" , this , b ) 
{}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::MoreRooFit::Division::Division 
( const Ostap::MoreRooFit::Division& right , 
  const char*               name  ) 
  : RooAbsReal ( right , name ) 
  , m_A ( "!A" , this , right.m_A ) 
  , m_B ( "!B" , this , right.m_B ) 
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::MoreRooFit::Division::~Division(){}
// ============================================================================
// cloning
// ============================================================================
Ostap::MoreRooFit::Division* 
Ostap::MoreRooFit::Division::clone ( const char* newname ) const 
{ return new Division ( *this , newname ) ; }
// ============================================================================
// the actual evaluation of the result 
// ============================================================================
Double_t Ostap::MoreRooFit::Division::evaluate () const 
{
  const double a = m_A ;
  const double b = m_B ;
  //
  static const std::string s_self { "Ostap::MoreRooFit::Division" } ;
  Ostap::Assert ( a + b != 0 , s_division_by_zero , s_self ) ;
  //
  return a / b ;
}
// ============================================================================



// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Fraction::Fraction
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : Division 
    ( name_  ( name  , "divide" , a , b ) , 
      title_ ( title , "/"      , a , b ) , a , b )
{}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::MoreRooFit::Fraction::Fraction
( const Ostap::MoreRooFit::Fraction& right , 
  const char*                        name  ) 
  : Division ( right , name ) 
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::MoreRooFit::Fraction::~Fraction(){}
// ============================================================================
// cloning
// ============================================================================
Ostap::MoreRooFit::Fraction* 
Ostap::MoreRooFit::Fraction::clone ( const char* newname ) const 
{ return new Fraction ( *this , newname ) ; }
// ============================================================================
// the actual evaluation of the result 
// ============================================================================
Double_t Ostap::MoreRooFit::Fraction::evaluate () const 
{
  const double a = m_A ;
  const double b = m_B ;
  //
  static const std::string s_self { "Ostap::MoreRooFit::Fraction" } ;
  Ostap::Assert ( a + b != 0 , s_division_by_zero , s_self ) ;
  //
  return a / ( a + b ) ;
}
// ============================================================================

 
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Asymmetry::Asymmetry 
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : Division 
    ( name_  ( name  , "asymmetry" , a , b ) , 
      title_ ( title , " asym "    , a , b ) , a , b )
{}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::MoreRooFit::Asymmetry::Asymmetry
( const Ostap::MoreRooFit::Asymmetry& right , 
  const char*                            name  ) 
  : Division ( right , name ) 
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::MoreRooFit::Asymmetry::~Asymmetry(){}
// ============================================================================
// cloning
// ============================================================================
Ostap::MoreRooFit::Asymmetry* 
Ostap::MoreRooFit::Asymmetry::clone ( const char* newname ) const 
{ return new Asymmetry ( *this , newname ) ; }
// ============================================================================
// the actual evaluation of the result 
// ============================================================================
Double_t Ostap::MoreRooFit::Asymmetry::evaluate () const 
{
  const double a = m_A ;
  const double b = m_B ;
  //
  static const std::string s_self { "Ostap::MoreRooFit::Asymmetry" } ;
  Ostap::Assert ( a + b != 0 , s_division_by_zero , s_self ) ;
  //
  return ( a - b ) / ( a + b ) ;
}


// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Power::Power
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : Division 
    ( name_  ( name  , "pow" , a , b ) , 
      title_ ( title , "** " , a , b ) , a , b )
{}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::MoreRooFit::Power::Power
( const Ostap::MoreRooFit::Power& right , 
  const char*               name  ) 
  : Division  ( right , name ) {}
// ============================================================================
// destructor 
// ============================================================================
Ostap::MoreRooFit::Power::~Power(){}
// ============================================================================
// cloning
// ============================================================================
Ostap::MoreRooFit::Power* 
Ostap::MoreRooFit::Power::clone ( const char* newname ) const 
{ return new Power ( *this , newname ) ; }
// ============================================================================
// the actual evaluation of the result 
// ============================================================================
Double_t Ostap::MoreRooFit::Power::evaluate () const 
{
  const double a = m_A ;
  const double b = m_B ;
  //
  const long double aL = a ;
  const long double bL = b ;
  //
  if      ( s_zero ( b )          ) { return 1.0 ; }        // RETURN 
  else if ( 0 < b && s_zero ( a ) ) { return 0.0 ; }        // RETURN
  else if ( 0 < b && Ostap::Math::isint ( b ) )
  {
    const int nb = Ostap::Math::round   ( b ) ;
    if      ( 0 == nb )             { return 1.0 ; }        // RETURN  
    else if ( 0 <  nb ) 
    {
      const unsigned long NB = nb ;
      return 0 == NB ? 1.0L : Ostap::Math::POW ( aL , NB ) ;  //   RETURN
    } 
  }
  //
  return std::pow ( aL , bL ) ;
}
// ============================================================================

// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Exp::Exp
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : Division 
    ( name_   ( name , "exp" , a , b ) , 
      title1_ ( name , "exp" , a , b ) , a , b )
{}
// ============================================================================
// cloning
// ============================================================================
Ostap::MoreRooFit::Exp* 
Ostap::MoreRooFit::Exp::clone ( const char* newname ) const 
{ return new Exp ( *this , newname ) ; }
// ============================================================================
// the actual evaluation of the result 
// ============================================================================
Double_t Ostap::MoreRooFit::Exp::evaluate () const 
{
  const long double a = m_A ;
  const long double b = m_B ;
  //
  return s_zero ( a ) || s_zero ( b ) ? 0 : std::exp ( a * b ) ;
}
// ============================================================================

// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Log::Log
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : Division 
    ( name_   ( name , "log" , a , b ) , 
      title1_ ( name , "log" , a , b ) , a , b )
{}
// ============================================================================
// cloning
// ============================================================================
Ostap::MoreRooFit::Log* 
Ostap::MoreRooFit::Log::clone ( const char* newname ) const 
{ return new Log ( *this , newname ) ; }
// ============================================================================
// the actual evaluation of the result 
// ============================================================================
Double_t Ostap::MoreRooFit::Log::evaluate () const 
{
  const long double a  = m_A ;
  const long double b  = m_B ;
  //
  const long double ab = a * b ;
  
  static const std::string s_self { "Ostap::MoreRooFit::Log" } ;
  Ostap::Assert ( ab > 0 , s_negative_log , s_self ) ;
  //
  return std::log ( ab ) ;
}
// ============================================================================


// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Erf::Erf
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : Division 
    ( name_   ( name , "erf" , a , b ) , 
      title1_ ( name , "erf" , a , b ) , a , b )
{}
// ============================================================================
// cloning
// ============================================================================
Ostap::MoreRooFit::Erf* 
Ostap::MoreRooFit::Erf::clone ( const char* newname ) const 
{ return new Erf ( *this , newname ) ; }
// ============================================================================
// the actual evaluation of the result 
// ============================================================================
Double_t Ostap::MoreRooFit::Erf::evaluate () const 
{
  const long double a = m_A ;
  const long double b = m_B ;
  //
  return std::erf ( a * b ) ;
}
// ============================================================================


// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Sin::Sin
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : Division 
    ( name_   ( name , "sin" , a , b ) , 
      title1_ ( name , "sin" , a , b ) , a , b )
{}
// ============================================================================
// cloning
// ============================================================================
Ostap::MoreRooFit::Sin* 
Ostap::MoreRooFit::Sin::clone ( const char* newname ) const 
{ return new Sin ( *this , newname ) ; }
// ============================================================================
// the actual evaluation of the result 
// ============================================================================
Double_t Ostap::MoreRooFit::Sin::evaluate () const 
{
  const long double a  = m_A ;
  const long double b  = m_B ;
  //
  const long double ab = a * b ;  
  //
  return std::sin( ab ) ;
}
// ============================================================================

// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Cos::Cos
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : Division 
    ( name_   ( name , "cos" , a , b ) , 
      title1_ ( name , "cos" , a , b ) , a , b )
{}
// ============================================================================
// cloning
// ============================================================================
Ostap::MoreRooFit::Cos* 
Ostap::MoreRooFit::Cos::clone ( const char* newname ) const 
{ return new Cos ( *this , newname ) ; }
// ============================================================================
// the actual evaluation of the result 
// ============================================================================
Double_t Ostap::MoreRooFit::Cos::evaluate () const 
{
  const long double a  = m_A ;
  const long double b  = m_B ;
  //
  const long double ab = a * b ;  
  //
  return std::cos ( ab ) ;
}
// ============================================================================



// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Tan::Tan
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : Division 
    ( name_   ( name , "tan" , a , b ) , 
      title1_ ( name , "tan" , a , b ) , a , b )
{}
// ============================================================================
// cloning
// ============================================================================
Ostap::MoreRooFit::Tan* 
Ostap::MoreRooFit::Tan::clone ( const char* newname ) const 
{ return new Tan ( *this , newname ) ; }
// ============================================================================
// the actual evaluation of the result 
// ============================================================================
Double_t Ostap::MoreRooFit::Tan::evaluate () const 
{
  const long double a  = m_A ;
  const long double b  = m_B ;
  //
  const long double ab = a * b ;  
  //
  return std::tan( ab ) ;
}
// ============================================================================


// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Tanh::Tanh
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : Division 
    ( name_   ( name , "tanh" , a , b ) , 
      title1_ ( name , "tanh" , a , b ) , a , b )
{}
// ============================================================================
// cloning
// ============================================================================
Ostap::MoreRooFit::Tanh* 
Ostap::MoreRooFit::Tanh::clone ( const char* newname ) const 
{ return new Tanh ( *this , newname ) ; }
// ============================================================================
// the actual evaluation of the result 
// ============================================================================
Double_t Ostap::MoreRooFit::Tanh::evaluate () const 
{
  const long double a  = m_A ;
  const long double b  = m_B ;
  //
  const long double ab = a * b ;  
  //
  return std::tanh( ab ) ;
}
// ============================================================================


// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Atan2::Atan2
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : Division 
    ( name_   ( name , "atan2" , a , b ) , 
      title1_ ( name , "atan2" , a , b ) , a , b )
{}
// ============================================================================
// cloning
// ============================================================================
Ostap::MoreRooFit::Atan2* 
Ostap::MoreRooFit::Atan2::clone ( const char* newname ) const 
{ return new Atan2 ( *this , newname ) ; }
// ============================================================================
// the actual evaluation of the result 
// ============================================================================
Double_t Ostap::MoreRooFit::Atan2::evaluate () const 
{
  const long double a  = m_A ;
  const long double b  = m_B ;
  //
  return std::atan2 ( a , b ) ;
}
// ============================================================================



// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Gamma::Gamma
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : Division 
    ( name_   ( name , "gamma" , a , b ) , 
      title1_ ( name , "gamma" , a , b ) , a , b )
{}
// ============================================================================
// cloning
// ============================================================================
Ostap::MoreRooFit::Gamma* 
Ostap::MoreRooFit::Gamma::clone ( const char* newname ) const 
{ return new Gamma ( *this , newname ) ; }
// ============================================================================
// the actual evaluation of the result 
// ============================================================================
Double_t Ostap::MoreRooFit::Gamma::evaluate () const 
{
  const long double a = m_A ;
  const long double b = m_B ;
  //
  return std::tgamma ( a * b ) ;
}
// ============================================================================

// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::LGamma::LGamma
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : Division 
    ( name_   ( name , "lgamma" , a , b ) , 
      title1_ ( name , "lgamma" , a , b ) , a , b )
{}
// ============================================================================
// cloning
// ============================================================================
Ostap::MoreRooFit::LGamma* 
Ostap::MoreRooFit::LGamma::clone ( const char* newname ) const 
{ return new LGamma ( *this , newname ) ; }
// ============================================================================
// the actual evaluation of the result 
// ============================================================================
Double_t Ostap::MoreRooFit::LGamma::evaluate () const 
{
  const long double a  = m_A ;
  const long double b  = m_B ;
  //
  return std::lgamma ( a * b ) ;
}
// ============================================================================

// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::IGamma::IGamma
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : Division 
    ( name_   ( name , "igamma" , a , b ) , 
      title1_ ( name , "igamma" , a , b ) , a , b )
{}
// ============================================================================
// cloning
// ============================================================================
Ostap::MoreRooFit::IGamma* 
Ostap::MoreRooFit::IGamma::clone ( const char* newname ) const 
{ return new IGamma ( *this , newname ) ; }
// ============================================================================
// the actual evaluation of the result 
// ============================================================================
Double_t Ostap::MoreRooFit::IGamma::evaluate () const 
{
  const long double a  = m_A ;
  const long double b  = m_B ;
  //
  return Ostap::Math::igamma ( a * b ) ;
}
// ============================================================================

// ============================================================================
// constructor with variable
// ============================================================================
Ostap::MoreRooFit::Id::Id
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        v     ) 
  : RooAbsReal 
    ( name2_  ( name  , "Id" , v ).c_str() ,
      title2_ ( title , "Id" , v ).c_str() )
  , m_V ( "!V" , "V" , this , v ) 
{}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::MoreRooFit::Id::Id
( const Ostap::MoreRooFit::Id& right , 
  const char*               name  ) 
  : RooAbsReal ( right , name ) 
  , m_V ( "!V" , this , right.m_V ) 
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::MoreRooFit::Id::~Id(){}
// ============================================================================
// cloning
// ============================================================================
Ostap::MoreRooFit::Id* 
Ostap::MoreRooFit::Id::clone ( const char* newname ) const 
{ return new Id( *this , newname ) ; }
// ============================================================================
// the actual evaluation of the result 
// ============================================================================
Double_t Ostap::MoreRooFit::Id::evaluate () const 
{
  const double v = m_V ;
  //
  return v ;
}
// ============================================================================


// ============================================================================
//                                                                      The END
// ============================================================================


