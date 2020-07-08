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
ClassImp(Ostap::MoreRooFit::Abs           )
ClassImp(Ostap::MoreRooFit::Exp           )
ClassImp(Ostap::MoreRooFit::Log           )
ClassImp(Ostap::MoreRooFit::Log10         )
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
ClassImp(Ostap::MoreRooFit::OneVar        )
ClassImp(Ostap::MoreRooFit::TwoVars       )
ClassImp(Ostap::MoreRooFit::FunOneVar     )
ClassImp(Ostap::MoreRooFit::FunTwoVars    )
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
  inline bool a_zero ( const double x ) { return s_zero ( x ) ; }
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
// constructor with variable
// ============================================================================
Ostap::MoreRooFit::OneVar::OneVar
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        x     ) 
  : RooAbsReal 
    ( name2_  ( name  , "one" , x ).c_str() ,
      title2_ ( title , "one" , x ).c_str() )
  , m_x ( "!x" , "x" , this , x ) 
{}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::MoreRooFit::OneVar::OneVar
( const Ostap::MoreRooFit::OneVar& right , 
  const char*               name  ) 
  : RooAbsReal ( right , name ) 
  , m_x ( "!x" , this , right.m_x ) 
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::MoreRooFit::OneVar::~OneVar(){}
// ============================================================================


// ============================================================================
// constructor with variable
// ============================================================================
Ostap::MoreRooFit::TwoVars::TwoVars
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        x     ,
  RooAbsReal&        y     ) 
  : OneVar ( name_   ( name  , "two" , x , y ) ,
             title1_ ( title , "two" , x , y ) , x )
  , m_y ( "!y" , "y" , this , y ) 
{}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::MoreRooFit::TwoVars::TwoVars
( const Ostap::MoreRooFit::TwoVars& right , 
  const char*                       name  ) 
  : OneVar ( right , name ) 
  , m_y    ( "!y" , this , right.m_y ) 
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::MoreRooFit::TwoVars::~TwoVars(){}
// ============================================================================


// ============================================================================
// copy constructor 
// ============================================================================
Ostap::MoreRooFit::FunOneVar::FunOneVar
( const Ostap::MoreRooFit::FunOneVar& right , 
  const char*                          name  ) 
  : OneVar ( right , name ) 
  , m_fun  ( right.m_fun ) 
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::MoreRooFit::FunOneVar::~FunOneVar () {}
// ============================================================================
// clone method 
// ============================================================================
Ostap::MoreRooFit::FunOneVar*
Ostap::MoreRooFit::FunOneVar::clone ( const char* newname ) const 
{ return new Ostap::MoreRooFit::FunOneVar(*this, newname ) ;}
// ============================================================================
// the actual evaluation of the result 
// ============================================================================
Double_t Ostap::MoreRooFit::FunOneVar::evaluate () const
{
  const long double x = m_x ;
  //
  return m_fun ( x ) ;
} ;
// ============================================================================



// ============================================================================
// copy constructor 
// ============================================================================
Ostap::MoreRooFit::FunTwoVars::FunTwoVars
( const Ostap::MoreRooFit::FunTwoVars& right , 
  const char*                          name  ) 
  : TwoVars ( right , name ) 
  , m_fun2  ( right.m_fun2 ) 
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::MoreRooFit::FunTwoVars::~FunTwoVars(){}
// ============================================================================
// the actual evaluation of the result 
// ============================================================================
Double_t Ostap::MoreRooFit::FunTwoVars::evaluate () const
{
  const double a = m_x ;
  const double b = m_y ;
  //
  return m_fun2 ( a , b ) ;
} ;
// ============================================================================
Ostap::MoreRooFit::FunTwoVars*
Ostap::MoreRooFit::FunTwoVars::clone ( const char* newname ) const 
{ return new Ostap::MoreRooFit::FunTwoVars( *this , newname ) ; }
// ============================================================================

// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Division::Division 
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : FunTwoVars ( name_  ( name  , "divide" , a , b ) ,
                 title_ ( title , "/"      , a , b ) , 
                 []( const double x , const double y ) -> double
                 { return x / y ; } ,
                 a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Fraction::Fraction
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : FunTwoVars ( name_  ( name  , "fraction" , a , b ) ,
                 title_ ( title , "frac"     , a , b ) , 
                 []( const double x , const double y ) -> double
                 { return x / ( x + y ) ; } ,
                 a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Asymmetry::Asymmetry
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : FunTwoVars ( name_  ( name  , "asymmetry" , a , b ) ,
                 title_ ( title , "asym"      , a , b ) , 
                 []( const double x , const double y ) -> double
                 { return ( x  - y ) / ( x + y ) ; } ,
                 a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Power::Power
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : FunTwoVars ( name_  ( name  , "pow" , a , b ) ,
                 title_ ( title , "**"  , a , b ) , 
                 []( const double x , const double y ) -> double
                 { 
                   if      ( a_zero ( y )                      ) { return 1.0 ; }
                   else if ( 0 < y && a_zero ( x )             ) { return 0.0 ; }
                   else if ( 0 < y && Ostap::Math::isint ( y ) ) 
                   { 
                     const int ny = Ostap::Math::round ( y ) ;
                     if      ( 0 == ny ) { return 1.0 ; }
                     return std::pow ( x , ny ) ;
                   }
                   return std::pow ( x , y ) ; } ,
                 a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Abs::Abs
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : FunTwoVars ( name_   ( name  , "abs" , a , b ) ,
                 title1_ ( title , "abs" , a , b ) , 
                 []( const double x , const double y ) -> double
                 { return std::abs ( x * y ) ; } ,
                 a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Exp::Exp
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : FunTwoVars ( name_   ( name  , "exp" , a , b ) ,
                 title1_ ( title , "exp" , a , b ) , 
                 []( const double x , const double y ) -> double
                 { return std::exp ( x * y ) ; } ,
                 a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Log::Log
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : FunTwoVars ( name_   ( name  , "log" , a , b ) ,
                 title1_ ( title , "log" , a , b ) , 
                 []( const double x , const double y ) -> double
                 { return std::log ( x * y ) ; } ,
                 a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Log10::Log10
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : FunTwoVars ( name_   ( name  , "log10" , a , b ) ,
                 title1_ ( title , "log10" , a , b ) , 
                 []( const double x , const double y ) -> double
                 { return std::log10 ( x * y ) ; } ,
                 a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Erf::Erf
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : FunTwoVars ( name_   ( name  , "erf" , a , b ) ,
                 title1_ ( title , "erf" , a , b ) , 
                 []( const double x , const double y ) -> double
                 { return std::erf ( x * y ) ; } ,
                 a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Sin::Sin
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : FunTwoVars ( name_   ( name  , "sin" , a , b ) ,
                 title1_ ( title , "sin" , a , b ) , 
                 []( const double x , const double y ) -> double
                 { return std::sin ( x * y ) ; } ,
                 a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Cos::Cos
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : FunTwoVars ( name_   ( name  , "cos" , a , b ) ,
                 title1_ ( title , "cos" , a , b ) , 
                 []( const double x , const double y ) -> double
                 { return std::sin ( x * y ) ; } ,
                 a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Tan::Tan
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : FunTwoVars ( name_   ( name  , "tan" , a , b ) ,
                 title1_ ( title , "tan" , a , b ) , 
                 []( const double x , const double y ) -> double
                 { return std::tan ( x * y ) ; } ,
                 a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Tanh::Tanh
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : FunTwoVars ( name_   ( name  , "tanh" , a , b ) ,
                 title1_ ( title , "tanh" , a , b ) , 
                 []( const double x , const double y ) -> double
                 { return std::tanh ( x * y ) ; } ,
                 a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Atan2::Atan2
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : FunTwoVars ( name_   ( name  , "atan2" , a , b ) ,
                 title1_ ( title , "atan2" , a , b ) , 
                 []( const double x , const double y ) -> double
                 { return std::atan2 ( x , y ) ; } ,
                 a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Gamma::Gamma
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : FunTwoVars ( name_   ( name  , "gamma" , a , b ) ,
                 title1_ ( title , "gamma" , a , b ) , 
                 []( const double x , const double y ) -> double
                 { return std::tgamma ( x * y ) ; } ,
                 a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::LGamma::LGamma
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : FunTwoVars ( name_   ( name  , "lgamma" , a , b ) ,
                 title1_ ( title , "lgamma" , a , b ) , 
                 []( const double x , const double y ) -> double
                 { return std::lgamma ( x * y ) ; } ,
                 a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::IGamma::IGamma
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : FunTwoVars ( name_   ( name  , "igamma" , a , b ) ,
                 title1_ ( title , "igamma" , a , b ) , 
                 []( const double x , const double y ) -> double
                 { return Ostap::Math::igamma ( x * y ) ; } ,
                 a , b )
{}
// ============================================================================



// ============================================================================
Ostap::MoreRooFit::Id::Id
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        v     ) 
  : OneVar 
    ( name2_  ( name  , "Id" , v ) ,
      title2_ ( title , "Id" , v ) , v ) 
{}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::MoreRooFit::Id::Id
( const Ostap::MoreRooFit::Id& right , 
  const char*               name  ) 
  : OneVar ( right , name ) 
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
{ return new Id ( *this , newname ) ; }
// ============================================================================
// the actual evaluation of the result 
// ============================================================================
Double_t Ostap::MoreRooFit::Id::evaluate () const 
{
  const double v = m_x ;
  return v ;
}
// ============================================================================
Double_t Ostap::MoreRooFit::Id::analyticalIntegral
( Int_t            code     ,
  const char*      range    ) const
{ return m_x.arg().analyticalIntegral    ( code , range ) ; }
//
Double_t Ostap::MoreRooFit::Id::analyticalIntegralWN
( Int_t            code     ,
  const RooArgSet* normset  ,
  const char*      range    ) const
{ return m_x.arg().analyticalIntegralWN  ( code , normset , range ) ; }
//  
Int_t    Ostap::MoreRooFit::Id::getAnalyticalIntegral
( RooArgSet&       allVars  ,
  RooArgSet&       analVars ,
  const char*      range    ) const
{ return m_x.arg().getAnalyticalIntegral   ( allVars , analVars , range ) ; }
//
Int_t    Ostap::MoreRooFit::Id::getAnalyticalIntegralWN
( RooArgSet&       allVars  ,
  RooArgSet&       analVars ,
  const RooArgSet* normset  ,
  const char*      range    ) const
{ return m_x.arg().getAnalyticalIntegralWN ( allVars , analVars , normset , range ) ; }
// ============================================================================

// ============================================================================
//                                                                      The END
// ============================================================================


