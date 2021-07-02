// ============================================================================
// Include files 
// ============================================================================
// ROOT/RooFit 
// ============================================================================
#include "RooAbsReal.h"
#include "RooArgList.h"
#include "RooAddition.h"
#include "RooAbsPdf.h"
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
ClassImp(Ostap::MoreRooFit::Erfc          )
ClassImp(Ostap::MoreRooFit::Sin           )
ClassImp(Ostap::MoreRooFit::Cos           )
ClassImp(Ostap::MoreRooFit::Tan           )
ClassImp(Ostap::MoreRooFit::Sinh          )
ClassImp(Ostap::MoreRooFit::Cosh          )
ClassImp(Ostap::MoreRooFit::Tanh          )
ClassImp(Ostap::MoreRooFit::Sech          )
ClassImp(Ostap::MoreRooFit::Atan2         )
ClassImp(Ostap::MoreRooFit::Gamma         )
ClassImp(Ostap::MoreRooFit::LGamma        )
ClassImp(Ostap::MoreRooFit::IGamma        )
ClassImp(Ostap::MoreRooFit::Id            )
ClassImp(Ostap::MoreRooFit::OneVar        )
ClassImp(Ostap::MoreRooFit::TwoVars       )
ClassImp(Ostap::MoreRooFit::FunOneVar     )
ClassImp(Ostap::MoreRooFit::FunTwoVars    )
ClassImp(Ostap::MoreRooFit::ProductPdf    )
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
// construct c1*a + c2*b 
// ============================================================================
Ostap::MoreRooFit::Addition::Addition
( const std::string& name  , 
  const std::string& title ,
  RooAbsReal&        a     ,
  RooAbsReal&        b     ,
  RooAbsReal&        c1    ,
  RooAbsReal&        c2    ) 
  : RooAddition ( name .c_str() , title.c_str() ,
                  RooArgList ( a  ,  b  ) , 
                  RooArgList ( c1 ,  c2 ) )
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
      title_ ( title , "-"        , a , b )  , a , b , 1 , -1 )    
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
  : TwoVars ( name_  ( name  , "divide" , a , b ) ,
              title_ ( title , "/"      , a , b ) , a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Fraction::Fraction
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : TwoVars ( name_  ( name  , "fraction" , a , b ) ,
              title_ ( title , "frac"     , a , b ) , a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Asymmetry::Asymmetry
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : TwoVars ( name_  ( name  , "asymmetry" , a , b ) ,
              title_ ( title , "asym"      , a , b ) , a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Power::Power
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : TwoVars ( name_  ( name  , "pow" , a , b ) ,
              title_ ( title , "**"  , a , b ) , a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Abs::Abs
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : TwoVars ( name_   ( name  , "abs" , a , b ) ,
              title1_ ( title , "abs" , a , b ) , a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Exp::Exp
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : TwoVars ( name_   ( name  , "exp" , a , b ) ,
              title1_ ( title , "exp" , a , b ) , a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Log::Log
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : TwoVars ( name_   ( name  , "log" , a , b ) ,
              title1_ ( title , "log" , a , b ) , a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Log10::Log10
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : TwoVars ( name_   ( name  , "log10" , a , b ) ,
              title1_ ( title , "log10" , a , b ) , a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Erf::Erf
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : TwoVars ( name_   ( name  , "erf" , a , b ) ,
              title1_ ( title , "erf" , a , b ) , a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Erfc::Erfc
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : TwoVars ( name_   ( name  , "erfc" , a , b ) ,
              title1_ ( title , "erfc" , a , b ) , a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Sin::Sin
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : TwoVars ( name_   ( name  , "sin" , a , b ) ,
              title1_ ( title , "sin" , a , b ) , a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Cos::Cos
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : TwoVars ( name_   ( name  , "cos" , a , b ) ,
              title1_ ( title , "cos" , a , b ) , a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Tan::Tan
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : TwoVars ( name_   ( name  , "tan" , a , b ) ,
              title1_ ( title , "tan" , a , b ) , a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Sinh::Sinh
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : TwoVars ( name_   ( name  , "sinh" , a , b ) ,
              title1_ ( title , "sinh" , a , b ) , a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Cosh::Cosh
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : TwoVars ( name_   ( name  , "cosh" , a , b ) ,
              title1_ ( title , "cosh" , a , b ) , a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Tanh::Tanh
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : TwoVars ( name_   ( name  , "tanh" , a , b ) ,
              title1_ ( title , "tanh" , a , b ) , a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Sech::Sech
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : TwoVars ( name_   ( name  , "sech" , a , b ) ,
              title1_ ( title , "sech" , a , b ) , a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Atan2::Atan2
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : TwoVars ( name_   ( name  , "atan2" , a , b ) ,
              title1_ ( title , "atan2" , a , b ) , a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Gamma::Gamma
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : TwoVars ( name_   ( name  , "gamma" , a , b ) ,
              title1_ ( title , "gamma" , a , b ) , a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::LGamma::LGamma
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : TwoVars ( name_   ( name  , "lgamma" , a , b ) ,
              title1_ ( title , "lgamma" , a , b ) , a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::IGamma::IGamma
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : TwoVars ( name_   ( name  , "igamma" , a , b ) ,
                 title1_ ( title , "igamma" , a , b ) , a , b )
{}
// ============================================================================

// ============================================================================
// the actual evaluation of the result 
// ============================================================================
Double_t Ostap::MoreRooFit::Division::evaluate () const 
{ const double a = m_x ; const double b = m_y ; return a / b ; }
// ============================================================================
Double_t Ostap::MoreRooFit::Fraction::evaluate () const 
{ const double a = m_x ; const double b = m_y ; return a / ( a + b ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::Asymmetry::evaluate () const 
{ const double a = m_x ; const double b = m_y ; return ( a - b ) / ( a + b ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::Power::evaluate () const 
{ 
  const double x = m_x ; 
  const double y = m_y ;
  if      ( a_zero ( y )                      ) { return 1.0 ; }
  else if ( 0 < y && a_zero ( x )             ) { return 0.0 ; }
  else if ( 0 < y && Ostap::Math::isint ( y ) ) 
  { 
    const int ny = Ostap::Math::round ( y ) ;
    if      ( 0 == ny ) { return 1.0 ; }
    return std::pow ( x , ny ) ;
  }
  return std::pow ( x , y ) ;
}
// ============================================================================
Double_t Ostap::MoreRooFit::Abs::evaluate () const 
{ const double a = m_x ; const double b = m_y ; return std::abs    ( a * b ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::Exp::evaluate () const 
{ const double a = m_x ; const double b = m_y ; return std::exp    ( a * b ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::Log::evaluate () const 
{ const double a = m_x ; const double b = m_y ; return std::log    ( a * b ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::Log10::evaluate () const 
{ const double a = m_x ; const double b = m_y ; return std::log10  ( a * b ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::Erf::evaluate () const 
{ const double a = m_x ; const double b = m_y ; return std::erf    ( a * b ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::Erfc::evaluate () const 
{ const double a = m_x ; const double b = m_y ; return std::erfc   ( a * b ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::Gamma::evaluate () const 
{ const double a = m_x ; const double b = m_y ; return std::tgamma ( a * b ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::LGamma::evaluate () const 
{ const double a = m_x ; const double b = m_y ; return std::lgamma ( a * b ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::IGamma::evaluate () const 
{ const double a = m_x ; const double b = m_y ; 
  return Ostap::Math::igamma ( a * b ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::Sin::evaluate () const 
{ const double a = m_x ; const double b = m_y ; return std::sin    ( a * b ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::Cos::evaluate () const 
{ const double a = m_x ; const double b = m_y ; return std::cos    ( a * b ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::Tan::evaluate () const 
{ const double a = m_x ; const double b = m_y ; return std::tan    ( a * b ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::Sinh::evaluate () const 
{ const double a = m_x ; const double b = m_y ; return std::sinh   ( a * b ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::Cosh::evaluate () const 
{ const double a = m_x ; const double b = m_y ; return std::cosh   ( a * b ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::Tanh::evaluate () const 
{ const double a = m_x ; const double b = m_y ; return std::tanh   ( a * b ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::Sech::evaluate () const 
{ 
  const double a = m_x ; const double b = m_y ; 
  return Ostap::Math::sech ( a * b ) ; 
}
// ============================================================================
Double_t Ostap::MoreRooFit::Atan2::evaluate () const 
{ const double a = m_x ; const double b = m_y ; return std::atan2  ( a , b ) ; }
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
/* constructor from name, title and two pdfs
 *  @param name  name 
 *  @param title name 
 *  @param pdf1 the first pdf 
 *  @param pdf2 the second pdf 
 */
// ============================================================================
Ostap::MoreRooFit::ProductPdf::ProductPdf 
( const char* name  , 
  const char* title , 
  RooAbsPdf&  pdf1  , 
  RooAbsPdf&  pdf2  )
  : RooAbsPdf  ( name , title ) 
    //
  , m_pdf1 ( "pdf1" , "The first PDF"  , this , pdf1 ) 
  , m_pdf2 ( "pdf2" , "The second PDF" , this , pdf2 ) 
{}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::MoreRooFit::ProductPdf::ProductPdf 
( const Ostap::MoreRooFit::ProductPdf& right , 
  const char*                       name  ) 
  : RooAbsPdf ( right , name ) 
    //
  , m_pdf1 ( "pdf1" , this , right.m_pdf1 ) 
  , m_pdf2 ( "pdf2" , this , right.m_pdf2 )
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::MoreRooFit::ProductPdf::~ProductPdf(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::MoreRooFit::ProductPdf*
Ostap::MoreRooFit::ProductPdf::clone ( const char* newname ) const 
{ return new Ostap::MoreRooFit::ProductPdf( *this , newname ) ; }
// ============================================================================
// the main method 
// ============================================================================
Double_t Ostap::MoreRooFit::ProductPdf::evaluate () const
{
  const double v1 = m_pdf1 ;
  const double v2 = m_pdf2 ;
  return v1 * v2 ;
}
// ============================================================================

// ============================================================================
//                                                                      The END
// ============================================================================


