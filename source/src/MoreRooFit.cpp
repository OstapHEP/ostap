// ============================================================================
// Include files 
// ============================================================================
// ROOT
// ============================================================================
#include "TClass.h"
#include "TObject.h"
#include "TNamed.h"
// ============================================================================
// ROOT/RooFit 
// ============================================================================
#include "RooAbsReal.h"
#include "RooFitResult.h"
#include "RooArgList.h"
#include "RooAddition.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooAbsDataStore.h"
#include "RooFitResult.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooAbsCategory.h"
#include "RooAbsRealLValue.h"
#include "RooAbsCategoryLValue.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Iterator.h"
#include "Ostap/Power.h"
#include "Ostap/MoreMath.h"
#include "Ostap/MoreVars.h"
#include "Ostap/MoreRooFit.h"
#include "Ostap/RooFitUtils.h"
#include "Ostap/Peaks.h"
#include "Ostap/ToStream.h"
// ============================================================================
// local
// ============================================================================
#include "Exception.h"
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
ClassImp(Ostap::MoreRooFit::Constant      )
ClassImp(Ostap::MoreRooFit::Addition      )
ClassImp(Ostap::MoreRooFit::Addition2     )
ClassImp(Ostap::MoreRooFit::Subtraction   )
ClassImp(Ostap::MoreRooFit::Division      )
ClassImp(Ostap::MoreRooFit::Combination   )
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
ClassImp(Ostap::MoreRooFit::Sigmoid       )
ClassImp(Ostap::MoreRooFit::Hypot         )
ClassImp(Ostap::MoreRooFit::AbsAplusB     )
ClassImp(Ostap::MoreRooFit::BesselJ       )
ClassImp(Ostap::MoreRooFit::BesselY       )
ClassImp(Ostap::MoreRooFit::BesselI       )
ClassImp(Ostap::MoreRooFit::BesselK       )
ClassImp(Ostap::MoreRooFit::Gamma         )
ClassImp(Ostap::MoreRooFit::LGamma        )
ClassImp(Ostap::MoreRooFit::IGamma        )
ClassImp(Ostap::MoreRooFit::Id            )
ClassImp(Ostap::MoreRooFit::MaxV          )
ClassImp(Ostap::MoreRooFit::MinV          )
ClassImp(Ostap::MoreRooFit::OneVar        )
ClassImp(Ostap::MoreRooFit::TwoVars       )
ClassImp(Ostap::MoreRooFit::NVars         )
ClassImp(Ostap::MoreRooFit::FunOneVar     )
ClassImp(Ostap::MoreRooFit::FunTwoVars    )
ClassImp(Ostap::MoreRooFit::ProductPdf    )
ClassImp(Ostap::MoreRooFit::WrapPdf       )
ClassImp(Ostap::MoreRooFit::AddDeps       )
ClassImp(Ostap::MoreRooFit::NVars         )
ClassImp(Ostap::MoreRooFit::Minimal       )
ClassImp(Ostap::MoreRooFit::Maximal       )
ClassImp(Ostap::MoreRooFit::Rank          )
ClassImp(Ostap::MoreRooFit::ABC           )
ClassImp(Ostap::MoreRooFit::Clamp         )
ClassImp(Ostap::MoreRooFit::CrystalBallN  )
// ============================================================================
namespace 
{
  //
  const std::string s_division_by_zero      { "division by zero"       } ;
  const std::string s_negative_log          { "log of negativenumber"  } ;
  //
  inline std::string  name_ 
  ( const std::string&  name , 
    const std::string&  oper , 
    const TNamed&       a    ,
    const TNamed&       b    )
  { 
    return 
      name.empty () ? ( oper + "_" + a.GetName() + "_" + b.GetName() ) : name ; 
  }  
  // ==========================================================================
  inline std::string  title_ 
  ( const std::string& title , 
    const std::string& oper  , 
    const TNamed&      a    ,
    const TNamed&      b    )
  { return 
      title.empty () ? std::string ( "(" ) + a.GetName() + oper + b.GetName() + ")" : title  ; }
  // ==========================================================================
  inline std::string  title1_ 
  ( const std::string& title , 
    const std::string& oper  , 
    const TNamed&      a    ,
    const TNamed&      b    )
  { return 
      title.empty () ? oper + "(" + a.GetName() + "," + b.GetName() + ")" : title  ; }
  // ==========================================================================
  inline std::string  name2_
  ( const std::string&  name , 
    const std::string&  oper , 
    const TNamed&       b    )
  { 
    return 
      name.empty () ? ( oper + "_" + b.GetName() ) : name ; 
  }  
  // ==========================================================================
  inline std::string  title2_ 
  ( const std::string& title , 
    const std::string&  oper  , 
    const TNamed&        b    )
  { return 
      title.empty () ? oper + "(" + b.GetName() + ")" : title  ; }
  // ==========================================================================
  inline bool a_zero ( const double x ) { return s_zero ( x ) ; }
  // ==========================================================================
}

// ============================================================================
// constructor with no variables 
// ============================================================================
Ostap::MoreRooFit::Constant::Constant
( const std::string& name  , 
  const std::string& title , 
  const double       value ) 
  : RooAbsReal ( name.c_str () , title.c_str() ) 
  , m_value    ( value ) 
  , m_vlst     ( "!vlst" , "variables" , this ) 
{
  // _fast = true ;
  setAttribute("Constant",true) ;
}
// ============================================================================
// constructor with one variables 
// ============================================================================
Ostap::MoreRooFit::Constant::Constant
( const std::string& name  , 
  const std::string& title , 
  const double       value ,
  RooAbsReal&        x     )
  : Constant ( name , title , value ) 
{
  m_vlst.add ( x ) ;
}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Constant::Constant
( const std::string& name  , 
  const std::string& title , 
  const double       value ,
  RooAbsReal&        x     ,
  RooAbsReal&        y     )
  : Constant ( name , title , value  ) 
{
  m_vlst.add ( x ) ;
  m_vlst.add ( y ) ;
}
// ============================================================================
// constructor with three variables 
// ============================================================================
Ostap::MoreRooFit::Constant::Constant
( const std::string& name  , 
  const std::string& title , 
  const double       value ,
  RooAbsReal&        x     ,
  RooAbsReal&        y     ,
  RooAbsReal&        z     )
  : Constant ( name , title , value ) 
{
  m_vlst.add ( x ) ;
  m_vlst.add ( y ) ;
  m_vlst.add ( z ) ;
}
// ============================================================================
// constructor with list of variables 
// ============================================================================
Ostap::MoreRooFit::Constant::Constant
( const std::string& name  , 
  const std::string& title , 
  const double       value ,
  const RooArgList&  v     )
  : Constant ( name , title , value ) 
{
  ::copy_real ( v , m_vlst , "Invalid var parameter" , "Ostap::MoreRooFit::Constant!" ) ;  
}
// ============================================================================
// copy 
// ============================================================================
Ostap::MoreRooFit::Constant::Constant
( const Ostap::MoreRooFit::Constant& right   , 
  const char*                        newname ) 
  : RooAbsReal ( right , newname ) 
  , m_value    ( right.m_value   )
  , m_vlst     ( "!vlst" , this , right.m_vlst )
{
  // _fast = true ;
}
// ============================================================================
Ostap::MoreRooFit::Constant*
Ostap::MoreRooFit::Constant::clone ( const char* newname ) const 
{ return new Ostap::MoreRooFit::Constant(*this , newname ) ; }
// ============================================================================
// REDEFINED METHODS 
// ============================================================================
#if ROOT_VERSION_CODE<=ROOT_VERSION(6,29,0)
#include "RunContext.h"
RooSpan<const double> 
Ostap::MoreRooFit::Constant::getValues 
( RooBatchCompute::RunContext& evalData   , 
  const                        RooArgSet* ) const 
{
  if ( evalData.spans.end() == evalData.spans.find ( this ) )
  { evalData.spans[this] = {&m_value, 1}; }
  return evalData.spans[this];
}
#endif
// ============================================================================
// write to the stream 
// ============================================================================
void Ostap::MoreRooFit::Constant::writeToStream(std::ostream& os, bool compact) const
{ os << m_value ; }
// ============================================================================
// the actual evaluation of the result 
// ================m============================================================
Double_t Ostap::MoreRooFit::Constant::evaluate () const
{ return m_value ; } 
// ============================================================================
// the actual evaluation of the result
// ============================================================================
double   Ostap::MoreRooFit::Constant::getValV 
( const RooArgSet* /* a */ ) const 
{ return m_value ; }



// ============================================================================
// construct x + y 
// ============================================================================
Ostap::MoreRooFit::Addition::Addition
( const std::string& name  , 
  const std::string& title ,
  RooAbsReal&        x     ,
  RooAbsReal&        y     ) 
  : RooAddition ( name_  ( name  , "add" , x , y ).c_str() , 
                  title_ ( title , "+"   , x , y ).c_str() ,
                  RooArgList ( x  , y  ) ) 
{}
// ============================================================================
// construct x + y + z 
// ============================================================================
Ostap::MoreRooFit::Addition::Addition
( const std::string& name  , 
  const std::string& title ,
  RooAbsReal&        x     ,
  RooAbsReal&        y     ,
  RooAbsReal&        z     ) 
  : RooAddition ( name_  ( name  , "add" , x , z ).c_str() , 
                  title_ ( title , "+"   , x , z ).c_str() ,
                  RooArgList ( x  , y  , z ) ) 
{}
// ============================================================================
// several variables 
// ============================================================================
Ostap::MoreRooFit::Addition::Addition
( const std::string& name  ,  
  const std::string& title ,
  const RooArgList&  vars  )
  : RooAddition ( name.c_str() , title.c_str() , vars )
{
  Ostap::Assert ( 2 <= size() ,
		  "Not enoght elements!"
		  "Ostap::MoreRooFit::Addition" ) ;
}

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
// construct c1*x + c2*y 
// ============================================================================
Ostap::MoreRooFit::Addition2::Addition2
( const std::string& name  , 
  const std::string& title ,
  RooAbsReal&        x     ,
  RooAbsReal&        y     , 
  RooAbsReal&        c1    ,
  RooAbsReal&        c2    )
  : RooAddition ( name_  ( name  , "add" , x , y ).c_str() , 
                  title_ ( title , "+"   , x , y ).c_str() ,
                  RooArgList ( x  , y  ) , 
                  RooArgList ( c1 , c2 ) )
  , m_a ( "!a" , "variables" , this )
  , m_c ( "!c" , "variables" , this )    
{
  m_a.add ( x  ) ;
  m_a.add ( y  ) ;
  m_c.add ( c1 ) ;
  m_c.add ( c2 ) ;
}
// ============================================================================
/// construct \f$ \sum_i a_ic_i\f$ 
Ostap::MoreRooFit::Addition2::Addition2
( const std::string& name  ,  
  const std::string& title ,
  const RooArgList&  a     ,
  const RooArgList&  c     )
  : RooAddition ( name.c_str() , title.c_str() , a , c )
  , m_a ( "!a" , "variables" , this )
  , m_c ( "!c" , "variables" , this )    
{
  ::copy_real ( a , m_a , "Invalid var parameter!" , "Ostap::MoreRooFit::Addition2" ) ;
  ::copy_real ( c , m_c , "Invalid var parameter!" , "Ostap::MoreRooFit::Addition2" ) ;
  //
  Ostap::Assert ( 2 <= m_a.getSize()
		  &&   m_a.getSize() ==  m_c.getSize()
		  &&   m_a.getSize() == _set.getSize() ,
		  "Invalid size of variables" "Ostap::MoreRooFit::Addition2" ) ;
}
// ============================================================================
//  copy constructor 
// ============================================================================
Ostap::MoreRooFit::Addition2::Addition2
( const Ostap::MoreRooFit::Addition2& right   , 
  const char*                         newname ) 
  : RooAddition ( right , newname )
  , m_a ( "!a" , this , right.m_a )
  , m_c ( "!a" , this , right.m_c )
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::MoreRooFit::Addition2::~Addition2(){}
// ============================================================================
Ostap::MoreRooFit::Addition2* 
Ostap::MoreRooFit::Addition2::clone ( const char* newname ) const 
{ return new Addition2 ( *this , newname ) ; }
// ============================================================================



// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Product::Product
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        x     , 
  RooAbsReal&        y     )
  : RooProduct 
    ( name_  ( name  , "mult" , x , y ).c_str() , 
      title_ ( title , "*"    , x , y ).c_str() , 
      RooArgList ( x , y ) ) 
{}
// ============================================================================
Ostap::MoreRooFit::Product::Product
( const std::string& name  , 
  const std::string& title ,
  const RooArgList&  vars  ) 
  : RooProduct ( name.c_str () , title.c_str() , vars )
{
  Ostap::Assert ( 2 <= vars.getSize ()
		  && size () == vars.getSize()
		  &&  0 == _compCSet.getSize()  , 
		  "Invalid size of varibles"    ,
		  "Ostap::MoreRooFit::Product"  ) ; 
}
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
  : Addition2 
    ( name_  ( name  , "subtract" , a , b )  ,
      title_ ( title , "-"        , a , b )  , a , b , 1 , -1 )    
{}
// ============================================================================
Ostap::MoreRooFit::Subtraction::Subtraction
( const Ostap::MoreRooFit::Subtraction& right   , 
  const char*                           newname ) 
  : Addition2 ( right , newname )
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
// constructor with variable
// ============================================================================
Ostap::MoreRooFit::OneVar::OneVar
( const std::string& name  , 
  const std::string& title , 
  const double       x     )
  : OneVar ( name , title , RooFit::RooConst ( x ) )
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
// constructor with two variables
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
// constructor with two variables
// ============================================================================
Ostap::MoreRooFit::TwoVars::TwoVars
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        x     ,
  const double       y     ) 
  : TwoVars ( name , title , x , RooFit::RooConst ( y ) )
{}
// ============================================================================
// constructor with two variables
// ============================================================================
Ostap::MoreRooFit::TwoVars::TwoVars
( const std::string& name  , 
  const std::string& title , 
  const double       x     ,
  RooAbsReal&        y     ) 
  : TwoVars ( name , title , RooFit::RooConst ( x ) , y )
{}
// ============================================================================
// constructor with two variables
// ============================================================================
Ostap::MoreRooFit::TwoVars::TwoVars
( const std::string& name  , 
  const std::string& title , 
  const double       x     ,
  const double       y     )
  : TwoVars ( name , title , RooFit::RooConst ( x ) , RooFit::RooConst ( y ) )	      
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
// constructor with wtwo variables 
// ============================================================================
Ostap::MoreRooFit::NVars::NVars
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a1    ,
  RooAbsReal&        a2    ) 
  : RooAbsReal ( name.c_str()  , title.c_str() ) 
  , m_vars     ( "!vars"  , "variables" , this ) 
{
  m_vars.add ( a1 ) ;
  m_vars.add ( a2 ) ;  
}  
// ============================================================================
// constructor with three variables 
// ============================================================================
Ostap::MoreRooFit::NVars::NVars
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a1    ,
  RooAbsReal&        a2    ,
  RooAbsReal&        a3    ) 
  : RooAbsReal ( name.c_str()  , title.c_str() ) 
  , m_vars     ( "!vars"  , "variables" , this ) 
{
  m_vars.add ( a1 ) ;
  m_vars.add ( a2 ) ;  
  m_vars.add ( a3 ) ;  
}  
// ============================================================================
// constructor with four variables 
// ============================================================================
Ostap::MoreRooFit::NVars::NVars
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a1    ,
  RooAbsReal&        a2    ,
  RooAbsReal&        a3    , 
  RooAbsReal&        a4    ) 
  : RooAbsReal ( name.c_str()  , title.c_str() ) 
  , m_vars     ( "!vars"  , "variables" , this ) 
{
  m_vars.add ( a1 ) ;
  m_vars.add ( a2 ) ;  
  m_vars.add ( a3 ) ;  
  m_vars.add ( a4 ) ;  
}  
// ============================================================================
// constructor with five variables 
// ============================================================================
Ostap::MoreRooFit::NVars::NVars
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a1    ,
  RooAbsReal&        a2    ,
  RooAbsReal&        a3    , 
  RooAbsReal&        a4    ,
  RooAbsReal&        a5    ) 
  : RooAbsReal ( name.c_str()  , title.c_str() ) 
  , m_vars     ( "!vars"  , "variables" , this ) 
{
  m_vars.add ( a1 ) ;
  m_vars.add ( a2 ) ;  
  m_vars.add ( a3 ) ;  
  m_vars.add ( a4 ) ;  
  m_vars.add ( a5 ) ;  
}
// ============================================================================
// constructor with list of variables 
// ============================================================================
Ostap::MoreRooFit::NVars::NVars
( const std::string& name  , 
  const std::string& title , 
  const RooArgList&  lst   )
  : RooAbsReal ( name.c_str()  , title.c_str() ) 
  , m_vars     ( "!vars"  , "variables" , this ) 
{
  ::copy_real   ( lst , m_vars, "Invalid var parameter!" ,
                  "Ostap::MoreRooFit::NVars" ) ;
}
// ============================================================================
// constructor with list of variables 
// ============================================================================
Ostap::MoreRooFit::NVars::NVars
( const std::string& name  , 
  const std::string& title , 
  const RooArgSet&   lst   )
  : RooAbsReal ( name.c_str()  , title.c_str() ) 
  , m_vars     ( "!vars"  , "variables" , this ) 
{
  ::copy_real   ( lst , m_vars, "Invalid var parameter!" ,
                  "Ostap::MoreRooFit::NVars" ) ;
}
// ============================================================================
// constructor with list of variables 
// ============================================================================
Ostap::MoreRooFit::NVars::NVars
( const std::string&      name  , 
  const std::string&      title , 
  const RooAbsCollection& lst   )
  : RooAbsReal ( name.c_str()  , title.c_str() ) 
  , m_vars     ( "!vars"  , "variables" , this ) 
{
  ::copy_real   ( lst , m_vars, "Invalid var parameter!" ,
                  "Ostap::MoreRooFit::NVars" ) ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::MoreRooFit::NVars::NVars
( const Ostap::MoreRooFit::NVars& right ,
  const char*                     name  ) 
  : RooAbsReal ( right , name )
  , m_vars     ( "!vars"  , this , right.m_vars )
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::MoreRooFit::NVars::~NVars(){}
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
Ostap::MoreRooFit::Combination::Combination
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     , 
  const double       alpha , 
  const double       beta  , 
  const double       gamma )
  : TwoVars ( name_  ( name  , "combination" , a , b ) ,
              title_ ( title , "comb"        , a , b ) , a , b )
  , m_alpha ( alpha ) 
  , m_beta  ( beta  ) 
  , m_gamma ( gamma )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Asymmetry::Asymmetry
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ,
  const double       scale ) 
  : TwoVars ( name_  ( name  , "asymmetry" , a , b ) ,
              title_ ( title , "asym"      , a , b ) , a , b )
  , m_scale ( scale ) 
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
Ostap::MoreRooFit::Sigmoid::Sigmoid
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : TwoVars ( name_   ( name  , "sigmoid" , a , b ) ,
              title1_ ( title , "sigmoid" , a , b ) , a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Hypot::Hypot
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : TwoVars ( name_   ( name  , "hypot" , a , b ) ,
              title1_ ( title , "hypot" , a , b ) , a , b )
{}

// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::AbsAplusB::AbsAplusB
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : TwoVars ( !name.empty  () ? name  : ( std::string("abs(") + a.GetName() + ")_plus_" + b.GetName() ) ,
	      !title.empty () ? title : ( std::string("abs(") + a.GetName() + " + "     + b.GetName() ) ,
	      a , b ) 
{}




// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::BesselJ::BesselJ
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : TwoVars ( name_   ( name  , "besselJ" , a , b ) ,
              title1_ ( title , "besselJ" , a , b ) , a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::BesselY::BesselY
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : TwoVars ( name_   ( name  , "besselY" , a , b ) ,
              title1_ ( title , "besselY" , a , b ) , a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::BesselI::BesselI
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : TwoVars ( name_   ( name  , "besselI" , a , b ) ,
              title1_ ( title , "besselI" , a , b ) , a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::BesselK::BesselK
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : TwoVars ( name_   ( name  , "besselK" , a , b ) ,
              title1_ ( title , "besselK" , a , b ) , a , b )
{}

// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::MaxV::MaxV
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : TwoVars ( name_   ( name  , "max" , a , b ) ,
              title1_ ( title , "max" , a , b ) , a , b )
{}
// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::MinV::MinV
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     , 
  RooAbsReal&        b     ) 
  : TwoVars ( name_   ( name  , "min" , a , b ) ,
              title1_ ( title , "min" , a , b ) , a , b )
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
Double_t Ostap::MoreRooFit::Combination::evaluate () const 
{ const double a = m_x ; const double b = m_y ; 
  return m_alpha * a * ( m_beta + m_gamma * b )  ; }
// ============================================================================
Double_t Ostap::MoreRooFit::Asymmetry::evaluate () const 
{ const double a = m_x ; const double b = m_y ; return m_scale * ( a - b ) / ( a + b ) ; }
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
Double_t Ostap::MoreRooFit::Sigmoid::evaluate () const 
{ const double a = m_x ; const double b = m_y ; return 0.5 * ( 1 + std::tanh ( a * b ) ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::Hypot::evaluate () const 
{ const double a = m_x ; const double b = m_y ; return std::hypot ( a , b ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::MaxV::evaluate () const 
{ const double a = m_x ; const double b = m_y ; return std::max  ( a , b ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::MinV::evaluate () const 
{ const double a = m_x ; const double b = m_y ; return std::min  ( a , b ) ; }
// ============================================================================


// ============================================================================
Double_t Ostap::MoreRooFit::AbsAplusB::evaluate () const 
{ const double a = m_x ; const double b = m_y ; return std::abs ( a ) + b ; }
// ============================================================================


// ============================================================================
Double_t Ostap::MoreRooFit::BesselJ::evaluate () const 
{ const double x = m_x ; const double nu = m_y ; return Ostap::Math::bessel_Jnu ( nu , x ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::BesselY::evaluate () const 
{ const double x = m_x ; const double nu = m_y ; return Ostap::Math::bessel_Ynu ( nu , x ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::BesselI::evaluate () const 
{ const double x = m_x ; const double nu = m_y ; return Ostap::Math::bessel_Inu ( nu , x ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::BesselK::evaluate () const 
{ const double x = m_x ; const double nu = m_y ; return Ostap::Math::bessel_Knu ( nu , x ) ; }
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
Ostap::MoreRooFit::Id::Id
( const std::string& name  , 
  RooAbsReal&        v     ) 
  : OneVar 
    ( name2_  ( name , "Id" , v ) ,
      title2_ ( name , "Id" , v ) , v ) 
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
Double_t Ostap::MoreRooFit::Id::analyticalIntegral
( Int_t            code     ,
  const char*      range    ) const
{ return m_x.arg().analyticalIntegral    ( code , range ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::Id::analyticalIntegralWN
( Int_t            code     ,
  const RooArgSet* normset  ,
  const char*      range    ) const
{ return m_x.arg().analyticalIntegralWN  ( code , normset , range ) ; }
// ============================================================================
Int_t    Ostap::MoreRooFit::Id::getAnalyticalIntegral
( RooArgSet&       allVars  ,
  RooArgSet&       analVars ,
  const char*      range    ) const
{ return m_x.arg().getAnalyticalIntegral   ( allVars , analVars , range ) ; }
// ============================================================================
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
  , m_pdf1 ( "!pdf1" , "The first PDF"  , this , pdf1 ) 
  , m_pdf2 ( "!pdf2" , "The second PDF" , this , pdf2 ) 
{}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::MoreRooFit::ProductPdf::ProductPdf 
( const Ostap::MoreRooFit::ProductPdf& right , 
  const char*                       name  ) 
  : RooAbsPdf ( right , name ) 
    //
  , m_pdf1 ( "!pdf1" , this , right.m_pdf1 ) 
  , m_pdf2 ( "!pdf2" , this , right.m_pdf2 )
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
// constructor 
// ============================================================================
Ostap::MoreRooFit::WrapPdf::WrapPdf 
( const char *name  , 
  const char *title , 
  RooAbsReal& func  ) 
  : RooAbsPdf ( name , title ) 
  , m_func ( "!func" , "a function" , this , func )
{}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::MoreRooFit::WrapPdf::WrapPdf 
( const Ostap::MoreRooFit::WrapPdf& right   , 
  const char*           newname ) 
  : RooAbsPdf ( right , newname ) 
  , m_func ( "!func" , this , right.m_func )
{}
// ============================================================================
/// cloning method 
// ============================================================================
Ostap::MoreRooFit::WrapPdf*
Ostap::MoreRooFit::WrapPdf::clone ( const char* newname ) const 
{ return new Ostap::MoreRooFit::WrapPdf ( *this , newname ) ; }
// ============================================================================
// Analytical Integration handling
// ============================================================================
bool Ostap::MoreRooFit::WrapPdf::forceAnalyticalInt 
( const RooAbsArg& dep ) const 
{ return m_func.arg().forceAnalyticalInt ( dep ); }
// ============================================================================
Int_t Ostap::MoreRooFit::WrapPdf::getAnalyticalIntegralWN 
( RooArgSet& allVars         , 
  RooArgSet& analVars        , 
  const RooArgSet* normSet   ,
  const char*      rangeName ) const 
{
  return m_func.arg().getAnalyticalIntegralWN ( allVars   ,  
                                                analVars  , 
                                                normSet   , 
                                                rangeName ) ;
}
// ============================================================================
Int_t Ostap::MoreRooFit::WrapPdf::getAnalyticalIntegral
( RooArgSet&       allVars   , 
  RooArgSet&       numVars   ,
  const char*      rangeName ) const 
{ return m_func.arg().getAnalyticalIntegral ( allVars , numVars , rangeName ); }
// ============================================================================
// Hints for optimized brute-force sampling
Int_t  Ostap::MoreRooFit::WrapPdf::getMaxVal
( const RooArgSet& vars) const 
{ return m_func.arg().getMaxVal ( vars ) ; }
// ============================================================================
double  Ostap::MoreRooFit::WrapPdf::maxVal(Int_t code) const
{ return m_func.arg().maxVal ( code ) ; }
// ============================================================================
Int_t Ostap::MoreRooFit::WrapPdf::minTrialSamples
( const RooArgSet& arGenObs ) const 
{ return m_func.arg().minTrialSamples ( arGenObs ) ; }
// Plotting and binning hints
bool  Ostap::MoreRooFit::WrapPdf::isBinnedDistribution 
( const RooArgSet& obs     ) const 
{ return m_func.arg().isBinnedDistribution ( obs ) ; }
// ============================================================================
std::list<double>* Ostap::MoreRooFit::WrapPdf::binBoundaries
( RooAbsRealLValue& obs , 
  double            xlo , 
  double            xhi ) const 
{ return m_func.arg().binBoundaries ( obs , xlo , xhi ) ; }
// ============================================================================
std::list<double>* Ostap::MoreRooFit::WrapPdf::plotSamplingHint
( RooAbsRealLValue& obs , 
  double            xlo , 
  double            xhi ) const 
{ return m_func.arg().plotSamplingHint ( obs , xlo , xhi ) ; }
// ============================================================================



// ============================================================================
// constructor with list of variables 
// ============================================================================
Ostap::MoreRooFit::AddDeps::AddDeps
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        x     ,
  const RooArgList&  v     ) 
  : OneVar ( name , title , x ) 
  , m_vlst ( "!vlst" , "variables" , this )
{
  ::copy_real ( v , m_vlst , 
                "Invalid var parameter" , 
                "Ostap::MoreRooFit::AddDeps!" ) ;  
}
// ============================================================================
// copy 
// ============================================================================
Ostap::MoreRooFit::AddDeps::AddDeps
( const AddDeps& right   , 
  const char*    newname ) 
  : OneVar ( right , newname ) 
  , m_vlst ( "!vlst" , this , right.m_vlst )
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::MoreRooFit::AddDeps::~AddDeps(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::MoreRooFit::AddDeps*
Ostap::MoreRooFit::AddDeps::clone ( const char* newname ) const 
{ return new Ostap::MoreRooFit::AddDeps( *this , newname ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::AddDeps::analyticalIntegral
( Int_t            code     ,
  const char*      range    ) const
{ return m_x.arg().analyticalIntegral    ( code , range ) ; }
// ============================================================================
Double_t Ostap::MoreRooFit::AddDeps::analyticalIntegralWN
( Int_t            code     ,
  const RooArgSet* normset  ,
  const char*      range    ) const
{ return m_x.arg().analyticalIntegralWN  ( code , normset , range ) ; }
// ============================================================================
Int_t    Ostap::MoreRooFit::AddDeps::getAnalyticalIntegral
( RooArgSet&       allVars  ,
  RooArgSet&       analVars ,
  const char*      range    ) const
{ return m_x.arg().getAnalyticalIntegral   ( allVars , analVars , range ) ; }
// ============================================================================
Int_t    Ostap::MoreRooFit::AddDeps::getAnalyticalIntegralWN
( RooArgSet&       allVars  ,
  RooArgSet&       analVars ,
  const RooArgSet* normset  ,
  const char*      range    ) const
{ return m_x.arg().getAnalyticalIntegralWN ( allVars , analVars , normset , range ) ; }
// ============================================================================




// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Minimal::Minimal
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a1    ,
  RooAbsReal&        a2    ) 
  : NVars ( name , title , a1 , a2 ) 
{}
// ============================================================================
// constructor with three variables 
// ============================================================================
Ostap::MoreRooFit::Minimal::Minimal
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a1    ,
  RooAbsReal&        a2    , 
  RooAbsReal&        a3    ) 
  : NVars ( name , title , a1 , a2 , a3 ) 
{}
// ============================================================================
// constructor with four variables 
// ============================================================================
Ostap::MoreRooFit::Minimal::Minimal
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a1    ,
  RooAbsReal&        a2    , 
  RooAbsReal&        a3    ,
  RooAbsReal&        a4    ) 
  : NVars ( name , title , a1 , a2 , a3 , a4 ) 
{}
// ============================================================================
// constructor with five variables 
// ============================================================================(
Ostap::MoreRooFit::Minimal::Minimal
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a1    ,
  RooAbsReal&        a2    , 
  RooAbsReal&        a3    ,
  RooAbsReal&        a4    , 
  RooAbsReal&        a5    ) 
  : NVars ( name , title , a1 , a2 , a3 , a4 , a5 ) 
{}
// ============================================================================
// constructor with list of variables 
// ============================================================================
Ostap::MoreRooFit::Minimal::Minimal
( const std::string& name  , 
  const std::string& title , 
  const RooArgList&  lst   )
  : NVars ( name , title , lst ) 
{
  Ostap::Assert ( 2 <= ::size( m_vars )        ,
                  "Invalid var parameter!"     ,
                  "Ostap::MoreRooFit::Minimal" ) ;
}
// ============================================================================
// constructor with list of variables 
// ============================================================================
Ostap::MoreRooFit::Minimal::Minimal
( const std::string& name  , 
  const std::string& title , 
  const RooArgSet&   lst   )
  : NVars ( name , title , lst ) 
{
  Ostap::Assert ( 2 <= ::size( m_vars )        ,
                  "Invalid var parameter!"     ,
                  "Ostap::MoreRooFit::Minimal" ) ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::MoreRooFit::Minimal::Minimal
( const Ostap::MoreRooFit::Minimal& right ,
  const char*                       name  ) 
  : NVars ( right , name ) 
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::MoreRooFit::Minimal::~Minimal(){}
// ============================================================================
// clone method 
// ============================================================================
Ostap::MoreRooFit::Minimal*
Ostap::MoreRooFit::Minimal::clone ( const char* newname ) const
{ return new Ostap::MoreRooFit::Minimal ( *this , newname ) ; }
// ============================================================================
// the actual evaluation of the result 
// ================m===========================================================
Double_t Ostap::MoreRooFit::Minimal::evaluate () const
{
  // ==========================================================================
  double value = s_INFINITY ;
  // ==========================================================================
  static const std::string s_message { "Invalid variable type!"      } ;
  static const std::string s_tag     { "Ostap::MoreRooFit::Minimal!" } ;
  // ==========================================================================
  unsigned ii = 0 ;
  for ( auto* c : this->m_vars ) 
  {
    const RooAbsReal* v = dynamic_cast<RooAbsReal*> ( c ) ;
    Ostap::Assert ( v != nullptr , s_message , s_tag , 510 ) ;
    value = std::min ( value , v->getVal() ) ;
  }
  // 
  return value ; 
} 
// ============================================================================



// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Maximal::Maximal
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a1    ,
  RooAbsReal&        a2    ) 
  : NVars ( name , title , a1 , a2 ) 
{}
// ============================================================================
// constructor with three variables 
// ============================================================================
Ostap::MoreRooFit::Maximal::Maximal
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a1    ,
  RooAbsReal&        a2    , 
  RooAbsReal&        a3    ) 
  : NVars ( name , title , a1 , a2 , a3 ) 
{}
// ============================================================================
// constructor with four variables 
// ============================================================================
Ostap::MoreRooFit::Maximal::Maximal
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a1    ,
  RooAbsReal&        a2    , 
  RooAbsReal&        a3    ,
  RooAbsReal&        a4    ) 
  : NVars ( name , title , a1 , a2 , a3 , a4 ) 
{}
// ============================================================================
// constructor with five variables 
// ============================================================================(
Ostap::MoreRooFit::Maximal::Maximal
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a1    ,
  RooAbsReal&        a2    , 
  RooAbsReal&        a3    ,
  RooAbsReal&        a4    , 
  RooAbsReal&        a5    ) 
  : NVars ( name , title , a1 , a2 , a3 , a4 , a5 ) 
{}
// ============================================================================
// constructor with list of variables 
// ============================================================================
Ostap::MoreRooFit::Maximal::Maximal
( const std::string& name  , 
  const std::string& title , 
  const RooArgList&  lst   )
  : NVars ( name , title , lst ) 
{
  Ostap::Assert ( 2 <= ::size( m_vars )        ,
                  "Invalid var parameter!"     ,
                  "Ostap::MoreRooFit::Minimal" ) ;
}
// ============================================================================
// constructor with list of variables 
// ============================================================================
Ostap::MoreRooFit::Maximal::Maximal
( const std::string& name  , 
  const std::string& title , 
  const RooArgSet&   lst   )
  : NVars ( name , title , lst ) 
{
  Ostap::Assert ( 2 <= ::size( m_vars )        ,
                  "Invalid var parameter!"     ,
                  "Ostap::MoreRooFit::Minimal" ) ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::MoreRooFit::Maximal::Maximal
( const Ostap::MoreRooFit::Maximal& right ,
  const char*                       name  ) 
  : NVars ( right , name ) 
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::MoreRooFit::Maximal::~Maximal (){}
// ============================================================================
// clone method 
// ============================================================================
Ostap::MoreRooFit::Maximal*
Ostap::MoreRooFit::Maximal::clone ( const char* newname ) const
{ return new Ostap::MoreRooFit::Maximal ( *this , newname ) ; }
// ============================================================================
// the actual evaluation of the result 
// ================m===========================================================
Double_t Ostap::MoreRooFit::Maximal::evaluate () const
{
  // ==========================================================================
  double value = - s_INFINITY ;
  // ==========================================================================
  static const std::string s_message { "Invalid variable type!"      } ;
  static const std::string s_tag     { "Ostap::MoreRooFit::Maximal!" } ;
  // ==========================================================================
  unsigned ii = 0 ;
  for ( auto* c : this->m_vars ) 
  {
    const RooAbsReal* v = dynamic_cast<RooAbsReal*> ( c ) ;
    Ostap::Assert ( v != nullptr , s_message , s_tag , 510 ) ;
    value = std::max ( value , v->getVal() ) ;
  }
  // 
  return value ; 
} 
// ============================================================================




// ============================================================================
// constructor with two variables 
// ============================================================================
Ostap::MoreRooFit::Rank::Rank 
( const std::string& name  , 
  const std::string& title ,
  const int          rank  , 
  RooAbsReal&        a1    ,
  RooAbsReal&        a2    )
  : Rank ( name , title , rank , RooArgList ( a1 , a2 ) )
{}
// ============================================================================
// constructor with three variables 
// ============================================================================
Ostap::MoreRooFit::Rank::Rank 
( const std::string& name  , 
  const std::string& title ,
  const int          rank  , 
  RooAbsReal&        a1    ,
  RooAbsReal&        a2    ,
  RooAbsReal&        a3    )
  : Rank ( name , title , rank , RooArgList ( a1 , a2 , a3 ) )
{}
// ============================================================================
// constructor with four variables 
// ============================================================================
Ostap::MoreRooFit::Rank::Rank 
( const std::string& name  , 
  const std::string& title ,
  const int          rank  , 
  RooAbsReal&        a1    ,
  RooAbsReal&        a2    ,
  RooAbsReal&        a3    ,
  RooAbsReal&        a4    )
  : Rank ( name , title , rank , RooArgList ( a1 , a2 , a3 , a4 ) )
{}
// ============================================================================
// constructor with five variables 
// ============================================================================
Ostap::MoreRooFit::Rank::Rank 
( const std::string& name  , 
  const std::string& title ,
  const int          rank  , 
  RooAbsReal&        a1    ,
  RooAbsReal&        a2    ,
  RooAbsReal&        a3    ,
  RooAbsReal&        a4    ,
  RooAbsReal&        a5    )
  : Rank ( name , title , rank , RooArgList ( a1 , a2 , a3 , a4 , a5 ) )
{}
// ============================================================================
// constructor with six variables 
// ============================================================================
Ostap::MoreRooFit::Rank::Rank 
( const std::string& name  , 
  const std::string& title ,
  const int          rank  , 
  RooAbsReal&        a1    ,
  RooAbsReal&        a2    ,
  RooAbsReal&        a3    ,
  RooAbsReal&        a4    ,
  RooAbsReal&        a5    ,
  RooAbsReal&        a6    )
  : Rank ( name , title , rank , RooArgList ( a1 , a2 , a3 , a4 , a5 , a6 ) )
{}
// ============================================================================
// constructor with many variables 
// ============================================================================
Ostap::MoreRooFit::Rank::Rank 
( const std::string&      name  , 
  const std::string&      title ,
  const int               rank  , 
  const RooAbsCollection& vars  )
  : NVars  ( name , title , vars ) 
  , m_rank ( rank )
  , m_aux  () 
{
  //
  const std::string s_INVALIDPAR  = "Invalid parameter!"                ;
  const std::string s_v8          = "Ostap::MoreRooFit::Rank"           ;
  const std::string s_NOTENOUGH   = "Vector of coefficients is short!"  ;
  const std::string s_INVALIDRANK = "Invalid rank!"                     ;
  //
  const std::size_t  NN = size () ;
  //
  ::copy_real   ( vars , m_vars , s_INVALIDPAR , s_v8 ) ;
  Ostap::Assert ( 1 <= NN , s_NOTENOUGH  , s_v8 , 510 ) ;
  //
  if ( m_rank < 0 ) { m_rank += NN ; }
  //
  Ostap::Assert ( 0 <= m_rank && m_rank < NN , s_INVALIDRANK , s_v8 ) ;
  //
  m_aux.resize ( NN ) ;
}
// =============================================================================
// copy constructor 
// =============================================================================
Ostap::MoreRooFit::Rank::Rank
( const Ostap::MoreRooFit::Rank& right , 
  const char*                    name  )
  : NVars  ( right , name   ) 
  , m_rank ( right.m_rank   )
  , m_aux  ( right.m_aux    ) 
{}
// =============================================================================
// destructor  
// ============================================================================
Ostap::MoreRooFit::Rank::~Rank(){}
// ============================================================================
// clone it!
// ============================================================================
Ostap::MoreRooFit::Rank*
Ostap::MoreRooFit::Rank::clone ( const char* newname ) const
{ return new Rank( *this , newname ) ; }
// ============================================================================
// Evaluate it!
// ============================================================================
Double_t 
Ostap::MoreRooFit::Rank::evaluate() const
{
  ::set_pars ( m_vars , m_aux ) ;
  std::sort ( m_aux.begin () , m_aux.end() ) ;
  return m_aux [ m_rank ] ;
}


// ============================================================================
// constructor with three variables 
// ============================================================================
Ostap::MoreRooFit::ABC::ABC
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        a     ,
  RooAbsReal&        b     , 
  RooAbsReal&        c     ) 
  : NVars ( name , title , a , b , c ) 
{}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::MoreRooFit::ABC::ABC
( const Ostap::MoreRooFit::ABC& right ,
  const char*                   name  ) 
  : NVars ( right , name ) 
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::MoreRooFit::ABC::~ABC (){}
// ============================================================================
// clone method 
// ============================================================================
Ostap::MoreRooFit::ABC*
Ostap::MoreRooFit::ABC::clone ( const char* newname ) const
{ return new Ostap::MoreRooFit::ABC ( *this , newname ) ; }
// ============================================================================
// the actual evaluation of the result 
// ================m===========================================================
Double_t Ostap::MoreRooFit::ABC::evaluate () const
{
  //
  Ostap::Assert ( ::size ( this->m_vars ) == 3    , 
                  "Invalid number of parameters!" ,
                  "Ostap::MoreRooFit::ABC"        ) ;
  //
  const RooAbsArg* aa = this->m_vars.at ( 0 ) ;
  const RooAbsArg* ab = this->m_vars.at ( 1 ) ;
  const RooAbsArg* ac = this->m_vars.at ( 2 ) ;
  //
  Ostap::Assert ( ( nullptr != aa ) && 
                  ( nullptr != ab ) && 
                  ( nullptr != ac )        , 
                  "Invalid parameter!"     ,
                  "Ostap::MoreRooFit::ABC" ) ;
  //
  const RooAbsReal* va = static_cast<const RooAbsReal*> ( aa ) ;
  const RooAbsReal* vb = static_cast<const RooAbsReal*> ( ab ) ;
  const RooAbsReal* vc = static_cast<const RooAbsReal*> ( ac ) ;
  //
  Ostap::Assert ( ( nullptr != va ) && 
                  ( nullptr != vb ) && 
                  ( nullptr != vc )        , 
                  "Invalid parameter!"     ,
                  "Ostap::MoreRooFit::ABC" ) ;
  //
  const double a = va -> getVal() ;
  if ( s_zero ( a ) ) { return  0 ; }
  //
  const double b = vb -> getVal() ;
  if ( s_zero ( b ) ) { return  0 ; }
  //
  const double c = vc -> getVal() ;
  //
  return a * ( b / c ) ;
}
// ============================================================================




// ============================================================================
// constructor with one variables 
// ============================================================================
Ostap::MoreRooFit::Clamp::Clamp
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        e     ,
  const double       a     ,
  const double       b     )  
  : OneVar ( name , title , e )
  , m_a ( std::min ( a , b ) ) 
  , m_b ( std::max ( a , b ) ) 
{
  Ostap::Assert ( m_a < m_b && !s_equal ( m_a , m_b ) ,
                  "Invalid parameters!"     ,
                  "Ostap::MoreRooFit::Clamp" ) ;
}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::MoreRooFit::Clamp::Clamp
( const Ostap::MoreRooFit::Clamp& right ,
  const char*                          name  ) 
  : OneVar ( right , name )
  , m_a ( right.m_a )
  , m_b ( right.m_b )    
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::MoreRooFit::Clamp::~Clamp(){}
// ============================================================================
// clone method 
// ============================================================================
Ostap::MoreRooFit::Clamp*
Ostap::MoreRooFit::Clamp::clone ( const char* newname ) const
{ return new Ostap::MoreRooFit::Clamp ( *this , newname ) ; }
// ============================================================================
// the actual evaluation of the result 
// ================m===========================================================
Double_t Ostap::MoreRooFit::Clamp::evaluate () const
{
  const double xx = m_x ;
  return
    xx <= m_a || s_equal ( xx , m_a ) ? m_a : 
    xx >= m_b || s_equal ( xx , m_b ) ? m_b : xx ;
}


// ============================================================================
// constructor with one variables 
// ============================================================================
Ostap::MoreRooFit::CrystalBallN::CrystalBallN
( const std::string& name  , 
  const std::string& title , 
  RooAbsReal&        n     )
  : OneVar ( name , title , n )
{}
// ============================================================================
// constructor with one variables 
// ============================================================================
Ostap::MoreRooFit::CrystalBallN::CrystalBallN
( const std::string& name  , 
  RooAbsReal&        n     )
  : CrystalBallN ( name , "N-parameter: n -> N transformation" , n )
{}
// ============================================================================
// "copy" constructor 
// ============================================================================
Ostap::MoreRooFit::CrystalBallN::CrystalBallN
( const Ostap::MoreRooFit::CrystalBallN& right ,
  const char*                            name  ) 
  : OneVar ( right , name )
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::MoreRooFit::CrystalBallN::~CrystalBallN(){}
// ============================================================================
// clone method 
// ============================================================================
Ostap::MoreRooFit::CrystalBallN*
Ostap::MoreRooFit::CrystalBallN::clone ( const char* newname ) const
{ return new Ostap::MoreRooFit::CrystalBallN ( *this , newname ) ; }
// ============================================================================
// the actual evaluation of the result 
// ================m===========================================================
Double_t Ostap::MoreRooFit::CrystalBallN::evaluate () const
{ return Ostap::Math::CrystalBall::N ( m_x ) ; }
// ============================================================================


// ============================================================================
/* reset RooabsData and ubnderlying store 
 *  @see RooAbsData
 *  @see RooAbsDataStore 
 *  @param data dataset to be reset/clean 
 */
// ============================================================================
RooAbsData*
Ostap::MoreRooFit::reset_data
( RooAbsData* data )
{
  if ( nullptr == data ) { return nullptr ; }
  //
  RooAbsDataStore* store = data->store() ;
  if ( store )
    {
      store -> resetCache   () ;
      store -> resetBuffers () ;
      store -> reset        () ; 
    }
  //
  data -> resetBuffers () ;
  data -> reset        () ;
  //
  return data ;
}
// ============================================================================
/*  delete  RooAbsData
 *  @see RooAbsData
 *  @see RooAbsDataStore 
 *  @param data dataset to be reset/clean 
 *  @return nullptr 
 */
// ============================================================================
RooAbsData*
Ostap::MoreRooFit::delete_data
( RooAbsData* data )
{
  if ( nullptr == data ) { return nullptr ; }
  reset_data ( data ) ;
  delete data ;
  return nullptr ;
}
// ========================================================================
// helper function to call RooAbsPdf::fitTo ( data , options ) 
// ============================================================================
RooFitResult* 
Ostap::MoreRooFit::fitTo
( RooAbsPdf&           model , 
  RooAbsData&          data  , 
  const RooLinkedList& opts  )
{ return model.fitTo ( data , opts ) ; }
// ============================================================================
// helper function to call RooAbsPdf::createNLL( data , options ) 
// ============================================================================
RooAbsReal* 
Ostap::MoreRooFit::createNLL
( RooAbsPdf&           model , 
  RooAbsData&          data  , 
  const RooLinkedList& opts  )
{ return model.createNLL ( data , opts ) ; }
// ============================================================================
/* assign RooAbsCollection
 *  @see RooAbsCollection
 *  @param  from  source collection  
 *  @param  to    destination collection
 */
// ===========================================================================
void Ostap::MoreRooFit::assign 
(       RooAbsCollection& to   ,
  const RooAbsCollection& from ) 
{ 
::assign ( to , from ) ; 
}
// ============================================================================
/*  print TNamed
 *  @see  TNamed
 */
// ============================================================================
std::ostream&
Ostap::Utils::toStream
( const TNamed& o ,
  std::ostream& s )
{
  const TClass* c = o.IsA() ;
  if (   nullptr !=  c          ) { s << c->GetName() ; }
  const char* n = o.GetName() ;
  if ( ( nullptr != n ) && n[0] ) { s << " Name:"  << n ; }
  const char* t = o.GetTitle() ;
  if ( ( nullptr != t ) && t[0] ) { s << " Title:" << t ; }
  return s ;
}
// ============================================================================
/** print RooAbsReal 
 *  @see  RooAbsReal
 */
// ============================================================================
std::ostream&
Ostap::Utils::toStream
( const RooAbsReal& o ,
  std::ostream&     s )
{
  // ==========================================================================
  const unsigned short s_precision = 4 ; // ====================================
  // ==========================================================================
  if ( nullptr != dynamic_cast<const RooConstVar*> ( &o ) )
    { return toStream ( o.getVal() , s , s_precision  ) ; }           // RETURN
  //
  const RooRealVar* rrv = dynamic_cast<const RooRealVar*>  ( &o ) ;
  if ( nullptr == rrv ) { o.printClassName ( s ) ; s << "/"     ; }
  //
  s << o.GetName() << " : " ;
  if ( o.isConstant () ) { return toStream ( o.getVal() , s , s_precision ) ; } // RETURN  
  //
  if ( rrv && rrv->hasAsymError () && ( rrv->getErrorLo() || rrv->getErrorHi() ) ) 
    {
      s << "( " ;
      toStream ( o.getVal()                , s , s_precision ) << " -/ " ;
      toStream ( rrv->getErrorLo ()        , s , s_precision ) << " +/ " ;
      return toStream ( rrv->getErrorHi () , s , s_precision ) << " )" ;  // RETURN 
    }
  //
  if ( rrv && rrv->hasError () && rrv->getError () )
    {
      s << "( " ;
      toStream ( o.getVal()              , s , s_precision ) << " +/- " ;
      return toStream ( rrv->getError () , s , s_precision ) << " )" ;     // RETURN 
    }
  //
  return toStream ( o.getVal() , s , s_precision ) ;                       // RETURN       
}
// ============================================================================
/*  print RooAbsCategory  
 *  @see  RooAbsCategory
 */
// ============================================================================
std::ostream&
Ostap::Utils::toStream
( const RooAbsCategory& o ,
  std::ostream&         s )
{
  // ==========================================================================
  s << o.GetName() << "::" ;
  toStream ( ::getLabel ( o ) , s ) << "/" ;
  return toStream ( ::getValue ( o )  , s ) ;
  // ==========================================================================
}
// ============================================================================
/*  print RooAbsArg 
 *  @see  RooAbsArg 
 */
// ============================================================================
std::ostream&
Ostap::Utils::toStream
( const RooAbsArg& o ,
  std::ostream&    s )
{
  const RooAbsReal*     r = dynamic_cast<const RooAbsReal*>    (&o) ;
  if ( nullptr != r ) { return toStream ( *r , s ) ; }
  const RooAbsCategory* c = dynamic_cast<const RooAbsCategory*>(&o) ;
  if ( nullptr != c ) { return toStream ( *c , s ) ; }
  o.printClassName ( s ) ;
  return s << o.GetName() ;  
}
// ============================================================================
/* print RooAbsCollection
 *  @see  RooAbsCollection
 */
// ============================================================================
std::ostream&
Ostap::Utils::toStream
( const RooAbsCollection& o ,
  std::ostream&           s )
{
  // ==========================================================================
  static const std::string  s_open   { "{ " } ;
  static const std::string  s_close  { " }" } ;
  static const std::string  l_open   { "[ " } ;
  static const std::string  l_close  { " ]" } ;
  static const std::string  delim    { ", " } ;
  //
  const RooArgSet* aset = dynamic_cast<const RooArgSet*>( &o ) ;
  //
  const std::string& open  = aset ? s_open  : l_open  ;
  const std::string& close = aset ? s_close : l_close ;
  //
  // ==========================================================================
  return toStream ( o.begin () , o.end () , s , open , close , delim ) ;
  // ==========================================================================
}
// ============================================================================


// ============================================================================
//                                                                      The END
// ============================================================================


