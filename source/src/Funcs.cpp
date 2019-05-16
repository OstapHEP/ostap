// ============================================================================
// Include files
// ============================================================================
// local
// ============================================================================
#include "Ostap/Funcs.h"
#include "Ostap/Formula.h"
#include "Ostap/Iterator.h"
#include "Ostap/StatusCode.h"
#include "Ostap/HFuncs.h"
// ============================================================================
// Root
// ============================================================================
#include "TTree.h"
#include "RooFormulaVar.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooAbsData.h"
// ============================================================================
//  Local
// ============================================================================
#include "Exception.h"
// ============================================================================
/** @file 
 *  Implementation file for classes from namespace Ostap::Functions
 *  @date 2018-03-31 
 *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
 */
// ============================================================================
ClassImp(Ostap::Functions::FuncFormula)
ClassImp(Ostap::Functions::FuncTH1)
ClassImp(Ostap::Functions::FuncTH2)
ClassImp(Ostap::Functions::FuncTH3)
// ============================================================================
/*  constructor from the formula expression 
 *  @param expression the  formula expression 
 *  @param tree       the tree 
 *  @param name       the name for the formula 
 */
// ============================================================================
Ostap::Functions::FuncFormula::FuncFormula
( const std::string& expression , 
  const TTree*       tree       ,
  const std::string& name       )
: Ostap::IFuncTree()
  , TObject      () 
  , m_tree       ( tree       ) 
  , m_formula    ( nullptr    )    
  , m_expression ( expression )  
  , m_name       ( name       )  
{
  if ( m_tree && !make_formula () )
  { throw Ostap::Exception ( "Invalid Formula '" + m_expression + "'" , 
                             "Ostap::Function::FuncFormula"           , 
                             Ostap::StatusCode(700)                   ) ; }
}
// ============================================================================
// copy
// ============================================================================
Ostap::Functions::FuncFormula::FuncFormula
( const Ostap::Functions::FuncFormula& right )  
  : Ostap::IFuncTree( right )
  , TObject         ( right ) 
  , m_tree          ( nullptr            )  // ATTENTION! 
  , m_formula       ( nullptr            )    
  , m_expression    ( right.m_expression )  
  , m_name          ( right.m_name       )  
{}
// ============================================================================

// ============================================================================
// destructor
// ============================================================================
Ostap::Functions::FuncFormula::~FuncFormula(){}
// ============================================================================
// notify 
// ============================================================================
bool Ostap::Functions::FuncFormula::notify() const
{ return ( m_formula &&  m_formula->ok() ) ? m_formula->Notify() : false ; }
// ============================================================================
// make formula  
// ============================================================================
bool Ostap::Functions::FuncFormula::make_formula () const
{
  m_formula.reset ( nullptr ) ;
  if ( nullptr == m_tree ) { return false ; }
  TTree* t  = const_cast<TTree*> ( m_tree ) ;
  m_formula = std::make_unique<Ostap::Formula> ( m_name , m_expression , t ) ; 
  return  ( m_formula && m_formula -> ok () ) ?  m_formula->Notify() : false ;  
}
// ============================================================================
//  evaluate the formula for  TTree
// ============================================================================
double Ostap::Functions::FuncFormula::operator() ( const TTree* tree ) const
{
  //
  if ( nullptr != tree && m_tree != tree )
  {
    m_tree = tree ;
    m_formula.reset( nullptr) ; 
  }
  //
  if ( m_formula &&  m_formula->ok() &&  m_formula->GetTree() != m_tree  )
  { m_formula.reset( nullptr) ; }
  //
  Ostap::Assert ( nullptr != m_tree                ,
                  "InvalidTree"                    ,
                  "Ostap::Function::FuncFormula"   , 
                  Ostap::StatusCode(701)           ) ;
  //
  if ( !m_formula || !m_formula->ok() ) { make_formula () ;}
  //
  Ostap::Assert ( m_formula && m_formula->ok()    ,
                  "Invalid Formula"               , 
                  "Ostap::Function::FuncFormula"  , 
                  Ostap::StatusCode(700)          ) ;  
  //
  return m_formula->evaluate() ;
}
// ===========================================================================
/* constructor from the formula expression 
 *  @param expression the formula expression 
 *  @param data       the data
 *  @param name       the name for the formula 
 */
// ===========================================================================
Ostap::Functions::FuncRooFormula::FuncRooFormula 
( const std::string& expression , 
  const RooAbsData*  data       ,
  const std::string& name       ) 
  : Ostap::IFuncData () 
  , m_data       ( data       )
  , m_formula    ( nullptr    )
  , m_expression ( expression ) 
  , m_name       ( name       )
{
  if ( m_data && !make_formula() )
  { throw Ostap::Exception ( "Invalid Formula '" + m_expression + "'" , 
                             "Ostap::Function::FuncRooFormula"        , 
                             Ostap::StatusCode(706)                   ) ; }
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Functions::FuncRooFormula::~FuncRooFormula(){}
// ============================================================================
// make formula  
// ============================================================================
bool Ostap::Functions::FuncRooFormula::make_formula () const
{
  m_formula.reset ( nullptr ) ;
  if  ( !m_data ) { return false ; }
  //
  const RooArgSet* varset  = m_data->get() ;
  if (  nullptr == varset ) 
  { throw Ostap::Exception ( "Invalid RooArgSet", 
                             "Ostap::Function::FuncRooFormula"  , 
                             Ostap::StatusCode(705)             ) ; }
  RooArgList varlst ;
  Ostap::Utils::Iterator iter ( *varset ) ;
  while ( RooAbsArg* a = iter.static_next<RooAbsArg>() ) { varlst.add ( *a ) ; }
  //
  m_formula = std::make_unique<RooFormulaVar> 
    ( m_name .c_str () , m_expression.c_str () , varlst ) ;
  //
  return m_formula && m_formula -> ok () ;
}

  // ============================================================================
//  evaluate the formula for  Data
// ============================================================================
double Ostap::Functions::FuncRooFormula::operator() ( const RooAbsData* data ) const
{
  //
  if ( nullptr != data  && data != m_data )
  { 
    m_data   = data ;
    m_formula.reset ( nullptr ) ;
  }  
  //
  Ostap::Assert ( nullptr != m_data                  ,  
                  "Invalid RooAbsData"               , 
                  "Ostap::Function::FuncRooFormula"  , 
                  Ostap::StatusCode(709)             ) ; 
  //
  if ( !m_formula || !m_formula->ok() ) { make_formula () ; }
  //
  Ostap::Assert  ( m_formula && m_formula->ok()      , 
                   "Invalid RooFormula"              , 
                   "Ostap::Function::FuncRooFormula" , 
                   Ostap::StatusCode(708)            ) ; 
  //
  return m_formula->getVal() ;
}
// ============================================================================



// ===========================================================================
// FuncTH 
// ===========================================================================
/*  constructor
 *  @param tree          (INPUT) the tree 
 *  @param interpolation (INPUT) interpolation type 
 *  @param edges         (INPUT) special tretament of edges?
 *  @param extrapolate   (INPUT) use extrapolation?
 *  @param density       (INPUT) use  density?
 */
// ===========================================================================
Ostap::Functions::FuncTH::FuncTH 
( const TTree*         tree          ,
  const bool           edges         ,
  const bool           extrapolate   , 
  const bool           density       )
  : Ostap::IFuncTree() 
  , TObject        () 
  , m_tree         ( tree          ) 
  , m_edges        ( edges         ) 
  , m_extrapolate  ( extrapolate   ) 
  , m_density      ( density       )   
{}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Functions::FuncTH::FuncTH
( const Ostap::Functions::FuncTH&  right ) 
  : Ostap::IFuncTree ( right               ) 
  , TObject          ( right               ) 
  , m_tree           ( nullptr             )  //  ATTENTION! Tree is not copied!
  , m_edges          ( right.m_edges       ) 
  , m_extrapolate    ( right.m_extrapolate ) 
  , m_density        ( right.m_density     ) 
{}
// ======================================================================
// destructor 
// ======================================================================
Ostap::Functions::FuncTH::~FuncTH(){}
// ======================================================================
/*  (protected) constructor without the historgam
 *  @param xvar          (INPUT) the expression/variable 
 *  @param tree          (INPUT) the tree 
 *  @param interpolation (INPUT) interpolation type 
 *  @param edges         (INPUT) special tretament of edges?
 *  @param extrapolate   (INPUT) use extrapolation?
 *  @param density       (INPUT) use  density?
 */
// ===========================================================================
Ostap::Functions::FuncTH1::FuncTH1
( const std::string&   xvar          , 
  const TTree*         tree          ,
  const Ostap::Math::HistoInterpolation::Type tx , 
  const bool           edges         ,
  const bool           extrapolate   , 
  const bool           density       )
  : FuncTH     ( tree , edges , extrapolate , density )
  , m_xvar     ( nullptr ) 
  , m_xvar_exp ( xvar    )
  , m_tx       ( tx      )
{
  Ostap::Assert ( nullptr == m_tree || make_xvar () ,
                  "Invalid Formula '" + m_xvar_exp + "'" , 
                  "Ostap::Function::FuncTH1"             ) ;
}
// ======================================================================
/*  constructor from the histogram 
 *  @param histo         (INPUT) the historgam 
 *  @param xvar          (INPUT) the expression/variable 
 *  @param tree          (INPUT) the tree 
 *  @param interpolation (INPUT) interpolation type 
 *  @param edges         (INPUT) special tretament of edges?
 *  @param extrapolate   (INPUT) use extrapolation?
 *  @param density       (INPUT) use  density?
 */
// ===========================================================================
Ostap::Functions::FuncTH1::FuncTH1
( const TH1F&          histo         , 
  const std::string&   xvar          , 
  const TTree*         tree          ,
  const Ostap::Math::HistoInterpolation::Type tx , 
  const bool           edges         ,
  const bool           extrapolate   , 
  const bool           density       )
  : FuncTH1 ( xvar , tree , tx , edges  , extrapolate , density ) 
{
  /// copy the histogram 
  histo.Copy ( m_h1 ) ;
  m_h1.SetDirectory ( nullptr ) ;
  ///
}
// ===========================================================================
/*  constructor from the formula expression 
 *  @param histo         (INPUT) the historgam 
 *  @param xvar          (INPUT) the expression/variable 
 *  @param tree          (INPUT) the tree 
 *  @param interpolation (INPUT) interpolation type 
 *  @param edges         (INPUT) special tretament of edges?
 *  @param extrapolate   (INPUT) use extrapolation?
 *  @param density       (INPUT) use  density?
 */
// ===========================================================================
Ostap::Functions::FuncTH1::FuncTH1
( const TH1D&          histo         , 
  const std::string&   xvar          , 
  const TTree*         tree          ,
  const Ostap::Math::HistoInterpolation::Type tx , 
  const bool           edges         ,
  const bool           extrapolate   , 
  const bool           density       )
  : FuncTH1( xvar , tree , tx , edges  , extrapolate , density ) 
{
  /// copy the histogram 
  histo.Copy ( m_h1 ) ;
  m_h1.SetDirectory ( nullptr ) ;
  ///
}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Functions::FuncTH1::FuncTH1
( const Ostap::Functions::FuncTH1&  right ) 
  : FuncTH     ( right            ) 
  , m_xvar     ( nullptr          ) 
  , m_xvar_exp ( right.m_xvar_exp ) 
  , m_tx       ( right.m_tx       ) 
  , m_h1       ( right.m_h1       ) 
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Functions::FuncTH1::~FuncTH1(){}
// ============================================================================
// notify 
// ============================================================================
bool Ostap::Functions::FuncTH1::notify () const
{ return ( m_xvar &&  m_xvar->ok() ) ? m_xvar->Notify() : false ; }
// ============================================================================
// make the formula
// ============================================================================
bool Ostap::Functions::FuncTH1::make_xvar() const
{
  m_xvar.reset ( nullptr ) ;
  if ( nullptr == m_tree ) { return false ; }
  m_xvar.reset ( nullptr ) ;
  TTree* t = const_cast<TTree*> ( m_tree ) ; 
  m_xvar   = std::make_unique<Ostap::Formula> ( "" , m_xvar_exp , t ) ;
  if ( m_tree && m_xvar && m_xvar->ok() ) { m_xvar->Notify() ; }
  return m_xvar && m_xvar->ok () ;
}
// ============================================================================
//  evaluate the formula for  TTree
// ============================================================================
double Ostap::Functions::FuncTH1::operator() ( const TTree* tree ) const
{
  //
  // the tree 
  if ( nullptr != tree  && tree != m_tree )
  { 
    m_tree = tree  ;
    m_xvar.reset ( nullptr ) ;
  }
  Ostap::Assert ( nullptr != m_tree , 
                  "Invalid Tree"    , 
                  "Ostap::Function::FuncTH1" ) ;
  //
  // check consistency
  if ( m_xvar && ( m_xvar -> GetTree() != m_tree ) ) { m_xvar.reset ( nullptr ) ; }
  //
  // the  axis 
  if ( !m_xvar || !m_xvar->ok() ) { make_xvar()  ; }
  Ostap::Assert ( m_xvar && m_xvar->ok()                 , 
                  "Invalid Formula '" + m_xvar_exp + "'" , 
                  "Ostap::Function::FuncTH1"             ) ;
  //
  // agree? 
  Ostap::Assert ( m_tree == m_xvar->GetTree()            , 
                  "mismatch in tree"                     ,
                  "Ostap::Function::FuncTH1"             ) ;
  //
  const double xvar = m_xvar->evaluate() ;
  //
  return Ostap::Math::HistoInterpolation::interpolate_1D
    ( m_h1  , xvar , m_tx , m_edges , m_extrapolate , m_density ).value() ;
}
// ===========================================================================
/*  (protected) constructor without histogram 
 *  @param histo         (INPUT) the historgam 
 *  @param xvar          (INPUT) the expression/variable 
 *  @param yvar          (INPUT) the expression/variable 
 *  @param tree          (INPUT) the tree 
 *  @param tx            (INPUT) interpolation type 
 *  @param ty            (INPUT) interpolation type 
 *  @param edges         (INPUT) special tretament of edges?
 *  @param extrapolate   (INPUT) use extrapolation?
 *  @param density       (INPUT) use  density?
 */
// ===========================================================================
Ostap::Functions::FuncTH2::FuncTH2 
( const std::string&   xvar                      , 
  const std::string&   yvar                      , 
  const TTree*         tree                      ,
  const Ostap::Math::HistoInterpolation::Type tx , 
  const Ostap::Math::HistoInterpolation::Type ty , 
  const bool           edges                     ,
  const bool           extrapolate               , 
  const bool           density                   ) 
  : FuncTH ( tree , edges , extrapolate , density ) 
  , m_xvar     ( nullptr ) 
  , m_yvar     ( nullptr ) 
  , m_xvar_exp ( xvar    )
  , m_yvar_exp ( yvar    )
  , m_tx       ( tx      )
  , m_ty       ( ty      )
{
  Ostap::Assert ( nullptr == m_tree || make_xvar () , 
                  "Invalid Formula '" + m_xvar_exp + "'" ,
                  "Ostap::Function::FuncTH2"             ) ;
  Ostap::Assert ( nullptr == m_tree || make_yvar () , 
                  "Invalid Formula '" + m_yvar_exp + "'" ,
                  "Ostap::Function::FuncTH2"             ) ;
} 
// ======================================================================
/*  constructor from the histogram 
 *  @param histo         (INPUT) the historgam 
 *  @param xvar          (INPUT) the expression/variable 
 *  @param tree          (INPUT) the tree 
 *  @param interpolation (INPUT) interpolation type 
 *  @param edges         (INPUT) special tretament of edges?
 *  @param extrapolate   (INPUT) use extrapolation?
 *  @param density       (INPUT) use  density?
 */
// ======================================================================
Ostap::Functions::FuncTH2::FuncTH2 
( const TH2F&          histo                     , 
  const std::string&   xvar                      , 
  const std::string&   yvar                      , 
  const TTree*         tree                      ,
  const Ostap::Math::HistoInterpolation::Type tx ,
  const Ostap::Math::HistoInterpolation::Type ty , 
  const bool           edges                     ,
  const bool           extrapolate               , 
  const bool           density                   )
  : FuncTH2 ( xvar ,  yvar , tree ,  tx , ty, edges , extrapolate , density ) 
{
  /// copy the histogram 
  histo.Copy ( m_h2 ) ;
  m_h2.SetDirectory ( nullptr ) ;
  ///
}
// ======================================================================
/*  constructor from the histogram 
 *  @param histo         (INPUT) the historgam 
 *  @param xvar          (INPUT) the expression/variable 
 *  @param tree          (INPUT) the tree 
 *  @param interpolation (INPUT) interpolation type 
 *  @param edges         (INPUT) special tretament of edges?
 *  @param extrapolate   (INPUT) use extrapolation?
 *  @param density       (INPUT) use  density?
 */
// ======================================================================
Ostap::Functions::FuncTH2::FuncTH2 
( const TH2D&          histo                     , 
  const std::string&   xvar                      , 
  const std::string&   yvar                      , 
  const TTree*         tree                      ,
  const Ostap::Math::HistoInterpolation::Type tx ,
  const Ostap::Math::HistoInterpolation::Type ty , 
  const bool           edges                     ,
  const bool           extrapolate               , 
  const bool           density                   )
  : FuncTH2 ( xvar ,  yvar , tree ,  tx , ty, edges , extrapolate , density ) 
{
  /// copy the histogram 
  histo.Copy ( m_h2 ) ;
  m_h2.SetDirectory ( nullptr ) ;
  ///
}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Functions::FuncTH2::FuncTH2
( const Ostap::Functions::FuncTH2&  right ) 
  : FuncTH     ( right      ) 
  , m_xvar     ( nullptr    ) 
  , m_yvar     ( nullptr    ) 
  , m_xvar_exp ( right.m_xvar_exp ) 
  , m_yvar_exp ( right.m_yvar_exp ) 
  , m_tx       ( right.m_tx ) 
  , m_ty       ( right.m_ty )
  , m_h2       ( right.m_h2 )
{}
// ============================================================================
//  destructor 
// ============================================================================
Ostap::Functions::FuncTH2::~FuncTH2(){}
// ============================================================================
// notify 
// ============================================================================
bool Ostap::Functions::FuncTH2::notify () const
{  
  const bool b1 = ( m_xvar && m_xvar->ok() ) ? m_xvar->Notify() : false ; 
  const bool b2 = ( m_yvar && m_yvar->ok() ) ? m_yvar->Notify() : false ; 
  return b1 && b2 ;
}
// ============================================================================
// make the formula
// ============================================================================
bool Ostap::Functions::FuncTH2::make_xvar() const
{
  m_xvar.reset ( nullptr ) ;
  if ( nullptr == m_tree ) { return false ; }
  m_xvar.reset ( nullptr ) ;
  TTree* t = const_cast<TTree*> ( m_tree ) ; 
  m_xvar = std::make_unique<Ostap::Formula> ( "" , m_xvar_exp , t ) ;
  return m_xvar && m_xvar->ok () ;
}
// ============================================================================
// make the formula
// ============================================================================
bool Ostap::Functions::FuncTH2::make_yvar() const
{
  m_yvar.reset ( nullptr ) ;
  if ( nullptr == m_tree ) { return false ; }
  m_yvar.reset ( nullptr ) ;
  TTree* t = const_cast<TTree*> ( m_tree ) ; 
  m_yvar = std::make_unique<Ostap::Formula> ( "" , m_yvar_exp , t ) ;
  return m_yvar && m_yvar->ok () ;
}
// ============================================================================
//  evaluate the formula for  TTree
// ============================================================================
double Ostap::Functions::FuncTH2::operator() ( const TTree* tree ) const
{
  //
  if ( nullptr != tree  && tree != m_tree )
  { 
    m_tree = tree  ;
    m_xvar.reset ( nullptr ) ;
    m_yvar.reset ( nullptr ) ;
  }
  //
  Ostap::Assert ( nullptr != m_tree          , 
                  "Invalid Tree"             , 
                  "Ostap::Function::FuncTH2" ) ;
  //
  if ( !m_xvar || !m_xvar->ok() ) { make_xvar()  ; }
  if ( !m_yvar || !m_yvar->ok() ) { make_yvar()  ; }
  //
  Ostap::Assert ( m_xvar && m_xvar->ok()  ,
                  "Invalid Formula '" + m_xvar_exp + "'" , 
                  "Ostap::Function::FuncTH2"             ) ;
  Ostap::Assert ( m_yvar && m_yvar->ok()  ,
                  "Invalid Formula '" + m_yvar_exp + "'" , 
                  "Ostap::Function::FuncTH2"             ) ;
  //
  const double xvar = m_xvar->evaluate() ;
  const double yvar = m_yvar->evaluate() ;
  //
  return Ostap::Math::HistoInterpolation::interpolate_2D
    ( m_h2  , xvar , yvar , m_tx , m_ty , m_edges , m_extrapolate , m_density ).value() ;
}
// ===========================================================================
/*  (protected) constructor without histogram 
 *  @param xvar          (INPUT) the expression/variable 
 *  @param yvar          (INPUT) the expression/variable 
 *  @param tree          (INPUT) the tree 
 *  @param tx            (INPUT) interpolation type 
 *  @param ty            (INPUT) interpolation type 
 *  @param edges         (INPUT) special tretament of edges?
 *  @param extrapolate   (INPUT) use extrapolation?
 *  @param density       (INPUT) use  density?
 */
// ===========================================================================
Ostap::Functions::FuncTH3::FuncTH3 
( const std::string&   xvar                      , 
  const std::string&   yvar                      , 
  const std::string&   zvar                      , 
  const TTree*         tree                      ,
  const Ostap::Math::HistoInterpolation::Type tx , 
  const Ostap::Math::HistoInterpolation::Type ty , 
  const Ostap::Math::HistoInterpolation::Type tz , 
  const bool           edges                     ,
  const bool           extrapolate               , 
  const bool           density                   ) 
  : FuncTH ( tree , edges , extrapolate , density ) 
  , m_xvar     ( nullptr ) 
  , m_yvar     ( nullptr ) 
  , m_zvar     ( nullptr ) 
  , m_xvar_exp ( xvar    )
  , m_yvar_exp ( yvar    )
  , m_zvar_exp ( zvar    )
  , m_tx       ( tx      )
  , m_ty       ( ty      )
  , m_tz       ( tz      )
{
  Ostap::Assert ( nullptr == m_tree || make_xvar () , 
                  "Invalid Formula '" + m_xvar_exp + "'" ,
                  "Ostap::Function::FuncTH3"             ) ;
  Ostap::Assert ( nullptr == m_tree || make_yvar () , 
                  "Invalid Formula '" + m_yvar_exp + "'" ,
                  "Ostap::Function::FuncTH3"             ) ;
  Ostap::Assert ( nullptr == m_tree || make_zvar () , 
                  "Invalid Formula '" + m_zvar_exp + "'" ,
                  "Ostap::Function::FuncTH3"             ) ;
} 
// ======================================================================
/*  constructor from the histogram 
 *  @param histo         (INPUT) the historgam 
 *  @param xvar          (INPUT) the expression/variable 
 *  @param yvar          (INPUT) the expression/variable 
 *  @param zvar          (INPUT) the expression/variable 
 *  @param tree          (INPUT) the tree 
 *  @param tx            (INPUT) interpolation type 
 *  @param ty            (INPUT) interpolation type 
 *  @param tz            (INPUT) interpolation type 
 *  @param edges         (INPUT) special tretament of edges?
 *  @param extrapolate   (INPUT) use extrapolation?
 *  @param density       (INPUT) use  density?
 */
// ======================================================================
Ostap::Functions::FuncTH3::FuncTH3 
( const TH3F&          histo                     , 
  const std::string&   xvar                      , 
  const std::string&   yvar                      , 
  const std::string&   zvar                      , 
  const TTree*         tree                      ,
  const Ostap::Math::HistoInterpolation::Type tx ,
  const Ostap::Math::HistoInterpolation::Type ty , 
  const Ostap::Math::HistoInterpolation::Type tz , 
  const bool           edges                     ,
  const bool           extrapolate               , 
  const bool           density                   )
  : FuncTH3 ( xvar ,  yvar , zvar , tree ,  tx , ty , tz , edges , extrapolate , density ) 
{
  /// copy the histogram 
  histo.Copy ( m_h3 ) ;
  m_h3.SetDirectory ( nullptr ) ;
  ///
}
// ======================================================================
/*  constructor from the histogram 
 *  @param histo         (INPUT) the historgam 
 *  @param xvar          (INPUT) the expression/variable 
 *  @param tree          (INPUT) the tree 
 *  @param interpolation (INPUT) interpolation type 
 *  @param edges         (INPUT) special tretament of edges?
 *  @param extrapolate   (INPUT) use extrapolation?
 *  @param density       (INPUT) use  density?
 */
// ======================================================================
Ostap::Functions::FuncTH3::FuncTH3
( const TH3D&          histo                     , 
  const std::string&   xvar                      , 
  const std::string&   yvar                      , 
  const std::string&   zvar                      , 
  const TTree*         tree                      ,
  const Ostap::Math::HistoInterpolation::Type tx ,
  const Ostap::Math::HistoInterpolation::Type ty , 
  const Ostap::Math::HistoInterpolation::Type tz , 
  const bool           edges                     ,
  const bool           extrapolate               , 
  const bool           density                   )
  : FuncTH3 ( xvar ,  yvar , zvar , tree ,  tx , ty , tz , edges , extrapolate , density ) 
{
  /// copy the histogram 
  histo.Copy ( m_h3 ) ;
  m_h3.SetDirectory ( nullptr ) ;
  ///
}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Functions::FuncTH3::FuncTH3
( const Ostap::Functions::FuncTH3&  right ) 
  : FuncTH     ( right      ) 
  , m_xvar     ( nullptr    ) 
  , m_yvar     ( nullptr    ) 
  , m_zvar     ( nullptr    ) 
  , m_xvar_exp ( right.m_xvar_exp ) 
  , m_yvar_exp ( right.m_yvar_exp ) 
  , m_zvar_exp ( right.m_zvar_exp ) 
  , m_tx       ( right.m_tx ) 
  , m_ty       ( right.m_ty )
  , m_tz       ( right.m_tz )
  , m_h3       ( right.m_h3 )
{}
// ============================================================================
//  destructor 
// ============================================================================
Ostap::Functions::FuncTH3::~FuncTH3(){}
// ============================================================================
// notify 
// ============================================================================
bool Ostap::Functions::FuncTH3::notify () const
{  
  const bool b1 = ( m_xvar && m_xvar->ok() ) ? m_xvar->Notify() : false ; 
  const bool b2 = ( m_yvar && m_yvar->ok() ) ? m_yvar->Notify() : false ; 
  const bool b3 = ( m_zvar && m_zvar->ok() ) ? m_zvar->Notify() : false ; 
  return b1 && b2 && b3 ;
}
// ============================================================================
// make the formula
// ============================================================================
bool Ostap::Functions::FuncTH3::make_xvar() const
{
  m_xvar.reset ( nullptr ) ;
  if ( nullptr == m_tree ) { return false ; }
  m_xvar.reset ( nullptr ) ;
  TTree* t = const_cast<TTree*> ( m_tree ) ; 
  m_xvar = std::make_unique<Ostap::Formula> ( "" , m_xvar_exp , t ) ;
  return m_xvar && m_xvar->ok () ;
}
// ============================================================================
// make the formula
// ============================================================================
bool Ostap::Functions::FuncTH3::make_yvar() const
{
  m_yvar.reset ( nullptr ) ;
  if ( nullptr == m_tree ) { return false ; }
  m_yvar.reset ( nullptr ) ;
  TTree* t = const_cast<TTree*> ( m_tree ) ; 
  m_yvar = std::make_unique<Ostap::Formula> ( "" , m_yvar_exp , t ) ;
  return m_yvar && m_yvar->ok () ;
}
// ============================================================================
// make the formula
// ============================================================================
bool Ostap::Functions::FuncTH3::make_zvar() const
{
  m_zvar.reset ( nullptr ) ;
  if ( nullptr == m_tree ) { return false ; }
  m_zvar.reset ( nullptr ) ;
  TTree* t = const_cast<TTree*> ( m_tree ) ; 
  m_zvar = std::make_unique<Ostap::Formula> ( "" , m_zvar_exp , t ) ;
  return m_zvar && m_zvar->ok () ;
}
// ============================================================================
//  evaluate the formula for  TTree
// ============================================================================
double Ostap::Functions::FuncTH3::operator() ( const TTree* tree ) const
{
  //
  if ( nullptr != tree  && tree != m_tree )
  { 
    m_tree = tree  ;
    m_xvar.reset ( nullptr ) ;
    m_yvar.reset ( nullptr ) ;
    m_zvar.reset ( nullptr ) ;
  }
  //
  Ostap::Assert ( nullptr != m_tree          , 
                  "Invalid Tree"             , 
                  "Ostap::Function::FuncTH3" ) ;
  //
  if ( !m_xvar || !m_xvar->ok() ) { make_xvar()  ; }
  if ( !m_yvar || !m_yvar->ok() ) { make_yvar()  ; }
  if ( !m_zvar || !m_zvar->ok() ) { make_zvar()  ; }
  //
  Ostap::Assert ( m_xvar && m_xvar->ok()  ,
                  "Invalid Formula '" + m_xvar_exp + "'" , 
                  "Ostap::Function::FuncTH3"             ) ;
  Ostap::Assert ( m_yvar && m_yvar->ok()  ,
                  "Invalid Formula '" + m_yvar_exp + "'" , 
                  "Ostap::Function::FuncTH3"             ) ;
  Ostap::Assert ( m_zvar && m_zvar->ok()  ,
                  "Invalid Formula '" + m_zvar_exp + "'" , 
                  "Ostap::Function::FuncTH3"             ) ;
  //
  const double xvar = m_xvar->evaluate() ;
  const double yvar = m_yvar->evaluate() ;
  const double zvar = m_yvar->evaluate() ;
  //
  return Ostap::Math::HistoInterpolation::interpolate_3D
    ( m_h3 , 
      xvar , yvar , zvar , 
      m_tx , m_ty , m_tz , m_edges , m_extrapolate , m_density ).value() ;
}
// ============================================================================
// The END 
// ============================================================================
