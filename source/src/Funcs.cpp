// ============================================================================
// Include files
// ============================================================================
// ROOT&RooFit
// ============================================================================
#include "RVersion.h"
#include "TTree.h"
#include "TChain.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooAbsData.h"
// ============================================================================
// Ostap 
// ============================================================================
#include "Ostap/Funcs.h"
#include "Ostap/Formula.h"
#include "Ostap/StatusCode.h"
#include "Ostap/HFuncs.h"
#include "Ostap/FormulaVar.h"
// ============================================================================
//  Local
// ============================================================================
#include "Exception.h"
#include "local_utils.h"
#include "local_roofit.h"
#include "status_codes.h"
// ============================================================================
/** @file 
 *  Implementation file for classes from namespace Ostap::Functions
 *  @date 2018-03-31 
 *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
 */
// ============================================================================
ClassImp(Ostap::Functions::FuncFormula)
ClassImp(Ostap::Functions::FuncTree)
ClassImp(Ostap::Functions::Func1D)
ClassImp(Ostap::Functions::Func2D)
ClassImp(Ostap::Functions::Func3D)
ClassImp(Ostap::Functions::FuncTH1)
ClassImp(Ostap::Functions::FuncTH2)
ClassImp(Ostap::Functions::FuncTH3)
ClassImp(Ostap::Functions::Expression)
ClassImp(Ostap::Functions::RooTreeFun)
// ============================================================================
Ostap::Functions::FunTree::~FunTree(){};
Ostap::Functions::Func1D::~Func1D(){};
Ostap::Functions::Func2D::~Func2D(){};
Ostap::Functions::Func3D::~Func3D(){};
//
Ostap::Functions::FuncRoo1D::~FuncRoo1D(){};
Ostap::Functions::FuncRoo2D::~FuncRoo2D(){};
Ostap::Functions::FuncRoo3D::~FuncRoo3D(){};
//
Ostap::Functions::FuncTH1::~FuncTH1(){};
Ostap::Functions::FuncTH2::~FuncTH2(){};
Ostap::Functions::FuncTH3::~FuncTH3(){};
//
Ostap::Functions::FuncRooTH1::~FuncRooTH1(){};
Ostap::Functions::FuncRooTH2::~FuncRooTH2(){};
Ostap::Functions::FuncRooTH3::~FuncRooTH3(){};
//
Ostap::Functions::RooTreeFun::~RooTreeFun() {}
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
  //
  if ( nullptr != m_tree )
    {
      const TChain* chain = dynamic_cast<const TChain*>( m_tree ) ;
      if ( chain )
        {
          const TTree* chain_tree = chain->GetTree() ;
          if ( chain_tree ) { m_tree = chain_tree ; }
        }
    }
  //
  Ostap::Assert ( !m_tree || make_formula ()               , 
                  "Invalid Formula '" + m_expression + "'" , 
                  "Ostap::Function::FuncFormula"           ,
                  INVALID_FORMULA , __FILE__ , __LINE__    ) ;
  //
}
// ============================================================================
// copy
// ============================================================================
Ostap::Functions::FuncFormula::FuncFormula
( const Ostap::Functions::FuncFormula& right )  
  : Ostap::IFuncTree ( right )
  , TObject          ( right ) 
  , m_tree           ( nullptr            )  // ATTENTION! 
  , m_formula        ( nullptr            )  // ATTENTION  
  , m_expression     ( right.m_expression )  
  , m_name           ( right.m_name       )  
{}
// ============================================================================
Ostap::Functions::FuncFormula*
Ostap::Functions::FuncFormula::Clone ( const char* /* newname */ ) const 
{ return new FuncFormula ( *this ) ; }
// ============================================================================
Ostap::Functions::FuncFormula*
Ostap::Functions::FuncFormula::clone ( const char* /* newname */ ) const 
{ return new FuncFormula ( *this ) ; }
// ============================================================================
// destructor
// ============================================================================
Ostap::Functions::FuncFormula::~FuncFormula(){}
// ============================================================================
// notify 
// ============================================================================
Bool_t Ostap::Functions::FuncFormula::Notify () 
{ 
  m_formula.reset ( nullptr ) ;
  return ( m_formula &&  m_formula->ok() ) ? m_formula->Notify() : false ; 
}
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
  // if ( nullptr != tree ) 
  // {
  //   const TChain* chain = dynamic_cast<const TChain*> ( tree ) ;
  //   if ( chain ) { tree = chain->GetTree() ; }
  // }
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
  Ostap::Assert ( nullptr != m_tree                  ,
                  "Invalid Tree"                     ,
                  "Ostap::Function::FuncFormula"     , 
                  INVALID_TREE , __FILE__ , __LINE__ ) ;
  //
  if ( !m_formula || !m_formula->ok() ) { make_formula () ;}
  //
  Ostap::Assert ( m_formula && m_formula->ok()    ,
                  "Invalid Formula"               , 
                  "Ostap::Function::FuncFormula"  , 
                  INVALID_FORMULA  , __FILE__ , __LINE__) ;  
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
  , m_name       ( !name.empty() ? name : Ostap::tmp_name ( "expr_" , expression ) ) 
{
  Ostap::Assert ( !m_data || make_formula () ,
                  "Invalid Formula '" + m_expression + "'" , 
                  "Ostap::Function::FuncRooFormula"        , 
                  INVALID_FORMULA , __FILE__ , __LINE__    ) ;
}
// ============================================================================
// copy
// ============================================================================
Ostap::Functions::FuncRooFormula::FuncRooFormula
( const Ostap::Functions::FuncRooFormula& right )  
  : Ostap::IFuncData ( right )
  , m_data       ( nullptr            ) // ATTENTION! 
  , m_formula    ( nullptr            ) // ATTENTION!
  , m_expression ( right.m_expression ) 
  , m_name       ( right.m_name       )
{}
// ============================================================================
// clone 
// ============================================================================
Ostap::Functions::FuncRooFormula*
Ostap::Functions::FuncRooFormula::clone ( const char* /* name */ ) const
{ return new FuncRooFormula( *this ) ; }
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
  Ostap::Assert ( nullptr != varset                     , 
                  "Invalid RooArgSet"                   , 
                  "Ostap::Function::FuncRooFormula"     ,
                  INVALID_ARGSET  , __FILE__ , __LINE__ ) ;
  //
  RooArgList varlst ;
  //
  ::copy ( *varset , varlst ) ;
  //
  m_formula = std::make_unique<Ostap::FormulaVar> ( m_name , m_expression , varlst , false ) ;
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
  Ostap::Assert ( nullptr != m_data                     ,  
                  "Invalid RooAbsData"                  , 
                  "Ostap::Function::FuncRooFormula"     , 
                  INVALID_ABSDATA , __FILE__ , __LINE__ ) ; 
  //
  if ( !m_formula || !m_formula->ok() ) { make_formula () ; }
  //
  Ostap::Assert  ( m_formula && m_formula->ok()           , 
                   "Invalid RooFormula"                   , 
                   "Ostap::Function::FuncRooFormula"      ,
                   INVALID_FORMULA  , __FILE__ , __LINE__ ) ;
  //
  return m_formula->getVal() ;
}
// ============================================================================


// ============================================================================
Ostap::Functions::FunTree::FunTree
( std::function<double(const TTree*)> fun  ,
  const TTree*                        tree )
  : TObject ()
  , m_fun   ( fun  )
  , m_tree  ( tree )
{}
// ============================================================================
Ostap::Functions::FunTree::FunTree
( const Ostap::Functions::FunTree& right )
  : TObject ( right        )
  , m_fun   ( right.m_fun  )
  , m_tree  ( right.m_tree )
{}
// ============================================================================
// clone :
// ===========================================================================
Ostap::Functions::FunTree* 
Ostap::Functions::FunTree::clone ( const char* /* newname */ ) const
{ return new FunTree ( *this ) ; }
// ===========================================================================
Ostap::Functions::FunTree* 
Ostap::Functions::FunTree::Clone ( const char* /* newname */ ) const
{ return new FunTree ( *this ) ; }
// ===========================================================================
// notify 
// ============================================================================
Bool_t Ostap::Functions::FunTree::Notify ()  { return true ; }
// ============================================================================
//  evaluate the function for  TTree
// ============================================================================
double Ostap::Functions::FunTree::operator() ( const TTree* tree ) const
{
  //
  if ( nullptr != tree ) 
    {
      const TChain* chain = dynamic_cast<const TChain*>( tree ) ;
      if ( chain ) { tree = chain->GetTree() ; }
    }
  //
  // the tree 
  if ( tree != m_tree ) {  m_tree = tree  ; }
  //
  Ostap::Assert ( nullptr != m_tree , 
                  "Invalid Tree"    , 
                  "Ostap::Function::FunTree" ) ;
  //
  return m_fun ( m_tree ) ;
}
// ============================================================================


Ostap::Functions::Func1D::Func1D
( std::function<double(double)> fun  , 
  const std::string&            x    ,
  const TTree*                  tree ) 
  : TObject () 
  , m_fun      ( fun     )
  , m_xvar_exp ( x       ) 
  , m_xvar     { nullptr }
  , m_tree     { tree    }
{}     
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Functions::Func1D::Func1D
( const Ostap::Functions::Func1D&  right ) 
  : Ostap::IFuncTree ( right            ) 
  , TObject          ( right            ) 
  , m_fun            ( right.m_fun      )
  , m_xvar_exp       ( right.m_xvar_exp ) 
  , m_xvar           ( nullptr )
  , m_tree           ( nullptr ) 
{}
// ===========================================================================
// clone :
// ===========================================================================
Ostap::Functions::Func1D* 
Ostap::Functions::Func1D::clone ( const char* /* newname */ ) const
{ return new Func1D ( *this ) ; }
// ===========================================================================
Ostap::Functions::Func1D* 
Ostap::Functions::Func1D::Clone ( const char* /* newname */ ) const
{ return new Func1D ( *this ) ; }
// ===========================================================================
// notify 
// ============================================================================
Bool_t Ostap::Functions::Func1D::Notify () 
{
  /// attention! here  we delete the variable instead of notify/reset 
  m_xvar.reset ( nullptr ) ;
  return ( m_xvar && m_xvar->ok() ) ? m_xvar->Notify() : false ; 
}
// ============================================================================
// make the formula
// ============================================================================
bool Ostap::Functions::Func1D::make_xvar() const
{
  m_xvar.reset ( nullptr ) ;
  if ( nullptr == m_tree ) { return false ; }
  m_xvar.reset ( nullptr ) ;
  TTree* t = const_cast<TTree*> ( m_tree ) ; 
  m_xvar   = std::make_unique<Ostap::Formula> ( m_xvar_exp , t ) ;
  if ( m_tree && m_xvar && m_xvar->ok() ) { m_xvar->Notify() ; }
  return m_xvar && m_xvar->ok () ;
}
// ============================================================================
//  evaluate the function for  TTree
// ============================================================================
double Ostap::Functions::Func1D::operator() ( const TTree* tree ) const
{
  //
  if ( nullptr != tree ) 
    {
      const TChain* chain = dynamic_cast<const TChain*>( tree ) ;
      if ( chain ) { tree = chain->GetTree() ; }
    }
  //
  // the tree 
  if ( tree != m_tree )
  { 
    m_tree = tree  ;
    m_xvar.reset ( nullptr ) ;
  }
  //
  Ostap::Assert ( nullptr != m_tree , 
                  "Invalid Tree"    , 
                  "Ostap::Function::Func1D" ) ;
  //
  // check consistency
  if ( m_xvar && ( m_xvar -> GetTree() != m_tree ) ) { m_xvar.reset ( nullptr ) ; }
  //
  // the  axis 
  if ( !m_xvar || !m_xvar->ok() ) { make_xvar()  ; }
  Ostap::Assert ( m_xvar && m_xvar->ok()                 , 
                  "Invalid Formula '" + m_xvar_exp + "'" , 
                  "Ostap::Function::Func1D"              ,
                  INVALID_FORMULA  , __FILE__ , __LINE__ ) ;
  //
  // agree? 
  Ostap::Assert ( m_tree == m_xvar->GetTree()            , 
                  "Mimatch in tree"                      ,
                  "Ostap::Function::Func1D"              ,
                  MISMATCH_TREE  , __FILE__ , __LINE__   ) ;
  //
  const double xvar = m_xvar->evaluate() ;
  //
  return m_fun ( xvar ) ;
}
// ============================================================================
Ostap::Functions::Func2D::Func2D
( std::function<double(double,double)> fun  , 
  const std::string&                   x    ,
  const std::string&                   y    ,
  const TTree*                         tree ) 
  : TObject () 
  , m_fun      ( fun     )
  , m_xvar_exp ( x       ) 
  , m_yvar_exp ( y       ) 
  , m_xvar     { nullptr }
  , m_yvar     { nullptr }
  , m_tree     { tree    }
{}     
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Functions::Func2D::Func2D
( const Ostap::Functions::Func2D&  right ) 
  : Ostap::IFuncTree ( right            ) 
  , TObject          ( right            ) 
  , m_fun            ( right.m_fun      )
  , m_xvar_exp       ( right.m_xvar_exp ) 
  , m_yvar_exp       ( right.m_yvar_exp ) 
  , m_xvar           ( nullptr )
  , m_yvar           ( nullptr )
  , m_tree           ( nullptr ) 
{}
// ===========================================================================
// clone :
// ===========================================================================
Ostap::Functions::Func2D* 
Ostap::Functions::Func2D::clone ( const char* /* newname */ ) const
{ return new Func2D ( *this ) ; }
// ===========================================================================
Ostap::Functions::Func2D* 
Ostap::Functions::Func2D::Clone ( const char* /* newname */ ) const
{ return new Func2D ( *this ) ; }
// ===========================================================================
// notify 
// ============================================================================
Bool_t Ostap::Functions::Func2D::Notify () 
{  
  m_xvar.reset ( nullptr ) ;
  m_yvar.reset ( nullptr ) ;
  //
  const bool b1 = ( m_xvar && m_xvar->ok() ) ? m_xvar->Notify() : false ; 
  const bool b2 = ( m_yvar && m_yvar->ok() ) ? m_yvar->Notify() : false ; 
  //
  return b1 && b2 ;
}
// ============================================================================
// make the formula
// ============================================================================
bool Ostap::Functions::Func2D::make_xvar() const
{
  m_xvar.reset ( nullptr ) ;
  if ( nullptr == m_tree ) { return false ; }
  m_xvar.reset ( nullptr ) ;
  TTree* t = const_cast<TTree*> ( m_tree ) ; 
  m_xvar   = std::make_unique<Ostap::Formula> ( m_xvar_exp , t ) ;
  if ( m_tree && m_xvar && m_xvar->ok() ) { m_xvar->Notify() ; }
  return m_xvar && m_xvar->ok () ;
}
// ============================================================================
// make the formula
// ============================================================================
bool Ostap::Functions::Func2D::make_yvar() const
{
  m_yvar.reset ( nullptr ) ;
  if ( nullptr == m_tree ) { return false ; }
  m_yvar.reset ( nullptr ) ;
  TTree* t = const_cast<TTree*> ( m_tree ) ; 
  m_yvar   = std::make_unique<Ostap::Formula> ( m_yvar_exp , t ) ;
  if ( m_tree && m_yvar && m_yvar->ok() ) { m_yvar->Notify() ; }
  return m_yvar && m_yvar->ok () ;
}
// ============================================================================
//  evaluate the function for  TTree
// ============================================================================
double Ostap::Functions::Func2D::operator() ( const TTree* tree ) const
{
  //
  if ( nullptr != tree ) 
  {
    const TChain* chain = dynamic_cast<const TChain*> ( tree ) ;
    if ( chain ) { tree = chain->GetTree() ; }
  }
  //
  // the tree 
  if ( tree != m_tree )
  { 
    m_tree = tree  ;
    m_xvar.reset ( nullptr ) ;
    m_yvar.reset ( nullptr ) ;
  }
  //
  Ostap::Assert ( nullptr != m_tree                  , 
                  "Invalid Tree"                     ,   
                  "Ostap::Function::Func2D"          ,
                  INVALID_TREE , __FILE__ , __LINE__ ) ;
  //
  // check consistency
  if ( m_xvar && ( m_xvar -> GetTree() != m_tree ) ) { m_xvar.reset ( nullptr ) ; }
  if ( m_yvar && ( m_yvar -> GetTree() != m_tree ) ) { m_yvar.reset ( nullptr ) ; }
  //
  // the  axis 
  if ( !m_xvar || !m_xvar->ok() ) { make_xvar()  ; }
  Ostap::Assert ( m_xvar && m_xvar->ok()                 , 
                  "Invalid Formula '" + m_xvar_exp + "'" , 
                  "Ostap::Function::Func2D"              ,
                  INVALID_FORMULA , __FILE__ , __LINE__  ) ;
  // the  axis 
  if ( !m_yvar || !m_yvar->ok() ) { make_yvar()  ; }
  Ostap::Assert ( m_yvar && m_yvar->ok()                 , 
                  "Invalid Formula '" + m_yvar_exp + "'" , 
                  "Ostap::Function::Func2D"              ,
                  INVALID_FORMULA , __FILE__ , __LINE__  ) ;
  //
  // agree? 
  Ostap::Assert ( m_tree == m_xvar->GetTree() && 
                  m_tree == m_yvar->GetTree()           ,
                  "mismatch in tree"                    ,
                  "Ostap::Function::Func2D"             ,
                  MISMATCH_TREE  , __FILE__ , __LINE__  ) ;
  //
  const double xvar = m_xvar->evaluate() ;
  const double yvar = m_yvar->evaluate() ;
  //
  return m_fun ( xvar , yvar ) ;
}
// ============================================================================

// ============================================================================
// Func3D 
// ============================================================================
Ostap::Functions::Func3D::Func3D
( std::function<double(double,double,double)> fun  , 
  const std::string&                          x    ,
  const std::string&                          y    ,
  const std::string&                          z    ,
  const TTree*                                tree ) 
  : TObject () 
  , m_fun      ( fun     )
  , m_xvar_exp ( x       ) 
  , m_yvar_exp ( y       ) 
  , m_zvar_exp ( z       ) 
  , m_xvar     { nullptr }
  , m_yvar     { nullptr }
  , m_zvar     { nullptr }
  , m_tree     { tree    }
{}   
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Functions::Func3D::Func3D
( const Ostap::Functions::Func3D&  right ) 
  : Ostap::IFuncTree ( right            ) 
  , TObject          ( right            ) 
  , m_fun            ( right.m_fun      )
  , m_xvar_exp       ( right.m_xvar_exp ) 
  , m_yvar_exp       ( right.m_yvar_exp ) 
  , m_zvar_exp       ( right.m_zvar_exp ) 
  , m_xvar           ( nullptr )
  , m_yvar           ( nullptr )
  , m_zvar           ( nullptr )
  , m_tree           ( nullptr ) 
{}
// ===========================================================================
// clone :
// ===========================================================================
Ostap::Functions::Func3D* 
Ostap::Functions::Func3D::clone ( const char* /* newname */ ) const
{ return new Func3D ( *this ) ; }
// ===========================================================================
Ostap::Functions::Func3D* 
Ostap::Functions::Func3D::Clone ( const char* /* newname */ ) const
{ return new Func3D ( *this ) ; }
// ===========================================================================
// notify 
// ============================================================================
Bool_t Ostap::Functions::Func3D::Notify () 
{  
  m_xvar.reset ( nullptr ) ;
  m_yvar.reset ( nullptr ) ;
  m_zvar.reset ( nullptr ) ;
  //
  const bool b1 = ( m_xvar && m_xvar->ok() ) ? m_xvar->Notify() : false ; 
  const bool b2 = ( m_yvar && m_yvar->ok() ) ? m_yvar->Notify() : false ; 
  const bool b3 = ( m_zvar && m_zvar->ok() ) ? m_zvar->Notify() : false ; 
  //
  return b1 && b2 && b3 ;
}
// ============================================================================
// make the formula
// ============================================================================
bool Ostap::Functions::Func3D::make_xvar() const
{
  m_xvar.reset ( nullptr ) ;
  if ( nullptr == m_tree ) { return false ; }
  m_xvar.reset ( nullptr ) ;
  TTree* t = const_cast<TTree*> ( m_tree ) ; 
  m_xvar   = std::make_unique<Ostap::Formula> ( m_xvar_exp , t ) ;
  if ( m_tree && m_xvar && m_xvar->ok() ) { m_xvar->Notify() ; }
  return m_xvar && m_xvar->ok () ;
}
// ============================================================================
// make the formula
// ============================================================================
bool Ostap::Functions::Func3D::make_yvar() const
{
  m_yvar.reset ( nullptr ) ;
  if ( nullptr == m_tree ) { return false ; }
  m_yvar.reset ( nullptr ) ;
  TTree* t = const_cast<TTree*> ( m_tree ) ; 
  m_yvar   = std::make_unique<Ostap::Formula> ( m_yvar_exp , t ) ;
  if ( m_tree && m_yvar && m_yvar->ok() ) { m_yvar->Notify() ; }
  return m_yvar && m_yvar->ok () ;
}
// ============================================================================
// make the formula
// ============================================================================
bool Ostap::Functions::Func3D::make_zvar() const
{
  m_zvar.reset ( nullptr ) ;
  if ( nullptr == m_tree ) { return false ; }
  m_zvar.reset ( nullptr ) ;
  TTree* t = const_cast<TTree*> ( m_tree ) ; 
  m_zvar   = std::make_unique<Ostap::Formula> ( m_zvar_exp , t ) ;
  if ( m_tree && m_zvar && m_zvar->ok() ) { m_zvar->Notify() ; }
  return m_zvar && m_zvar->ok () ;
}
// ============================================================================
//  evaluate the function for  TTree
// ============================================================================
double Ostap::Functions::Func3D::operator() ( const TTree* tree ) const
{
  //
  if ( nullptr != tree ) 
  {
    const TChain* chain = dynamic_cast<const TChain*> ( tree ) ;
    if ( chain ) { tree = chain->GetTree() ; }
  }
  // the tree 
  if ( tree != m_tree )
  { 
    m_tree = tree  ;
    m_xvar.reset ( nullptr ) ;
    m_yvar.reset ( nullptr ) ;
    m_zvar.reset ( nullptr ) ;
  }
  //
  Ostap::Assert ( nullptr != m_tree                  , 
                  "Invalid Tree"                     , 
                  "Ostap::Function::Func3D"          ,
                  INVALID_TREE , __FILE__ , __LINE__ ) ;
  //
  // check consistency
  if ( m_xvar && ( m_xvar -> GetTree() != m_tree ) ) { m_xvar.reset ( nullptr ) ; }
  if ( m_yvar && ( m_yvar -> GetTree() != m_tree ) ) { m_yvar.reset ( nullptr ) ; }
  if ( m_zvar && ( m_zvar -> GetTree() != m_tree ) ) { m_zvar.reset ( nullptr ) ; }
  //
  // the  axis 
  if ( !m_xvar || !m_xvar->ok() ) { make_xvar()  ; }
  Ostap::Assert ( m_xvar && m_xvar->ok()                 , 
                  "Invalid Formula '" + m_xvar_exp + "'" , 
                  "Ostap::Function::Func2D"              , 
                  INVALID_FORMULA , __FILE__ , __LINE__  ) ;
  // the  axis 
  if ( !m_yvar || !m_yvar->ok() ) { make_yvar()  ; }
  Ostap::Assert ( m_yvar && m_yvar->ok()                 , 
                  "Invalid Formula '" + m_yvar_exp + "'" , 
                  "Ostap::Function::Func2D"              ,
                  INVALID_FORMULA , __FILE__ , __LINE__  ) ;
  // the  axis 
  if ( !m_zvar || !m_zvar->ok() ) { make_zvar()  ; }
  Ostap::Assert ( m_zvar && m_zvar->ok()                 , 
                  "Invalid Formula '" + m_zvar_exp + "'" , 
                  "Ostap::Function::Func3D"              , 
                  INVALID_FORMULA , __FILE__ , __LINE__  ) ;
  //
  // agree? 
  Ostap::Assert ( m_tree == m_xvar->GetTree() && 
                  m_tree == m_yvar->GetTree() && 
                  m_tree == m_zvar->GetTree()          , 
                  "mismatch in tree"                   ,
                  "Ostap::Function::Func2D"            , 
                  MISMATCH_TREE  , __FILE__ , __LINE__ ) ;
  //
  const double xvar = m_xvar->evaluate() ;
  const double yvar = m_yvar->evaluate() ;
  const double zvar = m_zvar->evaluate() ;
  //
  return m_fun ( xvar , yvar , zvar ) ;
}
// ======================================================================
/*  constructor from the histogram 
 *  @param histo         (INPUT) the histogram 
 *  @param xvar          (INPUT) the expression/variable 
 *  @param tree          (INPUT) the tree 
 *  @param interpolation (INPUT) interpolation type 
 *  @param edges         (INPUT) special tretament of edges?
 *  @param extrapolate   (INPUT) use extrapolation?
 *  @param density       (INPUT) use  density?
 */
// ===========================================================================
Ostap::Functions::FuncTH1::FuncTH1
( const TH1&           histo         , 
  const std::string&   xvar          , 
  const TTree*         tree          ,
  const Ostap::Math::HistoInterpolation::Type tx , 
  const bool           edges         ,
  const bool           extrapolate   , 
  const bool           density       )
  : FuncTH1 ( Ostap::Math::Histo1D ( histo , tx , edges , extrapolate , density ) , 
	      xvar , 
	      tree ) 
{}
// ============================================================================
/*  constructor from the histogram 
 *  @param histo         (INPUT) the historgam 
 *  @param xvar          (INPUT) the expression/variable 
 *  @param tree          (INPUT) the tree 
 */
// ============================================================================
Ostap::Functions::FuncTH1::FuncTH1
( const Ostap::Math::Histo1D& histo , 
  const std::string&          xvar  , 
  const TTree*                tree  ) 
  : Func1D  ( histo , xvar , tree ) 
  , m_histo ( histo ) 
{}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Functions::FuncTH1::FuncTH1
( const Ostap::Functions::FuncTH1&  right ) 
  : Func1D  ( right         )
  , m_histo ( right.m_histo ) 
{}
// ===========================================================================
// clone :
// ===========================================================================
Ostap::Functions::FuncTH1* 
Ostap::Functions::FuncTH1::clone ( const char* /* newname */ ) const
{ return new FuncTH1 ( *this ) ; }
// ============================================================================
Ostap::Functions::FuncTH1* 
Ostap::Functions::FuncTH1::Clone ( const char* /* newname */ ) const
{ return new FuncTH1 ( *this ) ; }
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
( const TH2&           histo                     , 
  const std::string&   xvar                      , 
  const std::string&   yvar                      , 
  const TTree*         tree                      ,
  const Ostap::Math::HistoInterpolation::Type tx ,
  const Ostap::Math::HistoInterpolation::Type ty , 
  const bool           edges                     ,
  const bool           extrapolate               , 
  const bool           density                   )
  : FuncTH2 ( Ostap::Math::Histo2D ( histo , tx , ty , edges , extrapolate , density ) , 
	      xvar ,   
	      yvar , 
	      tree )
{}
// ======================================================================
/*  constructor from the histogram 
 *  @param histo         (INPUT) the historgam 
 *  @param xvar          (INPUT) the expression/variable 
 *  @param tree          (INPUT) the tree 
 */
// ======================================================================
Ostap::Functions::FuncTH2::FuncTH2 
( const Ostap::Math::Histo2D& histo , 
  const std::string&          xvar  , 
  const std::string&          yvar  , 
  const TTree*                tree  ) 
  : Func2D ( histo , xvar , yvar , tree )
  , m_histo ( histo ) 
{}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Functions::FuncTH2::FuncTH2
( const Ostap::Functions::FuncTH2&  right ) 
  : Func2D  ( right ) 
  , m_histo ( right.m_histo ) 
{}
// ===========================================================================
// clone :
// ===========================================================================
Ostap::Functions::FuncTH2* 
Ostap::Functions::FuncTH2::clone ( const char* /* newname */ ) const
{ return new FuncTH2 ( *this ) ; }
// ===========================================================================
Ostap::Functions::FuncTH2* 
Ostap::Functions::FuncTH2::Clone ( const char* /* newname */ ) const
{ return new FuncTH2 ( *this ) ; }
// ===========================================================================


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
( const TH3&           histo                     , 
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
: FuncTH3 ( Ostap::Math::Histo3D ( histo , tx , ty , tz , edges , extrapolate , density ) , 
	    xvar , 
	    yvar , 
	    zvar , 
	    tree )
{}
// ======================================================================
/*  constructor from the histogram 
 *  @param histo         (INPUT) the historgam 
 *  @param xvar          (INPUT) the expression/variable 
 *  @param tree          (INPUT) the tree 
 */
// ======================================================================
Ostap::Functions::FuncTH3::FuncTH3
( const Ostap::Math::Histo3D& histo , 
  const std::string&          xvar  , 
  const std::string&          yvar  , 
  const std::string&          zvar  , 
  const TTree*                tree  )
  : Func3D ( histo , xvar , yvar , zvar , tree )
  , m_histo ( histo ) 
{}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Functions::FuncTH3::FuncTH3
( const Ostap::Functions::FuncTH3&  right ) 
  : Func3D  ( right  )
  , m_histo ( right.m_histo ) 
{}
// ===========================================================================
// clone :
// ============================================================================
Ostap::Functions::FuncTH3* 
Ostap::Functions::FuncTH3::clone ( const char* /* newname */ ) const
{ return new FuncTH3 ( *this ) ; }
// ============================================================================
Ostap::Functions::FuncTH3* 
Ostap::Functions::FuncTH3::Clone ( const char* /* newname */ ) const
{ return new FuncTH3 ( *this ) ; }
// ============================================================================





// ============================================================================
// RooAbsData functions
// ============================================================================

// ============================================================================
// copy constructor
// ============================================================================
Ostap::Functions::FuncRoo1D::FuncRoo1D 
( const Ostap::Functions::FuncRoo1D& right ) 
  : Ostap::IFuncData (  right ) 
  , m_fun      ( right.m_fun      ) 
  , m_xvar_exp ( right.m_xvar_exp )
  , m_xvar     ( nullptr )   
  , m_data     ( nullptr ) 
{}
// ============================================================================
// IFuncData::clone
// ============================================================================
Ostap::Functions::FuncRoo1D*
Ostap::Functions::FuncRoo1D::clone ( const char* /* name */  ) const
{ return  new FuncRoo1D ( *this ) ; }
// ============================================================================
// make formula 
// ============================================================================
bool Ostap::Functions::FuncRoo1D::make_xvar () const 
{
  m_xvar.reset ( nullptr ) ;
  if ( nullptr == m_data ) { return false ; }
  m_xvar.reset ( nullptr ) ;
  // 
  const RooArgSet* varset  = m_data->get() ;
  Ostap::Assert ( nullptr != varset                     ,
                  "Invalid RooArgSet"                   , 
                  "Ostap::Function::FuncRoo1D"          ,
                  INVALID_ARGSET  , __FILE__ , __LINE__ ) ;
  //
  RooArgList varlst ;
  ::copy ( *varset , varlst ) ;
  //
  m_xvar = std::make_unique<Ostap::FormulaVar> ( m_xvar_exp , varlst , false ) ;
  //
  return m_xvar && m_xvar -> ok () ;
}
// ============================================================================
//  evaluate the formula for RooAbsData
// ============================================================================
double Ostap::Functions::FuncRoo1D::operator() 
  ( const RooAbsData* data ) const
{
  //
  if ( nullptr != data  && data != m_data )
  { 
    m_data   = data ;
    m_xvar.reset ( nullptr ) ;
  }  
  //
  Ostap::Assert ( nullptr != m_data                     ,  
                  "Invalid RooAbsData"                  , 
                  "Ostap::Function::FuncRoo1D"          ,
                  INVALID_ABSDATA , __FILE__ , __LINE__ ) ;
  //
  if ( !m_xvar || !m_xvar->ok() ) { make_xvar () ; }
  //
  Ostap::Assert  ( m_xvar && m_xvar->ok()        , 
                   "Invalid RooFormula"          , 
                   "Ostap::Function::FuncRoo1D"  ,
                   INVALID_FORMULA , __FILE__ , __LINE__ ) ;
  //
  const double x = m_xvar->getVal() ;
  //
  return m_fun ( x ) ;
}
// ======================================================================

// ======================================================================
// copy constructor
// ======================================================================
Ostap::Functions::FuncRoo2D::FuncRoo2D 
( const Ostap::Functions::FuncRoo2D& right ) 
  : Ostap::IFuncData (  right ) 
  , m_fun      ( right.m_fun      ) 
  , m_xvar_exp ( right.m_xvar_exp )
  , m_yvar_exp ( right.m_yvar_exp )
  , m_xvar     ( nullptr )   
  , m_yvar     ( nullptr )   
  , m_data     ( nullptr ) 
{}
// ============================================================================
// IFuncData::clone
// ============================================================================
Ostap::Functions::FuncRoo2D*
Ostap::Functions::FuncRoo2D::clone ( const char* /* name */  ) const
{ return  new FuncRoo2D ( *this ) ; }
// ============================================================================
// make formula 
// ============================================================================
bool Ostap::Functions::FuncRoo2D::make_xvar () const 
{
  m_xvar.reset ( nullptr ) ;
  if ( nullptr == m_data ) { return false ; }
  m_xvar.reset ( nullptr ) ;
  // 
  const RooArgSet* varset  = m_data->get() ;
  Ostap::Assert ( nullptr != varset                     ,
                  "Invalid RooArgSet"                   , 
                  "Ostap::Function::FuncRoo2D"          ,
                  INVALID_ARGSET  , __FILE__ , __LINE__ ) ;
  //
  RooArgList varlst ;
  ::copy ( *varset , varlst ) ;
  //
  m_xvar = std::make_unique<Ostap::FormulaVar> ( m_xvar_exp , varlst , false ) ;
  //
  return m_xvar && m_xvar -> ok () ;
}
// ============================================================================
// make formula 
// ============================================================================
bool Ostap::Functions::FuncRoo2D::make_yvar () const 
{
  m_yvar.reset ( nullptr ) ;
  if ( nullptr == m_data ) { return false ; }
  m_yvar.reset ( nullptr ) ;
  // 
  const RooArgSet* varset  = m_data->get() ;
  Ostap::Assert ( nullptr != varset                     ,
                  "Invalid RooArgSet"                   , 
                  "Ostap::Function::FuncRoo2D"          ,
                  INVALID_ARGSET  , __FILE__ , __LINE__ ) ;
  //
  RooArgList varlst ;
  ::copy ( *varset , varlst ) ;
  //
  m_yvar = std::make_unique<Ostap::FormulaVar> ( m_yvar_exp , varlst , false ) ;
  //
  return m_yvar && m_yvar -> ok () ;
}
// ============================================================================
//  evaluate the formula for RooAbsData
// ============================================================================
double Ostap::Functions::FuncRoo2D::operator() 
  ( const RooAbsData* data ) const
{
  //
  if ( nullptr != data  && data != m_data )
  { 
    m_data   = data ;
    m_xvar.reset ( nullptr ) ;
    m_yvar.reset ( nullptr ) ;
  }  
  //
  Ostap::Assert ( nullptr != m_data                     ,  
                  "Invalid RooAbsData"                  , 
                  "Ostap::Function::FuncRoo2D"          ,
                  INVALID_ABSDATA , __FILE__ , __LINE__ ) ;
  //
  if ( !m_xvar || !m_xvar->ok() ) { make_xvar () ; }
  Ostap::Assert  ( m_xvar && m_xvar->ok()                , 
                   "Invalid RooFormula"                  , 
                   "Ostap::Function::FuncRoo2D"          ,
                   INVALID_FORMULA , __FILE__ , __LINE__ ) ; 
  //
  if ( !m_yvar || !m_yvar->ok() ) { make_yvar () ; }
  Ostap::Assert  ( m_yvar && m_yvar->ok()                , 
                   "Invalid RooFormula"                  , 
                   "Ostap::Function::FuncRoo2D"          , 
                   INVALID_FORMULA , __FILE__ , __LINE__ ) ; 
  //
  const double x = m_xvar->getVal() ;
  const double y = m_yvar->getVal() ;
  //
  return m_fun ( x , y ) ;
}


// ======================================================================
// copy constructor
// ======================================================================
Ostap::Functions::FuncRoo3D::FuncRoo3D 
( const Ostap::Functions::FuncRoo3D& right ) 
  : Ostap::IFuncData (  right ) 
  , m_fun      ( right.m_fun      ) 
  , m_xvar_exp ( right.m_xvar_exp )
  , m_yvar_exp ( right.m_yvar_exp )
  , m_zvar_exp ( right.m_zvar_exp )
  , m_xvar     ( nullptr )   
  , m_yvar     ( nullptr )   
  , m_zvar     ( nullptr )   
  , m_data     ( nullptr ) 
{}
// ============================================================================
// IFuncData::clone
// ============================================================================
Ostap::Functions::FuncRoo3D*
Ostap::Functions::FuncRoo3D::clone ( const char* /* name */  ) const
{ return  new FuncRoo3D ( *this ) ; }
// ============================================================================
// make formula 
// ============================================================================
bool Ostap::Functions::FuncRoo3D::make_xvar () const 
{
  m_xvar.reset ( nullptr ) ;
  if ( nullptr == m_data ) { return false ; }
  m_xvar.reset ( nullptr ) ;
  // 
  const RooArgSet* varset  = m_data->get() ;  
  Ostap::Assert ( nullptr != varset                    ,
                  "Invalid RooArgSet"                  , 
                  "Ostap::Function::FuncRoo3D"         ,
                  INVALID_ARGSET , __FILE__ , __LINE__ ) ;
  //
  RooArgList varlst ;
  ::copy ( *varset , varlst ) ;
  //
  m_xvar = std::make_unique<Ostap::FormulaVar> ( m_xvar_exp , varlst , false ) ;
  //
  return m_xvar && m_xvar -> ok () ;
}
// ============================================================================
// make formula 
// ============================================================================
bool Ostap::Functions::FuncRoo3D::make_yvar () const 
{
  m_yvar.reset ( nullptr ) ;
  if ( nullptr == m_data ) { return false ; }
  m_yvar.reset ( nullptr ) ;
  // 
  const RooArgSet* varset  = m_data->get() ;
  Ostap::Assert ( nullptr != varset                    ,
                  "Invalid RooArgSet"                  , 
                  "Ostap::Function::FuncRoo3D"         ,
                  INVALID_ARGSET , __FILE__ , __LINE__ ) ;
  //
  RooArgList varlst ;
  ::copy ( *varset , varlst ) ;
  //
  m_yvar = std::make_unique<Ostap::FormulaVar> ( m_yvar_exp , varlst , false ) ;
  //
  return m_yvar && m_yvar -> ok () ;
}
// ============================================================================
// make formula 
// ============================================================================
bool Ostap::Functions::FuncRoo3D::make_zvar () const 
{
  m_zvar.reset ( nullptr ) ;
  if ( nullptr == m_data ) { return false ; }
  m_zvar.reset ( nullptr ) ;
  // 
  const RooArgSet* varset  = m_data->get() ;
  Ostap::Assert ( nullptr != varset                    ,
                  "Invalid RooArgSet"                  , 
                  "Ostap::Function::FuncRoo3D"         ,
                  INVALID_ARGSET , __FILE__ , __LINE__ ) ;
  //
  RooArgList varlst ;
  ::copy ( *varset , varlst ) ;
  //
  m_zvar = std::make_unique<Ostap::FormulaVar> ( m_zvar_exp , varlst , false ) ;
  //
  return m_zvar && m_zvar -> ok () ;
}
// ============================================================================
//  evaluate the formula for RooAbsData
// ============================================================================
double Ostap::Functions::FuncRoo3D::operator() 
  ( const RooAbsData* data ) const
{
  //
  if ( nullptr != data  && data != m_data )
  { 
    m_data   = data ;
    m_xvar.reset ( nullptr ) ;
    m_yvar.reset ( nullptr ) ;
    m_zvar.reset ( nullptr ) ;
  }  
  //
  Ostap::Assert ( nullptr != m_data                     ,  
                  "Invalid RooAbsData"                  , 
                  "Ostap::Function::FuncRoo2D"          ,
                  INVALID_ABSDATA , __FILE__ , __LINE__ ) ;
  //
  if ( !m_xvar || !m_xvar->ok() ) { make_xvar () ; }
  Ostap::Assert  ( m_xvar && m_xvar->ok()                , 
                   "Invalid RooFormula"                  , 
                   "Ostap::Function::FuncRoo2D"          ,
                   INVALID_FORMULA , __FILE__ , __LINE__ ) ;
  //
  if ( !m_yvar || !m_yvar->ok() ) { make_yvar () ; }
  Ostap::Assert  ( m_yvar && m_yvar->ok()                , 
                   "Invalid RooFormula"                  , 
                   "Ostap::Function::FuncRoo2D"          , 
                   INVALID_FORMULA , __FILE__ , __LINE__ ) ;
  //
  if ( !m_zvar || !m_zvar->ok() ) { make_zvar () ; }
  Ostap::Assert  ( m_zvar && m_zvar->ok()                , 
                   "Invalid RooFormula"                  , 
                   "Ostap::Function::FuncRoo2D"          , 
                   INVALID_FORMULA , __FILE__ , __LINE__ ) ;
  //
  const double x = m_xvar->getVal() ;
  const double y = m_yvar->getVal() ;
  const double z = m_zvar->getVal() ;
  //
  return m_fun ( x , y , z ) ;
}



 



// ======================================================================
/*  constructor from the histogram 
 *  @param histo         (INPUT) the histogram 
 *  @param xvar          (INPUT) the expression/variable 
 *  @param tree          (INPUT) the tree 
 *  @param interpolation (INPUT) interpolation type 
 *  @param edges         (INPUT) special tretament of edges?
 *  @param extrapolate   (INPUT) use extrapolation?
 *  @param density       (INPUT) use  density?
 */
// ===========================================================================
Ostap::Functions::FuncRooTH1::FuncRooTH1
( const TH1&           histo         , 
  const std::string&   xvar          , 
  const RooAbsData*    data          ,
  const Ostap::Math::HistoInterpolation::Type tx , 
  const bool           edges         ,
  const bool           extrapolate   , 
  const bool           density       )
  : FuncRoo1D ( Ostap::Math::Histo1D ( histo , tx , edges , extrapolate , density ) , 
                xvar , 
                data )
{}
// ============================================================================
/*  constructor from the histogram 
 *  @param histo         (INPUT) the historgam 
 *  @param xvar          (INPUT) the expression/variable 
 *  @param tree          (INPUT) the tree 
 */
// ============================================================================
Ostap::Functions::FuncRooTH1::FuncRooTH1
( const Ostap::Math::Histo1D& histo , 
  const std::string&          xvar  , 
  const RooAbsData*           data  )
  : FuncRoo1D ( histo , xvar , data )
  , m_histo   ( histo )
{}
// ============================================================================
/*  constructor from the histogram 
 *  @param histo         (INPUT) the histogram 
 *  @param xvar          (INPUT) the expression/variable 
 *  @param yvar          (INPUT) the expression/variable 
 *  @param tree          (INPUT) the tree 
 *  @param interpolation (INPUT) interpolation type 
 *  @param edges         (INPUT) special tretament of edges?
 *  @param extrapolate   (INPUT) use extrapolation?
 *  @param density       (INPUT) use  density?
 */
// ===========================================================================
Ostap::Functions::FuncRooTH2::FuncRooTH2
( const TH2&           histo         , 
  const std::string&   xvar          , 
  const std::string&   yvar          , 
  const RooAbsData*    data          ,
  const Ostap::Math::HistoInterpolation::Type tx , 
  const Ostap::Math::HistoInterpolation::Type ty , 
  const bool           edges         ,
  const bool           extrapolate   , 
  const bool           density       )
  : FuncRoo2D ( Ostap::Math::Histo2D ( histo , tx , ty , edges , extrapolate , density ) , 
                xvar , 
                yvar , 
                data ) 
{}
// ============================================================================
/*  constructor from the histogram 
 *  @param histo         (INPUT) the historgam 
 *  @param xvar          (INPUT) the expression/variable 
 *  @param tree          (INPUT) the tree 
 */
// ============================================================================
Ostap::Functions::FuncRooTH2::FuncRooTH2
( const Ostap::Math::Histo2D& histo , 
  const std::string&          xvar  , 
  const std::string&          yvar  , 
  const RooAbsData*           data  )
  : FuncRoo2D ( histo , xvar , yvar , data )
  , m_histo   ( histo )
{}
// ======================================================================
/*  constructor from the histogram 
 *  @param histo         (INPUT) the histogram 
 *  @param xvar          (INPUT) the expression/variable 
 *  @param yvar          (INPUT) the expression/variable 
 *  @param zvar          (INPUT) the expression/variable 
 *  @param tree          (INPUT) the tree 
 *  @param interpolation (INPUT) interpolation type 
 *  @param edges         (INPUT) special tretament of edges?
 *  @param extrapolate   (INPUT) use extrapolation?
 *  @param density       (INPUT) use  density?
 */
// ===========================================================================
Ostap::Functions::FuncRooTH3::FuncRooTH3
( const TH3&           histo         , 
  const std::string&   xvar          , 
  const std::string&   yvar          , 
  const std::string&   zvar          , 
  const RooAbsData*    data          ,
  const Ostap::Math::HistoInterpolation::Type tx , 
  const Ostap::Math::HistoInterpolation::Type ty , 
  const Ostap::Math::HistoInterpolation::Type tz , 
  const bool           edges         ,
  const bool           extrapolate   , 
  const bool           density       )
  : FuncRoo3D (  Ostap::Math::Histo3D ( histo , tx , ty , tz , edges , extrapolate , density ) , 
                 xvar , 
                 yvar , 
                 zvar , 
                 data ) 
{}
// ============================================================================
/*  constructor from the histogram 
 *  @param histo         (INPUT) the historgam 
 *  @param xvar          (INPUT) the expression/variable 
 *  @param tree          (INPUT) the tree 
 */
// ============================================================================
Ostap::Functions::FuncRooTH3::FuncRooTH3
( const Ostap::Math::Histo3D& histo , 
  const std::string&          xvar  , 
  const std::string&          yvar  , 
  const std::string&          zvar  , 
  const RooAbsData*           data  )
  : FuncRoo3D ( histo , xvar , yvar , zvar , data )
  , m_histo   ( histo )
{}
// ============================================================================

// ============================================================================
// IFuncData::clone
// ============================================================================
Ostap::Functions::FuncRooTH1*
Ostap::Functions::FuncRooTH1::clone ( const char* /* name */  ) const
{ return  new FuncRooTH1 ( *this ) ; }
// ============================================================================
Ostap::Functions::FuncRooTH2*
Ostap::Functions::FuncRooTH2::clone ( const char* /* name */  ) const
{ return  new FuncRooTH2 ( *this ) ; }
// ============================================================================
Ostap::Functions::FuncRooTH3*
Ostap::Functions::FuncRooTH3::clone ( const char* /* name */  ) const
{ return  new FuncRooTH3 ( *this ) ; }
// ============================================================================


// ============================================================================
/*  constructor from the formula expression 
 *  @param expression the formula expression 
 *  @param tree       the tree 
 *  @param name       the name for the formula 
 */
// ============================================================================
Ostap::Functions::Expression::Expression
( const std::string& expression , 
  const TTree*       tree       ,
  const std::string& name       ) 
  :  FuncFormula ( expression , tree    , name ) 
  ,  m_roofun    ( expression , nullptr , name )
{}
// ============================================================================
/*  constructor from the formula expression 
 *  @param expression the formula expression 
 *  @param data       the data 
 *  @param name       the name for the formula 
 */
// ============================================================================
Ostap::Functions::Expression::Expression
( const std::string& expression , 
  const RooAbsData*  data       ,
  const std::string& name       ) 
  :  FuncFormula ( expression , nullptr , name ) 
  ,  m_roofun    ( expression , data    , name )
{}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Functions::Expression::Expression 
( const Ostap::Functions::Expression& right )
  : FuncFormula ( right          ) 
  , m_roofun    ( right.m_roofun ) 
{}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Functions::Expression::~Expression(){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Functions::Expression* 
Ostap::Functions::Expression::clone ( const char* /* newname */ ) const 
{ return new Ostap::Functions::Expression (*this ) ; }
// ============================================================================
Ostap::Functions::Expression* 
Ostap::Functions::Expression::Clone ( const char* /* newname */ ) const 
{ return new Ostap::Functions::Expression (*this ) ; }
// ============================================================================
// evaluate the function from TTree
// ============================================================================
double Ostap::Functions::Expression::operator () 
  ( const TTree* tree ) const 
{
  // const Ostap::Functions::FuncFormula& formula = *this ;
  // return formula ( tree ) ; 
  return Ostap::Functions::FuncFormula::operator() ( tree ) ;
}
// ============================================================================
// evaluate the function from RooAbsData 
// ============================================================================
double Ostap::Functions::Expression::operator () 
  ( const RooAbsData* data ) const { return m_roofun ( data ) ; }
// ============================================================================



// ============================================================================
/*  full constructor 
 *  @param fun the functon 
 *  @param observables observables 
 *  @param normalization nornalization 
 *  @Param mapping  RooFit varibale <-> TTree branch mapping 
 *  @param tree input tree 
 */
// ============================================================================
Ostap::Functions::RooTreeFun::RooTreeFun
( const RooAbsReal&       fun           , 
  const RooAbsData&       observables   , 
  const RooAbsCollection* normalization ,
  const DCT&              mapping       ,
  const TTree*            tree          )
  : RooTreeFun ( fun , *observables.get() , normalization , mapping , tree )
{}
// ============================================================================
/* full constructor 
 *  @param fun the functon 
 *  @param observables observables 
 *  @Param mapping  RooFit varibale <-> TTree branch mapping 
 *  @param tree input tree 
 */
// ============================================================================
Ostap::Functions::RooTreeFun::RooTreeFun
( const RooAbsReal&       fun         , 
  const RooAbsData&       observables , 
  const DCT&              mapping     ,       
  const TTree*            tree        )
  : RooTreeFun ( fun , *observables.get() , mapping , tree )
{}
// ============================================================================
/*  full constructor 
 *  @param fun the functon 
 *  @param observables observables 
 *  @Param mapping  RooFit varibale <-> TTree branch mapping 
 */
// ============================================================================
Ostap::Functions::RooTreeFun::RooTreeFun
( const RooAbsReal&       fun           , 
  const RooAbsCollection& observables   , 
  const DCT&              mapping       ,
  const TTree*            tree          ) 
  : RooTreeFun ( fun , observables , nullptr , mapping , tree )
{}
// ============================================================================
/** full constructor 
 *  @param fun the functon 
 *  @param observables observables 
 *  @param normalization nornalization 
 *  @Param mapping  RooFit varibale <-> TTree branch mapping 
 *  @param tree input tree 
 */
// ============================================================================
Ostap::Functions::RooTreeFun::RooTreeFun
( const RooAbsReal&       fun           , 
  const RooAbsCollection& observables   , 
  const RooAbsCollection* normalization ,
  const DCT&              mapping       ,
  const TTree*            tree          )
  : Ostap::IFuncTree ()
  , Ostap::Trees::RooGetter ( mapping , tree )
  , m_fun  ( fun , observables , normalization ) 
{
  // ==========================================================================
#if ROOT_VERSION_CODE < ROOT_VERSION(6,18,0) // ===============================
  // ==========================================================================
  //
  Ostap::Utils::Iterator iter( this->observables () ) ; // only for ROOT < 6.18
  RooAbsArg* o = 0 ;
  while ( o = (RooAbsArg*) iter .next() )
    {
      // ======================================================================
#else // ======================================================================
      // ======================================================================
      for ( auto* o : this->observables() )
    {
      // ======================================================================
#endif// ======================================================================
      // ======================================================================
      if ( nullptr != o ) { add ( o->GetName() , o->GetName() ) ; } 
    }
  //
}
// ============================================================================
// copy constructor
// ============================================================================ 
Ostap::Functions::RooTreeFun::RooTreeFun
( const Ostap::Functions::RooTreeFun& right )
  : Ostap::IFuncTree        ( right )
  , Ostap::Trees::RooGetter ( right )
  , m_fun                   ( right.m_fun ) 
{}
// ============================================================================ \
// clone 
// ============================================================================
Ostap::Functions::RooTreeFun*
Ostap::Functions::RooTreeFun::clone ( const char* /* newname */ ) const
{ return new RooTreeFun ( *this ) ; }
// ============================================================================
Ostap::Functions::RooTreeFun*
Ostap::Functions::RooTreeFun::Clone ( const char* /* newname */ ) const
{ return new RooTreeFun ( *this ) ; }
// ============================================================================
// the main method
// ============================================================================
double Ostap::Functions::RooTreeFun::operator()
  ( const TTree* tree ) const
{
  const Ostap::StatusCode sc = assign ( m_fun.observables()  , tree ) ;
  Ostap::Assert ( sc.isSuccess ()      ,
                  "Invaild RooGetter!" ,
                  "Ostap::Functions::RooTreeFun:()" , sc , __FILE__ , __LINE__ ) ;
  //
  return m_fun.evaluate() ;
}




// ============================================================================
//                                                                      The END 
// ============================================================================
