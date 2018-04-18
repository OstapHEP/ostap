// ============================================================================
// Include files
// ============================================================================
// local
// ============================================================================
#include "Ostap/Funcs.h"
#include "Ostap/Formula.h"
#include "Ostap/Iterator.h"
#include "Ostap/StatusCode.h"
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
  TTree* t = const_cast<TTree*> ( m_tree ) ;
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



// ============================================================================
ClassImp(Ostap::Functions::FuncFormula) 
// ============================================================================
// The END 
// ============================================================================
