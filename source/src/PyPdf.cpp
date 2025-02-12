// ============================================================================
// Include files
// ============================================================================
// Python
// ============================================================================
#include "Python.h"
// ============================================================================
// STD&STL
// ============================================================================
#include <cstring>
// ============================================================================
// ROOT 
// ============================================================================
#include "TPython.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/PyPdf.h"
#include "Ostap/PyVar.h"
#include "Ostap/Iterator.h"
// ============================================================================
// Local
// ============================================================================
#include "CallPython.h"
#include "local_roofit.h"
#include "status_codes.h"
// ============================================================================
/** @file 
 *  Implementation file for class Ostap::Models::PyPdf 
 *  @see  Ostap::Models::PyPdf
 *  @see  Ostap::Models::PyPdf2
 *  @date 2018-06-06 
 *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  static char s_evaluate[] = "evaluate"                 ;
  static char s_clone   [] = "clone"                    ;
  static char s_getAI   [] = "get_analytical_integral"  ;
  static char s_AI      [] = "analytical_integral"      ;
  // ==========================================================================
}
// ============================================================================
/*  Standard constructor
 *  @param name      the name of PDF 
 *  @param title     the title  of PDF 
 *  @param variables all variables 
 */
// ============================================================================
Ostap::Models::PyPdf::PyPdf
( const char*       name      , 
  const char*       title     ,
  const RooArgList& variables )
  : RooAbsPdf ( name , title  ) 
  , m_varlist ( "!varlist" , "All variables(list)" , this ) 
{
  //
  ::copy_real ( variables , m_varlist ) ;
  //
}
// ============================================================================


// ============================================================================
// helper function to be redefined in python  
// ============================================================================
int    Ostap::Models::PyPdf::get_analytical_integral () const { return 0 ; }
// ============================================================================
// helper function to be redefined in python  
// ============================================================================
double Ostap::Models::PyPdf::analytical_integral     () const 
{
  Ostap::throwException ( "Method `analytical_integral' must be overriden!"  ,
                          "Ostap::Models::PyPdf"                             ,
                          UNDEFINED_METHOD ,  __FILE__ , __LINE__ ) ;
  return 0 ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::PyPdf::PyPdf
( const Ostap::Models::PyPdf& right , 
  const char*                 name  ) 
  : RooAbsPdf ( right , name ) 
  , m_varlist  ( "!varlist" , this , right.m_varlist ) 
{}
// ============================================================================


// ============================================================================
// virtual destructor 
// ============================================================================
Ostap::Models::PyPdf::~PyPdf() {}
// ============================================================================
Ostap::Models::PyPdf* 
Ostap::Models::PyPdf::clone ( const char* name ) const 
{  
  Ostap::throwException ( "Clone method MUST be overridden!"     ,  
                          "Ostap::Functions::PyPdf"              , 
                          UNDEFINED_METHOD , __FILE__ , __LINE__ ) ;
  return nullptr  ;  
  // ============================================================================
}
// ============================================================================
// declare the analysitical integrals 
// ============================================================================
Int_t Ostap::Models::PyPdf::getAnalyticalIntegral
( RooArgSet&  allVars    ,
  RooArgSet&  analVars   ,
  const char* rangeName  ) const 
{
  //
  ///
  m_allDeps   = &allVars  ;
  m_analDeps  = &analVars ;
  m_rangeName = rangeName ;
  m_intCode   = 0         ;
  //  
  // ==========================================================================
  const int code = get_analytical_integral() ;
  // ==========================================================================
  m_allDeps   = nullptr   ;
  m_analDeps  = nullptr   ;
  m_rangeName = nullptr   ;
  m_intCode   = 0         ;
  // ==========================================================================  
  return code ;
  // ==========================================================================
}
// ============================================================================
// get the integral 
// ============================================================================
Double_t Ostap::Models::PyPdf::analyticalIntegral 
( Int_t       code      , 
  const char* rangeName ) const 
{
  //
  m_intCode    = code      ;
  m_rangeName  = rangeName ;
  ///
  // ==========================================================================
  const double result =   analytical_integral() ;
  // ==========================================================================
  //
  m_intCode    = 0         ;
  m_rangeName  = nullptr   ;
  //
  return result ;
}
// ============================================================================
/// move the function from protected to public interface 
Bool_t Ostap::Models::PyPdf::matchArgs ( const RooArgSet& refVars ) const 
{ 
  return 
    nullptr != m_allDeps  && 
    nullptr != m_analDeps && 
    RooAbsReal::matchArgs ( *m_allDeps , *m_analDeps , refVars ) ;
}
// ============================================================================
// the actual evaluation of function
// ============================================================================
Double_t Ostap::Models::PyPdf::evaluate() const 
{ 
  Ostap::throwException ( "evaluate method must be overrided!"   , 
                          "Ostap::Functions::PyPdf"              ,
                          UNDEFINED_METHOD , __FILE__ , __LINE__ ) ;
  return -1 ;
}
// ============================================================================
// get a variable with index 
// ============================================================================
double Ostap::Models::PyPdf::variable ( const unsigned short index ) const 
{
  Ostap::Assert ( index < m_varlist.getSize() , 
                  "Invalid index"           , 
                  "PyPdf::variable(index)"  ,
                  INVALID_VARIABLE  , __FILE__ , __LINE__ ) ;
  //
  const RooAbsArg*  a = m_varlist.at ( index ) ;
  Ostap::Assert ( a , 
                  "Invalid element"        , 
                  "PyPdf::variable(index)" ,
                  INVALID_VARIABLE  , __FILE__ , __LINE__ ) ;
  //
  const RooAbsReal* v = dynamic_cast<const RooAbsReal*>( a ) ;
  Ostap::Assert ( v , 
                  "Invalid element type"    , 
                  "PyPdf::variable(index)"  , 
                  INVALID_VARIABLE  , __FILE__ , __LINE__ ) ;
  return v->getVal() ;  
}
// ============================================================================
// get a variable with name 
// ============================================================================
double Ostap::Models::PyPdf::variable ( const char* name  ) const 
{
  const RooAbsArg*  a = m_varlist.find ( name  ) ;
  Ostap::Assert ( a , 
                  "Invalid element"        , 
                  "PyPdf::variable(name)"  , 
                  INVALID_VARIABLE  , __FILE__ , __LINE__ ) ;
  const RooAbsReal* v = dynamic_cast<const RooAbsReal*>( a ) ;
  Ostap::Assert ( v , 
                  "Invalid element type"    , 
                  "PyPdf::variable(name)"  , 
                  INVALID_VARIABLE  , __FILE__ , __LINE__ ) ;
  return v->getVal() ;  
}
// ============================================================================


// ============================================================================
// PyPdf2 
// ============================================================================


// ============================================================================
/*  Standard constructor
 *  @param self      python-partner for this C++ instance 
 *  @param name      the name of PDF 
 *  @param title     the title  of PDF 
 *  @param variables all variables 
 */
// ============================================================================
Ostap::Models::PyPdf2::PyPdf2
( const char*       name      , 
  const char*       title     ,
  PyObject*         function  , 
  const RooArgList& variables )
  : RooAbsPdf  (   name , title  ) 
  , m_function ( function )
  , m_varlist  ( "!varlist" , "All variables(list)" , this ) 
{
  //
  Ostap::Assert ( m_function , 
                  "Invalid py-functon"   , 
                  "PyPdf2::consructor"   , 
                  Ostap::StatusCode(400) ) ;
  //
  ::copy_real ( variables , m_varlist ) ;
  //
  if ( m_function ) { Py_XINCREF ( m_function ) ; }
}
// ============================================================================
/*  Standard constructor
 *  @param name      the name of PDF 
 *  @param function  callable function
 *  @param variables all variables 
 *  @param title     the title  of PDF 
 */
// ============================================================================
Ostap::Models::PyPdf2::PyPdf2
( const std::string& name       , 
  PyObject*          function   ,
  const RooArgList&  variables  ,
  const std::string& title      )
  : PyPdf2 ( name.c_str() , 
             title.empty() ? name.c_str() : title.c_str() , 
             function     ,  
             variables    )
{}      
// ============================================================================
/*  Standard constructor
 *  @param name      the name of PDF 
 *  @param function  callable function
 *  @param variables all variables 
 *  @param title     the title  of PDF 
 */
// ============================================================================
Ostap::Models::PyPdf2::PyPdf2
( const std::string& name       , 
  const RooArgList&  variables  ,
  PyObject*          function   ,
  const std::string& title      ) 
  : PyPdf2 ( name , function , variables , title ) 
{}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::PyPdf2::PyPdf2
( const Ostap::Models::PyPdf2& right , 
  const char*                  name  ) 
  : RooAbsPdf  ( right , name     ) 
    //
  , m_function ( right.m_function ) 
  , m_varlist  ( "!varlist" , this , right.m_varlist ) 
{
  if ( m_function ) { Py_XINCREF ( m_function  ) ; }
}
// ============================================================================
// virtual destructor 
// ============================================================================
Ostap::Models::PyPdf2::~PyPdf2() 
{
  if ( m_function ) { Py_DECREF ( m_function  ) ; m_function  = nullptr ; }
}
// ============================================================================
Ostap::Models::PyPdf2* 
Ostap::Models::PyPdf2::clone ( const char* name ) const 
{ return new Ostap::Models::PyPdf2 ( *this , name ) ; }

// ============================================================================
// the actual evaluation of function
// ============================================================================
Double_t Ostap::Models::PyPdf2::evaluate() const 
{
  // 
  if  ( 0 == m_function || !PyCallable_Check( m_function ) ) 
  {
    PyErr_Print() ;
    Ostap::throwException ( "Function is not callable/invalid" ,
                            "PyPdf2::evaluate"                 ,
                            Ostap::StatusCode(500) )           ;
  }
  //
   //
  PyObject* arguments = PyTuple_New ( m_varlist.getSize() ) ;
  if ( !arguments )
  {
    PyErr_Print () ;
    Ostap::throwException ( "Can't create PyTuple"   ,
                            "PyPdf2::evaluate"     ,
                            Ostap::StatusCode(500) ) ;
    return 0 ;
  }
  //
  unsigned short index = 0 ;
  //
  for  ( auto* av : m_varlist )
  {
    if ( nullptr == av ) {continue  ; }
    RooAbsReal* v = static_cast<RooAbsReal*> ( av ) ;
    //  
    if ( nullptr == v  ) { continue ; }
    PyObject* pv =  PyFloat_FromDouble ( v->getVal()  ) ;
    if ( 0 != PyTuple_SetItem ( arguments , index , pv ) ) 
    {
      PyErr_Print () ;
      Py_XDECREF ( arguments ) ; arguments = nullptr ;
      Ostap::throwException ( "Can't fill PyTuple"   ,
                              "PyPdf2::evaluate"     ,
                              Ostap::StatusCode(500) ) ;
    }
    ++index ;
  }
  //
  PyObject* result = PyObject_CallObject ( m_function , arguments ) ;
  //
  Py_XDECREF ( arguments ) ;
  //
  return result_to_double ( result , "PyPdf2::evaluate" ) ;
}
// ============================================================================


// ============================================================================
// needed ? 
// ============================================================================
ClassImp(Ostap::Models::PyPdf)
ClassImp(Ostap::Models::PyPdf2)
// ============================================================================
//                                                                      The END
// ============================================================================
