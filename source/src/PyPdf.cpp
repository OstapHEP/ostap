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
 *  @see  Ostap::Models::PyPdfLite
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
( const std::string& name        , 
  const std::string& title       ,
  const RooArgList&  variables   ) 
  : RooAbsPdf   ( name.c_str() , title.c_str() )
  , m_varlist ( "!varlist" , "All variables(list)" , this ) 
{ ::copy_real ( variables , m_varlist ) ; }
// ============================================================================
/*  Standard constructor
 *  @param name      the name of PDF 
 *  @param title     the title  of PDF 
 *  @param observables observables 
 *  @param parameters  parameters 
 */
// ============================================================================
Ostap::Models::PyPdf::PyPdf
( const std::string& name        ,   
  const std::string& title       ,
  const RooArgList&  observables , 
  const RooArgList&  parameters  ) 
  : RooAbsPdf ( name.c_str() , title.c_str() )
  , m_varlist ( "!varlisy" , "All variables(list)" , this ) 
{
  ::copy_real ( observables , m_varlist ) ;
  ::copy_real ( parameters  , m_varlist ) ;
}
// ============================================================================
// helper function to be redefined in python  
// ============================================================================
int    Ostap::Models::PyPdf::get_analytical_integral () const { return 0 ; }
// ============================================================================
// helper function to be redefined in python  
// ============================================================================
double Ostap::Models::PyPdf::analytical_integral     () const 
{
  Ostap::throwException ( "Method `analytical_integral' *MUST* be overriden!" ,
                          "Ostap::Models::PyPdf"                              ,
                          UNDEFINED_METHOD ,  __FILE__ , __LINE__ ) ;
  return 0 ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::PyPdf::PyPdf
( const Ostap::Models::PyPdf& right , 
  const char*                 name  ) 
  : RooAbsPdf   ( right , name ) 
  , m_varlist ( "!varlist" , this , right.m_varlist ) 
{}
// ============================================================================
// virtual destructor 
// ============================================================================
Ostap::Models::PyPdf::~PyPdf() {}
// ============================================================================
Ostap::Models::PyPdf* 
Ostap::Models::PyPdf::clone
( const char* name ) const 
{  
  Ostap::throwException ( "Clone method *MUST* be overridden!"   ,  
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
// move the function from protected to public interface 
// ============================================================================
Bool_t Ostap::Models::PyPdf::match_args
( const RooArgSet& vars ) const 
{ 
  return 
    nullptr != m_allDeps  && 
    nullptr != m_analDeps && 
    RooAbsReal::matchArgs ( *m_allDeps , *m_analDeps , vars ) ;
}
// ============================================================================
// move the function from protected to public interface 
// ============================================================================
Bool_t Ostap::Models::PyPdf::match_arg
( const RooAbsArg& var ) const 
{ 
  return 
    nullptr != m_allDeps  && 
    nullptr != m_analDeps && 
    RooAbsReal::matchArgs ( *m_allDeps , *m_analDeps , var ) ;
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
double Ostap::Models::PyPdf::value ( const unsigned short index ) const 
{
  Ostap::Assert ( index < m_varlist.getSize() , 
                  "Invalid index"             , 
                  "PyPdf::value(index)"       ,
                  INVALID_VARIABLE  , __FILE__ , __LINE__ ) ;
  //
  const RooAbsArg*  a = m_varlist.at ( index ) ;
  Ostap::Assert ( a                        , 
                  "Invalid element"        , 
                  "PyPdf::value(index)"    ,
                  INVALID_VARIABLE  , __FILE__ , __LINE__ ) ;
  //
  const RooAbsReal* v = dynamic_cast<const RooAbsReal*>( a ) ;
  Ostap::Assert ( v                         , 
                  "Invalid element type"    , 
                  "PyPdf::value(index)"     , 
                  INVALID_VARIABLE  , __FILE__ , __LINE__ ) ;
  return v->getVal() ;  
}
// ============================================================================
// get a variable with name 
// ============================================================================
double Ostap::Models::PyPdf::value ( const char* name  ) const 
{
  const RooAbsArg*  a = m_varlist.find ( name  ) ;
  Ostap::Assert ( a                        , 
                  "Invalid element"        , 
                  "PyPdf::value(name)"     , 
                  INVALID_VARIABLE  , __FILE__ , __LINE__ ) ;
  const RooAbsReal* v = dynamic_cast<const RooAbsReal*>( a ) ;
  Ostap::Assert ( v                        , 
                  "Invalid element type"   , 
                  "PyPdf::value(name)"     , 
                  INVALID_VARIABLE  , __FILE__ , __LINE__ ) ;
  return v->getVal() ;  
}
// ============================================================================


// ============================================================================
// PyPdfLite 
// ============================================================================
/*  Standard constructor
 *  @param self      python-partner for this C++ instance 
 *  @param name      the name of PDF 
 *  @param title     the title  of PDF 
 *  @param variables all variables 
 */
// ============================================================================
Ostap::Models::PyPdfLite::PyPdfLite
( const char*       name      , 
  const char*       title     ,
  PyObject*         function  , 
  const RooArgList& variables )
  : RooAbsPdf  ( name     , title  ) 
  , m_function ( function )
  , m_varlist  ( "!varlist" , "All variables(list)" , this ) 
{
  //
  Ostap::Assert ( m_function && PyCallable_Check ( m_function ) , 
                  "Invalid py-function"                         , 
                  "Ostap::Models::PyPdfLite"                    ,
                  INVALID_CALLABLE , __FILE__ , __LINE__        ) ;
  //
  ::copy_real ( variables , m_varlist ) ;
  //
  if ( m_function ) { Py_XINCREF ( m_function ) ; }
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::PyPdfLite::PyPdfLite
( const Ostap::Models::PyPdfLite& right , 
  const char*                     name  ) 
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
Ostap::Models::PyPdfLite::~PyPdfLite() 
{
  if ( m_function ) { Py_DECREF ( m_function  ) ; m_function  = nullptr ; }
}
// ============================================================================
std::size_t Ostap::Models::PyPdfLite::numrefs   () const
{ return nullptr == m_function ? 0 : Py_REFCNT ( m_function ) ; }
// ============================================================================
/*  get the underlyaing function 
 *  @attention referenc odut is incremented!
 */
// ============================================================================
const PyObject*
Ostap::Models::PyPdfLite::function  () const
{
  if ( m_function ) { Py_XINCREF ( m_function  ) ; }  
  return m_function ;
}
// =============================================================================
Ostap::Models::PyPdfLite* 
Ostap::Models::PyPdfLite::clone ( const char* name ) const 
{ return new Ostap::Models::PyPdfLite ( *this , name ) ; }
// ============================================================================
// the actual evaluation of function
// ============================================================================
Double_t Ostap::Models::PyPdfLite::evaluate() const 
{
  // 
  if  ( 0 == m_function || !PyCallable_Check( m_function ) ) 
    {
      PyErr_Print() ;
      Ostap::throwException ( "Function is not callable/invalid"     ,
                              "Ostap::Models::PyPdfLite::evaluate"   ,
                              INVALID_CALLABLE , __FILE__ , __LINE__ ) ;
    }
  //
  //
  PyObject* arguments = PyTuple_New ( m_varlist.getSize() ) ;
  if ( !arguments )
    {
      PyErr_Print () ;
      Ostap::throwException ( "Can't create PyTuple"               ,
                              "Ostap::Models::PyPdfLite::evaluate" ,
                              ERROR_PYTHON , __FILE__ , __LINE__   ) ;
      return 0 ;
    }
  //
  unsigned short index = 0 ;
  //
  for  ( auto* av : m_varlist )
    {
      Ostap::Assert ( nullptr != av                         ,
                      "Invalid RooAbsArg"                   ,
                      "Ostap::Models::PyPdfLite::evaluate"  ,
                      INVALID_ABSARG , __FILE__ , __LINE__  ) ;
      //
      RooAbsReal* v = static_cast<RooAbsReal*> ( av ) ;
      Ostap::Assert ( nullptr != v                          ,
                      "Invalid RooAbsReal"                  ,
                      "Ostap::Models::PyPdfLite::evaluate"  ,
                      INVALID_ABSREAL , __FILE__ , __LINE__ ) ;
      //
      PyObject* pv =  PyFloat_FromDouble ( v->getVal()  ) ;
      if ( 0 != PyTuple_SetItem ( arguments , index , pv ) ) 
        {
          PyErr_Print () ;
          Py_XDECREF ( arguments ) ; arguments = nullptr ;
          Ostap::throwException ( "Can't fill PyTuple"                 ,
                                  "Ostap::Models::PyPdfLite::evaluate" ,
                                  ERROR_PYTHON , __FILE__ , __LINE__   ) ;
        }
      ++index ;
    }
  //
  PyObject* result = PyObject_CallObject ( m_function , arguments ) ;
  //
  Py_XDECREF ( arguments ) ;
  //
  return result_to_double ( result , "PyPdfLite::evaluate" ) ;
  }
// ============================================================================

// ============================================================================
// needed ? 
// ============================================================================
ClassImp(Ostap::Models::PyPdf)
ClassImp(Ostap::Models::PyPdfLite)
// ============================================================================
//                                                                      The END
// ============================================================================
