// ============================================================================
// Include files
// ============================================================================
// Python
// ============================================================================
#include "Python.h"
// ============================================================================
// ROOT 
// ============================================================================
#include "TPython.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/PyPdf.h"
#include "Ostap/Iterator.h"
// ============================================================================
// Local
// ============================================================================
#include "CallPython.h"
// ============================================================================
/** @file 
 *  Implementation file for class Ostap::Models::PyPdf
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
 *  @param self      python-partner for this C++ instance 
 *  @param name      the name of PDF 
 *  @param title     the title  of PDF 
 *  @param variables all variables 
 */
// ============================================================================
Ostap::Models::PyPdf::PyPdf
( PyObject*               self      , 
  const char*             name      , 
  const char*             title     ,
  const RooAbsCollection& variables )
  : RooAbsPdf (   name , title ) 
  , m_self ( self )
  , m_varlist ( "varlist" , "All variables(list)" , this ) 
  , m_varset  ( "varset"  , "All variables(set)"  , this ) 
{
  //
  Ostap::Assert ( m_self , 
                  "self* points to NULL" , 
                  "PyPdf::consructor"    , 
                  Ostap::StatusCode(400) ) ;
  //
  Ostap::Utils::Iterator it ( variables ) ;
  while ( RooAbsReal* v = it.static_next<RooAbsReal>() )
  {
    m_varlist.add ( *v ) ;
    m_varset .add ( *v ) ; 
  } 
  Py_XINCREF ( m_self ) ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::PyPdf::PyPdf
( const Ostap::Models::PyPdf& right , 
  const char*                 name  ) 
  : RooAbsPdf ( right , name ) 
    //
  , m_self     ( right.m_self ) 
  , m_varlist  ( "!varlist" , this , right.m_varlist ) 
  , m_varset   ( "!varset"  , this , right.m_varset  ) 
{
  Py_XINCREF ( m_self ) ;
}
// ============================================================================
// virtual destructor 
// ============================================================================
Ostap::Models::PyPdf::~PyPdf() { if ( m_self ) { Py_DECREF ( m_self ) ; } }
// ============================================================================
Ostap::Models::PyPdf* 
Ostap::Models::PyPdf::clone ( const char* name ) const 
{
  /// create the python clone  
  PyObject* method = PyObject_GetAttrString ( m_self , s_clone ) ;
  if ( !method ) 
  {
    PyErr_Print();
    Ostap::throwException ( "No method ``clone'' is found"  ,
                            "PyPdf::clone"                  ,
                            Ostap::StatusCode(500)          ) ;
  }
  if  ( !PyCallable_Check ( method ) ) 
  {
    PyErr_Print();
    Py_DECREF ( method ); 
    Ostap::throwException ( "Attribute ``clone'' is not callable" ,
                            "PyPdf::clone"              ,
                            Ostap::StatusCode(500) ) ;
  }
  /// create C++ clone 
  PyPdf*     cl     = new Ostap::Models::PyPdf ( *this , name ) ;
  /// create kwargs 
  PyObject*  kwargs = PyDict_New  (   ) ;
  ///
  /// set "name"-item 
  if  ( 0 != PyDict_SetItem ( kwargs                         ,   
#if defined (PY_MAJOR_VERSION)  and PY_MAJOR_VERSION < 3
                              PyString_FromString  ( "name" ) ,
                              PyString_FromString  ( ( name ? name : "" ) ) ) )
#else 
                              PyUnicode_FromString ( "name" ) ,
                              PyUnicode_FromString ( ( name ? name : "" ) ) ) )
#endif
  {
    PyErr_Print();
    Py_DECREF ( method ) ; 
    Py_DECREF ( kwargs ) ;
    Ostap::throwException ( "Can't set ``name'' item"        ,
                            "PyPdf::clone"                   ,
                            Ostap::StatusCode(500)           ) ;
  }
  /// create "pdf"-item 
  PyObject*  pycl = TPython::ObjectProxy_FromVoidPtr ( cl , cl->IsA()->GetName() , false ) ;  
  if ( !pycl ) 
  {
    PyErr_Print();
    Py_DECREF ( method ) ; 
    Py_DECREF ( kwargs ) ;
    Ostap::throwException ( "Can't pythonize PyPdf instance" ,
                            "PyPdf::clone"                   ,
                            Ostap::StatusCode(500)           ) ; 
  }
  if  ( 0 != PyDict_SetItem ( kwargs                             ,   
#if defined (PY_MAJOR_VERSION)  and PY_MAJOR_VERSION < 3
                              PyString_FromString  ( "pdf" )     ,
#else 
                              PyUnicode_FromString ( "pdf" )     ,
#endif 
                              pycl                               ) )
  {
    PyErr_Print();
    Py_DECREF ( method ) ; 
    Py_DECREF ( kwargs ) ;
    Py_DECREF ( pycl   ) ;
    Ostap::throwException ( "Can't set ``pdf'' item"        ,
                            "PyPdf::clone"                   ,
                            Ostap::StatusCode(500)           ) ;
  }
  /// create args 
  PyObject*  args   = PyTuple_New ( 0 ) ;
  // 
  // create python clone!
  PyObject* pyclone = PyObject_Call ( method , args , kwargs ) ;
  if ( !pyclone ) 
  {
    PyErr_Print();
    Py_DECREF ( method ) ; 
    Ostap::throwException ( "Can't create  python ``clone''" ,
                            "PyPdf::clone"                   ,
                            Ostap::StatusCode(500)           ) ;
  }
  //
  Py_INCREF ( pyclone ) ;
  ///
  PyObject *old = cl->m_self ;
  cl->m_self = pyclone   ;   // the most important line!!!
  //
  Py_DECREF  ( method  ) ;
  if ( old ) { Py_DECREF ( old ) ; }
  //
  return cl ;
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
  Ostap::Assert ( m_self                         , 
                  "self* points to NULL"         , 
                  "PyPdf::getAnalyticalIntegral" , 
                  Ostap::StatusCode(400)         ) ;
  ///
  m_allDeps   = &allVars  ;
  m_analDeps  = &analVars ;
  m_rangeName = rangeName ;
  m_intCode   = 0         ;
  //  
  if ( 1 != PyObject_HasAttrString ( m_self , s_getAI ) ||
       1 != PyObject_HasAttrString ( m_self , s_AI    ) )
  {
    //
    m_allDeps   = nullptr   ;
    m_analDeps  = nullptr   ;
    m_rangeName = nullptr   ;
    m_intCode   = 0         ;
    //
    return 0 ;
  }
  //
  /// create the python method  
  PyObject* getai = PyObject_GetAttrString ( m_self , s_getAI ) ;
  if ( !getai || !PyCallable_Check ( getai ) ) 
  {
    if ( nullptr == getai) { PyErr_Print ()        ; }
    else                   { Py_DECREF   ( getai ) ; }
    //
    m_allDeps   = nullptr   ;
    m_analDeps  = nullptr   ;
    m_rangeName = nullptr   ;
    m_intCode   = 0         ;
    //
    return 0 ;                                                            // RETURN
  }
  /// create args 
  PyObject*  args   = PyTuple_New ( 0 ) ;
  // 
  // call python method 
  PyObject* icode = PyObject_Call ( getai , args , nullptr ) ;
  //
  //
#if defined (PY_MAJOR_VERSION)  and PY_MAJOR_VERSION < 3  
  //
  if ( !icode || !PyInt_Check ( icode ) )   // integer value ? ) 
  {
    if ( nullptr == icode ) { PyErr_Print ()        ; }
    else                    { Py_DECREF   ( icode ) ; }
    //
    Py_DECREF ( getai ) ; 
    //
    Ostap::throwException ( "Can't get proper code"        ,
                            "PyPdf::getAnalyticalIntegral" ,
                            Ostap::StatusCode(500)         ) ;
  }
  //
  const Int_t code = PyInt_AS_LONG( icode );
  Py_DECREF ( icode ) ;
  //
#else 
  //
  if ( !icode || !PyLong_Check ( icode ) )   // integer value ? 
  {
    if ( nullptr == icode ) { PyErr_Print ()        ; }
    else                    { Py_DECREF   ( icode ) ; }
    //
    Py_DECREF ( getai ) ; 
    //
    Ostap::throwException ( "Can't get proper code"        ,
                            "PyPdf::getAnalyticalIntegral" ,
                            Ostap::StatusCode(500)         ) ;
  }
  //
  int overflow = 0 ;
  const long code = PyLong_AsLongAndOverflow ( icode , &overflow );
  if      ( -1 == overflow || 1 ==overflow ) 
  {
    Py_DECREF ( icode );      
    Ostap::throwException ( "Can't get proper code/overflow" ,
                            "PyPdf::getAnalyticalIntegral"   ,
                            Ostap::StatusCode(600) ) ;
  }
  else if ( -1 == code && PyErr_Occurred() ) 
  {
    PyErr_Print();
    Py_DECREF ( icode );      
    Ostap::throwException ( "Can't get proper code" ,
                            "PyPdf::getAnalyticalIntegral"   ,
                            Ostap::StatusCode(700) ) ;
  }
  Py_DECREF ( icode ) ;
  //
#endif 
  //

  ///
  m_allDeps   = nullptr   ;
  m_analDeps  = nullptr   ;
  m_rangeName = nullptr   ;
  m_intCode   = 0         ;
  //
  return code ;
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
  const double result = call_python ( m_self , s_AI ) ;
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
{ return call_python ( m_self , s_evaluate ) ; }
// ============================================================================
// needed ? 
// ============================================================================
ClassImp(Ostap::Models::PyPdf)
// ============================================================================
//                                                                      The END
// ============================================================================
