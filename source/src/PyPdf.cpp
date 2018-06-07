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
  static char s_evaluate[] = "evaluate" ;
  static char s_clone   [] = "clone"    ;
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
                              PyString_FromString ( "name" ) ,
                              PyString_FromString ( ( name ? name : "" ) ) ) )
  {
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
                              PyString_FromString ( "pdf" )      ,
                              pycl                               ) )
  {
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
// the actual evaluation of function
// ============================================================================
Double_t Ostap::Models::PyPdf::evaluate() const 
{ return call_python ( m_self , s_evaluate ) ; }
// ============================================================================
//                                                                      The END
// ============================================================================
