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
#include "RVersion.h"
#include "TPython.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/PyVar.h"
// ============================================================================
// Local
// ============================================================================
#include "CallPython.h"
#include "Exception.h"
#include "local_roofit.h"
// ============================================================================
/** @file 
 *  Implementation file for classes Ostap::Functions::PyVar 
 *  @see  Ostap::Functions::PyVar 
 *  @see  Ostap::Functions::PyVar2 
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
// Standard constructor
// ============================================================================
#if defined(OSTAP_OLD_PYROOT) && OSTAP_OLD_PYROOT
// ============================================================================
Ostap::Functions::PyVar::PyVar 
(  PyObject*         self      , 
   const char*       name      , 
   const char*       title     ,
   const RooArgList& variables ) 
  : RooAbsReal  ( name  , title ) 
  , m_self      (  self ) 
  , m_variables ( "variables", "The actual variables/parameters" , this ) 
{
  //
  ::copy_real ( variables , m_variables , "Variable is not RooAbsReal" , "Ostap::Functions::PyVar::PyVar" );
  //
  Py_XINCREF ( m_self ) ;  
}
// ============================================================================
#else 
// ============================================================================
Ostap::Functions::PyVar::PyVar 
( const char*       name      , 
  const char*       title     ,
  const RooArgList& variables ) 
  : RooAbsReal  ( name  , title ) 
  , m_variables ( "variables", "The actual variables/parameters" , this ) 
{
  //
  ::copy_real ( variables , m_variables , "Variable is not RooAbsReal" , "Ostap::Functions::PyVar::PyVar" );
  //
}
// ============================================================================
#endif  
// ============================================================================

// =============================================================================
// Copy constructor
// =============================================================================
#if defined(OSTAP_OLD_PYROOT) && OSTAP_OLD_PYROOT
// =============================================================================
Ostap::Functions::PyVar::PyVar 
( const Ostap::Functions::PyVar& right , const char* newname ) 
  : RooAbsReal  ( right , newname )
  , m_self      ( right.m_self )
  , m_variables ( "variables"  , this , right.m_variables )
{
  Py_XINCREF ( m_self ) ;  
}
// ============================================================================
#else
// ============================================================================
Ostap::Functions::PyVar::PyVar 
( const Ostap::Functions::PyVar& right , const char* newname ) 
  : RooAbsReal  ( right , newname )
  , m_variables ( "variables"  , this , right.m_variables )
{}
// ============================================================================
#endif  
// ============================================================================

// =============================================================================
// virtual destructor
// =============================================================================
Ostap::Functions::PyVar::~PyVar() 
{
#if defined(OSTAP_OLD_PYROOT) && OSTAP_OLD_PYROOT
  Py_DECREF ( m_self ) ; 
#endif 
}
// ============================================================================
//  Clone method 
// ============================================================================
Ostap::Functions::PyVar* 
Ostap::Functions::PyVar::clone ( const char* name ) const 
{
  // ==========================================================================
#if defined(OSTAP_OLD_PYROOT) && OSTAP_OLD_PYROOT
  // ==========================================================================
  /// create the python clone  
  PyObject* method = PyObject_GetAttrString ( m_self , s_clone ) ;
  if ( !method ) 
  {
    PyErr_Print();
    Ostap::throwException ( "No method ``clone'' is found"  ,
                            "PyVar::clone"                  ,
                            Ostap::StatusCode(500)          ) ;
  }
  if  ( !PyCallable_Check ( method ) ) 
  {
    PyErr_Print();
    Py_DECREF ( method ); 
    Ostap::throwException ( "Attribute ``clone'' is not callable" ,
                            "PyVar::clone"              ,
                            Ostap::StatusCode(500) ) ;
  }
  /// create C++ clone 
  PyVar*     cl     = new Ostap::Functions::PyVar ( *this , name ) ;
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
  /// create "pyvar"-item 
  PyObject*  pycl = 
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
       TPython::CPPInstance_FromVoidPtr ( cl , cl->IsA()->GetName() , false ) ;  
#else
       TPython::ObjectProxy_FromVoidPtr ( cl , cl->IsA()->GetName() , false ) ;  
#endif

  if ( !pycl ) 
  {
    PyErr_Print();
    Py_DECREF ( method ) ; 
    Py_DECREF ( kwargs ) ;
    Ostap::throwException ( "Can't pythonize PyVar instance" ,
                            "PyVar::clone"                   ,
                            Ostap::StatusCode(500)           ) ; 
  }
  if  ( 0 != PyDict_SetItem ( kwargs                           ,   
#if defined (PY_MAJOR_VERSION)  and PY_MAJOR_VERSION < 3
                            PyString_FromString  ( "pyvar" )   ,
#else 
                            PyUnicode_FromString ( "pyvar" )   ,
#endif 
                            pycl                               ) )
  {
    PyErr_Print();
    Py_DECREF ( method ) ; 
    Py_DECREF ( kwargs ) ;
    Py_DECREF ( pycl   ) ;
    Ostap::throwException ( "Can't set ``pyvar'' item"       ,
                            "PyVar::clone"                   ,
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
                            "PyVar::clone"                   ,
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
  // ==========================================================================
#else 
  // ==========================================================================  
  Ostap::throwException ( "clone method must be implemented!" , 
                        "Ostap::Functions::PyVar"  ) ;
  return nullptr ;
// ==========================================================================
#endif
}
// ============================================================================
// get a variable with index 
// ============================================================================
double Ostap::Functions::PyVar::variable ( const unsigned short index ) const 
{
  Ostap::Assert ( index < m_variables.getSize() , 
                  "Invalid index"          , 
                  "PyVar::variable(index)" , 
                  Ostap::StatusCode ( 800 ) ) ;
  //
  const RooAbsArg*  a = m_variables.at ( index ) ;
  Ostap::Assert ( a , 
                  "Invalid element"        , 
                  "PyVar::variable(index)" , 
                  Ostap::StatusCode ( 801 ) ) ;
  const RooAbsReal* v = static_cast<const RooAbsReal*>( a ) ;
  Ostap::Assert ( a , 
                  "Invalid element type"    , 
                  "PyVar::variable(index)"  , 
                  Ostap::StatusCode ( 802 ) ) ;
  return v->getVal() ;  
}
// ============================================================================
// get a variable with name 
// ============================================================================
double Ostap::Functions::PyVar::variable ( const char* name  ) const 
{
  const RooAbsArg*  a = m_variables.find ( name  ) ;
  Ostap::Assert ( a , 
                  "Invalid element"        , 
                  "PyVar::variable(name)"  , 
                  Ostap::StatusCode ( 803 ) ) ;
  const RooAbsReal* v = static_cast<const RooAbsReal*>( a ) ;
  Ostap::Assert ( a , 
                  "Invalid element type"    , 
                  "PyVar::variable(name)"  , 
                  Ostap::StatusCode ( 804 ) ) ;
  return v->getVal() ;  
}
// ============================================================================
// the actual evaluation of function
// ============================================================================
Double_t Ostap::Functions::PyVar::evaluate() const 
{ 
  // ==========================================================================
#if defined(OSTAP_OLD_PYROOT) && OSTAP_OLD_PYROOT
  // ==========================================================================
  return call_method ( m_self , s_evaluate ) ; 
  // ==========================================================================  
#else 
  // ==========================================================================  
  Ostap::throwException ( "evaluate method must be implemented!" , 
                          "Ostap::Functions::PyVar"  ) ;
  return -1000 ;
  // ==========================================================================
#endif
}
// ============================================================================



// ============================================================================
// PyVar2
// ============================================================================


// ============================================================================
// Standard constructor
// ============================================================================
Ostap::Functions::PyVar2::PyVar2 
( const char*       name      , 
  const char*       title     ,
  PyObject*         function  , 
  const RooArgList& variables ) 
  : RooAbsReal  ( name  , title ) 
  , m_function  ( function ) 
  , m_variables ( "variables", "The actual variables/parameters" , this ) 
{
  //
  ::copy_real ( variables , m_variables , "Variable is not RooAbsReal" , "Ostap::Functions::PyVar2::PyVar2" );
  //
  Py_XINCREF ( m_function ) ;
  m_arguments = PyTuple_New ( m_variables.getSize() ) ;
}
// =============================================================================
// Copy constructors 
// =============================================================================
Ostap::Functions::PyVar2::PyVar2 
( const Ostap::Functions::PyVar2& right , const char* newname ) 
  : RooAbsReal  ( right , newname  )
  , m_function  ( right.m_function )
  , m_variables ( "variables"  , this , right.m_variables )
{
  Py_XINCREF ( m_function ) ;  
  m_arguments = PyTuple_New ( m_variables.getSize() ) ;
}
// =============================================================================
// virtual destructor
// =============================================================================
Ostap::Functions::PyVar2::~PyVar2() 
{ 
  if ( m_function  ) { Py_DECREF ( m_function  ) ; }  
  if ( m_arguments ) { Py_DECREF ( m_arguments ) ; } 
}
// ============================================================================
//  Clone method 
// ============================================================================
Ostap::Functions::PyVar2* 
Ostap::Functions::PyVar2::clone ( const char* name ) const 
{ return new Ostap::Functions::PyVar2 ( *this , name ) ; }
// ============================================================================
// the actual evaluation of function
// ============================================================================
Double_t Ostap::Functions::PyVar2::evaluate() const 
{
  // 
  if  ( 0 == m_function || !PyCallable_Check( m_function ) ) 
  {
    PyErr_Print() ;
    Ostap::throwException ( "Function is not callable/invalid" ,
                            "PyVar2::evaluate"                 ,
                            Ostap::StatusCode(500) )           ;
  }
  //
 Ostap:Assert ( PySequence_Size ( m_arguments ) == m_variables.getSize() , 
                "Invalid argument/varlist  size!" ,
                "PyVar2::evaluate"                ,
                Ostap::StatusCode(500) ) ;
  //
  unsigned short index = 0 ;
  //
#if ROOT_VERSION_CODE < ROOT_VERSION(6,18,0)
  //
  Ostap::Utils::Iterator it ( m_variables ) ; // only for ROOT < 6.18 
  while ( RooAbsReal* v = it.static_next<RooAbsReal>() )
  {
    //
#else
    //
  for  ( auto* vv : m_variables )
  {
    RooAbsReal* v = static_cast<RooAbsReal*>( vv ) ;
    //
#endif
    //
    const double value = v->getVal() ;
    PyObject* pv =  PyFloat_FromDouble ( value ) ;
    if ( 0 != PyTuple_SetItem ( m_arguments , index , pv ) ) 
    {
      PyErr_Print () ;
      Ostap::throwException ( "Can't fill PyTuple"   ,
                              "PyVar2::evaluate"     ,
                              Ostap::StatusCode(500) ) ;
    }
    ++index ;
  }
  //
  PyObject* result = PyObject_CallObject ( m_function , m_arguments ) ;
  //
  return result_to_double ( result , "PyVar2::evaluate" ) ;
 }
// ============================================================================


// ============================================================================
ClassImp(Ostap::Functions::PyVar)
ClassImp(Ostap::Functions::PyVar2)
// ============================================================================

// ============================================================================
//                                                                      The END 
// ============================================================================
