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
// ============================================================================
/** @file 
 *  Implementation file for class Ostap::Models::PyPdf and Ostap::Functions::PyPdf 
 *  @see  Ostap::Models::PyPd
 *  @see  Ostap::Functions::PyPdf 
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
( PyObject*         self      , 
  const char*       name      , 
  const char*       title     ,
  const RooArgList& variables )
  : RooAbsPdf (   name , title  ) 
  , m_self ( self )
  , m_varlist ( "!varlist" , "All variables(list)" , this ) 
{
  //
  Ostap::Assert ( m_self , 
                  "self* points to NULL" , 
                  "PyPdf::consructor"    , 
                  Ostap::StatusCode(400) ) ;
  //
  Ostap::Utils::Iterator it ( variables ) ;
  while ( RooAbsReal* v = it.static_next<RooAbsReal>() ) { m_varlist.add ( *v ) ; } 
  //
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
{
  Py_XINCREF ( m_self ) ;
}
// ============================================================================
// virtual destructor 
// ============================================================================
Ostap::Models::PyPdf::~PyPdf() 
{
  if ( m_self ) 
  {
    const int err1 = PyObject_SetAttrString ( m_self , "pdf"   , Py_None ) ;
    if  ( 0 != err1 ) { PyErr_Print(); }
    const int err2 = PyObject_SetAttrString ( m_self , "pypdf" , Py_None ) ;
    if  ( 0 != err2 ) { PyErr_Print(); }
    Py_DECREF ( m_self ) ; 
  }
}
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
                              PyString_FromString  ( "pypdf" )   ,
#else 
                              PyUnicode_FromString ( "pypdf" )   ,
#endif 
                              pycl                               ) )
  {
    PyErr_Print();
    Py_DECREF ( method ) ; 
    Py_DECREF ( kwargs ) ;
    Py_DECREF ( pycl   ) ;
    Ostap::throwException ( "Can't set ``pypdf'' item"       ,
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
  const double result = call_method ( m_self , s_AI ) ;
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
{ return call_method ( m_self , s_evaluate ) ; }
// ============================================================================
// get a variable with index 
// ============================================================================
double Ostap::Models::PyPdf::variable ( const unsigned short index ) const 
{
  Ostap::Assert ( index < m_varlist.getSize() , 
                  "Invalid index"          , 
                  "PyPdf::variable(index)" , 
                  Ostap::StatusCode ( 805 ) ) ;
  //
  const RooAbsArg*  a = m_varlist.at ( index ) ;
  Ostap::Assert ( a , 
                  "Invalid element"        , 
                  "PyPdf::variable(index)" , 
                  Ostap::StatusCode ( 806 ) ) ;
  const RooAbsReal* v = static_cast<const RooAbsReal*>( a ) ;
  Ostap::Assert ( a , 
                  "Invalid element type"    , 
                  "PyPdf::variable(index)"  , 
                  Ostap::StatusCode ( 807 ) ) ;
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
                  Ostap::StatusCode ( 808 ) ) ;
  const RooAbsReal* v = static_cast<const RooAbsReal*>( a ) ;
  Ostap::Assert ( a , 
                  "Invalid element type"    , 
                  "PyPdf::variable(name)"  , 
                  Ostap::StatusCode ( 809 ) ) ;
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
                  "self* points to NULL" , 
                  "PyPdf2::consructor"    , 
                  Ostap::StatusCode(400) ) ;
  //
  Ostap::Utils::Iterator it ( variables ) ;
  while ( RooAbsReal* v = it.static_next<RooAbsReal>() ) { m_varlist.add ( *v ) ; } 
  //
  if ( m_function ) { Py_XINCREF ( m_function ) ; }
  m_arguments = PyTuple_New ( m_varlist.getSize() ) ;
}
// ============================================================================
// copy constructor
// ============================================================================
Ostap::Models::PyPdf2::PyPdf2
( const Ostap::Models::PyPdf2& right , 
  const char*                 name  ) 
  : RooAbsPdf  ( right , name     ) 
    //
  , m_function ( right.m_function ) 
  , m_varlist  ( "!varlist" , this , right.m_varlist ) 
{
  if ( m_function ) { Py_XINCREF ( m_function  ) ; }
  m_arguments = PyTuple_New ( m_varlist.getSize() ) ;
}
// ============================================================================
// virtual destructor 
// ============================================================================
Ostap::Models::PyPdf2::~PyPdf2() 
{
  if ( m_function  ) { Py_DECREF ( m_function  ) ; m_function  = nullptr ; }
  if ( m_arguments ) { Py_DECREF ( m_arguments ) ; m_arguments = nullptr ; }
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
 Ostap:Assert ( PySequence_Size ( m_arguments ) == m_varlist.getSize() , 
                "Invalid argument/varlist  size!" ,
                "PyPdf2::evaluate"                 ,
                Ostap::StatusCode(500) ) ;
  //
  Ostap::Utils::Iterator it ( m_varlist ) ;
  unsigned short index = 0 ;
  while ( RooAbsReal* v = it.static_next<RooAbsReal>() )
  {
    PyObject* pv =  PyFloat_FromDouble ( v->getVal()  ) ;
    if ( 0 != PyTuple_SetItem ( m_arguments , index , pv ) ) 
    {
      PyErr_Print () ;
      Ostap::throwException ( "Can't fill PyTuple"   ,
                              "PyPdf2::evaluate"     ,
                              Ostap::StatusCode(500) ) ;
    }
    ++index ;
  }
  //
  PyObject* result = PyObject_CallObject ( m_function , m_arguments ) ;
  //
  return result_to_double ( result , "PyPdf2::evaluate" ) ;
}
// ============================================================================


// ============================================================================
// PyVar
// ============================================================================


// ============================================================================
// Standard constructor
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
  Ostap::Utils::Iterator it ( variables ) ;
  while ( RooAbsReal* v = it.static_next<RooAbsReal>() ) { m_variables.add ( *v ) ; } 
  //
  Py_XINCREF ( m_self ) ;  
}
// =============================================================================
// Copy constructors 
// =============================================================================
Ostap::Functions::PyVar::PyVar 
( const Ostap::Functions::PyVar& right , const char* newname ) 
  : RooAbsReal  ( right , newname )
  , m_self      ( right.m_self )
  , m_variables ( "variables"  , this , right.m_variables )
{
  Py_XINCREF ( m_self ) ;  
}
// =============================================================================
// virtual destructor
// =============================================================================
Ostap::Functions::PyVar::~PyVar() { Py_DECREF ( m_self ) ; }
// ============================================================================
//  Clone method 
// ============================================================================
Ostap::Functions::PyVar* 
Ostap::Functions::PyVar::clone ( const char* name ) const 
{
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
  PyObject*  pycl = TPython::ObjectProxy_FromVoidPtr ( cl , cl->IsA()->GetName() , false) ;  
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
{ return call_method ( m_self , s_evaluate ) ; }
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
  Ostap::Utils::Iterator it ( variables ) ;
  while ( RooAbsReal* v = it.static_next<RooAbsReal>() ) { m_variables.add ( *v ) ; } 
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
  Ostap::Utils::Iterator it ( m_variables ) ;
  unsigned short index = 0 ;
  while ( RooAbsReal* v = it.static_next<RooAbsReal>() )
  {
    PyObject* pv =  PyFloat_FromDouble ( v->getVal()  ) ;
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
// needed ? 
// ============================================================================
ClassImp(Ostap::Models::PyPdf)
ClassImp(Ostap::Models::PyPdf2)
ClassImp(Ostap::Functions::PyVar)
ClassImp(Ostap::Functions::PyVar2)


// ============================================================================
//                                                                      The END
// ============================================================================
