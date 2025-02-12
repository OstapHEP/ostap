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
#include "status_codes.h"
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
// Copy constructor
// =============================================================================
Ostap::Functions::PyVar::PyVar 
( const Ostap::Functions::PyVar& right , const char* newname ) 
  : RooAbsReal  ( right , newname )
  , m_variables ( "variables"  , this , right.m_variables )
{}
// ============================================================================
// virtual destructor
// =============================================================================
Ostap::Functions::PyVar::~PyVar() {}
// ============================================================================
//  Clone method 
// ============================================================================
Ostap::Functions::PyVar* 
Ostap::Functions::PyVar::clone ( const char* name ) const 
{
  Ostap::throwException ( "clone method must be implemented!"    , 
                          "Ostap::Functions::PyVar"              ,
                          UNDEFINED_METHOD , __FILE__ , __LINE__ ) ;
  return nullptr ;
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
                  "Invalid element"                          , 
                  "PyVar::variable " + std::string ( name ) ,
                  INVALID_VARIABLE , __FILE__ , __LINE__     ) ;
  const RooAbsReal* v = dynamic_cast<const RooAbsReal*>( a ) ;
  Ostap::Assert ( v , 
                  "Invalid element type" , 
                  "PyVar::variable " + std::string ( name ) ,
                  Ostap::StatusCode ( 804 ) ) ;
  return v->getVal() ;  
}
// ============================================================================
// the actual evaluation of function
// ============================================================================
Double_t Ostap::Functions::PyVar::evaluate() const 
{ 
  Ostap::throwException ( "evaluate method must be implemented!" , 
                          "Ostap::Functions::PyVar"              ,
                          UNDEFINED_METHOD , __FILE__ , __LINE__ ) ;
  return -1000 ;
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
  for  ( auto* vv : m_variables )
  {
    RooAbsReal* v = static_cast<RooAbsReal*>( vv ) ;
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
