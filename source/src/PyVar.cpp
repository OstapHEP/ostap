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
 *  @see  Ostap::Functions::PyVarLite 
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
  // ==========================================================================}
  /// defautl (invalid) value 
  const double s_DEFAULT {-1000 } ;   /// default (invalid) value 
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
  , m_varlist ( "!varlist", "The actual variables/parameters" , this ) 
{
  ::copy_real ( variables                          ,
                m_varlist                          ,
                "Not all variables are RooAbsReal" ,
                "Ostap::Functions::PyVa"           ,
                __FILE__ , __LINE__                ) ;
}
// ============================================================================
// Copy constructor
// =============================================================================
Ostap::Functions::PyVar::PyVar 
( const Ostap::Functions::PyVar& right , const char* newname ) 
  : RooAbsReal  ( right , newname )
  , m_varlist ( "!varlist"  , this , right.m_varlist )
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
  Ostap::Assert( false                                  ,
                 "clone method must be implemented!"    , 
                 "Ostap::Functions::PyVar"              ,
                 UNDEFINED_METHOD , __FILE__ , __LINE__ ) ;
  return nullptr ;
}
// ============================================================================
// get a variable with index 
// ============================================================================
double Ostap::Functions::PyVar::value ( const unsigned short index ) const 
{
  Ostap::Assert ( index < m_varlist.getSize()           , 
                  "Invalid index"                         , 
                  "PyVar::value(index)"                   ,
                  INVALID_VARIABLE  , __FILE__ , __LINE__ ) ;
  //
  const RooAbsArg*  a = m_varlist.at ( index ) ;
  Ostap::Assert ( a , 
                  "Invalid element"                       , 
                  "PyVar::value(index)"                   , 
                  INVALID_VARIABLE  , __FILE__ , __LINE__ ) ;
  const RooAbsReal* v = static_cast<const RooAbsReal*>( a ) ;
  Ostap::Assert ( v , 
                  "Invalid element type"                  , 
                  "PyVar::value(index)"                   ,
                  INVALID_VARIABLE  , __FILE__ , __LINE__ ) ;
  return v->getVal() ;  
}
// ============================================================================
// get a variable with name 
// ============================================================================
double Ostap::Functions::PyVar::value ( const char* name  ) const 
{
  const RooAbsArg*  a = m_varlist.find ( name  ) ;
  Ostap::Assert ( a , 
                  "Invalid element"                          , 
                  "PyVar::value " + std::string ( name )     ,
                  INVALID_VARIABLE , __FILE__ , __LINE__     ) ;
  const RooAbsReal* v = dynamic_cast<const RooAbsReal*>( a ) ;
  Ostap::Assert ( v , 
                  "Invalid element type"                     , 
                  "PyVar::value" + std::string ( name )      ,
                  INVALID_VARIABLE , __FILE__ , __LINE__     ) ;
  return v->getVal() ;  
}
// ============================================================================
// the actual evaluation of function
// ============================================================================
Double_t Ostap::Functions::PyVar::evaluate() const 
{ 
  Ostap::Assert ( false ,
                  "evaluate method must be implemented!" , 
                  "Ostap::Functions::PyVar"              ,
                  UNDEFINED_METHOD , __FILE__ , __LINE__ ) ;
  return s_DEFAULT ;
}
// ============================================================================
std::vector<double>
Ostap::Functions::PyVar::get_values() const
{
  std::vector<double> values  {} ; values.reserve ( m_varlist.size() ) ;
  for ( const RooAbsArg* a : m_varlist )
    {
      // it is safe since itis checked in constructur 
      const RooAbsReal* v = static_cast<const RooAbsReal*> ( a ) ;
      values.push_back ( v->getVal() ) ;
    }
  return values ;
}
// ============================================================================
// PyVarLite
// ============================================================================
// Standard constructor
// ============================================================================
Ostap::Functions::PyVarLite::PyVarLite 
( const char*       name      , 
  const char*       title     ,
  PyObject*         function  , 
  const RooArgList& variables ) 
  : RooAbsReal  ( name  , title ) 
  , m_function  ( function ) 
  , m_varlist   ( "!varlist", "The actual variables/parameters" , this ) 
{
  //
  ::copy_real   ( variables                          ,
                  m_varlist                          ,
                  "Not all variables are RooAbsReal" ,
                  "Ostap::Functions::PyVarLite"      ,
                  __FILE__ , __LINE__                ) ; 
  //
  Ostap::Assert ( m_function && PyCallable_Check ( m_function ) , 
                  "Invalid py-function"                         , 
                  "Ostap::Functions::PyVarLite"                 ,
                  INVALID_CALLABLE , __FILE__ , __LINE__        ) ;
  //
  Py_XINCREF ( m_function ) ;
}
// =============================================================================
// Copy constructors 
// =============================================================================
Ostap::Functions::PyVarLite::PyVarLite 
( const Ostap::Functions::PyVarLite& right , const char* newname ) 
  : RooAbsReal ( right , newname  )
  , m_function ( right.m_function )
  , m_varlist  ( "!varlist"  , this , right.m_varlist )
{
  Py_XINCREF ( m_function ) ;  
}
// =============================================================================
// virtual destructor
// =============================================================================
Ostap::Functions::PyVarLite::~PyVarLite() 
{ 
  if ( m_function  ) { Py_DECREF ( m_function  ) ; }  
}
// ============================================================================
//  Clone method 
// ============================================================================
Ostap::Functions::PyVarLite* 
Ostap::Functions::PyVarLite::clone ( const char* name ) const 
{ return new Ostap::Functions::PyVarLite ( *this , name ) ; }
// ============================================================================
std::size_t Ostap::Functions::PyVarLite::numrefs   () const
{ return nullptr == m_function ? 0 : Py_REFCNT ( m_function ) ; }
// ============================================================================
/*  get the underlyaing function 
 *  @attention referenc odut is incremented!
 */
// ============================================================================
const PyObject*
Ostap::Functions::PyVarLite::function  () const
{
  if ( m_function ) { Py_XINCREF ( m_function  ) ; }  
  return m_function ;
}
// ============================================================================
// the actual evaluation of function
// ============================================================================
Double_t Ostap::Functions::PyVarLite::evaluate() const 
{
  // 
  if  ( 0 == m_function || !PyCallable_Check( m_function ) ) 
    {
      PyErr_Print() ;
      Ostap::Assert ( false ,
                      "Function is not callable/invalid"      ,
                      "?Ostap::Fuctions::PyVarLite::evaluate" ,
                      INVALID_CALLABLE , __FILE__ , __LINE__  ) ;
      return s_DEFAULT ;
    }
  //
  PyObject* arguments = PyTuple_New ( m_varlist.getSize() ) ;
  if ( !arguments )
    {
      PyErr_Print () ;
      Ostap::Assert ( false                                   ,
                      "Can't create PyTuple"                  ,
                      "?Ostap::Fuctions::PyVarLite::evaluate" ,
                      ERROR_PYTHON , __FILE__ , __LINE__      ) ;
      return s_DEFAULT ;
    }
  //
  unsigned short index = 0 ;
  //
  for  ( auto* av : m_varlist )
    {
      Ostap::Assert ( nullptr != av                           ,
                      "Invalid RooAbsArg"                     ,
                      "?Ostap::Fuctions::PyVarLite::evaluate" ,
                      INVALID_ABSARG , __FILE__ , __LINE__    ) ;
      //
      RooAbsReal* v = static_cast<RooAbsReal*> ( av ) ;
      Ostap::Assert ( nullptr != v                            ,
                      "Invalid RooAbsReal"                    ,
                      "?Ostap::Fuctions::PyVarLite::evaluate" ,
                      INVALID_ABSREAL , __FILE__ , __LINE__   ) ;
      //
      PyObject* pv =  PyFloat_FromDouble ( v->getVal()  ) ;
      if ( 0 != PyTuple_SetItem ( arguments , index , pv ) ) 
        {
          PyErr_Print () ;
          Py_XDECREF ( arguments ) ; arguments = nullptr ;
          Ostap::Assert ( false ,
                          "Can't fill PyTuple"                    ,
                          "?Ostap::Fuctions::PyVarLite::evaluate" ,
                          ERROR_PYTHON , __FILE__ , __LINE__      ) ;
        }
      ++index ;
    }
  //
  PyObject* result = PyObject_CallObject ( m_function , arguments ) ;
  //
  Py_XDECREF ( arguments ) ;
  //
  return result_to_double ( result , "PyVarLite::evaluate" ) ;
}
// ============================================================================
std::vector<double>
Ostap::Functions::PyVarLite::get_values() const
{
  std::vector<double> values  {} ; values.reserve ( m_varlist.size() ) ;
  for ( const RooAbsArg* a : m_varlist )
    {
      // it is safe since itis checked in constructur 
      const RooAbsReal* v = static_cast<const RooAbsReal*> ( a ) ;
      values.push_back ( v->getVal() ) ;
    }
  return values ;
}


// ============================================================================
ClassImp(Ostap::Functions::PyVar)
ClassImp(Ostap::Functions::PyVarLite)
// ============================================================================

// ============================================================================
//                                                                      The END 
// ============================================================================
