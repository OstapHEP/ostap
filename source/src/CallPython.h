// ============================================================================
#ifndef CALLPYTHON_H 
#define CALLPYTHON_H 1
// ============================================================================
// Incldue files 
// ============================================================================
// Python
// ============================================================================
#include "Python.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatusCode.h"
// ============================================================================
//   Local
// ============================================================================
#include "status_codes.h"
// ============================================================================
namespace 
{
  // ==========================================================================
  /// convert Python to double
  double result_to_double ( PyObject* r , const char* tag  )
  {
    // ========================================================================
    if  ( !r ) 
    {
      PyErr_Print () ;
      Ostap::Assert ( false                                  ,
                      "CallPython:invalid `result'"          ,
                      tag                                    ,
                      INVALID_PYOBJECT , __FILE__ , __LINE__ ) ;
    }
    // ========================================================================
    // float or integer ?
    // ========================================================================
    if      ( PyFloat_Check ( r ) )  // floating value? 
      {        const double result = PyFloat_AS_DOUBLE ( r ) ;
        Py_DECREF ( r ) ;
        return result ;                                    // RETURN 
      }
    // ========================================================================
    // int in python3 
    // ========================================================================
    else if ( PyLong_Check ( r ) )   // integer value ? 
      {
        int overflow = 0 ;
        const long result = PyLong_AsLongAndOverflow ( r , &overflow );
        if      ( -1 == overflow || 1 == overflow ) 
          {
            Py_DECREF ( r ) ;      
            Ostap::Assert ( false                                  ,
                            "CallPython:long overflow"             ,
                            tag                                    ,
                            INVALID_PYOBJECT , __FILE__ , __LINE__ ) ;
          }
        else if ( -1 == result && PyErr_Occurred() ) 
          {
            PyErr_Print();
            Py_DECREF ( r ) ;      
            Ostap::Assert ( false                                  ,
                            "CallPython:invalid conversion"        ,
                            tag                                    ,
                            INVALID_PYOBJECT , __FILE__ , __LINE__ ) ;
          }
        Py_DECREF ( r ) ;
        return result ;                                    //  RETURN
      } 
    // ========================================================================
    const double result = PyFloat_AsDouble ( r ) ;
    if ( PyErr_Occurred() ) 
      {
        PyErr_Print();
        Py_DECREF ( r ) ;      
        Ostap::Assert ( false                                  ,
                        "CallPython:invalid conversion"        ,
                        tag                                    ,
                        INVALID_PYOBJECT , __FILE__ , __LINE__ ) ;
      }
    //
    Py_DECREF ( r ) ; 
    return result   ;                                     // RETURN
  }
  // ==========================================================================
  /// call python Method  and convert resutl to double 
  double call_method ( PyObject* self , char* method )  
  {
    // check arguments
    Ostap::Assert ( self                                   ,
                    "CallPython:invalid `object'"          , 
                    "CallPython::call_method"              ,
                    INVALID_PYOBJECT , __FILE__ , __LINE__ ) ;
    // call Python
    PyObject* result = PyObject_CallMethod ( self , method , nullptr );    
    // error/exception ?
    return result_to_double ( result , method ) ;
  }
  // ==========================================================================
}
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // CALLPYTHON_H
// ============================================================================
