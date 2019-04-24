// ============================================================================
// Include files
// ============================================================================
// Python
// ============================================================================
#include "Python.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/BLOB.h"
#include "Ostap/PyBLOB.h"
// ============================================================================
/** @file
 *  Implementation file for functions form  file Ostap/PyBLOB.h
 *  @see Ostap::BLOB
 *  @date 2019-04-24 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 */
// ============================================================================
/*  convert blob to python bytes 
 *  @see Ostap::BLOB
 *  @param blob
 *  @return PyBytes object from the blob
 */
// ============================================================================
PyObject* Ostap::blob_to_bytes ( const Ostap::BLOB& blob ) 
{
#if defined (PY_MAJOR_VERSION) and PY_MAJOR_VERSION < 3
  return PyString_FromStringAndSize ( (const char*) blob.buffer () , blob.size () ) ;
#else
  return  PyBytes_FromStringAndSize ( (const char*) blob.buffer () , blob.size () ) ;
#endif  
}
// ============================================================================
/*  convert bytes to blob 
 *  @see   Ostap::BLOB
 *  @param blob the blib to be updated 
 *  @param bytes (INPUT) the bytes 
 *  @return True if conversion successul
 */
// ============================================================================
PyObject* Ostap::blob_from_bytes ( Ostap::BLOB& blob , PyObject* bytes ) 
{
  // 
  // check the arguments 
  //
#if defined (PY_MAJOR_VERSION) and PY_MAJOR_VERSION < 3
  if ( nullptr == bytes || !PyString_Check ( bytes ) )
#else
  if ( nullptr == bytes || !PyBytes_Check  ( bytes ) ) 
#endif 
    {
    PyErr_SetString( PyExc_TypeError, "Invalid bytes/string object" ) ;
    return NULL ;
  } 
  //
  // set the blob 
  //
#if defined (PY_MAJOR_VERSION) and PY_MAJOR_VERSION < 3
  blob.setBuffer ( PyString_Size ( bytes ) , PyString_AsString ( bytes ) ) ;
#else   
  blob.setBuffer ( PyBytes_Size  ( bytes ) , PyBytes_AsString  ( bytes ) ) ;
#endif 
  // 
  Py_INCREF ( Py_True );
  //
  return Py_True ;
}

// ==========================================================================


// ============================================================================
//                                                                      The END 
// ============================================================================
