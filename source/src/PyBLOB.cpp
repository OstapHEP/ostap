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
{ return  PyBytes_FromStringAndSize ( (const char*) blob.buffer () , blob.size () ) ; }
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
  if ( nullptr == bytes || !PyBytes_Check  ( bytes ) ) 
    {
      PyErr_SetString( PyExc_TypeError, "Invalid bytes/string object" ) ;
      return NULL ;
    } 
  //
  // set the blob 
  //
  blob.setBuffer ( PyBytes_Size  ( bytes ) , PyBytes_AsString  ( bytes ) ) ;
  // 
  Py_INCREF ( Py_True );
  //
  return Py_True ;
}

// ==========================================================================


// ============================================================================
//                                                                      The END 
// ============================================================================
