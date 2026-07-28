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
#include <iostream>
PyObject* Ostap::blob_to_bytes ( const Ostap::BLOB& blob ) 
{
  std::cerr << " blob_to_bytes : "
            <<  blob.GetName()
            << "size="
            <<  blob.size() << std::endl ;
  if ( blob.empty () || !blob.buffer () ) { return PyBytes_FromStringAndSize ( "" , 0 ) ; }
  return PyBytes_FromStringAndSize ( reinterpret_cast<const char*> ( blob.buffer () ) , blob.size () );
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
  if ( nullptr == bytes || !PyBytes_Check  ( bytes ) ) 
  {
    PyErr_SetString ( PyExc_TypeError, "Invalid bytes/string object" ) ;
    return NULL ;
  } 
  //
  std::cerr << " blob_from_bytes/1 : "
            <<  blob.GetName()
            << " bytes=" << PyBytes_Size ( bytes ) << std::endl ;
  
  // set the blob 
  //
  blob.setBuffer ( PyBytes_Size  ( bytes ) , PyBytes_AsString  ( bytes ) ) ;
  // 
  Py_INCREF ( Py_True );
  //
  std::cerr << " blob_from_bytes/2 : "
            <<  blob.GetName()
            <<  blob.size () << std::endl ;
  //
  return Py_True ;
}
// ============================================================================
//                                                                      The END 
// ============================================================================
