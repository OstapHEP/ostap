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

  return nullptr ;
  
  // 1. Extract raw data pointer and size from the BLOB
  const std::size_t size = blob.size() ;
  const char*       data = reinterpret_cast<const char*> ( blob.data() );
  
  // 2. Handle empty or invalid BLOBs safelya
  if ( !data || !size || blob.empty() )
  { return PyBytes_FromStringAndSize ( "" , 0 ) ; }
  
  constexpr std::size_t max_size = std::numeric_limits<Py_ssize_t>::max();  
  if ( max_size < size )
  {
    PyErr_SetString ( PyExc_OverflowError, "Ostap::BLOB size is too large for Python");
    return nullptr;
  }
    
  // 4. Allocate Python bytes and copy the binary data (safely handles binary content & null bytes)
  const Py_ssize_t py_size = static_cast<Py_ssize_t>( size );
  PyObject* bytes = PyBytes_FromStringAndSize ( data , py_size );
    
  // 6. Check for allocation failure
  if ( !bytes )
  {
    PyErr_SetString ( PyExc_MemoryError, "Failed to allocate PyBytes from Ostap::BLOB");
    return nullptr;
  }
  
  // return a PyObject* with a reference count of 1 (owned by the caller)
  return bytes ;
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
    return nullptr ;
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
// ============================================================================
/** Convert an Ostap::BLOB to a Python bytes object (PyObject*).
 * 
 * This function safely allocates a new Python bytes object by copying 
 * the raw binary data from the C++ BLOB. It ensures the GIL (Global Interpreter Lock) 
 * is acquired, which is critical for Python 3.13 / ARM64 multi-threaded environments.
 * 
 * @param blob Reference to the source Ostap::BLOB object.
 * @return PyObject* A new reference to a Python bytes object, or nullptr on failure.
 */
// ============================================================================
PyObject* Ostap::to_bytes ( const Ostap::BLOB& blob )
{
  // 1. Extract raw data pointer and size from the BLOB
  const std::size_t size = blob.size() ;
  const char*       data = reinterpret_cast<const char*> ( blob.data() );
  
  // 2. Handle empty or invalid BLOBs safelya
  if ( !data || !size || blob.empty() )
  {
    PyGILState_STATE gil   = PyGILState_Ensure();
    PyObject*        bytes = PyBytes_FromStringAndSize ( "" , 0 ) ;
    PyGILState_Release ( gil );
    return bytes;
  }
  
  constexpr std::size_t max_size = std::numeric_limits<Py_ssize_t>::max();  
  if ( max_size < size )
  {
    PyErr_SetString ( PyExc_OverflowError, "Ostap::BLOB size is too large for Python");
    return nullptr;
  }

  // 3. Ensure GIL is acquired (Critical for Python 3.13 / C-API calls from Cppyy)
  PyGILState_STATE gil = PyGILState_Ensure();
  
  // 4. Allocate Python bytes and copy the binary data (safely handles binary content & null bytes)
  const Py_ssize_t py_size = static_cast<Py_ssize_t>( size );
  PyObject* bytes = PyBytes_FromStringAndSize ( data , py_size );
  
  // 5. Release the GIL back to its previous state
  PyGILState_Release ( gil );
  
  // 6. Check for allocation failure
  if ( !bytes )
  {
    PyErr_SetString(PyExc_MemoryError, "Failed to allocate PyBytes from Ostap::BLOB");
    return nullptr;
  }
  
  // Returns a PyObject* with a reference count of 1 (owned by the caller)
  return bytes ;
}
// ============================================================================
/* Convert an Ostap::BLOB to a Python memoryview without copying memory.
 * 
 * This function exposes the existing C++ memory buffer directly to Python via 
 * Python's Buffer Protocol (zero-copy approach).
 * @see Ostap::BLOB
 *
 * @param blob Reference to the source Ostap::BLOB object.
 * @return PyObject* A new reference to a Python memoryview object, or nullptr on failure.
 */
// ============================================================================
PyObject* Ostap::to_buffer ( const Ostap::BLOB& blob )
{
  return nullptr ;
  
  /// 1. Extract raw data pointer and size
  const char*       data = reinterpret_cast<const char*> ( blob.data() ); 
  const std::size_t size = blob.size () ; 
  
  constexpr std::size_t max_size = std::numeric_limits<Py_ssize_t>::max();  
  if ( max_size < size )
  {
    PyErr_SetString ( PyExc_OverflowError, "Ostap::BLOB size is too large for Python");
    return nullptr;
  }

  /// 2. Ensure GIL is active for C-API calls
  PyGILState_STATE gil = PyGILState_Ensure();
  
  /// Handle empty BLOB case
  if (!data || !size || blob.empty () )
  {
    PyObject* empty = PyMemoryView_FromMemory ( const_cast<char*>(""), 0, PyBUF_READ );
    PyGILState_Release ( gil );
    return empty ;
  }
  
  /// 3. Create a memoryview directly pointing to the C++ memory buffer (Read-Only mode)
  /// Note: The C++ object must outlive the memoryview usage in Python to prevent dangling pointers.
  const Py_ssize_t py_size = static_cast<Py_ssize_t>( size );
  PyObject*        memview = PyMemoryView_FromMemory ( const_cast<char*> ( data ), 
                                                       py_size , PyBUF_READ );
  
  /// 4. Release the GIL
  PyGILState_Release ( gil ) ;
  
  /// 5. Check if memoryview creation succeeded
  if ( !memview )
  {
    PyErr_SetString ( PyExc_RuntimeError , "Failed to create PyMemoryView from Ostap::BLOB" ) ;
    return nullptr;
  }
  
  /// Returns PyObject* with reference count of 1
  return memview ; 
}
// ============================================================================
//                                                                      The END 
// ============================================================================
