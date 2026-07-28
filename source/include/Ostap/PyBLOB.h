// ============================================================================
#ifndef OSTAP_PYBLOB_H 
#define OSTAP_PYBLOB_H 1
// ============================================================================
// Include files
// ============================================================================
// Python
// ============================================================================
struct  _object ;
typedef _object PyObject ;
// ============================================================================
// Ostap
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  class BLOB ;
  // ==========================================================================
  /** convert blob to python bytes 
   *  @see Ostap::BLOB
   *  @param blob
   *  @return PyBytes object from the blob
   */
  PyObject* blob_to_bytes   ( const BLOB& blob ) ;
  // ==========================================================================
  /** convert bytes to blob 
   *  @see   Ostap::BLOB
   *  @param blob the blib to be updated 
   *  @param bytes (INPUT) the bytes 
   *  @return True if conversion successul
   */
  PyObject* blob_from_bytes ( BLOB& blob , PyObject* bytes ) ;
  // ==========================================================================
  /** Convert an Ostap::BLOB to a Python memoryview without copying memory.
   * 
   * This function exposes the existing C++ memory buffer directly to Python via 
   * Python's Buffer Protocol (zero-copy approach).
   * @see Ostap::BLOB
   *
   * @param blob Reference to the source Ostap::BLOB object.
   * @return PyObject* A new reference to a Python memoryview object, or nullptr on failure.
   */
  PyObject* to_buffer ( const Ostap::BLOB& blob )  ;
  // ==============================================================================
  /** Convert an Ostap::BLOB to a Python bytes object (PyObject*).
   * 
   * This function safely allocates a new Python bytes object by copying 
   * the raw binary data from the C++ BLOB. It ensures the GIL (Global Interpreter Lock) 
   * is acquired, which is critical for Python 3.13 / ARM64 multi-threaded environments.
   * 
   * @param blob Reference to the source Ostap::BLOB object.
   * @return PyObject* A new reference to a Python bytes object, or nullptr on failure.
   */
  PyObject* to_bytes ( const Ostap::BLOB& blob ) ;      
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                     The END 
// ============================================================================
#endif // OSTAP_PYBLOB_H
// ============================================================================
