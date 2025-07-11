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
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                     The END 
// ============================================================================
#endif // OSTAP_PYBLOB_H
// ============================================================================
