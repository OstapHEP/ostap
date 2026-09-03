// ============================================================================
#ifndef OSTAP_BLOB_H 
#define OSTAP_BLOB_H 1
// ============================================================================
// Include files
// ============================================================================
// ROOT 
// ============================================================================
#include "TNamed.h"
#include "TArrayC.h"
// ============================================================================
// Ostap:
// ============================================================================
#include "Ostap/Types.h"
#include "Ostap/Span.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  /** @class  BLOB Blob.h Ostap/Blob.h
   *  Trivial ROOT-based class to store blobs in ROOT file
   *  @author Vanya Belyaev
   *  @date   2019-03-27
   */
  class BLOB : public TNamed 
  {
  public:
    // ========================================================================
    ClassDefOverride(Ostap::BLOB,1) ;
    // ========================================================================
  public: 
    // ========================================================================
    /// Standard constructor
    BLOB
    ( const std::string& name   = "" , 
      const std::string& title  = "" ,
      const std::size_t  len    = 0  , 
      const void*        buffer = 0  ) ;
    /// destructor 
    virtual ~BLOB() ;
    // ========================================================================
  public: // gettters 
    // ========================================================================
    /// get the size of the buffer 
    inline std::size_t size   () const { return 0 >= m_data.GetSize  () ? 0 : m_data.GetSize () ; }
    /// get the buffer/data itself 
    inline const void* buffer () const { return      m_data.GetArray () ; }    
    /// get the buffer/data itself 
    inline const void* data   () const { return      m_data.GetArray () ; }    
    /// get the buffer/data itself 
    inline bool        empty  () const { return !m_data.GetArray() || 0 >= m_data.GetSize  () ; }    
    // ========================================================================
  public: // setters 
    // ========================================================================
    /// redefine the buffer 
    void setBuffer 
    ( const std::size_t size   , 
      const void*       buffer ) ;
    // ========================================================================
  public : // Bytes & span 
    // ========================================================================
#if defined(OSTAP_HAS_STD_BYTE) && OSTAP_HAS_STD_BYTE
    // ========================================================================
    /// get the buffer as const std::byte*
    inline const std::byte*     bytes     () const 
    { return reinterpret_cast<const std::byte*> ( m_data.GetArray() ) ; }
    /// set the buffer from bytes 
    inline void                 setBuffer ( const std::size_t size   ,
                                            const std::byte*  buffer )
    { this->setBuffer ( size , static_cast<const void*> ( buffer ) ) ; }    
    // ========================================================================
#if defined(OSTAP_HAS_STD_SPAN) && OSTAP_HAS_STD_SPAN
    // ========================================================================
    /// get std::span<const std::byte> view of the BLOB
    inline std::span<const std::byte> byte_span () const 
    { return std::span<const std::byte> ( bytes () , size() ) ; }
    /// set the buffer from span 
    inline void                 setBuffer ( const std::span<const std::byte>& buffer )
    { this->setBuffer ( buffer.size () , buffer.data () ) ; }     
    // ========================================================================   
#endif // OSTAP_HAS_STD_SPAN 
#endif // OSTAP_HAS_STD_BYTE 
    // ========================================================================
  private:
    // ========================================================================
    /// the data itself  
    TArrayC m_data {} ; // the data buffer 
    // ========================================================================
  };
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                     The END 
// ============================================================================
#endif // OSTAP_BLOB_H
// ============================================================================
