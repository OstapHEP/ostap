// ============================================================================
#ifndef OSTAP_BUFFER_H 
#define OSTAP_BUFFER_H 1
// ============================================================================
// Incldue files
// ============================================================================
// STD STL 
// ============================================================================
#include <cstddef>
#include <iterator>
#include <algorithm>
#include <map>
#include <type_traits>
// ============================================================================
// ROOT 
// ============================================================================
#include "RtypesCore.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Types.h"
#include "Ostap/Span.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Utils
  {
    // ========================================================================
#if defined(OSTAP_HAS_STD_SPAN) && OSTAP_HAS_STD_SPAN
    // ========================================================================
    /** @class Buffer 
     *  Helper class to add the content of buffer to TTree
     *  - actually it is span + default value 
     *  @date 2025-02-05
     */
    template <class DATA>
    class Buffer
    {
      // ======================================================================
      static_assert ( Ostap::is_numeric_or_byte<DATA>::value, "Ostap::Utils::Buffer can only be used with numeric types or std::byte!");
      // ======================================================================
    public :
      // ======================================================================
      /// data type 
      using TYPE = std::remove_cv_t<std::remove_reference_t<DATA>>;
      /// span type
      using SPAN = std::span<const TYPE> ;
      // ======================================================================
    public : // ===============================================================
      // ======================================================================
      Buffer
      ( const DATA*       data  = nullptr ,
        const std::size_t size  = 0       ,
        const DATA        value = DATA {} )
        : m_span  ( data  , size ) 
        , m_value ( value ) 
      {} ;                   
      // ======================================================================
    public: // ================================================================
      // ======================================================================
      /// create new buffer with offset 
      inline Buffer offset ( const std::size_t offset ) const
      { return   ( this -> size () <= offset ) ? 
          Buffer ( this -> data () + this->size () , 0                        , this -> value () ) :
          Buffer ( this -> data () + offset        , this -> size () - offset , this -> value () ) ; }      
      // ======================================================================
    public: // ================================================================
      // ======================================================================
      bool        empty () const { return m_span.empty () ; }
      std::size_t size  () const { return m_span.size  () ; }
      const DATA& value () const { return m_value         ; }
      typename SPAN::const_pointer   data  () const { return m_span.data  () ; }
      typename SPAN::iterator        begin () const { return m_span.begin () ; }
      typename SPAN::iterator        end   () const { return m_span.end   () ; }
      typename SPAN::const_reference operator[] ( const std::size_t index ) const
      { return  index < size () ? m_span [ index ] : m_value ; }
      // ======================================================================
    public: // ================================================================
      // ======================================================================
      /// set new default value 
      void setValue ( const DATA new_value ) { m_value = new_value ; }
      // ======================================================================
    public: // ===============================================================
      // ======================================================================
      /// swap two buffers 
      void swap ( Buffer& another )
      {
        std::swap ( m_span  , another.m_span  ) ;
        std::swap ( m_value , another.m_value ) ;
      }
      // ======================================================================
    public : // ===============================================================
      // ======================================================================
      /// conversion to span
      operator const SPAN& () const { return m_span ; } 
      // ======================================================================
    private: // ===============================================================
      // ======================================================================
      SPAN m_span  {   } ;
      DATA m_value { 0 } ;
      // ======================================================================
    } ;  // ===================================================================
    // ========================================================================
#else // ======================================================================
    // ========================================================================
    /** @class Buffer 
     *  Helper buffer class 
     *  - actually it is span + default value 
     *  @date 2025-02-05
     */
    template <typename DATA>
    class Buffer
    {
      // ======================================================================
      static_assert ( Ostap::is_numeric_or_byte<DATA>::value, "Ostap::Utils::Buffer can only be used with numeric types or std::byte!");
      // ======================================================================
    public:
      // ======================================================================
      /// data type 
      using TYPE = std::remove_cv_t<std::remove_reference_t<DATA>>;
      // ======================================================================      
    public: // ================================================================
      // ======================================================================
      Buffer
      ( const TYPE*       data  = nullptr ,
        const std::size_t size  = 0       ,
        const TYPE        value = TYPE {} )
        : m_data  ( data  )
        , m_size  ( size  )
        , m_value ( value ) 
      {} ;                   
      // ======================================================================
    public: // ================================================================
      // ======================================================================
      /// create new fuffer with offset 
      inline Buffer offset ( const std::size_t offset ) const
      { return   ( this -> size () <= offset ) ? 
          Buffer ( this -> data () + this->size () , 0                        , this -> value () ) :
          Buffer ( this -> data () + offset        , this -> size () - offset , this -> value () ) ; }            
      // ======================================================================
    public: // ================================================================
      // ======================================================================
      bool        empty () const { return m_size == 0     ; }
      std::size_t size  () const { return m_size          ; }
      const TYPE& value () const { return m_value         ; }
      const TYPE* data  () const { return m_data          ; }
      const TYPE* begin () const { return m_data          ; }
      const TYPE* end   () const { return m_data + m_size ; }      
      const TYPE& operator[] ( const std::size_t index ) const
      { return  index < m_size ? *(m_data+index) : m_value ; }
      // ======================================================================
    public: // ================================================================
      // ======================================================================
      /// set new default value 
      void setValue ( const DATA new_value ) { m_value = new_value ; }
      // ======================================================================
    public: // ===============================================================
      // ======================================================================
      /// swap two buffers 
      void swap ( Buffer& another )
      {
        std::swap ( m_data  , another.m_data  ) ;
        std::swap ( m_size  , another.m_size  ) ;
        std::swap ( m_value , another.m_value ) ;
      }
      // ======================================================================
    private: // ===============================================================
      // ======================================================================
      const TYPE* m_data  { nullptr } ;
      std::size_t m_size  { 0 } ;
      TYPE        m_value {} ;
      // ======================================================================
    } ;  // ===================================================================
    // ========================================================================
#endif // OSTAP_HAS_STD_SPAN ==================================================
    // ========================================================================
    /// swap two buffers
    template<class DATA>
    inline void swap
    ( Buffer<DATA>& a ,
      Buffer<DATA>& b ) { a.swap ( b ) ; }
    // ========================================================================
    template <class SCALAR>    
    inline
    Buffer<SCALAR>
    make_buffer
    ( const SCALAR*       data              ,
      const std::size_t size                ,
      const SCALAR        value = SCALAR {} )
    { return Buffer<SCALAR> ( data , size , value ) ; }
    // =========================================================================
    inline 
    Buffer<signed char>
    schar_buffer
    ( const void*       data       ,
      const std::size_t size       ,
      const signed char value = 0  )
    { return Buffer<signed char> ( static_cast<const signed char*> ( data ) , size , value ) ; }
    // =========================================================================
    inline
    Buffer<unsigned char>
    uchar_buffer
    ( const void*         data      ,
      const std::size_t   size      ,
      const unsigned char value = 0 )
    { return Buffer<unsigned char> ( static_cast<const unsigned char*> ( data ) , size , value ) ; }
    // =======================================================================
    inline
    Buffer<char>
    char_buffer
    ( const void*         data       ,
      const std::size_t   size       ,
      const char          value = 0  )
    { return Buffer<char> ( static_cast<const char*> ( data ) , size , value ) ; }
    // =======================================================================
#if defined(OSTAP_HAS_STD_BYTE) && OSTAP_HAS_STD_BYTE
    // ========================================================================
    inline
    Buffer<std::byte>
    byte_buffer
    ( const void*         data  ,
      const std::size_t   size  ,
      const std::byte     value = std::byte{} )
    { return Buffer<std::byte> ( static_cast<const std::byte*> ( data ) , size , value ) ; }
    // =========================================================================
#endif // OSTAP_HAS_STD_BYTE
    // =========================================================================    
    /** @class Buffers 
     *  a collection of several named buffers 
     */
    template <class DATA>
    class Buffers
    {
      // ======================================================================
    public:
      // ======================================================================
      typedef Ostap::Utils::Buffer<DATA>                              BUFFER  ;
      typedef Ostap::Dict<BUFFER>                                     BUFFERS ;
      typedef typename BUFFERS::const_iterator                 const_iterator ;
      // ======================================================================
    public :
      // ======================================================================
      /// add new buffer into the map 
      inline void add
      ( const std::string& name ,
        const BUFFER& buffer )
      {
        typedef typename BUFFERS::value_type value_type ;
        m_buffers.insert ( value_type ( name , buffer ) ) ;
      }
      // ======================================================================
    public :
      // ======================================================================
      bool           empty () const { return m_buffers.empty () ; }
      std::size_t    size  () const { return m_buffers.size  () ; }
      const_iterator begin () const { return m_buffers.begin () ; }
      const_iterator end   () const { return m_buffers.end   () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// create new buffers with offset 
      Buffers offset ( const std::size_t offset ) const
      {
        Buffers result {}  ;
        for ( const_iterator it = this->begin() ; this->end() != it ; ++it )
          { result.add ( it->first , it->second.offset ( offset ) ) ; } 
        return result;       
      }
      // ======================================================================
    private :
      // ====================================================================== 
      /// actual map of  buffers
      BUFFERS m_buffers {} ; //!        actual map of  buffers
      // ======================================================================
    } ; // ====================================================================
    // ========================================================================
  } //                                        The END of namespace Ostap::Utils  
  // ==========================================================================
} //                                                  The END of namespaceOstap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_BUFFER_H
// ============================================================================
