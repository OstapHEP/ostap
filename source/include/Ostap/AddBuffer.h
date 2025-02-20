// ============================================================================
#ifndef OSTAP_ADDBUFFER_H 
#define OSTAP_ADDBUFFER_H 1
// ============================================================================
// Include files 
// ============================================================================
// STD&STL
// ============================================================================ 
#include <string>
// ============================================================================
#if defined ( __cplusplus ) && ( 202002L <= __cplusplus ) // ==================
// ============================================================================
#include <span>  // ===========================================================
// ============================================================================
#endif // =====================================================================
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/IFuncs.h"
#include "Ostap/StatusCode.h"
#include "Ostap/ProgressBar.h"
// ============================================================================
// Forward declarations from ROOT 
// ============================================================================
class TTree ; 
// ============================================================================
namespace Ostap
{  
  // ==========================================================================
  namespace Trees
  { 
    // ========================================================================
    /** valid name for branch ?
     *  - not empty 
     *  - no blanks 
     *  - no special symbols 
     */
    bool valid_name_for_branch ( const std::string& name ) ;
    // ========================================================================
#if defined ( __cplusplus ) && ( 202002L <= __cplusplus ) && ( 1 > 2 ) // ==================
    // ========================================================================
    /** @class Buffer 
     *  Helper cladd to add the content of buffer to TTrer
     *  - actually it is span + defautl value 
     *  @date 2025-02-05
     */
    template <class DATA> 
    class Buffer
    {
      // ======================================================================
      typedef std::span<const DATA>  SPAN ;
      // ======================================================================
    public: // ================================================================
      // ======================================================================
      Buffer
      ( const DATA*       data  = nullptr    ,
        const std::size_t size  = 9          ,
        const DATA        value = DATA ( 0 ) )
        : m_span  ( data  , size ) 
        , m_value ( value ) 
      {} ;                   
      // ======================================================================
    public: // ================================================================
      // ======================================================================
      /// create new fuffer with offset 
      inline Buffer offset ( const std::size_t offset ) const ; 
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
     *  Helper class to add the content of buffer to TTrer
     *  - actually it is span + defautl value 
     *  @date 2025-02-05
     */
    template <typename DATA> 
    class Buffer
    {
      // ======================================================================
    public: // ================================================================
      // ======================================================================
      Buffer
      ( const DATA*       data  = nullptr    ,
        const std::size_t size  = 9          ,
        const DATA        value = DATA ( 0 ) )
        : m_data  ( data  )
        , m_size  ( size  )
        , m_value ( value ) 
      {} ;                   
      // ======================================================================
    public: // ================================================================
      // ======================================================================
      /// create new fuffer with offset 
      inline Buffer offset ( const std::size_t offset ) const ;
      // ======================================================================
    public: // ================================================================
      // ======================================================================
      bool        empty () const { return m_size == 0     ; }
      std::size_t size  () const { return m_size          ; }
      const DATA& value () const { return m_value         ; }
      const DATA* data  () const { return m_data          ; }
      const DATA* begin () const { return m_data          ; }
      const DATA* end   () const { return m_data + m_size ; }      
      const DATA& operator[] ( const std::size_t index ) const
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
      const DATA* m_data  { nullptr } ;
      std::size_t m_size  { 0 } ;
      DATA        m_value { 0 } ;
      // ======================================================================
    } ;  // ===================================================================
    // ========================================================================
#endif // =====================================================================
    // ========================================================================
    /// swap two buffers
    template<class DATA>
    inline void swap ( Buffer<DATA>& a , Buffer<DATA>& b ) { a.swap ( b ) ; }
    // ========================================================================
    /// create new fuffer with offset 
    template <class DATA> 
    inline
    Buffer<DATA>
    Buffer<DATA>::offset ( const std::size_t offset ) const 
    { return ( size() <= offset ) ? 
        Buffer ( data () + size   () , 0                , value () ) :
        Buffer ( data () + offset    , size () - offset , value () ) ;
    }
    // =========================================================================
    template <class DATA>
    inline
    Buffer<DATA>
    make_buffer
    ( const DATA*       data           ,
      const std::size_t size           ,
      const DATA        value = DATA() )
    { return Buffer<DATA> ( data , size , value ) ; }
    // =========================================================================
    inline
    Buffer<char>
    char_buffer
    ( const void*       data              ,
      const std::size_t size              ,
      const char        value = char( 0 ) )
    { return Buffer<char> ( static_cast<const  char*> ( data ) , size , value ) ; }
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
      typedef Ostap::Trees::Buffer<DATA>                              BUFFER  ;
      typedef std::map<std::string,BUFFER>                            BUFFERS ;
      typedef typename BUFFERS::const_iterator                 const_iterator ;
      // ======================================================================
    public :
      // ======================================================================
      /// add new buffer into the map 
      void add ( const std::string& name , const BUFFER& buffer )
      { m_buffers.insert ( BUFFERS::value_type ( name , buffer ) ) ; }
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
      /// actual map of  bufferss 
      BUFFERS m_buffers {} ; //!        actual map of  bufferss       
      // ======================================================================
    } ; // ====================================================================
    // ========================================================================
    /// add buffer to tree 
    Ostap::StatusCode
    add_buffer
    ( TTree*                            tree     ,
      const std::string&                name     ,
      const Buffer<Double_t>&           buffer   , 
      const Ostap::Utils::ProgressConf& progress = false ) ;                 
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree* tree                                ,
      const std::string&                name     ,
      const Buffer<Float_t>&            buffer   ,
      const Ostap::Utils::ProgressConf& progress = false ) ; 
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree*                            tree     ,
      const std::string&                name     ,
      const Buffer<Char_t>&             buffer   , 
      const Ostap::Utils::ProgressConf& progress = false ) ; 
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree*                            tree     ,
      const std::string&                name     ,
      const Buffer<UChar_t>&            buffer   , 
      const Ostap::Utils::ProgressConf& progress = false ) ; 
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree*                            tree     ,
      const std::string&                name     ,
      const Buffer<Short_t>&            buffer   , 
      const Ostap::Utils::ProgressConf& progress = false ) ;
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree*                            tree     ,
      const std::string&                name     ,
      const Buffer<UShort_t>&           buffer   , 
      const Ostap::Utils::ProgressConf& progress = false ) ; 
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree*                            tree     ,
      const std::string&                name     ,
      const Buffer<Int_t>        &      buffer   , 
      const Ostap::Utils::ProgressConf& progress = false ) ; 
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree*                            tree     ,
      const std::string&                name     ,
      const Buffer<UInt_t>&             buffer   , 
      const Ostap::Utils::ProgressConf& progress = false ) ; 
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree*                            tree     ,
      const std::string&                name     ,
      const Buffer<Long64_t>&           buffer   , 
      const Ostap::Utils::ProgressConf& progress = false ) ; 
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree*                            tree     ,
      const std::string&                name     ,
      const Buffer<ULong64_t>&          buffer   , 
      const Ostap::Utils::ProgressConf& progress = false ) ;
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree*                            tree     ,
      const std::string&                name     ,
      const Buffer<Long_t>&             buffer   , 
      const Ostap::Utils::ProgressConf& progress = false ) ; 
    // ========================================================================
    Ostap::StatusCode
    add_buffer
    ( TTree*                            tree     ,
      const std::string&                name     ,
      const Buffer<ULong_t>&            buffer   , 
      const Ostap::Utils::ProgressConf& progress = false ) ;
    // ========================================================================
    /** add <code>long dobule</code> buffer to TTree
     *  @attentoon actually doubel valeus are stored 
     */ 
    Ostap::StatusCode
    add_buffer
    ( TTree*                            tree     ,
      const std::string&                name     ,
      const Buffer<long double>&        buffer   , 
      const Ostap::Utils::ProgressConf& progress = false ) ;                 
    // ========================================================================
  } //                                        The end of namespace Ostap::Trees 
  // ==========================================================================
} //                                                 The end of namesapce Ostap 
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_ADDBUFFER_H
// ============================================================================
