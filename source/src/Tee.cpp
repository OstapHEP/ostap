// ============================================================================
// Include 
// ============================================================================
#include <memory>
#include <fstream>
#include <streambuf>
#include <iostream>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Tee.h"
// ============================================================================
/** @file  Tee.cpp
 *  Implementation file for class Ostap::Utils::Tee
 *  @see Gaudi::Utils::Tee
 *  @see Gaudi::Utils::Mute
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2013-07-07
 */
// ============================================================================
namespace std 
{
  // ==========================================================================
  /** @class basic_teebuf 
   *  implementation of "tee"-functionality 
   *  the actual code is copied from 
   *  http://wordaligned.org/articles/cpp-streambufs
   */  
  template <typename CHAR_TYPE, typename TRAITS = std::char_traits<CHAR_TYPE> >
  class basic_teebuf: public std::basic_streambuf<CHAR_TYPE,TRAITS>
  {
  public:
    // ========================================================================
    typedef typename TRAITS::int_type int_type;
    // ========================================================================
  public:
    // ========================================================================
    basic_teebuf 
    ( std::basic_streambuf<CHAR_TYPE,TRAITS> * sb1 ,
      std::basic_streambuf<CHAR_TYPE,TRAITS> * sb2 )
      : sb1(sb1)
      , sb2(sb2)
    {}
    // ========================================================================    
  private:    
    // ========================================================================
    virtual int sync()
    {
      int const r1 = sb1->pubsync();
      int const r2 = sb2->pubsync();
      return r1 == 0 && r2 == 0 ? 0 : -1;
    }
    //
    virtual int_type overflow ( int_type c )
    {
      //
      int_type const eof = TRAITS::eof();
      //
      if ( TRAITS::eq_int_type(c, eof)) { return TRAITS::not_eof(c); }
      else
      {
        //
        CHAR_TYPE const ch = TRAITS::to_char_type(c);
        int_type  const r1 = sb1->sputc(ch);
        int_type  const r2 = sb2->sputc(ch);
        return
          TRAITS::eq_int_type(r1, eof) ||
          TRAITS::eq_int_type(r2, eof) ? eof : c;
        //
      }
    }
    // ========================================================================
  private:
    // ========================================================================
    std::basic_streambuf<CHAR_TYPE,TRAITS> *sb1 ;
    std::basic_streambuf<CHAR_TYPE,TRAITS> *sb2 ;
    // ========================================================================
  };
  //
  // ==========================================================================  
  /** @class basic_teestream 
   *  implementation of "tee"-functionality 
   *  the actual code is copied from 
   *  http://wordaligned.org/articles/cpp-streambufs
   */  
  template < typename CHAR_TYPE, typename TRAITS = std::char_traits<CHAR_TYPE> >
  struct basic_teestream : public std::basic_ostream<CHAR_TYPE,TRAITS>
  {
    // ========================================================================
    typedef std::basic_ostream<CHAR_TYPE,TRAITS> stream_type ;
    typedef std::basic_teebuf <CHAR_TYPE,TRAITS> streambuff_type ;
    // ========================================================================    
    basic_teestream ( stream_type& first, stream_type& second )
      : stream_type(&stmbuf), stmbuf( first.rdbuf(), second.rdbuf() ) {}
    //
    basic_teestream ( streambuff_type* first, streambuff_type* second )
      : stream_type(&stmbuf), stmbuf(first,second ) {}
    //
    ~basic_teestream() { stmbuf.pubsync() ; }    
    // ========================================================================
  private: 
    // ========================================================================    
    streambuff_type stmbuf ;
    // ========================================================================
  };
  // ==========================================================================
  typedef basic_teebuf<char>    teebuf    ;
  typedef basic_teestream<char> teestream ;
  // ==========================================================================
}
// ============================================================================

// ============================================================================
// constructor from the filename 
// ============================================================================
Ostap::Utils::Tee::Tee  ( const std::string&  filename   )
  : m_file        ( new std::ofstream ( filename.c_str() ) ) 
  , m_own         ( true    ) 
  , m_buffer_cout ( nullptr ) 
  , m_buffer_cerr ( nullptr ) 
  , m_keep_cout   ( nullptr ) 
  , m_keep_cerr   ( nullptr ) 
{
  std::cout << std::flush ;
  std::cerr << std::flush ;
  m_keep_cout = std::cout.rdbuf() ;
  m_keep_cerr = std::cerr.rdbuf() ;
  m_buffer_cout.reset ( new std::teebuf ( m_keep_cout , m_file->rdbuf() ) ) ;
  m_buffer_cerr.reset ( new std::teebuf ( m_keep_cerr , m_file->rdbuf() ) ) ;
  std::cout.rdbuf ( m_buffer_cout.get() ) ;
  std::cerr.rdbuf ( m_buffer_cerr.get() ) ;
}
// ============================================================================
// constructor from the stream
// ============================================================================
// Ostap::Utils::Tee::Tee  ( std::ostream&  filestream   )
//   : m_file   ( &filestream ) 
//   , m_own    ( false       ) 
//   , m_buffer (   ) 
//   , m_keep   ( nullptr     ) 
// {
//   std::cout << std::flush ;
//   m_keep   = std::cout.rdbuf() ;
//   m_buffer.reset ( new std::teebuf ( m_keep , m_file->rdbuf() ) ) ;
//   std::cout.rdbuf ( m_buffer.get() ) ;
// }
// ============================================================================
// destructor 
// ============================================================================
Ostap::Utils::Tee::~Tee() { exit () ; }
// ============================================================================
/*  helper function to implement python's __enter__  
 *  the action is performed in constructor 
 */
// ============================================================================
void Ostap::Utils::Tee::enter (){}
// ============================================================================
//  helper function to implement python's __exit__
// ============================================================================
void Ostap::Utils::Tee::exit  ()  
{
  std::cout << std::flush ;
  std::cerr << std::flush ;
  if ( m_file ) { (*m_file) << std::flush ; }
  //
  // 1. restore std::cout 
  if ( nullptr != m_keep_cout ) 
  {
    std::cout.rdbuf ( m_keep_cout ) ;
    m_keep_cout = nullptr ;
  }
  // 2. restore std::cerr
  if ( nullptr != m_keep_cerr ) 
  {
    std::cerr.rdbuf ( m_keep_cerr ) ;
    m_keep_cerr = nullptr ;
  }
  // 3. delete/close  the file (if needed) 
  if ( m_file ) 
  {
    if ( m_own ) { m_file.reset   () ; }
    else         { m_file.release () ; }
  }
}
// ============================================================================
//                                                                      The END 
// ============================================================================


