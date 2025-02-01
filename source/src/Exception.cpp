// Include files 
// ============================================================================
// STD&STL
// ============================================================================
#include <sstream>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatusCode.h"
// ============================================================================
// Local 
// ============================================================================
#include "Exception.h"
// ============================================================================
/** Implementation file for class Ostap::Exception
 *  @see Ostap::Exception
 */ 
// ============================================================================
/*  Constructor (1)
 *  @param Message error message
 *  @param Tag "name tag", or exeption type
 *  @param Code status code
 */
// ============================================================================
Ostap::Exception::Exception 
( const std::string& message ,
  const std::string& tag     ,
  const StatusCode&  code    ,
  const char*        file    ,
  const std::size_t  line    ) 
  : m_message    ( message   )   
  , m_tag        ( tag       )
  , m_code       ( code      )
  , m_file       ( file ? file : "" )
  , m_line       ( line      )
  , m_what       ()
  , m_previous   () 
{
  m_what = toString() ;
}
// ============================================================================
/*  constructor from std:exception 
 *  @param exc exception 
 *  @param  file file name 
 *  @param line line numbr 
 */
// ============================================================================
Ostap::Exception::Exception 
( const std::exception& exc   ,
  const char*           file  ,
  const std::size_t     line  )
  : m_message    ( exc.what() )
  , m_tag        ( typeid ( exc ).name()      )
  , m_code       ( Ostap::StatusCode::FAILURE )
  , m_file       ( file ? file : "" )
  , m_line       ( line      )
  , m_what       ()
  , m_previous   () 
{
  m_what = toString() ;  
}
// ============================================================================
/*  Constructor (2)
 *  @param Message error message
 *  @param Tag "name tag", or exeption type
 *  @param Code status code
 *  @param previous "previous"  exception
 */
// ============================================================================
Ostap::Exception::Exception 
( const std::string&      message  ,
  const std::string&      tag      ,
  const Ostap::Exception& previous , 
  const StatusCode&       code     ,
  const char*             file     ,
  const std::size_t       line     )
  : m_message    ( message   )   
  , m_tag        ( tag       )
  , m_code       ( code      )
  , m_file       ( file ? file : "" )
  , m_line       ( line      )
  , m_what       ()
  , m_previous   ( previous.clone() )
{
  m_what = toString() ;
}
// ============================================================================
/** Constructor (3)
 *  @param message error message
 *  @param tag "name tag", or exeption type
 *  @param code status code
 *  @param previous "previous" exception (used to improve the error message)
 */
// ============================================================================
Ostap::Exception::Exception 
( const std::string&       message   ,
  const std::string&       tag       ,
  const std::exception&    previous  , 
  const Ostap::StatusCode& code      ,
  const char*              file      ,
  const std::size_t        line      ) 
  : Exception ( message , tag , Exception ( previous ) , code , file , line ) 
{}
// ============================================================================
// Copy constructor (deep copying!)
// ============================================================================
Ostap::Exception::Exception ( const Ostap::Exception& right ) 
  : std::exception( right )
  , m_message { right.m_message  }  
  , m_tag     { right.m_tag      }
  , m_code    { right.m_code     }
  , m_file    ( right.m_file     )
  , m_line    ( right.m_line     )
  , m_what    ( right.m_what     )
  , m_previous{ right.previous () ? right.previous()->clone() : nullptr }
{}
// ============================================================================
// clone operation
// ============================================================================
Ostap::Exception* Ostap::Exception::clone() const { return new Exception(*this); }
// ============================================================================
// ============================================================================
// method  for overloaded printout to std::ostream& and MsgStream&
// ============================================================================
std::ostream&
Ostap::Exception::fillStream
( std::ostream& os ) const 
{
  //
  const static std::string s_exception = " EXCEPTION : " ;
  const static std::string s_index     = " --- INDEX : " ;
  const static std::string s_tag       = " ---   TAG : " ;
  const static std::string s_code      = " ---  CDDE : " ;
  const static std::string s_file      = " ---  FILE : " ;
  const static std::string s_line      = " ---  line : " ;
  //
  os << s_exception << m_message ; 
  //
  const std::string newline = "\n" ; 
  const int ind = index() ;
  if ( ind    ) { os << newline << s_index  << ind ; }  
  if ( !m_tag.empty() ) { os << newline << s_tag << m_tag   ; }
  // 
  switch ( m_code )
    {
    case Ostap::StatusCode::SUCCESS     : os << newline << s_code << "SUCCESS"     ;  break ;
    case Ostap::StatusCode::FAILURE     : os << newline << s_code << "FAILURE"     ;  break ;
    case Ostap::StatusCode::RECOVERABLE : os << newline << s_code << "RECOVERABLE" ;  break ;
    default                             : os << newline << s_code << m_code        ;  break ;
    }
  //
  if ( !m_file.empty() ) { os << newline << s_file << m_file ; } 
  if (  m_line         ) { os << newline << s_line << m_line ; }
  //
  if ( m_previous ) { os << newline ; m_previous->fillStream ( os ) ; }
  //
  return os ;
}
// ============================================================================
// conversion to string 
// ============================================================================
std::string Ostap::Exception::toString  () const 
{
  std::ostringstream os ;
  fillStream ( os ) ;
  return os.str() ;
}
// ============================================================================
/*  throw the exception
 *  @param message the reason 
 *  @param tag     the tag 
 *  @param code    the code 
 */
// ============================================================================
Ostap::StatusCode Ostap::throwException
( const std::string&       message ,
  const std::string&       tag     , 
  const Ostap::StatusCode& code    ,
  const char*              file    ,
  const std::size_t        line    )
{
  throw Ostap::Exception ( message , tag , code , file , line ) ;
  return code ;
}
// ===========================================================================  


// ============================================================================
//  The END 
// ============================================================================
