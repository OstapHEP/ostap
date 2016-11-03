// Include files 
// ============================================================================
// STD&STL
// ============================================================================
#include <iostream>
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
  const StatusCode&  code    )
  : m_message    ( message   )   
  , m_tag        ( tag       )
  , m_code       ( code      ) 
{}
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
  const StatusCode&       code     ,
  const Ostap::Exception& previous )
  : m_message    ( message   )   
  , m_tag        ( tag       )
  , m_code       ( code      ) 
  , m_previous   ( previous.clone()  )
{}
// ============================================================================
/** Constructor (3)
 *  @param message error message
 *  @param tag "name tag", or exeption type
 *  @param code status code
 *  @param previous "previous" exception (used to improve the error message)
 */
// ============================================================================
Ostap::Exception::Exception 
( const std::string&    message   ,
  const std::string&    tag       ,
  const StatusCode&     code      ,
  const std::exception& previous  )
  : m_message    ( message   )   
  , m_tag        ( tag       )
  , m_code       ( code      ) 
{
  m_message += std::string(":exception(") + previous.what() + ")" ;
}
// ============================================================================
// Copy constructor (deep copying!)
// ============================================================================
Ostap::Exception::Exception ( const Ostap::Exception& right ) 
  : std::exception( right )
  , m_message {     right.message () } 
  , m_tag     {     right.tag     () }
  , m_code    {     right.code    () }
  , m_previous{     right.previous() ? right.previous()->clone() : nullptr }
{}
// ============================================================================
// clone operation
// ============================================================================
Ostap::Exception* Ostap::Exception::clone() const { return new Exception(*this); }
// ============================================================================
// assignment operator
// ============================================================================
Ostap::Exception& Ostap::Exception::operator=( const Ostap::Exception& right )
{
  if (&right == this ) { return *this ; }
  m_message  =   right.message() ;
  m_tag      =   right.tag    () ;
  m_code     =   right.code   () ;
  m_previous.reset( right.previous() ? right.previous()->clone() : nullptr );
  return *this;
}
// ============================================================================
// update the error message to be printed
// ============================================================================
void Ostap::Exception::setMessage ( const std::string& m ) { m_message = m ; }
// ============================================================================
// update name tag
// ============================================================================
void Ostap::Exception::setTag ( const std::string& t ) { m_tag = t ; }
// ============================================================================
//  update the status code for the exception
// ============================================================================
void Ostap::Exception::setCode( const Ostap::StatusCode& s ) { m_code = s ; }
// ============================================================================
// method  for overloaded printout to std::ostream& and MsgStream&
// ============================================================================
std::ostream& Ostap::Exception::fillStream ( std::ostream& os ) const 
{
  os << tag() << " \t " << message() ;
  switch( code() ) {
  case Ostap::StatusCode::SUCCESS : os << "\t StatusCode=SUCCESS"    ;  break ;
  case Ostap::StatusCode::FAILURE : os << "\t StatusCode=FAILURE"    ;  break ;
  default                         : os << "\t StatusCode=" << code() ;  break ;
  }
  return ( 0 != previous() ) ? previous()->fillStream( os << std::endl ) : os ;
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
/* throw the exception
 *  @param message the reason 
 *  @param tag     the tag 
 *  @param code    the code 
 */
// ============================================================================
void Ostap::throwException
( const std::string&       message ,
  const std::string&       tag     , 
  const Ostap::StatusCode& code    ) 
{ throw Ostap::Exception ( message , tag , code ) ; }
// ===========================================================================  


// ============================================================================
//  The END 
// ============================================================================
