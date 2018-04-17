#ifndef OSTAP_EXCEPTION_H
#define OSTAP_EXCEPTION_H
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <string>
#include <exception>
#include <memory>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatusCode.h"
// ============================================================================
namespace Ostap 
{
  // ==========================================================================
  /** @class Exception Exception.h 
   *  Define Ostap exception
   *  Actually it is a bit simplified vertion of GaudiException 
   *  class from Gaudi Project
   *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
   */
  class Exception: public std::exception 
  {
  public:
    // ========================================================================
    /** Constructor (1)
     *  @param message error message
     *  @param tag     "name tag", or exception type
     *  @param code     status code
     */
    Exception ( const std::string& message ,
                const std::string& tag     ,
                const StatusCode&  code    ) ;
    // ========================================================================
    /** Constructor (2)
     *  @param message error message
     *  @param tag "name tag", or exeption type
     *  @param code status code
     *  @param previous  "previous"  exception
     */
    Exception ( const std::string& message  ,
                const std::string& tag      ,
                const StatusCode&  code     ,
                const Exception&   previous ) ;
    // ========================================================================
    /** Constructor (3)
     *  @param message  error message
     *  @param tag      "name tag", or exeption type
     *  @param code     status code
     *  @param previous "previous" exception (used to improve the error message)
     */
    Exception ( const std::string&    message  ,
                const std::string&    tag      ,
                const StatusCode&     code     ,
                const std::exception& previous ) ;
    // ========================================================================
    /// Copy constructor (deep copying!)
    Exception( const Exception& right ) ;
    /// destructor (perform the deletion of "previous" field!)
    virtual ~Exception() throw() {}
    /// assignment operator
    Exception& operator=( const Exception& right ) ;
    /// clone operation
    virtual Exception* clone() const ;
    // ========================================================================
  public: // getters 
    // ========================================================================
    ///  error message to be printed
    const std::string& message   () const { return m_message ; } 
    ///  name tag for the exception, or exception type
    const std::string& tag       () const { return m_tag     ; }
    /// StatusCode for Exception
    const StatusCode&  code      () const { return m_code    ; }
    /// get the previous exception ( "previous" element in the linked list)
    const Exception*   previous  () const { return m_previous.get() ; }
    // ========================================================================
  public: // setters 
    // ========================================================================    
    /// update the error message to be printed
    virtual void setMessage ( const std::string& newMessage ) ;
    /// update name tag
    virtual void setTag     ( const std::string& newTag     ) ;
    ///  update the status code for the exception
    virtual void setCode    ( const StatusCode& newStatus  ) ; 
    // ========================================================================    
  public: // other 
    // ========================================================================
    /// methods  for overloaded printout to std::ostream& and MsgStream&
    virtual std::ostream& fillStream ( std::ostream& os ) const ;
    /// conversion to string 
    virtual std::string   toString   () const ;
    // ========================================================================
  public: // std::exception
    // ========================================================================
    /// method from std::exception
    virtual const char* what () const throw()  { return message().c_str() ; }
    // ========================================================================
  private:
    // ========================================================================
    /// error message 
    std::string       m_message ;  // error message
    /// exception tag 
    std::string       m_tag     ;  // exception tag
    /// status code 
    Ostap::StatusCode m_code    ;  // status code for exception
    /// "previous" element in the linked list
    std::unique_ptr<Exception> m_previous; 
    // ========================================================================
  } ;
  // ==========================================================================
  /// overloaded printout to std::ostream
  inline std::ostream& operator<< ( std::ostream& os , const Exception& e  ) 
  { return e.fillStream ( os ); }
  // ===========================================================================
  inline bool Assert
  ( const bool assertion                                          ,
    const std::string& message                               ,
    const std::string& tag     = "Ostap"                    ,
    const StatusCode&  sc      =  StatusCode::FAILURE )
  { return assertion ? true : throwException ( message , tag , sc ).isSuccess() ; }
  // ===========================================================================
  template < unsigned int N1, unsigned int N2>
  inline bool Assert
  ( const bool assertion      ,
    const char (&message)[N1] ,
    const char (&tag)    [N2] ,
    const StatusCode& sc = StatusCode::FAILURE )
  { return assertion ? true :
      throwException ( std::string ( message , message + N1 ) ,
                       std::string ( tag     , tag     + N2 ) , sc ).isSuccess() ; }
  // ===========================================================================
} //                                                      end of namespace Ostap
// =============================================================================
//                                                                       The END 
// =============================================================================
#endif  // OSTAP_EXCEPTION_H
// =============================================================================
