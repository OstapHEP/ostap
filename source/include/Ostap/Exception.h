// ============================================================================
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
     *  @param  file file name 
     *  @param line line numbr 
     */
    Exception
    ( const std::string& message                             ,
      const std::string& tag                                 ,
      const StatusCode&  code   = Ostap::StatusCode::FAILURE , 
      const char*        file   = nullptr                    ,
      const std::size_t  line   = 0                          ) ;
    // ========================================================================
    /** Constructor (2)
     *  @param message error message
     *  @param tag "name tag", or exeption type
     *  @param code status code
     *  @param previous  "previous"  exception
     *  @param  file file name 
     *  @param line line numbr 
     */
    Exception
    ( const std::string& message                             ,
      const std::string& tag                                 ,
      const Exception&   previous                            , 
      const StatusCode&  code   = Ostap::StatusCode::FAILURE , 
      const char*        file   = nullptr                    ,
      const std::size_t  line   = 0                          ) ;
    // ========================================================================
    /** Constructor (3)
     *  @param message  error message
     *  @param tag      "name tag", or exeption type
     *  @param code     status code
     *  @param previous "previous" exception (used to improve the error message)
     *  @param  file file name 
     *  @param line line numbr 
     */
    Exception
    ( const std::string&    message                             ,
      const std::string&    tag                                 ,
      const std::exception& previous                            , 
      const StatusCode&     code   = Ostap::StatusCode::FAILURE , 
      const char*           file   = nullptr                    ,
      const std::size_t     line   = 0                          ) ;
    // ========================================================================
    /** constructor from std:exception 
     *  @param exc exception 
     *  @param  file file name 
     *  @param line line numbr 
     */
    Exception
    ( const std::exception& exc             ,
      const char*           file  = nullptr ,
      const std::size_t     line  = 0       ) ;
    /// Copy constructor (deep copying!)
    Exception ( const Exception& right ) ;
    /// destructor (perform the deletion of "previous" field!)
    virtual ~Exception() throw() {}
    /// clone operation
    virtual Exception* clone() const ;
    // ========================================================================
  public: // standard getters 
    // ========================================================================
    /// const std::string&       message  () const { return m_message ; }
    const std::string&       tag      () const { return m_tag     ; }
    const Ostap::StatusCode& code     () const { return m_code    ; }
    /// previousescepintin the chain 
    const Exception*         previous () const { return m_previous.get() ; } 
    /// index of this exception in the chain 
    unsigned int             index    () const
    { return m_previous ? 1 + m_previous->index() : 0 ;  }
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
    //. const char* what () const throw () override { return m_what.c_str() ; }
    //  const char* what () const noexcept override { return m_what.c_str() ; }
    // ========================================================================
#if defined ( __cplusplus ) && ( 201103L <= __cplusplus ) // ==================
    // ========================================================================
    const char* what () const noexcept override { return m_what.c_str() ; }
    // ========================================================================
#else // ======================================================================
    const char* what () const throw () override { return m_what.c_str() ; }
    // ========================================================================
#endif  // ====================================================================
    // ========================================================================
  private:
    // ========================================================================
    /// error message 
    std::string                m_message  {}    ;  // error message
    /// exception tag 
    std::string                m_tag      {}    ;  // exception tag
    /// status code 
    Ostap::StatusCode          m_code { Ostap::StatusCode::FAILURE } ;  // status code for exception
    /// file name
    std::string                m_file     {}    ;
    /// line number 
    std::size_t                m_line     { 0 } ;
    /// "what"
    std::string                m_what     {   } ; 
    /// "previous" element in the linked list
    std::unique_ptr<Exception> m_previous {} ;
    // ========================================================================
  } ;
  // ==========================================================================
  /// overloaded printout to std::ostream
  inline std::ostream& operator<< ( std::ostream& os , const Exception& e  ) 
  { return e.fillStream ( os ); }
  // ===========================================================================
} //                                                      end of namespace Ostap
// =============================================================================
//                                                                       The END 
// =============================================================================
#endif  // OSTAP_EXCEPTION_H
// =============================================================================
