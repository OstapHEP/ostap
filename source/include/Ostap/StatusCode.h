// ============================================================================
#ifndef OSTAP_STATUSCODE_H
#define OSTAP_STATUSCODE_H
// ============================================================================
// Incldue files
// ============================================================================
// STD&STL
// ============================================================================
#include <ostream>
#include <memory>
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  /** @class StatusCode StatusCode.h
   *  Simplified version of StatusCode form Gaudi Project 
   */
  class StatusCode final 
  {
  public:
    // ========================================================================
    enum
    {
      FAILURE     = 0,
      SUCCESS     = 1,
      RECOVERABLE = 2
    };
    // ========================================================================
  public:
    // ========================================================================
    /// Constructor.
    StatusCode ( unsigned long code =  SUCCESS ) : m_code ( code ) {}
    // ========================================================================
  public: // check 
    // ========================================================================
    bool isSuccess        () const { return m_code == SUCCESS     ; }
    bool isFailure        () const { return !isSuccess()          ; }  // NB!!
    bool isRecoverable    () const { return m_code == RECOVERABLE ; }
    // ========================================================================
  public: // getters 
    // ========================================================================
    /// Get the status code by value.
    unsigned long getCode () const { return m_code    ; }
    unsigned long    code () const { return m_code    ; }
    operator unsigned long() const { return getCode() ; }
    // ========================================================================
  public: // setters 
    // ========================================================================
    /// Set the status code by value.
    void setCode( unsigned long value )  { m_code = value; }
    /// assignement from the value 
    StatusCode& operator=(unsigned long value) { setCode(value); return *this; }
    // ========================================================================    
  public: // comparison operators 
    // ========================================================================    
    bool operator< ( const StatusCode& b ) const { return m_code <  b.code() ; }
    bool operator<=( const StatusCode& b ) const { return m_code <= b.code() ; }
    bool operator> ( const StatusCode& b ) const { return m_code >  b.code() ; }
    bool operator>=( const StatusCode& b ) const { return m_code >= b.code() ; }
    bool operator!=( const StatusCode& b ) const { return m_code != b.code() ; }
    bool operator==( const StatusCode& b ) const { return m_code == b.code() ; }
    // ========================================================================    
  private :
    // ========================================================================
    /// the actual code
    unsigned long m_code ; // the actual code
    // ========================================================================
  };
  // ==========================================================================
  // printout 
  std::ostream& operator<< ( std::ostream& s , const StatusCode& sc ) ;
  // ===========================================================================
} //                                                      end of namespace Ostap
// =============================================================================
//                                                                       The END 
// =============================================================================
#endif  // OSTAP_STATUSCODES_H
// =============================================================================


