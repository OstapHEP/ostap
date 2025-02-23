// ============================================================================
#ifndef OSTAP_STATUSCODE_H
#define OSTAP_STATUSCODE_H
// ============================================================================
// Incldue files
// ============================================================================
// STD&STL
// ============================================================================
#include <string>
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  /** @class StatusCode Ostap/StatusCode.h
   *  Simplified version of the class StatusCode from Gaudi Project 
   */
  class StatusCode final 
  {
    // ========================================================================
  public: // ==================================================================
    // ========================================================================
    enum ErrorCodes 
      {
        SUCCESS     = 0 ,
        FAILURE     = 1 ,
        RECOVERABLE = 2 ,
      };
    // ========================================================================
  public: // ==================================================================
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
    unsigned long    code () const { return m_code    ; }
    // ========================================================================
  private : // ================================================================
    // ========================================================================
    /// the actual code
    unsigned long m_code ; // the actual code
    // ========================================================================
  }; // =======================================================================
  // ==========================================================================
  // printout 
  std::ostream& operator<< ( std::ostream& s , const StatusCode& sc ) ;
  // ===========================================================================
} //                                                      end of namespace Ostap
// =============================================================================
namespace Ostap
{
  // ===========================================================================
  /** throw the exception
   *  @param message the reason 
   *  @param tag     the tag/category  
   *  @param code    the code 
   *  @param file    filename 
   *  @param line    linenumber 
   */
  StatusCode throwException
  ( const std::string& message                              , 
    const std::string& tag     = "Ostap"                    , 
    const StatusCode&  code    = Ostap::StatusCode::FAILURE ,
    const char*        file    = nullptr                    ,
    const std::size_t  line    = 0                          ) ;
  // ==========================================================================
  /** Assert certain condition
   *  @param assertion condition 
   *  @param message   message 
   *  @param tag       tag/category 
   *  @param file      filename 
   *  @param line      linenumber 
   */ 
  inline bool Assert
  ( const bool         assertion                       ,
    const std::string& message                         ,
    const std::string& tag       = "Ostap"             ,
    const StatusCode&  sc        = StatusCode::FAILURE ,
    const char*        file      = nullptr             ,
    const std::size_t  line      =  0                  )
  { return assertion ? true : throwException ( message , tag , sc , file , line ).isSuccess() ; }
  // ===========================================================================
  /** Assert certain condition
   *  @param assertion condition 
   *  @param message   message 
   *  @param tag       tag/category 
   *  @param file      filename 
   *  @param line      linenumber 
   */ 
  inline bool Assert
  ( const bool         assertion                       ,
    const char*        message                         ,
    const char*        tag       = "Ostap"             ,
    const StatusCode&  sc        = StatusCode::FAILURE ,
    const char*        file      = nullptr             ,
    const std::size_t  line      =  0                  )
  { return assertion ? true :
      throwException ( std::string ( message ) ,
                       std::string ( tag     ) , sc , file , line ).isSuccess() ; }
  // ===========================================================================
  /** Assert certain condition
   *  @param assertion condition 
   *  @param message   message 
   *  @param tag       tag/category 
   *  @param file      filename 
   *  @param line      linenumber 
   */ 
  template < unsigned int N1,
             unsigned int N2>
  inline bool Assert
  ( const bool        assertion                              ,
    const char        (&message)[N1+1]                       ,
    const char        (&tag)    [N2+1]                       ,
    const StatusCode& sc               = StatusCode::FAILURE ,
    const char*       file             = nullptr             ,
    const std::size_t line             =  0                  )
  {
    return assertion ? true :
      throwException ( std::string ( message , message + N1 ) ,
                       std::string ( tag     , tag     + N2 ) , sc , file , line ).isSuccess() ;
  }
  // ===========================================================================
} //                                                      end of namespace Ostap
// =============================================================================
//                                                                       The END 
// =============================================================================
#endif  // OSTAP_STATUSCODES_H
// =============================================================================


