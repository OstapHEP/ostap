// ============================================================================
#ifndef OSTAP_ERROR2EXCEPTION_H 
#define OSTAP_ERROR2EXCEPTION_H 
// ============================================================================
// Include files
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Utils 
  {
    // ========================================================================
    /// use local error handler for ROOT 
    bool useErrorHandler ( const bool use = true ) ;
    // ========================================================================
    /** @class ErrorHandlerSentry Error2Exception.h Ostap/Error2Exception.h
     *  Simple error handler for ROOT that converts error messages into exceptions
     *  @author Vanya Belyaev Ivan.Belayev@itep.ru
     *  @date   2016-12-10
     */
    // ========================================================================
    class ErrorSentry
    {
    public:
      // =======================================================================
      /// constructor: make use of local error handler
      ErrorSentry  () ;
      /// destructor: stop using local error handler
      ~ErrorSentry () ;
      // =======================================================================
    private : 
      // =======================================================================
      /// is actually used? 
      bool m_previous ; // is actually used? 
      // ======================================================================      
    } ;
    // ========================================================================
    /** @class GslError
     *  helper class to manipulate with GSL error handlers 
     *  @author Vanya BELAYEV Ivan.Belyaev@itep.ru
     */
    class GslError
    {
      //=======================================================================
    public:
      //=======================================================================
      /// constructor: make use of Gsl Error Handler: print error to stderr 
      GslError () ; 
      /// destructor: stop using the error  handler 
      ~GslError() ;
      // ======================================================================
    private:
      // ======================================================================
      void*  m_previous ;
      // ======================================================================
    };
    // ========================================================================
    /** @class GslIgnore
     *  helper class to manipulate with GSL error handlers 
     *  @author Vanya BELAYEV Ivan.Belyaev@itep.ru
     */
    class GslIgnore
    {
      //=======================================================================
    public:
      //=======================================================================
      /// constructor: make use of Gsl Error Handler: ignore error 
      GslIgnore  () ; 
      /// destructor: stop using the error  handler 
      ~GslIgnore () ;
      // ======================================================================
    private:
      // ======================================================================
      void*  m_previous ;
      // ======================================================================
    };
    // ========================================================================
    /** @class GslException
     *  helper class to manipulate with GSL error handlers : 
     *  throw exception on error 
     *  @author Vanya BELAYEV Ivan.Belyaev@itep.ru
     */
    class GslException
    {
      //=======================================================================
    public:
      //=======================================================================
      /// constructor: make use of Gsl Error Handler: throw Exception 
      GslException  () ; 
      /// destructor: stop using the error  handler 
      ~GslException () ;
      // ======================================================================
    private:
      // ======================================================================
      void*  m_previous ;
      // ======================================================================
    };
    // ========================================================================
  }    
  // ==========================================================================
}
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_ERROR2EXCEPTION_H
// ============================================================================
