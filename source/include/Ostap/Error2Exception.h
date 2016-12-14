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
    /** @class ErrorSentry 
     *  Simple error handler for ROOT that converts error messages into exceptions
     *  @author Vanya Belyaev Ivan.Belayev@itep.ru
     *  @date   2016-12-10
     */
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
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     */
    class GslError
    {
      //=======================================================================
    public:
      //=======================================================================
      /// constructor: make use of Gsl Error Handler: print error to stderr 
      GslError  () ; 
      /// destructor: stop using the error  handler 
      ~GslError () ;
      // ======================================================================
    protected: 
      // ======================================================================
      typedef  void handler ( const char* , const char* , int , int ) ;
      GslError ( handler* h ) ;
      handler*  m_previous ;
      // ======================================================================
    };
    // ========================================================================
    /** @class GslIgnore
     *  helper class to manipulate with GSL error handlers 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     */
    class GslIgnore : public GslError 
    {
      //=======================================================================
    public:
      //=======================================================================
      /// constructor: make use of Gsl Error Handler: ignore error 
      GslIgnore  () ; // constructor: make use of Gsl Error Handler: ignore error 
      // ======================================================================
    };
    // ========================================================================
    /** @class GslException
     *  helper class to manipulate with GSL error handlers : 
     *  throw exception on error 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     */
    class GslException : public GslError
    {
      //=======================================================================
    public:
      //=======================================================================
      /// constructor: make use of Gsl Error Handler: throw Exception 
      GslException  () ; // constructor: make use of Gsl Error Handler: throw Exception 
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
