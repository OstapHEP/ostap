// ============================================================================
#ifndef OSTAP_ERROR2EXCEPTION_H 
#define OSTAP_ERROR2EXCEPTION_H 
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <vector>
#include <string>
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
    /** @class ErrorSentry Ostap/Utils.h
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
      GslError  ( const bool force = false ) ; 
      /// destructor: stop using the error  handler 
      ~GslError () ;
      // ======================================================================
    protected: 
      // ======================================================================
      typedef  void handler ( const char* , const char* , int , int ) ;
      GslError ( handler* h , const bool force = false ) ;
      handler*  m_previous ;
      bool      m_force    ;
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
      GslIgnore  ( const bool force = false ) ;
      // ======================================================================
    };
    // ========================================================================
    /** @class GslCount
     *  helper class to manipulate with GSL error handlers 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     */
    class GslCount : public GslError 
    {
      //=======================================================================
    public:
      //=======================================================================
      /// constructor: make use of Gsl Error Handler: counts errors 
      GslCount ( const bool force = false ) ;
      // ======================================================================
      /// get total number of errors 
      static std::size_t size  () ; // get total number of errors 
      // ======================================================================
      /// clear summary of errors 
      static std::size_t clear () ;
      // ======================================================================
      typedef std::vector<std::string> Row   ;
      typedef std::vector<Row>         Table ;
      // ======================================================================
      /// get all errors in a form of the  table 
      static Table table () ;
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
      GslException  ( const bool force = false ) ; 
      // ======================================================================
    };
    // ========================================================================
  } //                                        The end of namespace Ostap::Utils 
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
#endif // OSTAP_ERROR2EXCEPTION_H
// ============================================================================
//                                                                      The END 
// ============================================================================
