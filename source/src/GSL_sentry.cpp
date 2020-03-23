// $Id$
// ============================================================================
// Include files 
// ============================================================================
// STD & STL 
// ============================================================================
#include <iostream>
#include <sstream>
// ============================================================================
// Ostap 
// ============================================================================
#include "Ostap/Error2Exception.h"
// ============================================================================
// local
// ============================================================================
#include "GSL_sentry.h"
#include "Exception.h"
// ============================================================================
/** @file 
 *  Implementation file for class Ostap::Math::GSL::GSL_Error_Handler
 *  
 *  @see Ostap::Math::GSL::GSL_Error_Handler
 *  @date 2012-05-27 
 *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
 *
 *                    $Revision$
 *  Last modification $Date$
 *                 by $Author$
 */
// ============================================================================
namespace
{
  // ==========================================================================
  // GSL 
  // ==========================================================================
  void GSL_local_error
  ( const char * reason    ,
    const char * file      ,
    int          line      ,
    int          gsl_errno ) 
  {
    std::cerr 
      << " GSL_ERROR : "   
      << gsl_errno << "/'" << gsl_strerror ( gsl_errno ) << "'"
      << "\t reason '"     
      << reason    << "' "
      << "\t file/line '"  
      << file      << "'/" << line 
      << std::endl ;  
  }
  // ==========================================================================
  void GSL_ignore_error
  ( const char * /* reason    */ ,
    const char * /* file      */ ,
    int          /* line      */ ,
    int          /* gsl_errno */ ) {}
  // ==========================================================================
  void GSL_exception_error
  ( const char * reason    ,
    const char * file      ,
    int          line      ,
    int          gsl_errno ) 
  {
    std::string tag = "GSL/Error" ;
    std::ostringstream ss ;
    ss << gsl_strerror ( gsl_errno ) << "(" << gsl_errno << ") "
       << reason
       << " in "     << file   << " at line " << line ;
    Ostap::throwException ( tag + ": " + ss.str() , 
                            tag                   , 
                            Ostap::StatusCode ( 20000 + gsl_errno ) ) ;
  }
}
// ============================================================================
// constructor: make use of Gsl Error Handler: print error to stderr 
// ============================================================================
Ostap::Utils::GslError::GslError ( Ostap::Utils::GslError::handler* h )
  : m_previous ( gsl_set_error_handler ( h ) ) 
{ 
  static_assert( std::is_same<handler,gsl_error_handler_t>::value  ,
                 "``handler'' type is not ``gsl_error_handler_t''" ) ;
}
// ============================================================================
// constructor: make use of Gsl Error Handler: print error to stderr 
// ============================================================================
Ostap::Utils::GslError::GslError () : GslError ( &GSL_local_error ) {}
// ============================================================================
// destructor: stop using the error  handler 
// ============================================================================
Ostap::Utils::GslError::~GslError() { gsl_set_error_handler ( m_previous ) ; }
// ============================================================================
// constructor: make use of Gsl Error Handler: print error to stderr 
// ============================================================================
Ostap::Utils::GslIgnore::GslIgnore() : GslError( &GSL_ignore_error ) {}
// ============================================================================
// constructor: make use of Gsl Error Handler: print error to stderr 
// ============================================================================
Ostap::Utils::GslException::GslException() : GslError( &GSL_exception_error ) {}
// ============================================================================
// The END 
// ============================================================================
