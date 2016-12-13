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
// constructor: loc
// ============================================================================
Ostap::Math::GSL::GSL_Error_Handler::GSL_Error_Handler () 
  : m_old ( 0 ) 
{ m_old = gsl_set_error_handler ( &GSL_local_error ) ; }
// ============================================================================
// destructor/ unlock 
// ============================================================================
Ostap::Math::GSL::GSL_Error_Handler::~GSL_Error_Handler () 
{ gsl_set_error_handler ( m_old ) ; }

// ============================================================================
// constructor: make use of Gsl Error Handler: print error to stderr 
// ============================================================================
Ostap::Utils::GslError::GslError() 
  : m_previous(0) 
{ m_previous = (void*) gsl_set_error_handler ( &GSL_local_error ) ; }
// ============================================================================
// destructor: stop using the error  handler 
// ============================================================================
Ostap::Utils::GslError::~GslError() 
{ 
  if ( m_previous ) 
  { gsl_set_error_handler ( (gsl_error_handler_t*) m_previous ) ; } 
}

// ============================================================================
// constructor: make use of Gsl Error Handler: print error to stderr 
// ============================================================================
Ostap::Utils::GslIgnore::GslIgnore() 
  : m_previous (0) 
{ m_previous = (void*) gsl_set_error_handler ( &GSL_ignore_error ) ; }
// ============================================================================
// destructor: stop using the error  handler 
// ============================================================================
Ostap::Utils::GslIgnore::~GslIgnore() 
{ 
  if ( m_previous ) 
  { gsl_set_error_handler ( (gsl_error_handler_t*) m_previous ) ; } 
}

// ============================================================================
// constructor: make use of Gsl Error Handler: print error to stderr 
// ============================================================================
Ostap::Utils::GslException::GslException() 
  : m_previous (0) 
{ m_previous = (void*) gsl_set_error_handler ( &GSL_exception_error ) ; }
// ============================================================================
// destructor: stop using the error  handler 
// ============================================================================
Ostap::Utils::GslException::~GslException() 
{ 
  if ( m_previous ) 
  { gsl_set_error_handler ( (gsl_error_handler_t*) m_previous ) ; } 
}

// ============================================================================
// The END 
// ============================================================================
