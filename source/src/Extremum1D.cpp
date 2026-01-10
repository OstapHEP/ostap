// ============================================================================
// Include files 
// ============================================================================
// Ostap
// ============================================================================
// GSL
// ============================================================================
#include "gsl/gsl_errno.h"
#include "gsl/gsl_min.h"
// ============================================================================
// Local stuff 
// ============================================================================
#include "Extremum1D.h" 
#include "status_codes.h"
// ============================================================================
/** @file
 *  Implementation file for Ostap::Math::GSL::Extremum1D
 *  @see Ostap::Math::GSL::Extremum1D 
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch 
 *  @date 2019-05-14
 */
// ============================================================================
// constructor: allocate GSL minimizer 
// ============================================================================
Ostap::Math::GSL::Minimizer::Minimizer
( const gsl_function*            fun   , 
  const double                   guess , 
  const double                   low   , 
  const double                   high  ,  
  const gsl_min_fminimizer_type* mtype )
  : m_minimizer ( nullptr )
{
 
  Ostap::Assert  ( fun                                      , 
                    "Invalid GSL function"                  , 
                    "Ostap::Math::GSL::Minimizer"           , 
                    INVALID_FUNCTION , __FILE__ , __LINE__  ) ;  

  Ostap::Assert  ( low < high  , 
                   "Invalid low/high parameters!"            , 
                   "Ostap::Math::GSL::Minimizer"             , 
                   INVALID_PARAMETERS , __FILE__  , __LINE__ ) ; 
                   
  if ( !mtype ) { mtype = gsl_min_fminimizer_brent ; }
  m_minimizer = gsl_min_fminimizer_alloc ( mtype ) ;
               
  Ostap::Assert  (  m_minimizer                              , 
                    "Invalid GSL minimizer!"                 , 
                    "Ostap::Math::GsL::Minimizer"            , 
                    INVALID_MINIMIZER  , __FILE__ , __LINE__ ) ;

  // initialize it!
  const int status = gsl_min_fminimizer_set ( m_minimizer , 
                             const_cast<gsl_function*> ( fun ) , 
                             guess , low , high ) ; 

  Ostap::Assert ( GSL_SUCCESS == status         , 
                  "Cannot set GSL minimizer"    , 
                  "Ostap::Math::GSL::Minimizer" , 
                   ERROR_GSL + status , __FILE__ , __LINE__ ) ;
};

// ============================================================================
// destructor: deallocate the minimier 
// ============================================================================
Ostap::Math::GSL::Minimizer::~Minimizer() 
{
  if ( m_minimizer ) { gsl_min_fminimizer_free ( m_minimizer ) ; }
  m_minimizer = nullptr ;  
}
// ============================================================================
//                                                                      The END 
// ============================================================================
