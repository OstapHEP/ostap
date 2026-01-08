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
( const gsl_min_fminimizer_type* mtype )
  : m_minimizer ( nullptr )
{
  if ( !mtype ) { mtype = gsl_min_fminimizer_brent ; }
  m_minimizer = gsl_min_fminimizer_alloc ( mtype ) ;
}
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
