// ============================================================================
#ifndef OSTAP_GSL_UTILS_H 
#define OSTAP_GSL_UTILS_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <ostream>
// ============================================================================
// GSL 
// ============================================================================
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
// =============================================================================
/** @file Ostap/GSL_utils.h
 *  utilities for GSL 
 */
// =============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Utils 
  {
    // ========================================================================
    /** print GSL-vector to the stream 
     *  @param v the vector 
     *  @param s the stream 
     *  @return the stream 
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2012-05-28
     */
    std::ostream& toStream 
    ( const gsl_vector&  v , 
      std::ostream&      s ) ;
    // ========================================================================
    /** print GSL-matrix to the stream 
     *  @param m the matrix 
     *  @param s the stream 
     *  @return the stream 
     *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
     *  @date 2012-05-28
     */
    std::ostream& toStream 
    ( const gsl_matrix&  m , 
      std::ostream&      s ) ;
    // ========================================================================
  } //                                            end of namespace Ostap::Utils
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
/// print operator 
inline std::ostream& operator<<( std::ostream& s , const gsl_vector& v ) 
{ return Ostap::Utils::toStream ( v , s ) ; }
// ============================================================================
/// print operator 
inline std::ostream& operator<<( std::ostream& s , const gsl_matrix& m ) 
{ return Ostap::Utils::toStream ( m , s ) ; }
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_GSL_UTILS_H
// ============================================================================
