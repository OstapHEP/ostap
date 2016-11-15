// $Id$
// ============================================================================
// Include files 
// ============================================================================
// GSL 
// ============================================================================
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/ToStream.h"
#include "Ostap/GSL_utils.h"
// ============================================================================
// local
// ============================================================================
#include "format.h"
// ============================================================================
/** @file 
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2012-05-28
 *                    $Revision$
 *  Last modification $Date$
 *  by                $Author$
 */
// ============================================================================
/* print GSL-vector to the stream 
 *  @param v the vector 
 *  @param s the stream 
 *  @return the stream 
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2012-05-28
 */
// ============================================================================
std::ostream& 
Ostap::Utils::toStream 
( const gsl_vector&  v , 
  std::ostream&      s ) 
{
  s << "[ " ;
  for ( std::size_t i = 0 ; i < v.size ; ++i ) 
  {
    if ( 0 != i ) { s << " , " ; } 
    s << gsl_vector_get ( &v , i ) ;
  }
  s << "]" ;
  return s ;
}
// ============================================================================
/*  print GSL-matrix to the stream 
 *  @param m the matrix 
 *  @param s the stream 
 *  @return the stream 
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2012-05-28
 */
// ============================================================================
std::ostream& 
Ostap::Utils::toStream 
( const gsl_matrix&  m , 
  std::ostream&      s ) 
{
  //
  for ( std::size_t i = 0 ; i < m.size1 ; ++i ) 
  {
    s << " | " ;
    for ( std::size_t j = 0 ; j < m.size2 ; ++j ) 
    {
      if ( 0 != j ) { s << ", " ; }
      //    
      s << Ostap::format ( "%|11.5g|" , gsl_matrix_get ( &m , i , j ) ) ;
    }
    s << " | "<< std::endl ;
  }
  //
  return s ;  
}
// ============================================================================
// The END 
// ============================================================================
