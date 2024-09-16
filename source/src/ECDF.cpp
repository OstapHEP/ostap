// ============================================================================
// Incldue files 
// ============================================================================
// STD&STL
// ============================================================================
#include <memory>
#include <cstring>
// ============================================================================
// Ostap 
// ============================================================================
#include "Ostap/ECDF.h"
// ============================================================================
// Local
// ============================================================================
#include "Exception.h"
// ============================================================================
/** @file 
 *  Implementation file for class Ostap::Math::ECDF
 *  @see Ostap::Math::ECDF
 *  @date 2024-09-16 
 *  @author Vanya BELYAEV 
 */
// ============================================================================
// Standard constructor from  data
// ============================================================================
Ostap::Math::ECDF::ECDF
( const Ostap::Math::ECDF::Data&  data  )
  : m_data ( data ) 
{
  this->sort_me() ;
}
// ============================================================================
// sort container 
// ============================================================================
void Ostap::Math::ECDF::sort_me()
{
  Ostap::Assert ( 0 < m_data.size()        ,
                  "Empty data constainer!" ,
                  "Ostap::Math::ECDF"      ) ;
  // sort it! 
  std::sort ( m_data.begin() , m_data.end  () ) ;
}
// ============================================================================
// swap two objects 
// ============================================================================
void Ostap::Math::ECDF::swap ( Ostap::Math::ECDF& right )
{ std::swap ( m_data , right.m_data ) ; }
// ============================================================================
// the main method 
// ============================================================================
double Ostap::Math::ECDF::evaluate   ( const double x ) const
{
  return
    ( x < m_data.front () ) ? 0.0 : 
    ( x > m_data.back  () ) ? 1.0 :
    double ( std::upper_bound ( m_data.begin () ,
                                m_data.end   () , x ) - m_data.begin() ) / m_data.size () ;
}
// ============================================================================
// add a value to data container  
// ============================================================================
unsigned long
Ostap::Math::ECDF::add
( const double value  )
{
  auto where = std::upper_bound ( m_data.begin () , m_data.end   () , value ) ;
  m_data.insert ( where , value ) ;
  return m_data.size () ;
}
// ============================================================================
// add a value to data container  
// ============================================================================
unsigned long
Ostap::Math::ECDF::add
( const Ostap::Math::ECDF::Data& values )
{
  return  add ( values.begin() , values.end() ) ;
}
// ============================================================================
// add a value to data container  
// ============================================================================
unsigned long
Ostap::Math::ECDF::add
( const Ostap::Math::ECDF& values )
{
  Data tmp2  ( m_data.size() + values.size() ) ;
  /// merge two sorted continers 
  std::merge ( m_data.begin        () ,
               m_data.end          () ,
               values.m_data.begin () ,
               values.m_data.end   () ,
               tmp2.begin          () ) ;
  std::swap ( m_data , tmp2 ) ;
  return m_data.size () ;
}


// ============================================================================
//                                                                      The END 
// ============================================================================

