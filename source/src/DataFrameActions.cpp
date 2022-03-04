// ============================================================================
// Include files 
// ============================================================================
// STD&STL 
// ============================================================================
#include <memory>
// ============================================================================
// Ostap 
// ============================================================================
#include "Ostap/DataFrame.h"
#include "Ostap/DataFrameActions.h"
// ============================================================================
// Local
// ============================================================================
#include "OstapDataFrame.h"
// ============================================================================
/// ONLY starting from ROOT 6.16
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,16,0)
// ============================================================================
/** @file 
 *  Implementation file for objects from file Ostap/DataFrameActions.h 
 *  @date 2021-06-09 
 *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
 */
// ============================================================================
// constructor 
// ============================================================================
ROOT::Detail::RDF::StatVar::StatVar ()
  : m_result ( std::make_shared<Ostap::StatEntity>() ) 
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
  , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetThreadPoolSize     () ) : 1u )
#else 
  , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetImplicitMTPoolSize () ) : 1u )
#endif
  , m_slots  ( this->m_N ) 
{}
// ============================================================================
// Finalize 
// ============================================================================
void ROOT::Detail::RDF::StatVar::Finalize() 
{ 
  Result_t sum { m_slots [0] } ;
  for ( unsigned int i = 1 ; i < m_N ; ++i ) { sum += m_slots [ i ] ; }
  *m_result = sum ;
}
// ============================================================================
// constructor 
// ============================================================================
ROOT::Detail::RDF::WStatVar::WStatVar ()
  : m_result ( std::make_shared<Ostap::WStatEntity>() ) 
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
  , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetThreadPoolSize     () ) : 1u )
#else 
  , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetImplicitMTPoolSize () ) : 1u )
#endif
  , m_slots  ( this->m_N ) 
{}
// ============================================================================
// Finalize 
// ============================================================================
void ROOT::Detail::RDF::WStatVar::Finalize() 
{ 
  Result_t sum { m_slots [0] } ;
  for ( unsigned int i = 1 ; i < m_N ; ++i ) { sum += m_slots [ i ] ; }
  *m_result = sum ;
}
// ============================================================================

// ============================================================================
#endif // #if ROOT_VERSION_CODE >= ROOT_VERSION(6,16,0)
// ============================================================================

// ============================================================================
//                                                                      The END 
// ============================================================================
