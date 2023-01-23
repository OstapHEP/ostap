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
// constructor 
// ============================================================================
ROOT::Detail::RDF::LegendrePoly::LegendrePoly
( const unsigned short N    , 
  const double         xmin ,
  const double         xmax )
  : m_result ( std::make_shared<Result_t> ( N , xmin , xmax ) )
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
  , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetThreadPoolSize     () ) : 1u )
#else 
  , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetImplicitMTPoolSize () ) : 1u )
#endif
  , m_slots  ( this->m_N , *(this->m_result.get() ) ) 
{}
// ============================================================================
// constructor 
// ============================================================================
ROOT::Detail::RDF::LegendrePoly::LegendrePoly
( const Ostap::Math::LegendreSum& p )
  : LegendrePoly ( p.degree() , p.xmin () , p.xmax () )
{}
// ============================================================================
// Finalize 
// ============================================================================
void ROOT::Detail::RDF::LegendrePoly::Finalize() 
{ 
  (*m_result) *= 0.0 ;
  for ( unsigned int i = 0 ; i < m_N ; ++i ) { *m_result += m_slots [ i ] ; }
}
// ============================================================================
// constructor 
// ============================================================================
ROOT::Detail::RDF::ChebyshevPoly::ChebyshevPoly
( const unsigned short N    , 
  const double         xmin ,
  const double         xmax )
  : m_result ( std::make_shared<Result_t> ( N , xmin , xmax ) )
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
  , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetThreadPoolSize     () ) : 1u )
#else 
  , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetImplicitMTPoolSize () ) : 1u )
#endif
  , m_slots  ( this->m_N , *(this->m_result.get() ) ) 
{}
// ============================================================================
// constructor 
// ============================================================================
ROOT::Detail::RDF::ChebyshevPoly::ChebyshevPoly
( const Ostap::Math::ChebyshevSum& p )
  : ChebyshevPoly ( p.degree() , p.xmin () , p.xmax () )
{}
// ============================================================================
// Finalize 
// ============================================================================
void ROOT::Detail::RDF::ChebyshevPoly::Finalize() 
{ 
  (*m_result) *= 0.0 ;
  for ( unsigned int i = 0 ; i < m_N ; ++i ) { *m_result += m_slots [ i ] ; }
}
// ============================================================================
// constructor 
// ============================================================================
ROOT::Detail::RDF::BernsteinPoly::BernsteinPoly
( const unsigned short N    , 
  const double         xmin ,
  const double         xmax )
  : m_result ( std::make_shared<Result_t> ( N , xmin , xmax ) )
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
  , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetThreadPoolSize     () ) : 1u )
#else 
  , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetImplicitMTPoolSize () ) : 1u )
#endif
  , m_slots  ( this->m_N , *(this->m_result.get() ) ) 
{}
// ============================================================================
// constructor 
// ============================================================================
ROOT::Detail::RDF::BernsteinPoly::BernsteinPoly
( const Ostap::Math::Bernstein& p )
  : BernsteinPoly ( p.degree() , p.xmin () , p.xmax () )
{}
// ============================================================================
// Finalize 
// ============================================================================
void ROOT::Detail::RDF::BernsteinPoly::Finalize() 
{ 
  (*m_result) *= 0.0 ;
  for ( unsigned int i = 0 ; i < m_N ; ++i ) { *m_result += m_slots [ i ] ; }
}
// ============================================================================
// constructor 
// ============================================================================
ROOT::Detail::RDF::LegendrePoly2::LegendrePoly2
( const unsigned short NX   , 
  const unsigned short NY   , 
  const double         xmin ,
  const double         xmax ,
  const double         ymin ,
  const double         ymax )
  : m_result ( std::make_shared<Result_t> ( NX , NY , xmin , xmax , ymin , ymax ) )
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
  , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetThreadPoolSize     () ) : 1u )
#else 
  , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetImplicitMTPoolSize () ) : 1u )
#endif
  , m_slots  ( this->m_N , *(this->m_result.get() ) ) 
{}
// ============================================================================
// constructor 
// ============================================================================
ROOT::Detail::RDF::LegendrePoly2::LegendrePoly2
( const Ostap::Math::LegendreSum2& p )
  : LegendrePoly2 ( p.nX   () , p.nY   () , 
                    p.xmin () , p.xmax () , 
                    p.ymin () , p.ymax () )
{}
// ============================================================================
// Finalize 
// ============================================================================
void ROOT::Detail::RDF::LegendrePoly2::Finalize() 
{ 
  (*m_result) *= 0.0 ;
  for ( unsigned int i = 0 ; i < m_N ; ++i ) { *m_result += m_slots [ i ] ; }
}
// ============================================================================
// constructor 
// ============================================================================
ROOT::Detail::RDF::BernsteinPoly2::BernsteinPoly2
( const unsigned short NX   , 
  const unsigned short NY   , 
  const double         xmin ,
  const double         xmax ,
  const double         ymin ,
  const double         ymax )
  : m_result ( std::make_shared<Result_t> ( NX , NY , xmin , xmax , ymin , ymax ) )
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
  , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetThreadPoolSize     () ) : 1u )
#else 
  , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetImplicitMTPoolSize () ) : 1u )
#endif
  , m_slots  ( this->m_N , *(this->m_result.get() ) ) 
{}
// ============================================================================
// constructor 
// ============================================================================
ROOT::Detail::RDF::BernsteinPoly2::BernsteinPoly2
( const Ostap::Math::Bernstein2D& p )
  : BernsteinPoly2 ( p.nX   () , p.nY   () , 
                     p.xmin () , p.xmax () , 
                     p.ymin () , p.ymax () )
{}
// ============================================================================
// Finalize 
// ============================================================================
void ROOT::Detail::RDF::BernsteinPoly2::Finalize() 
{ 
  (*m_result) *= 0.0 ;
  for ( unsigned int i = 0 ; i < m_N ; ++i ) { *m_result += m_slots [ i ] ; }
}
// ============================================================================
// constructor 
// ============================================================================
ROOT::Detail::RDF::LegendrePoly3::LegendrePoly3
( const unsigned short NX   , 
  const unsigned short NY   , 
  const unsigned short NZ   , 
  const double         xmin ,
  const double         xmax ,
  const double         ymin ,
  const double         ymax ,
  const double         zmin ,
  const double         zmax )
  : m_result ( std::make_shared<Result_t> ( NX   , NY   , NZ , 
                                            xmin , xmax , 
                                            ymin , ymax , 
                                            zmin , zmax ) )
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
  , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetThreadPoolSize     () ) : 1u )
#else 
  , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetImplicitMTPoolSize () ) : 1u )
#endif
  , m_slots  ( this->m_N , *(this->m_result.get() ) ) 
{}
// ============================================================================
// constructor 
// ============================================================================
ROOT::Detail::RDF::LegendrePoly3::LegendrePoly3
( const Ostap::Math::LegendreSum3& p )
  : LegendrePoly3 ( p.nX   () , p.nY   () , p.nZ () ,  
                    p.xmin () , p.xmax () , 
                    p.ymin () , p.ymax () ,
                    p.zmin () , p.zmax () )
{}
// ============================================================================
// Finalize 
// ============================================================================
void ROOT::Detail::RDF::LegendrePoly3::Finalize() 
{ 
  (*m_result) *= 0.0 ;
  for ( unsigned int i = 0 ; i < m_N ; ++i ) { *m_result += m_slots [ i ] ; }
}
// ============================================================================
// constructor 
// ============================================================================
ROOT::Detail::RDF::BernsteinPoly3::BernsteinPoly3
( const unsigned short NX   , 
  const unsigned short NY   , 
  const unsigned short NZ   , 
  const double         xmin ,
  const double         xmax ,
  const double         ymin ,
  const double         ymax ,
  const double         zmin ,
  const double         zmax )
  : m_result ( std::make_shared<Result_t> ( NX   , NY   , NZ , 
                                            xmin , xmax , 
                                            ymin , ymax ,
                                            zmin , zmax ) )
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
  , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetThreadPoolSize     () ) : 1u )
#else 
  , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetImplicitMTPoolSize () ) : 1u )
#endif
  , m_slots  ( this->m_N , *(this->m_result.get() ) ) 
{}
// ============================================================================
// constructor 
// ============================================================================
ROOT::Detail::RDF::BernsteinPoly3::BernsteinPoly3
( const Ostap::Math::Bernstein3D& p )
  : BernsteinPoly3 ( p.nX   () , p.nY   () , p.nZ() , 
                     p.xmin () , p.xmax () , 
                     p.ymin () , p.ymax () ,
                     p.zmin () , p.zmax () )
{}
// ============================================================================
// Finalize 
// ============================================================================
void ROOT::Detail::RDF::BernsteinPoly3::Finalize() 
{ 
  (*m_result) *= 0.0 ;
  for ( unsigned int i = 0 ; i < m_N ; ++i ) { *m_result += m_slots [ i ] ; }
}
// ============================================================================
// constructor 
// ============================================================================
ROOT::Detail::RDF::LegendrePoly4::LegendrePoly4
( const unsigned short NX   , 
  const unsigned short NY   , 
  const unsigned short NZ   , 
  const unsigned short NU   , 
  const double         xmin ,
  const double         xmax ,
  const double         ymin ,
  const double         ymax ,
  const double         zmin ,
  const double         zmax ,
  const double         umin ,
  const double         umax )
  : m_result ( std::make_shared<Result_t> ( NX   , NY   , NZ , NU , 
                                            xmin , xmax , 
                                            ymin , ymax , 
                                            zmin , zmax ,
                                            umin , umax ) )
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,22,0)
  , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetThreadPoolSize     () ) : 1u )
#else 
  , m_N      ( ROOT::IsImplicitMTEnabled() ? std::max ( 1u , ROOT::GetImplicitMTPoolSize () ) : 1u )
#endif
  , m_slots  ( this->m_N , *(this->m_result.get() ) ) 
{}
// ============================================================================
// constructor 
// ============================================================================
ROOT::Detail::RDF::LegendrePoly4::LegendrePoly4
( const Ostap::Math::LegendreSum4& p )
  : LegendrePoly4 ( p.nX   () , p.nY   () , 
                    p.nZ   () , p.nU   () , 
                    p.xmin () , p.xmax () , 
                    p.ymin () , p.ymax () ,
                    p.zmin () , p.zmax () ,
                    p.umin () , p.umax () )
{}
// ============================================================================
// Finalize 
// ============================================================================
void ROOT::Detail::RDF::LegendrePoly4::Finalize() 
{ 
  (*m_result) *= 0.0 ;
  for ( unsigned int i = 0 ; i < m_N ; ++i ) { *m_result += m_slots [ i ] ; }
}





// ============================================================================
#endif // #if ROOT_VERSION_CODE >= ROOT_VERSION(6,16,0)
// ============================================================================

// ============================================================================
//                                                                      The END 
// ============================================================================
