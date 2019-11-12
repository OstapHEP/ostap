// ============================================================================
// Include files 
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Topics.h"
// ============================================================================
// ROOT & Roofit
// ============================================================================
#include "RooMsgService.h"
// ============================================================================
/** @file 
 *  Implementation file for classe:
 *  - Ostap::Utils::AddTopics
 *  - Ostap::Utils::RemoveTopics
 *  @date 2019-11-11 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 */
// ============================================================================
/*  remove topic from the stream 
 *  @see RooMgsService 
 *  @see RooFit::MsgTopic
 *  @see RooFit::MsgLevel 
 *  @return true if topic is removed 
 */
// ============================================================================
bool Ostap::Utils::remove_topic 
( const unsigned short   stream , 
  const unsigned  short  topics , 
  const RooFit::MsgLevel level  ) 
{
  RooMsgService& svc = RooMsgService::instance() ;
  bool removed = false ;
  if ( stream < svc.numStreams() && svc.getStreamStatus ( stream ) )
  {
    RooMsgService::StreamConfig& s = svc.getStream ( stream ) ;
    if ( s.minLevel <= level )
    { 
      for ( unsigned short j = 0 ; j <= 15 ; ++j )
      {
        RooFit::MsgTopic topic = (RooFit::MsgTopic) ( topics & ( 1 << j ) ) ;
        if ( s.topic & topic ) { s.removeTopic ( topic ) ; removed = true ; } 
      }
    }
  }
  return removed ;
}
// ============================================================================
/*  add topic topic from the stream 
 *  @see RooMgsService 
 *  @see RooFit::MsgTopic
 *  @return true if topic is added 
 */
// ============================================================================
bool Ostap::Utils::add_topic 
( const unsigned short stream , 
  const unsigned short topics ) 
{
  RooMsgService& svc = RooMsgService::instance() ;
  bool added = false ;  
  if ( stream < svc.numStreams() && svc.getStreamStatus ( stream ) )
  {
    RooMsgService::StreamConfig& s = svc.getStream ( stream ) ;
      for ( unsigned short j = 0 ; j <= 15 ; ++j )
      {
        RooFit::MsgTopic topic = (RooFit::MsgTopic) ( topics & ( 1 << j ) ) ;
        s.addTopic ( topic ) ; added = true ; 
      }   
  }
  return added ;
}
// ===========================================================================
Ostap::Utils::RemoveTopic::RemoveTopic 
( const unsigned short   topics  ,              
  const RooFit::MsgLevel level   , 
  const int              stream  ) 
  : m_topics  ( topics ) 
  , m_level   ( level  )
  , m_streams () 
{
  RooMsgService& svc = RooMsgService::instance() ;
  svc.saveState () ;
  const int num_streams =  svc.numStreams() ;
  for ( unsigned short i = 0 ; i < num_streams ; ++i )
  {
    if ( stream < 0 || stream == i ) 
    { if ( remove_topic ( i , m_topics , m_level ) ) { m_streams.insert ( i ) ; } }  
  }
}
// =============================================================================
// destructor 
// =============================================================================
Ostap::Utils::RemoveTopic::~RemoveTopic() { exit() ; }
// =============================================================================~
void Ostap::Utils::RemoveTopic::exit ()
{
  if ( m_streams.empty() ) { return ; }
  RooMsgService& svc = RooMsgService::instance() ;
  svc.restoreState () ;
  m_streams.clear() ;
}
// ============================================================================
Ostap::Utils::AddTopic::AddTopic
( const unsigned short   topics ,  
  const int              stream ) 
  : m_topics  ( topics        ) 
  , m_level   ( RooFit::INFO  )
  , m_streams () 
{ 
  RooMsgService& svc = RooMsgService::instance() ;
  svc.saveState () ;
  const int num_streams =  svc.numStreams() ;
  for ( unsigned short i = 0 ; i < num_streams ; ++i )
  {
    if ( stream < 0 || stream == i ) 
    { if  ( add_topic ( i , m_topics ) ) { m_streams.insert ( i ) ; } }
  }
}
// =============================================================================
// destructor 
// =============================================================================
Ostap::Utils::AddTopic::~AddTopic() { exit() ; }
// =============================================================================~
void Ostap::Utils::AddTopic::exit ()
{
  if ( m_streams.empty() ) { return ; }
  RooMsgService& svc = RooMsgService::instance() ;
  svc.restoreState () ;
  m_streams.clear() ;
}
// ============================================================================
//                                                                      The END 
// ============================================================================
