// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <iostream>
#include <string>
#if defined ( __cplusplus ) && ( 201103L <= __cplusplus ) 
#include <chrono>
#endif
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Math.h"
#include "Ostap/ProgressBar.h"
// ============================================================================
// Lcoal 
// ============================================================================
#include "format.h"
// ============================================================================
/** @file
 *  implementation file for class Ostap::Utils::ProgressBar
 *  @see Ostap::Utils::ProgressBar
 *  @see Ostap::Utils::ProgressConf
 *  @date 2022-09-03 
 *  @author Vanya Belyaev Ivan/Belyaev@itep.ru
 */
// ======================================================================
/** full constructor
 *  @param symbol     symbol to show 
 *  @param empty      empty symbol
 *  @param left       left/prefix 
 *  @param right      right/prefix 
 *  @param width      the width 
 *  @param use_timer  use the timer?
 */
// ======================================================================
Ostap::Utils::ProgressConf::ProgressConf
( const unsigned short     width     ,
  const std::string&       symbol    ,
  const std::string&       empty     ,
  const std::string&       left      ,
  const std::string&       right     ,
  const std::string&       what      ,
  const bool               use_timer ,
  const bool               atty      ) 
  : m_width     ( width     )
  , m_symbol    ( symbol    ) 
  , m_empty     ( empty     ) 
  , m_left      ( left      ) 
  , m_right     ( right     ) 
  , m_what      ( what      ) 
  , m_use_timer ( use_timer )
  , m_atty      ( atty      )
{
  // ==========================================================================
  if ( m_symbol.empty() ) { m_symbol = "#" ; }
  if ( m_empty .empty() ) { m_empty  = " " ; }                          \
  // ==========================================================================
  setWidth    ( width     ) ;
  setUseTimer ( use_timer ) ;
  // ==========================================================================
}
// ==========================================================================
/// use timer ? 
// ==========================================================================
void Ostap::Utils::ProgressConf::setUseTimer ( const bool    value   ) 
{ 
#if defined ( __cplusplus ) && ( 201103L <= __cplusplus ) 
  m_use_timer = value ; 
#else 
  m_use_timer = value && false ;  // always false 
#endif 
}
// ==========================================================================
// width 
// ==========================================================================
void Ostap::Utils::ProgressConf::setWidth    
( const unsigned short value ) { m_width  = value < 512 ? value : 512 ; }
// ==========================================================================
/*  Constructor from the configuration and maximal count
 *  @param maxcount   maximal count (==0 : no display bar)
 *  @param conf       configuration 
 */
// ==========================================================================
Ostap::Utils::ProgressBar::ProgressBar 
( const ProgressConf&      conf     , 
  const unsigned long long maxcount )
  : ProgressConf ( conf     )  
  , m_maxcount   ( maxcount ) 
  , m_start      ( 0        )
    // ==========================================================================
{
  m_wtot = ( what().length () + left().length() + right().length () + 
             ( width () + 2 ) * std::max ( symbol().length() , empty().length() ) + 10 ) ;
  // ==========================================================================
#if defined ( __cplusplus ) && ( 201103L <= __cplusplus ) 
  // ==========================================================================
  if ( use_timer() ) 
  {
    m_start = 
      std::chrono::system_clock::now().time_since_epoch().count() 
      * std::chrono::system_clock::period::num 
      / std::chrono::system_clock::period::den ;
  }
  // ==========================================================================
#endif 
  // ==========================================================================
}
// ==========================================================================
/*  Constructor from configuration and maximal conut 
 *  @param conf       configuration 
 *  @param maxcount   maximal count (==0 : no display bar)
 */
// ==========================================================================
Ostap::Utils::ProgressBar::ProgressBar 
( const unsigned long long maxcount ,  
  const ProgressConf&      conf     ) 
  : ProgressBar ( conf , maxcount ) {}
// ============================================================================
/** full constructor for print to std::cout
 *  @param maxcount   maximal count (==0 : no display bar)
 *  @param symbol     symbol to show 
 *  @param empty      empty symbol
 *  @param left       left/prefix 
 *  @param right      right/prefix 
 *  @param width      the width 
 */
// ============================================================================
Ostap::Utils::ProgressBar::ProgressBar 
( const unsigned long long maxcount  ,
  const unsigned short     width     ,
  const std::string&       symbol    ,
  const std::string&       empty     ,
  const std::string&       left      ,
  const std::string&       right     , 
  const std::string&       what      , 
  const bool               use_timer ,
  const bool               atty      )
  : ProgressBar ( ProgressConf ( width , symbol , empty , left , right , what , use_timer , atty ) , 
                  maxcount )
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Utils::ProgressBar::~ProgressBar()
{
  if ( enabled () ) { show_bar ( false ) ; std::cout << std::endl ; }
}
// ============================================================================
// show the bar
// ============================================================================
Ostap::Utils::ProgressBar& 
Ostap::Utils::ProgressBar::show_bar ( const bool show_eta ) 
{ 
  //
  const unsigned int w = width () ;
  if      ( 0 == w          ) { return *this ; }         // RETURN
  else if ( 0 == m_maxcount ) { return *this ; }         // RETURN
  
  //
  const double       fraction = 1.0 * m_count / m_maxcount ;
  const unsigned int rtics    = w * fraction  ;
  const unsigned int mtics    = std::min ( rtics , w ) ;
  //
  m_next_count = m_maxcount * double ( rtics + 1 ) / w ;
  //
  std::string line =left ()  ;
  line.reserve ( m_wtot ) ;
  //
  // ==========================================================================
#if defined ( __cplusplus ) && ( 201103L <= __cplusplus ) 
  // ==========================================================================  
  //
  const std::string& s1 = symbol() ;
  const std::string& s2 = empty () ;
  //
  const unsigned int  s_LEN = 256 ;
  static char s_buffer [ s_LEN ] ;
  // ==========================================================================
  if ( 3 < mtics && use_timer ()  ) 
  {
    //
    const unsigned long long duration = 
      std::chrono::system_clock::now().time_since_epoch().count() 
      * std::chrono::system_clock::period::num 
      / std::chrono::system_clock::period::den  - m_start ;
    //
    const unsigned long long time = 
      show_eta ? 
      duration * ( std::max ( 1 - fraction , 0. ) / fraction ) : 
      duration ;
    //
    const auto d_ = std::div ( time   , 60 * 60 * 24 ) ;
    const auto h_ = std::div ( d_.rem , 60 * 60      ) ;
    const auto m_ = std::div ( h_.rem , 60           ) ;
    //
    const unsigned short days    = d_ .quot ;
    const unsigned short hours   = h_ .quot ;
    const unsigned short minutes = m_ .quot ;
    const unsigned short seconds = m_ .rem  ;
    //
    const unsigned int used = 
      days && days < 100  ?
                     snprintf ( s_buffer , s_LEN , "%02d:%02d:%02d:%02ds " , days , hours , minutes , seconds ) :
      !days && hours ? 
                     snprintf ( s_buffer , s_LEN , "%02d:%02d:%02ds "             , hours , minutes , seconds ) :
      !days && !hours &&  minutes ? 
                     snprintf ( s_buffer , s_LEN , "%02d:%02ds "                          , minutes , seconds ) :
      !days && !hours && !minutes ? 
                     snprintf ( s_buffer , s_LEN , "%02ds "                                         , seconds ) : 0 ;
    
    if      (  show_eta && used && used + 4 < mtics ) 
    {
      line.append ( "ETA " ) ; 
      line.append ( std::string ( s_buffer , s_buffer + used ) ) ; 
      for ( unsigned short i = used + 4 ; i < mtics   ; ++i ) { line.append ( s1 ) ; }
    } 
    else if ( !show_eta && used && used < mtics ) 
    {
      { line.append ( std::string ( s_buffer , s_buffer + used ) ) ; }
      for ( unsigned short i = used ; i < mtics   ; ++i ) { line.append ( s1 ) ; }
    }
    else 
    { for ( unsigned short i = 0 ; i < mtics   ; ++i ) { line.append ( s1  ) ; } } 
  }
  else 
  { for ( unsigned short i = 0 ; i < mtics   ; ++i ) { line.append ( s1 ) ; } }
  // ==========================================================================
#else 
  // ==========================================================================
  for ( unsigned short i = 0     ; i < mtics   ; ++i ) { line.append ( s1 ) ; }
  // ==========================================================================
#endif 
  // ==========================================================================
  //
  for ( unsigned short i = mtics ; i < w ; ++i ) { line.append ( s2  ) ; }
  //
  line.append ( right () ) ;
  line.append ( " "      ) ;
  line.append ( Ostap::format ( "%4.1f" , 100 * fraction ) ) ;
  line.append ( "%"     ) ;  
  //
  // show it!
  std::cout  << what()  << line << "\r" ;
  if ( atty() ) { std::cout << std::flush ; }
  //
  return *this ;
}

      


// ============================================================================
//                                                                      The END 
// ============================================================================
