// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include  <iostream> 
#include  <mutex> 
#include  <cmath> 
// ============================================================================
// ROOT 
// ============================================================================
#include "RtypesCore.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/DataFrame.h"
#include "Ostap/DataFrameUtils.h"
#include "Ostap/ProgressBar.h"
// ============================================================================
// local 
// ============================================================================
#include "OstapDataFrame.h"
// ============================================================================
/** @file 
 *  Implementation file for DataFrame related functions from  
 *  namespace Ostap::Utils
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2019-11-08
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  /// local mutex from progress  bar 
  std::mutex s_mutex_bar {} ;
  // ==========================================================================
  /** @class DataFrameProgress 
   *  Helper class to show th eprogree bar for data frame 
   */
  class DataFrameProgress : public Ostap::Utils::ProgressConf 
  {
  public:
    // ========================================================================
    /** constructor with all parameters 
     *  @param  nchunks total number of chunks 
     *  @param  width   bar width 
     *  @param  symbol  symbol to use as progress 
     *  @param  blank   blank symbol
     *  @param  left    prefix 
     *  @param  right   suffix 
     */
    DataFrameProgress
    ( const unsigned short               nchunks  , 
      const Ostap::Utils::ProgressConf & progress )
      : ProgressConf ( progress ) 
      , m_nchunks    ( nchunks ) 
      , m_chunks     ( 0       )
    {
      if ( !m_nchunks ) { setWidth ( 0 ) ; }   // DISBALE IT!
    }
    // ========================================================================
    /// default move constructor
    DataFrameProgress (       DataFrameProgress&& ) = default ;
    /// disabled copy constructir 
    DataFrameProgress ( const DataFrameProgress&  ) = default ;
    // ========================================================================
  public:
    // ========================================================================
    // the main method 
    // ========================================================================
    void operator() ( unsigned int islot ,  ULong64_t& /* u */ ) 
    {
      std::lock_guard<std::mutex> lock ( s_mutex_bar ) ;
      ///
      /// increment number of processed chunks 
      ++m_chunks ; // increment number of processed chunks
      ///
      const unsigned int w  = width () ;
      if ( !w || !m_nchunks ) { return ; } // DISABLED! 
      //
      const double done     = m_nchunks <= m_chunks ? 100 : double ( m_chunks * 100 ) / m_nchunks  ;
      const unsigned int ns = m_nchunks <= m_chunks ? w   : int ( std::floor ( 0.01 * done * w ) ) ;  
      //
      const std::string& s1 = symbol () ;
      const std::string& s2 = empty  () ;
      //
      if ( m_chunks <= m_nchunks ) 
      {
        std::string bar      ;
        bar.reserve ( 1024 ) ;
        bar += left () ;
        for ( unsigned int i = 0  ; i < ns ; ++i ) { bar += s1 ; }
        for ( unsigned int i = ns ; i < w  ; ++i ) { bar += s2 ; }
        bar += right ();
        //
        std::cout << bar << ' ' << std::floor ( done ) <<  "%" ;
        if ( m_chunks < m_nchunks ) { std::cout << '\r'                    ; }
        else                        { std::cout << std::endl << std::flush ; } 
      }
    }
    // ========================================================================
  private :
    // ========================================================================
    /// number of currrently processes chunks       
    unsigned int   m_nchunks         ; // number of chunks 
    // ========================================================================
    /// number of processed chunks 
    unsigned int m_chunks { 0 }   ; // number of currently processed chunks 
    // ========================================================================
  } ; 
  // ==========================================================================
} //                                     The end of local anonymous namespace 
// ============================================================================
/*  helper function to create callable  for drawing of progress bar to the frame 
 *  @param  nchunks total number of chunks 
 *  @param  width   effective bar width (no left, rigtht & percentage) 
 *  @param  symbol symbol to use as progress 
 *  @param  blank  blank symbol
 *  @param  left   prefix 
 *  @param  right  suffix 
 *  The format of progress bar is:
 *  <code>left+(%*psymbol)+(N-%)*blank+right+percentage</code>
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2019-11-08
 */
// ============================================================================
std::function<void(unsigned int,ULong64_t&)>
Ostap::Utils::frame_progress
( const unsigned short nchunks , 
  const unsigned short width   , 
  const std::string&   symbol  ,
  const std::string&   blank   ,
  const std::string&   left    , 
  const std::string&   right   )
{
  //
  return DataFrameProgress 
    ( nchunks ,
      Ostap::Utils::ProgressConf ( width   , 
                                   symbol  , 
                                   blank   , 
                                   left    , 
                                   right   ) ) ; 
}
// ===========================================================================
/*  helper function to create callable  for drawing of progress bar to the frame 
 *  @param  nchunks  total number of chunks 
 *  @param  progress progress bar  configuration 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date   2019-11-08
 */
// ============================================================================
std::function<void(unsigned int,ULong64_t&)>
Ostap::Utils::frame_progress
( const unsigned short               nchunks  , 
  const Ostap::Utils::ProgressConf&  progress ) 
{ return DataFrameProgress ( nchunks , progress ) ; }
// ============================================================================

// ============================================================================
//                                                                      The END 
// ============================================================================
