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
  class DataFrameProgress
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
    ( const unsigned short nchunks        , 
      const unsigned short width          , 
      const std::string&   symbol  = "#"  , 
      const std::string&   blank   = " "  ,
      const std::string&   left    = "[ " ,
      const std::string&   right   = " ]" ) 
      : m_nchunks ( nchunks ) 
      , m_width   ( width < 10u  ? 10u : width ) 
      , m_symbol  ( symbol  ) 
      , m_blank   ( blank   ) 
      , m_left    ( left    ) 
      , m_right   ( right   ) 
      , m_chunks  ( 0       )
    {}
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
      const double done = 
        m_nchunks <= m_chunks ? 100  : double ( m_chunks * 100 ) / m_nchunks  ;
      //
      const unsigned short ns =
        m_nchunks <= m_chunks ? m_width :
        int ( std::floor ( 0.01 * done * m_width ) ) ;  
      //
      if ( m_chunks <= m_nchunks ) 
      {
        std::string bar      ;
        bar.reserve ( 1024 ) ;
        bar += m_left ;
        for ( unsigned short i = 0  ; i < ns       ; ++i ) { bar += m_symbol ; }
        for ( unsigned short i = ns ; i < m_width  ; ++i ) { bar += m_blank  ; }
        bar += m_right ;
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
    /// bar width 
      unsigned short m_width           ; // bar width 
    /// symbol 
    std::string    m_symbol { "#"  } ; // the symbol 
    /// blank symbol
    std::string    m_blank  { " "  } ; // blank symbol
    // left edge 
    std::string    m_left   { "[ " } ; // left edge 
    // right edge 
    std::string    m_right  { " ]" } ; // right edge 
    // ========================================================================
  private :
    // ========================================================================
    /// number of processed chunks 
    unsigned int m_chunks { 0 }   ; // number of currently processed chunks 
    // ========================================================================
  } ; 
  // ==========================================================================
} //                                     The end of local anonymous namespace 
// ============================================================================
/* helper function to create callable  for drawing of progress bar to the frame 
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
  return DataFrameProgress ( nchunks , 
                             width   , 
                             symbol  , 
                             blank   , 
                             left    , 
                             right   ) ;
}


// ============================================================================
//                                                                      The END 
// ============================================================================
