// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include  <iostream> 
#include  <mutex> 
#include  <cmath> 
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
}
// ============================================================================
/* helper function to create callable  for drawing of progress bar to the frame 
 *  @param  chunk   sise of each chunk  (parameter <code>everyN</code> for 
 *                                       <code>OnPartialResultSlot</code>)
 *  @param  nchunks total number of chunks 
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
( const unsigned long  chunk   , 
  const unsigned short nchunks , 
  const std::string&   symbol  ,
  const std::string&   blank   ,
  const std::string&   left    , 
  const std::string&   right   )
{
  //
  return [symbol,nchunks,left,right,blank] 
    ( unsigned int islot , ULong64_t& /* u */ ) -> void 
  {
    //
    static unsigned short istep { 0 } ;
    std::lock_guard<std::mutex> lock ( s_mutex_bar ) ;
    //
    istep        += 1      ;
    //
    if ( istep <= nchunks ) 
    {
      const int percent = 
        istep >= nchunks ? 100 : 
        int ( std::floor ( double ( istep * 100 ) / nchunks ) ) ; 
      //
      std::string show = "\r" + left ;
      for  ( unsigned short i = 0     ; i < istep    ; ++i ) { show += symbol ; } 
      for  ( unsigned short j = istep ; j < nchunks  ; ++j ) { show += blank  ; } 
      show += right ;
      std::cout << show << std::right << std::setw  ( 3 ) << percent << "% " ;
    }
    //
    if ( istep == nchunks ) { std::cout << std::endl ; }
    //
    std::cout << std::flush ;
    //
  } ;
  //
}
// ============================================================================
//                                                                      The END 
// ============================================================================
