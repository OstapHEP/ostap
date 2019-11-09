// ============================================================================
#ifndef OSTAP_DATAFRAMEUTILS_H 
#define OSTAP_DATAFRAMEUTILS_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <string>
#include <functional>
// ============================================================================
// Ostap 
// ============================================================================
#include "Ostap/DataFrame.h"
// ============================================================================
namespace Ostap 
{
  // ==========================================================================
  namespace Utils 
  {
    // =======================================================================
    /** helper function to create callable  for drawing of progress bar to the frame 
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
    std::function<void(unsigned int,ULong64_t&)>
    frame_progress
    ( const unsigned long  chunk          , 
      const unsigned short nchunks        , 
      const std::string&   symbol  = "#"  , 
      const std::string&   blank   = " "  ,
      const std::string&   left    = "[ " ,
      const std::string&   right   = " ]" ) ;
    // ========================================================================
  } //                                        The end of namespace Ostap::Utils 
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
#endif // OSTAP_DATAFRAMEUTILS_H
// ============================================================================
//                                                                      The END 
// ============================================================================
