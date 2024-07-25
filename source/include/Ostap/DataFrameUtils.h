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
    /// forward declaration 
    class ProgressConf ; // forward declaration 
    // =======================================================================
    /** helper function to create callable  for drawing of progress bar to the frame 
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
    std::function<void(unsigned int,ULong64_t&)>
    frame_progress
    ( const unsigned short nchunks        , 
      const unsigned short width   =  80  , 
      const std::string&   symbol  = "#"  , 
      const std::string&   blank   = " "  ,
      const std::string&   left    = "[ " ,
      const std::string&   right   = " ]" ) ;
    // ========================================================================
    /** helper function to create callable  for drawing of progress bar to the frame 
     *  @param  nchunks  total number of chunks 
     *  @param  progress progress bar  configuration 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2019-11-08
     */
    std::function<void(unsigned int,ULong64_t&)>
    frame_progress
    ( const unsigned short nchunks        , 
      const ProgressConf&  progress       ) ;
    // ========================================================================
    /** get pool size 	      
     *  @see ROOT::IsImplicitMTEnabled   () 
     *  @see ROOT::GetThreadPoolSize     () 
     *  @see ROOT::GetImplicitMTPoolSize () 
     */
    unsigned int mt_pool_size () ;
    // ========================================================================
  } //                                        The end of namespace Ostap::Utils 
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
#endif // OSTAP_DATAFRAMEUTILS_H
// ============================================================================
//                                                                      The END 
// ============================================================================
