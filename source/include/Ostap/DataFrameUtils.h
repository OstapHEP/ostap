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
// ROOT 
// ============================================================================
#include "RVersion.h"
// ============================================================================
// Ostap 
// ============================================================================
#include "Ostap/DataFrame.h"
// ============================================================================
#if ROOT_VERSION(6,32,0) <= ROOT_VERSION_CODE
// ============================================================================
namespace ROOT { namespace RDF { template <class T> class RResultPtr  ; } }
// ============================================================================
#endif
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
#if ROOT_VERSION(6,32,0) <= ROOT_VERSION_CODE
    // ========================================================================
    /** helper utility to add progress bar
     *  to "Count"-action 
     *  @see https://root-forum.cern.ch/t/problems-with-onpartialresultslot-in-new-root-version-6-32-02/60257/3
     *  @param count    input counter 
     *  @param nchunks 
     *  @param progress progress bar configuration      
     */
    ROOT::RDF::RResultPtr<ULong64_t>&
    add_progress_bar
    ( ROOT::RDF::RResultPtr<ULong64_t>& result      ,
      const unsigned short             nchunks     ,
      const unsigned long              how_often   , 
      const ProgressConf&              progress    ) ;  
    // ========================================================================
#endif
    // ========================================================================
  } //                                        The end of namespace Ostap::Utils 
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
#endif // OSTAP_DATAFRAMEUTILS_H
// ============================================================================
//                                                                      The END 
// ============================================================================
