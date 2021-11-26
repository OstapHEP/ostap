#ifndef OSTAP_HISTOHASH_H 
#define OSTAP_HISTOHASH_H 1
// ============================================================================
// Include files
// ============================================================================
// Forward declarations 
// ============================================================================
class TH1    ; // ROOT 
class TGraph ; // ROOT
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Utils
  {
    // ========================================================================
    /** get hash for the given graph 
     *  @param graph the graphs 
     *  @return hash value 
     */
    std::size_t hash_graph ( const TGraph* graph  ) ;
    // ========================================================================
    /** get hash for the given histogram 
     *  @param histo the historgam  
     *  @return hash value 
     */
    std::size_t hash_histo ( const TH1* hist      ) ;
    // ========================================================================
    /*  get hash for the given axis 
     *  @param axis the axis 
     *  @return hash value 
     */
    std::size_t hash_axis  ( const TAxis* axis )  ;
    // ========================================================================    
  } //                                        The end of namepsace Ostap::Utils 
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_HISTOHASH_H
// ============================================================================
