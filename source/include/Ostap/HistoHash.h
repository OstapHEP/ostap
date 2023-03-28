#ifndef OSTAP_HISTOHASH_H 
#define OSTAP_HISTOHASH_H 1
// ============================================================================
// Include files
// ============================================================================
// Forward declarations 
// ============================================================================
class TH1    ; // ROOT 
class TAxis  ; // ROOT 
class TGraph ; // ROOT
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    class Histo1D ;
    class Histo2D ;
    class Histo3D ;
    // ========================================================================
  } 
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
    /** get hash for Ostap::Math::Histo1D  object 
     *  @see Ostap::Math::Histo1D 
     */
    std::size_t hash_histo ( const Ostap::Math::Histo1D& h ) ;
    // ========================================================================    
    /** get hash for Ostap::Math::Histo2D  object 
     *  @see Ostap::Math::Histo2D 
     */
    std::size_t hash_histo ( const Ostap::Math::Histo2D& h ) ;
    // ========================================================================
    /** get hash for Ostap::Math::Histo3D  object 
     *  @see Ostap::Math::Histo3D 
     */
    std::size_t hash_histo ( const Ostap::Math::Histo3D& h ) ;
    // ========================================================================    
  } //                                        The end of namepsace Ostap::Utils 
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_HISTOHASH_H
// ============================================================================
