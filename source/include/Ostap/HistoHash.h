// ============================================================================
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
    /** same binning ?
     *  @param a the first axis 
     *  @param b the second axis 
     *  @return True if axes have the same binnings , False otherwise
     */
    bool same_binning
      ( const TAxis& a ,
	const TAxis& b ) ;
    // ========================================================================
    /** same binning ?
     *  @param  a the first histo 
     *  @param  b the second histo 
     *  @return True if histos have the same binnings , False otherwise
     */
    bool same_binning
      ( const TH1&    a ,
	const TH1&    b ) ;    
    // ========================================================================
    /** check if bins are uniform 
     *  @param axis axis to be checked 
     *  @param check actually check via loop over bin edges)
     */
    bool uniform_bins 
      ( const TAxis* axis         , 
	const bool   check = true ) ;
    // ========================================================================    
    /** check if bins of histogram are uniform 
     *  @param histo historgam to be checked 
     *  @param check actually check via loop over bin edges)
     */
    bool uniform_bins 
      ( const TH1& histo          , 
	const bool check = true ) ;
    // ========================================================================    
  } //                                        The end of namepsace Ostap::Utils 
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_HISTOHASH_H
// ============================================================================
