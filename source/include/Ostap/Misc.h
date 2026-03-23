// ============================================================================
#ifndef OSTAP_MISC_H 
#define OSTAP_MISC_H 1
// ============================================================================
/** @file Ostap/Misc.h
 *  collection of various C++ utilities 
 *  @author Vanya Belyaev
 *  @date   2018-03-23
 */
// ============================================================================
// forward declarations
// ============================================================================
class TVirtualPad ; // ROOT
class TCanvas     ; // ROOT 
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Utils 
  {
    // ========================================================================
    /** get the (current) pad
     *  @see TVirtualPad::Pad 
     */
    TVirtualPad* get_pad    () ;
    // ========================================================================
    /** get the (current) canvas
     *  - current canvas is a canvas associated with current pad 
     *  @see TCanvas 
     */
    TCanvas*     get_canvas () ;
    // ========================================================================
    /** get existing canvas by name 
     *  @param name   Canvas name
     *  @return canvas or nullptr 
     */
    TCanvas*     get_canvas
    ( const std::string& name ) ; 
    // ========================================================================
    /** check existing canvas by name 
     *  @param name   Canvas name
     *  @return canvas or nullptr 
     */
    bool         has_canvas
    ( const std::string& name )  ;
    // ========================================================================
    /** call for TVirtualPad::Update 
     *  @see TVirtualPad::Update 
     */
    TVirtualPad* pad_update
    ( TVirtualPad* pad = nullptr ) ;
    // ========================================================================
    /** call for TVirtualPad::UpdateAsync 
     *  @see TVirtualPad::Update
     *  @see TVirtualPad::UpdateAsync 
     */
    TVirtualPad* pad_update_async
    ( TVirtualPad* pad = nullptr ) ;
    // ========================================================================
    /** @class CanvasContext 
     *  @see TVirtualPad::TContext 
     */
    class CanvasContext
    {
    public:
      // ======================================================================
      /// default constructor, catch the current canvas 
      CanvasContext  () ;
      /// destructor
      ~CanvasContext () ;
      // =====================================================================
    public:
      // =====================================================================
      /// get the saved canvas 
      inline TCanvas* saved  () const { return m_saved  ; }
      // =====================================================================
    public:
      // =====================================================================
      /// CONTEXT MANAGER: (re-catch the curent canvas)
      bool        enter  () ;
      /// CONTEXT MANAGER: exit  
      bool        exit   () ; 
      /// Is the context active (keep some vaild canvas)
      inline bool active () const { return m_saved ; } 
      /// =====================================================================
    private:
      // ======================================================================
      /// saved canvas 
      TCanvas* m_saved       { nullptr } ; // saved canvas 
      // ======================================================================
    } ;
    // ========================================================================
  } //                                        The end of namespace Ostap::Utils 
  // ==========================================================================
} //                                                The  end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_MISC_H
// ============================================================================
