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
     *  @see TROOT::GetSelectedPad 
     *  @see TVirtualPad::Pad 
     */
    TVirtualPad* get_pad    () ;
    // ========================================================================
    /** get the (current) canvas 
     *  @see TCanvas 
     */
    TCanvas*     get_canvas () ;
    // ========================================================================
    /** call for TVirtualPad::Update 
     *  @see TVirtualPad::Update 
     */
    TVirtualPad* pad_update       ( TVirtualPad* pad = nullptr ) ;
    // ========================================================================
    /** call for TVirtualPad::UpdateAsync 
     *  @see TVirtualPad::Update
     *  @see TVirtualPad::UpdateAsync 
     */
    TVirtualPad* pad_update_async ( TVirtualPad* pad = nullptr ) ;
    // ========================================================================
  } //                                        The end of namespace Ostap::Utils 
  // ==========================================================================
} //                                                The  end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_MISC_H
// ============================================================================
