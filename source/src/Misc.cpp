// $Id$
// ============================================================================
// Include files 
// ============================================================================
// ROOT 
// ============================================================================
#include "RVersion.h"
#include "TVirtualPad.h"
#include "TROOT.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Misc.h"
// ============================================================================
/** get the (current) pad
 *  @see TROOT::GetSelectedPad 
 *  @see TVirtualPad::Pad 
 */
// ============================================================================
TVirtualPad* Ostap::Utils::get_pad ()
{
  TVirtualPad* pad = nullptr ; 
  const TROOT* groot = ROOT::GetROOT() ;
  // (1) current` pad?
  if ( nullptr != groot ) { pad = groot->GetSelectedPad() ; }
  // (2) global  pad? 
  if ( nullptr != pad   ) { return pad ; }
  return TVirtualPad::Pad() ;
}
// ============================================================================
/*  call for TVirtualPad::Update 
 *  @see TVirtualPad::Update 
 */
// ============================================================================
TVirtualPad* Ostap::Utils::pad_update  ( TVirtualPad* pad )
{
  //
  if ( nullptr == pad ) { pad = get_pad () ; }
  // call for Update 
  if ( nullptr != pad ) { pad -> Update () ; }
  //
  return pad ;
}
// ============================================================================
/*  call for TVirtualPad::UpdateAsync 
 *  @see TVirtualPad::UpdateAsync 
 */
// ============================================================================
TVirtualPad* Ostap::Utils::pad_update_async  ( TVirtualPad* pad )
{
  if ( nullptr == pad ) { pad = get_pad () ; }
  // call for Update 
  if ( nullptr != pad )
    {
#if ROOT_VERSION(6,30,0)<=ROOT_VERSION_CODE
      pad -> UpdateAsync () ;
#else
      pad -> Update      () ;
#endif      
    }
  //
  return pad ;
}
// ============================================================================


// ============================================================================
//                                                                      The END 
// ============================================================================
