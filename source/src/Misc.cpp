// $Id$
// ============================================================================
// Include files 
// ============================================================================
// ROOT 
// ============================================================================
#include "RVersion.h"
#include "TList.h"
#include "TVirtualPad.h"
#include "TCanvas.h"
#include "TROOT.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Misc.h"
// ============================================================================
/** get the (current) pad
 *  @see TVirtualPad::Pad 
 */
// ============================================================================
TVirtualPad* Ostap::Utils::get_pad ()
{
  // TVirtualPad* pad  = nullptr ; 
  // const TROOT* groot = ROOT::GetROOT() ;
  // (1) current` pad?
  // if ( nullptr != groot ) { pad = groot->GetSelectedPad() ; }
  // (2) global  pad? 
  // if ( nullptr != pad   ) { return pad ; }
  return TVirtualPad::Pad() ;
}
// ============================================================================
/* get the (current) canvas 
 *  @see TCanvas 
 */
// ============================================================================
TCanvas*  Ostap::Utils::get_canvas ()
{
  TVirtualPad* pad = get_pad() ;
  if ( nullptr == pad ) { return nullptr ; }
  return pad->GetCanvas() ; 
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
  if ( nullptr != pad && pad->IsModified() ) { pad -> Update () ; }
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
  if ( nullptr != pad && pad -> IsModified() )
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
// default constructor
// ============================================================================
Ostap::Utils::CanvasContext::CanvasContext ()
  : m_saved ( nullptr )
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Utils::CanvasContext::~CanvasContext()
{
  if ( m_saved ) { exit () ; }
  m_saved = nullptr ;
}
// ============================================================================
// CONTEXT MANAGER: (re)catch the current canvas 
// ============================================================================
bool Ostap::Utils::CanvasContext::enter ()
{
  m_saved = get_canvas () ;
  if ( m_saved  && m_saved -> IsModified() ) { m_saved->Update() ; } 
  return active () ;
} 
// ============================================================================
// CONTEXT MANAGER: exit  
// ============================================================================
bool Ostap::Utils::CanvasContext::exit  ()
{
  // nothing to do..
  if ( !m_saved ) { return active () ; }
  //
  // check the saved canvas is existing valid TCanvas and make at active
  // otherwise switch to the current_canvas if it is valid 
  //  
  TROOT* groot = ROOT::GetROOT ()  ;
  if ( groot )
  {
    /// get the list of all alived canvases
    TSeqCollection* all_canvases = groot->GetListOfCanvases() ;
    if ( all_canvases )
    {
      for ( const TObject* obj : *all_canvases )
      {
	if ( obj && obj == m_saved )
	{
	  m_saved -> cd ()   ;
	  if ( m_saved -> IsModified() ) { m_saved -> Update() ; } 
	  m_saved =  nullptr ;
	  return active ()   ;	               // RETURN
	  // ==================================================================
	} // valid canvas is found in  the list 
      }	// list of all canvas 
    } // valid collection of all canvases 
  } // valid TROOT 
  /// 
  /// otherwise switch to the current canvas if valid 
  m_saved      = nullptr ;
  TCanvas* cnv = get_canvas () ;
  if ( cnv )
  {
    cnv -> cd () ;
    if ( cnv -> IsModified() ) { cnv -> Update () ; }
  }
  //
  return active ();
}
// ============================================================================
//                                                                      The END 
// ============================================================================
