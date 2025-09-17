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
// default contructor
// ============================================================================
Ostap::Utils::PadContext::PadContext ()
  : PadContext ( ROOT::GetROOT() && !ROOT::GetROOT()->IsBatch() ? true : false )
{}
// ============================================================================
// contructor with interactive flag 
// ============================================================================
Ostap::Utils::PadContext::PadContext
( const bool interactive )
  : m_active      ( true        )
  , m_interactive ( interactive )
  , m_saved       ( TVirtualPad::Pad() ) 
{}
// ============================================================================
// full contructor
// ============================================================================
Ostap::Utils::PadContext::PadContext
( TVirtualPad* pad         ,
  const bool   interactive ,
  const bool   not_null    )
  : m_active      ( true        )
  , m_interactive ( interactive )
  , m_saved       ( TVirtualPad::Pad() ) 
{
  if       ( pad && m_interactive ) { pad->cd() ; }
  else if  ( pad || !not_null     ) { TVirtualPad::Pad() = pad ; }
}
// ============================================================================
// destructor
// ============================================================================
Ostap::Utils::PadContext::~PadContext()
{
  exit () ;
  m_saved = nullptr ;
}
// ============================================================================
// CONTEXT MANAGER: (fake) enter
// ============================================================================
const Ostap::Utils::PadContext&
Ostap::Utils::PadContext::enter () { return *this ; } 
// ============================================================================
// CONTEXT MANAGER: exit  
// ============================================================================
const Ostap::Utils::PadContext&
Ostap::Utils::PadContext::exit  ()
{
  /// avoid multiple exit
  if ( !m_active ) { return *this ; }
  //
  if      ( !m_interactive || !m_saved ) { TVirtualPad::Pad() = m_saved ; }
  else if ( m_saved ) { m_saved->cd () ; } 
  //
  m_active = false ;
  //
  return *this ;
}



// ============================================================================
//                                                                      The END 
// ============================================================================
