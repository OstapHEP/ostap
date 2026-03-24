// ============================================================================
// Include files 
// ============================================================================
// STD&STL
// ============================================================================
#include <string>
#include <limits>
// ============================================================================
// ROOT 
// ============================================================================
#include "TROOT.h"
#include "TObject.h"
#include "TNamed.h"
// ============================================================================
// ROOT/roofit  
// ============================================================================
#include "RooNameReg.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/RootID.h"
#include "Ostap/Names.h"
// ============================================================================
// local 
// ============================================================================
#include "format.h"
// ============================================================================
/** @file 
 *  Implementation file for function Ostap::Utils::rootID
 *  @date 2020-09-02 
 *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
 */
// ============================================================================
/** @fn rootID 
 *  Find "not-used" name in ROOT 
 *  @param prefix use this prefix 
 *  @param suffix use this suffix
 *  @see TROOT::FindObject 
 *  @see RooNameReg  
 *  @see TROOT::FindObject 
 *  @author Vanya Belyaev
 *  @date   2020-09-02
 */
// ============================================================================
std::string Ostap::Utils::rootID
( const std::string& prefix ,
  const std::string& suffix ) 
{
  //
  if ( !suffix.empty() || !prefix.empty() )
  {
    const std::string name { prefix + suffix } ;
    if ( !usedRootID ( name ) ) { return name ; } 
  }
  //
  if ( prefix.empty() ) { return rootID ( "root_" , suffix ) ; }
  //
  TROOT* root = ROOT::GetROOT() ;
  if ( !root ) { return prefix + "0000" + suffix ; }
  // 
  /// lock it...
  R__LOCKGUARD(gROOTMutex) ;    
  //
  static const unsigned long N = std::numeric_limits<unsigned long>::max ()  ; 
  for ( unsigned long label = 1001 ; label < N ; ++label )
  {
    const std::string tag { prefix +
			    Ostap::format ( ( label <   10000 ) ? "%04u" :
					    ( label <  100000 ) ? "%05u" :
					    ( label < 1000000 ) ? "%06u" : "%u" , label ) + suffix } ;      
    /// check ROOT & ROOT/RooFit 
    if      ( nullptr != RooNameReg::known ( tag.c_str () ) ) { continue ; }
    else if ( nullptr != root->FindObject  ( tag.c_str () ) ) { continue ; }    
    /// 
    return tag ;
  }
  //
  return prefix + "XXXX" + suffix ;
}
// ============================================================================
/*  @fn usedRootID 
 *  Is this name already used by ROOT/RooFit ? 
 *  @see TROOT::FindObject 
 *  @see RooNameReg 
 */
// ============================================================================
#include <iostream> 
bool Ostap::Utils::usedRootID
( const std::string& name )
{
  std::cerr  << "usedRootID/1  " << name << std::endl ;
  TROOT* root = ROOT::GetROOT() ;
  std::cerr  << "usedRootID/2  " << name << std::endl ;
  if ( !root ) { return false ; }
  //
  std::cerr  << "usedRootID/3  " << name << std::endl ;

  //
  if ( nullptr != RooNameReg::known ( name.c_str() ) ) { return true ; }

   std::cerr  << "usedRootID/4  " << name << std::endl ;
  /// lock it...
  R__LOCKGUARD(gROOTMutex) ;

  std::cerr  << "usedRootID/5  " << name << std::endl ;
  
  return  ( nullptr != root->FindObject  ( name.c_str() ) ) ;
  //
}
// ============================================================================
//                                                                     The END 
// ============================================================================

