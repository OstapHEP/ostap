// ============================================================================
// Include files 
// ============================================================================
// ROOT 
// ============================================================================
#include "TCut.h"
#include "TTree.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Formula.h"
#include "Ostap/PySelectorWithCuts.h"
// ============================================================================
// local
// ============================================================================
#include "Exception.h"
// ============================================================================
/** @file 
 *  Implementation file for class Analysis::SelectorWithCuts
 * 
 *  @date 2013-05-06 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  const std::string s_whitespaces = " \t\n\r\f\v" ;
  inline std::string strip ( const std::string& s                          ,
                             const std::string& whitespace = s_whitespaces )
  {
    const std::string::size_type p1 = s .find_first_not_of ( whitespace );
    if ( std::string::npos == p1 ) { return ""; }
    //
    const std::string::size_type p2 = s .find_last_not_of ( whitespace ) ; 
    const auto range = p2 + 1 - p1 ;
    //
    return s .substr ( p1 , range );
  }
  // ==========================================================================
}
// ============================================================================
// ClassImp(Ostap::SelectorWithCuts) ;
// ============================================================================
Ostap::SelectorWithCuts::SelectorWithCuts
(
#if defined(OSTAP_OLD_PYROOT) && OSTAP_OLD_PYROOT
  PyObject*   self          , 
#endif  
  const TCut& cuts          ,
  TTree*      tree          ) :
#if defined(OSTAP_OLD_PYROOT) && OSTAP_OLD_PYROOT
  SelectorWithCuts ( self , std::string ( cuts.GetTitle () ) , tree ) 
#else
  SelectorWithCuts (        std::string ( cuts.GetTitle () ) , tree ) 
#endif 
{}
// ============================================================================
// constructor 
// ============================================================================
Ostap::SelectorWithCuts::SelectorWithCuts
( 
#if defined(OSTAP_OLD_PYROOT) && OSTAP_OLD_PYROOT
  PyObject*          self , 
#endif
  const std::string& cuts ,     
  TTree*             tree )
#if defined(OSTAP_OLD_PYROOT) && OSTAP_OLD_PYROOT
  : Ostap::Selector ( self , tree ) 
#else
  : Ostap::Selector (        tree ) 
#endif  
  , fMyCuts            ( strip ( cuts ) ) 
  , fMyFormula         ()
  , m_good             ( 0              )            
{
  if ( tree && !fMyCuts.empty() )
  {
    Ostap::Assert ( make_formula ( tree )     , 
                    "Invalid formula "        , 
                    "Ostap::SelectorWithCuts" ) ;
  }
}
// ============================================================================
// virtual destructor 
// ============================================================================
Ostap::SelectorWithCuts::~SelectorWithCuts () {}
// ============================================================================
bool Ostap::SelectorWithCuts::make_formula ( TTree*   tree ) 
{
  if ( fMyCuts.empty() ) { return true ; }
  fMyFormula.reset() ;
  if ( nullptr == tree ) { return false ; }
  //
  fMyFormula.reset ( new Ostap::Formula ( fMyCuts , tree ) ) ; 
  return ok  () ;    
}
// ============================================================================
// reset formula unisy new tree 
// ============================================================================
void Ostap::SelectorWithCuts::reset_formula ( TTree* tree ) 
{
  //
  set_tree ( tree ) ;
  fMyFormula.reset() ;
  //
  if ( tree && !fMyCuts.empty() ) 
  { Ostap::Assert ( make_formula ( tree )           , 
                    "Invalid formula "              , 
                    "Ostap::SelectorWithCuts::Init" ) ; }
}   
// ============================================================================
// notify 
// ============================================================================
Bool_t Ostap::SelectorWithCuts::Notify() 
{
  if ( fMyFormula ) { fMyFormula->Notify() ; }
  return Ostap::Selector::Notify () ;
}
// ============================================================================
// init 
// ============================================================================
void Ostap::SelectorWithCuts::Init ( TTree* tree ) 
{
  reset_formula  ( tree  ) ;
  Ostap::Selector::Init ( tree ) ;
}
// ============================================================================
//  begin 
// ============================================================================
void Ostap::SelectorWithCuts::Begin ( TTree* tree ) 
{
  //
  reset_formula  ( tree  ) ;
  Ostap::Selector::Begin ( tree ) ; 
}
// ============================================================================
// slave begin 
// ============================================================================
void Ostap::SelectorWithCuts::SlaveBegin ( TTree* tree ) 
{
  //
  Ostap::Assert ( tree                  , 
                  "Invalid tree"        , 
                  "Ostap::SelectorWithCuts::SlaveBegin" ) ;
  reset_formula ( tree  ) ;
  //
  Ostap::Selector::SlaveBegin  ( tree ) ;
}
// ============================================================================
// check the entry 
// ============================================================================
bool Ostap::SelectorWithCuts::good_entry  ( Long64_t entry ) 
{
  //
  if ( Ostap::Selector::GetEntry ( entry ) < 0 ) 
  {
    Ostap::Selector::Abort ( "" , TSelector::kAbortFile ) ;
    return false ; 
  }
  //
  if ( !fMyCuts.empty() && fMyFormula && fMyFormula->GetNdim() && !fMyFormula ->evaluate() )
  { return false ; }
  //
  return true ;  
}
// ============================================================================
// process 
// ============================================================================
Bool_t Ostap::SelectorWithCuts::Process      ( Long64_t entry ) 
{ 
  //
  if ( !good_entry ( entry ) ) { return false ; }
  //
  // increment the event counter for good events 
  ++m_good  ;
  //
#if ROOT_VERSION_CODE < ROOT_VERSION(6,22,0)
  //
  return Ostap::Selector::Process ( entry ) ;
  //
#else
  // ==========================================================================
  // total event counter 
  increment_event () ;
  return process_entry () ;
  // ==========================================================================
#endif 
}
// ============================================================================
// is formula OK?
// ============================================================================
bool Ostap::SelectorWithCuts::ok () const // is formula OK ? 
{ return fMyCuts.empty() || ( fMyFormula && fMyFormula->ok () ) ; }
// ============================================================================
// process good entry 
// ============================================================================
bool Ostap::SelectorWithCuts::process_entry () { return true ; }
// ============================================================================
// new cuts 
// ============================================================================
void Ostap::SelectorWithCuts::set_cuts  ( const TCut& cuts ) 
{ set_cuts ( std::string ( cuts.GetTitle() ) ) ; }
// ============================================================================
// new cuts 
// ============================================================================
void Ostap::SelectorWithCuts::set_cuts  ( const std::string& cuts ) 
{
  fMyFormula.reset() ;
  fMyCuts = strip ( cuts ) ;
}
// ============================================================================



// ============================================================================
//                                                                     The END 
// ============================================================================
