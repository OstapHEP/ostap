// ============================================================================
// Include files 
// ============================================================================
// ROOT 
// ============================================================================
#include "RVersion.h"
#include "TCut.h"
#include "TTree.h"
#include "TChain.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Formula.h"
#include "Ostap/StatusCode.h"
#include "Ostap/PySelectorWithCuts.h"
#include "Ostap/ProgressConf.h"
#include "Ostap/ProgressBar.h"
// ============================================================================
// local
// ============================================================================
#include "status_codes.h"
// ============================================================================
/** @file 
 *  Implementation file for class SelectorWithCuts
 *  @date 2013-05-06 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  const std::string s_whitespaces = " \t\n\r\f\v" ;
  inline std::string strip
  ( const std::string& s                          ,
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
// Full constructor 
// ============================================================================
Ostap::SelectorWithCuts::SelectorWithCuts
(  const std::string&                cuts     ,     
   TTree*                            tree     ,
   const Ostap::Utils::ProgressConf& progress )   
  : Ostap::Selector    ( tree , progress ) 
  , m_cuts             ( strip ( cuts )  ) 
  , m_formula          ()
  , m_good             ( 0              )            
{
 Ostap:Assert ( !get_tree() || m_cuts.empty() || make_formula ( tree ) , 
                "Invalid formula"        ,
                "Ostap::SelectorWithCuts" ,
                INVALID_FORMULA , __FILE__ , __LINE__ ) ;
}
// ============================================================================
// Full constructor 
// ============================================================================
Ostap::SelectorWithCuts::SelectorWithCuts
(  const TCut&                       cuts     ,     
   TTree*                            tree     ,
   const Ostap::Utils::ProgressConf& progress )
  : SelectorWithCuts ( std::string ( cuts.GetTitle() ) , tree , progress )
{}
// ============================================================================
// Full constructor 
// ============================================================================
Ostap::SelectorWithCuts::SelectorWithCuts
(  const std::string&                cuts  ,     
   TTree*                            tree  )
  : SelectorWithCuts ( cuts , tree , false  )
{}
// ============================================================================
// Full constructor 
// ============================================================================
Ostap::SelectorWithCuts::SelectorWithCuts
(  const TCut&                       cuts     ,     
   TTree*                            tree     ) 
  : SelectorWithCuts ( std::string ( cuts.GetTitle() ) , tree )
{}
// ============================================================================
// Full constructor 
// ============================================================================
Ostap::SelectorWithCuts::SelectorWithCuts
(  const std::string&                cuts     ,     
   const Ostap::Utils::ProgressConf& progress , 
   TTree*                            tree     ) 
  : SelectorWithCuts ( cuts , tree , progress )
{}
// ============================================================================
// Full constructor 
// ============================================================================
Ostap::SelectorWithCuts::SelectorWithCuts
(  const TCut&                       cuts     ,     
   const Ostap::Utils::ProgressConf& progress , 
   TTree*                            tree     ) 
  : SelectorWithCuts ( cuts , tree , progress )
{}
// ============================================================================
// virtual destructor 
// ============================================================================
Ostap::SelectorWithCuts::~SelectorWithCuts () {}
// ============================================================================
bool Ostap::SelectorWithCuts::make_formula ( TTree*   tree ) 
{
  if ( m_cuts.empty()  ) { return true ; }
  m_formula.reset() ;
  if ( nullptr == tree ) { return false ; }
  m_formula = std::make_unique<Ostap::Formula> ( m_cuts , tree) ; 
  return ok  () ;    
}
// ============================================================================
// reset formula unisy new tree 
// ============================================================================
void Ostap::SelectorWithCuts::reset_formula ( TTree* tree ) 
{
  set_tree ( tree ) ;
  m_formula.reset() ;
  //
 Ostap:Assert ( !get_tree() || m_cuts.empty() || make_formula ( tree ) , 
                "Invalid formula "        ,
                "Ostap::SelectorWithCuts" ,
                INVALID_FORMULA , __FILE__ , __LINE__ ) ;
}   
// ============================================================================
// notify 
// ============================================================================
Bool_t Ostap::SelectorWithCuts::Notify() 
{
  if ( m_formula ) { m_formula->Notify() ; }
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
  Ostap::Assert ( tree                                   , 
                  "Invalid tree"                         , 
                  "Ostap::SelectorWithCuts::SlaveBegin"  ,
                  INVALID_TREE , __FILE__ , __LINE__     ) ;
  reset_formula ( tree  ) ;
  Ostap::Selector::SlaveBegin  ( tree ) ;
}
// ============================================================================
// check the entry 
// ============================================================================
bool Ostap::SelectorWithCuts::good_entry  ( Long64_t entry ) 
{
  if ( Ostap::Selector::GetEntry ( entry ) < 0 ) 
    {
      Ostap::Selector::Abort ( "" , TSelector::kAbortFile ) ;
      return false ; 
    }
  //
  if ( !m_cuts.empty() && m_formula && m_formula->GetNdim() && !m_formula ->evaluate() )
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
  if ( !good_entry ( entry ) ) 
    {
      increment_event () ;
      return false       ;   // RETURN 
    }
  //
  // increment the event counter for good events 
  ++m_good  ;
  //
  // total event counter 
  increment_event () ;
  return process_entry () ;
  // ==========================================================================
}
// ============================================================================
// terminate the slave 
// ============================================================================
void Ostap::SelectorWithCuts::SlaveTerminate () 
{ return Ostap::Selector::SlaveTerminate() ; }
// ============================================================================
// terminate
// ============================================================================
void Ostap::SelectorWithCuts::Terminate      ()
{ return Ostap::Selector::Terminate() ; }
// ============================================================================
// get entry 
// ============================================================================
Int_t Ostap::SelectorWithCuts::GetEntry       
( Long64_t entry  , 
  Int_t    getall )
{ return Ostap::Selector::GetEntry ( entry , getall ) ; }
// ============================================================================
// Version 
// ============================================================================
Int_t Ostap::SelectorWithCuts::Version () const                
{ return Ostap::Selector::Version () ; }
// ============================================================================
// is formula OK?
// ============================================================================
bool Ostap::SelectorWithCuts::ok () const // is formula OK ? 
{ return m_cuts.empty() || ( m_formula && m_formula->ok () ) ; }
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
  m_formula.reset() ;
  m_cuts = strip ( cuts ) ;
}


// ============================================================================
#if ROOT_VERSION_CODE < ROOT_VERSION(6,36,0)
// ============================================================================
 ClassImp(Ostap::SelectorWithCuts) ;
// ============================================================================
#endif 
// ============================================================================

// ============================================================================
//                                                                     The END 
// ============================================================================
