// $Id$
// ============================================================================
// Include files 
// ============================================================================
// ROOT 
// ============================================================================
#include "TTree.h"
#include "TChain.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/PySelector.h"
// ============================================================================
/** @file 
 * 
 *  Implementation file for class Ostap::PySelector
 *
 *  @see Ostap::PySelector
 *  @see TPySelector 
 * 
 *  @author Vanya Belyaev Ivan.Belyaev@cern.ch
 *  @date   2011-01-21
 */
// ============================================================================
ClassImp(Ostap::Selector) ;
// ============================================================================
// constructor 
// ============================================================================
Ostap::Selector::Selector
( TTree*    tree , 
  PyObject* self ) 
  : TPySelector ( tree, self )
{}
// ============================================================================
// destructor
// ============================================================================
Ostap::Selector::~Selector(){}
// ============================================================================
/*  helper function to use TTree::Process in python 
 * 
 *  @param tree      root-tree 
 *  @param selector  the selector 
 *  
 *  @see TTree 
 *  @see TTree::Process 
 *  @see TSelector 
 *
 *  @author Vanya Belyaev Ivan.Belyaev@cern.ch
 *  @date   2011-01-21
 */
// ============================================================================
long Ostap::Process::process
( TTree*             tree     ,
  TSelector*         selector )
{
  if ( 0 == tree || 0 == selector ) { return 0 ; }
  return tree->Process ( selector ) ;
}
// ============================================================================
/*  helper function to use TTree::Process in python 
 * 
 *  @param tree      root-tree 
 *  @param selector  the selector 
 *  @param events    events to be processed 
 *  
 *  @see TTree 
 *  @see TTree::Process 
 *  @see TSelector 
 *
 *  @author Vanya Belyaev Ivan.Belyaev@cern.ch
 *  @date   2013-02-10
 */
// ============================================================================
long Ostap::Process::process
( TTree*              tree      ,
  TSelector*          selector  , 
  const unsigned long events    ) 
{
  if ( 0 == tree || 0 == selector ) { return 0 ; }
  return tree->Process ( selector , "" , events ) ;
} 
// ============================================================================
/* helper function to use TChain::Process in python 
 * 
 *  @param chain     root-chain
 *  @param selector  the selector 
 *  
 *  @see TTree 
 *  @see TTree::Process 
 *  @see TSelector 
 *
 *  @author Vanya Belyaev Ivan.Belyaev@cern.ch
 *  @date   2011-01-21
 */
// ============================================================================
long Ostap::Process::process
( TChain*            chain    ,
  TSelector*         selector )
{
  if ( 0 == chain || 0 == selector ) { return 0 ; }
  return chain -> Process ( selector ) ;
}
// ============================================================================
/* helper function to use TChain::Process in python 
 * 
 *  @param chain     root-chain
 *  @param selector  the selector 
 *  @param events    events to be processed 
 *  
 *  @see TTree 
 *  @see TTree::Process 
 *  @see TSelector 
 *
 *  @author Vanya Belyaev Ivan.Belyaev@cern.ch
 *  @date   2011-01-21
 */
// ============================================================================
long Ostap::Process::process
( TChain*             chain    ,
  TSelector*          selector ,
  const unsigned long events   ) 
{
  if ( 0 == chain || 0 == selector ) { return 0 ; }
  return chain -> Process ( selector , "" , events ) ;
}
// ============================================================================
// The END 
// ============================================================================
