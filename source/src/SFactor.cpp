// $Id$
// ============================================================================
// Include files 
// ============================================================================
// ROOT 
// ============================================================================
#include "TTree.h"
#include "RooAbsData.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/SFactor.h"
// ============================================================================
/** @file 
 *  Implementation file for class Analysis::SFactor
 *  @see Ostap::SFactor
 *  
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2013-04-27
 *
 */
// ==========================================================================
/*  Get sum and sum of squares for the simple branch in Tree, e.g.
 *  s-factor from usage of s_weight 
 *  The direct summation in python is rather slow, thus C++ routine helps
 *  to speedup procedure drastically 
 * 
 *  @code 
 *
 *  tree  = ...
 *  sf = tree.sFactor ( "S_sw")
 *  sumw  = sf.value () 
 *  sumw2 = sf.cov2  () 
 * 
 *  scale = sumw/sumw2 ## use in fit! 
 *
 *  @endcode 
 *  
 *  Also it is a way to get the signal component (with right uncertainty)
 *
 *  @param  tree    (INPUT) the tree 
 *  @param  varname (INPUT) name for the simple variable 
 *  @return s-factor in a form of value +- sqrt(cov2)  
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2013-04-27
 */
// ============================================================================
Ostap::Math::ValueWithError
Ostap::SFactor::sFactor 
( TTree*             tree    ,  
  const std::string& varname ) 
{
  //
  typedef Ostap::Math::ValueWithError VE ;
  //
  if ( 0 == tree                                 ) { return VE ( 0 , -100 ) ; } // INVALID TREE
  if ( varname.empty()                           ) { return VE ( 0 , -200 ) ; } // invalid branch
  if ( 0 == tree->FindBranch ( varname.c_str() ) ) { return VE ( 0 , -300 ) ; } // non-exiting branch
  if ( 0 == tree->GetBranch  ( varname.c_str() ) ) { return VE ( 0 , -400 ) ; } // non-exiting branch
  //
  Double_t   value      ;
  TBranch*   branch = 0 ;
  //
  tree -> SetBranchAddress ( varname.c_str()  , &value , &branch ) ;
  //
  const bool status = tree -> GetBranchStatus ( varname.c_str() ) ;
  //
  tree -> SetBranchStatus ( varname.c_str() , true  ) ;
  //
  Long64_t nEntries = tree->GetEntries() ;
  //
  double sumw  = 0 ;
  double sumw2 = 0 ;
  //
  for ( Long64_t i = 0 ; i< nEntries ;  ++i) 
  {
    tree -> GetEntry ( i ) ;
    //
    sumw  +=         value  ;
    sumw2 += value * value  ;
    //
  }
  // recover the status 
  tree -> SetBranchStatus ( varname.c_str() , status ) ;
  //
  return VE ( sumw , sumw2 ) ;
}
// ============================================================================
/*  Get sum and sum of squares for the weights in sataset, e.g. 
 *  s-factor from usage of s_weight 
 *  The direct summation in python is rather slow, thus C++ routine helps
 *  to speedup procedure drastically 
 * 
 *  @code 
 *
 *  data  = ...
 *  sf    = data.sFactor ()
 *  sumw  = sf.value () 
 *  sumw2 = sf.cov2  () 
 * 
 *  scale = sumw/sumw2 ## use in fit! 
 *
 *  @endcode 
 *  
 *  Also it is a way to get the signal component (with right uncertainty)
 *
 *  @param  dataset (INPUT) the tree 
 *  @return s-factor in a form of value +- sqrt(cov2)  
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2013-04-27
 */
// ============================================================================
Ostap::Math::ValueWithError Ostap::SFactor::sFactor ( const RooAbsData* data ) 
{
  //
  typedef Ostap::Math::ValueWithError VE ;
  //
  if ( !data               ) { return VE ( -1 , -1 ) ; }
  if ( !data->isWeighted() ) { return VE (  1 ,  1 ) ; }  // non-weighted dataset 
  //
  long double sumw  = 0 ;
  long double sumw2 = 0 ;
  const unsigned long nEntries = data->numEntries() ;
  //
  for ( unsigned long entry = 0 ; entry < nEntries ; ++entry )   
  {
    //
    if ( 0 == data->get ( entry)  ) { break ; }           // BREAK
    //
    sumw  += data -> weight        () ;
    sumw2 += data -> weightSquared () ;
    //
  }
  return VE ( sumw , sumw2 ) ;
}
// ============================================================================
// The END 
// ============================================================================
