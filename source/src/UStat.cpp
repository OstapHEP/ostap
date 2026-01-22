// $Id$
// ============================================================================
// Include files 
// ============================================================================
// STD& STL 
// ============================================================================
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <memory>
// ============================================================================
// ROOT & RooFit 
// ============================================================================
#include "TIterator.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "TH1.h"
// ============================================================================
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGlobalFunc.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Math.h"
#include "Ostap/Power.h"
#include "Ostap/UStat.h"
#include "Ostap/Iterator.h"
#include "Ostap/ProgressBar.h"
// ============================================================================
// Local
// ============================================================================
#include "local_roofit.h"
#include "status_codes.h"
// ============================================================================
/** @file
 *  Implementation file for class Analysis::UStat
 *  @see Analysis::Ustat
 *  @date 2011-09-27 
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  double getDistance
  ( const RooArgSet* x , 
    const RooArgSet* y )
  {
    if ( !x || !y ) { return -1 ; } // RETURN 
    //
    double result = 0.;
    //
    for ( auto* xa : *x )
    {
      if ( !xa ) { continue ; } 
      const RooAbsArg*  ya   = y -> find ( *xa ) ;
      if ( !ya ) { continue ; }       
      const RooAbsReal* xv   = static_cast<const RooAbsReal*> ( xa ) ;
      if ( !xv ) { continue ; }
      const RooAbsReal* yv   = static_cast<const RooAbsReal*> ( ya ) ;
      if ( !yv ) { continue ; }
      const double      val  = xv -> getVal() - yv->getVal() ;
      result                += val * val ;
      //
    }
    //
    return std::sqrt ( result ) ;
  }
  // ==========================================================================
} //                                                 end of anonymous namespace  
// ============================================================================
/*  Calculate U-statistics 
 *  @param pdf   (input) PDF
 *  @param data  (input) data 
 *  @param hist  (update) the histogram with U-statistics 
 *  @param args  (input)  the arguments
 *  @param tStat (output,optional) value for T-statistics 
 */
// ============================================================================
Ostap::StatusCode Ostap::UStat::calculate
( const RooAbsPdf&  pdf   , 
  const RooDataSet& data  ,  
  double&           tStat ,
  TH1*              hist  ,
  const RooArgSet*  args  ) 
{
  /// make a fake progress bar 
  Ostap::Utils::ProgressConf progress { 0 } ;
  return calculate ( progress ,
                     pdf      ,
                     data     ,
                     tStat    ,
                     hist     ,
                     args     ) ;                  
}
// ============================================================================
/*  calculate U-statistics 
 *  @param pdf   (input) PDF
 *  @param data  (input) data 
 *  @param hist  (update) the histogram with U-statistics 
 *  @param args  (input)  the arguments
 *  @param tStat (output,optional) value for T-statistics 
 */
// ============================================================================
Ostap::StatusCode Ostap::UStat::calculate
( const Ostap::Utils::ProgressConf& progress , 
  const RooAbsPdf&                  pdf      , 
  const RooDataSet&                 data     ,  
  double&                           tStat    ,
  TH1*                              hist     ,
  const RooArgSet*                  args     ) 
{
  //
  if ( nullptr == args ) { args = pdf.getObservables ( data ) ; }
  if ( nullptr == args ) { return INVALID_OBSERVABLES ; }
  //
  RooArgSet rargs {} ;
  for ( auto *a : *args )
  {
    if ( !a ) { continue ; }
    RooAbsReal* r = dynamic_cast<RooAbsReal*> ( a ) ;
    if ( !r ) { continue ; }
    rargs.add ( *r , true ) ;
  }
  //
  const unsigned int dim    = rargs.getSize   () ;
  if ( 1 > dim   ) { return INVALID_ARGSET      ; }
  /// volume of n-ball 
  const double volume = Ostap::Math::nball_volume ( dim ) ;
  //
  // const RooDataSet*  cloned = (RooDataSet*)data.Clone() ;
  // std::unique_ptr<RooDataSet> cloned  { new RooDataSet ( data ) } ;
  std::unique_ptr<RooDataSet> cloned = std::make_unique<RooDataSet> ( data , nullptr ) ;
  //
  const unsigned int num    = data.numEntries () ;
  //
  typedef std::vector<double> StatU ;
  StatU ustat ( num , 0.0 ) ; 
  //
  const RooArgSet* event_x = 0 ;
  const RooArgSet* event_y = 0 ;
  //
  // get the observables 
  std::unique_ptr<RooArgSet> observables { pdf.getObservables ( data ) } ;
  //
  Ostap::Utils::ProgressBar bar ( num , progress ) ;
  for ( unsigned int i = 0 ; i < num ; ++i , ++bar ) 
  {
    //
    // 1. Get "Event"
    event_x = data . get(i) ;      
    if ( 0 == event_x || 0 == event_x->getSize() ) { return INVALID_ENTRY ; } // RETURN 
    //
    std::unique_ptr<RooArgSet> event_i ( event_x -> selectCommon ( rargs ) ) ;
    if ( !event_i || 0 == event_i->getSize() )     { return INVALID_ENTRY ; }             // RETURN 
    //
    // 2.Evaluate PDF 
    //
    ::assign ( *observables, *event_x ) ;
    const double pdfValue = pdf . getVal ( *observables ) ;
    //
    const Int_t xs = event_i -> getSize () ;
    //
    // minimal distance 
    double min_distance  = std::numeric_limits<double>::max () ;
    for ( unsigned int j = 0 ; j < num ; ++j ) 
    {
      if ( i == j ) { continue ; }
      //
      event_y = cloned -> get ( j ) ;      
      if ( 0 == event_y || 0  == event_y->getSize() )  { return INVALID_ENTRY ; }            // RETURN 
      //
      std::unique_ptr<RooArgSet> event_j ( event_y -> selectCommon ( rargs ) ) ;
      if ( !event_j     || xs != event_j->getSize() )  { return INVALID_ENTRY ; }            // RETURN 
      //
      const double distance = getDistance ( event_i.get() , event_j.get() ) ;
      if ( distance < 0                             )  { return INVALID_ENTRY ; }
      //
      if ( 0 == j || distance <= min_distance ) { min_distance = distance  ; }
      //
    }
    //
    // volume of n-ball: 
    const double vol   = volume * Ostap::Math::POW ( min_distance , dim ) ;
    //
    const double value = std::exp ( - vol * num * pdfValue ) ;
    //
    if ( hist ) { hist -> Fill ( value ) ; }
    //
    ustat [ i ] = value ; 
    //
  }
  //
  // calculate T-statistics from U-distribution
  //
  std::stable_sort ( ustat.begin() , ustat.end() ) ;
  //
  double tS = 0 ;
  double nD = ustat.size() ;
  for ( StatU::const_iterator u = ustat.begin() ; ustat.end() != u  ; ++u ) 
  {
    const double e = ( double ( u - ustat.begin() + 1 ) )  / nD ;
    const double d = (*u) - e ;
    //
    tS += d * d ;
  }
  // finally return the value:
  tStat = tS ;
  // 
 return Ostap::StatusCode::SUCCESS ;
}
// ============================================================================
//                                                                      The END 
// ============================================================================
