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
#include "Ostap/Power.h"
#include "Ostap/UStat.h"
#include "Ostap/Iterator.h"
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
  double getDistance ( const RooArgSet* x , 
                       const RooArgSet* y )
  {
    if ( 0 == x || 0 == y ) { return -1 ; } // RETURN 
    //
    double result = 0.;
    //
    Ostap::Utils::Iterator xIter ( *x ) ;
    Ostap::Utils::Iterator yIter ( *y ) ;
    //
    RooRealVar* xVar = 0 ;
    RooRealVar* yVar = 0 ;
    while ( ( xVar = (RooRealVar*)xIter -> Next() ) && 
            ( yVar = (RooRealVar*)yIter -> Next() )    ) 
    {
      const double val  = xVar->getVal() - yVar->getVal() ;
      result           += val*val ;
    }
    //
    return std::sqrt ( result ) ;
  }
  // ==========================================================================
  /// get the volume of n-ball with unit radius 
  double nBallVolume ( const unsigned int n )
  {
    return 
      0 == n ?  0.0 :  // 0-ball : nothing 
      1 == n ?  2.0 :  // 1-ball : interval 
      2 == n ? M_PI :  // 2-ball : circle 
      2 * M_PI / double ( n ) * nBallVolume ( n - 2 ) ;
  }
  // ==========================================================================
} //                                                 end of anonymous namespace  
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
( const RooAbsPdf&  pdf   , 
  const RooDataSet& data  ,  
  TH1&              hist  ,
  double&           tStat ,
  RooArgSet*        args  ) 
{
  //
  if ( 0 == args ) { args = pdf.getObservables ( data ) ; }
  if ( 0 == args ) { return Ostap::StatusCode( InvalidArgs ) ; }
  //
  const unsigned int dim    = args->getSize   () ;
  if ( 1 > dim   ) { return Ostap::StatusCode( InvalidDims ) ; }
  const double volume = nBallVolume ( dim ) ;
  //
  typedef std::vector<double> TStat ;
  TStat tstat ;
  //
  const RooDataSet*  cloned = (RooDataSet*)data.Clone() ;
  //
  const unsigned int num    = data.numEntries () ;
  //
  const RooArgSet * event_x = 0 ;
  const RooArgSet * event_y = 0 ;
  //
  for ( unsigned int i = 0 ; i < num ; ++i ) 
  {
    //
    // 1. Get "Event"
    event_x = data . get(i) ;      
    if ( 0 == event_x || 0 == event_x->getSize() ) 
    { return Ostap::StatusCode ( InvalidItem1 ) ; }             // RETURN 
    //
    std::unique_ptr<RooArgSet> event_i ( ( RooArgSet*)event_x->selectCommon( *args ) ) ;
    if ( !event_i || 0 == event_i->getSize() ) 
    { return Ostap::StatusCode ( InvalidItem2 ) ; }             // RETURN 
    //
    // 2.Evaluate PDF 
    Ostap::Utils::Iterator iter  ( *args ) ;
    RooRealVar * var = 0 ;
    while ( (var = (RooRealVar*)iter->Next() ) ) 
    { var->setVal( event_i->getRealValue( var->GetName() ) ); }
    //
    const double pdfValue = pdf . getVal( args ) ;
    //
    double min_distance  = 1.e+100 ;
    for ( unsigned int j = 0 ; j < num ; ++j ) 
    {
      if ( i == j ) { continue ; }
      //
      event_y = cloned->get(j) ;      
      if ( 0 == event_y || 0 == event_y->getSize() ) 
      { return Ostap::StatusCode ( InvalidItem1 ) ; }            // RETURN 
      //
      std::unique_ptr<RooArgSet> event_j ( ( RooArgSet*)event_y->selectCommon( *args ) ) ;
      if ( !event_j || 0 == event_j->getSize() )  
      { return Ostap::StatusCode ( InvalidItem2 ) ; }            // RETURN 
      //
      const double distance = getDistance ( event_i.get() , event_j.get() ) ;
      if ( 0 > distance ) { return Ostap::StatusCode( InvalidDist ) ; }  // RETURN 
      //
      if ( 0 == j || distance < min_distance ) 
      { min_distance = distance  ; }
      //
    }
    //
    // volume of n-ball: 
    const double val1 = volume * Ostap::Math::pow ( min_distance , dim ) ;
    //
    const double value = std::exp ( -val1 * num * pdfValue ) ;
    //
    hist.Fill ( value ) ;
    //
    tstat.push_back ( value ) ; 
    //
  } 
  delete cloned ;
  //
  // calculate T-statistics
  //
  std::sort ( tstat.begin() , tstat.end() ) ;
  double tS = 0 ;
  double nD = tstat.size() ;
  for ( TStat::const_iterator t = tstat.begin() ; tstat.end() != t  ; ++t ) 
  {
    const double e = ( double ( t - tstat.begin() + 1 ) )  / nD ;
    const double d = (*t) - e ;
    //
    tS += d * d ;
  }
  // finally return the value:
  tStat = tS ;
  //
  return Ostap::StatusCode::SUCCESS ;
}
// ============================================================================
// The END 
// ============================================================================
