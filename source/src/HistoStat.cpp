// $Id$
// ============================================================================
// Include files 
// ============================================================================
// STD & STL 
// ============================================================================
#include <cmath>
#include <limits>
#include <climits>
// ============================================================================
// ROOT 
// ============================================================================
#include "TH1.h"
#include "TAxis.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/HistoStat.h"
// ============================================================================
/** @file
 *  Implementation file for class Gaudi::Utils::HStats
 *  @author Vanya BELYAEV Ivan.BElayev@cern.ch
 *  @see Gaudi::Utils::HistoStats 
 *  @see Gaudi::Utils::HStats
 *  @date 2010-10-28 
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  /// define local "bad" value 
  const double s_bad      = -1 * std::numeric_limits<float>::max() ;
  /// define local "bad" value 
  const long   s_long_bad =      std::numeric_limits<int>::min()   ;
  // ==========================================================================
}
// ============================================================================
/*  get the moment of the certain around the specified  "value"
 *  @param histo histogram
 *  @param order the momentm order
 *  @param value central value 
 *  @return the evaluated moment
 */
// ============================================================================
double Ostap::Utils::HistoStat::moment        
( const TH1*         histo , 
  const unsigned int order ,  
  const double       value ) 
{
  if ( 0 == histo ) { return s_bad ; }                       // RETURN
  if ( 0 == order ) { return 1.0   ; }                       // RETURN
  if ( 1 == order ) { return mean( histo ) - value ; }       // RETURN
  if ( 2 == order ) 
  {
    const double _r =         rms  ( histo ) ;
    const double _d = value - mean ( histo ) ;
    return _r *_r + _d * _d ;                                // RETURN
  }
  const double n = nEff ( histo )  ;
  if ( 0 >= n     ) { return 0.0   ; }                       // RETURN 
  //
  // get the exis 
  const TAxis* axis = histo->GetXaxis() ;
  if ( 0 == axis ) { return s_bad ; }                        // RETURN 
  //
  // number of bins 
  const int nBins = axis->GetNbins () ;
  double result = 0 ;
  double weight = 0 ;
  // loop over all bins 
  for ( int i = 1 ; i <= nBins ; ++i ) 
  {
    const double xBin = axis  -> GetBinCenter  ( i ) ;     // bin center 
    const double yBin = histo -> GetBinContent ( i ) ;     // bin content 
    //
    weight += yBin ;
    result += yBin * std::pow ( xBin - value , order ) ;
  }    
  //
  if ( 0 != weight ) { result /= weight ; }
  //
  return result ;
}
// ============================================================================
/** evaluate the uncertanty for 'bin-by-bin'-moment
 *  @param histo histogram
 *  @param order the moment parameter 
 *  @param value central value 
 *  @return the evaluated uncertanty in the moment
 */
// ============================================================================
double Ostap::Utils::HistoStat::momentErr
( const TH1*         histo , 
  const unsigned int order ) 
{
  if ( 0 == histo ) { return s_bad ; }                   // RETURN 
  const double n = nEff ( histo ) ;
  if ( 0 >= n     ) { return 0.0   ; }                   // RETURN
  //
  const double a2o = moment ( histo , 2 * order ) ;      // a(2o)
  const double ao  = moment ( histo ,     order ) ;      // a(o) 
  double result = a2o - ao*ao ;
  result        /= n ;
  result = std::max ( 0.0 , result ) ;
  //
  return std:: sqrt ( result ) ;                            // RETURN  
}
// ============================================================================
/*  evaluate the central momentum (around the mean value) 
 *  @param histo histogram
 *  @param order the momentm order
 *  @param value central value 
 *  @return the evaluated central moment
 */
// ============================================================================
double Ostap::Utils::HistoStat::centralMoment 
( const TH1*         histo , 
  const unsigned int order ) 
{
  if ( 0 == histo ) { return s_bad ; }                        // RETURN
  if ( 0 == order ) { return 1.0   ; }                        // RETURN
  if ( 1 == order ) { return 0.0   ; }                        // RETURN
  if ( 2 == order ) 
  {
    const double sigma = rms ( histo )  ;
    return sigma * sigma ;                                    // RETURN
  }
  // delegate the actual evaluation to another method:
  return moment ( histo , order , mean ( histo ) ) ;
}
// ============================================================================
/*  evaluate the uncertanty for 'bin-by-bin'-central momentum 
 *  (around the mean value)  
 *  ( the uncertanty is calculated with O(1/n2) precision)
 *  @param histo histogram
 *  @param order the moment parameter 
 *  @param value central value 
 *  @return the evaluated uncertanty in the central moment
 */
// ============================================================================
double Ostap::Utils::HistoStat::centralMomentErr
( const TH1*         histo , 
  const unsigned int order ) 
{
  //
  if ( 0 == histo ) { return s_bad ; }                    // RETURN
  const double n    = nEff ( histo ) ;
  if ( 0 >= n     ) { return 0.0   ; }                    // RETURN
  //
  const double mu2  = centralMoment ( histo , 2             ) ; // mu(2)
  const double muo  = centralMoment ( histo ,     order     ) ; // mu(o)
  const double mu2o = centralMoment ( histo , 2 * order     ) ; // mu(2o)
  const double muom = centralMoment ( histo ,     order - 1 ) ; // mu(o-1)
  const double muop = centralMoment ( histo ,     order + 1 ) ; // mu(o+1)
  //
  double result  =  mu2o  ;
  result        -=  2.0   * order * muom * muop ;
  result        -=  muo   * muo   ;
  result        +=  order * order * mu2  * muom * muom ;
  result        /=  n     ;
  result = std::max ( 0.0 , result ) ;
  //
  return std:: sqrt ( result ) ;                            // RETURN
}
// ============================================================================
// get the skewness for the histogram 
// ============================================================================
double Ostap::Utils::HistoStat::skewness ( const TH1* histo ) 
{
  if ( 0 == histo ) { return s_bad ; }                      // RETURN
  //
  const double mu3 = centralMoment ( histo , 3 ) ;
  const double s3  = std::pow ( rms ( histo ) , 3 ) ;
  //
  return ( std::fabs(s3)>0 ? mu3/s3 : 0.0 );
}
// ============================================================================
// get the error in skewness 
// ============================================================================
double Ostap::Utils::HistoStat::skewnessErr ( const TH1* histo ) 
{
  if ( 0 == histo ) { return s_bad ; }                     // RETURN 
  //
  const double n = nEff ( histo ) ;
  if ( 2 > n      ) { return 0.0   ; }                     // RETURN
  //
  double result = 6 ;
  result *= ( n - 2 )  ;
  result /= ( n + 1 ) * ( n + 3 ) ;    
  //
  return std::sqrt ( result ) ;  
}
// ============================================================================
// get the kurtosis for the histogram 
// ============================================================================
double Ostap::Utils::HistoStat::kurtosis ( const TH1* histo ) 
{
  if ( 0 == histo ) { return s_bad ; }                    // RETURN 
  //
  const double mu4 = centralMoment ( histo , 4 ) ;
  const double s4  = std::pow ( rms ( histo ) , 4 ) ;
  //
  return ( std::fabs(s4)>0 ? mu4/s4 - 3.0 : 0.0 );
}
// ============================================================================
// get the error in kurtosis
// ============================================================================
double Ostap::Utils::HistoStat::kurtosisErr ( const TH1* histo ) 
{
  if ( 0 == histo ) { return s_bad ; }                    // RETURN 
  const double n = nEff ( histo ) ;
  if ( 3 > n      ) { return 0.0 ; }                      // RETURN   
  //
  double result = 24 * n ;
  result *= ( n - 2 ) * ( n - 3 ) ;
  result /= ( n + 1 ) * ( n + 1 ) ;
  result /= ( n + 3 ) * ( n + 5 ) ;
  //
  return std::sqrt ( result ) ;  
}
// ============================================================================
// get the effective entries 
// ============================================================================
double Ostap::Utils::HistoStat::nEff ( const TH1* histo ) 
{
  if ( 0 == histo ) { return s_bad ; }
  //
  return histo -> GetEffectiveEntries () ;
}
// ============================================================================
// get the mean value for the histogram 
// ============================================================================
double Ostap::Utils::HistoStat::mean ( const TH1* histo ) 
{
  if ( 0 == histo ) { return s_bad ; }
  //
  return histo -> GetMean () ;    
}
// ============================================================================
// get an error in the mean value 
// ============================================================================
double Ostap::Utils::HistoStat::meanErr ( const TH1* histo ) 
{
  if ( 0 == histo ) { return s_bad ; }
  //
  const double n = nEff ( histo ) ;
  //
  return 0 >= n ? 0.0 : rms ( histo ) / std::sqrt ( n ) ;
}
// ============================================================================
// get the rms value for the histogram 
// ============================================================================
double Ostap::Utils::HistoStat::rms ( const TH1* histo ) 
{
  if ( 0 == histo ) { return s_bad ; }
  //
  return histo -> GetRMS () ;    
}
// ============================================================================
// get the error in rms 
// ============================================================================
double Ostap::Utils::HistoStat::rmsErr ( const TH1* histo ) 
{
  if ( 0 == histo ) { return s_bad ; }  
  //
  const double n = nEff ( histo ) ;    
  if ( 1 >=  n ) { return 0.0 ; }
  //
  double result = 2 + kurtosis ( histo ) ;
  result += 2.0 /( n - 1 ) ;
  result /= 4.0 * n ;
  result = std::max ( result , 0.0 ) ;
  //
  return rms ( histo ) * std::sqrt ( result ) ;
}
// ============================================================================
// The END 
// ============================================================================
