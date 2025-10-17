//
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
#include "TH2.h"
#include "TH3.h"
#include "TAxis.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/HistoStat.h"
// ============================================================================
// Local
// ============================================================================
#include "status_codes.h"
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
/*  get the certain moment with respect to the specified  "value"
 *  @param histo histogram
 *  @param order the momentm order
 *  @param value central value 
 *  @return the evaluated moment
 */
// ============================================================================
double Ostap::Utils::HistoStat::moment        
( const TH1*           histo , 
  const unsigned short order ,  
  const double         value ) 
{
  if ( 0 == histo ) { return s_bad ; }                       // RETURN
  if ( 0 == order ) { return 1.0   ; }                       // RETURN
  if ( 1 == order ) { return mean ( histo ) - value ; }       // RETURN
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
/*  evaluate the uncertanty for 'bin-by-bin'-moment
 *  @param histo histogram
 *  @param order the moment parameter 
 *  @param value central value 
 *  @return the evaluated uncertanty in the moment
 */
// ============================================================================
double Ostap::Utils::HistoStat::momentErr
( const TH1*           histo , 
  const unsigned short order ) 
{
  if ( 0 == histo ) { return s_bad ; }                   // RETURN 
  const double n = nEff ( histo ) ;
  if ( 0 >= n     ) { return 0.0   ; }                   // RETURN
  //
  const double a2o = moment ( histo , 2 * order ) ;      // a(2o)
  const double ao  = moment ( histo ,     order ) ;      // a(o) 
  double result    = a2o - ao*ao ;
  result          /= n ;
  result = std::max ( 0.0 , result ) ;
  //
  return std:: sqrt ( result ) ;                         // RETURN  
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
( const TH1*           histo , 
  const unsigned short order ) 
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
( const TH1*           histo , 
  const unsigned short order ) 
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
// 2D-stuff 
// ============================================================================
/*  get the "bin-by-bin"-moment around the specified  "value"
 *  \f$ m(k_x,k_y; x , y  ) \equiv 
 *   \frac{ \sum_i (x_i - x)^{n_x} (y_i - y)^{n_y} N_i }
 *        { \sum_j N_j } \f$ 
 *  @param histo histogram
 *  @param kx    the moment parameter
 *  @param ky    the moment parameter
 *  @param xv    the central value
 *  @param yv    the central value
 *  @return the evaluated moment
 */
// ============================================================================
double Ostap::Utils::HistoStat::moment2        
( const TH2*                histo ,
  const unsigned short      kx    ,
  const unsigned short      ky    ,
  const double              xv    ,
  const double              yv    ) 
{
  if ( 0 == histo         ) { return s_bad ; }                // RETURN
  if ( 0 == kx && 0 == ky ) { return 1.0   ; }                // RETURN
  //
  // get the axis 
  const TAxis* xaxis = histo->GetXaxis() ;
  if ( 0 == xaxis ) { return s_bad ; }                        // RETURN 
  const TAxis* yaxis = histo->GetYaxis() ;
  if ( 0 == yaxis ) { return s_bad ; }                        // RETURN 
  //
  // number of bins 
  const int xBins = xaxis -> GetNbins () ;
  const int yBins = yaxis -> GetNbins () ;
  //
  double result = 0 ;
  double weight = 0 ;
  //
  // loop over all bins 
  for ( int i = 1 ; i <= xBins ; ++i ) 
  {
    const double xbin = xaxis -> GetBinCenter  ( i ) ;     // bin center 
    for ( int j = 1 ; j <= yBins ; ++j ) 
    {
      const double ybin = yaxis -> GetBinCenter  ( j ) ;      // bin center 
      // bin content 
      const double w    = histo -> GetBinContent ( i , j ) ;  // bin content
      //
      weight += w ;
      result += w 
        * std::pow ( xbin - xv , kx ) 
        * std::pow ( ybin - xv , ky ) ;
    }   
  }
  //
  if ( 0 != weight ) { result /= weight ; }
  //
  return result ;
}
// ============================================================================
/** get the "bin-by-bin"-central moment
 *  \f$ m(k_x,k_y; x , y  ) \equiv 
 *   \frac{ \sum_i (x_i - \mu_x)^{k_x} (y_i - \mu_y)^{k_y} N_i }
 *        { \sum_j N_j } \f$ 
 *  @param histo histogram
 *  @param kx    the moment parameter
 *  @param ky    the moment parameter
 *  @return the evaluated moment
 */
// ============================================================================
double Ostap::Utils::HistoStat::central_moment2
( const TH2*                histo   ,
  const unsigned short      kx      ,
  const unsigned short      ky      ) 
{
  if ( 0 == histo         ) { return s_bad ; }                // RETURN
  if ( 0 == kx && 0 == ky ) { return 1.0   ; }                // RETURN
  //
  const double xmean = moment2 ( histo , 1 , 0 ) ;
  const double ymean = moment2 ( histo , 0 , 1 ) ;
  //
  return moment2 ( histo , kx , ky , xmean , ymean ) ;
}
// ============================================================================
/** get the "bin-by-bin"-standartized moment
 *  \f$ m(k_x,k_y; x , y  ) \equiv 
 *   \frac{1}{\sigma_x^{k_x}\sigma_y^{k_y}}
 *   \frac{ \sum_i (x_i - \mu_x)^{k_x} (y_i - \mu_y)^{k_y} N_i }
 *        { \sum_j N_j } \f$ 
 *  @param histo histogram
 *  @param kx    the moment parameter
 *  @param ky    the moment parameter
 *  @return the evaluated moment
 */
// ============================================================================
double Ostap::Utils::HistoStat::std_moment2
( const TH2*                histo   ,
  const unsigned short      kx      ,
  const unsigned short      ky      ) 
{
  if ( 0 == histo         ) { return s_bad ; }                // RETURN
  if ( 0 == kx && 2 == ky ) { return 1.0   ; }                // RETURN
  if ( 2 == kx && 0 == ky ) { return 1.0   ; }                // RETURN
  //
  const double sigma2_x = central_moment2 ( histo , 2 , 0 ) ;
  const double sigma2_y = central_moment2 ( histo , 0 , 2 ) ;
  //
  return central_moment2 ( histo , kx , ky ) / 
    ( std::pow ( sigma2_x , 0.5 * kx ) * std::pow ( sigma2_y , 0.5 * ky ) ) ;
}
// ============================================================================
/*  get the "bin-by-bin"-moment around the specified  "value"
 *  \f$ m(k_x,k_y,k_z; x , y  , z ) \equiv 
 *   \frac{ \sum_i (x_i - x)^{k_x} (y_i - y)^{k_y} (z_i - z )^{k_z} N_i 
 *        { \sum_j N_j } \f$ 
 *  @param histo histogram
 *  @param kx    the moment parameter
 *  @param ky    the moment parameter
 *  @param kz    the moment parameter
 *  @param xv    the central value
 *  @param yv    the central value
 *  @param zv    the central value
 *  @return the evaluated moment
 */
// ============================================================================
double Ostap::Utils::HistoStat::moment3        
( const TH3*                histo ,
  const unsigned short      kx    ,
  const unsigned short      ky    ,
  const unsigned short      kz    ,
  const double              xv    ,
  const double              yv    ,
  const double              zv    ) 
  
{
  if ( 0 == histo                    ) { return s_bad ; }     // RETURN
  if ( 0 == kx && 0 == ky && 0 == kz ) { return 1.0   ; }     // RETURN
  //
  // get the axis 
  const TAxis* xaxis = histo->GetXaxis() ;
  if ( 0 == xaxis ) { return s_bad ; }                        // RETURN 
  const TAxis* yaxis = histo->GetYaxis() ;
  if ( 0 == yaxis ) { return s_bad ; }                        // RETURN 
  const TAxis* zaxis = histo->GetZaxis() ;
  if ( 0 == zaxis ) { return s_bad ; }                        // RETURN 
  //
  // number of bins 
  const int xBins = xaxis -> GetNbins () ;
  const int yBins = yaxis -> GetNbins () ;
  const int zBins = zaxis -> GetNbins () ;
  //
  double result = 0 ;
  double weight = 0 ;
  //
  // loop over all bins 
  for ( int i = 1 ; i <= xBins ; ++i ) 
  {
    const double xbin = xaxis -> GetBinCenter  ( i ) ;              // bin center 
    for ( int j = 1 ; j <= yBins ; ++j ) 
    {
      const double ybin = yaxis -> GetBinCenter  ( j ) ;            // bin center 
      for ( int k = 1 ; k <= zBins ; ++k ) 
      {
        const double zbin = zaxis -> GetBinCenter  ( k ) ;          // bin center
        // bin content 
        const double w    = histo -> GetBinContent ( i , j , k ) ;  // bin content
        //
        weight += w ;
        result += w
          * std::pow ( xbin - xv , kx ) 
          * std::pow ( ybin - xv , ky ) 
          * std::pow ( zbin - zv , kz ) ;
      }   
    }
  }  
  //
  if ( 0 != weight ) { result /= weight ; }
  //
  return result ;
}
// ============================================================================
/*  get the "bin-by-bin"-central moment
 *  \f$ m(k_x,k_y,k_z; x , y  , z ) \equiv 
 *   \frac{ \sum_i (x_i -\mu_x)^{k_x} 
 *                 (y_i - \mu_y)^{k_y} 
 *                 (z_i - \mu_zz )^{k_z} N_i 
 *        { \sum_j N_j } \f$ 
 *  @param histo histogram
 *  @param kx    the moment parameter
 *  @param ky    the moment parameter
 *  @param kz    the moment parameter
 *  @param xv    the central value
 *  @param yv    the central value
 *  @param zv    the central value
 *  @return the evaluated moment
 */
// ============================================================================
double Ostap::Utils::HistoStat::central_moment3
( const TH3*                histo   ,
  const unsigned short      kx      ,
  const unsigned short      ky      ,
  const unsigned short      kz      ) 
{
  if ( 0 == histo                    ) { return s_bad ; }     // RETURN
  if ( 0 == kx && 0 == ky && 0 == kz ) { return 1.0   ; }     // RETURN
  //
  const double xmean = moment3 ( histo , 1 , 0 , 0 ) ;
  const double ymean = moment3 ( histo , 0 , 1 , 0 ) ;
  const double zmean = moment3 ( histo , 0 , 0 , 1 ) ;
  //
  return moment3 ( histo , kx , ky , kz , xmean , ymean , zmean ) ;
}
// ============================================================================
/*  get the "bin-by-bin"-standartized moment
 *  \f$ m(k_x,k_y,k_z; x , y  , z ) \equiv 
 *   \frac{1}{\sigma_x^{k_x}\sigma_y^{k_y}\sigma_z^{k_z}}
 *   \frac{ \sum_i (x_i -\mu_x)^{k_x} 
 *                 (y_i - \mu_y)^{k_y} 
 *                 (z_i - \mu_zz )^{k_z} N_i 
 *        { \sum_j N_j } \f$ 
 *  @param histo histogram
 *  @param kx    the moment parameter
 *  @param ky    the moment parameter
 *  @param kz    the moment parameter
 *  @return the evaluated moment
 */
// ============================================================================
double Ostap::Utils::HistoStat::std_moment3
( const TH3*                histo   ,
  const unsigned short      kx      ,
  const unsigned short      ky      ,
  const unsigned short      kz      ) 
{
  if ( 0 == histo                    ) { return s_bad ; }                // RETURN
  if ( 2 == kx && 0 == ky && 0 == kz ) { return 1.0   ; }                // RETURN
  if ( 0 == kx && 2 == ky && 0 == kz ) { return 1.0   ; }                // RETURN
  if ( 0 == kx && 0 == ky && 2 == kz ) { return 1.0   ; }                // RETURN
  //
  const double sigma2_x = central_moment3 ( histo , 2 , 0 , 0 ) ;
  const double sigma2_y = central_moment3 ( histo , 0 , 2 , 0 ) ;
  const double sigma2_z = central_moment3 ( histo , 0 , 0 , 2 ) ;
  //
  return central_moment3 ( histo , kx , ky , kz ) / 
    ( std::pow ( sigma2_x , 0.5 * kx ) * 
      std::pow ( sigma2_y , 0.5 * ky ) *
      std::pow ( sigma2_z , 0.5 * kz ) ) ;
}
// =======================================================================
/* get Riemann sum for the histogram:
 *  \f[ R = \sum_{i} y_i \Delta x_i \f]
 *  where 
 *  - \f$ y_i \f$ is a content of bin 
 *  - \f$ \Delta x _i \f$ is a "volume" of bin
 */
// =======================================================================
double Ostap::Utils::HistoStat::riemann_sum3
( const TH3& histo )
{
  Ostap::Assert ( 3 == histo.GetDimension()               ,
		  "Invalid TH3 dimension != 3 "           ,
		  "Ostap::Utils::HistoStat::riemann_sum3" ,
		  INVALID_TH3 , __FILE__ , __LINE__       ) ;
  //
  const TAxis* xaxis = histo.GetXaxis() ;
  const TAxis* yaxis = histo.GetYaxis() ;
  const TAxis* zaxis = histo.GetZaxis() ;
  //
  Ostap::Assert ( xaxis && yaxis && zaxis                 ,
		  "Invalid TAxis*"                        ,
		  "Ostap::Utils::HistoStat::riemann_sum3" ,
		  INVALID_TH3 , __FILE__ , __LINE__       ) ;
  //
  const Int_t nx = xaxis->GetNbins() ;
  const Int_t ny = yaxis->GetNbins() ;
  const Int_t nz = zaxis->GetNbins() ;
  //
  double result = 0 ;
  for ( Int_t ix = 1 ; ix <= nx ; ++ix )
    {
      const double sx = xaxis->GetBinWidth ( ix ) ;
      for ( Int_t iy = 1 ; iy <= ny ; ++iy )
	{
	  const double sy = yaxis->GetBinWidth ( iy ) ;	  
	  for ( Int_t iz = 1 ; iz <= nz ; ++iz )
	    {
	      const double sz = zaxis->GetBinWidth ( iz ) ;	      
	      result += histo.GetBinContent ( ix , iy , iz ) * sx * sy * sz ;
	    }
	}
    }
  //
  return result ;
}
// =======================================================================
/* get Riemann sum for the histogram:
 *  \f[ R = \sum_{i} y_i \Delta x_i \f]
 *  where 
 *  - \f$ y_i \f$ is a content of bin 
 *  - \f$ \Delta x _i \f$ is a "volume" of bin
 */
// =======================================================================
double Ostap::Utils::HistoStat::riemann_sum2
( const TH2& histo )
{
  //
  if ( 3 <= histo.GetDimension() )
    {
      const TH3* h3 = dynamic_cast<const TH3*> ( &histo ) ;
      Ostap::Assert ( h3                                      ,
		      "Invalid TH3* cast"                     ,
		      "Ostap::Utils::HistoStat::riemann_sum2" ,
		      INVALID_TH3 , __FILE__ , __LINE__       ) ;    
      return riemann_sum3 ( *h3 ) ; 
    }
  //
  Ostap::Assert ( 2 == histo.GetDimension()               ,
		  "Invalid TH2 dimension != 2 "           ,
		  "Ostap::Utils::HistoStat::riemann_sum2" ,
		  INVALID_TH3 , __FILE__ , __LINE__       ) ;
  //
  const TAxis* xaxis = histo.GetXaxis() ;
  const TAxis* yaxis = histo.GetYaxis() ;
  //
  Ostap::Assert ( xaxis && yaxis                          ,
		  "Invalid TAxis*"                        ,
		  "Ostap::Utils::HistoStat::riemann_sum2" ,
		  INVALID_TH3 , __FILE__ , __LINE__       ) ;
  //
  const Int_t nx = xaxis->GetNbins() ;
  const Int_t ny = yaxis->GetNbins() ;
  //
  double result = 0 ;
  for ( Int_t ix = 1 ; ix <= nx ; ++ix )
    {
      const double sx = xaxis->GetBinWidth ( ix ) ;
      for ( Int_t iy = 1 ; iy <= ny ; ++iy )
	{
	  const double sy = yaxis->GetBinWidth ( iy ) ;	  
	  result += histo.GetBinContent ( ix , iy ) * sx * sy ;
	}
    }
  //
  return result ;
}
// =======================================================================
/* get Riemann sum for the histogram:
 *  \f[ R = \sum_{i} y_i \Delta x_i \f]
 *  where 
 *  - \f$ y_i \f$ is a content of bin 
 *  - \f$ \Delta x _i \f$ is a "volume" of bin
 */
// =======================================================================
double Ostap::Utils::HistoStat::riemann_sum1
( const TH1& histo )
{
  //
  if (      3 <= histo.GetDimension() )
    {
      const TH3* h3 = dynamic_cast<const TH3*> ( &histo ) ;
      Ostap::Assert ( h3                                      ,
		      "Invalid TH3* cast"                     ,
		      "Ostap::Utils::HistoStat::riemann_sum1" ,
		      INVALID_TH3 , __FILE__ , __LINE__       ) ;    
      return riemann_sum3 ( *h3 ) ; 
    }
  else if ( 2 <= histo.GetDimension() )
    {
      const TH2* h2 = dynamic_cast<const TH2*> ( &histo ) ;
      Ostap::Assert ( h2                                      ,
		      "Invalid TH2* cast"                     ,
		      "Ostap::Utils::HistoStat::riemann_sum1" ,
		      INVALID_TH2 , __FILE__ , __LINE__       ) ;    
      return riemann_sum2 ( *h2 ) ; 
    }
  //
  Ostap::Assert ( 1 == histo.GetDimension()               ,
		  "Invalid TH1 dimension != 1 "           ,
		  "Ostap::Utils::HistoStat::riemann_sum1" ,
		  INVALID_TH3 , __FILE__ , __LINE__       ) ;
  //
  const TAxis* xaxis = histo.GetXaxis() ;
  //
  Ostap::Assert ( xaxis                                   ,
		  "Invalid TAxis*"                        ,
		  "Ostap::Utils::HistoStat::riemann_sum1" ,
		  INVALID_TH3 , __FILE__ , __LINE__       ) ;
  //
  const Int_t nx = xaxis->GetNbins() ;
  //
  double result = 0 ;
  for ( Int_t ix = 1 ; ix <= nx ; ++ix )
    {
      const double sx = xaxis->GetBinWidth ( ix ) ;
      result += histo.GetBinContent ( ix ) * sx  ;
    }
  //
  return result ;
}
// =======================================================================
/* get Riemann sum for the histogram:
 *  \f[ R = \sum_{i} y_i \Delta x_i \f]
 *  where 
 *  - \f$ y_i \f$ is a content of bin 
 *  - \f$ \Delta x _i \f$ is a "volume" of bin
 */
// =======================================================================
double Ostap::Utils::HistoStat::riemann_sum
( const TH1& histo )
{
  //
  if (      3 <= histo.GetDimension() )
    {
      const TH3* h3 = dynamic_cast<const TH3*> ( &histo ) ;
      Ostap::Assert ( h3                                      ,
		      "Invalid TH3* cast"                     ,
		      "Ostap::Utils::HistoStat::riemann_sum"  ,
		      INVALID_TH3 , __FILE__ , __LINE__       ) ;    
      return riemann_sum3 ( *h3 ) ; 
    }
  else if ( 2 <= histo.GetDimension() )
    {
      const TH2* h2 = dynamic_cast<const TH2*> ( &histo ) ;
      Ostap::Assert ( h2                                      ,
		      "Invalid TH2* cast"                     ,
		      "Ostap::Utils::HistoStat::riemann_sum"  ,
		      INVALID_TH2 , __FILE__ , __LINE__       ) ;    
      return riemann_sum2 ( *h2 ) ; 
    }
  //
  return riemann_sum1 ( histo ) ;
}
// =======================================================================

// ============================================================================
//                                                                      The END 
// ============================================================================
