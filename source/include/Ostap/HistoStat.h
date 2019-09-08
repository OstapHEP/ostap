// ============================================================================
#ifndef OSTAP_HSTATS_H
#define OSTAP_HSTATS_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <utility>
// ============================================================================
// forward declarations
// ============================================================================
class TH1 ; // ROOT 
class TH2 ; // ROOT 
class TH3 ; // ROOT 
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Utils
  {
    // ========================================================================
    /** @class HistoStat   Ostap/HistoStat.h
     *  The collection of trivial functions to access  the
     *  statistical information for the histograms
     *  @see Ostap::Utils::HistoStats 
     *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
     *  @date 2010-10-28
     */
    class HistoStat
    {
    public :
      // ======================================================================
      /** get the "bin-by-bin"-moment around the specified  "value"
       *  @param histo histogram
       *  @param order the moment parameter
       *  @param value central value
       *  @return the evaluated moment
       */
      static double moment
      ( const TH1*                histo     ,
        const unsigned short      order     ,
        const double              value = 0 ) ;
      // ======================================================================
      /** evaluate the uncertanty for 'bin-by-bin'-moment
       *  @param histo histogram
       *  @param order the moment parameter
       *  @return the evaluated uncertanty in the moment
       */
      static double momentErr
      ( const TH1*                histo ,
        const unsigned short      order ) ;
      // ======================================================================
      /** evaluate the 'bin-by-bin'-central moment (around the mean value)
       *  @param histo histogram
       *  @param order the moment parameter
       *  @return the evaluated central moment
       */
      static double centralMoment
      ( const TH1*                histo ,
        const unsigned short      order ) ;
      // ======================================================================
      /** evaluate the uncertanty for 'bin-by-bin'-central moment
       *  (around the mean value)
       *  ( the uncertanty is calculated with O(1/n2) precision)
       *  @param histo histogram
       *  @param order the moment parameter
       *  @return the evaluated uncertanty in the central moment
       */
      static double centralMomentErr
      ( const TH1*                histo ,
        const unsigned short      order ) ;
      // ======================================================================
      /// get the skewness for the histogram
      static double skewness    ( const TH1* histo ) ;
      // ======================================================================
      /// get the error in skewness for the histogram
      static double skewnessErr ( const TH1* histo ) ;
      // ======================================================================
      /// get the kurtosis for the histogram
      static double kurtosis    ( const TH1* histo ) ;
      // ======================================================================
      /// get the error in kurtosis for the histogram
      static double kurtosisErr ( const TH1* histo ) ;
      // ======================================================================
      /// get the mean value for the histogram  (just for completeness)
      static double mean        ( const TH1* histo ) ;
      // ======================================================================
      /// get an error in the mean value
      static double meanErr     ( const TH1* histo ) ;
      // ======================================================================
      /// get the rms value for the histogram  (just for completeness)
      static double rms         ( const TH1* histo ) ;
      // ======================================================================
      /// get an error in the rms value
      static double rmsErr      ( const TH1* histo ) ;
      // ======================================================================
      /// get the effective entries   (just for completeness)
      static double nEff        ( const TH1* histo ) ;
      // ======================================================================
      // 2D histograms 
      // ======================================================================
      /** get the "bin-by-bin"-moment around the specified  "value"
       *  \f$ m(k_x,k_y; x , y  ) \equiv 
       *   \frac{ \sum_i (x_i - x)^{k_x} 
       *                 (y_i - y)^{k_y} N_i }
       *        { \sum_j N_j } \f$ 
       *  @param histo histogram
       *  @param kx    the moment parameter
       *  @param ky    the moment parameter
       *  @param xv    the central value
       *  @param yv    the central value
       *  @return the evaluated moment
       */
      static double moment2
      ( const TH2*                histo   ,
        const unsigned short      kx      ,
        const unsigned short      ky      ,
        const double              xv = 0  ,
        const double              yv = 0  ) ;
      // ======================================================================
      /** get the "bin-by-bin"-central moment
       *  \f$ m(k_x,k_y; x , y  ) \equiv 
       *   \frac{ \sum_i (x_i - \mu_x)^{k_x} 
       *                 (y_i - \mu_y)^{k_y} N_i }
       *        { \sum_j N_j } \f$ 
       *  @param histo histogram
       *  @param kx    the moment parameter
       *  @param ky    the moment parameter
       *  @return the evaluated moment
       */
      static double central_moment2
      ( const TH2*                histo   ,
        const unsigned short      kx      ,
        const unsigned short      ky      ) ;
      // ======================================================================
      /** get the "bin-by-bin"-standartized moment
       *  \f$ m(k_x,k_y; x , y  ) \equiv 
       *   \frac{1}{\sigma_x^{k_x}\sigma_y^{k_y}}
       *   \frac{ \sum_i (x_i - \mu_x)^{k_x} 
                         (y_i - \mu_y)^{k_y} N_i }
       *        { \sum_j N_j } \f$ 
       *  @param histo histogram
       *  @param kx    the moment parameter
       *  @param ky    the moment parameter
       *  @return the evaluated moment
       */
      static double std_moment2
      ( const TH2*                histo   ,
        const unsigned short      kx      ,
        const unsigned short      ky      ) ;
      // ======================================================================
      // 3D histograms 
      // ======================================================================
      /** get the "bin-by-bin"-moment around the specified  "value"
       *  \f$ m(k_x,k_y,k_z; x , y  , z ) \equiv 
       *   \frac{ \sum_i (x_i - x)^{k_x} 
       *                 (y_i - y)^{k_y} 
       *                 (z_i - z )^{k_z} N_i }
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
      static double moment3
      ( const TH3*                histo   ,
        const unsigned short      kx      ,
        const unsigned short      ky      ,
        const unsigned short      kz      ,
        const double              xv = 0  ,
        const double              yv = 0  ,
        const double              zv = 0  ) ;
      // ======================================================================
      /** get the "bin-by-bin"-central moment
       *  \f$ m(k_x,k_y,k_z; x , y  , z ) \equiv 
       *   \frac{ \sum_i (x_i -\mu_x)^{k_x} 
       *                 (y_i - \mu_y)^{k_y} 
       *                 (z_i - \mu_zz )^{k_z} N_i}
       *        { \sum_j N_j } \f$ 
       *  @param histo histogram
       *  @param kx    the moment parameter
       *  @param ky    the moment parameter
       *  @param kz    the moment parameter
       *  @return the evaluated moment
       */
      static double central_moment3
      ( const TH3*                histo   ,
        const unsigned short      kx      ,
        const unsigned short      ky      ,
        const unsigned short      kz      ) ;
      // ======================================================================
      /** get the "bin-by-bin"-standartized moment
       *  \f$ m(k_x,k_y,k_z; x , y  , z ) \equiv 
       *   \frac{1}{\sigma_x^{k_x}\sigma_y^{k_y}\sigma_z^{k_z}}
       *   \frac{ \sum_i (x_i -\mu_x)^{k_x} 
       *                 (y_i - \mu_y)^{k_y} 
       *                 (z_i - \mu_zz )^{k_z} N_i }
       *        { \sum_j N_j } \f$ 
       *  @param histo histogram
       *  @param kx    the moment parameter
       *  @param ky    the moment parameter
       *  @param kz    the moment parameter
       *  @return the evaluated moment
       */
      static double std_moment3
      ( const TH3*                histo   ,
        const unsigned short      kx      ,
        const unsigned short      ky      ,
        const unsigned short      kz      ) ;
      // ======================================================================
    } ;
    // ========================================================================
  } //                                        The end of namespace Ostap::Utils
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_HSTATS_H
// ============================================================================
