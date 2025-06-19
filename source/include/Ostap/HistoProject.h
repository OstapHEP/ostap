// ===========================================================================
#ifndef OSTAP_HISTOPROJECT_H 
#define OSTAP_HISTOPROJECT_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <limits>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatusCode.h"
// ============================================================================
// Forward declarations 
// =============================================================================
class TH1        ;  // ROOT 
class TH2        ;  // ROOT 
class TH3        ;  // ROOT 
class TTree      ;  // ROOT 
// =============================================================================
class RooAbsData ; // RooFit 
class RooAbsReal ; // RooFit 
// =============================================================================
#include "Ostap/Statistic.h"
// =============================================================================
/** @file Ostap/HistoProject 
 *  Collecton of useful counters for "projetcion" of 
 *  data sources (TTree/RooAbsData) into data accumulators:
 * - counters
 * - histograms 
 * - paramterizations 
 * @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 * @data 2025-06-19
 */
namespace Ostap
{  
  // ===========================================================================
  namespace  Utils
  {
    // =======================================================================
    /** @class TH1_Statistics 
     *  1D-histgram as statistics 
     *  @see Ostap::Math::Statistics
     */
    class TH1_Statistic : public Ostap::Math::Statistic
    {
      // ======================================================================
      /// constructor 
      TH1_Statistic  ( TH1* histo ) ;
      // ======================================================================
    public:
      // =====================================================================
      /** reset the histogram
       *  @see Ostap::Math::Statistic::reset 
       *  @see TH1::Reset 
       */
      void reset  () override ;
      // =====================================================================
      /** update the content 
       *  @see Ostap::Math::Statistic::update
       *  @see TH1::Fill
       */
      void update
      ( const double x ) override ; 
      // =====================================================================
    private:
      // ======================================================================
      TH1* m_histo { nullptr } ;
      // ======================================================================      
    } ;
    // =======================================================================
    /** @class TH1_@Statistic 
     *  1D-histgram as statistics 
     *  @see Ostap::Math::WStatistic
     */
    class TH1_WStatistic : public Ostap::Math::WStatistic
    {
      // ======================================================================
      TH1_WStatistic  ( TH1* histo ) ;
      // ======================================================================
    public:
      // =====================================================================
      /** reset the histogram
       *  @see Ostap::Math::WStatistic::reset 
       *  @see TH1::Reset 
       */
      void reset  () override ;
      // =====================================================================
      /** update the content 
       *  @see Ostap::Math::WStatistic::update
       *  @see TH1::Fill
       */
      void update
      ( const double x ,
	const double w ) override ; 
      // =====================================================================
    private:
      // ======================================================================
      TH1* m_histo { nullptr } ;
      // ======================================================================      
    } ;
    // =======================================================================
    /** @class TH2_Statistic 
     *  2D-histgram as statistics
     *  @see Ostap::Math::Statistics2
     */
    class TH2_Statistic : public Ostap::Math::Statistic2
    {
      // ======================================================================
      /// constructor 
      TH2_Statistic  ( TH2* histo ) ;
      // ======================================================================
    public:
      // =====================================================================
      /** reset the histogram
       *  @see Ostap::Math::Statistic2::reset 
       *  @see TH1::Reset 
       */
      void reset  () override ;
      // =====================================================================
      /** update the content 
       *  @see Ostap::Math::Statistic2::update
       *  @see TH1::Fill
       */
      void update
      ( const double x ,
	const double y ) override ; 
      // =====================================================================
    private:
      // ======================================================================
      TH2* m_histo { nullptr } ;
      // ======================================================================      
    } ;
    // =======================================================================
    /** @class TH2_WStatistic 
     *  2D-histgram as statistics 
     *  @see Ostap::Math::WStatistic2
     */
    class TH2_WStatistic : public Ostap::Math::WStatistic2
    {
      // ======================================================================
      TH2_WStatistic  ( TH2* histo ) ;
      // ======================================================================
    public:
      // =====================================================================
      /** reset the histogram
       *  @see Ostap::Math::WStatistic2::reset 
       *  @see TH1::Reset 
       */
      void reset  () override ;
      // =====================================================================
      /** update the content 
       *  @see Ostap::Math::WStatistic2::update
       *  @see TH1::Fill
       */
      void update
      ( const double x ,
	const double y , 
	const double w ) override ; 
      // =====================================================================
    private:
      // ======================================================================
      TH2* m_histo { nullptr } ;
      // ======================================================================      
    } ;
    // ========================================================================
    /** @class TH3_Statistic
     *  3D-histgram as statistics 
     *  @see Ostap::Math::Statistic3
     */
    class TH3_Statistic : public Ostap::Math::Statistic3
    {
      // ======================================================================
      /// constructor 
      TH3_Statistic  ( TH3* histo ) ;
      // ======================================================================
    public:
      // =====================================================================
      /** reset the histogram
       *  @see Ostap::Math::Statistic3::reset 
       *  @see TH1::Reset 
       */
      void reset  () override ;
      // =====================================================================
      /** update the content 
       *  @see Ostap::Math::Statistic3::update
       *  @see TH1::Fill
       */
      void update
      ( const double x ,
	const double y ,
	const double z ) override ; 
      // =====================================================================
    private:
      // ======================================================================
      TH3* m_histo { nullptr } ;
      // ======================================================================      
    } ;
    // =======================================================================
    /** @class TH3_WStatistic 
     *  3D-histgram as statistics 
     *  @see Ostap::Math::WStatistic3
     */
    class TH3_WStatistic : public Ostap::Math::WStatistic3
    {
      // ======================================================================
      TH3_WStatistic  ( TH3* histo ) ;
      // ======================================================================
    public:
      // =====================================================================
      /** reset the histogram
       *  @see Ostap::Math::WStatistic3::reset 
       *  @see TH1::Reset 
       */
      void reset  () override ;
      // =====================================================================
      /** update the content 
       *  @see Ostap::Math::WStatistics3::update
       *  @see TH1::Fill
       */
      void update
      ( const double x ,
	const double y , 
	const double z , 
	const double w ) override ; 
      // =====================================================================
    private:
      // ======================================================================
      TH3* m_histo { nullptr } ;
      // ======================================================================      
    } ;
    // ========================================================================
  } //                                        The end of namespace Ostap::Utils 
  // ===========================================================================
} // The end of namespace Ostap 
// =============================================================================
#endif // OSTAP_HPROJECT_H
// =============================================================================
//                                                                       The END 
// =============================================================================
