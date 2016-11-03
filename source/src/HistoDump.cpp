// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <cmath>
#include <vector>
#include <iostream>
#include <utility>
#include <sstream>
#include <numeric>
// ============================================================================
// ROOT
// ============================================================================
#include "TH1.h"
#include "TProfile.h"
#include "TAxis.h"
// ============================================================================
// Boost
// ============================================================================
#include "boost/format.hpp"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatusCode.h"
#include "Ostap/HistoDump.h"
// ============================================================================
namespace
{
  // ==========================================================================
  /** @struct Histo
   *  helper structure to keep the representation of the histogram
   *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
   *  @date 2009-09-18
   */
  struct Histo
  {
    // ========================================================================
    Histo() = default;
    // ========================================================================
    /** @struct Bin
     *  helper structure to keep the representation of bin
     *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
     *  @date 2009-09-18
     */
    struct Bin
    {
      // ======================================================================
      Bin ( const double h = 0  ,
            const double e = 0  ,
            const double l = -1 )
        : height ( h  )
        , error  ( e  )
        , lower ( l ) {}
      // ======================================================================
      /// bin height
      double height ;                                            //  bin height
      /// bin error
      double error  ;                                            //   bin error
      /// lower edge
      double lower  ;                                            //  lower edge
      // ======================================================================
      Bin& operator+= ( const Bin& right )
      {
        height          += right.height     ;
        const double e2  = error * error + right.error * right.error ;
        error            = std::sqrt ( e2 ) ;
        return *this ;
      }
      // ======================================================================
    } ;
    // ========================================================================
    /// get Y-max
    double maxY ( const bool withErr ) const
    {
      return std::accumulate( std::begin(bins), std::end(bins),
                              std::max ( under.height, over.height  ),
                              [&](double m, const Bin& b) {
                                return std::max ( m, withErr ? b.height + b.error : b.height ) ;
                              });
    }
    /// get Y-min
    double minY ( const bool withErr ) const
    {
      return std::accumulate( std::begin(bins), std::end(bins),
                              std::min ( under.height, over.height  ),
                              [&](double m, const Bin& b) {
                                return std::min ( m, withErr ? b.height - b.error : b.height ) ;
                              });
    }
    /// rebin the histogram
    Histo rebin ( const unsigned int bin ) const
    {
      // create new histogram
      Histo nh ;
      // copy overflow & underflow bins
      nh.under = under ;
      nh.over  = over  ;
      // rebin
      for ( unsigned int ibin = 0 ; ibin < bins.size() ; ++ibin )
      {
        const Bin& current = bins[ibin] ;
        if      ( nh.bins.empty() ) { nh.bins.push_back ( current ) ; }
        else if ( 0 == ibin % bin ) { nh.bins.push_back ( current ) ; }
        else                        { nh.bins.back()   += current   ; }
      }
      return nh ;
    }
    // find "null-bin", if any
    int nullBin () const
    {
      for ( auto ib = bins.cbegin() ; bins.cend() != ib + 1 ; ++ib )
      { if ( ib->lower <= 0 && 0 < (ib+1)->lower ) { return ib - bins.cbegin() ; } }
      return -1 ;
    }
    // ========================================================================
    typedef std::vector<Bin> Bins ;
    // ========================================================================
    /// histogram bins:
    Bins bins  ;                                              // histogram bins
    /// underflow bin
    Bin  under ;                                              //  underflow bin
    /// overflow bin
    Bin  over  ;                                              //   overflow bin
    // ========================================================================
  } ;
  // ==========================================================================
  /** "convert ROOT histogram into "local" histogram
   *   @param root (INPUT)  pointer to ROOT histogram
   *   @param hist (OUTPUT) local histogram
   *   @return StatusCode
   *   @author Vanya BELYAEV  Ivan.Belyaev@nikhef.nl
   *   @date 2009-09-19
   */
  Ostap::StatusCode _getHisto ( const TH1* root , Histo& hist )
  {
    // clear the histogram
    hist.bins.clear() ;
    //
    if ( !root      ) { return Ostap::StatusCode::FAILURE ; } // RETURN
    const TAxis* axis  = root->GetXaxis() ;
    if ( !axis      ) { return Ostap::StatusCode::FAILURE ; } // RETURN
    const int    nbins = axis->GetNbins () ;
    if ( 0 == nbins ) { return Ostap::StatusCode::FAILURE ; } // RETURN

    // underflow bin
    hist.under = Histo::Bin ( root -> GetBinContent ( 0         ) ,
                              root -> GetBinError   ( 0         ) ,
                              axis -> GetXmin       ()            ) ;
    // overflow bin
    hist.over  = Histo::Bin ( root -> GetBinContent ( nbins + 1 ) ,
                              root -> GetBinError   ( nbins + 1 ) ,
                              axis -> GetXmax       ()            ) ;
    //
    //
    for ( int ibin = 1 ; ibin <= nbins ; ++ibin )
    {
      // add to the local histo
      hist.bins.emplace_back( root -> GetBinContent ( ibin ) ,
                              root -> GetBinError   ( ibin ) ,
                              axis -> GetBinLowEdge ( ibin ) ) ;
    }
    return Ostap::StatusCode::SUCCESS ;
  }
  // ==========================================================================
  /** "convert ROOT 1D-profile into "local" histogram
   *   @param root (INPUT)  pointer to ROOT 1D-profile
   *   @param hist (OUTPUT) local histogram
   *   @return StatusCode
   *   @author Vanya BELYAEV  Ivan.Belyaev@nikhef.nl
   *   @date 2009-09-19
   */
  Ostap::StatusCode _getHisto ( const TProfile* root      ,
                         Histo&          hist      ,
                         const bool   /* spread */ )
  {
    const TH1* histo = root ;
    return _getHisto ( histo , hist ) ;
  }
  // ==========================================================================
  /** find the appropriate rebin factor
   *  (I suspect that there is some intelligent STD/STL routine for this purpose)
   *  @author Vanya BELYAEV  Ivan.BElyaev@nikhef.nl
   *  @date 2009-09-19
   */
  inline unsigned int rebin
  ( const unsigned int bins ,
    const unsigned int imax )
  {
    if ( 0 == imax  ) { return 1 ; }                                  // RETURN
    unsigned int ibin = 1 ;
    while ( bins > imax * ibin ) { ++ibin ; }
    return ibin ;                                                     // RETURN
  }
  // ==========================================================================
  /** decompose the double value into decimal significand
   *  (I suspect that there is some intelligent STD/STL routine for this purpose)
   *  @author Vanya BELYAEV  Ivan.BElyaev@nikhef.nl
   *  @date 2009-09-19
   */
  std::pair<double,int> decompose ( double v )
  {
    if      ( 0 == v ) { return { 0.0 , 0 } ;  }       // RETURN
    else if ( 1 == v ) { return { 1.0 , 0 } ;  }       // RETURN
    else if ( 0 >  v )
    {
      auto r = decompose ( -v ) ;
      return { -r.first , r.second } ;              // RETURN
    }
    else if ( 0.1 > v  )
    {
      int i       = 0 ;
      while ( 0.1 > v   ) { ++i ; v *= 10 ; }       // STUPID
      return { v , -i } ;
    }
    else if ( 1  < v  )
    {
      int i       = 0 ;
      while ( 1  <= v  ) { ++i ; v /= 10 ; }       // STUPID
      return { v , i } ;
    }
    return { v , 1 } ;
  }
  // ==========================================================================
  /** "efficient" pow-function
   *  @author Vanya BELYAEV  Ivan.Belyaev@nikhef.nl
   *  @date 2009-09-19
   */
  inline double _pow ( double x , unsigned long n )
  {
    double y = n % 2 ? x : 1;
    while ( n >>= 1 )
    {
      x = x * x;
      if ( n % 2) { y *= x; }
    }
    return y ;
  }
  // forward declaration to break cyclic dependency between rValMin and rValMax
  inline double rValMin ( double v );
  // ==========================================================================
  /** find the proper "round" value
   *  @author Vanya BELYAEV  Ivan.BElyaev@nikhef.nl
   *  @date 2009-09-19
   */
  inline double rValMax ( double v )
  {
    if      ( 0 == v  ) { return 0                   ; }             // RETURN
    else if ( 0 >  v  ) { return -1 * rValMin ( -v ) ; }             // RETURN
    // decompose the double value into decimal significand and mantissa
    std::pair<double,int> r = decompose ( v ) ;
    //
    const double f = std::ceil ( 20 * r.first ) / 2 ; // + 1 ;
    const int    l = r.second                   - 1 ;
    return 0 < l ? f * _pow ( 10 ,  l ) : f / _pow ( 10 , -l ) ;
  }
  // ==========================================================================
  /** find the proper "round" value
   *  @author Vanya BELYAEV  Ivan.BElyaev@nikhef.nl
   *  @date 2009-09-19
   */
  inline double rValMin ( double v )
  {
    if      ( 0 == v  ) { return 0                   ; }             // RETURN
    else if ( 0 >  v  ) { return -1 * rValMax ( -v ) ; }             // RETURN
    // decompose the double value into decimal significand and mantissa
    std::pair<double,int> r = decompose ( v ) ;
    const double f = std::floor ( 20 * r.first ) / 2 ; // - 1 ;
    const int    l = r.second - 1   ;
    return 0 < l ? f * _pow ( 10 ,  l ) : f / _pow ( 10 , -l ) ;
  }
  // ==========================================================================
  /** make the label for y-axis
   *  @author Vanya BELYAEV  Ivan.BElyaev@nikhef.nl
   *  @date 2009-09-19
   */
  inline std::string yLabel  ( double value )
  {
    boost::format fmt ( "%10.3g" ) ;
    fmt % value ;
    return fmt.str () ;
  }
  // ==========================================================================
  /** make the label for x-axis
   *  @author Vanya BELYAEV  Ivan.BElyaev@nikhef.nl
   *  @date 2009-09-19
   */
  inline std::string xLabel  ( const double value )
  {
    boost::format fmt ( "%9.3g" ) ;
    fmt % value ;
    return fmt.str () ;
  }
  // ==========================================================================
  /// get "correct" symbol
  char symbBin ( const Histo::Bin& bin    ,
                 const double      yLow   ,
                 const double      yHigh  ,
                 const bool        yNull  ,
                 const bool        errors )
  {
    if      ( errors && yLow  <= bin.height && bin.height < yHigh             ) { return '*' ; } // O
    else if ( errors && yHigh <  bin.height -  bin.error                      ) { return ' ' ; }
    else if ( errors && yLow  >= bin.height +  bin.error                      ) { return ' ' ; }
    else if ( errors                                                          ) { return 'I' ; }
    else if ( yLow  <= bin.height  && bin.height < yHigh                      ) { return '*' ; }
    else if ( 0 <= bin.height && yLow  <= bin.height  && 0 <  yHigh && !yNull ) { return '*' ; } // +
    else if ( 0 >  bin.height && yHigh >  bin.height  && 0 >= yLow  && !yNull ) { return '*' ; } // -
    //
    return ' ' ;
  }
  // ==========================================================================
  /** dump the local histogram in text form
   *  @author Vanya BELYAEV  Ivan.BElyaev@nikhef.nl
   *  @date 2009-09-19
   */
  std::ostream& dumpText
  ( const Histo&              histo  ,
    const std::size_t         width  ,
    const std::size_t         height ,
    const bool                errors ,
    std::ostream&             stream )
  {
    if (  40    > width  ) { return dumpText ( histo ,  40   , height , errors , stream ) ; }
    if ( 200    < width  ) { return dumpText ( histo , 200   , height , errors , stream ) ; }
    if ( 150    < height ) { return dumpText ( histo , width , 150    , errors , stream ) ; }
    if (  20    > height ) { return dumpText ( histo , width , 20     , errors , stream ) ; }
    if ( height > width  ) { return dumpText ( histo , width , width  , errors , stream ) ; }
    //
    const unsigned int nBins = histo.bins.size() ;
    if ( nBins > width  )
    {
      // rebin histogram
      Histo r = histo.rebin ( rebin ( nBins , width ) ) ;
      return dumpText ( r , width , height , errors , stream ) ;
    }
    //
    // get the Y-scale
    double yMax = std::max ( rValMax ( histo.maxY ( errors ) ) , 0.0 ) ;
    double yMin = std::min ( rValMin ( histo.minY ( errors ) ) , 0.0 ) ;

    if ( yMin == yMax ) { yMax = yMin + 1 ; }
    /// try to define the proper "Y-binning"
    std::pair<double,int> r = decompose ( yMax - yMin ) ;
    double _ny = std::ceil ( 10 * r.first ) ; //   1 <= ny < 10
    if ( 1 >=  _ny  ) { _ny = 10 ; }
    int yBins  = (int) std::max ( 1. , std::ceil ( height / _ny ) ) ;

    yBins     *= (int) _ny ;
    if ( 20     > yBins ) { yBins = 20    ; }
    const double yScale = ( yMax - yMin ) / yBins ;

    //
    const int ySkip =
      0 == yBins % 13 ? 13 :
      0 == yBins % 11 ? 11 :
      0 == yBins %  9 ?  9 :
      0 == yBins %  8 ?  8 :
      0 == yBins %  7 ?  7 :
      0 == yBins %  6 ?  6 :
      0 == yBins %  5 ?  5 :
      0 == yBins %  4 ?  4 : 10 ;

    const int xSkip =
      // 0 == nBins % 8 ? 8 :
      0 == nBins % 7 ? 7 :
      0 == nBins % 6 ? 6 :
      0 == nBins % 5 ? 5 :
      0 == nBins % 4 ? 4 : 10 ;

    int iNull = histo.nullBin() ;

    if ( 0 <= iNull ) { iNull %= xSkip ; }
    else              { iNull  = 0     ; }

    stream << std::endl ;

    for ( int yLine = -1  ; yLine < yBins ; ++yLine )
    {
      const double yHigh = yMax - yScale *   yLine       ;
      // const double yLow  = yHigh - yScale         ;
      const double yLow  = yMax - yScale * ( yLine + 1 ) ;
      //
      std::string line1 = " " ;

      const bool ynull  = ( yLow <= 0 && 0 < yHigh )    ;
      const bool yfirst = -1 == yLine || yBins -1  == yLine ;
      const bool ylab   =
        ( 0  == ( yLine + 1 ) % ySkip ) || yfirst || ynull ;

      if ( ylab ) { line1 += yLabel ( ( yLow <= 0 && 0 < yHigh ) ? 0.0 : yLow )  ; }
      else        { line1 += std::string( 10 , ' ' )                             ; }
      //
      line1 += " " ;
      //
      /// underflow bin:
      line1 += symbBin ( histo.under , yLow , yHigh , ynull , errors ) ;
      //
      line1 +=
        ynull ? "-+" :
        ylab  ? " +" : " |" ;
      //
      std::string line2 ;
      //
      for  ( auto ibin = histo.bins.cbegin() ;
             histo.bins.cend() != ibin ; ++ibin )
      {
        //char symb = ' ' ;
        const int i = ibin - histo.bins.cbegin () ;
        //
        const bool xnull =
          ibin->lower <= 0 && ( ibin + 1 ) != histo.bins.cend() && 0 < (ibin+1)->lower ;
        const bool xlab  =  iNull == i % xSkip ;
        //
        char symb = symbBin ( *ibin, yLow , yHigh , ynull , errors ) ;
        //
        if ( ' ' == symb )
        {
          if      ( ( ynull || yfirst ) && xlab ) { symb = '+' ; }
          else if (   ynull || yfirst           ) { symb = '-' ; }
          //
          else if ( ylab && xnull               ) { symb = '+' ; }
          else if (         xnull               ) { symb = '|' ; }
          //
          else if ( ylab || xlab                ) { symb = '.' ; }
          //
        }
        line2 += symb ;
        //
      }
      //
      std::string line3 =
        ynull ? "->" :
        ylab  ? "+ " : "| " ;
      //
      /// overflow bin:
      line3 += symbBin ( histo.over , yLow , yHigh , ynull , errors ) ;
      //
      stream << line1 << line2 << line3  << std::endl ;
      //
    }

    // get x-labels
    std::vector<std::string> xlabels ;
    for ( auto ib = histo.bins.cbegin() ; histo.bins.cend() != ib ; ++ib )
    {  xlabels.push_back ( xLabel ( ib->lower ) ) ; }
    // overflow& underflow  label
    const std::string oLabel = xLabel ( histo.over.lower  ) ;
    const std::string uLabel = xLabel ( histo.under.lower ) ;
    //
    static const std::string s_UNDERFLOW ( "UNDERFLOW" ) ;
    static const std::string s_OVERFLOW  ( " OVERFLOW" ) ;
    //
    //
    for ( unsigned int yLine = 0 ; yLine < 12 ; ++yLine )
    {
      std::string line = std::string ( 12 , ' ' ) ;
      //
      if ( yLine < s_UNDERFLOW.size() ) { line += s_UNDERFLOW[yLine] ; }
      else                              { line += ' '                ; }
      //
      line += ' ' ;
      //
      if ( uLabel.size() > yLine ) { line += uLabel[yLine] ; }
      else                         { line += ' '           ; }
      //
      for  ( auto ibin = histo.bins.cbegin() ; histo.bins.cend() != ibin ; ++ibin )
      {
        int ib  = ibin - histo.bins.cbegin() ;
        const bool xlab  =  ( iNull == ib % xSkip ) ;
        if ( xlab && yLine < xlabels[ib].size() ) { line += xlabels[ib][yLine] ; }
        else { line += ' ' ; }
      }
      //
      if ( oLabel.size() > yLine ) { line += oLabel[yLine] ; }
      else                         { line += ' '           ; }
      //
      line += ' ' ;
      //
      if ( yLine < s_OVERFLOW.size() ) { line += s_OVERFLOW[yLine] ; }
      else                             { line += ' '               ; }
      //
      stream << line << std::endl ;
    }
    //
    return stream ; // RETURN
  }
}
// ============================================================================
/*  dump the text representation of the histogram
 *  @param histo  (INPUT) the histogram
 *  @param stream (OUTUT) the stream
 *  @param width  (INPUT) the maximal column width
 *  @param errors (INPUT) print/plot errors
 *  @return the stream
 *  @author Vanya BELYAEV  Ivan.Belyaev@nikhef.nl
 *  @date 2009-09-19
 */
// ============================================================================
std::ostream& Ostap::Utils::Histos::histoDump
( const TH1*                histo  ,
  std::ostream&             stream ,
  const std::size_t         width  ,
  const std::size_t         height ,
  const bool                errors )
{
  const TProfile* profile = dynamic_cast<const TProfile*> ( histo ) ;
  if ( profile ) { return histoDump  ( profile , stream , width , height ) ; }
  //
  stream << std::endl ;
  if ( !histo     ) { return stream ; }  // RETURN
  Histo      hist ;
  Ostap::StatusCode sc   = _getHisto ( histo , hist ) ;
  if ( sc.isFailure() ) { return stream ; }  // RETURN
  //
  stream
    << boost::format ( " Histo Name  : \"%s\"")
    % histo -> GetName ()
    << std::endl
    << boost::format ( " Histo Title : \"%s\"")
    % histo -> GetTitle ()
    << std::endl
    << std::endl ;
  //
  stream
    << boost::format ( " Mean        : %11.5g +- %-10.4g ")
    % histo -> GetMean      ()
    % histo -> GetMeanError ()
    << std::endl
    << boost::format ( " Rms         : %11.5g +- %-10.4g ")
    % histo -> GetRMS      ()
    % histo -> GetRMSError ()
    << std::endl
    << boost::format ( " Skewness    : %11.5g            ")
    % histo -> GetSkewness ()
    << std::endl
    << boost::format ( " Kurtosis    : %11.5g            ")
    % histo -> GetKurtosis ()
    << std::endl
    << std::endl ;
  //
  stream
    << boost::format ( " Entries     :\n | %=11s | %=11s | %=11s | %=11s | %=11s |")
    % "All"
    % "Underflow"
    % "Overflow"
    % "#Equivalent"
    % "Integral"
    << std::endl
    << boost::format ( " | %=11.5g | %=11.5g | %=11.5g | %=11.5g | %=11.5g |")
    % histo -> GetEntries ()
    % histo -> GetBinContent ( 0                       )
    % histo -> GetBinContent ( histo->GetNbinsX() + 1  )
    % histo -> GetEffectiveEntries  ()
    % histo -> Integral             ()
    << std::endl
    << std::endl ;
  //
  return dumpText ( hist , width , height , errors , stream ) ;
}
// ============================================================================
/*  dump the text representation of the histogram
 *  @param histo  (INPUT) the histogram
 *  @param stream (OUTUT) the stream
 *  @param width  (INPUT) the maximal column width
 *  @param errors (INPUT) print/plot errors
 *  @return the stream
 *  @author Vanya BELYAEV  Ivan.Belyaev@nikhef.nl
 *  @date 2009-09-19
 */
// ============================================================================
std::ostream& Ostap::Utils::Histos::histoDump
( const TProfile*           histo  ,
  std::ostream&             stream ,
  const std::size_t         width  ,
  const std::size_t         height )
{
  stream << std::endl ;
  if ( ! histo     ) { return stream ; }  // RETURN
  Histo hist ;
  Ostap::StatusCode sc = _getHisto ( histo , hist , true ) ;
  if ( sc.isFailure() ) { return stream ; }  // RETURN
  //
  stream
    << boost::format ( " Profile Name  : \"%s\"")
    % histo -> GetName ()
    << std::endl
    << boost::format ( " Profile Title : \"%s\"")
    % histo -> GetTitle ()
    << std::endl
    << std::endl ;
  //
  stream
    << boost::format ( " Mean          : %11.5g ") % histo -> GetMean      ()
    << std::endl
    << boost::format ( " Rms           : %11.5g ") % histo -> GetRMS      ()
    << std::endl
    << std::endl ;
  //
  stream
    << boost::format ( " Entries       :\n | %=11s | %=11s | %=11s | %=11s |")
    % "All"
    % "Underflow"
    % "Overflow"
    % "Integral"
    << std::endl
    << boost::format ( " | %=11.5g | %=11.5g | %=11.5g | %=11.5g |")
    % histo -> GetEntries ()
    % histo -> GetBinContent ( 0                       )
    % histo -> GetBinContent ( histo->GetNbinsX() + 1  )
    % histo -> Integral             ()
    << std::endl
    << std::endl ;
  //
  return dumpText ( hist , width , height , true , stream ) ;
}
// ============================================================================
/*  dump the text representation of the histogram
 *  @param histo  (INPUT) the histogram
 *  @param width  (INPUT) the maximal column width
 *  @param errors (INPUT) print/plot errors
 *  @return string representation of the histogram
 *  @author Vanya BELYAEV  Ivan.Belyaev@nikhef.nl
 *  @date 2009-09-19
 */
// ============================================================================
std::string Ostap::Utils::Histos::histoDump
( const TH1*                histo  ,
  const std::size_t         width  ,
  const std::size_t         height ,
  const bool                errors )
{
  std::ostringstream stream ;
  histoDump ( histo , stream , width , height , errors );
  return stream.str() ;
}
// ============================================================================
/*  dump the text representation of the histogram
 *  @param histo  (INPUT) the histogram
 *  @param width  (INPUT) the maximal column width
 *  @param errors (INPUT) print/plot errors
 *  @return string representation of the histogram
 *  @author Vanya BELYAEV  Ivan.Belyaev@nikhef.nl
 *  @date 2009-09-19
 */
// ============================================================================
std::string Ostap::Utils::Histos::histoDump
( const TProfile*           histo  ,
  const std::size_t         width  ,
  const std::size_t         height )
{
  std::ostringstream stream ;
  histoDump ( histo , stream , width , height );
  return stream.str() ;
}
// ============================================================================

// ============================================================================
// The END
// ============================================================================
