// ============================================================================
#ifndef OSTAP_PLOT_H 
#define OSTAP_PLOT_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <vector>
#include <string>
#include <limits>
// ============================================================================
// ROOT 
// ============================================================================
#include "RooPlot.h"
// ============================================================================
/// forward decalrations
// ============================================================================
class RooCurve   ;  // ROOT/RooFit 
class RooEllipse ;  // ROOT/RooFit 
class RooHist    ;  // ROOT/RooFit 
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Utils 
  {
    // ========================================================================
    /// add a copy of curve to the plot 
    void add_copy
    ( RooPlot&           plot      ,
      const RooCurve*    curve     ,
      const char*        opts      ,
      const bool         invisible ) ;
    // ======================================================================
    ///  add a copy of ellipse to the plot 
    void add_copy
    ( RooPlot&           plot      ,
      const RooEllipse*  ellipse ,
      const char*        opts      ,
      const bool         invisible ) ;
    // ======================================================================
    ///  add a copy of hist to the plot 
    void add_copy
    ( RooPlot&           plot      ,
      const RooHist*     histo     ,
      const char*        opts      ,
      const bool         invisible ) ;
    // ======================================================================
    ///  add a copy of plotable to the plot 
    void add_copy
    ( RooPlot&           plot      ,
      const RooPlotable* plotable  ,
      const char*        opts      ,
      const bool         invisible ) ;
    // ======================================================================
    ///  add a copy of hist to the plot 
    void add_copy
    ( RooPlot&           plot      ,
      const TH1*         th1       ,
      const char*        opts      ,
      const bool         invisible ) ;
    ///  add a copy of object to the plot 
    void add_copy
    ( RooPlot&           plot      ,
      const TObject*     obj       ,
      const char*        opts      ,
      const bool         invisible ) ;
    // ======================================================================
    /** @fn copy_plot
     *  make a copy/clone of RooPlot object
     *  @see RooPlot
     *  Note that RooPlot has no copy constructor!
     */
    RooPlot*     copy_plot  ( const RooPlot*     plot ) ;
    // ========================================================================    
    /// helper copy method 
    RooPlotable* copy  ( const RooPlotable* right ) ;
    /// helper copy method 
    RooCurve*    copy  ( const RooCurve*    right ) ;
    /// helper copy method 
    RooEllipse*  copy  ( const RooEllipse*  right ) ;
    /// helper copy method 
    RooHist*     copy  ( const RooHist*     right ) ;
    /// helper copy method 
    TH1*         copy  ( const TH1*         right ) ;
    /// helper copy method 
    TObject*     copy  ( const TObject*     right ) ;
    // ========================================================================    
  } //                                        The end of namespace Ostap::Utils
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
#endif // OSTAP_PLOT_H 
// ============================================================================
//                                                                      The END
// ============================================================================
