
// ============================================================================
// Incldue files 
// ============================================================================
// STD&STL
// ============================================================================
#include <memory>
#include <cstring>
// ============================================================================
// ROOT
// ============================================================================
#include "RVersion.h"
#include "TH1.h"
#include "RooPlot.h"
#include "RooPlotable.h"
#include "RooCurve.h"
#include "RooEllipse.h"
#include "RooHist.h"
// ============================================================================
// Ostap 
// ============================================================================
#include "Ostap/Plot.h"
// ============================================================================
/** @file 
 *  @see Ostap::BLOB
 *  @date 2019-04-24 
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  inline RooCurve*    copy_curve   ( const RooCurve&   p ) { return new RooCurve    ( p ) ; }
  inline RooHist*     copy_histo   ( const RooHist&    p ) { return new RooHist     ( p ) ; }
  inline RooEllipse*  copy_ellipse ( const RooEllipse& p ) { return new RooEllipse ( p ) ; }
  inline TH1*         copy_th1     ( const TH1&        p ) 
  { 
    TObject* result = p.Clone () ;
    return dynamic_cast<TH1*> ( result ) ;
  }
  // ==========================================================================
  RooPlotable* copy_plotable  ( const RooPlotable& p )
  {
    const RooCurve*    curve    = dynamic_cast<const RooCurve*>    ( &p ) ;  
    if  ( curve    ) { return copy_curve    ( *curve    ) ; }
    const RooHist*     hist     = dynamic_cast<const RooHist*>     ( &p ) ;  
    if  ( hist     ) { return copy_histo    ( *hist     ) ; }
    const RooEllipse*  ellipse  = dynamic_cast<const RooEllipse*>  ( &p ) ;  
    if  ( ellipse  ) { return copy_ellipse  ( *ellipse  ) ; }
    const TObject*     object   = dynamic_cast<const TObject*>     ( &p ) ;
    if  ( object   )
    {
      TObject* copied = object->Clone ()  ; 
      if ( copied )
      {
        RooPlotable*result = dynamic_cast<RooPlotable*> ( copied ) ;
        if ( result ) { return result ; }
      }
    }
    return nullptr ; 
  }
  // ==========================================================================
  TObject* copy_object ( const TObject& p )
  {
    const RooCurve*    curve  = dynamic_cast<const RooCurve*>    ( &p ) ;  
    if ( curve    ) { return copy_curve    ( *curve    ) ; }
    const RooHist*     hist   = dynamic_cast<const RooHist*>     ( &p ) ;  
    if ( hist     ) { return copy_histo    ( *hist     ) ; }
    const RooEllipse* ellipse = dynamic_cast<const RooEllipse*>  ( &p ) ;  
    if ( ellipse  ) { return copy_ellipse  ( *ellipse  ) ; }
    const TH1*         th1    = dynamic_cast<const TH1*>         ( &p ) ;
    if ( th1      ) { return copy_th1     ( *th1       ) ; }
    const RooPlotable* plot   = dynamic_cast<const RooPlotable*> ( &p ) ; 
    if ( plot     ) 
    { 
      RooPlotable* result = copy_plotable ( *plot ) ; 
      return dynamic_cast<TObject*> ( result ) ; 
    }
    return p.Clone () ; 
  }
  // ==========================================================================
}
// ============================================================================
// add a copy of curve to the plot 
// ============================================================================
void Ostap::Utils::add_copy
( RooPlot&          plot      , 
  const RooCurve*   curve     ,
  const char*       opts      ,
  const bool        invisible )
{
  if ( curve ) { plot.addPlotable ( ::copy_curve ( *curve ) , opts , invisible ) ; } 
}
// ============================================================================
//  add a copy of ellipse to the plot 
// ============================================================================
void Ostap::Utils::add_copy
( RooPlot&          plot      , 
  const RooEllipse* ellipse ,
  const char*       opts      ,
  const bool        invisible )
{
  if ( ellipse ) { plot.addPlotable ( ::copy_ellipse ( *ellipse ) , opts , invisible ) ; }
}
// ==========================================================================
//  add a copy of hist to the plot 
// ==========================================================================
void Ostap::Utils::add_copy
( RooPlot&          plot      , 
  const RooHist*    histo     ,
  const char*       opts      ,
  const bool        invisible )
{
  if ( histo ) { plot.addPlotable ( ::copy_histo ( *histo ) , opts , invisible ) ; } 
}
// ==========================================================================
//  add a copy of plotable to the plot 
// ==========================================================================
void Ostap::Utils::add_copy
( RooPlot&           plot      , 
  const RooPlotable* plotable  ,
  const char*        opts      ,
  const bool         invisible )
{
  if ( plotable ) { plot.addPlotable ( ::copy_plotable ( *plotable ) , opts , invisible ) ; } 
}
// ===========================================================================
//  add a copy of hist to the plot 
// ===========================================================================
void Ostap::Utils::add_copy
( RooPlot&           plot      , 
  const TH1*         th1       ,
  const char*        opts      ,
  const bool         invisible ) 
{
  if ( th1 ) { plot.addTH1( ::copy_th1 ( *th1 ) , opts , invisible ) ; } 
}
// ============================================================================
//  add a copy of object to the plot 
// ============================================================================
void Ostap::Utils::add_copy
( RooPlot&           plot      , 
  const TObject*     obj       ,
  const char*        opts      ,
  const bool         invisible )
{
  if ( obj ) { plot.addObject ( ::copy_object ( *obj ) , opts , invisible ) ; }
}
// ============================================================================  
// copy the most important internal structure 
// ============================================================================
RooPlot* Ostap::Utils::copy_plot
( const RooPlot* plot )
{
  if ( !plot ) { return nullptr ; } 
  // ==========================================================================
#if ROOT_VERSION(6,32,0) <= ROOT_VERSION_CODE // ==============================
  // ==========================================================================
  auto newplot { std::make_unique<RooPlot> ( plot->GetXaxis()->GetXmin () ,
					     plot->GetXaxis()->GetXmin () ,
					     plot->GetNbinsX ()           ) } ;
  
  // ==========================================================================
#else // ======================================================================
  // ==========================================================================
  auto newplot { std::make_unique<RooPlot> ( plot->GetXaxis()->GetXmin () ,
					     plot->GetXaxis()->GetXmin () ) } ;
  // ==========================================================================
#endif // =====================================================================
  // ==========================================================================
  // copy the content 
  const unsigned int N = plot->numItems() ;
  for ( unsigned short index = 0 ; index < N ; ++index )
  {
    TObject* object = plot->getObject ( index ) ; 
    if ( !object ) { continue ; }
    const char* name      = plot->nameOf         ( index ) ;
    TString     dopts     = plot->getDrawOptions ( name  ) ;
    bool        invisible = plot->getInvisible   ( name  ) ;
    //
    const RooCurve*    curve = dynamic_cast<RooCurve*> ( object ) ;
    if ( curve    ) { add_copy ( *newplot , curve   , dopts , invisible ) ;  continue ; }
    //
    const RooHist*     histo = dynamic_cast<RooHist*> ( object ) ;
    if ( histo    ) { add_copy  ( *newplot , histo   , dopts , invisible ) ;  continue ; }
    //
    const RooEllipse*  ellipse = dynamic_cast<RooEllipse*> ( object ) ;
    if ( ellipse  ) { add_copy  ( *newplot , ellipse , dopts , invisible ) ;  continue ; }
    //
    const RooPlotable* plotable = dynamic_cast<RooPlotable*> ( object ) ;
    if ( plotable ) { add_copy  ( *newplot , plotable , dopts , invisible ) ;  continue ; }
    //
    const TH1*         th1      = dynamic_cast<TH1*> ( object ) ;
    if ( th1      ) { add_copy  ( *newplot , th1      , dopts , invisible ) ;  continue ; }
    //
    add_copy ( *newplot , object , dopts , invisible ) ;      
    //
  }
  //
  // ==========================================================================
  newplot->SetMinimum    ( plot->GetMinimum () ) ;
  newplot->SetMaximum    ( plot->GetMaximum () ) ;
  // ==========================================================================
#if ROOT_VERSION(6,39,0) <= ROOT_VERSION_CODE // ==============================
  // ==========================================================================  
  TH1*       hnew = newplot->hist () ;
  const TH1* hold =    plot->hist () ;
  /// copy all TH1-attributes 
  if ( hnew && hold ) { hold->Copy ( *hnew ) ; }
  // ==========================================================================  
#else // ======================================================================
  // ==========================================================================     
  newplot->SetAxisColor ( plot->GetXaxis() -> GetAxisColor () , "x" ) ;
  newplot->SetAxisColor ( plot->GetYaxis() -> GetAxisColor () , "y" ) ; 
  //		
  newplot->SetNdivisions ( plot->GetNdivisions ( "x" ) , "x") ;
  newplot->SetNdivisions ( plot->GetNdivisions ( "y" ) , "y") ;
  newplot->SetNdivisions ( plot->GetNdivisions ( "z" ) , "z") ;
  // ==========================================================================      
#endif // =====================================================================
  // ==========================================================================
  //  
  return newplot.release() ;
}
// ============================================================================   
RooHist*     Ostap::Utils::copy ( const RooHist*     p )
{ return p ? ::copy_histo    ( *p ) : nullptr ; } 
RooCurve*    Ostap::Utils::copy ( const RooCurve*    p )
{ return p ? ::copy_curve    ( *p ) : nullptr ; } 
RooEllipse*  Ostap::Utils::copy ( const RooEllipse*  p )
{ return p ? ::copy_ellipse  ( *p ) : nullptr ; } 
RooPlotable* Ostap::Utils::copy ( const RooPlotable* p )
{ return p ? ::copy_plotable ( *p ) : nullptr ; } 
TH1*         Ostap::Utils::copy ( const TH1*         p )
{ return p ? ::copy_th1      ( *p ) : nullptr ; } 
TObject*     Ostap::Utils::copy ( const TObject*     p )
{ return p ? ::copy_object   ( *p ) : nullptr ; }
// ============================================================================
//                                                                      The END 
// ============================================================================

