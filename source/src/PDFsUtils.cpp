// ============================================================================
// Include files 
// ============================================================================
// STD & ST:
// ============================================================================
#include <limits>
#include <cmath>
// ============================================================================
// ROOT 
// ============================================================================
#include "RVersion.h"
#include "RooAbsPdf.h"
#include "RooConstVar.h"
#include "RooRealVar.h"
#include "RooRealProxy.h"
#include "RooGlobalFunc.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/PDFsUtils.h"
// ============================================================================
// Local
// ============================================================================
// #include "local_roofit.h"
// #include "local_gsl.h"
// #include "status_codes.h"
// ============================================================================
/** @file 
 *  Implementation file for functions & classes from ffrom the file Ostap/PDFsUtils.h
 *  @see Ostap/PDFsUtils.h
 *  @author Vanya BELYAEV  Ivan.Belyaev@cern.ch
 *  @date   2011-11-30
 */
// ============================================================================
//  ShiftAndScalle 
// ============================================================================
#include <iostream>
Ostap::Models::ShiftAndScale::ShiftAndScale
( const char* name  ,
  const char* title ,
  RooAbsReal& x     ,
  RooAbsReal& scale ,
  RooAbsReal& shift )
  : RooAbsPdf   ( name     , title        ) 
  , m_x         ( "!x"     , "Observable"      , this , x     ) 
  , m_scale     ( "!scale" , "scale-parameter" , this , scale ) 
  , m_shift     ( "!shift" , "shift-parameter" , this , shift )
{
  std::cerr << " I AM CONSTRUCTOR (1)" << std::endl ;
  m_x.print ( std::cerr ) ;
  std::cerr << std::endl  ; 
}
// ============================================================================
//  ShiftAndScale 
// ============================================================================
Ostap::Models::ShiftAndScale::ShiftAndScale
( const char*  name  ,
  const char*  title ,
  RooAbsReal&  x     ,
  const double scale , 
  const double shift )
  : ShiftAndScale ( name  ,
		    title ,
		    x     ,
		    RooFit::RooConst ( scale ) ,
		    RooFit::RooConst ( shift ) )
{}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Models::ShiftAndScale::ShiftAndScale
( const Ostap::Models::ShiftAndScale& right ,
  const char*                         name  )
  : RooAbsPdf ( right , name ) 
  , m_x       ( "!x"     , this , right.m_x     )
  , m_scale   ( "!scale" , this , right.m_scale )
  , m_shift   ( "!shift" , this , right.m_shift )
{
  std::cerr << " I AM CONSTRUCTOR (2) " << std::endl  ;
}
// ============================================================================
// fake default constructor, needed just for proper (de)serialization
// ============================================================================
Ostap::Models::ShiftAndScale::ShiftAndScale  ()
{
  std::cerr << " I AM CONSTRUCTOR (3) " << std::endl ;
} 
// ============================================================================
// virtual destructor
// ============================================================================
Ostap::Models::ShiftAndScale::~ShiftAndScale () {}
// ============================================================================

// ============================================================================
#if ROOT_VERSION_CODE < ROOT_VERSION(6,36,0)
// ============================================================================
ClassImp(Ostap::Models::ShiftAndScale ) 
// ============================================================================
#endif
// ============================================================================
//                                                                      The END 
// ============================================================================
