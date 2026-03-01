// ============================================================================
#ifndef OSTAP_PDF_UTILS_H
#define OSTAP_PDF_UTILS_H 1
// ============================================================================
// Include files
// ============================================================================
// ROOT/RooFit 
// ============================================================================
#include "RooAbsPdf.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Models
  {
    // ========================================================================
    /** @class ShiftAndScale   Ostap/PDFsUtils.h
     *  Helper RooFit base class to keep "shift" and "scale" variables
     */
    class ShiftAndScale: public RooAbsPdf
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Models::ShiftAndScale, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from all parameters
      ShiftAndScale 
      ( const char*  name      ,
	const char*  title     ,
	RooAbsReal&  x         ,
	RooAbsReal&  scale     , 
	RooAbsReal&  shift     ) ;
      /// constructor from all parameters
      ShiftAndScale 
      ( const char*  name      ,
	const char*  title     ,
	RooAbsReal&  x         ,
	const double scale = 1 ,
	const double shift = 0 ) ;      
      /// "copy" constructor
      ShiftAndScale
      ( const ShiftAndScale& right ,
	const char*          name  = nullptr) ;
      /// virtual destructor
      virtual ~ShiftAndScale() ;
      // ======================================================================
    public: // some fake functionality
      // ======================================================================
      /// fake default constructor, needed just for proper (de)serialization
      ShiftAndScale () ;
      // ======================================================================
    public:
      // ======================================================================
      inline const RooAbsReal& x      () const { return m_x    .arg() ; }
      inline const RooAbsReal& xvar   () const { return m_x    .arg() ; }
      //
      inline const RooAbsReal& scale  () const { return m_scale.arg() ; }
      inline const RooAbsReal& shift  () const { return m_shift.arg() ; }
      // ======================================================================
    protected :
      // ======================================================================
      /// (OPTIONAL) explicit x -> t transformation
      inline double x2t ( const double x ) const
      {
	const double scale_ = m_scale  ;
	const double shift_ = m_shift  ;	
	return ( x - shift_ ) / scale_ ;
      }
      /// (OPTIONAL) explicit t -> x transformation
      inline double t2x ( const double t ) const
      {
	const double scale_ = m_scale  ;
	const double shift_ = m_shift  ;	
	return   t * scale_   + shift_ ;
      }
      // ======================================================================
    protected :
      // ======================================================================
      /// x observable
      RooRealProxy m_x     {} ; // x-observable
      /// scale  
      RooRealProxy m_scale {} ; // scale
      /// shift 
      RooRealProxy m_shift {} ; // shift 
      // ======================================================================
    } ; //                        The end of class Ostap::Models::ShiftAndScale
    // ========================================================================
  } //                                       The end of namespace Ostap::Models
  // ==========================================================================
} //                                                 The end of namesapce Ostap
// ============================================================================
#endif // OSTAP_PDF_UTILS_H 1
// ============================================================================
//                                                                      The END 
// ============================================================================
