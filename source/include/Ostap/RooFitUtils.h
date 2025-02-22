// ============================================================================
#ifndef OSTAP_ROOFITUTILS_H 
#define OSTAP_ROOFITUTILS_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <ostream>
// ============================================================================
// Ostap
// ============================================================================
// #include "Ostap/ToStream.h"
// ============================================================================
// Forward declarations 
// ============================================================================
class TNamed           ; // ROOT/RooFit
class RooAbsCollection ; // ROOT/RooFit
class RooPrintable     ; // ROOT/RooFit
class RooAbsArg        ; // ROOT/RooFit
class RooAbsReal       ; // ROOT/RooFit
class RooAbsCategory   ; // ROOT/RooFit
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Utils
  {
    // ========================================================================
    /** print RooAbsCollection
     *  @see  RooAbsCollection
     */
    std::ostream& toStream
      ( const RooAbsCollection& o ,
        std::ostream&           s ) ;
    // ========================================================================
    /** print TNamed
     *  @see  TNamed
     */
    std::ostream& toStream
      ( const TNamed& o ,
        std::ostream& s ) ;
    // ========================================================================
    /** print RooPrintable 
     *  @see  RooPrintable 
     */
    std::ostream& toStream
      ( const RooPrintable& o ,
        std::ostream&       s ) ;
    // ========================================================================
    /** print RooAbsArg 
     *  @see  RooAbsArg 
     */
    std::ostream& toStream
      ( const RooAbsArg& o ,
        std::ostream&    s ) ;
    // ========================================================================    
    /** print RooAbsReal 
     *  @see  RooAbsReal
     */
    std::ostream& toStream
      ( const RooAbsReal& o ,
        std::ostream&     s ) ;
    // ========================================================================
    /** print RooAbsCategory  
     *  @see  RooAbsCategory
     */
    std::ostream& toStream
      ( const RooAbsCategory& o ,
        std::ostream&         s ) ;
    // ========================================================================
  } //                                        The end of namespace Ostap::Utils
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_ROOFITUTILS_H
// ============================================================================
