// $Id$
// ============================================================================
#ifndef OSTAP_PRINTABLE_H 
#define OSTAP_PRINTABLE_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <string>
#include <iosfwd>
// ============================================================================
// forward declaration 
// ============================================================================
class RooPrintable ;  // RooFit 
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Utils
  {
    // ========================================================================
    /// convert  Printable object into string 
    std::string   toString        ( const RooPrintable& obj       , 
                                    const std::string&  opts = "" ) ;                   
    /// push Printable object into output stream 
    std::ostream& toStream        ( const RooPrintable& obj       , 
                                    std::ostream&       stream    , 
                                    const std::string&  opts = "" ) ;                   
    /// convert  Printable object into string 
    inline 
    std::string   to_string       ( const RooPrintable& obj       , 
                                    const std::string&  opts = "" ) 
    { return toString ( obj , opts ) ; }      
    /// convert  Printable object into string 
    inline
    std::string   print_printable ( const RooPrintable& obj       , 
                                    const std::string&  opts = "" ) 
    { return toString ( obj , opts ) ; }
    // ========================================================================
  } //                                            end of namespace Ostap::Utils 
  // ==========================================================================
} //                                                     end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_PRINTABLE_H
// ============================================================================
