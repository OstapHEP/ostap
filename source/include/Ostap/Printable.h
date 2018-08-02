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
/** @file Ostap/Printable.h
 *  utilities to deal with class RooPrintable
 *  @see  RooPrintable
 */   
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
    std::string   toString        ( const RooPrintable& obj       , 
                                    const std::string&  opts      ,
                                    const int           style     ) ;                   
    /// push Printable object into output stream 
    std::ostream& toStream        ( const RooPrintable& obj       , 
                                    std::ostream&       stream    , 
                                    const std::string&  opts = "" ) ;                   
    std::ostream& toStream        ( const RooPrintable& obj       , 
                                    std::ostream&       stream    , 
                                    const std::string&  opts      , 
                                    const int           style     ) ;                   
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
    /// convert  Printable object into string 
    inline
    std::string   print_printable ( const RooPrintable& obj       , 
                                    const std::string&  opts      ,
                                    const int           style     ) 
    { return toString ( obj , opts , style ) ; }

    // ==========================================================================
    /** helper function to print printable object into string (needed for python)
     *  @param object the object
     *  @param verbose flag 
     *  @param indent indent 
     *  @return string
     *  @see RooPrintable::printMultiline 
     */
    std::string  print_printable1 
    ( const RooPrintable&  object          , 
      const int            content         , 
      const bool           verbose = false , 
      std::string          indent  = ""    ) ;
    // ==========================================================================
    /** helper function to prin=t printable object into string (needed for python)
     *  @param object  the object
     *  @param content the content 
     *  @param style  the style  
     *  @param indent indent 
     *  @return string 
     *  @see RooPrintable::printSstream
     */
    std::string  print_printable2 
    ( const RooPrintable&  object       , 
      const int            content      , 
      const short          style        , 
      std::string          indent  = "" ) ;
    // ==========================================================================
    /** helper function to print printable object into tree 
     *  @param object  the object
     *  @param indent indent 
     *  @return string 
     */
    std::string  print_printable_tree 
    ( const RooPrintable&  object  , 
      std::string          indent  = "" ) ;
    // ========================================================================
  } //                                            end of namespace Ostap::Utils 
  // ==========================================================================
} //                                                     end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_PRINTABLE_H
// ============================================================================
