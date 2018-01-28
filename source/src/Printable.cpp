// Include files 
// ============================================================================
// STD&STL
// ============================================================================
#include <sstream>
// ============================================================================
// ROOT&RooFit 
// ============================================================================
#include "RooPrintable.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Printable.h"
// ============================================================================
// push Printable object into output stream 
// ============================================================================
std::ostream& Ostap::Utils::toStream 
( const RooPrintable& obj    , 
  std::ostream&       stream , 
  const std::string&  opts   ) 
{
  obj.printMultiline ( stream                                    , 
                       obj.defaultPrintContents ( opts.c_str() ) , 
                       obj.defaultPrintStyle    ( opts.c_str() ) ) ;
  return stream ;
}
// ============================================================================
std::ostream& Ostap::Utils::toStream 
( const RooPrintable& obj    , 
  std::ostream&       stream , 
  const std::string&  opts   , 
  const int           style  )
{
  obj.printMultiline ( stream                                    , 
                       obj.defaultPrintContents ( opts.c_str() ) , 
                       style    ) ;
  return stream ;
}
// ============================================================================
// convert  Printable object into string 
// ============================================================================
std::string   Ostap::Utils::toString 
( const RooPrintable& obj  , 
  const std::string&  opts ) 
{
  std::ostringstream ss ;
  toStream ( obj , ss , opts )  ;
  return ss.str() ;
}
// ============================================================================
// convert  Printable object into string 
// ============================================================================
std::string   Ostap::Utils::toString 
( const RooPrintable& obj   , 
  const std::string&  opts  ,
  const int           style )  
{
  std::ostringstream ss ;
  toStream ( obj , ss , opts , style )  ;
  return ss.str() ;
}
// ============================================================================
/*  helper function to print printable object into string (needed for python)
 *  @param object the object
 *  @param verbose flag 
 *  @param indent indent 
 *  @return string
 *  @see RooPrintable::printMultiline 
 */
// ==========================================================================
std::string Ostap::Utils::print_printable1 
( const RooPrintable&  object  , 
  const int            content , 
  const bool           verbose , 
  std::string          indent  ) 
{
  std::ostringstream s ;
  object.printMultiline ( s , content ,  verbose , indent ) ;
  return s.str() ;
}
// ==========================================================================
/*  helper function to print printable object into string (needed for python)
 *  @param object  the object
 *  @param content the content 
 *  @param style  the style  
 *  @param indent indent 
 *  @return string 
 *  @see RooPrritable::printSstream
 */
// ==========================================================================
std::string Ostap::Utils::print_printable2
( const RooPrintable&  object  , 
  const int            content , 
  const short          style   , 
  std::string          indent  ) 
{
  std::ostringstream s ;
  object.printStream 
    ( s , content ,  
      style == RooPrintable::kInline        ? RooPrintable::kInline        :
      style == RooPrintable::kSingleLine    ? RooPrintable::kSingleLine    :
      style == RooPrintable::kStandard      ? RooPrintable::kStandard      :
      style == RooPrintable::kVerbose       ? RooPrintable::kVerbose       :
      style == RooPrintable::kTreeStructure ? RooPrintable::kTreeStructure :
      RooPrintable::kStandard      ,
      indent ) ;
  return s.str() ;
}
// ==========================================================================
/*  helper function to print printable object into tree 
 *  @param object  the object
 *  @param indent indent 
 *  @return string 
 */
// ==========================================================================
std::string  Ostap::Utils::print_printable_tree 
( const RooPrintable&  object  , 
  std::string          indent  ) 
{
  std::ostringstream s ;
  object.printTree ( s , indent ) ;
  return s.str() ;
}
// ============================================================================
// The END 
// ============================================================================
