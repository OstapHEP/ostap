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
// The END 
// ============================================================================
