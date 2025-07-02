// ============================================================================
#ifndef OSTAP_NAMES_H 
#define OSTAP_NAMES_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <string>
// ============================================================================
// Forward declarations 
// ============================================================================
class TNamed ; // ROOT
// ============================================================================
namespace  Ostap
{
  // ==========================================================================
  /** Genrate some valid (ranodm) name
   *  @param prefix (INPUT) prefix 
   *  @param name   (INPUT) the base name 
   *  @param named  (INPUT) the TNNamed object 
   *  @param random (INPUT) use random generator 
   *  @return some random name 
   *  @see TNamed
   */  
  std::string tmp_name
  ( const std::string& prefix        , 
    const std::string& name          ,
    const TNamed*      named         ,
    const bool         random = true ) ;
  // ==========================================================================
  /** Genrate some valid (ranodm) name
   *  @param prefix (INPUT) prefix 
   *  @param name   (INPUT) the base name 
   *  @param name   (INPUT) the nase name 
   *  @param random (INPUT) use random generator 
   *  @return some random name 
   */  
  inline std::string tmp_name
  ( const std::string& prefix        , 
    const std::string& name          ,
    const bool         random = true ) 
  { return tmp_name ( prefix , name , nullptr , random ) ; }
  // ==========================================================================
  /** Is  the name "primitive" 
   *  the name is primitive if it could represent
   *  the variable name in container (e.g. TTree, RooabsData)
   *  - no whotespace symbols 
   *  - no symbols of operaions 
   *  @see TTree
   *  @see RooAbsData 
   */
  bool primitive    ( const std::string& name ) ;
  // ==========================================================================
  /** Trivial name for selection ? 
   *  - "1"
   *  - "1."
   *  - "1.0"
   *  - "true" , "True" , 'TRUE'
   *  - "yes"  . "Yes"  , "YES"
   *  - omly whitespaces 
   */
  bool trivial ( const std::string& selection ) ;
  // ==========================================================================
  /// remove all leading and yraling whotespaces
  std::string strip   ( const std::string& name ) ;
  /// convert to lower case
  std::string tolower ( const std::string& name ) ;
  /// convert to upper case 
  std::string toupper ( const std::string& name ) ;
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTYAP_NAMES_H
// ============================================================================
