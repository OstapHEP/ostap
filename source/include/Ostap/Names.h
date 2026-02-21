// ============================================================================
#ifndef OSTAP_NAMES_H 
#define OSTAP_NAMES_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================/
#include <string>
#include <typeinfo>
// ============================================================================
// Forward declarations 
// ============================================================================
class TNamed ; // ROOT
// ============================================================================
namespace  Ostap
{
  // ==========================================================================
  /** Generate some valid (ranodm) name
   *  @param prefix (INPUT) prefix 
   *  @param name   (INPUT) the base name 
   *  @param named  (INPUT) the TNNamed object 
   *  @param random (INPUT) use random generator 
   *  @param suffix (INPUT) suffix 
   *  @return some random name 
   *  @see TNamed
   */  
  std::string tmp_name
  ( const std::string& prefix        , 
    const std::string& name          ,
    const TNamed*      named         ,
    const std::string& suffix = ""   , 
    const bool         random = true ) ;
  // ==========================================================================
  /** Generate some valid (ranodm) name
   *  @param prefix (INPUT) prefix 
   *  @param name   (INPUT) the base name 
   *  @param name   (INPUT) the nase name 
   *  @param random (INPUT) use random generator 
   *  @return some random name 
   */  
  inline std::string tmp_name
  ( const std::string& prefix        , 
    const std::string& name          ,
    const std::string& suffix = ""   , 
    const bool         random = true ) 
  { return tmp_name ( prefix , name , nullptr , suffix , random ) ; }
  // ==========================================================================
  /** Is  the name "primitive" 
   *  the name is primitive if it could represent
   *  the variable name in container (e.g. TTree, RooAbsData)
   *  - no whitespace symbols 
   *  - no symbols of operations
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
   *  - "yes"  , "Yes"  , "YES"
   *  - only whitespaces 
   */
  bool trivial ( const std::string& selection ) ;
  // ==========================================================================
  /// remove all leading and trailing whitespaces
  std::string strip   ( const std::string& name ) ;
  /// convert to lower case
  std::string tolower ( const std::string& name ) ;
  /// convert to upper case 
  std::string toupper ( const std::string& name ) ;
  // ==========================================================================
  /** @fn class_name
   *  Get the de-mangled class name
   *  @see TClassEdit::DemangleName 
   *  @param  mangled mangled C++ class name
   *  @return demangled class name
   */
  std::string class_name
  ( const std::string& mangled ) ;
  // ==========================================================================
  /** @fn class_name
   *  Get the de-mangled class name
   *  @see TClassEdit::DemangleName 
   *  @param  mangled mangled C++ class name
   *  @return demangled class name
   */
  std::string class_name
  ( const char* mangled ) ;  
  // ==========================================================================
  /** @fn class_name
   *  Get the de-mangled class name
   *  @see TClassEdit::DemangleName 
   *  @param  into type-info object 
   *  @return demangled class name
   */
  std::string class_name
  ( const std::type_info& info  ) ;  
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_NAMES_H
// ============================================================================
