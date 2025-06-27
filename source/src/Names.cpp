// ============================================================================
// Incldue files
// ============================================================================
// STD&STL
// ============================================================================
#include <cctype>
#include <string>
#include <cstdlib>
#include <algorithm>
// ============================================================================
// ROOT 
// ============================================================================
#include "TNamed.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Names.h"
#include "Ostap/Hash.h"
// ============================================================================
// Local 
// ============================================================================
#include "local_hash.h"
// ============================================================================
namespace
{
  // ==========================================================================
  /** Any of from these symbols implies that the variable name 
   *  is NOT a primitive one
   */
  const std::string s_FORMULA = " */+-%|&^()[]!$?<>=\t\n\r\t\v";
  // ==========================================================================
  /// Good (non-whetespace symbol) 
  static const auto s_GOOD_SYMBOL =
    [] ( unsigned char c ) -> bool { return !std::isspace ( c ) ; } ;
  // ==========================================================================
  /// Predefined names for the trivial selections 
  const std::vector<std::string> s_TRIVIAL {
    ""     , "1"    , "1."   , "1.0"  ,
    "true" , "True" , "TRUE" , 
    "yes"  , "Yes"  , "YES"  } ;
  // ==========================================================================
} // ==========================================================================
// ===============================================================================
/*  Genrate some valid (ranodm) name
 *  @param prefix (INPUT) prefix 
 *  @param name   (INPUT) the base name 
 *  @param named  (INPUT) the TNNamed object 
 *  @param random (INPUT) use random generator 
 *  @return some random name 
 *  @see TNamed
 */  
// ===============================================================================
std::string Ostap::tmp_name 
( const std::string&  prefix ,
  const std::string&  name   ,
  const TNamed*       named  ,
  const bool          random ) 
{
  std::size_t hv = 
    nullptr == named || random ? 
    Ostap::Utils::hash_combiner ( prefix , name , random ) :
    Ostap::Utils::hash_combiner ( prefix , name , random , 
                                  std::string ( named -> GetName  () ) , 
                                  std::string ( named -> GetTitle () ) ) ;
  //
  if ( random ) { hv = Ostap::Utils::hash_combiner ( prefix , hv , std::rand() ) ; }
  //
  return prefix + std::to_string ( hv ) ;
}
// ============================================================================
/* Is  the name "primitive" 
 *  the name os primitive if it can corresponds to 
 *  the varible name in container (e.g. TTree, RooabsData)
 *  @see TTree
 *  @see RooAbsData 
 */
// ============================================================================
bool Ostap::primitive ( const std::string& name ) 
{ return std::string::npos == name.find_first_of ( s_FORMULA ) ; }
// ============================================================================
/*  Trivial name for selection ? 
 *  - "1"
 *  - "1."
 *  - "1.0"
 *  - "true" , "True" , 'TRUE'
 *  - "yes"  . "Yes"  , "YES"
 *  - omly whitespaces 
 */
// ============================================================================
bool Ostap::trivial ( const std::string& selection )
{
  return
    // empty string 
    selection.empty ()                                                        ? true :
    // only whotespaces 
    selection.end () == std::find_if ( selection.begin () ,
				       selection.end   () , s_GOOD_SYMBOL ) ? true :
    // predefined trivials 
    s_TRIVIAL.end () != std::find    ( s_TRIVIAL.begin () ,
				       s_TRIVIAL.end   () , strip ( selection ) ) ;				         
}
// ============================================================================
// remove leading and trailing spaces 
// ============================================================================
std::string Ostap::strip ( const std::string& s )
{
  //
  // position of the first non-empty (spaxe)  symbol 
  const std::string::size_type s1 =
    std::find_if ( s.begin  () , s.end  ()      , s_GOOD_SYMBOL ) - s.begin  () ;
  // position on of the last no-empty (space) symbol
  const std::string::size_type s2 =
    std::find_if ( s.rbegin () , s.rend () - s1 , s_GOOD_SYMBOL ) - s.rbegin () ;
  //
  return s.substr ( s1 , s.size() - s2 ) ;
}
// ===============================================================================
// convert to lower case
// ===============================================================================
std::string Ostap::tolower
( const std::string& name )
{
  std::string result { name } ;
  std::transform ( name   .begin () ,
		   name   .end   () ,
		   result .begin () ,
		   [] ( unsigned char c )
		   { return std::tolower ( c ) ; } ) ;
  return result ;
}
// ===============================================================================
// convert to upper case
// ===============================================================================
std::string Ostap::toupper
( const std::string& name )
{
  std::string result { name } ;
  std::transform ( name   .begin () ,
		   name   .end   () ,
		   result .begin () ,
		   [] ( unsigned char c )
		   { return std::toupper ( c ) ; } ) ;
  return result ;
}
// ============================================================================

// ============================================================================
//                                                                     The END 
// ============================================================================
