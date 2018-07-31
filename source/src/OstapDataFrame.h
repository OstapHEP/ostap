// ============================================================================
#ifndef OSTAPDATAFRAME_H 
#define OSTAPDATAFRAME_H 1
// ============================================================================
// STD&STL
// ============================================================================
#include <string>
// ============================================================================
// ROOTm ROOT::ROOT
// ============================================================================
#include "RVersion.h" // ROOT 
// ============================================================================
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,14,0)
#include "ROOT/RDataFrame.hxx"
#else 
#include "ROOT/TDataFrame.hxx"
#endif 
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/DataFrame.h"
// ============================================================================
namespace
{
  // ==========================================================================
  /** Any of from these symbols implies that the variable name 
   *  is not a primitive one
   */
  const std::string s_SYMBOLS = " */+-%|&^()[]!$?<>=";
  /// Is a variable name "primitive" (no operations) ?
  inline bool primitive ( const std::string& name )
  { return std::string::npos == name.find_first_of ( s_SYMBOLS ) ; }
  // ==========================================================================
  /// is selection/weight a trivial one? 
  inline bool trivial ( const std::string& selection ) 
  {
    return 
      ""     == selection || 
      "1"    == selection || 
      "1."   == selection || 
      "1.0"  == selection || 
      "true" == selection || 
      std::string::npos == selection.find_first_not_of (' ') ;
  }
  // ==========================================================================
} //                                             The end of anonymous namespace
// ============================================================================
namespace  Ostap
{
  // ==========================================================================
  std::string tmp_name ( std::string         prefix , 
                         const std::string&  name   ) ;
  // ==========================================================================  
} //                                                 The end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAPDATAFRAME_H
// ============================================================================
