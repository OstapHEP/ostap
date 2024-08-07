// ===============================================================================
// Include files 
// ===============================================================================
// STD& STL
// ===============================================================================
#include <string>
#include <cstdio>
// ===============================================================================
// ROOT 
// ===============================================================================
#include "TNamed.h" 
// ===============================================================================
// Ostap
// ===============================================================================
#include "Ostap/Hash.h"
// ===============================================================================
// local
// ===============================================================================
#include "OstapDataFrame.h"
#include "local_hash.h"
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
                                  std::string ( named->GetName  () ) , 
                                  std::string ( named->GetTitle () ) ) ;
  if ( random ) { hv = Ostap::Utils::hash_combiner ( prefix , hv , rand() ) ; }
  return prefix + std::to_string ( hv ) ;
}
// ==========================================================================  
//                                                                    The END 
// ==========================================================================  
