// ===========================================================================
// Include files  
// ===========================================================================
// STD&STL  
// ===========================================================================
#include <limits>
// ===========================================================================
// ROOT 
// ===========================================================================
#include "TVirtualTreePlayer.h"
// ===========================================================================
// Ostap
// ===========================================================================
#include "Ostap/Types.h" 
// ===========================================================================
namespace
{
  // =========================================================================
  static_assert ( std::numeric_limits<Ostap::EventIndex>::is_specialized     ,
                  "numeric_limits<Ostap::EventIndex>      is NOT specialized!" ) ;
  static_assert ( std::numeric_limits<Ostap::EventIndex>::is_integer         ,
                  "numeric_limits<Ostap::EventIndex>      is NOT integer!"     ) ;
  static_assert (TVirtualTreePlayer::kMaxEntries <= std::numeric_limits<Ostap::EventIndex>::max(),
		 "numeric_limits<Ostap::EventIndex>::max is too small"        ) ;
  // =========================================================================
  static_assert ( std::numeric_limits<Ostap::DataType>::is_specialized       ,
                  "numeric_limits<Ostap::DataType>      is NOT specialized!" ) ;
  // =========================================================================
} //                                            The end of anonymous namespace
// ===========================================================================
//                                                                     The END 
// ===========================================================================
