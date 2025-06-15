// ============================================================================
// Include files
// ============================================================================
// STD/STL
// ============================================================================
#include <limits>
// ============================================================================
// ROOT
// ============================================================================
#include "TVirtualTreePlayer.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Statistic.h"
// ============================================================================
/** @file
 *  Implementation file for classes Ostap::Math::Statistic and Ostap::Math::wStatistics
 *  @see Ostap::Math::Statistic
 *  @see Ostap::Math::WStatistic
 *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
 *  @date 2025-06-11
 */
// ===========================================================================
namespace
{
    // =======================================================================
    static_assert ( std::numeric_limits<Ostap::EventIndex>::is_specialized     ,
                  "numeric_limits<Ostap::EventIndex>      is NOT specialized!" ) ;
    static_assert ( std::numeric_limits<Ostap::EventIndex>::is_integer         ,
                  "numeric_limits<Ostap::EventIndex>      is NOT integer!"     ) ;
    static_assert (TVirtualTreePlayer::kMaxEntries <= std::numeric_limits<Ostap::EventIndex>::max(), 
                  "numeric_limits<Ostap::EventIndex>::max is too small"        ) ;
    // =======================================================================
}
// ===========================================================================
constexpr Ostap::EventIndex Ostap::FirstEvent { 0 } ;
constexpr Ostap::EventIndex Ostap::LastEvent  { TVirtualTreePlayer::kMaxEntries } ;
// ===========================================================================
// virtual destructor
// ===========================================================================
Ostap::Math::Statistic::~Statistic(){}
// ===========================================================================
// virtual destructor
// ===========================================================================
Ostap::Math::WStatistic::~WStatistic(){}
// ============================================================================
//                                                                      The END
// ============================================================================

