// ============================================================================
// Include files 
// ============================================================================
// STD&STL
// =============================================================================
#include <algorithm>
#include <cmath> 
// =============================================================================
// Ostap
// =============================================================================
#include "Ostap/Math.h"
#include "Ostap/Quantiles.h"
#include "Ostap/StatusCode.h"
// =============================================================================
// Local
// =============================================================================
#include "local_math.h"
#include "status_codes.h"
// =============================================================================
Ostap::Math::QBase::QBase
( const  bool check )
:  m_check ( check)
{}
// ==============================================================================
void Ostap::Math::QBase::throw_exception 
( const char* message , 
  const char* f , 
  const long  l ) const 
{
    Ostap::throwException 
    ( message , 
      "Ostap::Math::QBase" , 
       INVALID_DATA , f ? f : __FILE__  , 0 <= l ? l : __LINE__ ) ;
}
// =============================================================================
// constructor
// =============================================================================
Ostap::Math::HyndmanFan::HyndmanFan
( const Ostap::Math::HyndmanFan::QuantileType t     , 
  const bool                                 check ) 
: QBase ( check )
,  m_t   ( t ) 
{
    Ostap::Assert ( One <= t  && t <= Nine , 
        "Invalid QuantileType!" , 
        "Ostap::Math::HyndmanFan" , 
        INVALID_QUANTILE  , __FILE__ , __LINE__ ) ;
} 
// ==============================================================================
// constructor
// =============================================================================
Ostap::Math::ABQuantile::ABQuantile
( const double alpha , 
  const double beta  ,  
  const bool   check ) 
: QBase   ( check )
, m_alpha ( alpha )
, m_beta  ( beta  )
{
if ( s_zero  ( m_alpha    ) ) { m_alpha = 0 ; }
if ( s_equal ( m_beta , 1 ) ) { m_beta  = 1 ; }
//
  Ostap::Assert ( 0 <= m_alpha && m_alpha <= 1            ,
                  "Invalid alpha!"                        ,
                  "Ostap::Math::ABQuantile"               ,
                  INVALID_QUANTILE , __FILE__  , __LINE__ ) ;
  Ostap::Assert ( 0 <= m_beta && m_beta <= 1              ,
                  "Invalid beta!"                         ,
                  "Ostap::Math::ABQuantile"               ,
                  INVALID_QUANTILE , __FILE__  , __LINE__ ) ;

} 
// =============================================================================


// =============================================================================
//                                                                       The END
// ============================================================================= 