// ============================================================================
// Include files 
// ============================================================================
//STD&STL
// ============================================================================
// ROOT 
// ============================================================================
#include "TMatrixTSym.h"
// ============================================================================
// Ostap 
// ============================================================================
#include "Ostap/StatEntity.h"
#include "Ostap/WStatEntity.h"
#include "Ostap/Covariances.h"
#include "Ostap/StatusCode.h"
// ============================================================================
// Local 
// ============================================================================
#include "Exception.h"
#include "status_codes.h"
// ============================================================================

// ============================================================================
// Constructor from the dimensions
// ============================================================================
Ostap::Math::Covariances::Covariances
( const unsigned short N )
  : m_counters ( N )
  , m_cov2     ( N )
  , m_delta    ( N )
{
  Ostap::Assert ( 2 <= N                 , 
    "At least two variables are requred" , 
    "Ostap::Math::Covariances"           ,  
    INVALID_SIZE , __FILE__ , __LINE__   ) ; 
}
// ======================================================================
// construct for the content
// ======================================================================
Ostap::Math::Covariances::Covariances
( const Ostap::Math::Covariances::Counters&  counters  , 
  const Ostap::Math::Covariances::CovMatrix& cov2      )
  : m_counters ( counters        )
  , m_cov2     ( cov2            )
  , m_delta    ( counters.size() )
{
  Ostap::Assert ( 2  <= m_counters.size () &&
                m_counters.size() == m_cov2.GetNrows () &&
		            m_counters.size() == m_cov2.GetNcols ()         ,
		            "Invalid size if counters/covariance structure" ,
		            "Ostap::Math::NCovariance"                      ,
		            INVALID_TMATRIX  , __FILE__  , __LINE__ )  ;
}
// ======================================================================
// update the correlation counter 
// ======================================================================
Ostap::Math::Covariances&
Ostap::Math::Covariances::add 
( const std::vector<double>& input )
{
  // (1) check input 
  const std::size_t NN  = N() ;
  Ostap::Assert ( NN == input.size() , 
		  "Invalid size of input da ta "        ,
		  "Ostap::Math::NCovarianceL  Ladd"      ,
		  INVALID_DATA  , __FILE__  ,   __LINE__ ) ;

  // (2) skip infinitis 
  if ( !std::all_of ( input.begin() , input.end() , 
     [] ( const double v ) ->  bool { return std::isfinite ( v ) ; } ) ) { return *this  ; }

  // (2) number of entries 
  const double nn = n    () ;
  //
  /// update the matrix 
  if ( nn )
  {
    m_delta.resize ( NN ) ;
    for ( std::size_t i = 0 ; i < NN ; ++i )
	  {
      const double di =  input [ i ] - m_counters [ i ] .mean() ;
      m_delta [ i ] = di ;
	    for ( std::size_t j = 0 ; j <= i ; ++j )
      {
        const double dj = m_delta [ j ] ;  
        m_cov2 ( j , i ) += di * nn * dj / ( nn + 1 ) ;
        // Strange feature of ROOT: TMatrix needs to be symmetrized manually
        if ( i != j )  { m_cov2 ( i , j ) = m_cov2 ( j , i ) ; } 
	    }
	  }
  }
  /// update the counters 
  for ( std::size_t i = 0 ; i < NN ; ++i ) { m_counters [ i ] += input [ i ] ; }
  //
  return *this ;
}
// ============================================================================
// Constructor from the dimensions
// ============================================================================
Ostap::Math::WCovariances::WCovariances
( const unsigned short N )
  : m_counters ( N )
  , m_cov2     ( N )
  , m_delta    ( N )
{
  Ostap::Assert ( 2 <= N                 , 
    "At least two variables are requred" , 
    "Ostap::Math::WCovariances"          ,  
    INVALID_SIZE , __FILE__ , __LINE__   ) ; 
}
// ======================================================================
// construct for the content
// ======================================================================
Ostap::Math::WCovariances::WCovariances
( const Ostap::Math::WCovariances::Counters&  counters  , 
  const Ostap::Math::WCovariances::CovMatrix& cov2      )
  : m_counters ( counters        )
  , m_cov2     ( cov2            )
  , m_delta    ( counters.size() )
{
  Ostap::Assert ( 2  <= m_counters.size () &&
                m_counters.size() == m_cov2.GetNrows () &&
		            m_counters.size() == m_cov2.GetNcols ()         ,
		            "Invalid size if counters/covariance structure" ,
		            "Ostap::Math::WCovariances"                     ,
		            INVALID_TMATRIX  , __FILE__  , __LINE__ )  ;
}
// ======================================================================
// update the correlation counter 
// ======================================================================
Ostap::Math::WCovariances&
Ostap::Math::WCovariances::add 
( const std::vector<double>& input  , 
  const double               weight )
{
  // (1) check input 
  const std::size_t NN  = N() ;
  Ostap::Assert ( NN == input.size() , 
		  "Invalid size of input da ta "        ,
		  "Ostap::Math::NCovarianceL  Ladd"      ,
		  INVALID_DATA  , __FILE__  ,   __LINE__ ) ;
  // (2) skip zero or infinite weigts 
  if  ( !weight || !std::isfinite ( weight ) ) { return *this ;} 
  // (3) skip infinitis 
  if ( !std::all_of ( input.begin() , input.end() , 
     [] ( const double v ) ->  bool { return std::isfinite ( v ) ; } ) ) { return *this  ; }
   //
  // (4) number of entries 
  const double nn = n     () ;
  const double ww = sumw  () ; 
  const double tw = ww + weight ;  
  //
  /// update the matrix 
  if ( ww && tw )
  {
    const double nw = ww * weight / tw ;  
    m_delta.resize ( NN ) ;
    for ( std::size_t i = 0 ; i < NN ; ++i )
	  {
      const double di =  input [ i ] - m_counters [ i ] .mean() ;
      m_delta [ i ] = di ;
	    for ( std::size_t j = 0 ; j <= i ; ++j )
      {
        const double dj = m_delta [ j ] ;  
        m_cov2 ( j , i ) += di * nw * dj ;
        // Strange feature of ROOT: TMatrix needs to be symmetrized manually
        if ( i != j )  { m_cov2 ( i , j ) = m_cov2 ( j , i ) ; } 
	    }
	  }
  }
  /// update the counters 
  for ( std::size_t i = 0 ; i < NN ; ++i ) { m_counters [ i ].add ( input [ i ] , weight ) ; }
  //
  return *this ;
}

// ============================================================================
//                                                                      The END 
// ============================================================================
