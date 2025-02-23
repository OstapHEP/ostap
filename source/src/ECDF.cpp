// ============================================================================
// Incldue files 
// ============================================================================
// STD&STL
// ============================================================================
#include <memory>
#include <cstring>
#include <numeric>
// ============================================================================
// Ostap 
// ============================================================================
#include "Ostap/ECDF.h"
// ============================================================================
// Local
// ============================================================================
#include "Exception.h"
// ============================================================================
/** @file 
 *  Implementation file for class Ostap::Math::ECDF
 *  @see Ostap::Math::ECDF
 *  @date 2024-09-16 
 *  @author Vanya BELYAEV 
 */
// ============================================================================
// Standard constructor from  data
// ============================================================================
Ostap::Math::ECDF::ECDF
( const Ostap::Math::ECDF::Data&  data          ,
  const bool                      complementary )
  : Ostap::Math::ECDF::ECDF ( data.begin() , data.end() , complementary )
{}
// ============================================================================
// constructor to create complementary/oridnary ECDF
// ============================================================================
Ostap::Math::ECDF::ECDF
( const Ostap::Math::ECDF & right         , 
  const bool                complementary ) 
  : ECDF ( right ) 
{
  m_complementary = complementary ; 
}
// ============================================================================
// swap two objects 
// ============================================================================
void Ostap::Math::ECDF::swap ( Ostap::Math::ECDF& right )
{
  std::swap ( m_data          , right.m_data          ) ;
  std::swap ( m_complementary , right.m_complementary ) ;
}
// ============================================================================
// the main method 
// ============================================================================
double Ostap::Math::ECDF::evaluate   ( const double x ) const
{
  if      ( x < m_data.front () ) { return m_complementary ? 1.0 : 0.0 ; } 
  else if ( x > m_data.back  () ) { return m_complementary ? 0.0 : 1.0 ; } 
  //
  const double result = double ( rank ( x ) ) / m_data.size () ;
  return m_complementary ? ( 1 - result ) : result ; 
}
// ============================================================================
// the main method 
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::ECDF::estimate ( const double x ) const
{
  const std::size_t NN = m_data.size() ;
  //
  if      ( x < m_data.front () )
    { return Ostap::Math::binomEff ( m_complementary ? NN : 0u , NN ) ; }
  else if ( x > m_data.back  () )
    { return Ostap::Math::binomEff ( m_complementary ? 0u : NN , NN ) ; }
  //
  const std::size_t success  =
    std::upper_bound ( m_data.begin () , m_data.end   () , x ) - m_data.begin() ;
  //
  return Ostap::Math::binomEff ( m_complementary ? NN - success : success , NN ) ;
}
// ============================================================================
// add a value to data container  
// ============================================================================
Ostap::Math::ECDF::Data::size_type 
Ostap::Math::ECDF::add
( const double value  )
{
  auto where = std::upper_bound ( m_data.begin () , m_data.end   () , value ) ;
  m_data.insert ( where , value ) ;
  return m_data.size () ;
}
// ============================================================================
// add a value to data container  
// ============================================================================
Ostap::Math::ECDF::Data::size_type 
Ostap::Math::ECDF::add
( const Ostap::Math::ECDF::Data& values )
{ return add ( values.begin() , values.end() ) ; }
// ============================================================================
// add values to data container  
// ============================================================================
Ostap::Math::ECDF::Data::size_type 
Ostap::Math::ECDF::add
( const Ostap::Math::ECDF& values )
{ return add_sorted ( values.m_data.begin() , values.m_data.end () ) ; }
// ============================================================================
Ostap::Math::ECDF
Ostap::Math::ECDF::__add__  ( const Ostap::Math::ECDF&  x ) const 
{
  ECDF c { *this } ;
  c.add ( x ) ;
  return c ;
} 
// ============================================================================
// get ranks of the elements in the (pooled) sample
// ============================================================================
Ostap::Math::ECDF::Indices
Ostap::Math::ECDF::ranks
( const Ostap::Math::ECDF& sample ) const
{
  const Data::size_type N  = size() ;
  // fill outptut array with N
  Indices result ( sample.size () ,  N ) ;
  Data::size_type NS= sample.size() ;
  for ( Data::size_type i = 0 ; i < NS ; ++i )
    {
      Data::size_type r = rank ( sample.m_data [ i ] ) ;
      result [ i ] = r ;
      // try to be a bit more efficient, the rest of array is filled with N 
      if ( N <= r ) { break ; }
    }
  return result ;
}
// ============================================================================
/*  assuming that x comes from the same distribution 
 *  return transformed value  \f$ g = f(x) \f$, such that 
 *  \f$ g \f$  has Gaussian distribution
 */
// ============================================================================
double Ostap::Math::ECDF::gauss   ( const double x ) const
{ return Ostap::Math::probit ( uniform ( x ) ) ; }
// ============================================================================
/*  assuming that x comes from the same distribution 
 *  return transformed value  \f$ u = f(x) \f$, such that 
 *  \f$ u \f$  has uniform distribution for \f$ 0 \le  u \le 1 \f$ 
 */
// ============================================================================
double Ostap::Math::ECDF::uniform ( const double x ) const
{
  return 
    ( ( x < xmin () ) ? (       1.0 / size () ) :
      ( x > xmax () ) ? ( 1.0 - 1.0 / size () ) : 
      ( std::lower_bound ( m_data.begin () , m_data.end () , x ) - m_data.begin() ) * 1.0 / size () ) ;
}
// ============================================================================
// For weighted data 
// ============================================================================

// ============================================================================
/* Constructor from  data
 *  data must be non-empty!
 */ 
// ============================================================================
Ostap::Math::WECDF::WECDF
( const Ostap::Math::WECDF::Data& data          ,
  const bool                      complementary )
  : m_data          ( data )
  , m_sumw          ( 0    ) 
  , m_sumw2         ( 0    ) 
  , m_complementary ( complementary ) 
{
  //
  Ostap::Assert ( 0 < m_data.size()       ,
                  "Empty data container!" ,
                  "Ostap::Math::WECDF"    ) ;
  //
  if ( !std::is_sorted ( m_data.begin() ,
			 m_data.end  () ,			 
			 []  ( const Entry& e1 ,
			       const Entry& e2 ) -> bool
			 { return e1.first < e2.first ; } ) ) { this->sort_me() ; }
  //
  m_sumw  = calc_sumw  () ;
  m_sumw2 = calc_sumw2 () ;
  //
  Ostap::Assert ( 0 < m_sumw && 0 < m_sumw2      ,
		  "Non-positive sum of weights!" ,
		  "Ostap::Math::WECDF"           ) ;
}
// ============================================================================
/*  Constructor from data
 *  data must be non-empty!
 */ 
// ============================================================================
Ostap::Math::WECDF::WECDF
( const Ostap::Math::ECDF::Data&  data          ,
  const Ostap::Math::ECDF::Data&  weights       ,
  const bool                      complementary )
  : m_data          ( data.size() )
  , m_sumw          ( 0 ) 
  , m_sumw2         ( 0 ) 
  , m_complementary ( complementary )
{
  //
  Ostap::Assert ( 0 < m_data.size()       ,
                  "Empty data container!" ,
                  "Ostap::Math::WECDF"    ) ;
  //
  const unsigned long nn = m_data .size() ;
  const unsigned long nw = weights.size() ;
  for ( unsigned int i = 0 ; i < nn ; ++i )
    {
      Entry& entry = m_data [ i ] ;
      const double weight = ( i < nw ) ? weights [ i ] : 1.0 ;
      //
      entry.first  = data [ i ]       ;
      entry.second = weight           ;
      m_sumw       += weight          ; 
      m_sumw2      += weight * weight ; 
    }
  //
  if ( !std::is_sorted ( m_data.begin() , m_data.end() ,
			 []  ( const Entry& e1 ,
			       const Entry& e2 ) -> bool
			 { return e1.first < e2.first ; } ) ) { this->sort_me() ; }
  //
  Ostap::Assert ( 0 < m_sumw && 0 <= m_sumw2     ,
                  "Non-positive sum of weights!" ,
                  "Ostap::Math::WECDF"           ) ;
}
// ============================================================================
// Standard constructor from  data
// ============================================================================
Ostap::Math::WECDF::WECDF
( const Ostap::Math::ECDF::Data&  data          ,
  const bool                      complementary )
  : Ostap::Math::WECDF::WECDF ( data                       ,
				Ostap::Math::ECDF::Data ( data.size() , 1.0 ) ,
				complementary              )    
{}
// ============================================================================
Ostap::Math::WECDF::WECDF
( const Ostap::Math::WECDF& right         ,
  const bool                complementary )
  : WECDF ( right )
{
  m_complementary = complementary ;
}
// ============================================================================
Ostap::Math::WECDF::WECDF
( const Ostap::Math::ECDF& right         ,
  const bool               complementary )
  : WECDF ( right )
{
  m_complementary = complementary ;
}
// ============================================================================
Ostap::Math::WECDF::WECDF
( const Ostap::Math::ECDF& right            )
  : m_data          ( right.size()          )
  , m_sumw          ( right.size()          ) 
  , m_sumw2         ( right.size()          ) 
  , m_complementary ( right.complementary() ) 
{
  const unsigned long nn = m_data .size() ;
  const ECDF::Data&   rd = right.data()   ;
  for ( unsigned int i = 0 ; i < nn ; ++i )
    {
      Entry& entry = m_data [ i ] ;
      entry.first  = rd [ i ]  ;
      entry.second = 1.0   ;	
    }
}
// ============================================================================
// sort container 
// ============================================================================
void Ostap::Math::WECDF::sort_me()
{
  Ostap::Assert ( 0 < m_data.size()       ,
                  "Empty data container!" ,
                  "Ostap::Math::WECDF"    ) ;  
  // sort it! 
  std::stable_sort ( m_data.begin() ,
		     m_data.end  () ,
		     []  ( const Entry& e1 ,
			   const Entry& e2 ) -> bool
		     { return e1.first < e2.first ; } ) ;  				         		    
}
// ============================================================================
// swap two objects 
// ============================================================================
void Ostap::Math::WECDF::swap ( Ostap::Math::WECDF& right )
{
  std::swap ( m_data          , right.m_data          ) ;
  std::swap ( m_sumw          , right.m_sumw          ) ;
  std::swap ( m_sumw2         , right.m_sumw2         ) ;
  std::swap ( m_complementary , right.m_complementary ) ;
}
// ============================================================================
// add a value to data container  
// ============================================================================
Ostap::Math::WECDF::Data::size_type 
Ostap::Math::WECDF::add
( const Ostap::Math::WECDF::Entry& entry ) 
{
  auto where = std::upper_bound ( m_data.begin () ,
				  m_data.end   () ,
				  entry           ,
				  []  ( const Entry& e1 ,
					const Entry& e2 ) -> bool
				  { return e1.first < e2.first ; } ) ;  				       
  m_data.insert ( where , entry ) ;
  const double weight = entry.second ;
  m_sumw  += weight ;
  m_sumw2 += weight ;
  //
  Ostap::Assert ( 0 < m_sumw && 0 <= m_sumw2     ,
                  "Non-positive sum of weights!" ,
                  "Ostap::Math::WECDF::add"      ) ;
  //
  return m_data.size () ;
}
// ============================================================================
// add values to data container  
// ============================================================================
Ostap::Math::WECDF::Data::size_type 
Ostap::Math::WECDF::add
( const Ostap::Math::WECDF& values ) 
{
  Data tmp2  ( m_data.size() + values.size() ) ;
  /// merge two sorted containers 
  std::merge ( m_data.begin        () ,
               m_data.end          () ,
               values.m_data.begin () ,
               values.m_data.end   () ,
               tmp2.begin          () , 
	       []  ( const Entry& e1 , 
		     const Entry& e2 ) -> bool
	       { return e1.first < e2.first ; } ) ;
  //
  std::swap ( m_data , tmp2 ) ;
  m_sumw  += values.m_sumw  ;
  m_sumw2 += values.m_sumw2 ;    
  return m_data.size () ;
}
// ============================================================================
// add valuee to data container  
// ============================================================================
Ostap::Math::WECDF::Data::size_type 
Ostap::Math::WECDF::add
( const Ostap::Math::ECDF& values ) 
{ return add ( WECDF ( values ) ) ; }
// ============================================================================
// add values to data container  
// ============================================================================
Ostap::Math::WECDF::Data::size_type 
Ostap::Math::WECDF::add
( const Ostap::Math::WECDF::Data& values ) 
{
  Data tmp2     ( m_data.size() + values.size() ) ;
  const Data* input = &values ;
  //
  Data values2 ;  
  if ( !std::is_sorted ( values.begin() , values.end() ,
			 []  ( const Entry& e1 ,
			       const Entry& e2 ) -> bool
			 { return e1.first < e2.first ; } ) )
    {
      values2 = values ;
      std::stable_sort ( values2.begin() ,
			 values2.end  () ,
			 []  ( const Entry& e1 ,
			       const Entry& e2 ) -> bool
			 { return e1.first < e2.first ; } ) ;  				         		    
      //
      const Data* input = &values2 ;
    }
  
  /// merge two sorted containers 
  std::merge ( m_data.begin        () ,
	       m_data.end          () ,
	       input->begin        () ,
	       input->end          () ,
	       tmp2.begin          () , 
	       []  ( const Entry& e1 , 
		     const Entry& e2 ) -> bool
	       { return e1.first < e2.first ; } ) ;
  //
  std::swap ( m_data , tmp2 ) ;
  //
  m_sumw  = calc_sumw  () ;
  m_sumw2 = calc_sumw2 () ;
  //
  Ostap::Assert ( 0 < m_sumw && 0 < m_sumw2      ,
		  "Non-positive sum of weights!" ,
		  "Ostap::Math::WECDF::add"      ) ;
  //
  return m_data.size () ;
}
// ============================================================================
// the main method 
// ============================================================================
double Ostap::Math::WECDF::evaluate   ( const double x ) const
{
  if      ( x < xmin () ) { return m_complementary ? 1.0 : 0.0 ; } 
  else if ( x > xmax () ) { return m_complementary ? 0.0 : 1.0 ; } 
  //
  const Entry entry { x , 1.0 } ;
  auto found = std::upper_bound ( m_data.begin () ,
				  m_data.end   () ,
				  entry           ,
				  []  ( const Entry& e1 , 
					const Entry& e2 ) -> bool
				  { return e1.first < e2.first ; } ) ;
  const double wsum   = calc_sumw ( found - m_data.begin() ) ;
  const double result = wsum / m_sumw ;
  //
  return m_complementary ? ( 1 - result ) : result ; 
}
// ============================================================================
// the main method 
// ============================================================================
Ostap::Math::ValueWithError
Ostap::Math::WECDF::estimate ( const double x ) const
{
  const std::size_t NN = m_data.size() ;
  //
  typedef Ostap::Math::ValueWithError  VE ;
  //
  const VE all  { m_sumw , m_sumw2 } ;
  //
  if      ( x < xmin () )
    {
      const VE none { 0 , std::pow ( m_data.front().second , 2 ) } ;
      return m_complementary ?
	Ostap::Math::binomEff2 ( all  , none ) :
	Ostap::Math::binomEff2 ( none , all  ) ; }
  else if ( x > xmax () )
    {
      const VE none { 0 , std::pow ( m_data.back().second , 2 ) } ;
      return m_complementary ?
	Ostap::Math::binomEff2 ( none , all  ) :
	Ostap::Math::binomEff2 ( all  , none ) ;
    }
  //
  const Entry entry { x , 1.0 } ;
  auto found = std::upper_bound ( m_data.begin () ,
				  m_data.end   () ,
				  entry           ,
				  []  ( const Entry& e1 , 
					const Entry& e2 ) -> bool
				  { return e1.first < e2.first ; } ) ;
  
  const double wsum   = calc_sumw  ( found - m_data.begin() ) ;
  const double w2sum  = calc_sumw2 ( found - m_data.begin() ) ;
  //  
  const VE  acc { wsum          ,           w2sum } ;
  const VE  rej { m_sumw - wsum , m_sumw2 - w2sum } ;
  //
  return m_complementary ?
    Ostap::Math::binomEff2 ( acc , rej ) : 
    Ostap::Math::binomEff2 ( rej , acc ) ;
}
// ============================================================================
// get ranks of the elements in the (pooled) sample
// ============================================================================
Ostap::Math::WECDF::Indices
Ostap::Math::WECDF::ranks
( const Ostap::Math::ECDF& sample ) const
{
  const Data::size_type N  = size() ;
  // fill outptut array with N
  Indices result ( sample.size () ,  N ) ;
  Data::size_type NS = sample.size() ;
  for ( Data::size_type i = 0 ; i < NS ; ++i )
    {
      Data::size_type r = rank ( sample.data ( i ) ) ;
      result [ i ] = r ;
      // try to be a bit more efficient, the rest of array is filled with N 
      if ( N <= r ) { break ; }
    }
  return result ;
}
// ============================================================================
// get ranks of the elements in the (pooled) sample
// ============================================================================
Ostap::Math::WECDF::Indices
Ostap::Math::WECDF::ranks
( const Ostap::Math::WECDF& sample ) const
{
  const Data::size_type N  = size() ;
  // fill output array with N
  Indices result ( sample.size () ,  N ) ;
  Data::size_type NS = sample.size() ;
  for ( Data::size_type i = 0 ; i < NS ; ++i )
    {
      Data::size_type r = rank ( sample.data ( i ) ) ;
      result [ i ] = r ;
      // try to be a bit more efficient, the rest of array is filled with N 
      if ( N <= r ) { break ; }
    }
  return result ;
}
// ============================================================================




// ============================================================================
//                                                                      The END 
// ============================================================================

