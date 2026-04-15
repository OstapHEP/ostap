// ============================================================================
// Include files/
// ============================================================================
// STD&STL  
// ============================================================================
#include <limits>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/ToStream.h"
// ============================================================================
/*  the printtout of the strings.
 *  the string is printed a'la Python using the quotes
 *  @author Alexander MAZUROV Alexander.Mazurov@gmail.com
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2006-05-12
 */
// ============================================================================
std::ostream& Ostap::Utils::toStream
( const std::string& obj , std::ostream& s )
{
  auto c = ( std::string::npos == obj.find('\'') ? '\'' : '\"' );
  return s << c << obj << c;
}
// ===========================================================================
/*  the printout of boolean values "a'la Python"
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2006-09-09
 */
// ===========================================================================
std::ostream& Ostap::Utils::toStream
( const bool         obj , std::ostream& s )
{ return s << ( obj ? "True" : "False" ) ; }
// ===========================================================================
/*  The printout of float values with the reasonable precision
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2006-09-09
 */
// ===========================================================================
std::ostream& Ostap::Utils::toStream
( const float          obj   ,
  std::ostream&        s     ,
  const unsigned short prec  )
{
  /// full precision 
  if ( !prec ) { return toStream ( obj , s , std::numeric_limits<float>::max_digits10 ) ; }
  const int  p = s.precision() ;
  return s << std::setprecision ( prec )
	   << std::showpos 
	   << obj
	   << std::setprecision ( p    ) ;
}
// ============================================================================
/*  The printout of double values with the reasonable precision
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2006-09-09
 */
// ============================================================================
std::ostream& Ostap::Utils::toStream
( const double         obj  ,
  std::ostream&        s    ,
  const unsigned short prec ) 
{
  /// full precision 
  if ( !prec ) { return toStream ( obj , s , std::numeric_limits<double>::max_digits10 ) ; }
  const int p = s.precision() ;
  return s << std::setprecision ( prec )
	   << std::showpos 
	   << obj
	   << std::setprecision ( p    ) ;
}
// ============================================================================
/*  The printout of long double values with the reasonable precision
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2006-09-09
 */
// ============================================================================
std::ostream& Ostap::Utils::toStream
( const long double    obj  ,
  std::ostream&        s    ,
  const unsigned short prec ) 
{
  /// full precision 
  if ( !prec ) { return toStream ( obj , s , std::numeric_limits<long double>::max_digits10 ) ; }
  const int p = s.precision() ;
  return s << std::setprecision ( prec )
	   << std::showpos 
	   << obj
	   << std::setprecision ( p    ) ;
}
// ============================================================================
/*  the printout of complex values with the reasonable precision
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2006-09-09
 */
// ============================================================================
std::ostream& Ostap::Utils::toStream
( const std::complex<float>& obj  ,
  std::ostream&              s    ,
  const unsigned short       prec ) 
{
  s << "(" ;
  toStream ( obj.real() , s , prec ) ;
  toStream ( obj.imag() , s , prec ) << "j" ;
  return s << ")" ;
}
// ============================================================================
/*  the printout of complex values with the reasonable precision
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2006-09-09
 */
// ============================================================================
std::ostream& Ostap::Utils::toStream
( const std::complex<double>& obj  ,
  std::ostream&               s    ,
  const unsigned short        prec ) 
{
  s << "(" ;
  toStream ( obj.real () , s , prec ) ; 
  toStream ( obj.imag () , s , prec ) << "j" ;
  return s << ")" ;
}    
// ============================================================================
/*  the printout of complex values with the reasonable precision
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2006-09-09
 */
// ============================================================================
std::ostream& Ostap::Utils::toStream
( const std::complex<long double>& obj  ,
  std::ostream&               s    ,
  const unsigned short        prec ) 
{
  s << "(" ;
  toStream ( obj.real () , s , prec ) ; 
  toStream ( obj.imag () , s , prec ) << "j" ;
  return s << ")" ;
}    
// ============================================================================
//                                                                      The END 
// ============================================================================
