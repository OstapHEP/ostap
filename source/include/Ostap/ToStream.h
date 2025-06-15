// ============================================================================
#ifndef OSTAP_TOSTREAM_H
#define OSTAP_TOSTREAM_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <ostream>
#include <iomanip>
#include <vector>
#include <map>
#include <set>
#include <list>
#include <array>
#include <string>
#include <sstream>
#include <type_traits>
// ============================================================================
/** @file Ostap/ToStream.h
 *  implementation of various functions for streaming.
 *  this functionality is essential for usage of various types as property for
 *  the various Gaudi components
 *  @attention the implementation of the specific specializations must be done
 *                    before the inclusion of this file
 *  @todo ToStream.h : reimplement in terms of functors, to allow
 *                     easier specializations
 */
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Utils
  {
    // ========================================================================
    /** the generic implementation of the printout to the std::ostream
     *  @author Alexander MAZUROV Alexander.Mazurov@gmail.com
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-12
     */
    template<class TYPE>
    std::ostream& toStream ( const TYPE& obj, std::ostream& s ) ;
    // ========================================================================
    namespace details 
    {
      /// copied from GaudiKernel/SerializeSTL
      // ======================================================================
      struct IdentityOutputter 
      {
        template <typename T>
        std::ostream& operator() ( std::ostream& os , T&& t ) const { return os << std::forward<T> ( t ) ; }
      } ;
      // =======================================================================
      template <typename Stream, typename Iterator, typename Separator, typename OutputElement = IdentityOutputter>
      Stream& ostream_joiner
      ( Stream&       os                       ,
        Iterator      first                    ,
        Iterator      last                     ,
        Separator     sep                      ,
        OutputElement output = OutputElement{} ) 
      {
        if  (   first != last           ) { output ( os        , *first ) ; ++first ; }
        for ( ; first != last ; ++first ) { output ( os << sep , *first ) ; }
        return os;
      }
      // ======================================================================
      template <typename Stream, typename Container, typename Separator, typename OutputElement = IdentityOutputter>
      Stream& ostream_joiner(Stream& os, const Container& c, Separator sep, OutputElement output = OutputElement{}) 
      { return ostream_joiner( os, std::begin(c), std::end(c), sep, output ) ; }
      // =======================================================================     
      // helper function to print a tuple of any size
      template<class Tuple, std::size_t N>
      struct TuplePrinter 
      {
        static std::ostream& toStream ( const Tuple& t , std::ostream& s )
        {
          TuplePrinter<Tuple, N-1>::toStream ( t , s ) << " , ";
          return Ostap::Utils::toStream ( std::get<N-1>( t ) , s ) ;
        }
      };
      // =======================================================================
      template<class Tuple>
      struct TuplePrinter<Tuple, 1>
      {
        static std::ostream& toStream ( const Tuple& t, std::ostream& s)
        { return Ostap::Utils::toStream (  std::get<0> ( t ) , s ) ; }
      };
      // ======================================================================
    } //                             The end of namespace Ostap::Utils::details 
    // ========================================================================
    /** the helper function to print the sequence
     *  @param first (INPUT)  begin-iterator for the sequence
     *  @param last  (INPUT)  end-iterator for the sequence
     *  @param s     (UPDATE) the stream itself
     *  @param open  (INPUT)  "open"-symbol
     *  @param close (INPUT)  "close"-symbol
     *  @param delim (INPUT)  "delimiter"-symbol
     *  @return the stream
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2009-09-15
     */
    template <class ITERATOR>
    inline std::ostream& toStream
    ( ITERATOR           first ,                       // begin of the sequence
      ITERATOR           last  ,                       //   end of the sequence
      std::ostream&      s     ,                       //            the stream
      const std::string& open  ,                       //               opening
      const std::string& close ,                       //               closing
      const std::string& delim ) ;                     //             delimiter
    // ========================================================================
    /** the printtout of the strings.
     *  the string is printed a'la Python using the quotes
     *  @author Alexander MAZUROV Alexander.Mazurov@gmail.com
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-12
     */
    inline std::ostream& toStream
    ( const std::string& obj , std::ostream& s )
    {
      auto c = ( std::string::npos == obj.find('\'') ? '\'' : '\"' );
      return s << c << obj << c;
    }
    // ========================================================================
    /** the printout of boolean values "a'la Python"
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-09-09
     */
    inline std::ostream& toStream
    ( const bool         obj , std::ostream& s )
    { return s << ( obj ? "True" : "False" ) ; }
    // ========================================================================
    /** the printout of float values with the reasonable precision
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-09-09
     */
    inline std::ostream& toStream
    ( const float        obj , std::ostream& s , const int prec = 6 )
    {
      const int  p = s.precision() ;
      return s << std::setprecision (  prec ) << obj << std::setprecision ( p ) ;
    }
    // ========================================================================
    /** the printout of double values with the reasonable precision
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-09-09
     */
    inline std::ostream& toStream
    ( const double       obj , std::ostream& s , const int prec = 8 )
    {
      const int p = s.precision() ;
      return s << std::setprecision ( prec ) << obj << std::setprecision ( p ) ;
    }
    // ========================================================================
    /** the printout of long double values with the reasonable precision
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-09-09
     */
    inline std::ostream& toStream
    ( const long double  obj , std::ostream& s , const int prec = 10 )
    {
      const int p = s.precision() ;
      return s << std::setprecision ( prec ) << obj << std::setprecision ( p ) ;
    }
    // ========================================================================
    /** the partial template specialization of
     *  <c>std::pair<KTYPE,VTYPE></c> printout
     *  the pair is printed a'la Python tuple: " ( a , b )"
     *  @author Alexander MAZUROV Alexander.Mazurov@gmail.com
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-12
     */
    template<class KTYPE, class VTYPE>
    inline std::ostream& toStream
    ( const std::pair<KTYPE,VTYPE>& obj, std::ostream& s)
    { return toStream( obj.second, toStream( obj.first, s << "( " ) << " , "  ) << " )" ; }
    // ========================================================================
    /** the partial template specialization of <c>std::vector<TYPE,ALLOCATOR></c>
     *  printout. The vector is printed a'la Python list: "[ a, b, c ]"
     *  @author Alexander MAZUROV Alexander.Mazurov@gmail.com
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-12
     */
    template<class TYPE,class ALLOCATOR>
    inline std::ostream& toStream
    ( const std::vector<TYPE,ALLOCATOR>& obj, std::ostream& s)
    { return toStream ( obj.begin() , obj.end () , s , "[ " , " ]" , " , " ) ;  }
    // ========================================================================
    /** the partial template specialization of <c>std::list<TYPE,ALLOCATOR></c>
     *  printout. The vector is printed a'la Python list: "[ a, b, c ]"
     *  @author Alexander MAZUROV Alexander.Mazurov@gmail.com
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2007-04-08
     */
    template<class TYPE,class ALLOCATOR>
    inline std::ostream& toStream
    ( const std::list<TYPE,ALLOCATOR>& obj, std::ostream& s)
    { return toStream ( obj.begin() , obj.end () , s , "[ " , " ]" , " , " ) ; }
    // ========================================================================
    /** the partial template specialization of <c>std::set<TYPE,CMP,ALLOCATOR></c>
     *  printout. The vector is printed a'la Python list: "[ a, b, c ]"
     *  @author Alexander MAZUROV Alexander.Mazurov@gmail.com
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-12
     */
    template<class TYPE,class CMP,class ALLOCATOR>
    inline std::ostream& toStream
    ( const std::set<TYPE,CMP,ALLOCATOR>& obj, std::ostream& s)
    { return toStream ( obj.begin() , obj.end () , s , "[ " , " ]" , " , " ) ; }
    // =========================================================================
    /** the partial template specialization of <c>std::array<TYPE,N></c>
     *  printout. The vector is printed a'la Python list: "[ a, b, c ]"
     *  @author Alexander MAZUROV Alexander.Mazurov@gmail.com
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-12
     */
    template <class T, std::size_t N >
    inline std::ostream& toStream ( const std::array<T,N>& obj , std::ostream& s ) 
    { return toStream ( obj.begin() , obj.end()  , s , "( " , " )" , " , " ) ; }
    // ========================================================================
    /** the partial template specialization of
     *  <c>std::map<KTYPE,VTYPE,CMP,ALLOCATOR></c> printout
     *  the map is printed a'la Python dict: " ( a : b , c: d , e : f )"
     *  @author Alexander MAZUROV Alexander.Mazurov@gmail.com
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-12
     */
    template<class KTYPE, class VTYPE,class CMP,class ALLOCATOR>
    inline std::ostream& toStream
    ( const std::map<KTYPE,VTYPE,CMP,ALLOCATOR>& obj, std::ostream& s )
    {
      using Ostap::Utils::details::ostream_joiner;
      return ostream_joiner( s << "{ ", obj, " , ",
                             [](std::ostream& os, const std::pair<const KTYPE,VTYPE>& i)
                             -> std::ostream&
                             { return toStream( i.second, toStream( i.first, os ) << " : " ); }
                             ) << " }";
    }
    // ========================================================================
    /** the specialization for C-arrays, a'la python tuple
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2009-10-05
     */
    template <class TYPE, unsigned int N>
    std::ostream& toStream (       TYPE(&obj)[N] , std::ostream& s )
    { return toStream ( obj , obj + N , s , "( " , " )" , " , " ) ; }
    // ========================================================================
    /** the specialization for C-arrays, a'la python tuple
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2009-10-05
     */
    template <class TYPE, unsigned int N>
    std::ostream& toStream ( const TYPE(&obj)[N] , std::ostream& s )
    { return toStream ( obj , obj + N , s , "( " , " )" , " , " ) ; }
    // ========================================================================
    /** the specialization for C-string, a'la python tuple
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2009-10-05
     */
    template <unsigned int N>
    std::ostream& toStream (       char (&obj)[N] , std::ostream& s )
    { return toStream ( std::string ( obj , obj + N ) , s ) ; }
    // ========================================================================
    /** the specialization for C-string, a'la python tuple
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2009-10-05
     */
    template <unsigned int N>
    std::ostream& toStream ( const char (&obj)[N] , std::ostream& s )
    { return toStream ( std::string ( obj , obj + N ) , s ) ; }
    // ========================================================================
    /** the specialization for C-string, a'la python tuple
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2009-10-05
     */
    inline std::ostream& toStream ( const char* obj , std::ostream& s )
    { return toStream ( std::string ( obj ) , s ) ; }
    // ========================================================================    
    /** the generic implementation of the printout to the std::ostream
     *  @author Alexander MAZUROV Alexander.Mazurov@gmail.com
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-12
     */
    template<class TYPE>
    inline std::ostream& toStream ( const TYPE& obj, std::ostream& s )
    { return s << std::showpoint << std::boolalpha << obj ; }
    //
    // ========================================================================
    /** the helper function to print the sequence
     *  @param first (INPUT)  begin-iterator for the sequence
     *  @param last  (INPUT)  end-iterator for the sequence
     *  @param s     (UPDATE) the stream itself
     *  @param open  (INPUT)  "open"-symbol
     *  @param close (INPUT)  "close"-symbol
     *  @param delim (INPUT)  "delimiter"-symbol
     *  @return the stream
     *  @author Vanya BELYAEV Ivan.BElyaev@itep.ru
     *  @date 2009-09-15
     */
    template <class ITERATOR>
    inline std::ostream& toStream
    ( ITERATOR           first ,                       // begin of the sequence
      ITERATOR           last  ,                       //   end of the sequence
      std::ostream&      s     ,                       //            the stream
      const std::string& open  ,                       //               opening
      const std::string& close ,                       //               closing
      const std::string& delim )                       //             delimiter
    {
      using ref_t = typename std::iterator_traits<ITERATOR>::reference;
      using Ostap::Utils::details::ostream_joiner;
      return ostream_joiner( s << open, first, last, delim,
                             [](std::ostream& os, ref_t i ) -> std::ostream&
                             { return toStream( i, os ); } ) << close;
    }
    // ========================================================================
    /** the helper function to print the tuple
     *  @param tuple (INPUT) the tuple
     *  @param s     (UPDATE) the stream 
     *  @return the stream
     *  @author Aleander Mazurov alexander.mazurov@cern.ch
     *  @date 2015-03-21
     */
    template<typename... Args>
    inline std::ostream& toStream
    ( const std::tuple<Args...>& tuple ,
      std::ostream&              s     ) 
    { return Ostap::Utils::details::TuplePrinter<decltype(tuple), sizeof...(Args)>::toStream ( tuple , s << " ( " )<< " ) "; }
    // ========================================================================
    /// print generic pointer 
    // ========================================================================
    template <class TYPE>
    inline std::ostream& toStream
    ( const TYPE*   o ,
      std::ostream& s )
    {
      if ( nullptr == o ) { s << "nullptr" ; return o ;}
      return toStream ( *o , s ) ; 
    }
    // ========================================================================
    /** the generic implementation of the type conversion to the string
     *  @author Alexander MAZUROV Alexander.Mazurov@gmail.com
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-12
     *  @see  Ostap::Utils::toStream
     *  @todo need to be compared with boost::lexical_cast
     */
    template <class TYPE>
    inline std::string   toString ( const TYPE& obj )
    {
      std::ostringstream s;
      std::ios::fmtflags orig_flags = s.flags();
      s.setf ( std::ios::showpoint ) ; // to display correctly floats
      s.setf ( std::ios::boolalpha ) ; // for booleans
      toStream ( obj , s  ) ;
      s.flags( orig_flags ) ;
      return s.str();
    }
    // ========================================================================
    template <class TYPE>
    inline std::string to_string ( const TYPE& obj )
    { return toString ( obj ) ; }
    // ========================================================================
  } //                                            end of namespace Ostap::Utils
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END
// ============================================================================
#endif
// ============================================================================
