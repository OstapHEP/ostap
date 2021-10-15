// ============================================================================
#ifndef OSTAP_INTERPOLATION_H 
#define OSTAP_INTERPOLATION_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <vector>
#include <map>
#include <array>
#include <utility>
#include <iterator>
#include <algorithm>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Math.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math  
  {
    // ========================================================================
    /** @namespace Ostap::Math::Interpolation Ostap/Interpolation.h
     *  Collection of simple utilities for various types of interpolation
     *  - straightforward Lagrange interpolation 
     *  - Neville interpolation 
     *  - Barycentric Lagrange interpolation 
     *  @see https://en.wikipedia.org/wiki/Lagrange_polynomial
     *  @see https://en.wikipedia.org/wiki/Neville%27s_algorithm
     *  @see Jean-Paul Berrut and Lloyd N. Trefethen, 
     *       Barycentric Lagrange Interpolation, SIAM Rev., 46(3), 501–517.
     *       ISSN (print): 0036-1445
     *       ISSN (online): 1095-7200
     *  @see https://doi.org/10.1137/S0036144502417715
     *  @see https://en.wikipedia.org/wiki/Lagrange_polynomial
     *  @see https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
     *  
     *  None of this method should be applies for "long" sequence of 
     *  interpolation points ( e.g. >20), especially for uniform grid 
     *  (see https://en.wikipedia.org/wiki/Runge%27s_phenomenon)
     *  
     *  - Lagrange interpolation is numerically not very stable, and rather slow O(n^2)
     *  - Neville's algorithm has (a little bit) better numerical stability and a bit faster
     *  - Barycentric Lagrange interpolation is very efficienct: 
     *    O(n) for evaluation and O(n^2) for data-independendent initialization
     *
     *  Using simple Lagrange algorithm it is easy to get derivative with 
     *  respect to the data points, while using Neville's algorithm 
     *  one can easily calcualate the derivative with respect to 
     *  the argument.  
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2016-07-23
     */
    namespace Interpolation
    {
      // ======================================================================
      /// the actual type of "simple" data 
      typedef std::vector<std::pair<double,double> >           DATA    ;
      /// the actual type of "simple" data 
      typedef DATA                                             TABLE   ;
      /// the actual type of "simple" data 
      typedef std::vector<double>                              DATAVCT ;
      // ======================================================================
      /** @class Abscissas 
       *  Collection of interpolation abscissas 
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2016-07-23
       */
      class Abscissas
      {
        // ====================================================================
      public:
        // ====================================================================
        /// the actual vector of data 
        typedef std::vector<double>  Data ; // the actual vector of data 
        // ====================================================================
        enum Type { Uniform          = 0          , 
                    Chebyshev        = 1          , //   roots of T_n     ( x ) 
                    Chebyshev1       = Chebyshev  ,
                    GaussChebyshev   = Chebyshev  , 
                    Chebyshev2       = 2          , // extrema of T_{n-1} ( x ) 
                    Lobatto          = Chebyshev2 , 
                    ChebyshevLobatto = Chebyshev2 ,
                    GaussLobatto     = Chebyshev2 } ;
        // ====================================================================
      public:
        // ====================================================================
        /** create the abscissas from vector of abscissas 
         *  @param x       input vector of abscissas (to be sorted internally)
         *  @param sorted  indicate if input data is already sorted 
         *  @attention duplicated abscissas will be removed 
         */
        Abscissas ( const Data&         x     , 
                    const bool sorted = false ) ;
        // ===================================================================
        /** templated constructor from the arbitrary sequence
         *  @param begin start of the sequnce 
         *  @param end   end   of the sequnce 
         *  @param sorted  indicate if input data is already sorted 
         *  @attention duplicated abscissas will be removed 
         */
        template <class ITERATOR>
        Abscissas ( ITERATOR begin , 
                    ITERATOR end   , 
                    const bool sorted = false ) 
          :  m_x ( begin , end ) 
        { if ( !sorted ) { this->sort () ; } }
        // =================================================================
        /** templated constructor from the arbitrary sequence
         *  @param begin start of the sequnce 
         *  @param end   end   of the sequnce 
         *  @param fun   the function to be applied 
         *  @param sorted  indicate if input data is already sorted 
         *  @attention duplicated abscissas will be removed 
         */
        template <class ITERATOR, class FUNCTION>
        Abscissas ( ITERATOR   begin          , 
                    ITERATOR   end            , 
                    FUNCTION   fun            , 
                    const bool sorted = false ) 
          : m_x ( std::distance ( begin, end ) , 0.0 )  
        {
          std::transform ( begin , end , m_x.begin() , fun );
          if ( !sorted ) { this->sort () ; } 
        }
        // ====================================================================
        /** special constructor for the given interpolation type 
         *  @param n    number of interpolation points 
         *  @param low  low edge of the interval 
         *  @param high high edge of the interval 
         *  @param t    interpolation type 
         *  @see Ostap::Math::Interpolation::Abscissas::Type 
         */
        Abscissas ( const unsigned short n                ,
                    const double         low              ,   
                    const double         high             , 
                    const Type           t    = Chebyshev ) ;
        // ====================================================================
        /// default constructor  
        Abscissas ()                    = default ;
        // ====================================================================
      public:
        // ====================================================================
        unsigned int n       () const { return m_x.size  () ; }
        unsigned int size    () const { return m_x.size  () ; }
        bool         empty   () const { return m_x.empty () ; }
        // ====================================================================
      public:
        // ====================================================================
        /// get all abscissas 
        const Data&  x       () const { return m_x ; }
        /// get all abscissas (as cast-operator)
        operator const Data& () const { return m_x ; }
        // ====================================================================
        ///  get the absissas for the given index 
        double x    ( const unsigned short index ) const 
        { return index < m_x.size() ? m_x [ index ] : m_x.back() ; }
        // ====================================================================
      public:
        // ====================================================================
        /// minimal  abscissa 
        double xmin () const { return  m_x.front () ; }
        /// maximal abscissas 
        double xmax () const { return  m_x.back  () ; }
        // ====================================================================
      public: // expose (const) iterators 
        // ====================================================================
        Data::const_iterator begin () const { return  m_x.begin () ; }
        Data::const_iterator end   () const { return  m_x.end   () ; }        
        // ====================================================================
      public: // add & remove the point 
        // ====================================================================
        /** add one more abscissa in the list 
         *  @param xnew value to be added 
         *  @return -1 if point is NOT added or new index for added points  
         */
        int add     ( const double         xnew ) ;
        /** remove the point with the given index 
         *  @param index poitn with the index 
         *  @return true if point is really removed 
         */
        bool remove ( const unsigned short index ) ;
        // ====================================================================
      public: // make a slice fo the given  ascissas 
        // ====================================================================
        /// make a slice fo the given  ascissas 
        Abscissas slice ( const int i , const int j ) const ;
        // ====================================================================
      public:
        // ====================================================================
        friend void swap ( Abscissas&  a , Abscissas& b ) { swap (  a.m_x , b.m_x ) ; }
        // ====================================================================
      private:
        // ====================================================================
        /** sort abscissas and eliminate the duplicates  
         *  @return number of removed duplicates  
         */
        unsigned int sort ( ) ;
        // ====================================================================
      private:
        // ====================================================================
        /// (sorted) vector of abscissas 
        Data m_x {} ; // (sorted) vector of abscissas 
        // ====================================================================
      };
      // ======================================================================
      /** @class Table
       *  Interpolation table 
       *  it contains the ordered list of pairs (abscissa,value) 
       *  - duplicated entries are removed 
       */
      class Table 
      {
      public:
        // ====================================================================
        /** the simplest constructor 
         *  @param data   input data 
         *  @param sorted indicate if data already  sorted and duplicated removed  
         */
        Table ( const TABLE& data           , 
                const bool   sorted = false ) ;
        // ====================================================================
        /** simple (efficient) constructor from abscissas and y-list 
         *  @param x input vector of abscissas  
         *  @param y input vector of y 
         *  - if vector of y is longer  than vector x, extra values are ignored 
         *  - if vector of y is shorter than vector x, missing entries are assumed to be zero  
         */
        Table ( const Abscissas&       x , 
                const Abscissas::Data& y ) ;
        // ====================================================================
        /** simple contructor from x&y-lists 
         *  @param x input vector of x 
         *  @param y input vector of y 
         *  - if vector of y is longer  than vector x, extra values are ignored 
         *  - if vector of y is shorter than vector x, missing entries are assumed to be zero  
         *  @attention duplicated abscissas will be removed 
         */
        Table ( const Abscissas::Data& x , 
                const Abscissas::Data& y ) ;
        // ====================================================================
        /** simple constructor from map/dictionary { x : f(x) }
         *  @param data map/dictionary with data  { x : f(x) }
         *  It is elatively efficient: no sorting
         *  - "numerical" duplicates are removed
         */
        template <class KEY, class VALUE>
        Table ( const std::map<KEY,VALUE>& data ) 
          : m_table ( data.size() )  
        {
          TABLE::iterator i = m_table.begin() ;
          for ( const auto& item : data ) 
          { *i = std::make_pair ( item.first , item.second ) ; ++i ; }
          /// eliminate duplicates, no need to sort it  
          this->get_sorted ( true ) ;
        }
        // ====================================================================
        /** simple contructor from x&y-lists 
         *  @param xbegin start  of sequence of abscissas 
         *  @param xend   end    of sequence of abscissas 
         *  @param ybegin start  of sequence of values 
         *  @param yend   end    of sequence of values 
         *  @param sorted indicate if data already sorted and duplicates removed  
         *  - if vector of y is longer  than vector x, extra values are ignored 
         *  - if vector of y is shorter than vector x, missing entries are assumed to be zero  
         *  @attention duplicated abscissas will be removed 
         */
        template <class XITERATOR, class YITERATOR>
        Table ( XITERATOR  xbegin         ,
                XITERATOR  xend           ,
                YITERATOR  ybegin         ,
                YITERATOR  yend           , 
                const bool sorted = false )
          : m_table ( std::distance ( xbegin , xend ) ) 
        {  
          // ==================================================================
          /// 1) copy x into the table 
          std::transform 
            ( xbegin          , xend , 
              m_table.begin() , 
              [] ( double v ) { return std::make_pair ( v , 0.0 ) ; } );
          /// 2) copy y into the table 
          const unsigned int N = std::min ( std::distance ( xbegin , xend ) , 
                                            std::distance ( ybegin , yend ) ) ;
          std::transform
            ( m_table.begin () , m_table.begin () + N , 
              ybegin           , 
              m_table.begin () , 
              [] ( const auto& p , const double v ) 
              { return std::make_pair ( p.first , v ) ; } );
          /// 3) sort it & eliminate  duplicates 
          this->get_sorted() ;  
          // ==================================================================
        }
        // ====================================================================
        /** templated constructor 
         *  (very efficient: no sorting, no removal of duplicate...)
         *  @param fun       function object 
         *  @param abscissas interpolation abscissas 
         */
        template <class FUNCTION>
        Table ( const Abscissas& a   , 
                FUNCTION         fun )                
          : m_table ( a.size() ) 
        {
          // ==================================================================
          std::transform 
            ( a.begin () , a.end () ,  
              m_table.begin ()      , 
              [&fun]  ( const double x ) 
              { return std::make_pair ( x , fun ( x ) ) ; } ) ;
          // ==================================================================
        }     
        // ====================================================================
        /** templated constructor 
         *  @param begin start  of x-sequence 
         *  @param end   end    of x-sequence 
         *  @param fun   function object 
         *  @attention duplicated abscissas will be removed 
         */
        template <class ITERATOR, class FUNCTION>
        Table  ( FUNCTION fun   , 
                 ITERATOR begin , 
                 ITERATOR end   ) 
          : Table ( Abscissas ( begin, end ) , fun ) {}
        // ====================================================================
        /** templated constructor 
         *  @param fun  function object 
         *  @param n    number of interpolation points 
         *  @param low  low edge of interpolation region 
         *  @param high high edge of interpolation region 
         *  @param t    interpolation type 
         */
        template <class FUNCTION>
        Table ( FUNCTION              fun  , 
                const unsigned short  n    , 
                const double          low  ,   
                const double          high , 
                const Abscissas::Type t    ) 
          : Table ( Abscissas ( n , low , high , t ) , fun ) {}
        // ====================================================================
        /// default constructor
        Table () = default ;
        // ====================================================================
      public:
        // ====================================================================
        unsigned int n       () const { return m_table.size  () ; }
        unsigned int size    () const { return m_table.size  () ; }
        bool         empty   () const { return m_table.empty () ; }
        // ====================================================================
      public : // 
        // ====================================================================
        /// get the abscissas 
        Abscissas abscissas () const ;
        // ====================================================================
      public:
        // ====================================================================
        ///  get the absissas for the given index 
        double x    ( const unsigned short index ) const 
        { return index < m_table.size() ? m_table[ index ].first  : m_table.back().first  ; }
        ///  get the value    for the given index 
        double y    ( const unsigned short index ) const 
        { return index < m_table.size() ? m_table[ index ].second : m_table.back().second ; }
        // =====================================================================
      public:
        // ====================================================================
        /// minimal  abscissa 
        double xmin () const { return  m_table.front ().first ; }
        /// maximal abscissas 
        double xmax () const { return  m_table.back  ().first ; }
        // ====================================================================
      public: // expose begin/end iterators 
        // ====================================================================
        /// begin of the table 
        TABLE::const_iterator begin () const { return m_table.begin () ; }
        /// end of the table 
        TABLE::const_iterator end   () const { return m_table.end   () ; }
        // =====================================================================
      public: // show internal data 
        // =====================================================================
        const TABLE& table () const { return m_table ; }
        const TABLE& data  () const { return m_table ; }
        // =====================================================================
      public:
        // ====================================================================
        /** add the point (x,y) into interpolation table 
         *  @param x abscissas of the point to be added 
         *  @param y the value of function at x 
         *  @return the index of new point, or -1  if point if not added 
         */
        int add     ( const  double x ,   const  double y ) ;
        // ====================================================================
        /** remove the point with the  given index 
         *  @param index the point to be removed 
         *  @return   true if point is removed 
         */
        bool remove ( unsigned short index ) ;
        // ====================================================================
      public: // make interpolation 
        // ====================================================================
        /** interpolation using the straightforward Lagrange interpolant 
         *  https://en.wikipedia.org/wiki/Lagrange_polynomial
         *  - it is rather slow O(n^2)
         */
        double lagrange ( const double x ) const ;
        // ====================================================================
        /** interpolation using Neville's algorithm
         *  @see https://en.wikipedia.org/wiki/Neville%27s_algorithm
         *  - it is rather slow O(n^2)
         */
        double neville  ( const double x ) const ;
        // ====================================================================
      public: // interpolation with derivatives:
        // ====================================================================
        /** Interpolation using Neville's algorithm
         *  @see https://en.wikipedia.org/wiki/Neville%27s_algorithm
         *  - it is rather slow O(n^2)
         *  @return ( y(x) , dy/dx )
         */
        std::pair<double,double> 
        neville2  
        ( const double       x ) const ;
        // ====================================================================
        /** Simple lagrange interpolation 
         *  - it also evaluate the derivative wity respect to y_i 
         *  @param x interpolation point
         *  @param iy index of y_i
         *  @return ( y(x) , dy/d(y_i))
         */
        std::pair<double,double>
        lagrange2
        ( const double       x  , 
          const unsigned int iy ) const ;         
        // ====================================================================
      public: // make a slice for the given ascissas 
        // ====================================================================
        /// make a slice for the given range of points 
        Table slice ( const int i , const int j ) const ;
        // ====================================================================
      public:
        // ====================================================================
        /// swap two collections of points 
        friend void swap ( Table& a , Table& b ) 
        { swap ( a.m_table , b.m_table ) ; }
        // ====================================================================
      private:
        // ====================================================================
        ///  sort internal data and remove duplicated abscissas 
        void get_sorted ( bool sorted = false ) ;
        // ====================================================================
      private :
        // ====================================================================
        /// the actual table with data 
        TABLE m_table ;
        // ====================================================================
      } ;  
      // ======================================================================
      // The basic interpolation functions 
      // ======================================================================
      /** Very basic simple lagrange interpolation 
       *
       *  @param  xbegin INPUT start iterator for the sequence of abscissas 
       *  @param  xend   INPUT end   iterator for the sequence of abscissas 
       *  @param  ybegin INPUT start iterator for the sequence of values 
       *  @param  yend   INPUT end   iterator for the sequence of values 
       *  @param  x      INPUT evaluate the polynomial in this point
       *  @param  result result 
       *  @param  xvalue INPUT adapter for x-values  
       *  @param  yvalue INPUT adapter for y-values  
       *  @return the value of Largange interpolation polynomial at point x
       *
       *  - If y-vector is shorter than x-vector, it is assumed to be zero padded
       *  - If y-vector is longer  than x-vector, extra avalues are ignored       
       *  - If x-vector is empty, polynomial is zero
       *
       *  @see https://en.wikipedia.org/wiki/Lagrange_polynomial
       *  @warning it could be CPU inefficient
       *  @warning it should *NOT* be applied for very long sequence of points 
       *           (e.g. >20) due to bad numerical  stability and Runge phenomenon
       *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
       *  @date 2016-07-23
       */
      template <class XITERATOR, 
                class YITERATOR,
                class RESULT   , 
                class XADAPTER ,
                class YADAPTER >
      inline RESULT
      lagrange
      ( XITERATOR    xbegin , 
        XITERATOR    xend   , 
        YITERATOR    ybegin , 
        YITERATOR    yend   , 
        const double x      , 
        RESULT       result , 
        XADAPTER     xvalue ,   // adaptor to get y-value
        YADAPTER     yvalue ) ; // adaptor to get x-value
      // ======================================================================
      /** simple interpolation using Neville's algorithm 
       *
       *  In general it should be faster than largange algorithm, 
       *  but it includes a copy on input data, that could affect CPU performance
       *  Numerically it is more stable that Lagrange interpolation, 
       *  but anyhow it also should not be used for very high polynomial 
       *  degrees (e.g.>20), especially for uniform grid (Runge phenomenon)
       *
       *  - If y-vector is shorter than x-vector, it is assumed to be zero padded
       *  - If y-vector is longer  than x-vector, extra avalues are ignored       
       *  - If x-vector is empty, polynomial is zero
       *
       *  @param  xbegin INPUT start iterator for the sequence of abscissas 
       *  @param  xend   INPUT end   iterator for the sequence of abscissas 
       *  @param  ybegin INPUT start iterator for the sequence of values 
       *  @param  yend   INPUT end   iterator for the sequence of values 
       *  @param  x      INPUT evaluate the polynomial in this point
       *  @param  xvalue INPUT adapter for x-values  
       *  @param  yvalue INPUT adapter for y-values  
       *  @return the value of interpolation polynomial at point x
       *  @see https://en.wikipedia.org/wiki/Neville%27s_algorithm
       */
      template <class XITERATOR, 
                class YITERATOR, 
                class XADAPTER ,
                class YADAPTER >
      inline double
      neville
      ( XITERATOR    xbegin ,
        XITERATOR    xend   ,
        YITERATOR    ybegin ,
        YITERATOR    yend   ,
        const double x      , 
        XADAPTER     xvalue ,   // adaptor to get y-value
        YADAPTER     yvalue ) ; // adaptor to get x-value
      // ======================================================================
      /** simple interpolation using Neville's algorithm: 
       *    evaluate the interpolation polynomial and also the derivative 
       *
       *  In general it should be faster than largange algorithm, 
       *  but it includes a copy on input data, that could affect CPU performance
       *  Numerically it is more stable that Lagrange interpolation, 
       *  but anyhow it also should not be used for very high polynomial 
       *  degrees (e.g.>20), especially for uniform grid (Runge phenomenon)
       *
       *  - If y-vector is shorter than x-vector, it is assumed to be zero padded
       *  - If y-vector is longer  than x-vector, extra avalues are ignored       
       *  - If x-vector is empty, polynomial is zero
       *
       *  @param  xbegin INPUT start iterator for the sequence of abscissas 
       *  @param  xend   INPUT end   iterator for the sequence of abscissas 
       *  @param  ybegin INPUT start iterator for the sequence of values 
       *  @param  yend   INPUT end   iterator for the sequence of values 
       *  @param  x      INPUT evaluate the polynomial in this point
       *  @param  xvalue INPUT adapter for x-values  
       *  @param  yvalue INPUT adapter for y-values  
       *  @return the value of interpolation polynomial at point x
       *  @see https://en.wikipedia.org/wiki/Neville%27s_algorithm
       */
      template <class XITERATOR, 
                class YITERATOR, 
                class XADAPTER ,
                class YADAPTER >
      inline 
      std::pair<double,double> 
      neville2
      ( XITERATOR    xbegin ,
        XITERATOR    xend   ,
        YITERATOR    ybegin ,
        YITERATOR    yend   ,
        const double x      , 
        XADAPTER     xvalue ,   // adaptor to get y-value
        YADAPTER     yvalue ) ; // adaptor to get x-value
      // ======================================================================
      /** simple interpolation using Neville's algorithm 
       *
       *  In general it should be faster than largange algorithm, 
       *  but it modified input data!
       *
       *  Numerically it is more stable than Lagrange interpolation, 
       *  but anyhow it also should not be used for very high polynomial 
       *  degrees (e.g.>20), especially for uniform grid (Runge phenomenon)
       *
       *  @attention y-sequence must be "simple", convertible to doubles  
       *  @attention y-sequence is *MODIFIED*
       *
       *  @param  xbegin INPUT  start iterator for the sequence of abscissas 
       *  @param  xend   INPUT  end   iterator for the sequence of abscissas 
       *  @param  ybegin UPDATE start iterator for the sequence of values 
       *  @param  x      INPUT  evaluate the polynomial in this point
       *  @param  xvalue INPUT  adapter for x-values  
       *  @return the value of  interpolation polynomial at point x
       *  @see https://en.wikipedia.org/wiki/Neville%27s_algorithm
       */
      template <class XITERATOR, 
                class YITERATOR, 
                class XADAPTOR >
      inline double 
      neville
      ( XITERATOR    xbegin , 
        XITERATOR    xend   , 
        YITERATOR    ybegin , // NON-const!
        const double x      , 
        XADAPTOR     xvalue ) ;
      // ======================================================================      
      /** simple interpolation using Neville's algorithm with simultaneous 
       *  estimation of the derivative  
       *
       *  Numerically it is more stable than Lagrange interpolation, 
       *  but anyhow it also should not be used for very high polynomial 
       *  degrees (e.g.>20), especially for uniform grid (Runge phenomenon)
       *
       *  @attention y-sequence must be "simple", convertible to doubles  
       *  @attention y-sequence is *MODIFIED*
       *  @attention d-sequence is *MODIFIED*
       *
       *  @param  xbegin INPUT  start iterator for the sequence of abscissas 
       *  @param  xend   INPUT  end   iterator for the sequence of abscissas 
       *  @param  ybegin UPDATE start iterator for the sequence of values 
       *  @param  dbegin UPDATE start iterator
       *  @param  x      INPUT  evaluate the polynomial in this point
       *  @param  xvalue INPUT  adapter for x-values  
       *  @return the pair (function,derivative) at point x 
       *  @see https://en.wikipedia.org/wiki/Neville%27s_algorithm
       */
      template <class XITERATOR, 
                class YITERATOR, 
                class DITERATOR, 
                class XADAPTOR >
      inline std::pair<double,double>
      neville2
      ( XITERATOR    xbegin , 
        XITERATOR    xend   , 
        YITERATOR    ybegin , // NON-const!
        DITERATOR    dbegin , // NON-const!
        const double x      , 
        XADAPTOR     xvalue ) ;
      // ======================================================================      
      /** very simple lagrange interpolation 
       *  @param  xs INPUT sequence of abscissas 
       *  @param  ys INPUT sequence of values 
       *  @param  x      INPUT evaluate the polynomial in this point
       *  @return the value of Largange interpolation polynomial at point x
       *
       *  - If ys-vector is shorter than x-vector, it is assumed to be zero padded
       *  - If ys-vector is longer  than x-vector, extra avalues are ignored       
       *  - If xs-vector is empty, polynomial is zero
       *
       *  @warning it could be CPU inefficient
       *  @warning it should *NOT* be applied for very long sequence of points 
       *           (e.g. >20) due to bad numerical  stability
       *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
       *  @date 2016-07-23
       */
      double
      lagrange 
      ( const std::vector<double>& xs , 
        const std::vector<double>& ys , 
        const double               x  ) ;
      // ======================================================================      
      /** very simple lagrange interpolation 
       *  @param  data INPUT sequence of (x,y)
       *  @param  x      INPUT evaluate the polynomial in this point
       *  @return the value of Largange interpolation polynomial at point x
       *
       *  - If data-vector is empty, polynomial is zero
       *
       *  @warning it could be CPU inefficient
       *  @warning it should *NOT* be applied for very long sequence of points 
       *           (e.g. >20) due to bad numerical  stability
       *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
       *  @date 2016-07-23
       */
      double 
      lagrange
      ( const DATA&  data , 
        const double x    ) ;
      // ======================================================================      
      /** Simple lagrange interpolation 
       *  - it also evaluate the derivative wity respect to y_i 
       *  @param  xs INPUT sequence of abscissas 
       *  @param  ys INPUT sequence of values 
       *  @param  x  INPUT evaluate the polynomial in this point
       *  @param  iy INPUT index of y_i, the derivative shodul be calculated.
       *  @return the value of Largange interpolation polynomial at point x
       *
       *  - If ys-vector is shorter than x-vector, it is assumed to be zero padded
       *  - If ys-vector is longer  than x-vector, extra avalues are ignored       
       *  - If xs-vector is empty, polynomial is zero
       *
       *  @warning it could be CPU inefficient
       *  @warning it should *NOT* be applied for very long sequence of points 
       *           (e.g. >20) due to bad numerical  stability
       *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
       *  @date 2016-07-23
       */
      std::pair<double,double>
      lagrange2 
      ( const std::vector<double>& xs , 
        const std::vector<double>& ys , 
        const double               x  , 
        const unsigned int         iy ) ;
      // ======================================================================      
      /** Simple lagrange interpolation 
       *  - it also evaluate the derivative with respect to y_i 
       *  @param  data INPUT sequence of (x,y)
       *  @param  x    INPUT evaluate the polynomial in this point
       *  @param  iy   INPUT index of y_i, the derivative shodul be calculated.
       *  @return the value of Largange interpolation polynomial at point x
       *
       *  - If data-vector is empty, polynomial is zero
       *
       *  @warning it could be CPU inefficient
       *  @warning it should *NOT* be applied for very long sequence of points 
       *           (e.g. >20) due to bad numerical  stability
       *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
       *  @date 2016-07-23
       */
      std::pair<double,double>
      lagrange2
      ( const DATA&        data , 
        const double       x    , 
        const unsigned int iy   ) ;
      // ======================================================================      
      /** very simple Neville's interpolation 
       *  @param  xs INPUT sequence of abscissas 
       *  @param  ys INPUT sequence of values 
       *  @param  x  INPUT evaluate the polynomial in this point
       *  @return the value of Largange interpolation polynomial at point x
       *
       *  - If ys-vector is shorter than x-vector, it is assumed to be zero padded
       *  - If ys-vector is longer  than x-vector, extra avalues are ignored       
       *  - If xs-vector is empty, polynomial is zero
       *
       *  @warning it could be CPU inefficient
       *  @warning it should *NOT* be applied for very long sequence of points 
       *           (e.g. >20) due to bad numerical  stability
       *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
       *  @date 2016-07-23
       */
      double 
      neville
      ( const std::vector<double>& xs ,
        const std::vector<double>& ys , 
        const double               x  ) ;
      // ======================================================================      
      /** very simple Neville interpolation 
       *  @param  data INPUT sequence of (x,y)
       *  @param  x    INPUT evaluate the polynomial in this point
       *  @return the value of Largange interpolation polynomial at point x
       *
       *  - If data-vector is empty, polynomial is zero
       *
       *  @warning it could be CPU inefficient
       *  @warning it should *NOT* be applied for very long sequence of points 
       *           (e.g. >20) due to bad numerical  stability
       *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
       *  @date 2016-07-23
       */
      double 
      neville
      ( const DATA&  data , 
        const double x    ) ;
      // ======================================================================      
      /** very simple Neville's interpolation 
       *  -  it evalutes the polynomial and the derivative
       *  @param  xs INPUT sequence of abscissas 
       *  @param  ys INPUT sequence of values 
       *  @param  x  INPUT evaluate the polynomial in this point
       *  @return the value of Largange interpolation polynomial at point x
       *
       *  - If ys-vector is shorter than x-vector, it is assumed to be zero padded
       *  - If ys-vector is longer  than x-vector, extra avalues are ignored       
       *  - If xs-vector is empty, polynomial is zero
       *
       *  @warning it could be CPU inefficient
       *  @warning it should *NOT* be applied for very long sequence of points 
       *           (e.g. >20) due to bad numerical  stability
       *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
       *  @date 2016-07-23
       */
      std::pair<double,double> 
      neville2
      ( const std::vector<double>& xs ,
        const std::vector<double>& ys , 
        const double               x  ) ;
      // ======================================================================      
      /** very simple Neville interpolation 
       *  -  it evaluates the polynomial and the derivative
       *  @param  data INPUT sequence of (x,y)
       *  @param  x    INPUT evaluate the polynomial in this point
       *  @return the value of Largange interpolation polynomial at point x
       *
       *  - If data-vector is empty, polynomial is zero
       *
       *  @warning it could be CPU inefficient
       *  @warning it should *NOT* be applied for very long sequence of points 
       *           (e.g. >20) due to bad numerical  stability
       *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
       *  @date 2016-07-23
       */
      std::pair<double,double> 
      neville2 
      ( const DATA&  data , 
        const double x  ) ;
      // ======================================================================      
      /** @class Weights
       *  helper class to keep barycentric weigths for Lagrange interpolation
       *  @see Jean-Paul Berrut and Lloyd N. Trefethen, 
       *       Barycentric Lagrange Interpolation, SIAM Rev., 46(3), 501–517.
       *       ISSN (print): 0036-1445
       *       ISSN (online): 1095-7200
       *  @see https://doi.org/10.1137/S0036144502417715
       *  @see https://en.wikipedia.org/wiki/Lagrange_polynomial
       *  @see https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
       *  Barycentric weigths are calcualted for 
       *  - O(n) for <code>Chebyshev</code> and <code>Chebvyshev2</code>  abscissas 
       *  - O(n) for <code>Uniform</code> abscissas, but much slower 
       *  - O(n^2) for general case 
       *  @see Ostap::Math::Interpolation::Abscissas
       *  @see Ostap::Math::Interpolation::Abscissas::Type
       */
      class Weights
      {
      public:
        // ====================================================================
        ///  the actual type for vector of weights 
        typedef Abscissas::Data                                         Data  ;
        // ==================================================================== 
      public:
        // ====================================================================
        /// constructor from abscissas 
        Weights ( const Abscissas& a       ) ;
        /// constructor from abscissas 
        Weights ( const Abscissas::Data& a ) ; 
        /// constructor from abscissas
        Weights ( const unsigned short  n    ,
                  const double          xmin , 
                  const double          xmax , 
                  const Abscissas::Type t    ) ;
        // ====================================================================
        template <class ITERATOR>
        Weights 
        ( ITERATOR   begin          , 
          ITERATOR   end            ,
          const bool sorted = false )
          : m_a ( begin , end , sorted ) 
          , m_w () 
        { 
          m_w.resize( m_a.size() ) ;
          this->get_weights() ; 
        }
        // ====================================================================
        template <class ITERATOR, class  FUNCTION>
        Weights 
        ( ITERATOR   begin          , 
          ITERATOR   end            ,
          FUNCTION   fun            ,
          const bool sorted = false )
          : m_a ( begin , end , fun , sorted ) 
          , m_w () 
        { 
          m_w.resize( m_a.size() ) ;
          this->get_weights() ; 
        }
        // ====================================================================
        /// default constructor 
        Weights ( ) = default ;
        // ====================================================================
      public :
        // ====================================================================
        unsigned int n       () const { return m_a.size  () ; }
        unsigned int size    () const { return m_a.size  () ; }
        bool         empty   () const { return m_a.empty () ; }
        // ======================================================================
      public:
        // ====================================================================
        /// get abscissas
        const Abscissas&       abscissas () const { return m_a     ; }
        const Abscissas&       a         () const { return m_a     ; }
        const Abscissas::Data& x         () const { return m_a.x() ; }
        /// get weights         
        const Data&            weights   () const { return m_w     ; }
        const Data&            w         () const { return m_w     ; }
        /// get abscissa:
        double abscissa ( unsigned short index ) const { return m_a.x ( index ) ; }
        double x        ( unsigned short index ) const { return m_a.x ( index ) ; }
        /// get weight: 
        double weight   ( unsigned short index ) const 
        { return index < m_w.size() ? m_w[index] : m_w.back() ; }
        /// get weight: 
        double w        ( unsigned short index ) const { return weight ( index ) ; }
        // ====================================================================
      public:
        // ====================================================================
        /// minimal  abscissa 
        double xmin () const { return  m_a.xmin () ; }
        /// maximal abscissas 
        double xmax () const { return  m_a.xmax () ; }
        // ====================================================================
      public:
        // ====================================================================
        /** add the point x into  collection 
         *  @param x abscissas of the point to be added 
         *  @return the index of new point, or -1  if point if not added 
         */
        int add     ( const  double x ) ;
        // ====================================================================
        /** remove the point with the  given index 
         *  @param index the point to be removed 
         *  @return   true if point is removed 
         */
        bool remove ( unsigned short index ) ;
        // ====================================================================
      public: 
        // ====================================================================
        /// get the slice 
        Weights slice ( const int i , const int j ) const ;
        // ====================================================================
      public: 
        // ====================================================================
        /// swap two objects 
        friend void swap ( Weights& a , Weights& b ) 
        { swap ( a.m_a , b.m_a ) ; swap ( a.m_w , b.m_w ) ; }
        // ====================================================================
      private: 
        // ====================================================================
        /// calculate the barycentric weights 
        void get_weights () ; // calculate the barycentric weights 
        // ====================================================================
      private: 
        // ====================================================================
        /// the abscissas 
        Abscissas m_a {} ; // the abscissas        
        ///  the weights 
        Data      m_w {} ; // the weights 
        // ====================================================================        
      };
      // ======================================================================      
    } //                            end of namespace Ostap::Math::Interpolation 
    // ========================================================================
    /** @class Neville 
     *  Simple interpolation polynomial using Neville's algorithm 
     *  @see https://en.wikipedia.org/wiki/Neville%27s_algorithm
     *  @attention  it is not CPU efficient!
     */
    class Neville : public Interpolation::Table
    {
    public:
      // ====================================================================== 
      /** simple contructor from interpolation points  
       *  @param x input vector of abscissas  
       */
      Neville ( const Interpolation::Table& p ) ;
      // ====================================================================
      // other constructors from the base class 
      // ====================================================================
      /** the simplest constructor 
       *  @param data   input data 
       *  @param sorted indicate if data already  sorted and duplicated removed  
       */
      Neville ( const Interpolation::TABLE& data           , 
                const bool                  sorted = false ) 
        : Interpolation::Table ( data , sorted ) {}
      // ====================================================================== 
      /** simple contructor from abscissas and y-list 
       *  @param x input vector of abscissas  
       *  @param y input vector of y 
       *  - if vector of y is longer  than vector x, extra values are ignored 
       *  - if vector of y is shorter than vector x, missing entries are assumed to be zero  
       */
      Neville ( const Interpolation::Abscissas&       x , 
                const Interpolation::Abscissas::Data& y ) 
        : Interpolation::Table ( x , y ){}
      // ===================================================================
      /** simple contructor from x&y-lists 
       *  @param x input vector of x 
       *  @param y input vector of y 
       *  - if vector of y is longer  than vector x, extra values are ignored 
       *  - if vector of y is shorter than vector x, missing entries are assumed to be zero  
       *  @attention duplicated abscissas will be removed 
       */
      Neville ( const Interpolation::Abscissas::Data& x , 
                const Interpolation::Abscissas::Data& y ) 
        : Interpolation::Table ( x , y ) {}
      // ====================================================================
      /** simple constructor from map/dictionary { x : f(x) }
       *  @param data map/dictionary with data  { x : f(x) }
       *  It is elatively efficient: no sorting
       *  - "numerical" duplicates are removed
       */
      template <class KEY, class VALUE>
      Neville ( const std::map<KEY,VALUE>& data ) 
        : Interpolation::Table ( data ) {}
      // ====================================================================
      /** simple contructor from x&y-lists 
       *  @param xbegin start  of sequence of abscissas 
       *  @param xend   end    of sequence of abscissas 
       *  @param ybegin start  of sequence of values 
       *  @param yend   end    of sequence of values 
       *  @param sorted indicate if data already sorted and duplicates removed  
       *  - if vector of y is longer  than vector x, extra values are ignored 
       *  - if vector of y is shorter than vector x, missing entries are assumed to be zero  
       *  @attention duplicated abscissas will be removed 
       */
      template <class XITERATOR, class YITERATOR>
      Neville ( XITERATOR  xbegin         ,
                XITERATOR  xend           ,
                YITERATOR  ybegin         ,
                YITERATOR  yend           , 
                const bool sorted = false ) 
        : Interpolation::Table ( xbegin , xend , ybegin , yend , sorted ){}
      // ======================================================================
      /** templated constructor 
       *  (very efficient: no sorting, no removal of duplicate...)
       *  @param fun       function object 
       *  @param abscissas interpolation abscissas 
       */
      template <class FUNCTION>
      Neville ( FUNCTION                        fun , 
                const Interpolation::Abscissas& a   )
        : Interpolation::Table ( fun , a ){}
      // ======================================================================
      /** templated constructor 
       *  @param begin start  of x-sequence 
       *  @param end   end    of x-sequence 
       *  @param fun   function object 
       *  @attention duplicated abscissas will be removed 
       */
      template <class ITERATOR, class FUNCTION>
      Neville ( FUNCTION fun   , 
                ITERATOR begin , 
                ITERATOR end   ) 
        : Interpolation::Table ( fun , begin, end ) {}
      // ======================================================================
      /** templated constructor 
       *  @param fun  function object 
       *  @param n    number of interpolation points 
       *  @param low  low edge of interpolation region 
       *  @param high high edge of interpolation region 
       *  @param t    interpolation type 
       */
      template <class FUNCTION>
      Neville ( FUNCTION                             fun  , 
                const unsigned short                 n    , 
                const double                         low  ,   
                const double                         high , 
                const Interpolation::Abscissas::Type t    ) 
        : Interpolation::Table ( fun , n , low , high , t ) {}
      // ======================================================================
      /// default constructor
      Neville () = default ;
      // =====================================================================
    public:
      // ======================================================================
      /// the main method: get the value of interpolated polynomial 
      double evaluate    ( const  double x ) const { return neville  ( x ) ; }
      /// the main method: get the value of interpolated polynomial 
      double operator () ( const  double x ) const { return evaluate ( x ) ; }
      // ======================================================================
      /// get the derivative   (dy/dx) at point x 
      double derivative  ( const  double x ) const { return neville2( x ).second  ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Lagrange  
     *  https://en.wikipedia.org/wiki/Lagrange_polynomial
     *  @attention  it is not CPU efficient and numerically not stable 
     */
    class Lagrange : public Interpolation::Table 
    {
    public:
      // ====================================================================== 
      /** simple contructor from interpolation points  
       *  @param x input vector of abscissas  
       */
      Lagrange ( const Interpolation::Table& p ) ;
      // ====================================================================
      // other constructors from the base class 
      // ====================================================================
      /** the simplest constructor 
       *  @param data   input data 
       *  @param sorted indicate if data already  sorted and duplicated removed  
       */
      Lagrange ( const Interpolation::TABLE& data           , 
                 const bool                  sorted = false ) 
        : Interpolation::Table ( data , sorted ) {}
      // ====================================================================== 
      /** simple contructor from abscissas and y-list 
       *  @param x input vector of abscissas  
       *  @param y input vector of y 
       *  - if vector of y is longer  than vector x, extra values are ignored 
       *  - if vector of y is shorter than vector x, missing entries are assumed to be zero  
       */
      Lagrange ( const Interpolation::Abscissas&       x , 
                 const Interpolation::Abscissas::Data& y ) 
        : Interpolation::Table ( x , y ){}
      // ===================================================================
      /** simple contructor from x&y-lists 
       *  @param x input vector of x 
       *  @param y input vector of y 
       *  - if vector of y is longer  than vector x, extra values are ignored 
       *  - if vector of y is shorter than vector x, missing entries are assumed to be zero  
       *  @attention duplicated abscissas will be removed 
       */
      Lagrange ( const Interpolation::Abscissas::Data& x , 
                 const Interpolation::Abscissas::Data& y ) 
        : Interpolation::Table ( x , y ) {}
      // ====================================================================
      /** simple constructor from map/dictionary { x : f(x) }
       *  @param data map/dictionary with data  { x : f(x) }
       *  It is elatively efficient: no sorting
       *  - "numerical" duplicates are removed
       */
      template <class KEY, class VALUE>
      Lagrange ( const std::map<KEY,VALUE>& data ) 
        : Interpolation::Table ( data ) {}
      // ====================================================================
      /** simple contructor from x&y-lists 
       *  @param xbegin start  of sequence of abscissas 
       *  @param xend   end    of sequence of abscissas 
       *  @param ybegin start  of sequence of values 
       *  @param yend   end    of sequence of values 
       *  @param sorted indicate if data already sorted and duplicates removed  
       *  - if vector of y is longer  than vector x, extra values are ignored 
       *  - if vector of y is shorter than vector x, missing entries are assumed to be zero  
       *  @attention duplicated abscissas will be removed 
       */
      template <class XITERATOR, class YITERATOR>
      Lagrange ( XITERATOR  xbegin         ,
                 XITERATOR  xend           ,
                 YITERATOR  ybegin         ,
                 YITERATOR  yend           , 
                 const bool sorted = false ) 
        : Interpolation::Table ( xbegin , xend , ybegin , yend , sorted ){}
      // ======================================================================
      /** templated constructor 
       *  (very efficient: no sorting, no removal of duplicate...)
       *  @param fun       function object 
       *  @param abscissas interpolation abscissas 
       */
      template <class FUNCTION>
      Lagrange ( FUNCTION                        fun , 
                 const Interpolation::Abscissas& a   )
        : Interpolation::Table ( fun , a ){}
      // ======================================================================
      /** templated constructor 
       *  @param begin start  of x-sequence 
       *  @param end   end    of x-sequence 
       *  @param fun   function object 
       *  @attention duplicated abscissas will be removed 
       */
      template <class ITERATOR, class FUNCTION>
      Lagrange ( FUNCTION fun   , 
                 ITERATOR begin , 
                 ITERATOR end   ) 
        : Interpolation::Table ( fun , begin, end ) {}
      // ======================================================================
      /** templated constructor 
       *  @param fun  function object 
       *  @param n    number of interpolation points 
       *  @param low  low edge of interpolation region 
       *  @param high high edge of interpolation region 
       *  @param t    interpolation type 
       */
      template <class FUNCTION>
      Lagrange ( FUNCTION                             fun  , 
                 const unsigned short                 n    , 
                 const double                         low  ,   
                 const double                         high , 
                 const Interpolation::Abscissas::Type t    ) 
        : Interpolation::Table ( fun , n , low , high , t ) {}
      // ======================================================================
      /// default constructor
      Lagrange () = default ;
      // =====================================================================
    public:
      // ======================================================================
      /// the main method: get the value of interpolated polynomial 
      double evaluate    ( const  double x ) const { return lagrange ( x ) ; }
      /// the main method: get the value of interpolated polynomial 
      double operator () ( const  double x ) const { return evaluate ( x ) ; }
      // ======================================================================
      /// get the drivative with respect to i-th parameter (dy/d(y_i)) at point x 
      double derivative  ( const double       x , 
                           const unsigned int iy ) const 
      { return lagrange2 ( x , iy ).second  ; }
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Barycentric 
     *  Very efficient Barycentric Largange interpolation
     *  For intialization the barycentric weigths are calculated as:
     *  - O(n) for <code>Chebyshev</code> and <code>Chebvyshev2</code>  abscissas 
     *  - O(n) for <code>Uniform</code> abscissas, but much slower 
     *  - O(n^2) for general case 
     *  @see Ostap::Math::Interpolation::Abscissas
     *  @see Ostap::Math::Interpolation::Abscissas::Type
     *  For evaluation it takes O(n) - that is very fast!
     *  @see Jean-Paul Berrut and Lloyd N. Trefethen, 
     *       Barycentric Lagrange Interpolation, SIAM Rev., 46(3), 501–517.
     *       ISSN (print): 0036-1445
     *       ISSN (online): 1095-7200
     *  @see https://doi.org/10.1137/S0036144502417715
     *  @see https://en.wikipedia.org/wiki/Lagrange_polynomial
     *  @see https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
     */
    class Barycentric 
    {
      // ======================================================================
    public:
      // ======================================================================      
      typedef  Interpolation::Abscissas::Data                            Data ;
      // ======================================================================
    public:
      // ======================================================================
      /** templated constructor  (efficient, no sorting/no duplicate removal)
       *  @param fun       function object 
       *  @param abscissas interpolation abscissas 
       */
      template <class FUNCTION>
      Barycentric 
      ( FUNCTION                      fun , 
        const Interpolation::Weights& w   ) 
        : m_w ( w ) 
        , m_y ( w.size() , 0.0 ) 
      { std::transform ( this->x() . begin () , this->x() . end () , m_y.begin () , fun ) ; }
      // ======================================================================
      /** templated constructor (efficient, no sorting/no duplicate removal)
       *  @param fun       function object 
       *  @param abscissas interpolation abscissas 
       */
      template <class FUNCTION>
      Barycentric 
      ( FUNCTION                        fun , 
        const Interpolation::Abscissas& a   ) 
        : m_w ( a ) 
        , m_y ( a.size() , 0.0 ) 
      { std::transform ( a.begin () , a.end () , m_y.begin () , fun ) ; }
      // ======================================================================
      /** templated constructor (very efficient, no sorting/no duplicate removal)
       *  @param fun       function object 
       *  @param abscissas interpolation abscissas 
       */
      template <class FUNCTION>
      Barycentric 
      ( FUNCTION                             fun  , 
        const unsigned int                   n    , 
        const double                         low  ,   
        const double                         high , 
        const Interpolation::Abscissas::Type t    ) 
        : Barycentric ( fun , Interpolation::Abscissas  ( n , low , high , t ) ) 
      {}
      // ======================================================================
      /** templated constructor 
       *  @param fun function object 
       *  @param begin start of abscissas sequence 
       *  @param end   end   of abscissas sequence 
       *  @attention duplicated abscissas will be removed 
       */
      template <class FUNCTION, class ITERATOR>
      Barycentric 
      ( FUNCTION fun   , 
        ITERATOR begin ,
        ITERATOR end   )
        : m_w ( begin , end ) 
        , m_y () 
      { 
        m_y.resize     ( this -> size () ) ;
        std::transform ( this -> x () . begin () , 
                         this -> x () . end   () , m_y . begin () , fun ) ; 
      }
      // ======================================================================
      /** templated constructor 
       *  @param fun function object 
       *  @param x   interpolation abscissas 
       *  @attention duplicated abscissas will be removed 
       */
      template <class FUNCTION>
      Barycentric
      ( FUNCTION    fun , 
        const Data& x   ) 
        : Barycentric ( fun , x.begin() , x.end() ) 
      {}
      // ======================================================================
      /** simple constructor from the interpolation table 
       *  @param data interpolation table            
       */
      Barycentric 
      ( const Interpolation::Table& data ) ; 
      // ======================================================================
      /** simple constructor from the interpolation data 
       *  @param data interpolation data 
       *  @param sorted indicate if data already sorted and duplicates are removed 
       */
      Barycentric
      ( const Interpolation::TABLE& data           , 
        const bool                  sorted = false ) ;
      // ======================================================================
      /** simple constructor from map/dictionary { x : f(x) }
       *  @param data map/dictionary with data  { x : f(x) }
       *  It is elatively efficient: no sorting
       *  - "numerical" duplicates are removed
       */
      template <class KEY, class VALUE>
      Barycentric 
      ( const std::map<KEY,VALUE>& data ) 
        : Barycentric ( Interpolation::Table ( data ) ) {}
      // ====================================================================
      /** simple constructor from x&y-lists 
       *  @param x input vector of x 
       *  @param y input vector of y 
       *  - if vector of y is longer  than vector x, extra values are ignored 
       *  - if vector of y is shorter than vector x, missing entries are assumed to be zero  
       */
      Barycentric 
      ( const Interpolation::Abscissas& x , 
        const Data&                     y ) ; 
      // ===================================================================
      /** simple constructor from x&y-lists 
       *  @attention this is the most CPU efficient constructor!
       *  @param x input vector of 
       *  @param y input vector of y 
       *  - if vector of y is longer  than vector x, extra values are ignored 
       *  - if vector of y is shorter than vector x, missing entries are assumed to be zero  
       */
      Barycentric
      ( const Interpolation::Weights&   x , 
        const Data&                     y ) ;
      // ===================================================================
      /** simple constructor from x&y-lists 
       *  @param x input vector of x 
       *  @param x input vector of y 
       *  - if vector of y is longer  than vector x, extra values are ignored 
       *  - if vector of y is shorter than vector x, missing entries are assumed to be zero  
       *  @attention duplicated abscissas will be removed 
       */
      Barycentric 
      ( const Data& x , const Data& y ) ;
      // ======================================================================
      /// default constructor 
      Barycentric () = default ;
      // ======================================================================
    public:
      // ======================================================================
      unsigned int n       () const { return m_w.size  () ; }
      unsigned int size    () const { return m_w.size  () ; }
      bool         empty   () const { return m_w.empty () ; }
      // ======================================================================
    public:
      // ======================================================================
      const Interpolation::Abscissas&       abscissas () const { return m_w.abscissas () ; }
      const Interpolation::Weights::Data&   weights   () const { return m_w.w ()         ; }
      const Interpolation::Abscissas::Data& x         () const { return m_w.x ()         ; }
      const Interpolation::Weights::Data&   w         () const { return m_w.w ()         ; }
      const Data&                           y         () const { return m_y              ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get abscissa:
      double abscissa ( unsigned short index ) const { return m_w.x ( index ) ; }
      double x        ( unsigned short index ) const { return m_w.x ( index ) ; }
      /// get weights 
      double weight   ( unsigned short index ) const { return m_w.w ( index ) ; }
      double w        ( unsigned short index ) const { return m_w.w ( index ) ; }
      /// get function valeus 
      double y        ( unsigned short index ) const 
      { return index < m_y.size() ? m_y[index] : m_y.back() ; }
      // ======================================================================
    public:
      // ====================================================================
      /// minimal  abscissa 
      double xmin () const { return  m_w.xmin () ; }
      /// maximal abscissas 
      double xmax () const { return  m_w.xmax () ; }
      // ====================================================================
    public:
      // ======================================================================
      /// evaluate the interpolation polynomial 
      double evaluate   ( const double  x ) const ; //evaluate the polynomial 
      /// evaluate the interpolation polynomial 
      double operator() ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /** add the point (x,y)  into  collection 
       *  @param x abscissas of the point to be added 
       *  @param y the value of function at x 
       *  @return the index of new point, or -1  if point if not added 
       */
      int add     ( const  double x ,   const  double y ) ;
      // ======================================================================
      /** remove the point with the  given index 
       *  @param index the point to be removed 
       *  @return   true if point is removed 
       */
      bool remove ( const unsigned short index ) ;
      // ======================================================================
    public: 
      // ======================================================================
      /// get the slice
      Barycentric slice 
      ( const unsigned short i , 
        const unsigned short j ) const ;
      // ====================================================================
    private: 
      // ======================================================================
      /// abscissas and weigths 
      Interpolation::Weights m_w {} ; // abscissas and weigths 
      /// function  values 
      Data                   m_y {} ; // function values 
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Newton 
     *  Interpolation polynomial
     *  https://en.wikipedia.org/wiki/Newton_polynomial
     *  @attention it is efficient  and relatively stable numerically
     */
    class Newton
    {
    public:
      // ====================================================================== 
      /** simple contructor from interpolation points  
       *  @param x input vector of abscissas  
       */
      Newton ( const Interpolation::Table& p ) ;
      // ======================================================================
      /** the simplest constructor 
       *  @param data   input data 
       *  @param sorted indicate if data already  sorted and duplicated removed  
       */
      Newton
      ( const Interpolation::TABLE& data           , 
        const bool                  sorted = false ) 
        : m_table ( data        ) 
        , m_diffs ( data.size() ) 
      { this->get_differences() ; }
      // ======================================================================
      /** simple contructor from abscissas and y-list 
       *  @param x input vector of abscissas  
       *  @param y input vector of y 
       *  - if vector of y is longer  than vector x, extra values are ignored 
       *  - if vector of y is shorter than vector x, missing entries are assumed to be zero  
       */
      Newton
      ( const Interpolation::Abscissas&       x , 
        const Interpolation::Abscissas::Data& y ) 
        : m_table ( x , y    ) 
        , m_diffs ( x.size() ) 
      { this->get_differences() ; }      
      // ======================================================================
      /** simple contructor from x&y-lists 
       *  @param x input vector of x 
       *  @param y input vector of y 
       *  - if vector of y is longer  than vector x, extra values are ignored 
       *  - if vector of y is shorter than vector x, missing entries are assumed to be zero  
       *  @attention duplicated abscissas will be removed 
       */
      Newton
      ( const Interpolation::Abscissas::Data& x , 
        const Interpolation::Abscissas::Data& y ) 
        : m_table ( x , y    ) 
        , m_diffs ( x.size() ) 
      { this->get_differences() ; }
      // ======================================================================
      /** simple constructor from map/dictionary { x : f(x) }
       *  @param data map/dictionary with data  { x : f(x) }
       *  It is elatively efficient: no sorting
       *  - "numerical" duplicates are removed
       */
      template <class KEY, class VALUE>
      Newton ( const std::map<KEY,VALUE>& data ) 
        : m_table ( data        ) 
        , m_diffs ( data.size() ) 
      { this->get_differences() ; }
      // ======================================================================
      /** simple contructor from x&y-lists 
       *  @param xbegin start  of sequence of abscissas 
       *  @param xend   end    of sequence of abscissas 
       *  @param ybegin start  of sequence of values 
       *  @param yend   end    of sequence of values 
       *  @param sorted indicate if data already sorted and duplicates removed  
       *  - if vector of y is longer  than vector x, extra values are ignored 
       *  - if vector of y is shorter than vector x, missing entries are assumed to be zero  
       *  @attention duplicated abscissas will be removed 
       */
      template <class XITERATOR, class YITERATOR>
      Newton 
      ( XITERATOR  xbegin         ,
        XITERATOR  xend           ,
        YITERATOR  ybegin         ,
        YITERATOR  yend           , 
        const bool sorted = false ) 
        : m_table ( xbegin , xend , ybegin , yend , sorted )
        , m_diffs () 
      { this->get_differences() ; }
      // ======================================================================
      /** templated constructor 
       *  (very efficient: no sorting, no removal of duplicate...)
       *  @param fun       function object 
       *  @param abscissas interpolation abscissas 
       */
      template <class FUNCTION>
      Newton
      ( FUNCTION                        fun , 
        const Interpolation::Abscissas& a   )
        : m_table ( a  , fun )
        , m_diffs ( a.size() )
      { this->get_differences() ; }
      // ======================================================================
      /** templated constructor 
       *  @param begin start  of x-sequence 
       *  @param end   end    of x-sequence 
       *  @param fun   function object 
       *  @attention duplicated abscissas will be removed 
       */
      template <class ITERATOR, class FUNCTION>
      Newton
      ( FUNCTION fun   , 
        ITERATOR begin , 
        ITERATOR end   ) 
        : m_table ( fun , begin , end )
        , m_diffs ()
      { this->get_differences() ; }
      // ======================================================================
      /** templated constructor 
       *  @param fun  function object 
       *  @param n    number of interpolation points 
       *  @param low  low edge of interpolation region 
       *  @param high high edge of interpolation region 
       *  @param t    interpolation type 
       */
      template <class FUNCTION>
      Newton
      ( FUNCTION                             fun  , 
        const unsigned short                 n    , 
        const double                         low  ,   
        const double                         high , 
        const Interpolation::Abscissas::Type t    ) 
        : m_table ( fun , n , low , high , t )
        , m_diffs ()
      { this->get_differences() ; }
      /// default constructor
      Newton () = default ;
      // ======================================================================
    public:
      // ======================================================================
      /// the main method: get the value of interpolated polynomial 
      double evaluate   ( const double x ) const ; //evaluate the polynomial 
      /// evaluate the interpolation polynomial 
      double operator() ( const double x ) const { return evaluate ( x ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get abscissa:
      double abscissa ( unsigned short index ) const { return m_table.x ( index ) ; }
      /// get the abscissa
      double x        ( unsigned short index ) const { return m_table.x ( index ) ; }
      /// get the function value 
      double y        ( unsigned short index ) const { return m_table.x ( index ) ; }
      /// get the divided differece 
      double d        ( unsigned short index ) const 
      { return index < m_diffs.size() ? m_diffs [ index ] : 0.0 ; }
      // ======================================================================
    public:
      // ======================================================================
      /// minimal  abscissa 
      double xmin () const { return  m_table.xmin () ; }
      /// maximal abscissas 
      double xmax () const { return  m_table.xmax () ; }
      // ======================================================================
    public: // show some internal data 
      // ======================================================================
      const Interpolation::Table& table () const { return m_table         ; }
      const Interpolation::TABLE& data  () const { return m_table.data () ; }
      // ======================================================================
    public: // add & remove interpolation points 
      // ======================================================================
      /** add the point (x,y) into interpolation table 
       *  @param x abscissas of the point to be added 
       *  @param y the value of function at x 
       *  @return the index of new point, or -1  if point if not added
       */
      int add     ( const  double x ,   const  double y ) ;
      // ======================================================================
      /** remove the point with the  given index 
       *  @param index the point to be removed 
       *  @return   true if point is removed 
       */
      bool remove ( unsigned short index ) ;
      // ======================================================================
    public: // make a slice fo the given  ascissas 
      // ======================================================================
      /// make a slice for the given range of points 
      Newton slice ( const int i , const int j ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// swap two collections of points 
      friend void swap ( Newton& a , Newton& b ) 
      { 
        swap ( a.m_table , b.m_table ) ;
        swap ( a.m_diffs , b.m_diffs ) ;
      }
      // ======================================================================
    private: // get the divided differences 
      // ======================================================================
      /// get the divided differences 
      void get_differences () ; // get the divided differences 
      // ======================================================================
    private:
      // ======================================================================
      /// the interpolation table 
      Interpolation::Table m_table {} ; // the interpolation table 
      /// vector of divided differences 
      std::vector<double>  m_diffs {} ; // vector of divided differences 
      // ======================================================================
    } ;
    // ========================================================================
    namespace Interpolation 
    {
      // ======================================================================
      /** Very efficient Barycentric Lagrange Interpolation: O(n)
       *  @code 
       *  auto fun             = [] ( double x ) { return std::sin(x) ; } ;
       *  Weights  weights     = ... ;
       *  auto interpolant     = lagrange ( fun , weights ) ;
       *  const double   x     = ... ;
       *  const double   value = interpolant ( value ) 
       *  @endcode 
       *  @param func the function 
       *  @param weights precomputed barycentric weights 
       *  @see Jean-Paul Berrut and Lloyd N. Trefethen, 
       *       Barycentric Lagrange Interpolation, SIAM Rev., 46(3), 501–517.
       *       ISSN (print): 0036-1445
       *       ISSN (online): 1095-7200
       *  @see https://doi.org/10.1137/S0036144502417715
       *  @see https://en.wikipedia.org/wiki/Lagrange_polynomial
       *  @see https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
       *  @see Ostap::Math::Barycentric
       */
      template <class FUNCTION>
      inline 
      Ostap::Math::Barycentric
      lagrange ( FUNCTION         func    ,  
                 const Weights&   weights )
      { return Ostap::Math::Barycentric ( func , weights ) ; }
      // ======================================================================
      /** Very efficient Barycentric Lagrange Interpolation: O(n)
       *  @code 
       *  auto fun             = [] ( double x ) { return std::sin(x) ; } ;
       *  Abscissas abscissas  = ... ;
       *  auto interpolant     = lagrange ( fun , abscissas ) ;
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param abscissas interpolation abscissas  
       *  @see Jean-Paul Berrut and Lloyd N. Trefethen, 
       *       Barycentric Lagrange Interpolation, SIAM Rev., 46(3), 501–517.
       *       ISSN (print): 0036-1445
       *       ISSN (online): 1095-7200
       *  @see https://doi.org/10.1137/S0036144502417715
       *  @see https://en.wikipedia.org/wiki/Lagrange_polynomial
       *  @see https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
       *  @see Ostap::Math::Barycentric
       */
      template <class FUNCTION>
      inline  
      Ostap::Math::Barycentric
      lagrange ( FUNCTION         func      ,  
                 const Abscissas& abscissas )
      { return Ostap::Math::Barycentric ( func , abscissas ) ; }
      // ======================================================================
      /** Very efficient Barycentric Lagrange Interpolation: 
       *  @code 
       *  auto fun             = [] ( double x ) { return std::sin(x) ; } ;
       *  auto interpolant     = lagrange ( fun , 12 , 0. , 1. , Abscissas::Chebyshev ) 
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param N         number of interpolation absicssas
       *  @param low       low  edge of interpolating ranage 
       *  @param high      high edge of interpolating ranage 
       *  @param i         interpolation type       
       *  @see Jean-Paul Berrut and Lloyd N. Trefethen, 
       *       Barycentric Lagrange Interpolation, SIAM Rev., 46(3), 501–517.
       *       ISSN (print): 0036-1445
       *       ISSN (online): 1095-7200
       *  @see https://doi.org/10.1137/S0036144502417715
       *  @see https://en.wikipedia.org/wiki/Lagrange_polynomial
       *  @see https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
       *  @see Ostap::Math::Barycentric
       */
      template <class FUNCTION>
      inline 
      Ostap::Math::Barycentric
      lagrange ( FUNCTION              func ,  
                 const unsigned int    N    , 
                 const double          low  ,   
                 const double          high , 
                 const Abscissas::Type t    )        
      { return Ostap::Math::Barycentric ( func , N , low , high , t ) ; }      
      // ======================================================================
      /** Very efficient Barycentric Lagrange Interpolation: 
       *  @code 
       *  auto fun             = [] ( double x ) { return std::sin(x) ; } ;
       *  auto interpolant     = lagrange ( fun , {{ 0.0 , 0.1 , 0.2,  0.3, 0.7, 1.0 }} ) 
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param x         interpolation abscissas 
       *  @see Jean-Paul Berrut and Lloyd N. Trefethen, 
       *       Barycentric Lagrange Interpolation, SIAM Rev., 46(3), 501–517.
       *       ISSN (print): 0036-1445
       *       ISSN (online): 1095-7200
       *  @see https://doi.org/10.1137/S0036144502417715
       *  @see https://en.wikipedia.org/wiki/Lagrange_polynomial
       *  @see https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
       *  @see Ostap::Math::Barycentric
       */      
      template <class FUNCTION>
      inline 
      Ostap::Math::Barycentric
      lagrange  ( FUNCTION                               func ,  
                  const  Ostap::Math::Barycentric::Data& x    )
      { return Ostap::Math::Barycentric ( func , x ) ; }
      // ======================================================================
      /** Very efficient Barycentric Lagrange Interpolation: 
       *  @code 
       *  Table          table = ... ; // interpolation table 
       *  auto interpolant     = lagrange ( table ) 
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param table    interpolation table 
       *  @see Jean-Paul Berrut and Lloyd N. Trefethen, 
       *       Barycentric Lagrange Interpolation, SIAM Rev., 46(3), 501–517.
       *       ISSN (print): 0036-1445
       *       ISSN (online): 1095-7200
       *  @see https://doi.org/10.1137/S0036144502417715
       *  @see https://en.wikipedia.org/wiki/Lagrange_polynomial
       *  @see https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
       *  @see Ostap::Math::Barycentric
       */      
      inline 
      Ostap::Math::Barycentric
      lagrange ( const Table& data )   
      { return Ostap::Math::Barycentric ( data ) ; }
      // ======================================================================      
      /** Very efficient Barycentric Lagrange Interpolation: 
       *  @code 
       *  TABLE          table = ... ; // interpolation table 
       *  auto interpolant     = lagrange ( table ) 
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param table    interpolation table 
       *  @param sorted indicate if data is sorted and duplicates are removed 
       *  @see Jean-Paul Berrut and Lloyd N. Trefethen, 
       *       Barycentric Lagrange Interpolation, SIAM Rev., 46(3), 501–517.
       *       ISSN (print): 0036-1445
       *       ISSN (online): 1095-7200
       *  @see https://doi.org/10.1137/S0036144502417715
       *  @see https://en.wikipedia.org/wiki/Lagrange_polynomial
       *  @see https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
       *  @see Ostap::Math::Barycentric
       */
      inline 
      Ostap::Math::Barycentric
      lagrange ( const TABLE& data , const bool sorted = false )   
      { return Ostap::Math::Barycentric ( data , sorted ) ; }
      // ======================================================================
      /** Very efficient Barycentric Lagrange Interpolation: 
       *  @code
       *  typdef std::map<double,double> MAP ;
       *  MAP data ;
       *  data[1.0] = std::sin ( 1.0 ) ;
       *  data[1.5] = std::sin ( 1.5 ) ;
       *  data[2.0] = std::sin ( 2.0 ) ;
       *  data[2.5] = std::sin ( 2.5 ) ;
       *  auto interpolant     = lagrange ( data) 
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param data      interpolation table 
       *  @see Jean-Paul Berrut and Lloyd N. Trefethen, 
       *       Barycentric Lagrange Interpolation, SIAM Rev., 46(3), 501–517.
       *       ISSN (print): 0036-1445
       *       ISSN (online): 1095-7200
       *  @see https://doi.org/10.1137/S0036144502417715
       *  @see https://en.wikipedia.org/wiki/Lagrange_polynomial
       *  @see https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
       *  @see Ostap::Math::Barycentric
       */
      template <class KEY,class VALUE>
      inline 
      Ostap::Math::Barycentric
      lagrange ( const std::map<KEY,VALUE>& data )   
      { return Ostap::Math::Barycentric ( data ) ; }
      // ======================================================================
      /** Very efficient Barycentric Lagrange Interpolation: 
       *  @code
       *  typedef std::vector<double> DATA ;
       *  DATA data {{ std::sin( 0.1 ) , std::sin( 0.1 ) , std::sin( 0.3) }} ;
       *  Abscissas a          = ... ;
       *  auto interpolant     = lagrange ( a , data ) 
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param x     interpolation abscissas 
       *  @param y     interpolation data 
       *  @see Jean-Paul Berrut and Lloyd N. Trefethen, 
       *       Barycentric Lagrange Interpolation, SIAM Rev., 46(3), 501–517.
       *       ISSN (print): 0036-1445
       *       ISSN (online): 1095-7200
       *  @see https://doi.org/10.1137/S0036144502417715
       *  @see https://en.wikipedia.org/wiki/Lagrange_polynomial
       *  @see https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
       *  @see Ostap::Math::Barycentric
       */
      inline   
      Ostap::Math::Barycentric
      lagrange ( const Abscissas&                      x , 
                 const Ostap::Math::Barycentric::Data& y )
      { return Ostap::Math::Barycentric ( x , y ) ; }
      // ======================================================================
      /** Very efficient Barycentric Lagrange Interpolation: 
       *  @code
       *  typedef std::vector<double> DATA ;
       *  DATA data {{ std::sin( 0.1 ) , std::sin( 0.1 ) , std::sin( 0.3) }} ;
       *  Weights w          = ... ;
       *  auto interpolant     = lagrange ( w , data ) 
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param w     precomputed   interpolation weights 
       *  @param y     interpolation data 
       *  @see Jean-Paul Berrut and Lloyd N. Trefethen, 
       *       Barycentric Lagrange Interpolation, SIAM Rev., 46(3), 501–517.
       *       ISSN (print): 0036-1445
       *       ISSN (online): 1095-7200
       *  @see https://doi.org/10.1137/S0036144502417715
       *  @see https://en.wikipedia.org/wiki/Lagrange_polynomial
       *  @see https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
       *  @see Ostap::Math::Barycentric
       */
      inline 
      Ostap::Math::Barycentric
      lagrange ( const Weights&                        w , 
                 const Ostap::Math::Barycentric::Data& y )
      { return Ostap::Math::Barycentric ( w , y ) ; }
      // ======================================================================      
      /** Very efficient Barycentric Lagrange Interpolation: 
       *  @code
       *  typedef std::vector<double> DATA ;
       *  DATA xx {{           0.1   ,           0.2   ,            0.3   }} ;
       *  DATA yy {{ std::sin( 0.1 ) , std::sin( 0.2 ) , std::sin ( 0.3 ) }} ;
       *  Weights w          = ... ;
       *  auto interpolant     = lagrange ( xx , yy) 
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param x     interpolation abscissas 
       *  @param y     interpolation data 
       *  @see Jean-Paul Berrut and Lloyd N. Trefethen, 
       *       Barycentric Lagrange Interpolation, SIAM Rev., 46(3), 501–517.
       *       ISSN (print): 0036-1445
       *       ISSN (online): 1095-7200
       *  @see https://doi.org/10.1137/S0036144502417715
       *  @see https://en.wikipedia.org/wiki/Lagrange_polynomial
       *  @see https://people.maths.ox.ac.uk/trefethen/barycentric.pdf
       *  @see Ostap::Math::Barycentric
       */
      inline 
      Ostap::Math::Barycentric
      lagrange ( const Ostap::Math::Barycentric::Data& x , 
                 const Ostap::Math::Barycentric::Data& y )
      { return Ostap::Math::Barycentric ( x , y ) ; }
      // ======================================================================      
      //  Newton interpolation 
      // ======================================================================      
      /** Newton interpolation
       *  @code 
       *  auto fun             = [] ( double x ) { return std::sin(x) ; } ;
       *  Abscissas abscissas  = ... ;
       *  auto interpolant     = newton ( fun , abscissas ) ;
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param abscissas interpolation abscissas  
       *  @see https://en.wikipedia.org/wiki/Newton_polynomial
       *  @see Ostap::Math::Newton 
       */
      template <class FUNCTION>
      inline 
      Ostap::Math::Newton
      newton ( FUNCTION         func      ,  
               const Abscissas& abscissas )
      { return Ostap::Math::Newton ( func , abscissas ) ; }
      // =======================================================================
      /** Newton interpolation
       *  @code 
       *  auto fun             = [] ( double x ) { return std::sin(x) ; } ;
       *  auto interpolant     = newton ( fun , 12 , 0. , 1. , Abscissas::Chebyshev ) 
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param N         number of interpolation absicssas
       *  @param low       low  edge of interpolating ranage 
       *  @param high      high edge of interpolating ranage 
       *  @param i         interpolation type       
       *  @see https://en.wikipedia.org/wiki/Newton_polynomial
       *  @see Ostap::Math::Newton 
       */
      template <class FUNCTION>
      inline 
      Ostap::Math::Newton
      newton ( FUNCTION              func ,  
               const unsigned int    N    , 
               const double          low  ,   
               const double          high , 
               const Abscissas::Type t    )        
      { return Ostap::Math::Newton( func , N , low , high , t ) ; }
      // ======================================================================      
      /** Newton interpolation
       *  @code 
       *  auto fun             = [] ( double x ) { return std::sin(x) ; } ;
       *  auto interpolant     = newton ( fun , {{ 0.0 , 0.1 , 0.2,  0.3, 0.7, 1.0 }} ) 
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param x         interpolation abscissas 
       *  @see https://en.wikipedia.org/wiki/Newton_polynomial
       *  @see Ostap::Math::Newton 
       */
      template <class FUNCTION>
      inline 
      Ostap::Math::Newton
      newton ( FUNCTION               func ,  
               const Abscissas::Data& x    )
      { return Ostap::Math::Newton ( func , x ) ; }
      // ======================================================================      
      /** Newton interpolation
       *  @code 
       *  Table          table = ... ; // interpolation table 
       *  auto interpolant     = newton ( table ) 
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param table    interpolation table  
       *  @see https://en.wikipedia.org/wiki/Newton_polynomial
       *  @see Ostap::Math::Newton 
       */
      inline 
      Ostap::Math::Newton
      Newton  ( const Table& data )   
      { return Ostap::Math::Newton ( data ) ; }
      // ======================================================================      
      /** Newton interpolation
       *  @code 
       *  TABLE          table = ... ; // interpolation table 
       *  auto interpolant     = newton ( table ) 
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param table    interpolation table 
       *  @param sorted indicate if data is sorted and duplicates are removed 
       *  @see https://en.wikipedia.org/wiki/Newton_polynomial
       *  @see Ostap::Math::Newton 
       */
      inline 
      Ostap::Math::Newton
      newton ( const TABLE& data , const bool sorted = false )   
      { return Ostap::Math::Newton ( data , sorted ) ; }
      // ======================================================================
      /** Newton interpolation
       *  @code
       *  typdef std::map<double,double> MAP ;
       *  MAP data ;
       *  data[1.0] = std::sin ( 1.0 ) ;
       *  data[1.5] = std::sin ( 1.5 ) ;
       *  data[2.0] = std::sin ( 2.0 ) ;
       *  data[2.5] = std::sin ( 2.5 ) ;
       *  auto interpolant     = newton ( data) 
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param data      interpolation table 
       *  @see https://en.wikipedia.org/wiki/Newton_polynomial
       *  @see Ostap::Math::Newton 
       */
       template <class KEY,class VALUE>
      inline 
      Ostap::Math::Newton
      newton ( const std::map<KEY,VALUE>& data )   
      { return Ostap::Math::Newton ( data ) ; }
      // ======================================================================
      /** Newton interpolation
       *  @code
       *  typedef std::vector<double> DATA ;
       *  DATA data {{ std::sin( 0.1 ) , std::sin( 0.1 ) , std::sin( 0.3) }} ;
       *  Abscissas a          = ... ;
       *  auto interpolant     = lagrange ( a , data ) 
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param x     interpolation abscissas 
       *  @param y     interpolation data 
       *  @see https://en.wikipedia.org/wiki/Newton_polynomial
       *  @see Ostap::Math::Newton 
       */
      inline   
      Ostap::Math::Newton
      newton ( const Abscissas&       x , 
               const Abscissas::Data& y )
      { return Ostap::Math::Newton( x , y ) ; }
      // ======================================================================
      /** Newton interpolation
       *  @code
       *  typedef std::vector<double> DATA ;
       *  DATA xx {{           0.1   ,           0.2   ,            0.3   }} ;
       *  DATA yy {{ std::sin( 0.1 ) , std::sin( 0.2 ) , std::sin ( 0.3 ) }} ;
       *  Weights w          = ... ;
       *  auto interpolant     = lagrange ( xx , yy) 
       *  const double  x      = ... ;
       *  const double  value  = interpolant ( value ) 
       *  @endcode 
       *  @param func      the function 
       *  @param x     interpolation abscissas 
       *  @param y     interpolation data 
       *  @see https://en.wikipedia.org/wiki/Newton_polynomial
       *  @see Ostap::Math::Newton 
       */
      inline 
      Ostap::Math::Newton
      newton ( const Abscissas::Data& x , 
               const Abscissas::Data& y )
      { return Ostap::Math::Newton ( x , y ) ; }
      // ======================================================================
    }
    // ========================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
/*  Very simple lagrange interpolation 
 *
 *  @param  xbegin INPUT start iterator for the sequence of abscissas 
 *  @param  xend   INPUT end   iterator for the sequence of abscissas 
 *  @param  ybegin INPUT start iterator for the sequence of values 
 *  @param  yend   INPUT end   iterator for the sequence of values 
 *  @param  x      INPUT evaluate the polynomial in this point
 *  @param  xvalue INPUT adapter for x-values  
 *  @param  yvalue INPUT adapter for y-values  
 *  @return the value of Largange interpolation polynomial at point x
 *
 *  - If y-vector is shorter than x-vector, it is assumed to be zero padded
 *  - If y-vector is longer  than x-vector, extra avalues are ignored       
 *  - If x-vector is empty, polynomial is zero
 *
 *  @see https://en.wikipedia.org/wiki/Lagrange_polynomial
 *  @warning it could be CPU inefficient
 *  @warning it should *NOT* be applied for very long sequence of points 
 *           (e.g. >20) due to bad numerical  stability and Runge phenomenon
 *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
 *  @date 2016-07-23
 */
// ============================================================================
template <class XITERATOR, 
          class YITERATOR, 
          class RESULT   ,
          class XADAPTER ,
          class YADAPTER >
inline RESULT 
Ostap::Math::Interpolation::lagrange
( XITERATOR    xbegin , 
  XITERATOR    xend   , 
  YITERATOR    ybegin , 
  YITERATOR    yend   , 
  const double x      , 
  RESULT       result ,
  XADAPTER     xvalue ,   // adaptor to get y-value
  YADAPTER     yvalue )   // adaptor to get x-value
{
  /// 1)special case no data 
  if ( xbegin == xend || ybegin == yend ) { return result ; }    // RETURN
  const typename std::iterator_traits<XITERATOR>::difference_type nx = std::distance ( xbegin , xend ) ;
  /// 2) special case: constant function 
  if ( 1 == nx ) { return  result + yvalue ( *ybegin )    ; } // RETURN
  const typename std::iterator_traits<XITERATOR>::difference_type ny = std::distance ( ybegin , yend ) ;
  //
  /// redefine yend : skip extra y-values 
  if ( nx < ny ) { std::advance ( yend , nx - ny ) ; }
  //
  unsigned int i  =      0 ;
  XITERATOR    ix = xbegin ;
  YITERATOR    iy = ybegin ;
  //
  for ( ; iy != yend ; ++iy , ++ix , ++i ) 
  {
    const double xi = xvalue ( *ix ) ; // get values 
    //
    double       r = 1 ;
    unsigned int j = 0 ;
    for ( XITERATOR jx = xbegin ; jx != xend ; ++jx, ++j  )
    {
      const double xj = xvalue ( *jx ) ;  // get the value 
      if ( i == j ) { continue ; }
      r *= ( x - xj ) / ( xi - xj ) ;
    }
    //
    result += r * yvalue ( *iy ) ;
  }
  return result ;                            // RETURN 
}
// ============================================================================
/*  simple interpolation using Neville's algorithm 
 *
 *  In general it should be faster than largange algorithm, 
 *  but it includes a copy on input data, that could affect CPU performance
 *  Numerically it is more stable that Lagrange interpolation, 
 *  but anyhow it also should not be used for very high polynomial 
 *  degrees (e.g.>20), especially for uniform grid (Runge phenomenon)
 *
 *  - If y-vector is shorter than x-vector, it is assumed to be zero padded
 *  - If y-vector is longer  than x-vector, extra avalues are ignored       
 *  - If x-vector is empty, polynomial is zero
 *
 *  @param  xbegin INPUT start iterator for the sequence of abscissas 
 *  @param  xend   INPUT end   iterator for the sequence of abscissas 
 *  @param  ybegin INPUT start iterator for the sequence of values 
 *  @param  yend   INPUT end   iterator for the sequence of values 
 *  @param  x      INPUT evaluate the polynomial in this point
 *  @param  xvalue INPUT adapter for x-values  
 *  @param  yvalue INPUT adapter for y-values  
 *  @return the value of interpolation polynomial at point x
 *  @see https://en.wikipedia.org/wiki/Neville%27s_algorithm
 */
// ============================================================================
template <class XITERATOR, 
          class YITERATOR, 
          class XADAPTER ,
          class YADAPTER >
inline double 
Ostap::Math::Interpolation::neville
( XITERATOR    xbegin ,
  XITERATOR    xend   ,
  YITERATOR    ybegin ,
  YITERATOR    yend   ,
  const double x      , 
  XADAPTER     xvalue , // adaptor to get y-value
  YADAPTER     yvalue ) // adaptor to get x-value  
{
  /// 1)special case no data 
  if ( xbegin == xend || ybegin == yend ) { return                   0 ; } // RETURN
  const typename std::iterator_traits<XITERATOR>::difference_type NX = std::distance ( xbegin , xend ) ;  
  /// 2) special case: constant function 
  if ( 1 == NX                          ) { return  yvalue ( *ybegin ) ; } // RETURN
  const typename std::iterator_traits<YITERATOR>::difference_type NY = std::distance ( ybegin , yend ) ;  
  // temporary storage 
  std::vector<double> _y ( NX , 0 ) ;
  YITERATOR ylast = yend ;
  if ( NX < NY ) { std::advance ( ylast , NX - NY ) ; }
  std::transform ( ybegin , ylast , _y.begin() , yvalue ) ;
  // simple version of neville 
  return neville ( xbegin , xend , _y.begin() , x , xvalue ) ;
}
// ============================================================================
/*  simple interpolation using Neville's algorithm:
 *  - evaluate the interpolation polynomial and the derivative 
 *
 *  In general it should be faster than largange algorithm, 
 *  but it includes a copy on input data, that could affect CPU performance
 *  Numerically it is more stable that Lagrange interpolation, 
 *  but anyhow it also should not be used for very high polynomial 
 *  degrees (e.g.>20), especially for uniform grid (Runge phenomenon)
 *
 *  - If y-vector is shorter than x-vector, it is assumed to be zero padded
 *  - If y-vector is longer  than x-vector, extra avalues are ignored       
 *  - If x-vector is empty, polynomial is zero
 *
 *  @param  xbegin INPUT start iterator for the sequence of abscissas 
 *  @param  xend   INPUT end   iterator for the sequence of abscissas 
 *  @param  ybegin INPUT start iterator for the sequence of values 
 *  @param  yend   INPUT end   iterator for the sequence of values 
 *  @param  x      INPUT evaluate the polynomial in this point
 *  @param  xvalue INPUT adapter for x-values  
 *  @param  yvalue INPUT adapter for y-values  
 *  @return the value of interpolation polynomial and derivative at point x 
 *  @see https://en.wikipedia.org/wiki/Neville%27s_algorithm
 */
// ============================================================================
template <class XITERATOR, 
          class YITERATOR, 
          class XADAPTER ,
          class YADAPTER >
inline std::pair<double,double>
Ostap::Math::Interpolation::neville2
( XITERATOR    xbegin ,
  XITERATOR    xend   ,
  YITERATOR    ybegin ,
  YITERATOR    yend   ,
  const double x      , 
  XADAPTER     xvalue , // adaptor to get y-value
  YADAPTER     yvalue ) // adaptor to get x-value  
{
  /// 1)special case no data 
  if ( xbegin == xend || ybegin == yend ) { return std::make_pair ( 0                  , 0 ) ; } // RETURN
  const typename std::iterator_traits<XITERATOR>::difference_type NX = std::distance ( xbegin , xend ) ;  
  /// 2) special case: constant function 
  if ( 1 == NX                          ) { return std::make_pair ( yvalue ( *ybegin ) , 0 ) ; } // RETURN
  const typename std::iterator_traits<YITERATOR>::difference_type NY = std::distance ( ybegin , yend ) ;  
  // temporary storage 
  std::vector<long double> _y ( NX , 0 ) ;
  std::vector<long double> _d ( NX , 0 ) ;
  YITERATOR ylast = yend ;
  if ( NX < NY ) { std::advance ( ylast , NX - NY ) ; }
  std::transform  ( ybegin , ylast , _y.begin() , yvalue ) ;
  // simple version of neville 
  return neville2 ( xbegin , xend , _y.begin() , _d.begin() , x , xvalue ) ;
}
// ============================================================================
/* simple interpolation using Neville's algorithm: 
 *
 *
 *  In general it should be faster than largange algorithm, 
 *  but it modifies input data!
 *
 *  Numerically it is more stable than Lagrange interpolation, 
 *  but anyhow it also should not be used for very high polynomial 
 *  degrees (e.g.>20), especially for uniform grid (Runge phenomenon)
 *
 *  - If y-vector is shorter than x-vector, it is assumed to be zero padded
 *  - If y-vector is longer  than x-vector, extra avalues are ignored       
 *  - If x-vector is empty, polynomial is zero
 *
 *  @attention y-sequence must be "simple", convertible to doubles  
 *  @attention y-sequence is *MODIFIED*
 *
 *  @param  xbegin INPUT  start iterator for the sequence of abscissas 
 *  @param  xend   INPUT  end   iterator for the sequence of abscissas 
 *  @param  ybegin UPDATE start iterator for the sequence of values 
 *  @param  yend   UPDATE end   iterator for the sequence of values 
 *  @param  x      INPUT  evaluate the polynomial in this point
 *  @param  xvalue INPUT  adapter for x-values  
 *  @param  yvalue INPUT  adapter for y-values  
 *  @return the value of  interpolation polynomial at point x
 *  @see https://en.wikipedia.org/wiki/Neville%27s_algorithm
 */
// ============================================================================
template <class XITERATOR, 
          class YITERATOR, 
          class XADAPTOR >
inline double 
Ostap::Math::Interpolation::neville
( XITERATOR    xbegin , 
  XITERATOR    xend   , 
  YITERATOR    ybegin , // NON-CONST 
  const double x      , 
  XADAPTOR     xvalue )
{
  if ( xbegin == xend ) { return 0        ; } // RETURN 
  const typename std::iterator_traits<XITERATOR>::difference_type NX = std::distance ( xbegin , xend ) ;
  if ( 1 == NX        ) { return *ybegin  ; } // RETURN
  //
  for ( unsigned int k = 1 ; k < NX ; ++k ) 
  {
    for ( unsigned int i = 0 ; i < NX - k ; ++i ) 
    {
      const long double xj  = xvalue ( *(xbegin +(i + k) ) ) ;
      const long double xi  = xvalue ( *(xbegin + i      ) ) ;
      const long double yi  = *(ybegin + i      ) ;
      const long double yi1 = *(ybegin +(i + 1) ) ;
      *(ybegin+i) = ( ( x - xj ) * yi + ( xi - x ) * yi1 ) /( xi - xj ) ;
    }
  }
  return *ybegin ;
}
// ============================================================================
/* simple interpolation using Neville's algorithm with simultaneous 
 *  estimation of the derivative  
 *
 *  Numerically it is more stable than Lagrange interpolation, 
 *  but anyhow it also should not be used for very high polynomial 
 *  degrees (e.g.>20), especially for uniform grid (Runge phenomenon)
 *
 *  @attention y-sequence must be "simple", convertible to doubles  
 *  @attention y-sequence is *MODIFIED*
 *  @attention d-sequence is *MODIFIED*
 *
 *  @param  xbegin INPUT  start iterator for the sequence of abscissas 
 *  @param  xend   INPUT  end   iterator for the sequence of abscissas 
 *  @param  ybegin UPDATE start iterator for the sequence of values 
 *  @param  x      INPUT  evaluate the polynomial in this point
 *  @param  xvalue INPUT  adapter for x-values  
 *  @param  yvalue INPUT  adapter for y-values  
 *  @return the pair (function,derivative) at point x 
 *  @see https://en.wikipedia.org/wiki/Neville%27s_algorithm
 */
// ============================================================================
template <class XITERATOR, 
          class YITERATOR, 
          class DITERATOR, 
          class XADAPTOR >
inline std::pair<double,double>
Ostap::Math::Interpolation::neville2
( XITERATOR    xbegin , 
  XITERATOR    xend   , 
  YITERATOR    ybegin , // NON-const!
  DITERATOR    dbegin , // NON-const!
  const double x      , 
  XADAPTOR     xvalue ) 
{
  if ( xbegin == xend ) { return std::make_pair ( 0       , 0 ) ; } // RETURN 
  const typename std::iterator_traits<XITERATOR>::difference_type NX = std::distance ( xbegin , xend ) ;
  if ( 1 == NX        ) { return std::make_pair ( *ybegin , 0 ) ; } // RETURN
  //
  for ( unsigned int k = 1 ; k < NX ; ++k ) 
  {
    for ( unsigned int i = 0 ; i < NX - k ; ++i ) 
    {
      const long double xj  = xvalue ( *(xbegin +(i + k) ) ) ;
      const long double xi  = xvalue ( *(xbegin + i      ) ) ;
      const long double yi  = *(ybegin + i      ) ;
      const long double yi1 = *(ybegin +(i + 1) ) ;
      const long double di  = *(dbegin + i      ) ;
      const long double di1 = *(dbegin +(i + 1) ) ;
      *(ybegin+i) = ( ( x - xj ) * yi      + ( xi - x ) * yi1        ) / ( xi - xj ) ;
      *(dbegin+i) = ( ( x - xj ) * di + yi + ( xi - x ) * di1 - yi1  ) / ( xi - xj ) ;
    }
  }
  return std::make_pair ( *ybegin , *dbegin ) ;
}
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_INTERPOLATION_H
// ============================================================================
