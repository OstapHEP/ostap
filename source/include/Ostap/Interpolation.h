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
     *  @see Kai Hormann, "Barycentric interpolaiton", 
     *  @see https://www.inf.usi.ch/hormann/papers/Hormann.2014.BI.pdf            
     *
     * 
     *  None of this method should be applies for "long" sequence of 
     *  interpolation points ( e.g. >20), especially for uniform grid 
     *  (see https://en.wikipedia.org/wiki/Runge%27s_phenomenon)
     *  
     *  - Lagrange interpolation is numerically not very stable, and rather slow O(n^2)
     *  - Neville's algorithm has (a little bit) better numerical stability and a bit faster
     *  - True Barycentric Lagrange interpolation is very efficienct: 
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
        typedef std::vector<double>  Data           ; // the actual vector of data 
        typedef Data::const_iterator const_iterator ; // iterator type 
        typedef Data::const_iterator iterator       ; // iterator type
        // ====================================================================
        enum AType 
          { Generic          = -1          , 
            Uniform          =  0          , 
            Chebyshev        =  1          , //   roots of T_n     ( x ) 
            Chebyshev1       =  Chebyshev  ,
            GaussChebyshev   =  Chebyshev  , 
            Chebyshev2       =  2          , // extrema of T_{n-1} ( x ) 
            Lobatto          =  Chebyshev2 , 
            ChebyshevLobatto =  Chebyshev2 ,
            GaussLobatto     =  Chebyshev2 } ;
        // ====================================================================
      public:
        // ====================================================================
        /** create the abscissas from vector of abscissas 
         *  @param x       input vector of abscissas (to be sorted internally)
         *  @param sorted  indicate if input data is already sorted 
         *  @attention duplicated abscissas will be removed 
         */
        Abscissas
        ( const Data&         x     , 
          const bool sorted = false ) ;
        // ===================================================================
        /** templated constructor from the arbitrary sequence
         *  @param begin start of the sequnce 
         *  @param end   end   of the sequnce 
         *  @param sorted  indicate if input data is already sorted 
         *  @attention duplicated abscissas will be removed 
         */
        template <class ITERATOR,
                  typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                  typename = std::enable_if<std::is_convertible<value_type,long double>::value> >
        Abscissas 
        ( ITERATOR begin , 
          ITERATOR end   , 
          const bool sorted = false ) 
          :  m_x     ( begin , end ) 
          ,  m_atype ( Generic     )  
        { 
          if ( !sorted ) { this->sort () ; } 
          if ( !empty () )
          {
            m_xmin = m_x.front () ;
            m_xmax = m_x.back  () ;            
          } 
        }
        // =================================================================
        /** templated constructor from the arbitrary sequence
         *  @param begin start of the sequnce 
         *  @param end   end   of the sequnce 
         *  @param fun   the function to be applied 
         *  @param sorted  indicate if input data is already sorted 
         *  @attention duplicated abscissas will be removed 
         */
        template <class ITERATOR, class FUNCTION>
        Abscissas
        ( ITERATOR   begin          , 
          ITERATOR   end            , 
          FUNCTION   fun            , 
          const bool sorted = false ) 
          : m_x     ( std::distance ( begin, end ) , 0.0 )  
          , m_atype ( Generic     )
        {
          std::transform ( begin , end , m_x.begin() , fun );
          if ( !sorted   ) { this->sort () ; } 
          if ( !empty () )
          {
            m_xmin = m_x.front () ;
            m_xmax = m_x.back  () ;            
          }
        }
        // ====================================================================
        /** special constructor for the given interpolation type 
         *  @param n    number of interpolation points 
         *  @param low  low edge of the interval 
         *  @param high high edge of the interval 
         *  @param t    interpolation type 
         *  @see Ostap::Math::Interpolation::Abscissas::AType 
         */
        Abscissas 
        ( const unsigned short n                ,
          const double         low              ,   
          const double         high             , 
          const AType          t    = Chebyshev ) ;
        // ====================================================================
        /// default constructor  
        Abscissas ()                    = default ;
        // ====================================================================
      public:
        // ====================================================================
        /// number of interpolating points 
        unsigned int n       () const { return m_x.size  () ; }
        /// number of interpolating points 
        unsigned int size    () const { return m_x.size  () ; }
        /// no interpolation points ?
        bool         empty   () const { return m_x.empty () ; }
        // ====================================================================
        /// type of interpolation points 
        AType        atype   () const { return m_atype      ; }
        // ====================================================================
      public:
        // ====================================================================
        /// get all abscissas 
        const Data&    x     () const { return m_x ; }
        /// get all abscissas (as cast-operator)
        operator const Data& () const { return m_x ; }
        // ====================================================================
        ///  get the absissas for the given index 
        double x    ( const unsigned short index ) const 
        { 
          static const double s_nan = std::numeric_limits<double>::quiet_NaN () ;
          return m_x.empty() ? s_nan : index < m_x.size() ? m_x [ index ] : m_x.back()   ; 
        }
        // ====================================================================
        /// get abscissa for the given index 
        double operator[] ( const unsigned short index ) const 
        { return x ( index ) ; }
        // ====================================================================
      public:
        // ====================================================================
        /// minimal  abscissa 
        double xmin () const { return m_xmin ; }
        /// maximal abscissas 
        double xmax () const { return m_xmax ; }
        // ====================================================================
      public: // expose (const) iterators 
        // ====================================================================
        /// begin-iterator 
        Data::const_iterator begin () const { return  m_x.begin () ; }
        /// end-iterator 
        Data::const_iterator end   () const { return  m_x.end   () ; }        
        // ====================================================================
      public: // add & remove the point 
        // ====================================================================
        /** add one more abscissa in the list 
         *  @param xnew value to be added 
         *  @return -1 if point is NOT added or new index for added points  
         *  @attention it can change the type of abscissas! 
         */
        int add  
        ( const double         xnew ) ;
        /** remove the point with the given index 
         *  @param index poitn with the index 
         *  @return true if point is really removed 
         *  @attention it can change the type of abscissas! 
         */
        bool remove 
        ( const unsigned short index ) ;
        // ====================================================================
      public: // make a slice for the given  ascissas 
        // ====================================================================
        /// make a slice fo the given abscissas 
        Abscissas slice ( const int i , const int j ) const ;
        // ====================================================================
      public:
        // ====================================================================
        // efficient swap to two abscissas 
        void exchange ( Abscissas& right ) 
        { 
          swap ( m_x     , right.m_x     ) ; 
          std::swap ( m_atype , right.m_atype ) ;
          std::swap ( m_xmin  , right.m_xmin  ) ;
          std::swap ( m_xmax  , right.m_xmax  ) ;
        }
        // ====================================================================
      private:
        // ====================================================================
        /* sort abscissas and eliminate the duplicates  
         * @return number of removed duplicates  
         */
        unsigned int sort()  ;
        // ====================================================================
      private:
        // ====================================================================
        /// (sorted) vector of abscissas 
        Data   m_x     {}          ; // (sorted) vector of abscissas 
        /// type of abscissas 
        AType  m_atype { Generic } ; // type of abscissas
        /// minimal abssissa
        double m_xmin  { 0 }       ; // minimal abssissa
        /// maximal abssissa
        double m_xmax  { 0 }       ; // maximal abssissa
        //  ====================================================================
      };
      // ======================================================================
      /// swap interpolation abscissas
      inline void swap ( Abscissas& a , Abscissas& b ) { a.exchange ( b ) ; }
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
        /** simple and efficient constructor from abscissas and y-list 
         *  @param x input vector of abscissas  
         *  @param ybegin start-iterator for sequence of function values 
         *  @param yend   end-iterator for sequence of function values 
         *  - if vector of y is longer  than vector x, extra values are ignored 
         *  - if vector of y is shorter than vector x, missing entries are assumed to be zero  
         */
        template <class ITERATOR , 
                  typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                  typename = std::enable_if<std::is_convertible<value_type,long double>::value> >
        Table 
        ( const Abscissas& x      , 
          ITERATOR         ybegin , 
          ITERATOR         yend   ) 
          : m_abscissas ( x ) 
          , m_values    ( x.size() , 0.0 ) 
        {
          const std::size_t nx = x.size() ;
          const std::size_t ny = std::distance ( ybegin , yend ) ;
          const std::size_t N  = std::min ( nx , ny ) ;
          //
          std::copy ( ybegin , ybegin + N , m_values.begin () ) ;
        }
        // ====================================================================
        /** simple & efficient constructor from abscissas and y-list 
         *  @param x input vector of abscissas  
         *  @param y input vector of y 
         *  - if vector of y is longer  than vector x, extra values are ignored 
         *  - if vector of y is shorter than vector x, missing entries are assumed to be zero  
         */
        Table
        ( const Abscissas&       x , 
          const Abscissas::Data& y ) 
          : Table ( x , y.begin() , y.end() )
        {}
        // ==================================================================== 
        /** simple constructor from x&y-lists 
         *  @param x input vector of x 
         *  @param y input vector of y 
         *  @param sorted  indicate that x-list is already sorted 
         *  - if vector of y is longer  than vector x, extra values are ignored 
         *  - if vector of y is shorter than vector x, missing entries are assumed to be zero  
         *  @attention It is not very efficient: 
         *   - there is sorting and duplicates removal 
         */
        Table
        ( const Abscissas::Data& x              , 
          const Abscissas::Data& y              , 
          const bool             sorted = false ) ;
        // ====================================================================
        /** simple contructor from x&y-lists 
         *  @param xbegin start  of sequence of abscissas 
         *  @param xend   end    of sequence of abscissas 
         *  @param ybegin start  of sequence of values 
         *  @param yend   end    of sequence of values 
         *  @param sorted indicate if data already sorted and duplicates removed  
         *  - if vector of y is longer  than vector x, extra values are ignored 
         *  - if vector of y is shorter than vector x, missing entries are assumed to be zero  
         *  @attention It is not very efficient: 
         *   - there is sorting and duplicates removal 
         */
        template <class XITERATOR, 
                  class YITERATOR,
                  typename xvalue_type = typename std::iterator_traits<XITERATOR>::value_type    ,
                  typename yvalue_type = typename std::iterator_traits<YITERATOR>::value_type    ,
                  typename = std::enable_if<std::is_convertible<xvalue_type,long double>::value> ,
                  typename = std::enable_if<std::is_convertible<yvalue_type,long double>::value> >
        Table
        ( XITERATOR  xbegin         ,
          XITERATOR  xend           ,
          YITERATOR  ybegin         ,
          YITERATOR  yend           ,
          const bool sorted = false )
          : m_abscissas () 
          , m_values    () 
        {  
          // =================================================================
          /// helper intermediate table 
          const int NNN = std::distance ( xbegin , xend ) ;
          const unsigned int NN = std::max ( 0 , NNN ) ;
          TABLE table { NN } ;
          /// 1) copy x into the table 
          std::transform 
            ( xbegin        , 
              xend          , 
              table.begin() , 
              [] ( double v ) { return std::make_pair ( v , 0.0 ) ; } );
          /// 2) copy y into the table 
          const unsigned int N = std::min
            ( std::distance ( xbegin , xend ) , 
              std::distance ( ybegin , yend ) ) ;
          std::transform
            ( table.begin ()     , 
              table.begin () + N , 
              ybegin             , 
              table.begin ()     , 
              [] ( const auto& p , const double v ) 
              { return std::make_pair ( p.first , v ) ; } );
          /// 3) sort it & eliminate duplicates 
          this->get_sorted ( sorted , table ) ;  
          /// 4) fill abscissas 
          m_abscissas = Abscissas 
            ( table.begin () , 
              table.end   () ,  
              []( const auto& p )-> double { return p.first ; } , true ) ;
          /// 5) fill values 
          m_values.resize( table.size () ) ;
          std::transform 
            ( table    .begin () ,
              table    .end   () , 
              m_values .begin () , 
              []( const auto& p )-> double { return p.second ; }       ) ;
          // ==================================================================
        }
        // ======================================================================
        /** the simplest constructor 
         *  @param data   input data 
         *  @param sorted indicate if data already  sorted
         *  @attention It is not very efficient: 
         *   - there is data copy, sorting and duplicates removal 
         */
        Table 
        ( TABLE        data           , 
          const bool   sorted = false ) ;
        // ====================================================================
        /** simple constructor from map/dictionary { x : f(x) }
         *  @param data map/dictionary with data  { x : f(x) }
         *  It is relatively efficient: no sorting
         *  - "numerical" duplicates are removed
         */
        template <class KEY, class VALUE>
        Table
        ( const std::map<KEY,VALUE>& data ) 
          : m_abscissas () 
          , m_values    () 
        {
          /// 1) create helper intermediate table 
          TABLE table { data.size() } ;
          TABLE::iterator i = table.begin() ;
          for ( const auto& item : data ) 
          { *i = std::make_pair ( item.first , item.second ) ; ++i ; }
          /// 2) eliminate duplicates, no need to sort it  
          this->get_sorted ( true , table ) ;
          /// 3) fill abscissas 
          m_abscissas = Abscissas 
            ( table.begin () , 
              table.end   () ,  
              []( const auto& p )-> double { return p.first ; } , true ) ;
          /// 4) fill values 
          m_values.resize( table.size () ) ;
          std::transform 
            ( table.begin    () ,
              table.end      () , 
              m_values.begin () ,
              []( const auto& p )-> double { return p.second ; }       ) ;
          // ==================================================================
        }
        // ====================================================================
        /** templated constructor 
         *  (very efficient: no sorting, no removal of duplicate...)
         *  @param fun       function object 
         *  @param abscissas interpolation abscissas 
         */
        template <class FUNCTION>
        Table 
        ( const Abscissas& a     , 
          FUNCTION         fun   )                
          : m_abscissas ( a              )  
          , m_values    ( a.size() , 0.0 )
        { 
          std::transform ( a.begin () , a.end   () , m_values.begin () , 
                           [&fun] ( const double x ) -> double { return fun ( x ) ; } ) ; 
        }
        // ====================================================================
        /** templated constructor 
         *  (very efficient: no sorting, no removal of duplicate...)
         *  @param fun       function object 
         *  @param abscissas interpolation abscissas 
         */
        template <class FUNCTION>
        Table 
        ( FUNCTION         fun   ,
          const Abscissas& a     )
          : m_abscissas ( a              )  
          , m_values    ( a.size() , 0.0 )
        { 
          std::transform ( a.begin () , a.end   () , m_values.begin () , 
                           [&fun] ( const double x ) -> double { return fun ( x ) ; } ) ; 
        }
        // ====================================================================
        /** templated constructor 
         *  @param begin start  of x-sequence 
         *  @param end   end    of x-sequence 
         *  @param fun   function object 
         *  @attention duplicated abscissas will be removed 
         */
        template <class ITERATOR , 
                  class FUNCTION ,
                  typename value_type = typename std::iterator_traits<ITERATOR>::value_type,
                  typename = std::enable_if<std::is_convertible<value_type,long double>::value> >
        Table 
        ( ITERATOR begin , 
          ITERATOR end   ,
          FUNCTION fun   )  
          : Table ( Abscissas ( begin, end ) , fun ) 
        {}
        // ====================================================================
        /** templated constructor 
         *  @param fun  function object 
         *  @param n    number of interpolation points 
         *  @param low  low edge of interpolation region 
         *  @param high high edge of interpolation region 
         *  @param t    interpolation type 
         */
        template <class FUNCTION>
        Table 
        ( FUNCTION               fun  , 
          const unsigned short   n    , 
          const double           low  ,   
          const double           high , 
          const Abscissas::AType t    ) 
          : Table ( Abscissas ( n , low , high , t ) , fun ) 
        {}
        // ====================================================================
        /// default constructor
        Table () = default ;
        /// copy constructor 
        Table ( const Table&  ) = default ;
        /// move constructor 
        Table (       Table&& ) = default ;
        // ====================================================================
      public: // creator for interpolation table 
        // ====================================================================
        /// creator for the interpolation table 
        template <class FUNCTION, typename ...Args>
        static Table 
        create ( FUNCTION fun  , 
                 Args ... args ) { return Table ( fun , args... ) ; } 
        // ====================================================================
      public:
        // ====================================================================
        /// number of interpolation points 
        unsigned int size    () const { return m_abscissas.size  () ; }
        /// no interpolation points 
        bool         empty   () const { return m_abscissas.empty () ; }
        // ====================================================================
      public : // 
        // ====================================================================
        /// get the entry for given index 
        std::pair<double,double> operator[] ( const unsigned short index ) const 
        { return std::make_pair ( x ( index ) , y ( index ) ) ;  }
        // ====================================================================
      public:
        // ====================================================================
        ///  get the abscissas for the given index 
        double x    ( const unsigned short index ) const 
        { return m_abscissas.x ( index ) ; }
        ///  get the value    for the given index 
        double y    ( const unsigned short index ) const 
        { 
          static const double s_nan = std::numeric_limits<double>::quiet_NaN () ;
          return 
            m_values.empty()        ? s_nan :
            index < m_values.size() ? m_values[ index ] : m_values.back()  ;
        }
        // =====================================================================
      public: // get abscissas & data 
        // ====================================================================
        /// get abscissas      
        const Abscissas&       abscissas () const { return m_abscissas ; }
        /// get function values 
        const Abscissas::Data& values    () const { return m_values    ; }
        // =====================================================================
        /// get abscissas type 
        Abscissas::AType       atype     () const { return m_abscissas.atype()  ; }
        // ====================================================================
      public:
        // ====================================================================
        /// minimal  abscissa 
        double xmin () const { return m_abscissas.xmin () ; }
        /// maximal abscissas 
        double xmax () const { return m_abscissas.xmax () ; }
        // ====================================================================
      public: // iterators 
        // ====================================================================
        Abscissas::const_iterator       x_begin () const { return m_abscissas.begin () ; }
        Abscissas::const_iterator       x_end   () const { return m_abscissas.end   () ; }
        Abscissas::Data::const_iterator y_begin () const { return m_values   .begin () ; }
        Abscissas::Data::const_iterator y_end   () const { return m_values   .end   () ; }        
        // ====================================================================
      public:
        // ====================================================================
        /** add the point (x,y) into interpolation table 
         *  @param x abscissas of the point to be added 
         *  @param y the value of function at x 
         *  @return the index of new point, or -1  if point if not added 
         */
        int add 
        ( const  double x ,   
          const  double y ) ;
        // ====================================================================
        /** remove the point with the  given index 
         *  @param index the point to be removed 
         *  @return   true if point is removed 
         */
        bool remove 
        ( unsigned short index ) ;
        // ====================================================================
      public: // make interpolation 
        // ====================================================================
        /** interpolation using the straightforward Lagrange interpolant 
         *  https://en.wikipedia.org/wiki/Lagrange_polynomial
         *  - it is rather slow O(n^2)
         */
        double lagrange 
        ( const double x ) const ;
        // ====================================================================
        /** interpolation using Neville's algorithm
         *  @see https://en.wikipedia.org/wiki/Neville%27s_algorithm
         *  - it is rather slow O(n^2)
         */
        double neville 
        ( const double x ) const ;
        // ====================================================================
        /** interpolation using the 1st rational Berrut's interpolant 
         *  @see Kai Hormann, "Barycentric interpolation", 
         *  @see https://www.inf.usi.ch/hormann/papers/Hormann.2014.BI.pdf            
         */
        double berrut1st 
        ( const double x ) const ;
        // =====================================================================
        /** interpolation using the 2nd rational Berrut's interpolant 
         *  @see Kai Hormann, "Barycentric interpolation", 
         *  @see https://www.inf.usi.ch/hormann/papers/Hormann.2014.BI.pdf            
         */
        double berrut2nd  
        ( const double x ) const ;        
        // ====================================================================
      public: // interpolation with derivatives:
        // ====================================================================
        /** Interpolation using Neville's algorithm
         *  @see https://en.wikipedia.org/wiki/Neville%27s_algorithm
         *  - it is rather slow O(n^2)
         *  @return pair \f$(y(x),\frac{dy}{dx})\f$ 
         */
        std::pair<double,double> 
        neville2  
        ( const double       x ) const ;
        // ====================================================================
        /** Simple lagrange interpolation 
         *  - it also evaluates the derivative with respect to \f$ y_i\f$  
         *  @param x interpolation point
         *  @param iy index of y_i
         *  @return pair of \f$ (y(x),\frac{dy}{dy_i})\f$ 
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
      protected :
        // ====================================================================
        /// swap two interpolation tables 
        void exchange ( Table& right ) 
        { 
          swap ( m_abscissas , right.m_abscissas ) ;
          swap ( m_values    , right.m_values    ) ;
        }
        // ====================================================================
      public:
        // ====================================================================
        ///  sort table and remove duplicated abscissas 
        void get_sorted ( bool sorted , TABLE& table  ) ;
        // ====================================================================
      private :
        // ====================================================================
        /// interpolation absciccas 
        Abscissas       m_abscissas {} ; // interpolation absciccas 
        /// function values 
        Abscissas::Data m_values    {} ; // interpolation absciccas 
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
    } //                            end of namespace Ostap::Math::Interpolation 
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
