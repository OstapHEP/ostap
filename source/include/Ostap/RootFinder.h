// ============================================================================
#ifndef OSTAP_ROOTFINDER_H 
#define OSTAP_ROOTFINDER_H 1
// ============================================================================
// Include files 
// ============================================================================
// STD&STL
// ============================================================================
#include <cmath>
#include <functional>
#include <utility>
#include <string>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatusCode.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    /** @class Integrator Ostap/Integrator.h 
     *  simple numerical integrator for 1D&,2D&3D-cases 
     */
    class RootFinder 
    {
    public:
        // ===================================================================
        /// the actual type of function
        typedef std::function<double(double)>               function1 ;
        // ===================================================================
      public:
        // ===================================================================
        class Point
        { 
            // ===============================================================
        public:
            // ===============================================================
            Point
            ( const double x  ,  
              const double fx )
            : m_x  { x  }
            , m_fx ( fx )
            {}
            // ===============================================================
            Point () = delete ;
            // ===============================================================
        public:
            // ===============================================================
            /// x-value 
            inline double x  () const { return m_x  ; }
            /// function value
            inline double fx () const { return m_fx ; }
            // ===============================================================
        public:
            // ===============================================================
            inline void swap ( Point& other )
            {
                std::swap ( m_x  , other.m_x  ) ;
                std::swap ( m_fx , other.m_fx ) ; 
            }
            // ===============================================================
            inline bool operator<( const Point& other ) const
            {
              return ( m_x < other.m_x ) ? true  :
                     ( m_x > other.m_x ) ? false : m_fx < other.m_fx ; 
            }
            // ===============================================================
        private :
            // ===============================================================
            /// the point 
            double m_x  ;  // the point 
            /// funcnton value 
            double m_fx ;  // function value    
            // ===============================================================
        } ; 
        // ===================================================================
    public:
        // ===================================================================
        /** constructor from full configuration 
         *  @param max_calls  maximal number of function calls
         *  @param froot      consider x as root if  \f$ \left| f(x) \right| < f_{root} \f$ if \f$ 0 < f_{root}\f$ 
         *  @param atolerance absolute tolerance 
         *  @param rtolerance relative tolerance  
         */
        RootFinder 
        ( const unsigned short max_calls   =  100 ,
          const double         froot       = -1   ,
          const double         atolerance  =  0   , 
          const double         rtotlerance =  0   ) ;
        // ===================================================================
    public:
        // ===================================================================
        /// current number of function/derivative calls 
        inline std::size_t ncalls () const { return m_ncalls ; }
        // ===================================================================
    public:
        // ===================================================================
        /// maximal umber of function/derivatives calls 
        inline std::size_t max_calls  () const { return m_maxcalls    ; }
        inline double      froot      () const { return m_froot       ; }
        /// abolute tolerance 
        inline double      atolerance () const { return m_atolerance  ; }
        /// relative tolerance 
        inline double      rtolerance () const { return m_rtolerance  ; } 
        // ===================================================================
    public:
        // ===================================================================
        template <class FUNCTION>
        inline Ostap::StatusCode
        root 
        ( const FUNCTION& fun ,
          double&         r   ,
          double&         a   , 
          double&         b   ) const 
          {
            const auto cfun = std::cref ( fun ) ;
            return root ( cfun , r , a , b , function1 () , function1 () ) ;
          } 
        // =====================================================================
        template <class FUNCTION, 
                  class DERIVATIVE>
        inline Ostap::StatusCode
        root 
        ( const FUNCTION&   fun        ,
          const DERIVATIVE& derivative , 
          double&           r          ,
          double&           a          , 
          double&           b          ) const 
          {
            const auto cfun = std::cref ( fun        ) ;
            const auto cder = std::cref ( derivative ) ;
            return root ( cfun , r , a , b , cder , function1 () ) ;
          } 
        // =====================================================================
        template <class FUNCTION, 
                  class DERIVATIVE1,
                  class DERIVATIVE2>
        inline Ostap::StatusCode
        root 
        ( const FUNCTION&    fun         ,
          const DERIVATIVE1& derivative1 ,
          const DERIVATIVE2& derivative2 , 
          double&            r           ,
          double&            a           , 
          double&            b           ) const 
          {
            const auto cfun  = std::cref ( fun         ) ;
            const auto cder1 = std::cref ( derivative1 ) ;
            const auto cder2 = std::cref ( derivative2 ) ;  
            return root ( cfun , r , a , b , cder1 , cder2 ) ;
          } 
        // ===================================================================
    public:
        // ===================================================================
        // ===================================================================
    public:
        // ===================================================================
        /** find a root in [a,b]
         *  @param fun the function
         *  @param r   (update) the root 
         *  @param a   (update) the bracket interval
         *  @param b   (update) the bracket innteval
         *  @param derivative1 the first derivative (optional)
         *  @param derivative1 the second  derivative (optional)
         *  @return status code 
         */
        Ostap::StatusCode  
        root 
        ( function1 fun                        , 
          double&   r                          , 
          double&   a                          , 
          double&   b                          ,
          function1 derivative1 = function1 () , 
          function1 derivative2 = function1 () ) const ;
        // ====================================================================
    protected:
        // ====================================================================
        /** find a root throgh the sequence of steps/iterations
         *  @param fun the function
         *  @param derivative1 the first derivative 
         *  @param derivative1 the second derivative 
         *  @param r   (update) the root 
         *  @param a   (update) the bracket interval
         *  @param b   (update) the bracket innteval
         *  @return status code 
         */
        Ostap::StatusCode 
        root 
        ( function1 fun                   , 
          Point&    r                     , 
          Point&    a                     ,
          Point&    b                     , 
          function1 deriv1 = function1 () , 
          function1 deriv2 = function1 () ) const ;
        // =====================================================================
        /** Single step. It sequentially apply several methods 
         *  @param fun the function
         *  @param derivative1 the first derivative 
         *  @param derivative1 the second derivative 
         *  @param r   (update) the root 
         *  @param a   (update) the bracket interval
         *  @param b   (update) the bracket innteval
         *  @return status code 
         */
        Ostap::StatusCode 
        step 
        ( function1 fun                   , 
          Point&    r                     , 
          Point&    a                     ,
          Point&    b                     , 
          function1 deriv1 = function1 () , 
          function1 deriv2 = function1 () ) const ;
        // ====================================================================
    private :
        // ====================================================================
        ///  
        std::size_t          m_maxcalls   {  100   } ; 
        double               m_froot      { -1     } ; 
        double               m_atolerance {  1.e-9 } ; 
        double               m_rtolerance {  1.e-9 } ;        
        /// current number of calls 
        mutable std::size_t  m_ncalls     {  0     } ;
        // ====================================================================
    public:
        // ====================================================================
        /// status code if maximal @call is reached 
        static const Ostap::StatusCode  NumCallsLimit ;
        // ====================================================================
    };  
    // ========================================================================
    /// swap two points 
    inline void swap 
    ( RootFinder::Point& a , 
      RootFinder::Point& b ) { a.swap ( b ) ; }
    // ========================================================================
  } //                                          The end of namepace Ostap::Math
  // ==========================================================================
} // The end of namespace Ostap
// ============================================================================
#endif  //  OSTAP_ROOTFINDER_H
// ============================================================================