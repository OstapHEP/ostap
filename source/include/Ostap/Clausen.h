// ============================================================================
#ifndef OSTAP_CLAUSEN_H 
#define OSTAP_CLAUSEN_H 1
// ============================================================================
// Include files
// ============================================================================
// STD& STL 
// ============================================================================
#include <array>
#include <cmath>
#include <numeric>
// ============================================================================
// Ostap 
// ============================================================================ 
#include "Ostap/Clenshaw.h"
#include "Ostap/MakeArray.h"
#include "Ostap/Math.h"
// ============================================================================
/** @file Ostap/Clausen.h
 *  Helper structures to implement Clausen functions 
 *  @see https://en.wikipedia.org/wiki/Clausen_function
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2015-03-26
 */  
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    namespace Clausen
    {
      // ======================================================================      
      /** Limit the number of of terms for explicit Fourier series: 
       *  \f[ \sum_{k=1}^{L} \frac{ \sin kx }{k^N} \f]
       *  \f[ \sum_{k=1}^{L} \frac{ \cos kx }{k^N} \f]
       */
      template <unsigned int N>
      struct Len_max
      {
        enum : unsigned int { value =
			      N  <= 3 ? 65535u :
			      4  == N ? 11000u :
			      5  == N ?  1600u :
			      6  == N ?   600u :
			      7  == N ?   300u :
			      8  == N ?   200u :
			      9  == N ?    80u :
			      10 == N ?    50u :
			      11 == N ?    40u :
			      12 == N ?    35u :
			      13 == N ?    25u :
			      14 == N ?    20u :
			      15 == N ?    18u : 15u } ;
      };
      // ======================================================================
      /** S-sum: Explicit Fourrier' sine sum
       *  \f[ S_n ( x ) \sum_{k=1}^{K} \frac { \sin k x }{k^n}  \f] 
       *  @see https://en.wikipedia.org/wiki/Clausen_function
       *  - Explicit specializations are provided for small n
       *  - Explicit Fourrier' sine  sum for "large" n 
       */
      template <unsigned int N>
      struct S_
      {
        // =====================================================================
      public: 
        // =====================================================================
        enum { L = Len_max<N>::value } ;
        // =====================================================================
      public:
        // =====================================================================        
        inline double operator () ( const double x ) const
        { return Ostap::Math::Clenshaw::sine_sum ( s_SK.begin() , s_SK.end() , x ) ; }
        // ==================================================================        
      protected:
        // =====================================================================        
        /// the array of 1/k^{N} : 
        static const std::array<long double,L> s_SK ;
        // =====================================================================    
      } ;
      // =======================================================================
      /** C-sum: Explicit Fourrier' cosine sum
       *  \f[ S_n ( x ) \sum_{k=1}^{K} \frac { \sin k x }{k^n}  \f] 
       *  @see https://en.wikipedia.org/wiki/Clausen_function
       *  - Explicit specializations are provided for small n
       *  - Explicit Fourrier' sine  sum for "large" n 
       */
      template <unsigned int N>
      struct C_
      {
        // =====================================================================
      public: 
        // =====================================================================
        enum { L = Len_max<N>::value } ;
        // =====================================================================
      public:
        // =====================================================================        
        inline double operator () ( const double x ) const
        { return Ostap::Math::Clenshaw::cosine_sum ( s_CK.begin() , s_CK.end() , x ) ; }
        // ==================================================================        
      protected:
        // =====================================================================        
        /// the array of 1/k^{N} : 
        static const std::array<long double,L> s_CK ;
        // =====================================================================    
      } ;
      // =======================================================================
      template <unsigned int N>
      const std::array<long double,S_<N>::L>
      S_<N>::s_SK { Ostap::Math::make_array
          ( [] ( const std::size_t  k ) -> long double  
          { return std::pow ( 1.0L / ( k + 1 ) , N ) ; } ,
            std::make_index_sequence<S_<N>::L> () ) } ;      
      // =======================================================================
      template <unsigned int N>
      const std::array<long double,C_<N>::L>
      C_<N>::s_CK { Ostap::Math::make_array
          ( [] ( const std::size_t k ) -> long double  
          { return 0 == k ? 0.0 : std::pow ( 1.0L / k , N ) ; } ,
            std::make_index_sequence<C_<N>::L> () ) } ;      
      // ======================================================================
      // Specific cases: 
      double S0  ( const double x ) ; // explicit 
      double S1  ( const double x ) ; // explicit 
      double S2  ( const double x ) ; // explicit: "standard" Clausen function: use GSL 
      double S3  ( const double x ) ; // explicit
      double S4  ( const double x ) ; // explicit -> Polylogarithm 
      double S5  ( const double x ) ; // explicit
      double S6  ( const double x ) ; // explicit -> Polylogarithm       
      double S7  ( const double x ) ; // explicit -> Bernulli polynomials
      double S8  ( const double x ) ; // explicit -> polylogarithm            
      double S9  ( const double x ) ; // explicit -> Bernulli polynomials
      double S10 ( const double x ) ; // explicit -> polylogarithm            
      // ======================================================================
      // specific case for N=0
      template <> struct S_<0>  { inline double operator () ( const double x ) const { return S0 ( x ) ; } } ; 
      // =================================================================
      // specific case for N=1
      template <> struct S_<1>  { inline double operator () ( const double x ) const { return S1 ( x ) ; } } ; 
      // =================================================================
      // specific case for N=2
      template <> struct S_<2>  { inline double operator () ( const double x ) const { return S2 ( x ) ; } } ; 
      // =================================================================
      // specific case for N=3
      template <> struct S_<3>  { inline double operator () ( const double x ) const { return S3 ( x ) ; } } ; 
      // ======================================================================
      // specific case for N=4
      template <> struct S_<4>  { inline double operator () ( const double x ) const { return S4 ( x ) ; } } ;
      // ==========================================================================================
      // specific case for N=5
      template <> struct S_<5>  { inline double operator () ( const double x ) const { return S5 ( x ) ; } } ;
      // ==========================================================================================
      // specific case for N=6
      template <> struct S_<6>  { inline double operator () ( const double x ) const { return S6 ( x ) ; } } ;
      // ==========================================================================================
      // specific case for N=7
      template <> struct S_<7>  { inline double operator () ( const double x ) const { return S7 ( x ) ; } } ;
      // ==========================================================================================
      // specific case for N=8
      template <> struct S_<8>  { inline double operator () ( const double x ) const { return S8 ( x ) ; } } ;
      // ==========================================================================================      
      // specific case for N=9
      template <> struct S_<9>  { inline double operator () ( const double x ) const { return S9 ( x ) ; } } ;
      // specific case for N=10
      template <> struct S_<10> { inline double operator () ( const double x ) const { return S10 ( x ) ; } } ;
      // ==========================================================================================
      // Specific cases: 0, 1 , 2 & 4
      // ======================================================================
      double C0  ( const double x ) ; // explicit 
      double C1  ( const double x ) ; // explicit 
      double C2  ( const double x ) ; // explicit
      double C3  ( const double x ) ; // explicit -> polylogarithm       
      double C4  ( const double x ) ; // explicit
      double C5  ( const double x ) ; // explicit -> polylogarithm             
      double C6  ( const double x ) ; // explicit -> Bernulli polynomials
      double C7  ( const double x ) ; // explicit -> polylogarithm                   
      double C8  ( const double x ) ; // explicit -> Bernulli polynomials
      double C9  ( const double x ) ; // explicit -> polylogarithm                         
      double C10 ( const double x ) ; // explicit -> Bernulli polynomials 
      // ======================================================================
      // specific case for N=0
      template <> struct C_<0>  { inline double operator () ( const double x ) const { return C0  ( x ) ; } } ; 
      // specific case for N=1
      template <> struct C_<1>  { inline double operator () ( const double x ) const { return C1  ( x ) ; } } ; 
      // ======================================================================
      // specific case for N=2
      template <> struct C_<2>  { inline double operator () ( const double x ) const { return C2  ( x ) ; } } ; 
      // ======================================================================
      // specific case for N=3
      template <> struct C_<3>  { inline double operator () ( const double x ) const { return C3  ( x ) ; } } ; 
      // ======================================================================
      // specific case for N=4
      template <> struct C_<4>  { inline double operator () ( const double x ) const { return C4  ( x ) ; } } ;      
      // ======================================================================
      // specific case for N=5
      template <> struct C_<5>  { inline double operator () ( const double x ) const { return C5  ( x ) ; } } ;      
      // ======================================================================
      // specific case for N=6
      template <> struct C_<6>  { inline double operator () ( const double x ) const { return C6  ( x ) ; } } ;      
      // ======================================================================
      // specific case for N=7
      template <> struct C_<7>  { inline double operator () ( const double x ) const { return C7  ( x ) ; } } ;      
      // ======================================================================
      // specific case for N=8
      template <> struct C_<8>  { inline double operator () ( const double x ) const { return C8  ( x ) ; } } ;            
      // ======================================================================
      // specific case for N=9
      template <> struct C_<9>  { inline double operator () ( const double x ) const { return C9  ( x ) ; } } ;      
      // ======================================================================
      // specific case for N=10
      template <> struct C_<10> { inline double operator () ( const double x ) const { return C10 ( x ) ; } } ;      
      // ======================================================================
      template <unsigned int N, bool> struct SC_ ;      
      template <unsigned int N>       struct SC_<N,true>  : public S_<N> {} ;
      template <unsigned int N>       struct SC_<N,false> : public C_<N> {} ; 
      // ======================================================================
      /** Standard Clausen functions
       *  \f[ \begin{array}{lcl}
       *      Cl_{ 2m + 2 } ( x ) & = & \sum_k \frac{ \sin kx }{k^{ 2m + 2 }} \\ 
       *      Cl_{ 2m + 1 } ( x ) & = & \sum_k \frac{ \cos kx }{k^{ 2m + 1 }} \\ 
       *      Sl_{ 2m + 2 } ( x ) & = & \sum_k \frac{ \cos kx }{k^{ 2m + 2 }} \\ 
       *      Sl_{ 2m + 1 } ( x ) & = & \sum_k \frac{ \sin kx }{k^{ 2m + 1 }}
       *      \end{array}  \f] 
       */
      template <unsigned int N> struct Cl_ : public SC_<N,0==N%2> {} ;  
      template <unsigned int N> struct Sl_ : public SC_<N,1==N%2> {} ;       
      // ======================================================================
      /// \f$ \sum_{i=1} \frac{ \sin kx }{ k^n }\f$ 
      double S
      ( const unsigned int n ,
        const double       x ) ;
      // ======================================================================
      /// \f$ \sum_{i=1} \frac{ \cos kx }{ k^n }\f$ 
      double C
      ( const unsigned int n ,
        const double       x ) ;      
      // ======================================================================
      /** Generalized Clausen' function
       *  \f[ S_s(x) = \Im Li ( s , \mathrm{e}{ix} ) = \sum \frac{ \sin kx}{k^s}\f]
       */
      double S
      ( const double s ,
        const double x ) ;
      // ======================================================================
      /** Generalized Clausen' function
       *  \f[ C_s(x) = \Re Li ( s , \mathrm{e}^{ix}x ) = \sum \frac{ \cos kx}{k^s} \f] 
       */
      double C
      ( const double s ,
        const double x ) ;      
      // ======================================================================      
    } //                              The end of namespace Ostap::Math::Clausen
    // ========================================================================
    /** get the standard Clausen function \f$ Cl_2(x) \f$  using GSL
     *  @param x argument 
     *  @return value of clausen fnction \f$ Cl_2 \f$
     *  @see https://en.wikipedia.org/wiki/Clausen_function    
     */
    double clausen ( const double x ) ;
    // ========================================================================
    /** Standard Clausen functions
     *  \f[ \begin{array}{lcc}
     *      Cl_{2m+2} ( x ) & = \sum_k \frac{ \sin kx }{k^{2m+2}}& \\ 
     *      Cl_{2m+1} ( x ) & = \sum_k \frac{ \cos kx }{k^{2m+1}}& 
     *      \end{array}   \f] 
     * @see https://en.wikipedia.org/wiki/Clausen_function    
     */
    double Cl
    ( const unsigned int n , 
      const double       x ) ;
    // ========================================================================
    /** Standard Clausen functions, aka Gleisher-Clausen functions  
     *  \f[ \begin{array}{lcc}
     *      Sl_{2m+2} ( x ) & = \sum_k \frac{ \cos kx }{k^{2m+2}}& \\ 
     *      Sl_{2m+1} ( x ) & = \sum_k \frac{ \sin kx }{k^{2m+1}}& 
     *      \end{array}   \f] 
     * The function are related to Bernulli' polynomials 
     * @see https://en.wikipedia.org/wiki/Clausen_function    
     */
    double Sl
    ( const unsigned int  n , 
      const double        x ) ;    
    // ========================================================================
  } //                                         The end of namespace Ostap::Math 
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_CLAUSEN_H
// ============================================================================
