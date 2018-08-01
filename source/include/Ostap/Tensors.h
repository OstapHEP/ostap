// ===============================================================
#ifndef OSTAP_TENSORS_H
#define OSTAP_TENSORS_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <cmath>
#include <complex>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Vector4DTypes.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    /** @namespace Tensors
     *  Collection of general purpose tensors and operations with them
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2008-0725
     */
    namespace Tensors
    {
      // ======================================================================
      /** @enum Indices
       *  The list of Lorentz indices
       *  The numbers are in accordance to
       *  ROOT::Math::LorentzVector
       *  @see ROOT::Math::Lorentz::Vector
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2008-0725
       */
      enum {
        //
        X    = 0 ,
        Y    = 1 ,
        Z    = 2 ,
        T    = 3 ,
        //
        PX   = X ,
        PY   = Y ,
        PZ   = Z ,
        E    = T ,
        //
        LAST = 4
      } ;
      // ======================================================================
      /** @struct Delta_
       *
       *  (Compile-time) Kronecker delta: \f$ \delta^{\mu}_{\nu} \f$
       *
       *  @code
       *
       *    int d12 = Delta_<1,2>::value ;
       *
       *  @endcode
       *
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2008-0725
       */
      template <unsigned int I, unsigned int J>
      struct Delta_      { enum { value =  0 } ; } ;
      /// the proper template specialization for diagonal elements
      template <unsigned int I>
      struct Delta_<I,I> { enum { value =  1 } ; } ;
      // ======================================================================
      /** @struct  Delta
       *  Kronecker delta: \f$ \delta^i_j \f$
       *
       *  @code
       *
       *   Delta delta ;
       *
       *    const size_t i = ... ;
       *    const size_t j = ... ;
       *
       *    int ij = delta( i , j ) ;
       *
       *  @endcode
       *
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2008-0725
       */
      struct Delta
      {
        // ====================================================================
        /// Kroneker delta
        inline int delta
        (  const  size_t i ,
           const  size_t j ) const { return i == j ; }
        /// Kroneker delta (functional form)
        inline int operator()
        ( const size_t i ,
          const size_t j ) const { return i == j ; }
        // ====================================================================
      } ;
      // ======================================================================
      /** @struct G_
       *
       *  (Compile-time) Minkowski metric:
       * \f$  g_{\mu\nu} =
       *   \begin{pmatrix}
       *      -1 &  0  &  0  & 0 \\
       *       0 & -1  &  0  & 0 \\
       *       0 &  0  & -1  & 0 \\
       *       0 &  0  &  0  & 1
       * \end{pmatrix}
       * \f$
       *
       *  The metric has been choosen in the way (-,-,-,+) in
       *   accordance with ROOT::Math::LorentzVector
       *
       *  @code
       *
       *    const int gXX = G_<X,X>::value ;
       *
       *  @endcode
       *
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2008-0725
       */
      template <unsigned int I, unsigned int J>
      struct G_        : std::integral_constant<int, 0> {}; 
      /// the proper template specialzation for diagonal elements
      template <unsigned int I>
      struct G_<I,I>   : std::integral_constant<int,-1> {}; 
      /// the proper template specialization for time component
      template <>
      // struct G_<T,T>   : std::integral_constant<int, 1> {};   
      // struct G_<E,E>   : std::integral_constant<int, 1> {};   
      struct G_<3,3>   : std::integral_constant<int, 1> {};   
      // ======================================================================
      /** struct G
       *
       *  Minkowski metric:
       * \f$  g_{\mu\nu} =
       *   \begin{pmatrix}
       *      -1 &  0  &  0  & 0 \\
       *       0 & -1  &  0  & 0 \\
       *       0 &  0  & -1  & 0 \\
       *       0 &  0  &  0  & 1
       * \end{pmatrix}
       * \f$
       *
       *  The metric has been choosen in the way (-,-,-,+) in
       *  accordance with ROOT::Math::LorentzVector
       *
       *  @code
       *
       *    G g ;
       *
       *    int g00 = g(0,0) ;
       *
       *  @endcode
       *
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2008-0725
       */
      struct G
      {
        // ====================================================================
        /// the only one important function: get the metric
        inline int  operator ()
        ( const size_t i ,
          const size_t j ) const { return g( i, j ) ; }
        /// the only one important function: get the metric
        inline int  g
        ( const size_t i ,
          const size_t j ) const
        {
          return
            ( i    != j )  ?  0 :
            ( LAST <= i )  ?  0 :
            ( T    == i )  ?  1 : -1 ;
        }
        // ====================================================================
      };
      // ======================================================================
      /** @struct Epsilon_
       *
       *  (Compile-time) 4D Antisymmetric Levi-Civita tensor
       *   \f$ \epsilon_{\mu\nu\lambda\delta} \f$
       *
       *  Convention here: \f$ \epsilon_{0123} = \epsilon_{XYZT} = 1 \f$
       *
       *  @code
       *
       *   int result = Epsilon_<X,Y,Z,T>::value ;
       *
       *  @endcode
       *
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2008-0725
       */
      template <unsigned int I ,
                unsigned int J ,
                unsigned int K ,
                unsigned int L > struct Epsilon_ ;
      // ======================================================================
      /// stopping criteria for compile-time recursion
      template <unsigned int I, unsigned int J, unsigned int K>
      struct Epsilon_<I,I,J,K> : std::integral_constant<int, 0> {};
      /// stopping criteria for compile-time recursion
      template <unsigned int I, unsigned int J, unsigned int K>
      struct Epsilon_<I,J,I,K> : public Epsilon_<I,I,J,K> {} ;
      /// stopping criteria for compile-time recursion
      template <unsigned int I, unsigned int J, unsigned int K>
      struct Epsilon_<I,J,K,I> : public Epsilon_<I,I,J,K> {} ;
      /// stopping criteria for compile-time recursion
      template <unsigned int I, unsigned int J, unsigned int K>
      struct Epsilon_<J,I,K,I> : public Epsilon_<I,I,J,K> {} ;
      /// stopping criteria for compile-time recursion
      template <unsigned int I, unsigned int J, unsigned int K>
      struct Epsilon_<J,I,I,K> : public Epsilon_<I,I,J,K> {} ;
      /// stopping criteria for compile-time recursion
      template <unsigned int I, unsigned int J, unsigned int K>
      struct Epsilon_<J,K,I,I> : public Epsilon_<I,I,J,K> {} ;
      // ======================================================================
      /// stopping criteria for compile-time recursion
      template <unsigned int I, unsigned int J>
      struct Epsilon_<I,I,J,J> : std::integral_constant<int,0>{};
      /// stopping criteria for compile-time recursion
      template <unsigned int I, unsigned int J>
      struct Epsilon_<I,J,I,J> : public Epsilon_<I,I,J,J> {} ;
      /// stopping criteria for compile-time recursion
      template <unsigned int I, unsigned int J>
      struct Epsilon_<I,J,J,I> : public Epsilon_<I,I,J,J> {} ;
      // ======================================================================
      /// stopping criteria for compile-time recursion
      template <unsigned int I, unsigned int J>
      struct Epsilon_<I,I,I,J> : std::integral_constant<int,0>{};
      template <unsigned int I, unsigned int J>
      struct Epsilon_<I,I,J,I> : public Epsilon_<I,I,I,J> {} ;
      /// stopping criteria for compile-time recursion
      template <unsigned int I, unsigned int J>
      struct Epsilon_<I,J,I,I> : public Epsilon_<I,I,I,J> {} ;
      /// stopping criteria for compile-time recursion
      template <unsigned int I, unsigned int J>
      struct Epsilon_<J,I,I,I> : public Epsilon_<I,I,I,J> {} ;
      // ======================================================================
      /// stopping criteria for compile-time recursion
      template <unsigned int I>
      struct Epsilon_<I,I,I,I> : std::integral_constant<int,0>{};
      // ======================================================================
      /// stopping criteria for compile-time recursion
      template <>
      // struct Epsilon_<X,Y,Z,T> : std::integral_constant<int,1>{};
      // struct Epsilon_<X,Y,Z,E> : std::integral_constant<int,1>{};
      struct Epsilon_<0,1,2,3> : std::integral_constant<int,1>{};
      // ======================================================================
      namespace detail
      {
        // ====================================================================
        /// helper structure for conditional selection
        template <bool C,int VALUE,class TYPE>
        struct _Value : std::integral_constant<int,VALUE>{};
        // ====================================================================
        /// helper structure for conditional selection
        template <int VALUE,class TYPE>
        struct _Value<false,VALUE,TYPE> : std::integral_constant<int, -TYPE::value>{};
        // ====================================================================
      }
      //  =====================================================================
      /// the generic evaluation of Levi-Chivita elements
      template <unsigned int I , unsigned int J ,
                unsigned int K , unsigned int L >
      struct Epsilon_
      {
      private:
        // ====================================================================
        // helper types for permutations
        typedef Epsilon_<J,I,K,L>  _12 ; // permute the 1st and the 2nd index
        typedef Epsilon_<I,K,J,L>  _23 ; // permute the 2nd and the 3rd index
        typedef Epsilon_<I,J,L,K>  _34 ; // permute the 3rd and the 4th index
        // ====================================================================
      private:
        // ====================================================================
        /// valid indices?
        enum { valid = ( I != J ) && ( J != K ) && ( K != L )  && ( L != I ) } ;
        // ====================================================================
      public:
        // ====================================================================
        /** compile-time recursion here!
         *  The most important lines: the actual compile-evaluation
         *  of Levi-Civita symbols
         */
        enum { value = !valid ? 0 :
               detail::_Value<(I<J),
               detail::_Value<(J<K),
               detail::_Value<(K<L),(L<LAST)?1:0,_34>::value,_23>::value,_12>::value } ;
        // ====================================================================
      };
      // ======================================================================
      /** @struct Epsilon1_
       *  (Compile-time) evaluation of the tensor product of
       *  two Levi-Civita symbols:
       *
       *  The following identity has been used:
       * \f[
       *  \alpha^{IJK}_{LMN} =
       *  \epsilon^{IJK\kappa}
       *  \epsilon_{LMN\kappa} = -
       *  \begin{Vmatrix}
       *    \delta^{I}_{L} & \delta^{I}_{M} & \delta^{I}_{N} \  \
       *    \delta^{J}_{L} & \delta^{J}_{M} & \delta^{J}_{N}  \ \
       *    \delta^{K}_{L} & \delta^{K}_{M} & \delta^{K}_{N}
       *  \end{Vmatrix}
       * = \delta^{I}_{N}\delta^{J}_{M}\delta^{K}_{L}
       *  + \delta^{I}_{M}\delta^{J}_{L}\delta^{K}_{N}
       *  + \delta^{I}_{L}\delta^{J}_{N}\delta^{K}_{M}
       *  - \delta^{I}_{L}\delta^{J}_{M}\delta^{K}_{N}
       *  - \delta^{I}_{M}\delta^{J}_{N}\delta^{K}_{L}
       *  - \delta^{I}_{N}\delta^{J}_{L}\delta^{K}_{M}
       * \f]
       *
       *  @code
       *
       *  const int alpha_XYZXYE = Epsilon1_<X,Y,Z,X,Y,E>::value ;
       *
       *  @endcode
       *
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2008-07-28
       */
      template <unsigned int I,
                unsigned int J,
                unsigned int K,
                unsigned int L,
                unsigned int M,
                unsigned int N>
      struct Epsilon1_
      {
        enum
          {
            value = Delta_<I,N>::value * Delta_<J,M>::value * Delta_<K,L>::value
            +       Delta_<I,M>::value * Delta_<J,L>::value * Delta_<K,N>::value
            +       Delta_<I,L>::value * Delta_<J,N>::value * Delta_<K,M>::value
            //
            -       Delta_<I,L>::value * Delta_<J,M>::value * Delta_<K,N>::value
            -       Delta_<I,M>::value * Delta_<J,N>::value * Delta_<K,L>::value
            -       Delta_<I,N>::value * Delta_<J,L>::value * Delta_<K,M>::value
          } ;
      } ;
      // ======================================================================
      /** @struct Epsilon2_
       *  (Compile-time) evaluation of the tensor product of
       *  two Levi-Civita symbols:
       *
       * \f[
       *  \alpha^{IJ}_{KL} =
       *  \epsilon^{IJ\gamma\kappa}
       *  \epsilon_{KL\rho\kappa} = -
       *  \begin{Vmatrix}
       *    \delta^{I}_{K} & \delta^{I}_{L} \   \
       *    \delta^{J}_{K} & \delta^{J}_{L}  \  \
       *  \end{Vmatrix}
       * \f]
       *
       *  @code
       *
       *  const int alpha_XYYE = Epsilon2_<X,Y,Y,E>::value ;
       *
       *  @endcode
       *
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2008-07-28
       */
      template <unsigned int I,
                unsigned int J,
                unsigned int K,
                unsigned int L>
      struct Epsilon2_
      {
        enum
          { value = -2 * ( Delta_<I,K>::value * Delta_<J,L>::value -
                           Delta_<J,K>::value * Delta_<I,L>::value ) } ;
      } ;
      // ======================================================================
      namespace detail 
      {
        // ====================================================================
        template <class TYPE>
        struct _TYPES
        {
          typedef long double iTYPE ;
          typedef      double rTYPE ;
        };
        template <class TYPE> 
        struct _TYPES< std::complex<TYPE> >
        {
          typedef std::complex<double> iTYPE ;
          typedef std::complex<double> rTYPE ;
        } ;
        template <class COORDINATES>
        struct LVTYPES 
        {
          typedef typename _TYPES<typename COORDINATES::Scalar>::iTYPE iTYPE ;
          typedef typename _TYPES<typename COORDINATES::Scalar>::rTYPE rTYPE ;          
          typedef ROOT::Math::LorentzVector<COORDINATES>               vTYPE ;
        } ;  
        // ====================================================================
      }
      // ======================================================================
      /** @struct Epsilon
       *  Simple implementation of 4D Antisymmetric Levi-Civita symbols
       *   \f$ \epsilon_{\mu\nu\lambda\delta} \f$ and some related operations
       *
       *  Convention here: \f$ \epsilon_{0123} = \epsilon_{XYZT} = 1 \f$, that gives  
       *  - \f$ x == epsilon ( t , y , z ) \f$ 
       *  - \f$ y == epsilon ( t , z , x ) \f$ 
       *  - \f$ z == epsilon ( t , x , y ) \f$ 
       *  - \f$ t == epsilon ( x , y , z ) \f$ 
       * 
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2008-0725
       */
      struct Epsilon
      {
      public:
        // ====================================================================
        /** the major method for evaluation of elements of Levi-Civita tensor
         *   \f$ \epsilon_{\mu\nu\lambda\delta} \f$
         *
         *  Convention here: \f$ \epsilon_{0123} = \epsilon_{XYZT} = 1 \f$
         *
         *  @code
         *
         *   Epsilon e ;
         *
         *   int symbol = e(0,1,2,3) ;
         *
         *  @endcode
         *
         *  @attention the evaluation could be rather CPU expensive,
         *             please use the templated form if indices
         *             are known at compile-time
         *
         *  @param i the first  index
         *  @param j the second index
         *  @param k the third  index
         *  @param l the last   index
         *  @return the value of Levi-Civita symbol
         */
        inline int operator()
        ( const unsigned short i ,
          const unsigned short j ,
          const unsigned short k ,
          const unsigned short l ) const { return epsilon ( i , j , k , l ) ; }
        // ====================================================================
        /** evaluate the tensor product: e*v
         *
         *  \f$  t_{\mu\nu\lambda} =
         *  \epsilon_{\mu\nu\lambda\kappa}v^{\kappa} \f$
         *
         *  @code
         *
         *   const LorentzVector& v = ... ;
         *
         *   Epsilon e ;
         *
         *   const double t_XYE = e ( X , Y , E , v ) ;
         *
         *  @endcode
         *
         *  Convention here: \f$ \epsilon_{0123} = \epsilon_{xyzt} = 1 \f$
         *
         *  @attention the evaluation could be rather CPU expensive,
         *             please use the templated form if indices
         *             are known at compile-time
         *
         *  @see Epsilon::e_3
         *
         *  @param mu     the first  index
         *  @param nu     the second index
         *  @param lambda the third  index
         *  @param v the input vector
         *  @return the product (component)
         */
        inline double operator ()
        ( const unsigned short        mu     ,
          const unsigned short        nu     ,
          const unsigned short        lambda ,
          const Ostap::LorentzVector& v      ) const 
        { return epsilon ( mu , nu , lambda , v ) ; }
        // ====================================================================
        /** evaluate the tensor product: e*v1*v2
         *
         *  \f$  v_{\mu\nu} =
         *  \epsilon_{\mu\nu\lambda\kappa}v_1^{\lambda}v_2^{\kappa} \f$
         *
         *  @code
         *
         *   const LorentzVector& v1 = ... ;
         *   const LorentzVector& v2 = ... ;
         *
         *   Epsilon e ;
         *
         *   const double v_XY = e ( X , Y , v1 , v2 ) ;
         *
         *  @endcode
         *
         *  Convention here: \f$ \epsilon_{0123} = \epsilon_{XYZT} = 1 \f$
         *
         *  @attention the evaluation could be rather CPU expensive,
         *             please use the templated form if indices
         *             are known at compile-time
         *
         *  @see Epsilon::e_2
         *
         *  @param mu the first index
         *  @param nu the second index
         *  @param v1 the first  vector
         *  @param v2 the second vector
         *  @return the product (component)
         */
        inline double operator()
        ( const unsigned short        mu ,
          const unsigned short        nu ,
          const Ostap::LorentzVector& v1 ,
          const Ostap::LorentzVector& v2 ) const 
        { return epsilon ( mu , nu , v1 , v2 ) ; }
        // ======================================================================
        /** evaluate the e*v1*v2*v3 product
         *
         *  \f$  v_{\mu} =
         *  \epsilon_{\mu\nu\lambda\kappa}v_1^{\nu}v_2^{\lambda}v_3^{\kappa} \f$
         *
         *  @code
         *
         *   const LorentzVector& v1 = ... ;
         *   const LorentzVector& v2 = ... ;
         *   const LorentzVector& v3 = ... ;
         *
         *   Epsilon e ;
         *
         *   const double v_X = e ( X , v1 , v2 , v3 ) ;
         *
         *  @endcode
         *
         *  Convention here: \f$ \epsilon_{0123} = \epsilon_{XYZT} = 1 \f$
         *
         *  @param mu the first index
         *  @param v1 the first  vector
         *  @param v2 the second vector
         *  @param v3 the third  vector
         *  @return the product (component)
         */
        inline double operator()
        ( const unsigned short        mu ,
          const Ostap::LorentzVector& v1 ,
          const Ostap::LorentzVector& v2 ,
          const Ostap::LorentzVector& v3 ) const 
        { return epsilon ( mu , v1 , v2 , v3 ) ; }
        // ======================================================================
        /** evaluate the e*v1*v2*v3 product
         *
         *  \f$  v_{\mu} =
         *  \epsilon_{\mu\nu\lambda\kappa}v_1^{\nu}v_2^{\lambda}v_3^{\kappa} \f$
         *
         *  @code
         *
         *   const LorentzVector& v1 = ... ;
         *   const LorentzVector& v2 = ... ;
         *   const LorentzVector& v3 = ... ;
         *
         *   Epsilon e ;
         *
         *   const LorentzVector = e ( v1 , v2 , v3 ) ;
         *
         *  @endcode
         *
         *  Convention here: \f$ \epsilon_{0123} = \epsilon_{xyzt} = 1 \f$
         *
         *  The following identity holds numerically:
         *  @code
         *
         *   const LorentzVector& v1  = ... ;
         *   const LorentzVector& v2  = ... ;
         *   const LorentzVector& v3  = ... ;
         *
         *   Epsilon e ;
         *
         *   // NUMERICAL INDENTITY:
         *
         *   Assert( e ( v1 , v2 , v3 , v4 ) == v1.Dot( e ( v2 , v3 , v4 ) ) , ... ) ;
         *
         *  @endcode
         *
         *  @see Epsilon::e_1
         *
         *  @param v1 the first  vector
         *  @param v2 the second vector
         *  @param v3 the third  vector
         *  @return the product (vector)
         */
        inline LorentzVector operator()
        ( const Ostap::LorentzVector& v1 ,
          const Ostap::LorentzVector& v2 ,
          const Ostap::LorentzVector& v3 ) const
        { return epsilon ( v1 , v2 , v3 ) ; }
        // ======================================================================
        /** evaluate the e*v1*v2*v3*v4 product
         *
         *  \f$  r =
         *  \epsilon_{\mu\nu\lambda\kappa}
         *       v_1^{\mu}
         *       v_2^{\nu}
         *       v_3^{\lambda}
         *       v_4^{\kappa} \f$
         *
         *  @code
         *
         *   const LorentzVector& v1 = ... ;
         *   const LorentzVector& v2 = ... ;
         *   const LorentzVector& v3 = ... ;
         *   const LorentzVector& v4 = ... ;
         *
         *   Epsilon e ;
         *
         *   const double v_X = e ( v1 , v2 , v3 , v4 ) ;
         *
         *  @endcode
         *
         *  Convention here: \f$ \epsilon_{0123} = \epsilon_{XYZT} = 1 \f$
         *
         *  The following identity holds numerically:
         *  @code
         *
         *   const LorentzVector& v1  = ... ;
         *   const LorentzVector& v2  = ... ;
         *   const LorentzVector& v3  = ... ;
         *   const LorentzVector& v4  = ... ;
         *
         *   Epsilon e ;
         *
         *   // NUMERICAL INDENTITY:
         *
         *   Assert( e ( v1 , v2 , v3 , v4 ) == v1.Dot( e ( v2 , v3 , v4 ) ) , ... ) ;
         *
         *  @endcode
         *
         *  @param v1 the first  vector
         *  @param v2 the second vector
         *  @param v3 the third  vector
         *  @param v3 the fourth vector
         *  @return the product (vector)
         */
        inline double operator()
        ( const Ostap::LorentzVector& v1 ,
          const Ostap::LorentzVector& v2 ,
          const Ostap::LorentzVector& v3 ,
          const Ostap::LorentzVector& v4 ) const
        { return epsilon ( v1 , v2 , v3 , v4 ) ; }
        // ======================================================================
        /* evaluate the tensor product: (e*v1*v2*v3)*(e*v4*v5*v6)
         *
         *  \f$  r =
         *  \epsilon_{\mu\nu\lambda\kappa}
         *  \epsilon_{\mu\rho\delta\tau}
         *      v_1^{\nu}v_2^{\lambda}v_3^{\kappa}
         *      v_4^{\rho}v_5^{\delta}v_6^{\tau}   \f$
         *
         *  This expression typically appears in evaution of
         *  various "plane-angles".
         *
         *  @code
         *
         *   const LorentzVector& v1 = ... ;
         *   const LorentzVector& v2 = ... ;
         *   const LorentzVector& v3 = ... ;
         *   const LorentzVector& v4 = ... ;
         *   const LorentzVector& v5 = ... ;
         *   const LorentzVector& v6 = ... ;
         *
         *   Epsilon e ;
         *
         *   const double r = e ( v1 , v2 , v3 , v4 , v5 , v6 ) ;
         *
         *  @endcode
         *
         *  Convention here: \f$ \epsilon_{0123} = \epsilon_{XYZT} = 1 \f$
         *
         *  @param v1 the first  vector
         *  @param v2 the second vector
         *  @param v3 the third  vector
         *  @param v4 the fourth vector
         *  @param v5 the fith   vector
         *  @param v6 the sixth  vector
         *  @return the tensor product
         */
        double operator()
        ( const Ostap::LorentzVector& v1 ,
          const Ostap::LorentzVector& v2 ,
          const Ostap::LorentzVector& v3 ,
          const Ostap::LorentzVector& v4 ,
          const Ostap::LorentzVector& v5 ,
          const Ostap::LorentzVector& v6 ) const
        { return epsilon ( v1 , v2 , v3 , v4 , v5 , v6 ) ; }
        // =====================================================================
      public:
        // ======================================================================
        /** the major method for evaluation of elements of Levi-Civita tensor
         *  @code
         *
         *   Epsilon lcs ;
         *
         *   int symbol = lcs.epsilon (0,1,2,3) ;
         *
         *  @endcode
         *
         *  @attention The evaluation could  be rather time consuming,
         *             If the indices are known at compile time, try to use
         *             the templated methods - they are much efficient
         *
         *  Convention here: \f$ \epsilon_{0123} = \epsilon_{XYZT} = 1 \f$
         *
         *  @param i the first  index
         *  @param j the second index
         *  @param k the third  index
         *  @param l the last   index
         *  @return the value of Levi-Civita symbol
         */
        static constexpr inline int epsilon 
        ( const unsigned short i , 
          const unsigned short j ,
          const unsigned short k , 
          const unsigned short l ) 
        {
          /// the regular cases
          if (  i <  j && j <  k && k <  l && l <  4 ) { return 1 ; } // RETURN
          if (  i == j || j == k || k == l || l == i ) { return 0 ; } // RETURN
          if (  i >  3 || j >  3 || k >  3 || l >  3 ) { return 0 ; } // RETURN
          /// permutations are required:
          if ( i > j  ) { return -epsilon ( j , i , k , l ) ; } // RETURN
          if ( j > k  ) { return -epsilon ( i , k , j , l ) ; } // RETURN
          if ( k > l  ) { return -epsilon ( i , j , l , k ) ; } // RETURN
          /// here we can go only if some of number >=4, return 0..
          return 0 ;
        }
        // =====================================================================
      public: // various tensor products 
        // ======================================================================
        /** evaluate the tensor product e*v1
         *
         *  \f$  t_{\mu\nu\lambda} =
         *  \epsilon_{\mu\nu\lambda\kappa}v_1^{\kappa} \f$
         *
         *  @code
         *  const LorentzVector& v = ... ;
         *  const double t_XYE = Epsilon::epsilon(X,Y,E, v )
         *  @endcode
         *
         *  Convention here: \f$ \epsilon_{0123} = \epsilon_{XYZT} = 1 \f$
         *
         *  @see Epsilon::e_3
         *
         *  @param v the input vector
         *  @return the product (component)
         */
        template <class  COORDINATES>
        static inline typename detail::LVTYPES<COORDINATES>::rTYPE
        epsilon
        ( const unsigned short                          mu     ,
          const unsigned short                          nu     ,
          const unsigned short                          lambda ,
          const ROOT::Math::LorentzVector<COORDINATES>& v1     ) ;
        // ======================================================================
        /** evaluate the tensor e*v1*v2 product
         *
         *  \f$  v_{\mu\nu} =
         *  \epsilon_{\mu\nu\lambda\kappa}
         *       v_1^{\lambda}v_2^{\kappa} \f$
         *
         *  @code
         *  const LorentzVector& v1 = ... ;
         *  const LorentzVector& v2 = ... ;
         *  const double v_XY = Epsion::epsilon(X,Y,v1 , v2 )
         *  @endcode
         *
         *  Convention here: \f$ \epsilon_{0123} = \epsilon_{XYZT} = 1 \f$
         *
         *  @see Epsilon::e_2
         *
         *  @param mu the first index
         *  @param nu the second index
         *  @param v1 the first  vector
         *  @param v2 the second vector
         *  @return the product (component)
         */
        template <class  COORDINATES>
        static inline typename detail::LVTYPES<COORDINATES>::rTYPE
        epsilon
        ( const unsigned short                          mu ,
          const unsigned short                          nu ,
          const ROOT::Math::LorentzVector<COORDINATES>& v1 ,
          const ROOT::Math::LorentzVector<COORDINATES>& v2 ) ;
        // ======================================================================
        /** evaluate the e*v1*v2*v3 product
         *
         *  \f$  v_{\mu} =
         *  \epsilon_{\mu\nu\lambda\kappa}v_1^{\nu}v_2^{\lambda}v_3^{\kappa} \f$
         *
         *  @code
         *
         *   const LorentzVector& v1 = ... ;
         *   const LorentzVector& v2 = ... ;
         *   const LorentzVector& v3 = ... ;
         *
         *   const double v_X = Epsilon::epsilon ( Epsilon::X , v1 , v2 , v3 )
         *
         *  @endcode
         *
         *  Convention here: \f$ \epsilon_{0123} = \epsilon_{XYZT} = 1 \f$
         *
         *  @see Epsilon::e_1
         *
         *  @param mu the index of the result
         *  @param v1 the first  vector
         *  @param v2 the second vector
         *  @param v3 the third  vector
         *  @return the product (component)
         */
        template <class  COORDINATES>
        static inline typename detail::LVTYPES<COORDINATES>::rTYPE
        epsilon
        ( const unsigned short                          mu ,
          const ROOT::Math::LorentzVector<COORDINATES>& v1 ,
          const ROOT::Math::LorentzVector<COORDINATES>& v2 ,
          const ROOT::Math::LorentzVector<COORDINATES>& v3 ) ;
        // ======================================================================
        /** evaluate the e*v1*v2*v3 product
         *
         *  \f$  v_{\mu} =
         *  \epsilon_{\mu\nu\lambda\kappa}
         *   v_1^{\nu}v_2^{\lambda}v_3^{\kappa} \f$
         *
         *  @code
         *
         *   const LorentzVector& v1 = ... ;
         *   const LorentzVector& v2 = ... ;
         *   const LorentzVector& v3 = ... ;
         *
         *   Epsilon e ;
         *
         *   const LorentzVector v = Epsilon::epsilon ( v1 , v2 , v3 )
         *
         *  @endcode
         *
         *  Convention here: \f$ \epsilon_{0123} = \epsilon_{XYZT} = 1 \f$
         *
         *  @see Epsilon::e_1
         *
         *  @param v1 the first  vector
         *  @param v2 the second vector
         *  @param v3 the third  vector
         *  @return the product vector
         */
        template <class  COORDINATES>
        static inline typename detail::LVTYPES<COORDINATES>::vTYPE
        epsilon
        ( const ROOT::Math::LorentzVector<COORDINATES>& v1 ,
          const ROOT::Math::LorentzVector<COORDINATES>& v2 ,
          const ROOT::Math::LorentzVector<COORDINATES>& v3 ) ;
        // ======================================================================
        /** evaluate the tensor product: e*v1*v2*v3*v4
         *
         *  \f$  r =
         *  \epsilon_{\mu\nu\lambda\kappa}
         *    v_1^{\mu}
         *    v_2^{\nu}
         *    v_3^{\lambda}
         *    v_4^{\kappa} \f$
         *
         *  @code
         *
         *   const LorentzVector& v1 = ... ;
         *   const LorentzVector& v2 = ... ;
         *   const LorentzVector& v3 = ... ;
         *   const LorentzVector& v4 = ... ;
         *
         *   Epsilon e ;
         *
         *   const double r = e.epsilon ( v1 , v2 , v3 , v4 ) ;
         *
         *  @endcode
         *
         *  Convention here: \f$ \epsilon_{0123} = \epsilon_{XYZT} = 1 \f$
         *
         *  @param v1 the first  vector
         *  @param v2 the second vector
         *  @param v3 the third  vector
         *  @param v4 the fourth vector
         *  @return the tensor product
         */
        template <class  COORDINATES>
        static inline typename detail::LVTYPES<COORDINATES>::rTYPE
        epsilon
        ( const ROOT::Math::LorentzVector<COORDINATES>& v1 ,
          const ROOT::Math::LorentzVector<COORDINATES>& v2 ,
          const ROOT::Math::LorentzVector<COORDINATES>& v3 ,
          const ROOT::Math::LorentzVector<COORDINATES>& v4 ) ;
        // ======================================================================
        /** evaluate the tensor product: (e*a1*a2*a3)*(e*b1*b2*b3)
         *
         *  \f$  r =
         *  \epsilon_{\mu\nu\lambda\kappa}
         *  \epsilon_{\mu\rho\delta\tau}
         *          a_1^{\nu}a_2^{\lambda}a_3^{\kappa}
         *          b_1^{\rho}b_2^{\delta}b_3^{\tau}   \f$
         *
         *  This expression typically appears in evalution of
         *  various "plane-angles" (for this case a3=b3)
         *
         *  @code
         *
         *   const LorentzVector& a1 = ... ;
         *   const LorentzVector& a2 = ... ;
         *   const LorentzVector& a3 = ... ;
         *   const LorentzVector& b1 = ... ;
         *   const LorentzVector& b2 = ... ;
         *   const LorentzVector& b3 = ... ;
         *
         *   const double r = Epsilon::epsilon ( a1 , a2 , a3 , b1 , b2 , b3 ) ;
         *
         *  @endcode
         *
         *  The following identity has been used:
         * \f[
         *     \epsilon^{\alpha\beta\gamma\kappa}
         *     \epsilon_{\mu\mu\rho\kappa} = -
         *  \begin{Vmatrix}
         *   \delta^{\alpha}_{\mu} & \delta^{\alpha}_{\nu} & \delta^{\alpha}_{\rho} \ \
         *   \delta^{\beta}_{\mu}  & \delta^{\beta}_{\nu}  & \delta^{\beta}_{\rho}  \ \
         *   \delta^{\gamma}_{\mu} & \delta^{\gamma}_{\nu} & \delta^{\gamma}_{\rho}
         *  \end{Vmatrix}
         * \f]
         *
         *  Convention here: \f$ \epsilon_{0123} = \epsilon_{XYZT} = 1 \f$
         *
         *  @param a1 the first  vector
         *  @param a2 the second vector
         *  @param a3 the third  vector
         *  @param b1 the fourth vector
         *  @param b2 the fith   vector
         *  @param b3 the sixth  vector
         *  @return the tensor product
         */
        template <class  COORDINATES>
        static inline typename detail::LVTYPES<COORDINATES>::rTYPE
        epsilon
        ( const ROOT::Math::LorentzVector<COORDINATES>& a1 ,
          const ROOT::Math::LorentzVector<COORDINATES>& a2 ,
          const ROOT::Math::LorentzVector<COORDINATES>& a3 ,
          const ROOT::Math::LorentzVector<COORDINATES>& b1 ,
          const ROOT::Math::LorentzVector<COORDINATES>& b2 ,
          const ROOT::Math::LorentzVector<COORDINATES>& b3 ) ;        
        // ======================================================================
      public:
        // ======================================================================
        /** evaluate the magnitude of the "4-normal"
         *
         *  \f$ l^2 = L_{\mu}L^{\mu}
         *      = \left( \epsilon_{\mu\nu\lambda\delta}
         *     a^{\nu}b^{\lambda}c^{\delta} \right)^2
         *    =
         *     \left(ab\right)c^2 +
         *     \left(ac\right)b^2 +
         *     \left(bc\right)a^2 -
         *      a^2b^2c^2 - 2\left(ab\right)\left(bc\right)\left(ac\right)
         *   \f$
         *
         *  @attention For time-like input vectors,
         *             the 4-normal is a space-like vector,
         *             and therefore the result must be non-positive!
         *
         *  @param a the first  vector
         *  @param b the second vector
         *  @param c the third  vector
         *  @return the magnitude of 4-normal
         */
        template <class  COORDINATES>
        static inline typename detail::LVTYPES<COORDINATES>::rTYPE
        mag2
        ( const ROOT::Math::LorentzVector<COORDINATES>& a ,
          const ROOT::Math::LorentzVector<COORDINATES>& b ,
          const ROOT::Math::LorentzVector<COORDINATES>& c ) ;
        // ======================================================================
      public:
        // ======================================================================
        /** evaluate the tensor product e*v1
         *
         *  \f$  t_{\mu\nu\lambda} =
         *  \epsilon_{\mu\nu\lambda\kappa}v_1^{\kappa} \f$
         *
         *  @code
         *
         *   const LorentzVector& v = ... ;
         *   const double t_XYE = Epsilon::e_3<Epsilon::X,Epsilon::Y,Epsilon::E>( v )
         *
         *  @endcode
         *
         *  @param v the input vector
         *  @return the product (component)
         */
        template <unsigned int I, unsigned int J, unsigned int K, class  COORDINATES>
        static inline typename detail::LVTYPES<COORDINATES>::rTYPE
        e_3 ( const ROOT::Math::LorentzVector<COORDINATES>& v ) ;
        // ======================================================================
        /** evaluate the tensor e*v1*v2 product
         *
         *  \f$  v_{\mu\nu} =
         *  \epsilon_{\mu\nu\lambda\kappa}v_1^{\lambda}v_2^{\kappa} \f$
         *
         *  @code
         *
         *   const LorentzVector& v1 = ... ;
         *   const LorentzVector& v2 = ... ;
         *   const double v_XY = Epsilon::e_2<Epsilon::X,Epsilon::Y>( v1 , v2 )
         *
         *  @endcode
         *
         *  @param v1 the first  vector
         *  @param v2 the second vector
         *  @return the product (component)
         */
        template <unsigned int I, unsigned int J, class COORDINATES>
        static inline typename detail::LVTYPES<COORDINATES>::rTYPE 
        e_2 ( const ROOT::Math::LorentzVector<COORDINATES>& v1 ,
              const ROOT::Math::LorentzVector<COORDINATES>& v2 ) ;
        // ======================================================================
        /** evaluate the e*v1*v2*v3 product
         *
         *  \f$  v_{\mu} =
         *  \epsilon_{\mu\nu\lambda\kappa}v_1^{\nu}v_2^{\lambda}v_3^{\kappa} \f$
         *
         *  @code
         *
         *   const LorentzVector& v1 = ... ;
         *   const LorentzVector& v2 = ... ;
         *   const LorentzVector& v3 = ... ;
         *   const double v_X = Epsilon::e_1<Epsilon::X>( v1 , v2 , v3 )
         *
         *  @endcode
         *
         *  @param v1 the first  vector
         *  @param v2 the second vector
         *  @param v3 the third  vector
         *  @return the product (component)
         */
        template <unsigned int I, class COORDINATES>
        static inline typename detail::LVTYPES<COORDINATES>::rTYPE 
        e_1 ( const ROOT::Math::LorentzVector<COORDINATES>& v1 ,
              const ROOT::Math::LorentzVector<COORDINATES>& v2 ,
              const ROOT::Math::LorentzVector<COORDINATES>& v3 ) ;
        // ====================================================================
      };
      // ======================================================================
    } //                                  end of namespace Ostap::Math::Tensors
    // ========================================================================
  } //                                             end of namespace Ostap::Math
    // ========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Tensors.icpp"
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_TENSORS_H
// ============================================================================
