// ============================================================================
#ifndef OSTAP_DALITZ_H
#define OSTAP_DALITZ_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <cmath>
#include <array>
#include <utility>
#include <functional>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Kinematics.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Kinematics 
  {
    // ========================================================================
    /** @class Dalitz0 Ostap/Dalitz.h 
     *  Simple Kinematics of Dalitz plot 
     *  @see E.Byckling, K.Kajantie, "Particle kinematics", John Wiley & Sons,
     *              London, New York, Sydney, Toronto, 1973, p.89, eq. (5.23)     
     *  @see Section V.1
     *  @see https://userweb.jlab.org/~rafopar/Book/byckling_kajantie.pdf     
     */
    class Dalitz0
    {
    public :
      // ======================================================================
      /** constructor from three masses 
       *  - m1 : the mass of the first particle  \f$ m_1 \f$;
       *  - m2 : the mass of the second particle  \f$ m_2 \f$;
       *  - m3 : the mass of the third particle  \f$ m_3 \f$;
       */
      Dalitz0
      ( const double m1 = 0 , 
        const double m2 = 0 , 
        const double m3 = 0 ) ;
      // ======================================================================
    public:  // trivial getters 
      // ====================================================================== 
      ///   the first mass 
      double m1     () const  { return m_m1  ; }
      ///   the second mass 
      double m2     () const  { return m_m2  ; }
      ///   the third  mass 
      double m3     () const  { return m_m3  ; }
      // ======================================================================
      /// minimal value of s1: \f$ \left. s_1\right|_{min} =  (m_1+m_2)^2 \f$
      double s1_min () const { return m_cache [ 0 ] ; }
      /// minimal value of s2: \f$ \left. s_2\right|_{min} =  (m_2+m_3)^2 \f$
      double s2_min () const { return m_cache [ 1 ] ; }
      /// minimal value of s3: \f$ \left. s_3\right|_{min} =  (m_3+m_1)^2 \f$
      double s3_min () const { return m_cache [ 2 ] ; }
      // ======================================================================
      /// maximal value of s1: \f$ \left. s_1\right|_{max} =  (M-m_3)^2 \f$
      double s1_max ( const double M ) const { const double d = M - m3 () ; return d * d ; }
      /// maximal value of s2: \f$ \left. s_1\right|_{max} =  (M-m_1)^2 \f$
      double s2_max ( const double M ) const { const double d = M - m1 () ; return d * d ; }
      /// maximal value of s3: \f$ \left. s_3\right|_{max} =  (M-m_2)^2 \f$
      double s3_max ( const double M ) const { const double d = M - m2 () ; return d * d ; }
      // ======================================================================
      /// \f$ m_1^2\f$  precalculated 
      double m1sq   () const { return m_cache [ 3 ] ; }
      /// \f$ m_2^2\f$  precalculated 
      double m2sq   () const { return m_cache [ 4 ] ; }
      /// \f$ m_3^2\f$ precalculated 
      double m3sq   () const { return m_cache [ 5 ] ; }
      /// sum of squared masses  \f$  m_1^2+m_2^1+m_3^3\f$ 
      double summ2  () const { return m_cache [ 6 ] ; }
      /// sum of masses  \f$  m_1 + m_2 + m_3 \f$ 
      double summ   () const { return m_cache [ 7 ] ; }
      /// squared sum of masses  \f$  (m_1 + m_2 + m_3)^2 \f$ 
      double sqsumm () const { return m_cache [ 8 ] ; }
      ///  minimal value of s: \f$ s_{min} = \left(  m_1 +m_2+ m_3 \right)^2\f$
      double s_min  () const { return sqsumm () ; }      
      // ======================================================================
    public:
      // ======================================================================
      /// Is m1 equal to zero ? 
      bool m1_zero () const { return m_cacheb [ 0 ] ; }
      /// Is m1 equal to zero ? 
      bool m2_zero () const { return m_cacheb [ 1 ] ; }
      /// Is m1 equal to zero ? 
      bool m3_zero () const { return m_cacheb [ 2 ] ; }
      // ======================================================================
    public:  // limits for E1 , E2 and E3
      // ======================================================================
      double E1_min () const { return m_m1          ; }
      double E2_min () const { return m_m2          ; }
      double E3_min () const { return m_m3          ; }
      // ======================================================================
    public:
      // ======================================================================
      /// maximal value of momenta of the first  particle for the given s 
      double p1_max ( const double s ) const ;
      /// maximal value of momenta of the second particle for the given s 
      double p2_max ( const double s ) const ;
      /// maximal value of momenta of the third  particle for the given s 
      double p3_max ( const double s ) const ;
      // ======================================================================
    public:  // only two s_i are independent 
      // ======================================================================
      /** get s_1 for given s_2 and s_3 are defined 
       *  \f$ s_1 = s_{12} = \sum s -   s_2 - s_3  \f$ 
       */
      inline double s1 ( const double s , const double s2 , const double s3 ) const 
      { return s  + summ2 () - s2 - s3 ; }
      // ======================================================================
      /** get s_2 for given s_1 and s_3 are defined 
       *  \f$ s_2 = s_{23} = \sum s -   s_1 - s_3  \f$ 
       */
      inline double s2 ( const double s , const double s1 , const double s3 ) const 
      { return s + summ2 () - s1 - s3 ; }
      // ======================================================================
      /** get s_3 for given s_1 and s_2 are defined 
       *  \f$ s_3 = s_{31} = \sum s -   s_1 - s_2  \f$ 
       */
      inline double s3 ( const double s , const double s1 , const double s2 ) const 
      { return s + summ2  () - s1 - s2 ; }
      // ======================================================================
    public :  // geometry of Dalitz plot 
      // ======================================================================
      /** Is the point \f$ (s, s_1,s_2)\f$ "inside" the Dalitz plot?
       *  Get the sign of G-function 
       *  \f$ g(s_1,s_2) = G ( s_1, s_2 , s , m_2^2, m_1^2, m_3^2) \f$
       *  @see Ostap::Kinematics::G
       *  Physical region corresponds to \f$ g\le0 \f$  
       */
      bool   inside   ( const double s  ,
                        const double s1 ,
                        const double s2 ) const ;
      // =======================================================================
      /** get the measure of the distance form the point to the boundary of Dalitz plot. 
       *  Distance is defined as 
       *  \f$ d \equiv = \lambda ( P_1^2 , P_2^2 , P_3^2 ) \f$ 
       */
      double distance ( const double s ,
                        const double s1 ,
                        const double s2 ) const ;
      // ======================================================================
    public:  // more invariants  
      // ======================================================================
      /** get the product of 4-momenta
       *  \f$ \left(p_1p_2\right) = \frac{1}{2}\left( s_{12} - m_1^2 - m_2^2 \right) \f$
       */
      inline double p1p2
      ( const double /* s  */ ,
        const double    s1    ,
        const double /* s2 */ )  const 
      { return 0.5 * ( s1 - m1sq () - m2sq () ) ; }
      // ======================================================================
      /** get the product of 4-momenta
       *  \f$ \left(p_2p_3\right) = \frac{1}{2}\left( s_{23} - m_2^2 - m_3^2 \right) \f$      
       */
      inline double p2p3
      ( const double /* s  */ ,
        const double /* s1 */ ,
        const double    s2    ) const
      { return 0.5 * ( s2 - m2sq () -  m3sq () ) ; }      
      // ======================================================================
      /** get the product of 4-momenta
       *  \f$ \left(p_1p_3\right) = \frac{1}{2}\left( s_{13} - m_1^2 - m_3^2 \right) \f$      
       */
      inline double p1p3
      ( const double    s     ,
        const double    s1    ,
        const double    s2    ) const
      { return 0.5 * ( s3  ( s , s1 , s2 ) - m1sq () -  m3sq () ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /** get the product of 4-momenta
       *  \f$ \left(pp1\right) = \frac{1}{2}\left( s - s_{23} +  m_1^2 \right) \f$
       */
      inline double pp1
      ( const double    s     ,
        const double /* s1 */ ,
        const double    s2    ) const
      { return 0.5 * ( s - s2 + m1sq () ) ; }
      // ======================================================================
      /** get the product of 4-momenta
       *  \f$ \left(pp2\right) = \frac{1}{2}\left( s - s_{13} + m_2^2 \right) \f$
       */
      inline double pp2
      ( const double s  ,
        const double s1 ,
        const double s2 ) const
      { return 0.5 * ( s - s3  ( s , s1 , s2 )  + m2sq () ) ; }
      // ======================================================================
      /** get the product of 4-momenta
       *  \f$ \left(pp3\right) = \frac{1}{2}\left( s - s_{12} + m_3^2 \right) \f$
       */
      inline double pp3 
      ( const double    s     ,
        const double    s1    ,
        const double /* s2 */ ) const
      { return 0.5 * ( s - s1 + m3sq () ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /** get product of 4-momenta 
       *  \f$ pp_{12} = \frac{1}{2}\left(s + s_{12} - m^2_3\right)\f$
       */
      inline double pp12 
      ( const double    s     ,
        const double    s1    ,
        const double /* s2 */ ) const
      { return 0.5 * ( s + s1 - m3sq () ) ; }
      // ======================================================================
      /** get product of 4-momenta 
       *  \f$ pp_{23} = \frac{1}{2}\left(s + s_{23} - m^2_1\right)\f$
       */
      inline double pp23
      ( const double    s     ,
        const double /* s1 */ ,
        const double    s2    ) const
      { return 0.5 * ( s + s2 - m1sq () ) ; }
      // ======================================================================
      /** get product of 4-momenta 
       *  \f$ pp_{13} = \frac{1}{2}\left(s + s_{13} - m^2_2\right)\f$
       */
      inline double pp13
      ( const double    s     ,
        const double    s1    ,
        const double    s2    ) const
      { return 0.5 * ( s + s3 ( s , s1 , s2 )  - m2sq () ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /** get product of 4-momenta \f$ p_1p_{12} = m_1^2 +  p_1p_2\f$ 
       */
      inline double p1p12 
      ( const double    s     ,
        const double    s1    ,
        const double    s2    ) const
      { return m1sq () + p1p2 ( s , s1 , s2 ) ; }
      // =====================================================================
      /** get product of 4-momenta \f$ p_1p_{13} = m_1^2 +  p_1p_3\f$ 
       */
      inline double p1p13 
      ( const double    s     ,
        const double    s1    ,
        const double    s2    ) const
      { return m1sq () + p1p3 ( s , s1 , s2 ) ; }
      // ======================================================================
      /** get product of 4-momenta \f$ p_1p_{23} = p_1p_2 +  p_1p_3\f$ 
       */
      inline double p1p23 
      ( const double    s     ,
        const double    s1    ,
        const double    s2    ) const
      { return p1p2 ( s , s1 , s2 ) + p1p3 ( s , s1 , s2 ) ; }
      // ======================================================================
      /** get product of 4-momenta \f$ p_2p_{12} = p_1p_2 + m^2_2\f$ 
       */
      inline double p2p12 
      ( const double    s     ,
        const double    s1    ,
        const double    s2    ) const
      { return p1p2 ( s , s1 , s2 ) + m2sq () ; }
      // =====================================================================
      /** get product of 4-momenta \f$ p_2p_{13} = p_1p_2 + p2p3 \f$ 
       */
      inline double p2p13 
      ( const double    s     ,
        const double    s1    ,
        const double    s2    ) const
      { return p1p2 ( s , s1 , s2 ) + p2p3 ( s , s1 , s2 )  ; }
      // =====================================================================
      /** get product of 4-momenta \f$ p_2p_{23} = m^2_2 p_2p_3 \f$ 
       */
      inline double p2p23 
      ( const double    s     ,
        const double    s1    ,
        const double    s2    ) const
      { return m2sq () + p2p3 ( s , s1 , s2 )  ; }
      // =====================================================================      
      /** get product of 4-momenta \f$ p_3p_{12} = p_1p_3 + p2p3 \f$ 
       */
      inline double p3p12 
      ( const double    s     ,
        const double    s1    ,
        const double    s2    ) const
      { return p1p3 ( s , s1 , s2 ) + p2p3 ( s , s1 , s2 ) ; }
      // =====================================================================
      /** get product of 4-momenta \f$ p_3p_{13} = p_1p_3 + m^2_3 \f$ 
       */
      inline double p3p13 
      ( const double    s     ,
        const double    s1    ,
        const double    s2    ) const
      { return p1p3 ( s , s1 , s2 ) + m3sq () ; }
      // =====================================================================
      /** get product of 4-momenta \f$ p_3p_{23} = p_2p_3 + m^2_3 \f$ 
       */
      inline double p3p23 
      ( const double    s     ,
        const double    s1    ,
        const double    s2    ) const
      { return p2p3 ( s , s1 , s2 ) + m3sq () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// variable transformation \f$ (s,s_1,s_2) \leftrigtharrow (s,x_1,x_2) \f$
      // ======================================================================
      /// the first x-varible is just \f$ x_1 = \cos_{R23}(12) \f$ 
      double x1 ( const double    s    , const double    s1    , const double s2 ) const ;
      /// the second x-variable is  \f$  x_2 = s_2 \f$ 
      double x2 ( const double /* s */ , const double /* s1 */ , const double s2 ) const
      { return s2 ;}
      // ======================================================================
      /** (inverse) variable transformation  
       *   \f[ \begin{array{l} 
       *        s_1 = f_1 ( x_1 ,x_2  )  \\ 
       *        s_2 = f_2 ( x_1 ,x_2  )
       *       \end{array}\f]
       *   where 
       *   \f[ \begin{array{l} 
       *        x_1 = \cos_{R23)(12)  \\ 
       *        x_2 = s_2 
       *       \end{array} \f]
       *  @code
       *  Dalitz d = ... ;
       *  const double x1 = 0.1  ;
       *  const double x2 = 14.5 ;
       *  double s1, s2 ;
       *   std::tie ( s1, s2 ) = d.x2s ( s , x1 , x2 ) ; 
       *  @endcode  
       */
      std::pair<double,double> x2s
      ( const double s  ,
        const double x1 ,
        const double x2 ) const ;
      // =======================================================================
      /**  absolute value of the jacobian  
       *   \f$ J( s,s_1,s_2) = \left| \frac{\partial(s_1,s_2) }{\partial(x_1,x_2)} \right| \f$ 
       */
      double J
      ( const double s  ,
        const double s1 ,
        const double s2 ) const ;
      // ======================================================================
      /// variable transformation \f$ (s,s_1,s_2) \leftrigtharrow (s_2,y_1,y_2) \f$
      // ======================================================================
      /// the first  y-varible is just \f$ y_1 = s \f$ 
      double y1 ( const double s , const double /* s1 */ , const double /* s2 */ ) const
      { return s ; }
      /// the second y-varible is just \f$ y_2 = \cos_{R23}(12) \f$ 
      double y2 ( const double s  , const double   s1    , const double s2 ) const
      { return x1 ( s , s1  , s2 ) ; }
      // =======================================================================
      /** (inverse) variable transformation  
       *   \f[ \begin{array{l} 
       *        s   = f_1 ( y_1 ,y_2  )  \\ 
       *        s_1 = f_2 ( y_1 ,y_2  )
       *       \end{array}\f]
       *   where 
       *   \f[ \begin{array{l} 
       *        y_1 = s \\  
       *        y_2 = \cos_{R23)(12) 
       *       \end{array} \f]
       *  @code
       *  Dalitz d = ... ;
       *  const double y1 = 0.1  ;
       *  const double y2 = 0.5 ;
       *  double s, s1 ;
       *   std::tie ( s , s1 ) = d.y2s ( s2 , y1 , y2 ) ; 
       *  @endcode  
       */
      std::pair<double,double> y2s
      ( const double s2 ,
        const double y1 ,
        const double y2 ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// Dalitz plot boundaries \f$ s_1^{min/max} ( s, s_2 ) \f$ 
      std::pair<double,double>  s1_minmax_for_s_s2
      ( const double s  ,
        const double s2 ) const ;
      // ======================================================================
      /// Dalitz plot boundaries \f$ s_2^{min/max} ( s, s_1 ) \f$ 
      std::pair<double,double>  s2_minmax_for_s_s1
      ( const double s  ,
        const double s1 ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /** "transpose it", such that \f$ s_{i1} \f$ and \f$ s_{i2}\f$ 
       *  becomes  the main  variable
       */
      Dalitz0 transpose 
      ( const unsigned short i1 , 
        const unsigned short i2 ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag/hash value  
      std::size_t tag () const { return m_tag  ; }
      // ======================================================================
    private :
      // ======================================================================
      /// the first mass 
      double m_m1 ;                                      //      the first mass 
      /// the second mass 
      double m_m2 ;                                      //     the second mass 
      /// the third mass 
      double m_m3 ;                                      //      the third mass 
      /// tag/hash value 
      std::size_t          m_tag    ;                         // tag/hash value 
      /// the precalculated quantities 
      // std::array<double,9> m_cache  ;            // the precalculated quantities 
      // std::array<int,7>    m_cacheb ;            // the precalculated quantities
      double m_cache  [9] ;                     // the precalculated quantities 
      bool   m_cacheb [7] ;                     // the precalculated quantities
      // ======================================================================      
    } ;
    // ========================================================================
    /** @class Dalitz Ostap/Dalitz.h 
     *  Simple Kinematics of Dalitz plot 
     *  @see E.Byckling, K.Kajantie, "Particle kinematics", John Wiley & Sons,
     *              London, New York, Sydney, Toronto, 1973, p.89, eq. (5.23)     
     *  @see Section V.1
     *  @see https://userweb.jlab.org/~rafopar/Book/byckling_kajantie.pdf     
     */
    class Dalitz : public Dalitz0
    {
    public :
      // ======================================================================
      /** constructor from all masses 
       *  - M  : overall mass of the system, \f$\sqrt{s}\f$;
       *  - m1 : the mass of the first particle  \f$ m_1 \f$;
       *  - m2 : the mass of the second particle  \f$ m_2 \f$;
       *  - m3 : the mass of the third particle  \f$ m_3 \f$;
       */
      Dalitz ( const double M  = 1 , 
               const double m1 = 0 , 
               const double m2 = 0 , 
               const double m3 = 0 ) : Dalitz ( M , Dalitz0 ( m1 , m2 , m3 ) ) {}
      // ======================================================================
      /** constructor from all masses 
       *  @param m overlal mass of the system, \f$\sqrt{s}\f$;
       *  @param b individual masses 
       */
      Dalitz ( const double   M ,
               const Dalitz0& b ) ;
      // ======================================================================
      /** constructor from all masses 
       *  @param b individual masses 
       *  @param m overlal mass of the system, \f$\sqrt{s}\f$;
       */
      Dalitz ( const Dalitz0& b , 
               const double   M ) : Dalitz ( M , b ) {}
      // ======================================================================
    public:  // trivial getters 
      // ====================================================================== 
      /// \f$ s =  (p_1 + p^2 + p^3)^2 \f$
      double s      () const  { return S ()  ; }
      double sqs    () const  { return m_M   ; }
      ///   total mass 
      double M      () const  { return m_M   ; }
      // ======================================================================
      using Dalitz0::s1_max ;
      using Dalitz0::s2_max ;
      using Dalitz0::s3_max ;
      // ======================================================================
      /// maximal value of s1: \f$ \left. s_1\right|_{max} =  (M-m_3)^2 \f$
      double s1_max () const { return m_cache2 [ 0 ] ; }
      /// maximal value of s2: \f$ \left. s_1\right|_{max} =  (M-m_1)^2 \f$
      double s2_max () const { return m_cache2 [ 1 ] ; }
      /// maximal value of s3: \f$ \left. s_3\right|_{max} =  (M-m_2)^2 \f$
      double s3_max () const { return m_cache2 [ 2 ] ; }
      // ======================================================================
      /** Get the sum of all invariants 
       *    \f$ s_1 +  s_2 + s_3 = s_{12} + s_{23} + s_{31} = 
       *           s + m_1^2 + m_2^2 + m_3^2 \f$ 
       *  - The sum is precalculated 
       */
      double sums   () const { return m_cache2 [ 3 ] ; }
      // ======================================================================
      /// \f$ M^2\f$   precalculated 
      double Msq    () const { return m_cache2 [ 4 ] ; }
      /// \f$ s \f$    precalculated 
      double S      () const { return Msq ()      ; }
      // ======================================================================
    public:  // only two s_i are independent 
      // ======================================================================
      /** get s_1 for given s_2 and s_3 are defined 
       *  \f$ s_1 = s_{12} = \sum s -   s_2 - s_3  \f$ 
       */
      inline double s1 ( const double s2 , const double s3 ) const 
      { return sums () - s2 - s3 ; }
      // ======================================================================
      /** get s_2 for given s_1 and s_3 are defined 
       *  \f$ s_2 = s_{23} = \sum s -   s_1 - s_3  \f$ 
       */
      inline double s2 ( const double s1 , const double s3 ) const 
      { return sums () - s1 - s3 ; }
      // ======================================================================
      /** get s_3 for given s_1 and s_2 are defined 
       *  \f$ s_3 = s_{31} = \sum s -   s_1 - s_2  \f$ 
       */
      inline double s3 ( const double s1 , const double s2 ) const 
      { return sums () - s1 - s2 ; }
      // ======================================================================
    public: //invariants 
      // ======================================================================
      using Dalitz0::p1p2 ;
      using Dalitz0::p2p3 ;
      using Dalitz0::p1p3 ;
      using Dalitz0::pp1  ;
      using Dalitz0::pp2  ;
      using Dalitz0::pp3  ;
      // ======================================================================
      /** get the product of 4-momenta
       *  \f$ \left(p_1p_2\right) = \frac{1}{2}\left( s_{12} - m_1^2 - m_2^2 \right) \f$
       */
      inline double p1p2 ( const double s1 ,
                           const double s2 )  const
      { return p1p2 ( s () , s1 , s2 ) ; }
      // ======================================================================
      /** get the product of 4-momenta
       *  \f$ \left(p_2p_3\right) = \frac{1}{2}\left( s_{23} - m_2^2 - m_3^2 \right) \f$      
       */
      inline double p2p3 ( const double s1 ,
                           const double s2 ) const
      { return p2p3 ( s () , s1 , s2  ) ; }      
      // ======================================================================
      /** get the product of 4-momenta
       *  \f$ \left(p_1p_3\right) = \frac{1}{2}\left( s_{13} - m_1^2 - m_3^2 \right) \f$      
       */
      inline double p1p3 ( const double s1 ,
                           const double s2 ) const
      { return p1p3 ( s() , s1  , s2 ) ; }
      // ======================================================================
      /** get the product of 4-momenta
       *  \f$ \left(pp1\right) = \frac{1}{2}\left( s - s_{23} + m_1^2 \right) \f$
       */
      inline double pp1 ( const double s1 ,
                          const double s2 ) const
      { return pp1 ( s () , s1 , s2 )  ; }
      // ======================================================================
      /** get the product of 4-momenta
       *  \f$ \left(pp2\right) = \frac{1}{2}\left( s - s_{13} + m_2^2 \right) \f$
       */
      inline double pp2 ( const double s1 ,
                          const double s2 ) const
      { return pp2 ( s () , s1 , s2 )  ; }
      // ======================================================================
      /** get the product of 4-momenta
       *  \f$ \left(pp3\right) = \frac{1}{2}\left( s - s_{12} + m_3^2 \right) \f$
       */
      inline double pp3 ( const double s1 ,
                          const double s2 ) const
      { return pp3 ( s () , s1 , s2 )  ; }
      // ======================================================================
    public: 
      // ======================================================================
      /**  get the energies/momenta of particles in the rest-frame 
       *   @see Eq. (V.1.3) in E.Byckling, K.Kajantie, "Particle kinematics", 
       */
      // ======================================================================
      /// Energy of the 1st particle 
      inline double E1 ( const double /* s1 */ , const double s2        ) const 
      { return  ( s () + m1sq () - s2             ) / ( 2 *  sqs () ) ; }
      // ============-=========================================================
      /// Energy of the 2nd particle 
      inline double E2 ( const double s1       , const double s2        ) const 
      { return  ( s () + m2sq () - s3 ( s1 , s2 ) ) / ( 2 *  sqs () ) ; }
      // ======================================================================
      /// Energy of the 3rd particle 
      inline double E3 ( const double s1       , const double /* s2 */  ) const 
      { return  ( s () + m3sq () - s1             ) / ( 2 *  sqs () ) ; }
      // ======================================================================
      /// momentum of the 1st particle 
      double P1 ( const double s1  , const double s2 ) const ;
      /// momentum of the 2nd particle 
      double P2 ( const double s1  , const double s2 ) const ;
      /// momentum of the 3rd particle 
      double P3 ( const double s1  , const double s2 ) const ;
      // ======================================================================
      /// kinetic energy of 1st particle 
      inline double T1 ( const double s1 , const double s2        ) const 
      { return E1 ( s1 , s2 ) -  m1 () ; }
      /// kinetic energy of 2nd particle 
      inline double T2 ( const double s1 , const double s2        ) const 
      { return E2 ( s1 , s2 ) -  m2 () ; }
      /// kinetic energy of 3rd particle 
      inline double T3 ( const double s1 , const double s2        ) const 
      { return E3 ( s1 , s2 ) -  m3 () ; }          
      // ======================================================================
    public:  // limits for E1 , E2 and E3
      // ======================================================================
      double E1_max () const { return m_cache2 [ 5 ] ; }
      double E2_max () const { return m_cache2 [ 6 ] ; }
      double E3_max () const { return m_cache2 [ 7 ] ; }
      // ======================================================================
    public:
      // ======================================================================
      /**  \f$ \theta^{*}_{12}\f$ is angle between 
       *  \f$ p_1\f$  and \f$ p_2 \f$  in the rest frame:
       *   \f$ \cos \theta^{*}_{12} = 
       *   \left.\frac { p_1p_2}{P_1P_2}\right|_{\vec{P}=0}\f$
       */
      double cos_12  ( const double s1 , const double s2 ) const ;
      // =====================================================================
      /**  \f$ \theta^{*}_{23}\f$ is angle between 
       *  \f$ p_2\f$  and \f$ p_3 \f$  in the rest frame:
       *   \f$ \cos \theta^{*}_{23} = 
       *   \left.\frac { p_2p_3}{P_2P_3}\right|_{\vec{P}=0}\f$
       */
      double cos_23  ( const double s1 , const double s2 ) const ;
      // ======================================================================
      /**  \f$ \theta^{*}_{31}\f$ is angle between 
       *  \f$ p_3\f$  and \f$ p_1 \f$  in the rest frame:
       *   \f$ \cos \theta^{*}_{31} = 
       *   \left.\frac { p_3p_1}{P_3P_1}\right|_{\vec{P}=0}\f$
       */
      double cos_31  ( const double s1 , const double s2 ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /**  \f$ \theta^{*}_{12}\f$ is angle between 
       *  \f$ p_1\f$  and \f$ p_2 \f$  in the rest frame:
       *   \f$ \sin^2 \theta^{*}_{12} = 
       *   -4s \frac{ G(s_1,s_2,s, m_2^2, m_1^2,m_3^2) }
       *   {\lambda(s, m_1^2, s2) \lambda(s, m_2^2, s3 ) } \f$
       */
      double sin2_12 ( const double s1 , const double s2 ) const ;
      // ======================================================================
      /**  \f$ \theta^{*}_{23}\f$ is angle between 
       *  \f$ p_2\f$  and \f$ p_3 \f$  in the rest frame:
       *   \f$ \sin^2 \theta^{*}_{23} = 
       *   -4s \frac{ G(s_2,s_3,s, m_3^2, m_2^2,m_1^2) }
       *   {\lambda(s, m_2^2, s3 ) \lambda(s, m_3^2, s1 ) } \f$
       */
      double sin2_23 ( const double s1 , const double s2 ) const ;
      // ======================================================================
      /**  \f$ \theta^{*}_{31}\f$ is angle between 
       *  \f$ p_3\f$  and \f$ p_1 \f$  in the rest frame:
       *   \f$ \sin^2 \theta^{*}_{31} = 
       *   -4s \frac{ G(s_3,s_1,s, m_1^2, m_3^2,m_2^2) }
       *   {\lambda(s, m_3^2, s1 ) \lambda(s, m_1^2, s2 ) } \f$
       */
      double sin2_31 ( const double s1 , const double s2 ) const ;
      // ======================================================================
    public: // (1,2)-rest frame 
      // ======================================================================
      ///  full energy in (1,2)-rest frame 
      inline double E_R12  ( const double s1 , const double /* s2 */ ) const    
      { return ( s () + s1 - m3sq ()  )     / ( 2 * std::sqrt ( s1 ) ) ; }      
      /// energy of 1st  particle in (2,3)-rest frame 
      inline double E1_R12 ( const double s1 , const double /* s2 */ ) const    
      { return ( s1   + m1sq () - m2sq () ) / ( 2 * std::sqrt ( s1 ) ) ; }
      /// energy of 2nd  particle in (2,3)-rest frame 
      inline double E2_R12 ( const double s1 , const double /* s2 */ ) const    
      { return ( s1   + m2sq () - m1sq () ) / ( 2 * std::sqrt ( s1 ) ) ; }
      /// energy of 3rd  particle in (2,3)-rest frame 
      inline double E3_R12 ( const double s1 , const double /* s2 */ ) const    
      { return ( s () - s1 - m3sq ()  )     / ( 2 * std::sqrt ( s1 ) ) ; }
      /// total momentum in (1,2)-rest frame
      double P_R12   ( const double    s1 , const double s2 ) const ;
      /// momentum of 3rd particle in (1,2)-rest frame 
      double P3_R12  ( const double    s1    , const double s2 ) const 
      { return P_R12  ( s1 , s2 ) ; }
      /// momentum of 1st particle in (1,2)-rest frame 
      double P1_R12  ( const double    s1    , const double s2 ) const ;
      /// momentum of 2nd particle in (1,2)-rest frame 
      double P2_R12  ( const double    s1    , const double s2 ) const 
      { return P1_R12 ( s1 , s2 ) ; }
      /** cosine on the angle between 3rd and 1st particles in the  (1,2) rest frame
       *  \f$ \cos \theta_{31}^{R(1,2)}
       * - the same as <c>Ostap::Kinematics::cosThetaRest (p3,p1,p1+p2)</c>
       * - the same as <c>Ostap::Kinematics::cos_theta    (p3,p1,p1+p2)</c>
       * - the same as <c>Ostap::Kinematics::decayAngle   (p1,p1+p2,-p3)</c>
       */
      double cos_31_R12 ( const double    s1    , const double s2 ) const ;
      /** sine squared  on the angle between 3rd and 1st particles in the  (1,2) rest frame
       *  \f$ \sin^2 \theta_{31}^{R(1,2)}
       */
      double sin2_31_R12 ( const double   s1    , const double s2 ) const ;
      // ======================================================================
    public: // (2,3)-rest frame 
      // ======================================================================
      ///  full energy in (2,3)-rest frame 
      inline double E_R23  ( const double /* s1 */ , const double s2 ) const    
      { return ( s () + s2 - m1sq () )        / ( 2 * std::sqrt ( s2 ) ) ; }
      /// energy of 1st  particle in (2,3)-rest frame 
      inline double E1_R23 ( const double /* s1 */ , const double s2 ) const    
      { return ( s () - s2 - m1sq ()  )       / ( 2 * std::sqrt ( s2 ) ) ; }
      /// energy of 2nd  particle in (2,3)-rest frame 
      inline double E2_R23 ( const double /* s1 */ , const double s2 ) const    
      { return ( s2   +  m2sq () - m3sq ()  ) / ( 2 * std::sqrt ( s2 ) ) ; }
      /// energy of 3rd  particle in (2,3)-rest frame 
      inline double E3_R23 ( const double /* s1 */ , const double s2 ) const    
      { return ( s2   +  m3sq () - m2sq () ) / ( 2 * std::sqrt ( s2 ) ) ; }
      /// total momentum in (2,3)-rest frame
      double P_R23   ( const double /* s1 */ , const double s2 ) const ;
      /// momentum of 1st particle in (2,3)-rest frame 
      double P1_R23  ( const double    s1    , const double s2 ) const 
      { return P_R23  ( s1 , s2 ) ; }
      /// momentum of 2nd particle in (2,3)-rest frame 
      double P2_R23  ( const double    s1    , const double s2 ) const ;
      /// momentum of 3rd particle in (2,3)-rest frame 
      double P3_R23  ( const double    s1    , const double s2 ) const 
      { return P2_R23 ( s1 , s2 ) ; }
      // ======================================================================
      /** cosine on the angle between 1st and 2nd particles in the  (2,3) rest frame
       *  \f$ \cos \theta_{12}^{R(2,3)}
       * - the same as <c>Ostap::Kinematics::cosThetaRest (p1,p2,p2+p3)</c>
       * - the same as <c>Ostap::Kinematics::cos_theta    (p1,p2,p2+p3)</c>
       * - the same as <c>Ostap::Kinematics::decayAngle   (p1,p2+p3,-p2)</c>
       */
      double cos_12_R23  ( const double    s1    , const double s2 ) const ;
      /** sine squared  on the angle between 1st and 2nd particles in the  (2,3) rest frame
       *  \f$ \sin^2 \theta_{12}^{R(2,3)}
       */
      double sin2_12_R23 ( const double    s1    , const double s2 ) const ;
      // ======================================================================
    public: // (3,1)-rest frame 
      // ======================================================================
      ///  full energy in (3,1)-rest frame 
      inline double E_R31  ( const double s1 , const double s2 ) const    
      { 
        const double s3_ = s3 ( s1 , s2 ) ;
        return ( s () + s3_ - m2sq ()  )     / ( 2 * std::sqrt ( s3_) ) ; 
      }
      /// energy of 1st  particle in (3,1)-rest frame 
      inline double E1_R31 ( const double  s1 , const double s2 ) const    
      { 
        const double s3_ = s3 ( s1 , s2 ) ;
        return ( s3_ +  m1sq ()  - m3sq () ) / ( 2 * std::sqrt ( s3_ ) ) ;
      }
      /// energy of 2nd  particle in (3,1)-rest frame 
      inline double E2_R31 ( const double s1 , const double s2 ) const    
      { 
        const double s3_ = s3 ( s1 , s2 ) ;
        return ( s () - s3_ - m2sq ()  )    / ( 2 * std::sqrt ( s3_ ) ) ; 
      }
      /// energy of 3rd  particle in (2,3)-rest frame 
      inline double E3_R31 ( const double s1 , const double s2 ) const    
      { 
        const double s3_ = s3 ( s1 , s2 ) ;
        return ( s3_ +  m3sq ()  - m1sq () ) / ( 2 * std::sqrt ( s3_ ) ) ; 
      }
      /// total momentum in (3,1)-rest frame
      double P_R31   ( const double /* s1 */ , const double s2 ) const ;
      /// momentum of 2nd particle in (3,1)-rest frame 
      double P2_R31  ( const double    s1    , const double s2 ) const 
      { return P_R31  ( s1 , s2 ) ; }
      /// momentum of 1st particle in (2,3)-rest frame 
      double P3_R31  ( const double    s1    , const double s2 ) const ;
      /// momentum of 3rd particle in (2,3)-rest frame 
      double P1_R31  ( const double    s1    , const double s2 ) const 
      { return P3_R31 ( s1 , s2 ) ; }
      /** cosine on the angle between 2nd and 3rd particles in the  (3,1) rest frame
       *  \f$ \cos \theta_{23}^{R(3,1)}
       * - the same as <c>Ostap::Kinematics::cosThetaRest (p2,p3,p1+p3)</c>
       * - the same as <c>Ostap::Kinematics::cos_theta    (p2,p3,p1+p3)</c>
       * - the same as <c>Ostap::Kinematics::decayAngle   (p2,p1+p3,-p3)</c>
       */
      double cos_23_R31 ( const double    s1    , const double s2 ) const ;
      /** sine squared  on the angle between 2nd and 3rd particles in the  (3,1) rest frame
       *  \f$ \sin^2 \theta_{23}^{R(3,1)}
       */
      double sin2_23_R31 ( const double    s1    , const double s2 ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /** Dalitz plot density:
       *  \f$ R_3 = 
       *  \frac{1}{32s} \int {\mathrm{d}} s_1 {\mathrm{d}} s_2 {\mathrm{d}} Omega {\mathrm{d}} \phi_3 
       *  \Theta { -G \left( s_1, s_2 , s , m_2^2, m_1^2, m_3^2 \right) } \f$
       *  Here we take 
       *  \f$ \int {\mathrm{d}} \Omega = 4\pi\f$ and \f$ \int {\mathrm{d}} \phi_3 = 2\pi\f$. 
       */
      double density       ( const double s1  , const double s2  ) const ;
      // ======================================================================
      /** Dalitz density as  function of masses 
       *  \f$m_{12},m_{23}=\sqrt{s_1},\sqrt{s_2} \f$
       */
      double density_mass  ( const double m12 , const double m23 ) const 
      { return 
          m12 < ( m1 () + m2 () ) ? 0 : m12 > ( m_M  - m3 () ) ? 0 :
          m23 < ( m2 () + m3 () ) ? 0 : m23 > ( m_M  - m1 () ) ? 0 : 
          4 * m12 * m23 * density ( m12 * m12 , m23 * m23 ) ;
      }
      // ======================================================================
    public:   // 1-dimension 
      // ======================================================================
      /** Dalitz density in 1-dimension:
       *  \f$  \frac{{\mathrm{d}} R_3}{{\mathrm{d}} s_2} = 
       *  \frac{\pi^2}{2ss_2}
       *  \lambda^{1/2} ( s_2, s, m_1^2) 
       *  \lambda^{1/2} ( s_2, m_2^2, m_3^2) 
       *  \f$ 
       */
      double dRds2   ( const double s1 ) const ;
      // ======================================================================
      /** Dalitz density in 1-dimension:
       *  \f$  \frac{{\mathrm{d}} R_3}{{\mathrm{d}} s_3} = 
       *  \frac{\pi^2}{2ss_3}
       *  \lambda^{1/2} ( s_3, s, m_2^2) 
       *  \lambda^{1/2} ( s_3, m_3^2, m_1^2) 
       *  \f$ 
       */
      double dRds3   ( const double s3 ) const ;
      // ======================================================================
      /** Dalitz density in 1-dimension:
       *  \f$  \frac{{\mathrm{d}} R_3}{{\mathrm{d}} s_1} = 
       *  \frac{\pi^2}{2ss_1}
       *  \lambda^{1/2} ( s_1, s, m_3^2) 
       *  \lambda^{1/2} ( s_1, m_1^2, m_2^2) 
       *  \f$ 
       */
      double dRds1   ( const double s2 ) const ;
      // ======================================================================
    public: // Dalitz density as function of masses 
      // ======================================================================
      /// Dalitz density as function of \f$ m_12 = \sqrt{s_1}\f$ 
      inline double dRdm12 ( const double m12 ) const
      { return
          m12 < ( m1 () + m2 () ) ? 0 : m12 > ( m_M  - m3 () ) ? 0 : 
          2 * m12 * dRds1 ( m12 * m12 ) ; }
      /// Dalitz density as function of \f$ m_23 = \sqrt{s_2}\f$ 
      inline double dRdm23 ( const double m23 ) const
      { return 
          m23 < ( m2 () + m3 () ) ? 0 : m23 > ( m_M  - m1 () ) ? 0 :
          2 * m23 * dRds2 ( m23 * m23 ) ; }
      /// Dalitz density as function of \f$ m_31 = \sqrt{s_3}\f$ 
      inline double dRdm31 ( const double m31 ) const
      { return 
          m31 < ( m3 () + m1 () ) ? 0 : m31 > ( m_M  - m2 () ) ? 0 : 
          2 * m31 * dRds3 ( m31 * m31 ) ; }      
      // ======================================================================
    public :  // geometry of Dalitz plot 
      // ======================================================================
      using Dalitz0::inside   ;
      using Dalitz0::distance ;
      // ======================================================================
      /** Is the point \f$ (s_1,s_2)\f$ "inside" the Dalitz plot?
       *  Get the sign of G-function 
       *  \f$ g(s_1,s_2) = G ( s_1, s_2 , s , m_2^2, m_1^2, m_3^2) \f$
       *  @see Ostap::Kinematics::G
       *  Physical region corresponds to \f$ g\le0 \f$  
       */
      bool   inside   ( const double s1 , const double s2 ) const ;
      // =======================================================================
      /** get the measure of the distance form the point to the boundary of Dalitz plot. 
       *  Distance is defined as 
       *  \f$ d \equiv = \lambda ( P_1^2 , P_2^2 , P_3^2 ) \f$ 
       */
      double distance ( const double s1 , const double s2 ) const ;
      // ======================================================================
      /// Dalitz plot boundaries \f$ s_1^{min/max} ( s_2 ) \f$ 
      std::pair<double,double>  s1_minmax_for_s2 ( const double s2 ) const
      { return s1_minmax_for_s_s2 ( s () , s2 ) ; }
      // ======================================================================
      /// Dalitz plot boundaries \f$ s_2^{min/max} ( s_1 ) \f$ 
      std::pair<double,double>  s2_minmax_for_s1 ( const double s1 ) const 
      { return s2_minmax_for_s_s1 ( s () , s1 ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag/hash value  
      std::size_t tag () const { return m_tag2 ; }
      // ======================================================================
    private :
      // ======================================================================
      /// the total mass 
      double m_M ;                                       //      the total mass
      /// the precalcualted quantities 
      // std::array<double,8> m_cache2 ;           // the precalcualted quantities
      double m_cache2 [8] ;                     // the precalcualted quantities
      /// tag/hash value
      std::size_t          m_tag2   ;                         // tag/hash value
      // ======================================================================      
    } ;
    // ========================================================================
    /** Get a full integrated phase space over Dalitz plot 
     *  \f$  R(s) = \int \int R(s_1,s_2) {\mathrm{d}} s_1 {\mathrm{d}} s_2 =
     *  \int _{(m_2+m_3)^2}^{ (\sqrt{s}-m_1)^2}
     *  \frac{{\mathrm{d}} s_2}{s_2}
     *  \lambda^{1/2}(s_2,s,m_1^2)
     *  \lambda^{1/2}(s_2,m_2^2,m_3^2)\f$ 
     */
    double phase_space ( const Dalitz& dalitz ) ;    
    // ========================================================================
  } //                                   The end of namespace Ostap::Kinematics 
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_DALITZ_H
// ============================================================================
