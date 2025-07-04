// ============================================================================
#ifndef OSTAP_SVECTORWITHERROR_H 
#define OSTAP_SVECTORWITHERROR_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL 
// ============================================================================
#include <iosfwd>
// ============================================================================
// ROOT
// ============================================================================
#include "Math/SVector.h"
#include "Math/SMatrix.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/ValueWithError.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    /** @class SVectorWithError Ostap/SVectorWithError.h
     *  Simple class with represent SVector with 
     *  the associated covariance matrix
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     */
    template <unsigned int N, typename SCALAR=double>
    class SVectorWithError
    {
      // ======================================================================
    public:
      // ======================================================================
      /// the actual type of data
      typedef ROOT::Math::SVector<SCALAR,N>                             Value ;
      /// the actual type of covarinace matrix
      typedef ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> > Covariance ;
      /// self type 
      typedef SVectorWithError<N,SCALAR>                                Self  ;
      // ======================================================================
    public:
      // ======================================================================
      enum
        {
          /// vector size
          kSize = N // vector size 
        } ;  
      // ======================================================================
    public:
      // ======================================================================
      /// full constructor from vector and covariance matrix 
      SVectorWithError 
      ( const Value&      value = Value()       , 
        const Covariance& cov2  = Covariance () ) 
        : m_value ( value ) 
        , m_cov2  ( cov2  ) 
      {}
      /// full constructor from covariance matrix 
      SVectorWithError ( const Covariance& cov2 ) 
        : m_value (       ) 
        , m_cov2  ( cov2  ) 
      {}
      /// construct from expressions
      template <class B, class R>
      SVectorWithError 
      ( const Value&                            value , 
        const ROOT::Math::Expr<B,SCALAR,N,N,R>& cov2  ) 
        : m_value ( value )
        , m_cov2  ( cov2  ) 
      {}
      /// constructor from different scalar types 
      template <class B>
      SVectorWithError 
      ( const ROOT::Math::VecExpr<B,SCALAR,N>& value , 
        const Covariance&                      cov2  ) 
        : m_value ( value )
        , m_cov2  ( cov2  ) 
      {}
      /// construct from expressions
      template <class B1, class B2, class R>
      SVectorWithError 
      ( const ROOT::Math::VecExpr<B1,SCALAR,N>&  value , 
        const ROOT::Math::Expr<B2,SCALAR,N,N,R>& cov2  ) 
        : m_value ( value )
        , m_cov2  ( cov2  ) 
      {}
      /// construct from expressions
      template <class B>
      SVectorWithError 
      ( const ROOT::Math::VecExpr<B,SCALAR,N>& value ) 
        : m_value ( value )
        , m_cov2  (       ) 
      {}
      /// construct from expressions
      template <class B, class R>
      SVectorWithError 
      ( const ROOT::Math::Expr<B,SCALAR,N,N,R>& cov2 ) 
        : m_value (       )
        , m_cov2  ( cov2  ) 
      {}
      // ======================================================================
    public: // trivial accessors 
      // ======================================================================
      const  Value&      value       () const { return m_value      ; }
      const  Covariance& cov2        () const { return m_cov2       ; }
      const  Covariance& covariance  () const { return this->cov2() ; }
      // ======================================================================
      inline Value&      value       ()       { return m_value      ; }
      inline Covariance& cov2        ()       { return m_cov2       ; }      
      inline Covariance& covariance  ()       { return this->cov2() ; }
      // ======================================================================
      const  SCALAR& value
      ( unsigned int i ) const  { return m_value ( i )     ; }
      const  SCALAR& cov2
      ( unsigned int i , 
        unsigned int j ) const  { return m_cov2  ( i , j ) ; }
      // ======================================================================
      inline SCALAR& value
      ( unsigned int i )        { return m_value ( i )     ; }
      inline SCALAR& cov2
      ( unsigned int i , 
        unsigned int j )        { return m_cov2  ( i , j ) ; }      
      // ======================================================================
    public:  // finally it is just a vector 
      // ======================================================================
      const  SCALAR& operator()
        ( unsigned int i ) const { return m_value(i) ; }
      inline SCALAR& operator()
        ( unsigned int i )       { return m_value(i) ; }
      const  SCALAR& operator[]
      ( unsigned int i ) const   { return m_value[i] ; }
      inline SCALAR& operator[]
      ( unsigned int i )         { return m_value[i] ; }
      const  SCALAR& operator()
        ( unsigned int i ,
          unsigned int j ) const { return m_cov2(i,j) ; }
      inline SCALAR& operator()
        ( unsigned int i ,
          unsigned int j )       { return m_cov2(i,j) ; }      
      // ======================================================================
    public: // 
      // ======================================================================
      /// set value 
      void setValue
      ( const unsigned int i     ,
        const SCALAR       value ) { m_value ( i )   = value ; }
      void setCov2
      ( const unsigned int i     ,
        const unsigned int j     ,
        const SCALAR       value ) { m_cov2 ( i , j ) = value ; }
      /// set value 
      void set
      ( const unsigned int i     ,
        const SCALAR       value ) { setValue ( i , value ) ; }
      void set
      ( const unsigned int i     ,
        const unsigned int j     ,
        const SCALAR       value ) { setCov2 ( i , j , value ) ; }
      // ======================================================================
    public: // correlations 
      // ======================================================================
      /** get the correlation coefficient between "i" and "j"
       *  for invalid setup , return large negative value 
       *  @param i the first index 
       *  @param j the second index 
       *  @return correlation coefficient 
       */
      inline SCALAR  corr  ( unsigned int i  , unsigned int j  ) const ;
      /** get the full correlation matrix 
       *  @return false for invalid setup 
       */
      inline bool    corr  ( Covariance& corrm ) const ;
      // ======================================================================
    public: // setters 
      // ======================================================================
      void setValue      ( const Value&      v ) { m_value = v ; }
      void setCovariance ( const Covariance& c ) { m_cov2  = c ; }
      void setCov2       ( const Covariance& c ) { m_cov2  = c ; }
      // ======================================================================
      template <class B>
      void setValue       ( const ROOT::Math::VecExpr<B,SCALAR,N>&  v ) 
      { m_value = v ; }
      template <class B, class R>
      void setCovariance  ( const ROOT::Math::Expr<B,SCALAR,N,N,R>& c ) 
      { m_cov2 = c ; }
      template <class B, class R>
      void setCov2        ( const ROOT::Math::Expr<B,SCALAR,N,N,R>& c ) 
      { m_cov2 = c ; }
      // ======================================================================
    public: // cast:
      // ======================================================================
      operator const Value&      () const { return value () ; }
      operator       Value&      ()       { return value () ; }
      operator const Covariance& () const { return cov2  () ; }
      operator       Covariance& ()       { return cov2  () ; }
      // ======================================================================
    public: //  operators
      // ======================================================================
      Self& operator+= ( const Self&   right ) 
      { m_value += right.m_value ; m_cov2 += right.m_cov2 ; return *this ; }
      Self& operator-= ( const Self&   right ) 
      { m_value -= right.m_value ; m_cov2 += right.m_cov2 ; return *this ; }
      Self& operator+= ( const Value&  right ) 
      { m_value += right                                  ; return *this ; }
      Self& operator-= ( const Value&  right ) 
      { m_value -= right                                  ; return *this ; }
      Self& operator*= ( const double  s     ) 
      { m_value *= s             ; m_cov2 *= (s*s)        ; return *this ; }
      Self& operator/= ( const double  s     ) 
      { m_value /= s             ; m_cov2 /= (s*s)        ; return *this ; }
      // ======================================================================
      template <class B>
      Self&  operator+= ( const ROOT::Math::VecExpr<B,SCALAR,N>&  right ) 
      { m_value += right ; return *this ; }
      template <class B>
      Self&  operator-= ( const ROOT::Math::VecExpr<B,SCALAR,N>&  right ) 
      { m_value -= right ; return *this ; }
      // ======================================================================
      // unary- 
      Self operator-() const { return Self ( -m_value , m_cov2 ) ; }  // unary- 
      // ======================================================================
    public: //  chi2 distances
      // ======================================================================
      double chi2 ( const Self&  right ) const ;
      double chi2 ( const Value& right ) const ;      
      template <class B>
      double chi2 ( const ROOT::Math::VecExpr<B,SCALAR,N>&  right ) const ;
      // ======================================================================
    public:  // more functions 
      // ======================================================================
      /// calculate the weighted average for two vectors 
      inline Self mean    ( const Self& right ) const ;
      /// calculate the weighted average for two vectors 
      inline Self average ( const Self& right ) const { return mean ( right ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /** Get symmetrized Kullback-Leibler divergency,
       *  aka Jewffrey's divergency, for two objects 
       *  @see https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence
       *  @see Ostap::Math::kullback_leibler 
       *  @return symmetrised KL-divergency (-1 in case of error)
       */
      double kullback_leibler 
      ( const SVectorWithError& a ) const ;
      // ========================================================================
      /** Get asymmetric Kullback-Leibler divergency for two objects 
       *  @see https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence
       *  @see Ostap::Math::asymmetric_kullback_leibler 
       *  @return KL-divergency (-1 in case of error)
       */
      double asymmetric_kullback_leibler 
      ( const SVectorWithError& a ) const ;
      // ========================================================================
    public:
      // ========================================================================
      /** get Mahalanobis distance 
       *  https://en.wikipedia.org/wiki/Mahalanobis_distance  
       *  @return Mahalanobis distance (-1 in case of error)
       */
      double mahalanobis 
      ( const SVectorWithError& a ) const ;
      // ========================================================================
      /** get Mahalanobis distance 
       *  https://en.wikipedia.org/wiki/Mahalanobis_distance  
       *  @return Mahalanobis distance (-1 in case of error)
       */
      double mahalanobis 
      ( const Value& a ) const ;
      // ========================================================================
    public:
      // ========================================================================
      /** get the (unnornalized) weighted sum with set of weights 
       *  \f[ r = \sum_i v_i w_i \f]
       *  @param weigths (INPUT) vector of weights 
       */ 
      Ostap::Math::ValueWithError
      dot ( const SVectorWithError& weights ) const ;
      // ========================================================================
      /** get the (unnormalized) weighted sum with set of weights 
       *  \f[ r = \sum_i v_i w_i \f]
       *  @param weigths (INPUT) vector of weights 
       */ 
      Ostap::Math::ValueWithError
      dot ( const Value&            weights ) const ;
      // ========================================================================
    public:
      // ========================================================================
      /** get the (normalized) weighted sum with set of weights 
       *  \f[ r = \frac{\sum_i v_i w_i}{ \sum_i w_i } \f]
       *  @param weigths (INPUT) vector of weights 
       */ 
      Ostap::Math::ValueWithError
      weighted_sum ( const SVectorWithError& weights ) const ;
      // ========================================================================
      /** get the (normalized) weighted sum with set of weights 
       *  \f[ r = \frac{\sum_i v_i w_i}{ \sum_i w_i } \f]
       *  @param weigths (INPUT) vector of weights 
       */ 
      Ostap::Math::ValueWithError
      weighted_sum ( const Value&            weights ) const ;
      // ========================================================================
    public: //  helper functions for pythonizations
      // ======================================================================
      Self  __add__     ( const Self&  right ) const ;
      Self  __sub__     ( const Self&  right ) const ;      
      Self  __add__     ( const Value& right ) const ;
      Self  __sub__     ( const Value& right ) const ;      
      Self  __radd__    ( const Value& right ) const ;
      Self  __rsub__    ( const Value& right ) const ;
      // ======================================================================
      Self  __mul__     ( const double v     ) const { return (*this) *  v ; }
      Self  __truediv__ ( const double v     ) const { return (*this) /  v ; }
      Self  __div__     ( const double v     ) const { return (*this) /  v ; }
      Self  __rmul__    ( const double v     ) const { return (*this) *  v ; }
      // ======================================================================
      Self& __imul__    ( const double v     )       { return (*this) *= v ; }
      Self& __idiv__    ( const double v     )       { return (*this) /= v ; }
      Self& __iadd__    ( const double v     )       { return (*this) += v ; }
      Self& __isub__    ( const double v     )       { return (*this) -= v ; }
      Self& __iadd__    ( const Self&  v     )       { return (*this) += v ; }
      Self& __isub__    ( const Self&  v     )       { return (*this) -= v ; }
      Self& __iadd__    ( const Value& v     )       { return (*this) += v ; }
      Self& __isub__    ( const Value& v     )       { return (*this) -= v ; }
      // ======================================================================
    public:
      // ======================================================================
      /** "transform" vector withj uncertaintie susnig matrix 
       *   @param M (INPUT)
       *   @return transformed vector
       */
      template <unsigned int K, typename R> 
      SVectorWithError<K,SCALAR>
      __rmul__ 
      ( const ROOT::Math::SMatrix<SCALAR,K,N,R>& M ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// transform it 
      // template <unsigned int K, typename R>
      // Ostap::Math::SVectorWithError<K,SCALAR>
      // transform ( )
      // const ROOT::Math::SMatrix<SCALAR,K,N,R>& L ) const      
      // ======================================================================      
    public: //  printout 
      // ======================================================================
      /// printpout 
      std::ostream& fillStream ( std::ostream& s ) const ;         // printpout 
      /// conversion to string
      std::string   toString   () const ;               // conversion to string
      // ======================================================================
    private:
      // ======================================================================
      /// the data 
      Value      m_value  ;                                         // the data 
      /// the covariance matrix
      Covariance m_cov2   ;                           // the covarinance matrix  
      // ======================================================================
    } ;
    // ========================================================================
    /// disable null-size vectors
    template <class SCALAR>
    class SVectorWithError<0,SCALAR> {} ;
    // ========================================================================
    /** Get element of the vector as ValueWithError Object 
     *  @code
     *  SVectorWithError<5> vct = ... ;
     *  get<1> ( vct ) ;
     *  @endcode 
     *  @see Ostap::Math::ValueWithError 
     */
    template <unsigned int I  ,  
              unsigned int N  , 
              typename SCALAR ,
              typename std::enable_if<(I<N),bool>::type = true > 
    inline 
    ValueWithError
    get
    ( const SVectorWithError<N,SCALAR>& v ) 
    { return ValueWithError ( v.value ( I ) , v.cov2 ( I , I ) ) ; }
    // ===============================================================
    /** Get element of the vector as ValueWithError Object 
     *  @code
     *  SVectorWithError<5> vct = ... ;
     *  get ( vct , 1 ) ;
     *  @endcode 
     *  @see Ostap::Math::ValueWithError 
     */
     template <unsigned int N  , typename SCALAR >
     inline 
     ValueWithError
     get
     ( const SVectorWithError<N,SCALAR>& v , 
       const unsigned short              i )
     { return i < N ? ValueWithError ( v.value ( i ) , v.cov2 ( i , i ) ) : ValueWithError()  ; }
    // ========================================================================
    /** Get element of the vector as ValueWithError Object 
     *  @code
     *  SVectorWithError<5> vct = ... ;
     *  get ( 1 , vct ) ;
     *  @endcode 
     *  @see Ostap::Math::ValueWithError 
     */
     template <unsigned int N  , typename SCALAR >
     inline 
     ValueWithError
     get
     ( const unsigned short              i , 
       const SVectorWithError<N,SCALAR>& v )
     { return i < N ? ValueWithError ( v.value ( i ) , v.cov2 ( i , i ) ) : ValueWithError()  ; }
    // ========================================================================
    // /// specific case for N=1 
    // template <> 
    // class SVectorWithError<1,double> : public Ostap::Math::ValueWithError
    // {
    // public:
    //   // =======================================================================
    //   /// the actual type of data
    //   typedef ROOT::Math::SVector<double,1>                                    Value       ;
    //   /// the actual type of covarinace matrix
    //   typedef ROOT::Math::SMatrix<double,1,1,ROOT::Math::MatRepSym<double,1> > Covariance ;
    //   // ======================================================================
    //   enum {
    //     /// vector size
    //     kSize = 1 // vector size 
    //   } ;  
    //   // ======================================================================
    // public:
    //   // ======================================================================
    //   SVectorWithError 
    //   ( const double value = 0 , 
    //     const double cov2  = 0 ) 
    //     : ValueWithError ( value , cov2 ) 
    //   {}
    //   // ======================================================================
    //   SVectorWithError
    //   ( const Value&      value , 
    //     const Covariance& cov2  )
    //     : SVectorWithError ( value [ 0 ] , cov2 ( 0 , 0 ) ) 
    //   {}
    //   // ======================================================================
    //   /// constructor from expressions 
    //   template <class B>
    //   SVectorWithError 
    //   ( const ROOT::Math::VecExpr<B,double,1>& value , 
    //     const Covariance&                      cov2  ) 
    //     : SVectorWithError ( Value ( value ) , cov2  ) 
    //   {}
    //   /// constructor from expressions
    //   template <class B, class R>
    //   SVectorWithError 
    //   ( const Value&                            value , 
    //     const ROOT::Math::Expr<B,double,1,1,R>& cov2  ) 
    //     : SVectorWithError ( value , Covariance ( cov2 ) ) 
    //   {}
    //   /// constructor from expressions
    //   template <class B1, class B2, class R>
    //   SVectorWithError 
    //   ( const ROOT::Math::VecExpr<B1,double,1>&  value , 
    //     const ROOT::Math::Expr<B2,double,1,1,R>& cov2  ) 
    //     : SVectorWithError ( Value ( value ) , cov2  ) 
    //   {}
    //   // ======================================================================      
    // } ;  
    // ========================================================================
    /// printout 
    template <unsigned int N, class SCALAR> 
    inline std::ostream& operator<<
    ( std::ostream& s , const SVectorWithError<N,SCALAR>& vct ) 
    { return vct.fillStream ( s ) ; }
    // ========================================================================
    template <unsigned int N, class SCALAR> 
    inline double chi2 
    ( const SVectorWithError<N,SCALAR>&      v1 , 
      const SVectorWithError<N,SCALAR>&      v2 ) { return v1.chi2 ( v2 )  ; }
    // ========================================================================
    template <unsigned int N, class SCALAR> 
    inline double chi2 
    ( const ROOT::Math::SVector<SCALAR,N>&   v1 , 
      const SVectorWithError<N,SCALAR>&      v2 ) { return v2.chi2 ( v1 )  ; }
    // ========================================================================
    template <unsigned int N, class SCALAR> 
    inline double chi2 
    ( const SVectorWithError<N,SCALAR>&       v1 , 
      const ROOT::Math::SVector<SCALAR,N>&    v2 ) { return v1.chi2 ( v2 )  ; }
    // ========================================================================
    template <unsigned int N, class SCALAR, class B> 
    inline double chi2 
    (  const ROOT::Math::VecExpr<B,SCALAR,N>& v1 , 
       const SVectorWithError<N,SCALAR>&      v2 ) { return v2.chi2 ( v1 )  ; }
    // ========================================================================
    template <unsigned int N, class SCALAR, class B> 
    inline double chi2 
    ( const SVectorWithError<N,SCALAR>&       v1 , 
      const ROOT::Math::VecExpr<B,SCALAR,N>&  v2 ) { return v1.chi2 ( v2 )  ; }
    // ========================================================================
    template <unsigned int N, class SCALAR>
    inline 
    SVectorWithError<N,SCALAR>
    operator+ 
    ( const SVectorWithError<N,SCALAR>& v1 , 
      const SVectorWithError<N,SCALAR>& v2 ) 
    {
      SVectorWithError<N,SCALAR> tmp ( v1 ) ;
      return tmp += v2 ;
    }
    // ========================================================================
    template <unsigned int N, class SCALAR>
    inline 
    SVectorWithError<N,SCALAR>
    operator- 
    ( const SVectorWithError<N,SCALAR>& v1 , 
      const SVectorWithError<N,SCALAR>& v2 ) 
    {
      SVectorWithError<N,SCALAR> tmp ( v1 ) ;
      return tmp -= v2 ;
    }
    // ========================================================================    
    template <unsigned int N, class SCALAR>
    inline 
    SVectorWithError<N,SCALAR>
    operator+ 
    ( const SVectorWithError<N,SCALAR>&    v1 , 
      const ROOT::Math::SVector<SCALAR,N>& v2 ) 
    {
      SVectorWithError<N,SCALAR> tmp ( v1 ) ;
      return tmp += v2 ;
    }
    // ========================================================================
    template <unsigned int N, class SCALAR>
    inline 
    SVectorWithError<N,SCALAR>
    operator+ 
    ( const ROOT::Math::SVector<SCALAR,N>& v2 ,
      const SVectorWithError<N,SCALAR>&    v1 ) { return v1 + v2 ; }
    // ========================================================================
    template <unsigned int N, class SCALAR>
    inline 
    SVectorWithError<N,SCALAR>
    operator- 
    ( const SVectorWithError<N,SCALAR>&    v1 , 
      const ROOT::Math::SVector<SCALAR,N>& v2 )  
    {
      SVectorWithError<N,SCALAR> tmp ( v1 ) ;
      return tmp -= v2 ;
    }
    // ========================================================================
    template <unsigned int N, class SCALAR>
    inline 
    SVectorWithError<N,SCALAR>
    operator- 
    ( const ROOT::Math::SVector<SCALAR,N>& v2 ,
      const SVectorWithError<N,SCALAR>&    v1 ) 
    {
      return SVectorWithError<N,SCALAR> ( v2 - v1.value() , v1.cov2()  ) ;
    }
    // ========================================================================
    template <unsigned int N, class SCALAR, class B>
    inline 
    SVectorWithError<N,SCALAR>
    operator+ 
    ( const SVectorWithError<N,SCALAR>&      v1 , 
      const ROOT::Math::VecExpr<B,SCALAR,N>& v2 ) 
    {
      SVectorWithError<N,SCALAR> tmp ( v1 ) ;
      return tmp += v2 ;
    }
    // ========================================================================
    template <unsigned int N, class SCALAR, class B>
    inline 
    SVectorWithError<N,SCALAR>
    operator- 
    ( const SVectorWithError<N,SCALAR>&      v1 , 
      const ROOT::Math::VecExpr<B,SCALAR,N>& v2 ) 
    {
      SVectorWithError<N,SCALAR> tmp ( v1 ) ;
      return tmp -= v2 ;
    }
    // ========================================================================
    template <unsigned int N, class SCALAR, class B>
    inline 
    SVectorWithError<N,SCALAR>
    operator+ 
    ( const ROOT::Math::VecExpr<B,SCALAR,N>& v2 , 
      const SVectorWithError<N,SCALAR>&      v1 ) { return v1 + v2 ; }
    // ========================================================================
    template <unsigned int N, class SCALAR, class B>
    inline 
    SVectorWithError<N,SCALAR>
    operator- 
    ( const ROOT::Math::VecExpr<B,SCALAR,N>& v2 , 
      const SVectorWithError<N,SCALAR>&      v1 ) 
    {
      return SVectorWithError<N,SCALAR> ( v2 - v1.value() , v1.cov2()  ) ;
    }
    // ========================================================================
    template <unsigned int N, class SCALAR>
    inline 
    SVectorWithError<N,SCALAR>
    operator* 
    ( const SVectorWithError<N,SCALAR>&      v1 , 
      const SCALAR                           v2 ) 
    {
      SVectorWithError<N,SCALAR> tmp ( v1 ) ;
      return tmp *= v2 ;
    }
    // ========================================================================
    template <unsigned int N, class SCALAR>
    inline 
    SVectorWithError<N,SCALAR>
    operator/
    ( const SVectorWithError<N,SCALAR>&      v1 , 
      const SCALAR                           v2 ) 
    {
      SVectorWithError<N,SCALAR> tmp ( v1 ) ;
      return tmp /= v2 ;
    }
    // ========================================================================
    template <unsigned int N, class SCALAR>
    inline 
    SVectorWithError<N,SCALAR>
    operator* 
    ( const SCALAR                           v2 , 
      const SVectorWithError<N,SCALAR>&      v1 )
    {
      SVectorWithError<N,SCALAR> tmp ( v1 ) ;
      return tmp *= v2 ;
    }
    // ========================================================================
    /** "transform" vector with uncertainties using matrix 
     *   @param M (INPUT)
     *   @return transformed vector
     */
    template <unsigned int K, unsigned int N, typename SCALAR> 
    SVectorWithError<K,SCALAR>
    operator*
    ( const  ROOT::Math::SMatrix<SCALAR,K,N,ROOT::Math::MatRepStd<SCALAR,K,N> >& M , 
      const  SVectorWithError<N,SCALAR>&        v ) 
    { return SVectorWithError<K,SCALAR> ( M * v.value() , v.cov2().Similarity ( M ) ) ; }
    // ========================================================================
    /** "transform" vector with uncertainties using matrix 
     *   @param M (INPUT)
     *   @return transformed vector
     */
    template <unsigned int N, typename SCALAR> 
    SVectorWithError<N,SCALAR>
    operator*
    ( const  ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& M , 
      const  SVectorWithError<N,SCALAR>&        v ) 
    { return SVectorWithError<N,SCALAR> ( M * v.value() , v.cov2().Similarity ( M ) ) ; }
    // ========================================================================
    
    // ========================================================================
    /// evaluate the mean of a and b 
    template <unsigned int N, class SCALAR>
    inline 
    SVectorWithError<N,SCALAR> mean 
    ( const SVectorWithError<N,SCALAR>& v1 , 
      const SVectorWithError<N,SCALAR>& v2 ) { return v1.mean ( v2 ) ; }
    /// evaluate the mean of a and b 
    template <unsigned int N, class SCALAR>
    inline 
    SVectorWithError<N,SCALAR> average 
    ( const SVectorWithError<N,SCALAR>& v1 , 
      const SVectorWithError<N,SCALAR>& v2 ) { return v1.mean ( v2 ) ; }
    // ========================================================================
    /** Get symmetrized Kullback-Leibler divergency for two objects 
     *  @see https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence
     *  @see Ostap::Math::kullback_leibler 
     */
    template <unsigned int N, class SCALAR>
    inline double
    kullback_leibler 
    ( const SVectorWithError<N,SCALAR>& v1 , 
      const SVectorWithError<N,SCALAR>& v2 ) 
    { return v1.kullback_leibler  ( v2 ) ; }
    // ========================================================================
    /** Get asymmetric Kullback-Leibler divergency for two objects 
     *  @see https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence
     *  @see Ostap::Math::asymmetric_kullback_leibler 
     */
    template <unsigned int N, class SCALAR>
    inline double 
    asymmetric_kullback_leibler 
    ( const SVectorWithError<N,SCALAR>& v1 , 
      const SVectorWithError<N,SCALAR>& v2 ) 
    { return v1.asymmetric_kullback_leibler  ( v2 ) ; }
    // ========================================================================
    /**  get Mahalanobis distance 
     *  https://en.wikipedia.org/wiki/Mahalanobis_distance  
     *  @return Mahalanobis distance (-1 in case of error)
     */
    template <unsigned int N, class SCALAR>
    inline double
    mahalanobis  
    ( const Ostap::Math::SVectorWithError<N,SCALAR>& a , 
      const Ostap::Math::SVectorWithError<N,SCALAR>& b )
    { return a.mahalanobis ( b ) ; }
    // ========================================================================    
    /**  get Mahalanobis distance 
     *  https://en.wikipedia.org/wiki/Mahalanobis_distance  
     *  @return Mahalanobis distance (-1 in case of error)
     */
    template <unsigned int N, class SCALAR>
    inline double
    mahalanobis  
    ( const Ostap::Math::SVectorWithError<N,SCALAR>& a ,
      const ROOT::Math::SVector<SCALAR,N>&           b ) 
    { return a.mahalanobis ( b ) ; }
    // ========================================================================    
    /**  get Mahalanobis distance 
     *  https://en.wikipedia.org/wiki/Mahalanobis_distance  
     *  @return Mahalanobis distance (-1 in case of error)
     */
    template <unsigned int N, class SCALAR>
    inline double
    mahalanobis  
    ( const ROOT::Math::SVector<SCALAR,N>&           a ,  
      const Ostap::Math::SVectorWithError<N,SCALAR>& b )
    { return b.mahalanobis ( a ) ; }
    // ========================================================================
    /** get Cholesky decomposition for the covarance matrix 
     *  @param v (INPUT)  vector with uncrtainties 
     *  @param L (UPDATE) Cholesky decomposition of the covariance matrix 
     *  @return true if decomposition OK (matrix is positive definite) else false 
     */
    template <unsigned int N, class SCALAR>
    inline bool 
    cholesky 
    ( const Ostap::Math::SVectorWithError<N,SCALAR>&                      v , 
      ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepStd<SCALAR,N,N> >& L ) 
    { return cholesky ( v.cov2() , L ) ; }
    // ========================================================================
    /** "transform" vector with uncertainties using matrix 
     *   @param M (INPUT)
     *   @return transformed vector
     */
    template <unsigned int K, unsigned int N, typename SCALAR> 
    SVectorWithError<K,SCALAR>
    transform 
    ( const  ROOT::Math::SMatrix<SCALAR,K,N,ROOT::Math::MatRepStd<SCALAR,K,N> > & M , 
      const  Ostap::Math::SVectorWithError<N,SCALAR>&                             v ) 
    {
      return SVectorWithError<K,SCALAR> ( M * v.value() , v.cov2().Similarity ( M ) ) ;
    }
    // ========================================================================
    /** "transform" vector with uncertainties using matrix 
     *   @param M (INPUT)
     *   @return transformed vector
     */
    template <unsigned int N, typename SCALAR> 
    SVectorWithError<N,SCALAR>
    transform 
    ( const  ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> > & M , 
      const  Ostap::Math::SVectorWithError<N,SCALAR>&                           v ) 
    {
      return SVectorWithError<N,SCALAR> ( M * v.value() , v.cov2().Similarity ( M ) ) ;
    }
    // ========================================================================
    /** make unnormalized weighted sum, accounting the uncertainties  
     *  \f[ r  = \sum_i v_i w_i \f] 
     *  @param values  INPUT vector of values 
     *  @param weights INPUT vector of weigths 
     *  @return unnormalized weighted sum 
     */
    template <unsigned int N, typename SCALAR> 
    inline ValueWithError 
    dot 
    ( const SVectorWithError<N,SCALAR>& values  , 
      const SVectorWithError<N,SCALAR>& weights ) ;  
    // ========================================================================
    /** make unnormalized weighted sum, accounting the uncertainties  
     *  \f[ r  = \sum_i v_i w_i \f] 
     *  @param values  INPUT vector of values 
     *  @param weights INPUT vector of weigths 
     *  @return unnormalized weighted sum 
     */
    template <unsigned int N, typename SCALAR> 
    inline ValueWithError 
    dot 
    ( const SVectorWithError<N,SCALAR>&    values  , 
      const ROOT::Math::SVector<SCALAR,N>& weights ) ;
    // ========================================================================
    /** make unnormalized weighted sum, accounting the uncertainties  
     *  \f[ r  = \sum_i v_i w_i \f] 
     *  @param values  INPUT vector of values 
     *  @param weights INPUT vector of weigths 
     *  @return unnormalized weighted sum 
     */
    template <unsigned int N, typename SCALAR> 
    inline ValueWithError 
    dot 
    ( const ROOT::Math::SVector<SCALAR,N>& values  , 
      const SVectorWithError<N,SCALAR>&    weights ) ;
    // ========================================================================
    /** make normalized weighted sum, accounting the uncertainties  
     *  \f[ r  = \frac{ \sum_i v_i w_i } { \sum_i w_i }\f] 
     *  @param values  INPUT vector of values 
     *  @param weights INPUT vector of weigths 
     *  @return normalized weighted sum 
     */
    template <unsigned int N, typename SCALAR> 
    inline ValueWithError 
    weighted_sum
    ( const SVectorWithError<N,SCALAR>& values  , 
      const SVectorWithError<N,SCALAR>& weights ) ;
    // ========================================================================
    /** make normalized weighted sum, accounting the uncertainties  
     *  \f[ r  = \frac{ \sum_i v_i w_i } { \sum_i w_i }\f] 
     *  @param values  INPUT vector of values 
     *  @param weights INPUT vector of weigths 
     *  @return normalized weighted sum 
     */
    template <unsigned int N, typename SCALAR> 
    inline ValueWithError 
    weighted_sum
    ( const SVectorWithError<N,SCALAR>&    values  , 
      const ROOT::Math::SVector<SCALAR,N>& weights ) ;
    // ========================================================================
    /** make normalized weighted sum, accounting the uncertainties  
     *  \f[ r  = \frac{ \sum_i v_i w_i } { \sum_i w_i }\f] 
     *  @param values  INPUT vector of values 
     *  @param weights INPUT vector of weigths 
     *  @return normalized weighted sum 
     */
    template <unsigned int N, typename SCALAR> 
    inline ValueWithError 
    weighted_sum
    ( const ROOT::Math::SVector<SCALAR,N>& values  , 
      const SVectorWithError<N,SCALAR>&    weights ) ;
    // ========================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap 
// ============================================================================
#include "SVectorWithError.icpp"
// ============================================================================
namespace Ostap 
{
  // ==========================================================================
  namespace Math 
  {
    // =======================================================================
    /** get Hotelling's t-squared statistics 
     *  @see https://en.wikipedia.org/wiki/Hotelling%27s_T-squared_distribution#Two-sample_statistic
     *  \f[ t^2 = \frac{n_1 n_2}{n_1+n_2} \left(v_1-v_2\right)^T \Sigma^{-1} \left( v1-v2) \sim
     *   T^2 ( p , n_1 + n_2 -2 \f] 
     *  @author Vanya BELYUAEV Ivan.Belyaev@itep.ru
     *  @date 2023-03-07
     */
    template <unsigned int N, typename SCALAR>
    inline double hotelling 
    ( const Ostap::Math::SVectorWithError<N,SCALAR>& x  , 
      const unsigned long                            nx ,  
      const Ostap::Math::SVectorWithError<N,SCALAR>& y  , 
      const unsigned long                            ny )
    {
      return Ostap::Math::hotelling 
        ( x.value () , x.cov2 () , nx , 
          y.value () , y.cov2 () , ny ) ;      
    }
    // ========================================================================
    /// Are all elements are finite? 
    template <unsigned int N, typename SCALAR> 
    inline bool isfinite 
    ( const SVectorWithError<N,SCALAR>& vct )
    {
      return
	Ostap::Math::isfinite ( vct.value () ) &&
	Ostap::Math::isfinite ( vct.cov2  () ) ;	
    }
    // ========================================================================
  }
  // ==========================================================================
} 
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_SVECTORWITHERROR_H
// ============================================================================
