// ============================================================================
#ifndef OSTAP_COMPARISONS_H 
#define OSTAP_COMPARISONS_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL 
// ============================================================================
#include <type_traits> 
#include <functional> 
#include <algorithm> 
#include <cmath>
#include <array>
#include <vector>
// ============================================================================
// Ostap 
// ============================================================================
#include "Ostap/Param_t.h"
#include "Ostap/mULPS.h"
#include "Ostap/Lomont.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================    
    /** @struct abs_less 
     *  comparison by absolute value 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2007-08-17
     */
    template <class TYPE>
    struct abs_less 
    {
      // ======================================================================
      /// abs ( v1 ) < abs ( v2 ) ?
      constexpr TYPE operator() 
        ( param_t<TYPE> v1 ,
          param_t<TYPE> v2 ) const 
      { return m_eval ( std::abs ( v1 ) , std::abs ( v2 ) ) ; }
      // ======================================================================
      using CleanT = std::decay_t<TYPE> ; 
      /// evaluator: 
      std::less<CleanT> m_eval ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @struct abs_greater
     *  comparison by absolute value 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2007-08-17
     */
    template <class TYPE>
    struct abs_greater
    {
      // ======================================================================
      /// abs ( v1 ) > abs ( v2 ) ?
      constexpr TYPE operator()
      ( param_t<TYPE> v1 ,
        param_t<TYPE> v2 ) const 
      { return m_eval ( std::abs ( v2 ) , std::abs ( v1 ) ) ; }
      // ======================================================================
      /// actual evaluator: 
      using CleanT = std::decay_t<TYPE> ; 
      abs_less<CleanT> m_eval ;
      // ======================================================================
    } ;

    // ========================================================================
    /** compare two double numbers with relative precision 'epsilon'
     *
     *  Essentially it is a wrapper to gsl_fcmp function from GSL library
     *  See D.E.Knuth, "Seminumerical Algorithms", section 4.2.2
     *
     *  @param value1  (INPUT) the first value 
     *  @param value2  (INPUT) the second value 
     *  @param epsilon (INPUT) the (relative) precision 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2007-11-27
     */
    bool knuth_equal_to_double
    ( const double value1           ,
      const double value2           ,
      const double epsilon = 1.0e-8 ) ;
    // ========================================================================
    /** compare two double numbers with precision 'mULPS'
     *  @param value1 (INPUT) the first value 
     *  @param value2 (INPUT) the second value 
     *  @param mulps  (INPUT) the precision 
     *  @see Ostap::Math::lomont_compare_double 
     *  @see Ostap::Math::mULPS_double 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2007-11-27
     */
    inline bool equal_to_double
    ( const double       value1                 ,
      const double       value2                 ,
      const unsigned int mulps  = mULPS<double> ) 
    { return Ostap::Math::Lomont::compare_double ( value1 , value2 , mulps ) ; }
    // ========================================================================
    
    // ========================================================================
    /** @struct Equal_To
     *  helper structure for comparison of floating values
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2007-11-27
     */
    template <class TYPE>
    struct Equal_To 
    {
      // ======================================================================
      /// (fake) constructor
      constexpr Equal_To ( const unsigned int eps = mULPS<TYPE> ) {} 
      // ======================================================================        
      /// comparison
      constexpr bool operator()
      ( param_t<TYPE> v1 ,
        param_t<TYPE> v2 ) const
      { return m_cmp ( v1 , v2 ) ; }
      // ======================================================================
    private :      
      // ======================================================================
      using CleanT = std::decay_t<TYPE> ;
      /// comparison criteria 
      std::equal_to<CleanT> m_cmp {} ; // comparison criteria 
      // ======================================================================
    } ;
    // ========================================================================
    /// partial specialization for const-types
    template <class TYPE>
    struct Equal_To<const TYPE>: public Equal_To<TYPE> {} ;
    // ========================================================================
    /// partial specialization for references
    template <class TYPE>
    struct Equal_To<TYPE&>     : public Equal_To<TYPE> {} ;
    // ========================================================================
    /** explicit specialization for doubles
     *  @see Ostap::Math::mULPS_double 
     */
    template <>
    struct Equal_To<double>
    {
    public:
      // ======================================================================
      /// constructor
      constexpr Equal_To
      ( const unsigned int eps = mULPS<double> ) : m_cmp ( eps ) {}
      /// comparison:
      inline bool operator() ( const double v1 , const double v2 ) const
      { return m_cmp ( v1 , v2 ) ; }
      // ======================================================================
    private :
      // ======================================================================
      Equal_To ( const double /* eps */ ) ;
      // ======================================================================
    private :
      // ======================================================================
      /// evaluator 
      Lomont_<double> m_cmp ;                      // the evalautor 
      // ======================================================================
    };
    // ========================================================================
    /** explicit specialization for long doubles
     *  @see Ostap::Math::mULPS_double 
     */
    template <>
    struct Equal_To<long double>
    {
    public:
      // ======================================================================
      /// constructor
      constexpr Equal_To
      ( const unsigned int eps = mULPS<long double> )
        : m_cmp ( eps ) {}
      /// comparison:
      inline bool operator() 
      ( const long double v1 ,
        const long double v2 ) const
      { return m_cmp ( static_cast<double> ( v1 ) ,
                       static_cast<double> ( v2 ) ) ; }
      // ======================================================================
    private :
      // ======================================================================
      /// constructor
      Equal_To ( const long double /* eps */ ) ;
      // ======================================================================
    private :
      // ======================================================================
      /// the evaluator 
      Equal_To<double> m_cmp ;                                 // the evaluator 
      // ======================================================================
    };
    // ========================================================================
    /** explicit specialization for floats
     *  @see Ostap::Math::mULPS_float
     *  @see Ostap::Math::Lomont_
     *  @see Ostap::Math::Lomont_<float>
     */
    template <>
    struct Equal_To<float>
    {
    public:
      // ======================================================================
      /** constructor
       *  @see Ostap::Math::mULPS_float
       */
      constexpr Equal_To
      ( const unsigned short eps =  mULPS<float> )
      : m_cmp ( eps ) {}
      /// comparison:
      inline bool operator() ( const float v1 , const float v2 ) const
      { return m_cmp( v1 , v2 ) ; }
      // ======================================================================
    private :
      // ======================================================================
      /// the evaluator 
      Lomont_<float> m_cmp ;                       // the evaluator
      // ======================================================================
    } ;

    // ========================================================================
    /// specialization for complex values 
    template <class TYPE>
    struct Equal_To< std::complex<TYPE> > 
    {
    public: 
      // ======================================================================
      /// constructor
      constexpr Equal_To
        ( const unsigned int eps = mULPS<TYPE> )
        : m_equal ( eps )
        {}
      // ======================================================================
    public:
      // ======================================================================
      /// comparison:
      inline bool operator() 
      ( const std::complex<TYPE>& v1 ,
        const std::complex<TYPE>& v2 ) const
      { return
          m_equal ( v1.real () , v2.real () ) && 
          m_equal ( v1.imag () , v2.imag () ) ; }
      /// comparison:
      template <class TYPE2>
      inline bool operator() 
      ( const std::complex<TYPE> & v1 ,
        const std::complex<TYPE2>& v2 ) const
      { return
          m_equal ( v1.real () , v2.real () ) && 
          m_equal ( v1.imag () , v2.imag () ) ; }
      /// comparison:
      template <class TYPE2>
      inline bool operator() 
      ( const std::complex<TYPE2>& v1 ,
        const std::complex<TYPE> & v2 ) const
      { return
          m_equal ( v1.real () , v2.real () ) && 
          m_equal ( v1.imag () , v2.imag () ) ; }
      // =======================================================================
    private: 
      // =======================================================================
      /// comparison for real and imaginary parts 
      Equal_To<TYPE> m_equal ;
      // =======================================================================
    };
    // =========================================================================    
    namespace detail
    {
      /// Generic helper to compare any two sequence containers element-wise.
      template <class C1, class C2, class Comp>
      constexpr bool compare_containers
      ( const C1& c1 ,
        const C2& c2 , Comp&& cmp)
        {
          if constexpr ( std::is_same_v<std::decay_t<C1> , std::decay_t<C2>> )
          { if (&c1 == &c2) { return true; } }
          //
          return c1.size() == c2.size() && 
            std::equal ( c1.begin(), c1.end(), c2.begin(), std::forward<Comp> ( cmp ) ) ;
        }
      // ======================================================================
    } // namespace detail
    // ========================================================================    
    /** specialisation for vectors 
     *  @see Ostap::Math::mULPS_double
     *  @see Ostap::Math::Lomont_
     *  @see Ostap::Math::Lomont_<double>
     */
    template <class T, class Alloc>
    struct Equal_To<std::vector<T, Alloc> >
    {
    public:
      // ======================================================================
      using value_type = T;
      // ======================================================================
      /**
       * @brief Constructor with optional ULP threshold.
       * @see Ostap::Math::mULPS
       */
      constexpr explicit Equal_To
      ( const unsigned int eps = mULPS<T> ) 
        : m_cmp(eps) {}
      // ======================================================================
      /**
       * @brief Comparison operator for two vectors
       * @tparam Container1 First container type.
       * @tparam Container2 Second container type.
       * @param v1 First container.
       * @param v2 Second container.
       * @return True if vectors are equal within the given ULP threshold.
       */
      template <class Container1, class Container2>
      bool operator()
      ( const Container1& v1 ,
        const Container2& v2 ) const
      { return detail::compare_containers ( v1 , v2 , m_cmp ) ; }
      // ======================================================================
    private:
      // ======================================================================
      using CleanT = std::decay_t<T> ;
      Equal_To<CleanT> m_cmp; ///< Evaluator for individual vector elements
      // ======================================================================
    } ;
    // ========================================================================
    /// specializatinofor arrays 
    template <class T, std::size_t N>
    struct Equal_To<std::array<T, N>>
    {
      public:
      // =====================================================================
      constexpr explicit Equal_To
        ( const unsigned int eps = mULPS<T> ) 
        : m_cmp(eps) {}
      // =====================================================================      
      template <class Array1, class Array2>
      bool operator()
      ( param_t<Array1> a1 ,
        param_t<Array2> a2 ) const
      { return detail::compare_containers ( a1 , a2 , m_cmp ); }
      // =====================================================================
    private:
      // =====================================================================
      Equal_To<std::decay_t<T>> m_cmp;
      // =====================================================================
    };
    
    // ========================================================================
    template <class TYPE> struct    Zero ;
    template <class TYPE> struct NonZero ;
    // ========================================================================    
    
    /** @struct Zero
     *  helper structure for comparison of floating values
     *  @author Vanya BELYAEV Ivan.Belyaev@iep.ru
     *  @date 2007-11-27
     */
    template <class TYPE>
    struct Zero
    {
      // ======================================================================
      constexpr Zero ( const unsigned int eps  = mULPS<TYPE> )
        : m_cmp ( eps ) {}
      /// comparison
      inline bool operator() ( param_t<TYPE> v ) const
      { return !v || m_cmp ( v , 0 ) ; }
      // ======================================================================
    private:
      // ======================================================================
      // the comparizon criteria
      using CleanT = std::decay_t<TYPE> ;
      Equal_To<CleanT> m_cmp {} ;
      // ======================================================================
    } ;

     /// partial specialisation for the complex values 
    template <class TYPE>
    struct Zero< std::complex<TYPE> >
    {
      // ======================================================================
      constexpr Zero ( const unsigned int eps  = mULPS<TYPE> )
        : m_zero ( eps ) {}
      // ======================================================================
      /// comparison
      inline bool operator() ( const std::complex<TYPE>& v ) const 
      { return m_zero ( v.real() ) && m_zero ( v.imag () ) ; }
      // ======================================================================
    private:
      // ======================================================================
      // the comparison criteria 
      using CleanT = std::decay_t<TYPE> ;
      Zero<CleanT> m_zero {} ;
      // ======================================================================
    } ;
    
    // ========================================================================
    /** @struct NotZero
     *  helper structure for comparison of floating values
     *  @author Vanya BELYAEV Ivan.Belyaev@iep.ru
     *  @date 2007-11-27
     */
    template <class TYPE>
    struct NotZero 
    {
      // ======================================================================
      constexpr NotZero ( const unsigned int eps  = mULPS<TYPE> )
        : m_zero ( eps ) {}
      // ======================================================================      
      /// comparison
      inline bool operator() ( param_t<TYPE> v ) const
      { return !m_zero ( v ) ; }
      // ======================================================================
    private:
      // ======================================================================
      // the comparison criteria 
      using CleanT = std::decay_t<TYPE> ;
      Zero<CleanT> m_zero ;
      // ======================================================================
    } ;
    
    // ========================================================================
    /** specialisation for vectors 
     *  @see Ostap::Math::Zero
     *  @see Ostap::Math::Equal_To
     *  @see Ostap::Math::Lomont<float>
     */
    template < class TYPE, class ALLOCATOR>
    struct Zero< std::vector<TYPE,ALLOCATOR> > 
    {
    public:
      // ======================================================================
      constexpr Zero ( const unsigned int eps = mULPS<TYPE> )
        : m_zero ( eps ) {}
      // ======================================================================      
      ///  comparison
      inline bool operator () ( const std::vector<TYPE,ALLOCATOR>& v ) const
      {
        /// empty vector or all elements are zeros 
        return v.empty() || std::all_of ( v.begin() , v.end  () , m_zero )  ;
      }
      // ======================================================================
    private :
      // ======================================================================
      /// comparison criterion for elements
      using CleanT = std::decay_t<TYPE> ;
      Zero<CleanT> m_zero ;
      // ======================================================================
    } ;


    // ========================================================================
    /** specialisation for vectors 
     *  @see Ostap::Math::Zero
     *  @see Ostap::Math::Equal_To
     *  @see Ostap::Math::Lomont<float>
     */
    template < class TYPE, std::size_t N>
    struct Zero< std::array<TYPE,N> > 
    {
    public:
      // ======================================================================
      constexpr Zero ( const unsigned int eps = mULPS<TYPE> )
        : m_zero ( eps ) {}
      // ======================================================================      
      ///  comparison
      inline bool operator () ( const std::array<TYPE,N>& v ) const
      {
        /// all elements are zeros 
        return std::all_of ( v.begin() , v.end  () , m_zero )  ;
      }
      // ======================================================================
    private :
      // ======================================================================
      /// comparison criterion for elements
      using CleanT = std::decay_t<TYPE> ;
      Zero<CleanT> m_zero ;
      // ======================================================================
    } ;

    // ========================================================================
    /// Is value sufficiently  small ?
    template <class TYPE>
    struct Small 
    {
      // ======================================================================
      /// constructor with threshold 
      constexpr Small ( param_t<TYPE> a ) : m_a ( std::abs ( a ) ) {}
      /// the only one important method   
      constexpr  bool operator() ( param_t<TYPE> a ) const
      { return std::abs ( a ) <= m_a ; }
      // ======================================================================
    private: 
      // ======================================================================
      Small () = delete ;  // default constructor is disabled 
      // ======================================================================
    private :
      // ======================================================================
      using CleanT = std::decay_t<TYPE> ;
      CleanT m_a ;
      // ======================================================================
    } ;
   
    // ========================================================================
    /// specialization for complex values 
    template <class TYPE>
    struct Small< std::complex <TYPE> >  
    {
      // ======================================================================
      // constructor with threshold 
      constexpr Small ( param_t<TYPE> a ) : m_small ( a ) {}
      // the opnly one important method   
      constexpr bool operator () ( const std::complex<TYPE>& a ) const
      { return m_small ( a.real() ) && m_small ( a.imag() ); }
      // ======================================================================
    private :
      // ======================================================================
      /// default constructor is disabled 
      Small () = delete  ;  // default constructor is disabled 
      // ======================================================================
    private :
      // ======================================================================
      using CleanT = std::decay_t<TYPE> ;
      Small<CleanT> m_small  ;
      // ======================================================================
    } ;

    // ========================================================================
    /** specialization for vectors 
     *  vector is small, if empty or all elements are small 
     */ 
    template <class TYPE, class ALLOCATOR>
    struct Small< std::vector<TYPE,ALLOCATOR> > 
    {
      // ======================================================================
      /// constructor with threshold 
      constexpr Small ( param_t<TYPE> a ) : m_cmp ( a ) {}
      /// the only one important method   
      inline bool operator() ( const std::vector<TYPE,ALLOCATOR>& v ) const
      { return v.empty() || std::all_of ( v.begin() , v.end() , m_cmp) ; }
      // ======================================================================
    private :
      // ======================================================================
      /// default constructor is disabled 
      Small () ;  // default constructor is disabled 
      // ======================================================================
    private :
      // ======================================================================
      /// comparison 
      using CleanT = std::decay_t<TYPE> ;
      Small<CleanT> m_cmp ;
      // ======================================================================
    } ;

    // ========================================================================
    /** specialization for arrays 
     *  array all elements are small 
     */ 
    template <class TYPE, std::size_t N>
    struct Small< std::array<TYPE,N> > 
    {
      // ======================================================================
      /// constructor with threshold 
      constexpr Small ( param_t<TYPE> a ) : m_cmp ( a ) {}
      /// the only one important method   
      constexpr bool operator() ( const std::array<TYPE,N>& v ) const
      { return std::all_of ( v.begin() , v.end() , m_cmp) ; }
      // ======================================================================
    private :
      // ======================================================================
      /// default constructor is disabled 
      Small () ;  // default constructor is disabled 
      // ======================================================================
    private :
      // ======================================================================
      /// comparison 
      using CleanT = std::decay_t<TYPE> ;
      Small<CleanT> m_cmp ;
      // ======================================================================
    } ;

    // ========================================================================
    /** @struct MuchSmaller 
     *  Is a value of "a" tiny with respect to b ? 
     *  - if b is numerical zero, a is numerical zero also 
     *  - otherwise (a+b) is numerically equal to b 
     */
    template <class TYPE>
    struct MuchSmaller 
    {
    public:
      // ======================================================================
      constexpr MuchSmaller ( const unsigned int eps = mULPS<TYPE> )
        : m_zero  ( eps )
        , m_equal ( eps )
      {}
      // ======================================================================      
      /** Is a value of "a" tiny with respect to b ? 
       *  - if b is numerically zero, a is also zero  
       *  - otherwise (a+b) is numerically equal to b 
       */
      constexpr bool operator ()
        ( param_t<TYPE> a ,
          param_t<TYPE> b )  const 
      { return m_zero ( b ) ? m_zero ( a ) : m_equal ( a + b , b ) ; }
      // ======================================================================      
    private :
      // ======================================================================
      using CleanT = std::decay_t<TYPE> ;      
      /// zero ?
      Zero    <CleanT> m_zero  {} ; // zero ?
      /// eqiality ?
      Equal_To<CleanT> m_equal {} ; // equality ? 
      // ======================================================================
    } ;
    
    // ========================================================================
    /** @struct Tiny 
     *  Is a value of "a" tiny with respect to b ? 
     *  -  if b is numerical zero, a is numerical zero also 
     *  - otherwise (a+b) is numerically equal to b 
     */
    template <class TYPE>
    struct Tiny
    {
    public:
      // ======================================================================
      // constructor 
      constexpr Tiny
      ( param_t<TYPE>      b  ,
        const unsigned int eps = mULPS<TYPE> )
        : m_b   ( b   ) 
        , m_cmp ( eps )
      {}
      // ======================================================================
    public:
      // ======================================================================      
      /// Is a value of "a" tiny with respect to b ? 
      constexpr bool operator () ( param_t<TYPE> a )  const
      { return m_cmp ( a , m_b ) ; }
      // ======================================================================      
    private :
      // ======================================================================
      /// default constructor is disabled 
      Tiny () ; // default constructor is disabled 
      // ======================================================================      
    private :
      // ======================================================================
      /// the reference value
      using CleanT = std::decay_t<TYPE> ;            
      CleanT              m_b   ;
      /// smaller ? 
      MuchSmaller<CleanT> m_cmp ;
      // ======================================================================
    } ;
    
    // ========================================================================
    /** @struct LessOrEqual
     *  check if two values ar less or equal (numerically)
     *  @see Ostap::Math::Equal_To
     */
    template <class TYPE>
    struct LessOrEqual 
    {
      // ======================================================================
      constexpr LessOrEqual ( unsigned int eps = mULPS<TYPE> )
        : m_equal ( eps )
      {}
      // ====================================================================== 
      /// the only one method:  \f$ o_1 \le  o_2\f$  or \f$ o_1 \approx o_2\f$ 
      constexpr bool operator ()
        ( param_t<TYPE> o1 ,
          param_t<TYPE> o2 ) const 
      { return m_leq ( o1 , o2 ) || m_equal ( o1 , o2 ) ; }
      // ======================================================================
    private:
      // ======================================================================
      using CleanT = std::decay_t<TYPE> ;                  
      std::less_equal<CleanT> m_leq   ; // ordering criteria 
      Equal_To<CleanT>        m_equal ; // equality criteria 
      // ======================================================================
    } ;  

    // ========================================================================
    /** @struct GreaterOrEqual
     *  check if two values are greater or equal (numerically)
     *  @see Ostap::Math::Equal_To
     */
    template <class TYPE>
    struct GreaterOrEqual 
    {
      // ====================================================================== 
      constexpr GreaterOrEqual ( unsigned int eps = mULPS<TYPE> )
        : m_equal ( eps )
      {}
      // ====================================================================== 
      /// the only one method:  \f$ o_1 \ge o_2 \f$  or  \f$ o_1 \approx = o_2 \f$  
      constexpr  bool operator ()
        ( param_t<TYPE> o1 ,
          param_t<TYPE> o2 ) const 
      { return m_geq ( o1 , o2 ) || m_equal ( o1 , o2 ) ; }
      // ======================================================================
    private:
      // ======================================================================
      /// comparison
      using CleanT = std::decay_t<TYPE> ;                  
      std::greater_equal<CleanT>   m_geq   ; // ordering criteria
      /// equality
      Equal_To<CleanT>             m_equal ; // equality criteria 
      // ======================================================================
    } ;  

    // ========================================================================
    /** @struct  NumLess 
     *  "Numerically less"
     *  useful structure for sorting  
     *  @see Ostap::Math::Equal_To
     */
    template <class TYPE>
    struct NumLess 
    {
      // ====================================================================== 
      constexpr NumLess ( unsigned int eps = mULPS<TYPE> )
        : m_equal ( eps )
      {}
      // ======================================================================
      constexpr bool operator () ( const TYPE& o1 , const TYPE& o2 ) const 
      { return m_less ( o1 , o2 ) && !m_equal ( o1 , o2 ) ; }
      // ======================================================================      
    private:
      // ======================================================================
      using CleanT = std::decay_t<TYPE> ;                  
      /// comparion criteria for objects 
      std::less<CleanT>            m_less  ; // comparion criteria for objects
      /// equality criteria for  objects  
      Equal_To<CleanT>             m_equal ; // equality criteria for objects  
      // ======================================================================
    } ;

    // ========================================================================    
    template <class TYPE>
    struct IsReal
    {
    public:
      // ====================================================================== 
      constexpr IsReal ( unsigned int eps = mULPS<TYPE> )
        : m_ms    ( eps )
        {}
      // ======================================================================
      constexpr bool operator () ( const std::complex<TYPE>& z ) const 
      { return m_ms ( z.imag () , z.real() ) ;  }
      // ======================================================================
    private:
      // ======================================================================
      using CleanT = std::decay_t<TYPE> ;
      MuchSmaller<CleanT>   m_ms ;
      // =======================================================================
    } ;
    // ========================================================================
    template <class TYPE>
    struct IsImagine
    {
    public:
      // ====================================================================== 
      constexpr IsImagine ( unsigned int eps = mULPS<TYPE> )
        : m_ms    ( eps )
        {}
      // ======================================================================
      constexpr bool operator () ( const std::complex<TYPE>& z ) const 
      { return m_ms ( z.real () , z.imag () ) ;  }
      // ======================================================================
    private:
      // ======================================================================
      using CleanT = std::decay_t<TYPE> ;
      MuchSmaller<CleanT>   m_ms ;
      // =======================================================================
    } ;
    
    // ========================================================================
    /** Is the complex value actually real ?
     *  - imaginary part is exact zero
     *  - imaginary part is very small
     *  - imaginary part is negligible in comparison to the real part
     *
     *  @see Ostap::Math::Equal_To
     *  @see Ostap::Math::Zero 
     */
    template <class TYPE>
    inline bool isreal
    ( const std::complex<TYPE>& z )
    {
      using CleanT = std::decay<TYPE> ;
      constexpr static IsReal<CleanT> xreal{} ;
      return xreal ( z ) ;      
    }
    // ========================================================================
    /** Is the complex value pure imaginary? ?
     *  - real part is exact zero
     *  - real part is very small
     *  - real part is negligible in comparison to the real part
     *
     *  @see Ostap::Math::Equal_To
     *  @see Ostap::Math::Zero 
     */
    template <class TYPE>
    inline bool isimagine
    ( const std::complex<TYPE>& z )
    {
      using CleanT = std::decay<TYPE> ;
      constexpr static IsImagine<CleanT> ximag{} ;
      return ximag ( z ) ;      
    }
    
    // ========================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_COMPARISONS_H
// ============================================================================
