// $Id$ 
// ============================================================================
#ifndef OSTAP_MATH_H 
#define OSTAP_MATH_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL 
// ============================================================================
#include <cmath>
#include <algorithm>
#include <functional>
#include <vector>
#include <array>
#include <type_traits>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Lomont.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    /// Parameters for numerical calculations (M.Needham)
    static const double hiTolerance    = 1e-40;
    static const double lowTolerance   = 1e-20;
    static const double looseTolerance = 1e-5;
    static const double     sqrt_12 = 3.4641016151377546; // sqrt(12.)
    static const double inv_sqrt_12 = 0.2886751345948129; // 1./sqrt(12.)
    // ========================================================================
    /** @var mULPS_float
     *  "tolerance" parameter for "Lomont"-compare of floating point numbers.
     *  It corresponds to relative ("Knuth/GLS") tolerance of about ~6*10^-6
     *  for values in excess of 10^-37.
     *
     *  @see Gaudi::Math::lomont_compare_float 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2010-01-02
     */
    const unsigned short mULPS_float = 100 ;
    // ========================================================================
    /** @var mULPS_float_low
     *  "Low-tolerance" parameter for "Lomont"-compare of floating point numbers.
     *  It corresponds to relative ("Knuth/GLS") tolerance of about ~6*10^-5
     *  for values in excess of 10^-37.
     *
     *  @see Ostap::Math::lomont_compare_float 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2010-01-02
     */
    const unsigned short mULPS_float_low = 1000 ;
    // =========================================================================
    /** @var mULPS_double
     *  "tolerance" parameter for "Lomont"-compare of floating point numbers.
     *  It corresponds to relative ("Knuth/GLS") tolerance of about ~6*10^-13
     *  for values in excess of 10^-304.
     *  @see Ostap::Math::lomont_compare_double
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2010-01-02
     */
    const unsigned int mULPS_double = 1000 ;
    // ========================================================================
    namespace detail
    {
      // ======================================================================
      template <class T> struct param ;
      //
      template <class T> struct param<const T&> { typedef const T& param_type ; } ;
      template <class T> struct param<const T*> { typedef const T* param_type ; } ;
      template <class T> struct param<T&>: public param<const T&> {} ;
      template <class T> struct param<T*>: public param<const T*> {} ;
      //
      template <class T, bool small> 
      struct _param         { typedef const T& param_type ; } ;
      template <class T> 
      struct _param<T,true> { typedef const T  param_type ; } ;
      //
      template <class T> struct param 
      {
        typedef typename _param<T, 
                                (std::is_arithmetic<T>::value) ||
                                (std::is_pointer<T>   ::value) ||
                                (std::is_enum<T>      ::value) ||
                                (sizeof(T)<=sizeof(void*))>::param_type param_type ;
      };
      // ====================================================================
    }
    // ======================================================================
    /** @struct abs_less 
     *  comparison by absolute value 
     *  @author Vanya BELYAEV ibelyaev@physics.syr.edu
     *  @date 2007-08-17
     */
    template <class TYPE>
    struct abs_less : std::binary_function<TYPE,TYPE,TYPE>
    {
      inline TYPE operator() 
      ( typename detail::param<const TYPE>::param_type v1 ,
        typename detail::param<const TYPE>::param_type v2 ) const 
      { return m_eval ( std::fabs( v1 ) , std::fabs( v2 ) ) ; }
      /// evaluator: 
      std::less<TYPE> m_eval ;
    } ;
    // ========================================================================
    /** @struct abs_greater
     *  comparison by absolute value 
     *  @author Vanya BELYAEV ibelyaev@physics.syr.edu
     *  @date 2007-08-17
     */
    template <class TYPE>
    struct abs_greater : std::binary_function<TYPE,TYPE,TYPE>
    {
      inline TYPE operator() 
      ( typename detail::param<const TYPE>::param_type v1 ,
        typename detail::param<const TYPE>::param_type v2 ) const 
      { return m_eval ( std::fabs( v1 ) , std::fabs( v2 ) ) ; }
      /// evaluator: 
      std::greater<TYPE> m_eval ;
    } ;
    // ========================================================================
    /** return "min_by_abs"
     *  @author Vanya BELYAEV ibelyaev@physics.syr.edu
     *  @date 2007-08-17
     */        
    template <class TYPE> 
    inline TYPE absMin ( TYPE v1 , TYPE v2 ) 
    { return std::min ( std::fabs ( v1 ) , std::fabs ( v2 ) ) ; }
    // ========================================================================
    /** return  "max_by_abs"
     *  @author Vanya BELYAEV ibelyaev@physics.syr.edu
     *  @date 2007-08-17
     */
    template <class TYPE> 
    inline TYPE absMax ( TYPE v1 , TYPE v2 ) 
    { return std::max ( std::fabs ( v1 ) , std::fabs ( v2 ) ) ; }
    // ========================================================================
    /** compare two double numbers with relative precision 'epsilon'
     *
     *  Essentially it is a wrapper to gsl_fcmp function from GSL library
     *  See D.E.Knuth, "Seminumerical Algorithms", section 4.2.2
     *
     *  @param value1  (INPUT) the first value 
     *  @param value2  (INPUT) the second value 
     *  @param epsilon (INPUT) the (relative) precision 
     *  @author Vanya BELYAEV ibelyaev@physics.syr.edu
     *  @date 2007-11-27
     */
    bool knuth_equal_to_double
    ( const double value1           ,
      const double value2           ,
      const double epsilon = 1.0e-6 ) ;
    // ========================================================================
    /** compare two double numbers with precision 'mULPS'
     *  @param value1 (INPUT) the first value 
     *  @param value2 (INPUT) the second value 
     *  @param mULPS  (INPUT) the precision 
     *  @see Ostap::Math::lomont_compare_double 
     *  @see Ostap::Math::mULPS_double 
     *  @author Vanya BELYAEV ibelyaev@physics.syr.edu
     *  @date 2007-11-27
     */
    inline bool equal_to_double
    ( const double value1                      ,
      const double value2                      ,
      const unsigned int mULPS = mULPS_double  ) 
    { return lomont_compare_double ( value1 , value2 , mULPS ) ; }
    // ========================================================================
    /** @struct Equal_To
     *  helper structure for comparison of floating values
     *  @author Vanya BELYAEV ibelyaev@physics.syr.edu
     *  @date 2007-11-27
     */
    template <class TYPE>
    struct Equal_To : public std::binary_function<TYPE,TYPE,bool>
    {
      // ======================================================================
      typedef typename detail::param<const TYPE>::param_type T ;
      // ======================================================================
      /// comparison
      inline bool operator() ( T v1 , T v2 ) const
      {
        std::equal_to<TYPE> cmp ;
        return cmp ( v1 , v2 ) ;
      }
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
      Equal_To ( const unsigned int eps = mULPS_double ) : m_cmp ( eps ) {}
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
      Lomont<double> m_cmp ;                      // the evalautor 
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
      Equal_To ( const unsigned int eps = mULPS_double ) : m_cmp ( eps ) {}
      /// comparison:
      inline bool operator() 
      ( const long double v1 ,
        const long double v2 ) const
      { 
        using namespace std;
#ifdef __INTEL_COMPILER       // Disable ICC remark
#pragma warning(disable:2259) //  non-pointer conversion may lose significant bits
#pragma warning(push)
#endif
        return  m_cmp ( static_cast<double> ( v1 ) , 
                        static_cast<double> ( v2 ) ) ; 
#ifdef __INTEL_COMPILER         // End disable ICC remark
#pragma warning(pop)
#endif
      }
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
     *  @see Ostap::Math::Lomont
     *  @see Ostap::Math::Lomont<float>
     */
    template <>
    struct Equal_To<float>
    {
    public:
      // ======================================================================
      /** constructor
       *  @see Ostap::Math::mULPS_float
       */
      Equal_To ( const unsigned short eps =  mULPS_float ) : m_cmp ( eps ) {}
      /// comparison:
      inline bool operator() ( const float v1 , const float v2 ) const
      { return m_cmp( v1 , v2 ) ; }
      // ======================================================================
    private:
      // ======================================================================      
      /// constructor
      Equal_To ( const float /* eps */ ) ;
      // ======================================================================
    private :
      // ======================================================================
      /// the evaluator 
      Lomont<float> m_cmp ;                       // the evaluator
      // ======================================================================
    } ;
    // ========================================================================
    /** specialisation for vectors 
     *  @see Ostap::Math::mULPS_double
     *  @see Ostap::Math::Lomont
     *  @see Ostap::Math::Lomont<double>
     */
    template <>
    struct Equal_To<std::vector<double> > 
    {
    public:
      // ======================================================================
      /** constructor
       *  @see Ostap::Math::mULPS_double
       */
      Equal_To ( const unsigned int eps  = mULPS_double ) : m_cmp ( eps ) {}
      // ======================================================================
      /// comparison:
      inline bool operator() ( const std::vector<double>& v1 , 
                               const std::vector<double>& v2 ) const
      {
        return ( &v1 == &v2 ) || 
          ( v1.size() == v2.size() &&
            std::equal ( v1.begin () , v1.end () , v2.begin () , m_cmp ) ) ;
      }      
      /// comparison:
      inline bool operator() ( const std::vector<double>& v1 , 
                               const std::vector<float>&  v2 ) const
      {
        return v1.size() == v2.size() && 
          std::equal ( v1.begin () , v1.end () , v2.begin () , m_cmp ) ;
      }      
      /// comparison:
      inline bool operator() ( const std::vector<double>& v1 , 
                               const std::vector<int>&    v2 ) const
      {
        return v1.size() == v2.size() && 
          std::equal ( v1.begin () , v1.end () , v2.begin () , m_cmp ) ;
      }      
      /// comparison:
      inline bool operator() ( const std::vector<double>&       v1 , 
                               const std::vector<unsigned int>& v2 ) const
      {
        return v1.size() == v2.size() && 
          std::equal ( v1.begin () , v1.end () , v2.begin () , m_cmp ) ;
      }      
      /// comparison:
      inline bool operator() ( const std::vector<float>&  v1 , 
                               const std::vector<double>& v2 ) const
      {
        return v1.size() == v2.size() && 
          std::equal ( v1.begin () , v1.end () , v2.begin () , m_cmp ) ;
      }      
      /// comparison:
      inline bool operator() ( const std::vector<int>&    v1 , 
                               const std::vector<double>& v2 ) const
      {
        return v1.size() == v2.size() && 
          std::equal ( v1.begin () , v1.end () , v2.begin () , m_cmp ) ;
      }      
      /// comparison:
      inline bool operator() ( const std::vector<unsigned int>& v1 , 
                               const std::vector<double>&       v2 ) const
      {
        return v1.size() == v2.size() && 
          std::equal ( v1.begin () , v1.end () , v2.begin () , m_cmp ) ;
      }      
      // ======================================================================
    private:
      // ======================================================================
      /// the evaluator 
      Equal_To<double> m_cmp ;                                 // the evaluator 
      // ======================================================================
    } ;
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
    struct Zero : public std::unary_function<TYPE,bool>
    {
      // ======================================================================
      typedef typename detail::param<const TYPE>::param_type T ;
      /// comparison
      inline bool operator() ( T v ) const { return m_cmp ( v , 0 ) ; }
      // ======================================================================
    private:
      // ======================================================================
      // the comparizon criteria 
      Equal_To<TYPE> m_cmp ;
      // ======================================================================
    } ;
    // ========================================================================
    template <>
    struct Zero<double> : public std::unary_function<double,bool>
    {
      // ======================================================================
      /// comparison
      inline bool operator() ( const double  v ) const 
      { return !v || m_cmp ( v , 0 ) ; }
      // ======================================================================
    private:
      // ======================================================================
      // the comparizon criteria 
      Equal_To<double> m_cmp ;
      // ======================================================================
    } ;
    // ========================================================================
    template <>
    struct Zero<float> : public std::unary_function<float,bool>
    {
      // ======================================================================
      /// comparison
      inline bool operator() ( const float  v ) const 
      { return !v || m_cmp ( v , 0 ) ; }
      // ======================================================================
    private:
      // ======================================================================
      // the comparizon criteria 
      Equal_To<float> m_cmp ;
      // ======================================================================
    } ;
    // ========================================================================
    /// partial specialization for const-types
    template <class TYPE>
    struct Zero<const TYPE>    : public Zero<TYPE>     {} ;
    // ========================================================================
    /// partial specialization for references
    template <class TYPE>
    struct Zero<TYPE&>         : public Zero<TYPE>     {} ;
    // ========================================================================
    /** @struct NotZero
     *  helper structure for comparison of floating values
     *  @author Vanya BELYAEV Ivan.Belyaev@iep.ru
     *  @date 2007-11-27
     */
    template <class TYPE>
    struct NotZero : public std::unary_function<TYPE,bool>
    {
      // ======================================================================
      typedef typename detail::param<const TYPE>::param_type T ;
      /// comparison
      inline bool operator() ( T v ) const { return !m_zero ( v ) ; }
      // ======================================================================
    private:
      // ======================================================================
      // the comparison criteria 
      Zero<TYPE> m_zero ;
      // ======================================================================
    } ;
    /// partial specialization for const-types
    template <class TYPE>
    struct NotZero<const TYPE> : public NotZero<TYPE>  {} ;
    // ========================================================================
    /// partial specialization for references
    template <class TYPE>
    struct NotZero<TYPE&>      : public NotZero<TYPE>  {} ;
    // ========================================================================
    /** specialisation for vectors 
     *  @see Ostap::Math::Zero
     *  @see Ostap::Math::Equal_To
     *  @see Ostap::Math::Lomont<float>
     */
    template < class TYPE>
    struct Zero< std::vector<TYPE> > 
    {
    public:
      // ======================================================================
      inline bool operator () ( const std::vector<TYPE>& v ) const
      {
        /// empty vector or all elements are zeros 
        return v.empty() || ( v.end() == std::find_if ( v.begin() , v.end  () , m_nz ) ) ;
      }
      // ======================================================================
    private :
      // ======================================================================
      // comparison criteria for elements 
      NotZero<TYPE> m_nz ;
      // ======================================================================
    } ;
    // ========================================================================
    /// Is value sufficiently  small ?
    template <class TYPE>
    struct Small : public std::unary_function<TYPE,bool>
    {
      // ======================================================================
      typedef  TYPE   Inner ;
      // ======================================================================
      // constructor with threshold 
      Small ( const TYPE& a ) : m_a ( std::abs ( a ) ) {}
      // the opnly one important method   
      inline bool operator() ( const TYPE& a ) const
      { return std::abs ( a ) <= m_a ; }
      // ======================================================================
    private :
      // ======================================================================
      /// default constructor is disabled 
      Small () ;  // default constructor is disabled 
      // ======================================================================
    private :
      // ======================================================================
      TYPE m_a ;
      // ======================================================================
    } ;
    // ========================================================================
    /** specialization for vectors 
     *  vector is small, if empty or all elements are small 
     */ 
    template <class TYPE>
    struct Small<std::vector<TYPE> > : 
      public std::unary_function<std::vector<TYPE>,bool>
    {
      // ======================================================================
      typedef TYPE Inner ;
      // ======================================================================
      // constructor with threshold 
      Small ( const typename Small<TYPE>::Inner & a ) : m_cmp ( a ) {}
      // the only one important method   
      inline bool operator() ( const std::vector<TYPE>& v ) const
      {
        return 
          v.empty() || v.end() == std::find_if 
          ( v.begin() , v.end() , std::not1 ( m_cmp ) ) ;
      }
      // ======================================================================
    private :
      // ======================================================================
      /// default constructor is disabled 
      Small () ;  // default constructor is disabled 
      // ======================================================================
    private :
      // ======================================================================
      Small<TYPE> m_cmp ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @struct LessOrEqual
     *  check if two values ar less or equal (numerically)
     *  @see Ostap::Math::Equal_To
     */
    template <class TYPE>
    struct LessOrEqual : public std::binary_function<TYPE,TYPE,bool>
    {
      // ====================================================================== 
      /// the only one method:  o1 <= o2 or o1 \approx p2 
      inline bool operator () ( const TYPE& o1 , const TYPE& o2 ) const 
      { return m_leq ( o1 , o2 ) || m_equal ( o1 , o2 ) ; }
      // ======================================================================
    private:
      // ======================================================================
      std::less_equal<TYPE> m_leq   ; // ordering criteria 
      Equal_To<TYPE>        m_equal ; // equality criteria 
      // ======================================================================
    } ;  
    // ========================================================================
    /** @struct GreaterOrEqual
     *  check if two values are greater or equal (numerically)
     *  @see Ostap::Math::Equal_To
     */
    template <class TYPE>
    struct GreaterOrEqual : public std::binary_function<TYPE,TYPE,bool>
    {
      // ====================================================================== 
      /// the only one method:  o1 >= o2 or o1 \approx = o2 
      inline bool operator () ( const TYPE& o1 , const TYPE& o2 ) const 
      { return m_geq ( o1 , o2 ) || m_equal ( o1 , o2 ) ; }
      // ======================================================================
    private:
      // ======================================================================
      std::greater_equal<TYPE>   m_geq   ; // ordering criteria
      Equal_To<TYPE>             m_equal ; // equality criteria 
      // ======================================================================
    } ;  
    // ========================================================================
    /** @struct  NumLess 
     *  "Numerically less"
     *  useful structure for sorting  
     *  @see Ostap::Math::Equal_To
     */
    template <class TYPE>
    struct NumLess :  public std::binary_function<TYPE,TYPE,bool>
    {
      // ======================================================================
      inline bool operator () ( const TYPE& o1 , const TYPE& o2 ) const 
      { return m_less ( o1 , o2 ) && !m_equal ( o1 , o2 ) ; }
      // ======================================================================      
    private:
      // ======================================================================
      /// comparion criteria for objects 
      std::less<TYPE>            m_less  ; // comparion criteria for objects
      /// equality criteria for  objects  
      Equal_To<TYPE>             m_equal ; // equality criteria for objects  
      // ======================================================================
    } ;
    // ========================================================================
    /** round to nearest integer, rounds half integers to nearest even integer 
     *  @author Vanya BELYAEV Ivan.BElyaev
     */
    long round ( const double x ) ;
    // ========================================================================
    /** round to nearest integer, rounds half integers to nearest even integer 
     *  @author Vanya BELYAEV Ivan.BElyaev
     */
    inline long round ( const float  x ) { return round ( double ( x ) ) ; }
    // ========================================================================
    /** get mantissa and exponent 
     *  similar to std::frexp, but radix=10)
     *  @param x  INPUT  value 
     *  @param e  UPDATE exponent 
     *  @return  mantissa     (0.1<=m<1)
     *  @author Vanya BELYAEV Ivan.Belyaev       
     *  @date 2015-07-21
     */
    double frexp10 ( const double x , long& e ) ;
    // ========================================================================
    /** get mantissa and exponent 
     *  similar to std::frexp, but radix=10)
     *  @param x  INPUT  value 
     *  @param e  UPDATE exponent 
     *  @return  mantissa    (0.1<=m<1) 
     *  @author Vanya BELYAEV Ivan.Belyaev       
     *  @date 2015-07-21
     */
    float frexp10 ( const float x , long& e ) ;
    // ========================================================================
    /** get mantissa and exponent 
     *  similar to std::frexp, but radix=10)
     *  @param x  INPUT  value 
     *  @return   pair of mantissa (0.1<=m<1) and exponent 
     *  @author Vanya BELYAEV Ivan.Belyaev       
     *  @date 2015-07-21
     */
    std::pair<double,int>
    frexp10 ( const double x ) ;
    // ========================================================================
    /** round to N-significant digits 
     *  @param x  INPUT  input value 
     *  @param n  INPUT  number of significnat digits 
     *  @return rounded value 
     *  @author Vanya BELYAEV Ivan.Belyaev       
     *  @date 2015-07-21
     */
    double round_N ( const double x , const unsigned short n ) ;
    // ========================================================================
    /** round to N-significant digits 
     *  @param x  INPUT  input value 
     *  @param n  INPUT  number of significnat digits 
     *  @return rounded value 
     *  @author Vanya BELYAEV Ivan.Belyaev       
     *  @date 2015-07-21
     */
    float round_N ( const float x , const unsigned short n ) ;
    // ========================================================================
    /** is the value actually long ?
     *  @author Vanya BELYAEV Ivan.Belyaev       
     *  @date 2011-07-18
     */
    bool islong ( const double x ) ;
    // ========================================================================
    /** is the value actually long ?
     *  @author Vanya BELYAEV Ivan.Belyaev       
     *  @date 2011-07-18
     */
    bool islong ( const float  x ) ;
    // ========================================================================
    /** is the value actually int ?
     *  @author Vanya BELYAEV Ivan.Belyaev       
     *  @date 2011-07-18
     */
    bool isint  ( const double x ) ;
    // ========================================================================    
    /** is the value actually int ?
     *  @author Vanya BELYAEV Ivan.Belyaev       
     *  @date 2011-07-18
     */
    bool isint  ( const float  x ) ;
    // ========================================================================    
    /** check if the double value is actually equal to the integer value  
     *  @param val value to be compared with the integer 
     *  @param ref the reference integer number 
     *  @param mULPS the precision 
     *  @see Ostap::Math::lomont_compare_double 
     *  @see Ostap::Math::mULPS_double
     *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
     *  @date 2008-09-17
     */
    bool equal_to_int 
    ( const double       val                  , 
      const int          ref                  , 
      const unsigned int mULPS = mULPS_double ) ;
    // ========================================================================
    /** check if the double value is actually equal to the integer value  
     *  @param ref the reference integer  number 
     *  @param val value to be compared with the integer 
     *  @param mULPS the precision 
     *  @see Ostap::Math::lomont_compare_double 
     *  @see Ostap::Math::mULPS_double
     *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
     *  @date 2008-09-17
     */        
    inline bool equal_to_int 
    ( const int          ref                  , 
      const double       val                  , 
      const unsigned int mULPS = mULPS_double ) 
    { return equal_to_int ( val , ref , mULPS ) ; }
    // ========================================================================
    /** check if the double value is actually equal to the unsigned integer value  
     *  @param val value to be compared with the unsigned integer 
     *  @param ref the reference unsigned integer number 
     *  @param mULPS the precision 
     *  @see Ostap::Math::lomont_compare_double 
     *  @see Ostap::Math::mULPS_double
     *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
     *  @date 2008-09-17
     */        
    bool equal_to_uint 
    ( const double       val                  , 
      const unsigned int ref                  , 
      const unsigned int mULPS = mULPS_double ) ;
    // ========================================================================
    /** check if the double value is actually equal to the integer value  
     *  @param val value to be compared with the unsigned integer 
     *  @param ref the reference unsigned integer number 
     *  @param mULPS the precision 
     *  @see Ostap::Math::lomont_compare_double 
     *  @see Ostap::Math::mULPS_double
     *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
     *  @date 2008-09-17
     */        
    inline bool equal_to_uint 
    ( const unsigned int ref                  , 
      const double       val                  ,
      const unsigned int mULPS = mULPS_double ) 
    { return equal_to_uint ( val , ref , mULPS ) ; }
    // ========================================================================
    /// signed sqrt 
    inline double signed_sqrt ( const double value ) 
    {
      return 
        0 <= value ? std::sqrt ( value ) : -std::sqrt( std::abs ( value ) ) ;
    }  
    // ========================================================================
    /** make dot-multiplication of two sequences using std::fma 
     *  \f$ r = \sum_i  x_i y_i \f$
     *  @param begin  (INPUT) start of the first sequence 
     *  @param end    (INPUT) end   of the first sequence 
     *  @param begin2 (INPUT) start iterator of the second sequence
     *  @return   "dot" product of two sequences 
     *
     *  @see http://en.cppreference.com/w/cpp/numeric/math/fma
     *  "...the function std::fma evaluates faster 
     *   (in addition to being more precise) than the expression x*y+z for 
     *   float, double, and long double arguments, respectively. "
     */
    template <class ITERATOR1, class ITERATOR2>
    inline double dot_fma
    ( ITERATOR1 begin  , 
      ITERATOR1 end    ,
      ITERATOR2 begin2 ) 
    {
      long double dot = 0 ;
      for ( ; begin != end ; ++begin, ++begin2 ) 
      { dot = std::fma ( *begin , *begin2 , dot ) ; }
      return dot  ;  
    }
    // ========================================================================
    /** make dot-multiplication of two sequences using std::fma 
     *  \f$ r = \sum_i  x_i y_i \f$
     *  @param x     (INPUT) the first sequence 
     *  @param begin (INPUT) start iterator of the second sequence
     *  @return   "dot" product of two sequences 
     */
    template <unsigned int N, class TYPE, class ITERATOR>
    inline double dot_fma
    ( const std::array<TYPE,N>& x     , 
      ITERATOR                  begin )  
    { return dot_fma ( x.begin() , x.end() , begin ) ; }
    // ========================================================================
    /** make dot-multiplication of two sequences using std::fma 
     *  \f$ r = \sum_i  x_i y_i \f$
     *  @param x     (INPUT) the first sequence 
     *  @param begin (INPUT) start iterator of the second sequence
     *  @return   "dot" product of two sequences 
     */
    template <unsigned int N, class TYPE1, class TYPE2>
    inline double dot_fma
    ( const std::array<TYPE1,N>& x , 
      const std::array<TYPE2,N>& y )  
    { return dot_fma ( x.begin() , x.end() , y.begin() ) ; }
    // ========================================================================
    
    
    // ========================================================================
    /// simple scaling of elements of non-constant sequence        
    template <class ITERATOR, typename SCALAR>
    void scale ( ITERATOR first  ,
                 ITERATOR last   , 
                 SCALAR   factor )
    { for ( ; first != last ; ++first ) { (*first) *= factor ; } }
    // ========================================================================
    /// shift all elements of non-constant sequence        
    template <class ITERATOR, typename SCALAR>
    void shift ( ITERATOR first  ,
                 ITERATOR last   , 
                 SCALAR   factor )
    { for ( ; first != last ; ++first ) { (*first) += factor ; } }
    // ========================================================================
    /// scale all elements of vector 
    template <class TYPE , typename SCALAR>
    void scale ( std::vector<TYPE>& vct , SCALAR factor ) 
    { scale    ( vct.begin() , vct.end () , factor ) ; }
    // ========================================================================
    /// shift all elements of vector 
    template <class TYPE , typename SCALAR>
    void shift ( std::vector<TYPE>& vct , SCALAR factor ) 
    { shift    ( vct.begin() , vct.end () , factor ) ; }
    // ========================================================================
    template <class ITERATOR> 
    void negate ( ITERATOR first , ITERATOR last ) 
    { for ( ; first != last ; ++first ) { (*first) = -(*first) ; } }
    // ========================================================================
    template <class TYPE> 
    void negate ( std::vector<TYPE>& vct ) 
    { negate ( vct.begin() , vct.end() ) ; }
    // ========================================================================
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
// The END 
// ============================================================================
#endif // OSTAP_MATH_H
// ============================================================================
