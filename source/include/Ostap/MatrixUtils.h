// ============================================================================
#ifndef OSTAP_MATRIXUTILS_H
#define OSTAP_MATRIXUTILS_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL 
// ============================================================================
#include <algorithm>
#include <functional>
#include <utility>
#include <cmath>
#include <array>
// ============================================================================
// ROOT
// ============================================================================
#include "Math/SMatrix.h"
#include "Math/SVector.h"
#include "Math/Point3D.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Math.h"
// ============================================================================
/** @file Ostap/MatrixUtils.h
 *  The collection of functions for manipulation with matrices and vectors.
 *  In particular it includes
 *     - (re)setting all elements of matrices and vectors
 *     - set the diagonal matrix to be proportional to unit matrix 
 *     - efficient scaling of matrices&vectors
 *     - minimal&maximal elements of matrices and vectors
 *     - indices of the minimal&maximal elements of matrices and vectors 
 *     - minima&maximal by absolute value elements of matrices and vectors 
 *     - indices of minima&maximal by absolute value elements of 
 *       matrices and vectors 
 *     - the trace of the square matrices
 *     - find the minimal/maximal/abs.minimal&abs.maximal 
 *       diagonal elements of square matrices
 *     - count number of elements which satisfy some criteria
 *     - count number of diagonal elements which satisfy some criteria
 *     - check the presence of elements which satisfy some criteria
 *     - check the presence of diagonal elements which satisfy some criteria
 *     - efficient element-by-element "equality" for matrices 
 *     - few specific "updates" (in the spirit of BLAS)
 *
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 */
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    /// generic class 
    template <class T>  class Equal_To ;
    // ========================================================================
    /// specialisation for vectors  
    template <class T, unsigned int D>
    struct Equal_To<ROOT::Math::SVector<T,D> > 
    {
    public:
      // ======================================================================
      /** constructor
       *  @see Ostap::Math::mULPS_double
       */
      Equal_To ( const unsigned int eps  = mULPS_double ) : m_cmp ( eps ) {}
      // ======================================================================
      /// comparison:
      inline bool operator()
      ( const ROOT::Math::SVector<T,D>& v1 , 
        const ROOT::Math::SVector<T,D>& v2 ) const  
      {
        return &v1 == &v2 || 
          std::equal ( v1.begin() , v1.end() , v2.begin() , m_cmp ) ;
      }
      /// compare with another vector type (e.g. double and float)
      template <class T2>
      inline bool operator()
      ( const ROOT::Math::SVector<T ,D>& v1 , 
        const ROOT::Math::SVector<T2,D>& v2 ) const  
      { return std::equal ( v1.begin() , v1.end() , v2.begin() , m_cmp ) ; }
      /// compare with another vector type (e.g. double and float)
      template <class T2>
      inline bool operator()
      ( const ROOT::Math::SVector<T2,D>& v1 , 
        const ROOT::Math::SVector<T ,D>& v2 ) const  
      { return std::equal ( v1.begin() , v1.end() , v2.begin() , m_cmp ) ; }
      // ======================================================================
    private:
      // ======================================================================
      /// the evaluator 
      Equal_To<T> m_cmp ;                                 // the evaluator 
      // ======================================================================
    } ;
    // ========================================================================
    /// specialisation for matrices 
    template <class T, unsigned int D1, unsigned int D2, class R>
    struct Equal_To<ROOT::Math::SMatrix<T,D1,D2,R> > 
    {
    public:
      // ======================================================================
      /** constructor
       *  @see Ostap::Math::mULPS_double
       */
      Equal_To ( const unsigned int eps  = mULPS_double ) : m_cmp ( eps ) {}
      // ======================================================================
      /// comparison:
      inline bool operator()
      ( const ROOT::Math::SMatrix<T,D1,D2,R>& v1 , 
        const ROOT::Math::SMatrix<T,D1,D2,R>& v2 ) const  
      {
        return &v1 == &v2 || 
          std::equal ( v1.begin() , v1.end() , v2.begin() , m_cmp ) ;
      }
      /// compare with another matrix type (e.g. symmetric...)
      template <class T2, class R2>
      inline bool operator()
      ( const ROOT::Math::SMatrix<T ,D1,D2,R >& v1 , 
        const ROOT::Math::SMatrix<T2,D1,D2,R2>& v2 ) const  
      {
        for ( unsigned int i = 0 ; i < D1 ; ++i ) 
        { for ( unsigned int j = 0 ; j < D2 ; ++j ) 
          { if ( !m_cmp ( v1(i,j) , v2(i,j) ) ) { return false ; } } } // RETURN 
        return true ; 
      }
      /// compare with another matrix type (e.g. symmetric...)
      template <class T2, class R2>
      inline bool operator()
      ( const ROOT::Math::SMatrix<T2,D1,D2,R2>& v1 , 
        const ROOT::Math::SMatrix<T ,D1,D2,R >& v2 ) const  
      {
        for ( unsigned int i = 0 ; i < D1 ; ++i ) 
        { for ( unsigned int j = 0 ; j < D2 ; ++j ) 
          { if ( !m_cmp ( v1(i,j) , v2(i,j) ) ) { return false ; } } } // RETURN 
        return true ; 
      }
      // ======================================================================
    private:
      // ======================================================================
      /// the evaluator 
      Equal_To<T> m_cmp ;                                 // the evaluator 
      // ======================================================================
    } ;
    // ========================================================================
    /** set all elements of vector equal to some scalar value 
     * 
     *  @code 
     *
     *  Ostap::Vector3      vct   = ... ;
     * 
     *  // set all elements to be equal to 10.0 :
     *  Ostap::Math::setToScalar( vct  , 10.0 ) ;
     *  // set all elements to be eual to 0 :
     *  Ostap::Math::setToScalar( vct  ) ;
     * 
     *  @endcode 
     *  
     *  @param m      (input/output)vector  to be modified
     *  @param value  (input) new value for all vector elements 
     *  @return number of modified vector elemenets (for consistency)
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-24
     */
    template <class T, unsigned int D>
    inline size_t 
    setToScalar
    ( ROOT::Math::SVector<T,D>& m , const T& value = T() ) 
    { 
      std::fill ( m.begin() , m.end() , value ) ; 
      return D ;
    } 
    // ========================================================================
    /** set all elements of matrix equal to some scalar value 
     * 
     *  @code 
     *
     *  Ostap::Matrix4x3    mtrx1 = ... ;
     *  Ostap::SymMatrix5x5 mtrx2 = ... ; 
     * 
     *  // set all elements to be equal to 10.0 :
     *  Ostap::Math::setToScalar( mtrx1 , 10.0 ) ;
     *  // set all elements to be eual to 0 :
     *  Ostap::Math::setToScalar( mtrx2 , 10.0 ) ;
     * 
     *  @endcode 
     *  
     *  @param m      (input/output) matrix to be modified
     *  @param value  (input) new value for all matrix elements 
     *  @return number of modified matrix elemenets 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-24
     */
    template <class T, unsigned int D1, unsigned int D2, class R>
    inline size_t 
    setToScalar
    ( ROOT::Math::SMatrix<T,D1,D2,R>& m , const T& value = T() ) 
    { 
      std::fill ( m.begin() , m.end() , value ) ; 
      return m.end() - m.begin() ;
    } 
    // ========================================================================
    /** set square matrix to be proportional to unit matrix 
     *  
     *  @code
     * 
     *  Ostap::Matrix4x4    mtrx1 = ... ;
     *  Ostap::SymMatrix3x3 mtrx2 = ... ;
     * 
     *  // set mtrx1 to be equal unit matrix 
     *  Ostap::Math::setToUnit ( mtrx1 ) ; 
     *  // set mtrx2 to be equal unit matrix, multiplied by 2  
     *  Ostap::Math::setToUnit ( mtrx2 , 2.0 ) ;
     *  @endcode 
     * 
     *  @param[in,out] m  matrix to be modified 
     *  @param[in] value  value to be used as diagonal elements
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-24
     */  
    template <class T, unsigned int D, class R>
    inline size_t 
    setToUnit
    ( ROOT::Math::SMatrix<T,D,D,R>& m , const T& value = T(1) ) 
    { 
      /// nullify the matrix:
      std::fill ( m.begin() , m.end() , T(0.0) ) ; 
      /// set diagonal elements 
      for ( unsigned int i = 0 ; i < D ; ++i ) { m(i,i) = value ; }
      return m.end() - m.begin() ;
    }
    // ========================================================================
    /** efficient scale all elements of the matrix
     *
     *  @code 
     * 
     *  Ostap::Matrix4x3    mtrx1 = ...  ;
     *  Ostap::SymMatrix5x5 mtrx2 = ... ;
     *
     *  // multiply all elements by 100:
     *  Ostap::Math::scale ( mtrx1 , 100.0   ) ;
     *
     *  // divide all elements by 4 :
     *  Ostap::Math::scale ( mtrx2 , 1.0/4.0 ) ;
     *  
     *  @endcode  
     * 
     *  @param m      (input/output) matrix to be modified
     *  @param value  (input) scaling coefficient 
     *  @return number of modified matrix elemenets 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-24
     */
    template <class T, unsigned int D1, unsigned int D2, class R>
    inline size_t 
    scale 
    ( ROOT::Math::SMatrix<T,D1,D2,R>& m , const T& value ) 
    { 
      typedef typename ROOT::Math::SMatrix<T,D1,D2,R>::iterator iterator ;
      for ( iterator it = m.begin() ;  m.end() != it ; ++it ) { (*it) *= value ; }
      return m.end() - m.begin() ;
    } 
    // ========================================================================
    /** efficient scale all elements of the vector
     *
     *  @code 
     * 
     *  Ostap::Vector4      vct = ...  ;
     *
     *  // multiply all elements by 100:
     *  Ostap::Math::scale ( vct , 100.0   ) ;
     *
     *  @endcode  
     * 
     *  @param m      (input/output) vector to be modified
     *  @param value  (input) scaling coefficient 
     *  @return number of modified vector elemenets  (for consistency) 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-24
     */
    template <class T, unsigned int D>
    inline size_t 
    scale 
    ( ROOT::Math::SVector<T,D>& m , const T& value ) 
    { 
      typedef typename ROOT::Math::SVector<T,D>::iterator iterator ;
      for ( iterator it = m.begin() ;  m.end() != it ; ++it ) { (*it) *= value ; }
      return D ;
    } 
    // ========================================================================
    /** @struct _AbsCompare
     *  The trivial structure for comparison of "numbers" by the absolute value 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-25
     */
    template <class T>
    struct _AbsCompare 
    {
      inline bool operator() ( const T v1 , const T v2 ) const 
      { return std::abs( v1 ) < std::abs( v2 ) ; }
    } ;
    // ========================================================================
    /** find the maximal element in matrix 
     *  @param m (input) matrix to be studied
     *  @return the maximal element 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-24
     */
    template <class T,unsigned int D1,unsigned int D2, class R>
    inline T 
    max_element 
    ( const ROOT::Math::SMatrix<T,D1,D2,R>& m ) 
    { return *std::max_element ( m.begin() , m.end() ) ; }
    // ========================================================================
    /** find the minimal element in matrix 
     *  @param m (input) matrix to be studied
     *  @return the minimal element 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-24
     */
    template <class T,unsigned int D1,unsigned int D2, class R>
    inline T 
    min_element 
    ( const ROOT::Math::SMatrix<T,D1,D2,R>& m ) 
    { return *std::min_element ( m.begin() , m.end() ) ; }
    // ========================================================================
    /** find the maximal element in vector 
     *  @param m (input) vector to be studied
     *  @return the maximal element 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-24
     */
    template <class T,unsigned int D>
    inline T 
    max_element 
    ( const ROOT::Math::SVector<T,D>& m ) 
    { return *std::max_element ( m.begin() , m.end() ) ; }
    // ========================================================================
    /** find the minimal element in vector 
     *  @param m (input) vector to be studied
     *  @return the minimal element 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-24
     */
    template <class T,unsigned int D>
    inline T 
    min_element 
    ( const ROOT::Math::SVector<T,D>& m ) 
    { return *std::min_element ( m.begin() , m.end() ) ; }
    // ========================================================================
    /** find the element in matrix with the maximal absolute value  
     *
     *  @code
     *  
     *  const Ostap::Matrix4x3 mtrx1 = ... ;
     *  const Ostap::Matrix4x3 mtrx2 = ... ;
     * 
     *  // find the maximal absolute difference between elements : 
     *  const double diff  = Ostap::Math::maxabs_element ( mtrx1 - mtrx2 ) ;
     *  if ( diff < tolerance ) 
     *   {
     *      std::cout << " matrices are almost identical " << std::endl ;
     *   }
     *  
     *  @endcode 
     *  
     *  @param m (input) matrix to be studied
     *  @return the element with the maximal absolute value 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-24
     */
    template <class T,unsigned int D1,unsigned int D2, class R>
    inline T 
    maxabs_element 
    ( const ROOT::Math::SMatrix<T,D1,D2,R>& m ) 
    { return *std::max_element ( m.begin() , m.end()  , _AbsCompare<T>() ) ; }
    // ========================================================================
    /** find the element in matrix with the minimal absolute value  
     *  @param m (input) matrix to be studied
     *  @return the element with the minimal absolute value 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-24
     */
    template <class T,unsigned int D1,unsigned int D2, class R>
    inline T 
    minabs_element 
    ( const ROOT::Math::SMatrix<T,D1,D2,R>& m ) 
    { return *std::min_element ( m.begin() , m.end()  , _AbsCompare<T>() ) ; }
    // ========================================================================
    /** find an index of the  maximal element in matrix 
     *  @param m (input) matrix to be studied
     *  @param cmp comparison criteria
     *  @return the pair (i,j)-index of the maximal element 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-24
     */
    template <class T,unsigned int D1,unsigned int D2, class R, class CMP>
    inline std::pair<unsigned int,unsigned int> 
    ind_max_element 
    ( const ROOT::Math::SMatrix<T,D1,D2,R>& m , CMP cmp ) 
    { 
      std::pair<unsigned int,unsigned int> result (0,0) ;
      T tmp = m(0,0) ;
      for ( unsigned int i = 0 ; i < D1 ; ++i ) 
      {
        for ( unsigned int j = 0 ; j < D2 ; ++j ) 
        {
          const T val = m(i,j) ;
          if ( !cmp( tmp , val ) ) { continue ; }
          tmp = val ; 
          result.first  = i ; 
          result.second = j ;
        }
      }
      return result ;
    } 
    // ========================================================================
    /** find an index of the maximal element in symmetric matrix 
     *  @param m (input) symmetric matrix to be studied
     *  @param cmp comparison criteria
     *  @return the pair (i,j)-index of the maximal element 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-24
     */
    template <class T,unsigned int D, class CMP>
    inline std::pair<unsigned int,unsigned int> 
    ind_max_element 
    ( const ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >& m , CMP cmp ) 
    { 
      std::pair<unsigned int,unsigned int> result (0,0) ;
      T tmp = m(0,0) ;
      for ( unsigned int i = 0 ; i < D ; ++i ) 
      {
        for ( unsigned int j = i ; j < D ; ++j ) 
        {
          const T val = m(i,j) ;
          if ( !cmp ( tmp , val ) ) { continue ; }
          tmp = val ; 
          result.first  = i ; 
          result.second = j ;
        }
      }
      return result ;
    }     
    // ========================================================================
    /** find an index of the  maximal element in matrix 
     *  @param m (input) matrix to be studied
     *  @return the pair (i,j)-index of the maximal element 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-24
     */
    template <class T,unsigned int D1,unsigned int D2, class R>
    inline std::pair<unsigned int,unsigned int> 
    ind_max_element 
    ( const ROOT::Math::SMatrix<T,D1,D2,R>& m ) 
    { return ind_max_element ( m , std::less<T>() ) ; }
    // ========================================================================
    /** find an index of the minimal element in matrix 
     *  @param m (input) matrix to be studied
     *  @param cmp comparison criteria
     *  @return the pair (i,j)-index of the maximal element 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-24
     */
    template <class T,unsigned int D1,unsigned int D2, class R, class CMP>
    inline std::pair<unsigned int,unsigned int> 
    ind_min_element 
    ( const ROOT::Math::SMatrix<T,D1,D2,R>& m , CMP cmp ) 
    { 
      std::pair<unsigned int,unsigned int> result (0,0) ;
      T tmp = m(0,0) ;
      for ( unsigned int i = 0 ; i < D1 ; ++i ) 
      {
        for ( unsigned int j = 0 ; j < D2 ; ++j ) 
        {
          const T val = m(i,j) ;
          if ( !cmp( val , tmp ) ) { continue ; }
          tmp = val ; 
          result.first  = i ; 
          result.second = j ;
        }
      }
      return result ;
    } 
    // ========================================================================
    /** find an index of the minimal element in symmetric matrix 
     *  @param m (input) symmetric matrix to be studied
     *  @param cmp comparison criteria
     *  @return the pair (i,j)-index of the maximal element 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-24
     */
    template <class T,unsigned int D, class CMP>
    inline std::pair<unsigned int,unsigned int> 
    ind_min_element 
    ( const ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >& m , CMP cmp ) 
    { 
      std::pair<unsigned int,unsigned int> result (0,0) ;
      T tmp = m(0,0) ;
      for ( unsigned int i = 0 ; i < D ; ++i ) 
      {
        for ( unsigned int j = i ; j < D ; ++j ) 
        {
          const T val = m(i,j) ;
          if ( !cmp ( tmp , val ) ) { continue ; }
          tmp = val ; 
          result.first  = i ; 
          result.second = j ;
        }
      }
      return result ;
    } 
    // ========================================================================
    /** find an index of the minimal element in the matrix 
     *
     *  @code 
     *  
     *  const Ostap::Matrix4x4    mtrx1 = ... ;
     *  const Ostap::SymMatrix5x5 mtrx2 = ... ;
     * 
     *  typedef std::pair<unsigned int,unsigned int> PAIR ;
     *  
     *  // get the minimal element in mtrx1:
     *  PAIR index1 = Ostap::Math::ind_max_element ( mtrx1 ) ;
     *
     *  // get the minimal element in mtrx2:
     *  PAIR index2 = Ostap::Math::ind_max_element ( mtrx2 ) ;
     *
     *  @endcode 
     *
     *  @param m (input) matrix to be studied 
     *  @return the pair of indices for (the first) minimal element 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-24
     */
    template <class T, unsigned int D1, unsigned int D2, class R>
    inline std::pair<unsigned int,unsigned int> 
    ind_min_element 
    ( const ROOT::Math::SMatrix<T,D1,D2,R>& m ) 
    { return ind_min_element( m , std::less<T>() ) ; } 
    // ========================================================================
    /** find an index of the maximal element in the vector 
     *  @param m (input) vector to be studied
     *  @param cmp comparison criteria
     *  @return the index of the maximal element 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-24
     */
    template <class T,unsigned int D, class CMP>
    inline unsigned int 
    ind_max_element 
    ( const ROOT::Math::SVector<T,D>& m , CMP cmp )
    { return std::max_element( m.begin() , m.end() , cmp ) - m.begin() ; }    
    // ========================================================================
    /** find an index of the minimal element in the vector 
     *  @param m (input) vector to be studied
     *  @param cmp comparison criteria
     *  @return the index of the minimal  element 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-24
     */
    template <class T,unsigned int D, class CMP>
    inline unsigned int 
    ind_min_element 
    ( const ROOT::Math::SVector<T,D>& m , CMP cmp )
    { return std::min_element( m.begin() , m.end() , cmp ) - m.begin() ; }
    // ========================================================================
    /** find an index of the maximal element in the vector 
     *  @param m (input) vector to be studied
     *  @return the index of the maximal element 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-24
     */
    template <class T,unsigned int D>
    inline unsigned int 
    ind_max_element 
    ( const ROOT::Math::SVector<T,D>& m )
    { return std::max_element( m.begin() , m.end() ) - m.begin() ; }    
    // ========================================================================
    /** find an index of the minimal element in the vector 
     *  @param m (input) vector to be studied
     *  @return the index of the minimal element 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-24
     */
    template <class T,unsigned int D>
    inline unsigned int 
    ind_min_element 
    ( const ROOT::Math::SVector<T,D>& m )
    { return std::min_element ( m.begin() , m.end() ) - m.begin() ; }    
    // ========================================================================
    /** find an index of the element with the maximal absolute value 
     *
     *  @code 
     *  
     *  const Ostap::Matrix4x4    mtrx1 = ... ;
     *  const Ostap::SymMatrix5x5 mtrx2 = ... ;
     * 
     *  typedef std::pair<unsigned int,unsigned int> PAIR ;
     *  
     *  // get the element with the maximal absolute value in mtrx1:
     *  PAIR index1 = Ostap::Math::ind_maxabs_element ( mtrx1 ) ;
     *
     *  // get the element with the maximal absolute value in mtrx2:
     *  PAIR index2 = Ostap::Math::ind_maxabs_element ( mtrx2 ) ;
     *
     *  @endcode 
     *
     *  @param m (input) matrix to be studied 
     *  @return the pair of indices for the element with 
     *          the maximal absolute value 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-24
     */
    template <class T, unsigned int D1, unsigned int D2, class R>
    inline std::pair<unsigned int,unsigned int> 
    ind_maxabs_element 
    ( const ROOT::Math::SMatrix<T,D1,D2,R>& m ) 
    { return ind_max_element( m , _AbsCompare<T>() ) ; } 
    // ========================================================================
    /** find an index of the element with the minimal absolute value 
     *
     *  @code 
     *  
     *  const Ostap::Matrix4x4    mtrx1 = ... ;
     *  const Ostap::SymMatrix5x5 mtrx2 = ... ;
     * 
     *  typedef std::pair<unsigned int,unsigned int> PAIR ;
     *  
     *  // get the element with the minimal absolute value in mtrx1:
     *  PAIR index1 = Ostap::Math::ind_minabs_element ( mtrx1 ) ;
     *
     *  // get the element with the minimal absolute value in mtrx2:
     *  PAIR index2 = Ostap::Math::ind_minabs_element ( mtrx2 ) ;
     *
     *  @endcode 
     *
     *  @param m (input) matrix to be studied 
     *  @return the pair of indices for the element with 
     *          the maximal absolute value 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-24
     */    
    template <class T, unsigned int D1, unsigned int D2, class R>
    inline std::pair<unsigned int,unsigned int> 
    ind_minabs_element 
    ( const ROOT::Math::SMatrix<T,D1,D2,R>& m ) 
    { return ind_min_element( m , _AbsCompare<T>() ) ; } 
    // ========================================================================
    /** find an index of the element with maximal absolute value
     *  @param m (input) vector to be studied
     *  @return the index of the element with the maximal absolute value 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-24
     */
    template <class T,unsigned int D>
    inline unsigned int 
    ind_maxabs_element 
    ( const ROOT::Math::SVector<T,D>& m )
    { return ind_max_element ( m , _AbsCompare<T>() ) ; }    
    // ========================================================================
    /** find an index of the element with minimal absolute value 
     *  @param m (input) vector to be studied
     *  @return the index of the element with minimal absolute value 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-24
     */
    template <class T,unsigned int D>
    inline unsigned int 
    ind_minabs_element 
    ( const ROOT::Math::SVector<T,D>& m )
    { return ind_min_element ( m , _AbsCompare<T>() ) ; }
    // ========================================================================
    /** evaluate the trace (sum of diagonal elements) of the square matrix 
     *
     *  @code 
     *  
     *  const Ostap::Matrix4x4    mtrx1 = ... ;
     *  const Ostap::SymMatrix3x3 mtrx2 = ... ;
     * 
     *  // evaluate the trace of mtrx1:
     *  const double trace1 = Ostap::Math::trace ( mtrx1 ) ;
     *  // evaluate the trace of mtrx2:
     *  const double trace2 = Ostap::Math::trace ( mtrx2 ) ;
     *
     *  @endcode 
     *
     *  @param m (input) matrix to be studied 
     *  @return trace (sum of diagonal elements) of the matrix 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-24
     */
    template <class T, unsigned int D, class R>
    inline T 
    trace 
    ( const ROOT::Math::SMatrix<T,D,D,R>& m ) 
    {
      T result = m ( 0 , 0 ) ;
      for ( unsigned int i = 1 ; i < D ; ++i ) { result += m(i,i) ; }
      return result ;
    }
    // ========================================================================
    /** evaluate the trace (sum of diagonal elements) for matrix epxression 
     *  @param m (input) matrix expression
     *  @return trace (sum of diagonal elements) of the matrix 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-05-24
     */
    template <class A, class T, unsigned int D, class R>
    inline T
    trace
    ( const ROOT::Math::Expr<A,T,D,D,R>& m ) 
    {
      T result =  m ( 0 , 0 ) ;
      for ( unsigned int i = 1 ; i < D ; ++i ) { result += m ( i , i ) ; }
      return result ; 
    }
    // ========================================================================
    /** find the minimal diagonal element 
     *  @param   m (input) matrix to be studied 
     *  @param   cmp comparison criteria
     *  @return "min" diagonal element (in the sense of comparison criteria)
     */
    template <class T, unsigned int D, class R, class CMP>
    inline T 
    min_diagonal 
    ( const ROOT::Math::SMatrix<T,D,D,R>& m , CMP cmp ) 
    {
      T result = m(0,0);
      for ( unsigned int i = 1 ; i < D ; ++i ) 
      { 
        const T value = m(i,i) ;  
        if ( cmp ( value , result ) ) { result = value ; }
      }
      return result ;
    }
    // ========================================================================
    /** find the maximal diagonal element 
     *  @param   m (input) square matrix to be studied 
     *  @param   cmp comparison criteria
     *  @return "max" diagonal element (in the sense of comparison criteria)
     */
    template <class T, unsigned int D, class R, class CMP>
    inline T 
    max_diagonal 
    ( const ROOT::Math::SMatrix<T,D,D,R>& m , CMP cmp ) 
    {
      T result = m(0,0);
      for ( unsigned int i = 1 ; i < D ; ++i ) 
      { 
        const T value = m(i,i) ;  
        if ( cmp ( result , value ) ) { result = value ; }
      }
      return result ;
    } 
    // ========================================================================
    /** find the maximal diagonal element of square matrix 
     *
     *  @code 
     *
     *  const Ostap::Matrix4x4    mtrx1 = ... ;
     *  const Ostap::SymMatrix3x3 mtrx2 = ... ;
     *  
     *  // find the maximal diagonal element of mtrx1 
     *  const double m1 = Ostap::Math::max_diagonal ( mtrx1 ) ;
     *
     *  // find the maximal diagonal element of mtrx2 
     *  const double m2 = Ostap::Math::max_diagonal ( mtrx2 ) ;
     *
     *  @endcode
     *
     *  @param m (input) square matrix to be studied 
     *  @return the maximal diagonal element 
     *  @author Vanya BELYAEV ibelyaev@physics.cyr.edu
     *  @date 2006-05-24
     */
    template <class T, unsigned int D, class R>
    inline T 
    max_diagonal 
    ( const ROOT::Math::SMatrix<T,D,D,R>& m ) 
    { return max_diagonal( m , std::less<T>() ) ; }
    // ========================================================================
    /** find the maximal diagonal element of the square matrix 
     *
     *  @code 
     *
     *  const LHCb::Vertex* v = ..
     *  const Ostap::Matrix3x3& covariance = v->covMatrix() ;
     *  
     *  // find the minimal diagonal element of covariance matrix:
     *  const double m1 = Ostap::Math::min_diagonal ( covariance  ) ;
     *
     *  if ( 0 >= 0 ) 
     *  {
     *     std::err 
     *      << " Invlid covarinace matrix" 
     *      << " Non-positive elements on diagonal << m1 << std::endl ;
     *  }
     *
     *  @endcode
     *
     *  @param m (input) square matrix to be studied 
     *  @return the maximal diagonal element 
     *  @author Vanya BELYAEV ibelyaev@physics.cyr.edu
     *  @date 2006-05-24
     */
    template <class T, unsigned int D, class R>
    inline T 
    min_diagonal 
    ( const ROOT::Math::SMatrix<T,D,D,R>& m ) 
    { return min_diagonal( m , std::less<T>() ) ; }
    // ========================================================================
    /** find the diagonal element of square matrix with maximal absolute value 
     *
     *  @code 
     *
     *  const Ostap::Matrix4x4    mtrx1 = ... ;
     *  const Ostap::SymMatrix3x3 mtrx2 = ... ;
     *  
     *  // find the diagonal element of mtrx1 with the maximal absolute value 
     *  const double m1 = Ostap::Math::max_diagonal ( mtrx1 ) ;
     *
     *  // find the diagonal element of mtrx1 with the maximal absolute value 
     *  const double m2 = Ostap::Math::max_diagonal ( mtrx2 ) ;
     *
     *  @endcode
     *
     *  @param m (input) square matrix to be studied 
     *  @return the diagonal element withmaximal absolute value 
     *  @author Vanya BELYAEV ibelyaev@physics.cyr.edu
     *  @date 2006-05-24
     */
    template <class T, unsigned int D, class R>
    inline T 
    maxabs_diagonal 
    ( const ROOT::Math::SMatrix<T,D,D,R>& m ) 
    { return max_diagonal( m , _AbsCompare<T>() ) ; }
    // ========================================================================
    /** find the diagonal element of the square matrix with 
     *  the minimal absolute value 
     *
     *  @code 
     *
     *  const LHCb::Vertex* v = ..
     *  const Ostap::Matrix3x3& covariance = v->covMatrix() ;
     *  
     *  // find the diagonal element with minimal absolute value 
     *  const double m1 = Ostap::Math::min_diagonal ( covariance  ) ;
     *
     *  if ( 0.001 * Ostap::Units::micrometer ) 
     *  {
     *     std::err 
     *      << " Non-realistic element on diagonal << m1 << std::endl ;
     *  }
     *
     *  @endcode
     *
     *  @param m (input) square matrix to be studied 
     *  @return the diagonal element with the minimal absolute value 
     *  @author Vanya BELYAEV ibelyaev@physics.cyr.edu
     *  @date 2006-05-24
     */
    template <class T, unsigned int D, class R>
    inline T 
    minabs_diagonal 
    ( const ROOT::Math::SMatrix<T,D,D,R>& m ) 
    { return min_diagonal( m , _AbsCompare<T>() ) ; }
    // ========================================================================
    /** count the number of elements in matrix, which satisfy the certain criteria
     * 
     *  @code
     * 
     *  const Ostap::Matrix4x4 matrix = ... ;
     *  
     *  // number of NULL elements:
     *  const size_t nulls = 
     *      Ostap::Math::count_if ( matrix , 
     *       std::bind2nd( std::equal_to<double>() , 0.0 ) ) ;
     *
     *  // number of elements in excess of 100.0 
     *  const size_t large = 
     *      Ostap::Math::count_if ( matrix , 
     *       std::bind2nd( std::greater<double>() , 100.0 ) ) ;
     * 
     *  // number of elements which are less then 0.01 in absolute value 
     *  const size_t small = 
     *      Ostap::Math::count_if ( matrix , 
     *       std::bind2nd( _AbsCompare<double>() , 0.01 ) ) ;
     *  
     * 
     *  @endcode 
     *  
     *  @param m    (input) matrix to be studied 
     *  @param pred (input) predicate to be tested 
     *  @return number of elements for which the predicate is valid 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-04-24
     */
    template <class T, unsigned int D1, unsigned int D2, class R, class P>
    inline size_t 
    count_if 
    ( const ROOT::Math::SMatrix<T,D1,D2,R>& m , P pred )
    { return std::count_if ( m.begin() , m.end() , pred ) ; }
    // ========================================================================
    /** count the number of elements in matrix, which satisfy the certain criteria
     * 
     *  @code
     * 
     *  const Ostap::SymMatrix4x4 matrix = ... ;
     *  
     *  // number of NULL elements:
     *  const size_t nulls = 
     *      Ostap::Math::count_if ( matrix , 
     *       std::bind2nd( std::equal_to<double>() , 0.0 ) ) ;
     *
     *  // number of elements in excess of 100.0 
     *  const size_t large = 
     *      Ostap::Math::count_if ( matrix , 
     *       std::bind2nd( std::greater<double>() , 100.0 ) ) ;
     * 
     *  // number of elements which are less then 0.01 in absolute value 
     *  const size_t small = 
     *      Ostap::Math::count_if ( matrix , 
     *       std::bind2nd( Ostap::Math::_AbsCompare<double>() , 0.01 ) ) ;
     *  
     *  @endcode 
     *  
     *  If one needs check the presence of elements ("at least one")
     *  a bit more efficient algorithm 
     *   <c>Ostap::Math::check_if</c> should be used 
     *  @see Ostap::Math::check_if 
     *
     *  @param m    (input) symmetric matrix to be studied 
     *  @param pred (input) predicate to be tested 
     *  @return number of elements for which the predicate is valid 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-04-24
     */
    template <class T, unsigned int D, class P>
    inline size_t 
    count_if 
    ( const ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >& m , P pred )
    { 
      size_t result = 0 ;
      for ( unsigned int i = 0 ; i < D ; ++i ) 
      {
        if ( pred ( m ( i , i ) ) ) { result += 1 ; }
        for ( unsigned int j = i + 1 ; j < D ; ++j ) 
        {
          if ( pred (  m ( i , j ) ) ) { result +=2 ; }  // ATTENTION! 
        }
      }
      return result ;
    } 
    // ========================================================================
    /** count number of diagonal elements in matrix, 
     *  which satisfy certain criteria
     * 
     *  @code 
     * 
     *  const LHCb::Vertex* v = ... ;
     *  const Ostap::SymMatrix3x3& covariance = v->covMatrix() ;
     *
     *  // count number of VERY small (and negative) diagonal elements:
     *  const size_t bad = 
     *   Ostap::Math::cound_diagonal( covariance , 
     *   std::bind2nd( std::less<double>() , 0.01 * Ostap::Units::micrometer ) ;
     *  if ( 0 != bad ) 
     *   {
     *      std::cerr << " #bad diagonal elements is " << bad << std::endl ;  
     *   }
     *  
     *  @endcode 
     *  
     *  @param m    (input) square matrix to be studied 
     *  @param pred (input) predicate to be tested 
     *  @return number of diagonal elements for which the predicate is valid 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-04-24
     */
    template <class T, unsigned int D, class R, class P>
    inline size_t 
    count_diagonal
    ( const ROOT::Math::SMatrix<T,D,D,R>& m , P pred )
    {
      size_t result = 0 ;
      for ( unsigned int i = 0 ; i < D ; ++i ) 
      { if ( pred ( m ( i , i ) ) ) { result += 1 ; } }
      return result ;
    } 
    // ========================================================================
    /** check the presence of at least one element which satisfy the 
     *  criteria
     *
     *  @code 
     * 
     *  const LHCb::Vertex* v = ... ;
     *  const Ostap::SymMatrix3x3& covariance = v->covMatrix() ;
     *
     *  // check for "infinities" 
     *  const bool bad = 
     *   Ostap::Math::check_if( covariance , 
     *   std::bind2nd( std::greater<double>() , 1 * Ostap::Units::meter ) ;
     *  if ( bad ) 
     *   {
     *      std::cerr << " bad elements are detected " << std::endl ; 
     *   }
     *  
     *  @endcode 
     *
     *  In general this algorithm is faster than 
     *   <c>Ostap::Math::count_if</c>
     *  @see Ostap::Math::count_if 
     *
     *  @param m    (input)    matrix to be checked 
     *  @param pred (input) predicate tobe tested 
     *  @return true if at least one element is in the matrix 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-04-24
     */
    template <class T, unsigned int D1, unsigned int D2, class R, class P>
    inline bool
    check_if 
    ( const ROOT::Math::SMatrix<T,D1,D2,R>& m , P pred )
    { return m.end() != std::find_if ( m.begin() , m.end() , pred ) ; }
    // ========================================================================
    /** check the presence of at least one diagonal element which satisfy the 
     *  criteria
     *
     *  @code 
     * 
     *  const LHCb::Vertex* v = ... ;
     *  const Ostap::SymMatrix3x3& covariance = v->covMatrix() ;
     *
     *  // check for "almost nulls"
     *  const bool bad = 
     *   Ostap::Math::check_diagonal( covariance , 
     *   std::bind2nd( std::less<double>() , 0.001 * Ostap::Units::micrometer ) ;
     *  if ( bad ) 
     *   {
     *      std::cerr << " bad diagonal elements are detected " << std::endl ; 
     *   }
     *  
     *  @endcode 
     *
     *  @param m    (input) square matrix to be checked 
     *  @param pred (input) predicate to be tested 
     *  @return true if at least one element is in the matrix 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-04-24
     */
    template <class T, unsigned int D, class R, class P>
    inline bool
    check_diagonal 
    ( const ROOT::Math::SMatrix<T,D,D,R>& m , P pred )
    { 
      for ( unsigned int i = 0 ; i < D ; ++i ) 
      { if ( pred ( m ( i , i ) ) ) { return true ; } }
      return false ;
    } 
    // ========================================================================
    /** check the "equality" of the two matrices by checking 
     *  element-by-element: true == pred( m1(i,j) , m2(i,j) ) 
     *  
     *  It is an efficient way to compare the matrices in the different
     *  represenattions, e.g. symmetrical and non symmetrical matrices 
     *  of the same size. Also the "tolerance" coudl be introduced 
     *  for equality
     *
     *  @code 
     * 
     *  const Ostap::Matrix4x4 m1 = ... ;
     *  const Ostap::SymMatrix4x4 m2 = ... ;
     *  
     *   // comparison criteria:
     *   struct Equal 
     *   {
     *     Equal ( const double value ) : m_threshold ( value ) {} ;
     *     bool operator() ( const double v1 , const double v2 ) const 
     *     {
     *         return ::fabs( v1 , v2 ) < m_threshold ;
     *     }
     *    private: 
     *    double m_threshold ;
     *   } ;
     *  
     *   // "compare" the matrices 
     *   const bool eq = 
     *      Ostap::Math::equal_if ( m1 , m2 , Equal(0.001) ) ;
     *  
     *  @endcode 
     *  @param m1   (input) the first matrix to be checked 
     *  @param m2   (input) the second matrix to be checked 
     *  @param pred (input) predicate to be tested 
     *  @return true if at least once false == pred( m1(i,j) , m2(i,j) )
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-04-24
     */    
    template <class T1 , class T2 ,
              unsigned int D1 , unsigned int D2 ,
              class R1 , class R2 , class P>
    inline bool 
    equal_if  
    ( const ROOT::Math::SMatrix<T1,D1,D2,R1>& m1 , 
      const ROOT::Math::SMatrix<T2,D1,D2,R2>& m2 , P pred )
    { 
      for ( unsigned int i = 0 ; i < D1 ; ++i ) 
      {
        for ( unsigned int j = 0 ; j < D2 ; ++j ) 
        {
          if ( !pred( m1(i,j) , m2(i,j) ) ) { return false ; } // RETURN 
        }
      }
      return true ;
    } 
    // ========================================================================
    /** check the "equality" of the two matrices by checking 
     *  element-by-element: true == pred( m1(i,j) , m2(i,j) ) 
     *
     *  The specialization for matrices of the same representation.
     *  From the first principles it should be much more efficient. 
     *
     *  @code 
     * 
     *  const Ostap::Matrix4x4 m1 = ... ;
     *  const Ostap::Matrix4x4 m2 = ... ;
     *  
     *   // comparison criteria:
     *   struct Equal
     *   {
     *     Equal ( const double value ) : m_threshold ( value ) {} ;
     *     bool operator() ( const double v1 , const double v2 ) const 
     *     {
     *         return ::fabs( v1 , v2 ) < m_threshold ;
     *     }
     *    private: 
     *    double m_threshold ;
     *   } ;
     *  
     *   // "compare" the matrices 
     *   const bool eq = 
     *      Ostap::Math::equal_if ( m1 , m2 , Equal(0.001) ) ;
     *  
     *  @endcode 
     *  @param m1   (input) the first matrix to be checked 
     *  @param m2   (input) the second matrix to be checked 
     *  @param pred (input) predicate to be tested 
     *  @return true if at least once false == pred( m1(i,j) , m2(i,j) )
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2006-04-24
     */
    template <class T,unsigned int D1, unsigned int D2,class R,class P>
    inline bool 
    equal_if  
    ( const ROOT::Math::SMatrix<T,D1,D2,R>& m1 , 
      const ROOT::Math::SMatrix<T,D1,D2,R>& m2 , P pred )
    { return std::equal ( m1.begin() , m1.end() , m2.begin() , pred ) ; } 
    // ========================================================================
    // UPDATE
    // =========================================================================
    /** update the symmetric matrix according to the rule m +=  s*v*v^T
     *
     *  @param left the symmetric matrix to be updated 
     *  @param vect the vector 
     *  @param scale the scale factor
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2008-01-10
     */
    template <class T, class T2, unsigned int D>
    void update 
    ( ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > & left  , 
      const ROOT::Math::SVector<T2,D>&                         vect  , 
      const double                                             scale ) 
    {
      for ( unsigned int i = 0 ; i < D ; ++i ) 
      {
        for ( unsigned int j = i ; j < D ; ++j ) 
        { left ( i , j ) += scale * vect(i) * vect(j) ; }
      }
    }
    // =========================================================================
    /** update the symmetric matrix according to the rule m +=  s*v*v^T
     *
     *  @param left the symmetric matrix to be updated 
     *  @param vect the vector 
     *  @param scale the scale factor
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2008-01-10
     */
    template <class T, class B, class T2, unsigned int D>
    void update 
    ( ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > & left  , 
      const ROOT::Math::VecExpr<B,T2,D>&                       vect  , 
      const double                                             scale ) 
    {
      for ( unsigned int i = 0 ; i < D ; ++i ) 
      {
        for ( unsigned int j = i ; j < D ; ++j ) 
        { left ( i , j ) += scale * vect(i) * vect(j) ; }
      }
    }
    // =========================================================================
    /** update the matrix according to the rule m +=  s*v1*v2^T
     *
     *  @param left the matrix to be updated 
     *  @param vct1 the first vector 
     *  @param vct2 the second vector 
     *  @param scale the scale factor
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2008-01-10
     */
    template <class T, class R, class T2, class T3, unsigned int D1,unsigned int D2>
    void update 
    ( ROOT::Math::SMatrix<T,D1,D2,R>  & left        , 
      const ROOT::Math::SVector<T2,D1>& vct1        , 
      const ROOT::Math::SVector<T3,D2>& vct2        , 
      const double                      scale = 1.0 ) 
    {
      for ( unsigned int i = 0 ; i < D1 ; ++i ) 
      {
        for ( unsigned int j = 0 ; j < D2 ; ++j ) 
        { left ( i , j ) += scale * vct1(i) * vct2(j) ; }
      }
    }
    // =========================================================================
    /**  useful shortcut for product of vector, matrix and vector (v1^T*M*v2)
     *  @param vct1 the first vector 
     *  @param mtrx the matrix  
     *  @param vct2 the second vector 
     *  @return the product (v1^T*M*v2) 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2008-01-19
     */
    template<class T, class T1, class T2, class R, unsigned int D1, unsigned int D2>
    T mult 
    ( const ROOT::Math::SVector<T1,D1>&     vct1  ,
      const ROOT::Math::SMatrix<T,D1,D2,R>& mtrx  , 
      const ROOT::Math::SVector<T2,D2>&     vct2  ) 
    { 
      return ROOT::Math::Dot ( vct1 , mtrx * vct2 ) ; 
    }
    // =========================================================================
    /** update the symmetric matrix according to the rule m +=  scale * ( m + m^T )  
     *
     *  @param left the matrix  to be updated 
     *  @param right the matrix to be "symmetrized"  
     *  @param scale the scale factor to be applied 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2008-01-10
     */
    template <class T, class T2, class R, unsigned int D>
    void update 
    ( ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >&   left        , 
      const ROOT::Math::SMatrix<T2,D,D,R>&                      right       , 
      const double                                              scale = 1.0 ) 
    {
      for ( unsigned int i = 0 ; i < D ; ++i ) 
      {
        for ( unsigned int j = i ; j < D ; ++j ) 
        { left ( i , j ) += scale * ( right ( i , j ) + right ( j , i ) ) ; }
      }
    }
    // =========================================================================
    /** update the symmetric matrix according to the rule m +=  scale * ( m + m^T )  
     *
     *  @param left the matrix  to be updated 
     *  @param right the matrix to be "symmetrized"  
     *  @param scale the scale factor to be applied 
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date 2008-01-10
     */
    template <class T, class T2, class B, class R, unsigned int D>
    void update 
    ( ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >&   left        , 
      const ROOT::Math::Expr<B,T2,D,D,R>&                       right       , 
      const double                                              scale = 1.0 ) 
    {
      for ( unsigned int i = 0 ; i < D ; ++i ) 
      {
        for ( unsigned int j = i ; j < D ; ++j ) 
        { left ( i , j ) += scale * ( right ( i , j ) + right ( j , i ) ) ; }
      }
    }
    // ========================================================================
    /** inversion of symmetric positively defined matrices
     *  1) try fast method based on Cholesky's decomposiion 
     *  2) in case of failure, swith to the regular inversion
     *  @param what (INPUT) matrix to be inverted 
     *  @param result (UPDATE) th einverse matrix 
     *  @return problem flag     
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  - thanks to Manuel Schiller 
     *  @date 2015-02-11
     */
    template <class T, unsigned int D>
    inline int inverse 
    ( const ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >& what   , 
      ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >&       result ) 
    {
      int ifail = 0 ;     result = what.InverseChol ( ifail ) ;
      if ( 0 != ifail ) { result = what.Inverse     ( ifail ) ; }
      return ifail ;  
    }
    // ========================================================================
    /** Get the symmetrized Kullback-Leibler divergency for two objects 
     *  @see https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence
     *  \f[ f(v_1, C_1, v_2, C_2) = 
     *  (v_1-v_2)^{T} \left  ( C_1^{-1} + C_2^{-1} \right) (v_1 - v_2)  
     *   + Sp \left ( C_1 - C_2 \right ) 
     *    \times  \left ( C_2^{-1} - C_1^{-1} \right ) \f]
     *  @param v1 the first  data vector 
     *  @param c1 the first  covariance matrix 
     *  @param v2 the second data vector 
     *  @param c2 the second covariance matrix 
     *  @return Symmetrised Kullback-Leibler divergency, or -999 
     */
    template <unsigned int N, typename SCALAR>
    inline double kullback_leibler 
    ( const ROOT::Math::SVector<SCALAR,N>&                                    v1 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c1 , 
      const ROOT::Math::SVector<SCALAR,N>&                                    v2 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c2 )
    {
      /// the actual type of covariance matrix
      typedef typename ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> > COV ;
      ///
      double bad = -999 ;
      /// try to invert matrices 
      COV g1 { c1 } ;
      if  ( !g1.InvertChol () ) { return bad ; }
      COV g2 { c2 } ;
      if  ( !g2.InvertChol () ) { return bad ; }
      ///
      return ROOT::Math::Similarity ( g1 + g2 , v1 - v2 )
        + trace ( ( c1 - c2 ) * ( g2 - g1 ) ) ;
    }
    // ========================================================================
    /** Get the symmetrized Kullback-Leibler divergency for two objects 
     *  @see https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence
     *  \f[ f(v_1, C_1, v_2, C_2) = 
     *  (v_1-v_2)^{T} \left  ( C_1^{-1} + C_2^{-1} \right) (v_1 - v_2)  
     *   + Sp \left ( C_1 - C_2 \right ) 
     *    \times  \left ( C_2^{-1} - C_1^{-1} \right ) \f] 
     *  @param v1 the first  data vector 
     *  @param v2 the second data vector 
     *  @param c1 the first  covariance matrix 
     *  @param c2 the second covariance matrix 
     *  @return Symmetrised Kullback =Leibler divergency, or -999 
     */
    template <unsigned int N, typename SCALAR>
    inline double kullback_leibler 
    ( const ROOT::Math::SVector<SCALAR,N>&                                    v1 , 
      const ROOT::Math::SVector<SCALAR,N>&                                    v2 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c1 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c2 )
    { return kullback_leibler ( v1 , c1 , v2 , c2 ) ; }
    // ========================================================================
    /** Get the asymmetrical  Kullback-Leibler divergency for two objects 
     *  @see https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence
     *  @param v0 the first  data vector 
     *  @param c0 the first  covariance matrix 
     *  @param v1 the second data vector 
     *  @param c1 the second covariance matrix 
     *  @return asymmetric Kullback-Leibler divergency, or -999 
     */
    template <unsigned int N, typename SCALAR>
    inline double 
    asymmetric_kullback_leibler
    ( const ROOT::Math::SVector<SCALAR,N>&                                    v0 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c0 , 
      const ROOT::Math::SVector<SCALAR,N>&                                    v1 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c1 )
    {
      /// the actual type of covariance matrix
      typedef typename ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> > COV ;
      ///
      double bad = -999 ;
      ///
      SCALAR det0 = 1 ;
      if ( !c0.Det2 ( det0 ) ) { return bad ; }
      SCALAR det1 = 1 ;
      if ( !c1.Det2 ( det1 ) ) { return bad ; }
      //
      /// try to invert matrices 
      COV g1 { c1 } ;
      if  ( !g1.InvertChol () ) { return bad ; }
      ///
      return ( trace ( g1 * c0 ) - N 
               + ROOT::Math::Similarity ( g1 , v0 - v1 ) , 
               + std::log ( det1 / det0 ) ) / 2 ;
    }
    // ========================================================================
    /** Get the asymmetrical  Kullback-Leibler divergency for two objects 
     *  @see https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence
     *  @param v0 the first  data vector 
     *  @param c0 the first  covariance matrix 
     *  @param v1 the second data vector 
     *  @param c1 the second covariance matrix 
     *  @return asymmetric Kullback-Leibler divergency, or -999 
     */
    template <unsigned int N, typename SCALAR>
    inline double 
    asymmetric_kullback_leibler
    ( const ROOT::Math::SVector<SCALAR,N>&                                    v0 , 
      const ROOT::Math::SVector<SCALAR,N>&                                    v1 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c0 , 
      const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& c1 )
    { return asymmetric_kullback_leibkler ( v0 , c0 , v1 , c1 ) ; }
    // ========================================================================
    /*  get Cholesky decomposition for the covarance matrix 
     *  @param M (INPUT)  input symmetric positive definite matrix 
     *  @param L (OUTPUT) Cholesky decomposition of the covariance matrix 
     *  @return true if decomposition OK (matrix is positive definite) else false 
     */
    // ============================================================================
    template <unsigned int N, class SCALAR>
    inline bool 
    cholesky 
    ( const ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepSym<SCALAR,N> >& M , 
      ROOT::Math::SMatrix<SCALAR,N,N,ROOT::Math::MatRepStd<SCALAR,N,N> >    & L )
    {
      const ROOT::Math::CholeskyDecomp<SCALAR,N> decomp ( M ) ;
      return decomp.getL ( L ) ;
    }
    // ========================================================================

    // helper functions to allow proper operations in PyROOT
    // - need to avoid expressions  (no easy way to use them in PyROOT)
    // ========================================================================
    
    // ========================================================================
    // Vector operations 
    // ========================================================================
    template <class VECTOR>
    struct VctrOps ;
    // ========================================================================
    template <class T, unsigned int D>
    struct VctrOps < ROOT::Math::SVector<T,D> >
    {
      // a + b 
      static
      ROOT::Math::SVector<T,D> 
      add  ( const ROOT::Math::SVector<T,D>& a ,
             const ROOT::Math::SVector<T,D>& b ) { return a + b ; }
      // a + c 
      static
      ROOT::Math::SVector<T,D>
      add  ( const ROOT::Math::SVector<T,D>& a ,
             const double                    c ) { return a + c ; }
      // a - b 
      static
      ROOT::Math::SVector<T,D> 
      sub  ( const ROOT::Math::SVector<T,D>& a ,
             const ROOT::Math::SVector<T,D>& b ) { return a - b ; }
      // a - c 
      static
      ROOT::Math::SVector<T,D> 
      sub  ( const ROOT::Math::SVector<T,D>& a ,
             const double                    b ) { return a + b ; }
      // c - a  
      static
      ROOT::Math::SVector<T,D> 
      rsub ( const ROOT::Math::SVector<T,D>& a ,
             const double                    c ) { return c - a  ; }
    };    
    // ========================================================================
    // matrix operations 
    // ========================================================================
    template <class VECTOR>
    struct MtrxOps ;
    // generic matrices 
    template <class T, unsigned int D1, unsigned int D2>
    struct MtrxOps< ROOT::Math::SMatrix<T,D1,D2,ROOT::Math::MatRepStd<T,D1,D2> > > 
    {
      typedef ROOT::Math::SMatrix<T,D1,D2,ROOT::Math::MatRepStd<T,D1,D2> > MATRIX ;
      //
      // ======================================================================
      // m + m 
      static MATRIX add  ( const MATRIX& a ,
                           const MATRIX& b ) { return a + b ; }
      // m + c 
      static MATRIX add  ( const MATRIX& a ,
                           const double  c ) { return a + c ; }
      // ======================================================================
      // m - m 
      static MATRIX sub  ( const MATRIX& a ,
                           const MATRIX& b ) { return a - b ; }
      // m - c
      static MATRIX sub  ( const MATRIX& a ,
                           const double  c ) { return a - c ; }
      // ======================================================================
      // c - m
      static MATRIX rsub ( const MATRIX& a ,
                           const double  c ) { return c - a ; }
      // ======================================================================
    };  
    // square matrices 
    template <class T, unsigned int D>
    struct MtrxOps<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepStd<T,D,D> > > 
    {
      typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepStd<T,D,D> >      MATRIX ;
      typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >     SYMMATRIX ;
      //
      // ======================================================================
      // m + m
      static MATRIX add  ( const    MATRIX& a ,
                           const    MATRIX& b ) { return a + b ; }
      // m + s 
      static MATRIX add  ( const    MATRIX& a ,
                           const SYMMATRIX& b ) { return a + b ; }
      // m + c
      static MATRIX add  ( const    MATRIX& a ,
                           const    double  b ) { return a + b ; }
      // ======================================================================
      // m - m 
      static MATRIX sub  ( const    MATRIX& a ,
                           const    MATRIX& b ) { return a - b ; }
      // m - s 
      static MATRIX sub  ( const    MATRIX& a ,
                           const SYMMATRIX& b ) { return a - b ; }
      // m - c 
      static MATRIX sub  ( const    MATRIX& a ,
                           const    double  c ) { return a - c ; }
      // ======================================================================
      // s - m 
      static MATRIX rsub ( const    MATRIX& a ,
                           const SYMMATRIX& b ) { return b - a ; }
      // c - m  
      static MATRIX rsub ( const    MATRIX& a ,
                           const    double  c ) { return c - a ; }
      
      // ======================================================================
    };
    // symmetric matrices
    template <class T, unsigned int D>
    struct MtrxOps<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > > 
    {
      typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >      MATRIX ;
      typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepStd<T,D,D> > GENMATRIX ;
      //
      // ======================================================================
      // s + s 
      static MATRIX    add  ( const    MATRIX& a ,
                              const    MATRIX& b ) { return a + b ; }
      // s + m      
      static GENMATRIX add  ( const    MATRIX& a ,
                              const GENMATRIX& b ) { return a + b ; }
      // s + c
      static MATRIX    add  ( const    MATRIX& a ,
                              const     double c ) { return a + c ; }      
      // ======================================================================
      // s - s 
      static MATRIX    sub  ( const    MATRIX& a ,
                              const    MATRIX& b ) { return a - b ; }
      // s - m      
      static GENMATRIX sub  ( const    MATRIX& a ,
                              const GENMATRIX& b ) { return a - b ; }
      // s - c
      static MATRIX    sub  ( const    MATRIX& a ,
                              const     double c ) { return a - c ; }
      // ======================================================================      
      // m - s      
      static GENMATRIX rsub ( const    MATRIX& a ,
                              const GENMATRIX& b ) { return b - a ; }
      // c - s 
      static MATRIX    rsub ( const    MATRIX& a ,
                              const     double c ) { return c - a ; }
      // ======================================================================      
    };
    // ========================================================================
    /// multiplication
    template <class OBJ1, class OBJ2>
    struct MultiplyOp ;
    // ========================================================================
    // vector * vector  
    template <class T, unsigned int D>
    struct MultiplyOp < ROOT::Math::SVector<T,D>, ROOT::Math::SVector<T,D> >
    {
      // dot:
      static
      double 
      dot ( const ROOT::Math::SVector<T,D> & a , 
            const ROOT::Math::SVector<T,D> & b ) 
      {
        double result = 0 ;
        for ( unsigned short i = 0 ; i < D ; ++i ) { result +=  double( a[i] ) * b[i] ; }
        return result ;
      }
      // cross: 
      static
      ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepStd<T,D,D> >
      cross ( const ROOT::Math::SVector<T,D> & a , 
              const ROOT::Math::SVector<T,D> & b ) 
      {
        ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepStd<T,D,D> > result ;
        for ( unsigned short i = 0 ; i < D ; ++i ) 
        { for ( unsigned short j = 0 ; j < D ; ++j ) 
          { result(i,j) = a[i] * b[j] ; } }
        return result ;
      }
      // multiply
      static 
      double 
      multiply ( const ROOT::Math::SVector<T,D> & a , 
                 const ROOT::Math::SVector<T,D> & b ) { return dot ( a , b ) ; }
    } ;
    // ========================================================================    
    // cross: vector * vector  
    template <class T, unsigned int D1, unsigned int D2>
    struct MultiplyOp < ROOT::Math::SVector<T,D1>, ROOT::Math::SVector<T,D2> >
    {
      //
      static
      ROOT::Math::SMatrix<T,D1,D2,ROOT::Math::MatRepStd<T,D1,D2> >
      cross ( const ROOT::Math::SVector<T,D1> & a , 
              const ROOT::Math::SVector<T,D2> & b ) 
      {
        ROOT::Math::SMatrix<T,D1,D2,ROOT::Math::MatRepStd<T,D1,D2> > result ;
        for ( unsigned short i = 0 ; i < D1 ; ++i ) 
        { for ( unsigned short j = 0 ; j < D2 ; ++j ) 
          { result(i,j) = a[i] * b[j] ; } }
        return result ;
      }
    } ;
    // ========================================================================
    // vector * matrix 
    template <class T, unsigned int D, unsigned D2, class  R>
    struct MultiplyOp<ROOT::Math::SVector<T,D>,ROOT::Math::SMatrix<T,D,D2,R> >
    {
      static 
      ROOT::Math::SVector<T,D2> 
      multiply ( const ROOT::Math::SVector<T,D>      & a , 
                 const ROOT::Math::SMatrix<T,D,D2,R> & b ) { return a * b ; }
    } ;
    // =======================================================================
    // matrix * matrix  
    template <class T, unsigned int D1, unsigned D2, unsigned D3, class R1, class R2>
    struct MultiplyOp<ROOT::Math::SMatrix<T,D1,D2,R1>,ROOT::Math::SMatrix<T,D2,D3,R2> > 
    {
      static 
      ROOT::Math::SMatrix<T,D1,D3,ROOT::Math::MatRepStd<T,D1,D3> >  
      multiply  ( const ROOT::Math::SMatrix<T,D1,D2,R1> & a , 
                  const ROOT::Math::SMatrix<T,D2,D3,R2> & b ) { return a * b ; }
    } ;
    // =======================================================================
    // matrix * vector  
    template <class T, unsigned int D, unsigned D2, class  R>
    struct MultiplyOp <ROOT::Math::SMatrix<T,D,D2,R>,ROOT::Math::SVector<T,D2> >
    {
      static 
      ROOT::Math::SVector<T,D> 
      multiply  ( const ROOT::Math::SMatrix<T,D,D2,R> & a , 
                  const ROOT::Math::SVector<T,D2>     & b ) { return a * b ; }
    } ;
    // ========================================================================

    // ========================================================================
    template <class OBJ1, class OBJ2> 
    struct EqualityOp ;
    // ========================================================================
    // vector == vector 
    template <class T1,class T2, unsigned int D>
    struct EqualityOp< ROOT::Math::SVector<T1,D>,ROOT::Math::SVector<T2,D> >
    {
      static 
      bool 
      equal ( const ROOT::Math::SVector<T1,D>& v1 , 
              const ROOT::Math::SVector<T2,D>& v2 ) 
      {
        static const Equal_To<ROOT::Math::SVector<T1,D> > m_cmp{} ;
        return m_cmp ( v1 , v2 ) ;
      }
    } ;
    // ========================================================================
    // matrix == matrix 
    template <class T1,class T2, unsigned int D1,unsigned int D2, class R1, class R2>
    struct EqualityOp< ROOT::Math::SMatrix<T1,D1,D2,R1>,ROOT::Math::SMatrix<T2,D1,D2,R2> >
    {
      static 
      bool 
      equal ( const ROOT::Math::SMatrix<T1,D1,D2,R1>& v1 , 
              const ROOT::Math::SMatrix<T2,D1,D2,R2>& v2 ) 
      {
        static const Equal_To<ROOT::Math::SMatrix<T1,D1,D2,R1> > m_cmp{} ;
        return m_cmp ( v1 , v2 ) ;
      }
    } ;
    // ========================================================================
  } //                                                    end of namespace Math
  // ==========================================================================
} //                                               end of namespace Ostap::Math
// ============================================================================
// The END 
// ============================================================================
#endif // LHCBMATH_MATRIXUTILS_H
// ============================================================================
