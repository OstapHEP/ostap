#ifndef OSTAP_MATRIXUTILST_H 
#define OSTAP_MATRIXUTILST_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL 
// ============================================================================
#include <algorithm>
#include <numeric>
#include <functional>
#include <utility>
#include <cmath>
#include <array>
// ============================================================================
// ROOT
// ============================================================================
#include "TVectorT.h"
#include "TMatrixT.h"
#include "TMatrixTSym.h"
#include "TMatrixTUtils.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/MatrixUtils2.h"
// ============================================================================
/** @file Ostap/MatrixUtils2.h
 *  The collection of functions for manipulation with matrices and vectors.
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 */
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    /// Are all elements are finite? 
    template <class T> 
    inline bool isfinite 
    ( const TVectorT<T>& vct )
    {
      if ( !vct.IsValid() ) { return false ; }
      const T* begin =         vct.GetMatrixArray () ;
      const T* end   = begin + vct.GetNrows       () ;
      for ( const T* v = begin ; v != end ; ++v )
	{ if ( !std::isfinite ( *v ) ) { return false ; } }
      return true ;
    }
    /// Are all elements are finite? 
    template <class T> 
    inline bool isfinite 
    ( const TMatrixT<T>& mtrx )
    {
      if ( !mtrx.IsValid() ) { return false ; }
      const T* begin =         mtrx.GetMatrixArray () ;
      const T* end   = begin + mtrx.GetNrows() * mtrx.GetNcols() ;
      for ( const T* v = begin ; v != end ; ++v )
	{ if ( !std::isfinite ( *v ) ) { return false ; } }
      return true ;
    }
    /// Are all elements are finite? 
    template <class T> 
    inline bool isfinite 
    ( const TMatrixTSym<T>& mtrx )
    {
      if ( !mtrx.IsValid() ) { return false ; }
      const T* begin =         mtrx.GetMatrixArray () ;
      const T* end   = begin + mtrx.GetNrows() * mtrx.GetNcols() ;
      for ( const T* v = begin ; v != end ; ++v )
	{ if ( !std::isfinite ( *v ) ) { return false ; } }
      return true ;
    }
    // ========================================================================
    /// specialisation for vectors  
    template <class T>
    struct Equal_To<TVectorT<T> > 
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
      ( const TVectorT<T>& v1 , 
        const TVectorT<T>& v2 ) const  
      {
        return &v1 == &v2 ||
          ( v1.IsValid  () && v2.IsValid  () &&
            v1.GetNrows () == v2.GetNrows () && 
            std::equal ( v1.GetMatrixArray () , v1.GetMatrixArray () + v1.GetNrows() ,
                         v2.GetMatrixArray () , m_cmp ) ) ;                 
      }
      /// compare with another vector type (e.g. double and float)
      template <class T1, class T2>
      inline bool operator()
      ( const TVectorT<T1>&  v1 , 
        const TVectorT<T2>& v2 ) const
      {
        return ( v1.IsValid  () && v2.IsValid  () &&
                 v1.GetNrows () == v2.GetNrows () && 
                 std::equal ( v1.GetMatrixArray () , v1.GetMatrixArray () + v1.GetNrows() ,
                              v2.GetMatrixArray () , m_cmp ) ) ;                 
      }
      // ======================================================================
      template <class T1, class T2, unsigned int D>
      inline bool operator ()
      ( const TVectorT<T1>&               v1 ,
        const ROOT::Math::SVector<T2,D>&  v2 ) const
      {
        return v1.IsValid() && D == v1.GetNrows() &&
          std::equal ( v2.begin() , v2.end()  , v1.GetMatrixArray () , m_cmp );
      }
      // ======================================================================
      template <class T1, class T2, unsigned int D>
      inline bool operator ()
      ( const ROOT::Math::SVector<T2,D>& v1 , 
        const TVectorT<T1>&              v2 ) const
      { return (*this)( v2 , v1 ) ; }
      // ======================================================================
    private:
      // ======================================================================
      /// the evaluator 
      Equal_To<T> m_cmp ;                                 // the evaluator 
      // ======================================================================
    } ;
    // ========================================================================
    /// specialisation for matrices  
    template <class T>
    struct Equal_To<TMatrixT<T> > 
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
      ( const TMatrixT<T>& v1 , 
        const TMatrixT<T>& v2 ) const  
      {
        return &v1 == &v2 ||
          ( v1.IsValid  () && v2.IsValid  () &&
            v1.GetNrows () == v2.GetNrows () && 
            v1.GetNcols () == v2.GetNcols () &&
            std::equal ( v1.GetMatrixArray () ,
                         v1.GetMatrixArray () + v1.GetNrows()* v1.GetNcols()  ,
                         v2.GetMatrixArray () , m_cmp ) ) ;                 
      }
      // ======================================================================
      /// compare with another matrix type 
      template <class T1, class T2> 
      inline bool operator()
      ( const TMatrixT<T1>& v1 , 
        const TMatrixT<T2>& v2 ) const  
      {
        return ( v1.IsValid  () && v2.IsValid  () &&
                 v1.GetNrows () == v2.GetNrows () && 
                 v1.GetNcols () == v2.GetNcols () &&
                 std::equal ( v1.GetMatrixArray () ,
                              v1.GetMatrixArray () + v1.GetNrows()* v1.GetNcols()  ,
                              v2.GetMatrixArray () , m_cmp ) ) ;                 
      }
      // ======================================================================
      /// compare with another matrix type 
      template <class T1, class T2>
      inline bool operator()
      ( const TMatrixT<T1>&    v1 , 
        const TMatrixTSym<T2>& v2 ) const
      {
        if ( !v1.IsValid  () || !v2.IsValid  () ) { return false ; }
        if (  v1.GetNrows () !=  v2.GetNrows () ) { return false ; }
        if (  v1.GetNcols () !=  v2.GetNcols () ) { return false ; }
        //
        const unsigned long nr = v1.GetNrows () ;
        const unsigned long nc = v1.GetNcols () ;
        //
        for ( unsigned long i = 0 ; i < nr ; ++i )
        { for ( unsigned long j = 0 ; j < nc ; ++j )
          { if ( !m_cmp ( v1 ( i , j ) , v2 ( i , j ) ) )
            { return false ; } } }
        return true ;
      }
      // ======================================================================
      /// compare with another matrix type 
      template <class T1, class T2>      
      inline bool operator()
      ( const TMatrixTSym<T2>&  v1 , 
        const TMatrixT<T1>&     v2 ) const
      { return  (*this) ( v2 , v1 ) ; }
      // ======================================================================
      /// compare with another matrix type
      template <class T1, class T2, unsigned int D1 , unsigned int D2, class R1>
      inline bool operator()
      ( const TMatrixT<T1>&                     v1 , 
        const ROOT::Math::SMatrix<T2,D1,D2,R1>& v2 ) const
      {
        if ( !v1.IsValid() || v1.GetNrows() != D1 || v1.GetNcols() != D2 ) 
        { return false ; }
        //
        for ( unsigned long i = 0 ; i < D1 ; ++i )
        { for ( unsigned long j = 0 ; j < D2 ; ++j )
          { if ( !m_cmp ( v1 ( i , j ) , v2  ( i , j ) ) ) { return false ; } } }
        //
        return true ;
      }
      // ======================================================================
      /// compare with another matrix type
      template <class T1, class T2, unsigned int D1 , unsigned int D2, class R1>
      inline bool operator()
      ( const ROOT::Math::SMatrix<T2,D1,D2,R1>& v1 , 
        const TMatrixT<T1>&                     v2 ) const 
      { return  (*this) ( v2 , v1 ) ; }
      // ======================================================================
    private:
      // ======================================================================
      /// the evaluator 
      Equal_To<T> m_cmp ;                                 // the evaluator 
      // ======================================================================
    } ;
    // ========================================================================
    /// specialisation for matrices  
    template <class T>
    struct Equal_To<TMatrixTSym<T> > 
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
      ( const TMatrixTSym<T>& v1 , 
        const TMatrixTSym<T>& v2 ) const  
      {
        if ( &v1             ==  &v2            ) { return true  ; }
        if ( !v1.IsValid  () || !v2.IsValid  () ) { return false ; }
        if (  v1.GetNrows () !=  v2.GetNrows () ) { return false ; }
        //
        const unsigned long nc = v1.GetNcols () ;
        //
        for ( unsigned long i = 0 ; i < nc ; ++i )
        { for ( unsigned long j = i ; j < nc ; ++j ) // ATTENTION!!!
          { if ( !m_cmp ( v1 ( i , j ) , v2  ( i , j ) ) ) { return false ; } } }
        //
        return true ;
      }
      // ======================================================================
      /// compare with another matrix type (e.g. double and float
      template <class T1, class T2>
      inline bool operator()
      ( const TMatrixTSym<T1>& v1 , 
        const TMatrixTSym<T2>& v2 ) const  
      {
        if ( !v1.IsValid  () || !v2.IsValid  () ) { return false ; }
        if (  v1.GetNrows () !=  v2.GetNrows () ) { return false ; }
        //
        const unsigned long nc = v1.GetNcols () ;
        //
        for ( unsigned long i = 0 ; i < nc ; ++i )
        { for ( unsigned long j = i ; j < nc ; ++j ) // ATTENTION!!!
          { if ( !m_cmp ( v1 ( i , j ) , v2  ( i , j ) ) ) { return false ; } } }
        //
        return true ;
      }
      // ======================================================================
      /// compare with another matrix type
      template <class T1, class T2>
      inline bool operator()
      ( const TMatrixTSym<T1>& v1 , 
        const TMatrixT<T2>&    v2 ) const
      {
        if ( !v1.IsValid  () || !v2.IsValid  () ) { return false ; }
        if (  v1.GetNrows () !=  v2.GetNrows () ) { return false ; }
        if (  v1.GetNcols () !=  v2.GetNcols () ) { return false ; }
        //
        const unsigned long nc = v1.GetNcols () ;
        const unsigned long nr = v1.GetNrows () ;
        //
        for ( unsigned long i = 0 ; i < nr ; ++i )
        { for ( unsigned long j = 0 ; j < nc ; ++j )
          { if ( !m_cmp ( v1 ( i , j ) , v2  ( i , j ) ) ) { return false ; } } }
        //
        return true ;
      }
      // ======================================================================
      /// compare with another matrix type
      template <class T1, class T2>
      inline bool operator()
      ( const TMatrixT<T2>&    v1 ,
        const TMatrixTSym<T1>& v2 ) const  
      {
        if ( !v1.IsValid  () || !v2.IsValid  () ) { return false ; }
        if (  v1.GetNrows () !=  v2.GetNrows () ) { return false ; }
        if (  v1.GetNcols () !=  v2.GetNcols () ) { return false ; }
        //
        const unsigned long nc = v1.GetNcols () ;
        const unsigned long nr = v1.GetNrows () ;
        //
        for ( unsigned long i = 0 ; i < nr ; ++i )
        { for ( unsigned long j = 0 ; j < nc ; ++j )
          { if ( !m_cmp ( v1 ( i , j ) , v2  ( i , j ) ) ) { return false ; } } }
        //
        return true ;
      }
      // ======================================================================
      /// compare with another matrix type
      template <class T1, class T2, unsigned int D, class R1>
      inline bool operator()
      ( const TMatrixTSym<T1>&                v1 , 
        const ROOT::Math::SMatrix<T2,D,D,R1>& v2 ) const
      {
        if ( !v1.IsValid() || v1.GetNrows() != D || v1.GetNcols() != D ) { return false ; }
        //
        for ( unsigned long i = 0 ; i < D ; ++i )
        { for ( unsigned long j = 0 ; j < D ; ++j )
          { if ( !m_cmp ( v1 ( i , j ) , v2  ( i , j ) ) ) { return false ; } } }
        //
        return true ;
      }
      // ==========================================================================
      /// compare with another matrix type
      template <class T1, class T2, unsigned int D>
      inline bool operator()
      ( const TMatrixTSym<T1>&                                          v1 , 
        const ROOT::Math::SMatrix<T2,D,D,ROOT::Math::MatRepSym<T2,D> >& v2 ) const
      {
        if ( !v1.IsValid() || v1.GetNrows() != D || v1.GetNcols() != D ) { return false ; }
        //
        for ( unsigned long i = 0 ; i < D ; ++i )
        { for ( unsigned long j = i  ; j < D ; ++j ) // ATTENTION HERE 
          { if ( !m_cmp ( v1 ( i , j ) , v2  ( i , j ) ) ) { return false ; } } }
        //
        return true ;
      }
      // ======================================================================
      /// compare with another matrix type
      template <class T1, class T2, unsigned int D, class R1>
      inline bool operator()
      ( const ROOT::Math::SMatrix<T2,D,D,R1>& v1 , 
        const TMatrixTSym<T1>&                v2 ) const 
      { return  (*this) ( v2 , v1 ) ; }
      // ======================================================================
    private:
      // ======================================================================
      /// the evaluator 
      Equal_To<T> m_cmp ;                                 // the evaluator 
      // ======================================================================
    } ;
    // ========================================================================
    /// get element with maximal absolute value        
    template <class T>
    inline T 
    maxabs_element 
    ( const TMatrixT<T>&    m )
    {
      if ( !m.IsValid() ) { return 0 ; }
      T result = m ( 0 , 0 ) ;
      const unsigned int rows = m.GetNrows () ;
      const unsigned int cols = m.GetNcols () ;
      for ( unsigned int i = 0 ; i < rows ; ++i )
        {
          for ( unsigned int j = 0 ; j < cols ; ++j )
            {
              const double value = m ( i , j ) ;
              if ( std::abs ( result ) < std::abs ( value ) ) { result = value ; }
            }
        }
      return result ;
    }
    // ========================================================================
    /// get element with the maximal absolute value        
    template <class T>
    inline T 
    maxabs_element 
    ( const TMatrixTSym<T>& m )
    {
      if ( !m.IsValid() ) { return 0 ; }
      T result = m ( 0 , 0 ) ;
      const unsigned int rows = m.GetNrows () ;
      const unsigned int cols = m.GetNcols () ;
      for ( unsigned int i = 0 ; i < rows ; ++i )
        {
          for ( unsigned int j = i ; j < cols ; ++j )
            {
              const double value = m ( i , j ) ;
              if ( std::abs ( result ) < std::abs ( value ) ) { result = value ; }
            }
        }
      return result ;
    }
    // ========================================================================
    /// get element with the maximal absolute value        
    template <class T>
    inline T 
    maxabs_diagonal
    ( const TMatrixT<T>&    m )
    {
      if ( !m.IsValid() ) { return 0 ; }
      T result = m ( 0 , 0 ) ;
      const unsigned int rows = m.GetNrows () ;
      const unsigned int cols = m.GetNcols () ;
      const unsigned int D    = std::min ( rows , cols ) ; 
      for ( unsigned int i = 0 ; i < D ; ++i )
        {
          const double value = m ( i , i ) ;
          if ( std::abs ( result ) < std::abs ( value ) ) { result = value ; }
        }
      return result ;
    }
    // ========================================================================
    /// get diagonal element with the maximal absolute value        
    template <class T>
    inline T 
    maxabs_diagonal 
    ( const TMatrixTSym<T>& m )
    {
      if ( !m.IsValid() ) { return 0 ; }
      T result = m ( 0 , 0 ) ;
      const unsigned int rows = m.GetNrows () ;
      const unsigned int cols = m.GetNcols () ;
      const unsigned int D    = std::min ( rows , cols ) ; 
      for ( unsigned int i = 0 ; i < D ; ++i )
        {
          const double value = m ( i , i ) ;
          if ( std::abs ( result ) < std::abs ( value ) ) { result = value ; }
        }
      return result ;
    }
    // ========================================================================
    /// get element with maximal absolute value        
    template <class T>
    inline T 
    maxabs_element 
    ( const TVectorT<T>& v  )
    {
      if ( !v.IsValid() ) { return 0 ; }
      T result = v ( 0 ) ;
      const unsigned int rows = v.GetNrows () ;
      for ( unsigned int i = 1 ; i < rows ; ++i )
        {
          const double value = v ( i ) ;
          if ( std::abs ( result ) < std::abs ( value ) ) { result = value ; }
        }
      return result ;
    }
    // =======================================================================
    namespace  Ops
    {      
      // ======================================================================
      // converters 
      // ======================================================================
      template <class M> struct TM ;
      template <class M> struct SM ;
      // ======================================================================
      //  S --> T 
      // ======================================================================
      template <class T, unsigned int D1, unsigned int D2, class RR>
      struct TM<ROOT::Math::SMatrix<T,D1,D2,RR> >
      {
        typedef TMatrixT<T>  R ;
        // 
        static R operation  
        ( const ROOT::Math::SMatrix<T,D1,D2,RR>& m )
        {
          R result ( D1 , D2 ) ;
          for ( unsigned int i = 0 ; i < D1 ; ++i )
          { for ( unsigned int j = 0 ; j < D1 ; ++j )
            { result ( i , j ) = m ( i , j ) ; } }
          return result ;
        }
      } ;
      // ======================================================================
      template <class T, unsigned int D>
      struct TM<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > >
      {
        typedef TMatrixTSym<T>  R ;
        //
        static R operation
        ( const ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >& m )
        {
          R result ( D ) ;
          for ( unsigned int i = 0 ; i < D ; ++i )
          { for ( unsigned int j = 0 ; j < D ; ++j )
            { result ( i , j ) = m ( i , j ) ; } }
          return result ;
        }
      } ;
      // ======================================================================
      template <class T, unsigned int D>
      struct TM<ROOT::Math::SVector<T,D> >
      {
        typedef TVectorT<T>  R ;
        //
        static R operation 
        ( const ROOT::Math::SVector<T,D>& m )
        { return R ( D , m.begin() ) ; }
      } ;
      // ======================================================================
      
      // ======================================================================
      // T -> S
      // ======================================================================
      template <class T, unsigned int D1, unsigned int D2>
      struct SM<ROOT::Math::SMatrix<T,D1,D2> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2> R ;
        //
        static R operation 
        ( const TMatrixT<T>& m )
        {
          const T* start = m.GetMatrixArray() ;
          return R ( start , start + D1 * D2 ) ;
        }        
      } ;
      // ======================================================================
      template <class T, unsigned int D>
      struct SM<ROOT::Math::SMatrix<T,D,D> >
      {
        typedef ROOT::Math::SMatrix<T,D,D> R ;
        //
        static R operation 
        ( const TMatrixTSym<T>& m )
        {
          R result ;
          for ( unsigned int i = 0 ; i < D ; ++i  )
          { for ( unsigned int j = 0 ; j < D ; ++j )
            { result ( i , j ) = m ( i , j ) ; } }
          return result ;
        } 
        static R operation 
        ( const TMatrixT<T>& m )
        {
          const T* start = m.GetMatrixArray() ;
          return R ( start , start + D * D ) ;
        }        
      } ;  
      // ======================================================================      
      template <class T, unsigned int D>
      struct SM<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D>> >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > R ;
        //
        static R operation 
          ( const TMatrixTSym<T>& m )
        {
          R result ;
          for ( unsigned int i = 0 ; i < D ; ++i  )
          { for ( unsigned int j = i ; j < D ; ++j )
            { result ( i , j ) = 0.5 * ( m( i , j ) + m ( j, i ) ) ; } }
          return result ;
        } 
      } ;
      // ======================================================================
      template <class T, unsigned int D>
      struct SM<ROOT::Math::SVector<T,D> >
      {
        typedef ROOT::Math::SVector<T,D> R ;
        //
        static R operation
        ( const TVectorT<T>& m )
        {
          const T* start = m.GetMatrixArray() ;
          return R ( start , start + D ) ;
        }
      } ;
      // ======================================================================

      // ======================================================================
      // CHECKERS 
      // ======================================================================
      template <class T>
      struct CanAdd<TMatrixT<T>,TMatrixT<T> >
      {
        static bool operation
        ( const TMatrixT<T>& m1 , 
          const TMatrixT<T>& m2 )
        {
          return
            m1.IsValid  () && m2.IsValid  () &&
            m1.GetNrows () == m2.GetNrows () &&
            m1.GetNcols () == m2.GetNcols () ;
        }  
      } ;
      // ======================================================================
      template <class T>
      struct CanAdd<TMatrixT<T>,TMatrixTSym<T> >
      {
        static bool operation
        ( const TMatrixT<T>&    m1 , 
          const TMatrixTSym<T>& m2 )
        {
          return
            m1.IsValid  () && m2.IsValid  () &&
            m1.GetNrows () == m2.GetNrows () &&
            m1.GetNcols () == m2.GetNcols () ;
        }  
      } ;
      // ======================================================================
      template <class T>
      struct CanAdd<TMatrixTSym<T>,TMatrixT<T> >
      {
        static bool operation
        ( const TMatrixTSym<T>& m1 , 
          const TMatrixT<T>&    m2 )
        {
          return
            m1.IsValid  () && m2.IsValid  () &&
            m1.GetNrows () == m2.GetNrows () &&
            m1.GetNcols () == m2.GetNcols () ;
        }  
      } ;
      // ======================================================================
      template <class T>
      struct CanAdd<TMatrixTSym<T>,TMatrixTSym<T> >
      {
        static bool operation
        ( const TMatrixTSym<T>& m1 , 
          const TMatrixTSym<T>& m2 )
        {
          return
            m1.IsValid  () && m2.IsValid  () &&
            m1.GetNrows () == m2.GetNrows () &&
            m1.GetNcols () == m2.GetNcols () ;
        }  
      } ;
      // ======================================================================      
      template <class T,unsigned int D1, unsigned int D2, class R>
      struct CanAdd<ROOT::Math::SMatrix<T,D1,D2,R>,TMatrixT<T> >
      {
        static bool operation
        ( const ROOT::Math::SMatrix<T,D1,D2,R>& m1 , 
          const TMatrixT<T>&                    m2 )
        {
          return
            m2.IsValid        () &&
            D1 == m2.GetNrows () &&
            D2 == m2.GetNcols () ;
        }  
      } ;
      // ======================================================================
      template <class T,unsigned int D1, unsigned int D2, class R>
      struct CanAdd<TMatrixT<T>,ROOT::Math::SMatrix<T,D1,D2,R> >
      {
        static bool operation
        ( const TMatrixT<T>&                    m2 ,          
          const ROOT::Math::SMatrix<T,D1,D2,R>& m1 ) 
        {
          return
            m2.IsValid        () &&
            D1 == m2.GetNrows () &&
            D2 == m2.GetNcols () ;
        }  
      } ;
      // ======================================================================
      template <class T,unsigned int D1, unsigned int D2, class R>
      struct CanAdd<ROOT::Math::SMatrix<T,D1,D2,R>,TMatrixTSym<T> >
      {
        static bool operation
        ( const ROOT::Math::SMatrix<T,D1,D2,R>& m1 , 
          const TMatrixTSym<T>&                 m2 )
        {
          return
            m2.IsValid        () &&
            D1 == m2.GetNrows () &&
            D2 == m2.GetNcols () ;
        }  
      } ;
      // ======================================================================
      template <class T,unsigned int D1, unsigned int D2, class R>
      struct CanAdd<TMatrixTSym<T>,ROOT::Math::SMatrix<T,D1,D2,R> >
      {
        static bool operation
        ( const TMatrixTSym<T>&                 m2 ,          
          const ROOT::Math::SMatrix<T,D1,D2,R>& m1 ) 
        {
          return
            m2.IsValid        () &&
            D1 == m2.GetNrows () &&
            D2 == m2.GetNcols () ;
        }  
      } ;
      // ======================================================================
      template <class T>
      struct CanAdd<TVectorT<T>,TVectorT<T> >
      {
        static bool operation
        ( const TVectorT<T>& m1 ,  
          const TVectorT<T>& m2 ) 
        { return m1.IsValid() && m2.IsValid () && m1.GetNrows()== m2.GetNrows () ; }
      };
      // ======================================================================
      template <class T,unsigned int D>
      struct CanAdd<ROOT::Math::SVector<T,D>,TVectorT<T> >
      {
        static bool operation
        ( const ROOT::Math::SVector<T,D> & /* m1 */ ,  
          const TVectorT<T>&                  m2        ) 
        { return m2.IsValid () && D == m2.GetNrows () ; }
      };
      // ======================================================================
      template <class T,unsigned int D>
      struct CanAdd<TVectorT<T>,ROOT::Math::SVector<T,D> >
      {
        static bool operation
        ( const TVectorT<T>&                  m2     ,
          const ROOT::Math::SVector<T,D> & /* m1 */  )
        { return m2.IsValid () && D == m2.GetNrows () ; }
      };
      // ======================================================================
      
      // ======================================================================
      template <class T>
      struct CanMul<TMatrixT<T>,TMatrixT<T>>
      {
        static bool operation
        ( const TMatrixT<T>& m1 ,
          const TMatrixT<T>& m2 )
        { return m1.IsValid() &&  m2.IsValid() && m1.GetNcols() == m2.GetNrows() ; }
      };
      // ======================================================================
      template <class T>
      struct CanMul<TMatrixT<T>,TMatrixTSym<T>>
      {
        static bool operation
        ( const TMatrixT<T>   & m1 ,
          const TMatrixTSym<T>& m2 )
        { return m1.IsValid() &&  m2.IsValid() && m1.GetNcols() == m2.GetNrows() ; }
      };
      // ======================================================================      
      template <class T>
      struct CanMul<TMatrixTSym<T>,TMatrixT<T>>
      {
        static bool operation
        ( const TMatrixTSym<T> & m1 ,
          const TMatrixT<T>    & m2 )
        { return m1.IsValid() &&  m2.IsValid() && m1.GetNcols() == m2.GetNrows() ; }
      };
      // ======================================================================      
      template <class T>
      struct CanMul<TMatrixTSym<T>,TMatrixTSym<T>>
      {
        static bool operation
        ( const TMatrixTSym<T> & m1 ,
          const TMatrixTSym<T> & m2 )
        { return m1.IsValid() &&  m2.IsValid() && m1.GetNcols() == m2.GetNrows() ; }
      };
      // ======================================================================      
      template <class T>
      struct CanMul<TMatrixT<T>,TVectorT<T>>
      {
        static bool operation
        ( const TMatrixT<T>& m1 ,
          const TVectorT<T>& m2 )
        { return m1.IsValid() &&  m2.IsValid() && m1.GetNcols() == m2.GetNrows() ; }
      };
      // ======================================================================
      template <class T>
      struct CanMul<TMatrixTSym<T>,TVectorT<T>>
      {
        static bool operation
        ( const TMatrixTSym<T>& m1 ,
          const TVectorT<T>&    m2 )
        { return m1.IsValid() &&  m2.IsValid() && m1.GetNcols() == m2.GetNrows() ; }
      } ;
      // ======================================================================
      template <class T>
      struct CanMul<TVectorT<T>,TMatrixT<T>>
      {
        static bool operation
        ( const TVectorT<T>& m1 ,
          const TMatrixT<T>& m2 )
        { return m1.IsValid() &&  m2.IsValid() && m1.GetNrows () == m2.GetNrows() ; }
      };
      // ======================================================================
      template <class T>
      struct CanMul<TVectorT<T>,TMatrixTSym<T>>
      {
        static bool operation
        ( const TVectorT<T>&    m1 ,
          const TMatrixTSym<T>& m2 )
        { return m1.IsValid() &&  m2.IsValid() && m1.GetNrows () == m2.GetNrows() ; }
      };
      // ======================================================================
      
      template <class T>
      struct CanMul<TVectorT<T>,TVectorT<T>>
      {
        static bool operation
        ( const TVectorT<T>& m1 ,
          const TVectorT<T>& m2 )
        { return m1.IsValid() &&  m2.IsValid() && m1.GetNrows () == m2.GetNrows() ; }
      };
      // ======================================================================

      
      // ======================================================================
      template <class T,unsigned int D1, unsigned int D2, class R1>
      struct CanMul<ROOT::Math::SMatrix<T,D1,D2,R1>,TMatrixT<T>>
      {
        static bool operation
        ( const ROOT::Math::SMatrix<T,D1,D2,R1>& /*  m1 */ ,
          const TMatrixT<T>&                         m2 )
        { return m2.IsValid() &&  D2 == m2.GetNrows() ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D1, unsigned int D2, class R1>
      struct CanMul<ROOT::Math::SMatrix<T,D1,D2,R1>,TMatrixTSym<T>>
      {
        static bool operation
        ( const ROOT::Math::SMatrix<T,D1,D2,R1>& /*  m1 */ ,
          const TMatrixTSym<T>&                      m2 )
        { return m2.IsValid() &&  D2 == m2.GetNrows() ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D1, unsigned int D2, class R1>
      struct CanMul<ROOT::Math::SMatrix<T,D1,D2,R1>,TVectorT<T>>
      {
        static bool operation
        ( const ROOT::Math::SMatrix<T,D1,D2,R1>& /*  m1 */ ,
          const TVectorT<T>&                         m2 )
        { return m2.IsValid() &&  D2 == m2.GetNrows() ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct CanMul<ROOT::Math::SVector<T,D>,
                    TVectorT<T>>
      {
        static bool operation
        ( const ROOT::Math::SVector<T,D>& /*  m1 */ ,
          const TVectorT<T>&                 m2 )
        { return m2.IsValid() &&  D == m2.GetNrows() ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct CanMul<ROOT::Math::SVector<T,D>,
                    TMatrixT<T>>
      {
        static bool operation
          ( const ROOT::Math::SVector<T,D>& /*  m1 */ ,
            const TMatrixT<T>&                  m2 )
        { return m2.IsValid() &&  D == m2.GetNrows() ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct CanMul<ROOT::Math::SVector<T,D>,
                    TMatrixTSym<T>>
      {
        static bool operation
        ( const ROOT::Math::SVector<T,D>& /*  m1 */ ,
          const TMatrixTSym<T>&               m2 )
        { return m2.IsValid() &&  D == m2.GetNrows() ; }
      } ;

      // ======================================================================
      template <class T,unsigned int D1, unsigned int D2, class R1>
      struct CanMul<TMatrixT<T>,
                    ROOT::Math::SMatrix<T,D1,D2,R1> >
      {
        static bool operation
        ( const TMatrixT<T>&                         m2    ,
          const ROOT::Math::SMatrix<T,D1,D2,R1>& /*  m1 */ )
        { return m2.IsValid() &&  m2.GetNcols () == D1  ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D1, unsigned int D2, class R1>
      struct CanMul<TMatrixTSym<T>,
                    ROOT::Math::SMatrix<T,D1,D2,R1> >
      {
        static bool operation
        ( const TMatrixTSym<T>&                      m2    ,
          const ROOT::Math::SMatrix<T,D1,D2,R1>& /*  m1 */ )
        { return m2.IsValid() &&  m2.GetNcols () == D1  ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D1, unsigned int D2, class R1>
      struct CanMul<TVectorT<T>,
                    ROOT::Math::SMatrix<T,D1,D2,R1> >
      {
        static bool operation
        ( const TVectorT<T>&                         m2    ,
          const ROOT::Math::SMatrix<T,D1,D2,R1>& /*  m1 */ )
        { return m2.IsValid() &&  m2.GetNrows () == D1  ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct CanMul<TVectorT<T>,
                    ROOT::Math::SVector<T,D> >
      {
        static bool operation
        ( const TVectorT<T>&                  m2 , 
          const ROOT::Math::SVector<T,D>& /*  m1 */ ) 
        { return m2.IsValid() &&  m2.GetNrows () == D ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct CanMul<TMatrixT<T>,
                    ROOT::Math::SVector<T,D> >
      {
        static bool operation
        ( const TMatrixT<T>&                  m2    ,
          const ROOT::Math::SVector<T,D>& /*  m1 */ )
        { return m2.IsValid() && m2.GetNcols () == D ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct CanMul< TMatrixTSym<T>, 
                     ROOT::Math::SVector<T,D> > 
      {
        static bool operation
        ( const TMatrixTSym<T>&               m2    , 
          const ROOT::Math::SVector<T,D>& /*  m1 */ )
        { return m2.IsValid() && m2.GetNcols () == D ; }
      } ;

      // ======================================================================
      template <class T>
      struct CanIMul<TMatrixT<T>, TMatrixT<T> >
      {
        static bool operation
        ( const TMatrixT<T>& m1 ,
          const TMatrixT<T>& m2 )
        { return m1.IsValid() &&  m2.IsValid()
            && m1.GetNcols () == m2.GetNrows()
            && m2.GetNcols () == m2.GetNrows() ; }
      } ;      
      // ======================================================================
      template <class T>
      struct CanIMul<TMatrixT<T>, TMatrixTSym<T> >
      {
        static bool operation
        ( const TMatrixT<T>&    m1 ,
          const TMatrixTSym<T>& m2 )
        { return m1.IsValid() &&  m2.IsValid() 
            && m1.GetNcols () == m2.GetNrows() ; }
      } ;
      // ======================================================================
      template <class T, unsigned int D, class R1>
      struct CanIMul<TMatrixT<T>,
                     ROOT::Math::SMatrix<T,D,D,R1> >
      {
        static bool operation
        ( const TMatrixT<T>&                      m1    ,
          const ROOT::Math::SMatrix<T,D,D,R1>& /* m2 */ )
        { return m1.IsValid() && m1.GetNcols () == D ; }
      } ;
      // ======================================================================
      template <class T, unsigned int D1, unsigned D2>
      struct CanIMul<ROOT::Math::SMatrix<T,D1,D2>, 
                     TMatrixT<T> >
      {
        static bool operation
        ( const ROOT::Math::SMatrix<T,D1,D2>& /* m2 */ ,
          const TMatrixT<T>&                     m1    )
        { return m1.IsValid() && D2 == m1.GetNrows () && D2 == m1.GetNcols() ; }
      } ;
      // ======================================================================
      template <class T, unsigned int D1, unsigned D2>
      struct CanIMul<ROOT::Math::SMatrix<T,D1,D2>, 
                     TMatrixTSym<T> >
      {
        static bool operation
        ( const ROOT::Math::SMatrix<T,D1,D2>& /* m2 */ ,
          const TMatrixTSym<T>&                  m1    )
        { return m1.IsValid() && D2 == m1.GetNrows () ; }
      } ;
      // ======================================================================

      
      // add diagonal matrix 
      template <class T>
      struct CanAdd<TMatrixT<T> , double>
      {
        static bool operation
        ( const TMatrixT<T>&    m1    , 
          const double       /* m2 */ ) 
        { return m1.IsValid() && m1.GetNrows () == m1.GetNcols() ; }  
      } ;


      // ======================================================================
      template <class T>
      struct Add<TMatrixT<T> , double>
      {
        typedef TMatrixT<T>  M1 ;
        typedef TMatrixT<T>  R  ;
        //
        static R operation
        ( const M1& m1 , const double m2 ) 
        { 
          const Int_t D = m1.GetNrows() ;
          R result { D , D } ;
          for  ( unsigned int i = 0 ; i < D ; ++i ) { result ( i , i ) += m2 ; } 
          return result ;
        }
      } ;
      // =====================================================================
      template <class T>
      struct RAdd<TMatrixT<T>,double> 
        : public Add<TMatrixT<T> , double> {} ;      
      // ======================================================================      
      template <class T>
      struct IAdd<TMatrixT<T> , double>
      {
        typedef TMatrixT<T>  M1 ;
        static void operation
        ( M1& m1 , const double m2 ) 
        { 
          const Int_t D = m1.GetNrows() ;
          for  ( unsigned int i = 0 ; i < D ; ++i ) { m1 ( i , i ) += m2 ; } 
        }
      } ;
      // ======================================================================
      template <class T>
      struct Sub<TMatrixT<T> , double>
      {
        typedef TMatrixT<T>  M1 ;
        typedef TMatrixT<T>  R  ;
        //
        static R operation
        ( const M1& m1 , const double m2 ) 
        { 
          const Int_t D = m1.GetNrows() ;
          R result { D , D } ;
          for  ( Int_t i = 0 ; i < D ; ++i ) { result ( i , i ) -= m2 ; } 
          return result ;
        }
      } ;
      // =====================================================================
      template <class T>
      struct RSub<TMatrixT<T> , double>
      {
        typedef TMatrixT<T>  M1 ;
        typedef TMatrixT<T>  R  ;
        //
        static R operation
        ( const M1& m1 , const double m2 ) 
        { 
          const Int_t D = m1.GetNrows() ;
          R result { D , D } ; result *= -1 ;   // ATTENTION!!! 
          for  ( Int_t i = 0 ; i < D ; ++i ) { result ( i , i ) += m2 ; } 
          return result ;
        }
      } ;
      // ======================================================================
      template <class T>
      struct ISub<TMatrixT<T> , double>
      {
        typedef TMatrixT<T>  M1 ;
        static void operation
        ( M1& m1 , const double m2 ) 
        { 
          const Int_t D = m1.GetNrows() ;
          for  ( Int_t i = 0 ; i < D ; ++i ) { m1 ( i , i ) -= m2 ; } 
        }
      } ;


      // add diagonal matrix 
      template <class T>
      struct CanAdd<TMatrixTSym<T> , double>
      {
        static bool operation
        ( const TMatrixTSym<T>&    m1    , 
          const double          /* m2 */ ) { return m1.IsValid() ; }  
      } ;

      // ======================================================================
      template <class T>
      struct Add<TMatrixTSym<T> , double>
      {
        typedef TMatrixTSym<T>  M1 ;
        typedef TMatrixTSym<T>  R  ;
        //
        static R operation
        ( const M1& m1 , const double m2 ) 
        { 
          const Int_t D = m1.GetNrows() ;
          R result { D } ;
          for  ( Int_t i = 0 ; i < D ; ++i ) { result ( i , i ) += m2 ; } 
          return result ;
        }
      } ;
      // =====================================================================
      template <class T>
      struct RAdd<TMatrixTSym<T>,double> 
        : public Add<TMatrixTSym<T> , double> {} ;      
      // ======================================================================      
      template <class T>
      struct IAdd<TMatrixTSym<T> , double>
      {
        typedef TMatrixTSym<T>  M1 ;
        static void operation
        ( M1& m1 , const double m2 ) 
        { 
          const Int_t D = m1.GetNrows() ;
          for  ( Int_t  i = 0 ; i < D ; ++i ) { m1 ( i , i ) += m2 ; } 
        }
      } ;
      // ======================================================================
      template <class T>
      struct Sub<TMatrixTSym<T> , double>
      {
        typedef TMatrixTSym<T>  M1 ;
        typedef TMatrixTSym<T>  R  ;
        //
        static R operation
        ( const M1& m1 , const double m2 ) 
        { 
          const Int_t D = m1.GetNrows() ;
          R result { D } ;
          for  ( Int_t i = 0 ; i < D ; ++i ) { result ( i , i ) -= m2 ; } 
          return result ;
        }
      } ;
      // =====================================================================
      template <class T>
      struct RSub<TMatrixTSym<T> , double>
      {
        typedef TMatrixTSym<T>  M1 ;
        typedef TMatrixTSym<T>  R  ;
        //
        static R operation
        ( const M1& m1 , const double m2 ) 
        { 
          const Int_t D = m1.GetNrows() ;
          R result { D } ; result *= -1 ;   // ATTENTION!!! 
          for  ( Int_t i = 0 ; i < D ; ++i ) { result ( i , i ) += m2 ; } 
          return result ;
        }
      } ;
      // ======================================================================
      template <class T>
      struct ISub<TMatrixTSym<T> , double>
      {
        typedef TMatrixTSym<T>  M1 ;
        static void operation
        ( M1& m1 , const double m2 ) 
        { 
          const Int_t D = m1.GetNrows() ;
          for  ( Int_t i = 0 ; i < D ; ++i ) { m1 ( i , i ) -= m2 ; } 
        }
      } ;




      // =========================================================================      
      template <class T>
      struct CanPow<TMatrixT<T> >
      { 
        static bool operation 
        ( const TMatrixT<T>&      m1   ,
          const unsigned short /* n */ )
        { return m1.IsValid () && m1.GetNrows () == m1.GetNcols () ; } 
      } ; 
      
      template <class T>
      struct CanPow<TMatrixTSym<T> >
      { 
        static bool operation 
        ( const TMatrixTSym<T>&   m1   ,
          const unsigned short /* n */ )
        { return m1.IsValid () ; } 
      } ;

      template <class T>
      struct CanSym<TMatrixT<T> >
      { static bool operation ( const TMatrixT<T>& m1 )
        { return m1.IsValid () && m1.GetNrows () == m1.GetNcols () ; } } ;

      template <class T>
      struct CanSym<TMatrixTSym<T> >
      { static bool operation ( const TMatrixTSym<T>& m1 ) { return m1.IsValid () ; } } ;
      

      template <class T>
      struct CanASym<TMatrixT<T> >
      { static bool operation ( const TMatrixT<T>& m1 )
        { return m1.IsValid () && m1.GetNrows () == m1.GetNcols () ; } } ;

      template <class T>
      struct CanASym<TMatrixTSym<T> >
      { static bool operation ( const TMatrixTSym<T>& m1 ) { return m1.IsValid () ; } } ;
      
      
      // ======================================================================
      // add T + T 
      // ======================================================================
      template <class T>
      struct Add<TMatrixT<T>,TMatrixT<T> >
      {
        typedef TMatrixT<T> R ;
        //
        static  R operation 
        ( const TMatrixT<T>& m1 , 
          const TMatrixT<T>& m2 ) { return R ( m1 , R::kPlus ,  m2 ) ; }
      } ;
      // ======================================================================
      template <class T>
      struct Add<TMatrixT<T>,TMatrixTSym<T> >
      {
        typedef TMatrixT<T> R ;
        //
        static  R operation
        ( const TMatrixT<T>&    m1 , 
          const TMatrixTSym<T>& m2 ) { return R ( m1 , R::kPlus ,  m2 ) ; }
      } ;
      // ======================================================================
      template <class T>
      struct Add<TMatrixTSym<T>,TMatrixT<T> >
      {
        typedef TMatrixT<T> R ;
        //
        static  R operation 
        ( const TMatrixTSym<T>& m1 , 
          const TMatrixT<T>&    m2 ) { return R ( m2 , R::kPlus ,  m1 ) ; }
      } ;
      // ======================================================================
      template <class T>
      struct Add<TMatrixTSym<T>,TMatrixTSym<T> >
      {
        typedef TMatrixTSym<T> R ;
        //
        static  R operation 
        ( const TMatrixTSym<T>& m1 , 
          const TMatrixTSym<T>& m2 )
        { return R ( m1 , R::kPlus , m2 ); }
      } ;
      // ======================================================================


      // ======================================================================
      template <class T,unsigned int D1, unsigned int D2, class R1>
      struct Add<ROOT::Math::SMatrix<T,D1,D2,R1>,TMatrixT<T> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M1 ;
        typedef TMatrixT<T>                     M2 ;
        typedef ROOT::Math::SMatrix<T,D1,D2>    R  ;
        typedef SM<R>                           C1 ;
        //
        static R operation ( const M1& m1 , const M2& m2 )
        { return Add<M1,R>::operation ( m1 , C1::operation ( m2 ) ) ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D1, unsigned int D2, class R1>
      struct Add<TMatrixT<T>,ROOT::Math::SMatrix<T,D1,D2,R1> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M1 ;
        typedef TMatrixT<T>                     M2 ;
        typedef ROOT::Math::SMatrix<T,D1,D2>    R  ;
        typedef SM<R>                           C1 ;
        //
        static R operation ( const M2& m2 , const M1& m1 )
        { return Add<M1,R>::operation ( m1  , C1::operation ( m2 ) ) ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D1, unsigned int D2, class R1>
      struct Add<ROOT::Math::SMatrix<T,D1,D2,R1>,TMatrixTSym<T> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M1 ;
        typedef TMatrixTSym<T>                  M2 ;
        typedef ROOT::Math::SMatrix<T,D1,D2>    R  ;
        typedef SM<R>                           C1 ;
        //
        static R operation ( const M1& m1 , const M2& m2 )
        { return Add<M1,R>::operation ( m1 , C1::operation ( m2 ) ) ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D1, unsigned int D2, class R1>
      struct Add<TMatrixTSym<T>,ROOT::Math::SMatrix<T,D1,D2,R1> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M1 ;
        typedef TMatrixTSym<T>                  M2 ;
        typedef ROOT::Math::SMatrix<T,D1,D2>    R  ;
        typedef SM<R>                           C1 ;
        //
        static R operation ( const M2& m2 , const M1& m1 )
        { return Add<M1,R>::operation ( m1  , C1::operation ( m2 ) ) ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct Add<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >,
                 TMatrixTSym<T> >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef TMatrixTSym<T>                                         M2 ;
        typedef M1                                                     R  ;
        typedef SM<R>                                                  C1 ;
        //
        static M1 operation ( const M1& m1 , const M2& m2 )
        { return Add<M1,M1>::operation ( m1 , C1::operation ( m2 ) ) ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct Add<TMatrixTSym<T>,ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >   M1 ;
        typedef TMatrixTSym<T>                                           M2 ;
        typedef M1                                                       R  ;
        typedef SM<R>                                                    C1 ;
        //
        static M1 operation ( const M2& m2 , const M1& m1 )
        { return Add<M1,M1>::operation ( m1 , C1::operation ( m2 ) ) ; }
      } ;
      // ======================================================================

      // ======================================================================      
      template <class T>
      struct Add<TVectorT<T>, TVectorT<T> > 
      {
        typedef TVectorT<T>  M ;
        static TVectorT<T> operation ( const M& m1 , const M& m2 ) { return m1 + m2 ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct Add<ROOT::Math::SVector<T,D>, TVectorT<T> > 
      {
        typedef ROOT::Math::SVector<T,D>   M1 ;
        typedef TVectorT<T>                M2 ;
        typedef SM<M1>                     C1 ;
        typedef M1                         R  ;
        //
        static R operation ( const M1& m1 , const M2& m2 )
        { return Add<M1,M1>::operation ( m1 , C1::operation ( m2 ) ) ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct Add<TVectorT<T>,ROOT::Math::SVector<T,D> >
      {
        typedef ROOT::Math::SVector<T,D>    M1 ;
        typedef TVectorT<T>                 M2 ;
        typedef SM<M1>                       C1 ;
        typedef M1                           R  ;
        //
        static R operation ( const M2& m2 , const M1& m1 )
        { return Add<M1,M1>::operation ( C1::operation ( m2 ) , m1 ) ; }
      } ;
      // ======================================================================
      
      
      // ======================================================================
      // iadd
      // ======================================================================

      template <class T>
      struct IAdd<TMatrixT<T> , TMatrixT<T> >
      {
        typedef TMatrixT<T>                     M1 ;
        typedef TMatrixT<T>                     M2 ;
        // addition
        static void operation ( M1& m1 , const M2 & m2 ) { m1 += m2 ; }
      } ;

      template <class T>
      struct IAdd<TMatrixT<T> , TMatrixTSym<T> >
      {
        typedef TMatrixT<T>                     M1 ;
        typedef TMatrixTSym<T>                  M2 ;
        // addition
        static void operation ( M1& m1 , const M2 & m2 ) { m1 += m2 ; }
      } ;
      
      template <class T>
      struct IAdd<TMatrixTSym<T> , TMatrixTSym<T> >
      {
        typedef TMatrixTSym<T>                 M1 ;
        typedef TMatrixTSym<T>                 M2 ;
        // addition
        static void operation ( M1& m1 , const M2 & m2 ) { m1 += m2 ; }
      } ;

      
      template <class T,unsigned int D1,unsigned int D2, class R2>
      struct IAdd<TMatrixT<T>                     ,
                  ROOT::Math::SMatrix<T,D1,D2,R2> >
      {
        typedef TMatrixT<T>                     M1 ;
        typedef ROOT::Math::SMatrix<T,D1,D2,R2> M2 ;
        // addition
        static void operation ( M1& m1 , const M2 & m2 )
        {
          for ( unsigned int i = 0 ; i < D1 ; ++i )
          { for ( unsigned int j = 0 ; j < D2 ; ++j )
            { m1 ( i , j ) += m2 ( i , j ) ; } }
        } 
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct IAdd<TMatrixTSym<T> ,
                  ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > > 
      {
        typedef TMatrixTSym<T>                                          M1 ;
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >  M2 ;
        // addition
        static void operation ( M1& m1 , const M2 & m2 )
        {
          for ( unsigned int i = 0 ; i < D ; ++i )
          { for ( unsigned int j = i ; j < D ; ++j )
            { m1 ( i , j ) += m2 ( i , j ) ; } }
        } 
      } ;
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2>
      struct IAdd<ROOT::Math::SMatrix<T,D1,D2> ,
                  TMatrixT<T>                  >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2> M1 ;
        typedef TMatrixT<T>                  M2 ;
        // addition
        static void operation ( M1& m1 , const M2 & m2 )
        {
          for ( unsigned int i = 0 ; i < D1 ; ++i )
          { for ( unsigned int j = 0 ; j < D2 ; ++j )
            { m1 ( i , j ) += m2 ( i , j ) ; } }
        } 
      } ;
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2>
      struct IAdd<ROOT::Math::SMatrix<T,D1,D2> ,
                  TMatrixTSym<T>               >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2> M1 ;
        typedef TMatrixTSym<T>               M2 ;
        // addition
        static void operation ( M1& m1 , const M2 & m2 )
        {
          for ( unsigned int i = 0 ; i < D1 ; ++i )
          { for ( unsigned int j = 0 ; j < D2 ; ++j )
            { m1 ( i , j ) += m2 ( i , j ) ; } }
        } 
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct IAdd<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > ,
                  TMatrixTSym<T>                                         >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef TMatrixTSym<T>                                         M2 ;
        // addition
        static void operation ( M1& m1 , const M2 & m2 )
        {
          for ( unsigned int i = 0 ; i < D ; ++i )
          { for ( unsigned int j = 0 ; j < D ; ++j )
            { m1 ( i , j ) += m2 ( i , j ) ; } }
        } 
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct IAdd<ROOT::Math::SVector<T,D> ,
                  TVectorT<T>              >
      {
        typedef ROOT::Math::SVector<T,D> M1 ;
        typedef TVectorT<T>              M2 ;
        // addition
        static void operation ( M1& m1 , const M2 & m2 )
        { for ( unsigned int j = 0 ; j < D ; ++j ) { m1 ( j ) += m2 ( j ) ; } } 
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct IAdd<TVectorT<T>              ,
                  ROOT::Math::SVector<T,D> >
      {
        typedef TVectorT<T>              M1 ;
        typedef ROOT::Math::SVector<T,D> M2 ;
        // addition
        static void operation ( M1& m1 , const M2 & m2 )
        { for ( unsigned int j = 0 ; j < D ; ++j ) { m1 ( j ) += m2 ( j ) ; } } 
      } ;
      // ======================================================================
      
      // ======================================================================
      // sub
      // ======================================================================      
      template <class T>
      struct Sub<TMatrixT<T>,TMatrixT<T> >
      {
        typedef TMatrixT<T>  R ;
        //
        static R operation 
        ( const TMatrixT<T>& m1 , 
          const TMatrixT<T>& m2 ) { return R ( m1 , R::kMinus ,  m2 ) ; }
      } ;
      // ======================================================================
      template <class T>
      struct Sub<TMatrixT<T>,TMatrixTSym<T> >
      {
        typedef TMatrixT<T>  R ;
        //
        static R operation  
        ( const TMatrixT<T>&    m1 , 
          const TMatrixTSym<T>& m2 ) { return R ( m1 , R::kMinus ,  m2 )  ; }
      } ;
      // ======================================================================
      template <class T>
      struct Sub<TMatrixTSym<T>,TMatrixT<T> >
      {
        typedef TMatrixT<T>  R ;
        //
        static R operation
        ( const TMatrixTSym<T>& m1 , 
          const TMatrixT<T>&    m2 ) { return R ( m1 , R::kMinus ,  m2 )  ; }
      } ;
      // ======================================================================
      template <class T>
      struct Sub<TMatrixTSym<T>,TMatrixTSym<T> >
      {
        typedef TMatrixTSym<T>  R ;
        //        
        static R operation 
        ( const TMatrixTSym<T>& m1 , 
          const TMatrixTSym<T>& m2 ) { return R ( m1 , R::kMinus , m2 ); }
      } ;
      // ======================================================================
      template <class T,unsigned int D1, unsigned int D2, class R1>
      struct Sub<ROOT::Math::SMatrix<T,D1,D2,R1>,TMatrixT<T> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M1 ;
        typedef TMatrixT<T>                     M2 ;
        typedef ROOT::Math::SMatrix<T,D1,D2>    R  ;
        typedef SM<R>                           C1 ;
        //
        static R operation ( const M1& m1 , const M2& m2 )
        { return Sub<M1,R>::operation ( m1 , C1::operation ( m2 ) ) ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D1, unsigned int D2, class R1>
      struct Sub<TMatrixT<T>,ROOT::Math::SMatrix<T,D1,D2,R1> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M1 ;
        typedef TMatrixT<T>                     M2 ;
        typedef ROOT::Math::SMatrix<T,D1,D2>    R  ;
        typedef SM<R>                           C1 ;
        //
        static R operation ( const M2& m2 , const M1& m1 )
        { return Sub<R,M1>::operation ( C1::operation ( m2 ) , m1 ) ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D1, unsigned int D2, class R1>
      struct Sub<ROOT::Math::SMatrix<T,D1,D2,R1>,TMatrixTSym<T> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M1 ;
        typedef TMatrixTSym<T>                  M2 ;
        typedef ROOT::Math::SMatrix<T,D1,D2>    R  ;
        typedef SM<R>                           C1 ;
        //
        static R operation ( const M1& m1 , const M2& m2 )
        { return Sub<M1,R>::operation ( m1 , C1::operation ( m2 ) ) ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D1, unsigned int D2, class R1>
      struct Sub<TMatrixTSym<T>,ROOT::Math::SMatrix<T,D1,D2,R1> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M1 ;
        typedef TMatrixTSym<T>                  M2 ;
        typedef ROOT::Math::SMatrix<T,D1,D2>    R  ;
        typedef SM<R>                           C1 ;
        //
        static R operation ( const M2& m2 , const M1& m1 )
        { return Sub<R,M1>::operation ( C1::operation ( m2 ) , m1 ) ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct Sub<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >,
                 TMatrixTSym<T> >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef TMatrixTSym<T>                                         M2 ;
        typedef M1                                                     R  ;
        typedef SM<R>                                                  C1 ;
        //
        static M1 operation ( const M1& m1 , const M2& m2 )
        { return Sub<M1,R>::operation ( m1 , C1::operation ( m2 ) ) ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct Sub<TMatrixTSym<T>,ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef TMatrixTSym<T>                                         M2 ;
        typedef M1                                                     R  ;
        typedef SM<R>                                                  C1 ;
        //
        static M1 operation ( const M2& m2 , const M1& m1 )
        { return Sub<R,M1>::operation ( C1::operation ( m2 ) , m1 ) ; }
      } ;
      // ======================================================================      
      template <class T>
      struct Sub<TVectorT<T>, TVectorT<T> > 
      {
        typedef TVectorT<T>  M1 ;
        typedef TVectorT<T>  M2 ;
        typedef TVectorT<T>  R ;
        //
        static R operation ( const M1& m1 , const M1& m2 ) { return m1 - m2 ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct Sub<ROOT::Math::SVector<T,D>, TVectorT<T> > 
      {
        typedef ROOT::Math::SVector<T,D>    M1 ;
        typedef TVectorT<T>                 M2 ;
        typedef SM<M1>                      C1 ;
        typedef M1                          R  ;
        //
        static R operation ( const M1& m1 , const M2& m2 )
        { return Sub<M1,M1>::operation ( m1 , C1::operation ( m2 ) ) ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct Sub<TVectorT<T>,ROOT::Math::SVector<T,D> >
      {
        typedef ROOT::Math::SVector<T,D>    M1 ;
        typedef TVectorT<T>                 M2 ;
        typedef SM<M1>                       C1 ;
        typedef M1                           R  ;
        //
        static R operation  ( const M2& m2 , const M1& m1 )
        { return Sub<M1,M1>::operation ( C1::operation ( m2 ) , m1 ) ; }
      } ;
      // ======================================================================      
      
      // ======================================================================
      // isub
      // ======================================================================

      template <class T>
      struct ISub<TMatrixT<T> , TMatrixT<T> >
      {
        typedef TMatrixT<T>                     M1 ;
        typedef TMatrixT<T>                     M2 ;
        // addition
        static void operation ( M1& m1 , const M2 & m2 ) { m1 -= m2 ; }
      } ;

      template <class T>
      struct ISub<TMatrixT<T> , TMatrixTSym<T> >
      {
        typedef TMatrixT<T>                     M1 ;
        typedef TMatrixTSym<T>                  M2 ;
        // addition
        static void operation ( M1& m1 , const M2 & m2 ) { m1 -= m2 ; }
      } ;
      
      template <class T>
      struct ISub<TMatrixTSym<T> , TMatrixTSym<T> >
      {
        typedef TMatrixTSym<T>                 M1 ;
        typedef TMatrixTSym<T>                 M2 ;
        // addition
        static void operation ( M1& m1 , const M2 & m2 ) { m1 -= m2; }
      } ;

      
      template <class T,unsigned int D1,unsigned int D2, class R2>
      struct ISub<TMatrixT<T>                     ,
                  ROOT::Math::SMatrix<T,D1,D2,R2> >
      {
        typedef TMatrixT<T>                     M1 ;
        typedef ROOT::Math::SMatrix<T,D1,D2,R2> M2 ;
        // addition
        static void operation ( M1& m1 , const M2 & m2 )
        {
          for ( unsigned int i = 0 ; i < D1 ; ++i )
          { for ( unsigned int j = 0 ; j < D2 ; ++j )
            { m1 ( i , j ) -= m2 ( i , j ) ; } }
        } 
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct ISub<TMatrixTSym<T> ,
                  ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > > 
      {
        typedef TMatrixTSym<T>                                          M1 ;
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >  M2 ;
        // addition
        static void operation ( M1& m1 , const M2 & m2 )
        {
          for ( unsigned int i = 0 ; i < D ; ++i )
          { for ( unsigned int j = i ; j < D ; ++j )
            { m1 ( i , j ) -= m2 ( i , j ) ; } }
        } 
      } ;
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2>
      struct ISub<ROOT::Math::SMatrix<T,D1,D2> ,
                  TMatrixT<T>                  >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2> M1 ;
        typedef TMatrixT<T>                  M2 ;
        // addition
        static void operation ( M1& m1 , const M2 & m2 )
        {
          for ( unsigned int i = 0 ; i < D1 ; ++i )
          { for ( unsigned int j = 0 ; j < D2 ; ++j )
            { m1 ( i , j ) -= m2 ( i , j ) ; } }
        } 
      } ;
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2>
      struct ISub<ROOT::Math::SMatrix<T,D1,D2> ,
                  TMatrixTSym<T>               >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2> M1 ;
        typedef TMatrixTSym<T>               M2 ;
        // addition
        static void operation ( M1& m1 , const M2 & m2 )
        {
          for ( unsigned int i = 0 ; i < D1 ; ++i )
          { for ( unsigned int j = 0 ; j < D2 ; ++j )
            { m1 ( i , j ) -= m2 ( i , j ) ; } }
        } 
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct ISub<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > ,
                  TMatrixTSym<T>                                         >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef TMatrixTSym<T>                                         M2 ;
        // addition
        static void operation ( M1& m1 , const M2 & m2 )
        {
          for ( unsigned int i = 0 ; i < D ; ++i )
          { for ( unsigned int j = i ; j < D ; ++j )
            { m1 ( i , j ) -= m2 ( i , j ) ; } }
        } 
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct ISub<ROOT::Math::SVector<T,D> ,
                  TVectorT<T>              >
      {
        typedef ROOT::Math::SVector<T,D> M1 ;
        typedef TVectorT<T>              M2 ;
        // addition
        static void operation ( M1& m1 , const M2 & m2 )
        { for ( unsigned int j = 0 ; j < D ; ++j ) { m1 ( j ) -= m2 ( j ) ; } } 
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct ISub<TVectorT<T>              ,
                  ROOT::Math::SVector<T,D> >
      {
        typedef TVectorT<T>              M1 ;
        typedef ROOT::Math::SVector<T,D> M2 ;
        // addition
        static void operation ( M1& m1 , const M2 & m2 )
        { for ( unsigned int j = 0 ; j < D ; ++j ) { m1 ( j ) -= m2 ( j ) ; } } 
      } ;
      // ======================================================================
      
      // ======================================================================
      // mul
      // =======================================================================
      template <class T>
      struct Mul<TMatrixT<T>,TMatrixT<T> >
      {
        typedef TMatrixT<T> R ;
        static  R operation 
        ( const TMatrixT<T>& m1 , 
          const TMatrixT<T>& m2 ) { return m1 * m2 ; }
      } ;
      // =======================================================================      
      template <class T>
      struct Mul<TMatrixT<T>,TMatrixTSym<T> >
      {
        typedef TMatrixT<T> R ;
        static  R operation 
        ( const TMatrixT<T>&    m1 , 
          const TMatrixTSym<T>& m2 ) { return m1 * m2 ; }
      } ;
      // =======================================================================
      template <class T>
      struct Mul<TMatrixTSym<T>,TMatrixT<T> >
      {
        typedef TMatrixT<T> R ;
        static  R operation 
        ( const TMatrixTSym<T>& m1 , 
          const TMatrixT<T>&    m2 ) { return m1 * m2 ; }
      } ;
      // =======================================================================      
      template <class T>
      struct Mul<TMatrixTSym<T>,TMatrixTSym<T> >
      {
        typedef TMatrixT<T> R ;
        static  R operation 
        ( const TMatrixTSym<T>& m1 , 
          const TMatrixTSym<T>& m2 ) { return m1 * m2 ; }
      } ;
      // =======================================================================
      template <class T>
      struct Mul<TMatrixT<T>,TVectorT<T> >
      {
        typedef TVectorT<T> R ;
        static  R operation  
        ( const TMatrixT<T>& m1 , 
          const TVectorT<T>& m2 ) { return m1 * m2 ; }
      } ;
      // =======================================================================      
      template <class T>
      struct Mul<TMatrixTSym<T>,TVectorT<T> >
      {
        typedef TVectorT<T> R ;
        static  R operation 
        ( const TMatrixTSym<T>& m1 , 
          const TVectorT<T>&    m2 ) { return m1 * m2 ; }
      } ;
      // ======================================================================
      template <class T>
      struct Mul<TVectorT<T>,TMatrixT<T> >
      {
        typedef TVectorT<T> M1 ;
        typedef TMatrixT<T> M2 ;
        typedef TVectorT<T> R ;
        static  R operation 
        ( const TVectorT<T>& m1 , 
          const TMatrixT<T>& m2 ) { return  M2 ( M2::kTransposed , m2 ) * m1 ; }
      } ;
      // =======================================================================
      template <class T>
      struct Mul<TVectorT<T>,TMatrixTSym<T>>
      {
        typedef TVectorT<T> M1 ;
        typedef TMatrixT<T> M2 ;
        typedef TVectorT<T> R ;
        static  R operation 
        ( const TVectorT<T>&    m1 , 
          const TMatrixTSym<T>& m2 )  { return  M2 ( M2::kTransposed , m2 ) * m1 ; }
      } ;
      // ======================================================================      
      template <class T, unsigned int D1, unsigned int D2, class R1>
      struct Mul<ROOT::Math::SMatrix<T,D1,D2,R1> , 
                 TMatrixT<T> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1>   M1 ;
        typedef TMatrixT<T>                       M2 ;
        typedef TM<M1>                            C  ;
        typedef typename C::R                     NM ;
        typedef Mul<NM,M2>                        M  ;
        typedef typename M::R                     R  ;
        //
        static  R operation 
        ( const M1& m1 , 
          const M2& m2 )
        { return M::operation ( C::operation ( m1 )  , m2 ) ; }
      } ;
      // ======================================================================
      template <class T, unsigned int D1, unsigned int D2, class R1>
      struct Mul<ROOT::Math::SMatrix<T,D1,D2,R1> , 
                 TMatrixTSym<T> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1>   M1 ;
        typedef TMatrixTSym<T>                    M2 ;
        typedef TM<M1>                            C  ;
        typedef typename C::R                     NM ;
        typedef Mul<NM,M2>                        M  ;
        typedef typename M::R                     R  ;
        //        
        static  R operation 
        ( const M1& m1 , 
          const M2& m2 )
        { return M::operation ( C::operation ( m1 )  , m2 ) ; }
      } ;
      // ======================================================================
      template <class T, unsigned int D1, unsigned int D2, class R1>
      struct Mul<ROOT::Math::SMatrix<T,D1,D2,R1> , 
                 TVectorT<T> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1>   M1 ;
        typedef TVectorT<T>                       M2 ;
        typedef ROOT::Math::SVector<T,D2>         M3 ;
        typedef SM<M3>                            C  ;
        typedef typename C::R                     NM ;
        typedef Mul<M1,NM>                        M  ;
        typedef typename M::R                     R  ;
        //
        static  R operation 
        ( const M1& m1 , 
          const M2& m2 )
        { return M::operation ( m1 , C::operation ( m2 ) ) ; }
      } ;
      // ======================================================================




      // ======================================================================
      template <class T, unsigned int D1, unsigned int D2, class R1>
      struct Mul<TMatrixT<T>,
                 ROOT::Math::SMatrix<T,D1,D2,R1> > 
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1>   M1 ;
        typedef TMatrixT<T>                       M2 ;
        typedef TM<M1>                            C  ;
        typedef typename C::R                     NM ;
        typedef Mul<M2,NM>                        M  ;
        typedef typename M::R                     R  ;
        //
        static  R operation
        ( const M2& m2 , 
          const M1& m1 )
        { return M::operation ( m2 , C::operation ( m1 ) ) ; }
      } ;
      // ======================================================================
      template <class T, unsigned int D1, unsigned int D2, class R1>
      struct Mul<TMatrixTSym<T>,
                 ROOT::Math::SMatrix<T,D1,D2,R1> > 
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1>   M1 ;
        typedef TMatrixTSym<T>                    M2 ;
        typedef TM<M1>                            C  ;
        typedef typename C::R                     NM ;
        typedef Mul<M2,NM>                        M  ;
        typedef typename M::R                     R  ;
        //
        static  R operation 
        ( const M2& m2 , 
          const M1& m1 )
        { return M::operation ( m2 , C::operation ( m1 ) ) ; }
      } ;
      // ======================================================================
      template <class T, unsigned int D>
      struct Mul<TMatrixT<T>,
                 ROOT::Math::SVector<T,D> > 
      {
        typedef ROOT::Math::SVector<T,D>   M1 ;
        typedef TMatrixT<T>                M2 ;
        typedef TM<M1>                     C  ;
        typedef typename C::R              NM ;
        typedef Mul<M2,NM>                 M  ;
        typedef typename M::R              R  ;
        //
        static  R operation
        ( const M2& m2 , 
          const M1& m1 )
        { return M::operation ( m2 , C::operation ( m1 ) ) ; }
      } ;

      // ======================================================================
      template <class T, unsigned int D>
      struct Mul<TMatrixTSym<T> ,
                 ROOT::Math::SVector<T,D> > 
      {
        typedef ROOT::Math::SVector<T,D>   M1 ;
        typedef TMatrixTSym<T>             M2 ;
        typedef TM<M1>                     C  ;
        typedef typename C::R              NM ;
        typedef Mul<M2,NM>                 M  ;
        typedef typename M::R              R  ;
        //
        static  R operation
        ( const M2& m2 , 
          const M1& m1 )
        { return M::operation ( m2 , C::operation ( m1 ) ) ; }
      } ;
      // ======================================================================
      template <class T, unsigned int D>
      struct Mul<ROOT::Math::SVector<T,D> ,
                 TMatrixT<T> > 
      {
        typedef ROOT::Math::SVector<T,D>   M1 ;
        typedef TMatrixT<T>                M2 ;
        typedef TM<M1>                     C  ;
        typedef typename C::R              NM ;
        typedef Mul<NM,M2>                 M  ;
        typedef typename M::R              R  ;
        //
        static  R operation 
        ( const M1& m2 , 
          const M2& m1 )
        { return M::operation ( C::operation ( m2 ) , m1 ) ; }
      } ;
      // ======================================================================      
      template <class T, unsigned int D>
      struct Mul<ROOT::Math::SVector<T,D> ,
                 TMatrixTSym<T> > 
      {
        typedef ROOT::Math::SVector<T,D>   M1 ;
        typedef TMatrixTSym<T>             M2 ;
        typedef TM<M1>                     C  ;
        typedef typename C::R              NM ;
        typedef Mul<NM,M2>                 M  ;
        typedef typename M::R              R  ;
        //
        static  R operation
        ( const M1& m2 , 
          const M2& m1 )
        { return M::operation ( C::operation ( m2 ) , m1 ) ; }
      } ;
      // ======================================================================

      /// new (1)
      template <class T, unsigned int D>
      struct Mul<ROOT::Math::SVector<T,D> ,
                 TVectorT<T> > 
      {
        typedef ROOT::Math::SVector<T,D>   M1 ;
        typedef TVectorT<T>                M2 ;
        typedef SM<M1>                     C  ;
        typedef typename C::R              NM ;
        typedef Mul<M1,NM>                 M  ;
        typedef typename M::R              R  ;
        //
        static  R operation
        ( const M1& m1 , 
          const M2& m2 )
        { return M::operation ( m1 , C::operation ( m2 ) ) ; }
      } ;
      // ======================================================================

      /// new (2)
      template <class T, unsigned int D>
      struct Mul<TVectorT<T>              ,
                 ROOT::Math::SVector<T,D> >
      {
        typedef TVectorT<T>                M1 ;
        typedef ROOT::Math::SVector<T,D>   M2 ;
        typedef SM<M2>                     C  ;
        typedef typename C::R              NM ;
        typedef Mul<NM,M2>                 M  ;
        typedef typename M::R              R  ;
        //
        static  R operation
        ( const M1& m1 , 
          const M2& m2 )
        { return M::operation ( C::operation ( m1 ) , m2) ; }
      } ;
      // ======================================================================

      /// new (3)
      template <class T, unsigned int D>
      struct RMul<ROOT::Math::SVector<T,D> ,
                  TVectorT<T> > 
      {
        typedef ROOT::Math::SVector<T,D>   M1 ;
        typedef TVectorT<T>                M2 ;
        typedef SM<M1>                     C  ;
        typedef typename C::R              NM ;
        typedef Mul<M1,NM>                 M  ;
        typedef typename M::R              R  ;
        //
        static  R operation
        ( const M1& m1 , 
          const M2& m2 )
        { return M::operation ( C::operation ( m2 ) , m1 ) ; }
      } ;
      // ======================================================================

      /// new (4)
      template <class T, unsigned int D1, unsigned D2, class R1>
      struct RMul<ROOT::Math::SMatrix<T,D1,D2,R1> ,
                  TVectorT<T> > 
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M1 ;
        typedef TVectorT<T>                     M2 ;
        typedef ROOT::Math::SVector<T,D1>       P1 ;
        typedef SM<P1>                          C  ;
        typedef typename C::R                   NM ;
        typedef Mul<NM,M1>                      M  ;
        typedef typename M::R                   R  ;
        //
        static  R operation
        ( const M1& m1 , 
          const M2& m2 )
        { return M::operation ( C::operation ( m2 ) , m1 ) ; }
      } ;
      // ======================================================================




      // ======================================================================      
      template <class T>
      struct Mul<TVectorT<T> ,
                 TVectorT<T> > 
      {
        typedef TVectorT<T> M1 ;
        typedef TVectorT<T> M2 ;
        typedef double      R  ;
        //
        static  R operation 
        ( const M1& m1 , 
          const M2& m2 ) { return ::Dot ( m1 , m2 ) ; }
      } ;
      // ======================================================================

      // ======================================================================
      template <class T>
      struct IMul<TMatrixT<T>, TMatrixT<T> >
      {
        typedef TMatrixT<T> M1 ;
        typedef TMatrixT<T> M2 ;
        //
        static void operation ( M1& m1 , const M2& m2 ) { m1 *= m2 ; }
      } ;      
      // ======================================================================
      template <class T>
      struct IMul<TMatrixT<T>, TMatrixTSym<T> >
      {
        typedef TMatrixT<T>    M1 ;
        typedef TMatrixTSym<T> M2 ;
        //
        static void operation ( M1& m1 , const M2& m2 ) { m1 *= m2 ; }
      } ;
      // ======================================================================
      template <class T, unsigned int D, class R1>
      struct IMul<TMatrixT<T>,
                  ROOT::Math::SMatrix<T,D,D,R1> >
      {
        //
        typedef TMatrixT<T>                   M1 ;
        typedef ROOT::Math::SMatrix<T,D,D,R1> M2 ;
        typedef TM<M2>                        C  ;
        //
        static void operation ( M1& m1 , const M2& m2 ) { m1 *= C::operation ( m2 ) ; }
      } ;
      // ======================================================================
      template <class T, unsigned int D1, unsigned D2>
      struct IMul<ROOT::Math::SMatrix<T,D1,D2>, 
                  TMatrixT<T> >
      {
        //
        typedef ROOT::Math::SMatrix<T,D1,D2> M1 ;
        typedef TMatrixT<T>                  M2 ;
        typedef SM<M1>                        C  ;
        //
        static void operation ( M1& m1 , const M2& m2 ) { m1 *= C::operation ( m2 ) ; }
      } ;
      // ======================================================================
      template <class T, unsigned int D1, unsigned D2>
      struct IMul<ROOT::Math::SMatrix<T,D1,D2>, 
                  TMatrixTSym<T> >
      {
        //
        typedef ROOT::Math::SMatrix<T,D1,D2> M1 ;
        typedef TMatrixTSym<T>               M2 ;
        typedef SM<M1>                        C  ;
        //
        static void operation ( M1& m1 , const M2& m2 ) { m1 *= C::operation ( m2 ) ; }
      } ;
      // ======================================================================
      
      
      // ======================================================================
      /// can they be compared ?
      // ======================================================================
      template <class T1, class T2> 
      struct CanEq < TVectorT<T1> , 
                     TVectorT<T2> >
      {
        typedef TVectorT<T1> M1 ;
        typedef TVectorT<T2> M2 ;
        // equality 
        static bool operation ( const M1& m1 , 
                                const M2& m2 ) 
        { return m1.IsValid() && m2.IsValid() && m1.GetNrows() == m2.GetNrows() ; }
      } ;
      // =========================================================================
      template <class T1, unsigned int D, class T2> 
      struct CanEq < ROOT::Math::SVector<T1,D> , 
                     TVectorT<T2>              >
      {
        typedef ROOT::Math::SVector<T1,D> M1 ;
        typedef TVectorT<T2>              M2 ;
        // equality 
        static bool operation ( const M1& m1 , 
                                const M2& m2 ) 
        { return m2.IsValid() && D == m2.GetNrows() ; }
      } ;
      
      // =========================================================================
      template <class T1, unsigned int D, class T2> 
      struct CanEq < TVectorT<T2>              ,
                     ROOT::Math::SVector<T1,D> >
      {
        typedef TVectorT<T2>              M1 ;
        typedef ROOT::Math::SVector<T1,D> M2 ;
        // equality 
        static bool operation ( const M1& m1 , 
                                const M2& m2 ) 
        { return m1.IsValid() && D == m1.GetNrows() ; }
      } ;

      // ======================================================================
      // EQUALITY
      // ======================================================================
      template <class T>
      struct Eq <TVectorT<T> ,
                 TVectorT<T> >
      {
        typedef TVectorT<T>  M1 ;
        typedef TVectorT<T>  M2 ;
        typedef bool          R  ;
        // addition
        static R operation ( const M1& m1 , const M2& m2 )
        {
          static const Ostap::Math::Equal_To<M1> s_cmp ;
          return s_cmp ( m1 , m2 ) ;
        }
      };
      // ======================================================================
      template <class T, unsigned int D>
      struct Eq <TVectorT<T>,
                 ROOT::Math::SVector<T,D> >
      {
        typedef TVectorT<T>               M1 ;
        typedef ROOT::Math::SVector<T,D>  M2 ;
        typedef bool                      R  ;
        // addition
        static R operation ( const M1& m1 , const M2& m2 )
        {
          static const Ostap::Math::Equal_To<M1> s_cmp ;
          return s_cmp ( m1 , m2 ) ;
        }
      };
      // ======================================================================
      template <class T, unsigned int D>
      struct Eq <ROOT::Math::SVector<T,D>,
                 TVectorT<T> >
      {
        typedef ROOT::Math::SVector<T,D>  M1 ;
        typedef TVectorT<T>               M2 ;
        typedef bool                      R  ;
        // addition
        static R operation ( const M1& m1 , const M2& m2 )
        {
          static const Ostap::Math::Equal_To<M2> s_cmp ;
          return s_cmp ( m2 , m1 ) ;
        }
      };
      // ============================================================================

      
      // ======================================================================
      /// can they be compared ?
      // ======================================================================
      template <class T1, class T2> 
      struct CanEq < TMatrixT<T1> , 
                     TMatrixT<T2> >
      {
        typedef TMatrixT<T1> M1 ;
        typedef TMatrixT<T2> M2 ;
        // equality 
        static bool operation 
        ( const M1& m1 , 
          const M2& m2 ) 
        { return m1.IsValid() && m2.IsValid() 
            && m1.GetNrows() == m2.GetNrows() 
            && m1.GetNcols() == m2.GetNcols() ; }
      } ;
      // ======================================================================
      template <class T1, class T2> 
      struct CanEq < TMatrixT<T1>    , 
                     TMatrixTSym<T2> >
      {
        typedef TMatrixT<T1>    M1 ;
        typedef TMatrixTSym<T2> M2 ;
        // equality 
        static bool operation 
        ( const M1& m1 , 
          const M2& m2 ) 
        { return m1.IsValid() && m2.IsValid() 
            && m1.GetNrows() == m2.GetNrows() 
            && m1.GetNcols() == m2.GetNcols() ; }
      } ;
      // ======================================================================
      template <class T1, class T2> 
      struct CanEq < TMatrixTSym<T1>    , 
                     TMatrixT<T2>       >
      {
        typedef TMatrixTSym<T1> M1 ;
        typedef TMatrixT<T2>    M2 ;
        // equality 
        static bool operation ( const M1& m1 , 
                                const M2& m2 ) 
        { return m1.IsValid() && m2.IsValid() 
            && m1.GetNrows() == m2.GetNrows() 
            && m1.GetNcols() == m2.GetNcols() ; }
      } ;
      // ======================================================================
      template <class T1, class T2> 
      struct CanEq < TMatrixTSym<T1>    , 
                     TMatrixTSym<T2> >
      {
        typedef TMatrixTSym<T1> M1 ;
        typedef TMatrixTSym<T2> M2 ;
        // equality 
        static bool operation ( const M1& m1 , 
                                const M2& m2 ) 
        { return m1.IsValid() && m2.IsValid() 
            && m1.GetNrows() == m2.GetNrows() 
            && m1.GetNcols() == m2.GetNcols() ; }
      } ;
      // ======================================================================

      template <class T1, class T2>
      struct Eq <TMatrixT<T1> ,
                 TMatrixT<T2> >
      {
        typedef TMatrixT<T1>  M1 ;
        typedef TMatrixT<T2>  M2 ;
        typedef bool          R  ;
        // addition
        static R operation ( const M1& m1 , const M2& m2 )
        {
          static const Ostap::Math::Equal_To<M1> s_cmp ;
          return s_cmp ( m1  , m2 ) ;
        }
      };
      // ======================================================================
      template <class T1, class T2>
      struct Eq <TMatrixT<T1>    ,
                 TMatrixTSym<T2> >
      {
        typedef TMatrixT<T1>     M1 ;
        typedef TMatrixTSym<T2>  M2 ;
        typedef bool             R  ;
        // addition
        static R operation ( const M1& m1 , const M2& m2 )
        {
          static const Ostap::Math::Equal_To<M1> s_cmp ;
          return s_cmp ( m1  , m2 ) ;
        }
      };
      // ======================================================================      
      template <class T1, class T2>
      struct Eq <TMatrixTSym<T1> ,
                 TMatrixT<T2>    >
      {
        typedef TMatrixTSym<T1>  M1 ;
        typedef TMatrixT<T2>     M2 ;
        typedef bool             R  ;
        // addition
        static R operation ( const M1& m1 , const M2& m2 )
        {
          static const Ostap::Math::Equal_To<M2> s_cmp ;
          return s_cmp ( m2  , m1 ) ;
        }
      };
      // ======================================================================
      template <class T1, class T2>
      struct Eq <TMatrixTSym<T1> ,
                 TMatrixTSym<T2> >
      {
        typedef TMatrixTSym<T1>  M1 ;
        typedef TMatrixTSym<T2>  M2 ;
        typedef bool             R  ;
        // addition
        static R operation ( const M1& m1 , const M2& m2 )
        {
          static const Ostap::Math::Equal_To<M1> s_cmp ;
          return s_cmp ( m1  , m2 ) ;
        }
      };
      // ======================================================================

      
      // ======================================================================
      template <class T1, class T2, unsigned int D1, unsigned int D2, class R1>
      struct CanEq < ROOT::Math::SMatrix<T1,D1,D2,R1>, 
                     TMatrixT<T2>  >
      {
        typedef ROOT::Math::SMatrix<T1,D1,D2,R1>  M1 ;
        typedef TMatrixT<T2>                      M2 ;
        //
        static bool operation ( const M1& /* m1 */ , 
                                const M2&    m2    ) 
        { return m2.IsValid() && m2.GetNrows() == D1 && m2.GetNcols() == D2 ; }
      } ;
      //  =====================================================================
      template <class T1, class T2, unsigned int D1, unsigned int D2, class R1>
      struct CanEq < TMatrixT<T2>  , 
                     ROOT::Math::SMatrix<T1,D1,D2,R1> >
      {
        typedef ROOT::Math::SMatrix<T1,D1,D2,R1>  M1 ;
        typedef TMatrixT<T2>                      M2 ;
        //
        static bool operation 
          ( const M2&    m2    , 
            const M1& /* m1 */  ) 
        { return m2.IsValid() && m2.GetNrows() == D1 && m2.GetNcols() == D2 ; }
      } ;
      //  =====================================================================
      template <class T1, class T2, unsigned int D, class R1>
      struct CanEq < ROOT::Math::SMatrix<T1,D,D,R1>, 
                     TMatrixTSym<T2>  >
      {
        typedef ROOT::Math::SMatrix<T1,D,D,R1>  M1 ;
        typedef TMatrixTSym<T2>                   M2 ;
        static bool operation ( const M1& /* m1 */ , 
                                const M2&    m2    ) 
        { return m2.IsValid() && m2.GetNrows() == D  ; }
      } ;
      //  =====================================================================
      template <class T1, class T2, unsigned int D, class R1>
      struct CanEq < TMatrixTSym<T2> , 
                     ROOT::Math::SMatrix<T1,D,D,R1> >
      {
        typedef ROOT::Math::SMatrix<T1,D,D,R1>  M1 ;
        typedef TMatrixTSym<T2>                   M2 ;
        static bool operation 
        ( const M2&    m2    , 
          const M1& /* m1 */ ) 
        { return m2.IsValid() && m2.GetNrows() == D  ; }
      } ;
      
      
      // ======================================================================
      template <class T1, class T2, unsigned int D1, unsigned int D2, class R1>
      struct Eq <ROOT::Math::SMatrix<T1,D1,D2,R1>,
                 TMatrixT<T2> >
      {
        typedef ROOT::Math::SMatrix<T1,D1,D2,R1>  M1 ;
        typedef TMatrixT<T2>                      M2 ;
        typedef bool                              R  ;
        // addition
        static R operation ( const M1& m1 , const M2& m2 )
        {
          static const Ostap::Math::Equal_To<M2> s_cmp ;
          return s_cmp ( m2 , m1 ) ;
        }
      };
      // ======================================================================
      template <class T1 , class T2, unsigned int D, class R1>
      struct Eq <ROOT::Math::SMatrix<T1,D,D,R1>,
                 TMatrixTSym<T2> >
      {
        typedef ROOT::Math::SMatrix<T1,D,D,R1>  M1 ;
        typedef TMatrixTSym<T2>                 M2 ;
        typedef bool                            R  ;
        // addition
        static R operation ( const M1& m1 , const M2& m2 )
        {
          static const Ostap::Math::Equal_To<M2> s_cmp ;
          return s_cmp ( m2 , m1 ) ;
        }
      };
      // ======================================================================
      template <class T1, class T2, unsigned int D1, unsigned int D2, class R1>
      struct Eq <TMatrixT<T2> ,
                 ROOT::Math::SMatrix<T1,D1,D2,R1> >
      {
        typedef ROOT::Math::SMatrix<T2,D1,D2,R1>  M2 ;
        typedef TMatrixT<T2>                      M1 ;
        typedef bool                              R  ;
        // addition
        static R operation ( const M1& m1 , const M2& m2 )
        {
          static const Ostap::Math::Equal_To<M1> s_cmp ;
          return s_cmp ( m1 , m2 ) ;
        }
      };
      // ======================================================================
      template <class T1, class T2, unsigned int D, class R1>
      struct Eq <TMatrixTSym<T2>,
                 ROOT::Math::SMatrix<T1,D,D,R1> >
      {
        typedef ROOT::Math::SMatrix<T1,D,D,R1>  M2 ;
        typedef TMatrixTSym<T2>                 M1 ;
        typedef bool                            R  ;
        // addition
        static R operation ( const M1& m1 , const M2& m2 )
        {
          static const Ostap::Math::Equal_To<M1> s_cmp ;
          return s_cmp ( m1 , m2 ) ;
        }
      };
      // ======================================================================


      
      // ======================================================================
      template <class T>
      struct CanEq < TMatrixT<T> , double>
      {
        typedef TMatrixT<T>     M1 ;
        typedef double          M2 ;
        static bool operation 
        ( const M1&    m1    , 
          const M2  /* m2 */ ) 
        { return m1.IsValid() && m1.GetNrows() == m1.GetNcols (); }
      } ;
      // ======================================================================
      template <class T>
      struct CanEq < TMatrixTSym<T> , double>
      {
        typedef TMatrixTSym<T>  M1 ;
        typedef double          M2 ;
        static bool operation 
        ( const M1&    m1    , 
          const M2  /* m2 */ ) 
        { return m1.IsValid() ; }
      } ;


      // ======================================================================
      template <class T>
      struct Eq <TMatrixT<T>, double>
      {
        typedef TMatrixT<T>     M1 ;
        typedef double          M2 ;
        //
        static bool operation 
        ( const M1&    m1    , 
          const M2     m2    ) 
        { 
          //
          static const Ostap::Math::Equal_To<T>  s_cmp  ;
          static const Ostap::Math::Zero<T>      s_zero ;
          //
          const Int_t D = m1.GetNrows() ;
          for ( Int_t i = 0 ; i < D ; ++i ) 
          {
            if ( !s_cmp  ( m1 ( i , i ) , m2 )        ) { return false ; }            
            for ( Int_t j = 0 ; j < D ;  ++j ) 
            { 
              if ( i != j && !s_zero ( m1 ( i , j ) ) ) { return false ; }
            }
          }
          return true ;
        }
      } ;
      // =================================================
      template <class T>
      struct Eq <TMatrixTSym<T>, double>
      {
        typedef TMatrixTSym<T>  M1 ;
        typedef double          M2 ;
        //
        static bool operation 
        ( const M1&    m1    , 
          const M2     m2    ) 
        { 
          //
          static const Ostap::Math::Equal_To<T>  s_cmp  ;
          static const Ostap::Math::Zero<T>      s_zero ;
          //
          const Int_t D = m1.GetNrows() ;
          for ( Int_t i = 0 ; i < D ; ++i ) 
          {
            if ( !s_cmp  ( m1 ( i , i ) , m2 )  ) { return false ; }            
            for ( Int_t j = i + 1 ; j < D ; ++j ) 
            { 
              if ( !s_zero ( m1 ( i , j )    )  ) { return false ; }
            }
          }
          return true ;
        }
      } ;
      
      
      // ======================================================================
      template <class T>
      struct Pow<TMatrixT<T> >
      {
        //
        typedef TMatrixT<T> M ;
        typedef TMatrixT<T> R ;
        //
        static R operation ( const M& m , const unsigned short n )
        {
          //
          if      ( 0 == n ) { return M ( M::kUnit , m ) ; }
          else if ( 1 == n ) { return         m ; }
          else if ( 2 == n ) { return     m * m ; }
          else if ( 3 == n ) { return m * m * m ; }
          //
          R r = operation ( m , n / 2 )  ;
          if ( 0 == n / 2 ) { return r * r ; }
          //
          return r * r * m ;
        }
        // 
      } ;
      // ======================================================================
      template <class T>
      struct Pow<TMatrixTSym<T> >
      {
        //
        typedef TMatrixTSym<T> M ;
        typedef TMatrixT<T>    R ;
        //
        static R operation ( const M& m , const unsigned short n )
        {
          //
          if      ( 0 == n ) { return M ( M::kUnit , m ) ; }
          else if ( 1 == n ) { return         m ; }
          else if ( 2 == n ) { return     m * m ; }
          else if ( 3 == n ) { return m * m * m ; }
          //
          R r = operation ( m , n / 2 )  ;
          if ( 0 == n / 2 ) { return r * r ; }
          //
          return r * r * m ;
        }
      } ;

      // ======================================================================
      template <class T>
      struct CanInvert< TMatrixT<T> >  
      {
        typedef TMatrixT<T>     M1 ;
        typedef bool            R  ;
        //
        static R operation ( const M1& m1 ) 
        { return m1.IsValid() && m1.GetNrows() == m1.GetNcols() ; } 
      } ;
      // ======================================================================
      template <class T>
      struct CanInvert< TMatrixTSym<T> >  
      {
        typedef TMatrixTSym<T>  M1 ;
        typedef bool            R  ;
        //
        static R operation ( const M1& m1 ) { return m1.IsValid() ; } 
      } ;
      
      // ======================================================================
      template <class T>
      struct Invert< TMatrixT<T> >  
      {
        typedef TMatrixT<T>     M1 ;
        typedef TMatrixT<T>     R  ;
        //
        static R operation ( const M1& m1 , int& flag  )
        { 
          if ( !m1.IsValid() ) { flag  = 1 ; return m1 ; }
          //
          M1 result { m1 } ;
          Double_t det = 1 ;
          result.Invert ( &det ) ;
          flag = ( 0 == det ) ;
          return result ;
        }
      } ;      
      // ======================================================================
      template <class T>
      struct Invert< TMatrixTSym<T> >  
      {
        typedef TMatrixTSym<T>  M1 ;
        typedef TMatrixTSym<T>  R  ;
        //
        static R operation ( const M1& m1 , int& flag  )
        { 
          if ( !m1.IsValid() ) { flag  = 1 ; return m1 ; }
          //
          M1 result { m1 } ;
          Double_t det = 1 ;
          result.Invert ( &det ) ;
          flag = ( 0 == det ) ;
          return result ;
        }
      } ;      
      // ======================================================================

      // ======================================================================
      template <class T1,unsigned int D, class T2>
      struct CanDot<ROOT::Math::SVector<T1,D> , 
                    TVectorT<T2>              >
      {
        static bool operation 
        ( const ROOT::Math::SVector<T1,D>& /* m1 */ ,
          const TVectorT<T2>&                 m2    )  
        { return m2.IsValid() && D == m2.GetNrows () ; }
        
      } ;
      // ======================================================================
      template <class T1,unsigned int D, class T2>
      struct CanDot<TVectorT<T1>, 
                    ROOT::Math::SVector<T2,D> >
      {
        static bool operation 
        ( const TVectorT<T1>&                 m1    , 
          const ROOT::Math::SVector<T2,D>& /* m2 */ )
        { return m1.IsValid() && D == m1.GetNrows () ; }
      } ;
      // ======================================================================
      template <class T1, class T2>
      struct CanDot<TVectorT<T1>, 
                    TVectorT<T2>>  
      {
        static bool operation 
          ( const TVectorT<T1>& m1 ,
            const TVectorT<T2>& m2 )
        { return m1.IsValid() && m2.IsValid() && m1.GetNrows () == m2.GetNrows () ; }
      } ;
      // ======================================================================
      
      // ======================================================================
      // DOT
      // ======================================================================
      template <class T1,unsigned int D, class T2>
      struct Dot<ROOT::Math::SVector<T1,D> ,
                 TVectorT<T2> >
      {
        typedef ROOT::Math::SVector<T1,D>  M1 ;
        typedef TVectorT<T1>               M2 ;
        typedef double                     R  ;
        // multiplication
        static double operation 
        ( const M1 & m1 , 
          const M2 & m2 )
        { return std::inner_product ( m1.begin () , 
                                      m1.end   ()   , 
                                      m2.GetMatrixArray() , 0.0 ) ; }
      } ;            
      // ======================================================================
      template <class T1,unsigned int D,class T2>
      struct Dot< TVectorT<T1> , 
                  ROOT::Math::SVector<T2,D> >
      {
        typedef TVectorT<T1>               M1 ;
        typedef ROOT::Math::SVector<T2,D>  M2 ;
        typedef double                     R  ;
        // multiplication
        static double operation 
        ( const M1 & m1 , 
          const M2 & m2 )
        { return std::inner_product ( m2.begin () , 
                                      m2.end   ()   , 
                                      m1.GetMatrixArray() , 0.0 ) ; }
      } ;            
      // ======================================================================
      template <class T1,class T2>
      struct Dot< TVectorT<T1> , 
                  TVectorT<T2> >
      {
        typedef TVectorT<T1>               M1 ;
        typedef TVectorT<T2>               M2 ;
        typedef double                     R  ;
        // multiplication
        static double operation 
        ( const M1 & m1 , 
          const M2 & m2 )
        { return std::inner_product ( m1.GetMatrixArray() , 
                                      m1.GetMatrixArray() + m1.GetNrows () , 
                                      m2.GetMatrixArray() , 0.0 ) ; }
      } ;            
      // ======================================================================

      
      // ======================================================================
      template <class T1,unsigned int D, class T2>
      struct CanCross<ROOT::Math::SVector<T1,D> , 
                      TVectorT<T2>              >
      {
        static bool operation 
        ( const ROOT::Math::SVector<T1,D>& /* m1 */ ,
          const TVectorT<T2>&                 m2    )  
        { return m2.IsValid() ; }
      } ;
      // =======================================================================
      template <class T1,unsigned int D, class T2>
      struct CanCross< TVectorT<T2>  , 
                       ROOT::Math::SVector<T1,D> >
      {
        static bool operation 
        ( const TVectorT<T2>&                 m2    , 
          const ROOT::Math::SVector<T1,D>& /* m1 */ )
        { return m2.IsValid() ; }
      } ;
      // =======================================================================
      template <class T1>
      struct CanCross< TVectorT<T1> , 
                       TVectorT<T1> >
      {
        static bool operation 
        ( const TVectorT<T1>& m1    ,
          const TVectorT<T1>& m2    ) 
        { return m1.IsValid() && m2.IsValid() ; }
      } ;
      // =======================================================================
      template <class T>
      struct Cross< TVectorT<T> ,
                    TVectorT<T> > 
      {
        typedef TVectorT<T>  M1 ;
        typedef TVectorT<T>  M2 ;
        typedef TMatrixT<T>  R  ;
        //
        static  R operation
        ( const M1& m1 , 
          const M2& m2 )
        { return ::OuterProduct ( m1 , m2 ) ; }
      } ;
      // ======================================================================
      template <class T, unsigned int D>
      struct Cross<ROOT::Math::SVector<T,D> ,
                   TVectorT<T> > 
      {
        typedef ROOT::Math::SVector<T,D>   M1 ;
        typedef TVectorT<T>                M2 ;
        typedef TM<M1>                     C  ;
        typedef typename C::R              NM ;
        typedef Cross<NM,M2>               M  ;
        typedef typename M::R              R  ;
        //
        static  R operation
        ( const M1& m1 , 
          const M2& m2 )
        { return M::operation ( C::operation ( m1 ) , m2  ) ; }
      } ;
      // ======================================================================
      template <class T, unsigned int D>
      struct Cross< TVectorT<T> , 
                    ROOT::Math::SVector<T,D> >
      {
        typedef ROOT::Math::SVector<T,D>   M1 ;
        typedef TVectorT<T>                M2 ;
        typedef TM<M1>                     C  ;
        typedef typename C::R              NM ;
        typedef Cross<M2,NM>               M  ;
        typedef typename M::R              R  ;
        //
        static  R operation
        ( const M2& m2 , 
          const M1& m1 )
        { return M::operation ( m2 , C::operation ( m1 ) ) ; }
      } ;
      // ======================================================================

      
      template <class T,unsigned int D>
      struct CanSim<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D>> ,
                    TVectorT<T> >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef TVectorT<T>                                            M2 ;
        // check 
        static bool operation 
          ( const M1& /* m1 */ , const M2& m2 ) 
        { return m2.IsValid() && D == m2.GetNrows() ; }
      } ;
      // ======================================================================      
      template <class T,unsigned int D>
      struct CanSim<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D>> ,
                    TMatrixTSym<T>  >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef TMatrixTSym<T>                                         M2 ;
        // check 
        static bool operation 
          ( const M1& /* m1 */ , const M2& m2 ) 
        { return m2.IsValid() && D == m2.GetNrows() ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct CanSim<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D>> ,
                    TMatrixT<T>  >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef TMatrixT<T>                                            M2 ;
        // check 
        static bool operation 
          ( const M1& /* m1 */ , const M2& m2 ) 
        { return m2.IsValid() && D == m2.GetNcols() ; }
      } ;
      
      // ===================================================================
      template <class T>
      struct CanSim< TMatrixTSym<T> , 
                     TVectorT<T>    > 
      {
        typedef TMatrixTSym<T>  M1 ;
        typedef TVectorT<T>     M2 ;
        // check 
        static bool operation 
        ( const M1& m1 , const M2& m2 ) 
        { return m1.IsValid() && m2.IsValid() && m1.GetNcols() == m2.GetNrows() ; }
      } ;
      // ======================================================================
      template <class T>
      struct CanSim< TMatrixTSym<T> , 
                     TMatrixT<T>    > 
      {
        typedef TMatrixTSym<T>  M1 ;
        typedef TMatrixT<T>     M2 ;
        // check 
        static bool operation 
        ( const M1& m1 , const M2& m2 ) 
        { return m1.IsValid() && m2.IsValid() && m1.GetNcols() == m2.GetNcols() ; }
      } ;
      // ======================================================================
      template <class T>
      struct CanSim< TMatrixTSym<T> , 
                     TMatrixTSym<T> > 
      {
        typedef TMatrixTSym<T>  M1 ;
        typedef TMatrixTSym<T>  M2 ;
        // check 
        static bool operation 
        ( const M1& m1 , const M2& m2 ) 
        { return m1.IsValid() && m2.IsValid() && m1.GetNcols() == m2.GetNcols() ; }
      } ;
      // ======================================================================


      
      // ======================================================================
      template <class T, unsigned int D>
      struct CanSim< TMatrixTSym<T> , 
                     ROOT::Math::SVector<T,D> > 
      {
        typedef TMatrixTSym<T>           M1 ;
        typedef ROOT::Math::SVector<T,D> M2 ;
        // check 
        static bool operation 
        ( const M1&  m1 , const M2& /* m2 */ ) 
        { return m1.IsValid() && D == m1.GetNcols() ; }
      } ;
      // ======================================================================
      template <class T, unsigned int D1, unsigned int D2, class R1>
      struct CanSim< TMatrixTSym<T> , 
                     ROOT::Math::SMatrix<T,D1,D2,R1> > 
      {
        typedef TMatrixTSym<T>                   M1 ;
        typedef ROOT::Math::SMatrix<T,D1,D2,R1>  M2 ;
        // check 
        static bool operation 
        ( const M1&  m1 , const M2& /* m2 */ ) 
        { return m1.IsValid() && D2 == m1.GetNcols() ; }
      } ;
      // ======================================================================
      
      
      // ======================================================================
      template <class T, unsigned int D>
      struct Sim< TMatrixTSym<T> , 
                  ROOT::Math::SVector<T,D> > 
      {
        typedef TMatrixTSym<T>                                             M1 ;
        typedef ROOT::Math::SVector<T,D>                                   M2 ;
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >     M3 ;
        typedef SM<M3>                                                     C  ;
        typedef typename C::R                                              NM ;
        typedef Sim<NM,M2>                                                 O  ;
        typedef typename O::R                                              R  ;        
        // operation 
        static R operation 
        ( const M1&  m1 , const M2& m2 ) 
        { return O::operation ( C::operation ( m1 ) , m2 ) ; }
      } ;
      
      // ======================================================================
      template <class T, unsigned int D1, unsigned D2, class R1>
      struct Sim< TMatrixTSym<T> , 
                  ROOT::Math::SMatrix<T,D1,D2,R1> > 
      {
        typedef TMatrixTSym<T>                                                M1 ;
        typedef ROOT::Math::SMatrix<T,D1,D2,R1>                               M2 ;
        typedef ROOT::Math::SMatrix<T,D2,D2,ROOT::Math::MatRepSym<T,D2> >     M3 ;
        typedef SM<M3>                                                        C  ;
        typedef typename C::R                                                 NM ;
        typedef Sim<NM,M2>                                                    O  ;
        typedef typename O::R                                                 R  ;        
        // operation 
        static R operation 
        ( const M1&  m1 , const M2& m2 ) 
        { return O::operation ( C::operation ( m1 ) , m2 ) ; }
      } ;

      // ======================================================================
      template <class T>
      struct Sim< TMatrixTSym<T> , 
                  TVectorT<T>    > 
      {
        typedef TMatrixTSym<T>  M1 ;
        typedef TVectorT<T>     M2 ;
        typedef double          R  ;
        // operation 
        static R operation 
        ( const M1& m1 , const M2& m2 ) { return m1.Similarity ( m2 ) ; }
      } ;
      // ======================================================================
      template <class T>
      struct Sim< TMatrixTSym<T> , 
                  TMatrixT<T>    > 
      {
        typedef TMatrixTSym<T>  M1 ;
        typedef TMatrixT<T>     M2 ;
        typedef TMatrixTSym<T>  R  ;
        // operation 
        static R operation 
        ( const M1& m1 , 
          const M2& m2 ) 
        { M1 res { m1 } ; res.Similarity ( m2 ) ; return res ; }
      } ;
     // =======================================================================
      template <class T>
      struct Sim< TMatrixTSym<T> , 
                  TMatrixTSym<T> > 
      {
        typedef TMatrixTSym<T>  M1 ;
        typedef TMatrixTSym<T>  M2 ;
        typedef TMatrixTSym<T>  R  ;
        // operation 
        static R operation 
        ( const M1& m1 , 
          const M2& m2 ) 
        { M1 res { m1 } ; res.Similarity ( m2 ) ; return res ; }
      } ;

      // =====================================================================
      template <class T,unsigned int D> 
      struct Sim<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > ,
                 TVectorT<T> >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef TVectorT<T>                                            M2 ;
        typedef ROOT::Math::SVector<T,D>                               M3 ;
        typedef SM<M3>                                                 C  ;
        typedef Sim<M1,M3>                                             O  ;
        typedef typename O::R                                          R  ;
        // similarity 
        static R operation ( const M1& m1 , const M2& m2 )
        { return O::operation ( m1 , C::operation ( m2 ) ) ; }
      } ;
      // ======================================================================      
      template <class T,unsigned int D> 
      struct Sim<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > ,
                 TMatrixTSym<T> >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef TMatrixTSym<T>                                         M2 ; 
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M3 ;
        typedef SM<M3>                                                 C  ;
        typedef Sim<M1,M3>                                             O  ;
        typedef typename O::R                                          R  ;
        // similarity 
        static R operation ( const M1& m1 , const M2& m2 )
        { return O::operation ( m1 , C::operation ( m2 ) ) ; }
      } ;
      
      // ======================================================================      
      template <class T,unsigned int D> 
      struct Sim<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > ,
                 TMatrixT<T> >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef TMatrixT<T>                                            M2 ;
        typedef TMatrixTSym<T>                                         M3 ;
        typedef TM<M1>                                                 C  ;
        typedef Sim<M3,M2>                                             O  ;
        typedef typename O::R                                          R  ;
        // similarity 
        static R operation ( const M1& m1 , const M2& m2 )
        { return O::operation ( C::operation ( m1 ) , m2 ) ; }
      } ;



      // ======================================================================
      template <class T> 
      struct CanSimT< TMatrixTSym<T> , 
                      TMatrixT<T>    > 
      {
        typedef TMatrixTSym<T>   M1 ;
        typedef TMatrixT<T>      M2 ;
        // check 
        static bool operation 
        ( const M1& m1  , const M2& m2 ) 
        { return m1.IsValid() && m2.IsValid() && m1.GetNrows() == m2.GetNrows() ; }
      } ;
      // ======================================================================
      template <class T>
      struct CanSimT< TMatrixTSym<T> , 
                      TMatrixTSym<T> > 
      {
        typedef TMatrixTSym<T>   M1 ;
        typedef TMatrixTSym<T>   M2 ;
        // check 
        static bool operation 
        ( const M1& m1  , const M2& m2 ) 
        { return m1.IsValid() && m2.IsValid() && m1.GetNrows() == m2.GetNrows() ; }
      } ;
      
      
      
      // ======================================================================
      template <class T>
      struct SimT< TMatrixTSym<T> , 
                   TMatrixT<T>    > 
      {
        typedef TMatrixTSym<T>   M1 ;
        typedef TMatrixT<T>      M2 ;
        typedef TMatrixTSym<T>   R  ;
        // check 
        static R operation 
        ( const M1& m1  , const M2& m2 ) 
        { M1 result { m1 } ; result.SimilarityT ( m2 ) ; return result ; }
      } ;      
      // ======================================================================
      template <class T>
      struct SimT< TMatrixTSym<T> , 
                   TMatrixTSym<T> > 
      {
        typedef TMatrixTSym<T>   M1 ;
        typedef TMatrixTSym<T>   M2 ;
        typedef TMatrixTSym<T>   R  ;
        // check 
        static R operation 
        ( const M1& m1  , const M2& m2 ) 
        { M1 result { m1 } ; result.Similarity ( m2 ) ; return result ; }
      } ;
      
      // ==================================================================
      template <class T,unsigned int D>
      struct CanSimT<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > ,
                     TMatrixT<T>  >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef TMatrixT<T>                                            M2 ;
        // check 
        static bool operation 
        ( const M1& /* m1 */ , const M2& m2 ) 
        { return m2.IsValid() && D == m2.GetNrows() ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct CanSimT<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D>> ,
                     TMatrixTSym<T>  >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef TMatrixTSym<T>                                         M2 ;
        // check 
        static bool operation 
        ( const M1& /* m1 */ , const M2& m2 ) 
        { return m2.IsValid() && D == m2.GetNrows() ; }
      } ;


      // ===================================================================
      template <class T,unsigned int D> 
      struct SimT<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > ,
                  TMatrixT<T> >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef TMatrixT<T>                                            M2 ;
        typedef TMatrixTSym<T>                                         M3 ;
        typedef TM<M1>                                                 C  ;
        typedef SimT<M3,M2>                                            O  ;
        typedef typename O::R                                          R  ;
        // similarity 
        static R operation ( const M1& m1 , const M2& m2 )
        { return O::operation ( C::operation ( m1 ) , m2 ) ; }
      } ;
      // =====================================================================
      template <class T,unsigned int D> 
      struct SimT<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > ,
                  TMatrixTSym<T> >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef TMatrixTSym<T>                                         M2 ;
        typedef TMatrixTSym<T>                                         M3 ;
        typedef TM<M1>                                                 C  ;
        typedef SimT<M3,M2>                                            O  ;
        typedef typename O::R                                          R  ;
        // similarity 
        static R operation ( const M1& m1 , const M2& m2 )
        { return O::operation ( C::operation ( m1 ) , m2 ) ; }
      } ;

      
      // ======================================================================
      template <class T , unsigned int D1, unsigned int D2, class R1>
      struct CanSimT< TMatrixTSym<T> , 
                      ROOT::Math::SMatrix<T,D1,D2,R1> >
      {
        typedef TMatrixTSym<T>                  M1 ;
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M2 ;
        // check 
        static bool operation 
        ( const M1& m1  , const M2& /* m2 */ ) 
        { return m1.IsValid() && D1 == m1.GetNrows() ; }
      } ;

      
      // ======================================================================
      template <class T , unsigned int D1, unsigned int D2, class R1>
      struct SimT< TMatrixTSym<T> , 
                   ROOT::Math::SMatrix<T,D1,D2,R1> >
      {
        typedef TMatrixTSym<T>                                            M1 ;
        typedef ROOT::Math::SMatrix<T,D1,D2,R1>                           M2 ;
        typedef ROOT::Math::SMatrix<T,D1,D1,ROOT::Math::MatRepSym<T,D1> > M3 ;
        typedef SM<M3>                                                    C  ;
        typedef SimT<M3,M2>                                               O  ;
        typedef typename O::R                                             R  ;
        // 
        static R operation 
        ( const M1& m1  , const M2& m2 ) 
        { return O::operation ( C::operation ( m1 ) , m2 ) ; }
      } ;
      
      
      
      
      // ======================================================================
      template <class T>
      struct Sym<TMatrixT<T> >
      {
        //
        typedef TMatrixT<T>    M ;
        typedef TMatrixTSym<T> R ;
        //
        static R operation ( const M& m )
        {
          //
          const unsigned int nr = m.GetNrows () ;
          const unsigned int nc = m.GetNcols () ;
          //
          R r ( nr ) ;
          //
          for ( unsigned int i = 0 ; i < nr ; ++i )
          {
            r ( i , i ) = m ( i , i ) ;
            for ( unsigned int j = i + 1 ; j < nr ; ++j )
            {
              const T v = 0.5 * ( m ( i , j ) + m ( j , i ) ) ;
              r ( i , j ) = v ;
              r ( j , i ) = v ;
            }
          }
          return r ;
        }
      } ;
      
      
      // ======================================================================
      template <class T>
      struct Sym<TMatrixTSym<T> >
      {
        //
        typedef TMatrixTSym<T> M ;
        typedef TMatrixTSym<T> R ;
        //
        static R operation ( const M& m )
        {
          //
          const unsigned int nr = m.GetNrows () ;
          const unsigned int nc = m.GetNcols () ;
          //
          R r ( nr ) ;
          //
          for ( unsigned int i = 0 ; i < nr ; ++i )
          {
            r ( i , i ) = m ( i , i ) ;
            for ( unsigned int j = i + 1 ; j < nr ; ++j )
            {
              const T v = 0.5 * ( m ( i , j ) + m ( j , i ) ) ;
              r ( i , j ) = v ;
              r ( j , i ) = v ;
            }
          }
          return r ;
        }
      } ;
      
      // ======================================================================
      template <class T>
      struct ASym<TMatrixT<T> >
      {
        //
        typedef TMatrixT<T> M ;
        typedef TMatrixT<T> R  ;
        //
        static R operation ( const M& m )
        {
          //
          const unsigned int nr = m.GetNrows () ;
          const unsigned int nc = m.GetNcols () ;
          //
          R r ( nr , nc ) ;
          //
          for ( unsigned int i = 0 ; i < nr ; ++i )
          {
            r ( i , i ) = 0 ;
            for ( unsigned int j = i + 1 ; j < nr ; ++j )
            {
              const T v = 0.5 * ( m ( i , j ) -  m ( j , i ) ) ;
              r ( i , j ) =  v ;
              r ( j , i ) = -v ; 
            }
          }
          return r ;
        }
      } ;
      
      // ======================================================================
      template <class T>
      struct ASym<TMatrixTSym<T> >
      {
        //
        typedef TMatrixTSym<T> M ;
        typedef TMatrixT<T>    R  ;
        //
        static R operation ( const M& m )
        {
          //
          const unsigned int nr = m.GetNrows () ;
          const unsigned int nc = m.GetNcols () ;
          //
          R r ( nr , nc ) ;
          //
          for ( unsigned int i = 0 ; i < nr ; ++i )
          {
            r ( i , i ) = 0 ;
            for ( unsigned int j = i + 1 ; j < nr ; ++j )
            {
              const T v = 0.5 * ( m ( i , j ) -  m ( j , i ) ) ;
              r ( i , j ) =  v ;
              r ( j , i ) = -v ; 
            }
          }
          return r ;
        }
      } ;

      // ======================================================================      
    } //                                  The end of namespace Ostap::Math::Ops
    // ========================================================================

    
    
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                           The end of namespace Ostap::Math
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_MATRIXUTILST_H
// ============================================================================
