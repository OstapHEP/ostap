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
      template <class T2>
      inline bool operator()
      ( const TVectorT<T>&  v1 , 
        const TVectorT<T2>& v2 ) const
      {
        return ( v1.IsValid  () && v2.IsValid  () &&
                 v1.GetNrows () == v2.GetNrows () && 
                 std::equal ( v1.GetMatrixArray () , v1.GetMatrixArray () + v1.GetNrows() ,
                              v2.GetMatrixArray () , m_cmp ) ) ;                 
      }
      /// compare with another vector type (e.g. double and float)
      template <class T2>
      inline bool operator()
      ( const TVectorT<T2>& v1 , 
        const TVectorT<T>&  v2 ) const { return (*this)( v2  , v1 ) ; }
      // ======================================================================
      template <class T2, unsigned int D>
      inline bool operator ()
      ( const TVectorT<T>&               v1 ,
        const ROOT::Math::SVector<T,D>&  v2 ) const
      {
        return v1.IsValid() && D == v1.GetNrows() &&
          std::equal ( v2.begin() , v2.end()  , v1.GetMatrixArray () , m_cmp );
      }
      // ======================================================================
      template <class T2, unsigned int D>
      inline bool operator ()
      ( const ROOT::Math::SVector<T,D>&  v1 , 
        const TVectorT<T>&               v2 ) const
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
      /// compare with another matrix type (e.g. double and float)
      inline bool operator()
      ( const TMatrixT<T>&    v1 , 
        const TMatrixTSym<T>& v2 ) const
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
        //
        return true ;
      }
      /// compare with another matrix type (e.g. double and float)
      inline bool operator()
      ( const TMatrixTSym<T>&  v1 , 
        const TMatrixT<T>&     v2 ) const { return  (*this) ( v2 , v1 ) ; }
      // ======================================================================
      template <unsigned int D1 , unsigned int D2, class R1>
      inline bool operator()
      ( const TMatrixT<T>&                     v1 , 
        const ROOT::Math::SMatrix<T,D1,D2,R1>& v2 ) const
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
      template <unsigned int D1 , unsigned int D2, class R1>
      inline bool operator()
      ( const ROOT::Math::SMatrix<T,D1,D2,R1>& v1 , 
        const TMatrixT<T>&                     v2 ) const
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
      /// compare with another matrix type (e.g. double and float)
      inline bool operator()
      ( const TMatrixTSym<T>& v1 , 
        const TMatrixT<T>&    v2 ) const
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
      /// compare with another matrix type (e.g. double and float)
      inline bool operator()
      ( const TMatrixT<T>&    v1 , 
        const TMatrixTSym<T>& v2 ) const { return (*this)( v2  , v1 ) ; }
      // ======================================================================
      template <unsigned int D1 , unsigned int D2, class R1>
      inline bool operator()
      ( const TMatrixTSym<T>&                  v1 , 
        const ROOT::Math::SMatrix<T,D1,D2,R1>& v2 ) const
      {
        if ( !v1.IsValid() || v1.GetNrows() != D1 || v1.GetNcols() != D2 ) { return false ; }
        //
        for ( unsigned long i = 0 ; i < D1 ; ++i )
        { for ( unsigned long j = 0 ; j < D2 ; ++j )
          { if ( !m_cmp ( v1 ( i , j ) , v2  ( i , j ) ) ) { return false ; } } }
        //
        return true ;
      }
      // ======================================================================
      template <unsigned int D1 , unsigned int D2, class R1>
      inline bool operator()
      ( const ROOT::Math::SMatrix<T,D1,D2,R1>& v1 , 
        const TMatrixTSym<T>&                  v2 ) const { return  (*this) ( v2 , v1 ) ; }
      // ======================================================================
    private:
      // ======================================================================
      /// the evaluator 
      Equal_To<T> m_cmp ;                                 // the evaluator 
      // ======================================================================
    } ;
    // ========================================================================



    
    // ========================================================================    
    namespace  Ops
    {      
      // ======================================================================
      // converters 
      // ======================================================================
      template <class M>
      struct TM ;
      template <class M>
      struct SM ;

      
      // ======================================================================
      //  S --> T 
      // ======================================================================
      template <class T, unsigned int D1, unsigned int D2>
      struct TM<ROOT::Math::SMatrix<T,D1,D2> >
      {
        typedef TMatrixT<T>  R ;
        // 
        static R transform 
        ( const ROOT::Math::SMatrix<T,D1,D2>& m )
        {
          //
          R result ( D1 , D2 ) ;
          for ( unsigned int i = 0 ; i < D1 ; ++i )
          { for ( unsigned int j = 0 ; j < D1 ; ++j )
            { result ( i , j ) = m ( i , j ) ; } }
          //
          return result ;
        }
      } ;
      // ======================================================================
      template <class T, unsigned int D>
      struct TM<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > >
      {
        typedef TMatrixTSym<T>  R ;
        //
        static R transform
        ( const ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >& m )
        {
          //
          R result ( D ) ;
          for ( unsigned int i = 0 ; i < D ; ++i )
          { for ( unsigned int j = 0 ; j < D ; ++j )
            { result ( i , j ) = m ( i , j ) ; } }
          //
          return result ;
        }
      } ;
      // ======================================================================
      template <class T, unsigned int D>
      struct TM<ROOT::Math::SVector<T,D> >
      {
        typedef TVectorT<T>  R ;
        //
        static R transform ( const ROOT::Math::SVector<T,D>& m )
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
        static R transform ( const TMatrixT<T>& m )
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
        static R transform ( const TMatrixTSym<T>& m )
        {
          R result ;
          for ( unsigned int i = 0 ; i < D ; ++i  )
          { for ( unsigned int j = 0 ; j < D ; ++j )
            { result ( i , j ) = m ( i , j ) ; } }
          return result ;
        } 
        static R transform ( const TMatrixT<T>& m )
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
        static R transform ( const TMatrixTSym<T>& m )
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
        static R transform ( const TVectorT<T>& m )
        {
          const T* start = m.GetMatrixArray() ;
          return R ( start , start + D ) ;
        }
      } ;
      // ======================================================================
      


      // ======================================================================
      // CHECKERS 
      // ======================================================================
      template <class T1,class T2>
      struct CanAdd<TMatrixT<T1>,TMatrixT<T2> >
      {
        static bool ok
        ( const TMatrixT<T1>& m1 , 
          const TMatrixT<T2>& m2 )
        {
          return
            m1.IsValid  () && m2.IsValid  () &&
            m1.GetNrows () == m2.GetNrows () &&
            m1.GetNcols () == m2.GetNcols () ;
        }  
      } ;
      // ======================================================================
      template <class T1,class T2>
      struct CanAdd<TMatrixT<T1>,TMatrixTSym<T2> >
      {
        static bool ok
        ( const TMatrixT<T1>&    m1 , 
          const TMatrixTSym<T2>& m2 )
        {
          return
            m1.IsValid  () && m2.IsValid  () &&
            m1.GetNrows () == m2.GetNrows () &&
            m1.GetNcols () == m2.GetNcols () ;
        }  
      } ;
      // ======================================================================
      template <class T1,class T2>
      struct CanAdd<TMatrixTSym<T1>,TMatrixT<T2> >
      {
        static bool ok
        ( const TMatrixTSym<T1>& m1 , 
          const TMatrixT<T2>&    m2 )
        {
          return
            m1.IsValid  () && m2.IsValid  () &&
            m1.GetNrows () == m2.GetNrows () &&
            m1.GetNcols () == m2.GetNcols () ;
        }  
      } ;
      // ======================================================================
      template <class T1,class T2>
      struct CanAdd<TMatrixTSym<T1>,TMatrixTSym<T2> >
      {
        static bool ok
        ( const TMatrixTSym<T1>& m1 , 
          const TMatrixTSym<T2>& m2 )
        {
          return
            m1.IsValid  () && m2.IsValid  () &&
            m1.GetNrows () == m2.GetNrows () &&
            m1.GetNcols () == m2.GetNcols () ;
        }  
      } ;
      // ======================================================================      
      template <class T1,unsigned int D1, unsigned int D2, class R, class T2>
      struct CanAdd<ROOT::Math::SMatrix<T1,D1,D2,R>,TMatrixT<T2> >
      {
        static bool ok
        ( const ROOT::Math::SMatrix<T1,D1,D2,R>& m1 , 
          const TMatrixT<T2>&                    m2 )
        {
          return
            m2.IsValid        () &&
            D1 == m2.GetNrows () &&
            D2 == m2.GetNcols () ;
        }  
      } ;
      // ======================================================================
      template <class T1,unsigned int D1, unsigned int D2, class R, class T2>
      struct CanAdd<TMatrixT<T2>,ROOT::Math::SMatrix<T1,D1,D2,R> >
      {
        static bool ok
        ( const TMatrixT<T2>&                    m2 ,          
          const ROOT::Math::SMatrix<T1,D1,D2,R>& m1 ) 
        {
          return
            m2.IsValid        () &&
            D1 == m2.GetNrows () &&
            D2 == m2.GetNcols () ;
        }  
      } ;
      // ======================================================================
      template <class T1,unsigned int D1, unsigned int D2, class R, class T2>
      struct CanAdd<ROOT::Math::SMatrix<T1,D1,D2,R>,TMatrixTSym<T2> >
      {
        static bool ok
        ( const ROOT::Math::SMatrix<T1,D1,D2,R>& m1 , 
          const TMatrixTSym<T2>&                 m2 )
        {
          return
            m2.IsValid        () &&
            D1 == m2.GetNrows () &&
            D2 == m2.GetNcols () ;
        }  
      } ;
      // ======================================================================
      template <class T1,unsigned int D1, unsigned int D2, class R, class T2>
      struct CanAdd<TMatrixTSym<T2>,ROOT::Math::SMatrix<T1,D1,D2,R> >
      {
        static bool ok
        ( const TMatrixTSym<T2>&                 m2 ,          
          const ROOT::Math::SMatrix<T1,D1,D2,R>& m1 ) 
        {
          return
            m2.IsValid        () &&
            D1 == m2.GetNrows () &&
            D2 == m2.GetNcols () ;
        }  
      } ;
      // ======================================================================
      template <class T1,class T2>
      struct CanAdd<TVectorT<T1>,TVectorT<T2> >
      {
        static bool ok
        ( const TVectorT<T1>& m1 ,  
          const TVectorT<T2>& m2 ) 
        { return m1.IsValid() && m2.IsValid () && m1.GetNrows()== m2.GetNrows () ; }
      };
      // ======================================================================
      template <class T1,unsigned int D,class T2>
      struct CanAdd<ROOT::Math::SVector<T1,D>,TVectorT<T2> >
      {
        static bool ok
        ( const ROOT::Math::SVector<T1,D> & /* m1 */ ,  
          const TVectorT<T2>&               m2        ) 
        { return m2.IsValid () && D == m2.GetNrows () ; }
      };
      // ======================================================================
      template <class T1,unsigned int D,class T2>
      struct CanAdd<TVectorT<T2>,ROOT::Math::SVector<T1,D> >
      {
        static bool ok
        ( const TVectorT<T2>&                  m2     ,
          const ROOT::Math::SVector<T1,D> & /* m1 */  )
        { return m2.IsValid () && D == m2.GetNrows () ; }
      };
      // ======================================================================


      
      // ======================================================================
      template <class T1,class T2>
      struct CanMul<TMatrixT<T1>,TMatrixT<T2>>
      {
        static bool ok
        ( const TMatrixT<T1>& m1 ,
          const TMatrixT<T2>& m2 )
        { return m1.IsValid() &&  m2.IsValid() && m1.GetNcols() == m2.GetNrows() ; }
      };
      // ======================================================================
      template <class T1,class T2>
      struct CanMul<TMatrixT<T1>,TMatrixTSym<T2>>
      {
        static bool ok
        ( const TMatrixT<T1>   & m1 ,
          const TMatrixTSym<T2>& m2 )
        { return m1.IsValid() &&  m2.IsValid() && m1.GetNcols() == m2.GetNrows() ; }
      };
      // ======================================================================      
      template <class T1,class T2>
      struct CanMul<TMatrixTSym<T1>,TMatrixT<T2>>
      {
        static bool ok
        ( const TMatrixTSym<T1> & m1 ,
          const TMatrixT<T2>    & m2 )
        { return m1.IsValid() &&  m2.IsValid() && m1.GetNcols() == m2.GetNrows() ; }
      };
      // ======================================================================      
      template <class T1,class T2>
      struct CanMul<TMatrixTSym<T1>,TMatrixTSym<T2>>
      {
        static bool ok
        ( const TMatrixTSym<T1> & m1 ,
          const TMatrixTSym<T2> & m2 )
        { return m1.IsValid() &&  m2.IsValid() && m1.GetNcols() == m2.GetNrows() ; }
      };
      // ======================================================================      
      template <class T1,class T2>
      struct CanMul<TMatrixT<T1>,TVectorT<T2>>
      {
        static bool ok
        ( const TMatrixT<T1>& m1 ,
          const TVectorT<T2>& m2 )
        { return m1.IsValid() &&  m2.IsValid() && m1.GetNcols() == m2.GetNrows() ; }
      };
      // ======================================================================
      template <class T1,class T2>
      struct CanMul<TMatrixTSym<T1>,TVectorT<T2>>
      {
        static bool ok
        ( const TMatrixTSym<T1>& m1 ,
          const TVectorT<T2>&    m2 )
        { return m1.IsValid() &&  m2.IsValid() && m1.GetNcols() == m2.GetNrows() ; }
      } ;
      // ======================================================================
      template <class T1,class T2>
      struct CanMul<TVectorT<T2>,TMatrixT<T1>>
      {
        static bool ok
        ( const TVectorT<T2>& m1 ,
          const TMatrixT<T1>& m2 )
        { return m1.IsValid() &&  m2.IsValid() && m1.GetNrows () == m2.GetNrows() ; }
      };
      // ======================================================================
      template <class T1,class T2>
      struct CanMul<TVectorT<T2>,TMatrixTSym<T1>>
      {
        static bool ok
        ( const TVectorT<T2>& m1 ,
          const TMatrixT<T1>& m2 )
        { return m1.IsValid() &&  m2.IsValid() && m1.GetNrows () == m2.GetNrows() ; }
      };
      // ======================================================================
      
      template <class T1,class T2>
      struct CanMul<TVectorT<T2>,TVectorT<T1>>
      {
        static bool ok
        ( const TVectorT<T2>& m1 ,
          const TVectorT<T1>& m2 )
        { return m1.IsValid() &&  m2.IsValid() && m1.GetNrows () == m2.GetNrows() ; }
      };
      // ======================================================================

      
      // ======================================================================
      template <class T1,unsigned int D1, unsigned int D2, class R1,
                class T2>
      struct CanMul<ROOT::Math::SMatrix<T1,D1,D2,R1>,TMatrixT<T2>>
      {
        static bool ok
        ( const ROOT::Math::SMatrix<T1,D1,D2,R1>& /*  m1 */ ,
          const TMatrixT<T2>&                         m2 )
        { return m2.IsValid() &&  D2 == m2.GetNrows() ; }
      } ;
      // ======================================================================
      template <class T1,unsigned int D1, unsigned int D2, class R1,
                class T2>
      struct CanMul<ROOT::Math::SMatrix<T1,D1,D2,R1>,TMatrixTSym<T2>>
      {
        static bool ok
        ( const ROOT::Math::SMatrix<T1,D1,D2,R1>& /*  m1 */ ,
          const TMatrixTSym<T2>&                      m2 )
        { return m2.IsValid() &&  D2 == m2.GetNrows() ; }
      } ;
      // ======================================================================
      template <class T1,unsigned int D1, unsigned int D2, class R1,
                class T2>
      struct CanMul<ROOT::Math::SMatrix<T1,D1,D2,R1>,TVectorT<T2>>
      {
        static bool ok
        ( const ROOT::Math::SMatrix<T1,D1,D2,R1>& /*  m1 */ ,
          const TVectorT<T2>&                         m2 )
        { return m2.IsValid() &&  D2 == m2.GetNrows() ; }
      } ;
      // ======================================================================
      template <class T1,unsigned int D,
                class T2>
      struct CanMul<ROOT::Math::SVector<T1,D>,
                    TVectorT<T2>>
      {
        static bool ok
        ( const ROOT::Math::SVector<T1,D>& /*  m1 */ ,
          const TVectorT<T2>&                 m2 )
        { return m2.IsValid() &&  D == m2.GetNrows() ; }
      } ;
      // ======================================================================
      template <class T1,unsigned int D,
                class T2>
      struct CanMul<ROOT::Math::SVector<T1,D>,
                    TMatrixT<T2>>
      {
        static bool ok
        ( const ROOT::Math::SVector<T1,D>& /*  m1 */ ,
          const TMatrixT<T2>&                  m2 )
        { return m2.IsValid() &&  D == m2.GetNrows() ; }
      } ;
      // ======================================================================
      template <class T1,unsigned int D,
                class T2>
      struct CanMul<ROOT::Math::SVector<T1,D>,
                    TMatrixTSym<T2>>
      {
        static bool ok
        ( const ROOT::Math::SVector<T1,D>& /*  m1 */ ,
          const TMatrixTSym<T2>&               m2 )
        { return m2.IsValid() &&  D == m2.GetNrows() ; }
      } ;

      
      // ======================================================================
      template <class T1,unsigned int D1, unsigned int D2, class R1,
                class T2>
      struct CanMul<TMatrixT<T2>,
                    ROOT::Math::SMatrix<T1,D1,D2,R1> >
      {
        static bool ok
        ( const TMatrixT<T2>&                         m2    ,
          const ROOT::Math::SMatrix<T1,D1,D2,R1>& /*  m1 */ )
        { return m2.IsValid() &&  m2.GetNcols == D1  ; }
      } ;
      // ======================================================================
      template <class T1,unsigned int D1, unsigned int D2, class R1,
                class T2>
      struct CanMul<TMatrixTSym<T2>,
                    ROOT::Math::SMatrix<T1,D1,D2,R1> >
      {
        static bool ok
        ( const TMatrixTSym<T2>&                      m2    ,
          const ROOT::Math::SMatrix<T1,D1,D2,R1>& /*  m1 */ )
        { return m2.IsValid() &&  m2.GetNcols == D1  ; }
      } ;
      // ======================================================================
      template <class T1,unsigned int D1, unsigned int D2, class R1,
                class T2>
      struct CanMul<TVectorT<T2>,
                    ROOT::Math::SMatrix<T1,D1,D2,R1> >
      {
        static bool ok
        ( const TVectorT<T2>&                      m2    ,
          const ROOT::Math::SMatrix<T1,D1,D2,R1>& /*  m1 */ )
        { return m2.IsValid() &&  m2.GetNrows == D1  ; }
      } ;
      // ======================================================================
      template <class T1,unsigned int D,
                class T2>
      struct CanMul<TVectorT<T2>,
                    ROOT::Math::SVector<T1,D> >
      {
        static bool ok
        ( const TVectorT<T2>&                  m2 , 
          const ROOT::Math::SVector<T1,D>& /*  m1 */ ) 
        { return m2.IsValid() &&  m2.GetNrows() == D ; }
      } ;
      // ======================================================================
      template <class T1,unsigned int D,
                class T2>
      struct CanMul<TMatrixT<T2>,
                    ROOT::Math::SVector<T1,D> >
      {
        static bool ok
        ( const TMatrixT<T2>&                  m2    ,
          const ROOT::Math::SVector<T1,D>& /*  m1 */ )
        { return m2.IsValid() && m2.GetNcols() == D ; }
      } ;
      // ======================================================================
      template <class T1,unsigned int D,
                class T2>
      struct CanMul< TMatrixTSym<T2>, 
                     ROOT::Math::SVector<T1,D> > 
              {
        static bool ok
        ( const TMatrixTSym<T2>&               m2    , 
          const ROOT::Math::SVector<T1,D>& /*  m1 */ )
        { return m2.IsValid() && m2.GetNcols() == D ; }
      } ;

      
      
      
      // ======================================================================
      template <class T1, class T2>
      struct CanIMul<TMatrixT<T1>, TMatrixT<T2> >
      {
        static bool ok
        ( const TMatrixT<T1>& m1 ,
          const TMatrixT<T2>& m2 )
        { return m1.IsValid() &&  m2.IsValid()
            && m1.GetNcols () == m2.GetNrows()
            && m2.GetNcols () == m2.GetNrows() ; }
      } ;      
      // ======================================================================
      template <class T1, class T2>
      struct CanIMul<TMatrixT<T1>, TMatrixTSym<T2> >
      {
        static bool ok
        ( const TMatrixT<T1>&    m1 ,
          const TMatrixTSym<T2>& m2 )
        { return m1.IsValid() &&  m2.IsValid() && m1.GetNcols () == m2.GetNrows() ; }
      } ;
      // ======================================================================
      template <class T1, unsigned int D, class R1,
                class T2>
      struct CanIMul<TMatrixT<T2>,
                     ROOT::Math::SMatrix<T1,D,D,R1> >
      {
        static bool ok
        ( const TMatrixT<T2>&                      m1    ,
          const ROOT::Math::SMatrix<T1,D,D,R1>& /* m2 */ )
        { return m1.IsValid() && m1.GetNcols () == D ; }
      } ;
      // ======================================================================
      template <class T1, unsigned int D1, unsigned D2,
                class T2>
      struct CanIMul<ROOT::Math::SMatrix<T1,D1,D2>, 
                     TMatrixT<T2> >
      {
        static bool ok
        ( const ROOT::Math::SMatrix<T1,D1,D2>& /* m2 */ ,
          const TMatrixT<T2>&                     m1    )
        { return m1.IsValid() && D2 == m1.GeNrows () ; }
      } ;

      
      template <class T>
      struct CanPow<TMatrixT<T> >
      { static bool ok ( const TMatrixT<T>& m1 )
        { return m1.IsValid () && m1.GetNrows () == m1.GetNcols () ; } } ;
      
      
      template <class T>
      struct CanPow<TMatrixTSym<T> >
      { static bool ok ( const TMatrixTSym<T>& m1 ) { return m1.IsValid () ; } } ;

      template <class T>
      struct CanSym<TMatrixT<T> >
      { static bool ok ( const TMatrixT<T>& m1 )
        { return m1.IsValid () && m1.GetNrows () == m1.GetNcols () ; } } ;

      template <class T>
      struct CanSym<TMatrixTSym<T> >
      { static bool ok ( const TMatrixTSym<T>& m1 ) { return m1.IsValid () ; } } ;
      

      template <class T>
      struct CanASym<TMatrixT<T> >
      { static bool ok ( const TMatrixT<T>& m1 )
        { return m1.IsValid () && m1.GetNrows () == m1.GetNcols () ; } } ;

      template <class T>
      struct CanASym<TMatrixTSym<T> >
      { static bool ok ( const TMatrixTSym<T>& m1 ) { return m1.IsValid () ; } } ;
      
      
      // ======================================================================
      // add
      // ======================================================================
      template <class T>
      struct Add<TMatrixT<T>,TMatrixT<T> >
      {
        typedef TMatrixT<T> R ;
        static  R add 
        ( const TMatrixT<T>& m1 , 
          const TMatrixT<T>& m2 ) { return m1 + m2 ; }
      } ;
      // ======================================================================
      template <class T>
      struct Add<TMatrixT<T>,TMatrixTSym<T> >
      {
        typedef TMatrixT<T> R ;
        static  R add 
        ( const TMatrixT<T>&    m1 , 
          const TMatrixTSym<T>& m2 ) { return m1 + m2 ; }
      } ;
      // ======================================================================
      template <class T>
      struct Add<TMatrixTSym<T>,TMatrixT<T> >
      {
        typedef TMatrixT<T> R ;
        static  R add 
        ( const TMatrixTSym<T>& m1 , 
          const TMatrixT<T>&   m2 ) { return m1 + m2 ; }
      } ;
      // ======================================================================
      template <class T>
      struct Add<TMatrixTSym<T>,TMatrixTSym<T> >
      {
        typedef TMatrixTSym<T> R ;
        static  R add 
        ( const TMatrixTSym<T>& m1 , 
          const TMatrixTSym<T>& m2 )
        { return TMatrixTSym<T> ( m1 , TMatrixTSym<T>::kPlus , m2 ); }
      } ;
      // ======================================================================
      template <class T,unsigned int D1, unsigned int D2, class R1>
      struct Add<ROOT::Math::SMatrix<T,D1,D2,R1>,TMatrixT<T> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M1 ;
        typedef TMatrixT<T>                     M2 ;
        typedef ROOT::Math::SMatrix<T,D1,D2>    R  ;
        typedef SM<R>                           C1 ;
        //
        static R add ( const M1& m1 , const M2& m2 )
        { return Add<M1,R>::add ( m1 , C1::transform ( m2 ) ) ; }
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
        static R add ( const M2& m2 , const M1& m1 )
        { return Add<M1,R>::add ( m1  , C1::transform ( m2 ) ) ; }
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
        static M1 add ( const M1& m1 , const M2& m2 )
        { return Add<M1,M1>::add ( m1 , C1::transform ( m2 ) ) ; }
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
        static M1 add ( const M2& m2 , const M1& m1 )
        { return Add<M1,M1>::add ( m1 , C1::transform ( m2 ) ) ; }
      } ;
      // ======================================================================      
      template <class T>
      struct Add<TVectorT<T>, TVectorT<T> > 
      {
        typedef TVectorT<T>  M ;
        static TVectorT<T> add ( const M& m1 , const M& m2 ) { return m1 + m2 ; }
      } ;
      // ======================================================================
      template <class T1,unsigned int D, class T2>
      struct Add<ROOT::Math::SVector<T1,D>, TVectorT<T2> > 
      {
        typedef ROOT::Math::SVector<T1,D>   M1 ;
        typedef TVectorT<T2>                M2 ;
        typedef SM<M1>                      C1 ;
        typedef M1                          R  ;
        //
        static R add ( const M1& m1 , const M2& m2 )
        { return Add<M1,M1>::add ( m1 , C1::transform ( m2 ) ) ; }
      } ;
      // ======================================================================
      template <class T1,unsigned int D,class T2>
      struct Add<TVectorT<T2>,ROOT::Math::SVector<T1,D> >
      {
        typedef ROOT::Math::SVector<T1,D>    M1 ;
        typedef TVectorT<T2>                 M2 ;
        typedef SM<M1>                       C1 ;
        typedef M1                           R  ;
        //
        static R add ( const M2& m2 , const M1& m1 )
        { return Add<M1,M1>::add ( C1::transform ( m2 ) , m1 ) ; }
      } ;
      // ======================================================================
      
      
      // ======================================================================
      // iadd
      // ======================================================================      
      template <class T,unsigned int D1,unsigned int D2, class R2>
      struct IAdd<TMatrixT<T>                     ,
                  ROOT::Math::SMatrix<T,D1,D2,R2> >
      {
        typedef TMatrixT<T>                     M1 ;
        typedef ROOT::Math::SMatrix<T,D1,D2,R2> M2 ;
        // addition
        static M1& iadd ( M1& m1 , const M2 & m2 )
        {
          for ( unsigned int i = 0 ; i < D1 ; ++i )
          { for ( unsigned int j = 0 ; j < D2 ; ++j )
            { m1 ( i , j ) += m2 ( i , j ) ; } }
          //
          return m1 ;
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
        static M1& iadd ( M1& m1 , const M2 & m2 )
        {
          for ( unsigned int i = 0 ; i < D ; ++i )
          { for ( unsigned int j = i ; j < D ; ++j )
            { m1 ( i , j ) += m2 ( i , j ) ; } }
          //
          return m1 ;
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
        static M1& iadd ( M1& m1 , const M2 & m2 )
        {
          for ( unsigned int i = 0 ; i < D1 ; ++i )
          { for ( unsigned int j = 0 ; j < D2 ; ++j )
            { m1 ( i , j ) += m2 ( i , j ) ; } }
          //
          return m1 ;
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
        static M1& iadd ( M1& m1 , const M2 & m2 )
        {
          for ( unsigned int i = 0 ; i < D1 ; ++i )
          { for ( unsigned int j = 0 ; j < D2 ; ++j )
            { m1 ( i , j ) += m2 ( i , j ) ; } }
          //
          return m1 ;
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
        static M1& iadd ( M1& m1 , const M2 & m2 )
        {
          for ( unsigned int i = 0 ; i < D ; ++i )
          { for ( unsigned int j = i ; j < D ; ++j )
            { m1 ( i , j ) += m2 ( i , j ) ; } }
          //
          return m1 ;
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
        static M1& iadd ( M1& m1 , const M2 & m2 )
        {
          //
          for ( unsigned int j = 0 ; j < D ; ++j )
          { m1 ( j ) += m2 ( j ) ; }
          //
          return m1 ;
        } 
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct IAdd<TVectorT<T>              ,
                  ROOT::Math::SVector<T,D> >
      {
        typedef TVectorT<T>              M1 ;
        typedef ROOT::Math::SVector<T,D> M2 ;
        // addition
        static M1& iadd ( M1& m1 , const M2 & m2 )
        {
          //
          for ( unsigned int j = 0 ; j < D ; ++j )
          { m1 ( j ) += m2 ( j ) ; }
          //
          return m1 ;
        } 
      } ;
      // ======================================================================
      
      // ======================================================================
      // sub
      // ======================================================================      
      template <class T>
      struct Sub<TMatrixT<T>,TMatrixT<T> >
      {
        static TMatrixT<T> sub 
        ( const TMatrixT<T>& m1 , 
          const TMatrixT<T>& m2 ) { return m1 - m2 ; }
      } ;
      // ======================================================================
      template <class T>
      struct Sub<TMatrixT<T>,TMatrixTSym<T> >
      {
        static TMatrixT<T> sub 
        ( const TMatrixT<T>&    m1 , 
          const TMatrixTSym<T>& m2 ) { return m1 - m2 ; }
      } ;
      // ======================================================================
      template <class T>
      struct Sub<TMatrixTSym<T>,TMatrixT<T> >
      {
        static TMatrixT<T> sub
        ( const TMatrixTSym<T>& m1 , 
          const TMatrixT<T>&    m2 ) { return m1 - m2 ; }
      } ;
      // ======================================================================
      template <class T>
      struct Sub<TMatrixTSym<T>,TMatrixTSym<T> >
      {
        static TMatrixTSym<T> sub
        ( const TMatrixTSym<T>& m1 , 
          const TMatrixTSym<T>& m2 )
        { return TMatrixTSym<T> ( m1 , TMatrixTSym<T>::kMinus , m2 ); }
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
        static R sub ( const M1& m1 , const M2& m2 )
        { return Sub<M1,R>::sub ( m1 , C1::transform ( m2 ) ) ; }
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
        static R sub ( const M2& m2 , const M1& m1 )
        { return Sub<R,M1>::sub ( C1::transform ( m2 ) , m1 ) ; }
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
        static M1 sub ( const M1& m1 , const M2& m2 )
        { return Sub<M1,R>::sub ( m1 , C1::transform ( m2 ) ) ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct Sub<TMatrixTSym<T>,ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef TMatrixTSym<T>                                         M2 ;
        typedef M1                                                     R  ;
        typedef SM<R>                                                 C1 ;
        //
        static M1 sub ( const M2& m2 , const M1& m1 )
        { return Sub<R,M1>::sub ( C1::transform ( m2 ) , m1 ) ; }
      } ;
      // ======================================================================      
      template <class T>
      struct Sub<TVectorT<T>, TVectorT<T> > 
      {
        typedef TVectorT<T>  M1 ;
        typedef TVectorT<T>  M2 ;
        typedef TVectorT<T>  R ;
        //
        static R sub ( const M1& m1 , const M1& m2 ) { return m1 - m2 ; }
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
        static R sub ( const M1& m1 , const M2& m2 )
        { return Sub<M1,M1>::sub ( m1 , C1::transform ( m2 ) ) ; }
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
        static R sub ( const M2& m2 , const M1& m1 )
        { return Sub<M1,M1>::sub ( C1::transform ( m2 ) , m1 ) ; }
      } ;
      // ======================================================================
      
      
      
      // ======================================================================
      // isub
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2, class R2>
      struct ISub<TMatrixT<T>                     ,
                  ROOT::Math::SMatrix<T,D1,D2,R2> >
      {
        typedef TMatrixT<T>                     M1 ;
        typedef ROOT::Math::SMatrix<T,D1,D2,R2> M2 ;
        // addition
        static M1& isub ( M1& m1 , const M2 & m2 )
        {
          for ( unsigned int i = 0 ; i < D1 ; ++i )
          { for ( unsigned int j = 0 ; j < D2 ; ++j )
            { m1 ( i , j ) -= m2 ( i , j ) ; } }
          //
          return m1 ;
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
        static M1& isub ( M1& m1 , const M2 & m2 )
        {
          for ( unsigned int i = 0 ; i < D ; ++i )
          { for ( unsigned int j = i ; j < D ; ++j )
            { m1 ( i , j ) -= m2 ( i , j ) ; } }
          //
          return m1 ;
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
        static M1& isub ( M1& m1 , const M2 & m2 )
        {
          for ( unsigned int i = 0 ; i < D1 ; ++i )
          { for ( unsigned int j = 0 ; j < D2 ; ++j )
            { m1 ( i , j ) -= m2 ( i , j ) ; } }
          //
          return m1 ;
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
        static M1& isub ( M1& m1 , const M2 & m2 )
        {
          for ( unsigned int i = 0 ; i < D1 ; ++i )
          { for ( unsigned int j = 0 ; j < D2 ; ++j )
            { m1 ( i , j ) -= m2 ( i , j ) ; } }
          //
          return m1 ;
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
        static M1& isub ( M1& m1 , const M2 & m2 )
        {
          for ( unsigned int i = 0 ; i < D ; ++i )
          { for ( unsigned int j = i ; j < D ; ++j )
            { m1 ( i , j ) -= m2 ( i , j ) ; } }
          //
          return m1 ;
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
        static M1& isub ( M1& m1 , const M2 & m2 )
        {
          //
          for ( unsigned int j = 0 ; j < D ; ++j ) { m1 ( j ) -= m2 ( j ) ; }
          //
          return m1 ;
        } 
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct ISub<TVectorT<T>              ,
                  ROOT::Math::SVector<T,D> >
      {
        typedef TVectorT<T>              M1 ;
        typedef ROOT::Math::SVector<T,D> M2 ;
        // addition
        static M1& isub ( M1& m1 , const M2 & m2 )
        {
          //
          for ( unsigned int j = 0 ; j < D ; ++j ) { m1 ( j ) -= m2 ( j ) ; }
          //
          return m1 ;
        } 
      } ;
      // ======================================================================
      
      // ======================================================================
      // mul
      // =======================================================================
      template <class T>
      struct Mul<TMatrixT<T>,TMatrixT<T> >
      {
        typedef TMatrixT<T> R ;
        static  R mul
        ( const TMatrixT<T>& m1 , 
          const TMatrixT<T>& m2 ) { return m1 * m2 ; }
      } ;
      // =======================================================================      
      template <class T>
      struct Mul<TMatrixT<T>,TMatrixTSym<T> >
      {
        typedef TMatrixT<T> R ;
        static  R mul
        ( const TMatrixT<T>&    m1 , 
          const TMatrixTSym<T>& m2 ) { return m1 * m2 ; }
      } ;
      // =======================================================================
      template <class T>
      struct Mul<TMatrixTSym<T>,TMatrixT<T> >
      {
        typedef TMatrixT<T> R ;
        static  R mul
        ( const TMatrixTSym<T>& m1 , 
          const TMatrixT<T>&    m2 ) { return m1 * m2 ; }
      } ;
      // =======================================================================      
      template <class T>
      struct Mul<TMatrixTSym<T>,TMatrixTSym<T> >
      {
        typedef TMatrixT<T> R ;
        static  R mul
        ( const TMatrixTSym<T>& m1 , 
          const TMatrixTSym<T>& m2 ) { return m1 * m2 ; }
      } ;
      // =======================================================================
      template <class T>
      struct Mul<TMatrixT<T>,TVectorT<T> >
      {
        typedef TVectorT<T> R ;
        static  R mul
        ( const TMatrixT<T>& m1 , 
          const TVectorT<T>& m2 ) { return m1 * m2 ; }
      } ;
      // =======================================================================      
      template <class T>
      struct Mul<TMatrixTSym<T>,TVectorT<T> >
      {
        typedef TVectorT<T> R ;
        static  R mul
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
        static  R mul
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
        static  R mul
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
        static  R mul
        ( const M1& m1 , 
          const M2& m2 )
        { return M::mul ( C::transform ( m1 )  , m2 ) ; }
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
        static  R mul
        ( const M1& m1 , 
          const M2& m2 )
        { return M::mul ( C::transform ( m1 )  , m2 ) ; }
      } ;
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
        static  R mul
        ( const M2& m2 , 
          const M1& m1 )
        { return M::mul ( m2 , C::transform ( m1 ) ) ; }
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
        static  R mul
        ( const M2& m2 , 
          const M1& m1 )
        { return M::mul ( m2 , C::transform ( m1 ) ) ; }
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
        static  R mul
        ( const M2& m2 , 
          const M1& m1 )
        { return M::mul ( m2 , C::transform ( m1 ) ) ; }
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
        static  R mul
        ( const M2& m2 , 
          const M1& m1 )
        { return M::mul ( m2 , C::transform ( m1 ) ) ; }
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
        static  R mul
        ( const M1& m2 , 
          const M2& m1 )
        { return M::mul ( m2 , C::transform ( m1 ) ) ; }
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
        static  R mul
        ( const M1& m2 , 
          const M2& m1 )
        { return M::mul ( m2 , C::transform ( m1 ) ) ; }
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
        static  R mul
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
        static M1& imult ( M1& m1 , const M2& m2 ) { m1 *= m2  ; return m1 ; }
      } ;      
      // ======================================================================
      template <class T>
      struct IMul<TMatrixT<T>, TMatrixTSym<T> >
      {
        typedef TMatrixT<T>    M1 ;
        typedef TMatrixTSym<T> M2 ;
        //
        static M1& imult ( M1& m1 , const M2& m2 ) { m1 *= m2  ; return m1 ; }
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
        static M1& imult ( M1& m1 , const M2& m2 ) { m1 *= C::transform ( m2 ) ; return m1 ; }
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
        static M1& imult ( M1& m1 , const M2& m2 ) { m1 *= C::transform ( m2 ) ; return m1 ; }
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
        static M1& imult ( M1& m1 , const M2& m2 ) { m1 *= C::transform ( m2 ) ; return m1 ; }
      } ;
      // ======================================================================

      
      
      // ======================================================================
      // EQUALITY
      // ======================================================================
      template <class T1, class T2>
      struct Eq <TVectorT<T1> ,
                 TVectorT<T2> >
      {
        typedef TVectorT<T1>  M1 ;
        typedef TVectorT<T2>  M2 ;
        typedef bool          R  ;
        // addition
        static R eq  ( const M1& m1 , const M2& m2 )
        {
          static const Ostap::Math::Equal_To<M1> s_cmp ;
          return s_cmp ( m1 , m2 ) ;
        }
      };
      // ======================================================================
      template <class T1, unsigned int D,class T2>
      struct Eq <TVectorT<T2>,
                 ROOT::Math::SVector<T1,D> >
      {
        typedef TVectorT<T2>               M1 ;
        typedef ROOT::Math::SVector<T1,D>  M2 ;
        typedef bool                       R  ;
        // addition
        static R eq  ( const M1& m1 , const M2& m2 )
        {
          static const Ostap::Math::Equal_To<M1> s_cmp ;
          return s_cmp ( m1 , m2 ) ;
        }
      };
      // ======================================================================
      template <class T1, unsigned int D,class T2>
      struct Eq <ROOT::Math::SVector<T1,D>,
                 TVectorT<T2> >
      {
        typedef ROOT::Math::SVector<T1,D>  M1 ;
        typedef TVectorT<T1>               M2 ;
        typedef bool                       R  ;
        // addition
        static R eq  ( const M1& m1 , const M2& m2 )
        {
          static const Ostap::Math::Equal_To<M2> s_cmp ;
          return s_cmp ( m2 , m1 ) ;
        }
      };
      // ============================================================================
      template <class T>
      struct Eq <TMatrixT<T> ,
                 TMatrixT<T> >
      {
        typedef TMatrixT<T>  M1 ;
        typedef TMatrixT<T>  M2 ;
        typedef bool         R  ;
        // addition
        static R eq  ( const M1& m1 , const M2& m2 )
        {
          static const Ostap::Math::Equal_To<M1> s_cmp ;
          return s_cmp ( m1  , m2 ) ;
        }
      };
      // ======================================================================
      template <class T>
      struct Eq <TMatrixT<T>    ,
                 TMatrixTSym<T> >
      {
        typedef TMatrixT<T>     M1 ;
        typedef TMatrixTSym<T>  M2 ;
        typedef bool             R  ;
        // addition
        static R eq  ( const M1& m1 , const M2& m2 )
        {
          static const Ostap::Math::Equal_To<M1> s_cmp ;
          return s_cmp ( m1  , m2 ) ;
        }
      };
      // ======================================================================      
      template <class T>
      struct Eq <TMatrixTSym<T> ,
                 TMatrixT<T>    >
      {
        typedef TMatrixTSym<T>  M1 ;
        typedef TMatrixT<T>     M2 ;
        typedef bool             R  ;
        // addition
        static R eq  ( const M1& m1 , const M2& m2 )
        {
          static const Ostap::Math::Equal_To<M2> s_cmp ;
          return s_cmp ( m2  , m1 ) ;
        }
      };
      // ======================================================================
      template <class T>
      struct Eq <TMatrixTSym<T> ,
                 TMatrixTSym<T> >
      {
        typedef TMatrixTSym<T>  M1 ;
        typedef TMatrixTSym<T>  M2 ;
        typedef bool             R  ;
        // addition
        static R eq  ( const M1& m1 , const M2& m2 )
        {
          static const Ostap::Math::Equal_To<M1> s_cmp ;
          return s_cmp ( m1  , m2 ) ;
        }
      };
      // ======================================================================
      template <class T1, unsigned int D1, unsigned int D2, class R1,
                class T2>
      struct Eq <ROOT::Math::SMatrix<T1,D1,D2,R1>,TMatrixT<T2> >
      {
        typedef ROOT::Math::SMatrix<T1,D1,D2,R1>  M1 ;
        typedef TMatrixT<T2>                      M2 ;
        typedef bool                              R  ;
        // addition
        static R eq  ( const M1& m1 , const M2& m2 )
        {
          static const Ostap::Math::Equal_To<M2> s_cmp ;
          return s_cmp ( m2 , m1 ) ;
        }
      };
      // ======================================================================
      template <class T1, unsigned int D, class R1,
                class T2>
      struct Eq <ROOT::Math::SMatrix<T1,D,D,R1>,TMatrixTSym<T2> >
      {
        typedef ROOT::Math::SMatrix<T1,D,D,R1>  M1 ;
        typedef TMatrixTSym<T2>                 M2 ;
        typedef bool                            R  ;
        // addition
        static R eq  ( const M1& m1 , const M2& m2 )
        {
          static const Ostap::Math::Equal_To<M2> s_cmp ;
          return s_cmp ( m2 , m1 ) ;
        }
      };
      // ======================================================================
      template <class T1, unsigned int D1, unsigned int D2, class R1,
                class T2>
      struct Eq <TMatrixT<T2> ,
                 ROOT::Math::SMatrix<T1,D1,D2,R1> >
      {
        typedef ROOT::Math::SMatrix<T1,D1,D2,R1>  M2 ;
        typedef TMatrixT<T2>                      M1 ;
        typedef bool                              R  ;
        // addition
        static R eq  ( const M1& m1 , const M2& m2 )
        {
          static const Ostap::Math::Equal_To<M1> s_cmp ;
          return s_cmp ( m1 , m2 ) ;
        }
      };
      // ======================================================================
      template <class T1, unsigned int D, class R1,
                class T2>
      struct Eq <TMatrixTSym<T2>,
                 ROOT::Math::SMatrix<T1,D,D,R1> >
      {
        typedef ROOT::Math::SMatrix<T1,D,D,R1>  M2 ;
        typedef TMatrixTSym<T2>                 M1 ;
        typedef bool                            R  ;
        // addition
        static R eq  ( const M1& m1 , const M2& m2 )
        {
          static const Ostap::Math::Equal_To<M1> s_cmp ;
          return s_cmp ( m1 , m2 ) ;
        }
      };
      // ======================================================================

      
      // ======================================================================
      template <class T>
      struct Pow<TMatrixT<T> >
      {
        //
        typedef TMatrixT<T> M ;
        typedef TMatrixT<T> R ;
        //
        static R pow ( const M& m , const unsigned short n )
        {
          //
          if      ( 0 == n ) { return M ( M::kUnit , m ) ; }
          else if ( 1 == n ) { return         m ; }
          else if ( 2 == n ) { return     m * m ; }
          else if ( 3 == n ) { return m * m * m ; }
          //
          R r = pow ( m , n / 2 )  ;
          if ( 0 == n / 2 ) { return r * r ; }
          //
          return r * r * m ;
        }
      } ;


      // ======================================================================
      template <class T>
      struct Pow<TMatrixTSym<T> >
      {
        //
        typedef TMatrixTSym<T> M ;
        typedef TMatrixT<T>    R ;
        //
        static R pow ( const M& m , const unsigned short n )
        {
          //
          if      ( 0 == n ) { return M ( M::kUnit , m ) ; }
          else if ( 1 == n ) { return         m ; }
          else if ( 2 == n ) { return     m * m ; }
          else if ( 3 == n ) { return m * m * m ; }
          //
          R r = pow ( m , n / 2 )  ;
          if ( 0 == n / 2 ) { return r * r ; }
          //
          return r * r * m ;
        }
      } ;
      
      
      // ======================================================================
      template <class T>
      struct Sym<TMatrixT<T> >
      {
        //
        typedef TMatrixT<T>    M ;
        typedef TMatrixTSym<T> R ;
        //
        static R sym ( const M& m )
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
        static R sym ( const M& m )
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
        static R asym ( const M& m )
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
        static R asym ( const M& m )
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
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                           The end of namespace Ostap::Math
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_MATRIXUTILST_H
// ============================================================================
