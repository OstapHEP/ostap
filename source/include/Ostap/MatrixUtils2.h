// ============================================================================
#ifndef OSTAP_MATRIXUTILS2_H
#define OSTAP_MATRIXUTILS2_H 1
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
#include "Math/SMatrix.h"
#include "Math/SVector.h"
#include "Math/Point3D.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Math.h"
#include "Ostap/EigenSystem.h"
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
    namespace  Ops
    {
      // ======================================================================
      // CHECKERS
      // ======================================================================
      
      template <class M1, class M2>
      struct CanAdd   { static bool operation ( const M1& /* m1 */ , const M2& /* m2 */ ) { return false ; } } ;
      
      template <class M1, class M2>
      struct CanMul   { static bool operation ( const M1& /* m1 */ , const M2& /* m2 */ ) { return false ; } } ;
      template <class M1, class M2>
      struct CanIMul  { static bool operation ( const M1& /* m1 */ , const M2& /* m2 */ ) { return false ; } } ;
      
      template <class M1, class M2>
      struct CanDiv   { static bool operation ( const M1& /* m1 */ , const M2& /* m2 */ ) { return false ; } } ;
      template <class M1, class M2>
      struct CanIDiv  { static bool operation ( const M1& /* m1 */ , const M2& /* m2 */ ) { return false ; } } ;
      
      template <class M1, class M2>
      struct CanDot   { static bool operation ( const M1& /* m1 */ , const M2& /* m2 */ ) { return false ; } } ;
      template <class M1, class M2>
      struct CanCross { static bool operation ( const M1& /* m1 */ , const M2& /* m2 */ ) { return false ; } } ;
      template <class M1, class M2>
      struct CanSim   { static bool operation ( const M1& /* m1 */ , const M2& /* m2 */ ) { return false ; } } ;
      template <class M1, class M2>
      struct CanSimT  { static bool operation ( const M1& /* m1 */ , const M2& /* m2 */ ) { return false ; } } ;

      template <class M1, class M2>
      struct CanEq    { static bool operation ( const M1& /* m1 */ , const M2& /* m2 */ ) { return false ; } } ;
      
      template <class M1>
      struct CanPow   { static bool operation ( const M1& /* m1 */ , const double /* p */ ) { return false ; } } ;

      template <class M1>
      struct CanSym   { static bool operation ( const M1& /* m1 */ ) { return false ; } } ;

      template <class M1>
      struct CanASym  { static bool operation ( const M1& /* m1 */ ) { return false ; } } ;

      template <class M1>
      struct CanInvert { static bool operation ( const M1& /* m1 */ ) { return false ; } } ;

      template <class M1, class M2>
      struct CanRMul
      { static bool operation ( const M1& m1 , const M2& m2 ) { return CanMul<M2,M1>::operation ( m2 , m1 ) ; } } ;
      
      // ======================================================================
      // partial specializations with scalar/double 
      // ======================================================================
      template <class M1>
      struct CanMul<M1,double>  
      { static bool operation ( const M1&    /* m1 */ , const double /* m2 */ ) { return true ; } } ;
      template <class M1>
      struct CanMul<double,M1> 
      { static bool operation ( const double /* m2 */ , const M1&    /* m1 */ ) { return true ; } } ;
      template <class M1>
      struct CanRMul<M1,double> 
      { static bool operation ( const M1&    /* m1 */ , const double /* m2 */ ) { return true ; } } ;
      template <class M1>
      struct CanIMul<M1,double> 
      { static bool operation ( const M1&    /* m1 */ , const double /* m2 */ ) { return true ; } } ;      
      template <class M1>
      struct CanDiv<M1,double> 
      { static bool operation ( const M1&    /* m1 */ , const double /* m2 */ ) { return true ; } } ;
      template <class M1>
      struct CanIDiv<M1,double> 
      { static bool operation ( const M1&    /* m1 */ , const double /* m2 */ ) { return true ; } } ;
      
      
      template <class T, unsigned int D, class R1>
      struct CanInvert< ROOT::Math::SMatrix<T,D,D,R1> > 
      {
        typedef ROOT::Math::SMatrix<T,D,D,R1> M1 ;
        typedef bool                          R  ;
        //
        static R operation ( const M1& /* m1 */ ) { return true ; } 
      } ;
      
      
      // ======================================================================
      // new cases  with "almost" scalar
      // ======================================================================      
      template <class M1, class T, class R1>
      struct CanMul<M1, ROOT::Math::SMatrix<T,1,1,R1> > 
      {
        typedef ROOT::Math::SMatrix<T,1,1,R1> M2 ;
        static bool operation ( const M1& /* m1 */ , const M2& /* m2 */ ) { return true ; } 
      } ;
      // ======================================================================
      template <class M1, class T, class R1>
      struct CanRMul<M1, ROOT::Math::SMatrix<T,1,1,R1> > 
      {
        typedef ROOT::Math::SMatrix<T,1,1,R1> M2 ;
        static bool operation ( const M1& /* m1 */ , const M2& /* m2 */ ) { return true ; } 
      } ;
      // ======================================================================
      template <class M1, class T, class R1>
      struct CanIMul<M1, ROOT::Math::SMatrix<T,1,1,R1> > 
      {
        typedef ROOT::Math::SMatrix<T,1,1,R1> M2 ;
        static bool operation ( const M1& /* m1 */ , const M2& /* m2 */ ) { return true ; } 
      } ;
      // ======================================================================
      template <class M1, class T, class R1>
      struct CanMul<ROOT::Math::SMatrix<T,1,1,R1> , M1 > 
      {
        typedef ROOT::Math::SMatrix<T,1,1,R1> M2 ;
        static bool operation ( const M2& /* m1 */ , 
                                const M1& /* m2 */ ) { return true ; } 
      } ;      
      // ======================================================================
      template <class M1, class T, class R1>
      struct CanRMul<ROOT::Math::SMatrix<T,1,1,R1> , M1 > 
      {
        typedef ROOT::Math::SMatrix<T,1,1,R1> M2 ;
        static bool operation ( const M2& /* m1 */ , const M1& /* m2 */ ) { return true ; } 
      } ;
      // ======================================================================

      
      // ======================================================================
      // Can be added 
      // ======================================================================
      
      template <class T,unsigned int D1,unsigned int D2,class R1,
                class R2>
      struct CanAdd<ROOT::Math::SMatrix<T,D1,D2,R1> ,
                    ROOT::Math::SMatrix<T,D1,D2,R2> >
      {
        static bool operation
        ( const ROOT::Math::SMatrix<T,D1,D2,R1>& /* m1 */ , 
          const ROOT::Math::SMatrix<T,D1,D2,R2>& /* m2 */ ) { return true ; }  
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct CanAdd<ROOT::Math::SVector<T,D> ,
                    ROOT::Math::SVector<T,D> >
      {
        static bool operation
        ( const ROOT::Math::SVector<T,D>& /* m1 */ , 
          const ROOT::Math::SVector<T,D>& /* m2 */ ) { return true ; }  
      } ;
      // ======================================================================
      template <class T,unsigned int D,class R1>
      struct CanAdd<ROOT::Math::SMatrix<T,D,D,R1> , double>
      {
        static bool operation
        ( const ROOT::Math::SMatrix<T,D,D,R1>& /* m1 */ , 
          const double                         /* m2 */ ) { return true ; }  
      } ;
      
      // ======================================================================
      // Can be multiplied
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2,class R1,
                unsigned int D3,class R2>
      struct CanMul<ROOT::Math::SMatrix<T,D1,D2,R1> ,
                    ROOT::Math::SMatrix<T,D2,D3,R2> >
      {
        static bool operation
        ( const ROOT::Math::SMatrix<T,D1,D2,R1>& /* m1 */ , 
          const ROOT::Math::SMatrix<T,D2,D3,R2>& /* m2 */ ) { return true ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2,class R1>
      struct CanMul<ROOT::Math::SMatrix<T,D1,D2,R1> ,
                    ROOT::Math::SVector<T,D2> >
      {
        static bool operation
        ( const ROOT::Math::SMatrix<T,D1,D2,R1>& /* m1 */ , 
          const ROOT::Math::SVector<T,D2>&       /* m2 */ ) { return true ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2,class R1>
      struct CanMul<ROOT::Math::SVector<T,D1>       , 
                    ROOT::Math::SMatrix<T,D1,D2,R1> >
      {
        static bool operation
        ( const ROOT::Math::SVector<T,D1>&       /* m2 */ ,
          const ROOT::Math::SMatrix<T,D1,D2,R1>& /* m1 */ )  { return true ; } 
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct CanMul<ROOT::Math::SVector<T,D>       , 
                    ROOT::Math::SVector<T,D> >
      {
        static bool operation
        ( const ROOT::Math::SVector<T,D>& /* m2 */ ,
          const ROOT::Math::SVector<T,D>& /* m1 */ )  { return true ; } 
      } ;
      // ========================================================================
      
      // ========================================================================
      template <class T,unsigned int D1, unsigned int D2,
                class R2>
      struct CanIMul<ROOT::Math::SMatrix<T,D1,D2>    ,
                     ROOT::Math::SMatrix<T,D2,D2,R2> >
      {
        static bool operation
        ( const ROOT::Math::SMatrix<T,D1,D2>&    /* m1 */ , 
          const ROOT::Math::SMatrix<T,D2,D2,R2>& /* m2 */ ) { return true ; }
      } ;
      // ======================================================================
      
      // ======================================================================
      template <class T,unsigned int D>
      struct CanDot<ROOT::Math::SVector<T,D> , 
                    ROOT::Math::SVector<T,D> >
      {
        static bool operation 
        ( const ROOT::Math::SVector<T,D>& /* m2 */ ,
          const ROOT::Math::SVector<T,D>& /* m1 */ )  { return true ; } 
      } ;
      // ======================================================================
      
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2>
      struct CanCross<ROOT::Math::SVector<T,D1> , 
                      ROOT::Math::SVector<T,D2> >
      {
        static bool operation
        ( const ROOT::Math::SVector<T,D1>& /* m2 */ ,
          const ROOT::Math::SVector<T,D2>& /* m1 */ )  { return true ; } 
      } ;
      // ======================================================================
      
      // ======================================================================
      template <class T,unsigned int D,unsigned int D2, class R2>
      struct CanSim<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > ,
                    ROOT::Math::SMatrix<T,D2,D,R2> >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef ROOT::Math::SMatrix<T,D2,D,R2>                         M2 ;
        // check
        static bool operation
        ( const M1&  m1 ,
          const M2&  m2 ) { return true ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct CanSim<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D>> ,
                    ROOT::Math::SVector<T,D> >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef ROOT::Math::SVector<T,D>                               M2 ;
        // check 
        static bool operation 
          ( const M1& /* m1 */ , const M2& /* m2 */ ) { return true ; }
      } ;
      // =====================================================================
      template <class T,unsigned int D> 
      struct CanSim<ROOT::Math::SVector<T,D> ,
                    ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > >
      {
        typedef ROOT::Math::SVector<T,D>                               M1 ;
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M2 ;
        // check 
        static bool operation 
        ( const M1& /* m1 */ , const M2& /* m2 */ ) { return true ; }
      } ;
      // ======================================================================
      
      // ======================================================================
      template <class T,unsigned int D,unsigned int D2, class R2> 
      struct CanSimT<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > ,
                     ROOT::Math::SMatrix<T,D,D2,R2> >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >    M1 ;
        typedef ROOT::Math::SMatrix<T,D,D2,R2>                            M2 ;
        // check 
        static bool operation ( const M1& /* m1 */ , const M2& /* m2 */ ) { return true ; }
      } ;
      // ======================================================================
      
      template <class T, unsigned int D, class R1>
      struct CanPow<ROOT::Math::SMatrix<T,D,D,R1> >
      { 
        static bool operation 
        ( const ROOT::Math::SMatrix<T,D,D,R1>& /* m1 */ ,
          const unsigned short                 /* p  */ ) { return true ; } 
      } ; 

      template <class T, unsigned int D, class R1>
      struct CanSym<ROOT::Math::SMatrix<T,D,D,R1> >
      { static bool operation ( const ROOT::Math::SMatrix<T,D,D,R1>& /* m1 */ ) { return true ; } } ;
      
      template <class T, unsigned int D, class R1>
      struct CanASym<ROOT::Math::SMatrix<T,D,D,R1> >
      { static bool operation ( const ROOT::Math::SMatrix<T,D,D,R1>& /* m1 */ ) { return true ; } } ;

      // ======================================================================
      // Operations
      // ======================================================================
      template <class T1,class T2>
      struct  Add  ;
      template <class T1,class T2>
      struct IAdd ;
      template <class T1,class T2>
      struct  Sub  ;
      template <class T1,class T2>
      struct ISub ;
      template <class T1,class T2>
      struct  Mul  ;
      template <class T1,class T2>
      struct IMul ;
      template <class T1,class T2>
      struct  Div  ;
      template <class T1,class T2>
      struct IDiv ;
      
      template <class T1,class T2>
      struct RAdd ;
      template <class T1,class T2>
      struct RSub ;
      template <class T1,class T2>
      struct RDiv ;
      
      template <class T1,class T2>
      struct Dot   ;
      template <class T1,class T2>
      struct Cross ;
      template <class T1,class T2>
      struct Sim   ;
      template <class T1,class T2>
      struct SimT  ;


      template <class T1,class T2>
      struct Eq ;

      template <class T1>
      struct Pow ;

      template <class T1>
      struct Sym  ;
      
      template <class T1>
      struct ASym ;

      template <class T1>
      struct Invert ;

      // ======================================================================
      // "Right" operations
      // ======================================================================
      template <class M1, class M2>
      struct RAdd
      {
        typedef Add<M2,M1>    A ;
        typedef typename A::R R ;
        static R operation ( const M1& m1 , const  M2& m2 ) { return A::operation ( m2 , m1 ) ; }
      } ;
      // ======================================================================
      template <class M1, class M2>
      struct RSub
      {
        typedef Sub<M2,M1>    S ;
        typedef typename S::R R ;
        static R operation ( const M1& m1 , const  M2& m2 ) { return S::operation ( m2 , m1 ) ; }
      } ;
      // ======================================================================
      template <class M1, class M2>
      struct RMul
      {
        typedef Mul<M2,M1>    M ;
        typedef typename M::R R ;
        static R operation ( const M1& m1 , const  M2& m2 ) { return M::operation ( m2 , m1 ) ; }
      } ;
      
      // ======================================================================
      // scaling
      // ======================================================================
      template <class M1>
      struct IMul<M1,double>
      { static void operation ( M1& m1 , const double m2 ) 
        {  m1 *= m2 ; } } ;
      // ======================================================================
      template <class M1>
      struct IDiv<M1,double>
      { static void operation ( M1& m1 , const double m2 )
        { IMul<M1,double>::operation ( m1 , 1 / m2 ) ; } } ;
      // ======================================================================
      template <class M1>
      struct Mul<M1,double>
      {
        typedef M1 R ;
        static R operation ( const M1& m1 , const double m2 ) { return m1 * m2  ; }
      } ;
      // ======================================================================
      template <class M1>
      struct RMul<M1,double>
      {
        typedef M1 R ;
        static R operation ( const M1& m1 , const double m2 ) 
        { return Mul<M1,double>::operation ( m1 , m2 ) ; }
      } ;
      // ======================================================================
      template <class M1>
      struct Div<M1,double>
      {
        typedef M1 R ;
        static R operation ( const M1& m1 , const double m2 )
        { return Mul<M1,double>::operation ( m1 , 1 / m2 ) ; }
      } ;
      // ======================================================================

      // ======================================================================
      // Trivia 
      // ======================================================================
      template <class M>
      struct IAdd<M,M>     
      { static void operation ( M& m1 , const M& m2 ) { m1 += m2 ; } } ;
      // ======================================================================
      template <class M>
      struct ISub<M,M>
      { static void operation ( M& m1 , const M& m2 ) { m1 -= m2 ; } } ;
      

      // ======================================================================
      /// "almost constants"
      // ======================================================================
      template <class M1, class T, class R1>
      struct Mul<M1, ROOT::Math::SMatrix<T,1,1,R1> > 
      {
        typedef ROOT::Math::SMatrix<T,1,1,R1> M2 ;
        typedef Mul<M1,double>                O  ;
        typedef typename O::R                 R  ;
        //
        static R operation ( const M1& m1 , const M2& m2 )   
        { return O::operation (  m1 , m2 ( 0 , 0 ) ) ; }
      } ;
      // ======================================================================
      template <class M1, class T, class R1>
      struct RMul<M1, ROOT::Math::SMatrix<T,1,1,R1> > 
      {
        typedef ROOT::Math::SMatrix<T,1,1,R1> M2 ;
        typedef Mul<M1,double>                O  ;
        typedef typename O::R                 R  ;
        //
        static R operation ( const M1& m1 , const M2& m2 )
        { return O::operation ( m1 , m2 ( 0 , 0 ) ) ; }
      } ;      
      // ======================================================================
      template <class M1, class T, class R1>
      struct IMul<M1, ROOT::Math::SMatrix<T,1,1,R1> > 
      {
        typedef ROOT::Math::SMatrix<T,1,1,R1> M2 ;
        typedef IMul<M1,double>               O  ;
        static void operation ( M1& m1 , const M2& m2 ) { O::operation ( m1 , m2 ( 0 , 0 ) ) ; }
      } ;
      // ======================================================================
      template <class M1, class T, class R1>
      struct Mul<ROOT::Math::SMatrix<T,1,1,R1>, M1 > 
      {
        typedef ROOT::Math::SMatrix<T,1,1,R1> M2 ;
        typedef Mul<M1,double>                O  ;
        typedef typename O::R                 R  ;
        //
        static R operation ( const M2& m2 , const M1& m1 ) 
        { return O::operation ( m1 , m2 ( 0 , 0 ) ) ; } 
      } ;
      // ======================================================================
      template <class M1, class T, class R1>
      struct RMul<ROOT::Math::SMatrix<T,1,1,R1>, M1> 
      {
        typedef ROOT::Math::SMatrix<T,1,1,R1> M2 ;
        typedef RMul<M1,double>               O  ;
        typedef typename O::R                 R  ;
        //
        static R operation ( const M2& m2 , const M1& m1 ) 
        { return O::operation ( m1 , m2 ( 0 , 0 ) ) ; } 
      } ;
      // ======================================================================

      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2>
      struct IAdd<ROOT::Math::SMatrix<T,D1,D2> ,
                  ROOT::Math::SMatrix<T,D1,D2> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2> M ;
        static void operation ( M& m1 , const M& m2 ) { m1 += m2 ; }
      }; 
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2>
      struct ISub<ROOT::Math::SMatrix<T,D1,D2> ,
                  ROOT::Math::SMatrix<T,D1,D2> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2> M ;
        static void operation ( M& m1 , const M& m2 ) { m1 -= m2 ; }
      }; 
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2,class R>
      struct IAdd<ROOT::Math::SMatrix<T,D1,D2,R> ,
                  ROOT::Math::SMatrix<T,D1,D2,R> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R> M ;
        static void operation ( M& m1 , const M& m2 ) { m1 += m2 ; }
      }; 
      // =====================================================================
      template <class T,unsigned int D1,unsigned int D2,class R>
      struct ISub<ROOT::Math::SMatrix<T,D1,D2,R> ,
                  ROOT::Math::SMatrix<T,D1,D2,R> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R> M ;
        static void operation ( M& m1 , const M& m2 ) { m1 -= m2 ; }
      };
      // ======================================================================
      template <class T,unsigned int D>
      struct IAdd<ROOT::Math::SVector<T,D> ,
                  ROOT::Math::SVector<T,D> >
      {
        typedef ROOT::Math::SVector<T,D> M ;
        static void operation ( M& m1 , const M& m2 ) { m1 += m2 ; }
      }; 
      // ======================================================================
      template <class T,unsigned int D>
      struct ISub<ROOT::Math::SVector<T,D> ,
                  ROOT::Math::SVector<T,D> >
      {
        typedef ROOT::Math::SVector<T,D> M ;
        static void operation ( M& m1 , const M& m2 ) { m1 -= m2 ; }
      };

      // ======================================================================
      // ADD
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2,class R1,
                class R2>
      struct Add<ROOT::Math::SMatrix<T,D1,D2,R1> ,
                 ROOT::Math::SMatrix<T,D1,D2,R2> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M1 ;
        typedef ROOT::Math::SMatrix<T,D1,D2,R2> M2 ;
        typedef ROOT::Math::SMatrix<T,D1,D2>    R  ;
        // addition
        static R operation ( const M1& m1 , const M2& m2 ) { return m1 + m2 ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct Add<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > ,
                 ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M2 ;
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > R  ;
        // addition
        static R operation ( const M1 & m1 , const M2 & m2 ) { return m1 + m2 ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct Add<ROOT::Math::SVector<T,D> ,
                 ROOT::Math::SVector<T,D> >
      {
        typedef ROOT::Math::SVector<T,D> M1 ;
        typedef ROOT::Math::SVector<T,D> M2 ;
        typedef ROOT::Math::SVector<T,D> R  ;
        // addition
        static R operation ( const M1 & m1 , const M2 & m2 ) { return m1 + m2 ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D,class R1>
      struct Add<ROOT::Math::SMatrix<T,D,D,R1>, double> 
      {
        typedef ROOT::Math::SMatrix<T,D,D,R1> M1 ;
        typedef double                        M2 ;
        typedef ROOT::Math::SMatrix<T,D,D,R1> R  ;
        // addition
        static R operation ( const M1& m1 , const double m2 ) 
        { 
          R result  { m1 } ;
          for  ( unsigned int i = 0 ; i < D ; ++i ) { result ( i , i )  += m2 ; }
          return result ; 
        }
      } ;
      // ======================================================================
      template <class T,unsigned int D,class R1>
      struct RAdd<ROOT::Math::SMatrix<T,D,D,R1>, double> 
        : public Add<ROOT::Math::SMatrix<T,D,D,R1>, double> {} ;
      

      // ======================================================================
      // IADD
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2, class R2>
      struct IAdd<ROOT::Math::SMatrix<T,D1,D2>    ,
                  ROOT::Math::SMatrix<T,D1,D2,R2> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2>    M1 ;
        typedef ROOT::Math::SMatrix<T,D1,D2,R2> M2 ;
        // addition
        static void operation ( M1 & m1 , const M2 & m2 ) { m1 += m2 ; }
      } ;
      // ======================================================================      
      template <class T,unsigned int D>
      struct IAdd<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > ,
                 ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > >
      {
       typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
       typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M2 ;
       // addition
        static void operation ( M1 & m1 , const M2 & m2 ) { m1 += m2 ; }
      } ;
      // ======================================================================      
      template <class T,unsigned int D,class R1>
      struct IAdd<ROOT::Math::SMatrix<T,D,D,R1> , double> 
      {
        typedef ROOT::Math::SMatrix<T,D,D,R1> M1 ;
        typedef double                        M2 ;
        // addition
        static void operation ( M1& m1 , const double m2 ) 
        { for  ( unsigned int i = 0 ; i < D ; ++i ) { m1 ( i , i )  += m2 ; } }
      } ;

      
      // ======================================================================
      // SUB
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2,class R1,
                class R2>
      struct Sub<ROOT::Math::SMatrix<T,D1,D2,R1> ,
                 ROOT::Math::SMatrix<T,D1,D2,R2> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M1 ;
        typedef ROOT::Math::SMatrix<T,D1,D2,R2> M2 ;
        typedef ROOT::Math::SMatrix<T,D1,D2>    R  ;
        // subtraction 
        static R operation ( const M1& m1 , const M2& m2 ) { return m1 - m2 ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct Sub<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > ,
                 ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M2 ;
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > R  ;
        // subtraction 
        static R operation ( const M1 & m1 , const M2 & m2 ) { return m1 - m2 ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct Sub<ROOT::Math::SVector<T,D> ,
                 ROOT::Math::SVector<T,D> >
      {
        typedef ROOT::Math::SVector<T,D> M1 ;
        typedef ROOT::Math::SVector<T,D> M2 ;
        typedef ROOT::Math::SVector<T,D> R  ;
        // subtraction
        static R operation ( const M1 & m1 , const M2 & m2 ) { return m1 - m2 ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D,class R1>
      struct Sub<ROOT::Math::SMatrix<T,D,D,R1> , double> 
      {
        typedef ROOT::Math::SMatrix<T,D,D,R1> M1 ;
        typedef double                        M2 ;
        typedef ROOT::Math::SMatrix<T,D,D,R1> R  ;
        // addition
        static R operation ( const M1& m1 , const double m2 ) 
        { 
          R result  { m1 } ;
          for  ( unsigned int i = 0 ; i < D ; ++i ) { result ( i , i ) -= m2 ; }
          return result ; 
        }
      } ;
      // ======================================================================
      template <class T,unsigned int D,class R1>
      struct RSub<ROOT::Math::SMatrix<T,D,D,R1> , double> 
      {
        typedef ROOT::Math::SMatrix<T,D,D,R1> M1 ;
        typedef double                        M2 ;
        typedef ROOT::Math::SMatrix<T,D,D,R1> R  ;
        // addition
        static R operation ( const M1& m1 , const double m2 ) 
        { 
          R result  { m1 } ; result *= -1 ; // ATTENTION! 
          for  ( unsigned int i = 0 ; i < D ; ++i ) { result ( i , i ) += m2 ; }
          return result ; 
        }
      } ;
      
      // ======================================================================
      // ISUB
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2, class R2>
      struct ISub<ROOT::Math::SMatrix<T,D1,D2>    ,
                  ROOT::Math::SMatrix<T,D1,D2,R2> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2>    M1 ;
        typedef ROOT::Math::SMatrix<T,D1,D2,R2> M2 ;
        // subtraction
        static void operation ( M1 & m1 , const M2 & m2 ) { m1 -= m2 ; }
      } ;
      // ======================================================================

      template <class T,unsigned int D>
      struct ISub<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > ,
                  ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M2 ;
        // addition
        static void operation ( M1 & m1 , const M2 & m2 ) { m1 -= m2 ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D,class R1>
      struct ISub<ROOT::Math::SMatrix<T,D,D,R1> , double> 
      {
        typedef ROOT::Math::SMatrix<T,D,D,R1> M1 ;
        typedef double                        M2 ;
        // addition
        static void operation ( M1& m1 , const double m2 ) 
        { for  ( unsigned int i = 0 ; i < D ; ++i ) { m1 ( i , i ) -= m2 ; } }
      } ;

      
      // ======================================================================
      // MUL 
      // ======================================================================      
      template <class T,unsigned int D1,unsigned int D2,class R1,
                unsigned int D3,class R2>
      struct Mul<ROOT::Math::SMatrix<T,D1,D2,R1> ,
                 ROOT::Math::SMatrix<T,D2,D3,R2> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M1 ;
        typedef ROOT::Math::SMatrix<T,D2,D3,R2> M2 ;
        typedef ROOT::Math::SMatrix<T,D1,D3>     R  ;
        // multiplication
        static R operation ( const M1 & m1 , const M2 & m2 ) { return m1 * m2 ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2,class R1>
      struct Mul<ROOT::Math::SMatrix<T,D1,D2,R1> ,
                 ROOT::Math::SVector<T,D2>       >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M1 ;
        typedef ROOT::Math::SVector<T,D2>       M2 ;
        typedef ROOT::Math::SVector<T,D1>       R  ;
        // multiplication
        static R operation ( const M1 & m1 , const M2 & m2 ) { return m1 * m2 ; }
      } ;
      // ===================================================================
      template <class T,unsigned int D1,unsigned int D2,class R1>
      struct Mul<ROOT::Math::SVector<T,D1>       ,
                 ROOT::Math::SMatrix<T,D1,D2,R1> >
      {
        typedef ROOT::Math::SVector<T,D1>       M1 ;
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M2 ;
        typedef ROOT::Math::SVector<T,D2>       R  ;
        // multiplication
        static R operation ( const M1& m1 , const M2 & m2 ) { return m1 * m2 ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2,class R1>
      struct Mul<ROOT::Math::SMatrix<T,D1,D2,R1>,double>
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M1 ;
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> R  ;
        // multiplication
        static R operation ( const M1& m1 , const double m2 ) { return m1 * m2  ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2,class R1>
      struct Mul<double,ROOT::Math::SMatrix<T,D1,D2,R1>>
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M1 ;
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> R  ;
        // multiplication
        static R operation ( const double m2 , const M1& m1) { return m1 * m2  ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct Mul<ROOT::Math::SVector<T,D> ,
                 ROOT::Math::SVector<T,D> >
      {
        typedef ROOT::Math::SVector<T,D>       M1 ;
        typedef ROOT::Math::SVector<T,D>       M2 ;
        typedef double                         R  ;
        // multiplication
        static double operation ( const M1 & m1 , const M2 & m2 )
        { return std::inner_product ( m1.begin() , m1.end() , m2.begin() , 0.0 ); }
      } ;
      // ===============================================================
      template <class T,unsigned int D>
      struct Mul<ROOT::Math::SVector<T,D>,double>
      {
        typedef ROOT::Math::SVector<T,D> M1 ;
        typedef ROOT::Math::SVector<T,D> R  ;
        // multiplication
        static R operation ( const M1& m1 , const double m2 ) { return m1 * m2  ; }
      } ;
      // ===============================================================      
      template <class T,unsigned int D>
      struct Mul<double,ROOT::Math::SVector<T,D>>
      {
        typedef ROOT::Math::SVector<T,D> M1 ;
        typedef ROOT::Math::SVector<T,D> R  ;
        // multiplication
        static R operation ( const double m2,  const M1& m1 ) { return m1 * m2  ; }
      } ;
      // ======================================================================

      
      // ======================================================================
      // IMUL 
      // ======================================================================      
      template <class T,unsigned int D1,unsigned int D2,
                class R2>
      struct IMul<ROOT::Math::SMatrix<T,D1,D2>    ,
                  ROOT::Math::SMatrix<T,D2,D2,R2> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2>    M1 ;
        typedef ROOT::Math::SMatrix<T,D2,D2,R2> M2 ;
        // in-place multiplication
        static void operation ( M1& m1 , const M2&m2 ) { m1 *= m2 ; }
      } ;
      // ======================================================================


      // ======================================================================
      // DOT
      // ======================================================================
      template <class T1,unsigned int D, class T2>
      struct Dot<ROOT::Math::SVector<T1,D> ,
                 ROOT::Math::SVector<T2,D> >
      {
        typedef ROOT::Math::SVector<T1,D> M1 ;
        typedef ROOT::Math::SVector<T2,D> M2 ;
        typedef double                    R  ;
        // multiplication
        static double operation ( const M1 & m1 , const M2 & m2 )
        { return std::inner_product ( m1.begin() , m1.end() , m2.begin() , 0.0 ) ; }
      } ;            
      // ======================================================================


      // ======================================================================
      // CROSS
      // ======================================================================
      template <class T,unsigned int D1,unsigned D2>
      struct Cross<ROOT::Math::SVector<T,D1> ,
                   ROOT::Math::SVector<T,D2> >
      {
        typedef ROOT::Math::SVector<T,D1>       M1 ;
        typedef ROOT::Math::SVector<T,D2>       M2 ;
        typedef ROOT::Math::SMatrix<T,D1,D2>    R  ;
        // multiplication
        static R operation ( const M1 & m1 , const M2 & m2 )
        {
          R r ;
          for ( unsigned int i = 0 ; i < D1 ; ++i ) 
          { for ( unsigned int j = 0 ; j < D2 ; ++j ) 
            { r(i,j) = m1[i] * m2[j] ; } }
          return r;
        }
      } ;      
      // ======================================================================

      
      // ======================================================================
      // SIM 
      // ======================================================================
      template <class T,unsigned int D,
                unsigned int D2, class R2> 
      struct Sim<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > ,
                 ROOT::Math::SMatrix<T,D2,D,R2> >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >    M1 ;
        typedef ROOT::Math::SMatrix<T,D2,D,R2>                            M2 ;
        typedef ROOT::Math::SMatrix<T,D2,D2,ROOT::Math::MatRepSym<T,D2> > R  ;        
        // similarity
        static R operation ( const M1& A , const M2& U  )
        { return ROOT::Math::Similarity ( U , A ) ; }
      } ;
      // ======================================================================      
      template <class T,unsigned int D> 
      struct Sim<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > ,
                 ROOT::Math::SVector<T,D> >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef ROOT::Math::SVector<T,D>                               M2 ;
        typedef double                                                 R  ;
        // similarity 
        static double operation ( const M1& A , const M2& V )
        { return ROOT::Math::Similarity ( A , V ) ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D> 
      struct Sim<ROOT::Math::SVector<T,D> ,
                 ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > >
      {
        typedef ROOT::Math::SVector<T,D>                               M1 ;
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M2 ;
        typedef double                                                 R  ;
        // check 
        static double operation ( const M1& V , const M2& A  )
        { return ROOT::Math::Similarity ( A , V ) ; }
      } ;
      // ======================================================================
      
      // ======================================================================
      // SIMT 
      // ======================================================================
      template <class T,unsigned int D,
                unsigned int D2, class R2> 
      struct SimT<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > ,
                  ROOT::Math::SMatrix<T,D,D2,R2> >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >    M1 ;
        typedef ROOT::Math::SMatrix<T,D,D2,R2>                            M2 ;
        typedef ROOT::Math::SMatrix<T,D2,D2,ROOT::Math::MatRepSym<T,D2> > R  ;        
        // similarity
        static R operation ( const M1& A , const M2& U  )
        { return ROOT::Math::SimilarityT ( U , A ) ; }
      } ;
      // ======================================================================

      
      // ======================================================================
      template <class T, unsigned int D, class R1>
      struct Invert<ROOT::Math::SMatrix<T,D,D,R1> > 
      {
        typedef ROOT::Math::SMatrix<T,D,D,R1> M1 ;
        typedef ROOT::Math::SMatrix<T,D,D,R1> R  ;
        //
        static R operation ( const M1& m1 , int& flag  )
        { return m1.Inverse ( flag ) ; }
      } ;      
      // ======================================================================
      
      // ======================================================================
      /// can be compared ?
      // ======================================================================
      template <class T1, unsigned D, class T2> 
      struct CanEq < ROOT::Math::SVector<T1,D> , 
                     ROOT::Math::SVector<T2,D> >
      {
        typedef ROOT::Math::SVector<T1,D> M1 ;
        typedef ROOT::Math::SVector<T2,D> M2 ;
        // equality 
        static bool operation 
        ( const M1& /* m1 */ , 
          const M2& /* m2 */ ) { return true ; }
      } ;
      // ======================================================================
      template <class T, unsigned D> 
      struct CanEq < ROOT::Math::SVector<T,D> , 
                     ROOT::Math::SVector<T,D> >
      {
        typedef ROOT::Math::SVector<T,D> M1 ;
        typedef ROOT::Math::SVector<T,D> M2 ;
        // equality 
        static bool operation 
        ( const M1& /* m1 */ , 
          const M2& /* m2 */ ) { return true ; }
      } ;
      // ======================================================================
      template <class T, unsigned D1, unsigned D2, class R1> 
      struct CanEq < ROOT::Math::SMatrix<T,D1,D2,R1> , 
                     ROOT::Math::SMatrix<T,D1,D2,R1> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M1 ;
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M2 ;
        // equality 
        static bool operation 
        ( const M1& /* m1 */ , 
          const M2& /* m2 */ ) { return true ; }
      } ;
      // ======================================================================
      template <class T, unsigned D1, unsigned D2, class R1, class R2> 
      struct CanEq < ROOT::Math::SMatrix<T,D1,D2,R1> , 
                     ROOT::Math::SMatrix<T,D1,D2,R2> >
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M1 ;
        typedef ROOT::Math::SMatrix<T,D1,D2,R2> M2 ;
        // equality 
        static bool operation 
        ( const M1& /* m1 */ , 
          const M2& /* m2 */ ) { return true ; }
      } ;
      // ======================================================================
      template <class T1, unsigned D1, unsigned D2, class R1, class T2, class R2> 
      struct CanEq < ROOT::Math::SMatrix<T1,D1,D2,R1> , 
                     ROOT::Math::SMatrix<T2,D1,D2,R2> >
      {
        typedef ROOT::Math::SMatrix<T1,D1,D2,R1> M1 ;
        typedef ROOT::Math::SMatrix<T2,D1,D2,R2> M2 ;
        // equality 
        static bool operation 
        ( const M1& /* m1 */ , 
          const M2& /* m2 */ ) { return true ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D,class R1>
      struct CanEq<ROOT::Math::SMatrix<T,D,D,R1> , double >
      {
        typedef ROOT::Math::SMatrix<T,D,D,R1> M1 ;
        typedef double                        M2 ;
        typedef bool                          R  ;
        // 
        static R operation ( const M1& /* m1 */ , const M2 /* m2 */ ) { return true ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct CanEq<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > , double >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >  M1 ;
        typedef double                                                  M2 ;
        typedef bool                                                    R  ;
        // 
        static R operation ( const M1& m1 , const M2 m2 ) { return true ; }
      } ;

      // ======================================================================
      // Equality
      // ======================================================================
      template <class T1,unsigned int D1,unsigned int D2,class R1,
                class T2, class R2>
      struct Eq<ROOT::Math::SMatrix<T1,D1,D2,R1> ,
                ROOT::Math::SMatrix<T2,D1,D2,R2> >
      {
        typedef ROOT::Math::SMatrix<T1,D1,D2,R1> M1 ;
        typedef ROOT::Math::SMatrix<T2,D1,D2,R2> M2 ;
        typedef bool                      R  ;
        // addition
        static R operation ( const M1& m1 , const M2& m2 )
        {
          static const Ostap::Math::Equal_To<M1> s_cmp ;
          return s_cmp ( m1 ,  m2 ) ;
        }
      } ;
      // ======================================================================
      template <class T1,unsigned int D,
                class T2>
      struct Eq<ROOT::Math::SVector<T1,D> ,
                ROOT::Math::SVector<T2,D> >
      {
        typedef ROOT::Math::SVector<T1,D> M1 ;
        typedef ROOT::Math::SVector<T2,D> M2 ;
        typedef bool                      R  ;
        // addition
        static R operation ( const M1& m1 , const M2& m2 )
        {
          static const Ostap::Math::Equal_To<M1> s_cmp ;
          return s_cmp ( m1 ,  m2 ) ;
        }
      } ;
      // ======================================================================

      
      template <class T,unsigned int D,class R1>
      struct Eq<ROOT::Math::SMatrix<T,D,D,R1> , double >
      {
        typedef ROOT::Math::SMatrix<T,D,D,R1> M1 ;
        typedef double                        M2 ;
        typedef bool                          R  ;
        // 
        static R operation ( const M1& m1 , const M2 m2 )
        {
          static const Ostap::Math::Equal_To<T>  s_cmp  ;
          static const Ostap::Math::Zero<T>      s_zero ;
          for ( unsigned int i = 0 ; i < D ; ++i ) 
          {
            if ( !s_cmp  ( m1 ( i , i ) , m2 )        ) { return false ; }            
            for ( unsigned int j = 0 ; j < D ; ++j ) 
            { 
              if ( i != j && !s_zero ( m1 ( i , j ) ) ) { return false ; }
            }
          }
          return true ;
        }
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct Eq<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > , double >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef double                                                 M2 ;
        typedef bool                                                   R  ;
        //
        static R operation ( const M1& m1 , const M2 m2 )
        {
          static const Ostap::Math::Equal_To<T>  s_cmp  ;
          static const Ostap::Math::Zero<T>      s_zero ;
          for ( unsigned int i = 0 ; i < D ; ++i ) 
          {
            if ( !s_cmp  ( m1 ( i , i ) , m2 )   ) { return false ; }            
            for ( unsigned int j = i +1  ; j < D ; ++j ) 
            {
              if ( !s_zero ( m1 ( i , j )      ) ) { return false ; }
            }
          }
          return true ;
        }
      } ;

      // =======================================================================
      // EXTRA
      // =======================================================================
      template <class T>
      struct Eigen ;
      
      
      template <class T,unsigned int D> 
      struct Eigen<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef ROOT::Math::SVector<T,D>                               M2 ;
        typedef ROOT::Math::SMatrix<T,D,D>                             M3 ;
        // get eigen values 
        static Ostap::StatusCode operation 
        ( const M1&  m             , 
          M2&        values        ,
          const bool sorted = true )
        {
          Ostap::Math::GSL::EigenSystem eigen {} ;
          return eigen.eigenValues  ( m , values , sorted ) ;
        }
        // get eigen values and eigenvectors 
        static Ostap::StatusCode operation 
        ( const M1&  m                 , 
          M2&        values            , 
          M3&        vectors           ,
          const bool sorted     = true ,
          const bool ascending  = true )
        {
          Ostap::Math::GSL::EigenSystem eigen {} ;
          return eigen.eigenVectors ( m , values , vectors , sorted , ascending ) ;
        }
      } ;
      
      // ======================================================================
      template <class T, unsigned int D, class R1>
      struct Pow<ROOT::Math::SMatrix<T,D,D,R1> >
      {
        //
        typedef ROOT::Math::SMatrix<T,D,D,R1> M ;
        typedef ROOT::Math::SMatrix<T,D,D>    R ;
        //
        static R operation ( const M& m , const unsigned short n )
        {
          //
          if      ( 0 == n ) { return M ( ROOT::Math::SMatrixIdentity () ) ; }
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
      template <class T, class R1>
      struct Pow<ROOT::Math::SMatrix<T,1,1,R1> >
      {
        //
        typedef ROOT::Math::SMatrix<T,1,1,R1> M ;
        typedef double   R ;
        //
        static R operation ( const M& m , const int    n )
        { return 0 == n ? 1.0 : std::pow ( m ( 0 , 0 )  , n ) ; }
      } ;
      // ======================================================================
      template <class T, unsigned int D>
      struct Sym<ROOT::Math::SMatrix<T,D,D> >
      {
        //
        typedef ROOT::Math::SMatrix<T,D,D>                             M ;
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > R ;
        //
        static R operation ( const M& m )
        {
          //
          R r ;
          for ( unsigned int i = 0 ; i < D ; ++i )
          { r ( i , i ) = m ( i , i ) ;
            for ( unsigned int j = i + 1  ; j < D ; ++j )
            { r ( i , j ) = 0.5 * ( m ( i , j ) + m ( j , i ) ) ; } }
          //
          return r ;
        }
      } ;
      // ======================================================================
      template <class T, unsigned int D>
      struct Sym<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > >
      {
        //
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M ;
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > R ;
        //
        static const R& operation ( const M& m ) { return m ; }
      } ;
      // ======================================================================
      template <class T, unsigned int D>
      struct ASym<ROOT::Math::SMatrix<T,D,D> >
      {
        //
        typedef ROOT::Math::SMatrix<T,D,D> M ;
        typedef ROOT::Math::SMatrix<T,D,D> R ;
        //
        static R operation ( const M& m )
        {
          //
          R r ;
          for ( unsigned int i = 0 ; i < D ; ++i )
          {
            r ( i , i ) = 0 ;
            for ( unsigned int j = i + 1 ; j < D ; ++j )
            {
              const T v = 0.5 * ( m ( i , j ) - m ( j , i ) ) ;
              r ( i , j ) =  v ;
              r ( j , i ) = -v ;          
            } 
          }
          //
          return r ;
        }
      } ;
      // ======================================================================
      template <class T, unsigned int D>
      struct ASym<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > >
      {
        //
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M ;
        typedef ROOT::Math::SMatrix<T,D,D> R ;
        //
        static R operation ( const M& /* m */ ) { return R () ; }  
      } ;      
      // ======================================================================
    } //                                  The end of namespace Ostap::Math::Ops
    // ========================================================================
  } //                                        The end of namespace  Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // LHCBMATH_MATRIXUTILS2_H
// ============================================================================
