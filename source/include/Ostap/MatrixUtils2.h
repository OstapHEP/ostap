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
      struct CanAdd   { static bool ok ( const M1& /* m1 */ , const M2& /* m2 */ ) { return false ; } } ;
      
      template <class M1, class M2>
      struct CanMul   { static bool ok ( const M1& /* m1 */ , const M2& /* m2 */ ) { return false ; } } ;
      template <class M1, class M2>
      struct CanIMul  { static bool ok ( const M1& /* m1 */ , const M2& /* m2 */ ) { return false ; } } ;
      
      template <class M1, class M2>
      struct CanDiv   { static bool ok ( const M1& /* m1 */ , const M2& /* m2 */ ) { return false ; } } ;
      template <class M1, class M2>
      struct CanIDiv  { static bool ok ( const M1& /* m1 */ , const M2& /* m2 */ ) { return false ; } } ;
      
      template <class M1, class M2>
      struct CanDot   { static bool ok ( const M1& /* m1 */ , const M2& /* m2 */ ) { return false ; } } ;
      template <class M1, class M2>
      struct CanCross { static bool ok ( const M1& /* m1 */ , const M2& /* m2 */ ) { return false ; } } ;
      template <class M1, class M2>
      struct CanSim   { static bool ok ( const M1& /* m1 */ , const M2& /* m2 */ ) { return false ; } } ;
      template <class M1, class M2>
      struct CanSimT  { static bool ok ( const M1& /* m1 */ , const M2& /* m2 */ ) { return false ; } } ;

      template <class M1>
      struct CanPow   { static bool ok ( const M1& /* m1 */ ) { return false ; } } ;
      
      template <class M1>
      struct CanSym   { static bool ok ( const M1& /* m1 */ ) { return false ; } } ;
      template <class M1>
      struct CanASym  { static bool ok ( const M1& /* m1 */ ) { return false ; } } ;
      
      template <class M1, class M2>
      struct CanRMul
      { static bool ok ( const M1& m1 , const M2& m2 ) { return CanMul<M2,M1>::ok ( m2 , m1 ) ; } } ;
      
      // ======================================================================
      // partial specializations with scalar/double 
      // ======================================================================
      template <class M1>
      struct CanMul<M1,double>  { static bool ok ( const M1&    /* m1 */ , const double /* m2 */ ) { return true ; } } ;
      template <class M1>
      struct CanMul<double,M1>  { static bool ok ( const double /* m2 */ , const M1&    /* m1 */ ) { return true ; } } ;
      template <class M1>
      struct CanRMul<M1,double> { static bool ok ( const M1&    /* m1 */ , const double /* m2 */ ) { return true ; } } ;
      template <class M1>
      struct CanIMul<M1,double> { static bool ok ( const M1&    /* m1 */ , const double /* m2 */ ) { return true ; } } ;      
      template <class M1>
      struct CanDiv<M1,double>  { static bool ok ( const M1&    /* m1 */ , const double /* m2 */ ) { return true ; } } ;
      template <class M1>
      struct CanIDiv<M1,double> { static bool ok ( const M1&    /* m1 */ , const double /* m2 */ ) { return true ; } } ;
      
      // ======================================================================
      // Can be added 
      // ======================================================================
      
      template <class T,unsigned int D1,unsigned int D2,class R1,
                class R2>
      struct CanAdd<ROOT::Math::SMatrix<T,D1,D2,R1> ,
                    ROOT::Math::SMatrix<T,D1,D2,R2> >
      {
        static bool ok
        ( const ROOT::Math::SMatrix<T,D1,D2,R1>& /* m1 */ , 
          const ROOT::Math::SMatrix<T,D1,D2,R2>& /* m2 */ ) { return true ; }  
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct CanAdd<ROOT::Math::SVector<T,D> ,
                    ROOT::Math::SVector<T,D> >
      {
        static bool ok
        ( const ROOT::Math::SVector<T,D>& /* m1 */ , 
          const ROOT::Math::SVector<T,D>& /* m2 */ ) { return true ; }  
      } ;
      // ======================================================================
      
      // ======================================================================
      // Can be multiplied
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2,class R1,
                unsigned int D3,class R2>
      struct CanMul<ROOT::Math::SMatrix<T,D1,D2,R1> ,
                    ROOT::Math::SMatrix<T,D2,D3,R2> >
      {
        static bool ok
        ( const ROOT::Math::SMatrix<T,D1,D2,R1>& /* m1 */ , 
          const ROOT::Math::SMatrix<T,D2,D3,R2>& /* m2 */ ) { return true ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2,class R1>
      struct CanMul<ROOT::Math::SMatrix<T,D1,D2,R1> ,
                    ROOT::Math::SVector<T,D2> >
      {
        static bool ok
        ( const ROOT::Math::SMatrix<T,D1,D2,R1>& /* m1 */ , 
          const ROOT::Math::SVector<T,D2>&       /* m2 */ ) { return true ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2,class R1>
      struct CanMul<ROOT::Math::SVector<T,D1>       , 
                    ROOT::Math::SMatrix<T,D1,D2,R1> >
      {
        static bool ok
        ( const ROOT::Math::SVector<T,D1>&       /* m2 */ ,
          const ROOT::Math::SMatrix<T,D1,D2,R1>& /* m1 */ )  { return true ; } 
      } ;
      // ======================================================================
      template <class T,unsigned int D>
      struct CanMul<ROOT::Math::SVector<T,D>       , 
                    ROOT::Math::SVector<T,D> >
      {
        static bool ok
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
        static bool ok ( const ROOT::Math::SMatrix<T,D1,D2>&    /* m1 */ , 
                         const ROOT::Math::SMatrix<T,D2,D2,R2>& /* m2 */ ) { return true ; }
      } ;
      // ======================================================================
      
      // ======================================================================
      template <class T,unsigned int D>
      struct CanDot<ROOT::Math::SVector<T,D> , 
                    ROOT::Math::SVector<T,D> >
      {
        static bool ok
        ( const ROOT::Math::SVector<T,D>& /* m2 */ ,
          const ROOT::Math::SVector<T,D>& /* m1 */ )  { return true ; } 
      } ;
      // ======================================================================
      
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2>
      struct CanCross<ROOT::Math::SVector<T,D1> , 
                      ROOT::Math::SVector<T,D2> >
      {
        static bool ok
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
        static bool ok
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
        static bool ok ( const M1& /* m1 */ , const M2& /* m2 */ ) { return true ; }
      } ;
      // =====================================================================
      template <class T,unsigned int D> 
      struct CanSim<ROOT::Math::SVector<T,D> ,
                    ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > >
      {
        typedef ROOT::Math::SVector<T,D>                               M1 ;
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M2 ;
        // check 
        static bool ok ( const M1& /* m1 */ , const M2& /* m2 */ ) { return true ; }
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
        static bool ok ( const M1& /* m1 */ , const M2& /* m2 */ ) { return true ; }
      } ;
      // ======================================================================


      template <class T, unsigned int D, class R1>
      struct CanPow<ROOT::Math::SMatrix<T,D,D,R1> >
      { static bool ok ( const ROOT::Math::SMatrix<T,D,D,R1>& /* m1 */ ) { return true ; } } ;
      
      template <class T, unsigned int D, class R1>
      struct CanSym<ROOT::Math::SMatrix<T,D,D,R1> >
      { static bool ok ( const ROOT::Math::SMatrix<T,D,D,R1>& /* m1 */ ) { return true ; } } ;
      
      template <class T, unsigned int D, class R1>
      struct CanASym<ROOT::Math::SMatrix<T,D,D,R1> >
      { static bool ok ( const ROOT::Math::SMatrix<T,D,D,R1>& /* m1 */ ) { return true ; } } ;
      
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

      // ======================================================================
      // "Right" operations
      // ======================================================================
      template <class M1, class M2>
      struct RAdd
      {
        typedef Add<M2,M1>    A ;
        typedef typename A::R R ;
        static R radd ( const M1& m1 , const  M2& m2 ) { return A::add ( m2 , m1 ) ; }
      } ;
      // ======================================================================
      template <class M1, class M2>
      struct RSub
      {
        typedef Sub<M2,M1>    S ;
        typedef typename S::R R ;
        static R rsub ( const M1& m1 , const  M2& m2 ) { return S::sub ( m2 , m1 ) ; }
      } ;
      // ======================================================================
      template <class M1, class M2>
      struct RMul
      {
        typedef Mul<M2,M1>    M ;
        typedef typename M::R R ;
        static R rmul ( const M1& m1 , const  M2& m2 ) { return M::mul ( m2 , m1 ) ; }
      } ;
      
      // ======================================================================
      // scaling
      // ======================================================================
      template <class M1>
      struct IMul<M1,double>
      { static void imul ( M1& m1 , const double m2 ) {  m1 *= m2 ; } } ;
      // ======================================================================
      template <class M1>
      struct IDiv<M1,double>
      { static void  idiv ( M1& m1 , const double m2 ) { IMul<M1,double>::imul ( m1 , 1 / m2 ) ; } } ;
      // ======================================================================
      template <class M1>
      struct Mul<M1,double>
      {
        typedef M1 R ;
        static R mul ( const M1& m1 , const double m2 ) { return m1 * m2  ; }
      } ;
      // ======================================================================
      template <class M1>
      struct RMul<M1,double>
      {
        typedef M1 R ;
        static R rmul ( const M1& m1 , const double m2 ) { return Mul<M1,double>::mul ( m1 , m2 ) ; }
      } ;
      // ======================================================================
      template <class M1>
      struct Div<M1,double>
      {
        typedef M1 R ;
        static R div ( const M1& m1 , const double m2 ) { return Mul<M1,double>::mul ( m1 , 1 / m2 ) ; }
      } ;
      // ======================================================================
      template <class M>
      struct IAdd<M,M>     
      { static void iadd ( const M& m1 , const M& m2 ) { m1 += m2 ; } } ;
      // ======================================================================
      template <class M>
      struct ISub<M,M>
      { static void isub ( const M& m1 , const M& m2 ) { m1 -= m2 ; } } ;
      
      
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
        static R add ( const M1& m1 , const M2& m2 ) { return m1 + m2 ; }
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
        static R add ( const M1 & m1 , const M2 & m2 ) { return m1 + m2 ; }
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
        static R add ( const M1 & m1 , const M2 & m2 ) { return m1 + m2 ; }
      } ;
      // ======================================================================
      


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
        static void iadd ( M1 & m1 , const M2 & m2 ) { m1 += m2 ; }
      } ;
      // ======================================================================
      
      template <class T,unsigned int D>
      struct IAdd<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > ,
                 ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > >
      {
       typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
       typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M2 ;
       // addition
        static void iadd ( M1 & m1 , const M2 & m2 ) { m1 += m2 ; }
      } ;
      // ======================================================================
      
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
        static R sub ( const M1& m1 , const M2& m2 ) { return m1 - m2 ; }
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
        static R sub ( const M1 & m1 , const M2 & m2 ) { return m1 - m2 ; }
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
        static R sub ( const M1 & m1 , const M2 & m2 ) { return m1 - m2 ; }
      } ;
      // ======================================================================
      
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
        static void isub ( M1 & m1 , const M2 & m2 ) { m1 -= m2 ; }
      } ;
      // ======================================================================

      template <class T,unsigned int D>
      struct ISub<ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > ,
                  ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > >
      {
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M1 ;
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > M2 ;
        // addition
        static void isub ( M1 & m1 , const M2 & m2 ) { m1 -= m2 ; }
      } ;
      // ======================================================================
      
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
        static R mul ( const M1 & m1 , const M2 & m2 ) { return m1 * m2 ; }
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
        static R mul ( const M1 & m1 , const M2 & m2 ) { return m1 * m2 ; }
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
        static R mul ( const M1& m1 , const M2 & m2 ) { return m1 * m2 ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2,class R1>
      struct Mul<ROOT::Math::SMatrix<T,D1,D2,R1>,double>
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M1 ;
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> R  ;
        // multiplication
        static R mul ( const M1& m1 , const double m2 ) { return m1 * m2  ; }
      } ;
      // ======================================================================
      template <class T,unsigned int D1,unsigned int D2,class R1>
      struct Mul<double,ROOT::Math::SMatrix<T,D1,D2,R1>>
      {
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> M1 ;
        typedef ROOT::Math::SMatrix<T,D1,D2,R1> R  ;
        // multiplication
        static R mul ( const double m2 , const M1& m1) { return m1 * m2  ; }
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
        static double mul ( const M1 & m1 , const M2 & m2 )
        { return std::inner_product ( m1.begin() , m1.end() , m2.begin() , 0.0 ); }
      } ;
      // ===============================================================
      template <class T,unsigned int D>
      struct Mul<ROOT::Math::SVector<T,D>,double>
      {
        typedef ROOT::Math::SVector<T,D> M1 ;
        typedef ROOT::Math::SVector<T,D> R  ;
        // multiplication
        static R mul ( const M1& m1 , const double m2 ) { return m1 * m2  ; }
      } ;
      // ===============================================================      
      template <class T,unsigned int D>
      struct Mul<double,ROOT::Math::SVector<T,D>>
      {
        typedef ROOT::Math::SVector<T,D> M1 ;
        typedef ROOT::Math::SVector<T,D> R  ;
        // multiplication
        static R mul ( const double m2,  const M1& m1 ) { return m1 * m2  ; }
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
        static void imul ( M1& m1 , const M2&m2 ) { m1 *= m2 ; }
      } ;
      // ======================================================================

      
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
        static R eq  ( const M1& m1 , const M2& m2 )
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
        static R eq  ( const M1& m1 , const M2& m2 )
        {
          static const Ostap::Math::Equal_To<M1> s_cmp ;
          return s_cmp ( m1 ,  m2 ) ;
        }
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
        static double dot ( const M1 & m1 , const M2 & m2 )
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
        static R cross ( const M1 & m1 , const M2 & m2 )
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
        static R sim ( const M1& A , const M2& U  )
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
        static double sim ( const M1& A , const M2& V )
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
        static double sim ( const M1& V , const M2& A  )
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
        static R simt ( const M1& A , const M2& U  )
        { return ROOT::Math::SimilarityT ( U , A ) ; }
      } ;
      // ======================================================================




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
        // eigen values 
        static Ostap::StatusCode values  ( const M1& m , M2& v  , const bool sorted = true )
        {
          Ostap::Math::GSL::EigenSystem eigen {} ;
          return eigen.eigenValues  ( m , v      , sorted ) ;
        }
        // eigen vectors 
        static Ostap::StatusCode vectors ( const M1& m , M2& v  , M3& vs , const bool sorted = true )
        {
          Ostap::Math::GSL::EigenSystem eigen {} ;
          return eigen.eigenVectors ( m , v , vs , sorted ) ;
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
        static R pow ( const M& m , const unsigned short n )
        {
          //
          if      ( 0 == n ) { return M ( ROOT::Math::SMatrixIdentity () ) ; }
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
      template <class T, unsigned int D>
      struct Sym<ROOT::Math::SMatrix<T,D,D> >
      {
        //
        typedef ROOT::Math::SMatrix<T,D,D>                             M ;
        typedef ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> > R ;
        //
        static R sym ( const M& m )
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
        static R sym ( const M& m ) { return m ; }
      } ;
      
      // ======================================================================
      template <class T, unsigned int D>
      struct ASym<ROOT::Math::SMatrix<T,D,D> >
      {
        //
        typedef ROOT::Math::SMatrix<T,D,D> M ;
        typedef ROOT::Math::SMatrix<T,D,D> R ;
        //
        static R asym ( const M& m )
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
        static R asym ( const M& /* m */ ) { return R () ; }  
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
