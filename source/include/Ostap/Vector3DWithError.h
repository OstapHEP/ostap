// ============================================================================
#ifndef LHCBMATH_VECTOR3DWITHERROR_H 
#define LHCBMATH_VECTOR3DWITHERROR_H 1
// ============================================================================
// Include files 
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/SymmetricMatrixTypes.h"
#include "Ostap/GenericVectorTypes.h"
#include "Ostap/Vector3DTypes.h"
#include "Ostap/Point3DTypes.h"
// ============================================================================
#include "Ostap/Point3DWithError.h"
#include "Ostap/SVectorWithError.h"
// ============================================================================
/** @file Ostap/Vector3DWithError.h
 *  Collection of useful objects with associated "covarinaces".
 *  The concept has been stollen from Wouter Hulsbergen's lines 
 *  @author Vanya BELYAEV Ivane.B2lyaev@nikhef.nl
 *  @date 20090603
 */
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    /// forward declarations 
    class Point3DWithError  ;
    class Vector3DWithError ;
    // ========================================================================
    class Vector3DWithError : public Ostap::XYZVector 
    {
    public:
      // ======================================================================
      /// the actual type of the Vector 
      typedef Ostap::XYZVector                                     Vector3D   ;
      typedef Ostap::XYZVector                                     Value      ;
      /// the actual type of covariance 
      typedef Ostap::SymMatrix3x3                                  Covariance ;
      // ======================================================================
    public:
      // ======================================================================
      /// the actual type of generic 3-vector 
      typedef Ostap::Vector3                                       Vector     ;
      /// the actual type of vector with errors  
      typedef Ostap::Math::SVectorWithError<3,double>              VectorE    ;
      // ======================================================================
    public:
      // ======================================================================
      /// the actual type of the point 
      typedef Ostap::XYZPoint                                      Point3D    ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from 3d-vector and covariance matrix 
      Vector3DWithError 
      ( const Vector3D&   vct    = Vector3D   () , 
        const Covariance& matrix = Covariance () ) ;
      // ======================================================================
      /// constructor from 3d-vector and covariance matrix 
      Vector3DWithError 
      ( const Covariance& matrix                 , 
        const Vector3D&   vct    = Vector3D   () ) ;
      // ======================================================================
      /// constructor from generic vector and covariance matrix 
      Vector3DWithError 
      ( const Vector&     vct                    , 
        const Covariance& matrix = Covariance () ) ;
      // ======================================================================
      /// constructor from generic vector with errors 
      Vector3DWithError 
      ( const VectorE&    vct                    ) ;
      // ======================================================================
    public: // trivial accessors 
      // ======================================================================
      const  Vector3D&   vector3d   ( ) const { return *this       ; }
      const  Vector3D&   vector3D   ( ) const { return vector3d () ; }
      const  Covariance& covariance ( ) const { return cov2  ()    ; }      
      inline Vector3D&   vector3d   ( )       { return *this       ; }
      // ======================================================================      
      const  Vector3D&   value      ( ) const { return vector3d () ; }
      const  Covariance& cov2       ( ) const { return m_cov2      ; }
      // ======================================================================      
    public: // setters 
      // ======================================================================
      void setVector3d   ( const Vector3D&   v ) { vector3d() =  v   ; }
      void setVector3D   ( const Vector3D&   v ) { setVector3d ( v ) ; }
      void setVector     ( const Vector3D&   v ) { setVector3d ( v ) ; }
      // ======================================================================
      void setValue      ( const Vector3D&   v ) { setVector3d ( v ) ; }
      void setCovariance ( const Covariance& c ) { m_cov2      = c   ; }      
      // ======================================================================
      void setValue      ( const VectorE&    v ) ;
      void setValue      ( const Vector&     v ) ;
      // ======================================================================
    public: // finally it is just a point + covariance 
      // ======================================================================
      operator const Covariance& () const { return cov2  () ; }        
      operator       Covariance& ()       { return m_cov2   ; }        
      // ======================================================================
    public: // useful accessors to covarinace matrix 
      // ======================================================================
      /// access to elemens of covariance matrix 
      double cov2 ( unsigned int i , unsigned int j ) const 
      { return m_cov2 ( i , j ) ; }
      // ======================================================================
    public: // operators 
      // ======================================================================
      Vector3DWithError& operator+= ( const Vector3DWithError& right ) ;
      Vector3DWithError& operator+= ( const Vector3D&          right ) ;
      Vector3DWithError& operator+= ( const VectorE&           right ) ;
      Vector3DWithError& operator+= ( const Vector&            right ) ;
      Vector3DWithError& operator-= ( const Vector3DWithError& right ) ;
      Vector3DWithError& operator-= ( const Vector3D&          right ) ;
      Vector3DWithError& operator-= ( const VectorE&           right ) ;
      Vector3DWithError& operator-= ( const Vector&            right ) ;
      // ======================================================================
    public: // scaling
      // ======================================================================
      /// *= 
      Vector3DWithError& operator*= ( const double v ) ;                  // *= 
      /// /= 
      Vector3DWithError& operator/= ( const double v ) ;                  // /= 
      // ======================================================================
    public: // unary -
      // ======================================================================
      Vector3DWithError operator-() const ;                           // unary-
      // ======================================================================
    public:
      // ======================================================================
      /// get generic vector 
      void    asVector ( Vector&  data ) const ;
      /// get generic vector 
      void    asVector ( VectorE& data ) const ;
      /// convert to generic vector with errors:
      VectorE asVector () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate chi2 distance 
      double chi2 ( const Vector3DWithError& right )  const ;
      /// evaluate chi2 distance 
      double chi2 ( const Vector3D&          right )  const ;      
      /// evaluate chi2 distance 
      double chi2 ( const VectorE&           right )  const ;      
      /// evaluate chi2 distance 
      double chi2 ( const Vector&            right )  const ;      
      // ======================================================================
   public: // more functions 
      // ======================================================================
      Vector3DWithError  mean ( const Vector3DWithError& right ) const ;
      Vector3DWithError  mean ( const VectorE&           right ) const ;
      // ======================================================================
    public: // helper operators for Python
      // ======================================================================
      Vector3DWithError __add__  ( const Vector3DWithError& right ) const ;
      Vector3DWithError __add__  ( const Vector3D&          right ) const ;
      Point3DWithError  __add__  ( const Point3DWithError&  right ) const ;
      //
      Vector3DWithError __sub__  ( const Vector3DWithError& right ) const ;
      Vector3DWithError __sub__  ( const Vector3D&          right ) const ;
      // ======================================================================
      Vector3DWithError __radd__ ( const Vector3D&          right ) const 
      { return __add__ ( right ) ; }
      Vector3DWithError __rsub__ ( const Vector3D&          right ) const ;
      // ======================================================================
      Vector3DWithError& __imul__      ( const double v ) ;
      Vector3DWithError& __itruediv__  ( const double v ) ;
      Vector3DWithError& __idiv__      ( const double v ) { return __itruediv__ ( v ) ; }
      // ======================================================================
      Vector3DWithError  __mul__       ( const double v ) const ;
      Vector3DWithError  __truediv__   ( const double v ) const ;
      Vector3DWithError  __div__       ( const double v ) const { return __truediv__ ( v ) ; }
      Vector3DWithError  __rmul__      ( const double v ) const { return __mul__     ( v ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// printout 
      std::ostream& fillStream ( std::ostream& s ) const ;          // printout 
      /// conversion to the string 
      std::string   toString   () const ;           // conversion to the string 
      // ======================================================================
    private:
      // ======================================================================
      /// the covariance 
      Covariance m_cov2     ;                                 // the covariance 
      // ======================================================================
    } ;
    // ========================================================================
    inline Vector3DWithError operator+ 
    ( const Vector3DWithError& a , 
      const Vector3DWithError& b ) { return a.__add__  ( b ) ; }
    inline Vector3DWithError operator+ 
    ( const Vector3DWithError& a , 
      const Ostap::XYZVector&  b ) { return a.__add__  ( b ) ; }
    inline Vector3DWithError operator-
    ( const Vector3DWithError& a , 
      const Vector3DWithError& b ) { return a.__sub__  ( b ) ; }
    inline Vector3DWithError operator-
    ( const Vector3DWithError& a , 
      const Ostap::XYZVector&  b ) { return a.__sub__  ( b ) ; }
    inline Vector3DWithError operator-
    ( const Ostap::XYZVector&  b , 
      const Vector3DWithError& a ) { return a.__rsub__ ( b ) ; }   
    // ========================================================================
    Vector3DWithError operator- 
    ( const Ostap::XYZPoint&  b ,
      const Point3DWithError& a ) ;
    Vector3DWithError operator- 
    ( const Point3DWithError& a , 
      const Point3DWithError& b ) ;
    Vector3DWithError operator- 
    ( const Point3DWithError& a , 
      const Ostap::XYZPoint&  b ) ;
    // ========================================================================
    inline 
    Vector3DWithError operator* 
    ( const Vector3DWithError& a , 
      const double             b ) { return a.__mul__ ( b ) ; }
    inline 
    Vector3DWithError operator/ 
    ( const Vector3DWithError& a , 
      const double             b ) { return a.__truediv__ ( b ) ; }
    inline 
    Vector3DWithError operator* 
    ( const double             b , 
      const Vector3DWithError& a ) { return a.__mul__ ( b ) ; }
    // ========================================================================
    inline double chi2 
    ( const Vector3DWithError& a , 
      const Vector3DWithError& b ) { return a.chi2 ( b ) ; }  
    inline double chi2 
    ( const Vector3DWithError& a , 
      const Ostap::XYZVector&  b ) { return a.chi2 ( b ) ; }  
    inline double chi2 
    ( const Ostap::XYZVector&  b ,
      const Vector3DWithError& a ) { return a.chi2 ( b ) ; }  
    inline double chi2 
    ( const Vector3DWithError&           a , 
      const Vector3DWithError::VectorE&  b ) { return a.chi2 ( b ) ; }  
    inline double chi2 
    ( const Vector3DWithError::VectorE&  b ,
      const Vector3DWithError&           a ) { return a.chi2 ( b ) ; }  
    inline double chi2 
    ( const Vector3DWithError&           a , 
      const Vector3DWithError::Vector&   b ) { return a.chi2 ( b ) ; }  
    inline double chi2 
    ( const Vector3DWithError::Vector&   b ,
      const Vector3DWithError&           a ) { return a.chi2 ( b ) ; }  
    // ========================================================================
    inline Vector3DWithError mean 
    ( const Vector3DWithError& v1 , 
      const Vector3DWithError& v2 ) { return v1.mean ( v2 ) ; }
    // ========================================================================
  } //                                             end of namespace Ostap::Math 
  // ==========================================================================
} //                                                    end of namespace  Ostap
// ============================================================================
// The END 
// ============================================================================
#endif // LHCBMATH_POINTWITHERROR_H
// ============================================================================
