// $Id$
// ============================================================================
#ifndef LHCBMATH_LORENTZVECTORWITHERROR_H 
#define LHCBMATH_LORENTZVECTORWITHERROR_H 1
// ============================================================================
// Include files 
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/SymmetricMatrixTypes.h"
#include "Ostap/Vector4DTypes.h"
#include "Ostap/GenericVectorTypes.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/SVectorWithError.h"
#include "Ostap/ValueWithError.h"
// ============================================================================
/** @file 
 *  Collection of useful objects with associated "covarinaces".
 *  The concept has been stollen from Wouter Hulsbergen's lines 
 *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
 *  @date 20090603
 */
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    /** @class LorentzVectorWithError LHCbMath/LorentzVectorWithError.h
     *  
     *
     *  @author Vanya BELYAEV Ivan.Belyaev@nikhef.nl
     *  @date   2009-06-12
     */
    class LorentzVectorWithError : public Ostap::LorentzVector
    {
    public:
      // ======================================================================
      /// the actual type of the Vector 
      typedef Ostap::LorentzVector                                 Vector4D   ;
      typedef Ostap::LorentzVector                                 Value      ;
      /// the actual type of covariance 
      typedef Ostap::SymMatrix4x4                                  Covariance ;
      // ======================================================================
    public:
      // ======================================================================
      /// the actual type of generic 4-vector 
      typedef Ostap::Vector4                                       Vector     ;
      /// the actual type of vector with errors  
      typedef Ostap::Math::SVectorWithError<4,double>              VectorE    ;
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from lorentz vector and covariance
      LorentzVectorWithError 
      ( const Vector4D&   value  = Vector4D   () , 
        const Covariance& cov2   = Covariance () ) ;  
      // ======================================================================
      /// constructor from lorentz vector and covariance
      LorentzVectorWithError 
      ( const Covariance& cov2                   , 
        const Vector4D&   value  = Vector4D   () ) ;
      /// constructor from generic vector and covariance
      LorentzVectorWithError 
      ( const Vector&     value                  , 
        const Covariance& cov2   = Covariance () ) ;  
      /// constructor from generic vector and covariance
      LorentzVectorWithError 
      ( const VectorE&     value                 ) ;
      // ======================================================================
    public: // trivial accessors 
      // ======================================================================
      const  Vector4D&   vector4d   ( ) const { return *this       ; }
      const  Vector4D&   vector4D   ( ) const { return vector4d () ; }
      const  Covariance& covariance ( ) const { return cov2  ()    ; }      
      inline Vector4D&   vector4d   ( )       { return *this       ; }
      // ======================================================================      
      const  Vector4D&   value      ( ) const { return vector4d () ; }
      const  Covariance& cov2       ( ) const { return m_cov2      ; }
      // ======================================================================      
    public: // setters 
      // ======================================================================
      void setVector4d   ( const Vector4D&   v ) { vector4d() =  v   ; }
      void setVector4D   ( const Vector4D&   v ) { setVector4d ( v ) ; }
      void setVector     ( const Vector4D&   v ) { setVector4d ( v ) ; }
      // ======================================================================
      void setValue      ( const Vector4D&   v ) { setVector4d ( v ) ; }
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
      LorentzVectorWithError& operator+= ( const LorentzVectorWithError& right ) ;
      LorentzVectorWithError& operator+= ( const Vector4D&               right ) ;
      LorentzVectorWithError& operator+= ( const VectorE&                right ) ;
      LorentzVectorWithError& operator+= ( const Vector&                 right ) ;
      LorentzVectorWithError& operator-= ( const LorentzVectorWithError& right ) ;
      LorentzVectorWithError& operator-= ( const Vector4D&               right ) ;
      LorentzVectorWithError& operator-= ( const VectorE&                right ) ;
      LorentzVectorWithError& operator-= ( const Vector&                 right ) ;
      // ======================================================================
    public: // scaling
      // ======================================================================
      /// *= 
      LorentzVectorWithError& operator*= ( const double v ) ;            // *= 
      /// /= 
      LorentzVectorWithError& operator/= ( const double v ) ;            // /= 
      // ====================================================================== 
    public:
      // ======================================================================
      /// get generic vector 
      void    asVector  ( VectorE& data ) const ;
      /// get generic vector 
      void    asVector  ( Vector&  data ) const ;
      /// convert to generic vector with errors:
      VectorE asVector  () const ;
      /// convert to generic vector
      Vector  asVector4 () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate chi2 distance 
      double chi2 ( const LorentzVectorWithError& right )  const ;
      /// evaluate chi2 distance 
      double chi2 ( const Vector4D&               right )  const ;      
      /// evaluate chi2 distance 
      double chi2 ( const VectorE&                right )  const ;      
      /// evaluate chi2 distance 
      double chi2 ( const Vector&                 right )  const ;      
      // ======================================================================
    public: // more functions 
      // ======================================================================
      LorentzVectorWithError  mean ( const LorentzVectorWithError& right ) const ;
      LorentzVectorWithError  mean ( const VectorE&                right ) const ;
      // ======================================================================
    public: // useful accessors
      // ======================================================================
      /// get the mass with error 
      inline Ostap::Math::ValueWithError invariantMass           () const ;
      /// get the momentum with error 
      inline Ostap::Math::ValueWithError scalarMomentum          () const ; 
      /// get the transverse momentum with error
      inline Ostap::Math::ValueWithError transverseMomentum      () const ;      
      /// get the transverse mass     with error
      inline Ostap::Math::ValueWithError transverseMass          () const ;      
      /// get the transverse energy with error
      inline Ostap::Math::ValueWithError transverseEnergy        () const ;      
      /// get the kinetic energy with error
      inline Ostap::Math::ValueWithError kineticEnergy           () const ;
      /// get the transverse kinetic energy with error
      inline Ostap::Math::ValueWithError transverseKineticEnergy () const ;
      /// get the rapidity with error
      inline Ostap::Math::ValueWithError rapidity                () const ;
      /// get the pseudorapidity with error
      inline Ostap::Math::ValueWithError pseudorapidity          () const ;
      /// get the phi-angle with error
      inline Ostap::Math::ValueWithError phi                     () const ;
      /// get the theta-angle with error
      inline Ostap::Math::ValueWithError theta                   () const ;
      // ======================================================================
    public: // useful accessors
      // ======================================================================
      /// get the mass with error 
      inline Ostap::Math::ValueWithError mass () const 
      { return invariantMass           () ; }
      /// get the mass with error 
      inline Ostap::Math::ValueWithError m    () const 
      { return invariantMass           () ; }
      /// get the scalar momentum with error 
      inline Ostap::Math::ValueWithError p    () const
      { return scalarMomentum          () ; }
      /// get the transverse momentum with error 
      inline Ostap::Math::ValueWithError pt   () const
      { return transverseMomentum      () ; }
      /// get the transverse momentum with error 
      inline Ostap::Math::ValueWithError pT   () const
      { return transverseMomentum      () ; }
      /// get the transverse energy with error 
      inline Ostap::Math::ValueWithError et   () const
      { return transverseEnergy        () ; }
      /// get the transverse energy  with error 
      inline Ostap::Math::ValueWithError eT   () const
      { return transverseEnergy        () ; }
      /// get the transverse mass with error 
      inline Ostap::Math::ValueWithError mt   () const 
      { return transverseMass          () ; }
      /// get the transverse mass  with error 
      inline Ostap::Math::ValueWithError mT   () const
      { return transverseMass          () ; }
      /// get the kinetic energy with error
      inline Ostap::Math::ValueWithError eK   () const 
      { return kineticEnergy           () ; }
      /// get the transverse kinetic energy with error
      inline Ostap::Math::ValueWithError eTk   () const 
      { return transverseKineticEnergy () ; }
      // ======================================================================
    public:
      // ======================================================================
      /** calculate dispersion for some function <code>F(p)</code> using 
       *  the vector of derivatives <code>dFdP</code>
       *  @code
       *  LorentzVectorWithError v = ... ;
       *  // squared uncertainty in pt2 
       *  Ostap::Vector4 dFdp 
       *  dFdP[0] = 2 * v.Px() ; // d(pt^2)/d(px) 
       *  dFdP[1] = 2 * v.Py() ; // d(pt^2)/d(py) 
       *  dFdP[2] = 2 * v.Pz() ; // d(pt^2)/d(pz) 
       *  dFdP[3] = 2 * v.E () ; // d(pt^2)/d(e) 
       *  const  double sigma2 = v.dispersion ( dFdP ) ;
       *  // get a final result: 
       *  const ValueWithError result ( v.Perp2() , sigma2 ) ;       
       *  @endcode 
       */
      double dispersion ( const Vector& dFdp ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /** evaluate chi2-mass
       *  @see Ostap::Math::chi2mass
       */
      double chi2mass ( const double m0 ) const ;         // evaluate chi2-mass
      // ======================================================================
    public: // helper operators for Python
      // ======================================================================
      LorentzVectorWithError __add__   ( const LorentzVectorWithError& right ) const ;
      LorentzVectorWithError __sub__   ( const LorentzVectorWithError& right ) const ;
      //
      LorentzVectorWithError __add__   ( const Vector4D&               right ) const ;
      LorentzVectorWithError __add__   ( const VectorE&                right ) const ;
      LorentzVectorWithError __add__   ( const Vector&                 right ) const ;
      LorentzVectorWithError __sub__   ( const Vector4D&               right ) const ;
      LorentzVectorWithError __sub__   ( const VectorE&                right ) const ;
      LorentzVectorWithError __sub__   ( const Vector&                 right ) const ;
      // ======================================================================
      LorentzVectorWithError __radd__  ( const Vector4D&               right ) const 
      { return __add__ ( right ) ; }  
      LorentzVectorWithError __rsub__  ( const Vector4D&               right ) const ;
      // ======================================================================
      LorentzVectorWithError& __imul__ ( const double v ) ;
      LorentzVectorWithError& __idiv__ ( const double v ) ;
      LorentzVectorWithError  __mul__  ( const double v ) const ;
      LorentzVectorWithError  __div__  ( const double v ) const ;
      LorentzVectorWithError  __rmul__ ( const double v ) const 
      { return __mul__ ( v ) ; }
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
    };
    // ========================================================================
    inline LorentzVectorWithError
    operator+ 
    ( const LorentzVectorWithError& a ,
      const LorentzVectorWithError& b ){ return a.__add__ ( b ) ; }
    inline LorentzVectorWithError
    operator- 
    ( const LorentzVectorWithError& a ,
      const LorentzVectorWithError& b ){ return a.__sub__ ( b ) ; }
    inline LorentzVectorWithError
    operator+ 
    ( const LorentzVectorWithError& a ,
      const Ostap::LorentzVector&   b ){ return a.__add__ ( b ) ; }
    inline LorentzVectorWithError
    operator- 
    ( const LorentzVectorWithError& a ,
      const Ostap::LorentzVector&   b ){ return a.__sub__ ( b ) ; }
    inline LorentzVectorWithError
    operator+ 
    ( const Ostap::LorentzVector&   b ,
      const LorentzVectorWithError& a ) { return a + b ; }
    inline LorentzVectorWithError
    operator- 
    ( const Ostap::LorentzVector&   b , 
      const LorentzVectorWithError& a ) { return a.__rsub__ ( b ) ; }
    // ========================================================================
    inline LorentzVectorWithError
    operator* 
    ( const LorentzVectorWithError& a ,
      const double                  b ) { return a.__mul__ ( b ) ; }
    inline LorentzVectorWithError
    operator/ 
    ( const LorentzVectorWithError& a ,
      const double                  b ) { return a.__div__ ( b ) ; }
    inline LorentzVectorWithError
    operator* 
    ( const double                  b , 
      const LorentzVectorWithError& a ) { return a.__mul__ ( b ) ; }
    // ========================================================================
    inline double chi2 
    ( const LorentzVectorWithError& a , 
      const LorentzVectorWithError& b ) { return a.chi2  ( b ) ; }
    inline double chi2 
    ( const LorentzVectorWithError& a , 
      const Ostap::LorentzVector&   b ) { return a.chi2  ( b ) ; }
    inline double chi2 
    ( const Ostap::LorentzVector&   b ,
      const LorentzVectorWithError& a ) { return a.chi2  ( b ) ; }
    inline double chi2 
    ( const LorentzVectorWithError&          a , 
      const LorentzVectorWithError::VectorE& b ) { return a.chi2  ( b ) ; }
    inline double chi2 
    ( const LorentzVectorWithError::VectorE& b , 
      const LorentzVectorWithError&          a ) { return a.chi2  ( b ) ; }
    inline double chi2 
    ( const LorentzVectorWithError&          a , 
      const LorentzVectorWithError::Vector&  b ) { return a.chi2  ( b ) ; }
    inline double chi2 
    ( const LorentzVectorWithError::Vector&  b , 
      const LorentzVectorWithError&          a ) { return a.chi2  ( b ) ; }    
    // ========================================================================
    namespace Kinematics 
    {
      // ======================================================================
      // set of helper functions 
      // ======================================================================
      /** calculate the mass with uncertainty 
       *  @param mom 4-momentum 
       *  @param cov covariance 
       *  @return mass with uncertainty 
       */
      Ostap::Math::ValueWithError mass 
      ( const Ostap::LorentzVector& mom , 
        const Ostap::SymMatrix4x4&  cov ) ;
      // =====================================================================
      /** calculate the scalar momentum (p) with uncertainty 
       *  @param mom 4-momentum 
       *  @param cov covariance 
       *  @return p with uncertainty 
       */
      Ostap::Math::ValueWithError momentum 
      ( const Ostap::LorentzVector& mom , 
        const Ostap::SymMatrix4x4&  cov ) ;
      // =====================================================================
      /** calculate the rapidity (y) with uncertainty 
       *  @param mom 4-momentum 
       *  @param cov covariance 
       *  @return y with uncertainty 
       */
      Ostap::Math::ValueWithError rapidity 
      ( const Ostap::LorentzVector& mom , 
        const Ostap::SymMatrix4x4&  cov ) ;
      // ======================================================================
      /** calculate the pseudorapidity (eta) with uncertainty 
       *  @param mom 4-momentum 
       *  @param cov covariance 
       *  @return eta with uncertainty 
       */
      Ostap::Math::ValueWithError pseudorapidity 
      ( const Ostap::LorentzVector& mom , 
        const Ostap::SymMatrix4x4&  cov ) ;
      // ======================================================================
      /** calculate the phi (asymuthal angle) with uncertainy 
       *  \f$ \tan \phi = \frac{p_y}{p_x} \f$ 
       *  @param mom 4-momentum 
       *  @param cov covariance 
       *  @return phi with uncertainty 
       */
      Ostap::Math::ValueWithError phi
      ( const Ostap::LorentzVector& mom , 
        const Ostap::SymMatrix4x4&  cov ) ;
      // ======================================================================
      /** calculate the theta (polar angle) with uncertainy 
       *  \f$ \tan \theta = \frac{p_T}{p_z} \f$ 
       *  @param mom 4-momentum 
       *  @param cov covariance 
       *  @return theta with uncertainty 
       */
      Ostap::Math::ValueWithError theta
      ( const Ostap::LorentzVector& mom , 
        const Ostap::SymMatrix4x4&  cov ) ;
      // ======================================================================
      /** calculate the transverse momentum (pt) with uncertainty 
       *  @param mom 4-momentum 
       *  @param cov covariance 
       *  @return pt with uncertainty 
       */
      Ostap::Math::ValueWithError transverseMomentum  
      ( const Ostap::LorentzVector& mom , 
        const Ostap::SymMatrix4x4&  cov ) ;
      // ======================================================================
      /** calculate the squared transverse mass  (mT2) with uncertainty 
       *  \f$ m_T^2 = e^2 - p_z^2  = m^2 + p_T^2  \f$ 
       *  @param mom 4-momentum 
       *  @param cov covariance 
       *  @return mT2 with uncertainty 
       */
      Ostap::Math::ValueWithError transverseMass2
      ( const Ostap::LorentzVector& mom , 
        const Ostap::SymMatrix4x4&  cov ) ;
      // ======================================================================
      /** calculate the transverse mass  (mT) with uncertainty 
       *  \f$ m_T^2 = \sqrt{e^2 - p_z^2}  = \sqrt{m^2 + p_T^2}  \f$ 
       *  @param mom 4-momentum 
       *  @param cov covariance 
       *  @return mT with uncertainty 
       */
      Ostap::Math::ValueWithError transverseMass
      ( const Ostap::LorentzVector& mom , 
        const Ostap::SymMatrix4x4&  cov ) ;
      // ======================================================================
      /** calculate the squared transverse energy  (eT2) with uncertainty 
       *  \f$ e_T^2 = \frac{e^2p_t^2}{p^2} \f$ 
       *  @param mom 4-momentum 
       *  @param cov covariance 
       *  @return eT2 with uncertainty 
       */
      Ostap::Math::ValueWithError transverseEnergy2
      ( const Ostap::LorentzVector& mom , 
        const Ostap::SymMatrix4x4&  cov ) ;
      // ======================================================================
      /** calculate the transverse energy  (eT) with uncertainty 
       *  \f$ e_T^2 = \frac{e p_t}{p} \f$ 
       *  @param mom 4-momentum 
       *  @param cov covariance 
       *  @return eT with uncertainty 
       */
      Ostap::Math::ValueWithError transverseEnergy
      ( const Ostap::LorentzVector& mom , 
        const Ostap::SymMatrix4x4&  cov ) ;
      // ======================================================================
      /** calculate the kinetic energy  (eK) with uncertainty 
       *  \f$ e_K = e - m \f$ 
       *  @param mom 4-momentum 
       *  @param cov covariance 
       *  @return eK with uncertainty 
       */
      Ostap::Math::ValueWithError kineticEnergy
      ( const Ostap::LorentzVector& mom , 
        const Ostap::SymMatrix4x4&  cov ) ;
      // ======================================================================
      /** calculate the transverse kinetic energy  (eTK) with uncertainty 
       *  \f$ e_{T,k} = m_T - m \f$ 
       *  @param mom 4-momentum 
       *  @param cov covariance 
       *  @return eTK with uncertainty 
       */
      Ostap::Math::ValueWithError transverseKineticEnergy
      ( const Ostap::LorentzVector& mom , 
        const Ostap::SymMatrix4x4&  cov ) ;
      // ======================================================================
    } //                               end of namespace Ostap::Math::Kinematics 
    // ========================================================================
    // get the mass with error 
    // ========================================================================
    inline Ostap::Math::ValueWithError 
    Ostap::Math::LorentzVectorWithError::invariantMass      () const 
    { return Ostap::Math::Kinematics::mass         ( *this , cov2() ) ; }
    // ========================================================================
    // get the scalar momentum with error 
    // ========================================================================
    inline Ostap::Math::ValueWithError 
    Ostap::Math::LorentzVectorWithError::scalarMomentum     () const 
    { return Ostap::Math::Kinematics::momentum           ( *this , cov2() ) ; }
    // ========================================================================
    // get the transverse momentum with error 
    // ========================================================================
    inline Ostap::Math::ValueWithError 
    Ostap::Math::LorentzVectorWithError::transverseMomentum () const 
    { return Ostap::Math::Kinematics::transverseMomentum ( *this , cov2() ) ; }
    // ========================================================================
    // get the transverse mass with error 
    // ========================================================================
    inline Ostap::Math::ValueWithError 
    Ostap::Math::LorentzVectorWithError::transverseMass     () const 
    { return Ostap::Math::Kinematics::transverseMass     ( *this , cov2() ) ; }
    // ========================================================================
    // get the transverse energy with error 
    // ========================================================================
    inline Ostap::Math::ValueWithError 
    Ostap::Math::LorentzVectorWithError::transverseEnergy   () const 
    { return Ostap::Math::Kinematics::transverseEnergy ( *this , cov2() ) ; }
    // ========================================================================
    // get the kinetic energy with error 
    // ========================================================================
    inline Ostap::Math::ValueWithError 
    Ostap::Math::LorentzVectorWithError::kineticEnergy      () const 
    { return Ostap::Math::Kinematics::kineticEnergy    ( *this , cov2() ) ; }
    // ========================================================================
    // get the transverse kinetic energy with error 
    // ========================================================================
    inline Ostap::Math::ValueWithError 
    Ostap::Math::LorentzVectorWithError::transverseKineticEnergy () const 
    { return Ostap::Math::Kinematics::transverseKineticEnergy ( *this , cov2() ) ; }
    // ========================================================================
    // get the rapidity with error 
    // ========================================================================
    inline Ostap::Math::ValueWithError 
    Ostap::Math::LorentzVectorWithError::rapidity           () const 
    { return Ostap::Math::Kinematics::rapidity         ( *this , cov2() ) ; }
    // ========================================================================
    // get the pseudorapidity with error 
    // ========================================================================
    inline Ostap::Math::ValueWithError 
    Ostap::Math::LorentzVectorWithError::pseudorapidity     () const 
    { return Ostap::Math::Kinematics::pseudorapidity   ( *this , cov2() ) ; }
    // ========================================================================
    // get the phi-angle with error 
    // ========================================================================
    inline Ostap::Math::ValueWithError 
    Ostap::Math::LorentzVectorWithError::phi                () const 
    { return Ostap::Math::Kinematics::phi              ( *this , cov2() ) ; }
    // ========================================================================
    // get the theta-angle with error 
    // ========================================================================
    inline Ostap::Math::ValueWithError 
    Ostap::Math::LorentzVectorWithError::theta              () const 
    { return Ostap::Math::Kinematics::theta            ( *this , cov2() ) ; }
    // ========================================================================
    inline  Ostap::Math::LorentzVectorWithError mean 
    ( const Ostap::Math::LorentzVectorWithError& v1 , 
      const Ostap::Math::LorentzVectorWithError& v2 ) { return v1.mean ( v2 ) ; }
    // ========================================================================
  } //                                             end of namespace Ostap::Math 
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
// The END 
// ============================================================================
#endif // OSTAP_LORENTZVECTORWITHERROR_H
// ============================================================================
