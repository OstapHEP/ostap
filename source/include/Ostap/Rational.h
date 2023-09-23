// ============================================================================
#ifndef OSTAP_RATIONAL_H 
#define OSTAP_RATIONAL_H 1
// ============================================================================
// Include files
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Math.h"
#include "Ostap/Interpolation.h"
#include "Ostap/Interpolants.h"
#include "Ostap/Parameters.h"
#include "Ostap/Workspace.h"
#include "Ostap/Bernstein1D.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    /** @class Rational Ostap/Rational.h
     *  A simple pole-free rational function at interval \f$ x_{min} \le x \le x_{max}\f$
     *  \f[ F(x) = \frac{p(x)}{q(x)} \f]
     *  Actually internally it uses 
     *  the Floater-Hormann rational barycentric interpolant 
     *  and parameters are the function valeus at Chebyshev's nodes  
     *
     *  @see Ostap::Math::FloaterHormann
     *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
     *  @date   2023-09-14
     */
    class Rational : public Parameters 
    {
    public: 
      // ======================================================================
      /** constructor  
       *  @param n degree of numerator 
       *  @param d degree of denumerator is \f$ \max ( n - d , 0 ) \f$
       *  @param xmin  low-edge  of interval
       *  @param xmax  high-edge of interval
       */
      Rational
      ( const unsigned short n = 3 , 
        const unsigned short d = 1 , 
        const double   xmin    = 0 , 
        const double   xmax    = 1 ) ;  
      // =====================================================================
      /** constructor  
       *  @param pars  vector of parameters 
       *  @param d     degree defect,  
       *  @param xmin  low-edge  of interval
       *  @param xmax  high-edge of interval
       */
      Rational
      ( const std::vector<double>& pars      , 
        const unsigned short       d         , 
        const double               xmin  = 0 , 
        const double               xmax  = 1 ) ;  
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the rational function 
      double  evaluate ( const double x ) const ;
      /// evaluate the rational function 
      inline  double operator () ( const double x ) const 
      { return evaluate ( x ) ; }
      // ======================================================================
    public: // getters 
      // ======================================================================
      /// get xmin 
      double   xmin    () const { return m_abscissas.xmin () ; }
      /// get xmax
      double   xmax    () const { return m_abscissas.xmax () ; }
      /// get the n
      unsigned short n () const { return npars      ()  ; }
      /// get the d 
      unsigned short d () const { return m_weights.d () ; }
      // ======================================================================
    public: // integrals 
      // ======================================================================
      // get the integral between xmin and xmax 
      double integral  () const ;
      // ======================================================================
      // get the integral
      double integral 
      (  const double xlow  , 
         const double xhigh ) const ;
      // ======================================================================
    public: // some useful operations with Rational function 
      // ======================================================================
      /// scale rational function by some value 
      Rational& scale ( const double value ) ;
      /// add a constant to rational function 
      Rational& add   ( const double value ) ;
      /// scale as operator 
      Rational& operator *= ( const double value ) { return scale (       value ) ; }
      /// scale as operator 
      Rational& operator /= ( const double value ) { return scale ( 1.0 / value ) ; }
      /// add as operator 
      Rational& operator += ( const double value ) { return add   (       value ) ; }
      /// add as operator 
      Rational& operator -= ( const double value ) { return add   (     - value ) ; }
      /// negate it 
      Rational  operator -  () const ;  
      // ======================================================================
    public: // sme python-related stuff 
      // ======================================================================
      Rational  __add__     ( const double value ) const ;
      Rational  __sub__     ( const double value ) const ;
      Rational  __mul__     ( const double value ) const ;
      Rational  __div__     ( const double value ) const ;
      Rational  __truediv__ ( const double value ) const ;
      Rational  __radd__    ( const double value ) const ;
      Rational  __rsub__    ( const double value ) const ;
      Rational  __rmul__    ( const double value ) const ;
      Rational  __neg__     ()                     const ;
      Rational& __iadd__    ( const double value ) { return (*this) += value ; }        
      Rational& __imul__    ( const double value ) { return (*this) *= value ; }        
      Rational& __isub__    ( const double value ) { return (*this) /= value ; }        
      // ======================================================================
    public:
      // ======================================================================
      /// unique tag 
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// abscissas 
      Ostap::Math::Interpolation::Abscissas m_abscissas {} ; // abscissas 
      /// weights 
      Ostap::Math::FloaterHormann::Weights  m_weights   {} ; // weights 
      /// integration workspace 
      Ostap::Math::WorkSpace                m_workspace {} ; // workspace 
      // ======================================================================
    } ;  
    // ========================================================================
    /// scale Rational function 
    inline Rational operator* ( const Rational& a , const double    b ) 
    { return Rational ( a ) *= b ; }
    /// scale Rational function 
    inline Rational operator* ( const double    b , const Rational& a ) 
    { return Rational ( a ) *= b ; }
    /// scale Rational function 
    inline Rational operator/ ( const Rational& a , const double    b ) 
    { return Rational ( a ) /= b ; }
    /// add Rational function 
    inline Rational operator+ ( const Rational& a , const double    b ) 
    { return Rational ( a ) += b ; }
    /// add Rational function 
    inline Rational operator+ ( const double    b , const Rational& a ) 
    { return Rational ( a ) += b ; }
    /// subtract Rational function 
    inline Rational operator- ( const Rational& a , const double    b ) 
    { return Rational ( a ) -= b ; }
    /// subtract Rational function 
    inline Rational operator- ( const double    b , const Rational& a ) 
    { return  ( -a ) += b ; }
    // ========================================================================
    /** @class RationalBernstein
     *  Rational fnuction as ratio of Bernstein polyhomial and 
     *  positive Bernstein polynomial
     *  \f[ R ( x ) = \frac{B(x)}{P(x) \frac{1}{ x_{max} - x_{min} } \f]
     *  @see Ostap::Math::Bernstein
     *  @see Ostap::Math::Positive 
     *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
     *  @date   2023-09-14
     */
    class RationalBernstein 
    {
    public:
      // ======================================================================
      RationalBernstein
      ( const unsigned short       p      = 3 ,
        const unsigned short       q      = 0 , 
        const double               xmin   = 0 , 
        const double               xmax   = 1 ) ;
      // ======================================================================
      RationalBernstein
      ( const std::vector<double>& p ,
        const std::vector<double>& q ,     
        const double               xmin   = 0 , 
        const double               xmax   = 1 ) ;
      // ======================================================================
      RationalBernstein
      ( const std::vector<double>& a          ,
        const unsigned short       p          ,      
        const double               xmin   = 0 , 
        const double               xmax   = 1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the rational function
      double evaluate ( const double x ) const ;
      /// evaluate the rational function
      inline double operator()  ( const double x ) const { return evaluate ( x ) ; }      
      // ======================================================================
    public:
      // ======================================================================
      double         xmin  () const { return m_p.xmin   () ; }
      double         xmax  () const { return m_p.xmax   () ; }
      /// number of parameters 
      unsigned short npars () const { return m_p.npars() + m_q.npars() ; }
      /// get parameter 
      inline double  par   ( const unsigned short index ) const 
      {
        const unsigned short np = m_p.npars() ;
        return index < np ? m_p.par ( index ) : m_q.par ( index - np ) ; 
      }
      /// set parameter 
      inline double  setPar
      ( const unsigned short index , 
        const double         value )
      {
        const unsigned short np = m_p.npars() ;
        return index < np ? m_p.setPar ( index , value ) : m_q.setPar ( index - np , value ) ; 
      }
      /// all parameters (by value!!!)
      std::vector<double> pars () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// degree of numerator 
      unsigned short   pdegree () const { return m_p.degree () ; }
      /// degree of denomerator 
      unsigned short   qdegree () const { return m_q.degree () ; }
      /// #npars of numerator 
      unsigned short   pnpars  () const { return m_p.npars  () ; }
      /// #npars of denomerator 
      unsigned short   qnpars  () const { return m_q.npars  () ; }
      /// get the numerator 
      const Ostap::Math::Bernstein& numerator   () const { return m_p ; }
      /// get the denominator
      const Ostap::Math::Positive&  denominator () const { return m_q ; }      
      // ======================================================================
      /// numerator parameters 
      const std::vector<double>&    ppars       () const { return m_p.pars() ; }
      // ======================================================================
      // denominator parameters (by value) 
      std::vector<double>           qpars       () const { return m_q.pars() ; }
      // ======================================================================
    public: // integrals 
      // ======================================================================
      // get the integral betweek xmin and xmax 
      double integral  () const ;
      // ======================================================================
      // get the integral
      double integral 
      (  const double xlow  , 
         const double xhigh ) const ;
      // ======================================================================
    public: // some useful operations with Rational function 
      // ======================================================================
      /// scale rational function by some value 
      RationalBernstein& scale ( const double value ) ;
      /// add a constant to rational function 
      RationalBernstein& add   ( const double value ) ;
      /// scale as operator 
      RationalBernstein& operator *= ( const double value ) { return scale (       value ) ; }
      /// scale as operator 
      RationalBernstein& operator /= ( const double value ) { return scale ( 1.0 / value ) ; }
      /// add as operator 
      RationalBernstein& operator += ( const double value ) { return add   (       value ) ; }
      /// add as operator 
      RationalBernstein& operator -= ( const double value ) { return add   (     - value ) ; }
      /// negate it 
      RationalBernstein  operator -  () const ;  
      /// multiply with Bernstein  polynomial
      RationalBernstein& operator *=  ( const Bernstein& right ) ;
      /// add      with with Bernstein polynomial
      RationalBernstein& operator +=  ( const Bernstein& right ) ;
      /// subtarct      with with Bernstein polynomial
      RationalBernstein& operator -=  ( const Bernstein& right ) ;
      // ======================================================================
    public: // sme python-related stuff 
      // ======================================================================
      RationalBernstein  __add__     ( const double value ) const ;
      RationalBernstein  __sub__     ( const double value ) const ;
      RationalBernstein  __mul__     ( const double value ) const ;
      RationalBernstein  __div__     ( const double value ) const ;
      RationalBernstein  __truediv__ ( const double value ) const ;
      RationalBernstein  __radd__    ( const double value ) const ;
      RationalBernstein  __rsub__    ( const double value ) const ;
      RationalBernstein  __rmul__    ( const double value ) const ;
      RationalBernstein  __neg__     ()                     const ;
      RationalBernstein& __iadd__    ( const double value ) { return (*this) += value ; }        
      RationalBernstein& __imul__    ( const double value ) { return (*this) *= value ; }        
      RationalBernstein& __isub__    ( const double value ) { return (*this) /= value ; }        
      // ======================================================================
      RationalBernstein  __add__     ( const Bernstein& right ) const ;
      RationalBernstein  __sub__     ( const Bernstein& right ) const ;
      RationalBernstein  __mul__     ( const Bernstein& right ) const ;
      RationalBernstein  __radd__    ( const Bernstein& right ) const ;
      RationalBernstein  __rsub__    ( const Bernstein& right ) const ;
      RationalBernstein  __rmul__    ( const Bernstein& right ) const ;
      RationalBernstein& __imul__    ( const Bernstein& right ) { return (*this) *= right ; }        
      RationalBernstein& __isub__    ( const Bernstein& right ) { return (*this) -= right ; }        
      // ====================================================================== 
    public:
      // ======================================================================
      /// unique tag 
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// numerator 
      Ostap::Math::Bernstein  m_p         {} ; // numerator 
      /// denominator 
      Ostap::Math::Positive   m_q         {} ; // denominator
      /// integration workspace 
      Ostap::Math::WorkSpace  m_workspace {} ; // workspace 
      // ======================================================================
    } ;
    // ========================================================================
    /// scale Rational function 
    inline RationalBernstein operator* ( const RationalBernstein& a , const double    b ) 
    { return RationalBernstein ( a ) *= b ; }
    /// scale Rational function 
    inline RationalBernstein operator* ( const double    b , const RationalBernstein& a ) 
    { return RationalBernstein ( a ) *= b ; }
    /// scale Rational function 
    inline RationalBernstein operator/ ( const RationalBernstein& a , const double    b ) 
    { return RationalBernstein ( a ) /= b ; }
    /// add Rational function 
    inline RationalBernstein operator+ ( const RationalBernstein& a , const double    b ) 
    { return RationalBernstein ( a ) += b ; }
    /// add Rational function 
    inline RationalBernstein operator+ ( const double    b , const RationalBernstein& a ) 
    { return RationalBernstein  ( a ) += b ; }
    /// subtract Rational function 
    inline RationalBernstein operator- ( const RationalBernstein& a , const double    b ) 
    { return RationalBernstein ( a ) -= b ; }
    /// subtract Rational function 
    inline RationalBernstein operator- ( const double    b , const RationalBernstein& a ) 
    { return  ( -a ) += b ; }
    /// multiply with Bernstein polynomial 
    inline RationalBernstein operator* ( const RationalBernstein& a , const Bernstein& b  ) 
    { return RationalBernstein ( a ) *= b ; }
    /// multiply with Bernstein polynomial 
    inline RationalBernstein operator* ( const Bernstein& b , const RationalBernstein& a  ) 
    { return RationalBernstein ( a ) *= b ; }
    /// add with Bernstein polynomial 
    inline RationalBernstein operator+ ( const RationalBernstein& a , const Bernstein& b  ) 
    { return RationalBernstein ( a ) += b ; }
    /// add with Bernstein polynomial 
    inline RationalBernstein operator+ ( const Bernstein& b , const RationalBernstein& a  ) 
    { return RationalBernstein ( a ) += b ; }
    /// subtract with Bernstein polynomial 
    inline RationalBernstein operator- ( const RationalBernstein& a , const Bernstein& b  ) 
    { return RationalBernstein ( a ) -= b ; }
    /// subtract with Bernstein polynomial 
    inline RationalBernstein operator- ( const Bernstein& b , const RationalBernstein& a  ) 
    { return ( -a ) += b ; }
    // ========================================================================
    /** @class RationalPositive 
     *  Rational fnuction as ratio of two positive Bernstein polyhomials 
     *  \f[ R ( x ) = \frac{B(x)}{P(x) \f]
     *  @see Ostap::Math::Positive 
     *  @author Vanya BELYAEV Ivan.Belyaev@cern.ch
     *  @date   2023-09-14
     */
    class RationalPositive 
    {
    public:
      // ======================================================================
      RationalPositive 
      ( const unsigned short       p      = 3 ,
        const unsigned short       q      = 0 , 
        const double               xmin   = 0 , 
        const double               xmax   = 1 ) ;
      // ======================================================================
      RationalPositive
      ( const std::vector<double>& p ,
        const std::vector<double>& q ,     
        const double               xmin   = 0 , 
        const double               xmax   = 1 ) ;
      // ======================================================================
      RationalPositive 
      ( const std::vector<double>& a          ,
        const unsigned short       p          ,      
        const double               xmin   = 0 , 
        const double               xmax   = 1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the rational function
      double evaluate ( const double x ) const ;
      /// evaluate the rational function
      inline double operator()  ( const double x ) const { return evaluate ( x ) ; }      
      // ======================================================================
    public:
      // ======================================================================
      double         xmin  () const { return m_p.xmin   () ; }
      double         xmax  () const { return m_p.xmax   () ; }
      /// number of parameters 
      unsigned short npars () const { return m_p.npars() + m_q.npars() ; }
      /// get parameter 
      inline double  par   ( const unsigned short index ) const 
      {
        const unsigned short np = m_p.npars() ;
        return index < np ? m_p.par ( index ) : m_q.par ( index - np ) ; 
      }
      /// set parameter 
      inline double  setPar
      ( const unsigned short index , 
        const double         value )
      {
        const unsigned short np = m_p.npars() ;
        return index < np ? m_p.setPar ( index , value ) : m_q.setPar ( index - np , value ) ; 
      }
      /// all parameters (by value!!!)
      std::vector<double> pars () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// degree of numerator 
      unsigned short   pdegree () const { return m_p.degree () ; }
      /// degree of denomerator 
      unsigned short   qdegree () const { return m_q.degree () ; }
      /// #npars of numerator 
      unsigned short   pnpars  () const { return m_p.npars  () ; }
      /// #npars of denomerator 
      unsigned short   qnpars  () const { return m_q.npars  () ; }
      /// get the numerator 
      const Ostap::Math::Positive&  numerator   () const { return m_p ; }
      /// get the denominator
      const Ostap::Math::Positive&  denominator () const { return m_q ; }      
      // ======================================================================
      /// numerator parameters (by value)
      std::vector<double>           ppars       () const { return m_p.pars() ; }
      // ======================================================================
      // denominator parameters (by value) 
      std::vector<double>           qpars       () const { return m_q.pars() ; }
      // ======================================================================
    public: // integrals 
      // ======================================================================
      // get the integral betweek xmin and xmax 
      double integral  () const ;
      // ======================================================================
      // get the integral
      double integral 
      (  const double xlow  , 
         const double xhigh ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// unique tag 
      std::size_t tag () const ;
      // ======================================================================
    private:
      // ======================================================================
      /// numerator 
      Ostap::Math::Positive   m_p         {} ; // numerator 
      /// denominator 
      Ostap::Math::Positive   m_q         {} ; // denominator
      /// integration workspace 
      Ostap::Math::WorkSpace  m_workspace {} ; // workspace 
      // ======================================================================
    } ;
    // ========================================================================
  } //                                         The end of namespace Ostap::Math 
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_RATIONAL_H
// ============================================================================
