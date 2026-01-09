// ============================================================================
#ifndef OSTAP_MODELS2D_H 
#define OSTAP_MODELS2D_H 1
// ============================================================================
// include files
// ============================================================================
// STD & STL
// ============================================================================
#include <cmath>
#include <vector>
#include <complex>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Models.h"
#include "Ostap/Bernstein2D.h"
#include "Ostap/Integrator.h"
// ============================================================================
/** @file Ostap/Models2D.h
 *  collection of 2D-models
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2010-04-19
 */
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    /** @class PS2DPol
     *  The 2D-function, that represent a cross-product of two phase-space factors,
     *  \f$ \Phi_{k,l}(x)\f$ and \f$ \Phi_{m,n}(y)\f$,  
     *  modulated by the 2D-positive polynomial.
     *
     *  The function is:
     *  \f[ f(x,y) = 
     *      \Phi_{k,l}(x;x_{low}, x_{high})
     *      \Phi_{m,n}(y;y_{low}, y_{high})
     *      P_{N,M}(x,y) \f]
     *  where  
     *  - \f$ \Phi_{k,l}(x;x_{low},x_{high}) \f$ is a phase-space function for x-axis
     *  - \f$ \Phi_{m,n}(y;y_{low},y_{high}) \f$ is a phase-space function for y-axis
     *  - \f$ P_{N,M}(x,y) \f$ is 2D positive Bernstein polynomial 
     *  @see Ostap::Math::PhaseSpaceNL 
     *  @see Ostap::Math::Positive2D
     */
    class  PS2DPol : public Ostap::Math::PolyFactor2D 
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the order
      PS2DPol
      ( const PhaseSpaceNL&   psx  = PhaseSpaceNL () ,
        const PhaseSpaceNL&   psy  = PhaseSpaceNL () ,
        const unsigned short  Nx   = 1               ,
        const unsigned short  Ny   = 1               ) ;
      /// constructor from the order
      PS2DPol
      ( const PhaseSpaceNL&   psx        ,
        const PhaseSpaceNL&   psy        ,
        const unsigned short  Nx         ,
        const unsigned short  Ny         ,
        const double          xmin       ,
        const double          xmax       ,
        const double          ymin       ,
        const double          ymax       ) ;
      // ======================================================================
      /// constructor from the polynomial and phase spaces
      PS2DPol 
      ( const Ostap::Math::Positive2D& pol ,
        const PhaseSpaceNL&            psx ,
        const PhaseSpaceNL&            psy ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      inline double operator ()
      ( const double x ,
	const double y ) const { return evaluate ( x , y ) ; } 
      /// get the value      
      double evaluate
      ( const double x ,
	const double y ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /** get the integral over 2D-region
       *  \f[ \int_{x_low}^{x_high}\int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y\f]
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       *  @param ylow  low  edge in y
       *  @param yhigh high edge in y
       */
      double integral 
      ( const double xlow  ,
        const double xhigh ,
        const double ylow  , 
        const double yhigh ) const ;
      // ======================================================================
      /** integral over x-dimension
       *  \f[ \int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param y     variable
       *  @param xlow  low  edge in y
       *  @param xhigh high edge in y
       */
      double integrateX 
      ( const double y     ,
        const double xlow  , 
        const double xhigh ) const ;
      // ======================================================================
      /** integral over x-dimension
       *  \f[ \int_{x_low}^{x_high} \mathcal{B}(x,y) \mathrm{d}x\f]
       *  @param x     variable
       *  @param ylow  low  edge in x
       *  @param yhigh high edge in x
       */
      double integrateY 
      ( const double x     ,
        const double ylow  ,
        const double yhigh ) const ;
      // ======================================================================
    public:
      // ======================================================================
      const Ostap::Math::PhaseSpaceNL& psX         () const { return m_psx      ; }
      const Ostap::Math::PhaseSpaceNL& psY         () const { return m_psy      ; }
      const Ostap::Math::PhaseSpaceNL& phasespaceX () const { return psX ()     ; }
      const Ostap::Math::PhaseSpaceNL& phasespaceY () const { return psY ()     ; }
      // ====================================== ===============================
    public:
      // ====================================== ===============================
      /// get the tag 
      std::size_t tag () const ;
      // ====================================== ===============================
    private: // helper function to make the calculations
      // ======================================================================
      /// helper function to make calculations
      double calculate
      ( const std::vector<double>& fx , 
	const std::vector<double>& fy ) const ;
      // ======================================================================
    private:
      // ======================================================================
      /// Phase space factors 
      Ostap::Math::PhaseSpaceNL m_psx      ;
      Ostap::Math::PhaseSpaceNL m_psy      ;
      // ====================================================================== 
    private:
      // ======================================================================
      /// workspace for numerical integration
      Ostap::Math::WorkSpace m_workspace    ;
      // ======================================================================
   };
    // ========================================================================
    /** @class PS2DPolSym
     *  The symmetric 2D-function, that represent a cross-product 
     *  \f$ \Phi_{k,l}(x)\f$ and \f$ \Phi_{m,n}(y)\f$,  
     *  modulated by the 2D-positive symmetric polynomial.
     *  It is a "symmetrised" version of class Ostap::Math::PS2DPol
     *
     *  The function is:
     *  \f[ f(x,y) = 
     *      \Phi_{k,l}(x;x_{low}, x_{high})
     *      \Phi_{k,l}(y;y_{low}, y_{high})
     *       P_{N,N}(x,y) \f]
     *  where  
     *  - \f$ \Phi_{k,l}(x;x_{low},x_{high}) \f$ is a phase-space function,
     *        \f$ y_{low}=x_{low}\f$ and \f$y_{high}=x_{high}\f$
     *  - \f$ P_{N,N}(x,y) \f$ is 2D positive symmetric  Bernstein polynomial 
     *
     *  Clearly the function is symmetric: \f$f(x,y) = f(y,x) \f$ 
     *  @see Ostap::Math::PhaseSpaceNL 
     *  @see Ostap::Math::Positive2DSym
     *  @see Ostap::Math::PS2DPol
     */
    class  PS2DPolSym : public Ostap::Math::PolyFactor2DSym
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the order
      PS2DPolSym 
      ( const PhaseSpaceNL&   ps   = PhaseSpaceNL() ,
        const unsigned short  N    =  1             ) ;
      /// constructor from the order
      PS2DPolSym 
      ( const PhaseSpaceNL&   ps         ,
        const unsigned short  N          ,
        const double          xmin       ,
        const double          xmax       ) ;
      /// constructor from the polynomial and phase spaces
      PS2DPolSym 
      ( const Ostap::Math::Positive2DSym& pol ,
        const PhaseSpaceNL&               ps ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      inline double operator ()
      ( const double x ,
	const double y ) const { return evaluate ( x , y ) ; }
      /// get the value
      double evaluate 
      ( const double x ,
	const double y ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /** get the integral over 2D-region
       *  \f[ \int_{x_low}^{x_high}\int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y\f]
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       *  @param ylow  low  edge in y
       *  @param yhigh high edge in y
       */
      double integral 
      ( const double xlow  , 
        const double xhigh ,
        const double ylow  , 
        const double yhigh ) const ;
      // ======================================================================
      /** integral over x-dimension
       *  \f[ \int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param y     variable
       *  @param xlow  low  edge in y
       *  @param xhigh high edge in y
       */
      double integrateX 
      ( const double y     ,
        const double xlow  , 
        const double xhigh ) const ;
      // ======================================================================
      /** integral over x-dimension
       *  \f[ \int_{x_low}^{x_high} \mathcal{B}(x,y) \mathrm{d}x\f]
       *  @param x     variable
       *  @param ylow  low  edge in x
       *  @param yhigh high edge in x
       */
      double integrateY
      ( const double x     ,
        const double ylow  ,
        const double yhigh ) const ;
      // ======================================================================
    public:
      // ======================================================================
      const Ostap::Math::PhaseSpaceNL&  psX         () const { return m_ps       ; }
      const Ostap::Math::PhaseSpaceNL&  psY         () const { return m_ps       ; }
      const Ostap::Math::PhaseSpaceNL&  phasespaceX () const { return psX()      ; }
      const Ostap::Math::PhaseSpaceNL&  phasespaceY () const { return psY()      ; }
      // ====================================== ===============================
    public:
      // ====================================== ===============================
      /// get the tag 
      std::size_t tag () const ;
      // ====================================== ===============================
    private: // helper functions to make the calculations
      // ======================================================================
      /// helper function to make calculations
      double calculate
      ( const std::vector<double>& fx , 
	const std::vector<double>& fy ) const ;
      // ======================================================================
    private:
      // ======================================================================
      /// Phase space factor 
      Ostap::Math::PhaseSpaceNL m_ps        ;  // Phase space factor 
      // ======================================================================
    private:
      // ======================================================================
      /// workspace for numerical integration
      Ostap::Math::WorkSpace m_workspace    ;
      // ======================================================================
    };
    // ========================================================================
    /** @class PS2DPol2
     *  
     *  The 2D-function, that represent non-factorizeable "product" of  
     *  phase-space functions modulated by the 2D-positive polynomial.
     *  The  function is useful to describe e.g. 2D-distributions of 
     *  \f$ m_{23}\f$ vs \f$m_{45}\f$ from 5-body decays. 
     *
     *  The function is:
     *  \f[ f(x,y) = \frac{1}{2}
     *   \left( \Phi_{k,n}(x;x_{low},x_{high}) \Phi_{l,m-1}(y,y_{low},m_{max}-x) 
     *        + \Phi_{l,m}(y;y_{low},y_{high}) \Phi_{k,n-1}(x,x_{low},m_{max}-y) 
     *   \right) P_{N^{x},N^{y}}(x,y) \f]
     *  where 
     *  - \f$ \Phi_{i,j}(z;z_{low},z_{high}\f$ are normalized phase space functions
     *    for mass of \f$i\f$-particles from \f$j\f$-body decays
     *  - \f$ P_{N^{x},N^{y}}(x,y) \f$ is 2D positive Bernstein polynomial
     *  - \f$m_{max}\f$ is a maximal allowed mass for \f$x+y\f$
     *  
     *  @see Ostap::Math::PhaseSpaceNL 
     *  @see Ostap::Math::Positive2D
     */
    class  PS2DPol2 : public Ostap::Math::PolyFactor2D
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the order
      PS2DPol2
      ( const PhaseSpaceNL&   psx  = PhaseSpaceNL () ,
        const PhaseSpaceNL&   psy  = PhaseSpaceNL () ,
        const double          mmax =  -1 ,  
        const unsigned short  Nx   =   1 ,
        const unsigned short  Ny   =   1 ) ;
      /// constructor from the order
      PS2DPol2
      ( const PhaseSpaceNL&   psx        ,
        const PhaseSpaceNL&   psy        ,
        const double          mmax       ,  
        const unsigned short  Nx         ,
        const unsigned short  Ny         ,
        const double          xmin       ,
        const double          xmax       ,
        const double          ymin       ,
        const double          ymax       ) ;
      // ======================================================================
      /// constructor from the polynomial and phase spaces
      PS2DPol2
      ( const Ostap::Math::Positive2D& pol  ,
        const PhaseSpaceNL&            psx  ,
        const PhaseSpaceNL&            psy  ,
        const double                   mmax = -1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      inline double operator ()
      ( const double x ,
	const double y ) const { return evaluate ( x , y ) ; } 
      /// get the value
      double evaluate
      ( const double x ,
	const double y ) const ;
      // ======================================================================
    public :
      // ======================================================================
      /// gett the maximal value for m(x) + m(y) 
      inline double mmax () const { return m_mmax ; }  
      // ====================================== ===============================      
    public:
      // ======================================================================
      /** get the integral over 2D-region
       *  \f[ \int_{x_low}^{x_high}\int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y\f]
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       *  @param ylow  low  edge in y
       *  @param yhigh high edge in y
       */
      double integral 
      ( const double xlow  ,
        const double xhigh ,
        const double ylow  , 
        const double yhigh ) const ;
      // ======================================================================
      /** integral over x-dimension
       *  \f[ \int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param y     variable
       *  @param xlow  low  edge in y
       *  @param xhigh high edge in y
       */
      double integrateX 
      ( const double y     ,
        const double xlow  , 
        const double xhigh ) const ;
      // ======================================================================
      /** integral over x-dimension
       *  \f[ \int_{x_low}^{x_high} \mathcal{B}(x,y) \mathrm{d}x\f]
       *  @param x     variable
       *  @param ylow  low  edge in x
       *  @param yhigh high edge in x
       */
      double integrateY 
      ( const double x     ,
        const double ylow  ,
        const double yhigh ) const ;
      // ======================================================================
    public:
      // ======================================================================
      const Ostap::Math::PhaseSpaceNL& psX         () const { return m_psx      ; }
      const Ostap::Math::PhaseSpaceNL& psY         () const { return m_psy      ; }
      const Ostap::Math::PhaseSpaceNL& phasespaceX () const { return psX ()     ; }
      const Ostap::Math::PhaseSpaceNL& phasespaceY () const { return psY ()     ; }
      // ====================================== ===============================
    public:
      // ====================================== ===============================
      /// get the tag
      std::size_t tag () const ;
      // ====================================== ===============================
    private:
      // ======================================================================
      /// Phase space factors 
      Ostap::Math::PhaseSpaceNL m_psx      ;
      Ostap::Math::PhaseSpaceNL m_psy      ;
      /// Max-mass 
      double                    m_mmax     ;
      // ======================================================================
    private:
      // ======================================================================
      /// Phase space
      mutable Ostap::Math::PhaseSpaceNL m_psx_aux ;
      mutable Ostap::Math::PhaseSpaceNL m_psy_aux ;
      // ======================================================================
      /// workspace
      Ostap::Math::WorkSpace m_workspace   ;
      // ======================================================================
    };
    // ========================================================================
    /** @class PS2DPol2Sym
     *
     *  The symmetric 2D-function, that represent non-factorizeable "product" of  
     *  phase-space functions modulated by the 2D-positive polynomial.
     *  It is a  symmetrised version of class Ostap::Math::PS2DPol2.
     *  The  function is useful to describe e.g. 2D-distributions of 
     *  \f$ m_{23}\f$ vs \f$m_{45}\f$ from 5-body decays. 
     *
     *  The function is:
     *  \f[ f(x,y) = \frac{1}{2}
     *   \left( \Phi_{k,n}(x;x_{low},x_{high}) \Phi_{k,n-1}(y,y_{low},m_{max}-x) 
     *        + \Phi_{k,n}(y;y_{low},y_{high}) \Phi_{k,n-1}(x,x_{low},m_{max}-y) 
     *        \right)
     *        P_{N,N}(x,y) \f]
     *  where 
     *  - \f$ \Phi_{i,j}(x;x_{low},x_{high}\f$ are normalized phase space function,
     *    for mass of \f$i\f$-particles from \f$j\f$-body decays;
     *  - \f$ y_{low}=x_{low}\f$ and \f$y_{high}=x_{high}\f$
     *  - \f$ P_{N,N}(x,y) \f$ is 2D positive symmertic Bernstein polynomial
     *  - \f$m_{max}\f$ is a maximal allowed mass for \f$x+y\f$
     *
     *  Clearly the function is symmetric \f$f(x,y) = f(y,x) \f$ 
     *  @see Ostap::Math::PhaseSpaceNL 
     *  @see Ostap::Math::Positive2DSym
     *  @see Ostap::Math::PS2DPol2
     */
    class  PS2DPol2Sym : public Ostap::Math::PolyFactor2DSym
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the order
      PS2DPol2Sym
      ( const PhaseSpaceNL&   ps   = PhaseSpaceNL() ,
        const double          mmax =  -1            ,  
        const unsigned short  N    =   1            ) ;
      /// constructor from the order
      PS2DPol2Sym
      ( const PhaseSpaceNL&   ps         ,
        const double          mmax       ,  
        const unsigned short  N          ,
        const double          xmin       ,
        const double          xmax       ) ;
      // ======================================================================
      /// constructor from the polynomial and phase spaces
      PS2DPol2Sym 
      ( const Ostap::Math::Positive2DSym& pol       ,
        const PhaseSpaceNL&               ps        ,
        const double                      mmax = -1 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      inline double operator ()
      ( const double x ,
	const double y ) const { return evaluate ( x , y ) ; }
      /// get the value
      double evaluate
      ( const double x ,
	const double y ) const ;
      // ======================================================================
    public :
      // ======================================================================
      /// gett the maximal value for m(x) + m(y) 
      inline double mmax () const { return m_mmax ; }  
      // ====================================== ===============================      
    public:
      // ======================================================================
      /** get the integral over 2D-region
       *  \f[ \int_{x_low}^{x_high}\int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y\f]
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       *  @param ylow  low  edge in y
       *  @param yhigh high edge in y
       */
      double integral 
      ( const double xlow  , 
        const double xhigh ,
        const double ylow  , 
        const double yhigh ) const ;
      // ======================================================================
      /** integral over x-dimension
       *  \f[ \int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param y     variable
       *  @param xlow  low  edge in y
       *  @param xhigh high edge in y
       */
      double integrateX
      ( const double y     ,
        const double xlow  , 
        const double xhigh ) const ;
      // ======================================================================
      /** integral over x-dimension
       *  \f[ \int_{x_low}^{x_high} \mathcal{B}(x,y) \mathrm{d}x\f]
       *  @param x     variable
       *  @param ylow  low  edge in x
       *  @param yhigh high edge in x
       */
      double integrateY 
      ( const double x     ,
        const double ylow  ,
        const double yhigh ) const ;
      // ======================================================================
    public:
      // ======================================================================
      const Ostap::Math::PhaseSpaceNL&  psX         () const { return m_ps       ; }
      const Ostap::Math::PhaseSpaceNL&  psY         () const { return m_ps       ; }
      const Ostap::Math::PhaseSpaceNL&  phasespaceX () const { return psX()      ; }
      const Ostap::Math::PhaseSpaceNL&  phasespaceY () const { return psY()      ; }
      // ====================================== ===============================
    public:
      // ====================================== ===============================
      /// get the tag 
      std::size_t tag () const ;
      // ====================================== ===============================
    private:
      // ======================================================================
      /// Phase space
      Ostap::Math::PhaseSpaceNL m_ps        ;
      // ======================================================================
      /// Max-mass 
      double                    m_mmax           ;
      // ======================================================================
    private:
      // ======================================================================
      /// Phase space
      mutable Ostap::Math::PhaseSpaceNL m_psx_aux ;
      mutable Ostap::Math::PhaseSpaceNL m_psy_aux ;
      // ======================================================================
    private:
      // ======================================================================
      /// workspace for numerical integration
      Ostap::Math::WorkSpace m_workspace    ;
      // ======================================================================
    };
    // ========================================================================
    /** @class PS2DPol3
     *  
     *  The 2D-function, that represent non-factorizeable "product" of  
     *  two modulated phase-space functions.
     *  It can be considered as a simpler alternative for class Ostap::Math::PS2DPol2
     *  The  function is useful to describe e.g. 2D-distributions of 
     *  \f$ m_{23}\f$ vs \f$m_{45}\f$ from 5-body decays. 
     *
     *  The function is:
     *  \f[ f(x,y) = \frac{1}{2}
     *   \left( \Phi^{(N^{x})}_{k,n}(x;x_{low},x_{high}) \Phi_{l,m-1}(y,y_{low},m_{max}-x) 
     *        + \Phi^{(N^{y})}_{l,m}(y;y_{low},y_{high}) \Phi_{k,n-1}(x,x_{low},m_{max}-y) 
     *    \right) \f]
     *  where 
     *  - \f$\Phi_{i,j}(z;z_{low},z_{high}\f$ are normalized phase space functions
     *    for mass of \f$i\f$-particles from \f$j\f$-body decays
     *  - \f$\Phi^{(N)}_{i,j}(z;z_{low},z_{high}\f$ are normalized phase space functions
     *    for mass of \f$i\f$-particles from \f$j\f$-body decays, modulated by
     *    1D positive bernstein polynomial of degree \f$N\f$
     *  - \f$m_{max}\f$ is a maximal allowed mass for \f$x+y\f$
     *  
     *  @see Ostap::Math::PhaseSpacePol
     *  @see Ostap::Math::PhaseSpaceNL 
     *  @see Ostap::Math::Positive1D
     *  @see Ostap::Math::PS2DPol2
     */
    class  PS2DPol3 
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the order
      PS2DPol3 
      ( const PhaseSpaceNL&   psx  = PhaseSpaceNL () ,
        const PhaseSpaceNL&   psy  = PhaseSpaceNL () ,
        const double          mmax =  -1 ,  
        const unsigned short  Nx   =   1 ,
        const unsigned short  Ny   =   1 ) ;
      /// constructor from the order
      PS2DPol3
      ( const PhaseSpacePol&  psx        ,
        const PhaseSpacePol&  psy        ,
        const double          mmax =  -1 ) ;
      /// constructor from the order
      PS2DPol3
      ( const PhaseSpaceNL&   psx        ,
        const PhaseSpaceNL&   psy        ,
        const double          mmax       ,  
        const unsigned short  Nx         ,
        const unsigned short  Ny         ,
        const double          xmin       ,
        const double          xmax       ,
        const double          ymin       ,
        const double          ymax       ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      inline double operator ()
      ( const double x ,
	const double y ) const { return evaluate ( x , y ) ; } 
      /// get the value
      double evaluate
      ( const double x ,
	const double y ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get number of parameters
      inline std::size_t npars () const { return m_psx.npars () + m_psy.npars() ; }
      /// set k-parameter
      inline bool setPar       ( const unsigned int k , const double value )
      { 
        const unsigned int npx = m_psx.npars() ;
        return 
          k < npx     ? 
          m_psx.setPar ( k       , value ) :
          m_psy.setPar ( k - npx , value ) ;
      }
      /// set k-parameter
      inline bool setParameter ( const unsigned int k , const double value )
      { return setPar   ( k , value ) ; }
      /// get the parameter value
      inline double  par       ( const unsigned int k ) const
      {
        const unsigned int npx = m_psx.npars() ;
        return 
          k < npx     ? 
          m_psx.par ( k       ) :
          m_psy.par ( k - npx ) ;
      }
      /// get the parameter value
      inline double  parameter ( const unsigned int k ) const { return par ( k ) ; }
      /// get parameters/phases 
      std::vector<double> pars() const ;
      /// x-parameters 
      // const std::vector<double>& xpars() const { return m_psx.pars() ; }
      /// y-parameters 
      // const std::vector<double>& ypars() const { return m_psy.pars() ; }      
      /// get nX & nY
      inline unsigned short nX () const { return m_psx.degree () ; }
      inline unsigned short nY () const { return m_psy.degree () ; }
      // ======================================================================
    public:
      // ======================================================================
      inline double xmin () const { return m_psx.xmin () ; }
      inline double xmax () const { return m_psx.xmax () ; }
      inline double ymin () const { return m_psy.xmin () ; }
      inline double ymax () const { return m_psy.xmax () ; }
      inline double mmax () const { return m_mmax        ; }
      // ======================================================================
    public:
      // ======================================================================
      /** get the integral over 2D-region
       *  \f[ \int_{x_low}^{x_high}\int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y\f]
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       *  @param ylow  low  edge in y
       *  @param yhigh high edge in y
       */
      double integral
      ( const double xlow  ,
        const double xhigh ,
        const double ylow  , 
        const double yhigh ) const ;
      // ======================================================================
      /** integral over x-dimension
       *  \f[ \int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param y     variable
       *  @param xlow  low  edge in y
       *  @param xhigh high edge in y
       */
      double integrateX 
      ( const double y     ,
        const double xlow  , 
        const double xhigh ) const ;
      // ======================================================================
      /** integral over x-dimension
       *  \f[ \int_{x_low}^{x_high} \mathcal{B}(x,y) \mathrm{d}x\f]
       *  @param x     variable
       *  @param ylow  low  edge in x
       *  @param yhigh high edge in x
       */
      double integrateY 
      ( const double x     ,
        const double ylow  ,
        const double yhigh ) const ;
      // ======================================================================
    public:
      // ======================================================================
      const Ostap::Math::PhaseSpacePol& psX         () const { return m_psx  ; }
      const Ostap::Math::PhaseSpacePol& psY         () const { return m_psy  ; }
      const Ostap::Math::PhaseSpacePol& phasespaceX () const { return psX () ; }
      const Ostap::Math::PhaseSpacePol& phasespaceY () const { return psY () ; }
      // ====================================== ===============================
    public:
      // ====================================== ===============================
      /// get the tag 
      std::size_t tag () const ;
      // ====================================== ===============================
    private:
      // ======================================================================
      /// Phase space
      Ostap::Math::PhaseSpacePol m_psx      ;
      Ostap::Math::PhaseSpacePol m_psy      ;
      /// Max-mass 
      double                     m_mmax     ;
      // ======================================================================
    private:
      // ======================================================================
      /// Phase space
      mutable Ostap::Math::PhaseSpaceNL m_psx_aux ;
      mutable Ostap::Math::PhaseSpaceNL m_psy_aux ;
      // ======================================================================
      /// workspace
      Ostap::Math::WorkSpace m_workspace   ;
      // ======================================================================
    };
    // ========================================================================
    /** @class PS2DPol3Sym
     *  
     *  The symmetric 2D-function, that represent non-factorizeable "product" of  
     *  two modulated phase-space functions.
     *  It is a  symmetrized version of Ostap::Math::PS2DPol3
     *  The  function is useful to describe e.g. 2D-distributions of 
     *  \f$ m_{23}\f$ vs \f$m_{45}\f$ from 5-body decays. 
     *
     *  The function is:
     *  \f[ f(x,y) = \frac{1}{2}
     *   \left( \Phi^{(N)}_{k,n}(x;x_{low},x_{high}) \Phi_{k,n-1}(y,y_{low},m_{max}-x) 
     *        + \Phi^{(N)}_{k,n}(y;y_{low},y_{high}) \Phi_{k,n-1}(x,x_{low},m_{max}-y) 
     *   \right) \f]
     *  where 
     *  - \f$ \Phi_{i,j}(z;z_{low},z_{high}\f$ are normalized phase space functions
     *    for mass of \f$i\f$-particles from \f$j\f$-body decays
     *  - \f$\Phi^{(N)}_{i,j}(z;z_{low},z_{high}\f$ are normalized phase space functions
     *    for mass of \f$i\f$-particles from \f$j\f$-body decays, modulated by
     *    1D positive benrstein polynomial of degree \f$N\f$
     *  - \f$m_{max}\f$ is a maximal allowed mass for \f$x+y\f$
     *
     *  Clearly the function is symmetric:  \f$f(y,x)=f(x,y)\f$
     *  @see Ostap::Math::PhaseSpacePol
     *  @see Ostap::Math::PhaseSpaceNL
     *  @see Ostap::Math::Positive
     *  @see Ostap::Math::PS2DPol3
     */
    class  PS2DPol3Sym
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the order
      PS2DPol3Sym 
      ( const PhaseSpaceNL&   ps   = PhaseSpaceNL () ,
        const double          mmax =  -1 ,  
        const unsigned short  N    =   1 ) ;
      /// constructor from the order
      PS2DPol3Sym
      ( const PhaseSpacePol&  ps         ,
        const double          mmax =  -1 ) ;
      /// constructor from the order
      PS2DPol3Sym
      ( const PhaseSpaceNL&   ps         ,
        const double          mmax       ,  
        const unsigned short  N          ,
        const double          xmin       ,
        const double          xmax       );
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      inline double operator ()
      ( const double x ,
	const double y ) const { return evaluate ( x , y ) ; }
      /// get the value
      double evaluate
      ( const double x ,
	const double y ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get number of parameters
      inline std::size_t npars () const { return m_ps.npars ()  ; }
      /// set k-parameter
      inline bool setPar       ( const unsigned int k , const double value )
      { 
        const unsigned int np = m_ps.npars() ;
        return k < np ? m_ps.setPar ( k , value ) : m_ps.setPar ( k - np , value ) ;
      }
      /// set k-parameter
      inline bool setParameter ( const unsigned int k , const double value )
      { return setPar   ( k , value ) ; }
      /// get the parameter value
      inline double  par       ( const unsigned int k ) const
      {
        const unsigned int np = m_ps.npars() ;
        return k < np ? m_ps.par ( k  ) : m_ps.par ( k - np ) ;
      }
      /// get the parameter value
      inline double  parameter ( const unsigned int k ) const { return par ( k ) ; }
      /// get parameters/phases 
      // const std::vector<double>& pars () const { return m_ps.pars() ; }
      /// get nX & nY
      inline unsigned short n  () const { return m_ps.degree () ; }
      inline unsigned short nX () const { return m_ps.degree () ; }
      inline unsigned short nY () const { return m_ps.degree () ; }
      // ======================================================================
    public:
      // ======================================================================
      inline double xmin () const { return m_ps.xmin () ; }
      inline double xmax () const { return m_ps.xmax () ; }
      inline double ymin () const { return m_ps.xmin () ; }
      inline double ymax () const { return m_ps.xmax () ; }
      inline double mmax () const { return m_mmax       ; }
      // ======================================================================
    public:
      // ======================================================================
      /** get the integral over 2D-region
       *  \f[ \int_{x_low}^{x_high}\int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y\f]
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       *  @param ylow  low  edge in y
       *  @param yhigh high edge in y
       */
      double integral
      ( const double xlow  ,
        const double xhigh ,
        const double ylow  , 
        const double yhigh ) const ;
      // ======================================================================
      /** integral over x-dimension
       *  \f[ \int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param y     variable
       *  @param xlow  low  edge in y
       *  @param xhigh high edge in y
       */
      double integrateX 
      ( const double y     ,
        const double xlow  , 
        const double xhigh ) const ;
      // ======================================================================
      /** integral over x-dimension
       *  \f[ \int_{x_low}^{x_high} \mathcal{B}(x,y) \mathrm{d}x\f]
       *  @param x     variable
       *  @param ylow  low  edge in x
       *  @param yhigh high edge in x
       */
      double integrateY
      ( const double x     ,
        const double ylow  ,
        const double yhigh ) const ;
      // ======================================================================
    public:
      // ======================================================================
      const Ostap::Math::PhaseSpacePol& ps          () const { return m_ps   ; }
      const Ostap::Math::PhaseSpacePol& psX         () const { return m_ps   ; }
      const Ostap::Math::PhaseSpacePol& psY         () const { return m_ps   ; }
      const Ostap::Math::PhaseSpacePol& phasespace  () const { return ps  () ; }
      const Ostap::Math::PhaseSpacePol& phasespaceX () const { return psX () ; }
      const Ostap::Math::PhaseSpacePol& phasespaceY () const { return psY () ; }
      // ====================================== ===============================
    public:
      // ====================================== ===============================
      /// get the tag 
      std::size_t tag () const ;
      // ====================================== ===============================
    private:
      // ======================================================================
      /// Phase space
      Ostap::Math::PhaseSpacePol m_ps  ;
      /// Max-mass 
      double                    m_mmax ;
      // ======================================================================
    private:
      // ======================================================================
      /// Phase space
      mutable Ostap::Math::PhaseSpaceNL m_psx_aux ;
      mutable Ostap::Math::PhaseSpaceNL m_psy_aux ;
      // ======================================================================
      /// workspace
      Ostap::Math::WorkSpace m_workspace   ;
      // ======================================================================
    };
    // ========================================================================
    /** @class ExpoPS2DPol
     *  The 2D-function:
     *  \f$ f(x,y) = exp(tau*x)*Ps(y)*P_{pos}(x,y) \f$, where
     *  \f$Ps\f$ denotes phase-space function and
     * \f$P_{pos}\f$ denotes the positive polynomial
     */
    class  ExpoPS2DPol : public Ostap::Math::PolyFactor2D 
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the order
      ExpoPS2DPol 
      ( const PhaseSpaceNL&   psy  = PhaseSpaceNL() ,
        const double          xmin = 0 ,
        const double          xmax = 1 ,
        const unsigned short  Nx   = 1 ,
        const unsigned short  Ny   = 1 ,
        const double          tau  = 0 ) ;
      /// constructor from the order
      ExpoPS2DPol 
      ( const PhaseSpaceNL&   psy      ,
        const double          xmin     ,
        const double          xmax     ,
        const unsigned short  Nx       ,
        const unsigned short  Ny       ,
        const double          ymin     ,
        const double          ymax     ,
        const double          tau  = 0 ) ;
      // constructor from components
      ExpoPS2DPol 
      ( const Ostap::Math::Positive2D& pol     , 
        const PhaseSpaceNL&            ps      ,        
        const double                   tau = 0 ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      inline double operator ()
      ( const double x ,
	const double y ) const { return evaluate ( x , y ) ; } 
      /// get the value
      double evaluate 
      ( const double x ,
	const double y ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get number of parameters
      inline std::size_t npars () const { return m_positive.npars () ; }
      /// set k-parameter
      inline  bool setPar       ( const unsigned int k , const double value )
      { return m_positive.setPar ( k , value ) ;}
      /// set k-parameter
      inline bool setParameter ( const unsigned int k , const double value )
      { return setPar   ( k , value ) ; }
      /// get the parameter value
      inline double  par       ( const unsigned int k ) const
      { return m_positive.par ( k ) ; }
      /// get the parameter value
      inline double  parameter ( const unsigned int k ) const { return par ( k ) ; }
      /// get nX & nY
      inline unsigned short nX () const { return m_positive.nX () ; }
      inline unsigned short nY () const { return m_positive.nY () ; }
      // ======================================================================
    public:
      // ======================================================================
      inline double xmin () const { return m_positive.xmin () ; }
      inline double xmax () const { return m_positive.xmax () ; }
      inline double ymin () const { return m_positive.ymin () ; }
      inline double ymax () const { return m_positive.ymax () ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get tau
      inline double  tau () const { return m_tau ;}
      /// set tau
      bool           setTau ( const double val ) ;
      // ======================================================================
    public:
      // ======================================================================
      /** get the integral over 2D-region
       *  \f[ \int_{x_low}^{x_high}\int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y\f]
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       *  @param ylow  low  edge in y
       *  @param yhigh high edge in y
       */
      double integral
      ( const double xlow  ,
        const double xhigh ,
        const double ylow  , 
        const double yhigh ) const ;
      // ======================================================================
      /** integral over x-dimension
       *  \f[ \int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param y     variable
       *  @param xlow  low  edge in y
       *  @param xhigh high edge in y
       */
      double integrateX 
      ( const double y     ,
        const double xlow  ,
        const double xhigh ) const ;
      // ======================================================================
      /** integral over x-dimension
       *  \f[ \int_{x_low}^{x_high} \mathcal{B}(x,y) \mathrm{d}x\f]
       *  @param x     variable
       *  @param ylow  low  edge in x
       *  @param yhigh high edge in x
       */
      double integrateY 
      ( const double x     ,
        const double ylow  ,
        const double yhigh ) const ;
      // ======================================================================
    public:
      // ======================================================================
      const Ostap::Math::PhaseSpaceNL& psY         () const { return m_psy      ; }
      const Ostap::Math::PhaseSpaceNL& phasespaceY () const { return psY ()     ; }
      // ====================================== ===============================
    public:
      // ====================================== ===============================
      /// get the tag 
      std::size_t tag () const ;
      // ====================================== ===============================
    private: // helper function to make the calculations
      // ======================================================================
      /// helper function to make calculations
      double calculate
      ( const std::vector<double>& fx , 
	const std::vector<double>& fy ) const ;
      // ======================================================================
    private:
      // ======================================================================
      /// Phase space
      Ostap::Math::PhaseSpaceNL m_psy      ;
      /// exponential
      double                    m_tau      ;
      // ======================================================================
    private:
      // ======================================================================
      /// workspace
      Ostap::Math::WorkSpace m_workspace   ;
      // ======================================================================
    };
    // ========================================================================
    /** @class Expo2DPol
     *  The 2D-function:
     *  \f$ f(x,y) = exp(x)*expo(y)*P_{pos}(x,y) \f$, where
     * \f$P_{pos}\f$ denotes the positive polynomial
     */
    class  Expo2DPol : public Ostap::Math::PolyFactor2D 
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the order
      Expo2DPol
      ( const double          xmin = 0  ,
        const double          xmax = 1  ,
        const double          ymin = 0  ,
        const double          ymax = 1  ,
        const unsigned short  Nx   =  1 ,
        const unsigned short  Ny   =  1 ,
        const double          taux =  0 ,
        const double          tauy =  0 ) ;            
      /// constructor from polynomial 
      Expo2DPol
      ( const Ostap::Math::Positive2D& pol       , 
        const double                   taux =  0 ,
        const double                   tauy =  0 ) ;            
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      inline double operator ()
      ( const double x ,
	const double y ) const { return evaluate ( x , y ) ; }
      /// get the value
      double evaluate
      ( const double x ,
	const double y ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get tau
      inline double  tauX    () const { return m_tauX ;}
      inline double  tauY    () const { return m_tauY ;}
      /// set tau
      bool           setTauX ( const double val ) ;
      bool           setTauY ( const double val ) ;
      // ======================================================================
    public:
      // ======================================================================
      /** get the integral over 2D-region
       *  \f[ \int_{x_low}^{x_high}\int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y\f]
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       *  @param ylow  low  edge in y
       *  @param yhigh high edge in y
       */
      double integral
      ( const double xlow  ,
        const double xhigh ,
        const double ylow  , 
        const double yhigh ) const ;
      // ======================================================================
      /** integral over x-dimension
       *  \f[ \int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param y     variable
       *  @param xlow  low  edge in y
       *  @param xhigh high edge in y
       */
      double integrateX
      ( const double y     ,
        const double xlow  ,
        const double xhigh ) const ;
      // ======================================================================
      /** integral over x-dimension
       *  \f[ \int_{x_low}^{x_high} \mathcal{B}(x,y) \mathrm{d}x\f]
       *  @param x     variable
       *  @param ylow  low  edge in x
       *  @param yhigh high edge in x
       */
      double integrateY 
      ( const double x     ,
        const double ylow  , 
        const double yhigh ) const ;
      // ======================================================================
    public:
      // ====================================== ===============================
      /// get the tag 
      std::size_t tag () const ;
      // ====================================== ===============================
    private: // helper function to make the calculations
      // ======================================================================
      /// helper function to make calculations
      double calculate
      ( const std::vector<double>& fx , 
	const std::vector<double>& fy ) const ;
      // ======================================================================
    private:
      // ======================================================================
      /// exponential factors  
      double                    m_tauX     ;
      double                    m_tauY     ;
      // ======================================================================
    };
    // ========================================================================
    /** @class Expo2DPolSym
     *  The 2D-function:
     *  \f$ f(x,y) = exp(x)*expo(y)*P_{sym}(x,y) \f$, where
     * \f$P_{pos}\f$ denotes the symmetric positive polynomial
     */
    class  Expo2DPolSym : public Ostap::Math::PolyFactor2DSym 
    {
      // ======================================================================
    public:
      // ======================================================================
      /// constructor from the order
      Expo2DPolSym
      ( const double          xmin = 0 ,
        const double          xmax = 1 ,
        const unsigned short  N    = 1 ,
        const double          tau  = 0 ) ;
      // ======================================================================
      /// constructor from polynomial 
      Expo2DPolSym
      ( const Ostap::Math::Positive2DSym& pol     , 
        const double                      tau = 0 ) ;            
      // ======================================================================
    public:
      // ======================================================================
      /// get the value
      inline double operator ()
      ( const double x ,
	const double y ) const { return evaluate ( x , y ) ; }
      /// get the value
      double evaluate
      ( const double x ,
	const double y ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get tau
      inline double  tau     () const { return m_tau  ;}
      /// set tau
      bool           setTau  ( const double val ) ;
      // ======================================================================
    public:
      // ======================================================================
      /** get the integral over 2D-region
       *  \f[ \int_{x_low}^{x_high}\int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y\f]
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       *  @param ylow  low  edge in y
       *  @param yhigh high edge in y
       */
      double integral 
      ( const double xlow , const double xhigh ,
        const double ylow , const double yhigh ) const ;
      // ======================================================================
      /** integral over x-dimension
       *  \f[ \int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param y     variable
       *  @param xlow  low  edge in y
       *  @param xhigh high edge in y
       */
      double integrateX 
      ( const double y    ,
        const double xlow , const double xhigh ) const ;
      // ======================================================================
      /** integral over y-dimension
       *  \f[ \int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param x     variable
       *  @param ylow  low  edge in x
       *  @param yhigh high edge in x
       */
      double integrateY 
      ( const double x    ,
        const double ylow , const double yhigh ) const ;
      // ======================================================================
    public:
      // ====================================== ===============================
      /// get the tag 
      std::size_t tag () const ;
      // ====================================== ===============================
    private: // helper function to make the calculations
      // ======================================================================
      /// helper function to make calculations
      double calculate 
      ( const std::vector<double>& fx , 
        const std::vector<double>& fy ) const ;
      // ======================================================================
    private:
      // ======================================================================
      /// exponential
      double                     m_tau      ;
      // ======================================================================
    }; 
    // ========================================================================
    /** @class Gauss2D 
     *  Simple 2D rotated Gaussian function 
     *  \f[ f(x,y) = \frac{1}{2\pi\sqrt{\sigma_x}\sqrt{\sigma_y}}
     *   \mathrm{e}^{\left( -\frac{1}{2} \delta x^{\prime2}}{\sigma_x^2} 
     *   -\frac{1}{2} \delta y^{\prime2}}{\sigma_y^2} 
     *   \right) }\f]
     * - where 
     *   \f[ \left( \begin{array}{c} \delta x^prime \\ \delta y^\prime \end{array} \right) = 
     *    \left( \begin{array}{cc} \cos \theta & \sin\theta \		\
     *                             -\sin\theta & \cos\theta \end{array} \right) \times   
     *     \left( \begin{array}{c} x - \mu_x \\ y - \mu_y\end{array} \right)\f] 
     * 
     *  @date 2022-06-22
     */
    class Gauss2D 
    {
      // ======================================================================
    public: 
      // ======================================================================
      /// constructor 
      Gauss2D 
      ( const double muX    = 0 , 
        const double muY    = 0 , 
        const double sigmaX = 1 , 
        const double sigmaY = 1 ,
        const double theta  = 0 ) ;
      // =====================================================================
    public:
      // ======================================================================
      /// get the value
      double operator () ( const double x , const double y ) const ;
      // ======================================================================
    public:  // getters 
      // ======================================================================      
      /// mux                                    
      inline double muX    () const { return m_muX    ; }
      /// mux                                   
      inline double muY    () const { return m_muY    ; }
      /// sigmax                                    
      inline double sigmaX () const { return m_sigmaX ; }
      /// sigmay                                   
      inline double sigmaY () const { return m_sigmaY ; }
      /// theta 
      inline double theta  () const { return m_theta  ; }
      // ======================================================================
    public:  // setters 
      // ======================================================================      
      bool setMuX      ( const double value ) ;
      bool setMuY      ( const double value ) ;
      bool setSigmaX   ( const double value ) ;
      bool setSigmaY   ( const double value ) ;
      bool setTheta    ( const double value ) ;
      // ======================================================================      
    public:
      // ======================================================================      
      /// get the integral over the whole 2D-region
      double integral () const ;
      // ======================================================================      
      /** get the integral over 2D-region
       *  \f[ \int_{x_low}^{x_high}\int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}x\mathrm{d}y\f]
       *  @param xlow  low  edge in x
       *  @param xhigh high edge in x
       *  @param ylow  low  edge in y
       *  @param yhigh high edge in y
       */
      double integral 
      ( const double xlow  ,
        const double xhigh ,
        const double ylow  , 
        const double yhigh ) const ;
      // ======================================================================
      /** integral over x-dimension
       *  \f[ \int_{y_low}^{y_high} \mathcal{B}(x,y) \mathrm{d}y\f]
       *  @param y     variable
       *  @param xlow  low  edge in y
       *  @param xhigh high edge in y
       */
      double integrateX 
      ( const double y     ,
        const double xlow  , 
        const double xhigh ) const ;
      // ======================================================================
      /** integral over x-dimension
       *  \f[ \int_{x_low}^{x_high} \mathcal{B}(x,y) \mathrm{d}x\f]
       *  @param x     variable
       *  @param ylow  low  edge in x
       *  @param yhigh high edge in x
       */
      double integrateY 
      ( const double x     ,
        const double ylow  ,
        const double yhigh ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the tag 
      std::size_t tag () const ;
      // ======================================================================
    private: 
      // ======================================================================
      /// mux 
      double m_muX      { 0.0 } ; // mux
      /// muy 
      double m_muY      { 0.0 } ; // muy 
      /// sigmax 
      double m_sigmaX   { 1.0 } ; // sigmax 
      /// sigmay 
      double m_sigmaY   { 1.0 } ; // sigmay 
      /// theta 
      double m_theta    { 0.0 } ; // theta 
      /// sin(theta) 
      double m_sintheta { 0.0 } ; // sin(theta)
      /// cos(theta) 
      double m_costheta { 1.0 } ; // sin(theta)      
      // ======================================================================
    private:
      // ======================================================================
      /// workspace for numerical integration
      Ostap::Math::WorkSpace m_workspace {} ;
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Tsallis2 
     *  2D particle density distribution as function of pt and rapidity 
     *  @see L. Marques, J. Cleymans, A. Deppman, 
     *       "Description of High-Energy pp Collisions 
     *        Using Tsallis Thermodynamics: 
     *        Transverse Momentum and Rapidity Distributions", 
     *        Phys. Rev. D 91, 054025, 	arXiv:1501.00953 
     *  @see https://arxiv.org/abs/1501.00953
     *  @see https://doi.org/10.1103/PhysRevD.91.054025
     *  @see Ostap::Math::Tsallis
     */ 
    class Tsallis2 
    {
      // ======================================================================
    public:
      // ======================================================================
      /** constructor from all parameters 
       *  @param mass particle mass (needed to calculate transverse mass)
       *  @param q    q-parameter of Tsallis  (q=1 corresponds to Boltzman statistics)
       *  @param T    the temperature
       *  @param mu   chemical potential 
       */
      Tsallis2 
      ( const double mass = 1     ,   // mass
        const double T    = 0.050 ,   // temperature
        const double q    = 1.1   ,   // q=1 -> Boltzman statistics 
        const double mu   = 0     ) ; // chemical potential 
      // ======================================================================
    public:
      // ======================================================================
      /** evaluate Tsallis function
       *  @param pt transverse momentum of the particle 
       *  @param y  rapidity of the particle 
       */
      double evaluate 
      ( const double pt , const double y  ) const ;
      /** evaluate Tsallis functuon
       *  @param pt transverse momentum of the particle 
       *  @param y  rapidity of the particle 
       */
      double operator() 
      ( const double pt , const double y  ) const 
      { return evaluate ( pt , y ) ; }
      // ======================================================================
     public:
      // ======================================================================
      /// particle mass 
      double mass () const { return m_mass ; }
      /// q-parameter  
      double q    () const { return m_q    ; }
      /// temperature 
      double T    () const { return m_T    ; }
      /// chemical potential 
      double mu   () const { return m_mu   ; }
      // ======================================================================
    public:
      // ======================================================================
      /// update mass-parameter
      bool setMass ( const double value ) ; // update mass-parameter
      bool setM    ( const double value ) { return setMass ( value ) ; }
      /// update q-parameter
      bool setQ    ( const double value ) ; // update q-parameter
      /// update temperature 
      bool setT    ( const double value ) ; // update temperature 
      /// update chemiclas potyential 
      bool setMu   ( const double value ) ; // update chemical potential 
      // ======================================================================
    public:
      // ======================================================================
      /// get the transverse mass for the given pt 
      inline double mT   ( const double pt  ) const
      { return std::hypot ( pt , m_mass  ) ; }
      // ======================================================================
    public:
      // ======================================================================
      /// get the integral between ptlow-pthigh and ylow-yhigh
      double integral 
      ( const double ptlow  ,
        const double pthigh , 
        const double ylow   , 
        const double yhigh  ) const ;
      /// get the integral between ylow-yhigh for given pt 
      double integrate_y  
      ( const double pt ,
        const double ylow   , 
        const double yhigh  ) const ;
      /// get the integral between ptlow-pthigh for given rapidity 
      double integrate_pt  
      ( const double y      , 
        const double ptlow  ,
        const double pthigh ) const ;
      // ======================================================================
    public:
      // ======================================================================
      // get the tag
      std::size_t tag () const ;
      // ======================================================================
    private :
      // ======================================================================
      /// mass of the particle 
      double m_mass  { 1.0   } ; // mass of the particle 
      /// temperature 
      double m_T     { 0.050 } ; // temperature
      /// q-parameter 
      double m_q     { 1.1   } ; // q-parameter 
      /// chemical potential 
      double m_mu    { 0     } ; // chemical potential
      // ======================================================================
    private :
      // ======================================================================
      /// workspace
      Ostap::Math::WorkSpace m_workspace {} ; // workspace
      // ======================================================================
    } ;
    // ========================================================================    
  } //                                             end of namespace Ostap::Math
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
//                                                                      The END
// ============================================================================
#endif // OSTAP_MODELS2D_H
// ============================================================================
