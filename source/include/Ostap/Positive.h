#ifndef OSTAP_POSITIVE_H 
#define OSTAP_POSITIVE_H 1
// =============================================================================
// Include files
// =============================================================================
// Ostap
// =============================================================================
#include "Ostap/NSphere.h"
#include "Ostap/Workspace.h"
// =============================================================================
/** @file  Ostap/Positive.h
 *  Implemenetation of positive polynomials
 *  - A non-negative polynomial on the [a,b] interval 
 *  @see S.Karlin, L.S.Shapley, "Geometry of Moment Space", 
 *       Mem. Am. Math. Soc., 12, 1953 
 *  - A non-negative polynomial on the [ x_0,+infinityL interval
 *  @see S.Karlin, W.J.Studden, "Tchebysheff systems: with application
 *       in analysis and statistics", 
 *
 * @see https://bookstore.ams.org/view?ProductCode=MEMO/1/12 
 * @see https://www.scirp.org/(S(351jmbntvnsjt1aadkposzje))/reference/ReferencesPapers.aspx?ReferenceID=1762140
 */  
namespace Ostap 
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    /** @class KarlinShapley 
     *  A non-negative polynomial on the [a,b] interval 
     *  @see S.Karlin, L.S.Shapley, "Geometry of Moment Space", 
     *       Mem. Am. Math. Soc., 12, 1953  
     */
    class KarlinShapley 
    {
      // ======================================================================
    public: 
      // ======================================================================
      /** constructor from the order and interval 
       *  @param   N polynomial degree
       *  @param   xmin minimal x 
       *  @param   xmax maximal x 
       */  
      KarlinShapley
      ( const unsigned short N     = 0 , 
        const double         xmin  = 0 , 
        const double         xmax  = 1 ) ;
      // =======================================================================
      /** constructor from parameters 
       *  @param   pars parameters 
       *  @param   xmin minimal x 
       *  @param   xmax maximal x 
       */  
      KarlinShapley
      ( const std::vector<double>& pars      , 
        const double               xmin  = 0 , 
        const double               xmax  = 1 ) ;
      // =========================================================================
      /** constructor from the scale and phases 
       *  @param   A        global  scale 
       *  @param   phases1  phases of the 1st sphere 
       *  @param   phases2 phases  of the 2nd sphere 
       *  @param   xmin minimal x 
       *  @param   xmax maximal x 
       */  
      KarlinShapley
      ( const double               A         , 
        const std::vector<double>& phases1   , 
        const std::vector<double>& phases2   , 
        const double               xmin  = 0 , 
        const double               xmax  = 1 ) ;
      // ==========================================================================
      /** constructor from the scale, phase and phases 
       *  @param   A        global  scale 
       *  @param   phi      phase of the 1st sphere 
       *  @param   phases2  phases  of the 2nd sphere 
       *  @param   xmin minimal x 
       *  @param   xmax maximal x 
       */  
      KarlinShapley
      ( const double               A         , 
        const double               phi       , 
        const std::vector<double>& phases2   , 
        const double               xmin  = 0 , 
        const double               xmax  = 1 ) ;
      // ======================================================================
    public: 
      // ======================================================================
      /// evaluate Karlin-Shapley polynomial 
      inline double operator() ( const double x ) const 
      { return evaluate ( x ) ; }
      /// evaluate Karlin-Shapley polynomial 
      double evaluate (  const double x ) const ;
      // ======================================================================
    public: 
      // ======================================================================
      /// get degree of polynomial
      unsigned short N      () const { return m_sphere1.npars () + m_sphere2.npars () ; }
      /// get degree of polynomial
      unsigned short degree () const { return m_sphere1.npars () + m_sphere2.npars () ; }
      /// get degree of polynomial
      unsigned short order  () const { return m_sphere1.npars () + m_sphere2.npars () ; }
      /// get minimal x 
      double         xmin   () const { return m_xmin  ; }
      /// get maximal  x 
      double         xmax   () const { return m_xmax  ; }
      // ======================================================================
    public: 
      // ======================================================================
      /// number of parameters:  #sphere parameters + 2 
      inline unsigned short npars () const 
      { return 1 + m_sphere1.npars() + m_sphere2.npars() ; }
      /// global coefficient 
      double A () const { return m_A ; }
      /// coeffcienct at alpha-polynomial 
      double alpha () const { return m_A * m_sphere1.x2 ( 0 ) ; }
      /// coeffcienct at beta-polynomial
      double beta  () const { return m_A * m_sphere1.x2 ( 1 ) ; }
      /// get first phase parameters 
      const std::vector<double>& phases1 () const { return m_sphere1.pars () ; }
      /// get seconf phase-parameters 
      const std::vector<double>& phases2 () const { return m_sphere2.pars () ; }
      /** get parameter (phase)  by index 
       *  - k = 0 : parameter A 
       *  - k = 1 : phase from 1st sphere (if exists) 
       *  - 1 < k : phase from the sphere 
       */
      inline double   par ( const unsigned short k ) const 
      {
        const unsigned short n0 = 1  ;
        const unsigned short n1 = n0 + m_sphere1.npars () ;
        const unsigned short n2 = n1 + m_sphere2.npars () ;
        return 
          k < n0 ? m_A : 
          k < n1 ? m_sphere1.par ( k - 1 ) : m_sphere2.par ( k - n1  ) ;
      }
      /// 
      // ======================================================================
    public:
      // ======================================================================
      /** set parameter (phase)  by index 
       *  - k = 0 : parameter A 
       *  - k = 1 : phase from the 1st sphere 
       *  - 1 < k : phase from the 2nd sphere 
       */
      inline bool setPar
      ( const unsigned short k      , 
        const double         value  )
      {
        //
        const unsigned short n0 = 1  ;
        const unsigned short n1 = n0 + m_sphere1.npars () ;
        const unsigned short n2 = n1 + m_sphere2.npars () ;
        //
        if       ( 0 == k  ) { return setA ( value ) ; }
        else if  ( k <  n1 ) { return m_sphere1.setPar ( k - 1  , value ) ; }
        //
        const bool updated = m_sphere2.setPar ( k - n1 , value ) ;
        if ( updated ) { updateRoots() ; }
        //
        return updated ;
      }
      // ======================================================================
      /// set parameter alpha 
      bool setA ( const double value ) ;
      // ======================================================================
    public: 
      // ======================================================================
      /// get the internal variable from external
      inline double t ( const double x ) const 
      { return ( x - m_xmin ) / ( m_xmax - m_xmin ) ; }
      /// get the external variable from internal  
      inline double x ( const double t ) const 
      { return m_xmin + t * ( m_xmax - m_xmin ) ; }
      // ======================================================================
    public: 
      // ======================================================================
      /// get Karlin-Shapley t-roots 
      const std::vector<double>& troots () const { return m_troots ; }
      // ======================================================================
    public:  // numerical intergation
      // ======================================================================
      /** get (numerical) integral 
       * \f$ f = \int^{x_{max}}_{x_{min}} P(x) dx \f$ 
       */
      double integral 
      ( const double xmin , 
        const double xmax ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the unique tag 
      std::size_t tag () const ; // get the unique tag      
      // ======================================================================
    public: // swap two polynomials 
      // ======================================================================
      /// swap two polynomials
      void swap ( KarlinShapley& right ) ; // swap two polynomials
      // ======================================================================
    private:
      // ======================================================================
      /// update internal Karlin-Shapley t-roots 
      void updateRoots () ; // update internal Karlin-Shapley t-roots 
      // ======================================================================
    private:
      // ======================================================================
      /// minimal x 
      double                      m_xmin   { 0 } ; // minimal x 
      /// maximal  x 
      double                      m_xmax   { 1 } ; // maximal x 
      /// A 
      double                      m_A      { 1 } ; // global coefficient
      // ======================================================================
    private: 
      // ======================================================================
      /// 1st n-sphere 
      Ostap::Math::NSphere    m_sphere1   {} ; // n-sphere
      /// 2nd n-sphere 
      Ostap::Math::NSphere    m_sphere2   {} ; // n-sphere
      /// Karlin-Shapley t-roots
      std::vector<double>     m_troots    {} ; // Karlin-Shapley t-roots 
      /// integration workspace for numerical integration  
      Ostap::Math::WorkSpace  m_workspace {} ;
      // ======================================================================
    } ;
    // ========================================================================
    /// swap two Karlin-Shapley polynomials 
    inline void swap 
    ( KarlinShapley& a , 
      KarlinShapley& b ) { a.swap ( b ) ; }
    // ========================================================================
    /** @class KarlinStudden  
     *  A non-negative polynomial on the [ x_0,+infinityL interval
     *  @see S.Karlin, W.J.Studden, "Tchebysheff systems: with application
     *       in analysis and statistics", 
     * @see https://www.scirp.org/(S(351jmbntvnsjt1aadkposzje))/reference/ReferencesPapers.aspx?ReferenceID=1762140
     */
    class KarlinStudden 
    {
      // ======================================================================
    public:
      // ======================================================================
      /** constuctor from the order and interval 
       *  @param   N     polymnomial degree
       *  @param   xmin the point x_0 
       *  @param   scale global scale factor (for better numerical perfoamnce)
       */  
      KarlinStudden 
      ( const unsigned short N     = 0 , 
        const double         xmin  = 0 , 
        const double         scale = 1 ) ;
      // ======================================================================
      /** constructor from all parameters 
       *  @param   pars parameters 
       *  @param   xmin minimal x 
       *  @param   scale global scale factor (for better numerical perfoamnce)
       */  
      KarlinStudden
      ( const std::vector<double>& pars      , 
        const double               xmin  = 0 , 
        const double               scale = 1 ) ;
      // =========================================================================
      /** constructor from the scale and phases 
       *  @param   A        global  scale 
       *  @param   phases1  phases of the 1st sphere 
       *  @param   phases2 phases  of the 2nd sphere 
       *  @param   xmin minimal x 
       *  @param   scale global scale factor (for better numerical perfoamnce)
       */  
      KarlinStudden
      ( const double               A         , 
        const std::vector<double>& phases1   , 
        const std::vector<double>& phases2   , 
        const double               xmin  = 0 , 
        const double               scale = 1 ) ;
      // ==========================================================================
      /** constructor from the scale, phase and phases 
       *  @param   A        global  scale 
       *  @param   phi      phase of the 1st sphere 
       *  @param   phases2  phases  of the 2nd sphere 
       *  @param   xmin minimal x 
       *  @param   scale global scale factor (for better numerical perfoamnce)
       */  
      KarlinStudden
      ( const double               A          , 
        const double               phi        , 
        const std::vector<double>& phases2    , 
        const double               xmin   = 0 , 
        const double               scale  = 1 ) ;
      // ======================================================================
    public: 
      // ======================================================================
      /// evaluate Karlin-Studden  polynomial 
      inline double operator() ( const double x ) const 
      { return evaluate ( x ) ; }
      /// evaluate Karlin-Studden polynomial 
      double evaluate (  const double x ) const ;
      // ======================================================================
    public: 
      // ======================================================================
      /// get degree of polynomial
      unsigned short N      () const { return m_sphere1.npars () + m_sphere2.npars () ; }
      /// get degree of polynomial
      unsigned short degree () const { return m_sphere1.npars () + m_sphere2.npars () ; }
      /// get degree of polynomial
      unsigned short order  () const { return m_sphere1.npars () + m_sphere2.npars () ; }
      /// get minimal x 
      double         xmin   () const { return m_xmin  ; }
      /// get the scale 
      double         scale  () const { return m_scale ; }
      // ======================================================================
    public: 
      // ======================================================================
      /// number of parameters:  #sphere parameters + 2 
      inline unsigned short npars () const 
      { return 1 + m_sphere1.npars() + m_sphere2.npars() ; }
      /// global coefficient 
      double A () const { return m_A ; }
      /// coeffcienct at alpha-polynomial 
      double alpha () const { return m_A * m_sphere1.x2 ( 0 ) ; }
      /// coeffcienct at beta-polynomial
      double beta  () const { return m_A * m_sphere1.x2 ( 1 ) ; }
      /// get first phase parameters 
      const std::vector<double>& phases1 () const { return m_sphere1.pars () ; }
      /// get seconf phase-parameters 
      const std::vector<double>& phases2 () const { return m_sphere2.pars () ; }
      /** get parameter (phase)  by index 
       *  - k = 0 : parameter A 
       *  - k = 1 : phase from 1st sphere (if exists) 
       *  - 1 < k : phase from the sphere 
       */
      inline double   par ( const unsigned short k ) const 
      {
        const unsigned short n0 = 1  ;
        const unsigned short n1 = n0 + m_sphere1.npars () ;
        const unsigned short n2 = n1 + m_sphere2.npars () ;
        return 
          k < n0 ? m_A : 
          k < n1 ? m_sphere1.par ( k - 1 ) : m_sphere2.par ( k - n1  ) ;
      }
      /// 
      // ======================================================================
    public:
      // ======================================================================
      /** set parameter (phase)  by index 
       *  - k = 0 : parameter A 
       *  - k = 1 : phase from the 1st sphere 
       *  - 1 < k : phase from the 2nd sphere 
       */
      inline bool setPar
      ( const unsigned short k      , 
        const double         value  )
      {
        //
        const unsigned short n0 = 1  ;
        const unsigned short n1 = n0 + m_sphere1.npars () ;
        const unsigned short n2 = n1 + m_sphere2.npars () ;
        //
        if       ( 0 == k  ) { return setA ( value ) ; }
        else if  ( k <  n1 ) { return m_sphere1.setPar ( k - 1  , value ) ; }
        //
        const bool updated = m_sphere2.setPar ( k - n1 , value ) ;
        if ( updated ) { updateRoots() ; }
        //
        return updated ;
      }
      // ======================================================================
      /// set parameter alpha 
      bool setA ( const double value ) ;
      // ======================================================================
    public: 
      // ======================================================================
      /// get the internal variable from external
      inline double t ( const double x ) const 
      { return ( x - m_xmin ) / m_scale ; }
      /// get the external variable from internal  
      inline double x ( const double t ) const 
      { return   t * m_scale + m_xmin   ; }
      // ======================================================================
    public: 
      // ======================================================================
      /// get Karlin-Studden t-roots 
      const std::vector<double>& troots () const { return m_troots ; }
      // ======================================================================
    public:  // numerical intergation
      // ======================================================================
      /** get (numerical) integral 
       * \f$ f = \int^{x_{max}}_{x_{min}} P(x) dx \f$ 
       */
      double integral 
      ( const double xmin , 
        const double xmax ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the unique tag 
      std::size_t tag() const ; // get the unique tag      
      // ======================================================================
    public: // swap two polynomials 
      // ======================================================================
      /// swap two polynomials
      void swap ( KarlinStudden& right ) ; // swap two polynomials
      // ======================================================================
    private:
      // ======================================================================
      /// update internal Karlin-Studden  t-roots 
      void updateRoots () ; // update internal Karlin-Studden t-roots 
      // ======================================================================
    private:
      // ======================================================================
      /// minimal x 
      double                  m_xmin      { 0 } ; // minimal x 
      /// scale parameter
      double                  m_scale     { 1 } ; // scale parameter
      /// A 
      double                  m_A         { 1 } ; // global coefficient
      // ======================================================================
    private: 
      // ======================================================================
      /// 1st n-sphere 
      Ostap::Math::NSphere    m_sphere1   {} ; // n-sphere
      /// 2nd n-sphere 
      Ostap::Math::NSphere    m_sphere2   {} ; // n-sphere
      /// Karlin-Shapley t-roots
      std::vector<double>     m_troots    {} ; // Karlin-Shapley t-roots 
      /// integration workspace for numerical integration  
      Ostap::Math::WorkSpace  m_workspace {} ;
      // ======================================================================      
    } ;
    // ========================================================================
    /// swap two Karlin-Studdenpolynomials 
    inline void swap 
    ( KarlinStudden& a , 
      KarlinStudden& b ) { a.swap ( b ) ; }
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_POSITIVE_H
// ============================================================================
