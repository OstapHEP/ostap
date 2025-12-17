// ============================================================================
#ifndef OSTAP_HISTOINTERPOLATORS_H 
#define OSTAP_HISTOINTERPOLATORS_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL 
// ============================================================================
#include <array> 
// ============================================================================
// ROOT
// ============================================================================
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/HistoInterpolation.h"
#include "Ostap/Workspace.h"
// ============================================================================
namespace Ostap 
{
  // ==========================================================================
  namespace Math 
  {
    // ========================================================================
    /** @class HistoInterpolator HistoInterpolators.h Ostap/HistoInterpolators.h
     *  Simple function object/histogram interpolator 
     *  @see Ostap::Math::HistoInterpolation
     *  @see Ostap::Math::HistoInterpolation::interpolate_1D
     *  @see Ostap::Math::HistoInterpolation::interpolate_2D
     *  @see Ostap::Math::HistoInterpolation::interpolate_3D
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2019-11-16
     */
    class HistoInterpolator 
    {
    public:
      // ======================================================================
      /** constructor with full  specification 
       *  @see Ostap::Math::HistoInterpolation
       *  @see Ostap::Math::HistoInterpolation::interpolate_1D
       *  @see Ostap::Math::HistoInterpolation::interpolate_2D
       *  @see Ostap::Math::HistoInterpolation::interpolate_3D
       */
      HistoInterpolator
      ( const bool edges       = true  , 
        const bool extrapolate = false , 
        const bool density     = false ) ;
      // ======================================================================
    public:
      // ======================================================================
      inline bool  edges       () const { return m_edges       ; }
      inline bool  extrapolate () const { return m_extrapolate ; }
      inline bool  density     () const { return m_density     ; }
      // ======================================================================
    protected : 
      // ======================================================================
      /// special treatment of edges 
      bool                                  m_edges       { true  } ;
      /// extrapolate?
      bool                                  m_extrapolate { false } ;
      /// density                           
      bool                                  m_density     { false } ;
      // ======================================================================
    };
    // ========================================================================
    /** @class Histo1D HistoInterpolators.h Ostap/HistoInterpolators.h
     *  Simple function object/histogram interpolator 
     *  @see Ostap::Math::HistoInterpolation
     *  @see Ostap::Math::HistoInterpolation::interpolate_1D
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2019-11-16
     */
    class Histo1D : public HistoInterpolator 
    {
    public:
      // ======================================================================
      /** constructor with full  specification 
       *  @see Ostap::Math::HistoInterpolation
       *  @see Ostap::Math::HistoInterpolation::interpolate_1D
       */
      Histo1D 
      ( const TH1& histo ,
        const Ostap::Math::HistoInterpolation::Type t = Ostap::Math::HistoInterpolation::Default ,
        const bool edges       = true  , 
        const bool extrapolate = false , 
        const bool density     = false ) ;
      // ======================================================================
      /// constructor from the histogram and predefined configuration 
      Histo1D 
      ( const TH1&     histo ,
        const Histo1D& conf  ) ;
      // ======================================================================
      /// default  constructor 
      Histo1D () ;
      // ======================================================================
    public:
      // ======================================================================
      double xmin () const ;
      double xmax () const ;      
      // ======================================================================
    public:
      // ======================================================================
      inline double operator ()
      ( const double x ) const 
      {
        return Ostap::Math::HistoInterpolation::interpolate_1D
          ( m_h , x , m_t , edges () , extrapolate () , density () );                    
      }
      // ======================================================================
    public: // integrals 
      // ======================================================================
      /// integral over whole histogram range 
      double integral ()    const ;
      /// integral between low and high
      double integral 
      ( const double low  , 
        const double high ) const ;
      // ======================================================================
    public:
      // ======================================================================
      /** get a random number from this distribution
       *  - if maximum value is non-positive: generate uniform distribution
       *  - negative content is interpreted as zero
       *  @param rng  random generator
       *  @see TH1::GetRandom 
       *  @attention It can be rather inefficieint (e.g. for historgams 
       *  with large  number of empty vins 
       */
      double random
      ( TRandom* rng = nullptr ) const ;
      /// get a random number 
      std::size_t random
      ( double&  result           ,
        TRandom* rng    = nullptr ) const ;
      // ======================================================================
    public:
      // ======================================================================
      const TH1D&                           h () const  { return m_h ; }
      Ostap::Math::HistoInterpolation::Type t () const  { return m_t ; }      
      // ======================================================================
    public:
      // ======================================================================      
      /// get unique tag 
      std::size_t                           tag () const { return m_tag ; }
      // ======================================================================
    private :
      // ======================================================================
      // the histogram itself 
      TH1D                                  m_h {}  ; // the histogram itself 
      // ======================================================================
      /// interpolation type 
      Ostap::Math::HistoInterpolation::Type m_t { Ostap::Math::HistoInterpolation::Default };
      // ======================================================================
      // unique tag 
      std::size_t                           m_tag       { 0 } ;
      /// integration workspace 
      Ostap::Math::WorkSpace                m_workspace {   } ;
      // ======================================================================
    };
    // ========================================================================
    /** @class Histo2D HistoInterpolators.h Ostap/HistoInterpolators.h
     *  Simple function object/histogram interpolator 
     *  @see Ostap::Math::HistoInterpolation
     *  @see Ostap::Math::HistoInterpolation::interpolate_2D
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2019-11-16
     */
    class Histo2D : public HistoInterpolator 
    {
    public:
      // ======================================================================
      /** constructor with full  specification 
       *  @see Ostap::Math::HistoInterpolation
       *  @see Ostap::Math::HistoInterpolation::interpolate_1D
       */
      Histo2D 
      ( const TH2& histo , 
        const Ostap::Math::HistoInterpolation::Type tx = 
        Ostap::Math::HistoInterpolation::Default ,
        const Ostap::Math::HistoInterpolation::Type ty = 
        Ostap::Math::HistoInterpolation::Default ,
        const bool edges       = true  , 
        const bool extrapolate = false , 
        const bool density     = false ) ;
      // ======================================================================
      /// default  constructor 
      Histo2D () ;
      // ======================================================================
    public:
      // ======================================================================
      double xmin () const ;
      double xmax () const ;      
      double ymin () const ;
      double ymax () const ;      
      // ======================================================================
    public:
      // ======================================================================
      inline double operator ()
      ( const double x , 
        const double y ) const 
      {
        return Ostap::Math::HistoInterpolation::interpolate_2D
          ( m_h , x , y , m_tx , m_ty , 
            edges () , extrapolate () , density () );                    
      }
      // ======================================================================
    public:
      // ======================================================================
      /** get a random number from this distribution
       *  - if maximum value is non-positive: generate uniform distribution
       *  - negative content is interpreted as zero
       *  @param rng  random generator
       *  @see TH1::GetRandom 
       *  @attention It can be rather inefficieint (e.g. for historgams 
       *  with large  number of empty vins 
       */
      std::array<double,2> random
      ( TRandom* rnd = nullptr ) const ;
      /// get a random number 
      std::size_t random
      ( double&  xr            ,
        double&  yr            ,
        TRandom* rng = nullptr ) const ;
      // ======================================================================
    public:
      // ======================================================================
      const TH2D&                           h  () const  { return m_h  ; }
      Ostap::Math::HistoInterpolation::Type tx () const  { return m_tx ; }      
      Ostap::Math::HistoInterpolation::Type ty () const  { return m_ty ; }      
      // ======================================================================
    public: 
      // ======================================================================
      /// integral over whole histogram range 
      double integral ()    const ;
      /// integral between xmin, xmax, ymin, ymax 
      double integral 
      ( const double xmin  , 
        const double xmax  ,
        const double ymin  , 
        const double ymax  ) const ;
      double integrateX
      ( const double y     ) const ;
      double integrateX 
      ( const double y     , 
        const double xmin  , 
        const double xmax  ) const ;
      double integrateY 
      ( const double x     ) const ;
      double integrateY 
      ( const double x     , 
        const double ymin  , 
        const double ymax  ) const ;
      // ======================================================================
    public:
      // ======================================================================      
      /// get unique tag 
      std::size_t                           tag () const { return m_tag ; }
      // ======================================================================
    private :
      // ======================================================================
      // the histogram itself 
      TH2D                                  m_h {}  ;// the histogram itself 
      // ======================================================================
      /// interpolation type 
      Ostap::Math::HistoInterpolation::Type m_tx { Ostap::Math::HistoInterpolation::Default };
      /// interpolation type 
      Ostap::Math::HistoInterpolation::Type m_ty { Ostap::Math::HistoInterpolation::Default };
      // ======================================================================
      // unique tag 
      std::size_t                           m_tag       { 0 } ;
      /// integration workspace 
      Ostap::Math::WorkSpace                m_workspace {   } ;
      // ======================================================================
    };
    // ========================================================================
    /** @class Histo3D HistoInterpolators.h Ostap/HistoInterpolators.h
     *  Simple function object/histogram interpolator 
     *  @see Ostap::Math::HistoInterpolation
     *  @see Ostap::Math::HistoInterpolation::interpolate_3D
     *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
     *  @date   2019-11-16
     */
    class Histo3D : public HistoInterpolator 
    {
    public:
      // ======================================================================
      /** constructor with full  specification 
       *  @see Ostap::Math::HistoInterpolation
       *  @see Ostap::Math::HistoInterpolation::interpolate_3D
       */
      Histo3D 
      ( const TH3& histo , 
        const Ostap::Math::HistoInterpolation::Type tx = 
        Ostap::Math::HistoInterpolation::Default ,
        const Ostap::Math::HistoInterpolation::Type ty = 
        Ostap::Math::HistoInterpolation::Default ,
        const Ostap::Math::HistoInterpolation::Type tz = 
        Ostap::Math::HistoInterpolation::Default ,
        const bool edges       = true  , 
        const bool extrapolate = false , 
        const bool density     = false ) ;
      // ======================================================================
      /// default  constructor 
      Histo3D () ;
      // ======================================================================
    public:
      // ======================================================================
      double xmin () const ;
      double xmax () const ;      
      double ymin () const ;
      double ymax () const ;      
      double zmin () const ;
      double zmax () const ;      
      // ======================================================================
    public:
      // ======================================================================
      inline double operator () 
      ( const double x , 
        const double y ,
        const double z ) const 
      {
        return Ostap::Math::HistoInterpolation::interpolate_3D
          ( m_h ,  x  ,  y  ,  z  , m_tx  ,  m_ty  ,  m_tz , 
            edges () , extrapolate () , density () );                    
      }
      // ======================================================================
    public:
      // ======================================================================
      /** get a random number from this distribution
       *  - if maximum value is non-positive: generate uniform distribution
       *  - negative content is interpreted as zero
       *  @param rng  random generator
       *  @see TH1::GetRandom 
       *  @attention It can be rather inefficieint (e.g. for historgams 
       *  with large  number of empty vins 
       */
      std::array<double,3> random
      ( TRandom* rnd = nullptr ) const ;
      /// get a random number 
      std::size_t random
      ( double&  xr            ,
        double&  yr            ,
        double&  zr            ,
        TRandom* rng = nullptr ) const ;
      // ======================================================================
    public:
      // ======================================================================
      const TH3D&                           h  () const  { return m_h  ; }
      Ostap::Math::HistoInterpolation::Type tx () const  { return m_tx ; }      
      Ostap::Math::HistoInterpolation::Type ty () const  { return m_ty ; }      
      Ostap::Math::HistoInterpolation::Type tz () const  { return m_tz ; }      
      // ======================================================================
    public: // integrals 
      // ======================================================================
      /// integral over whole histogram range 
      double integral ()    const ;
      /// integral between xmin, xmax, ymin, ymax 
      double integral 
      ( const double xmin  , 
        const double xmax  ,
        const double ymin  , 
        const double ymax  ,
        const double zmin  , 
        const double zmax  ) const ;
      // ======================================================================
      double integrateXY
      ( const double z     ) const ;
      double integrateXY 
      ( const double z     , 
        const double xmin  , 
        const double xmax  ,
        const double ymin  , 
        const double ymax  ) const ;
      // ======================================================================
      double integrateXZ
      ( const double y     ) const ;
      double integrateXZ 
      ( const double y     , 
        const double xmin  , 
        const double xmax  ,
        const double zmin  , 
        const double zmax  ) const ;
      // ======================================================================
      double integrateYZ
      ( const double x     ) const ;
      double integrateYZ 
      ( const double y     , 
        const double ymin  , 
        const double ymax  ,
        const double zmin  , 
        const double zmax  ) const ;
      // ======================================================================
      double integrateX 
      ( const double y     , 
        const double z     ) const ;
      double integrateX 
      ( const double y     , 
        const double z     , 
        const double xmin  , 
        const double xmax  ) const ;
      // ======================================================================
      double integrateY 
      ( const double x     , 
        const double z     ) const ;
      double integrateY 
      ( const double x     , 
        const double z     , 
        const double ymin  , 
        const double ymax  ) const ;
      // ======================================================================
      double integrateZ 
      ( const double x     , 
        const double y     ) const ;
      double integrateZ 
      ( const double x     , 
        const double y     , 
        const double zmin  , 
        const double zmax  ) const ;
      // ======================================================================
    public:
      // ======================================================================      
      /// get unique tag 
      std::size_t                           tag () const { return m_tag ; }
      // ======================================================================
    private :
      // ======================================================================
      // the histogram itself 
      TH3D                                  m_h {}  ;// the histogram itself 
      // ======================================================================
      /// interpolation type 
      Ostap::Math::HistoInterpolation::Type m_tx { Ostap::Math::HistoInterpolation::Default };
      /// interpolation type 
      Ostap::Math::HistoInterpolation::Type m_ty { Ostap::Math::HistoInterpolation::Default };
      /// interpolation type 
      Ostap::Math::HistoInterpolation::Type m_tz { Ostap::Math::HistoInterpolation::Default };
      // ======================================================================
      // unique tag 
      std::size_t                           m_tag       { 0 } ;
      /// integration workspace 
      Ostap::Math::WorkSpace                m_workspace {   } ;
      // ======================================================================
    };
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
//                                                                     The END 
// ============================================================================
#endif // OSTAP_HISTOINTERPOLATORS_H
// ============================================================================
