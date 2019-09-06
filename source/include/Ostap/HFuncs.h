// ============================================================================
#ifndef OSTAP_HFUNCS_H 
#define OSTAP_HFUNCS_H 1
// ============================================================================
// Include files
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/IFuncs.h"
#include "Ostap/HistoInterpolation.h"
// ============================================================================
// ROOT
// ============================================================================
#include  "TObject.h"
#include  "TH1D.h"
#include  "TH2D.h"
#include  "TH3D.h"
// ============================================================================
namespace Ostap 
{
  // ==========================================================================
  /// from Ostap 
  class  Formula ;
  // ==========================================================================
  namespace Functions 
  {
    // ========================================================================
    /** @class FuncTH Ostap/HFuncs.h
     *  simple implementation of Tree-function based on TH
     */
    class FuncTH : public TObject , public Ostap::IFuncTree
    {
    public:
      // ======================================================================
      ClassDef(Ostap::Functions::FuncTH,1) ;
      // ======================================================================
    public :
      // ======================================================================
      /** constructor
       *  @param tree          (INPUT) the tree 
       *  @param edges         (INPUT) special tretament of edges?
       *  @param extrapolate   (INPUT) use extrapolation?
       *  @param density       (INPUT) use  density?
       */
      FuncTH ( const TTree*         tree          = nullptr ,
               const bool           edges         = true     ,
               const bool           extrapolate   = false    , 
               const bool           density       = false    );
      // ======================================================================
      /// copy constructor
      FuncTH ( const FuncTH& right ) ;
      // ======================================================================
      /// destructor 
      virtual ~FuncTH () ;
      // ======================================================================
    public :
      // ======================================================================
      /// get the histogram 
      virtual const TH1* histo() const = 0 ;
      // ======================================================================
    protected:
      // ======================================================================
      mutable const TTree* m_tree           { nullptr } ; //!
      /// special treatment of edges?
      bool                 m_edges          { true    } ; // special treatment for edges?
      /// extrapolate?
      bool                 m_extrapolate    { false   } ; // extrapolate ?
      /// density?
      bool                 m_density        { false   } ; // density ?
      // ======================================================================
    } ;
    // ========================================================================
    /** @class FuncTH1
     *  simple implementation of Tree-function based on TH1
     */
    class FuncTH1 : public FuncTH
    {
    public:
      // ======================================================================
      ClassDef(Ostap::Functions::FuncTH1,1) ;
      // ======================================================================
    public :
      // ======================================================================
      /** constructor from the histogram 
       *  @param histo         (INPUT) the historgam 
       *  @param xvar          (INPUT) the expression/variable 
       *  @param tree          (INPUT) the tree 
       *  @param tx            (INPUT) interpolation type 
       *  @param edges         (INPUT) special tretament of edges?
       *  @param extrapolate   (INPUT) use extrapolation?
       *  @param density       (INPUT) use  density?
       */
      FuncTH1 ( const TH1F&          histo                    , 
                const std::string&   xvar                     , 
                const TTree*         tree          =  nullptr ,
                const Ostap::Math::HistoInterpolation::Type 
                tx = Ostap::Math::HistoInterpolation::Cubic ,
                const bool           edges         = true     ,
                const bool           extrapolate   = false    , 
                const bool           density       = false    );
      // ======================================================================
      /** constructor from the histogram 
       *  @param histo         (INPUT) the historgam 
       *  @param xvar          (INPUT) the expression/variable 
       *  @param tree          (INPUT) the tree 
       *  @param tx            (INPUT) interpolation type 
       *  @param edges         (INPUT) special tretament of edges?
       *  @param extrapolate   (INPUT) use extrapolation?
       *  @param density       (INPUT) use  density?
       */
      FuncTH1 ( const TH1D&          histo                    , 
                const std::string&   xvar                     , 
                const TTree*         tree          =  nullptr ,
                const Ostap::Math::HistoInterpolation::Type 
                tx = Ostap::Math::HistoInterpolation::Cubic ,
                const bool           edges         = true     ,
                const bool           extrapolate   = false    , 
                const bool           density       = false    );
      // ======================================================================
      // copy constructor
      FuncTH1 ( const FuncTH1& right ) ;
      // ======================================================================
      /// default constructor, needed for serialization 
      FuncTH1 () = default ;
      /// destructor 
      virtual ~FuncTH1 () ;
      // ======================================================================
    public:
      // ======================================================================
      FuncTH1* Clone ( const char* newname = "" ) const override ;
      // ======================================================================
    protected : // private constructor without histogram 
      // ======================================================================
      /** constructor without histogram 
       *  @param xvar          (INPUT) the expression/variable 
       *  @param tree          (INPUT) the tree 
       *  @param tx            (INPUT) interpolation type 
       *  @param edges         (INPUT) special tretament of edges?
       *  @param extrapolate   (INPUT) use extrapolation?
       *  @param density       (INPUT) use  density?
       */
      FuncTH1 ( const std::string&   xvar                     , 
                const TTree*         tree          =  nullptr ,
                const Ostap::Math::HistoInterpolation::Type 
                tx = Ostap::Math::HistoInterpolation::Cubic ,
                const bool           edges         = true     ,
                const bool           extrapolate   = false    , 
                const bool           density       = false    );
      // ======================================================================
    public:
      // ======================================================================
      ///  evaluate the formula for  TTree
      double operator () ( const TTree* tree ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the histogram 
      const TH1D* histo() const override { return &m_h1 ;}
      // ======================================================================
    public:
      // ======================================================================
      Bool_t Notify () override  ; 
      // ======================================================================
    private: 
      // ======================================================================
      /// make x-var 
      bool make_xvar   () const ;
      // ======================================================================
    protected:
      // ======================================================================
      mutable std::unique_ptr<Ostap::Formula> m_xvar { nullptr } ; //!
      // ======================================================================
      /// the  expression for x-variable   
      std::string    m_xvar_exp       {       }  ; // the  expression for x-variable
      /// interpolation type 
      Ostap::Math::HistoInterpolation::Type 
      m_tx { Ostap::Math::HistoInterpolation::Cubic } ; // the interpolation type
      // ======================================================================
    private:
      // ======================================================================
      ///  the histogram 
      TH1D           m_h1             {       } ; // the histogram
      // ======================================================================
    } ;
    // ========================================================================
    /** @class FuncTH2
     *  simple implementation of Tree-function based on TH2
     */
    class FuncTH2 : public FuncTH
    {
    public:
      // ======================================================================
      ClassDef(Ostap::Functions::FuncTH2,1) ;
      // ======================================================================
    public :
      // ======================================================================
      /** constructor from the histogram 
       *  @param histo         (INPUT) the historgam 
       *  @param xvar          (INPUT) the expression/variable 
       *  @param yvar          (INPUT) the expression/variable 
       *  @param tree          (INPUT) the tree 
       *  @param tx            (INPUT) interpolation type 
       *  @param ty            (INPUT) interpolation type 
       *  @param edges         (INPUT) special tretament of edges?
       *  @param extrapolate   (INPUT) use extrapolation?
       *  @param density       (INPUT) use  density?
       */
      FuncTH2 ( const TH2F&          histo                     , 
                const std::string&   xvar                      , 
                const std::string&   yvar                      , 
                const TTree*         tree           =  nullptr ,
                const Ostap::Math::HistoInterpolation::Type 
                tx = Ostap::Math::HistoInterpolation::Cubic ,
                const Ostap::Math::HistoInterpolation::Type 
                ty = Ostap::Math::HistoInterpolation::Cubic , 
                const bool           edges          = true     ,
                const bool           extrapolate    = false    , 
                const bool           density        = false    );
      // ======================================================================
      /** constructor from the histogram 
       *  @param histo         (INPUT) the historgam 
       *  @param xvar          (INPUT) the expression/variable 
       *  @param yvar          (INPUT) the expression/variable 
       *  @param tree          (INPUT) the tree 
       *  @param tx            (INPUT) interpolation type 
       *  @param ty            (INPUT) interpolation type 
       *  @param edges         (INPUT) special tretament of edges?
       *  @param extrapolate   (INPUT) use extrapolation?
       *  @param density       (INPUT) use  density?
       */
      FuncTH2 ( const TH2D&          histo                    , 
                const std::string&   xvar                     , 
                const std::string&   yvar                     , 
                const TTree*         tree          =  nullptr ,
                const Ostap::Math::HistoInterpolation::Type 
                tx = Ostap::Math::HistoInterpolation::Cubic , 
                const Ostap::Math::HistoInterpolation::Type 
                ty = Ostap::Math::HistoInterpolation::Cubic , 
                const bool           edges         = true     ,
                const bool           extrapolate   = false    , 
                const bool           density       = false    );
      // ======================================================================
      // copy constructor
      FuncTH2 ( const FuncTH2& right ) ;
      // ======================================================================
      /// default constructor, needed for serialization 
      FuncTH2 () = default ;
      /// destructor 
      virtual ~FuncTH2 () ;
      // ======================================================================
    public:
      // ======================================================================
      FuncTH2* Clone ( const char* newname = "" ) const override ;
      // ======================================================================
    protected : // private constructor without histogram 
      // ======================================================================
      /** constructor without histogram 
       *  @param histo         (INPUT) the historgam 
       *  @param xvar          (INPUT) the expression/variable 
       *  @param yvar          (INPUT) the expression/variable 
       *  @param tree          (INPUT) the tree 
       *  @param tx            (INPUT) interpolation type 
       *  @param ty            (INPUT) interpolation type 
       *  @param edges         (INPUT) special tretament of edges?
       *  @param extrapolate   (INPUT) use extrapolation?
       *  @param density       (INPUT) use  density?
       */
      FuncTH2 ( const std::string&   xvar                     , 
                const std::string&   yvar                     , 
                const TTree*         tree          =  nullptr ,
                const Ostap::Math::HistoInterpolation::Type 
                tx = Ostap::Math::HistoInterpolation::Cubic ,
                const Ostap::Math::HistoInterpolation::Type 
                ty = Ostap::Math::HistoInterpolation::Cubic , 
                const bool           edges         = true     ,
                const bool           extrapolate   = false    , 
                const bool           density       = false    );
      // ======================================================================
    public:
      // ======================================================================
      ///  evaluate the formula for  TTree
      double operator() ( const TTree* tree ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the histogram 
      const TH2D* histo() const override { return &m_h2 ;}
      // ======================================================================
    public:
      // ======================================================================
      Bool_t Notify () override ;
      // ======================================================================
    protected:
      // ======================================================================
      /// make x-var 
      bool make_xvar   () const ;
      /// make y-var 
      bool make_yvar   () const ;
      // ======================================================================
    protected:
      // ======================================================================
      mutable std::unique_ptr<Ostap::Formula> m_xvar { nullptr } ; //!
      mutable std::unique_ptr<Ostap::Formula> m_yvar { nullptr } ; //!
      // ======================================================================
      /// the  expression for x-variable   
      std::string    m_xvar_exp       {       } ; // the  expression for x-variable
      /// the  expression for y-variable   
      std::string    m_yvar_exp       {       } ; // the  expression for y-variable
      /// interpolation type 
      Ostap::Math::HistoInterpolation::Type 
      m_tx { Ostap::Math::HistoInterpolation::Cubic } ; // the interpolation type
      /// interpolation type 
      Ostap::Math::HistoInterpolation::Type 
      m_ty { Ostap::Math::HistoInterpolation::Cubic } ; // the interpolation type
      // ======================================================================
    private:
      // ======================================================================
      ///  the histogram 
      TH2D           m_h2             {       } ; // the histogram
      // ======================================================================
    } ;
    // ========================================================================
    /** @class FuncTH3
     *  simple implementation of Tree-function based on TH3
     */
    class FuncTH3 : public FuncTH
    {
    public:
      // ======================================================================
      ClassDef(Ostap::Functions::FuncTH3,1) ;
      // ======================================================================
    public :
      // ======================================================================
      /** constructor from the histogram 
       *  @param histo         (INPUT) the historgam 
       *  @param xvar          (INPUT) the expression/variable 
       *  @param yvar          (INPUT) the expression/variable 
       *  @param zvar          (INPUT) the expression/variable 
       *  @param tree          (INPUT) the tree 
       *  @param tx            (INPUT) interpolation type 
       *  @param ty            (INPUT) interpolation type 
       *  @param tz            (INPUT) interpolation type 
       *  @param edges         (INPUT) special tretament of edges?
       *  @param extrapolate   (INPUT) use extrapolation?
       *  @param density       (INPUT) use  density?
       */
      FuncTH3 ( const TH3F&          histo                     , 
                const std::string&   xvar                      , 
                const std::string&   yvar                      , 
                const std::string&   zvar                      , 
                const TTree*         tree           =  nullptr ,
                const Ostap::Math::HistoInterpolation::Type 
                tx = Ostap::Math::HistoInterpolation::Cubic ,
                const Ostap::Math::HistoInterpolation::Type 
                ty = Ostap::Math::HistoInterpolation::Cubic , 
                const Ostap::Math::HistoInterpolation::Type 
                tz = Ostap::Math::HistoInterpolation::Cubic , 
                const bool           edges          = true     ,
                const bool           extrapolate    = false    , 
                const bool           density        = false    );
      // ======================================================================
      /** constructor from the histogram 
       *  @param histo         (INPUT) the historgam 
       *  @param xvar          (INPUT) the expression/variable 
       *  @param yvar          (INPUT) the expression/variable 
       *  @param zvar          (INPUT) the expression/variable 
       *  @param tree          (INPUT) the tree 
       *  @param tx            (INPUT) interpolation type 
       *  @param ty            (INPUT) interpolation type 
       *  @param tz            (INPUT) interpolation type 
       *  @param edges         (INPUT) special tretament of edges?
       *  @param extrapolate   (INPUT) use extrapolation?
       *  @param density       (INPUT) use  density?
       */
      FuncTH3 ( const TH3D&          histo                    , 
                const std::string&   xvar                     , 
                const std::string&   yvar                     , 
                const std::string&   zvar                     , 
                const TTree*         tree          =  nullptr ,
                const Ostap::Math::HistoInterpolation::Type 
                tx = Ostap::Math::HistoInterpolation::Cubic , 
                const Ostap::Math::HistoInterpolation::Type 
                ty = Ostap::Math::HistoInterpolation::Cubic , 
                const Ostap::Math::HistoInterpolation::Type 
                tz = Ostap::Math::HistoInterpolation::Cubic , 
                const bool           edges         = true     ,
                const bool           extrapolate   = false    , 
                const bool           density       = false    );
      // ======================================================================
      // copy contructor
      // ======================================================================
      FuncTH3 ( const FuncTH3& right ) ;
      // ======================================================================
      /// default constructor, needed for serialization 
      FuncTH3 () = default ;
      /// destructor 
      virtual ~FuncTH3 () ;
      // ======================================================================
    public:
      // ======================================================================
      FuncTH3* Clone ( const char* newname = "" ) const override ;
      // ======================================================================
    protected : // private constructor without histogram 
      // ======================================================================
      /** constructor without histogram 
       *  @param xvar          (INPUT) the expression/variable 
       *  @param yvar          (INPUT) the expression/variable 
       *  @param zvar          (INPUT) the expression/variable 
       *  @param tree          (INPUT) the tree 
       *  @param tx            (INPUT) interpolation type 
       *  @param ty            (INPUT) interpolation type 
       *  @param tz            (INPUT) interpolation type 
       *  @param edges         (INPUT) special tretament of edges?
       *  @param extrapolate   (INPUT) use extrapolation?
       *  @param density       (INPUT) use  density?
       */
      FuncTH3 ( const std::string&   xvar                     , 
                const std::string&   yvar                     , 
                const std::string&   zvar                     , 
                const TTree*         tree          =  nullptr ,
                const Ostap::Math::HistoInterpolation::Type 
                tx = Ostap::Math::HistoInterpolation::Cubic ,
                const Ostap::Math::HistoInterpolation::Type 
                ty = Ostap::Math::HistoInterpolation::Cubic , 
                const Ostap::Math::HistoInterpolation::Type 
                tz = Ostap::Math::HistoInterpolation::Cubic , 
                const bool           edges         = true     ,
                const bool           extrapolate   = false    , 
                const bool           density       = false    );
      // ======================================================================
    public:
      // ======================================================================
      ///  evaluate the formula for  TTree
      double operator() ( const TTree* tree ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the histogram 
      const TH3D* histo() const override { return &m_h3 ; }
      // ======================================================================
    public:
      // ======================================================================
      Bool_t Notify () override ; 
      // ======================================================================
    protected:
      // ======================================================================
      /// make x-var 
      bool make_xvar   () const ;
      /// make y-var 
      bool make_yvar   () const ;
      /// make z-var 
      bool make_zvar   () const ;
      // ======================================================================
    protected:
      // ======================================================================
      mutable std::unique_ptr<Ostap::Formula> m_xvar { nullptr } ; //!
      mutable std::unique_ptr<Ostap::Formula> m_yvar { nullptr } ; //!
      mutable std::unique_ptr<Ostap::Formula> m_zvar { nullptr } ; //!
      // ======================================================================
      /// the  expression for x-variable   
      std::string    m_xvar_exp       {       } ; // the  expression for x-variable
      /// the  expression for y-variable   
      std::string    m_yvar_exp       {       } ; // the  expression for y-variable
      /// the  expression for z-variable   
      std::string    m_zvar_exp       {       } ; // the  expression for z-variable
      /// interpolation type 
      Ostap::Math::HistoInterpolation::Type 
      m_tx { Ostap::Math::HistoInterpolation::Cubic } ; // the interpolation type
      /// interpolation type 
      Ostap::Math::HistoInterpolation::Type 
      m_ty { Ostap::Math::HistoInterpolation::Cubic } ; // the interpolation type
      /// interpolation type 
      Ostap::Math::HistoInterpolation::Type 
      m_tz { Ostap::Math::HistoInterpolation::Cubic } ; // the interpolation type
      // ======================================================================
    private:
      // ======================================================================
      ///  the histogram 
      TH3D           m_h3             {       } ; // the histogram
      // ======================================================================
    } ; 
    // ========================================================================
  } //                                    The end of namespace Ostap::Functions
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_HFUNCS_H
// ============================================================================
