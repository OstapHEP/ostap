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
#include "Ostap/HistoInterpolators.h"
// ============================================================================
// ROOT
// ============================================================================
#include  "TObject.h"
#include  "TH1D.h"
#include  "TH2D.h"
#include  "TH3D.h"
// ============================================================================
// RooFit 
// ============================================================================
#include  "RooFormulaVar.h"
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
    /** @class FuncTH1
     *  simple implementation of Tree-function based on TH1
     */
    class FuncTH1 : public TObject , public Ostap::IFuncTree
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
      FuncTH1 ( const TH1&           histo                     , 
                const std::string&   xvar                      , 
                const TTree*         tree          =  nullptr  ,
                const Ostap::Math::HistoInterpolation::Type tx = 
                Ostap::Math::HistoInterpolation::Default       ,
                const bool           edges         = true      ,
                const bool           extrapolate   = false     , 
                const bool           density       = false     ) ;
      // ======================================================================
      /** constructor from the histogram 
       *  @param histo         (INPUT) the historgam 
       *  @param xvar          (INPUT) the expression/variable 
       *  @param tree          (INPUT) the tree 
       */
      FuncTH1 ( const Ostap::Math::Histo1D& histo            , 
                const std::string&          xvar             , 
                const TTree*                tree  =  nullptr ) ;
      // ======================================================================
      // copy constructor
      FuncTH1 ( const FuncTH1& right ) ;
      // ======================================================================
      /// default constructor, needed for serialization 
      FuncTH1 () = default ;
      // ======================================================================
    public:
      // ======================================================================
      FuncTH1* Clone ( const char* newname = "" ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      ///  evaluate the formula for  TTree
      double operator () ( const TTree* tree ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the histogram 
      const Ostap::Math::Histo1D& histo () const { return m_histo     ; }
      /// get the histogram 
      const TH1D&                 h     () const { return m_histo.h() ; }
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
    private :
      // ======================================================================
      /// the histogram object 
      Ostap::Math::Histo1D  m_histo    ; /// the histogram object 
      /// expression for x-axis 
      std::string           m_xvar_exp ; /// expression for x-axis 
      /// the actual function for x-axis 
      mutable std::unique_ptr<Ostap::Formula> m_xvar { nullptr } ; //!
      /// the tree itself 
      mutable const TTree*                    m_tree { nullptr } ; //!
      // ======================================================================
    } ;
    // ========================================================================
    /** @class FuncTH2
     *  simple implementation of Tree-function based on TH2
     */
    class FuncTH2 : public TObject , public Ostap::IFuncTree
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
      FuncTH2 ( const TH2&           histo                     , 
                const std::string&   xvar                      , 
                const std::string&   yvar                      , 
                const TTree*         tree           =  nullptr ,
                const Ostap::Math::HistoInterpolation::Type tx = 
                Ostap::Math::HistoInterpolation::Default       ,
                const Ostap::Math::HistoInterpolation::Type ty = 
                Ostap::Math::HistoInterpolation::Default       , 
                const bool           edges          = true     ,
                const bool           extrapolate    = false    , 
                const bool           density        = false    );
      // ======================================================================
      /** constructor from the histogram 
       *  @param histo         (INPUT) the historgam 
       *  @param xvar          (INPUT) the expression/variable 
       *  @param yvar          (INPUT) the expression/variable 
       *  @param tree          (INPUT) the tree 
       */
      FuncTH2 ( const Ostap::Math::Histo2D& histo            , 
                const std::string&          xvar             , 
                const std::string&          yvar             , 
                const TTree*                tree  =  nullptr ) ;
      // ======================================================================
      // copy constructor
      FuncTH2 ( const FuncTH2& right ) ;
      // ======================================================================
      /// default constructor, needed for serialization 
      FuncTH2 () = default ;
      // ======================================================================
    public:
      // ======================================================================
      FuncTH2* Clone ( const char* newname = "" ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      ///  evaluate the formula for  TTree
      double operator() ( const TTree* tree ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the histogram 
      const Ostap::Math::Histo2D& histo () const { return m_histo     ; }
      /// get the histogram 
      const TH2D&                 h     () const { return m_histo.h() ; }
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
    private: 
      // ======================================================================
      /// the histogram object 
      Ostap::Math::Histo2D  m_histo    ; /// the histogram object 
      /// expression for x-axis 
      std::string           m_xvar_exp ; /// expression for x-axis 
      /// expression for y-axis 
      std::string           m_yvar_exp ; /// expression for y-axis 
      /// the actual function for x-axis 
      mutable std::unique_ptr<Ostap::Formula> m_xvar { nullptr } ; //!
      /// the actual function for y-axis 
      mutable std::unique_ptr<Ostap::Formula> m_yvar { nullptr } ; //!
      /// the tree itself 
      mutable const TTree*                    m_tree { nullptr } ; //!
      // ======================================================================
    } ;
    // ========================================================================
    /** @class FuncTH3
     *  simple implementation of Tree-function based on TH3
     */
    class FuncTH3 : public TObject , public Ostap::IFuncTree
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
      FuncTH3 ( const TH3&           histo                     , 
                const std::string&   xvar                      , 
                const std::string&   yvar                      , 
                const std::string&   zvar                      , 
                const TTree*         tree           =  nullptr ,
                const Ostap::Math::HistoInterpolation::Type tx = 
                Ostap::Math::HistoInterpolation::Default       ,
                const Ostap::Math::HistoInterpolation::Type ty = 
                Ostap::Math::HistoInterpolation::Default       , 
                const Ostap::Math::HistoInterpolation::Type tz =
                Ostap::Math::HistoInterpolation::Default       , 
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
       */
      FuncTH3 ( const Ostap::Math::Histo3D& histo            , 
                const std::string&          xvar             , 
                const std::string&          yvar             , 
                const std::string&          zvar             , 
                const TTree*                tree  =  nullptr ) ;
      // ======================================================================
      // copy contructor
      // ======================================================================
      FuncTH3 ( const FuncTH3& right ) ;
      // ======================================================================
      /// default constructor, needed for serialization 
      FuncTH3 () = default ;
      // ======================================================================
    public:
      // ======================================================================
      FuncTH3* Clone ( const char* newname = "" ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      ///  evaluate the formula for  TTree
      double operator() ( const TTree* tree ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the histogram 
      const Ostap::Math::Histo3D& histo () const { return m_histo     ; }
      /// get the histogram 
      const TH3D&                 h     () const { return m_histo.h() ; }
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
    private: 
      // ======================================================================
      /// the histogram object 
      Ostap::Math::Histo3D  m_histo    ; /// the histogram object 
      /// expression for x-axis 
      std::string           m_xvar_exp ; /// expression for x-axis 
      /// expression for y-axis 
      std::string           m_yvar_exp ; /// expression for y-axis 
      /// expression for z-axis 
      std::string           m_zvar_exp ; /// expression for z-axis 
      /// the actual function for x-axis 
      mutable std::unique_ptr<Ostap::Formula> m_xvar { nullptr } ; //!
      /// the actual function for y-axis 
      mutable std::unique_ptr<Ostap::Formula> m_yvar { nullptr } ; //!
      /// the actual function for z-axis 
      mutable std::unique_ptr<Ostap::Formula> m_zvar { nullptr } ; //!
      /// the tree itself 
      mutable const TTree*                    m_tree { nullptr } ; //!
      // ======================================================================
    } ; 
    // ========================================================================
    // 'RooAbsData'-functions
    // ========================================================================
    /** @class FuncRooTH1
     *  simple implementation of 'RooAbsData'-function based on TH1
     */
    class FuncRooTH1 : public Ostap::IFuncData
    {
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
      FuncRooTH1 ( const TH1&           histo                     , 
                   const std::string&   xvar                      , 
                   const RooAbsData*    data       =  nullptr     ,
                   const Ostap::Math::HistoInterpolation::Type tx = 
                   Ostap::Math::HistoInterpolation::Default       ,
                   const bool           edges         = true      ,
                   const bool           extrapolate   = false     , 
                   const bool           density       = false     ) ;
      // ======================================================================
      /** constructor from the histogram 
       *  @param histo         (INPUT) the historgam 
       *  @param xvar          (INPUT) the expression/variable 
       *  @param tree          (INPUT) the tree 
       */
      FuncRooTH1 ( const Ostap::Math::Histo1D& histo            , 
                   const std::string&          xvar             , 
                   const RooAbsData*           data  =  nullptr ) ;
      // ======================================================================
      // copy constructor
      FuncRooTH1 ( const FuncRooTH1& right ) ;
      // ======================================================================
      /// default constructor, needed for serialization 
      FuncRooTH1 () = default ;
      // ======================================================================
    public:
      // ======================================================================
      ///  evaluate the formula for  data
      double operator () ( const RooAbsData* data ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the histogram 
      const Ostap::Math::Histo1D& histo () const { return m_histo     ; }
      /// get the histogram 
      const TH1D&                 h     () const { return m_histo.h() ; }
      // ======================================================================
    private: 
      // ======================================================================
      /// make x-var 
      bool make_xvar   () const ;
      // ======================================================================
    private :
      // ======================================================================
      /// the histogram object 
      Ostap::Math::Histo1D  m_histo    ; /// the histogram object 
      /// expression for x-axis 
      std::string           m_xvar_exp ; /// expression for x-axis 
      /// the actual function for x-axis 
      mutable std::unique_ptr<RooFormulaVar>  m_xvar { nullptr } ; //!
      /// the tree itself 
      mutable const RooAbsData*               m_data { nullptr } ; //!
      // ======================================================================
    } ;
    // ========================================================================
    /** @class FuncRooTH2
     *  simple implementation of 'RooAbsData'-function based on TH1
     */
    class FuncRooTH2 : public Ostap::IFuncData
    {
    public :
      // ======================================================================
      /** constructor from the histogram 
       *  @param histo         (INPUT) the historgam 
       *  @param xvar          (INPUT) the expression/variable 
       *  @param yvar          (INPUT) the expression/variable 
       *  @param data          (INPUT) data
       *  @param tx            (INPUT) interpolation type 
       *  @param ty            (INPUT) interpolation type 
       *  @param edges         (INPUT) special tretament of edges?
       *  @param extrapolate   (INPUT) use extrapolation?
       *  @param density       (INPUT) use  density?
       */
      FuncRooTH2 ( const TH2&           histo                     , 
                   const std::string&   xvar                      , 
                   const std::string&   yvar                      , 
                   const RooAbsData*    data       =  nullptr     ,
                   const Ostap::Math::HistoInterpolation::Type tx = 
                   Ostap::Math::HistoInterpolation::Default       ,
                   const Ostap::Math::HistoInterpolation::Type ty = 
                   Ostap::Math::HistoInterpolation::Default       ,
                   const bool           edges         = true      ,
                   const bool           extrapolate   = false     , 
                   const bool           density       = false     ) ;
      // ======================================================================
      /** constructor from the histogram 
       *  @param histo         (INPUT) the historgam 
       *  @param xvar          (INPUT) the expression/variable 
       *  @param yvar          (INPUT) the expression/variable 
       *  @param data          (INPUT) data
       */
      FuncRooTH2 ( const Ostap::Math::Histo2D& histo            , 
                   const std::string&          xvar             , 
                   const std::string&          yvar             , 
                   const RooAbsData*           data  =  nullptr ) ;
      // ======================================================================
      // copy constructor
      FuncRooTH2 ( const FuncRooTH2& right ) ;
      // ======================================================================
      /// default constructor, needed for serialization 
      FuncRooTH2 () = default ;
      // ======================================================================
    public:
      // ======================================================================
      ///  evaluate the formula for  data
      double operator () ( const RooAbsData* data ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the histogram 
      const Ostap::Math::Histo2D& histo () const { return m_histo     ; }
      /// get the histogram 
      const TH2D&                 h     () const { return m_histo.h() ; }
      // ======================================================================
    private: 
      // ======================================================================
      /// make x-var 
      bool make_xvar   () const ;
      /// make y-var 
      bool make_yvar   () const ;
      // ======================================================================
    private :
      // ======================================================================
      /// the histogram object 
      Ostap::Math::Histo2D  m_histo    ; /// the histogram object 
      /// expression for x-axis 
      std::string           m_xvar_exp ; /// expression for x-axis 
      /// expression for y-axis 
      std::string           m_yvar_exp ; /// expression for y-axis 
      /// the actual function for x-axis 
      mutable std::unique_ptr<RooFormulaVar>  m_xvar { nullptr } ; //!
      /// the actual function for y-axis 
      mutable std::unique_ptr<RooFormulaVar>  m_yvar { nullptr } ; //!
      /// the tree itself 
      mutable const RooAbsData*               m_data { nullptr } ; //!
      // ======================================================================
    } ;
    // ========================================================================
    /** @class FuncRooTH3
     *  simple implementation of 'RooAbsData'-function based on TH3
     */
    class FuncRooTH3 : public Ostap::IFuncData
    {
    public :
      // ======================================================================
      /** constructor from the histogram 
       *  @param histo         (INPUT) the historgam 
       *  @param xvar          (INPUT) the expression/variable 
       *  @param yvar          (INPUT) the expression/variable 
       *  @param zvar          (INPUT) the expression/variable 
       *  @param data          (INPUT) data
       *  @param tx            (INPUT) interpolation type 
       *  @param ty            (INPUT) interpolation type 
       *  @param tz            (INPUT) interpolation type 
       *  @param edges         (INPUT) special tretament of edges?
       *  @param extrapolate   (INPUT) use extrapolation?
       *  @param density       (INPUT) use  density?
       */
      FuncRooTH3 ( const TH3&           histo                     , 
                   const std::string&   xvar                      , 
                   const std::string&   yvar                      , 
                   const std::string&   zvar                      , 
                   const RooAbsData*    data       =  nullptr     ,
                   const Ostap::Math::HistoInterpolation::Type tx = 
                   Ostap::Math::HistoInterpolation::Default       ,
                   const Ostap::Math::HistoInterpolation::Type ty = 
                   Ostap::Math::HistoInterpolation::Default       ,
                   const Ostap::Math::HistoInterpolation::Type tz = 
                   Ostap::Math::HistoInterpolation::Default       ,
                   const bool           edges         = true      ,
                   const bool           extrapolate   = false     , 
                   const bool           density       = false     ) ;
      // ======================================================================
      /** constructor from the histogram 
       *  @param histo         (INPUT) the historgam 
       *  @param xvar          (INPUT) the expression/variable 
       *  @param yvar          (INPUT) the expression/variable 
       *  @param zvar          (INPUT) the expression/variable 
       *  @param data          (INPUT) data
       */
      FuncRooTH3 ( const Ostap::Math::Histo3D& histo            , 
                   const std::string&          xvar             , 
                   const std::string&          yvar             , 
                   const std::string&          zvar             , 
                   const RooAbsData*           data  =  nullptr ) ;
      // ======================================================================
      // copy constructor
      FuncRooTH3 ( const FuncRooTH3& right ) ;
      // ======================================================================
      /// default constructor, needed for serialization 
      FuncRooTH3 () = default ;
      // ======================================================================
    public:
      // ======================================================================
      ///  evaluate the formula for  data
      double operator () ( const RooAbsData* data ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the histogram 
      const Ostap::Math::Histo3D& histo () const { return m_histo     ; }
      /// get the histogram 
      const TH3D&                 h     () const { return m_histo.h() ; }
      // ======================================================================
    private: 
      // ======================================================================
      /// make x-var 
      bool make_xvar   () const ;
      /// make y-var 
      bool make_yvar   () const ;
      /// make y-var 
      bool make_zvar   () const ;
      // ======================================================================
    private :
      // ======================================================================
      /// the histogram object 
      Ostap::Math::Histo3D  m_histo    ; /// the histogram object 
      /// expression for x-axis 
      std::string           m_xvar_exp ; /// expression for x-axis 
      /// expression for y-axis 
      std::string           m_yvar_exp ; /// expression for y-axis 
      /// expression for z-axis 
      std::string           m_zvar_exp ; /// expression for z-axis 
      /// the actual function for x-axis 
      mutable std::unique_ptr<RooFormulaVar>  m_xvar { nullptr } ; //!
      /// the actual function for y-axis 
      mutable std::unique_ptr<RooFormulaVar>  m_yvar { nullptr } ; //!
      /// the actual function for z-axis 
      mutable std::unique_ptr<RooFormulaVar>  m_zvar { nullptr } ; //!
      /// the tree itself 
      mutable const RooAbsData*               m_data { nullptr } ; //!
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
