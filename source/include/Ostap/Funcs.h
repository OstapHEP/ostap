// ============================================================================
#ifndef OSTAP_FUNCS_H 
#define OSTAP_FUNCS_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <string>
#include <memory>
#include <functional>
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/IFuncs.h"
#include "Ostap/Formula.h"
#include "Ostap/TreeGetter.h"
#include "Ostap/RooFun.h"
// ============================================================================
// ROOT
// ============================================================================
#include  "TObject.h"
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
    /** @class FuncFormula Ostap/Funcs.h
     *  simple implementation of TTree-function based on Ostap::Formula
     */
    class FuncFormula : public Ostap::IFuncTree, public TObject
    {
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Functions::FuncFormula,2) ;
      // ======================================================================
    public :
      // ======================================================================
      /** constructor from the formula expression 
       *  @param expression the  formula expression 
       *  @param tree       the tree 
       *  @param name       the name for the formula 
       */
      FuncFormula 
      ( const std::string& expression            , 
        const TTree*       tree       =  nullptr ,
        const std::string& name       = ""       ) ;
      // ======================================================================
      /// copy constructor 
      FuncFormula ( const  FuncFormula& right )  ;
      // ======================================================================
      virtual ~FuncFormula() ;
      // ======================================================================
    public:
      // ======================================================================
      /// default constructor, needed for (de) serialization 
      FuncFormula () = default ;
      // ======================================================================
    public: // Clone&clone 
      // ======================================================================
      /// IFuncTree::clone 
      FuncFormula* clone ( const char* newname = "" ) const override ;
      /// TObject::Clone 
      FuncFormula* Clone ( const char* newname = "" ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      ///  evaluate the formula for  TTree
      double operator() ( const TTree* tree ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the expressoon
      const std::string& expression () const { return m_expression ; }
      /// function name 
      const std::string& fun_name   () const { return m_name ; }
      // ======================================================================
    public:
      // ======================================================================
      bool ok () const { return m_formula && m_formula->ok () ;  }
      // ======================================================================
    public:
      // ======================================================================
      Bool_t Notify   () override ; 
      // ======================================================================
    private:
      // ======================================================================
      /// make formula 
      bool make_formula() const ;
      // ======================================================================
    private:
      // ======================================================================
      mutable const TTree*                    m_tree    { nullptr } ; //! 
      mutable std::unique_ptr<Ostap::Formula> m_formula { nullptr } ; //!
      // ======================================================================
      /// the  expression itself 
      std::string m_expression {} ; // the  expression itself 
      /// the name  
      std::string m_name       {} ; // the name  
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Func1D 
     *  Generic 1D-function 
     */
    class Func1D : public TObject , public Ostap::IFuncTree 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Functions::Func1D,2) ;
      // ======================================================================
    public :
      // ======================================================================
      template <class FUNCTION>
      Func1D 
      ( FUNCTION           fun             , 
        const std::string& x               ,
        const TTree*       tree =  nullptr ) 
        : TObject () 
        , m_fun      ( fun     )
        , m_xvar_exp ( x       ) 
        , m_xvar     { nullptr }
        , m_tree     { tree    }
      {}     
      // ======================================================================
      Func1D 
      ( std::function<double(double)> fun  , 
        const std::string& x               ,
        const TTree*       tree =  nullptr ) 
        : TObject () 
        , m_fun      ( fun     )
        , m_xvar_exp ( x       ) 
        , m_xvar     { nullptr }
        , m_tree     { tree    }
      {}     
      // ======================================================================
      /// copy constructor
      Func1D ( const Func1D& right ) ;
      // ======================================================================
      /// default constructor, needed for serialization 
      Func1D () = default ;
      // ======================================================================
      /// virtual destructor 
      virtual ~Func1D() ;
      // ======================================================================
    public: // Clone&clone 
      // ======================================================================
      /// IFuncTree::clone 
      Func1D* clone ( const char* newname = "" ) const override ;
      /// TObject::Clone 
      Func1D* Clone ( const char* newname = "" ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      template <class FUNCTION>
      static inline Func1D 
      create 
      ( FUNCTION           fun             , 
        const std::string& x               ,
        const TTree*       tree =  nullptr ) 
      { return Func1D ( fun , x , tree ) ; } 
      // ======================================================================
    public:
      // ======================================================================
      ///  evaluate the function for TTree
      double operator () ( const TTree* tree ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      const std::string& x() const { return m_xvar_exp ; }
      // ======================================================================
    public:
      // ======================================================================
      Bool_t Notify   () override ; 
      // ======================================================================
    private:
      // ======================================================================
      /// make formula 
      bool make_xvar () const ;
      // ======================================================================
    public:
      // ======================================================================
      //  evaluate the function 
      double func ( const double x ) const { return m_fun  ( x  ) ; }
      // ======================================================================
    protected :
      // ======================================================================
      /// the function  itself 
      std::function<double(double)> m_fun      {} ; /// the function 
      /// expression for x-axis 
      std::string                   m_xvar_exp {} ; /// expression for x-axis 
      /// the actual function for x-axis 
      mutable std::unique_ptr<Ostap::Formula> m_xvar { nullptr } ; //!
      /// the tree itself 
      mutable const TTree*                    m_tree { nullptr } ; //!
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Func2D 
     *  Generic 2D-function 
     */
    class Func2D : public TObject , public Ostap::IFuncTree 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Functions::Func2D,2) ;
      // ======================================================================
    public :
      // ======================================================================
      template <class FUNCTION>
      Func2D 
      ( FUNCTION           fun             , 
        const std::string& x               ,
        const std::string& y               ,
        const TTree*       tree =  nullptr ) 
        : TObject () 
        , m_fun      ( fun     )
        , m_xvar_exp ( x       ) 
        , m_yvar_exp ( y       ) 
        , m_xvar     { nullptr }
        , m_yvar     { nullptr }
        , m_tree     { tree    }
      {}     
      // ======================================================================
      /// copy constructor
      Func2D ( const Func2D& right ) ;
      // ======================================================================
      /// default constructor, needed for serialization 
      Func2D () = default ;
      // ======================================================================
      /// virtual destructor 
      virtual ~Func2D() ;
      // ======================================================================
    public: // Clone&clone 
      // ======================================================================
      /// IFuncTree::clone 
      Func2D* clone ( const char* newname = "" ) const override ;
      /// TObject::Clone 
      Func2D* Clone ( const char* newname = "" ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      template <class FUNCTION>
      static inline Func2D 
      create
      ( FUNCTION           fun             , 
        const std::string& x               ,
        const std::string& y               ,
        const TTree*       tree =  nullptr ) 
      { return Func2D ( fun , x , y ,  tree ) ; } 
      // ======================================================================
    public:
      // ======================================================================
      ///  evaluate the function for TTree
      double operator () ( const TTree* tree ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      const std::string& x () const { return m_xvar_exp ; }
      const std::string& y () const { return m_yvar_exp ; }
      // ======================================================================
    public:
      // ======================================================================
      Bool_t Notify   () override ; 
      // ======================================================================
    private:
      // ======================================================================
      /// make xvar
      bool make_xvar () const ;
      /// make yvar
      bool make_yvar () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the function 
      double func 
      ( const double x , 
        const double y ) const { return m_fun  ( x , y ) ; }
      // ======================================================================
    protected :
      // ======================================================================
      /// the function  itself 
      std::function<double(double,double)> m_fun {} ; /// the function 
      /// expression for x-axis 
      std::string                   m_xvar_exp {} ; /// expression for x-axis 
      /// expression for y-axis 
      std::string                   m_yvar_exp {} ; /// expression for y-axis 
      /// the actual function for x-axis 
      mutable std::unique_ptr<Ostap::Formula> m_xvar { nullptr } ; //!
      /// the actual function for y-axis 
      mutable std::unique_ptr<Ostap::Formula> m_yvar { nullptr } ; //!
      /// the tree itself 
      mutable const TTree*                    m_tree { nullptr } ; //!
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Func3D 
     *  Generic 3D-function 
     */
    class Func3D : public TObject , public Ostap::IFuncTree 
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Functions::Func3D,2) ;
      // ======================================================================
    public :
      // ======================================================================
      template <class FUNCTION>
      Func3D 
      ( FUNCTION           fun             , 
        const std::string& x               ,
        const std::string& y               ,
        const std::string& z               ,
        const TTree*       tree =  nullptr ) 
        : TObject () 
        , m_fun      ( fun     )
        , m_xvar_exp ( x       ) 
        , m_yvar_exp ( y       ) 
        , m_zvar_exp ( z       ) 
        , m_xvar     { nullptr }
        , m_yvar     { nullptr }
        , m_zvar     { nullptr }
        , m_tree     { tree    }
      {}     
      // ======================================================================
      /// copy constructor
      Func3D ( const Func3D& right ) ;
      // ======================================================================
      /// default constructor, needed for serialization 
      Func3D () = default ;
      // ======================================================================
      /// virtual destructor 
      virtual ~Func3D() ;
      // ======================================================================
    public: // Clone&clone 
      // ======================================================================
      /// IFuncTree::clone 
      Func3D* clone ( const char* newname = "" ) const override ;
      /// TObject::Clone 
      Func3D* Clone ( const char* newname = "" ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      template <class FUNCTION>
      static inline Func3D 
      create
      ( FUNCTION           fun             , 
        const std::string& x               ,
        const std::string& y               ,
        const std::string& z               ,
        const TTree*       tree =  nullptr ) 
      { return Func3D ( fun , x , y , z , tree ) ; } 
      // ======================================================================
    public:
      // ======================================================================
      ///  evaluate the function for TTree
      double operator () ( const TTree* tree ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      const std::string& x () const { return m_xvar_exp ; }
      const std::string& y () const { return m_yvar_exp ; }
      const std::string& z () const { return m_zvar_exp ; }
      // ======================================================================
    public:
      // ======================================================================
      Bool_t Notify   () override ; 
      // ======================================================================
    private:
      // ======================================================================
      /// make xvar
      bool make_xvar () const ;
      /// make yvar
      bool make_yvar () const ;
      /// make yvar
      bool make_zvar () const ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the function 
      double func
      ( const double x , 
        const double y , 
        const double z ) const { return m_fun  ( x , y , z ) ; }
      // ======================================================================
    protected :
      // ======================================================================
      /// the function  itself 
      std::function<double(double,double,double)> m_fun {} ; /// the function 
      /// expression for x-axis 
      std::string                   m_xvar_exp {} ; /// expression for x-axis 
      /// expression for y-axis 
      std::string                   m_yvar_exp {} ; /// expression for y-axis 
      /// expression for z-axis 
      std::string                   m_zvar_exp {} ; /// expression for z-axis 
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
    /** @class RooTreeFun 
     *  Specil "tree-function" that evalated RooFit-function 
     */
    class RooTreeFun : public Ostap::IFuncTree, public Ostap::Trees::RooGetter 
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Functions::RooTreeFun,1) ;
      // ======================================================================
    public:
      // ======================================================================
      /** full constructor 
       *  @param fun the functon 
       *  @param observables observables 
       *  @param normalization nornalization 
       *  @Param mapping  RooFit varibale <-> TTree branch mapping 
       *  @param tree input tree 
       */
      RooTreeFun
      ( const RooAbsReal&       fun           , 
        const RooAbsCollection& observables   , 
        const RooAbsCollection* normalization = nullptr ,
        const DCT&              mapping       = DCT()   ,
        const TTree*            tree          = nullptr ) ;
      // ======================================================================
      /** full constructor 
       *  @param fun the functon 
       *  @param observables observables 
       *  @Param mapping  RooFit varibale <-> TTree branch mapping 
       *  @param tree input tree 
       */
      RooTreeFun
      ( const RooAbsReal&       fun                     , 
        const RooAbsCollection& observables             , 
        const DCT&              mapping                 ,       
        const TTree*            tree          = nullptr ) ;
      // ======================================================================
      /** full constructor 
       *  @param fun the functon 
       *  @param observables observables 
       *  @param normalization nornalization 
       *  @Param mapping  RooFit varibale <-> TTree branch mapping 
       *  @param tree input tree 
       */
      RooTreeFun
      ( const RooAbsReal&       fun           , 
        const RooAbsData&       observables   , 
        const RooAbsCollection* normalization = nullptr ,
        const DCT&              mapping       = DCT()   ,
        const TTree*            tree          = nullptr ) ;
      // ======================================================================
      /** full constructor 
       *  @param fun the functon 
       *  @param observables observables 
       *  @Param mapping  RooFit varibale <-> TTree branch mapping 
       *  @param tree input tree 
       */
      RooTreeFun
      ( const RooAbsReal&       fun                     , 
        const RooAbsData&       observables             , 
        const DCT&              mapping                 ,       
        const TTree*            tree          = nullptr ) ;      
      // ======================================================================
      /// copy constructor 
      RooTreeFun  ( const RooTreeFun&  right ) ;
      /// move constructor 
      RooTreeFun  (       RooTreeFun&& right ) = default ;
      /// virtual constructor 
      virtual ~RooTreeFun() ;
      // ======================================================================
    public: // Clone&clone 
      // ======================================================================
      // IFunTree::clone 
      RooTreeFun* Clone ( const char* newname = "" ) const override ;
      // Tobject::clone 
      RooTreeFun* clone ( const char* newname = "" ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      ///  evaluate the formula for TTree
      double operator() ( const TTree* tree ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the fnuiction
      const RooAbsReal& function      () const { return m_fun.fun           () ; }
      /// get the observables 
      const RooArgSet&  observables   () const { return m_fun.observables   () ; }
      /// get the parameters 
      const RooArgSet&  parameters    () const { return m_fun.parameters    () ; }
      /// get the normalization 
      const RooArgSet*  normalization () const { return m_fun.normalization () ; }
      // ======================================================================
    private:
      // ======================================================================
      /// RooFun object
      Ostap::MoreRooFit::RooFun   m_fun {} ; //Roofun object
      // ======================================================================
    } ;
    // ========================================================================
    /** @class FuncRooFormula
     *  simple implementation of 'RooAbsData'-function based on RooFormulaVar
     */
    class FuncRooFormula : public Ostap::IFuncData
    {
    public :
      // ======================================================================
      /** constructor from the formula expression 
       *  @param expression the formula expression 
       *  @param data       the data
       *  @param name       the name for the formula 
       */
      FuncRooFormula 
      ( const std::string& expression            , 
        const RooAbsData*  data       =  nullptr ,
        const std::string& name       = ""       ) ;
      // ======================================================================
      /// copy contructor
      FuncRooFormula ( const FuncRooFormula&  right ) ;
      // ======================================================================
      /// default constructor, needed for serialization 
      FuncRooFormula () = default ;
      // ======================================================================
      /// destructor 
      virtual ~FuncRooFormula() ;
      // ======================================================================
    public:
      // ======================================================================
      // IFuncData::clone 
      FuncRooFormula* clone ( const char* name = "" ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// expression
      const std::string& expression () const { return m_expression ; }
      /// function name 
      const std::string& fun_name   () const { return m_name       ; }
      // ======================================================================
    public:
      // ======================================================================
      ///  evaluate the formula for  data
      double operator () ( const RooAbsData* data ) const override ;
      // ======================================================================
    private:
      // ======================================================================
      /// make formula 
      bool make_formula() const ;
      // ======================================================================
    private:
      // ======================================================================
      mutable const RooAbsData*              m_data    { nullptr } ;
      mutable std::unique_ptr<RooFormulaVar> m_formula { nullptr } ;
      // ======================================================================
      /// the  expression itself 
      std::string m_expression {} ; // the  expression itself 
      /// the name  
      std::string m_name       {} ; // the name  
      // ======================================================================
    } ;
    // ========================================================================
    /** @class FuncRoo1D 
     *  Generic 1D-function 
     */
    class FuncRoo1D : public Ostap::IFuncData
    {
    public :
      // ======================================================================
      template <class FUNCTION>
      FuncRoo1D 
      ( FUNCTION           fun              , 
        const std::string& x                ,
        const RooAbsData*  data  =  nullptr ) 
        : Ostap::IFuncData () 
        , m_fun      ( fun     )
        , m_xvar_exp ( x       ) 
        , m_xvar     { nullptr }
        , m_data     { data    }
      {}     
      // ======================================================================
      /// copy constructor
      FuncRoo1D ( const FuncRoo1D& right ) ;
      // ======================================================================
      /// default constructor, needed for serialization 
      FuncRoo1D () = default ;
      // ======================================================================
      /// virtual destructor
      virtual ~FuncRoo1D() ; 
      // ======================================================================
    public:
      // ======================================================================
      template <class FUNCTION>
      static inline FuncRoo1D 
      create 
      ( FUNCTION           fun             , 
        const std::string& x               ,
        const RooAbsData*  data =  nullptr ) 
      { return FuncRoo1D ( fun , x , data ) ; } 
      // ======================================================================
    public:
      // ======================================================================
      // IFuncData::clone 
      FuncRoo1D* clone ( const char* name = "" ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      ///  evaluate the function for TTree
      double operator () ( const RooAbsData* tree ) const override ;
      // ======================================================================
    private:
      // ======================================================================
      /// make formula 
      bool make_xvar () const ;
      // ======================================================================
    public:
      // ======================================================================
      const std::string& x () const { return m_xvar_exp ; }
      // ======================================================================
    public:
      // ======================================================================
      //  evaluate the function 
      double func ( const double x ) const { return m_fun  ( x  ) ; }
      // ======================================================================
    protected :
      // ======================================================================
      /// the function  itself 
      std::function<double(double)> m_fun      {} ; /// the function 
      /// expression for x-axis 
      std::string                   m_xvar_exp {} ; /// expression for x-axis 
      /// the actual function for x-axis 
      mutable std::unique_ptr<RooFormulaVar> m_xvar { nullptr } ; //!
      /// the tree itself 
      mutable const RooAbsData*              m_data { nullptr } ; //!
      // ======================================================================
    } ;
    // ========================================================================
    /** @class FuncRoo2D 
     *  Generic 2D-function 
     */
    class FuncRoo2D : public Ostap::IFuncData
    {
    public :
      // ======================================================================
      template <class FUNCTION>
      FuncRoo2D
      ( FUNCTION           fun              , 
        const std::string& x                ,
        const std::string& y                ,
        const RooAbsData*  data  =  nullptr ) 
        : Ostap::IFuncData () 
        , m_fun      ( fun     )
        , m_xvar_exp ( x       ) 
        , m_yvar_exp ( y       ) 
        , m_xvar     { nullptr }
        , m_yvar     { nullptr }
        , m_data     { data    }
      {}     
      // ======================================================================
      /// copy constructor
      FuncRoo2D ( const FuncRoo2D& right ) ;
      // ======================================================================
      /// default constructor, needed for serialization 
      FuncRoo2D () = default ;
      // ======================================================================
      /// virtual destructor
      virtual ~FuncRoo2D() ; 
      // ======================================================================
    public:
      // ======================================================================
      // IFuncData::clone 
      FuncRoo2D* clone ( const char* name = "" ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      template <class FUNCTION>
      static inline FuncRoo2D 
      create
      ( FUNCTION           fun             , 
        const std::string& x               ,
        const std::string& y               ,
        const RooAbsData*  data =  nullptr ) 
      { return FuncRoo2D ( fun , x , y , data ) ; } 
      // ======================================================================
    public:
      // ======================================================================
      ///  evaluate the function for TTree
      double operator () ( const RooAbsData* tree ) const override ;
      // ======================================================================
    private:
      // ======================================================================
      /// make formula 
      bool make_xvar () const ;
      /// make formula 
      bool make_yvar () const ;
      // ======================================================================
    public:
      // ======================================================================
      const std::string& x () const { return m_xvar_exp ; }
      const std::string& y () const { return m_yvar_exp ; }
      // ======================================================================
    public:
      // ======================================================================
      //  evaluate the function 
      double func 
      ( const double x , 
        const double y ) const { return m_fun  ( x , y) ; }
      // ======================================================================
    protected :
      // ======================================================================
      /// the function  itself 
      std::function<double(double,double)> m_fun  {} ; /// the function 
      /// expression for x-axis 
      std::string                   m_xvar_exp {} ; /// expression for x-axis 
      /// expression for y-axis 
      std::string                   m_yvar_exp {} ; /// expression for x-axis 
      /// the actual function for x-axis 
      mutable std::unique_ptr<RooFormulaVar> m_xvar { nullptr } ; //!
      /// the actual function for y-axis 
      mutable std::unique_ptr<RooFormulaVar> m_yvar { nullptr } ; //!
      /// the tree itself 
      mutable const RooAbsData*              m_data { nullptr } ; //!
      // ======================================================================
    } ;
    // ========================================================================
    /** @class FuncRoo3D 
     *  Generic 3D-function 
     */
    class FuncRoo3D : public Ostap::IFuncData
    {
    public :
      // ======================================================================
      template <class FUNCTION>
      FuncRoo3D
      ( FUNCTION           fun              , 
        const std::string& x                ,
        const std::string& y                ,
        const std::string& z                ,
        const RooAbsData*  data  =  nullptr ) 
        : Ostap::IFuncData () 
        , m_fun      ( fun     )
        , m_xvar_exp ( x       ) 
        , m_yvar_exp ( y       ) 
        , m_zvar_exp ( z       ) 
        , m_xvar     { nullptr }
        , m_yvar     { nullptr }
        , m_zvar     { nullptr }
        , m_data     { data    }
      {}     
      // ======================================================================
      /// copy constructor
      FuncRoo3D ( const FuncRoo3D& right ) ;
      // ======================================================================
      /// default constructor, needed for serialization 
      FuncRoo3D () = default ;
      // ======================================================================
      /// virtual destructor
      virtual ~FuncRoo3D() ; 
      // ======================================================================
    public:
      // ======================================================================
      template <class FUNCTION>
      static inline FuncRoo3D 
      create
      ( FUNCTION           fun             , 
        const std::string& x               ,
        const std::string& y               ,
        const std::string& z               ,
        const RooAbsData*  data =  nullptr ) 
      { return FuncRoo3D ( fun , x , y , z , data ) ; } 
      // ======================================================================
    public:
      // ======================================================================
      // IFuncData::clone 
      FuncRoo3D* clone ( const char* name = "" ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      ///  evaluate the function for TTree
      double operator () ( const RooAbsData* tree ) const override ;
      // ======================================================================
    private:
      // ======================================================================
      /// make formula 
      bool make_xvar () const ;
      /// make formula 
      bool make_yvar () const ;
      /// make formula 
      bool make_zvar () const ;
      // ======================================================================
    public:
      // ======================================================================
      const std::string& x () const { return m_xvar_exp ; }
      const std::string& y () const { return m_yvar_exp ; }
      const std::string& z () const { return m_zvar_exp ; }
      // ======================================================================
    public:
      // ======================================================================
      //  evaluate the function 
      double func 
      ( const double x , 
        const double y , 
        const double z ) const { return m_fun ( x , y , z ) ; }
      // ======================================================================
    protected :
      // ======================================================================
      /// the function  itself 
      std::function<double(double,double,double)> m_fun  {} ; /// the function 
      /// expression for x-axis 
      std::string                   m_xvar_exp {} ; /// expression for x-axis 
      /// expression for y-axis 
      std::string                   m_yvar_exp {} ; /// expression for x-axis 
      /// expression for z-axis 
      std::string                   m_zvar_exp {} ; /// expression for z-axis 
      /// the actual function for x-axis 
      mutable std::unique_ptr<RooFormulaVar> m_xvar { nullptr } ; //!
      /// the actual function for y-axis 
      mutable std::unique_ptr<RooFormulaVar> m_yvar { nullptr } ; //!
      /// the actual function for z-axis 
      mutable std::unique_ptr<RooFormulaVar> m_zvar { nullptr } ; //!
      /// the tree itself 
      mutable const RooAbsData*              m_data { nullptr } ; //!
      // ======================================================================
    } ;
    // ========================================================================
    /** @class Expression
     *  "Universal" formula that is simultaneously 
     *  - Ostap::IFuncTree 
     *  - Ostap::IFuncData
     *  @see Ostap::Functions::FuncFormula
     *  @see Ostap::Functions::FuncRooFormula
     *  @see Ostap::IFuncTree 
     *  @see Ostap::IFuncData
     *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
     *  @date 2023-03-03
     */
    class Expression: public FuncFormula, public Ostap::IFuncData 
    {
      // ======================================================================
    public :
      // ======================================================================
      ClassDefOverride(Ostap::Functions::Expression,1) ;
      // ======================================================================
    public :
      // ======================================================================
      /** constructor from the formula expression 
       *  @param expression the  formula expression 
       *  @param tree       the tree 
       *  @param name       the name for the formula 
       */
      Expression
      ( const std::string& expression            , 
        const TTree*       tree       =  nullptr ,
        const std::string& name       = ""       ) ;
      // ======================================================================
      /** constructor from the formula expression 
       *  @param expression the formula expression 
       *  @param data       the data 
       *  @param name       the name for the formula 
       */
      Expression
      ( const std::string& expression            , 
        const RooAbsData*  data                  ,
        const std::string& name       = ""       ) ;
      // ======================================================================
      /// copy constructor 
      Expression ( const Expression& right )  ;
      // ======================================================================
      /// default constructor, needed for serialization 
      Expression () = default ;
      /// destructor 
      virtual ~Expression() ;
      // ======================================================================
    public:
      // ======================================================================
    public: // Clone&clone 
      // ======================================================================
      /// IFuncTree::clone, IFuncData::clone 
      Expression* clone ( const char* newname = "" ) const override ;
      /// TObject::Clone 
      Expression* Clone ( const char* newname = "" ) const override ;
      // ======================================================================
    public:
      // ======================================================================
      /// evaluate the function from TTree 
      double operator ()  ( const TTree*      tree ) const override ;
      /// evaluate the function from RooAbsData 
      double operator ()  ( const RooAbsData* data ) const override ;
      // ======================================================================
    private:      
      // ======================================================================
      /// Ostap::IFuncData  for Roo-stuff 
      FuncRooFormula m_roofun{} ; // IFuncData  for Roo-stuff 
      // ======================================================================
    } ;
    // ========================================================================    
  } //                                   The END of  namespace Ostap::Functions
  // ==========================================================================
} //                                                The END of  namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_FUNCS_H
// ============================================================================

