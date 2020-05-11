// ============================================================================
#ifndef OSTAP_MOREROOFIT2_H 
#define OSTAP_MOREROOFIT2_H 1
// ============================================================================
// Include files
// ============================================================================
// ROOT/RooFit 
// ============================================================================
#include "RVersion.h"
#include "RooAbsReal.h"
#include "RooRealProxy.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/PhaseSpace.h"
// ============================================================================
namespace Ostap 
{
  // =========================================================================
  namespace  MoreRooFit
  {
    // =======================================================================
    /** @class M2Q
     *  get the momenum in 2-body system with mass m
     *  @see Ostap::Math::PhaseSpace2 
     *  @see Ostap::MoreRooFit::Q2M
     */
    class M2Q : public RooAbsReal
    {
      // ========================================================================
      ClassDef(Ostap::MoreRooFit::M2Q, 1 ) ;
      // ========================================================================
    public:
      // ========================================================================
      M2Q () = default ;
      /// constructor 
      M2Q ( const std::string& name       , 
            const std::string& title      , 
            RooAbsReal&        m          ,
            const double       m1         ,
            const double       m2         ) ;
      /// constructor 
      M2Q ( const std::string& name       , 
            const std::string& title      , 
            RooAbsReal&        m          ,
            const double       m1         )
        : M2Q ( name , title , m , m1 , m1 )
      {} ;
      /// constructor 
      M2Q ( RooAbsReal&        m          ,
            const double       m1         ,
            const double       m2         ,
            const std::string& name  = "" , 
            const std::string& title = "" ) 
        : M2Q ( name , title , m , m1 , m2 )
      {} ;
      /// copy 
      M2Q ( const M2Q&    right       , 
            const char*        newname = 0 ) ;
      /// destructor
      virtual ~M2Q () ;
      /// clone 
      M2Q* clone ( const char* newname ) const override ;
      // ======================================================================
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;    
      // ======================================================================
    public:
      // ======================================================================
      const Ostap::Math::PhaseSpace2& phasespace  ()  const { return m_ps ; }
      // ======================================================================
    protected: 
      // ======================================================================
      Ostap::Math::PhaseSpace2 m_ps {} ;
      // ======================================================================
      /// mass-variable
      RooRealProxy m_m {} ; // mass-variable
      // ======================================================================
    } ;   
    // ========================================================================
    /** @class Q2M
     *  get the mass from the momeutm in 2-body system 
     *  f[ m = \sqrt{m+1^2 + q^2} + \sqrt{m_2^2+q^2} \f]
     *  @see Ostap""MoreRooFit::M2Q
     */
    class Q2M : public M2Q 
    {
      // ========================================================================
      ClassDef(Ostap::MoreRooFit::Q2M, 1 ) ;
      // ========================================================================
    public:
      // ========================================================================
      Q2M () = default ;
      /// constructor
      Q2M ( const std::string& name       , 
            const std::string& title      , 
            RooAbsReal&        q          ,
            const double       m1         ,
            const double       m2         ) ;
      /// constructor 
      Q2M ( const std::string& name       , 
            const std::string& title      , 
            RooAbsReal&        q          ,
            const double       m1         )
        : Q2M ( name , title , q , m1 , m1 )
      {} ;
      /// constructor 
      Q2M ( RooAbsReal&        q          ,
            const double       m1         ,
            const double       m2         ,
            const std::string& name  = "" , 
            const std::string& title = ""  ) 
        : Q2M ( name , title , q , m1 , m2 )
      {} ;
      /// copy 
      Q2M ( const Q2M&    right       , 
            const char*        newname = 0 ) ;
      /// destructor
      virtual ~Q2M () ;
      /// clone 
      Q2M* clone ( const char* newname ) const override ;
      // ======================================================================
    protected:
      // ======================================================================
      // the actual evaluation of the result 
      Double_t evaluate () const override ;    
      // ======================================================================
    } ;   
    // ========================================================================    
  } //                                   The end of namepsace Ostap::MoreRooFit 
  // ==========================================================================
} //                                                 The end of namepsace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_MOREROOFIT2_H
// ============================================================================
