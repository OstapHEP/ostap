// ============================================================================
#ifndef OSTAP_FITRESULT_H 
#define OSTAP_FITRESULT_H 1
// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <vector>
#include <string>
#include <limits>
// ============================================================================
// ROOT 
// ============================================================================
#include "RooAbsReal.h"
#include "RooFitResult.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Utils 
  {
    // ========================================================================
    /** @class  FitResults  Ostap/FitResult.h
     *  small extension of the Class RooFitResult 
     *  @see RooFitResult
     *  @author Vanya Belyaev
     *  @date   2019-03-27
     */
    class FitResults : public RooFitResult
    {
    public:
      // ======================================================================
      ClassDefOverride(Ostap::Utils::FitResults, 1) ;
      // ======================================================================
    public:
      // ======================================================================
      typedef std::vector<std::pair<std::string,int> > History ;
      // ======================================================================
    public:
      // ======================================================================
      /// Constructor from  RooFitResult 
      FitResults
      ( const RooFitResult& right              , 
        const char*         new_name = nullptr ) ;
      // ======================================================================
      /// full constructor #1
      FitResults
      ( const std::string&         name      ,
        const std::string&         title     ,
        const RooArgList&          constvars , // setConstParList 
        const RooArgList&          initvars  , // setInitParLits 
        const RooArgList&          finalvars , // setFinalParList
        const int                  status    , // setStatus 
        const int                  covqual   , // setCovQual 
        const double               minnll    , // setMinNLL     
        const double               edm       , // setEDM 
        const int                  numinvnll , // setNumInvalidNLL
        //
        const TMatrixDSym&         v         , // setCovarianceMatrix
        const History&             history   = History () ) ; // setStatusHistory 
      // ======================================================================
      /// full constructor #2
      FitResults 
      ( const std::string&         name      ,
        const std::string&         title     ,
        const RooArgList&          constvars , // setConstParList 
        const RooArgList&          initvars  , // setInitParLits 
        const RooArgList&          finalvars , // setFinalParList 
        const int                  status    , // setStatus 
        const int                  covqual   , // setCovQual 
        const double               minnll    , // setMinNLL     
        const double               edm       , // setEDM 
        const int                  numinvnll , // setNumInvalidNLL
        //
        const std::vector<double>& globalcc  , // fillCorrMatrix 
        const TMatrixDSym&         corrs     , // fillCorrMatrix 
        const TMatrixDSym&         covs      , // fillCorrMatrix
        const History&             history   = History () ) ; // setStatusHistory 
      /// copy constructor
      FitResults ( const FitResults& right ) ;
      /// default constructor 
      FitResults() = default ;  
      /// destructor 
      virtual ~FitResults () ;
      // clones 
      FitResults* Clone ( const char* newname = nullptr ) const override ;
      FitResults* clone () const override ;
      // ======================================================================
    public:
      // ======================================================================
      // check the status of FitResut 
      bool ok () const 
      {
        if ( !_finalPars || _finalPars->getSize() <= 0 ) { return false ; }
        const std::size_t size = _finalPars->getSize()  ;
        for ( std::size_t i = 0 ; i < size ; ++i )
          {
            const RooAbsArg*  a = _finalPars->at ( i ) ;
            if ( nullptr == a ) { return false ; }             // RETURN
            const RooAbsReal* r = dynamic_cast<const RooAbsReal*>( a ) ; //
            if ( nullptr == r ) { return false ; }             // RETURN            
          }
        return true ;
      }
      // ======================================================================
    public:
      // ======================================================================
      /// get all values 
      void get_values ( std::vector<double>& out ) const
      {
        // bad value: quite NaN 
        static const double s_NaN { std::numeric_limits<double>::quiet_NaN() } ;
        //
        Int_t size = _finalPars ? _finalPars->getSize () : 0 ;
        if ( size <= 0 ) { out.clear() ; return ; }                 // RETURN
        //
        out.resize ( size ) ;
        for ( Int_t i = 0 ; i < size ; ++i )
          {
            const RooAbsArg*  a = _finalPars->at ( i ) ;
            if ( nullptr ==  a ) { out [ i ] = s_NaN ; continue ;  } 
            const RooAbsReal* r = static_cast<const RooAbsReal*>( a ) ;
            out [ i ] = r ? r->getVal () : s_NaN ;
          }
      }
      // ======================================================================
      /// get vector of global correlations coeffficients 
      std::vector<double> global_cc () const ;
      // ======================================================================
      /// add label/status pair to the history 
      void add_to_history 
      ( const std::string& label  , 
        const int          status ) ;
      // ======================================================================        
    }; //                                The end of the class Ostap::FitResults
    // ========================================================================
    /** calculate the global correlation coefficient 
     * \f$ \rho_k = \sqrt{    1 - \left[ C_{kk} V_{kk}\right]^{-1} } \f$
     *  where \f$ C \f$ is covarinace matrix and \f$ V = C^{-1}\$ is inverse
     *  @code
     *  @param index parameter index 
     *  @param r     fit result 
     *  @return global correlation coefficients (or -1 in case of failure) 
     */
    double global_cc 
    ( const RooFitResult&  r     , 
      const unsigned short index ) ;
    // ========================================================================
    /** calculate the global correlation coefficients 
     * \f$ \rho_k = \sqrt{    1 - \left[ C_{kk} V_{kk}\right]^{-1} } \f$
     *  where \f$ C \f$ is covarinace matrix and \f$ V = C^{-1}\$ is inverse
     *  @code
     *  @param r     fit result 
     *  @return vector global correlation coefficients (empty in case of failure)
     */
    std::vector<double> global_cc 
    ( const RooFitResult& r ) ;
    // ========================================================================
  } //                                        The end of namespace Ostap::Utils
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                     The END 
// ============================================================================
#endif // OSTAP_FITRESULT_H
// ============================================================================
