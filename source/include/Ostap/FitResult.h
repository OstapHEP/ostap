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
// ============================================================================
// ROOT 
// ============================================================================
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
      ( const RooFitResult& right , 
        const char* new_name = nullptr ) ;
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
  } //                                        The end of namespace Ostap::Utils
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                     The END 
// ============================================================================
#endif // OSTAP_FITRESULT_H
// ============================================================================
