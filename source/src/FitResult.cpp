// ============================================================================
// Include files 
// ============================================================================
// STD&STL 
// ============================================================================
// ROOT 
// ============================================================================
#include "TVector.h"
#include "RooFitResult.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/FitResult.h"
// ============================================================================
/** @file 
 *  Implementation file for class Ostap::FitResult
 *  @see Ostap::FitResult 
 *  @date 2022-08-08 
 */
// ============================================================================
// Constructor from RooFitResult 
// ============================================================================
Ostap::Utils::FitResults::FitResults
( const RooFitResult& right   , 
  const char*         newname ) 
  : RooFitResult ( right )
{
  if ( nullptr != newname ) { SetName ( newname ) ; }
}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Utils::FitResults::FitResults
( const Ostap::Utils::FitResults& right ) 
  : RooFitResult ( right )
{}
// ============================================================================
// full constructor #1
// ============================================================================
Ostap::Utils::FitResults::FitResults 
( const std::string&                       name      ,
  const std::string&                       title     ,
  const RooArgList&                        constvars , // setConstParList 
  const RooArgList&                        initvars  , // setInitParLits 
  const RooArgList&                        finalvars , // setFinalParList 
  const int                                status    , // setStatus 
  const int                                covqual   , // setCovQual 
  const double                             minnll    , // setMinNLL     
  const double                             edm       , // setEDM 
  const int                                numinvnll , // setNumInvalidNLL
  const TMatrixDSym&                       v         , // setCovarianceMatrix 
  const Ostap::Utils::FitResults::History& history   ) // setStatusHistory 
: RooFitResult ( name.c_str() , title.c_str() ) 
{
  setConstParList     ( constvars ) ;
  setInitParList      ( initvars  ) ;
  setNumInvalidNLL    ( numinvnll ) ;
  setStatus           ( status    ) ;
  setCovQual          ( covqual   ) ;
  setMinNLL           ( minnll    ) ;
  setFinalParList     ( finalvars ) ;
  setCovarianceMatrix ( const_cast<TMatrixDSym&> ( v ) ) ;
  setStatusHistory    ( const_cast<History&> ( history ) ) ;
  //
  _GC = new TVectorD ( _CM->GetNcols() ) ;
  fillLegacyCorrMatrix() ;
  //
}
// ============================================================================
// full constructor #2
// ============================================================================
Ostap::Utils::FitResults::FitResults
( const std::string&                       name      ,
  const std::string&                       title     ,
  const RooArgList&                        constvars , // setConstParList 
  const RooArgList&                        initvars  , // setInitParLits 
  const RooArgList&                        finalvars , // setFinalParList 
  const int                                status    , // setStatus 
  const int                                covqual   , // setCovQual 
  const double                             minnll    , // setMinNLL     
  const double                             edm       , // setEDM 
  const int                                numinvnll , // setNumInvalidNLL
  const std::vector<double>&               globalcc  , // fillCorrMatrix 
  const TMatrixDSym&                       corrs     , // fillCorrMatrix 
  const TMatrixDSym&                       covs      , // fillCorrMatrix 
  const Ostap::Utils::FitResults::History& history   ) // setStatusHistory 
{
  setConstParList     ( constvars ) ;
  setInitParList      ( initvars  ) ;
  setNumInvalidNLL    ( numinvnll ) ;
  setStatus           ( status    ) ;
  setCovQual          ( covqual   ) ;
  setMinNLL           ( minnll    ) ;
  setFinalParList     ( finalvars ) ;
  fillCorrMatrix      ( globalcc , corrs , covs ) ;
  setStatusHistory    ( const_cast<History&> ( history ) ) ;
}
// ============================================================================
// destructor 
// ============================================================================
Ostap::Utils::FitResults::~FitResults (){}
// ============================================================================
// clone 
// ============================================================================
Ostap::Utils::FitResults* 
Ostap::Utils::FitResults::Clone ( const char* newname ) const
{ return new Ostap::Utils::FitResults ( *this , newname ) ; }
// ============================================================================
Ostap::Utils::FitResults* 
Ostap::Utils::FitResults::clone () const { return Clone() ; }
// ============================================================================
std::vector<double> 
Ostap::Utils::FitResults::global_cc () const 
{
  if ( nullptr == _GC ) { return std::vector<double>() ; }
  return std::vector<double>( _GC->GetMatrixArray() , 
                              _GC->GetMatrixArray() + 
                              _GC->GetNrows()        ) ;
}
// ============================================================================
// add label/status pair to the history 
// ============================================================================
void Ostap::Utils::FitResults::add_to_history 
( const std::string& label  , 
  const int          status )
{ _statusHistory.push_back ( std::make_pair ( label , status ) ) ; }
// ===========================================================================



// ============================================================================
ClassImp(Ostap::Utils::FitResults    )
// ============================================================================
//                                                          The END 
// ============================================================================
