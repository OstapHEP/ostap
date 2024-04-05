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
#include "Ostap/Math.h"
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
  //
  if ( _CM && !_GC ) { _GC = new TVectorD ( _CM->GetNcols() ) ; }
  fillLegacyCorrMatrix() ;
  //
}
// ============================================================================
// copy constructor 
// ============================================================================
Ostap::Utils::FitResults::FitResults
( const Ostap::Utils::FitResults& right ) 
  : RooFitResult ( right )
{
  if ( _CM && !_GC ) { _GC = new TVectorD ( _CM->GetNcols() ) ; }
  fillLegacyCorrMatrix() ;
}
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
  if ( _CM && !_GC ) { _GC = new TVectorD ( _CM->GetNcols() ) ; }
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
  // use (pre)calcualted values 
  if ( nullptr != _GC ) 
  { 
    return std::vector<double>( _GC->GetMatrixArray () , 
                                _GC->GetMatrixArray () + 
                                _GC->GetNrows       () ) ; 
  }
  // recalculate: 
  return Ostap::Utils::global_cc ( *this ) ;
}
// ============================================================================
// add label/status pair to the history 
// ============================================================================
void Ostap::Utils::FitResults::add_to_history 
( const std::string& label  , 
  const int          status )
{ _statusHistory.push_back ( std::make_pair ( label , status ) ) ; }
// ============================================================================
/*  calculate the global correlation coefficients 
 * \f$ \rho_k = \sqrt{    1 - \left[ C_{kk} V_{kk}\right]^{-1} } \f$
 *  where \f$ C \f$ is covarinace matrix and \f$ V = C^{-1}\$ is inverse
 *  @code
 *  @param r     fit result 
 *  @return vector global correlation coefficients (empty in case of failure)
 */
// ============================================================================
std::vector<double>
Ostap::Utils::global_cc ( const RooFitResult& r ) 
{
  /// zero for doubles  
  const static Ostap::Math::Zero<double> s_zero {}      ; // zero for doubles
  ///
  /// 1) try to invert covariance matrix 
  TMatrixTSym<double> cinv { r.covarianceMatrix() } ;
  double det = 0 ;
  cinv.Invert ( &det  ) ;
  if ( s_zero ( det ) ) { return std::vector<double>() ; }
  ///
  const int ncols = cinv.GetNcols() ;
  std::vector<double> result ( ncols, 0.0 ) ;
  //
  const  TMatrixTSym<double>& cm = r.covarianceMatrix() ;
  for ( int index = 0 ; index < ncols ; ++index ) 
  {
    const double cv = cinv ( index , index ) * cm ( index , index ) ;
    if ( s_zero ( cv ) ) { return std::vector<double>() ; }
    const double rho2 = 1.0 - 1.0 / cv ;
    if ( rho2 < 0      ) { return std::vector<double>() ; }
    result [ index ] = std::sqrt ( rho2 ) ;
  }
  //
  return result ;
}
// ============================================================================
/*  calculate the global correlation coefficient 
 * \f$ \rho_k = \sqrt{    1 - \left[ C_{kk} V_{kk}\right]^{-1} } \f$
 *  where \f$ C \f$ is covarinace matrix and \f$ V = C^{-1}\$ is inverse
 *  @code
 *  @param r     fit result 
 *  @return global correlation coefficients (or -1 in case of failure) 
 */
// ============================================================================
double Ostap::Utils::global_cc 
( const RooFitResult&  r     , 
  const unsigned short index ) 
{
  const std::vector<double> result { Ostap::Utils::global_cc ( r ) } ;
  if ( result.size() <= index ) { return -1 ; }
  return result [ index ] ;
}
// ============================================================================
ClassImp(Ostap::Utils::FitResults    )
// ============================================================================
//                                                                      The END 
// ============================================================================
