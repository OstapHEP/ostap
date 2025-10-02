// Include files
// ============================================================================
// STD&STL
// ============================================================================
#include <functional>
#include <algorithm>
// ============================================================================
// ROOT 
// ============================================================================
#include "TH1.h"
#include "TGraph.h"
#include "TAxis.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/Hash.h"
#include "Ostap/HistoHash.h"
#include "Ostap/HistoInterpolators.h"
#include "Ostap/StatusCode.h"
// ============================================================================
// Local 
// ============================================================================
#include "local_hash.h"
#include "local_math.h"
#include "status_codes.h"
// ============================================================================
/** @file 
 *  implementation file for function from the file Ostap/HistoHash.h 
 *  @date  2021-04-09 
 *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
 */
// ============================================================================
/* get hash for the given graph 
 *  @param graph the graphs 
 *  @return hash value 
 */
// ============================================================================
std::size_t Ostap::Utils::hash_graph 
( const TGraph* graph )
{
  if ( nullptr == graph ) { return 0 ; }
  const Int_t N = graph->GetN() ;
  //
  static const std::string s_what { "TGraph" } ;
  std::size_t seed =Ostap::Utils::hash_combiner ( s_what , N ) ;
  //
  //
  const Double_t*  X    = graph->GetX       () ;
  if ( X    ) { seed = Ostap::Utils::hash_combiner 
      ( seed , Ostap::Utils::hash_range ( X    , X    + N ) ) ; }
  //
  const Double_t*  Y    = graph->GetY       () ;
  if ( Y    ) { seed = Ostap::Utils::hash_combiner 
      ( seed , Ostap::Utils::hash_range ( Y    , Y    + N ) ) ; }
  //
  const Double_t*  EX   = graph->GetEX     () ;
  if ( EX   ) { seed = Ostap::Utils::hash_combiner 
      ( seed , Ostap::Utils::hash_range ( EX   , EX   + N ) ) ; }
  //
  const Double_t*  EY   = graph->GetEY      () ;
  if ( EY   ) { seed = Ostap::Utils::hash_combiner 
      ( seed , Ostap::Utils::hash_range ( EY   , EY   + N ) ) ; }
  //
  const Double_t*  EXh  = graph->GetEXhigh  () ;
  if ( EXh  ) { seed = Ostap::Utils::hash_combiner 
      ( seed , Ostap::Utils::hash_range ( EXh  , EXh  + N ) ) ; }
  //
  const Double_t*  EXl  = graph->GetEXlow   () ;
  if ( EXl  ) { seed = Ostap::Utils::hash_combiner 
      ( seed , Ostap::Utils::hash_range ( EXl  , EXl  + N ) ) ; }
  //  
  const Double_t*  EYh  = graph->GetEYhigh  () ;
  if ( EYh  ) { seed = Ostap::Utils::hash_combiner 
      ( seed , Ostap::Utils::hash_range ( EYh  , EYh  + N ) ) ; }
  //
  const Double_t*  EYl  = graph->GetEYlow   () ;
  if ( EYl  ) { seed = Ostap::Utils::hash_combiner
      ( seed , Ostap::Utils::hash_range ( EYl  , EYl  + N ) ) ; }
  //
  const Double_t*  EXhd = graph->GetEXhighd () ;
  if ( EXhd ) { seed = Ostap::Utils::hash_combiner 
      ( seed , Ostap::Utils::hash_range ( EXhd , EXhd + N ) ) ; }
  //
  const Double_t*  EXld = graph->GetEXlowd  () ;
  if ( EXld ) { seed = Ostap::Utils::hash_combiner 
      ( seed , Ostap::Utils::hash_range ( EXld , EXld + N ) ) ; }
  //  
  const Double_t*  EYhd = graph->GetEYhighd () ;
  if ( EYhd ) { seed = Ostap::Utils::hash_combiner 
      ( seed , Ostap::Utils::hash_range ( EYhd , EYhd + N ) ) ; }
  //
  const Double_t*  EYld = graph->GetEYlowd  () ;
  if ( EYld ) { seed = Ostap::Utils::hash_combiner 
      ( seed , Ostap::Utils::hash_range ( EYld , EYld + N ) ) ; }
  //
  return seed ;
}
// ============================================================================
/*  get hash for the given axis 
 *  @param axis the axis 
 *  @return hash value 
 */
// ============================================================================
std::size_t Ostap::Utils::hash_axis 
( const TAxis* axis )  
{
  if ( nullptr == axis ) { return 0 ; }
  //
  static const std::string s_what { "TAxis" } ;
  std::size_t seed = Ostap::Utils::hash_combiner 
    ( s_what              ,
      axis -> GetNbins () , 
      axis -> GetXmin  () , 
      axis -> GetXmax  () ) ;
  //
  const TArrayD* A = axis->GetXbins() ;
  if ( A && A->GetSize() ) 
  {
    const Double_t* a = A->GetArray() ;
    if ( a ) {  seed = Ostap::Utils::hash_combiner ( seed , Ostap::Utils::hash_range ( a , a + A ->GetSize () ) ) ; }
  }
  //
  return seed ;
}
// ============================================================================
/*  get hash for the given histogram 
 *  @param histo the historgam  
 *  @return hash value 
 */
// ============================================================================
std::size_t Ostap::Utils::hash_histo 
( const TH1* histo )  
{
  if ( nullptr == histo ) { return 0 ; }
  //
  static const std::string s_what { "TH1" } ;
  //
  const Int_t NX = histo -> GetNbinsX () ;
  const Int_t NY = histo -> GetNbinsY () ;
  const Int_t NZ = histo -> GetNbinsZ () ;
  //
  std::size_t seed = Ostap::Utils::hash_combiner ( s_what                          ,
						   histo -> GetDimension ()        , 
                                                   histo -> GetEntries   ()        , 
                                                   NX , NY , NZ                    ,
                                                   hash_axis ( histo -> GetXaxis() ) , 
                                                   hash_axis ( histo -> GetYaxis() ) , 
                                                   hash_axis ( histo -> GetZaxis() ) ) ;
  //
  for ( int ix = 1 ; ix <= NX ; ++ix ) 
  { for ( int iy = 1 ; iy <= NY ; ++iy ) 
    { for ( int iz = 1 ; iz <= NZ ; ++iz ) 
      { seed = Ostap::Utils::hash_combiner 
          ( seed , 
            histo -> GetBinContent ( ix , iy , iz ) ,
            histo -> GetBinError   ( ix , iy , iz ) ) ;  
      }
    }
  }
  //
  // add sumw2 information 
  const TArrayD* sumw2 = histo->GetSumw2() ;
  if ( sumw2 && sumw2->GetSize() ) 
  {
    const Double_t* a = sumw2->GetArray() ;
    if ( a ) { seed = Ostap::Utils::hash_combiner ( seed , Ostap::Utils::hash_range ( a , a + sumw2 ->GetSize () ) ) ; }
  }
  //
  // // add statistics: 
  // double stats [7];
  // histo -> GetStats ( stats ) ;
  // seed = Ostap::Utils::hash_combiner ( seed , Ostap::Utils::hash_range ( stats , stats + 7 ) ) ;
  //
  return seed ;
}
// ============================================================================
/*  get hash for Ostap::Math::Histo1D  object 
 *  @see Ostap::Math::Histo1D 
 */
// ============================================================================
std::size_t Ostap::Utils::hash_histo 
( const Ostap::Math::Histo1D& h ) 
{ return Ostap::Utils::hash_combiner 
    ( hash_histo ( &h.h() ) , 
      h.t           () ,
      h.edges       () , 
      h.extrapolate () , 
      h.density     () ) ; }
// ============================================================================
/*  get hash for Ostap::Math::Histo2D  object 
 *  @see Ostap::Math::Histo2D 
 */
// ============================================================================
std::size_t Ostap::Utils::hash_histo 
( const Ostap::Math::Histo2D& h ) 
{ return Ostap::Utils::hash_combiner 
    ( hash_histo ( &h.h() ) , 
      h.tx          () ,
      h.ty          () ,
      h.edges       () , 
      h.extrapolate () , 
      h.density     () ) ; }
// ============================================================================
/*  get hash for Ostap::Math::Histo3D  object 
 *  @see Ostap::Math::Histo3D 
 */
// ============================================================================
std::size_t Ostap::Utils::hash_histo 
( const Ostap::Math::Histo3D& h ) 
{ return Ostap::Utils::hash_combiner 
    ( hash_histo ( &h.h() ) , 
      h.tx          () ,
      h.ty          () ,
      h.tz          () ,
      h.edges       () , 
      h.extrapolate () , 
      h.density     () ) ; }
// ============================================================================
/* same binning ?
 *  @param a the first axis 
 *  @param b the second axis 
 *  @return True if axes have the sam ebinnigs , False otherwise
 */
// ============================================================================
bool Ostap::Utils::same_binning
( const TAxis& a ,
  const TAxis& b )
{
  if ( &a == &b     ) { return true ; }                              // RETURN
  //
  const Int_t a_nb = a.GetNbins() ;
  const Int_t b_nb = b.GetNbins() ;
  if ( a_nb != b_nb ) { return false ; }                             // RETURN
  //
  const double amin = a.GetXmin () ;
  const double bmin = b.GetXmin () ;
  if ( !s_equal ( amin , bmin )     ) { return false ; }             // RETURN
  //
  const double amax = a.GetXmax () ;
  const double bmax = b.GetXmax () ;
  if ( !s_equal ( amax , bmax )     ) { return false ; }             // RETURN 
  //
  const TArrayD* abins = a.GetXbins () ;
  const TArrayD* bbins = b.GetXbins () ;
  //
  if ( abins && bbins )
    {
      //
      Ostap::Assert ( abins->GetSize() == bbins->GetSize() ,
		      "Inconsistent TAxis bins!"   ,
		      "Ostap::Utils::same_binning" , 
		      INVALID_TAXISBINS , __FILE__ , __LINE__ ) ;
      //
      return std::equal ( abins->GetArray() ,
			  abins->GetArray() + abins->GetSize() ,
			  bbins->GetArray() , s_equal ) ;               // RETURN   
    }
  else if ( abins )
    {
      Ostap::Assert ( abins->GetSize() == b_nb + 1 ,
		      "Inconsistent TAxis bins!"   ,
		      "Ostap::Utils::same_binning" , 
		      INVALID_TAXISBINS , __FILE__ , __LINE__ ) ;
      //            
      const unsigned long NA = abins->GetSize  () ;
      const Double_t*     xa = abins->GetArray () ;
      for ( Int_t         ia = 0 ; ia < NA ; ++ia , ++xa  )
	{
	  const double xb =  ( ia * bmax + ( b_nb - ia ) * bmin ) / b_nb ; 
	  if ( !s_equal ( *xa , xb ) ) { return false ; }            // RETURN 
	}
    }
  else if ( bbins )
    {
      Ostap::Assert ( bbins->GetSize() == a_nb + 1 ,
		      "Inconsistent TAxis bins!"   ,
		      "Ostap::Utils::same_binning" , 
		      INVALID_TAXISBINS , __FILE__ , __LINE__ ) ;
      //            
      const unsigned long NB = bbins->GetSize  () ;
      const Double_t*     xb = bbins->GetArray () ;
      for ( Int_t         ib = 0 ; ib < NB ; ++ib , ++xb  )
	{
	  const double xa =  ( ib * amax + ( a_nb - ib ) * amin ) / a_nb ; 
	  if ( !s_equal ( *xb , xa ) ) { return false ; }         // RETURN 
	}
    }
  //
  return true ;
}
// ============================================================================
/*  same binning ?
 *  @param  a the first histo 
 *  @param  b the second histo 
 *  @return True if histos have the same binnings , False otherwise
 */
// ============================================================================
bool Ostap::Utils::same_binning
( const TH1&    a ,
  const TH1&    b )
{
  //
  if ( &a == &b     ) { return true  ; }                              // RETURN
  //
  const Int_t adim = a.GetDimension () ;
  const Int_t bdim = b.GetDimension () ;
  if ( adim != bdim ) { return false ; }                             // RETURN 
  //
  if ( 3 <= adim && !same_binning ( *a.GetZaxis() , *b.GetZaxis() ) ) { return false ; }
  if ( 2 <= adim && !same_binning ( *a.GetYaxis() , *b.GetYaxis() ) ) { return false ; }
  if ( 1 <= adim && !same_binning ( *a.GetXaxis() , *b.GetXaxis() ) ) { return false ; }
  //
  return true ;
}
// ============================================================================
//                                                                      The END 
// ============================================================================
