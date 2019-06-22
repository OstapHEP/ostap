// ============================================================================
// Include files
// ============================================================================
// STD&STL
// ============================================================================
// ============================================================================
// ROOT&RooFit 
// ============================================================================
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "RooAbsArg.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooDataSet.h"
#include "RooFormulaVar.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/IFuncs.h"
#include "Ostap/AddVars.h"
// ============================================================================
/** @file
 *  Implementation fiel for functions from file Ostap/AddVars.h
 *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 *  @date 2019-06-22
 */
// ============================================================================
/*  add new variable to dataset
 *  @param  dataset input    dataset
 *  @param  name    variable name 
 *  @param  func    rule to  calculate new variable
 *  @return the added variable 
 */
// ============================================================================
const RooAbsReal* 
Ostap::Functions::add_var 
( RooDataSet&             dataset , 
  const std::string&      name    , 
  const Ostap::IFuncData& func    ) 
{  
  //
  RooRealVar var         { name.c_str() , ""        , 0.0 } ;
  RooArgSet  varset      { var          ,"one var" } ;
  RooDataSet new_dataset { "" , ""      , varset   } ;
  //
  // loop over events in the input data set 
  const unsigned long nEntries = dataset.numEntries() ;
  for ( unsigned long entry = 0 ; entry < nEntries ; ++entry )   
  {
    //
    if ( 0 == dataset.get( entry)  ) { break ; }                    // BREAK
    //
    // calculate the function 
    const double value = func ( &dataset ) ;
    // 
    var.setVal ( value ) ;
    //
    new_dataset.add ( varset ) ;
  }
  // 
  // merge  two  datasets 
  dataset.merge ( &new_dataset ) ;
  //
  const RooArgSet*  vars = dataset.get(0) ;
  if ( nullptr == vars ) { return nullptr ; }
  //
  const RooAbsArg*  nvar = vars->find ( name.c_str() );
  if  ( nullptr == nvar ) { return nullptr ; }
  //
  return dynamic_cast<const RooAbsReal*> ( nvar ) ;   
}
// ============================================================================
/*  add new variable to dataset
 *  @param  dataset input dataset
 *  @param  name    variable name 
 *  @param  formula formula 
 *  @return the added variable 
 */
// ============================================================================
const RooAbsReal* 
Ostap::Functions::add_var 
( RooDataSet&             dataset , 
  const std::string&      name    , 
  const std::string&      formula ) 
{
  //
  const RooArgSet* vars =  dataset.get(0);
  if ( nullptr == vars ) { return nullptr ; }
  //
  RooArgList lst { *vars } ;
  auto var = std::make_unique<RooFormulaVar>( name.c_str() , formula.c_str() , lst ) ;
  if ( !var ||  !var->ok() ) { return nullptr ; }
  //
  dataset.addColumn ( *var ) ;
  //
  vars = dataset.get(0) ;
  //
  if ( nullptr == vars ) { return nullptr ; }
  //
  const RooAbsArg*  nvar = vars->find ( name.c_str() );
  if  ( nullptr == nvar ) { return nullptr ; }
  //
  return dynamic_cast<const RooAbsReal*> ( nvar ) ;   
}
// ============================================================================
/*  add new variable to dataset, sampled from 1D-histogram
 *  @param  dataset input    dataset
 *  @param  name    variable name 
 *  @param  histo   histogram to be sampled 
 *  @return the added variable 
 *  @see TH1::GetRandom 
 */
// ============================================================================
const RooAbsReal* 
Ostap::Functions::add_var 
( RooDataSet&             dataset , 
  const std::string&      name    , 
  const TH1&              histo   ) 
{
  //
  const TH1* h1 = &histo ;
  if ( nullptr != dynamic_cast<const TH2*> ( h1 ) ) { return nullptr ; }
  //
  //
  RooRealVar var         { name.c_str() , ""        , 0.0 } ;
  RooArgSet  varset      { var          ,"one var" } ;
  RooDataSet new_dataset { "" , ""      , varset   } ;
  //
  // loop over events in the input data set 
  const unsigned long nEntries = dataset.numEntries() ;
  for ( unsigned long entry = 0 ; entry < nEntries ; ++entry )   
  {
    //
    if ( 0 == dataset.get( entry)  ) { break ; }                    // BREAK
    //
    // calculate the function 
    const double value = histo.GetRandom() ;
    // 
    var.setVal ( value ) ;
    //
    new_dataset.add ( varset ) ;
  }
  // 
  // merge  two  datasets 
  dataset.merge ( &new_dataset ) ;
  //
  const RooArgSet*  vars = dataset.get( 0 ) ;
  if ( nullptr == vars ) { return nullptr ; }
  //
  const RooAbsArg*  nvar = vars->find ( name.c_str() );
  if  ( nullptr == nvar ) { return nullptr ; }
  //
  return dynamic_cast<const RooAbsReal*> ( nvar ) ;   
}
// ============================================================================
/*  add new variables to dataset, sampled from 2D-histogram
 *  @param  dataset input    dataset
 *  @param  namex   variable name 
 *  @param  namey   variable name 
 *  @param  histo   histogram to be sampled 
 *  @return the added variable 
 *  @see TH2::GetRandom2
 */
// ============================================================================
const RooAbsReal* 
Ostap::Functions::add_var 
( RooDataSet&             dataset , 
  const std::string&      namex   , 
  const std::string&      namey   , 
  const TH2&              histo   ) 
{
  //
  const TH2* h = &histo ;
  if ( nullptr != dynamic_cast<const TH3*> ( h ) ) { return nullptr ; }
  //
  RooRealVar varx        { namex.c_str() , "" , 0.0 } ;
  RooRealVar vary        { namey.c_str() , "" , 0.0 } ;
  
  RooArgSet  varset      { varx, vary    ,"one var" } ;
  RooDataSet new_dataset { "" , ""       , varset   } ;
  //
  double value_x = 0 ;
  double value_y = 0 ;
  //
  TH2* h2 = const_cast<TH2*> ( h ) ;
  //
  // loop over events in the input data set 
  const unsigned long nEntries = dataset.numEntries() ;
  for ( unsigned long entry = 0 ; entry < nEntries ; ++entry )   
  {
    //
    if ( 0 == dataset.get( entry)  ) { break ; }                    // BREAK
    //
    // calculate the function 
    h2->GetRandom2 ( value_x , value_y ) ;
    // 
    varx.setVal ( value_x ) ;
    vary.setVal ( value_y ) ;
    //
    new_dataset.add ( varset ) ;
  }
  // 
  // merge  two  datasets 
  dataset.merge ( &new_dataset ) ;
  //
  const RooArgSet*  vars = dataset.get( 0 ) ;
  if ( nullptr == vars ) { return nullptr ; }
  //
  const RooAbsArg*  nvar = vars->find ( namey.c_str() );
  if  ( nullptr == nvar ) { return nullptr ; }
  //
  return dynamic_cast<const RooAbsReal*> ( nvar ) ;   
}
// ============================================================================
/*  add new variables to dataset, sampled from 2D-histogram
 *  @param  dataset input    dataset
 *  @param  namex   variable name 
 *  @param  namey   variable name 
 *  @param  namez   variable name 
 *  @param  histo   histogram to be sampled 
 *  @return the added variable 
 *  @see TH2::GetRandom2
 */
// ============================================================================
const RooAbsReal* 
Ostap::Functions::add_var 
( RooDataSet&             dataset , 
  const std::string&      namex   , 
  const std::string&      namey   , 
  const std::string&      namez   , 
  const TH3&              histo   ) 
{
  //
  RooRealVar varx        { namex.c_str() , "" , 0.0 } ;
  RooRealVar vary        { namey.c_str() , "" , 0.0 } ;
  RooRealVar varz        { namez.c_str() , "" , 0.0 } ;
  //
  RooArgSet  varset      { varx, vary , varz ,"one var" } ;
  RooDataSet new_dataset { "" , ""       , varset   } ;
  //
  double value_x = 0 ;
  double value_y = 0 ;
  double value_z = 0 ;
  //
  TH3* h3 = const_cast<TH3*> ( &histo) ;
  //
  // loop over events in the input data set 
  const unsigned long nEntries = dataset.numEntries() ;
  for ( unsigned long entry = 0 ; entry < nEntries ; ++entry )   
  {
    //
    if ( 0 == dataset.get( entry)  ) { break ; }                    // BREAK
    //
    // calculate the function 
    h3->GetRandom3 ( value_x , value_y , value_z ) ;
    // 
    varx.setVal ( value_x ) ;
    vary.setVal ( value_y ) ;
    varz.setVal ( value_z ) ;
    //
    new_dataset.add ( varset ) ;
  }
  // 
  // merge  two  datasets 
  dataset.merge ( &new_dataset ) ;
  //
  const RooArgSet*  vars = dataset.get( 0 ) ;
  if ( nullptr == vars ) { return nullptr ; }
  //
  const RooAbsArg*  nvar = vars->find ( namez.c_str() );
  if  ( nullptr == nvar ) { return nullptr ; }
  //
  return dynamic_cast<const RooAbsReal*> ( nvar ) ;   
}
// ============================================================================
//                                                                      The END 
// ============================================================================
