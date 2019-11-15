// Include files 
// ============================================================================
// Incldue files 
// ============================================================================
// STD &STL
// ============================================================================
#include <memory>
#include <cmath>
#include <limits>
// ============================================================================
// ROOT/RooFit
// ============================================================================
#include "RooRealVar.h"
#include "RooAbsData.h"
#include "RooDataSet.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/GetWeight.h"
// ============================================================================
/** @dile 
 * implementation file for function Ostap::Utils::getWeight
 * @see Ostap::Utils::getWeight
 * @date 2018-12-13 
 * @author Vanya BELYAEV Ivan.Belyaev@itep.ru
 */
// ============================================================================
namespace 
{
  // ==========================================================================
  /// helper class to get the protected info 
  class AuxDataSet : public RooDataSet
  {
  public:
    // ========================================================================
    ///  constructor from the dataset 
    AuxDataSet ( const RooDataSet& data ) : RooDataSet ( data ) {};
    ///   destructor 
    ~AuxDataSet(){}
  public:
    // =======================================================================
    /// get the weight variable 
    std::string wgtVar() const 
    { return RooDataSet::_wgtVar ? RooDataSet::_wgtVar->GetName() : "" ; }
    // ========================================================================
  } ;
  // ==========================================================================
  static_assert ( std::numeric_limits<double>::is_specialized        ,
                  "std::numeric_limits<double> is not specialized"   ) ;
  // ==========================================================================
  constexpr double s_VMIN = -0.9 * std::numeric_limits<double>::max ()  ;
  constexpr double s_VMAX =  0.9 * std::numeric_limits<double>::max ()  ;
  // ==========================================================================
}
// ============================================================================
/*  get the weight variable from data set (if possible)
 *  @param data (INPUT) dataset
 *  @return weigth variable, when possible 
 */
// ============================================================================
std::string Ostap::Utils::getWeight ( const RooAbsData* data ) 
{
  if ( nullptr == data ) { return "" ; }
  const RooDataSet* ds = dynamic_cast<const RooDataSet*>( data ) ;
  if ( nullptr == ds   ) { return ""  ; }
  //
  std::unique_ptr<RooAbsData> cloned { ds->emptyClone() } ;
  
  AuxDataSet aux { *dynamic_cast<RooDataSet*>( cloned.get() ) } ;
  //
  return aux.wgtVar() ;
}
// ============================================================================
/*  make un unweighted dataset from weighted one 
 *  @param weighted_data   (INPUT) input dataset 
 *  @return unweighted dataset, if and when possible
 */
// ============================================================================
std::pair<RooDataSet*,std::string> 
Ostap::Utils::unweight 
( const RooAbsData* weighted_data , 
  std::string       weight_var    ) 
{
  if ( nullptr ==  weighted_data     ) { return std::make_pair ( nullptr, "" ) ; } //  RETURN
  if ( ! weighted_data->isWeighted() ) { return std::make_pair ( nullptr, "" ) ; } //  RETURN
  //
  if ( weight_var.empty () ) { weight_var  = getWeight ( weighted_data ) ; }
  if ( weight_var.empty ()            ) { return std::make_pair ( nullptr, "" ) ; } //  RETURN 
  //
  static const std::string prefix { "unweighted_" } ;
  const std::string name  = prefix + weighted_data->GetName  () ;
  const std::string title = prefix + weighted_data->GetTitle () ;
  //
  RooRealVar weight ( weight_var.c_str() , "weight-variable" , s_VMIN , s_VMAX ) ;
  //
  RooArgSet  vset   ( *weighted_data->get() ) ;
  vset.add ( weight ) ;
  auto newdata = std::make_unique<RooDataSet> ( name .c_str ()      , 
                                                title.c_str ()      , 
                                                vset                ) ;
  //
  const unsigned long nEntries = weighted_data->numEntries () ;
  // loop 
  for ( unsigned long entry = 0 ; entry < nEntries ; ++entry )   
  {  
    if ( nullptr == weighted_data->get ( entry)  ) { break ; }  // BREAK
    Double_t low , high ;
    weighted_data->weightError ( low , high ) ;
    weight.setVal       ( weighted_data->weight () ) ;
    weight.setAsymError ( low , high ) ;
    newdata->add ( vset ) ;
  }
  //
  return std::make_pair ( newdata.release() , weight_var ) ;
}
// ============================================================================
//                                                                      The END 
// ============================================================================
