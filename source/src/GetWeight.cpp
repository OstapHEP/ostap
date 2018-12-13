// Include files 
// ============================================================================
// Incldue files 
// ============================================================================
// STD &STL
// ============================================================================
#include <memory>
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
//                                                                      The END 
// ============================================================================
