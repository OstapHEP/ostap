// ============================================================================
#ifndef OSTAP_VECTOR3DTYPES_H 
#define OSTAP_VECTOR3DTYPES_H 1
// ============================================================================
// Include files
// ============================================================================
// ROOT 
// ============================================================================
#include "Math/Vector3D.h"
// ============================================================================
/** @file Vector3DTypes.h
 *
 *  3D vector typedefs
 *
 *  @author Juan PALACIOS
 *  @date   2005-11-21
 */
namespace Ostap
{
  // ==========================================================================
  typedef ROOT::Math::XYZVector            XYZVector;        ///<  Cartesian 3D vector (double)
  typedef ROOT::Math::Polar3DVector        Polar3DVector;    ///<  Polar 3D vector (double)
  typedef ROOT::Math::RhoEtaPhiVector      RhoEtaPhiVector;  ///<  RhoEtaPhi 3D vector (double)
  typedef ROOT::Math::RhoZPhiVector        RhoZPhiVector;    ///<  RhoZPhi 3D vector (double)
  // ==========================================================================
}
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_VECTOR3DTYPES_H
// ============================================================================
