// ============================================================================
#ifndef OSTAP_STATUS_CODES_H
#define OSTAP_STATUS_CODES_H
// ============================================================================
// Incldue files
// ============================================================================
// Ostap 
// ============================================================================
#include "Ostap/StatusCode.h"
// ============================================================================
namespace 
{
  // ==========================================================================
  enum { // =================================================================== 
    // ========================================================================
    INVALID_TREE          = 750 ,
    INVALID_FORMULA       = 751 ,
    INVALID_FORMULAE      = 752 ,
    CANNOT_CREATE_BRANCH  = 753 , 
    CANNOT_CREATE_FORMULA = 754 , 
    // ========================================================================
    INVALID_TREEFUNCTION  = 755 ,
    INVALID_TH2           = 756 , 
    INVALID_TH1           = 757 , 
    INVALID_BUFFER        = 758 ,
    MISMATCH_TREE         = 759 , 
    //
    INVALID_ABSDATA       = 760 ,
    INVALID_ARGSET        = 761 ,
    INVALID_OBSERVABLE    = 762 , 
    INVALID_OBSERVABLES   = 763 ,
    //
    // =========================================================================
  }; // ========================================================================
  // ===========================================================================
} // the end of anonymous namespace 
// =============================================================================
#endif  // OSTAP_STATUS_CODES_H
// =============================================================================
//                                                                       The END 
// =============================================================================
