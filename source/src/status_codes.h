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
    INVALID_TREE          = 710 ,
    INVALID_FORMULA       = 711 ,
    INVALID_FORMULAE      = 712 ,
    INVALID_BRANCH_NAME   = 713 , 
    INVALID_BRANCH        = 714 , 
    CANNOT_CREATE_BRANCH  = 715 , 
    CANNOT_CREATE_FORMULA = 716 , 
    INVALID_NAME          = 717 , 
    // ========================================================================
    INVALID_TREEFUNCTION  = 720 ,
    INVALID_TH3           = 721 , 
    INVALID_TH2           = 722 , 
    INVALID_TH1           = 723 , 
    INVALID_BUFFER        = 724 ,
    MISMATCH_TREE         = 725  , 
    //
    INVALID_ABSDATA       = 730 ,
    INVALID_ARGSET        = 731 ,
    INVALID_ABSARG        = 732 ,   
    INVALID_ABSREAL       = 733 ,   
    INVALID_OBSERVABLE    = 734 , 
    INVALID_OBSERVABLES   = 735 ,
    INVALID_ABSPDF        = 736 ,
    //
    INVALID_PDF           = 740 ,
    INVALID_FITRESULT     = 741 ,
    //
    INVALID_PYSELF        = 750 ,
    //
    UNDEFINED_METHOD      = 760 ,
    INVALID_VARIABLE      = 761 ,
    INVALID_CALLABLE      = 762 ,
    INVALID_PYOBJECT      = 763 ,
    //
    ERROR_PYTHON          = 770 ,
    //
    INVALID_KNOTS         = 780 ,
    INVALID_PARS          = 781 , 
    INVALID_RANGE         = 782 , 
    INVALID_ORDER         = 783 , 
    INVALID_DATA          = 784 ,
    //
    INVALID_TMATRIX       = 790 ,
    INVALID_TVECTOR       = 791 ,
    INVALID_GMATRIX       = 792 ,
    INVALID_GVECTOR       = 793 ,
    INVALID_PERMUTATION   = 794 ,
    INVALID_ROWINDEX      = 795 ,
    INVALID_COLINDEX      = 796 ,
    INVALID_SCALE         = 797 ,
    // =========================================================================
    // =========================================================================
    ERROR_ROOT            = 1000000 ,
    ERROR_GSL             = 2000000 ,    
    // =========================================================================
  }; // ========================================================================
  // ===========================================================================
} // the end of anonymous namespace 
// =============================================================================
#endif  // OSTAP_STATUS_CODES_H
// =============================================================================
//                                                                       The END 
// =============================================================================
