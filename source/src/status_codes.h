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
    INVALID_TREE             = 710 ,
    INVALID_DATA             = 711 ,
    INVALID_FORMULA          = 712 ,
    INVALID_FORMULAE         = 713 ,
    INVALID_BRANCH_NAME      = 714 , 
    INVALID_BRANCH           = 715 , 
    CANNOT_CREATE_BRANCH     = 716 , 
    CANNOT_CREATE_FORMULA    = 717 , 
    INVALID_NAME             = 718 , 
    INVALID_HISTO            = 719 , 
    // ========================================================================
    INVALID_TREEFUNCTION     = 720 ,
    INVALID_TH3              = 721 , 
    INVALID_TH2              = 722 , 
    INVALID_TH1              = 723 ,
    INVALID_XAXIS            = 724 ,
    INVALID_YAXIS            = 725 ,
    INVALID_ZAXIS            = 726 ,
    INVALID_BUFFER           = 727 ,
    MISMATCH_TREE            = 728  , 
    //
    INVALID_ABSDATA          = 730 ,
    INVALID_ARGSET           = 731 ,
    INVALID_ARGLIST          = 732 ,
    INVALID_ABSARG           = 733 ,   
    INVALID_ABSREAL          = 734 ,   
    INVALID_OBSERVABLE       = 735 , 
    INVALID_OBSERVABLES      = 736 ,
    INVALID_ABSPDF           = 737 ,
    //
    INVALID_PDF              = 740 ,
    INVALID_FITRESULT        = 741 ,
    //
    INVALID_PYSELF           = 750 ,
    //
    UNDEFINED_METHOD         = 760 ,
    INVALID_VARIABLE         = 761 ,
    INVALID_CALLABLE         = 762 ,
    INVALID_PYOBJECT         = 763 ,
    //
    ERROR_PYTHON             = 770 ,
    //
    INVALID_KNOTS            = 780 ,
    INVALID_PARS             = 781 , 
    INVALID_RANGE            = 782 , 
    INVALID_ORDER            = 783 , 
    //
    INVALID_TMATRIX          = 790 ,
    INVALID_TVECTOR          = 791 ,
    INVALID_GMATRIX          = 792 ,
    INVALID_GVECTOR          = 793 ,
    INVALID_PERMUTATION      = 794 ,
    INVALID_ROWINDEX         = 795 ,
    INVALID_COLINDEX         = 796 ,
    INVALID_SCALE            = 797 ,
    //
    INVALID_CACHE            = 810 ,     
    INVALID_CHEBYSHEV        = 815 ,
    //
    INVALID_INTEGRATION_CODE = 820 ,
    INVALID_MAXVAL_CODE      = 821 ,
    //
    INVALID_PARAMETER        = 825 , 
    INVALID_PARAMETERS       = 826 , 
    INVALID_MINMAX           = 827 ,
    //
    INVALID_ENTRY            = 830 ,
    INVALID_EVENT            = 831 ,    
    //
    INVALID_KERNEL           = 835 , 
    INVALID_ECDF             = 836 , 
    INVALID_WECDF            = 837 , 
    INVALID_SMOOTH           = 838 ,
    // 
    INVALID_QUANTILE         = 840 ,
    INVALID_QUANTILE_INDEX   = 841 ,
    INVALID_PROBABILITY      = 842 ,
    INVALID_PROBABILITIES    = 843 ,
    INVALID_SELECTION        = 844 ,
    INVALID_SELECTIONS       = 845 ,  
    INVALID_COUNTER          = 846 ,
    INVALID_COUNTERS         = 847 ,
    //
    INVALID_FORMULA_CALL     = 850 ,
    //
    INVALID_WEIGHT           = 861 ,
    INVALID_SUMWEIGHT        = 862 ,
    //
    INVALID_SIZE             = 870 ,    
    INVALID_SUMW             = 871 ,    
    INVALID_SUMW2            = 872 ,
    //
    INVALID_DATA_WEIGHT      = 873 ,    
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
