// ============================================================================
#ifndef OSTAP_GSL_HELPERS_H 
#define OSTAP_GSL_HELPERS_H 1
// ============================================================================
// Include files
// ============================================================================
// GSL 
// ============================================================================
#include "gsl/gsl_linalg.h"
// ============================================================================
namespace Ostap 
{
  // ==========================================================================
  /** @class GSL_Matrix
   *  Internal class to  hold GSL-Matrix
   */
  class GSL_Matrix
  {
  public : 
    // ========================================================================
    struct Zero     {} ;
    struct Identity {} ;
    // ========================================================================        
  public : 
    // ========================================================================
    /// allocate GSL-matrix 
    GSL_Matrix  ( const unsigned short N1    , 
                  const unsigned short N2    ) ;
    /// allocate GSL-matrix and initialize all elements to value 
    GSL_Matrix  ( const unsigned short N1    , 
                  const unsigned short N2    , 
                  const double         value ) ;
    /// allocate GSL-matrix and initialize all elements to zero 
    GSL_Matrix  ( const unsigned short N1    , 
                  const unsigned short N2    , 
                  const Zero           /* zero */ ) ;
    /// allocate "identity" GSL-matrix
    GSL_Matrix  ( const unsigned short N1    , 
                  const unsigned short N2    , 
                  const Identity       /* id  */ ) ;
    /// allocate square GSL-matrix and initialize all elements to zero 
    GSL_Matrix  ( const unsigned short N1    , 
                  const Zero           zero  ) 
      : GSL_Matrix ( N1 ,  N1 , zero ) {}
    /// allocate square identity GSL-matrix 
    GSL_Matrix  ( const unsigned short N1    , 
                  const Identity       id    ) 
      : GSL_Matrix ( N1 ,  N1 , id   ) {}
    /// copy constructor 
    GSL_Matrix  ( const GSL_Matrix&  right ) ;
    /// move constructor 
    GSL_Matrix  (       GSL_Matrix&& right ) ;
    ///  destructor: free  GSL-matrix 
    ~GSL_Matrix () ;
    // ========================================================================
    /// no default constructor 
    GSL_Matrix  () = delete ;
    /// No copy assignement! 
    GSL_Matrix& operator= ( const GSL_Matrix&  ) = delete ;
    /// No move assignement! 
    GSL_Matrix& operator= (       GSL_Matrix&& ) = delete ;
    // ========================================================================
  public:
    // ========================================================================
    // get the matrix
    inline       gsl_matrix* matrix  ()       { return m_matrix ; }
    // get the matrix
    inline const gsl_matrix* matrix  () const { return m_matrix ; }
    // ========================================================================
    /// get the mattix element 
    double      get    ( const unsigned short n1 , 
                         const unsigned short n2 ) const 
    { return gsl_matrix_get ( m_matrix , n1 , n2 ) ; }
    /// set the matrix element 
    void        set    ( const unsigned short n1 , 
                         const unsigned short n2 , 
                         const double   value    )
    {        gsl_matrix_set ( m_matrix , n1 , n2 , value ) ; }
    // ========================================================================
  private:
    // ========================================================================
    /// the  actual pointer to GSL-matrix 
    gsl_matrix* m_matrix { nullptr } ; // the  actual pointer to GSL-matrix
    // ========================================================================
  };
  // ==========================================================================
  /** @class GSL_Vector
   *  Internal class to  hold GSL-Vector
   */
  class GSL_Vector
  {
  public: 
    // ========================================================================
    typedef GSL_Matrix::Zero Zero ;
    // ========================================================================
  public: 
    // ========================================================================
    /// allocate vector 
    GSL_Vector ( const short unsigned N     ) ;
    /// allocate vector 
    GSL_Vector ( const short unsigned N     , 
                 const double         value ) ;
    /// allocate vector 
    GSL_Vector ( const unsigned short N     , 
                 const Zero     /* zero */  ) ;
    /// copy constructor 
    GSL_Vector  ( const GSL_Vector&  right ) ;
    /// move constructor 
    GSL_Vector  (       GSL_Vector&& right ) ;
    /// destructor : free GSL-Vector
    ~GSL_Vector () ;
    // =======================================================================
    /// no default constructor 
    GSL_Vector() = delete ;
    /// No copy assignement! 
    GSL_Vector& operator= ( const GSL_Vector&  ) = delete ;
    /// No move assignement! 
    GSL_Vector& operator= (       GSL_Vector&& ) = delete ;
    // ========================================================================
  public:
    // ========================================================================
    // get the vector
    inline       gsl_vector* vector ()       { return m_vector ; }
    // get the vector
    inline const gsl_vector* vector () const { return m_vector ; }
    // ========================================================================
    /// get the vector element 
    double      get    ( const unsigned short n ) const 
    { return gsl_vector_get ( m_vector , n ) ; }
    /// set the vector element 
    void        set    ( const unsigned short n , 
                         const double   value    )
    {        gsl_vector_set ( m_vector , n , value ) ; }
    // ========================================================================
  private:
    // ========================================================================
    /// the  actual pointer to GSL-vector 
    gsl_vector* m_vector { nullptr } ; // the  actual pointer to GSL-vector
    // ========================================================================
  };
  // ==========================================================================
  /** @class   GSL_Permutation
   *  Internal class to keep GSL-permuation
   */
  class  GSL_Permutation
  {
  public:
    // ========================================================================
    /// constructor: allocate the permutation 
    GSL_Permutation ( const unsigned short N ) ;
    /// destructor: free permutation 
    ~GSL_Permutation() ;
    // ========================================================================
    GSL_Permutation () = delete ;
    GSL_Permutation ( const GSL_Permutation&  )  = delete ;
    GSL_Permutation (       GSL_Permutation&& )  = delete ;
    /// No copy assignement! 
    GSL_Permutation& operator= ( const GSL_Permutation&  ) = delete ;
    /// No move assignement! 
    GSL_Permutation& operator= (       GSL_Permutation&& ) = delete ;
    // =======================================================================
  public:
    // =======================================================================
    // get the permutation 
    inline       gsl_permutation* permutation ()       { return m_permutation ; }
    // get the permutation 
    inline const gsl_permutation* permutation () const { return m_permutation ; }
    // ========================================================================
  private :
    // ========================================================================
    /// the  actual pointer to GSL-permutation
    gsl_permutation* m_permutation { nullptr } ; // the  actual pointer to GSL-vector
    // ========================================================================    
  };
  // ==========================================================================
} //                                                 The end of namespace Ostap 
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_GSL_HELPERS_H
// ============================================================================
