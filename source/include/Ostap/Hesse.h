// $Id$
// ============================================================================
#ifndef LHCBMATH_HESSE_H 
#define LHCBMATH_HESSE_H 1
// ============================================================================
// Include files
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatusCode.h"
// ============================================================================
// GSL 
// ============================================================================
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    namespace GSL
    {
      // ======================================================================
      /** @class Hesse
       *  evaluate the hessian for the function
       *  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
       *  @date 2012-05-27
       */      
      class Hesse 
      {
      public:
        // ====================================================================
        /// the actual type of function to be used for hessian calculation 
        typedef double (*function) ( const gsl_vector  * x, void* params ) ;
        // ====================================================================
      public:
        // ====================================================================
        /*  constructor with all parameters 
         *  @param f      the function to be used 
         *  @param x      the point for hessian to be evaluated 
         *  @param params the parameters for the function 
         *  @param h the step-size (guess)
         */
        Hesse ( function          f      ,
                const gsl_vector* x      ,
                void*             params , 
                const double      h      ) ;
        /// destructor 
        ~Hesse() ;                                                // destructor 
        // ====================================================================
      private:
        // ====================================================================
        /// the default constructor is disabled 
        Hesse () ;                         // default constrictor is disabled
        /// the copy constructor is disabled 
        Hesse ( const Hesse& ) ;           // the copy constructor is disabled
        /// the assignement operator is disabled 
        Hesse& operator=( const Hesse& ) ; // no assignement
        // ====================================================================
      public:
        // ====================================================================
        Ostap::StatusCode calcHesse () ;
        Ostap::StatusCode calcCov2  () ;
        // ====================================================================
      public:
        // ====================================================================
        /// size of the problem
        std::size_t       size  () const { return m_x->size ; }
        /// get the hesse matrix 
        const gsl_matrix* hesse () const { return m_hesse ; }
        /// get the inverse hesse ("covariance") matrix 
        const gsl_matrix* cov2  () const { return m_cov2  ; }
        // ====================================================================
      private:
        // ====================================================================
        /// the function 
        function          m_func   ; // the function 
        /// the point 
        const gsl_vector* m_x      ; // the point 
        /// parameters 
        void*             m_params ; // the parameters 
        /// step-size 
        double            m_h      ; // the step-size 
        // ====================================================================
      private:
        // ====================================================================
        /// the actual hessian
        gsl_matrix* m_hesse ;                        // the actual hesse matrix 
        /// the axuxillary matrix 
        gsl_matrix* m_aux   ;                        // the axuxillary matrix 
        /// the inverse hesse matrix 
        gsl_matrix* m_cov2  ;                       // the inverse hesse matrix 
        /// 
        // ====================================================================
      private:
        // ====================================================================
        /// helper vector
        gsl_vector* m_a ;  // helper vector
        /// helper vector
        gsl_vector* m_b ;  // helper vector
        // ====================================================================
      } ;
      // ======================================================================
      /** invert the matrix using LU decomposition
       *  @param matrix (UPDATE) the matrix to be inverted 
       *  @param result (UPDATE) the result 
       *  @return status code 
       *  @attention the input matrix will be screwed up!
       *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
       *  @date 2012-05-28
       */
      Ostap::StatusCode invert_LU_1
      ( gsl_matrix* matrix , 
        gsl_matrix* result ) ;
      // ======================================================================
      /** invert the matrix using LU decomposition
       *  @param matrix (INPUT) the matrix to be inverted 
       *  @param result (UPDATE) the result 
       *  @return status code 
       *  @attention the input matrix will be preserved 
       *  @author Vanya BELYAEV  Ivan.Belyaev@itep.ru
       *  @date 2012-05-28
       */
      Ostap::StatusCode invert_LU_2 
      ( const gsl_matrix* matrix ,
        gsl_matrix*       result ) ;
      // ======================================================================
    } //                                      end of namespace Gaudi::Math::GSL 
  } //                                             end of namespace Gaudi::Math
  // ==========================================================================
} //                                                     end of namespace Gaudi
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // LHCBMATH_HESSE_H
// ============================================================================
