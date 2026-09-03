// ============================================================================
#ifndef OSTAP_EIGENSYSTEM_H 
#define OSTAP_EIGENSYSTEM_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL 
// ============================================================================
#include <vector>
#include <array>
// ============================================================================
// ROOT 
// ============================================================================
#include "Math/SMatrix.h"
#include "Math/SVector.h"
// ============================================================================
// Ostap
// ============================================================================
#include "Ostap/StatusCode.h"
#include "Ostap/LinAlg.h"
// ============================================================================
/** @file Ostap/EigenSystem.h
 *  Helper class with allows to find eigenvalues and eigenvectors 
 *  @author Vanya BELYAEV
 *  @date   2006-05-24
 */
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
      // Forward declarations 
      // ======================================================================
      class Matrix ;
      class Vector ;      
      // ======================================================================
      /** @class EigenSystem Ostap/EigenSystem.h
       *  Helper class to find eigenvalues and eigenvectors
       *  @author Vanya BELYAEV
       *  @date   2006-05-24
       */
      class EigenSystem
      {
        // ====================================================================
      public:
        // ====================================================================
        EigenSystem  ( const unsigned short N  = 0 ) ;
        /// copy constructor
        EigenSystem  ( const EigenSystem&  ) ;
        /// move constructor
        EigenSystem  (       EigenSystem&& ) ;        
        /// destructor 
        ~EigenSystem () ;
        // ====================================================================
      public:
        // ====================================================================
        /** Get the eigenvalues of symmetrical matrix
         *  @param matrix    (INPUT)  input matrix
         *  @param values    (UPDATE) output vector of eigenvalues
         *  @param sorted    (INPUT)  get the eigenvalues sorted ?
         *  @param ascending (INPUT)  sorting order  
         *  @return Status code 
         */
        Ostap::StatusCode
        eigenValues
        ( const Matrix& matrix            ,
          Vector&       values            ,
          const bool    sorted    = false , 
          const bool    ascending = true  ) const ;
        // ====================================================================
        /** Get the eigenvalues and eigenvectors of symmetrical matrix
         *  @param matrix    (INPUT)  input matrix
         *  @param values    (UPDATE) output vector of eigenvalues
         *  @param vectors   (UPDATE) output matrix where each column is eigen-vector
         *  @param sorted    (INPUT)  get the eigenvaleus sorted ?
         *  @param ascending (INPUT)  sorting order 
         *  @return Status code 
         */
        Ostap::StatusCode
        eigenVectors
        ( const Matrix& matrix            ,
          Vector&       values            ,
          Matrix&       vectors           ,
          const bool    sorted    = false ,
          const bool    ascending = true  ) const ;
        // =====================================================================
        /** Get the eigenvalues and eigenvectors of symmetrical matrix
         *  @param matrix    (INPUT)  input matrix
         *  @param values    (UPDATE) output vector of eigenvalues
         *  @param vectors   (UPDATE) output matrix where each column is eigen-vector
         *  @param sorted    (INPUT)  get the eigenvaleus sorted ?
         *  @param ascending (INPUT)  sorting order 
         *  @return Status code 
         */
        Ostap::StatusCode
        eigenVectors
        ( const Matrix&        matrix            ,
          Vector&              values            ,
          std::vector<Vector>& vectors           ,
          const bool           sorted    = false ,
          const bool           ascending = true  ) const ;
        // =====================================================================
      public : // symmetrix SMatrix 
        // =====================================================================
        /** Get the eigenvalues of symmetrical matrix
         *  @param  matrix    (INPUT)  input matrix
         *  @param  values    (UPDATE) output vector of eigenvalues
         *  @param  sorted    (INPUT)  get the eigenvalues sorted ?
         *  @param  ascending (INPUT)  sorting order  
         *  @return Status code 
         */
        template <class T,unsigned int D>
        inline Ostap::StatusCode
        eigenValues
        ( const ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >& matrix ,
          ROOT::Math::SVector<T,D>&                                     values , 
          const bool           sorted    = false ,
          const bool           ascending = true  ) const 
        {          
          /// copy S-matrix into GSL-matrix 
          const Matrix m { D , D , matrix.Array() , D * ( D + 1 ) / 2 } ;
          Vector v { D } ;
          // calculate the eigenvalues 
          const Ostap::StatusCode sc = this -> eigenValues ( m , v , sorted , ascending ) ;
          if ( sc.isFailure() ) { return sc ; }
          /// copy back
          std::copy ( v.begin () , v.end () , values.begin () ) ;
          return Ostap::StatusCode::SUCCESS ;
        }
        // =====================================================================
        /** Get the eigenvalues and eigenvectors of symmetrical matrix
         *  @param matrix    (INPUT)  input matrix
         *  @param values    (UPDATE) output vector of eigenvalues
         *  @param vectors   (UPDATE) output matrix where each column is eigen-vector
         *  @param sorted    (INPUT)  get the eigenvaleus sorted ?
         *  @param ascending (INPUT)  sorting order 
         *  @return Status code 
         */
        template <class T,unsigned int D>
        Ostap::StatusCode
        eigenVectors
        ( const ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >& matrix  ,
          ROOT::Math::SVector<T,D>&                                     values  , 
          ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepStd<T,D,D> >&     vectors , 
          const bool           sorted    = false ,
          const bool           ascending = true  ) const 
        {
          /// copy S-matrix into GSL-matrix 
          const Matrix m { D , D , matrix.Array() , D * ( D + 1 ) / 2 } ;
          /// vector of eigenvalues  
          Vector v { D } ;
          /// matrix of eigenvectors 
          Matrix b { D , D } ;
          // calculate the eigenvalues & eigenvectors 
          const Ostap::StatusCode sc = this -> eigenVectors ( m , v , b , sorted , ascending ) ;
          if ( sc.isFailure() ) { return sc ; }
          /// copy back
          std::copy ( v.begin () , v.end () , values.begin  () ) ;
          std::copy ( b.begin () , b.end () , vectors.begin () ) ;
          return Ostap::StatusCode::SUCCESS ;          
        }
        // =====================================================================
        /** Get the eigenvalues and eigenvectors of symmetrical matrix
         *  @param matrix    (INPUT)  input matrix
         *  @param values    (UPDATE) output vector of eigenvalues
         *  @param vectors   (UPDATE) output matrix where each column is eigen-vector
         *  @param sorted    (INPUT)  get the eigenvaleus sorted ?
         *  @param ascending (INPUT)  sorting order 
         *  @return Status code 
         */
        template <class T,unsigned int D>
        Ostap::StatusCode
        eigenVectors
        ( const ROOT::Math::SMatrix<T,D,D,ROOT::Math::MatRepSym<T,D> >& matrix  ,
          ROOT::Math::SVector<T,D>&                                     values  , 
          // std::vector<ROOT::Math::SVector<T,D>> &                       vectors ,
          std::array<ROOT::Math::SVector<T,D>,D> &                      vectors ,
          const bool           sorted    = false ,
          const bool           ascending = true  ) const 
        {
          /// copy S-matrix into GSL-matrix 
          const Matrix m { D , D , matrix.Array() , D * ( D + 1 ) / 2 } ;
          /// vector of eigenvalues  
          Vector v { D } ;
          /// vector of eigenvectors 
          std::vector<Vector> b { D } ;
          // calculate the eigenvalues & eigenvectors 
          const Ostap::StatusCode sc = this -> eigenVectors ( m , v , b , sorted , ascending ) ;
          if ( sc.isFailure() ) { return sc ; }
          /// copy back
          std::copy ( v.begin () , v.end () , values.begin  () ) ;
          std::copy ( b.begin () , b.end () , vectors.begin () ) ;
          return Ostap::StatusCode::SUCCESS ;          
        }
        // =====================================================================
      public :  
        // =====================================================================
        /// assignement operator
        EigenSystem& operator= ( const EigenSystem&  right ) ;
        /// assignement operator
        EigenSystem& operator= (       EigenSystem&& right ) ; 
        // =====================================================================
        /// swap two objects
        void swap ( EigenSystem& right ) ;
        // =====================================================================
      public : 
        // =====================================================================
        /// release allocated workspaces 
        void release () ;
        // =====================================================================
      private : 
        // =====================================================================
        void release1 () const ;
        void release2 () const ;
        // =====================================================================
      private :
        // =====================================================================
        /// SYMM&SYMMV  workspaces 
        // mutable gsl_eigen_symm_workspace*  m_ws_symm  { nullptr } ;
        // mutable gsl_eigen_symmv_workspace* m_ws_symmv { nullptr } ;
        mutable void* m_ws_symm  { nullptr } ;
        mutable void* m_ws_symmv { nullptr } ;        
        // =====================================================================        
      } ;
      // =======================================================================
      /// swap two objects
      inline void swap ( EigenSystem& a , EigenSystem& b ) { a.swap ( b ) ; } 
      // =======================================================================

      // =======================================================================
    } //                                                   end of namespace GSL
    // ========================================================================
  } //                                                    end of namespace Math 
  // ==========================================================================
} //                                                     end of namespace Ostap
// ============================================================================
#endif // LHCBMATH_EIGENSYSTEM_H
// ============================================================================
//                                                                      The END 
// ============================================================================
