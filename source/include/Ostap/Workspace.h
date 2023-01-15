// ============================================================================
#ifndef OSTAP_WORKSPACE_H 
#define OSTAP_WORKSPACE_H 1
// ============================================================================
// Include files
// ============================================================================
// ROOT 
// ============================================================================
#include "RVersion.h"
// ============================================================================
namespace Ostap
{
  // ==========================================================================
  namespace Math
  {
    // ========================================================================
    /** @class WorkSpace Ostap/Workspace.h
     *  helper utility to keep the integration workspace for GSL integration
     *  @author Vanya Belyaev Ivan.Belyaev@itep.ru
     *  @date 2011-12-03
     */
    class WorkSpace
    {
    public:
      // ======================================================================
      /** constructor
       *  @param size    size of the main integration worksoace 
       *  @param cquad   size of integration workspace for CQUAD   integrator 
       *  @param romberg size of integration workspace for Romberg integrator 
       */
      WorkSpace  
      ( const std::size_t    size         = 0 , 
        const unsigned short size_cquad   = 0 , 
        const unsigned short size_romberg = 0 ) ;
      /// (fictive) copy constructor
      WorkSpace  ( const WorkSpace&  right ) ;
      /// move constructor 
      WorkSpace  (       WorkSpace&& right ) ;
      /// destructor
      ~WorkSpace () ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the main integration workspace
      void* workspace         () const ; // get the main integration workspace
      /// get the integration workspace for CQUAD inetgrator 
      void* workspace_cquad   () const ; // get CQUAD integration workspace
      /// get the integration workspace for Romberg inetgrator 
      void* workspace_romberg () const ; // get Romberg integration workspace
      // ======================================================================
    public:
      // ======================================================================
      /// get the size of main allocated workspace 
      inline std::size_t size         () const { return m_size         ; }
      // ======================================================================
      /// get the size of allocated workspace for CQUAD integrator  
      inline std::size_t size_cquad   () const { return m_size_cquad   ; }
      // ======================================================================
      /// get the size of allocated workspace for Romberg integrator  
      inline std::size_t size_romberg () const { return m_size_romberg ; }
      // ======================================================================
    public:
      // ======================================================================
      /// resize the main integration workspace 
      std::size_t resize         ( const std::size_t newsize ) ;
      /// resize the integration workspace for CQUAD integrator 
      std::size_t resize_cquad   ( const std::size_t newsize ) ;
      /// resize the integration workspace for Romberg integrator 
      std::size_t resize_romberg ( const std::size_t newsize ) ;
      // ======================================================================
    public:
      // ======================================================================
      /// (fictive) assignement operator
      WorkSpace& operator= ( const WorkSpace&  right ) ;
      /// move      assignement operator
      WorkSpace& operator= (       WorkSpace&& right ) ;
      // ======================================================================
    public:
      // ======================================================================
      void swap ( WorkSpace& right ) ;
      // ======================================================================
    private:
      // ======================================================================
      /// the actual GSL-workspace
      /// mutable char*  m_workspace ;  /// the actual GSL-workspace
      // char* here to please dictionary generator...
      /// mutable char* m_workspace  { nullptr } ; //! the actual GSL-workspace
      mutable void* m_workspace          { nullptr } ; //! the actual GSL-workspace
      mutable void* m_workspace_cquad    { nullptr } ; //! the actual GSL-workspace
      mutable void* m_workspace_romberg  { nullptr } ; //! the actual GSL-workspace
      // ======================================================================
      /// size of the main allocated workspace 
      std::size_t   m_size         { 0 } ; // size of main workspace 
      /// size of the workspace for CQUAD integrator 
      std::size_t   m_size_cquad   { 0 } ; // size of CQUAD workspace 
      /// size of the workspace for Romberg integrator 
      std::size_t   m_size_romberg { 0 } ; // size of Romberg workspace 
      // ======================================================================
    } ;
    // ========================================================================
    /// swap two integration workspaces 
    inline void swap ( WorkSpace& a , WorkSpace& b ) { a.swap ( b ) ; }
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_WORKSPACE_H
// ============================================================================
