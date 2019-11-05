// ============================================================================
#ifndef OSTAP_WORKSPACE_H 
#define OSTAP_WORKSPACE_H 1
// ============================================================================
// Include files
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
      /// constructor
      WorkSpace ( const std::size_t size  = 0 );
      /// (fictive) copy constructor
      WorkSpace ( const WorkSpace&  right );
      /// move constructor 
      WorkSpace (       WorkSpace&& right );
      /// destructor
      ~WorkSpace () ;
      // ======================================================================
    public:
      // ======================================================================
      /// get the integration workspace
      void* workspace () const ;               // get the integration workspace
      // ======================================================================
      /// get the size of allocated workspace 
      // ======================================================================
      std::size_t size () const { return m_size ; }
      // ======================================================================
    public:
      // ======================================================================
      /// resize the workspace 
      std::size_t resize  ( const std::size_t newsize ) ;
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
      // mutable char*  m_workspace ;  /// the actual GSL-workspace
      // char* here to please dictionary generator...
      mutable char* m_workspace  { nullptr } ; /// the actual GSL-workspace
      // ======================================================================
      /// size of the allocated workspace 
      std::size_t   m_size       { 0 } ;   /// size of the allocated workspace 
      // ======================================================================
    } ;
    // ========================================================================
  } //                                         The end of namespace Ostap::Math
  // ==========================================================================
} //                                                 The end of namespace Ostap
// ============================================================================
//                                                                      The END 
// ============================================================================
#endif // OSTAP_WORKSPACE_H
// ============================================================================
