# @date   2010-04-30
# =============================================================================
""" Helper module to check "unpickeability" for different objects 
"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2010-04-30"
__version__ = "$Revision$" 
# =============================================================================
__all__ = (
    'DillChecker'    , ## check pickle-ability of objects 
    )
# =============================================================================
from ostap.io.checker import PickleChecker 
# ============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.parallel.dill_checker' )
else                      : logger = getLogger ( __name__                      )
# =============================================================================
DILL_COMMAND = """import sys, dill
with open('%s','rb') as f : dill.load ( f )"""
# =============================================================================
try : # =======================================================================
    # =========================================================================
    import dill 
    # =========================================================================
    ## @class DillChecker
    #  Check if the object can be properly pickled/unpickled
    class DillChecker(PickleChecker) :
        """ Check if the object can be properly pickled/unpickled
        """
        # =====================================================================
        ## check if the object can be properly pickled/unpickled 
        def pickles ( self , *objects ) :
            """ Check of the object can be properly pickled/unpickled
            """
            return self._pickles ( *objects               ,
                                   fun_dumps = dill.dumps ,
                                   fun_loads = dill.loads ) 
        # =====================================================================
        ## Check pickling of an object across another (sub) process
        def pickles_process ( self , *objects , fast = False  ) :
            """ Check pickling of an object across another (sub)process
            """
            return self._pickles_process ( *objects                ,
                                           fun_dump = dill.dump    ,
                                           fun_load = dill.load    ,
                                           command  = DILL_COMMAND ,
                                           fast     = fast         ) ;
        # =========================================================================
        ## add new type into the list of "known-types"
        def add ( self , *ntypes ) :
            """ Add new type into the list of "known-types
            """
            for ntype in ntypes :
                if ntype in self : continue 
                self.EXTRA_TYPES.add ( ntype )
                
    # =========================================================================
except ImportError : # ========================================================
    # =========================================================================
    DillChecker = PickleChecker
    # =========================================================================

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
#                                                                       The END 
# =============================================================================
