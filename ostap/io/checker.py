# @date   2010-04-30
# =============================================================================
""" Helper module to check "(unpickeability" for different objects 
"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2010-04-30"
__version__ = "$Revision$" 
# =============================================================================
__all__ = (
    'PickleChecker'    , ## check pickle-ability of objects 
    )
# =============================================================================
import pickle, array, sys
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.io.checker' )
else                      : logger = getLogger ( __name__               )
# =============================================================================
## the basic pickle-types that for sure always can be pickled/unpickled 
PICKLE_TYPES = frozenset ( [ type ( None ) , bool, int, float, str, bytes, bytearray, array.array ] )
# =============================================================================
PICKLE_COMMAND = """import sys, pickle
with open('%s','rb') as f : pickle.load ( f )"""
# =============================================================================
## @class PickleChecker
#  Check if the object/type can be pickled/unpickled
class PickleChecker ( object ) :
    """ Check if the object/type can be pickled/unnpickeld
    """
    MORE_TYPES          = [] 
    EXTRA_TYPES         = [] 
    NONPICKLEABLE_TYPES = [] 
    # ========================================================================
    ## if the object type is already knowns 
    def known ( self , *objtypes ) : return all ( ( o in self ) for o in objtypes ) 
    # ========================================================================
    ## type is knows for (un)pickle problems 
    def notgood ( self , *objtypes ) :
        return any ( ( o in self.NONPICKLEABLE_TYPES ) for o in objtypes )
    # =========================================================================
    ## check if the type is 'known'
    def __contains__ ( self , objtype ) :
        return objtype in PICKLE_TYPES or \
            objtype in self.MORE_TYPES or \
            objtype in self.EXTRA_TYPES
    # =========================================================================
    def unpack ( self , *objects , **kwargs  ) :
        unpacked = [] 
        for obj in objects :
            if   isinstance ( obj , ( list , tuple , set , frozenset ) ) :
                unpacked += list ( self.unpack (  *obj ) ) 
            elif isinstance ( obj ,   dict ) :                
                unpacked += list ( self.unpack ( *obj.keys   () ) )
                unpacked += list ( self.unpack ( *obj.values () ) )
            else :
                unpacked.append ( obj ) 

        ## unpack keyword stuff 
        for key , value in kwargs.items() : unpacked += [ key , value ] 
        
        return tuple ( unpacked )
    
    def _pickles ( self , *objects , fun_dumps , fun_loads ) :
        # =====================================================================
        objects  = self.unpack ( *objects ) 
        if not objects                : return True
        objtypes = set ( type ( o ) for o in objects ) 
        if self.notgood ( *objtypes ) : return False        
        ## check if the object can be properly pickled/unpickled
        if self.known   ( *objtypes ) : return True
        # =====================================================================
        try: # ================================================================
            # =================================================================
            from   ostap.core.core     import rootException
            with rootException() : fun_loads ( fun_dumps ( objects ) )
            return True 
        except ( pickle.PicklingError   ,
                 pickle.UnpicklingError ,
                 AttributeError         ,
                 TypeError              ) : # =================================
            # =================================================================
            return False
        except Exception : # ===================================================
            # ==================================================================
            logger.always ( 'EXCEPTION/2' , exc_info = True ) 
            return False
        except : # ============================================================
            # ==================================================================
            logger.always ( 'EXCEPTION/3' , exc_info = True ) 
            return False
        
    # =========================================================================
    ## Check pickling of an object across another (sub) process
    def _pickles_process ( self , *objects          , 
                           fun_dump                 ,
                           fun_load                 , 
                           command = PICKLE_COMMAND ,
                           fast    = False          ) :
        """ Check pickling of an object across another (sub)process
        """
        objects = self.unpack ( *objects ) 
        if not objects               : return True
        objtypes = set ( type ( o ) for o in objects ) 
        if self.notgood ( *objtypes ) : return False        
        ## check if the object can be properly pickled/unpickled
        if self.known   ( *objtypes ) : return True
        # =====================================================================
        from   ostap.utils.cleanup import CleanUp 
        tmpfile = CleanUp.tempfile ( suffix = '.pkl')
        # =====================================================================
        try : # ===============================================================
            # =================================================================
            from   ostap.core.core     import rootException
            with rootException () : 
                with open ( tmpfile , 'wb' ) as f : fun_dump ( objects , f )
            # =================================================================
        except ( pickle.PicklingError ,
                 AttributeError       ,
                 TypeError            ,
                 OSError              ) : # ===================================
            # =================================================================            
            CleanUp.remove_file ( tmpfile ) 
            return False
        except Exception : # ===================================================
            # ==================================================================
            logger.always ( 'EXCEPTION/2.1' , ext_info = True ) 
            CleanUp.remove_file ( tmpfile ) 
            return False
        except : # ============================================================
            # =================================================================
            logger.always ( 'EXCEPTION/3.1' , exc_info = True ) 
            CleanUp.remove_file ( tmpfile ) 
            return False
            
        # =====================================================================
        if ( fast and fun_load ) or not command : # ===========================
            # =================================================================
            try : # ===========================================================
                # =============================================================
                with open ( tmpfile , 'rb' ) as f :
                    fun_load ( f )
                    return True                                        # RETURN
            except ( pickle.UnpicklingError , 
                     AttributeError         ,
                     TypeError              ,
                     OSError                ) : # =============================
                # =============================================================
                return False                                           # RETURN
            else  : # =========================================================
                # =============================================================
                CleanUp.remove_file ( tmpfile ) 
            
        # =====================================================================
        ## the check-command 
        cmd = command % tmpfile
        # =====================================================================
        import subprocess 
        # =====================================================================
        try : # ===============================================================
            # =================================================================
            error = subprocess.run ( [ sys.executable , '-c' , cmd ] , 
                                     check  = True               ,
                                     stdout = subprocess.DEVNULL ,
                                     stderr = subprocess.DEVNULL ,                                      
                                     shell  = False              ).returncode
            # =================================================================
            return True if not error else False                        # RETURN
        except subprocess.CalledProcessError as e : # =========================
            # =================================================================
            return False                                               # RETURN
        else : # ==============================================================
            # =================================================================
            CleanUp.remove_file ( tmpfile ) 

    # =========================================================================
    ## check if all object can be properly pickled/unpickled 
    def pickles ( self , *objects ) :
        """ Check if all objects can be properly pickled/unpickled
        """
        return self._pickles  ( *objects                 ,
                                fun_dumps = pickle.dumps ,
                                fun_loads = pickle.loads )
    
    # =========================================================================
    ## Check pickling of all objects across another (sub) process
    def pickles_process ( self , *objects , fast = False   ) :
        """ Check pickling of all objects across another (sub)process
        """
        return self._pickles_process ( *objects                  ,
                                       fun_dump = pickle.dump    ,
                                       fun_load = pickle.load    ,
                                       command  = PICKLE_COMMAND ,
                                       fast     = fast           )
    
    # =========================================================================
    ## check if all arguments are (un)pickle-able
    def pickles_all ( self , *args , **kwargs ) :
        """ Check if all arguments are (un)pickle-able
        """
        objects = self.unpack ( *args , **kwargs ) 
        return self.pickles ( *objects )
                              
    # =========================================================================
    ## Check pickling of all objects across another (sub) process
    def pickles_process_all ( self , *args , **kwargs ) :
        """ Check pickling of all objects across another (sub) process
        """
        objects = self.unpack ( *args , **kwargs ) 
        return self.pickles_process ( *objects  )
    
    # =========================================================================
    ## add new type into the list of "known-types"
    def add ( self , *ntypes ) :
        """ Add new type into th elist of "known-types
        """
        for ntype in ntypes :
            if ntype in self : continue 
            self.MORE_TYPES.append ( ntype )
        
    # =========================================================================
    ## add new type into the list of "non-picleable" types 
    def add_nonpickleable ( self , *ntypes ) :
        """ Add new type into the list of "non-pickleable" types 
        """
        for ntype in ntypes :        
            self.NONPICKLEABLE_TYPES.add ( ntype ) 
        
    # =========================================================================
    ## Format and return the pickling table 
    def pickling_table ( self , *args , **kwargs ) :
        """ Format and return the pickling table 
        """        
        from ostap.core.ostap_types import string_types 
        from ostap.utils.basic      import typename, prntrf,  loop_items
    
        rows  = [  ( 'Argument' , 'type' , 'value' , "pickles'" , 'pickles"' ) ]

        from ostap.logger.symbols import checked_yes, checked_no 
        good, bad = checked_yes, checked_no 
        for i , a in enumerate ( args ) :
            pickable = self.pickles         ( a )
            process  = self.pickles_process ( a )
            row      = '#%d' % i , typename ( a ) , prntrf ( a ) , \
                good if pickable else bad , \
                good if process  else bad
            rows.append ( row )

        if 'prefix' in kwargs : kwargs [ '*prefix' ] = kwargs.pop ( 'prefix' )
        if 'title'  in kwargs : kwargs [ '*title'  ] = kwargs.pop ( 'title'  ) 
        
        for key , v in loop_items  ( kwargs ) :
            pickable = self.pickles         ( v )
            process  = self.pickles_process ( v )            
            row      = '%s' % key , typename ( v ) , prntrf  ( v ) , \
                good if pickable else bad , \
                good if process  else bad
            rows.append ( row  )
            
        prefix = kwargs.get ( '*prefix' , '' )
        if not prefix or not isinstance ( prefix , string_types ) : prefix = ''
        title  = kwargs.get ( '*ttile' , '(Un)Pickle-able?' )
        if not title  or not isinstance ( title  , string_types ) : title = '(Un)Pickle-able?'
        
        import ostap.logger.table as T        
        return T.table ( rows , title = title , prefix = prefix , alignment = 'lwwcc' )

    # ==========================================================================
    ## get table of known pickleable types 
    def table_pickleable  ( self , title = 'Pickleable Types' , prefix = '' , style = '' ) :
        """ Get table of known pickleable types 
        """
        from ostap.utils.basic import typename 
        extra = set ( self.MORE_TYPES ) | set ( self.EXTRA_TYPES ) | set ( PICKLE_TYPES ) 
        rows  = [ ( '#' , 'Type' ) ]
        extra = sorted ( typename ( e ) for e in extra ) 
        for i , e in enumerate ( extra , start = 1 ) :
            row = '%d' % i , e  
            rows.append ( row ) 
        import ostap.logger.table as T
        return T.table ( rows , title = title , prefix = prefix, alignment = 'll' , style = style )

    # ==========================================================================
    ## get table of known pickleable types 
    def table_nonpickleable  ( self , title = 'Non-pickleable Types' , prefix = '' , style = '' ) :
        """ Get table of known non-pickleable types 
        """
        from ostap.utils.basic import typename 
        rows  = [ ( '#' , 'Type' ) ]
        extra = sorted  ( typename ( t ) for t in self.NONPICKLEABLE_TYPES ) 
        for i , e in enumerate ( extra , start = 1 ) :
            row = '%d' % i , e  
            rows.append ( row ) 
        import ostap.logger.table as T
        return T.table ( rows , title = title , prefix = prefix, alignment = 'll' , style = style )
    
# =============================================================================
import atexit
@atexit.register
def _report_ () :
    # =========================================================================

    return

"""
    checker = PickleChecker()
    
    if checker.MORE_TYPES or checker.EXTRA_TYPES :
        title = 'Known pickleable types'
        table = checker.table_pickleable ( title = title , prefix = '# ' ) 
        logger.info ( '%s:\n%s' % ( title , table ) )
        
    if checker.NONPICKLEABLE_TYPES :
        title = 'Non-pickleable types'
        table = checker.table_nonpickleable ( title = title , prefix = '# ' ) 
        logger.info ( '%s:\n%s' % ( title , table ) ) 
"""

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
