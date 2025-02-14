#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file pickling.py
# Helper module to define pickling for various databases 
# @author Vanya BELYAEV Ivan.Belyaev@cern.ch
# @date   2010-04-30
# =============================================================================
""" Helper module to define pickling for various databases
"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2010-04-30"
__version__ = "$Revision$" 
# =============================================================================
__all__ = (
    'DEFAULT_PROTOCOL' ,
    'HIGHEST_PROTOCOL' ,
    'PROTOCOL'         ,
    'Pickler'          , 
    'Unpickler'        ,
    'BytesIO'          ,
    'dumps'            ,
    'loads'            ,
    'PickleChecker'    , ## check pickle-ability of objects 
    )
# =============================================================================
from   pickle         import ( Pickler, Unpickler, 
                               DEFAULT_PROTOCOL, HIGHEST_PROTOCOL,
                               dumps, loads, dump , load , 
                               PicklingError, UnpicklingError )
from   io             import  BytesIO 
import sys, array 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.io.pickling' )
else                      : logger = getLogger ( __name__               )
# =============================================================================
## the basic pickle-types that for sure always can be pickled/unpickled 
PICKLE_TYPES = type ( None ) , bool, int, float, str, bytes, bytearray, array.array 
# =============================================================================
PICKLE_COMMAND = """import sys, pickle
with open('%s','rb') as f : pickle.load ( f )"""
# =============================================================================
## @class PickleChecker
#  Check if the object/type can be pickled/unpickled
class PickleChecker ( object ) :
    """ Check if the object/type can be pickled/uunpickeld
    """
    MORE_TYPES          = set ()
    EXTRA_TYPES         = set ()
    NONPICKLEABLE_TYPES = set () 
    # ========================================================================
    ## if the object type is already knowns 
    def known ( self , *objtypes ) :
        return all ( ( o in PICKLE_TYPES    ) or
                     ( o in self.MORE_TYPES ) for o in objtypes ) 
    # ========================================================================
    ## type is knows for (un)pickle problems 
    def notgood ( self , *objtypes ) :
        return any ( ( o in self.NONPICKLEABLE_TYPES ) for o in objtypes )
    # =========================================================================
    def _pickles ( self , *objects , fun_dumps , fun_loads ) :
        # =====================================================================
        if not objects               : return True
        if self.notgood ( *objects ) : return False 
        ## check if the object can be properly pickled/unpickled
        if self.known ( *( type ( o ) for o in objects ) ) : return True
        from   ostap.core.core import rootException 
        # =====================================================================
        print ( 'EXCEPTION/0' ) 
        try: # ================================================================
            # =================================================================
            ## return fun_loads ( dumps ( objects ) ) == objects # ===============
            with rootException() : fun_loads ( fun_dumps ( objects ) )
            return True 
        except ( PicklingError, UnpicklingError, AttributeError, TypeError ) :
            # =================================================================
            print ( 'EXCEPTION/1' ) 
            CU.CleanUp.remove_file ( tmpfile ) 
            return False
        except Exception : # ===================================================
            # ==================================================================
            print ( 'EXCEPTION/2' ) 
            CU.CleanUp.remove_file ( tmpfile ) 
            return False
        except : # ============================================================
            # ==================================================================
            print ( 'EXCEPTION/3' ) 
            CU.CleanUp.remove_file ( tmpfile ) 
            return False
        
    # =========================================================================
    ## check if the type is 'known'
    def __contains__ ( self , objtype ) :
        return objtype in PICKLE_TYPES or \
            objtype in self.MORE_TYPES or \
            objtype in self.EXTRA_TYPES
    
    # =========================================================================
    ## Check pickling of an object across another (sub) process
    def _pickles_process ( self , *objects          ,
                           fun_dump                 ,
                           command = PICKLE_COMMAND ,
                           fast    = False          ) :
        """ Check pickling of an object across another (sub)process
        """
        if not objects               : return True
        if self.notgood ( *objects ) : return False 
        # =====================================================================
        ## check if the object can be properly pickled/unpickled
        if self.known ( *(type ( o ) for o in objects ) ) : return True
        # =====================================================================
        import ostap.utils.cleanup as     CU
        from   ostap.core.core     import rootException 
        tmpfile = CU.CleanUp.tempfile ( suffix = '.pkl')
        # =====================================================================
        print ( 'EXCEPTION/0.1' ) 
        try : # ===============================================================
            # =================================================================
            with rootException () : 
                with open ( tmpfile , 'wb' ) as f : fun_dump ( objects , f )
            # =================================================================
        except ( PicklingError  ,
                 AttributeError ,
                 TypeError      ,
                 OSError        ) : # ===
            # =================================================================            
            print ( 'EXCEPTION/1.1' ) 
            CU.CleanUp.remove_file ( tmpfile ) 
            return False
        except Exception : # ===================================================
            # ==================================================================
            print ( 'EXCEPTION/2.1' ) 
            CU.CleanUp.remove_file ( tmpfile ) 
            return False
        except : # ============================================================
            # ==================================================================
            print ( 'EXCEPTION/3.1' ) 
            CU.CleanUp.remove_file ( tmpfile ) 
            return False
            
        # =====================================================================
        if fast or not command : # ============================================
            # =================================================================
            try : # ===========================================================
                # =============================================================
                with open ( tmpfile , 'rb' ) as f :
                    ## return load ( f ) == objects                    # RETURN
                    load ( f )
                    return True                                        # RETURN
            except ( UnpicklingError , 
                     AttributeError  ,
                     TypeError       ,
                     OSError         ) : # ====================================
                # =============================================================
                return False                                           # RETURN
            else  : # =========================================================
                # =============================================================
                CU.CleanUp.remove_file ( tmpfile ) 
            
        # =====================================================================
        ## the check-command 
        cmd = command % tmpfile
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
            CU.CleanUp.remove_file ( tmpfile ) 

    # =========================================================================
    ## check if the object can be properly pickled/unpickled 
    def pickles ( self , *objects ) :
        """ Check of the object can be properly pickled/unpickled
        """
        return self._pickles  ( *objects          ,
                                fun_dumps = dumps ,
                                fun_loads = loads )
    
    # =========================================================================
    ## Check pickling of an object across another (sub) process
    def pickles_process ( self , *objects , fast = False   ) :
        """ Check pickling of an object across another (sub)process
        """
        return self._pickles_process ( *objects                  ,
                                       fun_dump = dump           ,
                                       command  = PICKLE_COMMAND ,
                                       fast     = fast           )
    
    # =========================================================================
    ## check if all arguments are (un)pickle-able
    def pickles_all ( self , *args , **kwargs ) :
        """ Check if all arguments are (un)pickle-able
        """
        objects = args + tuple ( v for v in kwargs.values() )
        return self.pickles ( *objects )
                              
    # =========================================================================
    ## Check pickling of all objects across another (sub) process
    def pickles_process_all ( self , *args , **kwargs ) :
        """ Check pickling of all objects across another (sub) process
        """
        objects = args + tuple ( v for v in kwargs.values() )
        return self.pickles_process ( *objects  )
    
    # =========================================================================
    ## add new type into the list of "known-types"
    def add ( self , ntype ) :
        """ Add new type into th elist of "known-types
        """
        if ntype in self : return 
        self.MORE_TYPES.add ( ntype )
        
    # =========================================================================
    ## add new type into the list of "non-picleable" types 
    def add_nonpickleable ( self , ntype ) :
        """ Add new type into the list of "non-pickleable" types 
        """
        self.NONPICKLEABLE_TYPES.add ( ntype ) 
        
    # =========================================================================
    ## Format and return the pickling table 
    def pickling_table ( self , *args , **kwargs ) :
        """ Format and return the pickling table 
        """        
        from ostap.core.ostap_types import string_types 
        from ostap.utils.basic      import typename, prntrf,  loop_items
    
        rows  = [  ( 'Argument' , 'type' , 'value' , 'pickles' , "pickles'" ) ]
        
        good, bad = '\u2714' , '\u2715' 
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
    
# =============================================================================
## helper function to get the protocol 
def get_protocol ( p ) :
    """ helper function to get the protocol"""    
    if   p.lower () in ( 'default' , 'def'  ) : return DEFAULT_PROTOCOL            
    elif p.lower () in ( 'highest' , 'high' ) : return HIGHEST_PROTOCOL
    elif p.lower () in ( 'compat'  , 'compatible' , 'backward_compatible' ) :
        return min ( 2 , HIGHEST_PROTOCOL )  
    else :
        return int ( p )
# =============================================================================    
## pickle protocol to be used 
PROTOCOL = None
#  (1) pickup protocol from the environment variable
# =============================================================================
try : # =======================================================================
    # =========================================================================
    from ostap.utils.basic import get_env, OSTAP_PROTOCOL  
    pe = get_env ( OSTAP_PROTOCOL  , '' )
    pp = get_protocol   ( pe )
    if pp < 0 : pp = HIGHEST_PROTOCOL 
    if 0  <= pp <= HIGHEST_PROTOCOL :
        PROTOCOL = pp
        logger.debug ( "Protocol %s is picked from 'OSTAP_PROTOCOL=%s' environment" % ( PROTOCOL , pe ) )
    # =========================================================================
except : # ====================================================================
    # =========================================================================
    pass
# =============================================================================
#  (2) take protocol from the configuration files
# =============================================================================
if PROTOCOL is None : # =======================================================
    # =========================================================================
    try : # ===================================================================
        # =====================================================================
        import ostap.core.config as OCC
        pe = OCC.general.get ( 'Protocol' , '' )
        pp = get_protocol ( pe )
        if pp < 0 : pp = HIGHEST_PROTOCOL 
        if 0  <= pp <= HIGHEST_PROTOCOL :
            PROTOCOL = pp 
            logger.debug ( "Protocol %s is picked from 'General: protocol=%s' section" %  ( PROTOCOL , pe ) )
        # =====================================================================
    except : # ================================================================
        # =====================================================================
        pass 
# =============================================================================
#  (3) use default protocol
# =============================================================================
if PROTOCOL is None : # =======================================================
    # =========================================================================
    PROTOCOL = DEFAULT_PROTOCOL 
    logger.debug ( "Default protocol %s is used" % PROTOCOL  )

# =============================================================================
import atexit
@atexit.register
def _report_ () :
    # =========================================================================
    if PickleChecker.MORE_TYPES or PickleChecker.EXTRA_TYPES :
        from ostap.utils.basic import typename 
        extra = PickleChecker.MORE_TYPES | PickleChecker.EXTRA_TYPES | set ( PICKLE_TYPES ) 
        rows = [ ( '#' , 'Type' ) ]
        extra = sorted ( typename ( e ) for e in extra ) 
        for i , e in enumerate ( extra , start = 1 ) :
            row = '%d' % i , e  
            rows.append ( row ) 
        import ostap.logger.table as T
        title = 'Pickle types' 
        logger.info ( '%s:\n%s' % ( title , T.table ( rows , title = title , prefix = '# ' , alignment = 'll' ) ) ) 

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    logger.info ( 'Pickling protocol: %s' % PROTOCOL )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
