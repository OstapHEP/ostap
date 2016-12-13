#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file roo_file.py
#  Module with decoration of TFile objects for efficient use in python
#  It provides TFile (well, actually any TDirectory) with python-like protocol
#  @code
#  rfile = ...
#
#  obj   = rfile['A/B/C/myobj']     ## READ  object from the file/directory
#  rfile['A/B/C/myobj2'] = object2  ## WRITE object to the file/directory 
#
#  obj  = rfile.A.B.C.myobj        ## another way to access to the object in file
#
#  obj  = rfile.get ( 'A/B/C/q' , None )   ## one more way to get object 
#
#  for obj in rfile : print obj     ## loop over all objects in file
#  for key,obj in rfile.iteritems() : print key, obj             ## another loop
#  for key,obj in rfile.iteritems( ROOT.TH1 ) : print key, obj   ## advanced loop
#  for key in rfile.keys()     : print k   ## get all keys and loop over them 
#  for key in rfile.iterkeys() : print k   ## loop over all keys in the file 
#
#  del  rfile['A/B']                       ## delete the object from the file
#
#  if 'A/MyHisto' in rfile          : print 'OK!' ## check presence of the key
#  if rfile.has_key ( 'A/MyHisto' ) : print 'OK!' ## check presence of the key
#
#  with ROOT.TFile('aa.root') as rfile : rfile.ls() ## context manager protocol
#  @endcode 
#
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
"""Decoration of TFile objects for efficient use in python

It provides TFile (well, actually any TDirectory) with python-like protocol

>>> rfile = ...

>>> obj   = rfile['A/B/C/myobj']     ## READ  object form the file/directory
>>> rfile['A/B/C/myobj2'] = object2  ## WRITE object to the file/directory 

>>> obj  = rfile.A.B.C.myobj              ## another way to access to the object
>>> obj  = rfile.get ( 'A/B/C/q' , None ) ## one more way to get object 

>>> for obj in rfile : print obj            ## loop over all objects in file
>>> for key,obj in rfile.iteritems() : print key, obj             ## another loop
>>> for key,obj in rfile.iteritems( ROOT.TH1 ) : print key, obj   ## advanced loop
>>> for k in rfile.keys()     : print k   ## get all keys and loop over them 
>>> for k in rfile.iterkeys() : print k   ## loop over all keys in the file

>>> del  rfile['A/B']                       ## delete the object from the file
>>> rfile.rm ( 'A/B' )                      ## delete the object from the file

>>> if 'A/MyHisto' in rfile          : print 'OK!' ## check presence of the key
>>> if rfile.has_key ( 'A/MyHisto' ) : print 'OK!' ## check presence of the key

>>> with ROOT.TFile('aa.root') as rfile : rfile.ls() ## context manager protocol
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'ROOTCWD'      , ## context manager to keep Current Directory
    'open_mode'    , ## decode open-mode for ROOT-files
    'open'         , ## just for completness 
    ) 
# =============================================================================
import ROOT, os , cppyy              ## attention here!!
cpp = cppyy.gbl
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.io.root_file' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'Some useful decorations for TFile objects')
# ==============================================================================
## context manager to preserve current directory (rather confusing stuff in ROOT)
from ostap.core.core import ROOTCWD
# ===============================================================================
## write the (T)object to ROOT-file/directory
#  @code
#  histo1 = ...
#  histo2 = ...
#  rfile['MyHisto'      ] = histo1 
#  rfile['MyDir/MyHisto'] = histo2
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@iep.ru
#  @date 2015-07-30
def _rd_setitem_ ( rdir , name , tobj ) :
    """Write the object to ROOT-file/directory
    >>> rfile['myhisto'] = h1
    >>> rfile['MyDir/MyHisto'] = h2
    """
    ##
    if   not rdir              :
        raise IOError("TDirectory is invalid")
    elif not rdir.IsWritable() :
        raise IOError("TDirectory '%s' is not writable" % rdir.GetPath() ) 
    ##
    with ROOTCWD() :
        ##
        dirname , sep , fname = name.partition('/')
        ##
        while sep :
            rdir.cd() 
            subdir = rdir.GetDirectory(dirname)
            if not subdir :
                rdir.cd() 
                subdir = rdir.mkdir ( dirname  , rdir.GetPath() + '/' + dirname )
                rdir.cd() 
                subdir = rdir.GetDirectory( dirname )
                if not subdir :
                    raise KeyError("TDirectory[]: can't create '%s%s'" % ( rdir.GetPath() , dirname ) )
            rdir   = subdir
            dirname, sep , fname = fname.partition ('/') 

        rdir.cd()
        return  rdir.WriteTObject( tobj , dirname , 'WriteDelete' )
        ## DO NOT DELETE OBJECT IN ROOT 
        ## ROOT.SetOwnership( tobj , False )
        ## return val 

# =============================================================================
## add named object into directory
#  @code
#  rdir   = ...    ## some writable ROOT directory
#  tnamed = ...    ## some named object
#  tnamed >> rdir  ## write it 
#  @endcore 
#  @author Vanya BELYAEV Ivan.Belyaev@iep.ru
#  @date 2016-06-03
def _rd_rrshift_ ( rdir , tnamed ) :
    """Add named object into directory
    
    rdir   = ...    ## some writable ROOT directory
    tnamed = ...    ## some named object
    tnamed >> rdir  ## write it
    """
    
    if   hasattr ( tnamed , 'GetName' ) : name = tnamed.GetName ()
    elif hasattr ( tnamed , 'name'    ) : name = tnamed.name    ()
    else : name = tnamed.__class__.__name__
    
    return _rd_setitem_ ( rdir , name , tnamed )


# =============================================================================
## add named object into directory
#  @code
#  rdir   = ...    ## some writable ROOT directory
#  tnamed = ...    ## some named object
#  tnamed >> rdir  ## write it 
#  @endcore 
#  @author Vanya BELYAEV Ivan.Belyaev@iep.ru
#  @date 2016-06-03
def _tnamed_rshift_ ( tnamed , rdir ) :
    """Add named object into directory
    
    rdir   = ...    ## some writable ROOT directory
    tnamed = ...    ## some named object
    tnamed >> rdir  ## write it
    """
    return _rd_setitem_ ( rdir , tnamed.GetName() , tnamed )

                 
# ===============================================================================
## get the object from ROOT-file/directory
#  @code
#  h1 = rfile['MyHisto'      ] 
#  h2 = rfile['MyDir/MyHisto'] 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@iep.ru
#  @date 2015-07-30
def _rd_getitem_ ( rdir , name ) :
    """Get the object from ROOT-file/directory
    >>> t = f['A/tup']
    >>> h = f['histo']
    """
    ##
    if not rdir : raise IOError("TDirectory is invalid")
    ##
    with ROOTCWD() :
        ##
        dirname , sep , fname = name.partition('/')
        ##
        while sep :
            rdir.cd()
            subdir = rdir.GetDirectory(dirname)
            if not subdir :
                raise KeyError("TDirectory[]: unknown directory name '%s%s'" % (rdir.GetPath() , dirname  ) ) 
            rdir = subdir
            dirname, sep , fname = fname.partition ('/') 

        rdir.cd()
        obj = rdir.Get ( dirname ) 
        if not obj : raise KeyError ("TDirectory[]: unknown object '%s%s'" % ( rdir.GetPath(), dirname ) )
        return obj 

# ===============================================================================
## get the object from ROOT-file/directory 
#  @code
#  h1 = rfile.MyHisto
#  h2 = rfile.MyDir.MyHisto
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@iep.ru
#  @date 2015-07-30
def _rd_getattr_ ( rdir , name ) :
    """
    Get the object from ROOT-file
    >>> t = f.A.tup
    >>> h = f.myhisto
    """
    ##
    if not rdir : return None
    ##
    with ROOTCWD() :
        ##
        dirname , sep , fname = name.partition('/')
        if sep : raise AttributeError('TDirectory[]: invalid attribute %s.%s' % rdir.GetPath() , dirname )

        rdir.cd()
        obj = rdir.Get ( dirname ) 
        if not obj : raise AttributeError ( 'TDirectory[]: unknown attribute %s.%s' % ( rdir.GetPath(), dirname ) )
        return obj 
    
# ===============================================================================
## check the existence of key in the root-file/directory
#  @code
#  if not 'myhisto'   in rfile : break 
#  if 'MyDir/MyHisto' in rfile : print 'OK!'
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@iep.ru
#  @date 2015-07-30
def _rd_contains_ ( rdir , name ) :
    """Check the existence of key in the root-file/directory
    >>> if not 'myhisto'  in rfile : break 
    >>> if 'MyDir/MyHisto' in rfile : print 'OK!'
    """
    ##
    if not rdir : return False 
    ##
    with ROOTCWD() :
        ##
        dirname , sep , fname = name.partition('/')
        ##
        while sep :
            rdir.cd()
            subdir = rdir.GetDirectory(dirname)
            if not subdir : return False                      ## RETURN 
            rdir = subdir
            dirname, sep , fname = fname.partition ('/') 

        rdir.cd()
        return rdir.FindKey  ( dirname )                       ## RETURN 

# =============================================================================
## delete object from ROOT file/directory
#  @code
#  del rfile['h1']
#  del rfile['mydir/h1']
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@iep.ru
#  @date 2015-07-30
def _rd_delitem_ ( rdir , name , cycle=';*') :
    """Delete an object from ROOT file/directory
    >>> del rfile['h1']
    >>> del rfile['mydir/h1']
    """
    ##
    if   not rdir              :
        raise IOError("TDirectory is invalid")
    elif not rdir.IsWritable() :
        raise IOError("TDirectory '%s' is not writable" % rdir.GetPath() ) 
    ##
    with ROOTCWD () :
        ##
        dirname , sep , fname = name.partition('/')
        ##
        while sep :
            rdir.cd()
            subdir = rdir.GetDirectory(dirname)
            if not subdir :
                raise KeyError("TDirectory, can't delete %s" % name ) 
            rdir = subdir
            dirname, sep , fname = fname.partition ('/') 

        rdir.cd()
        if not rdir.FindKey( dirname ) :
            raise KeyError("TDirectory, can't delete %s" % dirname ) 
        return rdir.Delete( dirname + cycle ) 

# =============================================================================a
## get all keys from ROOT file/directory
#  @code
#  keys = rfile.keys()
#  keys = rfile.keys( recursive = False )
#  keys = rfile.keys( no_dir    = False )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@iep.ru
#  @date 2015-07-30
def _rd_keys_ ( rdir , recursive = True , no_dir = True ) :
    """Get all keys from ROOT file/directory
    >>> keys = rfile.keys() 
    """
    _res = []
    if not rdir : return _res 
    ##
    with ROOTCWD() :
        ##
        rdir.cd() 
        _lst = rdir.GetListOfKeys()
        for i in _lst :
            inam = i.GetName()

            idir = rdir.GetDirectory ( inam )            
            if not idir or not no_dir : _res.append ( inam )
            
            if recursive and idir and not idir is rdir :
                ikeys = _rd_keys_ ( idir , recursive , no_dir )
                for k in ikeys : _res.append ( inam + '/' + k )
                
        return _res
    
# =============================================================================a
## Iterate over the content of ROOT file/directory 
#  @code
#  for key,obj in rfile.iteritems() : print key,obj
#  @endcode
#  The scan can be limited only by objects of certain type, e.g. histograms: 
#  @code
#  for key,hist in rfile.iteritems( ROOT.TH1 ) : print key,hist
#  @endcode
#  More complicated rules can be specified:
#  @code 
#  for key,obj in rfile.iteritems ( lambda k,o : k[0]=='M' ) : print key,obj
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@iep.ru
#  @date 2015-07-30
def _rd_iteritems_ ( rdir , fun = lambda k,t,o : True , recursive = True , no_dir = True  ) :
    """Iterate over the content of ROOT directory/file:
    >>> for key,obj  in rfile.iteritems()           : print key , obj
    >>> for key,hist in rfile.iteritems( ROOT.TH1 ) : print key , hist
    >>> for key,obj  in rfile.iteritems( lambda name,tkey,obj : name[0]=='M' ) : print key,obj
    """
    ##
    if isinstance ( fun , type ) and issubclass ( fun , ( ROOT.TObject, cpp.TObject) ) : 
        tobj = fun 
        fun  = lambda k,t,o : isinstance ( o , tobj )
    ##
    with ROOTCWD() :
        ##
        rdir.cd() 
        _lst = rdir.GetListOfKeys()
        for i in _lst :
            inam   = i.GetName()
            idir   = rdir.GetDirectory ( inam ) 
            if not idir or not no_dir : 
                obj = rdir.Get ( inam )
                if fun ( inam , i , obj ) : yield inam , obj
            if recursive and idir and not idir is rdir :
                for k, o in _rd_iteritems_ ( idir , fun , recursive , no_dir ) :
                    yield k,o
                    
# =============================================================================a
## Iterate over the keys in ROOT file/directory 
#  @code
#  for k in rfile.iterkeys() : print k
#  @endcode
#  The scan can be limited only by objects of certain type, e.g. histograms: 
#  @code
#  for k in rfile.iterkeys ( ROOT.TH1 ) : print k
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@iep.ru
#  @date 2015-07-30
def _rd_iterkeys_ ( rdir , typ = None , recursive = True , no_dir = True ) :
    """Iterate over the keys in ROOT file/directory 
    >>> for key,obj  in rfile.iteritems()           : print key , obj
    >>> for key,hist in rfile.iteritems( ROOT.TH1 ) : print key , hist
    """
    ##
    with ROOTCWD() :
        ##
        rdir.cd() 
        _lst = rdir.GetListOfKeys()
        for i in _lst :
            inam   = i.GetName()
            idir   = rdir.GetDirectory ( inam ) 
            if not idir or not no_dir : 
                if typ is None  or isinstance ( rdir.Get ( inam ) , typ ) : yield inam 
            if recursive and idir  and not idir is rdir :
                for k in _rd_iterkeys_ ( idir , typ , recursive , no_dir ) :
                    yield k

# =============================================================================a
## Iterate over the content of ROOT file/directory 
#  @code
#  for  obj  in rfile.itervalues() : print obj
#  @endcode
#  The scan can be limited only by objcets of certain type, e.g. histograms: 
#  @code
#  for  hist in rfile.itervalues( ROOT.TH1 ) : print hist 
#  @endcode
#  More complicated rules can be specified:
#  @code 
#  for  obj in rfile.itervalues( lambda k,o : k[0]=='M' ) : print obj
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@iep.ru
#  @date 2015-07-30
def _rd_itervalues_ ( rdir , fun = lambda k,t,o : True , recursive = True , no_dir = True ) :
    """Iterate over the content of ROOT directory/file:
    >>> for obj  in rfile.itervalues ()             : print obj
    >>> for hist in rfile.itervalues ( ROOT.TH1 )   : print hist
    >>> for obj  in rfile.itervalues ( lambda k,t,o : k[0]=='M' ) : print obj
    """
    ##
    with ROOTCWD() :
        ##
        rdir.cd()
        for key , obj in _rd_iteritems_ ( rdir , fun , recursive , no_dir ) :
            yield obj 

# =============================================================================a
## Iterate (recursive) over the content in ROOT file/directory
#  @code
#  for obj in rfile : print obj
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@iep.ru
#  @date 2015-07-30
def _rd_iter_ ( rdir ) :
    """Iterate (recursive) over the content in ROOT file/directory
    >>> for obj  in rfile : print obj
    """
    ##
    with ROOTCWD() :
        ##
        rdir.cd()
        for obj in _rd_itervalues_ ( rdir , recursive = True ) :
            yield obj

# =============================================================================
## delete object from ROOT file/directory
#  @code
#  rfile.rm('h1')
#  rfile.rm('mydir/h1')
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@iep.ru
#  @date 2015-07-30
def _rd_rm_ ( rdir , name , cycle=';*') :
    """Delete an object from ROOT file/directory
    >>> rfile.rm( 'h1')
    >>> rfile.rm('mydir/h1')
    """
    ##
    try :
        return _rd_delitem_ ( rdir , name , cycle )
    except KeyError :
        return
            
# ===============================================================================
## get the object from ROOT-file 
#  @code
#  h1 = rfile.get( 'MyHisto' , default_value ) 
#  h2 = rfile.get( 'MyDir/MyHisto' )
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@iep.ru
#  @date 2015-07-30
def _rd_get_ ( rdir , name , default = None ) :
    """
    Get the object from ROOT-file
    >>> h1 = rfile.get( 'MyHisto' , default_value ) 
    >>> h2 = rfile.get( 'MyDir/MyHisto')
    """
    ##
    if not rdir : return default 
    ##
    try :
        return _rd_getitem_ ( rdir , name )
    except KeyError :
        return default
    
# =============================================================================
## use ROOT-file with context-manager
#  @code
#  with ROOT.TFile('ququ.root') as f : f.ls() 
#  @endcode
def _rf_enter_ ( self      ) :
    """Use ROOT-file with the context manager
    >>> with ROOT.TFile('ququ.root') as f : f.ls() 
    """
    return self

# =============================================================================
## use ROOT-file with context-manager 
#  @code
#  with ROOT.TFile('ququ.root') as f : f.ls() 
#  @endcode
def _rf_exit_  ( self , *_ ) :
    """Use ROOT-file with the context manager    
    >>> with ROOT.TFile('ququ.root') as f : f.ls()
    """
    try :
        if self and self.IsOpen() : self.Close()
    except:
        pass
            

# =============================================================================
## the basic protocol:

ROOT.TDirectory.__getitem__  = _rd_getitem_ 
ROOT.TDirectory.__setitem__  = _rd_setitem_
ROOT.TDirectory.__contains__ = _rd_contains_
ROOT.TDirectory.__getattr__  = _rd_getattr_
ROOT.TDirectory.__delitem__  = _rd_delitem_
ROOT.TDirectory.__iter__     = _rd_iter_

# =============================================================================
## the extended protocol

ROOT.TDirectory.get          = _rd_get_
ROOT.TDirectory.keys         = _rd_keys_
ROOT.TDirectory.has_key      = _rd_contains_
ROOT.TDirectory.iteritems    = _rd_iteritems_
ROOT.TDirectory.iterkeys     = _rd_iterkeys_
ROOT.TDirectory.itervalues   = _rd_itervalues_

# =============================================================================
## some extra stuff 

ROOT.TDirectory.__rrshift__  = _rd_rrshift_
ROOT.TNamed    .__rshift__   = _tnamed_rshift_ 

ROOT.TDirectory.rm           = _rd_rm_

if hasattr ( ROOT.TFile , '__enter__' ) and hasattr ( ROOT.TFile , '__exit__' ) : pass
else :
    ROOT.TFile.__enter__ = _rf_enter_
    ROOT.TFile.__exit__  = _rf_exit_

    
# =============================================================================
_modes = {
    ##
    'n'        : 'NEW'      ,
    'c'        : 'RECREATE' ,
    'r'        : 'READ'     ,
    'u'        : 'UPDATE'   ,
    'w'        : 'WRITE'    ,
    'a'        : 'UPDATE'   ,
    ##
    '+'        : 'UPDATE'   ,
    'rw'       : 'UPDATE'   ,
    'w+'       : 'UPDATE'   ,
    'new'      : 'NEW'      ,
    'create'   : 'RECREATE' ,
    'recreate' : 'RECREATE' ,
    'read'     : 'READ'     ,
    'write'    : 'WRITE'    ,
    'update'   : 'UPDATE'   ,
    'append'   : 'APPEND'   ,
    ## 
    }

# =============================================================================
## define open modes for ROOT-file 
def open_mode ( mode ) :
    """ Define open-mode for ROOT-file
    >>> m = open_mode ( 'n' )
    """
    return _modes.get ( mode.lower() , mode.upper() )

# =============================================================================
## create ROOT.TFile without making it a current working directory 
#  @code
#  print ROOT.gROOT.CurrentDirectory()
#  f = ROOT.TFile('test_file.root','recreate')
#  print ROOT.gROOT.CurrentDirectory()
#  @endcode
#  @attention  IOError exception is raised for invalid file/open_mode 
def _rf_new_init_ ( rfile , fname , mode = '' , *args ) :
    """Open/create ROOT-file without making it a current working directory
    >>> print ROOT.gROOT.CurrentDirectory()
    >>> f = ROOT.TFile('test_file.root','recreate')
    >>> print ROOT.gROOT.CurrentDirectory()
    Attention:  IOError exception is raised for invalid file/open_mode
    """
    with ROOTCWD() :
        logger.debug ("Open  ROOT file %s/'%s'" % ( fname , mode ) )
        result = rfile._old_init_ ( fname , open_mode ( mode ) , *args )
        if not rfile or not rfile.IsOpen() :
            ## logger.error ( "Can't open ROOT file %s/'%s'" % ( fname , mode ) )
            raise IOError   ( "Can't open ROOT file %s/'%s'" % ( fname , mode ) )
        return result
    
# =============================================================================
## create ROOT.TFile without making it a current working directory 
#  @code
#  print ROOT.gROOT.CurrentDirectory()
#  f = ROOT.TFile.Open('test_file.root','recreate')
#  print ROOT.gROOT.CurrentDirectory()
#  @endcode
#  @attention  No exceptions are raised for invalid file/open_mode 
def _rf_new_open_ ( fname , mode = '' , *args ) :
    """Open/create ROOT-file without making it a current working directory
    >>> print ROOT.gROOT.CurrentDirectory()
    >>> f = ROOT.TFile.Open('test_file.root','recreate')
    >>> print ROOT.gROOT.CurrentDirectory()
    ATTENTION: No exceptions are raised for invalid file/open_mode 
    """
    with ROOTCWD() :
        logger.debug ( "Open  ROOT file %s/'%s'" % ( fname , mode ) )
        fopen = ROOT.TFile._old_open_ ( fname , open_mode ( mode ) , *args )
        if not fopen or not fopen.IsOpen() :
            logger.error ( "Can't open ROOT file %s/'%s'" % ( fname , mode ) )
            ## raise IOError   ( "Can't open ROOT file %s/'%s'" % ( fname , mode ) )
        #
        return fopen
        
# =============================================================================
## create ROOT.TFile without making it a current working directory 
#  @code
#  print ROOT.gROOT.CurrentDirectory()
#  f = ROOT.TFile.Open('test_file.root','recreate')
#  print ROOT.gROOT.CurrentDirectory()
#  @endcode
def _rf_new_close_ ( rfile , options = '' ) :
    """
    Open/create ROOT-file without making it a current working directory
    >>> print ROOT.gROOT.CurrentDirectory()
    >>> f = ROOT.TFile.Open('test_file.root','recreate')
    >>> print ROOT.gROOT.CurrentDirectory()
    """
    if rfile and rfile.IsOpen() :
        with ROOTCWD() :
            ##rfile .cd()
            logger.debug ( "Close ROOT file %s" % rfile.GetName() ) 
            rfile ._old_close_ ( options )
            
# =============================================================================
## another name, just for convinince
open = _rf_new_open_

# =============================================================================
if hasattr ( ROOT.TFile , '_new_init_' ) and hasattr ( ROOT.TFile , '_old_init_' ) : pass
else :
    
    _rf_new_init_.__doc__  += '\n' + ROOT.TFile.__init__.__doc__
    
    ROOT.TFile._old_init_   = ROOT.TFile.__init__
    ROOT.TFile._new_init_   = _rf_new_init_ 
    ROOT.TFile.__init__     = _rf_new_init_ 
    
# =============================================================================
if hasattr ( ROOT.TFile , '_new_open_' ) and hasattr ( ROOT.TFile , '_old_open_' ) : pass
else :

    _rf_new_open_.__doc__  += '\n' + ROOT.TFile.Open.__doc__
    _rf_new_open_           = staticmethod( _rf_new_open_ )
    
    ROOT.TFile._old_open_   = ROOT.TFile.Open
    ROOT.TFile._new_open_   = _rf_new_open_ 
    ROOT.TFile.Open         = _rf_new_open_ 
    ROOT.TFile.open         = _rf_new_open_ 

# =============================================================================
if hasattr ( ROOT.TFile , '_new_close_' ) and hasattr ( ROOT.TFile , '_old_close_' ) : pass
else :

    _rf_new_close_.__doc__  += '\n' + ROOT.TFile.Close.__doc__
    
    ROOT.TFile._old_close_   = ROOT.TFile.Close
    ROOT.TFile._new_close_   = _rf_new_close_ 
    ROOT.TFile.Close         = _rf_new_close_
    
if not hasattr ( ROOT.TFile , 'close' ) :
    ROOT.TFile.close = ROOT.TFile._new_close_

    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

    
# =============================================================================
# The END 
# =============================================================================
