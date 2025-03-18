#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/io/root_file.py
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
    'REOPEN'       , ## context manager to <code>ROOT.TFileReOpen('UPDATE')</code>
    'RootFiles'    , ## utility to collect an ddeal wth ROOT files 
    ) 
# =============================================================================
from   ostap.core.core   import ROOTCWD, valid_pointer
from   ostap.io.files    import Files
import ROOT, sys, os   
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.io.root_file' )
else                       : logger = getLogger( __name__ )
# =============================================================================
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
    """ Write the object to ROOT-file/directory
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
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@iep.ru
#  @date 2016-06-03
def _rd_rrshift_ ( rdir , tnamed ) :
    """ Add named object into directory
    
    >>> rdir   = ...    ## some writable ROOT directory
    >>> tnamed = ...    ## some named object
    >>> tnamed >> rdir  ## write it
    """
    
    if   hasattr ( tnamed , 'GetName' ) : name = tnamed.GetName ()
    elif hasattr ( tnamed , 'name'    ) : name = tnamed.name    ()
    else                                : name = typename ( tnamed ) 
    
    return _rd_setitem_ ( rdir , name , tnamed )


# =============================================================================
## add named object into directory
#  @code
#  rdir   = ...    ## some writable ROOT directory
#  tnamed = ...    ## some named object
#  tnamed >> rdir  ## write it 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@iep.ru
#  @date 2016-06-03
def _tnamed_rshift_ ( tnamed , rdir ) :
    """ Add named object into directory
    
    >>>  rdir   = ...    ## some writable ROOT directory
    >>> tnamed = ...    ## some named object
    >>> tnamed >> rdir  ## write it
    """
    return _rd_setitem_ ( rdir , tnamed.GetName() , tnamed )

                 
# ===============================================================================
## get the (key,object) pair from ROOT-file/directory
#  @code
#  t = f.get_key_object ('A/tup') 
#  h = f.get_key_object ('histo') 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@iep.ru
#  @date 2015-07-30
def _rd_key_object_ ( rdir , name ) :
    """ Get the (key.,object) pair from ROOT-file/directory
    >>> t = f.get_key_object ('A/tup') 
    >>> h = f.get_key_object ('histo') 
    """
    ##
    if not rdir : raise IOError ( "TDirectory is invalid" )
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
        
        return rdir.FindKey ( dirname ) , obj 
                 
# ===============================================================================
## get the object from ROOT-file/directory
#  @code
#  h1 = rfile['MyHisto'      ] 
#  h2 = rfile['MyDir/MyHisto'] 
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@iep.ru
#  @date 2015-07-30
def _rd_getitem_ ( rdir , name ) :
    """ Get the object from ROOT-file/directory
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
## get the key from ROOT-file/directory
#  @author Vanya BELYAEV Ivan.Belyaev@iep.ru
#  @date 2015-07-30
def _rd_key_ ( rdir , name ) :
    """ Get the key from ROOT-file/directory
    """
    ##
    key , obj = _rd_key_object ( rdir , name )
    
    return key 

# ===============================================================================
## get the object from ROOT-file/directory 
#  @code
#  h1 = rfile.MyHisto
#  h2 = rfile.MyDir.MyHisto
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@iep.ru
#  @date 2015-07-30
def _rd_getattr_ ( rdir , name ) :
    """ Get the object from ROOT-file
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
    """ Check the existence of key in the root-file/directory
    >>> if not 'myhisto'  in rfile : break 
    >>> if 'MyDir/MyHisto' in rfile : print 'OK!'
    """
    ##
    if not rdir : return False 
    ##
    name = str ( name ) 
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
    """ Delete an object from ROOT file/directory
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
    """ Get all keys from ROOT file/directory
    >>> keys = rfile.keys() 
    """
    _res = []
    if not rdir : return _res 
    ##
    with ROOTCWD() :
        ##
        rdir.cd() 
        _lst = tuple ( k.GetName() for k in rdir.GetListOfKeys() ) 
        if not _lst :  return _res
        
        for inam in _lst :

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
    """ Iterate over the content of ROOT directory/file:
    >>> for key,obj  in rfile.iteritems()           : print key , obj
    >>> for key,hist in rfile.iteritems( ROOT.TH1 ) : print key , hist
    >>> for key,obj  in rfile.iteritems( lambda name,tkey,obj : name[0]=='M' ) : print key,obj
    """
    ##
    if isinstance ( fun , type ) and issubclass ( fun , ( ROOT.TObject , ) ) : 
        tobj = fun 
        fun  = lambda k,t,o : isinstance ( o , tobj )
    ##
    with ROOTCWD() :
        ##
        rdir.cd() 
        _lst = rdir.GetListOfKeys()
        if _lst : 
            for i in _lst :
                inam   = i.GetName()
                idir   = rdir.GetDirectory ( inam ) 
                if not idir or not no_dir : 
                    obj = rdir.Get ( inam )
                    if fun ( inam , i , obj ) : yield inam , obj
                if recursive and idir and not idir is rdir :
                    for k, o in _rd_iteritems_ ( idir , fun , recursive , no_dir ) :
                        yield k , o

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
    """ Iterate over the keys in ROOT file/directory 
    >>> for key,obj  in rfile.iteritems()           : print key , obj
    >>> for key,hist in rfile.iteritems( ROOT.TH1 ) : print key , hist
    """
    ##
    with ROOTCWD() :
        ##
        rdir.cd() 
        _lst = rdir.GetListOfKeys()
        if _lst : 
            for i in _lst :
                inam   = i.GetName()
                idir   = rdir.GetDirectory ( inam ) 
                if not idir or not no_dir : 
                    if typ is None  or isinstance ( rdir.Get ( inam ) , typ ) : yield inam 
                if recursive and idir  and not idir is rdir :
                    for k in _rd_iterkeys_ ( idir , typ , recursive , no_dir ) :
                        yield k



# =============================================================================a
## iterator oiver keyname/keys pairs  from ROOT file/directory
#  @code
#  for kname, key in rfile.ikeyskeys() :
#    print kname , key.GetClassName()
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@iep.ru
#  @date 2015-07-30
def _rd_ikeyskeys_ ( rdir , recursive = True , no_dir = True ) :
    """ Iterator over  keyname/key pairs  from ROOT file/directory
    >>> for kname, key in rfile.ikeyskeys() :
    ...    print kname , key.GetClassName()
    """
    ##
    with ROOTCWD() :
        ##
        rdir.cd() 
        klst = rdir.GetListOfKeys()
        
        for key in klst :
            
            kname = key.GetName()

            kdir  = rdir.GetDirectory ( kname )            
            if not kdir or not no_dir :
                yield kname, key 
            
            if recursive and kdir and not kdir is rdir :
                for kn,kk in _rd_ikeyskeys_ ( kdir , recursive , no_dir ) :
                    yield kname + '/' + kn , kk 

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
    """ Iterate over the content of ROOT directory/file:
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
    """ Delete an object from ROOT file/directory
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
    """ Get the object from ROOT-file
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
## ``walk'' through the directory and its content
#  @code
#  rfile = ...
#  for r,d,o in rfile.walk() :
#    print 'root dir: %s, %d subdirs, %d objects' % ( r , len ( d ) , len( o ) )
#  @endcode
def _rd_walk_ ( rdir , topdown = True ) :
    """``walk'' through the directory and its content
    >>> rfile = ...
    >>> for r,d,o in rfile.walk() :
    ...  print 'root dir: %s, %d subdirs, %d objects' % ( r , len ( d ) , len( o ) )
    """
    with ROOTCWD() :
        ##
        rdir.cd()
        _subdirs = []
        _objects = []
        _idirs   = [] 
        _lst = rdir.GetListOfKeys()
        if _lst : 
            for i in _lst :
                inam   = i.GetName()
                idir   = rdir.GetDirectory ( inam ) 
                if idir :
                    _subdirs.append  ( inam )
                    _idirs.append    ( idir ) 
                else    :
                    _objects.append  ( inam )
                        
        _subdirs = tuple ( _subdirs )
        _objects = tuple ( _objects )
        
    if topdown : 
        yield rdir.GetName(), _subdirs, _objects 

    for isub in _idirs :
        for jj in _rd_walk_ ( isub , topdown ) : yield jj

    if not topdown :
        yield rdir.GetName(), _subdirs, _objects 
        
# =============================================================================
## use ROOT-file with context-manager
#  @code
#  with ROOT.TFile('ququ.root') as f : f.ls() 
#  @endcode
def _rf_enter_ ( self      ) :
    """ Use ROOT-file with the context manager
    >>> with ROOT.TFile('ququ.root') as f : f.ls() 
    """
    return self

# =============================================================================
## use ROOT-file with context-manager 
#  @code
#  with ROOT.TFile('ququ.root') as f : f.ls() 
#  @endcode
def _rf_exit_  ( self , *_ ) :
    """ Use ROOT-file with the context manager    
    >>> with ROOT.TFile('ququ.root') as f : f.ls()
    """
    try :
        if self and self.IsOpen() : self.Close()
    except:
        pass


# =============================================================================
## create the graphical representation of the directory structure
#  @code
#  rdir = ...
#  for o in rdir.make_graph () : print o.showme() 
#  @endocode
def _rd_make_tree_ ( rdir , parent = None , last = False ) :
    """ Create the graphical representation of the directory structure
    >>> rdir = ...
    >>> for o in rdir.make_graph () : print o.showme() 
    """
    from ostap.logger.utils import DisplayTree

    ## display function for each node 
    def display ( item , isdir ) :        
        name = item.GetName()
        if isdir  : 
            if isinstance ( item , ROOT.TFile ) :
                name = os.path.basename ( name )
            return name + '/'
        return '%s %% %s' % ( name , item .GetClassName() ) 
    
    ## create the root of the tree 
    root = DisplayTree ( rdir , parent , display , isdir = True , last = last ) 
                            
    yield root

    ## get all keys/objects and all keys/subdirectories
    
    dirs = []
    objs = []
    
    for key in rdir.GetListOfKeys() :
        kname = key.GetName()
        kdir  = rdir.GetDirectory ( kname ) 
        if kdir : dirs.append ( kdir )
        else    : objs.append ( key  ) 

    objs.sort ( key = lambda s : s.GetName() )
    dirs.sort ( key = lambda s : s.GetName() )

    ndirs = len ( dirs )
    for i , d  in enumerate ( dirs ) :
        last = ndirs == i + 1
        for l in _rd_make_tree_ ( d , parent = root , last = last ) :
            yield l
        
    nobjs = len ( objs ) 
    for i , o in enumerate ( objs ) :
        last = nobjs == i + 1
        yield DisplayTree ( o , parent = root , isdir = False , last = last )

# =============================================================================
## Show the tree directory structure
#  @code
#  rdir = ...
#  print rdir.show_tree()
#  @endcode
def _rd_show_tree_ ( rdir , prefix = '' ) :
    """ Show the tree directory structure
    >>> rdir = ...
    >>> print rdir.show_tree()
    """
    lines = [ prefix + line.showme () for line in _rd_make_tree_ ( rdir ) ]
    return '\n'.join ( lines ) 

# =============================================================================
## Show the tree directory structure
#  @code
#  rdir = ...
#  rdir.ls_tree()
#  @endcode
def _rd_ls_tree_ ( rdir ) :
    """ Show the tree directory structure
    >>> rdir = ...
    >>> rdir.ls_tree()
    """
    logger.info ( "Directory %s\n%s" % ( rdir.GetName () ,
                                         _rd_show_tree_ ( rdir , prefix = '# ' ) ) ) 

# =============================================================================
## Show the content of the directory as a table
#  @code
#  rdir = ...
#  rdir.ls_table ()
#  @endcode  
def _rd_table_ ( rdir , prefix = '# ' ) :
    """ Show the content of the directory as a table
    >>> rdir = ...
    >>> rdir.ls_table ()
    """
    
    lines   = [  ( n , k.GetClassName() , k.GetObjlen () ) for ( n , k ) in _rd_ikeyskeys_ ( rdir ) ]
    lines.sort()
    maxkey  = 5
    maxtype = 5
    for l in lines :
        maxkey  = max ( maxkey  , len ( l[0] ) )
        maxtype = max ( maxtype , len ( l[1] ) )
        
    fmt_type = '%%-%ds' % ( maxtype + 2 )
    fmt_key  = '%%-%ds' % ( maxkey  + 2 )
    
    table = [ ( fmt_key % 'Key' , fmt_type % 'type' , 'Size' ) ]
    for line in lines :
        
        size = line [ 2 ] 
        if 1024 * 1024 * 1024 <= size :
            size , _  = divmod ( size , 1024 * 1024 * 1024 )
            size =  '%s GB' % size
        elif 1024 * 1024      <= size :
            size , _  = divmod ( size , 1024 * 1024  )
            size =  '%s MB' % size
        elif 1024             <= size :
            size , _  = divmod ( size , 1024 )
            size =  '%s kB' % size
        else :
            size = '%3d  B' % size

        row = fmt_key % line[0] , fmt_type % line[1] , size 
        table.append ( row ) 

    name = rdir.GetName()
    if isinstance ( rdir , ROOT.TFile ) :
        name =  os.path.basename ( name )

    if maxkey + maxtype < len ( name ) :
        name = '<.>' + name [ -maxkey - maxtype : ]
        
    import ostap.logger.table as T
    return T.table ( table , title = '%s' % name , prefix = '# ' , alignment = 'llr' )

# =============================================================================
## Show the content of the directory as a table
#  @code
#  rdir = ...
#  rdir.ls_table ()
#  @endcode  
def _rd_ls_table_ ( rdir , prefix = '# ' ) :
    """ Show the content of the directory as a table
    >>> rdir = ...
    >>> rdir.ls_table ()
    """

    table = _rd_table_ ( rdir , prefix = prefix )

    logger.info ( 'Directory %s:\n%s' % ( rdir.GetName() , table ) )


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

ROOT.TDirectory.get            = _rd_get_
ROOT.TDirectory.keys           = _rd_keys_
ROOT.TDirectory.has_key        = _rd_contains_
ROOT.TDirectory.iteritems      = _rd_iteritems_
ROOT.TDirectory.iterkeys       = _rd_iterkeys_
ROOT.TDirectory.itervalues     = _rd_itervalues_
ROOT.TDirectory.ikeyskeys      = _rd_ikeyskeys_
ROOT.TDirectory.get_key        = _rd_key_
ROOT.TDirectory.get_key_object = _rd_key_object_

# =============================================================================
## some extra stuff 

ROOT.TDirectory.walk         = _rd_walk_
ROOT.TDirectory.__rrshift__  = _rd_rrshift_
ROOT.TNamed    .__rshift__   = _tnamed_rshift_ 

ROOT.TDirectory.show_tree    = _rd_show_tree_
ROOT.TDirectory.ls_tree      = _rd_ls_tree_
ROOT.TDirectory.ls_table     = _rd_ls_table_
ROOT.TDirectory.as_table     = _rd_table_

ROOT.TDirectory.rm           = _rd_rm_

if hasattr ( ROOT.TFile , '__enter__' ) and hasattr ( ROOT.TFile , '__exit__' ) : pass
else :
    ROOT.TFile.__enter__ = _rf_enter_
    ROOT.TFile.__exit__  = _rf_exit_

ROOT.TDirectory.make_tree = _rd_make_tree_
ROOT.TDirectory.show_tree = _rd_show_tree_
ROOT.TDirectory.ls_tree   = _rd_ls_tree_

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
    '+w'       : 'UPDATE'   ,
    'new'      : 'NEW'      ,
    'create'   : 'RECREATE' ,
    'recreate' : 'RECREATE' ,
    'read'     : 'READ'     ,
    'write'    : 'WRITE'    ,
    'update'   : 'UPDATE'   ,
    'append'   : 'APPEND'   ,
    '++'       : 'APPEND'   ,
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
    """ Open/create ROOT-file without making it a current working directory
    >>> print ROOT.gROOT.CurrentDirectory()
    >>> f = ROOT.TFile('test_file.root','recreate')
    >>> print ROOT.gROOT.CurrentDirectory()
    Attention:  IOError exception is raised for invalid file/open_mode
    """
    with ROOTCWD () :
        logger.debug ("Open  ROOT file %s/'%s'" % ( fname , mode ) )
        rinit = rfile._old_init_ ( fname , open_mode ( mode ) , *args )
        if not rfile or not rfile.IsOpen() or rfile.IsZombie () :
            raise IOError   ( "Can't open ROOT file %s/'%s'" % ( fname , mode ) )
        return rinit
    
# ============================================================================
## create ROOT.TFile without making it a current working directory 
#  @code
#  print ROOT.gROOT.CurrentDirectory()
#  f = ROOT.TFile.Open('test_file.root','recreate')
#  print ROOT.gROOT.CurrentDirectory()
#  @endcode
#  @attention  No exceptions are raised for invalid file/open_mode, unless specified 
def _rf_new_open_ ( fname , mode = '' , args = () , exception = False ) :
    """ Open/create ROOT-file without making it a current working directory
    >>> print ROOT.gROOT.CurrentDirectory()
    >>> f = ROOT.TFile.Open('test_file.root','recreate')
    >>> print ROOT.gROOT.CurrentDirectory()
    ATTENTION: No exceptions are raised for invalid file/open_mode, unless specified 
    """
    with ROOTCWD() :
        logger.debug ( "Open  ROOT file %s/'%s'" % ( fname , mode ) )
        try : 
            fopen = ROOT.TFile._old_open_ ( fname , open_mode ( mode ) , *args )
        except ( OSError, IOError ) :
            if exception : raise
            else : 
                logger.error  ( "Can't open ROOT file %s/'%s'" % ( fname , mode ) )
                return None 
            
        if not fopen or not fopen.IsOpen() :
            if  exception : 
                raise IOError ( "Can't open ROOT file %s/'%s'" % ( fname , mode ) )
            else :
                logger.error  ( "Can't open ROOT file %s/'%s'" % ( fname , mode ) )
            
        return fopen
        
# =============================================================================
## Close  ROOT.TFile without making it a current working directory 
#  @code
#  print ROOT.gROOT.CurrentDirectory()
#  f = ROOT.TFile.Open('test_file.root','recreate')
#  f.Close()
#  print ROOT.gROOT.CurrentDirectory()
#  @endcode
def _rf_new_close_ ( rfile , options = '' ) :
    """ Close ROOT-file without making it a current working directory
    >>> print ROOT.gROOT.CurrentDirectory()
    >>> f = ROOT.TFile.Open('test_file.root','recreate')
    >>> print ROOT.gROOT.CurrentDirectory()    
    >>> f.Close()
    >>> print ROOT.gROOT.CurrentDirectory()
    """
    if rfile and rfile.IsOpen() :
        ## with ROOTCWD() :
        logger.debug ( "Close ROOT file %s" % rfile.GetName() ) 
        return rfile ._old_close_ ( options )
            
# =============================================================================
## another name, just for convinince
open = _rf_new_open_

# =============================================================================
if hasattr ( ROOT.TFile , '_new_init_' ) and hasattr ( ROOT.TFile , '_old_init_' ) : pass
else :

    if ROOT.TFile.__init__.__doc__ :
        _rf_new_init_.__doc__  += '\n' + ROOT.TFile.__init__.__doc__

    ROOT.TFile._old_init_   = ROOT.TFile.__init__
    ROOT.TFile._new_init_   = _rf_new_init_ 
    ROOT.TFile.__init__     = _rf_new_init_ 
    
# =============================================================================
if hasattr ( ROOT.TFile , '_new_open_' ) and hasattr ( ROOT.TFile , '_old_open_' ) : pass
else :

    if ROOT.TFile.Open.__doc__ : 
        _rf_new_open_.__doc__  += '\n' + ROOT.TFile.Open.__doc__
        
    _rf_new_open_           = staticmethod ( _rf_new_open_ )
    
    ROOT.TFile._old_open_   = ROOT.TFile.Open
    ROOT.TFile._new_open_   = _rf_new_open_ 
    ROOT.TFile.Open         = _rf_new_open_ 
    ROOT.TFile.open         = _rf_new_open_ 

# =============================================================================
if hasattr ( ROOT.TFile , '_new_close_' ) and hasattr ( ROOT.TFile , '_old_close_' ) : pass
else :

    if ROOT.TFile.Close.__doc__ : 
        _rf_new_close_.__doc__  += '\n' + ROOT.TFile.Close.__doc__
        
    ROOT.TFile._old_close_   = ROOT.TFile.Close
    ROOT.TFile._new_close_   = _rf_new_close_ 
    ROOT.TFile.Close         = _rf_new_close_
    
if not hasattr ( ROOT.TFile , 'close' ) :
    ROOT.TFile.close = ROOT.TFile._new_close_


# =============================================================================
## `topdir': get the top directory for the given directory"""
def top_dir ( rdir ) :
    """`topdir': get the top directory for the given directory"""

    if not rdir : return None

    with ROOTCWD()  :

        top = rdir 
        if   isinstance ( rdir , ROOT.TDirectory ) : top = rdir
        elif hasattr    ( rdir , 'GetDirectory'  ) : top = rdir.GetDirectory()
        
        while top :
            moth = top.GetMotherDir()
            if not moth : return top  
            top  = moth
        else :
            return None 

ROOT.TDirectory.topdir = property ( top_dir , None , None )

# ==============================================================================
## get a full path of the directory
#  @code
#  rdir = ...
#  path = rdirt.full_path 
#  @endcode
def full_path ( rdir ) :
    """ Get a full path of the directory
    >>> rdir = ...
    >>> path = rdir.full_path 
    """
    assert rdir , 'Invalid ROOT TDirectory'
    
    path = rdir.GetPath()
    n    = path.find (':/')
    if n <= 0 :
        logger.error('Invalid directory path')
        return ''
    return path[n+2:] 

ROOT.TDirectory.full_path = property ( full_path , None , None )
    
# ==============================================================================
## Trivial context manager to treat TFile.ReOpen for 'UPDATE' mode
#  @code
#  tdir = ...
#  with REOPEN ( tdir ) as rfile : 
#     ...
#     rfile.Write("", ROOT.TFile.kOverwrite )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date 2019-05-16
class REOPEN(object) :
    """ Trivial context manager to treat TFile.ReOpen for 'UPDATE' mode
    >>> tdir = ...
    >>> with REOPEN ( tdir ) as rfile : 
    ...     ...
    ...     rfile.Write("", ROOT.TFile.kOverwrite )
    """
    def __init__ ( self , rfile ) :
        
        self.tfile = top_dir ( rfile ) 
        assert rfile and self.tfile and isinstance ( self.tfile , ROOT.TFile ),\
               'Invalid file directory %s' % ( rfile.GetName() if rfile else 'INVALID' )

    ##  enter the context: try to ReOpen fiel with <code>ReOpen('UPDATE')</code>
    def __enter__ ( self ) :
        """ Enter the context: make a try to reopen the file with `ReOpen('UPDATE')`
        """
        
        self.read = self.tfile.GetOption() == 'READ'
        if self.read :
            with ROOTCWD() : 
                r = self.tfile.ReOpen ('UPDATE' ) 
                if r < 0 : logger.error ("Can't reopen the file for UPDATE!")
                
        return self.tfile

    ## exit the context : return the mode to <code>READ</code> is needed 
    def __exit__ (  self , *_ ) :
        """ Exit the context: return mode to READ if needed 
        """
        if self.read :
            with ROOTCWD() : 
                r = self.tfile.ReOpen ('READ' ) 
                if r < 0 : logger.error ("Can't reopen the file for READ!")            


from ostap.core.core import _rd_valid_

# ==============================================================================
## length of the directory :  (recursive) number of keys
#  @code
#  rdir = ...
#  len(rdir) 
#  @endcode  
def  _rd_len_ ( rdir ) :
    """ Length of the directory : (recursive) number of keys
    >>> rdir = ...
    >>> len(rdir) 
    """
    
    if not rdir : return 0
    
    with ROOTCWD() :
        ##
        rdir.cd() 
        lst = rdir.GetListOfKeys()

        nkeys = lst.GetSize() if valid_pointer ( lst  ) else 0

        for item in lst :

            inam = item.GetName()
            idir = rdir.GetDirectory ( inam )            

            if idir and not idir is rdir :
                nkeys += _rd_len_ ( idir  ) 
    return nkeys


ROOT.TDirectory.__len__ = _rd_len_

# ========================================================================


# ========================================================================
## copy ROOT file (possible with protocol)
def copy_file ( source , destination , progress = True ) :
    """ Copy ROOT file (possible with protocol)
    """

    ## (1) trivial copy as plain file 
    if os.path.exists ( source ) and os.path.isfile ( source ) : 
        from ostap.utils.basic import copy_file as cpfile
        return cpfile ( source , destination , progress = progress )

    ## (2) check existance with ROOT 
    rf = ROOT.TFile.Open ( source , 'r' , exception = False )
    if not rf : raise IOError ( "``%s'' is not a readable ROOT file!" )

    ## (3) check the destination directory
    destination = os.path.realpath ( destination ) 
    destdir     = os.path.dirname  ( destination )
    if not os.path.exists ( destdir ) : 
        os.makedirs ( destdir )

    ## (4) the main line; use ROOT.TFile::Cp  
    with rf : rf.Cp ( destination , progress ) 
    
    assert os.path.exists ( destination ) and os.path.isfile ( destination ) , \
           "copy_file: Missing destination %s file!" % destination 
    return destination

## alias 
copy_root_file = copy_file

# =============================================================================
## @class RootFiles
#  Some specialzation of gheneric class Files for ROOT files
class RootFiles(Files) :
    """ Specialzation of gheneric class Files for ROOT files
    - it is aware of the ROOT protocols as file names. 
    """
    ## protocols for remote (ROOT) files 
    root_protocols = (
        'root:'  ,
        'xroot:' ,
        'http:'  ,
        'https:' ,    
    )
    # =========================================================================
    ## Has ROOT protocol in th efile name? 
    @staticmethod 
    def has_protocol ( fname ) :
        """ Has ROOT protocol in th efile name?
        """
        for p in RootFiles.root_protocols :
            if p and fname.startswith ( p ) : return p
            return False
        
    # =========================================================================
    ## Strip protocol from the file name    
    @staticmethod 
    def strip_protocol ( fname ) :
        """ Strip ROOT protocol from the file name"""
        for p in RootFiles.root_protocols :
            if p and fname.startswith ( p ) : return fname.replace ( p , '' ) 
        return fname

    # =========================================================================
    ## Get the list of files from the patterns 
    def the_files ( self ) :
        """ Get the list of files from the patterns"""

        # 1. get the regular files
        files = set ( Files.the_files ( self ) )

        # 2. add explicit ROOT patterns 
        for pattern in  self.patterns :
            if RootFiles.has_protocol ( pattern ) :
                files.add ( pattern )
        
        return tuple ( sorted ( files ) ) 

    # =========================================================================
    ## get a common path (prefix) for all files in collection
    #  - protocols are ignored 
    @property 
    def commonpath ( self ) :
        """'commonpath': common path (prefix) for all files in collection
        - protocols are ignored 
        """
        ##
        if any ( RootFiles.has_protocol ( f ) for f in self.files ) :
            files = []
            for f in self.files :
                files .append ( RootFiles.strip_protocol ( f ) ) 
        else : files = self.files
        
        if not files : return '' 

        from ostap.utils.basic import commonpath 
        cp = commonpath ( os.path.abspath ( f ) for f in files )
        return cp if os.path.isdir ( cp ) else os.path.dirname ( cp )
    
    # =========================================================================
    ## copy all the files to new directory
    #  - new directory will be created (if needed)
    #  - common path (prefix) for all files will be replaced by new directory
    def copy_files ( self , new_dir = None , parallel = False , also_bad = False ) :
        """ Copy all the files to new directory
        - new directory will be created (if needed)
        - common path (prefix) for all files will be replaced by new directory
        """
        
        ## use the temporary directory
        if new_dir is None or new_dir == '' :
            import ostap.utils.cleanup as CU
            new_dir = CU.CleanUp.tempdir()
        
        ## create directory if needed 
        if not os.path.exists ( new_dir ) : os.makedirs ( new_dir )
        
        from ostap.utils.basic  import writeable

        assert writeable ( new_dir ), \
            "copy_files: the destination directory `%s' is not writable!" % new_dir 
        
        nd = os.path.abspath  ( new_dir )
        nd = os.path.normpath ( nd      ) 

        cp = self.commonpath
        
        files_to_copy = set ( self.files )
        if also_bad and self.bad_files :
            files_to_copy |= set ( self.bad_files )

        for_copy = [] 
        for f in files_to_copy :
            fs   = os.path.abspath ( RootFiles.strip_protocol ( f ) ) 
            nf   = fs.replace ( cp , nd ) 
            nf   = os.path.abspath ( nf )
            pair = f , nf
            for_copy.append ( pair )

        nfiles = len ( for_copy ) 

        ## parallel copy 
        from   ostap.utils.basic      import numcpu
        if parallel and 1 < nfiles and 1 < numcpu () :
            
            from   ostap.utils.utils  import which
            if which ( 'parallel' ) : from ostap.utils.parallel_copy    import copy_files as parallel_copy
            else                    : from ostap.parallel.parallel_copy import copy_files as parallel_copy
            
            copied = parallel_copy ( for_copy , maxfiles = 1 , copier = copy_root_file , progress = not self.silent )
            copied = [ f [ 1 ] for f in copied ]

        ## sequential copy
        else :
            
            copied = []
            from ostap.utils.progress_bar import progress_bar
            for f, nf in progress_bar ( for_copy , silent = self.silent or nfiles <=1 ) :
                copied.append ( copy_root_file ( f , nf , progress = ( 1 == nfiles ) and self.verbose ) ) 
                
        if not self.silent :
            logger.info ( "copy_files: #%d files are copied to '%s'" %  ( len ( copied ) , nd ) )
            
        return self.clone ( files = copied )

    # ==========================================================================
    ## Get the size of the file
    #  @param fname input file name
    #  @return ROOT file sizeor -1 for non-existingt/invalid files
    def get_file_size ( self , fname ) :
        """ Get the size of the file
        - fname : input file name
        - return file size or -1 for non-existing/invalid files
        """
        
        ## (1) try as a regular file: 
        fsize = Files.get_file_size ( self , fname )
        if 0 <= fsize : return fsize                      ## RETURN

        ## (2) try to open it as ROOT file
        # ====================================================================+
        try : # ==============================================================+
            # ================================================================+
            with ROOT.TFile.Open ( fname , 'r' , exception = True ) as r :
                return r.GetSize()                        ## RETURN
            # ================================================================+
        except ( OSError , IOError ) : # =====================================+
            # ================================================================+
            pass

        return -1

    # =========================================================================
    ## merge all files using <code>hadd</code> script from ROOT
    #  @param output  name of the output merged file, if None,
    #                 the temporary name will be generated,
    #                 that will be deleted at the end of the session
    #  @param opts   options for command <code>hadd</code>
    #  @return the name of the merged file
    # OPTIONS:
    # -a                                   Append to the output
    # -k                                   Skip corrupt or non-existent files, do not exit
    # -T                                   Do not merge Trees
    # -O                                   Re-optimize basket size when merging TTree
    # -v                                   Explicitly set the verbosity level: 0 request no output, 99 is the default
    # -j                                   Parallelize the execution in multiple processes
    # -dbg                                 Parallelize the execution in multiple processes in debug mode (Does not delete partial files stored inside working directory)
    # -d                                   Carry out the partial multiprocess execution in the specified directory
    # -n                                   Open at most 'maxopenedfiles' at once (use 0 to request to use the system maximum)
    # -cachesize                           Resize the prefetching cache use to speed up I/O operations(use 0 to disable)
    # -experimental-io-features            Used with an argument provided, enables the corresponding experimental feature for output trees
    # -f                                   Gives the ability to specify the compression level of the target file(by default 4) 
    # -fk                                  Sets the target file to contain the baskets with the same compression
    #                                      as the input files (unless -O is specified). Compresses the meta data
    #                                      using the compression level specified in the first input or the
    #                                      compression setting after fk (for example 206 when using -fk206)
    # -ff                                  The compression level use is the one specified in the first input
    # -f0                                  Do not compress the target file
    # -f6                                  Use compression level 6. (See TFile::SetCompressionSettings for the support range of value.)    
    def hadd ( self , output = None , opts = "-ff -O" , **kwargs ) :
        """ Merge all files using <code>hadd</code> script from ROOT
        - `output`  name of the output merged file
        - `opts`   options for command <code>hadd</code>
        It returns the name of the merged file
        
        If no output file name is specified, the temporary name
        will be generate and the temporary file will be deleted
        at the end of the session

        OPTIONS:
        # -a                                   Append to the output
        # -k                                   Skip corrupt or non-existent files, do not exit
        # -T                                   Do not merge Trees
        # -O                                   Re-optimize basket size when merging TTree
        # -v                                   Explicitly set the verbosity level: 0 request no output, 99 is the default
        # -j                                   Parallelize the execution in multiple processes
        # -dbg                                 Parallelize the execution in multiple processes in debug mode (Does not delete partial files stored inside working directory)
        # -d                                   Carry out the partial multiprocess execution in the specified directory
        # -n                                   Open at most 'maxopenedfiles' at once (use 0 to request to use the system maximum)
        # -cachesize                           Resize the prefetching cache use to speed up I/O operations(use 0 to disable)
        # -experimental-io-features            Used with an argument provided, enables the corresponding experimental feature for output trees
        # -f                                   Gives the ability to specify the compression level of the target file(by default 4) 
        # -fk                                  Sets the target file to contain the baskets with the same compression
        #                                      as the input files (unless -O is specified). Compresses the meta data
        #                                      using the compression level specified in the first input or the
        #                                      compression setting after fk (for example 206 when using -fk206)
        # -ff                                  The compression level use is the one specified in the first input
        # -f0                                  Do not compress the target file
        # -f6                                  Use compression level 6. (See TFile::SetCompressionSettings for the support range of value.)                            
        """
        from ostap.utils.utils import hadd as hadd_
        return hadd_ ( self.files , output = output , opts = opts  , **kwargs )

    # =========================================================================
    ## Split data into severals chunks of smaller size and merge (using <code>hadd</code>) each chunk
    #  @code
    #  data   = ..
    #  merged = data.merge_chunks ( 10 ) 
    #  @endcode 
    def __merge_chunks (  self , chunks , opts = '-ff -O' , parallel = True ) :
        """ Split data into severals chunks of smaller size and merge (using <code>hadd</code>) each chunk        
        >>> data   = ..
        >>> merged = data.merge_chunks ( 10 ) 
        """

        import ostap.utils.cleanup as  CU
        cu     = CU.CleanUp()
        tmpdir = cu.tmpdir

        if parallel and chunks and 2 <= len ( chunks[0].files ) :
            
            from   ostap.parallel.parallel import WorkManager
            from   ostap.utils.utils       import hadd2       as hadd_
            
            if     '-j' in opts : opts = opts.replace('-j','')
            if not '-v' in opts : opts = '%s -v 0' % opts

            pargs  = [ ( c.files , None , tmpdir , opts ) for c in chunks ] 
            
            wm     = WorkManager ( silent = True , progress = not self.silent )
            
            merged = [ o for o in wm.iexecute ( hadd_ , pargs , progress = not self.silent ) ]

            return self.clone ( files = merged ) 
            
        if not '-j' in opts : opts = '%s -j' % opts

        output = []

        if self.silent and not '-v' in opts : opts = '%s -v 0' % opts
        for c in progress_bar ( chunks , silent = self.silent or not chunks ) :
            output.append (  c.hadd ( opts = opts , dir = tmpdir ) ) 
                    
        return self.clone ( output )

    # =========================================================================
    ## Split data into severals chunks of smaller size and merge (using <code>hadd</code>) each chunk
    #  @code
    #  data   = ..
    #  merged = data.merge_chunks ( 10 ) 
    #  @endcode 
    def merge_chunks (  self , chunk_size , opts = '-ff -O' , parallel = True ) :
        """ Split data into severals chunks of smaller size and merge (using <code>hadd</code>) each chunk        
        >>> data   = ..
        >>> merged = data.merge_chunks ( 10 ) 
        """
        chunks = self.split_chunks ( chunk_size )
        return self.__merge_chunks ( chunks , opts = opts , parallel = parallel )
        
    # =========================================================================
    ## Split data into severalgroup and merge (using  <code>hadd</code>) each group
    #  @code
    #  data   = ..
    #  merged = data.merge_groups ( 10 ) 
    #  @endcode 
    def merge_groups (  self , groups , opts = '-ff -O' , parallel = True ) :
        """ Split data into severals groups and merge (using <code>hadd</code>) each group
        >>> data   = ..
        >>> merged = data.merge_groups ( 10 ) 
        """
        chunks = self.split_groups ( groups )  
        return self.__merge_chunks ( chunks , opts = opts , parallel = parallel )

# =============================================================================
_decorated_classes_ = (
    ROOT.TFile       ,
    ROOT.TDirectory     
    )

_new_methods_   = (
    #
    ROOT.TFile.__init__          ,
    ROOT.TFile.Open              ,
    ROOT.TFile.__enter__         ,
    ROOT.TFile.__exit__          ,
    #
    ROOT.TDirectory.__getitem__  , 
    ROOT.TDirectory.__setitem__  , 
    ROOT.TDirectory.__contains__ , 
    ROOT.TDirectory.__getattr__  , 
    ROOT.TDirectory.__delitem__  , 
    ROOT.TDirectory.__iter__     , 
    #
    ROOT.TDirectory.get          , 
    ROOT.TDirectory.keys         ,
    ROOT.TDirectory.has_key      ,
    ROOT.TDirectory.iteritems    ,
    ROOT.TDirectory.iterkeys     ,
    ROOT.TDirectory.itervalues   ,
    #
    ROOT.TDirectory.__rrshift__  , 
    ROOT.TNamed    .__rshift__   , 
    )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
