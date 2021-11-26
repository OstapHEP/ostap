#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file
#  Module with helper functions to convert objects into HepDATA format
#
#  Supported types: 
#  - TH1D, TH1D
#  - TGraphErrors
#  - TGraphAsymmErrors
#
#  @see http://hepdata.cedar.ac.uk/resource/sample.input
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
#  
# =============================================================================
""" Helper functions to convert the objects into HepDATA format
Supported types: 
- TH1D, TH1D
- TGraphErrors
- TGraphAsymmErrors
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'HepDataFile'  , ## the HepDATA file
    'HepData'      , ## HepDATA record 
    ) 
# =============================================================================
import ROOT
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.utils.hepdata' )
else                       : logger = getLogger( __name__ )
# =============================================================================
logger.debug ( 'HepDATA format and conversion routines')
# =============================================================================
from sys import version_info as python_version
if python_version.major > 2 : items_loop = lambda d : d.    items () 
else                        : items_loop = lambda d : d.iteritems () 
    
# =============================================================================
## fields required for each dataset
dataset_fields = ( 'location'  , ## e.g. Figure 5a
                   'dscomment' , 
                   'reackey'   , ## e.g. P P --> Z0 Z0 X 
                   'obskey'    , ## e.g. SIG or DSIG/DPT
                   'qual'      , ## see http://hepdata.cedar.ac.uk/resource/sample.input
                   'xheader'   , ## e.g. PT IN GEV 
                   'yheader'   ) ## e.g. DSIG/DPT IN PB/GEV 
## fields required for each hepdata file 
hepfile_fields = ( 'author'     ,
                   'reference'  , ## format :  "REFERENCE : YEAR"
                   'doi'        , 
                   'status'     ,
                   'experiment' , ## e.g. CERN-LHC-LHCb
                   'detector'   , ## e.g. LHCb 
                   ## 'spiresId'   , ## id in SPIRES  (http://www.slac.stanford.edu/spires/) 
                   'inspireId'  , ## id in INSPIRE (http://inspirehep.net/)
                   'cdsId'      , ## id in CDS     (http://cds.cern.ch/)
                   'durhamId'   , ## will be added after submission 
                   'title'      ,
                   'comment'    )

from collections import defaultdict
METAINFO = defaultdict(list)

# =============================================================================
## Helper base class for HepDATA formatting 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
class HepDataBase(object) :
    """
    Helper base class for HepDATA formatting 
    """
    def __init__ ( self , metainfo = defaultdict(list) , **kwargs ) :
        ##
        from copy import deepcopy 
        self.meta = deepcopy(metainfo)
        for k , v in items_loop ( kwargs ) :
            if isinstance ( v , (list,tuple) ) :
                for i in v :
                    if not i in self.meta[k] : self.meta[k].append ( i )
            else :
                if not v in self.meta[k] : self.meta[k].append ( v )
                
    ## add meta-information
    def add ( self , key , info) : self.meta[key].append ( info ) 
        
    ## get missing keys 
    def missing ( self , keys ) :
        m = set()
        for key in keys : 
            if not key in self.meta : m.add ( key )
        return m 
    
    def format ( self , header , keys ) :
        "Make a proper valid representation of dataset"
        #
        lines = []
        if header : lines = [ '*%s:' % header ] 
        #
        for key in keys :
            items = self.meta[key]
            for i in items :
                lines.append ( '*%s: %s' % ( key , i ) )
                
        ## the actual format of the data 
        self.formatData( lines ) 
            
        return '\n'.join ( lines )

# ============================================================================-
## @class HepDataFile
#  Helper class to prepare the file with HepDATA information
#  @code 
#  metainfo = {...}
#  hf = HepDataFile(  *metainfo )
#  print hf
#  @endcode
#  Or use it as context manager 
#  @code
#  with HepDataFile(  *metainfo ) as hf :
#  ...            <do something>
#  @endcode 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2015-08-29
class HepDataFile(HepDataBase) :
    """Helper class to prepare the file with HepDATA information
    >>> metainfo = {...}
    >>> hf = HepDataFile(  *metainfo )
    >>> print hf
    Or as context manager 
    >>> with HepDataFile(  *metainfo ) as hf :
    ...       <do something> 
    """
    def __init__ ( self ,
                   filename = 'HepDATA.txt'     ,
                   datasets = []                ,
                   metainfo = defaultdict(list) , **kwargs ) :
        
        HepDataBase.__init__ ( self , metainfo , **kwargs ) 
        self.filename = filename 
        self.datasets = datasets
        m = self.missing ( hepfile_fields ) 
        if m : logger.warning('HepFILE missing keys: %s' % m )
        
    def append      ( self , dataset ) : self.datasets.append ( dataset )
    def add         ( self , dataset ) : self.append ( dataset ) 
    def __iadd__    ( self , dataset ) :
        self.append ( dataset )
        return self 
    def __lshift__  ( self , dataset ) :
        self += dataset
        return self
    def formatData ( self , lines ) : pass 
    def __enter__ ( self ) :
        
        result = self.format( '' , hepfile_fields ) + '\n'
        self.file = open ( self.filename , 'w' )
        self.file.write ( result + '\n' )
        return self
        
    def __exit__  ( self , *_ ) :
        for ds in self.datasets :
            self.file.write( '\n%s\n' % ds )
        self.file.write ('\n*e\n' )
        self.file.close ()

    def __str__  ( self ) :
        result = self.format( '' , hepfile_fields ) + '\n' 
        for ds in self.datasets :
            result += str( ds ) + '\n'
        return result + '\n*e'
    
# =============================================================================
## Convert simple object (presumably historgam or graph) to HepDATA format
#  @code
#  histo    = ...
#  metadata = { ... }
#  ds = HepData ( histo , **metadata  )
#  @endcode
#  @attention The object must have ``toHepDATA'' method 
#  Currently followong types are supported 
#  - ROOT.TH1F
#  - ROOT.TH1D
#  - ROOT.TGraphErrors 
#  - ROOT.TGraphAsymmErrors 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
class HepData(HepDataBase) :
    """ Convert simple object (presumably historgam or graph) to HepDATA format
    >>> histo    = ...
    >>> metadata = { ... }
    >>> ds = HepData ( histo , **metadata )
    The object must have ``toHepDATA'' method
    Currently following types are supported 
    - ROOT.TH1F
    - ROOT.TH1D
    - ROOT.TGraphErrors 
    - ROOT.TGraphAsymmErrors 
    """
    def __init__ ( self  ,
                   histo ,
                   ##
                   syst1 = '' , ## the first  systematic 
                   syst2 = '' , ## the second systematic 
                   syst3 = '' , ## the third  systematic 
                   ##
                   metainfo = defaultdict(list) , **kwargs ) :
        ##
        HepDataBase.__init__ ( self , metainfo , **kwargs )
        self.histo = histo 
        m = self.missing ( dataset_fields ) 
        if m : logger.warning('HepData missing keys: %s' % m )
        from copy import deepcopy 
        self.syst1  = deepcopy ( syst1 )
        self.syst2  = deepcopy ( syst2 )
        self.syst3  = deepcopy ( syst3 )

        self.result =  self.format ('dataset',dataset_fields)+'\n'
    ## the actual output :-) 
    def __str__ ( self ) : return self.result ## self.format ('dataset',dataset_fields)+'\n'

    ## the most important line: the proper delegation
    def formatData ( self , the_lines ) :
        """The most important line: the proper delegation"""
        args       = (self.syst1,self.syst2,self.syst3)
        ## delegate to the specific method of the histogram/graph
        result     = self.histo.toHepDATA ( *args )      ## invoke the specific method 
        the_lines += result.split('\n')

# ==============================================
## helper function to decode the systematic uncertainties
#  Systematic could be 
#   - just a string
#   - an object with index:   obj  [ibin]
#   - a kind of function:     func (ibin)
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07  
def get_syst ( syst , *index ) :
    """Helper function to decode the systematic uncertainties
    Systematic could be 
    - just a string
    - an object with index:   obj  [ibin]
    - a kind of function:     func (ibin)
    """
    if   isinstance ( syst , str )   : return syst
    elif syst and hasattr  ( syst , '__getitem__' )  : return str ( syst [  index ] )
    elif syst and callable ( syst )                  : return str ( syst ( *index ) )
    elif syst : raise AttributeError("Invalid systematic %s/%s" % ( syst , type( syst ) ) )
    
    return ''
    
# =============================================================================
## Dump TH1 in HepData format with the optional specification of
#  systematic uncertainties
#  @code
#  h1 = ...
#  print h1.toHepDATA()
#  print h1.toHepDATA( syst1 ='0.01:detector' )
#  print h1.toHepDATA( syst1 ='0.01:detector' ,
#                      syst2 = [ 0.02 ,0.03 , .... , 0.8] )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07  
def _h1_hepdata_ ( histo      ,
                   syst1 = '' , 
                   syst2 = '' ,
                   syst3 = '' ) :
    """
    Dump -histogram in HepData -compatible format with \
    optional specification of up to three systematic uncertainnties 
    >>> h1 = ...
    >>> print h1.toHepDATA()
    >>> print h1.toHepDATA( syst1 ='0.01:detector' )
    >>> print h1.toHepDATA( syst1 ='0.01:detector' ,
    ...                     syst2 = [ 0.02 ,0.03 , .... , 0.8] )
    >>> print h1.toHepDATA( syst1 ='0.01:detector' ,
    ...                     syst2 = [ 0.02 ,0.03 , .... , 0.8]   ,
    ...                     syst3 = lambda i : '%s:last' % h2[i] ) 
    
    """    
    ## the basic format 
    fmt   = '   %s TO %s ; %s +- %s ; '
    lines = [ '*data: x : y ' ]
    
    index = 0 

    for item in iteritems( histo ) :
        
        i = item[0] ## bin-number 
        x = item[1] ## x
        y = item[2] ## y 

        line = fmt %  ( x.value () - x.error () ,
                        x.value () + x.error () ,
                        y.value () , y.error () )
        
        #
        ## up to three systematic uncertainties 
        #

        s1   = get_syst ( syst1 , index ) 
        s2   = get_syst ( syst2 , index ) 
        s3   = get_syst ( syst3 , index )

        index += 1 

        if s2 and not s1 :
            raise AttributeError ( "HepData: syst2 is specified without syst1" )
        if s3 and not s2 :
            raise AttributeError ( "HepData: syst3 is specified without syst2" )
        
        if s1 :
            line +=        " (DSYS=%s" % s1
            if s2 : line += ",DSYS=%s" % s2
            if s3 : line += ",DSYS=%s" % s3
            line += ") ; "
            
        lines.append (line)
        
    lines.append ('*dataend:')
    return '\n'.join(lines)+'\n'

# =============================================================================
## Dump TGraphAsymmErrore in HepData format with the optional specification of
#  systematic uncertainties
#  @code
#  graph = ...
#  print graph.toHepDATA()
#  print graph.toHepDATA ( syst1 ='0.01:detector' )
#  print graph.toHepDATA ( syst1 ='0.01:detector' ,
#                          syst2 = [ 0.02 ,0.03 , .... , 0.8] )
#  @endcode
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07  
def _tgae_hepdata_ ( graph      ,
                     syst1 = '' , 
                     syst2 = '' ,
                     syst3 = '' ) :
    """
    Dump -histogram in HepData -compatible format with \
    optional specification of up to three systematic uncertainnties 
    >>> h1 = ...
    >>> print h1.toHepDATA()
    >>> print h1.toHepDATA( syst1 ='0.01:detector' )
    >>> print h1.toHepDATA( syst1 ='0.01:detector' ,
    ...                     syst2 = [ 0.02 ,0.03 , .... , 0.8] )
    >>> print h1.toHepDATA( syst1 ='0.01:detector' ,
    ...                     syst2 = [ 0.02 ,0.03 , .... , 0.8]   ,
    ...                     syst3 = lambda i : '%s:last' % h2[i] ) 
    
    """
    
    ## the basic format 
    fmt   = '   %s TO %s ; %s +%s -%s '
    lines = [ '*data: x : y ' ]

    index = 0 
    for item in items_loop ( graph ) :
        
        i   = item[0] ## bin-number 
        x   = item[1] ## x
        exl = item[2]
        exh = item[3]
        y   = item[4] ## x
        eyl = item[5]
        eyh = item[6]
        
        line = fmt %  ( x  - abs ( exl ) ,
                        x  + abs ( exh ) ,
                        y  , abs ( eyh ) , abs ( eyl ) )
        #
        ## up to three systematic uncertainties 
        #

        s1   = get_syst ( syst1 , index ) 
        s2   = get_syst ( syst2 , index ) 
        s3   = get_syst ( syst3 , index ) 

        index += 1
        
        if s2 and not s1 :
            raise AttributeError ( "HepData: syst2 is specified without syst1" )
        if s3 and not s2 :
            raise AttributeError ( "HepData: syst3 is specified without syst2" )
        
        if s1 :
            line +=        " (DSYS=%s" % s1
            if s2 : line += ",DSYS=%s" % s2
            if s3 : line += ",DSYS=%s" % s3
            line += ") ; "

        line += " ; "
        lines.append (line)
        
    lines.append ('*dataend:')
    return '\n'.join(lines)+'\n'


# =============================================================================
for t in ( ROOT.TH1D         ,
           ROOT.TH1F         ,
           ROOT.TGraphErrors ) : 
    
    t . toHepDATA = _h1_hepdata_
    t . toHepData = _h1_hepdata_

ROOT.TGraphAsymmErrors .toHepDATA = _tgae_hepdata_  ## different one 
ROOT.TGraphAsymmErrors .toHepData = _tgae_hepdata_  ## different one 

# =============================================================================
if '__main__' == __name__ :
        
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
# The END 
# =============================================================================
