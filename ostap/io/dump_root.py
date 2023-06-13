#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
# @file dump_root.py
#
# Helper module to solve the problem with evolution of ROOT serialization for 
# pickle-based databases:
# - <code>zipshelve</code>,
# - <code>bz2shelve</code>,
# - <code>lzshelve</code>,
# - <code>rootshelve</code>,
# - <code>rootshelve</code>,
# - <code>sqliteshelve</code>.
#
# If database with some pickled ROOT-objects is created using "old" verison
# of ROOT, and if the serialization methods for these objects are updated
# reading of the objects with newer verison of ostaap/ROOT
# can fail with  following error  messaged from ROOT, e.g. this one for class <code>ROOT.TH1D</code>:
# @code
# Error in <TBufferFile::ReadClassBuffer>: Could not find the StreamerInfo for version 2 of the class TH1D, object skipped at offset 19
# Error in <TBufferFile::CheckByteCount>: object of class TH1D read too few bytes: 2 instead of 878
# @endcode 
#
# The solution is rather simple:
#
# - one needs to prepare a plain ROOT file with the object of the problematic type.
#   @attention: The  file needs to be   created using "old" version of ROOT.
# @code
# h = ROOT.TH1D( ... ) ## instance of  problematic class
# import ostap.io.dump_file as DF
# root_file_name = DF.dump_root ( [ h ] ) 
# @endcode 
#
# - within new ostap/python  session (with newer version of ROOT), one need  to open this
#   file and read all objects. It activates the "old" streamer machinery
# @code
# import ostap.io.dump_file as DF
# DF.read_file ( root_file_name ) 
# @endcode
#  or, more explicit:
# @code 
# import ostap.io.root_file 
# rootfile = ROOT.TFile ( root_file_name ,  'r' )
# for k in root_file :
#   print k, root_file[k]
# @endcode 
#
# - then one needs to clone the pickle-based database into new database and use the new database 
#  @code
#  old_db = DBASE.open   ( 'old_dbase.db' , 'r')
#  new_db = old_db.clone ( 'new_dbase.db' ) 
#  @endcode 
#
# @author Vanya BELYAEV Ivan.Belyaev@cern.ch
# @date   2019-10-04 
# =============================================================================
""" Helper module to solve the problem with evolution of ROOT serialization for
all pickle-based databases:
- zipshelve,
- bz2shelve,
- lzshelve,
- rootshelve,
- rootshelve,
- sqliteshelve

If database with some pickled ROOT-objects is created using "old" verison
of ROOT, and if the serialization methods for these objects are updated
reading of the objects with newer verison of ostaap/ROOT
#can fail with  following error  messaged from ROOT, e.g. this one for class <code>ROOT.TH1D</code>:

> Error in <TBufferFile::ReadClassBuffer>: Could not find the StreamerInfo for version 2 of the class TH1D, object skipped at offset 19
> Error in <TBufferFile::CheckByteCount>: object of class TH1D read too few bytes: 2 instead of 878


The solution is rather simple:

1. one needs to prepare a plain ROOT file with the object of the problematic type.
Attention: The  file needs to be created using ``old'' version of ROOT.

>>> h = ROOT.TH1D( ... ) ## instance of  problematic class
>>> import ostap.io.dump_file as DF
>>> root_file_name = DF.dump_root ( [ h ] ) 

2. within new ostap/python  session (with newer version of ROOT), one need  to open this
file and read all objects. It activates the ``old'' ROOT streamer machinery

>>> import ostap.io.dump_file as DF
>>> DF.read_file ( root_file_name ) 

or, more explicit:

>>> import ostap.io.root_file 
>>> with ROOT.TFile ( root_file_name ,  'r' )   as rf 
>>> ...  for k in rf :
>>> ...  print k, rf[k]

3. then one needs to clone the pickle-based database into new database and use the new database 

>>> old_db = DBASE.open   ( 'old_dbase.db' , 'r')
>>> new_db = old_db.clone ( 'new_dbase.db' ) 

"""
# =============================================================================
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2015-07-31"
__version__ = "$Revision$" 
# =============================================================================
__all__ = (
    'dump_root' , ## dump  ROOT objects to the file 
    'read_root' , ## real all ROOT objects from the file 
    )
# =============================================================================
import ostap.io.root_file
from   ostap.core.meta_info import root_version_int 
import ROOT
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' == __name__ : logger = getLogger ( 'ostap.io.dump_root' )
else                      : logger = getLogger ( __name__             )
# =============================================================================
# dump list of objects into ROOT-file 
def dump_root ( objects , rfile  = '' ) :

    if not rfile :
        rfile = 'ROOT_Objects_%s.root' % root_version_int

    if isinstance ( objects , ROOT.TObject ) : objects = [ objects ]
    
    import ostap.io.root_file
    
    from   ostap.core.core import ROOTCWD
    
    with ROOTCWD(), ROOT.TFile ( rfile , 'RECREATE' ) as f :
        f.cd () 
        for o in objects :
            o.Write()
            logger.info ( 'Dump object of type %s' %   type ( o ) ) 
        f.ls    ()
            
    return rfile 

# =============================================================================
## read all objects from  ROOT file 
def read_root ( fname ) :    
    """Read all objects from  ROOT file """
    
    import ostap.io.root_file
    with ROOT.TFile.Open ( fname , 'READ' ) as f :
        
        f.ls()
        keys = f.GetListOfKeys()
        
        for k in keys :
            key = "%s;%d" % ( k.GetName() , k.GetCycle() )            
            obj = f.Get( key )            
            logger.info ( 'Read key/object  %s/%s' % ( key , type  (  obj ) ) )
                   
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
    
    v1   = ROOT.RooRealVar    ( 'v1' , '',1)
    v2   = ROOT.RooRealVar    ( 'v2' , '',1)
    vars = ROOT.RooArgSet     (  v1  ,    v2 )
    varl = ROOT.RooArgList    (  v1  ,    v2 )
    v3   = ROOT.RooFormulaVar ( 'v3' , 'v1*v2', varl )
    
    objects = [
        
        ROOT.TH1F ( 'h1f' , '' , 10 , 0 , 1 ) , 
        ROOT.TH1D ( 'h1d' , '' , 10 , 0 , 1 ) , 
        ROOT.TH2F ( 'h2f' , '' ,  3 , 0 , 1 , 3 , 0 , 1 ) , 
        ROOT.TH2D ( 'h2d' , '' ,  3 , 0 , 1 , 3 , 0 , 1 ) , 
        ROOT.TH3F ( 'h3f' , '' ,  2 , 0 , 1 , 2 , 0 , 1 , 2 , 0 , 1 ) , 
        ROOT.TH3D ( 'h3d' , '' ,  2 , 0 , 1 , 2 , 0 , 1 , 2 , 0 , 1 ) ,
        ROOT.TProfile ( 'prof1','',10,0,1) , 
        ROOT.TF1      (  'f1'  , 'sin(x)/x' ) , 
        ROOT.TF2      (  'f2'  , 'sin(x)*sin(y)/x/y') ,    
        ROOT.TGraph            ( 5 ) ,
        ROOT.TGraphErrors      ( 5 ) ,
        ROOT.TGraphAsymmErrors ( 5 ) ,
        ROOT.TTree  ('tt','') ,
        ROOT.TBox   ( 0 , 0 , 1 , 1 ) ,
        ROOT.TLine  ( 0 , 0 , 1 , 1 ) ,
        ROOT.TAxis  ( 2 , 0 , 1 ) ,
        ROOT.TPaletteAxis  ()  , 
        ROOT.TFitResult    ()  , 
        ROOT.TCut   ( 'cut' ) ,
        ##
        ROOT.RooRealVar    () ,
        ROOT.RooArgSet     () ,
        ROOT.RooArgList    () ,
        v1 , v2 , vars , varl , v3 , 
        ROOT.RooDataSet    ( 'ds', '', vars ) ,
        ROOT.TTreeFormula  () ,
        ROOT.TPaletteAxis  () , 
        ROOT.TGaxis        ()
        ]

    file_name = dump_root ( objects )
    read_root ( file_name ) 
    
# =============================================================================
##                                                                      The END 
# =============================================================================
