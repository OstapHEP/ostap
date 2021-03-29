#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/fitting/roocmdarg.py
#  Module with decoration of RooCmdArg 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2011-06-07
# =============================================================================
"""Decoration of RooCmdArg obejcts 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2011-06-07"
__all__     = (
    'match_arg'      , ## check the argument name mathhing 
    'check_arg'      , ## check the presense of argument in the list
    'nontrivial_arg' , ## check the presense of nontrivial arguments 
    ) 
# =============================================================================
import ROOT
from   ostap.core.ostap_types       import string_types, integer_types 
from   ostap.utils.utils            import chunked
import ostap.fitting.roocollections 
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger( 'ostap.fitting.roocmdarg' )
else                       : logger = getLogger( __name__ )
# =============================================================================
## print RooCmdArg object 
def _rca_print_ ( self ) :
    """Print RooCmdArg object 
    """
    name = self.GetName()

    ## RooAbsReal::plotOn arguments
    if   'DrawOption'           == name : return "DrawOptions('%s')"    %   self.getString ( 0 )
    elif 'SliceVars'            == name : return "Slice({.})"
    elif 'SliceCat'             == name : return "Slice({.})"
    elif 'Project'              == name : return "Project({.})"
    elif 'ProjData'             == name : return "ProjWData({.})"
    elif 'ProjData'             == name : return "ProjWData({.})"
    elif 'Asymmetry'            == name : return "Asymmetry({.})"
    elif 'Precision'            == name : return "Precision(%s)"         %   self.getDouble(0) 
    elif 'ShiftToZero'          == name : return "ShiftToZero()"
    elif 'Normalization'        == name : return "Normalizarion(%s)"     %   self.getDouble(0)
    elif 'RangeWithName'        == name : return "Range('%s',%s)"        % ( self.getString(0) ,
                                                                             self.getBool  ( ) ) 
    elif 'Range'                == name : return "Range(%s,%s,%s)"       % ( self.getDouble(0) ,
                                                                             self.getDouble(1) ,
                                                                             self.getBool  ( ) )
    elif 'NormRange'            == name : return "NormRange('%s')"       %   self.getString(0)
    elif 'VLines'               == name : return "VLines()" 
    elif 'LineColor'            == name : return 'LineColor(%d)'         %   self.getInt    ( 0 ) 
    elif 'LineStyle'            == name : return 'LineStyle(%d)'         %   self.getInt    ( 0 ) 
    elif 'LineWidth'            == name : return 'LineWidth(%d)'         %   self.getInt    ( 0 )     
    elif 'FillColor'            == name : return 'FillColor(%d)'         %   self.getInt    ( 0 ) 
    elif 'FillStyle'            == name : return 'FillStyle(%d)'         %   self.getInt    ( 0 ) 
    elif 'ProjectionRange'      == name : return "ProjectionRange('%s')" %   self.getString ( 0 )
    elif 'Name'                 == name : return "Name('%s')"            %   self.getString ( 0 )
    elif 'Invisible'            == name : return "Invisible()"           
    elif 'AddTo'                == name : return "AddTo('%s',%s,%s)"     % ( self.getString ( 0 ) ,
                                                                             self.getDouble ( 0 ) ,
                                                                             self.getDouble ( 1 ) ) 
    elif 'EvalErrorValue'       == name : return "EvalErrorValue(%s)"    %   self.getDouble ( 0 )
    elif 'MoveToBack'           == name : return "MoveToBack()"
    elif 'VisualizeError'       == name : return "VisializeError({.})" 
    elif 'VisualizeErrorData'   == name : return "VisializeError({.})" 
    elif 'ShowProgress'         == name : return "ShowProgress()"

    
    ## RooAbsPdf::plotOn arguments
    if   'SelectCompSet'        == name : return 'Components({.})'
    elif 'Normalization'        == name : return 'Normalization(%s,%d)'  % (  self.getDouble ( 0 ) ,
                                                                              self.getInt    ( 0 ) )  

    ## RooAbsData::plotOn arguments
    if   'CutSpec'              == name : return "Cut('%s')"             %    self.getString ( 0 )
    elif 'CutVar'               == name : return "Cut({.})"
    elif 'Binning'              == name : return "Binning({.})"
    elif 'BinningName'          == name : return "Binning('%s')"         %   self.getString ( 0 ) 
    elif 'BinningSpec'          == name : return 'Binning(%d,%s,%s)'     % ( self.getInt    ( 0 ) , self.getDouble ( 0 ) , self.getDouble ( 1 ) )
    elif 'MarkerColor'          == name : return 'MarkerColor(%d)'       %   self.getInt    ( 0 ) 
    elif 'MarkerStyle'          == name : return 'MarkerStyle(%d)'       %   self.getInt    ( 0 ) 
    elif 'MarkerSize'           == name : return 'MarkerSize (%s)'       %   self.getDouble ( 0 ) 
    elif 'CutRange'             == name : return "CutRange('%s')"        %   self.getString ( 0 ) 
    elif 'XErrorSize'           == name : return "XErrorSize(%s)"        %   self.getDouble ( 0 )
    elif 'RefreshNorm'          == name : return "RefreshNorm()"  
    elif 'Efficiency'           == name : return "Efficiency({.})"
    elif 'Rescale'              == name : return "Rescale(%s)"           %   self.getDouble ( 0 ) 


    
    ## RooDataHist::ctor arguments
    if   'Weight'               == name : return "Weight(%s)"            %   self.getDouble ( 0 )
    elif 'IndexCat'             == name : return "Index({.})"
    elif 'ImportHistoSlice'     == name : return "Import({.})"
    elif 'ImportDataHistSlice'  == name : return "Import({.})"
    elif 'ImportHisto'          == name : return "Import({.})"
    elif 'ImportDataHistSliceMany' == name : return "Import({.})"
    elif 'ImportHistoSliceMany'    == name : return "Import({.})"
    
    ## RooDataSet::ctor arguments
    if   'WeightVarName'        == name : return "WeightVar('%s',%s)"    % ( self.getString ( 0 ) ,
                                                                             self.getBool   (   ) ) 
    elif 'WeightVar'            == name : return "WeightVat({.})"
    elif 'LinkDataSlice'        == name : return "Link({.})"
    elif 'ImportDataSlice'      == name : return "Import({.})"
    elif 'ImportData'           == name : return "Import({.})"
    elif 'ImportTree'           == name : return "Import({.})"
    elif 'ImportFromFile'       == name : return "ImportFromFile('%s','%s')" % (  self.getString ( 0 ) ,
                                                                                  self.getString ( 1 ) )
    elif 'StoreError'           == name : return "StoreError({.})"
    elif 'StoreAsymError'       == name : return "StoreAsymError({.})"
    elif 'OwnLinked'            == name : return "OwnLinked()"
    elif 'ImportDataSliceMany'  == name : return "Import({.})"
    elif 'LinkDataSliceMany'    == name : return "Link({.})"


    ## RooChi2Var::ctor arguments
    if   'Extended'             == name : return "Extended(%s)"      %   self.getBool(   )
    elif 'DataError'            == name : return "DataError(%s)"     %   self.getInt ( 0 )
    elif 'NumCPU'               == name : return 'NumCPU(%d,%d)'     % ( self.getInt ( 0 ) ,
                                                                         self.getInt ( 1 ) ) 
    
    ## RooAbsCollection::printLatex arguments
    if   'Columns'              == name : return "Columns(%s)"       %   self.getInt   ( 0 ) 
    elif 'OutputFile'           == name : return "OutputFile('%s')"  %   self.getString( 0 )  
    elif 'Sibling'              == name : return "Sibling({.})"
    elif 'Format'               == name : return "Format('%s',%s)"   % ( self.getString( 0 ) ,
                                                                         self.etInt    ( 0 ) )    
    elif 'FormatArgs'           == name : return "Format({.})"
    
    
    ##  RooAbsRealLValue::frame arguments
    if   'Title'                == name : return "Title('%s')"        %   self.getString ( 0 )
    elif 'Bins'                 == name : return "Bins(%d)"           %   self.getInt    ( 0 )
    elif 'AutoRange'            == name : return "AutoRange({.})" 
    
    
    ##  RooAbsData::reduce arguments
    if   'SelectVars'           == name  : return "SelectVars({.})"
    elif 'EventRange'           == name  : return "EventRange(%d,%d)" % (  self.getInt ( 0 ) ,
                                                                           self.getInt ( 1 ) )

    ## RooAbsPdf::fitTo arguments
    if   'FitOptions'           == name : return "FitOptions('%s')"   % ( self.getString ( 0 ) )
    elif 'Optimize'             == name : return "Optimize(%s)"       %   self.getBool   (   )
    elif 'Verbose'              == name : return 'Verbose(%s)'        %   self.getBool   (   ) 
    elif 'Save'                 == name : return 'Save(%s)'           %   self.getBool   (   ) 
    elif 'Timer'                == name : return 'Timer(%s)'          %   self.getBool   (   ) 
    elif 'PrintLevel'           == name : return 'PrintLevel(%d)'     %   self.getInt    ( 0 ) 
    elif 'Warnings'             == name : return 'Warnings(%s)'       %   self.getBool   (   ) 
    elif 'Strategy'             == name : return 'Strategy(%d)'       %   self.getInt    ( 0 ) 
    elif 'InitialHesse'         == name : return 'InitialHesse(%s)'   %   self.getBool   (   ) 
    elif 'Hesse'                == name : return 'Hesse(%s)'          %   self.getBool   (   ) 
    elif 'Minos'                == name : return 'Minos(%s,{.})'      %   self.getBool   (   )
    elif 'ProjectedObservables' == name : return "ProjectedObservables({.})" 
    elif 'SplitRange'           == name : return "SplitRange(%s)"     %   self.getBool   (   )
    elif 'SumCoefRange'         == name : return "SumCoefRange('%s')" %   self.getString ( 0 ) 
    elif 'Constrain'            == name : return "Constrain({.})" 
    elif 'GlobalObservables'    == name : return "GlobalObservables({.})" 
    elif 'GlobalObservablesTag' == name : return "GlobalObservablesTag('%s')"  % self.getString ( 0 )
    elif 'Constrained'          == name : return "Constrained()"
    elif 'ExternalConstraints'  == name : return "ExternalConstraints({.})"
    elif 'PrintEvalErrors'      == name : return 'PrintEvalErrors(%d)' %  self.getInt  ( 0 ) 
    elif 'EvalErrorWall'        == name : return 'EvalErrorWall(%s)'   %  self.getBool (   ) 
    elif 'SumW2Error'           == name : return 'SumW2Error(%s)'      %   self.getBool () 
    elif 'CloneData'            == name : return 'CloneData(%s)'       %   self.getBool () 
    elif 'Integrate'            == name : return 'Integrate(%s)'       %   self.getBool () 
    elif 'Minimizer'            == name : return "Minimizer('%s','%s')"% ( self.getString ( 0 ) ,
                                                                           self.getString ( 1 ) )                                                                           
    elif 'OffsetLikelihood'     == name : return 'Offset(%s)'          %   self.getBool () 
    elif 'BatchMode'            == name : return 'BatchMode(%s)'       %   self.getBool () 
    elif 'AsymptoticError'      == name : return 'AsymptoticError(%s)' %   self.getBool () 
    
    
    ## RooAbsPdf::paramOn arguments
    if   'Label'                == name : return "Label('%s')"         %    self.getString ( 0 )
    elif 'Layout'               == name : return "Layout('%s')"        %  ( self.getDouble ( 0 ) ,
                                                                            self.getDouble ( 1 ) ,
                                                                            self.getInt    ( 0 ) / 10000.0 )
    elif 'Parameters'           == name : return "Parameters({.})"
    elif 'ShowConstants'        == name : return "ShowConstants(%s)"   %    self.getBool   ()
    
    ##  RooTreeData::statOn arguments
    if   'What'                 == name : return "What('%s')"          %    self.getString ( 0 )
    
    ## RooProdPdf::ctor arguments
    if   'Conditional'          == name : return "Conditional({.})"
    
    ## RooAbsPdf::generate arguments
    if   'PrototypeData'        == name : return "ProtoData({.})"
    elif 'NumEvents'            == name : return "NumEvents(%s)"     % self.getInt    ( 0 )
    elif 'NumEventsD'           == name : return "NumEvents(%s)"     % self.getDouble ( 0 )
    elif 'ExpectedData'         == name : return "ExpectedData(%s)"  % self.getBool   (   )
    elif 'Asimov'               == name : return "Asimov(%s)"        % self.getBool   (   )
    elif 'AutoBinned'           == name : return "AutoBinned(%s)"    % self.getBool   (   )
    elif 'GenBinned'            == name : return "GenBinned('%s')"   % self.getString ( 0  )


    ## RooAbsRealLValue::createHistogram arguments
    if   'YVar'                 == name : return "YVar({.})"
    elif 'ZVar'                 == name : return "ZVar({.})"
    elif 'AxisLabel'            == name : return "AxisLabel('%s')"     % self.getString ( 0 )
    elif 'Scaling'              == name : return "Scaling(%s)"         % self.getBool   (   )
   
    ## RooAbsReal::createHistogram arguments
    if   'IntrinsicBinning'     == name : return "IntrinsicBinning(%s)" % self.getBool  (   )
    
    ## RooAbsData::createHistogram arguments
    if   'AutoRangeData'        == name :
        if 1 ==   self.getInt ( 0 ) : 
            return  "AutoSymBining(%s,%s)" % ( self.getInt ( 1 ) , self.getDouble ( 0 ) )
        else :
            return  "AutoBining(%s,%s)"    % ( self.getInt ( 1 ) , self.getDouble ( 0 ) )
    
    ## RooAbsReal::fillHistogram arguments
    if   'IntOnb'              == name : return "IntObs({.})"

    ## RooAbsReal::createIntegral arguments
    if   'NormSet'             == name : return "NormSet({.})"
    elif 'NumIntConfig'        == name : return "NumIntConfig({.})"

    ## RooMCStudy::ctor arguments
    if   'Silence'             == name : return "Silence(%s)"        % self.getBool()
    elif 'FitModel'            == name : return "FitModel({.})"
    elif 'FitOptArgs'          == name : return "FitOptions({.})"
    elif 'Binned'              == name : return "Binned(%s)"         % self.getBool()
    elif 'BootStrapData'       == name : return "BootStrapData({.})"

    ## RooMCStudy::plot* arguments
    if   'FrameArg'            == name : return "Frame({.})"
    elif 'FrameBins'           == name : return "FrameBins(%s)"     %   self.getInt ( 0 )
    elif 'FrameRange'          == name : return "FrameRange(%s,%s)" % ( self.getDouble ( 0 ) , 
                                                                        self.getDouble ( 1 ) ) 
    elif 'FitGauss'            == name : return "FitGauss(%s)"      %   self.getBool   (   )


    ## RooRealVar::format arguments
    if   'ShowName'            == name : return "ShowName(%s)"        %   self.getBool   (   ) 
    elif 'ShowValue'           == name : return "ShowValue(%s)"       %   self.getBool   (   ) 
    elif 'ShowError'           == name : return "ShowError(%s)"       %   self.getBool   (   ) 
    elif 'ShowAsymError'       == name : return "ShowAsymError(%s)"   %   self.getBool   (   ) 
    elif 'ShowUnit'            == name : return "ShowUnit(%s)"        %   self.getBool   (   ) 
    elif 'AutoPrecision'       == name : return "AutoPrecision(%s)"   %   self.getInt    ( 0 ) 
    elif 'FixedPrecision'      == name : return "FixedPrecision(%s)"  %   self.getInt    ( 0 ) 
    elif 'TLatexStyle'         == name : return "TLatexStyle(%s)"     %   self.getBool   (   )
    elif 'LatexStyle'          == name : return "LatexStyle(%s)"      %   self.getBool   (   )
    elif 'LatexTableStyle'     == name : return "LatexTableStyle(%s)" %   self.getBool   (   )
    elif 'VerbatimName'        == name : return "VerbatimName(%s)"    %   self.getBool   (   )

    ## RooMsgService::addReportingStream arguments
    if   'Topic'               == name : return "Topic(%s)"           %   self.getInt    ( 0 )
    elif 'ObjectName'          == name : return "ObjectName('%s')"    %   self.getString ( 0 )
    elif 'ClassName'           == name : return "ClassName('%s')"     %   self.getString ( 0 )
    elif 'BaseClassName'       == name : return "BaseClassName('%s')" %   self.getString ( 0 )
    elif 'LabelName'           == name : return "TagName('%s')"       %   self.getString ( 0 )
    elif 'OutputStream'        == name : return "OutputStream({.})"
    elif 'Prefix'              == name : return "Prefix(%s)"          %   self.getBool   (   )
    elif 'Color'               == name : return "Color(%s)"           %   self.getBool   (   )
    

    ## RooWorkspace::import() arguments
    if   'RenameConflictNodes' == name : return "RenameConflictNodes('%s',%s)"     % ( self.getString ( 0 ) ,
                                                                                       self.getBool   (   ) )
    elif 'RecycleConflictNodes'== name : return "RecycleConflictNodes(%s)"         %   self.getBool   (   ) 
    elif 'RenameAllNodes'      == name : return "RenameAllNodes('%s')"             %   self.getString ( 0 ) 
    elif 'RenameAllVariables'  == name : return "RenameAllVariables('%s','%s')"    % ( self.getString ( 0 ) ,
                                                                                       self.getString ( 1 ) )
    elif 'RenameVariable'      == name : return "RenameVariable('%s','%s')"        % ( self.getString ( 0 ) ,
                                                                                       self.getString ( 1 ) )                                                                                       
    elif 'Rename'              == name : return "Renam('%s')"                      %   self.getString ( 0 ) 
    elif 'Embedded'            == name : return "Embedded(%s)"                     %   self.getBool   (   ) 
    elif 'NoRecursion'         == name : return "NoRecursion(%s)"                  %   self.getBool   (   ) 


    ## RooSimCloneTool::build() arguments
    if   'SplitParam'             == name : return "SplitParam({.})"
    elif 'SplitParamConstrained'  == name : return "SplitParamConstrained({.})"
    elif 'Restrict'               == name : return "Restrict('%s','%s')"    %  ( self.getString ( 0 ) ,
                                                                                 self.getString ( 1 ) ) 
    
    ## RooAbsPdf::createCdf() arguments
    if   'SupNormSet'             == name : return "SupNormSet({.})"
    elif 'ScanParameters'         == name : return "ScanParameters(%s,%s)"  % (  self.getInt    ( 0 ) ,
                                                                                 self.getInt    ( 1 ) ) 
    elif 'ScanNumCdf'             == name : return "ScanNumCdf()"
    elif 'ScanAllCdf'             == name : return "ScanAllCdf()"
    elif 'ScanNoCdf'              == name : return "ScanNoCdf()"


    ## ? 
    if   'MultiArg'               == name :
        return "MultiArg(%s)" % [ i for i in self.subArgs() if i.name ]
    
    return name

# ============================================================================
## flatten the list of arguments/commands  
def flat_args ( *args ) :
    """Flatten the list of arguments/commands
    """
    
    if not args : return ()
    
    flat = []
    for arg in args :

        if arg.name != "MultiArg" : flat.append ( arg )
        else :
            lst = [ i for i in arg.subArgs() if i.name ]
            flat = flat + list ( flat_args ( *lst ) ) 

    return tuple ( flat ) 


# =============================================================================
## check if the argument name matches the pattern
#  @code
#  ok =  mathch_arg ( "sumw2" , ... ) 
#  @endcode
#  - case-insenstive mathch 
#  - name starts
#  - <code>fnmatch</code> matching
#  - <code>regex</code> matching 
def match_arg ( pattern , arg ) :
    """Check if the argument name matches the pattern
    >>> ok =  mathch_arg ( 'sumw2' , ... ) 
    - case-insenstive mathhicng 
    - name starts
    - fnmatch matching
    - regex matching 
    """
    
    pl = pattern      .lower ()
    al = arg.GetName().lower () 

    if al == pl                    : return True ## 
    elif al.startswith ( pl )      : return True

    import fnmatch
    if fnmatch.fnmatch ( al , pl ) : return True
    
    ## regex
    import re
    try :
        expr = re.compile ( pattern , flags = re.I ) 
        return expr.match ( arg.GetName() ) 
    except :
        pass
    
    return False 

# =============================================================================
## Check the presense of the arg in the list
#  @code
#  arg =  check_arg ( "sumw2" , ... ) 
#  @endcode
#  - case-insenstive mathch 
#  - name starts
#  - <code>fnmatch</code> matching
#  - <code>regex</code> matching 
def check_arg  ( pattern , *args ) :
    """Check the presense of the arg in the list
    >>> arg =  check_arg ( "sumw2" , ... ) 
    - case-insenstive mathch 
    - name starts
    - <code>fnmatch</code> matching
    - <code>regex</code> matching 
    """

    flat = flat_args ( *args )
    
    for arg in flat :
        if  match_arg ( pattern , arg ) : return arg 
                
    return None

# =============================================================================
## check at least one command  different form the trivial commands  
def nontrivial_arg ( trivia , *args ) :
    """Check at least one command  different form the trivial commands
    """

    if not args : return False
    
    if isinstance ( trivia , string_types ) :
        trivia = trivia ,

    flat = flat_args ( *args )
    
    for arg in flat :

        for pattern in trivia :
            
            if match_arg ( pattern , arg ) : break 
            
        else :

            return True

    return False 

# ==============================================================================
## merge arguments into smaller chunks
#  @code
#  args = ...
#  margs = merge_args ( 3 , *args ) 
#  @endcode 
def merge_args ( num , *args ) : 
    """merge arguments into smaller chunks
    args = ...
    margs = merge_args ( 3 , *args ) 
    """
    assert isinstance ( num , integer_types ) and 1 <= num ,\
           "merge_args: invalid chunk size ``%s''" % num
    
    if len ( args ) < num : return tuple ( args ) 

    lst   = flat_args ( *args ) 
    keep  = [ l for l in lst ]
    
    while num < len ( lst ) : 
        
        nlst = chunked ( lst , 4 )
        ll   = [ ROOT.RooFit.MultiArg ( *l ) if 1 < len ( l ) else l [ 0 ] for l in nlst ] 
        for l in ll : keep.append ( l )  
        
        lst = tuple ( ll )

    return lst 

# =============================================================================
def _rca_bool_ ( self ) :
    """Get boolean value"""
    return True if self.getInt ( 0 ) else False


ROOT.RooCmdArg .__str__  = _rca_print_
ROOT.RooCmdArg .__repr__ = _rca_print_
ROOT.RooCmdArg .getBool  = _rca_bool_ 

_new_methods_ = [
    ROOT.RooCmdArg.__repr__ , 
    ROOT.RooCmdArg.__str__  , 
    ROOT.RooCmdArg.getBool  , 
    ]

# =============================================================================
_decorated_classes_ = (
    ROOT.RooCmdArg , 
    )

_new_methods_ = tuple ( _new_methods_ )

# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )
    
# =============================================================================
##                                                                      The END 
# =============================================================================
