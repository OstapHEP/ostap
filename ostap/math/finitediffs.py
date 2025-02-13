#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file  ostap/math/finitediffs.py 
#  Finite difference rules for numerrical differetiation 
#  @author Vanya BELYAEV Ivan.Belyaev@itep.ru
#  @date   2022-01-22
# =============================================================================
"""Finite difference rules for numerical differentiation 
"""
# =============================================================================
__version__ = "$Revision$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@itep.ru"
__date__    = "2022-01-22"
__all__     = (
    ##
    'Rule'              , ## abstract bse class for rules
    'darray'            , ## helper creator of array of doubles 
    'the_dot'           , ## dot-function
    ##
    'ForwardOpen'       , ## forward  (open) rule for differentiation 
    'BackwardOpen'      , ## backward (open) rule for differentiation 
    'CentralRule'       , ## central rule for differentiation 
    ##
    'ForwardOpenD1'     , ## forward  (open) rule for the 1st derivative 
    'ForwardOpenD2'     , ## forward  (open) rule for the 2nd derivative 
    'ForwardOpenD3'     , ## forward  (open) rule for the 3rd derivative 
    'ForwardOpenD4'     , ## forward  (open) rule for the 4th derivative 
    'ForwardOpenD5'     , ## forward  (open) rule for the 5th derivative 
    'ForwardOpenD6'     , ## forward  (open) rule for the 6th derivative 
    ## 
    'BackwardOpenD1'    , ## backward (open) rule for the 1st derivative 
    'BackwardOpenD2'    , ## backward (open) rule for the 2nd derivative 
    'BackwardOpenD3'    , ## backward (open) rule for the 3rd derivative 
    'BackwardOpenD4'    , ## backward (open) rule for the 4rh derivative 
    'BackwardOpenD5'    , ## backward (open) rule for the 5th derivative 
    'BackwardOpenD6'    , ## backward (open) rule for the 6th derivative 
    ##
    'CentralRuleD1'     , ## central rule for the 1st derivative 
    'CentralRuleD2'     , ## central rule for the 2nd derivative 
    'CentralRuleD3'     , ## central rule for the 3rd derivative 
    'CentralRuleD4'     , ## central rule for the 4th derivative 
    'CentralRuleD5'     , ## central rule for the 5th derivative 
    'CentralRuleD6'     , ## central rule for the 6th derivative 
    ##
    'Richardson'        , ## rule for (iterative) Richardson's extrapolation
    ##
    'DerivativeFD'      , ## generic class for evalaution of derivatives 
    'Derivative1'       , ## evaluate 1st derivative 
    'Derivative2'       , ## evaluate 2nd derivative 
    'Derivative3'       , ## evaluate 3rd derivative 
    'Derivative4'       , ## evaluate 4th derivative 
    'Derivative5'       , ## evaluate 5th derivative 
    'Derivative6'       , ## evaluate 6th derivative 
    ) 
# =============================================================================
from   collections       import namedtuple
from   ostap.math.base   import Ostap, iszero , isequal
from   ostap.math.ve     import VE
from   ostap.utils.utils import classprop, memoize 
import ROOT, math, abc, array, sys, bisect   
# =============================================================================
# logging 
# =============================================================================
from ostap.logger.logger import getLogger
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.finitediffs' )
else                       : logger = getLogger ( __name__                 )
# =============================================================================
## epsilon 
epsilon       = sys.float_info.epsilon
# =============================================================================
## Helper function \f$   f(d,n) = \epsilon^{\frac{1}{n+d}}\f$
@memoize 
def eps ( D , N ) :
    """ Helper function epsilon^(1/(n+d))
    """
    return pow ( epsilon , 1.0 / ( D + N ) )
# =============================================================================
_next_double_ = Ostap.Math.next_double
# =============================================================================
##  get "small" interval around x 
def delta ( x , ulps = 1000  ) :
    """Get ``small'' interval around x"""
    if abs ( x ) < epsilon : x = math.copysign ( 2*epsilon , x ) 
    n1 = _next_double_ ( x ,  ulps )
    n2 = _next_double_ ( x , -ulps )
    return max ( abs ( n1 - x ) , abs ( n2 - x ) )

# =============================================================================
## Four versions for dot & array settings 
## (1) use dot_fma from ostap            ## 17.7s
darray   = lambda x : array.array ( 'd' , x )
dot_fun  = Ostap.Math.dot_fma_
# 
# (2) use array and dot based on Kahan summation  ## 17.5s 
# darray  = lambda x : array.array ( 'd' , x )
# dot_fun = Ostap.Math.dot_kahan
#
# (3) use numpy:
# darray   = lambda x     : numpy.array ( x , dtype=float) 
# dot_fun  = lambda n,x,y,sx=0,sx=0 : numpy.dot(x[sx:sx+n],y[sy:sy+n])
#
# (4) generic python: use <code>math.fsum</code>, array can be arbitrary 
# dot_fun = lambda n,x,y,sx=0,sy=0 : math.fsum ( ( (i*j) for i,j,_ in zip(x[sx:],y[sy:],range(n)) ) )
# 
# =============================================================================
## dot-function
#  @param n number of elements to loop from
#  @param x the first  array
#  @param y the second array 
#  @param sx sx skip firtst <code>sx</code> elements from x 
#  @param sy sy skip firtst <code>sx</code> elements from x 
def the_dot ( n , x , y , sx = 0 , sy = 0 ) :
    """``dot''-function
    - n  : number of elements to loop from
    - x  : the first  array
    - y  : the second array 
    - sx : skip first <code>sx</code> elements from x 
    - sy : skip first <code>sx</code> elements from x 
    """
    return dot_fun ( n , x , y , sx , sy ) 
# =============================================================================
## Generator to evaluate a function seqeuntially for all points in the stencil
#  @code
#  point   = 2.0
#  step    = 0.1
#  stencil = ( -2, -1 , 0 , 1 , 2 )
#  for x , y  in fun_vals ( lambda s : s , point , step , stencil ) :
#  ... print ( x , y ) 
#  @endcode 
def fun_vals ( fun , point , step , stencil , args = () , kwargs = {} ) :
    """Generator to evaluate a function sequntially for all points in stencil
    >>> point   = 2.0
    >>> step    = 0.1
    >>> stencil = ( -2, -1 , 0 , 1 , 2 )
    >>> for x , y  in fun_vals ( lambda s : s , point , step , stencil ) :
    >>> ...  print ( fv ) 
    """
    point = float ( point )
    step  = float ( step  )
    
    for s in stencil :
        x = point + step * s 
        yield x , fun ( x , *args , **kwargs )
        
# =============================================================================
## calculate an expression \f$  \sum_i f ( x + s_i h  ) c_i \f$, where
#  @param fun the function
#  @param point the point 
#  @param step  the step
#  @param stencil the stencil  (s_i) 
#  @param coeffs  coefficients (c_i) 
def calc_dot ( fun , point , step , stencil , coeffs , args = () , kwargs = {} ) :
    """Calculate an expression   sum_i f ( x + s_i h  ) c_i 
    - fun the function
    - point the point 
    - step  the step
    - stencil the stencil  (s_i) 
    - coeffs  coefficients (c_i) 
    """
    
    point  = float ( point )
    step   = float ( step  )
    
    result = 0.0    
    for s , c in zip ( stencil , coeffs ) :
        
        if c :
            x       = point + s * step
            result += c * fun ( x , *args , **kwargs )
            
    return result


if (3,4) <= sys.version_info :
    # =========================================================================
    ## get hmax: select minimal from positive, else -1 
    def get_hmax ( *h ) :
        """get hmax : select minimal from positive, else -1"""
        return min ( ( i for i in h if 0 < i ) , default = -1 )
else :
    # =========================================================================
    ## get hmax: select minimal from positive, else -1 
    def get_hmax ( *h ) :
        """get hmax : select minimal from positive, else -1"""
        pos = [ i for i in h if 0 < i ]
        if   not pos          : return -1
        elif 1 == len ( pos ) : return pos [ 0 ] 
        return min ( pos )
    
# =============================================================================
## @class Rule
#  Define the rule for numerical differentiation (abstract class) 
class Rule(object):
    """Rule/algorithm
    Abstract base class for numeruical differentiation rule 
    """
    __metaclass__ = abc.ABCMeta

    # =========================================================================
    ## initialize the rule
    #  @param D          derivative order
    #  @param N          the order of the leading term of systematic (truncation) error
    #  @param stencil    the actual stencol for the rule
    #  @param interval   the interval of points
    #  @param with_error estmate errors ?
    #  @param hmax       maximal allowed step size 
    def __init__ ( self               ,
                   D                  ,
                   N                  ,
                   stencil    = ()    ,
                   interval   = ()    ,
                   with_error = False ,
                   max_step   = -1    ) :
        
        """Initialize the rule
        - D derivative order
        - N the order of the leading term of systematic (truncation) error
        - stencil    : the actual stencil for the rule 
        - interval   : the interval of points
        - with_error : estimate uncertainty?
        - max_step   : maximal allowed step size 
        """
        
        assert isinstance ( D , int ) and 1 <= D , \
               'Invalid derivative order D=%s/%s' % ( D , type ( D ) )
        assert isinstance ( N , int ) and 1 <= N , \
               'Invalid order of truncation error N=%s/%s' % ( N , type ( N ) )
                
        self.__D        = D
        self.__N        = N
        self.__interval = ()
        self.__with_error    = True if with_error    else False
        self.__max_step      = -1   if max_step <= 0 else float ( max_step ) 
        
        if stencil and 2 == len ( interval ) :

            self.__stencil  = stencil
            self.__interval = tuple ( interval ) 
            
        elif stencil and not interval : 
            
            l  = len ( stencil )
            
            assert D + N <= l , 'Invalid stencil lengh L=%d' % l 
            
            self.__stencil  = stencil
            self.__interval = stencil[0] , stencil[-1]
            
        elif interval and 2 == len ( interval ) : 
            
            self.__stencil  = ()
            self.__interval == tuple ( interval )
            
        else :

            raise TypeError ( "Stencil/intevals mismathch %s/%s" % ( stencil , interval ) )

        ## useful constant 
        self.__eps_DN = eps ( self.D , self.N )
        
    # =========================================================================
    ## The main method: evaluate the numerical derivative with the given step size 
    #  @param func   the function <code>func(x,*args,**kwargs)</code>
    #  @param x      point where derivative is evaluated
    #  @param h      the step 
    #  @param args   additional positional arguments
    #  @param kwargs additional kewword arguments
    @abc.abstractmethod 
    def __call__ ( self , func , x , h , args = () , kwargs = {} ) :
        """The main method: evalaute th enumerical derivarive 
        - func     the function with signature `func(x,*args,**kwargs)`
        - x        point where derivative is evaluated
        - h        the step 
        - args     additional positional arguments for the function 
        - kwargs   additional keyword arguments for th efunction 
        """
        pass

    # =========================================================================
    ## Calculate the optimal step (and f(x)) for numerical differentiation
    #  @code
    #  rule = ...
    #  hopt , f0 = rule.optimal_step ( fun , x = 0 , h = 0.1 ) 
    #  @endcode
    @abc.abstractmethod 
    def optimal_step ( self , fun , x , h , hmax , args = () , kwargs = {}) :
        """Calculate the optimal step (and f(x)) for numerical differentiation
        >>> rule = ...
        >>> hopt , f0 = rule.optimal_step ( fun , x = 0 , h = 0.1 ) 
        """
        pass 
    
    # =====================================================================================
    ## Get the value of the 1st derivative using the adaptive rule with the optimal step 
    @abc.abstractmethod 
    def derivative ( self , fun , x , h = 0 , args = () , kwargs = {} ) :
        """Get the value of the 1st derivative using the adaptive rule with the optimal step
        """
        pass
    
    # ===========================================================================
    ## the actual derivative order
    @property 
    def D  ( self ) :
        """``D'' : the actual derivative order"""
        return self.__D
    
    # =========================================================================
    ## get the order of the leading term of systematic (truncation) error 
    @property
    def N ( self ) :
        """Get the order of the leading term of systematic(truncation) error"""
        return self.__N 
    
    # ===========================================================================
    ## Get the actual stencil for the given rule
    #  @attention: it can be empty! 
    @property
    def stencil ( self ) :
        """``stencil'': The actual stencil used by the given rule (can be empty!)""" 
        return self.__stencil 

    # =========================================================================
    ## get the interval where function is evaluated
    #  (in units of the absolute value of the step size) 
    #  @code
    #  rule = ...
    #  left , right = rule.interval 
    #  @endcode
    @property
    def interval ( self ) :
        """get the interval where function is evaluated
        (in units of the absolute value of the step size)
        >>> rule = ...
        >>> left , right = rule.interval
        """
        return self.__interval 

    # ===========================================================================
    ## estimate the uncertainties?
    @property
    def with_error ( self ) :
        """``with_error'' : make error estimate?"""
        return self.__with_error
    @with_error.setter 
    def with_error ( self , value ) :
        self.__with_error = True if value else False

    # =========================================================================
    ## get the maximul allowed step size
    @property
    def max_step ( self ) :
        """``max_step'':  (if positive) maximum allowed step size"""
        return self.__max_step 
    @max_step.setter
    def max_step ( self , value ) :
        value = float ( value ) 
        self.__max_step = -1 if value <= 0 else value 

    # =========================================================================
    ## guess ``rule-of-thomb'' step size
    def h_guess ( self , x ) :
        """Guess ``rule-of-thomb'' step size
        """
        return  ( 1.0 + abs ( x ) ) * self.__eps_DN

    # =========================================================================
    ## Adjust the step size    
    def adjust_step ( self , x , h , hmax = -1 ) :
        """Adjust the step size    
        """
        x     = float ( x )  
        hh    = ( x + float ( h ) ) - x 
        
        dx    = delta ( x )
        
        ## too small ?
        h_rt      = self.h_guess ( x )  
        if iszero ( hh ) or abs ( hh ) <= dx or abs ( hh ) < 0.01 * h_rt : 
            hh    = 2 * h_rt 
            m , e = math.frexp    ( hh )
            hh    = math.copysign ( 2 ** ( e - 1 )  , m ) 
            
        ## get the current hmax 
        h_max = get_hmax ( hmax , self.max_step ) 
        
        if  dx < h_max <= abs ( hh ) :
            hh    = math.copysign ( h_max , hh )            
            m , e = math.frexp    ( hh )
            hh    = math.copysign ( 2 ** ( e - 1 )  , m )

        return hh 
 
# ============================================================================
## Context manager to switch on/off estimation of uncertainties for the given rule 
class WithError(object) :
    """Context manager to switch on/off estimation of uncertainties for the given rule
    """    
    def __init__ ( self , rule , value )  :
        self.__rule     = rule
        self.__previous = rule.with_error
        self.__value    = True if value else False
    def __enter__ ( self ) :
        self.__previous        = self.__rule.with_error
        self.__rule.with_error = self.__value
        return self
    def __exit__  ( self, *_ ) :
        self.__rule.with_error = self.__previous 

# ============================================================================
## Context manager to change the maximum allowed step size
class MaxStep (object ) :
    """Context manager to change the maximum allowed step size
    """
    def __init__ ( self , rule , value )  :
        self.__rule     = rule
        self.__previous = rule.max_step
        self.__value    = float ( value ) 
    def __enter__ ( self ) :
        self.__previous        = self.__rule.max_step 
        self.__rule.max_step   = self.__value 
        return self
    def __exit__  ( self, *_ ) :
        self.__rule.maX_step   = self.__previous 
        
        
# =============================================================================
## the basic configuration for the rule 
RuleConf = namedtuple ( 'RuleConf'  , ( 'order'        ,   ## 0
                                        'coefficients' ,   ## 1 
                                        'scale'        ,   ## 2 
                                        'accuracy'     ,   ## 3 
                                        'stencil'      ,   ## 4 
                                        'eroundoff'    ,   ## 5 
                                        'etruncation'  ,   ## 6  
                                        't_optimal'    ) ) ## 7 

# =============================================================================
## @class DiffRule
#  Helper intermediate base class to reduce code duplicaiton
class DiffRule(Rule) :
    """Helper intermediate base class to reduce code duplicaiton
    """
    def __init__ ( self                ,
                   config              ,
                   rule_err    = None  ,
                   rule_high   = None  ,
                   max_step    = -1    ) :

        assert isinstance ( config , RuleConf ) , \
               "Invalid type of ``config'' %s/%s" % ( config , type ( config ) )
        
        assert rule_err  is None or isinstance ( rule_err  , Rule ) , \
               "Invalid type for ``rule_err'' parameter %s" % rule_err

        assert 1 <= config.order , \
               "Invalid derivative order %s" % config.order 
        
        assert len ( config.stencil ) == len ( config.coefficients ) , \
               "Invalid configuration! len(stencil)!=len(coefficients!"
        
        Rule.__init__ ( self                                               ,
                        D          =         config.order                  ,
                        N          =         config.accuracy               ,
                        stencil    =         config.stencil                ,
                        interval   = ( min ( config.stencil [  0 ] , 0 ) ,
                                       max ( config.stencil [ -1 ] , 0 ) ) ,
                        with_error = True if rule_err else False           ,
                        max_step   = max_step                              )
        
        self.__config    = config
        self.__rule_err  = rule_err
        self.__rule_high = rule_high
        self.__eps_DN    = pow ( epsilon , 1.0 / ( self.D + self.N ) )

        
        ## the space for function values 
        self.__fv  = darray ( len ( self.coefficients ) * [ 0.0 ] )

    @property
    def funvals ( self ) :
        """``funvals'' : th eactual array of function values"""
        return self.__fv
    
    @property
    def config       ( self ) :
        """``config'' : th eactual configuration of the rule"""
        return self.__config
    @property
    def coefficients ( self ) :
        """``coefficients'' : the finite diffence coefficients (unscaled)"""
        return self.__config.coefficients
    @property
    def scale        ( self ) :
        """``scale'' : scale-factor for finite difference coefficients"""
        return self.__config.scale 
    @property
    def c_roundoff   ( self ) :
        """``c_roundoff'' : coefficient for caculation of the roundoff error"""
        return self.__config.eroundoff 
    @property
    def c_truncation  ( self ) :
        """``c_truncation'' : coefficient for caculation of the truncation error"""
        return self.__config.etruncation 
    @property
    def t_optimal     ( self ) :
        """``t_optimal'' : ``optimal'' value for the t-parameter of Richardson's extrapolation"""
        return self.__config.t_optimal
    @property
    def rule_err  ( self ) :
        """``rule_err''  : helper rule for evaluation of errors"""
        return self.__rule_err 
    @property
    def rule_high ( self ) :
        """``rule_high'' : helper rule for evaluation of errors"""
        return self.__rule_high

    # =========================================================================
    ## get an estimate for the roundoff error 
    def roundoff_error ( self , h , fmax = 1.0 ) :
        """Get an estimate for the roundoff error
        """
        return  epsilon * self.c_roundoff * abs ( fmax ) / pow ( h , self.D )

    # ==========================================================================
    ## Calculate the optimal step and the represenative value for abs(f) for numerical differentiation
    #  @code
    #  rule = ...
    #  hopt , fval = rule.optimal_step ( fun , x = 0 , h = 0.1 ) 
    #  @endcode
    def optimal_step ( self , fun , x , h , hmax , args = () , kwargs = {}) :
        """Calculate the optimal step and the representative value of abs(f) for numerical differentiation
        >>> rule = ...
        >>> hopt , fval = rule.optimal_step ( fun , x = 0 , h = 0.1 ) 
        """
        
        ## adjust inital step size 
        h = self.adjust_step ( x , h , hmax )
        
        ## (1) calculate a value of C using Richardson's extrapolation
        t  = self.t_optimal 
        h1 = h
        h2 = h * t 

        with WithError ( self , False ) :
            rule  = self 
            v1    = rule ( fun , x , h1 , args = args , kwargs = kwargs )
            fmax1 = max  ( self.funvals , key = abs )            
            v2    = rule ( fun , x , h2 , args = args , kwargs = kwargs )
            fmax2 = max  ( self.funvals , key = abs )
            fmax  = max  ( fmax1 , fmax2 ) 
    
        n  = self.N
        d  = self.D
        
        dv = v2 - v1
        if not dv :  dv = delta ( abs ( v1 ) ) 
        
        ## related to truncation 
        cn   = abs ( dv / ( ( 1.0 - pow ( t , n ) ) * pow ( abs ( h ) , n ) ) )
        
        ## roundoff error 
        fe   = epsilon * self.c_roundoff * abs ( fmax )
        
        ## final expression 
        hopt =   2 * d / ( n * cn ) * fe
        hopt = 1.5 * pow ( hopt , 1.0 / ( n + d ) )

        ## finbal adjustment 
        hopt = self.adjust_step ( x , hopt , hmax )
        
        return hopt , fmax 
    
    # =====================================================================================
    ## Get the value of the ``D'' derivative using the adaptive rule with the optimal step 
    def derivative ( self , fun , x , h = 0 , args = () , kwargs = {} ) :
        """Get the value of the ``D'' derivative using the adaptive rule with the optimal step
        """

        hopt , fval = self.optimal_step ( fun , x , h , self.max_step , args = args , kwargs = kwargs )
        
        return self ( fun , x , hopt , args = args , kwargs = kwargs )

# =============================================================================
## @class ForwardOpen
# Rule for the forward (open) finite differences 
class ForwardOpen (DiffRule) :
    """Rule for the forward (open) finite differences
    """
    # =========================================================================
    ## table of forward (open) rules for the finite differences 
    __FORWARD_RULES = (
        # =====================================================================
        ## 1st derivative
        # =====================================================================
        ( RuleConf ( 1 , darray ( [-1, 1] ) , 1 , 1 , darray  ( range ( 1 , 3 ) ) , 0.4082482905 , 0.66667  , 0.4142 ) ,
          RuleConf ( 1 , darray ( [-5, 8, -3] ) , 2 , 2 , darray  ( range ( 1 , 4 ) ) , 1.428869017 , -0.54545  , 0.5000 ) ,
          RuleConf ( 1 , darray ( [-26, 57, -42, 11] ) , 6 , 3 , darray  ( range ( 1 , 5 ) ) , 3.667297925 , 0.48000  , 0.5604 ) ,
          RuleConf ( 1 , darray ( [-77, 214, -234, 122, -25] ) , 12 , 4 , darray  ( range ( 1 , 6 ) ) , 8.402146441 , -0.43796  , 0.6058 ) ,
          RuleConf ( 1 , darray ( [-522, 1755, -2540, 1980, -810, 137] ) , 60 , 5 , darray  ( range ( 1 , 7 ) ) , 18.25702427 , 0.40816  , 0.6415 ) ,
          RuleConf ( 1 , darray ( [-669, 2637, -4745, 4920, -3015, 1019, -147] ) , 60 , 6 , darray  ( range ( 1 , 8 ) ) , 38.57200758 , -0.38567  , 0.6703 ) ,
          RuleConf ( 1 , darray ( [-5772, 26082, -56084, 72555, -59220, 30002, -8652, 1089] ) , 420 , 7 , darray  ( range ( 1 , 9 ) ) , 80.17366242 , 0.36794  , 0.6943 ) ,
          RuleConf ( 1 , darray ( [-13827, 70428, -176092, 272958, -278250, 187852, -81228, 20442, -2283] ) , 840 , 8 , darray  ( range ( 1 , 10 ) ) , 164.956598 , -0.35349  , 0.7145 ) ,
          RuleConf ( 1 , darray ( [-48610, 275445, -784920, 1417710, -1733004, 1461810, -842520, 317970, -71010, 7129] ) , 2520 , 9 , darray  ( range ( 1 , 11 ) ) , 337.1160279 , 0.34142  , 0.7319 ) ,
          RuleConf ( 1 , darray ( [-55991, 349255, -1117065, 2303430, -3283014, 3321822, -2392530, 1203690, -403155, 80939, -7381] ) , 2520 , 10 , darray  ( range ( 1 , 12 ) ) , 685.7317128 , -0.33114  , 0.7471 ) ,
          RuleConf ( 1 , darray ( [-699612, 4762626, -16891820, 39150045, -63737784, 75214524, -64992312, 40865220, -18247020, 5494434, -1002012, 83711] ) , 27720 , 11 , darray  ( range ( 1 , 13 ) ) , 1390.14087 , 0.32225  , 0.7603 ) ,
          RuleConf ( 1 , darray ( [-785633, 5794878, -22569206, 58074665, -106318179, 143343156, -144475716, 108993852, -60827415, 24419054, -6679398, 1115963, -86021] ) , 27720 , 12 , darray  ( range ( 1 , 14 ) ) , 2811.039273 , -0.31445  , 0.7721 ) ,
          RuleConf ( 1 , darray ( [-11359222, 90231323, -382787132, 1082724643, -2201521322, 3338354019, -3844708296, 3383444064, -2265649386, 1136832697, -414586172, 103894973, -16016182, 1145993] ) , 360360 , 13 , darray  ( range ( 1 , 15 ) ) , 5673.300608 , 0.30754  , 0.7827 ) ,
          RuleConf ( 1 , darray ( [-12530955, 106635585, -489414835, 1509235455, -3374426055, 5684163485, -7363422495, 7404831720, -5784363585, 3482642163, -1587490905, 530405785, -122643885, 17550255, -1171733] ) , 360360 , 14 , darray  ( range ( 1 , 16 ) ) , 11432.68878 , -0.30137  , 0.7922 ) ,
          RuleConf ( 1 , darray ( [-13726712, 124571940, -614969320, 2053304890, -5006634360, 9275021756, -13348186280, 15099528015, -13479059880, 9467405948, -5178349176, 2162614090, -666713320, 143104740, -19108088, 1195757] ) , 360360 , 15 , darray  ( range ( 1 , 17 ) ) , 23011.1514 , 0.29579  , 0.8008 ) ,
          RuleConf ( 1 , darray ( [-29889983, 288128824, -1522325720, 5471082820, -14447806100, 29192933224, -46208337032, 58073290990, -58316634090, 46809046856, -29868662824, 14968117892, -5767964020, 1650682520, -330603256, 41376458, -2436559] ) , 720720 , 16 , darray  ( range ( 1 , 18 ) ) , 46270.66031 , -0.29074  , 0.8086 ) ,
          RuleConf ( 1 , darray ( [-550271934, 5614607799, -31610879568, 121665119580, -345911194440, 757055940732, -1307093881392, 1806827899734, -2015860220660, 1820231237682, -1327349220912, 776010156012, -358831464264, 128360093580, -34276966992, 6434742114, -757839294, 42142223] ) , 12252240 , 17 , darray  ( range ( 1 , 19 ) ) , 92966.55294 , 0.28611  , 0.8158 ) ,
          RuleConf ( 1 , darray ( [-593094837, 6385420053, -38162783727, 156608608428, -476949277620, 1123962573636, -2102058252684, 3169623964806, -3889704810134, 3902280781542, -3201193810386, 2138806221084, -1153795835556, 495266726484, -165315050172, 41378230962, -7309743453, 812954477, -42822903] ) , 12252240 , 18 , darray  ( range ( 1 , 20 ) ) , 186664.0505 , -0.28187  , 0.8224 ) ,
          RuleConf ( 1 , darray ( [-12094689300, 137014841550, -866319635700, 3775848447825, -12263175825552, 30958707551400, -62347083656400, 101837669491350, -136326612632600, 150437160809364, -137116508357400, 103059539440650, -63536935035600, 31818044658600, -12744404605584, 3987325939050, -939170013300, 156672879950, -16505495700, 825887397] ) , 232792560 , 19 , darray  ( range ( 1 , 21 ) ) , 374588.3479 , 0.27795  , 0.8284 ) ,
          RuleConf ( 1 , darray ( [-12932216325, 153765382050, -1025449770450, 4730629256325, -16320994261677, 43943726547000, -94809631145400, 166762764469350, -241829891971850, 291108199928364, -291854651388300, 243730578559650, -169040214374850, 96743139636600, -45206952094584, 16972344934650, -4996988449425, 1111453688450, -175635630450, 17576427897, -837527025] ) , 232792560 , 20 , darray  ( range ( 1 , 22 ) ) , 751354.2731 , -0.27432  , 0.8341 ) ,
          RuleConf ( 1 , darray ( [-13780828710, 171586242135, -1203658371300, 5859283728375, -21399939385902, 61212139969365, -140858733605040, 265439412597150, -414514026195500, 540540838251414, -591173817375960, 543049744547310, -418472852697900, 269427273860250, -143883600222384, 63021447394290, -22265401871790, 6190398812675, -1304290102500, 195785028747, -18658387110, 848612385] ) , 232792560 , 21 , darray  ( range ( 1 , 23 ) ) , 1506476.362 , 0.27094  , 0.8393 ) ,
          RuleConf ( 1 , darray ( [-14640022575, 190488507165, -1402132154115, 7182442280475, -27684942508377, 83838151210275, -204965765454285, 411969771109710, -689258448406550, 967921050579714, -1146768093402750, 1149152591121990, -974067128724690, 696807486188550, -418628022433434, 209551805906850, -86372433721035, 28816410053585, -7589293224975, 1518943580847, -217132169925, 19750877415, -859193865] ) , 232792560 , 22 , darray  ( range ( 1 , 24 ) ) , 3019475.917 , -0.26779  , 0.8441 ) ,
          ) ,        
        # =====================================================================
        ## 2nd derivative
        # =====================================================================
        ( RuleConf ( 2 , darray ( [1, -2, 1] ) , 1 , 1 , darray  ( range ( 1 , 4 ) ) , 0.7071067812 , 0.50000  , 0.5961 ) ,
          RuleConf ( 2 , darray ( [3, -8, 7, -2] ) , 1 , 2 , darray  ( range ( 1 , 5 ) ) , 3.240370349 , -0.34286  , 0.6436 ) ,
          RuleConf ( 2 , darray ( [71, -236, 294, -164, 35] ) , 12 , 3 , darray  ( range ( 1 , 6 ) ) , 10.07190583 , 0.26667  , 0.6792 ) ,
          RuleConf ( 2 , darray ( [116, -461, 744, -614, 260, -45] ) , 12 , 4 , darray  ( range ( 1 , 7 ) ) , 26.63772301 , -0.22167  , 0.7071 ) ,
          RuleConf ( 2 , darray ( [2552, -11787, 23340, -25450, 16080, -5547, 812] ) , 180 , 5 , darray  ( range ( 1 , 8 ) ) , 64.70690341 , 0.19190  , 0.7298 ) ,
          RuleConf ( 2 , darray ( [3490, -18353, 43038, -58280, 48910, -25245, 7378, -938] ) , 180 , 6 , darray  ( range ( 1 , 9 ) ) , 149.4387081 , -0.17067  , 0.7486 ) ,
          RuleConf ( 2 , darray ( [127251, -750132, 2031932, -3285576, 3436650, -2360596, 1033452, -262512, 29531] ) , 5040 , 7 , darray  ( range ( 1 , 10 ) ) , 334.0934056 , 0.15472  , 0.7646 ) ,
          RuleConf ( 2 , darray ( [159826, -1043307, 3204632, -6021876, 7541100, -6465046, 3769752, -1435212, 322706, -32575] ) , 5040 , 8 , darray  ( range ( 1 , 11 ) ) , 730.5606649 , -0.14227  , 0.7784 ) ,
          RuleConf ( 2 , darray ( [976263, -6987865, 23994145, -51365340, 74903430, -76962746, 56046690, -28432020, 9584515, -1934205, 177133] ) , 25200 , 9 , darray  ( range ( 1 , 12 ) ) , 1572.4727 , 0.13225  , 0.7903 ) ,
          RuleConf ( 2 , darray ( [1166816, -9083948, 34474560, -82806585, 137785920, -164998232, 144082176, -91314510, 41025760, -12414620, 2273216, -190553] ) , 25200 , 10 , darray  ( range ( 1 , 13 ) ) , 3345.2627 , -0.12399  , 0.8009 ) ,
          RuleConf ( 2 , darray ( [45211732, -380251932, 1580309544, -4208114185, 7866803340, -10756730424, 10951798704, -8325167598, 4673718060, -1885179340, 517665192, -86769897, 6706804] ) , 831600 , 11 , darray  ( range ( 1 , 14 ) ) , 7053.445091 , 0.11706  , 0.8103 ) ,
          RuleConf ( 2 , darray ( [52315556, -472601644, 2134407816, -6239807849, 12946037500, -19899351912, 23141960688, -20515329582, 13816339548, -6964413500, 2549358856, -640868169, 99056516, -7103824] ) , 831600 , 12 , darray  ( range ( 1 , 15 ) ) , 14768.60568 , -0.11115  , 0.8187 ) ,
          RuleConf ( 2 , darray ( [5441543370, -52538338440, 256186438690, -815643823995, 1859598014274, -3173858227540, 4150444227930, -4203495912330, 3301812704190, -1996778832048, 913500257670, -306140313115, 70969470390, -10178036820, 680827774] ) , 75675600 , 13 , darray  ( range ( 1 , 16 ) ) , 30750.6277 , 0.10604  , 0.8263 ) ,
          RuleConf ( 2 , darray ( [6155179668, -63242882910, 331118249980, -1140348339585, 2833711561044, -5316908030434, 7722193899420, -8795745489960, 7894062281820, -5568528503538, 3056550060564, -1280253859885, 395673985980, -85109848110, 11385372244, -713636298] ) , 75675600 , 14 , darray  ( range ( 1 , 17 ) ) , 63737.72149 , -0.10157  , 0.8332 ) ,
          RuleConf ( 2 , darray ( [27600818349, -300653126472, 1682084961160, -6230249177460, 16758627656316, -34284707510872, 54753413811096, -69275322264720, 69930131970270, -56366454319032, 36090838455672, -18138090828676, 7006477356060, -2009295211560, 403153450216, -50536140024, 2980099677] ) , 302702400 , 15 , darray  ( range ( 1 , 18 ) ) , 131616.3567 , 0.09763  , 0.8395 ) ,
          RuleConf ( 2 , darray ( [30701312706, -353361530541, 2103752193712, -8338585340220, 24137804225976, -53470566591988, 93125131973328, -129573736519656, 145303149788940, -131739472137702, 96389252710608, -56509808990908, 26192336437176, -9388471781220, 2511489612976, -472203372576, 55688503746, -3100494357] ) , 302702400 , 16 , darray  ( range ( 1 , 19 ) ) , 270931.2609 , -0.09412  , 0.8453 ) ,
          RuleConf ( 2 , darray ( [1729792071433, -20973890279277, 132387205763643, -559112353067652, 1732944893211396, -4132366153713924, 7794344121938556, -11828196090443304, 14587871990154606, -14693614580043542, 12093263239159674, -8101935786477156, 4380771549594804, -1884179318364756, 630002847948396, -157926872717808, 27935957575377, -3110577433893, 164025123427] ) , 15437822400 , 17 , darray  ( range ( 1 , 20 ) ) , 556228.1373 , 0.09097  , 0.8506 ) ,
          ) , 
        # =====================================================================
        ## 3rd derivative
        # =====================================================================
        ( RuleConf ( 3 , darray ( [-1, 3, -3, 1] ) , 1 , 1 , darray  ( range ( 1 , 5 ) ) , 1.290994449 , 0.40000  , 0.6925 ) ,
          RuleConf ( 3 , darray ( [-7, 26, -36, 22, -5] ) , 2 , 2 , darray  ( range ( 1 , 6 ) ) , 7.260050505 , -0.23529  , 0.7221 ) ,
          RuleConf ( 3 , darray ( [-31, 137, -242, 214, -95, 17] ) , 4 , 3 , darray  ( range ( 1 , 7 ) ) , 26.35929627 , 0.16327  , 0.7454 ) ,
          RuleConf ( 3 , darray ( [-111, 568, -1219, 1408, -925, 328, -49] ) , 8 , 4 , darray  ( range ( 1 , 8 ) ) , 78.80107471 , -0.12410  , 0.7644 ) ,
          RuleConf ( 3 , darray ( [-2632, 15289, -38592, 54965, -47720, 25227, -7504, 967] ) , 120 , 5 , darray  ( range ( 1 , 9 ) ) , 211.3944508 , 0.09988  , 0.7802 ) ,
          RuleConf ( 3 , darray ( [-7667, 49802, -144468, 244498, -263650, 185022, -82292, 21158, -2403] ) , 240 , 6 , darray  ( range ( 1 , 10 ) ) , 530.001486 , -0.08357  , 0.7937 ) ,
          RuleConf ( 3 , darray ( [-663941, 4765806, -15614604, 30600654, -39405870, 34452306, -20381676, 7846074, -1779669, 180920] ) , 15120 , 7 , darray  ( range ( 1 , 11 ) ) , 1269.790623 , 0.07192  , 0.8053 ) ,
          RuleConf ( 3 , darray ( [-1748357, 13736362, -50150583, 111658308, -167111490, 174864312, -129063102, 66149148, -22480713, 4566590, -420475] ) , 30240 , 8 , darray  ( range ( 1 , 12 ) ) , 2945.993497 , -0.06320  , 0.8155 ) ,
          RuleConf ( 3 , darray ( [-11134014, 94996329, -382325510, 953009325, -1624993020, 1979531358, -1750525308, 1120181310, -507121350, 154405545, -28416894, 2392229] ) , 151200 , 9 , darray  ( range ( 1 , 13 ) ) , 6675.108392 , 0.05646  , 0.8245 ) ,
          RuleConf ( 3 , darray ( [-27624145, 254266062, -1118154742, 3084364390, -5901263955, 8201107380, -8450102724, 6482407284, -3665520615, 1487156830, -410337510, 69057862, -5356117] ) , 302400 , 10 , darray  ( range ( 1 , 14 ) ) , 14855.17186 , -0.05109  , 0.8325 ) ,
          RuleConf ( 3 , darray ( [-368973778, 3643333061, -17378140436, 52548948628, -111466254350, 174006412701, -204676771992, 183032122152, -124114958286, 62911075975, -23134652948, 5838074756, -905323666, 65108183] ) , 3326400 , 11 , darray  ( range ( 1 , 15 ) ) , 32599.20004 , 0.04672  , 0.8397 ) ,
          RuleConf ( 3 , darray ( [-440170953, 4640093511, -23857083361, 78464720328, -182734626525, 316543157051, -418481888517, 427380826752, -337920074811, 205447820325, -94403025123, 31753846456, -7384266591, 1061868633, -71197175] ) , 3326400 , 12 , darray  ( range ( 1 , 16 ) ) , 70742.57225 , -0.04310  , 0.8462 ) ,
          RuleConf ( 3 , darray ( [-235395809336, 2638012933320, -14542365629960, 51680149452295, -131080360178040, 249486567698368, -366174978008840, 420442771686795, -379738129553640, 269244476981480, -148412807671128, 62384105246645, -19338543001960, 4170542928720, -559165100440, 35118025721] ) , 1513512000 , 13 , darray  ( range ( 1 , 17 ) ) , 152128.8859 , 0.04005  , 0.8521 ) ,
          RuleConf ( 3 , darray ( [-546379942349, 6485439045472, -38155330101160, 145689760163710, -399731469448220, 829142933217872, -1337661252023096, 1705615966238470, -1732297984830270, 1403219376827840, -902136911347672, 454938008314426, -176247835096060, 50670547116560, -10188929042120, 1279649230274, -75588323677] ) , 3027024000 , 14 , darray  ( range ( 1 , 18 ) ) , 324705.2551 , -0.03744  , 0.8575 ) ,
          RuleConf ( 3 , darray ( [-627227265456, 7859843538291, -49150566043712, 200665939876470, -592148098442880, 1329426168603988, -2338227722795328, 3277934706023406, -3697696409561440, 3368617801559010, -2474455651132608, 1455504479086658, -676531070482176, 243087176111220, -65165108754880, 12274885172826, -1449992816496, 80847323107] ) , 3027024000 , 15 , darray  ( range ( 1 , 19 ) ) , 688725.1933 , 0.03519  , 0.8624 ) ,
          RuleConf ( 3 , darray ( [-713242079158, 9408110184927, -62310832540118, 270854027857302, -855353428371000, 2066401092402724, -3935006724359256, 6015270137275854, -7461532627533556, 7550658043750250, -6238291869104724, 4192839910339106, -2273310072046104, 980062099909956, -328370438683000, 82462973153658, -14610259312902, 1629113969743, -86014813702] ) , 3027024000 , 16 , darray  ( range ( 1 , 20 ) ) , 1453113.347 , -0.03323  , 0.8670 ) ,
          ) , 
        # =====================================================================
        ## 4th derivative
        # =====================================================================
        ( RuleConf ( 4 , darray ( [1, -4, 6, -4, 1] ) , 1 , 1 , darray  ( range ( 1 , 6 ) ) , 2.415229458 , 0.33333  , 0.7519 ) ,
          RuleConf ( 4 , darray ( [4, -19, 36, -34, 16, -3] ) , 1 , 2 , darray  ( range ( 1 , 7 ) ) , 16.05718946 , -0.17143  , 0.7721 ) ,
          RuleConf ( 4 , darray ( [59, -324, 741, -904, 621, -228, 35] ) , 6 , 3 , darray  ( range ( 1 , 8 ) ) , 66.55553237 , 0.10714  , 0.7885 ) ,
          RuleConf ( 4 , darray ( [115, -716, 1917, -2864, 2581, -1404, 427, -56] ) , 6 , 4 , darray  ( range ( 1 , 9 ) ) , 221.6390334 , -0.07484  , 0.8022 ) ,
          RuleConf ( 4 , darray ( [7807, -54296, 166476, -294152, 327730, -235752, 106876, -27896, 3207] ) , 240 , 5 , darray  ( range ( 1 , 10 ) ) , 650.460098 , 0.05614  , 0.8140 ) ,
          RuleConf ( 4 , darray ( [12082, -92771, 320376, -653252, 866380, -774402, 465976, -181796, 41682, -4275] ) , 240 , 6 , darray  ( range ( 1 , 11 ) ) , 1759.72136 , -0.04425  , 0.8241 ) ,
          RuleConf ( 4 , darray ( [1102859, -9261503, 35559873, -82158036, 126337470, -134893962, 101112018, -52456308, 18002151, -3686255, 341693] ) , 15120 , 7 , darray  ( range ( 1 , 12 ) ) , 4500.759925 , 0.03616  , 0.8330 ) ,
          RuleConf ( 4 , darray ( [1521002, -13861076, 58557738, -151151631, 264324660, -328076028, 294294084, -190443498, 86995746, -26684120, 4941266, -418143] ) , 15120 , 8 , darray  ( range ( 1 , 13 ) ) , 11052.99526 , -0.03037  , 0.8409 ) ,
          RuleConf ( 4 , darray ( [60566579, -595070508, 2742542394, -7820583110, 15323316705, -21672003888, 22630166076, -17543027988, 10003449285, -4086557780, 1134048234, -191782518, 14936519] ) , 453600 , 9 , darray  ( range ( 1 , 14 ) ) , 26328.37983 , 0.02606  , 0.8479 ) ,
          RuleConf ( 4 , darray ( [77975152, -821381957, 4100411088, -12799434988, 27770446400, -44076837339, 52503277344, -47416139256, 32408282736, -16533687475, 6112900112, -1549651212, 241247968, -17408573] ) , 453600 , 10 , darray  ( range ( 1 , 15 ) ) , 61253.02629 , -0.02274  , 0.8543 ) ,
          RuleConf ( 4 , darray ( [1077124482, -12106770867, 65069722678, -220654587708, 525092118210, -924079626349, 1236387674214, -1274550815736, 1015342733526, -621104977845, 286859109042, -96906966172, 22618928358, -3263063643, 219397810] ) , 4989600 , 11 , darray  ( range ( 1 , 16 ) ) , 139870.5745 , 0.02013  , 0.8600 ) ,
          RuleConf ( 4 , darray ( [1325001162, -15824921067, 91096774078, -333438477108, 863443786410, -1668453296389, 2477010457614, -2869637251536, 2610429169326, -1861727761245, 1031232779082, -435258634372, 135402817758, -29290115043, 3937548010, -247876680] ) , 4989600 , 12 , darray  ( range ( 1 , 17 ) ) , 314617.6275 , -0.01802  , 0.8652 ) ,
          ) ,
        # =====================================================================
        ## 5th derivative
        # =====================================================================
        ( RuleConf ( 5 , darray ( [-1, 5, -10, 10, -5, 1] ) , 1 , 1 , darray  ( range ( 1 , 7 ) ) , 4.582575695 , 0.28571  , 0.7922 ) ,
          RuleConf ( 5 , darray ( [-9, 52, -125, 160, -115, 44, -7] ) , 2 , 2 , darray  ( range ( 1 , 8 ) ) , 35.12477758 , -0.13043  , 0.8067 ) ,
          RuleConf ( 5 , darray ( [-73, 478, -1341, 2090, -1955, 1098, -343, 46] ) , 6 , 3 , darray  ( range ( 1 , 9 ) ) , 163.4947332 , 0.07407  , 0.8189 ) ,
          RuleConf ( 5 , darray ( [-154, 1126, -3609, 6626, -7625, 5634, -2611, 694, -81] ) , 6 , 4 , darray  ( range ( 1 , 10 ) ) , 599.7847762 , -0.04779  , 0.8293 ) ,
          RuleConf ( 5 , darray ( [-6709, 54141, -195084, 412116, -562638, 514854, -315756, 125124, -29061, 3013] ) , 144 , 5 , darray  ( range ( 1 , 11 ) ) , 1911.006797 , 0.03352  , 0.8384 ) ,
          RuleConf ( 5 , darray ( [-22009, 194192, -776763, 1855152, -2929386, 3194640, -2435622, 1281168, -444717, 91936, -8591] ) , 288 , 6 , darray  ( range ( 1 , 12 ) ) , 5549.232957 , -0.02494  , 0.8463 ) ,
          RuleConf ( 5 , darray ( [-704726, 6745939, -29651558, 78976797, -141554316, 179139534, -163200156, 106941738, -49357662, 15270191, -2848318, 242537] ) , 6048 , 7 , darray  ( range ( 1 , 13 ) ) , 15096.74656 , 0.01937  , 0.8534 ) ,
          RuleConf ( 5 , darray ( [-2033907, 20985338, -100517146, 295333694, -592213857, 852847428, -903396732, 708451836, -407820549, 167920482, -46910666, 7978534, -624455] ) , 12096 , 8 , darray  ( range ( 1 , 14 ) ) , 39146.76078 , -0.01555  , 0.8597 ) ,
          RuleConf ( 5 , darray ( [-42173356, 466421833, -2417607768, 7766124196, -17223504820, 27805245957, -33567663696, 30643490256, -21129842772, 10859104195, -4039778776, 1029528588, -161008588, 11664751] ) , 181440 , 9 , darray  ( range ( 1 , 15 ) ) , 97867.48788 , 0.01282  , 0.8654 ) ,
          RuleConf ( 5 , darray ( [-56325046, 664545493, -3705411558, 12917339356, -31389346510, 56136929337, -76065188766, 79212090336, -63627367842, 39190787575, -18205620466, 6180743748, -1448812378, 209788411, -14151690] ) , 181440 , 10 , darray  ( range ( 1 , 16 ) ) , 237772.9644 , -0.01079  , 0.8706 ) ,
          RuleConf ( 5 , darray ( [-2413491412, 30251474679, -180528895284, 678690225518, -1793102515140, 3518477644803, -5286749523748, 6183911073978, -5669615231676, 4069894284445, -2266744452060, 961218623994, -300228835244, 65173331433, -8788479180, 554764894] ) , 5987520 , 11 , darray  ( range ( 1 , 17 ) ) , 564614.1903 , 0.00924  , 0.8753 ) ,
          RuleConf ( 5 , darray ( [-3061210061, 40614973063, -258255133164, 1041412668958, -2971950456320, 6347712703635, -10473680464940, 13593812418538, -14005754244306, 11479795629005, -7453675393252, 3790453682826, -1479076776424, 427895774873, -86514717060, 10918263278, -647718649] ) , 5987520 , 12 , darray  ( range ( 1 , 18 ) ) , 1315969.908 , -0.00803  , 0.8796 ) ,
          ) ,
        # =====================================================================
        ## 6th derivative
        # =====================================================================
        ( RuleConf ( 6 , darray ( [1, -6, 15, -20, 15, -6, 1] ) , 1 , 1 , darray  ( range ( 1 , 8 ) ) , 8.774964387 , 0.25  , 0.8212 ) ,
          RuleConf ( 6 , darray ( [5, -34, 99, -160, 155, -90, 29, -4] ) , 1 , 2 , darray  ( range ( 1 , 9 ) ) , 76.13803255 , -0.1025641  , 0.8322 ) ,
          RuleConf ( 6 , darray ( [59, -448, 1488, -2824, 3350, -2544, 1208, -328, 39] ) , 4 , 3 , darray  ( range ( 1 , 10 ) ) , 393.0126032 , 0.053333333  , 0.8416 ) ,
          RuleConf ( 6 , darray ( [134, -1123, 4188, -9124, 12800, -11994, 7508, -3028, 714, -75] ) , 4 , 4 , darray  ( range ( 1 , 11 ) ) , 1574.232799 , -0.031944629  , 0.8498 ) ,
          RuleConf ( 6 , darray ( [15553, -142510, 589365, -1449000, 2345730, -2612916, 2028210, -1083240, 380925, -79630, 7513] ) , 240 , 5 , darray  ( range ( 1 , 12 ) ) , 5410.764452 , 0.021019443  , 0.8570 ) ,
          RuleConf ( 6 , darray ( [26971, -268108, 1217355, -3332970, 6113670, -7888032, 7303326, -4851180, 2264895, -707620, 133111, -11418] ) , 240 , 6 , darray  ( range ( 1 , 13 ) ) , 16786.98341 , -0.01478721  , 0.8634 ) ,
          RuleConf ( 6 , darray ( [10886713, -116643468, 576714846, -1739713060, 3565205235, -5227080696, 5619617556, -4461793992, 2595313935, -1078124860, 303485358, -51957588, 4090021] ) , 60480 , 7 , darray  ( range ( 1 , 14 ) ) , 48414.65439 , 0.010934591  , 0.8691 ) ,
          RuleConf ( 6 , darray ( [16417784, -188547391, 1008138384, -3321599366, 7519921000, -12345569073, 15110935392, -13953111828, 9713802312, -5032840625, 1885371664, -483381126, 75993944, -5531071] ) , 60480 , 8 , darray  ( range ( 1 , 15 ) ) , 132237.1236 , -0.0084029597  , 0.8743 ) ,
          RuleConf ( 6 , darray ( [23615248, -289311887, 1663107608, -5941476262, 14724582464, -26754892001, 36724919784, -38654808276, 31327786704, -19442163553, 9090033128, -3103258022, 730963168, -106295567, 7197464] ) , 60480 , 9 , darray  ( range ( 1 , 16 ) ) , 346366.0812 , 0.0066575991  , 0.8790 ) ,
          RuleConf ( 6 , darray ( [32699604, -425577227, 2616964988, -10074858242, 27124728404, -54035213069, 82192121564, -97112639136, 89785617564, -64909365333, 36370354196, -15503403962, 4864345148, -1060152947, 143462804, -9084356] ) , 60480 , 10 , darray  ( range ( 1 , 17 ) ) , 877688.447 , -0.005406885  , 0.8833 ) ,
          RuleConf ( 6 , darray ( [5792865167, -79800472988, 522621471096, -2156731053784, 6267725888308, -13582076298660, 22673311697960, -29710227868112, 30854480958378, -25459395726116, 16624838405384, -8495877496536, 3329355298516, -966789954844, 196119182808, -24823414016, 1476517439] ) , 7983360 , 11 , darray  ( range ( 1 , 18 ) ) , 2165434.499 , 0.0044816945  , 0.8873 ) ,
          RuleConf ( 6 , darray ( [7574191382, -110083018643, 764881836336, -3368032879984, 10507282280008, -24604922917080, 44719004934800, -64353460097432, 74158521245028, -68763436012766, 51268070634704, -30541570733376, 14352201916936, -5206346346544, 1407421009008, -267083779256, 31759063094, -1781326215] ) , 7983360 , 12 , darray  ( range ( 1 , 19 ) ) , 5226676.2 , -0.0037788837  , 0.8909 ) ,
          ) , 
        # =====================================================================
        )
    # =========================================================================
    ## get the overal configuriation for all known forward rules  
    @classprop
    def RULES ( cls ) :
        """Get the overal configuration for forward rules
        """
        return cls.__FORWARD_RULES
    # =========================================================================
    ## get the maximal derivative order (inclusive) 
    @classprop
    def DMAX ( cls ) :
        """get the maximal derivative order (inclusive)"""
        return len ( cls.RULES ) 

    # =========================================================================
    ## get maximal rule index for th egive derivative order 
    @classmethod
    def IMAX ( cls , D ) :
        """get maximal rule index for the given derivative order"""
        
        assert isinstance ( D , int ) and 1 <= D <= cls.DMAX , \
               'Invalid value of parameter ``D'' '        
        return len ( cls.RULES [ D - 1 ] ) 


    # =========================================================================
    def __init__ ( self , D , I = 3 , with_error = False , max_step = -1 ) :
        
        assert isinstance ( D , int ) and 1 <= D <= self.DMAX , \
               'Invalid value of parameter ``D'' '
        
        imax = len ( self.RULES [ D - 1 ]  ) 
        assert isinstance ( D , int ) and 0 <= I <  imax , \
               'Invalid value of parameter ``I'' must be 0<%s<%d ' % ( I , imax )

        ## the row from the difference table
        config = self.RULES [ D - 1 ] [ I ]

        assert D == config.order , 'Invalid configuration!'
        
        assert all ( ( 0 < s for s in config.stencil ) ) , \
               'For forward (open) differences all stencil component must be positive!' 
            
        if not with_error : rule_err = None
        elif 0 < I        : rule_err = ForwardOpen ( D = D , I = I - 1 , with_error = False , max_step = max_step )
        else              : rule_err = ForwardOpen ( D = D , I = 0     , with_error = False , max_step = max_step )
        
        # =====================================================================
        ## initialize the base class 
        DiffRule.__init__ ( self                ,
                            config              ,
                            rule_err = rule_err ,
                            max_step = -1       )
        ## index 
        self.__I = I

 
    @property
    def I  ( self ) :
        """``I'': index/row of the rule in the table"""
        return self.__I
    

    # ====================================================================================
    ## The main method: evaluate the numerical derivative 
    def __call__ ( self , func , x , h , args = () , kwargs = {}) :
        """The main method: evaluate the numerical derivative 
        - func   the function with signature `func(x,*args,**kwargs)`
        - x      point where derivative is evaluated
        - h      step size 
        - args   additional positional arguments
        - kwargs additional keyword arguments
        """

        fv = self.funvals 
        l  = len ( fv )

        for i in range ( l ) :
            fv [ i ] = func ( x + ( i + 1 ) * h , *args , **kwargs ) ## open rule! 

        d = self.D
        
        ## calculate derivatives
        coeff   = self.coefficients 
        scale   = self.scale         
        result  = the_dot ( l , fv , coeff )
        result /= scale * ( h ** d )
                
        if ( not self.with_error ) or ( not self.rule_err ) : return result       ## RETURN
        
        ## make an estimate for the roundoff error
        fmax  = max ( fv , key = abs )
        eroff = self.roundoff_error ( h , fmax ) 

        ## for this case use Richardson's extrapolation to get an estimate and the error 
        if 0 == self.I :
            ## make an error estimate using smaller steps 
            t     = self.t_optimal  
            with WithError ( self.rule_err , False ) : 
                result2 = self.rule_err ( func , x , t * h , args = args , kwargs = kwargs )
            tn    = t ** self.N 
            r      = ( result2 - tn * result ) / ( 1.0 - tn ) 
            etrunc =   r - result            
            return VE ( result , etrunc * etrunc + eroff * eroff ) 
            
        ## make an estimate for the truncation error from comparison with the lower order rule 
        
        coeff2   = self.rule_err.coefficients 
        scale2   = self.rule_err.scale 
        result2  = the_dot ( l - 1 , fv , coeff2 , 0 , 0  )
        result2 /= ( scale2 * ( h ** d ) )

        etrunc = result2 - result

        return VE ( result , etrunc * etrunc + eroff * eroff ) 

        

# =============================================================================
## @class BackwardOpen
# Rule for the backward (open) finite differences 
class BackwardOpen (DiffRule) :
    """Rule for the backward (open) finite differences
    """
    # =========================================================================
    ## table of backward (open) rules for the finite differences 
    __BACKWARD_RULES = (
        # =====================================================================
        ## 1st derivative
        # ====================================================================
        ( RuleConf ( 1 , darray ( [3, -8, 5] ) , 2 , 2 , darray  ( range ( -3 , 0 ) ) , 1.428869017 , -0.54545455  , 0.5000 ) ,
          RuleConf ( 1 , darray ( [-11, 42, -57, 26] ) , 6 , 3 , darray  ( range ( -4 , 0 ) ) , 3.667297925 , -0.48  , 0.5604 ) ,
          RuleConf ( 1 , darray ( [25, -122, 234, -214, 77] ) , 12 , 4 , darray  ( range ( -5 , 0 ) ) , 8.402146441 , -0.4379562  , 0.6058 ) ,
          RuleConf ( 1 , darray ( [-137, 810, -1980, 2540, -1755, 522] ) , 60 , 5 , darray  ( range ( -6 , 0 ) ) , 18.25702427 , -0.40816327  , 0.6415 ) ,
          RuleConf ( 1 , darray ( [147, -1019, 3015, -4920, 4745, -2637, 669] ) , 60 , 6 , darray  ( range ( -7 , 0 ) ) , 38.57200758 , -0.38567493  , 0.6703 ) ,
          RuleConf ( 1 , darray ( [-1089, 8652, -30002, 59220, -72555, 56084, -26082, 5772] ) , 420 , 7 , darray  ( range ( -8 , 0 ) ) , 80.17366242 , -0.36793693  , 0.6943 ) ,
          RuleConf ( 1 , darray ( [2283, -20442, 81228, -187852, 278250, -272958, 176092, -70428, 13827] ) , 840 , 8 , darray  ( range ( -9 , 0 ) ) , 164.956598 , -0.35348576  , 0.7145 ) ,
          RuleConf ( 1 , darray ( [-7129, 71010, -317970, 842520, -1461810, 1733004, -1417710, 784920, -275445, 48610] ) , 2520 , 9 , darray  ( range ( -10 , 0 ) ) , 337.1160279 , -0.34141715  , 0.7319 ) ,
          RuleConf ( 1 , darray ( [7381, -80939, 403155, -1203690, 2392530, -3321822, 3283014, -2303430, 1117065, -349255, 55991] ) , 2520 , 10 , darray  ( range ( -11 , 0 ) ) , 685.7317128 , -0.33113928  , 0.7471 ) ,
          RuleConf ( 1 , darray ( [-83711, 1002012, -5494434, 18247020, -40865220, 64992312, -75214524, 63737784, -39150045, 16891820, -4762626, 699612] ) , 27720 , 11 , darray  ( range ( -12 , 0 ) ) , 1390.14087 , -0.32224689  , 0.7603 ) ,
          RuleConf ( 1 , darray ( [86021, -1115963, 6679398, -24419054, 60827415, -108993852, 144475716, -143343156, 106318179, -58074665, 22569206, -5794878, 785633] ) , 27720 , 12 , darray  ( range ( -13 , 0 ) ) , 2811.039273 , -0.31445218  , 0.7721 ) ,
          RuleConf ( 1 , darray ( [-1145993, 16016182, -103894973, 414586172, -1136832697, 2265649386, -3383444064, 3844708296, -3338354019, 2201521322, -1082724643, 382787132, -90231323, 11359222] ) , 360360 , 13 , darray  ( range ( -14 , 0 ) ) , 5673.300608 , -0.30754447  , 0.7827 ) ,
          RuleConf ( 1 , darray ( [1171733, -17550255, 122643885, -530405785, 1587490905, -3482642163, 5784363585, -7404831720, 7363422495, -5684163485, 3374426055, -1509235455, 489414835, -106635585, 12530955] ) , 360360 , 14 , darray  ( range ( -15 , 0 ) ) , 11432.68878 , -0.30136558  , 0.7922 ) ,
          RuleConf ( 1 , darray ( [-1195757, 19108088, -143104740, 666713320, -2162614090, 5178349176, -9467405948, 13479059880, -15099528015, 13348186280, -9275021756, 5006634360, -2053304890, 614969320, -124571940, 13726712] ) , 360360 , 15 , darray  ( range ( -16 , 0 ) ) , 23011.1514 , -0.29579419  , 0.8008 ) ,
          RuleConf ( 1 , darray ( [2436559, -41376458, 330603256, -1650682520, 5767964020, -14968117892, 29868662824, -46809046856, 58316634090, -58073290990, 46208337032, -29192933224, 14447806100, -5471082820, 1522325720, -288128824, 29889983] ) , 720720 , 16 , darray  ( range ( -17 , 0 ) ) , 46270.66031 , -0.29073549  , 0.8086 ) ,
          RuleConf ( 1 , darray ( [-42142223, 757839294, -6434742114, 34276966992, -128360093580, 358831464264, -776010156012, 1327349220912, -1820231237682, 2015860220660, -1806827899734, 1307093881392, -757055940732, 345911194440, -121665119580, 31610879568, -5614607799, 550271934] ) , 12252240 , 17 , darray  ( range ( -18 , 0 ) ) , 92966.55294 , -0.28611418  , 0.8158 ) ,
          RuleConf ( 1 , darray ( [42822903, -812954477, 7309743453, -41378230962, 165315050172, -495266726484, 1153795835556, -2138806221084, 3201193810386, -3902280781542, 3889704810134, -3169623964806, 2102058252684, -1123962573636, 476949277620, -156608608428, 38162783727, -6385420053, 593094837] ) , 12252240 , 18 , darray  ( range ( -19 , 0 ) ) , 186664.0505 , -0.28186962  , 0.8224 ) ,
          RuleConf ( 1 , darray ( [-825887397, 16505495700, -156672879950, 939170013300, -3987325939050, 12744404605584, -31818044658600, 63536935035600, -103059539440650, 137116508357400, -150437160809364, 136326612632600, -101837669491350, 62347083656400, -30958707551400, 12263175825552, -3775848447825, 866319635700, -137014841550, 12094689300] ) , 232792560 , 19 , darray  ( range ( -20 , 0 ) ) , 374588.3479 , -0.27795237  , 0.8284 ) ,
          RuleConf ( 1 , darray ( [837527025, -17576427897, 175635630450, -1111453688450, 4996988449425, -16972344934650, 45206952094584, -96743139636600, 169040214374850, -243730578559650, 291854651388300, -291108199928364, 241829891971850, -166762764469350, 94809631145400, -43943726547000, 16320994261677, -4730629256325, 1025449770450, -153765382050, 12932216325] ) , 232792560 , 20 , darray  ( range ( -21 , 0 ) ) , 751354.2731 , -0.27432123  , 0.8341 ) ,
          RuleConf ( 1 , darray ( [-848612385, 18658387110, -195785028747, 1304290102500, -6190398812675, 22265401871790, -63021447394290, 143883600222384, -269427273860250, 418472852697900, -543049744547310, 591173817375960, -540540838251414, 414514026195500, -265439412597150, 140858733605040, -61212139969365, 21399939385902, -5859283728375, 1203658371300, -171586242135, 13780828710] ) , 232792560 , 21 , darray  ( range ( -22 , 0 ) ) , 1506476.362 , -0.27094297  , 0.8393 ) ,
          RuleConf ( 1 , darray ( [859193865, -19750877415, 217132169925, -1518943580847, 7589293224975, -28816410053585, 86372433721035, -209551805906850, 418628022433434, -696807486188550, 974067128724690, -1149152591121990, 1146768093402750, -967921050579714, 689258448406550, -411969771109710, 204965765454285, -83838151210275, 27684942508377, -7182442280475, 1402132154115, -190488507165, 14640022575] ) , 232792560 , 22 , darray  ( range ( -23 , 0 ) ) , 3019475.917 , -0.26778991  , 0.8441 ) ,
          ) ,
        # =====================================================================
        ## 2nd derivative
        # ====================================================================
        ( RuleConf ( 2 , darray ( [-2, 7, -8, 3] ) , 1 , 2 , darray  ( range ( -4 , 0 ) ) , 3.240370349 , -0.34285714  , 0.6436 ) ,
          RuleConf ( 2 , darray ( [35, -164, 294, -236, 71] ) , 12 , 3 , darray  ( range ( -5 , 0 ) ) , 10.07190583 , -0.26666667  , 0.6792 ) ,
          RuleConf ( 2 , darray ( [-45, 260, -614, 744, -461, 116] ) , 12 , 4 , darray  ( range ( -6 , 0 ) ) , 26.63772301 , -0.22167488  , 0.7071 ) ,
          RuleConf ( 2 , darray ( [812, -5547, 16080, -25450, 23340, -11787, 2552] ) , 180 , 5 , darray  ( range ( -7 , 0 ) ) , 64.70690341 , -0.19189765  , 0.7298 ) ,
          RuleConf ( 2 , darray ( [-938, 7378, -25245, 48910, -58280, 43038, -18353, 3490] ) , 180 , 6 , darray  ( range ( -8 , 0 ) ) , 149.4387081 , -0.17066811  , 0.7486 ) ,
          RuleConf ( 2 , darray ( [29531, -262512, 1033452, -2360596, 3436650, -3285576, 2031932, -750132, 127251] ) , 5040 , 7 , darray  ( range ( -9 , 0 ) ) , 334.0934056 , -0.15471988  , 0.7646 ) ,
          RuleConf ( 2 , darray ( [-32575, 322706, -1435212, 3769752, -6465046, 7541100, -6021876, 3204632, -1043307, 159826] ) , 5040 , 8 , darray  ( range ( -10 , 0 ) ) , 730.5606649 , -0.14226598  , 0.7784 ) ,
          RuleConf ( 2 , darray ( [177133, -1934205, 9584515, -28432020, 56046690, -76962746, 74903430, -51365340, 23994145, -6987865, 976263] ) , 25200 , 9 , darray  ( range ( -11 , 0 ) ) , 1572.4727 , -0.13224667  , 0.7903 ) ,
          RuleConf ( 2 , darray ( [-190553, 2273216, -12414620, 41025760, -91314510, 144082176, -164998232, 137785920, -82806585, 34474560, -9083948, 1166816] ) , 25200 , 10 , darray  ( range ( -12 , 0 ) ) , 3345.2627 , -0.12399348  , 0.8009 ) ,
          RuleConf ( 2 , darray ( [6706804, -86769897, 517665192, -1885179340, 4673718060, -8325167598, 10951798704, -10756730424, 7866803340, -4208114185, 1580309544, -380251932, 45211732] ) , 831600 , 11 , darray  ( range ( -13 , 0 ) ) , 7053.445091 , -0.11706371  , 0.8103 ) ,
          RuleConf ( 2 , darray ( [-7103824, 99056516, -640868169, 2549358856, -6964413500, 13816339548, -20515329582, 23141960688, -19899351912, 12946037500, -6239807849, 2134407816, -472601644, 52315556] ) , 831600 , 12 , darray  ( range ( -14 , 0 ) ) , 14768.60568 , -0.11115234  , 0.8187 ) ,
          RuleConf ( 2 , darray ( [680827774, -10178036820, 70969470390, -306140313115, 913500257670, -1996778832048, 3301812704190, -4203495912330, 4150444227930, -3173858227540, 1859598014274, -815643823995, 256186438690, -52538338440, 5441543370] ) , 75675600 , 13 , darray  ( range ( -15 , 0 ) ) , 30750.6277 , -0.10604225  , 0.8263 ) ,
          RuleConf ( 2 , darray ( [-713636298, 11385372244, -85109848110, 395673985980, -1280253859885, 3056550060564, -5568528503538, 7894062281820, -8795745489960, 7722193899420, -5316908030434, 2833711561044, -1140348339585, 331118249980, -63242882910, 6155179668] ) , 75675600 , 14 , darray  ( range ( -16 , 0 ) ) , 63737.72149 , -0.10157459  , 0.8332 ) ,
          RuleConf ( 2 , darray ( [2980099677, -50536140024, 403153450216, -2009295211560, 7006477356060, -18138090828676, 36090838455672, -56366454319032, 69930131970270, -69275322264720, 54753413811096, -34284707510872, 16758627656316, -6230249177460, 1682084961160, -300653126472, 27600818349] ) , 302702400 , 15 , darray  ( range ( -17 , 0 ) ) , 131616.3567 , -0.097630366  , 0.8395 ) ,
          RuleConf ( 2 , darray ( [-3100494357, 55688503746, -472203372576, 2511489612976, -9388471781220, 26192336437176, -56509808990908, 96389252710608, -131739472137702, 145303149788940, -129573736519656, 93125131973328, -53470566591988, 24137804225976, -8338585340220, 2103752193712, -353361530541, 30701312706] ) , 302702400 , 16 , darray  ( range ( -18 , 0 ) ) , 270931.2609 , -0.094118645  , 0.8453 ) ,
          RuleConf ( 2 , darray ( [164025123427, -3110577433893, 27935957575377, -157926872717808, 630002847948396, -1884179318364756, 4380771549594804, -8101935786477156, 12093263239159674, -14693614580043542, 14587871990154606, -11828196090443304, 7794344121938556, -4132366153713924, 1732944893211396, -559112353067652, 132387205763643, -20973890279277, 1729792071433] ) , 15437822400 , 17 , darray  ( range ( -19 , 0 ) ) , 556228.1373 , -0.090968689  , 0.8506 ) ,
          ) ,
        # =====================================================================
        ## 3rd derivative
        # ====================================================================
        ( RuleConf ( 3 , darray ( [5, -22, 36, -26, 7] ) , 2 , 2 , darray  ( range ( -5 , 0 ) ) , 7.260050505 , -0.23529412  , 0.7221 ) ,
          RuleConf ( 3 , darray ( [-17, 95, -214, 242, -137, 31] ) , 4 , 3 , darray  ( range ( -6 , 0 ) ) , 26.35929627 , -0.16326531  , 0.7454 ) ,
          RuleConf ( 3 , darray ( [49, -328, 925, -1408, 1219, -568, 111] ) , 8 , 4 , darray  ( range ( -7 , 0 ) ) , 78.80107471 , -0.12409514  , 0.7644 ) ,
          RuleConf ( 3 , darray ( [-967, 7504, -25227, 47720, -54965, 38592, -15289, 2632] ) , 120 , 5 , darray  ( range ( -8 , 0 ) ) , 211.3944508 , -0.099875156  , 0.7802 ) ,
          RuleConf ( 3 , darray ( [2403, -21158, 82292, -185022, 263650, -244498, 144468, -49802, 7667] ) , 240 , 6 , darray  ( range ( -9 , 0 ) ) , 530.001486 , -0.08357285  , 0.7937 ) ,
          RuleConf ( 3 , darray ( [-180920, 1779669, -7846074, 20381676, -34452306, 39405870, -30600654, 15614604, -4765806, 663941] ) , 15120 , 7 , darray  ( range ( -10 , 0 ) ) , 1269.790623 , -0.071918663  , 0.8053 ) ,
          RuleConf ( 3 , darray ( [420475, -4566590, 22480713, -66149148, 129063102, -174864312, 167111490, -111658308, 50150583, -13736362, 1748357] ) , 30240 , 8 , darray  ( range ( -11 , 0 ) ) , 2945.993497 , -0.063204651  , 0.8155 ) ,
          RuleConf ( 3 , darray ( [-2392229, 28416894, -154405545, 507121350, -1120181310, 1750525308, -1979531358, 1624993020, -953009325, 382325510, -94996329, 11134014] ) , 151200 , 9 , darray  ( range ( -12 , 0 ) ) , 6675.108392 , -0.056458811  , 0.8245 ) ,
          RuleConf ( 3 , darray ( [5356117, -69057862, 410337510, -1487156830, 3665520615, -6482407284, 8450102724, -8201107380, 5901263955, -3084364390, 1118154742, -254266062, 27624145] ) , 302400 , 10 , darray  ( range ( -13 , 0 ) ) , 14855.17186 , -0.051090352  , 0.8325 ) ,
          RuleConf ( 3 , darray ( [-65108183, 905323666, -5838074756, 23134652948, -62911075975, 124114958286, -183032122152, 204676771992, -174006412701, 111466254350, -52548948628, 17378140436, -3643333061, 368973778] ) , 3326400 , 11 , darray  ( range ( -14 , 0 ) ) , 32599.20004 , -0.046720955  , 0.8397 ) ,
          RuleConf ( 3 , darray ( [71197175, -1061868633, 7384266591, -31753846456, 94403025123, -205447820325, 337920074811, -427380826752, 418481888517, -316543157051, 182734626525, -78464720328, 23857083361, -4640093511, 440170953] ) , 3326400 , 12 , darray  ( range ( -15 , 0 ) ) , 70742.57225 , -0.043097867  , 0.8462 ) ,
          RuleConf ( 3 , darray ( [-35118025721, 559165100440, -4170542928720, 19338543001960, -62384105246645, 148412807671128, -269244476981480, 379738129553640, -420442771686795, 366174978008840, -249486567698368, 131080360178040, -51680149452295, 14542365629960, -2638012933320, 235395809336] ) , 1513512000 , 13 , darray  ( range ( -16 , 0 ) ) , 152128.8859 , -0.040046185  , 0.8521 ) ,
          RuleConf ( 3 , darray ( [75588323677, -1279649230274, 10188929042120, -50670547116560, 176247835096060, -454938008314426, 902136911347672, -1403219376827840, 1732297984830270, -1705615966238470, 1337661252023096, -829142933217872, 399731469448220, -145689760163710, 38155330101160, -6485439045472, 546379942349] ) , 3027024000 , 14 , darray  ( range ( -17 , 0 ) ) , 324705.2551 , -0.03744124  , 0.8575 ) ,
          RuleConf ( 3 , darray ( [-80847323107, 1449992816496, -12274885172826, 65165108754880, -243087176111220, 676531070482176, -1455504479086658, 2474455651132608, -3368617801559010, 3697696409561440, -3277934706023406, 2338227722795328, -1329426168603988, 592148098442880, -200665939876470, 49150566043712, -7859843538291, 627227265456] ) , 3027024000 , 15 , darray  ( range ( -18 , 0 ) ) , 688725.1933 , -0.035191892  , 0.8624 ) ,
          RuleConf ( 3 , darray ( [86014813702, -1629113969743, 14610259312902, -82462973153658, 328370438683000, -980062099909956, 2273310072046104, -4192839910339106, 6238291869104724, -7550658043750250, 7461532627533556, -6015270137275854, 3935006724359256, -2066401092402724, 855353428371000, -270854027857302, 62310832540118, -9408110184927, 713242079158] ) , 3027024000 , 16 , darray  ( range ( -19 , 0 ) ) , 1453113.347 , -0.033230045  , 0.8670 ) ,
          ) ,        
        # =====================================================================
        ## 4th derivative
        # ====================================================================
        ( RuleConf ( 4 , darray ( [-3, 16, -34, 36, -19, 4] ) , 1 , 2 , darray  ( range ( -6 , 0 ) ) , 16.05718946 , -0.17142857  , 0.7721 ) ,
          RuleConf ( 4 , darray ( [35, -228, 621, -904, 741, -324, 59] ) , 6 , 3 , darray  ( range ( -7 , 0 ) ) , 66.55553237 , -0.10714286  , 0.7885 ) ,
          RuleConf ( 4 , darray ( [-56, 427, -1404, 2581, -2864, 1917, -716, 115] ) , 6 , 4 , darray  ( range ( -8 , 0 ) ) , 221.6390334 , -0.074836296  , 0.8022 ) ,
          RuleConf ( 4 , darray ( [3207, -27896, 106876, -235752, 327730, -294152, 166476, -54296, 7807] ) , 240 , 5 , darray  ( range ( -9 , 0 ) ) , 650.460098 , -0.056140351  , 0.8140 ) ,
          RuleConf ( 4 , darray ( [-4275, 41682, -181796, 465976, -774402, 866380, -653252, 320376, -92771, 12082] ) , 240 , 6 , darray  ( range ( -10 , 0 ) ) , 1759.72136 , -0.044250248  , 0.8241 ) ,
          RuleConf ( 4 , darray ( [341693, -3686255, 18002151, -52456308, 101112018, -134893962, 126337470, -82158036, 35559873, -9261503, 1102859] ) , 15120 , 7 , darray  ( range ( -11 , 0 ) ) , 4500.759925 , -0.036159878  , 0.8330 ) ,
          RuleConf ( 4 , darray ( [-418143, 4941266, -26684120, 86995746, -190443498, 294294084, -328076028, 264324660, -151151631, 58557738, -13861076, 1521002] ) , 15120 , 8 , darray  ( range ( -12 , 0 ) ) , 11052.99526 , -0.030368522  , 0.8409 ) ,
          RuleConf ( 4 , darray ( [14936519, -191782518, 1134048234, -4086557780, 10003449285, -17543027988, 22630166076, -21672003888, 15323316705, -7820583110, 2742542394, -595070508, 60566579] ) , 453600 , 9 , darray  ( range ( -13 , 0 ) ) , 26328.37983 , -0.026056128  , 0.8479 ) ,
          RuleConf ( 4 , darray ( [-17408573, 241247968, -1549651212, 6112900112, -16533687475, 32408282736, -47416139256, 52503277344, -44076837339, 27770446400, -12799434988, 4100411088, -821381957, 77975152] ) , 453600 , 10 , darray  ( range ( -14 , 0 ) ) , 61253.02629 , -0.022742251  , 0.8543 ) ,
          RuleConf ( 4 , darray ( [219397810, -3263063643, 22618928358, -96906966172, 286859109042, -621104977845, 1015342733526, -1274550815736, 1236387674214, -924079626349, 525092118210, -220654587708, 65069722678, -12106770867, 1077124482] ) , 4989600 , 11 , darray  ( range ( -15 , 0 ) ) , 139870.5745 , -0.020129364  , 0.8600 ) ,
          RuleConf ( 4 , darray ( [-247876680, 3937548010, -29290115043, 135402817758, -435258634372, 1031232779082, -1861727761245, 2610429169326, -2869637251536, 2477010457614, -1668453296389, 863443786410, -333438477108, 91096774078, -15824921067, 1325001162] ) , 4989600 , 12 , darray  ( range ( -16 , 0 ) ) , 314617.6275 , -0.018024702  , 0.8652 ) ,
          ) ,
        # =====================================================================
        ## 5th derivative
        # ====================================================================
        ( RuleConf ( 5 , darray ( [7, -44, 115, -160, 125, -52, 9] ) , 2 , 2 , darray  ( range ( -7 , 0 ) ) , 35.12477758 , -0.13043478  , 0.8067 ) ,
          RuleConf ( 5 , darray ( [-46, 343, -1098, 1955, -2090, 1341, -478, 73] ) , 6 , 3 , darray  ( range ( -8 , 0 ) ) , 163.4947332 , -0.074074074  , 0.8189 ) ,
          RuleConf ( 5 , darray ( [81, -694, 2611, -5634, 7625, -6626, 3609, -1126, 154] ) , 6 , 4 , darray  ( range ( -9 , 0 ) ) , 599.7847762 , -0.047792897  , 0.8293 ) ,
          RuleConf ( 5 , darray ( [-3013, 29061, -125124, 315756, -514854, 562638, -412116, 195084, -54141, 6709] ) , 144 , 5 , darray  ( range ( -10 , 0 ) ) , 1911.006797 , -0.033523455  , 0.8384 ) ,
          RuleConf ( 5 , darray ( [8591, -91936, 444717, -1281168, 2435622, -3194640, 2929386, -1855152, 776763, -194192, 22009] ) , 288 , 6 , darray  ( range ( -11 , 0 ) ) , 5549.232957 , -0.024936401  , 0.8463 ) ,
          RuleConf ( 5 , darray ( [-242537, 2848318, -15270191, 49357662, -106941738, 163200156, -179139534, 141554316, -78976797, 29651558, -6745939, 704726] ) , 6048 , 7 , darray  ( range ( -12 , 0 ) ) , 15096.74656 , -0.019370491  , 0.8534 ) ,
          RuleConf ( 5 , darray ( [624455, -7978534, 46910666, -167920482, 407820549, -708451836, 903396732, -852847428, 592213857, -295333694, 100517146, -20985338, 2033907] ) , 12096 , 8 , darray  ( range ( -13 , 0 ) ) , 39146.76078 , -0.015554554  , 0.8597 ) ,
          RuleConf ( 5 , darray ( [-11664751, 161008588, -1029528588, 4039778776, -10859104195, 21129842772, -30643490256, 33567663696, -27805245957, 17223504820, -7766124196, 2417607768, -466421833, 42173356] ) , 181440 , 9 , darray  ( range ( -14 , 0 ) ) , 97867.48788 , -0.012821084  , 0.8654 ) ,
          RuleConf ( 5 , darray ( [14151690, -209788411, 1448812378, -6180743748, 18205620466, -39190787575, 63627367842, -79212090336, 76065188766, -56136929337, 31389346510, -12917339356, 3705411558, -664545493, 56325046] ) , 181440 , 10 , darray  ( range ( -15 , 0 ) ) , 237772.9644 , -0.010792896  , 0.8706 ) ,
          RuleConf ( 5 , darray ( [-554764894, 8788479180, -65173331433, 300228835244, -961218623994, 2266744452060, -4069894284445, 5669615231676, -6183911073978, 5286749523748, -3518477644803, 1793102515140, -678690225518, 180528895284, -30251474679, 2413491412] ) , 5987520 , 11 , darray  ( range ( -16 , 0 ) ) , 564614.1903 , -0.0092440136  , 0.8753 ) ,
          RuleConf ( 5 , darray ( [647718649, -10918263278, 86514717060, -427895774873, 1479076776424, -3790453682826, 7453675393252, -11479795629005, 14005754244306, -13593812418538, 10473680464940, -6347712703635, 2971950456320, -1041412668958, 258255133164, -40614973063, 3061210061] ) , 5987520 , 12 , darray  ( range ( -17 , 0 ) ) , 1315969.908 , -0.0080324127  , 0.8796 ) ,
          ) , 
        # =====================================================================
        ## 6th derivative
        # ====================================================================
        ( RuleConf ( 6 , darray ( [-4, 29, -90, 155, -160, 99, -34, 5] ) , 1 , 2 , darray  ( range ( -8 , 0 ) ) , 76.13803255 , -0.1025641  , 0.8322 ) ,
          RuleConf ( 6 , darray ( [39, -328, 1208, -2544, 3350, -2824, 1488, -448, 59] ) , 4 , 3 , darray  ( range ( -9 , 0 ) ) , 393.0126032 , -0.053333333  , 0.8416 ) ,
          RuleConf ( 6 , darray ( [-75, 714, -3028, 7508, -11994, 12800, -9124, 4188, -1123, 134] ) , 4 , 4 , darray  ( range ( -10 , 0 ) ) , 1574.232799 , -0.031944629  , 0.8498 ) ,
          RuleConf ( 6 , darray ( [7513, -79630, 380925, -1083240, 2028210, -2612916, 2345730, -1449000, 589365, -142510, 15553] ) , 240 , 5 , darray  ( range ( -11 , 0 ) ) , 5410.764452 , -0.021019443  , 0.8570 ) ,
          RuleConf ( 6 , darray ( [-11418, 133111, -707620, 2264895, -4851180, 7303326, -7888032, 6113670, -3332970, 1217355, -268108, 26971] ) , 240 , 6 , darray  ( range ( -12 , 0 ) ) , 16786.98341 , -0.01478721  , 0.8634 ) ,
          RuleConf ( 6 , darray ( [4090021, -51957588, 303485358, -1078124860, 2595313935, -4461793992, 5619617556, -5227080696, 3565205235, -1739713060, 576714846, -116643468, 10886713] ) , 60480 , 7 , darray  ( range ( -13 , 0 ) ) , 48414.65439 , -0.010934591  , 0.8691 ) ,
          RuleConf ( 6 , darray ( [-5531071, 75993944, -483381126, 1885371664, -5032840625, 9713802312, -13953111828, 15110935392, -12345569073, 7519921000, -3321599366, 1008138384, -188547391, 16417784] ) , 60480 , 8 , darray  ( range ( -14 , 0 ) ) , 132237.1236 , -0.0084029597  , 0.8743 ) ,
          RuleConf ( 6 , darray ( [7197464, -106295567, 730963168, -3103258022, 9090033128, -19442163553, 31327786704, -38654808276, 36724919784, -26754892001, 14724582464, -5941476262, 1663107608, -289311887, 23615248] ) , 60480 , 9 , darray  ( range ( -15 , 0 ) ) , 346366.0812 , -0.0066575991  , 0.8790 ) ,
          RuleConf ( 6 , darray ( [-9084356, 143462804, -1060152947, 4864345148, -15503403962, 36370354196, -64909365333, 89785617564, -97112639136, 82192121564, -54035213069, 27124728404, -10074858242, 2616964988, -425577227, 32699604] ) , 60480 , 10 , darray  ( range ( -16 , 0 ) ) , 877688.447 , -0.005406885  , 0.8833 ) ,
          RuleConf ( 6 , darray ( [1476517439, -24823414016, 196119182808, -966789954844, 3329355298516, -8495877496536, 16624838405384, -25459395726116, 30854480958378, -29710227868112, 22673311697960, -13582076298660, 6267725888308, -2156731053784, 522621471096, -79800472988, 5792865167] ) , 7983360 , 11 , darray  ( range ( -17 , 0 ) ) , 2165434.499 , -0.0044816945  , 0.8873 ) ,
          RuleConf ( 6 , darray ( [-1781326215, 31759063094, -267083779256, 1407421009008, -5206346346544, 14352201916936, -30541570733376, 51268070634704, -68763436012766, 74158521245028, -64353460097432, 44719004934800, -24604922917080, 10507282280008, -3368032879984, 764881836336, -110083018643, 7574191382] ) , 7983360 , 12 , darray  ( range ( -18 , 0 ) ) , 5226676.2 , -0.0037788837  , 0.8909 ) ,
          ) , 
        # ====================================================================
        )
        
    
    # =========================================================================
    ## get the overall configurition for all known backward rules  
    @classprop
    def RULES ( cls ) :
        """Get the overal configuration for all known backward rules
        """
        return cls.__BACKWARD_RULES
    # =========================================================================
    ## get the maximal derivative order (inclusive) 
    @classprop
    def DMAX ( cls ) :
        """get the maximal derivative order (inclusive)"""
        return len ( cls.RULES ) 
        
    
    # =========================================================================
    ## get maximal rule index for the given derivative order 
    @classmethod
    def IMAX ( cls , D ) :
        """get maximal rule index for the given derivative order"""
        
        assert isinstance ( D , int ) and 1 <= D <= cls.DMAX , \
               'Invalid value of parameter ``D'' '        
        return len ( cls.RULES [ D - 1 ] ) 


    # =========================================================================
    def __init__ ( self , D , I = 3 , with_error = False , max_step = -1 ) :
        
        assert isinstance ( D , int ) and 1 <= D <= self.DMAX , \
               'Invalid value of parameter ``D'' '
        
        imax = len ( self.RULES [ D - 1 ]  ) 
        assert isinstance ( D , int ) and 0 <= I <  imax , \
               'Invalid value of parameter ``I'' must be 0<%s<%d ' % ( I , imax )
        
        ## the row from the difference table
        config = self.RULES [ D - 1 ] [ I ]
        
        assert D == config.order , 'Invalid configuration!'
        
        assert all ( ( s < 0 for s in config.stencil ) ) , \
               'For backward (open) differences all stencil component must be negative!' 
        
        if not with_error : rule_err = None
        elif 0 < I        : rule_err = BackwardOpen ( D = D , I = I - 1 , with_error = False , max_step = max_step )
        else              : rule_err = BackwardOpen ( D = D , I = 0     , with_error = False , max_step = max_step )
        
        # =====================================================================
        ## initialize the base class 
        DiffRule.__init__ ( self                ,
                            config              ,
                            rule_err = rule_err ,
                            max_step = max_step )
        ## index 
        self.__I = I

            
    @property
    def I  ( self ) :
        """``I'': index/row of the rule in the table"""
        return self.__I


    # ====================================================================================
    ## The main method: evaluate the numerical derivative 
    def __call__ ( self , func , x , h , args = () , kwargs = {}) :
        """The main method: evaluate the numerical derivative 
        - func   the function with signature `func(x,*args,**kwargs)`
        - x      point where derivative is evaluated
        - h      step size 
        - args   additional positional arguments
        - kwargs additional keyword arguments
        """
        
        fv = self.funvals 
        l  = len ( fv )

        for i in range ( l ) :
            ii = i + 1 
            fv [ l - ii ] = func ( x - ii * h , *args , **kwargs ) ## open rule! 
            
        d = self.D
        
        ## calculate derivatives
        coeff   = self.coefficients 
        scale   = self.scale         
        result  = the_dot ( l , fv , coeff )
        result /= scale * ( h ** d )
                
        if ( not self.with_error ) or ( not self.rule_err ) : return result       ## RETURN
        
        ## make an estimate for the roundoff error
        fmax  = max ( fv , key = abs )
        eroff = self.roundoff_error ( h , fmax ) 

        ## for this case use Richardson's extrapolation to get an estimate and the error 
        if 0 == self.I :
            ## make an error estimate using smaller steps 
            t     = self.t_optimal  
            with WithError ( self.rule_err , False ) : 
                result2 = self.rule_err ( func , x , t * h , args = args , kwargs = kwargs )
            tn    = t ** self.N 
            r      = ( result2 - tn * result ) / ( 1.0 - tn ) 
            etrunc =   r - result            
            return VE ( result , etrunc * etrunc + eroff * eroff ) 
            
        ## make an estimate for the truncation error from comparison with the lower order rule 
        
        coeff2   = self.rule_err.coefficients 
        scale2   = self.rule_err.scale 
        result2  = the_dot ( l - 1 , fv , coeff2 , 1 , 0  )
        result2 /= ( scale2 * ( h ** d ) )

        etrunc = result2 - result

        return VE ( result , etrunc * etrunc + eroff * eroff ) 

        
# =============================================================================
## @class Central
# Rule for the central finite differences 
class CentralRule(DiffRule) :
    """Rule for the central finite differences
    """
    # =========================================================================
    __CENTRAL_RULES = (
        # =====================================================================
        ## 1st derivative
        # =====================================================================
        ( RuleConf ( 1 , darray ( [-1, 0, 1] ) , 2 , 2 , darray ( range ( -1 , 2 ) ) , 0.2041241452 , 6.0 , 0.5000 ) , 
          RuleConf ( 1 , darray ( [1, -8, 0, 8, -1] ) , 12 , 4 , darray ( range ( -2 , 3 ) ) , 0.2742835786 , -30.0 , 0.6058 ) , 
          RuleConf ( 1 , darray ( [-1, 9, -45, 0, 45, -9, 1] ) , 60 , 6 , darray ( range ( -3 , 4 ) ) , 0.3123240245 , 140.0 , 0.6703 ) , 
          RuleConf ( 1 , darray ( [3, -32, 168, -672, 0, 672, -168, 32, -3] ) , 840 , 8 , darray ( range ( -4 , 5 ) ) , 0.3370123643 , -630.0 , 0.7145 ) , 
          RuleConf ( 1 , darray ( [-2, 25, -150, 600, -2100, 0, 2100, -600, 150, -25, 2] ) , 2520 , 10 , darray ( range ( -5 , 6 ) ) , 0.3546772993 , 2772.0 , 0.7471 ) , 
          RuleConf ( 1 , darray ( [5, -72, 495, -2200, 7425, -23760, 0, 23760, -7425, 2200, -495, 72, -5] ) , 27720 , 12 , darray ( range ( -6 , 7 ) ) , 0.368118142 , -12012.0 , 0.7721 ) , 
          RuleConf ( 1 , darray ( [-15, 245, -1911, 9555, -35035, 105105, -315315, 0, 315315, -105105, 35035, -9555, 1911, -245, 15] ) , 360360 , 14 , darray ( range ( -7 , 8 ) ) , 0.3787871531 , 51480.0 , 0.7922 ) , 
          RuleConf ( 1 , darray ( [7, -128, 1120, -6272, 25480, -81536, 224224, -640640, 0, 640640, -224224, 81536, -25480, 6272, -1120, 128, -7] ) , 720720 , 16 , darray ( range ( -8 , 9 ) ) , 0.3875221123 , -218790.0 , 0.8086 ) , 
          RuleConf ( 1 , darray ( [-28, 567, -5508, 34272, -154224, 539784, -1559376, 4009824, -11027016, 0, 11027016, -4009824, 1559376, -539784, 154224, -34272, 5508, -567, 28] ) , 12252240 , 18 , darray ( range ( -9 , 10 ) ) , 0.3948445227 , 923780.0 , 0.8224 ) , 
          RuleConf ( 1 , darray ( [126, -2800, 29925, -205200, 1017450, -3907008, 12209400, -32558400, 79361100, -211629600, 0, 211629600, -79361100, 32558400, -12209400, 3907008, -1017450, 205200, -29925, 2800, -126] ) , 232792560 , 20 , darray ( range ( -10 , 11 ) ) , 0.401098173 , -3879876.0 , 0.8341 ) , 
          RuleConf ( 1 , darray ( [-30, 726, -8470, 63525, -344850, 1448370, -4924458, 14069880, -35174700, 82074300, -213393180, 0, 213393180, -82074300, 35174700, -14069880, 4924458, -1448370, 344850, -63525, 8470, -726, 30] ) , 232792560 , 22 , darray ( range ( -11 , 12 ) ) , 0.4065200913 , 16224936.0 , 0.8441 ) , 
          RuleConf ( 1 , darray ( [165, -4320, 54648, -445280, 2629935, -12022560, 44416680, -137057184, 364058145, -862952640, 1941643440, -4942365120, 0, 4942365120, -1941643440, 862952640, -364058145, 137057184, -44416680, 12022560, -2629935, 445280, -54648, 4320, -165] ) , 5354228880 , 24 , darray ( range ( -12 , 13 ) ) , 0.4112796986 , -67603900.0 , 0.8530 ) , 
          RuleConf ( 1 , darray ( [-198, 5577, -76050, 669240, -4275700, 21164715, -84658860, 282196200, -804259170, 2010647925, -4557468630, 9943567920, -24858919800, 0, 24858919800, -9943567920, 4557468630, -2010647925, 804259170, -282196200, 84658860, -21164715, 4275700, -669240, 76050, -5577, 198] ) , 26771144400 , 26 , darray ( range ( -13 , 14 ) ) , 0.4155017341 , 280816200.0 , 0.8607 ) , 
          RuleConf ( 1 , darray ( [143, -4312, 63063, -596232, 4099095, -21861840, 94279185, -338635440, 1037071035, -2765522760, 6568116555, -14330436120, 30452176755, -74959204320, 0, 74959204320, -30452176755, 14330436120, -6568116555, 2765522760, -1037071035, 338635440, -94279185, 21861840, -4099095, 596232, -63063, 4312, -143] ) , 80313433200 , 28 , darray ( range ( -14 , 15 ) ) , 0.4192803455 , -1163381400.0 , 0.8676 ) , 
          RuleConf ( 1 , darray ( [-1001, 32175, -502425, 5080075, -37407825, 213972759, -990614625, 3820942125, -12554524125, 35803642875, -90225180045, 205057227375, -432898591125, 899097073875, -2183521465125, 0, 2183521465125, -899097073875, 432898591125, -205057227375, 90225180045, -35803642875, 12554524125, -3820942125, 990614625, -213972759, 37407825, -5080075, 502425, -32175, 1001] ) , 2329089562800 , 30 , darray ( range ( -15 , 16 ) ) , 0.4226881038 , 4808643120.0 , 0.8738 ) , 
          RuleConf ( 1 , darray ( [15015, -512512, 8511360, -91660800, 719919200, -4398051840, 21770356608, -89845916160, 315864549000, -962634816000, 2583070089600, -6199368215040, 13561117970400, -27817677888000, 56628844272000, -135909226252800, 0, 135909226252800, -56628844272000, 27817677888000, -13561117970400, 6199368215040, -2583070089600, 962634816000, -315864549000, 89845916160, -21770356608, 4398051840, -719919200, 91660800, -8511360, 512512, -15015] ) , 144403552893600 , 32 , darray ( range ( -16 , 17 ) ) , 0.425781971 , -19835652870.0 , 0.8794 ) , 
          RuleConf ( 1 , darray ( [-3640, 131495, -2314312, 26449280, -220749760, 1434873440, -7565696320, 33289063808, -124833989280, 405710465160, -1159172757600, 2950621564800, -6786429599040, 14355908767200, -28711817534400, 57423635068800, -136381133288400, 0, 136381133288400, -57423635068800, 28711817534400, -14355908767200, 6786429599040, -2950621564800, 1159172757600, -405710465160, 124833989280, -33289063808, 7565696320, -1434873440, 220749760, -26449280, 2314312, -131495, 3640] ) , 144403552893600 , 34 , darray ( range ( -17 , 18 ) ) , 0.4286073668 , 81676217699.8 , 0.8844 ) , 
          RuleConf ( 1 , darray ( [884, -33696, 626535, -7574112, 66949740, -461438208, 2582772192, -12075298560, 48150253008, -166445319040, 505577656584, -1365716267136, 3319449260400, -7353856823040, 15101670261600, -29532155178240, 58141430507160, -136803365899200, 0, 136803365899200, -58141430507160, 29532155178240, -15101670261600, 7353856823040, -3319449260400, 1365716267136, -505577656584, 166445319040, -48150253008, 12075298560, -2582772192, 461438208, -66949740, 7574112, -626535, 33696, -884] ) , 144403552893600 , 36 , darray ( range ( -18 , 19 ) ) , 0.4312010116 , -335780006098.8 , 0.8890 ) , 
          RuleConf ( 1 , darray ( [-7956, 319124, -6251076, 79701219, -743878044, 5419682892, -32101198668, 158977364832, -672040678608, 2464149154896, -7940036165776, 22737376292904, -58467539038896, 136424257757424, -292337695194480, 584675390388960, -1120627831578840, 2175336378947160, -5075784884210040, 0, 5075784884210040, -2175336378947160, 1120627831578840, -584675390388960, 292337695194480, -136424257757424, 58467539038896, -22737376292904, 7940036165776, -2464149154896, 672040678608, -158977364832, 32101198668, -5419682892, 743878044, -79701219, 6251076, -319124, 7956] ) , 5342931457063200 , 38 , darray ( range ( -19 , 20 ) ) , 0.4335929578 , 1378465288198.6 , 0.8933 ) , 
          RuleConf ( 1 , darray ( [1938, -81600, 1679600, -22526400, 221392275, -1700292672, 10626829200, -55586491200, 248402132550, -963499180800, 3285532206528, -9956158201600, 27068305110600, -66629674118400, 149916766766400, -311826874874112, 609036864988500, -1146422334096000, 2197309473684000, -5088506149584000, 0, 5088506149584000, -2197309473684000, 1146422334096000, -609036864988500, 311826874874112, -149916766766400, 66629674118400, -27068305110600, 9956158201600, -3285532206528, 963499180800, -248402132550, 55586491200, -10626829200, 1700292672, -221392275, 22526400, -1679600, 81600, -1938] ) , 5342931457063200 , 40 , darray ( range ( -20 , 21 ) ) , 0.4358080697 , -5651707681799.0 , 0.8971 ) , 
          ) , 
        # =====================================================================
        ## 2nd derivative
        # =====================================================================
        ( RuleConf ( 2 , darray ( [-1, 16, -30, 16, -1] ) , 12 , 4 , darray ( range ( -2 , 3 ) ) , 0.9045921938 , -90.0 , 0.7071 ) , 
          RuleConf ( 2 , darray ( [2, -27, 270, -490, 270, -27, 2] ) , 180 , 6 , darray ( range ( -3 , 4 ) ) , 0.9981541606 , 560.0 , 0.7486 ) , 
          RuleConf ( 2 , darray ( [-9, 128, -1008, 8064, -14350, 8064, -1008, 128, -9] ) , 5040 , 8 , darray ( range ( -4 , 5 ) ) , 1.053089532 , -3150.0 , 0.7784 ) , 
          RuleConf ( 2 , darray ( [8, -125, 1000, -6000, 42000, -73766, 42000, -6000, 1000, -125, 8] ) , 25200 , 10 , darray ( range ( -5 , 6 ) ) , 1.089371493 , 16632.0 , 0.8009 ) , 
          RuleConf ( 2 , darray ( [-50, 864, -7425, 44000, -222750, 1425600, -2480478, 1425600, -222750, 44000, -7425, 864, -50] ) , 831600 , 12 , darray ( range ( -6 , 7 ) ) , 1.115189436 , -84084.0 , 0.8187 ) , 
          RuleConf ( 2 , darray ( [900, -17150, 160524, -1003275, 4904900, -22072050, 132432300, -228812298, 132432300, -22072050, 4904900, -1003275, 160524, -17150, 900] ) , 75675600 , 14 , darray ( range ( -7 , 8 ) ) , 1.134534513 , 411840.0 , 0.8332 ) , 
          RuleConf ( 2 , darray ( [-735, 15360, -156800, 1053696, -5350800, 22830080, -94174080, 538137600, -924708642, 538137600, -94174080, 22830080, -5350800, 1053696, -156800, 15360, -735] ) , 302702400 , 16 , darray ( range ( -8 , 9 ) ) , 1.149589301 , -1969110.0 , 0.8453 ) , 
          RuleConf ( 2 , darray ( [7840, -178605, 1982880, -14394240, 77728896, -340063920, 1309875840, -5052378240, 27788080320, -47541321542, 27788080320, -5052378240, 1309875840, -340063920, 77728896, -14394240, 1982880, -178605, 7840] ) , 15437822400 , 18 , darray ( range ( -9 , 10 ) ) , 1.161650227 , 9237800.0 , 0.8555 ) , 
          RuleConf ( 2 , darray ( [-31752, 784000, -9426375, 73872000, -427329000, 1969132032, -7691922000, 27349056000, -99994986000, 533306592000, -909151481810, 533306592000, -99994986000, 27349056000, -7691922000, 1969132032, -427329000, 73872000, -9426375, 784000, -31752] ) , 293318625600 , 20 , darray ( range ( -10 , 11 ) ) , 1.171536941 , -42678636.0 , 0.8643 ) , 
          RuleConf ( 2 , darray ( [75600, -2012472, 26087600, -220114125, 1365606000, -6691469400, 27301195152, -97504268400, 325014228000, -1137549798000, 5915258949600, -10053996959110, 5915258949600, -1137549798000, 325014228000, -97504268400, 27301195152, -6691469400, 1365606000, -220114125, 26087600, -2012472, 75600] ) , 3226504881600 , 22 , darray ( range ( -11 , 12 ) ) , 1.179793595 , 194699232.0 , 0.8720 ) , 
          RuleConf ( 2 , darray ( [-381150, 10886400, -151484256, 1371462400, -9112724775, 47609337600, -205205061600, 759845028096, -2522922944850, 7973682393600, -26911178078400, 137002361126400, -232272619118930, 137002361126400, -26911178078400, 7973682393600, -2522922944850, 759845028096, -205205061600, 47609337600, -9112724775, 1371462400, -151484256, 10886400, -381150] ) , 74209612276800 , 24 , darray ( range ( -12 , 13 ) ) , 1.186795929 , -878850700.0 , 0.8787 ) , 
          RuleConf ( 2 , darray ( [1097712, -33495462, 498279600, -4823346528, 34239805600, -190672917435, 871647622560, -3389740754400, 11592913380048, -36227854312650, 109488626367120, -358326413565120, 1791632067825600, -3030960911973290, 1791632067825600, -358326413565120, 109488626367120, -36227854312650, 11592913380048, -3389740754400, 871647622560, -190672917435, 34239805600, -4823346528, 498279600, -33495462, 1097712] ) , 964724959598400 , 26 , darray ( range ( -13 , 14 ) ) , 1.192811921 , 3931426800.0 , 0.8847 ) , 
          RuleConf ( 2 , darray ( [-245388, 7968576, -126252126, 1302170688, -9847665828, 58356538240, -283120392555, 1162196830080, -4152432424140, 13287783757248, -39448108029330, 114758132448960, -365791547181060, 1800819924583680, -3040805044214090, 1800819924583680, -365791547181060, 114758132448960, -39448108029330, 13287783757248, -4152432424140, 1162196830080, -283120392555, 58356538240, -9847665828, 1302170688, -126252126, 7968576, -245388] ) , 964724959598400 , 28 , darray ( range ( -14 , 15 ) ) , 1.198037945 , -17450721000.0 , 0.8900 ) , 
          ) ,
        # =====================================================================
        ## 3rd derivative
        # =====================================================================        
        ( RuleConf ( 3 , darray ( [1, -8, 13, 0, -13, 8, -1] ) , 8 , 4 , darray ( range ( -3 , 4 ) ) , 0.7806247498 , -17.1 , 0.7644 ) , 
          RuleConf ( 3 , darray ( [-7, 72, -338, 488, 0, -488, 338, -72, 7] ) , 240 , 6 , darray ( range ( -4 , 5 ) ) , 1.017242835 , 73.8 , 0.7937 ) , 
          RuleConf ( 3 , darray ( [205, -2522, 14607, -52428, 70098, 0, -70098, 52428, -14607, 2522, -205] ) , 30240 , 8 , darray ( range ( -5 , 6 ) ) , 1.198577403 , -315.7 , 0.8155 ) , 
          RuleConf ( 3 , darray ( [-479, 6840, -46296, 198760, -603315, 764208, 0, -764208, 603315, -198760, 46296, -6840, 479] ) , 302400 , 10 , darray ( range ( -6 , 7 ) ) , 1.343055916 , 1342.4 , 0.8325 ) , 
          RuleConf ( 3 , darray ( [1239, -20137, 155775, -766968, 2717891, -7345173, 8937819, 0, -8937819, 7345173, -2717891, 766968, -155775, 20137, -1239] ) , 3326400 , 12 , darray ( range ( -7 , 8 ) ) , 1.461649255 , -5675.4 , 0.8462 ) , 
          RuleConf ( 3 , darray ( [-266681, 4861024, -42325960, 235093600, -940620590, 2910104288, -7218002792, 8514769120, 0, -8514769120, 7218002792, -2910104288, 940620590, -235093600, 42325960, -4861024, 266681] ) , 3027024000 , 14 , darray ( range ( -8 , 9 ) ) , 1.561275009 , 23873.6 , 0.8575 ) , 
          RuleConf ( 3 , darray ( [63397, -1281033, 12405267, -76813928, 342868500, -1182036366, 3302404924, -7666346376, 8823005334, 0, -8823005334, 7666346376, -3302404924, 1182036366, -342868500, 76813928, -12405267, 1281033, -63397] ) , 3027024000 , 16 , darray ( range ( -9 , 10 ) ) , 1.646520741 , -99991.3 , 0.8670 ) , 
          RuleConf ( 3 , darray ( [-514639, 11419000, -121780250, 832461000, -4107729125, 15647010528, -48168199500, 124250212000, -273621591750, 308626058000, 0, -308626058000, 273621591750, -124250212000, 48168199500, -15647010528, 4107729125, -832461000, 121780250, -11419000, 514639] ) , 102918816000 , 18 , darray ( range ( -10 , 11 ) ) , 1.720559157 , 417253.5 , 0.8752 ) , 
          RuleConf ( 3 , darray ( [49208225, -1189505461, 13856535525, -103703531750, 561216226375, -2345810864775, 7912054151547, -22270808882100, 53867283890250, -113625406977250, 126034551856850, 0, -126034551856850, 113625406977250, -53867283890250, 22270808882100, -7912054151547, 2345810864775, -561216226375, 103703531750, -13856535525, 1189505461, -49208225] ) , 41064607584000 , 20 , darray ( range ( -11 , 12 ) ) , 1.785662061 , -1735622.7 , 0.8823 ) , 
          RuleConf ( 3 , darray ( [-260258575, 6808269600, -86028592392, 699916298400, -4125149443800, 18799608088800, -69122720605400, 211597080434784, -553762006877475, 1270152527547200, -2584831235461200, 2826897047553600, 0, -2826897047553600, 2584831235461200, -1270152527547200, 553762006877475, -211597080434784, 69122720605400, -18799608088800, 4125149443800, -699916298400, 86028592392, -6808269600, 260258575] ) , 903421366848000 , 22 , darray ( range ( -12 , 13 ) ) , 1.843504073 , 7199670.8 , 0.8885 ) , 
          ) , 
        # =====================================================================
        ## 4th derivative
        # =====================================================================        
        ( RuleConf ( 4 , darray ( [7, -96, 676, -1952, 2730, -1952, 676, -96, 7] ) , 240 , 6 , darray ( range ( -4 , 5 ) ) , 4.812152765 , 184.4 , 0.8241 ) , 
          RuleConf ( 4 , darray ( [-82, 1261, -9738, 52428, -140196, 192654, -140196, 52428, -9738, 1261, -82] ) , 15120 , 8 , darray ( range ( -5 , 6 ) ) , 5.471048086 , -947.0 , 0.8409 ) , 
          RuleConf ( 4 , darray ( [479, -8208, 69444, -397520, 1809945, -4585248, 6222216, -4585248, 1809945, -397520, 69444, -8208, 479] ) , 453600 , 10 , darray ( range ( -6 , 7 ) ) , 5.957910147 , 4698.3 , 0.8543 ) , 
          RuleConf ( 4 , darray ( [-1062, 20137, -186930, 1150452, -5435782, 22035519, -53626914, 72089160, -53626914, 22035519, -5435782, 1150452, -186930, 20137, -1062] ) , 4989600 , 12 , darray ( range ( -7 , 8 ) ) , 6.332845365 , -22701.5 , 0.8652 ) , 
          RuleConf ( 4 , darray ( [800043, -16666368, 169303840, -1128449280, 5643723540, -23280834304, 86616033504, -204354458880, 272701095810, -204354458880, 86616033504, -23280834304, 5643723540, -1128449280, 169303840, -16666368, 800043] ) , 18162144000 , 14 , darray ( range ( -8 , 9 ) ) , 6.630858144 , 107431.0 , 0.8744 ) , 
          RuleConf ( 4 , darray ( [-507176, 11529297, -127597032, 921767136, -4937306400, 21276654588, -79257718176, 275988469536, -635256384048, 842762184550, -635256384048, 275988469536, -79257718176, 21276654588, -4937306400, 921767136, -127597032, 11529297, -507176] ) , 54486432000 , 16 , darray ( range ( -9 , 10 ) ) , 6.873689694 , -499956.4 , 0.8822 ) , 
          RuleConf ( 4 , darray ( [9263502, -228380000, 2740055625, -21406140000, 123231873750, -563292379008, 2167568977500, -7455012720000, 24625943257500, -55552690440000, 73346273262262, -55552690440000, 24625943257500, -7455012720000, 2167568977500, -563292379008, 123231873750, -21406140000, 2740055625, -228380000, 9263502] ) , 4631346720000 , 18 , darray ( range ( -10 , 11 ) ) , 7.07555705 , 2294894.2 , 0.8890 ) , 
          RuleConf ( 4 , darray ( [-268408500, 7137032766, -92376903500, 777776488125, -4810424797500, 23458108647750, -94944649818564, 334062133231500, -1077345677805000, 3408762209317500, -7562073111411000, 9944398288852846, -7562073111411000, 3408762209317500, -1077345677805000, 334062133231500, -94944649818564, 23458108647750, -4810424797500, 777776488125, -92376903500, 7137032766, -268408500] ) , 615969113760000 , 20 , darray ( range ( -11 , 12 ) ) , 7.246153472 , -10413736.0 , 0.8949 ) , 
          ) , 
        # =====================================================================
        ## 5th derivative
        # =====================================================================              
        ( RuleConf ( 5 , darray ( [-13, 152, -783, 1872, -1938, 0, 1938, -1872, 783, -152, 13] ) , 288 , 6 , darray ( range ( -5 , 6 ) ) , 3.98337948 , 43.5 , 0.8463 ) , 
          RuleConf ( 5 , darray ( [139, -1936, 12500, -48176, 101559, -99744, 0, 99744, -101559, 48176, -12500, 1936, -139] ) , 12096 , 8 , darray ( range ( -6 , 7 ) ) , 5.089983867 , -175.1 , 0.8597 ) , 
          RuleConf ( 5 , darray ( [-518, 8301, -62710, 295244, -944862, 1819681, -1718382, 0, 1718382, -1819681, 944862, -295244, 62710, -8301, 518] ) , 181440 , 10 , darray ( range ( -7 , 8 ) ) , 6.057602156 , 712.6 , 0.8706 ) , 
          RuleConf ( 5 , darray ( [4201, -75908, 652023, -3539780, 13565962, -38061684, 68459875, -62714036, 0, 62714036, -68459875, 38061684, -13565962, 3539780, -652023, 75908, -4201] ) , 5987520 , 12 , darray ( range ( -8 , 9 ) ) , 6.908238146 , -2914.3 , 0.8796 ) , 
          RuleConf ( 5 , darray ( [-3739217, 75119112, -721271943, 4407497768, -19241468100, 63619040016, -161682804556, 275637687624, -246459164094, 0, 246459164094, -275637687624, 161682804556, -63619040016, 19241468100, -4407497768, 721271943, -75119112, 3739217] ) , 21794572800 , 14 , darray ( range ( -9 , 10 ) ) , 7.661614351 , 11944.9 , 0.8872 ) , 
          RuleConf ( 5 , darray ( [1824595, -40321144, 427576664, -2898570696, 14119093201, -52627196640, 155526600912, -365798390432, 597244221678, -523564225808, 0, 523564225808, -597244221678, 365798390432, -155526600912, 52627196640, -14119093201, 2898570696, -427576664, 40321144, -1824595] ) , 43589145600 , 16 , darray ( range ( -10 , 11 ) ) , 8.333977317 , -48997.9 , 0.8938 ) , 
          RuleConf ( 5 , darray ( [-113425697, 2733785665, -31719348453, 236068829960, -1267132147015, 5229615477963, -17266767656955, 46693491257712, -103170444595530, 162555496564570, -140176720604882, 0, 140176720604882, -162555496564570, 103170444595530, -46693491257712, 17266767656955, -5229615477963, 1267132147015, -236068829960, 31719348453, -2733785665, 113425697] ) , 11115232128000 , 18 , darray ( range ( -11 , 12 ) ) , 8.938399421 , 200998.8 , 0.8996 ) , 
          ) ,          
        # =====================================================================
        ## 6th derivative
        # =====================================================================        
        ( RuleConf ( 6 , darray ( [-695, 11616, -93750, 481760, -1523385, 2992320, -3735732, 2992320, -1523385, 481760, -93750, 11616, -695] ) , 60480 , 8 , darray ( range ( -6 , 7 ) ) , 29.02837162 , -408.6 , 0.8743 ) , 
          RuleConf ( 6 , darray ( [148, -2767, 25084, -147622, 629908, -1819681, 3436764, -4243668, 3436764, -1819681, 629908, -147622, 25084, -2767, 148] ) , 60480 , 10 , darray ( range ( -7 , 8 ) ) , 33.44293272 , 1900.3 , 0.8833 ) , 
          RuleConf ( 6 , darray ( [-4201, 86752, -869364, 5663648, -27131924, 101497824, -273839500, 501712288, -614231046, 501712288, -273839500, 101497824, -27131924, 5663648, -869364, 86752, -4201] ) , 7983360 , 12 , darray ( range ( -8 , 9 ) ) , 37.10243085 , -8743.0 , 0.8909 ) , 
          RuleConf ( 6 , darray ( [3739217, -84509001, 927349641, -6611246652, 34634642580, -143142840036, 485048413668, -1240369594308, 2218132476846, -2697076863910, 2218132476846, -1240369594308, 485048413668, -143142840036, 34634642580, -6611246652, 927349641, -84509001, 3739217] ) , 32691859200 , 14 , darray ( range ( -9 , 10 ) ) , 40.18024111 , 39816.3 , 0.8974 ) , 
          RuleConf ( 6 , darray ( [-3284271, 80642288, -962047494, 7453467504, -42357279603, 189457907904, -699869704104, 2194790342592, -5375197995102, 9424156064544, -11395096228516, 9424156064544, -5375197995102, 2194790342592, -699869704104, 189457907904, -42357279603, 7453467504, -962047494, 80642288, -3284271] ) , 130767436800 , 16 , darray ( range ( -10 , 11 ) ) , 42.80320126 , -179658.8 , 0.9030 ) , 
          RuleConf ( 6 , darray ( [61868562, -1640271399, 21146232302, -177051622470, 1086113268870, -5229615477963, 20720121188346, -70040236886568, 206340889191060, -487666489693710, 841060323629292, -1012227242852644, 841060323629292, -487666489693710, 206340889191060, -70040236886568, 20720121188346, -5229615477963, 1086113268870, -177051622470, 21146232302, -1640271399, 61868562] ) , 11115232128000 , 18 , darray ( range ( -11 , 12 ) ) , 45.06485914 , 803995.1 , 0.9080 ) , 
          ) , 
        # =====================================================================
        )
    # =========================================================================
        
    # =========================================================================
    ## get the overal configruiation for all known forward rules  
    @classprop
    def RULES ( cls ) :
        """Get the overal configuration for forward rules
        """
        return cls.__CENTRAL_RULES

    # =========================================================================
    ## get the maximal derivative order (inclusive) 
    @classprop
    def DMAX ( cls ) :
        """get the maximal derivative order (inclusive)"""
        return len ( cls.RULES ) 
        
    # =========================================================================
    ## get maximal rule index for the given derivative order 
    @classmethod
    def IMAX ( cls , D ) :
        """get maximal rule index for the given derivative order"""
        
        assert isinstance ( D , int ) and 1 <= D <= cls.DMAX , \
               'Invalid value of parameter ``D'' '        
        return len ( cls.RULES [ D - 1 ] ) 

    # =========================================================================
    def __init__ ( self , D , I = 3 , with_error = False , max_step = -1 ) :
        
        assert isinstance ( D , int ) and 1 <= D <= self.DMAX , \
               'Invalid value of parameter ``D'' '
        
        imax = len ( self.RULES [ D - 1 ]  ) 
        assert isinstance ( D , int ) and 0 <= I <  imax , \
               'Invalid value of parameter ``I'' must be 0<%s<%d ' % ( I , imax )
        
        ## the row from the difference table
        config = self.RULES [ D - 1 ] [ I ]

        assert config.order == D , 'Invalid cofnfigurtaion!'
        
        assert 1 == len ( config.stencil  ) % 2 , \
               'Stencil length must be odd!'
        
        if not with_error : rule_err = None
        elif 0 < I        : rule_err = CentralRule ( D = D , I = I - 1 , with_error = False , max_step = max_step )
        else              : rule_err = CentralRule ( D = D , I = 0     , with_error = False , max_step = max_step )
        
        # =====================================================================
        ## initialize the base class 
        DiffRule.__init__ ( self                ,
                            config              ,
                            rule_err = rule_err ,
                            max_step = max_step )
                            
        ## index 
        self.__I = I
        
    @property
    def I  ( self ) :
        """``I'': index/row of the rule in the global table"""
        return self.__I

    # ====================================================================================
    ## The main method: evaluate the numerical derivative 
    def __call__ ( self , func , x , h , args = () , kwargs = {}) :
        """The main method: evaluate the numerical derivative 
        - func   the function with signature `func(x,*args,**kwargs)`
        - x      point where derivative is evaluated
        - h      step size 
        - args   addtional positional arguments
        - kwargs addtional keyword arguments
        """

        fv   = self.funvals 
        l    = len ( fv )
        m    = l // 2 

        d    = self.D

        for j in range ( 1 , m + 1 ) :
            
            fv [ m - j ] = func ( x - j * h , *args , **kwargs ) ## 
            fv [ m + j ] = func ( x + j * h , *args , **kwargs ) ##

        if ( 0 == self.D % 2 ) :
            fv [ m     ] = func ( x , *args , **kwargs ) ##
        else : 
            fv [ m     ] = 0  ## not used 

            
        ## calculate derivatives
        coeff   = self.coefficients 
        scale   = self.scale         
        result  = the_dot ( l , fv , coeff )
        result /= scale * ( h ** d )
                
        if ( not self.with_error ) or ( not self.rule_err ) : return result       ## RETURN
        
        ## make an estimate for the roundoff error
        fmax  = max ( fv , key = abs )
        eroff = self.roundoff_error ( h , fmax ) 

        ## for this case use Richardson's extrapolation to get an estimate and the error 
        if 0 == self.I :
            ## make an error estimate using smaller steps 
            t     = self.t_optimal  
            with WithError ( self.rule_err , False ) : 
                result2 = self.rule_err ( func , x , t * h , args = args , kwargs = kwargs )
            tn    = t ** self.N 
            r      = ( result2 - tn * result ) / ( 1.0 - tn ) 
            etrunc =   r - result            
            return VE ( result , etrunc * etrunc + eroff * eroff ) 
            
        ## make an estimate for the truncation error from comparison with the lower order rule 
        
        coeff2   = self.rule_err.coefficients 
        scale2   = self.rule_err.scale 
        result2  = the_dot ( l - 2 , fv , coeff2  , 1 , 0  )
        result2 /= ( scale2 * ( h ** d ) )

        etrunc = result2 - result

        return VE ( result , etrunc * etrunc + eroff * eroff ) 



# =============================================================================
## @class ForwardOpenD1
#  Forward open rule for the 1st derivative 
class ForwardOpenD1 (ForwardOpen) :
    """Forward open rule for the 1st derivative"""
    def __init__ ( self , I = 3 , with_error = False , maX_step = -1 ) :
        ForwardOpen.__init__ ( self , 1 , I , with_error = with_error , max_step = max_step )
# =============================================================================
## @class ForwardOpenD2
#  Forward open rule for the 2nd derivative     
class ForwardOpenD2 (ForwardOpen) :
    """Forward open rule for the 2nd derivative"""
    def __init__ ( self , I = 3 , with_error = False , max_step = -1 ) :
        ForwardOpen.__init__ ( self , 2 , I , with_error = with_error , max_step = max_step )
# =============================================================================
## @class ForwardOpenD3
#  Forward open rule for the 3rd derivative     
class ForwardOpenD3 (ForwardOpen) :
    """Forward open rule for the 3rd derivative"""
    def __init__ ( self , I = 3 , with_error = False , max_step = -1 ) :
        ForwardOpen.__init__ ( self , 3 , I , with_error = with_error , max_step = max_step )
# =============================================================================
## @class ForwardOpenD4
#  Forward open rule for the 4th derivative     
class ForwardOpenD4 (ForwardOpen) :
    """Forward open rule for the 4th derivative"""
    def __init__ ( self , I = 3 , with_error = False , max_step = -1 ) :
        ForwardOpen.__init__ ( self , 4 , I , with_error = with_error , max_step = max_step )
# =============================================================================
## @class ForwardOpenD5
#  Forward open rule for the 5th derivative     
class ForwardOpenD5 (ForwardOpen) :
    """Forward open rule for the 5th derivative"""
    def __init__ ( self , I = 3 , with_error = False , max_step = -1 ) : 
        ForwardOpen.__init__ ( self , 5 , I , with_error = with_error , max_step = max_step )
# =============================================================================
## @class ForwardOpenD6
#  Forward open rule for the 6th derivative     
class ForwardOpenD6 (ForwardOpen) :
    """Forward open rule for the 6th derivative"""
    def __init__ ( self , I = 3 , with_error = False , max_step = -1 ) :
        ForwardOpen.__init__ ( self , 6 , I , with_error = with_error , max_step = max_step )
# =============================================================================


# =============================================================================
## @class BackwardOpenD1
#  Backward open rule for the 1st derivative 
class BackwardOpenD1 (BackwardOpen) :
    """Backward open rule for the 1st derivative"""
    def __init__ ( self , I = 3 , with_error = False , max_step = -1 ) :
        BackwardOpen.__init__ ( self , 1 , I , with_error = with_error , max_step = max_step )
# =============================================================================
## @class BackwardOpenD1
#  Backward open rule for the 1st derivative 
class BackwardOpenD2 (BackwardOpen) :
    """Backward open rule for the 2nd derivative"""
    def __init__ ( self , I = 3 , with_error = False , max_step = -1 ) :
        BackwardOpen.__init__ ( self , 2 , I , with_error = with_error , max_step = max_step )
# =============================================================================
## @class BackwardOpenD1
#  Backward open rule for the 1st derivative 
class BackwardOpenD3 (BackwardOpen) :
    """Backward open rule for the 3rd derivative"""
    def __init__ ( self , I = 3 , with_error = False , max_step = -1 ) :
        BackwardOpen.__init__ ( self , 3 , I , with_error = with_error , max_step = max_step )
# =============================================================================
## @class BackwardOpenD1
#  Backward open rule for the 1st derivative 
class BackwardOpenD4 (BackwardOpen) :
    """Backward open rule for the 4th derivative"""
    def __init__ ( self , I = 3 , with_error = False , maX_step = -1 ) :
        BackwardOpen.__init__ ( self , 4 , I , with_error = with_error , max_step = max_step )
# =============================================================================
## @class BackwardOpenD1
#  Backward open rule for the 1st derivative 
class BackwardOpenD5 (BackwardOpen) :
    """Backward open rule for the 5th derivative"""
    def __init__ ( self , I = 3 , with_error = False , max_step = -1 ) :
        BackwardOpen.__init__ ( self , 5 , I , with_error = with_error , max_step = max_step )
# =============================================================================
## @class BackwardOpenD6
#  Backward open rule for the 6th derivative 
class BackwardOpenD6 (BackwardOpen) :
    """Backward open rule for the 1st derivative"""
    def __init__ ( self , I = 3 , with_error = False , max_step = -1 ) :
        BackwardOpen.__init__ ( self , 6 , I , with_error = with_error , max_step = max_step )
# =============================================================================

# =============================================================================
## @class CentralRuleD1
#  Central rule for the 1st derivative 
class CentralRuleD1 (CentralRule) :
    """Central rule for the 1st derivative"""
    def __init__ ( self , I = 2 , with_error = False , max_step = -1 ) :
        CentralRule.__init__ ( self , 1 , I , with_error = with_error , max_step = max_step )
# =============================================================================
## @class CentralRuleD2
#  Central rule for the 2nd derivative 
class CentralRuleD2 (CentralRule) :
    """Central rule for the 2nd derivative"""
    def __init__ ( self , I = 2 , with_error = False , max_step = -1 ) :
        CentralRule.__init__ ( self , 2 , I , with_error = with_error , max_step = max_step )
# =============================================================================
## @class CentralRuleD3
#  Central rule for the 3rd derivative 
class CentralRuleD3 (CentralRule) :
    """Central rule for the 3rd derivative"""
    def __init__ ( self , I = 2 , with_error = False , maX_step = -1 ) :
        CentralRule.__init__ ( self , 3 , I , with_error = with_error , max_step = max_step )
# =============================================================================
## @class CentralRuleD4
#  Central rule for the 4th derivative 
class CentralRuleD4 (CentralRule) :
    """Central rule for the 4th derivative"""
    def __init__ ( self , I = 2 , with_error = False , max_step = -1 ) :
        CentralRule.__init__ ( self , 4 , I , with_error = with_error , max_step = max_step )
# =============================================================================
## @class CentralRuleD5
#  Central rule for the 5th derivative 
class CentralRuleD5 (CentralRule) :
    """Central rule for the 5th derivative"""
    def __init__ ( self , I = 2 , with_error = False , max_step = -1 ) :
        CentralRule.__init__ ( self , 5 , I , with_error = with_error , max_step = max_step )
# =============================================================================
## @class CentralRuleD6
#  Central rule for the 6th derivative 
class CentralRuleD6 (CentralRule) :
    """Central rule for the 6th derivative"""
    def __init__ ( self , I = 2 , with_error = False , max_step = -1 ) :
        CentralRule.__init__ ( self , 6 , I , with_error = with_error , max_step = max_step )

# =============================================================================
## Get an optimal value for paramters for Richardson's extrapolation
#  It is a (numerical)solution of equation \f$ f(t) = 0 \f$, where
#  \f$ f(t) = n t^{n+d} + (n+d) t^n - d \f$ 
@memoize 
def t_opt ( d , n ) :
    """Get an optimal value for paramters for Richardson's extrapolation
    It is a (numerical) solution of equation `f(t) = 0`, where
    `f(t) = n t^{n+d} + (n+d) t^n - d` 
    """
    fun = lambda t :  n*pow ( t , n + d ) + ( n + d ) * pow ( t , n ) - d
    from ostap.math.rootfinder import solve
    return solve ( fun , 1.e-15 , 1-1.e-15 )

# =============================================================================
## Richardson's extrapolation for numerical differentiation rule
#  - it is a Rule itself
#  - it can be chained, iterated, combined with other rules
#  @see https://en.wikipedia.org/wiki/Richardson_extrapolation 
class Richardson(Rule) :
    """Richarsson's extrapolation for numerical differentiation rule
    - it is a Rule itself
    - it can be chained, iterated, combined with other rules 
    - see https://en.wikipedia.org/wiki/Richardson_extrapolation 
    """
    
    def __init__ ( self , rule , t = 0.5 , level = 2 , with_error = False , max_step = -1 ) :

        assert isinstance ( rule , Rule ) , \
               "``rule'' must be Rule! got %s/%s" % ( rule , type ( rule ) )

        assert isinstance ( level , int ) and 1 <= level, \
               "Invalid ``level'' %s/%s" % ( level , type ( level ) ) 

        assert 0 < t < 1 , 'Invalid value of ``t'' parameter: it must be 0<%s<1' % t 
        
        Rule.__init__ ( self                        ,
                        D          = rule.D         ,
                        N          = rule.N + level , ## ATTENTION! 
                        stencil    = rule.stencil   , 
                        interval   = rule.interval  ,
                        with_error = with_error     ,
                        max_step   = max_step       ) 
        
        self.__rule  = rule 
        self.__t     = float ( t )
        self.__level = int ( level )

    @property
    def rule ( self ) :
        """``rule'': underlying rule for Richardson's extrapolation"""
        return self.__rule

    @property
    def t    ( self )  :
        """``t'' : t-parameter for Richardson's  extrapolation: 0<t<1"""
        return self.__t
    
    @property
    def level ( self ) :
        """``level'' : number of levels for (iterative) Richardson's extrapolation)"""
        return self.__level 
        
    # =========================================================================
    ## The main method: evaluate the numerical derivative with the given step 
    #  @param func   the function <code>func(x,*args,**kwargs)</code>
    #  @param x      point where derivative is evaluated
    #  @param h      the step 
    #  @param args   additional positional arguments
    #  @param kwargs additional kewword arguments
    def __call__ ( self , func , x , h , args = () , kwargs = {} ) :
        """The main method: evalaute th enumerical derivarive 
        - func     the function with signature `func(x,*args,**kwargs)`
        - x        point where derivative is evaluated
        - h        the step 
        - args     additional positional arguments
        - kwargs   additional keyword arguments
        """

        h0     = h * 1.0

        t      = self.t
        rule   = self.__rule
        
        with WithError ( rule , False ) : 
            values = [
                float ( rule ( func , x , h0 * ( t ** i ) , args = args , kwargs = kwargs ) ) for i in range ( self.level + 1 ) 
                ]
        
        tn   = t ** self.__rule.N
        dt   = 1.0 - tn

        ## iterate 
        while 2 < len ( values ) :
            values = [  ( values [ i + 1 ] - tn * values [ i ] ) / dt for i in range ( len ( values ) - 1 ) ]

            tn    *= t        ## each iteration improves ``effective n'' by one unit 
            dt     = 1.0 - tn            

            
        ## after iterations 
        v1 = values [ 0 ]
        v2 = values [ 1 ]

        ## final result (the last iteration) 
        v  = float ( v2 - tn * v1 ) / ( 1.0 - tn ) ## the last step 
        
        if not self.with_error : return v          ## RETURN 

        ## get an error estimate from the best estimate from the previous step 
        e = float ( v - v2 )
        
        return VE ( v , e * e ) 
    
# =============================================================================
## @class DerivativeFD
#  Very flexible class to calculate derivatives using combniation of
#  various rule. If <code>singular</code> points are specified, in their vinicity
#  special <code>forward</code> and <code>backward</code> are used
#  Optionally one can aqpply Richardson's iterative extrapolation
#  to improve results 
#  @code
#  central    = CentralRule  ( 1 , 3 )
#  forward    = ForwardOpen  ( 1 , 3 )
#  backward   = BackwardOpen ( 1 , 3 )
#  derivative = DerivativeFD  ( fun , central    = central         , 
#                                     forward    = forward         ,
#                                     backward   = backward        ,
#                                     singular   = [ -1 , 0 , -1 ] ,
#                                     with_error = True            ,
#                                     Richardson = 2               )
#  x     = 0.3
#  value = derivative ( x ) 
#  @endcode 
class DerivativeFD(object) :
    """Very flexible class to calculate derivatives using combniation of
    various rule. If <code>singular</code> points are specified, in their vinicity
    special <code>forward</code> and <code>backward</code> are used

    Optionally one can aqpply Richardson's iterative extrapolation
    to improve results
    
    >>> central    = CentralRule  ( 1 , 3 )
    >>> forward    = ForwardOpen  ( 1 , 3 )
    >>> backward   = BackwardOpen ( 1 , 3 )
    >>> derivative = DerivativeFD  ( fun , central    = central         , 
    ...                                    forward    = forward         ,
    ...                                    backward   = backward        ,
    ...                                    singular   = [ -1 , 0 , -1 ] ,
    ...                                    with_error = True            ,
    ...                                    Richardson = 2               )
    >>> x     = 0.3
    >>> value = derivative ( x ) 
    """
    def __init__ ( self               ,
                   fun                ,   ## the function
                   central            ,   ## central rule  
                   forward    = None  ,   ## forward rule 
                   backward   = None  ,   ## backward rule 
                   singular   = ()    ,   ## singular points
                   with_error = False ,   ## estimate the error?  
                   max_step   = -1    ,   ## max_step (if needed)
                   step       =  0    ,   ## proposed initial steo
                   richardson =  0    ) : ## use Richardson's extrapolation?
        
        assert callable  ( fun ) , "Invalid type for ``fun''!"
        assert isinstance ( central , Rule ) , \
               "Invalid type for the ``central'' rule %s" % central
        assert isinstance ( richardson , int ) and 0 <= richardson , \
               "Invalid type/valeue for ~`richardson'' %s/%s" % ( richardson , type ( richardson ) )
        
        points = tuple ( sorted ( set ( singular ) ) )
        
        if points :
            
            assert forward  and isinstance ( forward  , Rule ) and 0 <= forward.interval[0]   ,\
                   "Invalid ``forward'' rule %s" % forward
            assert backward and isinstance ( backward , Rule ) and backward.interval[-1] <= 0 ,\
                   "Invalid ``backward'' rule %s" % backward 
            
            assert central.D == forward .D , 'Mismatch for central/forward derivative orders!'
            assert central.D == backward.D , 'Mismatch for central/forward derivative orders!'
            
            
        self.__points     = points
        self.__central    = central
        self.__forward    = forward  if points else None 
        self.__backward   = backward if points else None 
        self.__fun        = fun
        self.__with_error = True if with_error else False 
        self.__max_step   = float ( max_step ) 
        self.__richardson = int ( richardson ) 
        self.__step       = max ( float ( step ) , 0 )
        
        if self.richardson : 
            if not hasattr ( central  , 't_optimal' ) :
                central .t_optimal = t_opt ( central .D , central .N )
            if forward  and not hasattr ( forward  , 't_optimal' ) :
                forward .t_optimal = t_opt ( forward .D , forward .N )
            if backward and not hasattr ( backward , 't_optimal' ) :
                backward.t_optimal = t_opt ( backward.D , backward.N )
            
    @property
    def D ( self ) :
        """``D'' : derivative order"""
        return self.__central.D
    
    @property
    def central ( self ) :
        """``central'' : central differentiation rule"""
        return self.__central
    @property
    def forward ( self ) :
        """``forward'' : forward differentiation rule"""
        return self.__forward
    @property
    def central ( self ) :
        """``backward'' : backward differentiation rule"""
        return self.__backward

    @property
    def step ( self ) :
        """``step'' : proposed intial step"""
        return self.__step
    @property
    def h0   ( self ) :
        """``h0'': proposed initial step (same as ``step''"""
        return self.__step 
    
    @property
    def points  ( self ) :
        """``points'' : list of singular points"""
        return self.__points
    
    @property
    def with_error ( self ) :
        return self.__with_error
    @with_error.setter
    def with_error ( self , value ) :
        self.__with_error = True if value else False 
        
    @property
    def max_step ( self ) :
        """``max_step'' : maximal value of the step"""
        return self.__max_step
        
    @max_step.setter
    def max_step ( self , value ) :
        self.__max_step = float ( value ) if 0 < value else -1.0

    @property
    def richardson ( self ) :
        """``richardson'': apply Richardson (iterative) extrapolation to the final result?"""
        return self.__richardson

    ## find an approriate rule ot be use and corresponding hmax 
    def find_rule ( self , x ) :
        """Find an approriate rule ot be use and corresponding hmax-parameter"""
        
        x      = float ( x ) 
        points = self.__points
        hmx    = self.max_step

        ## simple case 
        if not points : return self.__central, hmx 
        
        imax   = len ( points )
        
        i_c    = self.__central.interval
        
        h_c    = self.__central.h_guess ( x )
        n_p    = abs ( max ( i_c , key = abs ) ) 
        
        region = h_c * n_p

        ip   = bisect.bisect_left ( points  , x )

        ## below the minimal point 
        if   0    == ip :

            xmax  = points [ 0 ]
            dist1 = float ( xmax - x ) 

            ## use backward or cerntral rules 
            if dist1 < 1.1 * region : return self.__backward , hmx
            else                    : return self.__central  , get_hmax ( 0.95 * dist1 / n_p , hmx )  

        ## above maximal point 
        elif imax == ip :
            
            xmin  = points [-1]
            dist2 = float ( x - xmin ) 
            
            ## use forward or central rules 
            if dist2 < 1.1 * region : return self.__forward , hmx 
            else                    : return self.__central , get_hmax ( 0.95 * dist2 / n_p , hmx ) 
                            
        ## within some interval 
        xmin  = points [ ip -1 ]
        xmax  = points [ ip    ]

        dist1 = float ( xmax - x    )
        dist2 = float ( x    - xmin )

        ## use backward rule 
        if dist1   < 1.1 * region and dist1 <= dist2 :
            rule = self.__backward 
            return rule , get_hmax ( 0.95 * dist2 / abs ( rule.interval[0] )  , hmx )                 
            ## use forward rule 
        elif dist2 < 1.1 * region and dist2 <= dist1 :
            rule = self.__forward 
            return rule, get_hmax ( 0.95 * dist1 / abs ( rule.interval[1] )   , hmx )
        
        ## use central rule 
        return self.__central , get_hmax ( 0.95 * min ( dist1 , dist2 ) / n_p , hmx ) 

    ## the main method
    def __call__  ( self , x , *args , **kwargs ) :

        x      = float ( x ) 
        ## find a proper rule to use 
        rule , hmax = self.find_rule ( x )

        ## the function itself 
        fun    = self.__fun
                    
        ## get the optimal step
        hopt , fval = rule.optimal_step ( fun , x , self.step , hmax , args = args , kwargs = kwargs )
        
        if not self.richardson :
            with WithError ( rule , self.with_error ) : 
                return rule ( fun , x , hopt , args = args , kwargs = kwargs )
        
        # =====================================================================
        ## make Richardson's (iterative) extrapolation
        # =====================================================================
        
        ## get Richardson's parameter
        t = rule.t_optimal
        
        with WithError ( rule , False ) : 
            values = [
                float ( rule ( fun , x , hopt * ( t ** i ) , args = args , kwargs = kwargs ) ) for i in range ( self.richardson + 1 ) 
                ]
            
        tn   = t ** rule.N
        dt   = 1.0 - tn

        ## iterate 
        while 2 < len ( values ) :
            values = [ ( values [ i + 1 ] - tn * values [ i ] ) / dt for i in range ( len ( values ) - 1 ) ]
            
            tn    *= t        ## each iteration improves ``effective n'' by one unit 
            dt     = 1.0 - tn            
            
        ## after iterations 
        v1 = values [ 0 ]
        v2 = values [ 1 ]

        ## final result (the last iteration) 
        v  = float ( v2 - tn * v1 ) / ( 1.0 - tn ) ## the last step 
        
        if not self.with_error : return v          ## RETURN 

        ## get an error estimate from the best estimate from the previous step 
        e = float ( v - v2 )
        
        return VE ( v , e * e ) 


# =============================================================================
## @class Derivative1
#  Very flexible class to calculate 1st derivative using combniation of
#  various rule.
#  - If <code>singular</code> points are specified, in their vinicity
#  special <code>forward</code> and <code>backward</code> are used
#  - Optionally one can apply Richardson's iterative extrapolation
#  to improve results 
#  @code
#  derivative = Derivative1 ( fun )
#  x     = 0.3
#  value = derivative ( x ) 
#  @endcode 
class Derivative1(DerivativeFD) :
    """Very flexible class to calculate 1st derivative using combniation of
    various rule.
    - If `singular` points are specified, in their vinicity
    special `forward` and `backward` rules  are used
    - Optionally one can apply Richardson's iterative extrapolation
    to improve results
    
    >>> derivative = Derivative1 ( fun )
    >>> x     = 0.3
    >>> value = derivative ( x )
    """
    @classprop 
    def IMAX ( cls ) :
        """``IMAX'' : maximal rule index"""
        return CentralRule.IMAX(1)
    
    def __init__ ( self               ,
                   fun                ,
                   step       = 0     ,
                   calc       = 3     , ## rule index 
                   with_error = False ,
                   max_step   = -1    ,
                   central    = None  ,
                   forward    = None  ,
                   backward   = None  ,
                   singular   = ()    ,
                   richardson = 0     ) :
        
        D  = 1
        IC = min ( calc , CentralRule .IMAX ( D ) - 1 ) 
        IF = min ( calc , ForwardOpen .IMAX ( D ) - 1 ) 
        IB = min ( calc , BackwardOpen.IMAX ( D ) - 1 ) 
                   
        if not central :
            central = CentralRule( D , I = IC  , max_step = max_step )
            
        assert isinstance ( central , Rule ) and D == central.D , \
               'Invalid central finite difference!'

        points = tuple ( sorted ( set ( singular ) ) )
        if points :
            if not forward  : forward  = ForwardOpen  ( D , I = IF , max_step = max_step )
            if not backward : backward = BackwardOpen ( D , I = IB , max_step = max_step )
            
        ## initialize the base clas
        DerivativeFD.__init__ ( self                    ,
                                fun        = fun        ,
                                central    = central    ,
                                forward    = forward    , 
                                backward   = backward   , 
                                singular   = points     ,
                                with_error = with_error ,
                                max_step   = max_step   ,
                                step       = step       ,
                                richardson = richardson )
                                

# =============================================================================
## @class Derivative2
#  Very flexible class to calculate 2nd derivative using combniation of
#  various rule.
#  - If <code>singular</code> points are specified, in their vinicity
#  special <code>forward</code> and <code>backward</code> are used
#  - Optionally one can aqpply Richardson's iterative extrapolation
#  to improve results 
#  @code
#  derivative = Derivative1 ( fun )
#  x     = 0.3
#  value = derivative ( x ) 
#  @endcode 
class Derivative2(DerivativeFD) :
    """Very flexible class to calculate 2nd derivative using combniation of
    various rule.
    - If `singular` points are specified, in their vinicity
    special `forward` and `backward` rules  are used
    - Optionally one can aqpply Richardson's iterative extrapolation
    to improve results
    
    >>> derivative = Derivative2 ( fun )
    >>> x     = 0.3
    >>> value = derivative ( x )
    """
    @classprop 
    def IMAX ( cls ) :
        """``IMAX'' : maximal rule index"""
        return CentralRule.IMAX(2)
    
    def __init__ ( self               ,
                   fun                ,
                   step       = 0     ,
                   calc       = 3     , ## rule index 
                   with_error = False ,
                   max_step   = -1    ,
                   central    = None  ,
                   forward    = None  ,
                   backward   = None  ,
                   singular   = ()    ,
                   richardson = 0     ) :

        D  = 2
        IC = min ( calc , CentralRule .IMAX ( D ) - 1 ) 
        IF = min ( calc , ForwardOpen .IMAX ( D ) - 1 ) 
        IB = min ( calc , BackwardOpen.IMAX ( D ) - 1 ) 
                   
        if not central :
            central = CentralRule( D , I = IC  , max_step = max_step )
            
        assert isinstance ( central , Rule ) and D == central.D , \
               'Invalid central finite difference!'

        points = tuple ( sorted ( set ( singular ) ) )
        if points :
            if not forward  : forward  = ForwardOpen  ( D , I = IF , max_step = max_step )
            if not backward : backward = BackwardOpen ( D , I = IB , max_step = max_step )
            
        ## initialize the base clas
        DerivativeFD.__init__ ( self                    ,
                                fun        = fun        ,
                                central    = central    ,
                                forward    = forward    , 
                                backward   = backward   , 
                                singular   = points     ,
                                with_error = with_error ,
                                max_step   = max_step   ,
                                step       = step       ,
                                richardson = richardson )
                                
                                
# =============================================================================
## @class Derivative3
#  Very flexible class to calculate 3rd derivative using combniation of
#  various rule.
#  - If <code>singular</code> points are specified, in their vinicity
#  special <code>forward</code> and <code>backward</code> are used
#  - Optionally one can apply Richardson's iterative extrapolation
#  to improve results 
#  @code
#  derivative = Derivative3 ( fun )
#  x     = 0.3
#  value = derivative ( x ) 
#  @endcode 
class Derivative3(DerivativeFD) :
    """Very flexible class to calculate 3rd derivative using combniation of
    various rule.
    - If `singular` points are specified, in their vinicity
    special `forward` and `backward` rules  are used
    - Optionally one can apply Richardson's iterative extrapolation
    to improve results
    
    >>> derivative = Derivative3 ( fun )
    >>> x     = 0.3
    >>> value = derivative ( x )
    """
    @classprop 
    def IMAX ( cls ) :
        """``IMAX'' : maximal rule index"""
        return CentralRule.IMAX(3)
    
    def __init__ ( self               ,
                   fun                ,
                   step       = 0     ,
                   calc       = 3     , ## rule index 
                   with_error = False ,
                   max_step   = -1    ,
                   central    = None  ,
                   forward    = None  ,
                   backward   = None  ,
                   singular   = ()    ,
                   richardson = 0     ) :

        D  = 3
        IC = min ( calc , CentralRule .IMAX ( D ) - 1 ) 
        IF = min ( calc , ForwardOpen .IMAX ( D ) - 1 ) 
        IB = min ( calc , BackwardOpen.IMAX ( D ) - 1 ) 
                   
        if not central :
            central = CentralRule( D , I = IC  , max_step = max_step )
            
        assert isinstance ( central , Rule ) and D == central.D , \
               'Invalid central finite difference!'

        points = tuple ( sorted ( set ( singular ) ) )
        if points :
            if not forward  : forward  = ForwardOpen  ( D , I = IF , max_step = max_step )
            if not backward : backward = BackwardOpen ( D , I = IB , max_step = max_step )
            
        ## initialize the base clas
        DerivativeFD.__init__ ( self                    ,
                                fun        = fun        ,
                                central    = central    ,
                                forward    = forward    , 
                                backward   = backward   , 
                                singular   = points     ,
                                with_error = with_error ,
                                max_step   = max_step   ,
                                step       = step       ,
                                richardson = richardson )
                                

# =============================================================================
## @class Derivative4
#  Very flexible class to calculate 4th derivative using combniation of
#  various rule.
#  - If <code>singular</code> points are specified, in their vinicity
#  special <code>forward</code> and <code>backward</code> are used
#  - Optionally one can apply Richardson's iterative extrapolation
#  to improve results 
#  @code
#  derivative = Derivative4 ( fun )
#  x     = 0.3
#  value = derivative ( x ) 
#  @endcode 
class Derivative4(DerivativeFD) :
    """Very flexible class to calculate 4th derivative using combniation of
    various rule.
    - If `singular` points are specified, in their vinicity
    special `forward` and `backward` rules  are used
    - Optionally one can apply Richardson's iterative extrapolation
    to improve results
    
    >>> derivative = Derivative4 ( fun )
    >>> x     = 0.3
    >>> value = derivative ( x )
    """
    @classprop 
    def IMAX ( cls ) :
        """``IMAX'' : maximal rule index"""
        return CentralRule.IMAX(4)
    
    def __init__ ( self               ,
                   fun                ,
                   step       = 0     ,
                   calc       = 2     , ## rule index 
                   with_error = False ,
                   max_step   = -1    ,
                   central    = None  ,
                   forward    = None  ,
                   backward   = None  ,
                   singular   = ()    ,
                   richardson = 0     ) :

        D  = 4
        IC = min ( calc , CentralRule .IMAX ( D ) - 1 ) 
        IF = min ( calc , ForwardOpen .IMAX ( D ) - 1 ) 
        IB = min ( calc , BackwardOpen.IMAX ( D ) - 1 ) 
                   
        if not central :
            central = CentralRule( D , I = IC  , max_step = max_step )
            
        assert isinstance ( central , Rule ) and D == central.D , \
               'Invalid central finite difference!'

        points = tuple ( sorted ( set ( singular ) ) )
        if points :
            if not forward  : forward  = ForwardOpen  ( D , I = IF , max_step = max_step )
            if not backward : backward = BackwardOpen ( D , I = IB , max_step = max_step )
            
        ## initialize the base clas
        DerivativeFD.__init__ ( self                    ,
                                fun        = fun        ,
                                central    = central    ,
                                forward    = forward    , 
                                backward   = backward   , 
                                singular   = points     ,
                                with_error = with_error ,
                                max_step   = max_step   ,
                                step       = step       ,
                                richardson = richardson )
                                

# =============================================================================
## @class Derivative5
#  Very flexible class to calculate 5th derivative using combniation of
#  various rule.
#  - If <code>singular</code> points are specified, in their vinicity
#  special <code>forward</code> and <code>backward</code> are used
#  - Optionally one can apply Richardson's iterative extrapolation
#  to improve results 
#  @code
#  derivative = Derivative5 ( fun )
#  x     = 0.3
#  value = derivative ( x ) 
#  @endcode 
class Derivative5(DerivativeFD) :
    """Very flexible class to calculate 5th derivative using combniation of
    various rule.
    - If `singular` points are specified, in their vinicity
    special `forward` and `backward` rules  are used
    - Optionally one can apply Richardson's iterative extrapolation
    to improve results
    
    >>> derivative = Derivative5 ( fun )
    >>> x     = 0.3
    >>> value = derivative ( x )
    """
    @classprop 
    def IMAX ( cls ) :
        """``IMAX'' : maximal rule index"""
        return CentralRule.IMAX(5)
    
    def __init__ ( self               ,
                   fun                ,
                   step       = 0     ,
                   calc       = 2     , ## rule index 
                   with_error = False ,
                   max_step   = -1    ,
                   central    = None  ,
                   forward    = None  ,
                   backward   = None  ,
                   singular   = ()    ,
                   richardson = 0     ) :

        D  = 5
        IC = min ( calc , CentralRule .IMAX ( D ) - 1 ) 
        IF = min ( calc , ForwardOpen .IMAX ( D ) - 1 ) 
        IB = min ( calc , BackwardOpen.IMAX ( D ) - 1 ) 
                   
        if not central :
            central = CentralRule( D , I = IC  , max_step = max_step )
            
        assert isinstance ( central , Rule ) and D == central.D , \
               'Invalid central finite difference!'

        points = tuple ( sorted ( set ( singular ) ) )
        if points :
            if not forward  : forward  = ForwardOpen  ( D , I = IF , max_step = max_step )
            if not backward : backward = BackwardOpen ( D , I = IB , max_step = max_step )
            
        ## initialize the base clas
        DerivativeFD.__init__ ( self                    ,
                                fun        = fun        ,
                                central    = central    ,
                                forward    = forward    , 
                                backward   = backward   , 
                                singular   = points     ,
                                with_error = with_error ,
                                max_step   = max_step   ,
                                step       = step       ,
                                richardson = richardson )
                                


# =============================================================================
## @class Derivative6
#  Very flexible class to calculate 6th derivative using combniation of
#  various rule.
#  - If <code>singular</code> points are specified, in their vinicity
#  special <code>forward</code> and <code>backward</code> are used
#  - Optionally one can apply Richardson's iterative extrapolation
#  to improve results 
#  @code
#  derivative = Derivative6 ( fun )
#  x     = 0.3
#  value = derivative ( x ) 
#  @endcode 
class Derivative6(DerivativeFD) :
    """Very flexible class to calculate 6th derivative using combniation of
    various rule.
    - If `singular` points are specified, in their vinicity
    special `forward` and `backward` rules  are used
    - Optionally one can apply Richardson's iterative extrapolation
    to improve results
    
    >>> derivative = Derivative6 ( fun )
    >>> x     = 0.3
    >>> value = derivative ( x )
    """
    @classprop 
    def IMAX ( cls ) :
        """``IMAX'' : maximal rule index"""
        return CentralRule.IMAX(6)

    def __init__ ( self               ,
                   fun                ,
                   step       = 0     ,
                   calc       = 2     , ## rule index 
                   with_error = False ,
                   max_step   = -1    ,
                   central    = None  ,
                   forward    = None  ,
                   backward   = None  ,
                   singular   = ()    ,
                   richardson = 0     ) :

        D  = 6
        IC = min ( calc , CentralRule .IMAX ( D ) - 1 ) 
        IF = min ( calc , ForwardOpen .IMAX ( D ) - 1 ) 
        IB = min ( calc , BackwardOpen.IMAX ( D ) - 1 ) 
                   
        if not central :
            central = CentralRule( D , I = IC  , max_step = max_step )
            
        assert isinstance ( central , Rule ) and D == central.D , \
               'Invalid central finite difference!'

        points = tuple ( sorted ( set ( singular ) ) )
        if points :
            if not forward  : forward  = ForwardOpen  ( D , I = IF , max_step = max_step )
            if not backward : backward = BackwardOpen ( D , I = IB , max_step = max_step )
            
        ## initialize the base clas
        DerivativeFD.__init__ ( self                    ,
                                fun        = fun        ,
                                central    = central    ,
                                forward    = forward    , 
                                backward   = backward   , 
                                singular   = points     ,
                                with_error = with_error ,
                                max_step   = max_step   ,
                                step       = step       ,
                                richardson = richardson )
                                
        
# =============================================================================

if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )


# =============================================================================
##                                                                      The END 
# =============================================================================
