#!/usr/bin/env python
# -*- coding: utf-8 -*-
# =============================================================================
## @file ostap/math/minimize.py
#  Module with some useful utilities for minimization of scalar functiom
#  - a kind of replacement for scipy.minimize.minimize_scalar when scipy is not accessible
#  - the actual code is copied from scipy.minimize  0.18.11
#
#  The main entry point is a function <code>minimizescalar</code>.
#  - a copy from scipy 0.18.11 
# =============================================================================
""" Module with some useful utilities for minimization of scalar functiom
- a kind of replacement for scipy.minimize.minimize_scalar when scipy is not accessible
- the actual code is copied from scipy.minimize  0.18.11

The main entry point is a function <code>minimize_scalar</code>.
- a copy from scipy 0.18.11
"""
# =============================================================================
__version__ = "$Revision:$"
__author__  = "Vanya BELYAEV Ivan.Belyaev@cern.ch"
__date__    = "2018-10-05"
__all__     = (
    'scalar_minimize' , ## local copy of minimize_scalar from scipy 
    'minimize_scalar' , ## the main entry
    ## helper functions:
    'sp_minimum_1D'     , 
    'sp_maximum_1D'     , 
    'sp_minimum_2D'     , 
    'sp_maximum_2D'     , 
    'sp_minimum_3D'     , 
    'sp_maximum_3D'     , 
    )
# =============================================================================
# logging 
# =============================================================================
from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'ostap.math.minimize' )
else                       : logger = getLogger ( __name__              )
# =============================================================================
import math, warnings
from math import sqrt
try :
    import numpy
    import numpy  as np 
    _epsilon = math.sqrt(numpy.finfo(float).eps)
except ImportError :
    class numpy(object) :
        @staticmethod 
        def abs  ( value ) : return abs ( value ) 
        @staticmethod 
        def size ( value ) : return 1 
    import sys 
    _epsilon = sys.float_info.epsilon*0.5
    np =  numpy
# =============================================================================

class OptimizeWarning(UserWarning):
    pass

def is_array_scalar(x):
    """Test whether `x` is either a scalar or an array scalar.

    """
    return np.size(x) == 1

def _check_unknown_options(unknown_options):
    if unknown_options:
        msg = ", ".join(map(str, unknown_options.keys()))
        # Stack level 4: this is called from _minimize_*, which is
        # called from another function in Scipy. Level 4 is the first
        # level in user code.
        warnings.warn("Unknown solver options: %s" % msg, OptimizeWarning, 4)

class OptimizeResult(dict):
    """ Represents the optimization result.

    Attributes
    ----------
    x : ndarray
        The solution of the optimization.
    success : bool
        Whether or not the optimizer exited successfully.
    status : int
        Termination status of the optimizer. Its value depends on the
        underlying solver. Refer to `message` for details.
    message : str
        Description of the cause of the termination.
    fun, jac, hess: ndarray
        Values of objective function, its Jacobian and its Hessian (if
        available). The Hessians may be approximations, see the documentation
        of the function in question.
    hess_inv : object
        Inverse of the objective function's Hessian; may be an approximation.
        Not available for all solvers. The type of this attribute may be
        either np.ndarray or scipy.sparse.linalg.LinearOperator.
    nfev, njev, nhev : int
        Number of evaluations of the objective functions and of its
        Jacobian and Hessian.
    nit : int
        Number of iterations performed by the optimizer.
    maxcv : float
        The maximum constraint violation.

    Notes
    -----
    There may be additional attributes not listed above depending of the
    specific solver. Since this class is essentially a subclass of dict
    with attribute accessors, one can see which attributes are available
    using the `keys()` method.
    """
    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError(name)

    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __repr__(self):
        if self.keys():
            m = max(map(len, list(self.keys()))) + 1
            return '\n'.join([k.rjust(m) + ': ' + repr(v)
                              for k, v in sorted(self.items())])
        else:
            return self.__class__.__name__ + "()"

    def __dir__(self):
        return list(self.keys())



def scalar_minimize (fun, bracket=None, bounds=None, args=(),
                     method='brent', tol=None, options=None):
    """Minimization of scalar function of one variable.

    Parameters
    ----------
    fun : callable
        Objective function.
        Scalar function, must return a scalar.
    bracket : sequence, optional
        For methods 'brent' and 'golden', `bracket` defines the bracketing
        interval and can either have three items `(a, b, c)` so that `a < b
        < c` and `fun(b) < fun(a), fun(c)` or two items `a` and `c` which
        are assumed to be a starting interval for a downhill bracket search
        (see `bracket`); it doesn't always mean that the obtained solution
        will satisfy `a <= x <= c`.
    bounds : sequence, optional
        For method 'bounded', `bounds` is mandatory and must have two items
        corresponding to the optimization bounds.
    args : tuple, optional
        Extra arguments passed to the objective function.
    method : str or callable, optional
        Type of solver.  Should be one of

            - 'Brent'     :ref:`(see here) <optimize.minimize_scalar-brent>`
            - 'Bounded'   :ref:`(see here) <optimize.minimize_scalar-bounded>`
            - 'Golden'    :ref:`(see here) <optimize.minimize_scalar-golden>`
            - custom - a callable object (added in version 0.14.0),
              see below
    tol : float, optional
        Tolerance for termination. For detailed control, use solver-specific
        options.
    options : dict, optional
        A dictionary of solver options.
            maxiter : int
                Maximum number of iterations to perform.
            disp : bool
                Set to True to print convergence messages.

        See :func:`show_options()` for solver-specific options.

    Returns
    -------
    res : OptimizeResult
        The optimization result represented as a ``OptimizeResult`` object.
        Important attributes are: ``x`` the solution array, ``success`` a
        Boolean flag indicating if the optimizer exited successfully and
        ``message`` which describes the cause of the termination. See
        `OptimizeResult` for a description of other attributes.

    See also
    --------
    minimize : Interface to minimization algorithms for scalar multivariate
        functions
    show_options : Additional options accepted by the solvers

    Notes
    -----
    This section describes the available solvers that can be selected by the
    'method' parameter. The default method is *Brent*.
    Method :ref:`Brent <optimize.minimize_scalar-brent>` uses Brent's
    algorithm to find a local minimum.  The algorithm uses inverse
    parabolic interpolation when possible to speed up convergence of
    the golden section method.

    Method :ref:`Golden <optimize.minimize_scalar-golden>` uses the
    golden section search technique. It uses analog of the bisection
    method to decrease the bracketed interval. It is usually
    preferable to use the *Brent* method.

    Method :ref:`Bounded <optimize.minimize_scalar-bounded>` can
    perform bounded minimization. It uses the Brent method to find a
    local minimum in the interval x1 < xopt < x2.

    **Custom minimizers**

    It may be useful to pass a custom minimization method, for example
    when using some library frontend to minimize_scalar.  You can simply
    pass a callable as the ``method`` parameter.

    The callable is called as ``method(fun, args, **kwargs, **options)``
    where ``kwargs`` corresponds to any other parameters passed to `minimize`
    (such as `bracket`, `tol`, etc.), except the `options` dict, which has
    its contents also passed as `method` parameters pair by pair.  The method
    shall return an ``OptimizeResult`` object.

    The provided `method` callable must be able to accept (and possibly ignore)
    arbitrary parameters; the set of parameters accepted by `minimize` may
    expand in future versions and then these parameters will be passed to
    the method.  You can find an example in the scipy.optimize tutorial.

    .. versionadded:: 0.11.0

    Examples
    --------
    Consider the problem of minimizing the following function.

    >>> def f(x):
    ...     return (x - 2) * x * (x + 2)**2

    Using the *Brent* method, we find the local minimum as:

    >>> from scipy.optimize import minimize_scalar
    >>> res = minimize_scalar(f)
    >>> res.x
    1.28077640403

    Using the *Bounded* method, we find a local minimum with specified
    bounds as:

    >>> res = minimize_scalar(f, bounds=(-3, -1), method='bounded')
    >>> res.x
    -2.0000002026

    """
    if not isinstance(args, tuple):
        args = (args,)

    if callable(method):
        meth = "_custom"
    else:
        meth = method.lower()
    if options is None:
        options = {}
        
    if tol is not None:
        options = dict(options)
        if meth == 'bounded' and 'xatol' not in options:
            warn("Method 'bounded' does not support relative tolerance in x; "
                 "defaulting to absolute tolerance.", RuntimeWarning)
            options['xatol'] = tol
        elif meth == '_custom':
            options.setdefault('tol', tol)
        else:
            options.setdefault('xtol', tol)

    if meth == '_custom':
        return method(fun, args=args, bracket=bracket, bounds=bounds, **options)
    elif meth == 'brent':
        return _minimize_scalar_brent(fun, bracket, args, **options)
    elif meth == 'bounded':
        if bounds is None:
            raise ValueError('The `bounds` parameter is mandatory for '
                             'method `bounded`.')
        return _minimize_scalar_bounded(fun, bounds, args, **options)
    elif meth == 'golden':
        return _minimize_scalar_golden(fun, bracket, args, **options)
    else:
        raise ValueError('Unknown solver %s' % method)


def _minimize_scalar_brent(func, brack=None, args=(),
                           xtol=1.48e-8, maxiter=500,
                           **unknown_options):
    """
    Options
    -------
    maxiter : int
        Maximum number of iterations to perform.
    xtol : float
        Relative error in solution `xopt` acceptable for convergence.

    Notes
    -----
    Uses inverse parabolic interpolation when possible to speed up
    convergence of golden section method.

    """
    _check_unknown_options(unknown_options)
    tol = xtol
    if tol < 0:
        raise ValueError('tolerance should be >= 0, got %r' % tol)

    brent = Brent(func=func, args=args, tol=tol,
                  full_output=True, maxiter=maxiter)
    brent.set_bracket(brack)
    brent.optimize()
    x, fval, nit, nfev = brent.get_result(full_output=True)
    return OptimizeResult(fun=fval, x=x, nit=nit, nfev=nfev,
                          success=nit < maxiter)

class Brent:
    #need to rethink design of __init__
    def __init__(self, func, args=(), tol=1.48e-8, maxiter=500,
                 full_output=0):
        self.func = func
        self.args = args
        self.tol = tol
        self.maxiter = maxiter
        self._mintol = 1.0e-11
        self._cg = 0.3819660
        self.xmin = None
        self.fval = None
        self.iter = 0
        self.funcalls = 0

    # need to rethink design of set_bracket (new options, etc)
    def set_bracket(self, brack=None):
        self.brack = brack

    def get_bracket_info(self):
        #set up
        func = self.func
        args = self.args
        brack = self.brack
        ### BEGIN core bracket_info code ###
        ### carefully DOCUMENT any CHANGES in core ##
        if brack is None:
            xa, xb, xc, fa, fb, fc, funcalls = bracket(func, args=args)
        elif len(brack) == 2:
            xa, xb, xc, fa, fb, fc, funcalls = bracket(func, xa=brack[0],
                                                       xb=brack[1], args=args)
        elif len(brack) == 3:
            xa, xb, xc = brack
            if (xa > xc):  # swap so xa < xc can be assumed
                xc, xa = xa, xc
            if not ((xa < xb) and (xb < xc)):
                raise ValueError("Not a bracketing interval.")
            fa = func(*((xa,) + args))
            fb = func(*((xb,) + args))
            fc = func(*((xc,) + args))
            if not ((fb < fa) and (fb < fc)):
                raise ValueError("Not a bracketing interval.")
            funcalls = 3
        else:
            raise ValueError("Bracketing interval must be "
                             "length 2 or 3 sequence.")
        ### END core bracket_info code ###

        return xa, xb, xc, fa, fb, fc, funcalls
    def optimize(self):
        # set up for optimization
        func = self.func
        xa, xb, xc, fa, fb, fc, funcalls = self.get_bracket_info()
        _mintol = self._mintol
        _cg = self._cg
        #################################
        #BEGIN CORE ALGORITHM
        #################################
        x = w = v = xb
        fw = fv = fx = func(*((x,) + self.args))
        if (xa < xc):
            a = xa
            b = xc
        else:
            a = xc
            b = xa
        deltax = 0.0
        funcalls = 1
        iter = 0
        while (iter < self.maxiter):
            tol1 = self.tol * numpy.abs(x) + _mintol
            tol2 = 2.0 * tol1
            xmid = 0.5 * (a + b)
            # check for convergence
            if numpy.abs(x - xmid) < (tol2 - 0.5 * (b - a)):
                break
            # XXX In the first iteration, rat is only bound in the true case
            # of this conditional. This used to cause an UnboundLocalError
            # (gh-4140). It should be set before the if (but to what?).
            if (numpy.abs(deltax) <= tol1):
                if (x >= xmid):
                    deltax = a - x       # do a golden section step
                else:
                    deltax = b - x
                rat = _cg * deltax
            else:                              # do a parabolic step
                tmp1 = (x - w) * (fx - fv)
                tmp2 = (x - v) * (fx - fw)
                p = (x - v) * tmp2 - (x - w) * tmp1
                tmp2 = 2.0 * (tmp2 - tmp1)
                if (tmp2 > 0.0):
                    p = -p
                tmp2 = numpy.abs(tmp2)
                dx_temp = deltax
                deltax = rat
                # check parabolic fit
                if ((p > tmp2 * (a - x)) and (p < tmp2 * (b - x)) and
                        (numpy.abs(p) < numpy.abs(0.5 * tmp2 * dx_temp))):
                    rat = p * 1.0 / tmp2        # if parabolic step is useful.
                    u = x + rat
                    if ((u - a) < tol2 or (b - u) < tol2):
                        if xmid - x >= 0:
                            rat = tol1
                        else:
                            rat = -tol1
                else:
                    if (x >= xmid):
                        deltax = a - x  # if it's not do a golden section step
                    else:
                        deltax = b - x
                    rat = _cg * deltax
            if (numpy.abs(rat) < tol1):            # update by at least tol1
                if rat >= 0:
                    u = x + tol1
                else:
                    u = x - tol1
            else:
                u = x + rat
            fu = func(*((u,) + self.args))      # calculate new output value
            funcalls += 1

            if (fu > fx):                 # if it's bigger than current
                if (u < x):
                    a = u
                else:
                    b = u
                if (fu <= fw) or (w == x):
                    v = w
                    w = u
                    fv = fw
                    fw = fu
                elif (fu <= fv) or (v == x) or (v == w):
                    v = u
                    fv = fu
            else:
                if (u >= x):
                    a = x
                else:
                    b = x
                v = w
                w = x
                x = u
                fv = fw
                fw = fx
                fx = fu

            iter += 1
        #################################
        #END CORE ALGORITHM
        #################################

        self.xmin = x
        self.fval = fx
        self.iter = iter
        self.funcalls = funcalls

    def get_result(self, full_output=False):
        if full_output:
            return self.xmin, self.fval, self.iter, self.funcalls
        else:
            return self.xmin

def _minimize_scalar_bounded(func, bounds, args=(),
                             xatol=1e-5, maxiter=500, disp=0,
                             **unknown_options):
    """
    Options
    -------
    maxiter : int
        Maximum number of iterations to perform.
    disp : bool
        Set to True to print convergence messages.
    xatol : float
        Absolute error in solution `xopt` acceptable for convergence.

    """
    _check_unknown_options(unknown_options)
    maxfun = maxiter
    # Test bounds are of correct form
    if len(bounds) != 2:
        raise ValueError('bounds must have two elements.')
    x1, x2 = bounds

    if not (is_array_scalar(x1) and is_array_scalar(x2)):
        raise ValueError("Optimisation bounds must be scalars"
                         " or array scalars.")
    if x1 > x2:
        raise ValueError("The lower bound exceeds the upper bound.")

    flag = 0
    header = ' Func-count     x          f(x)          Procedure'
    step = '       initial'

    sqrt_eps = sqrt(2.2e-16)
    golden_mean = 0.5 * (3.0 - sqrt(5.0))
    a, b = x1, x2
    fulc = a + golden_mean * (b - a)
    nfc, xf = fulc, fulc
    rat = e = 0.0
    x = xf
    fx = func(x, *args)
    num = 1
    fmin_data = (1, xf, fx)

    ffulc = fnfc = fx
    xm = 0.5 * (a + b)
    tol1 = sqrt_eps * numpy.abs(xf) + xatol / 3.0
    tol2 = 2.0 * tol1

    if disp > 2:
        print(" ")
        print(header)
        print("%5.0f   %12.6g %12.6g %s" % (fmin_data + (step,)))
    while (numpy.abs(xf - xm) > (tol2 - 0.5 * (b - a))):
        golden = 1
        # Check for parabolic fit
        if numpy.abs(e) > tol1:
            golden = 0
            r = (xf - nfc) * (fx - ffulc)
            q = (xf - fulc) * (fx - fnfc)
            p = (xf - fulc) * q - (xf - nfc) * r
            q = 2.0 * (q - r)
            if q > 0.0:
                p = -p
            q = numpy.abs(q)
            r = e
            e = rat

            # Check for acceptability of parabola
            if ((numpy.abs(p) < numpy.abs(0.5*q*r)) and (p > q*(a - xf)) and
                    (p < q * (b - xf))):
                rat = (p + 0.0) / q
                x = xf + rat
                step = '       parabolic'

                if ((x - a) < tol2) or ((b - x) < tol2):
                    si = numpy.sign(xm - xf) + ((xm - xf) == 0)
                    rat = tol1 * si
            else:      # do a golden section step
                golden = 1

        if golden:  # Do a golden-section step
            if xf >= xm:
                e = a - xf
            else:
                e = b - xf
            rat = golden_mean*e
            step = '       golden'

        si = numpy.sign(rat) + (rat == 0)
        x = xf + si * numpy.max([numpy.abs(rat), tol1])
        fu = func(x, *args)
        num += 1
        fmin_data = (num, x, fu)
        if disp > 2:
            print("%5.0f   %12.6g %12.6g %s" % (fmin_data + (step,)))

        if fu <= fx:
            if x >= xf:
                a = xf
            else:
                b = xf
            fulc, ffulc = nfc, fnfc
            nfc, fnfc = xf, fx
            xf, fx = x, fu
        else:
            if x < xf:
                a = x
            else:
                b = x
            if (fu <= fnfc) or (nfc == xf):
                fulc, ffulc = nfc, fnfc
                nfc, fnfc = x, fu
            elif (fu <= ffulc) or (fulc == xf) or (fulc == nfc):
                fulc, ffulc = x, fu

        xm = 0.5 * (a + b)
        tol1 = sqrt_eps * numpy.abs(xf) + xatol / 3.0
        tol2 = 2.0 * tol1

        if num >= maxfun:
            flag = 1
            break

    fval = fx
    if disp > 0:
        _endprint(x, flag, fval, maxfun, xatol, disp)

    result = OptimizeResult(fun=fval, status=flag, success=(flag == 0),
                            message={0: 'Solution found.',
                                     1: 'Maximum number of function calls '
                                        'reached.'}.get(flag, ''),
                            x=xf, nfev=num)

    return result


def _minimize_scalar_golden(func, brack=None, args=(),
                            xtol=_epsilon, **unknown_options):
    """
    Options
    -------
    maxiter : int
        Maximum number of iterations to perform.
    xtol : float
        Relative error in solution `xopt` acceptable for convergence.

    """
    _check_unknown_options(unknown_options)
    tol = xtol
    if brack is None:
        xa, xb, xc, fa, fb, fc, funcalls = bracket(func, args=args)
    elif len(brack) == 2:
        xa, xb, xc, fa, fb, fc, funcalls = bracket(func, xa=brack[0],
                                                   xb=brack[1], args=args)
    elif len(brack) == 3:
        xa, xb, xc = brack
        if (xa > xc):  # swap so xa < xc can be assumed
            xc, xa = xa, xc
        if not ((xa < xb) and (xb < xc)):
            raise ValueError("Not a bracketing interval.")
        fa = func(*((xa,) + args))
        fb = func(*((xb,) + args))
        fc = func(*((xc,) + args))
        if not ((fb < fa) and (fb < fc)):
            raise ValueError("Not a bracketing interval.")
        funcalls = 3
    else:
        raise ValueError("Bracketing interval must be length 2 or 3 sequence.")

    _gR = 0.61803399
    _gC = 1.0 - _gR
    x3 = xc
    x0 = xa
    if (numpy.abs(xc - xb) > numpy.abs(xb - xa)):
        x1 = xb
        x2 = xb + _gC * (xc - xb)
    else:
        x2 = xb
        x1 = xb - _gC * (xb - xa)
    f1 = func(*((x1,) + args))
    f2 = func(*((x2,) + args))
    funcalls += 2
    while (numpy.abs(x3 - x0) > tol * (numpy.abs(x1) + numpy.abs(x2))):
        if (f2 < f1):
            x0 = x1
            x1 = x2
            x2 = _gR * x1 + _gC * x3
            f1 = f2
            f2 = func(*((x2,) + args))
        else:
            x3 = x2
            x2 = x1
            x1 = _gR * x2 + _gC * x0
            f2 = f1
            f1 = func(*((x1,) + args))
        funcalls += 1
    if (f1 < f2):
        xmin = x1
        fval = f1
    else:
        xmin = x2
        fval = f2

    return OptimizeResult(fun=fval, nfev=funcalls, x=xmin)

def bracket(func, xa=0.0, xb=1.0, args=(), grow_limit=110.0, maxiter=1000):
    """
    Bracket the minimum of the function.

    Given a function and distinct initial points, search in the
    downhill direction (as defined by the initital points) and return
    new points xa, xb, xc that bracket the minimum of the function
    f(xa) > f(xb) < f(xc). It doesn't always mean that obtained
    solution will satisfy xa<=x<=xb

    Parameters
    ----------
    func : callable f(x,*args)
        Objective function to minimize.
    xa, xb : float, optional
        Bracketing interval. Defaults `xa` to 0.0, and `xb` to 1.0.
    args : tuple, optional
        Additional arguments (if present), passed to `func`.
    grow_limit : float, optional
        Maximum grow limit.  Defaults to 110.0
    maxiter : int, optional
        Maximum number of iterations to perform. Defaults to 1000.

    Returns
    -------
    xa, xb, xc : float
        Bracket.
    fa, fb, fc : float
        Objective function values in bracket.
    funcalls : int
        Number of function evaluations made.

    """
    _gold = 1.618034
    _verysmall_num = 1e-21
    fa = func(*(xa,) + args)
    fb = func(*(xb,) + args)
    if (fa < fb):                      # Switch so fa > fb
        xa, xb = xb, xa
        fa, fb = fb, fa
    xc = xb + _gold * (xb - xa)
    fc = func(*((xc,) + args))
    funcalls = 3
    iter = 0
    while (fc < fb):
        tmp1 = (xb - xa) * (fb - fc)
        tmp2 = (xb - xc) * (fb - fa)
        val = tmp2 - tmp1
        if numpy.abs(val) < _verysmall_num:
            denom = 2.0 * _verysmall_num
        else:
            denom = 2.0 * val
        w = xb - ((xb - xc) * tmp2 - (xb - xa) * tmp1) / denom
        wlim = xb + grow_limit * (xc - xb)
        if iter > maxiter:
            raise RuntimeError("Too many iterations.")
        iter += 1
        if (w - xc) * (xb - w) > 0.0:
            fw = func(*((w,) + args))
            funcalls += 1
            if (fw < fc):
                xa = xb
                xb = w
                fa = fb
                fb = fw
                return xa, xb, xc, fa, fb, fc, funcalls
            elif (fw > fb):
                xc = w
                fc = fw
                return xa, xb, xc, fa, fb, fc, funcalls
            w = xc + _gold * (xc - xb)
            fw = func(*((w,) + args))
            funcalls += 1
        elif (w - wlim)*(wlim - xc) >= 0.0:
            w = wlim
            fw = func(*((w,) + args))
            funcalls += 1
        elif (w - wlim)*(xc - w) > 0.0:
            fw = func(*((w,) + args))
            funcalls += 1
            if (fw < fc):
                xb = xc
                xc = w
                w = xc + _gold * (xc - xb)
                fb = fc
                fc = fw
                fw = func(*((w,) + args))
                funcalls += 1
        else:
            w = xc + _gold * (xc - xb)
            fw = func(*((w,) + args))
            funcalls += 1
        xa = xb
        xb = xc
        xc = w
        fa = fb
        fb = fc
        fc = fw
    return xa, xb, xc, fa, fb, fc, funcalls


try :
    from scipy.optimize import minimize_scalar as ms 
    minimize_scalar = ms  
    scipy_available = True 
except ImportError :
    minimize_scalar = scalar_minimize 
    scipy_available = False 

# =============================================================================


if not scipy_available :
    
    sp_minimum_1D = None
    sp_maximum_1D = None
    sp_minimum_2D = None
    sp_maximum_2D = None
    sp_minimum_3D = None
    sp_maximum_3D = None
    
else :
        
    # =========================================================================
    ## get a minimum for 1D-function
    #  @code
    #  model = ...
    #  x = model.minimum() 
    #  @endcode 
    def sp_minimum_1D ( fun , xmin , xmax , x0 = None , *args ) :
        """Get a minimum for 1D-function
        >>> model = ...
        >>> x = model.minimum () 
        >>>
        """
        if x0 == None : x0 = 0.5 * ( xmin + xmax )
        
        import numpy as np
        x0     = np.array ( [ x0 ] )
        
        bounds = [ ( xmin , xmax ) ]
        
        import scipy.optimize as spo
        res    = spo.minimize ( fun , x0 = x0 , bounds = bounds )
        if not res.success :
            logger.error ( "Can't minimize the function: %s" % res.message )
        return res.x[0]
        
    # =========================================================================
    ## get a maximum for 1D-function
    #  @code
    #  model = ...
    #  x = model.maximum() 
    #  @endcode 
    def sp_maximum_1D ( fun , xmin , xmax , x0 = None , *args ) :
        """Get a maximum for 1D-function
        >>> model = ...
        >>> x = model.maximum () 
        >>>
        """
        funmin = lambda x , *a : -1.0 * ( float ( fun ( x , *a ) ) )
        return sp_minimum_1D ( funmin , xmin ,  xmax , x0 , *args )
    
    # =========================================================================
    ## get a minimum for 2D-function
    #  @code
    #  model2 = ...
    #  x , y = model2.minimum () 
    #  @endcode 
    def sp_minimum_2D ( fun  ,
                        xmin , xmax ,
                        ymin , ymax , x0 = () , *args ) :
        """Get a maximum for 2D-function
        >>> model2 = ...
        >>> x , y = model2.maximum() 
        >>>
        """
        if not x0 :  x0 = 0.5 * ( xmin + xmax ) , 0.5 * ( ymin + ymax ) 
        import numpy as np
        x0     = np.array ( *x0  )
        
        bounds = [ ( xmin , xmax ) , ( ymin , ymax ) ]
        
        import scipy.optimize as spo
        res    = spo.minimize ( fun , x0 = x0 , bounds = bounds )
        if not res.success :
            logger.error ( "Can't minimize the function: %s" % res.message )
        return res.x[0] , res.x[1] 
            
    # =========================================================================
    ## get a maximum for 2D-function
    #  @code
    #  model2 = ...
    #  x , y = model2.maximum() 
    #  @endcode 
    def sp_maximum_2D ( fun ,
                        xmin , xmax ,
                        ymin , ymax , x0 = () , *args ) :
        """Get a maximum for 2D-function
        >>> model2 = ...
        >>> x , y = model2.maximum () 
        >>>
        """
        funmin = lambda x , y , *a : -1.0 * ( float ( fun ( x , y , *a ) ) )
        return sp_minimum_2D  ( funmin ,
                                xmin , xmax ,
                                ymin , ymax , x0 , *args )
    
    # =========================================================================
    ## get a minimum for 3D-function
    #  @code
    #  model3 = ...
    #  x , y , z = model2.minimum () 
    #  @endcode 
    def sp_minimum_3D ( fun  ,
                        xmin , xmax ,
                        ymin , ymax ,
                        zmin , zmax , x0 = () , *args ) :
        """Get a minimum for 3D-function
        >>> model3 = ...
        >>> x , y , z = model3.minimum() 
        >>>
        """
        if not x0 :  x0 = 0.5 * ( xmin + xmax ) , 0.5 * ( ymin + ymax ) , 0.5 * ( zmin + zmax ) 
        import numpy as np
        x0     = np.array ( *x0  )
        
        bounds = [ ( xmin , xmax ) , ( ymin , ymax ) , ( zmin , zmax ) ]
        
        import scipy.optimize as spo
        res    = spo.minimize ( fun , x0 = x0 , bounds = bounds )
        if not res.success :
            logger.error ( "Can't minimize the function: %s" % res.message )
        return res.x[0] , res.x[1] , res.x[2] 
        
    # =========================================================================
    ## get a maximum for 3D-function
    #  @code
    #  model3 = ...
    #  x , y , z = model3.maximum() 
    #  @endcode 
    def sp_maximum_3D ( fun ,
                        xmin , xmax ,
                        ymin , ymax ,
                        zmin , zmax , x0 = () , *args ) :
        """Get a maximum for 3D-function
        >>> model3 = ...
        >>> x, y , z  = model3.maximum () 
        >>>
        """
        funmin = lambda x , y , z , *a : -1.0 * ( float ( fun ( x , y , z , *a ) ) )
        return sp_minimum_3D  ( funmin ,
                                xmin ,  xmax ,
                                ymin ,  ymax ,
                                zmin ,  zmax , x0 , *args )
    
# =============================================================================
if '__main__' == __name__ :
    
    from ostap.utils.docme import docme
    docme ( __name__ , logger = logger )

# =============================================================================
# The END 
# =============================================================================
