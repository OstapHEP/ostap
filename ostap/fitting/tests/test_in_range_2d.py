import sys
from ostap.core.pyrouts import *
import ROOT, random, time
import ostap.fitting.roofit 
import ostap.fitting.models as     Models 
from   ostap.core.core      import cpp, VE, dsID
from   ostap.logger.utils   import rooSilent
from   builtins             import range
from ostap.fitting.background import make_bkg 

from   ostap.logger.logger import getLogger 
if '__main__' ==  __name__ : logger = getLogger ( 'in_range' )
else                       : logger = getLogger ( __name__              )



## make simple test mass 
m_x     = ROOT.RooRealVar ( 'm_x' , 'Some test mass(X)' , 0 , 5 )
m_y     = ROOT.RooRealVar ( 'm_y' , 'Some test mass(Y)' , 6 , 10 )

## book very simple data set
varset  = ROOT.RooArgSet  ( m_x , m_y)
dataset = ROOT.RooDataSet ( dsID() , 'Test Data set-1' , varset )  



m1 = VE(3,0.10**2)
m2 = VE(7,0.10**2)

## fill it with three gausissians, 5k events each
N_ss = 5000
N_sb = 1000
N_bs = 500
N_bb = 100


random.seed(0)

## fill it : 5000 events  Gauss * Gauss *Gauss
for i in range(0,N_ss) : 
        m_x.value = m1.gauss() 
        m_y.value = m2.gauss() 
        dataset.add ( varset  )


## fill it : 500 events  Gauss * const * Gauss  
for i in range(0,N_sb) : 
        m_x.value = m1.gauss() 
        m_y.value = random.uniform ( *m_y.minmax() )  

        dataset.add ( varset  )

## fill it : 500 events  const * Gauss * Gauss
for i in range(0,N_bs) : 
        m_x.value = random.uniform ( *m_x.minmax() ) 
        m_y.value = m2.gauss() 
        dataset.add ( varset  )

## fill it : 1000 events  const * const *Gauss
for i in range(0,N_bb) :

        m_x.value  = random.uniform ( *m_x.minmax() )
        m_y.value  = random.uniform ( *m_y.minmax() )
        dataset.add ( varset  )



logger.info ('Dataset: %s' % dataset )  





## various fit components
signal_x1 = Models.Gauss_pdf ( 'G1x'  , xvar = m_x  , mean = m1.value() , sigma = m1.error() )  
signal_y1 = Models.Gauss_pdf ( name='G1y'  , xvar = m_y  , mean = m2.value() , sigma = m2.error()  ) 

bkg_x= make_bkg ( -1 , 'Bx' , m_x )
bkg_y= make_bkg ( -1 , name= 'By' , xvar =m_y )





model = Models.Fit2D (
        name    = 'fit_comp', 
        signal_x    = signal_x1, 
        signal_y    = signal_y1,
        bkg_1x  = bkg_x ,
        bkg_1y  = bkg_y ,
        )
with rooSilent() : 
        ## components
        model.SS.setVal ( 5000 )
        model.SB.setVal ( 1000 )
        model.BS.setVal ( 500 )
        model.BB.setVal ( 100 )

                
        r = model.fitTo ( dataset , ncpu=8 )


def draw_x() :
        dataset.m_y.setRange ( 'fit' , 8,10. )
        model.yvar.setRange ( 'fit' , 8,10. )
        model.draw1(dataset,nbins=200,in_range=(6,8))
        time.sleep (2)
        model.draw1(dataset,nbins=200, in_range='fit')
        time.sleep (2)

def draw_y() :
        dataset.m_x.setRange ( 'fit2' , 0,2.5 )
        model.xvar.setRange ( 'fit2' , 0,2.5 )
        model.draw2(dataset,nbins=200, in_range=(2.5,5))
        time.sleep (2)
        model.draw2(dataset,nbins=200, in_range='fit2')
        time.sleep (2)


if '__main__' == __name__ :

    ## draw x projections
    draw_x() 
    ## draw y projections
    draw_y() 
