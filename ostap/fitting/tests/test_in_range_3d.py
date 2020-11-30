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
if '__main__' ==  __name__ : logger = getLogger ( 'test_in_range_3d' )
else                       : logger = getLogger ( __name__              )



## make simple test mass 
m_x     = ROOT.RooRealVar ( 'm_x' , 'Some test mass(X)' , 0 , 5 )
m_y     = ROOT.RooRealVar ( 'm_y' , 'Some test mass(Y)' , 6 , 10 )
m_z     = ROOT.RooRealVar ( 'm_z' , 'Some test mass(z)' , 10 , 15 )

## book very simple data set
varset  = ROOT.RooArgSet  ( m_x , m_y,m_z )
dataset = ROOT.RooDataSet ( dsID() , 'Test Data set-1' , varset )  



m1 = VE(3,0.10**2)
m2 = VE(7,0.10**2)
m3 = VE(12,0.10**2)

## fill it with three gausissians, 5k events each
N_sss = 5000
N_ssb = 5000
N_sbs = 5000
N_sbb = 1000

N_bss = 500
N_bsb =  100
N_bbs =  100
N_bbb = 250

random.seed(0)

## fill it : 5000 events  Gauss * Gauss *Gauss
for i in range(0,N_sss) : 
        m_x.value = m1.gauss() 
        m_y.value = m2.gauss() 
        m_z.value = m3.gauss() 
        dataset.add ( varset  )


## fill it : 500 events  Gauss * const * Gauss  
for i in range(0,N_ssb) : 
        m_x.value = m1.gauss() 
        m_y.value = random.uniform ( *m_y.minmax() )  
        m_z.value = m3.gauss() 

        dataset.add ( varset  )

## fill it : 500 events  const * Gauss * Gauss
for i in range(0,N_sbs) : 
        m_x.value = random.uniform ( *m_x.minmax() ) 
        m_y.value = m2.gauss() 
        m_z.value = m3.gauss() 
        dataset.add ( varset  )

## fill it : 1000 events  const * const *Gauss
for i in range(0,N_sbb) :

        m_x.value  = random.uniform ( *m_x.minmax() )
        m_y.value  = random.uniform ( *m_y.minmax() )
        m_z.value = m3.gauss() 
        dataset.add ( varset  )
## fill it : 500 events    Gauss * Gauss * const 
for i in range(0,N_bss) : 
        m_x.value = m1.gauss() 
        m_y.value = m2.gauss() 
        m_z.value = random.uniform ( *m_z.minmax() )
        dataset.add ( varset  )
## fill it : 100 events  Gauss * const * const  
for i in range(0,N_bsb) : 
        m_x.value = m1.gauss() 
        m_y.value = random.uniform ( *m_y.minmax() )  
        m_z.value = random.uniform ( *m_y.minmax() )
        dataset.add ( varset  )

## fill it : 100 events  const * Gauss * const
for i in range(0,N_bbs) : 
        m_x.value = random.uniform ( *m_x.minmax() ) 
        m_y.value = m2.gauss() 
        m_z.value = random.uniform ( *m_y.minmax() )
        dataset.add ( varset  )

## fill it : 250 events  const * const * const
for i in range(0,N_bbb) :
        m_x.value  = random.uniform ( *m_x.minmax() )
        m_y.value  = random.uniform ( *m_y.minmax() )
        m_z.value = random.uniform ( *m_y.minmax() )
        dataset.add ( varset  )


logger.info ('Dataset: %s' % dataset )  





## various fit components
signal_x1 = Models.Gauss_pdf ( 'G1x'  , xvar = m_x  , mean = m1.value() , sigma = m1.error() )  
signal_y1 = Models.Gauss_pdf ( name='G1y'  , xvar = m_y  , mean = m2.value() , sigma = m2.error()  ) 
signal_z1 = Models.Gauss_pdf ( name='G1z'  , xvar = m_z  , mean = m3.value() , sigma = m3.error()  )

bkg_x= make_bkg ( -1 , 'Bx' , m_x )
bkg_y= make_bkg ( -1 , name= 'By' , xvar =m_y )
bkg_z= make_bkg ( -1 ,name='Bz' , xvar =m_z )





model = Models.Fit3D (
        name    = 'fit_comp', 
        signal_x    = signal_x1, 
        signal_y    = signal_y1,
        signal_z    = signal_z1,
        bkg_1x  = bkg_x ,
        bkg_1y  = bkg_y ,
        bkg_1z  = bkg_z ,
        )
with rooSilent() : 
        ## components
        model.SSS.setVal ( 5000 )
        model.SBS.setVal ( 5000 )
        model.SSB.setVal ( 5000 )
        model.SBB.setVal ( 1000 )
        model.BSS.setVal ( 500 )
        model.BBS.setVal ( 100 )
        model.BSB.setVal ( 100 )
        model.BBB.setVal ( 250 )

                
        r = model.fitTo ( dataset , ncpu=8 )


def draw_x() :
        dataset.m_y.setRange ( 'fit' , 6,8. )
        model.yvar.setRange ( 'fit' , 6,8. )
        model.draw1(dataset,nbins=200, in_range3=(11,12),in_range2=(8,10))
        time.sleep (1)
        model.draw1(dataset,nbins=200, in_range3=(11,12),in_range2='fit')
        time.sleep (1)  
        model.draw1(dataset,nbins=200, in_range3=(11,12))
        time.sleep (1)
        model.draw1(dataset,nbins=200, in_range2='fit')
        time.sleep (1)

def draw_y() :
        dataset.m_x.setRange ( 'fit2' , 2.5,3. )
        model.xvar.setRange ( 'fit2' , 2.5,3. )
        model.draw2(dataset,nbins=200, in_range3=(11,12),in_range1=(0,3))
        time.sleep (1)
        model.draw2(dataset,nbins=200, in_range3=(11,12),in_range1='fit2')
        time.sleep (1)
        model.draw2(dataset,nbins=200, in_range3=(11,12))
        time.sleep (1)
        model.draw2(dataset,nbins=200, in_range1='fit2')
        time.sleep (1)

def draw_z() :
        dataset.m_x.setRange ( 'fit3' , 2.5,3. )
        model.xvar.setRange ( 'fit3' , 2.5,3. )
        model.draw3(dataset,nbins=200, in_range2=(6,8),in_range1=(0,3))
        time.sleep (1)
        model.draw3(dataset,nbins=200, in_range2=(6,8),in_range1='fit3')
        time.sleep (1)
        model.draw3(dataset,nbins=200, in_range2=(6,8))
        time.sleep (1)
        model.draw3(dataset,nbins=200, in_range1='fit3')
        time.sleep (1)

if '__main__' == __name__ :

    ## draw x projections
    draw_x() 
    ## draw y projections
    draw_y() 
    ## draw z projections
    draw_z() 
