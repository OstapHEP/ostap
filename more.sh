export LCG_DIR=/cvmfs/lhcb.cern.ch/lib/lcg/releases/LCG_88


export OSTAP_DIR=$HOME/cmtuser/ostap

## at CERN-LCG scipy is in pyanalysis 
export PATH=$LCG_DIR/pyanalysis/2.0/x86_64-slc6-gcc62-opt/bin:$PATH
export PYTHONPATH=$LCG_DIR/pyanalysis/2.0/x86_64-slc6-gcc62-opt/lib/python2.7/site-packages:$PYTHONPATH

## at CERN-LCG ipython is in pytools  
export PATH=$LCG_DIR/pytools/2.0/x86_64-slc6-gcc62-opt/bin:$PATH
export PYTHONPATH=$LCG_DIR/pytools/2.0/x86_64-slc6-gcc62-opt/lib/python2.7/site-packages:$PYTHONPATH

export PATH=$PATH:$OSTAP_DIR/scripts
export PYTHONPATH=$OSTAP_DIR:$PYTHONPATH

export LD_LIBRARY_PATH=$OSTAP_DIR/build/lib:$LD_LIBRARY_PATH
