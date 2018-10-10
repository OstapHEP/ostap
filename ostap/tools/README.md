# Tools

* [ostap.tools](README.md)

Collection of high-level analysis tools:

  - [tmva.py](tmva.py): `Trainer` and `Reader`, the code to simplify communications with [`TMVA`](https://root.cern.ch/tmva), *The Toolkit for Multivariate Data Analysis with ROOT*.
     - small modification of the original code by *Albert PUIG*
  - [chopping.py](chopping.py): *chopping* utulity for `TMVA`
     - splits sample into several non-overlapping categories and performs separate `TMVA` training/testing/using for subsamples.
     - it is needed to  avoid of usage of the same events in `TMVA` training and using 
     - `Trainer` and `Reader` classes, the interface is alsmost identical to *simple* `TMVA` case 
  - [reweight.py](reweight.py) : tool for relatively easy *reweighting*
     - getting as the input the desired distributions for e.g. *data* sample, performs *reweighting* of *simulated* sample to achive the same shapes  
     - input *data* distributions are specified as collection of 1D or 2D distributions
 