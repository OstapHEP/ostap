# Utils 

* [ostap.utils](README.md)

Collection of various small utilities 

  - [basic.py](basic.py): small  set of very basic tiny utilities 
     - `isatty` : is `sys.stdout` attached to terminal or not  ?
     - `with_ipython` : do we run `ipython` or *plain* `python` ?
     - `terminal_size` : get the terminal console size
  - [docme.py](docme.py): internal helper utilty for *self-documenting* of ostap modules
  - [gsl.py](gsl.py): set of utilities to play with [`GSL`](https://www.gnu.org/software/gsl/doc/html/index.html) error handlers 
  - [hepdata.py](hepdata.py): functions to convert objects into [`HepDATA`](https://www.hepdata.net) format. Supported types are 
     - `TH1D`, `TH1F`, ...
     - `TGraphErrors` 
     - `TGraphAsymmErrors`
  - [pdg_format.py](pdg_format.py): set of utilities for rounding according to [PDG](http://pdg.lbl.gov) prescription. 
     - Quote: *The basic rule states that if the three highest order digits of the error lie between 100 and 354, we round to two significant digits. If they lie between 355 and 949, we round to one significant digit. Finally, if they lie between 950 and 999, we round up to 1000 and keep two significant digits. In all cases, the central value is given with a precision that matches that of the error.*
     - see [here](http://pdg.lbl.gov/2010/reviews/rpp2010-rev-rpp-intro.pdf)
     - see section 5.3 of [doi:10.1088/0954-3899/33/1/001](http://iopscience.iop.org/issue/0954-3899/33/1)
  - [memory.py](memory.py): useful utilities for memory profiling 
     - It is recommended to install `psutil` module
     - see also [here](http://fa.bianp.net/blog/2013/different-ways-to-get-memory-consumption-or-lessons-learned-from-memory_profiler)
  - [progress_bar.py](progress_bar.py): progress bar and related utilities 
  - [timing.py](timing.py): simple utilities to measure the timing performance 
  - [utils.py](utils.py): an unclassified collection of very different type of utilities
     - *profiler*
     - context managers to keep/preserve/protect various object  *canvas*, `sys.argv`, batch, `ImplicitMT`, etc...
     - ...

 