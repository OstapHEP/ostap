# Frequently Asked Questions 

## How to use `MINOS`?

The simplest way is to rely on `minos` keyword for `PDF.fitTo` command:
```python
dataset = ...
pdf     = ... 
result, frame = pdf.fitTo ( dataset , draw = True , silent = True , minos = True ) 
``` 
Since in this way `MINOS` is invoked for all parameters, it   can be rather CPU expensive, e.g. for large statistics fits or for the fits with large number of parameters. 

Alternatively one can use more flexible solution:
```python
dataset = ...
pdf     = ... 
result, frame = pdf.fitTo ( dataset , draw = True , silent = True , minos = True ) 
minuit  = pdf.minuit ( dataset )
minuit.migrad ()  
minuit.minos  ( 'A' , 'B' , 'C')
``` 
In this case `MINOS` is invoked only for parameters listed explicitely, and therefore it could be much more CPU efficient
