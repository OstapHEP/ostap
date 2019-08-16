# Visualization of fit results 


## Drawing options 
Consider the simple case of 1D fit with subsequent visualization of results 
```python
result , _ = pdf.fitTo ( dataset )
frame = pdf.draw ( dataset , nbins = 100 , ... ) 
```
A flexible control on the vizualiation options is provided via canse-insensitive (an underscrore-blind) keyword arguments: 

 - `data_options`                 : data options to draw the data points 
 - `signal_options`               : draw options for 'signal '      component(s)
 - `background_options`           : draw options for 'background'   component(s)
 - `background2D_options`         : draw options for 'background2D' component(s)
 - `crossterm1_options`           : draw options for 'crossterm1'   component(s)
 - `crossterm2_options`           : draw options for 'crossterm2'   component(s)    
 - `component_options`            : draw options for 'other'        component(s)
 - `total_fit_options`            : draw options for the total fit curve
All these options are the tuple/lists of `ROOT.RooFit.RooCmdArg` objects, that are directly transferred to 
corresponding `ROOT.RooXXX.plotOn`-method: 
```python
frame = pdf.draw ( dataset , nbins = 100 , data_options = ( ROOT.RooFit.MarkerStyle ( 20   ) , 
                                                            ROOT.RooFit.DrawOption  ( "zp" ) , 
                                                            ROOT.RooFit.MarkerSize  ( 0.5  ) ) )
```
The options can be stored in the pdf and later reused for all subsequent visualiztaions:
```python
pdf.draw_options['data_options'] = ( ROOT.RooFit.MarkerStyle ( 20   ) , 
                                     ROOT.RooFit.DrawOption  ( "zp" ) , 
                                     ROOT.RooFit.MarkerSize  ( 0.5  ) )
frame = pdf.draw ( dataset , nbins = 100 )
```

## Drawing styles 

In addition to the drawing options described above one can have drawing styles, defined  for certain fit components, e.g. `signal_style`, `background_style`,... The overall vizualization of fit components is a combination of options and styles. Styles can be specified using `ROOT.RooFit.RooCmdArg` object, but also in a parametric way using `Line`, `Area`, `Style` and `Styles`objects defiend in [fitdraw.py](fitdraw.py). It is convinient when the fit contains the several components for the same category. E.g. for the fit what constains  three "signal" components, they can be visualized as three solid areas :
```python
frame = pdf.draw ( dataset , nbins = 100 ,
   signal_style = ( Area ( fillcolor = 2 , fillstyle = 1001 ) , 
                    Area ( fillcolor = 4 , fillstyle = 1001 ) , 
                    Area ( fillcolor = 8 , fillstyle = 1001 ) ) ) 
```    
Or, e.g. for two background components 
```python
frame = pdf.draw ( dataset , nbins = 100 ,
   background_style = ( Line ( linecolor = 2 , linestyle = 2 ) , 
                        Line ( linecolor = 4 , linestyle = 3 ) ) )
```    
The styles also can be stored and reused later:
```python
pdf.draw_options['background_style'] = ( Line ( linecolor = 2 , linestyle = 2 ) , 
                                         Line ( linecolor = 4 , linestyle = 3 ) )
frame = pdf.draw ( dataset , nbins = 100 )
```    

Several styles are predefined 
 - `signal_style`       , style used for visualisation of "signal"        component(s)
 - `background_style`   , style used for visualisation of "background"    component(s)
 - `background2D_style` , style used for visualization of "background-2D" component(s)
 - `crossterm1_style`   , style used for visualization of "crossterm1"    component(s)
 - `crossterm2_style`   , style used for visualization of "crossterm2"    component(s)
 - `component_style`    , style used for visualization of other, non-classified  component(s)

The default styles are :
```python
default_signal_style  = (
    Style ( linecolor = ROOT.kRed      , linewidth = 2 , fillcolor = ROOT.kRed     , fillstyle = 1001 ) ,
    Style ( linecolor = ROOT.kBlue     , linewidth = 2 , fillcolor = ROOT.kBlue    , fillstyle = 1001 ) ,
    Style ( linecolor =    8           , linewidth = 2 , fillcolor =    8          , fillstyle = 1001 ) ,    
    Style ( linecolor = ROOT.kMagenta  , linewidth = 2 , fillcolor = ROOT.kMagenta , fillstyle = 1001 ) ,
    Style ( linecolor = ROOT.kCyan     , linewidth = 2 , fillcolor = ROOT.kCyan    , fillstyle = 1001 ) ,
    Style ( linecolor = ROOT.kYellow   , linewidth = 2 , fillcolor = ROOT.kYellow  , fillstyle = 1001 ) ,
    ) 

default_background_style = (
    Line  ( linecolor = ROOT.kBlue         , linestyle =  7 ) ,
    Line  ( linecolor = ROOT.kBlue    -  9 , linestyle = 11 ) ,
    Line  ( linecolor = ROOT.kBlue    +  3 , linestyle = 12 ) ,
    Line  ( linecolor = ROOT.kBlue    -  2 , linestyle = 13 ) ,
    Line  ( linecolor = ROOT.kBlue    - 10 , linestyle = 14 ) ,
    )

default_component_style  = (
    Style ( linecolor = ROOT.kOrange , fillcolor = ROOT.kOrange , fillstyle = 3345 ) ,
    Style ( linecolor = ROOT.kOrange , fillcolor = ROOT.kOrange , fillstyle = 3354 ) ,
    Style ( linecolor = ROOT.kOrange , fillcolor = ROOT.kOrange , fillstyle = 3305 ) ,
    Style ( linecolor = ROOT.kOrange , fillcolor = ROOT.kOrange , fillstyle = 3395 ) ,
    Style ( linecolor = ROOT.kOrange , fillcolor = ROOT.kOrange , fillstyle = 3422 ) ,
    Style ( linecolor = ROOT.kOrange , fillcolor = ROOT.kOrange , fillstyle = 3477 ) ,
    Style ( linecolor = ROOT.kOrange , fillcolor = ROOT.kOrange , fillstyle = 3544 ) ,
    Style ( linecolor = ROOT.kOrange , fillcolor = ROOT.kOrange , fillstyle = 3590 ) ,
    Style ( linecolor = ROOT.kOrange , fillcolor = ROOT.kOrange , fillstyle = 3509 ) ,
    )

default_crossterm1_style = (
    Line  ( linecolor = ROOT.kMagenta +  3 , linestyle = 11 ) ,
    Line  ( linecolor = ROOT.kMagenta - 10 , linestyle = 12 ) , 
    Line  ( linecolor = ROOT.kMagenta +  3 , linestyle = 13 ) , 
    Line  ( linecolor = ROOT.kMagenta -  3 , linestyle = 14 ) ,
    )

default_crossterm2_style = (
    Line  ( linecolor = ROOT.kGreen   +  1 , linestyle = 14 ) ,
    Line  ( linecolor = ROOT.kGreen   -  1 , linestyle = 13 ) ,
    Line  ( linecolor = ROOT.kGreen   - 10 , linestyle = 12 ) ,
    Line  ( linecolor = ROOT.kGreen   +  1 , linestyle = 11 ) ,
    ) 

default_background2D_style = default_background_style 

```

## Configuration file
All styles and options can redefined via the section `FitDraw` in the 
configuration files, e.g. 
```
[Fit Draw]
signal_style  = (
    Style ( linecolor = ROOT.kRed      , linewidth = 2 , fillcolor = ROOT.kRed     , fillstyle = 1001 ) ,
    Style ( linecolor = ROOT.kBlue     , linewidth = 2 , fillcolor = ROOT.kBlue    , fillstyle = 1001 ) ,
    Style ( linecolor =    8           , linewidth = 2 , fillcolor =    8          , fillstyle = 1001 ) ,    
    Style ( linecolor = ROOT.kMagenta  , linewidth = 2 , fillcolor = ROOT.kMagenta , fillstyle = 1001 ) ,
    Style ( linecolor = ROOT.kCyan     , linewidth = 2 , fillcolor = ROOT.kCyan    , fillstyle = 1001 ) ,
    Style ( linecolor = ROOT.kYellow   , linewidth = 2 , fillcolor = ROOT.kYellow  , fillstyle = 1001 ) ,
    ) 

Data_options = ROOT.RooFit.MarkerStyle ( 20   ) ,
               ROOT.RooFit.DrawOption  ( "zp" ) 

```


