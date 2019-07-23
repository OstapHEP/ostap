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
