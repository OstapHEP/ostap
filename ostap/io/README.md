# Input/Output 

* [ostap.io](README.md)


Several useful classes and functions for *input/output* operations 

- Very useful class `ZipShelf`, essentially the *zipped* version of the standard `shelve` module
   - it allows to store all  *pickable* objects, in particualr almost all `ROOT`-classes 
   - it stores them in *zipped* format, therefore the resulting database is rather compact 
   - it have very simple and intutitive `dict`-like/`shelve` interface and easy-to-use 
```
import zipshelve  ## import the ZipShelve module 
with zipshelve.open ('a_db', 'n')  as db :  ## create new DB
   abcde = ...
   db['some_key'] =  abcde              ## add information to DB 
   ...
   abcd = db['some_key']                ## get information from DB 
```
- `SQLiteShelf` - very similar to `ZipShelve` but uses SQLite as storage backend 
- `RootShelf` - very similar to `ZipShelve` but uses `ROOT.TFile` as storage backend 
- `RootOnlyShelf` - very similar to the previos one, also uses `ROOT.TFile` as storage backend. but allows to store only objects storable in `ROOT.TFile` 
 
Also  it provided a useful pythoniized decorations for `ROOT.TFile`/`ROOT.TDirectory`:
```
>>> rfile = ...

>>> obj   = rfile['A/B/C/myobj']     ## READ  object form the file/directory
>>> rfile['A/B/C/myobj2'] = object2  ## WRITE object to the file/directory 

>>> obj  = rfile.A.B.C.myobj              ## another way to access to the object
>>> obj  = rfile.get ( 'A/B/C/q' , None ) ## one more way to get object 

>>> for obj in rfile : print obj            ## loop over all objects in file
>>> for key,obj in rfile.iteritems() : print key, obj             ## another loop
>>> for key,obj in rfile.iteritems( ROOT.TH1 ) : print key, obj   ## advanced loop
>>> for k in rfile.keys()     : print k   ## get all keys and loop over them 
>>> for k in rfile.iterkeys() : print k   ## loop over all keys in the file

>>> del  rfile['A/B']                       ## delete the object from the file
>>> rfile.rm ( 'A/B' )                      ## delete the object from the file

>>> if 'A/MyHisto' in rfile          : print 'OK!' ## check presence of the key
>>> if rfile.has_key ( 'A/MyHisto' ) : print 'OK!' ## check presence of the key

>>> with ROOT.TFile('aa.root') as rfile : rfile.ls() ## context manager protocol
```