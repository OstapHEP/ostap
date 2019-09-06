# Paralell Pathos 

* [ostap.mp_pathos](PARALELL.md)

*Parallelisation* functionality based on [`pathos`](https://github.com/uqfoundation/pathos)-project 

It works both in multiprocessor and cluster regimes.
For cluster regime the password-less ssh-connection is requred. 


```python
servers = 'lxplus701.cern.ch', 'lxplus702.cern.ch', 'lxplus703.cern.ch', 
# wm    = WorkManager ( servers = servers )             ## with local CPU /"autodetect" 
# wm    = WorkManager ()                                ## only local CPU  (no remotes)
wm      = WorkManager ( ncpus = 0 , servers = servers ) ## without local CPU  
```

## Preconfigured  list of servers 

If `servers` is set to be 'auto' or 'congif',  the list of servers 
is picked up form the configrutaion file and/or `OSTAP_PPSERVERS` environment variable

### List of servers from configration file

If configuration   file constains a section `Parallel:<domain>`, where `<domain>` is the
domain of local host, and this section contains the property `ppservers`, ti will be used as 
default list of servers for the given domain
e.g.
```
[Parallel:CERN.CH]
ppservers = lxplus707 , lxplus706 , lxplus700  
```
If no such section/setting exists,  the property `ppserver` is looked into `Parallel` section:
```
[Parallel]
ppservers = lxplus707.cern.ch , lxplus706.cern.ch , lxplus700.cern.ch  
```

### List of servers from environment variable 

The environment variable `OSTAP_PPSERVERS` can be used to provide the list of remote servers:
```shell
export OSTAP_PPSERVERS="lxplus701,lxplus702,lxplus703" 
ostap 
```
Or: 
```shell
OSTAP_PPSERVERS="lxplus701,lxplus702,lxplus703" ostap 
```


The remote `ppserver` processes will be automatically launched using `pathos` 
ssh-tunneling  functionality

## Configuration of remote  servers

If remote hosts provides the proper  environment, and, in paricular `ppserver` is in  the path, no actionnns is required. Otherwise three possible (non-exclusive)  actions  can be specified to   setup the proper environment 
  
```python
wm      = WorkManager ( ... , 
    profile     = 'remote_script.sh' , 
    environment = """#!/bin/bash
        source /cvmfs/sft.cern.ch/lcg/views/LCG_95/${CMTCONFIG}/setup.sh
        source ${HOME}/cmtuser/ostap/build/INSTALL/${CMTCONFIG}/thisostap.sh
       """ , 
    script = 'local_sctipt.sh'           
   )
```

### `profile`

If `profile` is provdied, it should be the name of the existing script on remote host. The script wil lbe sourced to   provide the correct environment 

### `environment`

If `environment` is provided, it should be the (multiline) string, that will be dumped into  temporary file, transferred to the remote host and sources to provide the correct enbvironment 
 
### `script`

If `script` is provided, it should be the name of    exisitng script on local host. It will transferred to the remote host and sources to provide the correct enbvironment 
 


All three properties for each remote host can be specified via the configrutaion file into three sections
 - `Parallel:<remote host>` , e.g. `Parallel:lxplus701.cern.ch`
 - `Parallel:<domain>`, e.g.  `Parallel:CERN.CH` 
 - `Parallel`
e.g. 
```
[Paralell:lxplus701.cern.ch]
 script = my_local_script_for_701.sh

[Paralell:lxplus707.cern.ch]
 environment = #!/bin/bash
     source /cvmfs/sft.cern.ch/lcg/views/LCG_96/${CMTCONFIG}/setup.sh
     export PATH=$(IFS=':';t=($PATH);n=${#t[*]};a=();for ((i=0;i<n;i++)); do p="${t[i]%%*LCG_95*}"; [ "${p}" ] && a[i]="${p}"; done;echo "${a[*]}");
     export PYTHONPATH=$(IFS=':';t=($PYTHONPATH);n=${#t[*]};a=();for ((i=0;i<n;i++)); do p="${t[i]%%*LCG_95*}"; [ "${p}" ] && a[i]="${p}"; done;echo "${a[*]}");
     source ${HOME}/cmtuser/ostap/build/INSTALL/${CMTCONFIG}/thisostap.sh

[Paralell:CERN.CH]
 profile = ostap_profile.sh

[Parallel] 
 script = ostap_profile.sh

```

## Examples 
```python

def my_fun ( x ) : return x**2 

wm = WorkManager ( servers = "config" ) 

items = range ( 100 )

## get (unordered) list of squares for numbers  from  11  to 100: 
result1 =  wm.process ( my_fun , items , merger = TaskMerger ( lambda  a,b : a+[b] , init = [] ) )

## get sum of squares for numbers  from  11  to 100: 
result2 =  wm.process ( my_fun , items , merger = TaskMerger () )    
```




  





