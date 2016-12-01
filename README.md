[![build status](https://gitlab.cern.ch/amazurov/ostap/badges/master/build.svg)](https://gitlab.cern.ch/amazurov/ostap/commits/master)

# Ostap

## Setup

You need to provide ROOTSYS environment for building the library. Library should be built with the same compiler version as was used to build ROOT.

```
source /path/to/root/bin/thisroot.sh
```
or, if you are working at the LHCb environment
```
lb-run ROOT bash
```


## Example

Build and run python/example.py:
```
./scripts/bootstrap
```