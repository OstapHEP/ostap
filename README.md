[![build status](https://gitlab.cern.ch/amazurov/ostap/badges/master/build.svg)](https://gitlab.cern.ch/amazurov/ostap/commits/master)

# Ostap

[![Join the chat at https://gitter.im/OstapHEP/ostap](https://badges.gitter.im/OstapHEP/ostap.svg)](https://gitter.im/OstapHEP/ostap?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

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