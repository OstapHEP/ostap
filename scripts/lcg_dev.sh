#---Get the location this script (thisdir)
thisdir=$(cd "$1"; pwd)
# Note: readlink -f is not working on OSX
if [ "LINUX" = "OSX" ]; then
  thisdir=$(python -c "import os,sys; print(os.path.realpath(os.path.expanduser(sys.argv[1])))" ${thisdir})
else
  thisdir=$(readlink -f ${thisdir})
fi



#  First the compiler
if [ "$COMPILER" != "native" ] && [ -e /cvmfs/sft-nightlies.cern.ch/lcg/contrib/gcc/9.2.0/x86_64-centos7/setup.sh ]; then
    source /cvmfs/sft-nightlies.cern.ch/lcg/contrib/gcc/9.2.0/x86_64-centos7/setup.sh
fi

#  LCG version
LCG_VERSION=$2; export LCG_VERSION

#  then the rest...
if [ -z "${PATH}" ]; then
    PATH=${thisdir}/bin; export PATH
else
    PATH=${thisdir}/bin:$PATH; export PATH
fi
if [ -d ${thisdir}/scripts ]; then
    PATH=${thisdir}/scripts:$PATH; export PATH
fi

if [ -z "${LD_LIBRARY_PATH}" ]; then
    LD_LIBRARY_PATH=${thisdir}/lib; export LD_LIBRARY_PATH
else
    LD_LIBRARY_PATH=${thisdir}/lib:$LD_LIBRARY_PATH; export LD_LIBRARY_PATH
fi
if [ -d ${thisdir}/lib64 ]; then
    LD_LIBRARY_PATH=${thisdir}/lib64:$LD_LIBRARY_PATH; export LD_LIBRARY_PATH
fi

if [ -x ${thisdir}/bin/python ]; then
  PYTHON_VERSION=`expr $(readlink ${thisdir}/bin/python) : '.*/Python/\([0-9].[0-9]\).*'`
else
  PYTHON_VERSION=`python -c "import sys; print('{0}.{1}'.format(*sys.version_info))"`
fi
PY_PATHS=${thisdir}/lib/python$PYTHON_VERSION/site-packages

if [ -z "${PYTHONPATH}" ]; then
    PYTHONPATH=${thisdir}/lib:$PY_PATHS; export PYTHONPATH
else
    PYTHONPATH=${thisdir}/lib:$PY_PATHS:$PYTHONPATH; export PYTHONPATH
fi
if [ -d ${thisdir}/python ]; then
    PYTHONPATH=${thisdir}/python:$PYTHONPATH; export PYTHONPATH
fi

if [ -z "${MANPATH}" ]; then
    MANPATH=${thisdir}/man:${thisdir}/share/man; export MANPATH
else
    MANPATH=${thisdir}/man:${thisdir}/share/man:$MANPATH; export MANPATH
fi
if [ -z "${CMAKE_PREFIX_PATH}" ]; then
    CMAKE_PREFIX_PATH=${thisdir}; export CMAKE_PREFIX_PATH
else
    CMAKE_PREFIX_PATH=${thisdir}:$CMAKE_PREFIX_PATH; export CMAKE_PREFIX_PATH
fi
if [ -z "${CPLUS_INCLUDE_PATH}" ]; then
    CPLUS_INCLUDE_PATH=${thisdir}/include; export CPLUS_INCLUDE_PATH
else
    CPLUS_INCLUDE_PATH=${thisdir}/include:$CPLUS_INCLUDE_PATH; export CPLUS_INCLUDE_PATH
fi
if [ -z "${C_INCLUDE_PATH}" ]; then
    C_INCLUDE_PATH=${thisdir}/include; export C_INCLUDE_PATH
else
    C_INCLUDE_PATH=${thisdir}/include:$C_INCLUDE_PATH; export C_INCLUDE_PATH
fi
#if [ -z "${LIBRARY_PATH}" ]; then
#    LIBRARY_PATH=${thisdir}/lib; export LIBRARY_PATH
#else
#    LIBRARY_PATH=${thisdir}/lib:$LIBRARY_PATH; export LIBRARY_PATH
#fi
#if [ -d ${thisdir}/lib64 ]; then
#    export LIBRARY_PATH=${thisdir}/lib64:$LIBRARY_PATH
#fi

#---check for compiler variables
if [ -z "${CXX}" ]; then
    export FC=`command -v gfortran`
    export CC=`command -v gcc`
    export CXX=`command -v g++`
fi

#---Figure out the CMAKE_CXX_STANDARD (using Vc as a victim)
if [ -f $thisdir/include/Vc/Vc ]; then
    vc_home=$(dirname $(dirname $(dirname $(readlink $thisdir/include/Vc/Vc))))
    std=$(cat $vc_home/logs/Vc*configure.cmake | egrep -o "CMAKE_CXX_STANDARD=[0-9]+" | egrep -o "[0-9]+")
    export CMAKE_CXX_STANDARD=$std
fi

#---then ROOT
if [ -x $thisdir/bin/root ]; then
    if [ -x $thisdir/bin/python ]; then
        PYTHON_INCLUDE_PATH=$(dirname $(dirname $(readlink $thisdir/bin/python)))/include/$(\ls $(dirname $(dirname $(readlink $thisdir/bin/python)))/include)
    fi
    ROOTSYS=$(dirname $(dirname $(readlink $thisdir/bin/root))); export ROOTSYS
    if [ -z "${ROOT_INCLUDE_PATH}" ]; then
        ROOT_INCLUDE_PATH=${thisdir}/include:$PYTHON_INCLUDE_PATH; export ROOT_INCLUDE_PATH
    else
        ROOT_INCLUDE_PATH=${thisdir}/include:$PYTHON_INCLUDE_PATH:$ROOT_INCLUDE_PATH; export ROOT_INCLUDE_PATH
    fi
    if [ -d $thisdir/targets/x86_64-linux/include ]; then
        ROOT_INCLUDE_PATH=${thisdir}/targets/x86_64-linux/include:$ROOT_INCLUDE_PATH; export ROOT_INCLUDE_PATH
    fi
    if [ -z "${JUPYTER_PATH}" ]; then
        JUPYTER_PATH=${thisdir}/etc/notebook; export JUPYTER_PATH
    else
        JUPYTER_PATH=${thisdir}/etc/notebook:$JUPYTER_PATH; export JUPYTER_PATH
    fi
    export CPPYY_BACKEND_LIBRARY=$ROOTSYS/lib/libcppyy_backend${PYTHON_VERSION/./_}
    export CLING_STANDARD_PCH=none
fi

#---then Gaudi
if [ -x $thisdir/include/Gaudi ]; then
    jsoninc=$(dirname $(dirname $(readlink ${thisdir}/include/nlohmann/json.hpp))) # see https://github.com/root-project/root/issues/7950
    ROOT_INCLUDE_PATH=${jsoninc}:${thisdir}/src/cpp:$ROOT_INCLUDE_PATH; export ROOT_INCLUDE_PATH
fi
if [ -x $thisdir/scripts/gaudirun.py ]; then
    Gaudi_DIR=$(dirname $(dirname $(readlink $thisdir/scripts/gaudirun.py)));
    export CMAKE_PREFIX_PATH=$Gaudi_DIR:$CMAKE_PREFIX_PATH
fi

#---then PYTHON
if [ -x $thisdir/bin/python ]; then
    PYTHONHOME=$(dirname $(dirname $(readlink $thisdir/bin/python))); export PYTHONHOME
elif [ -x $thisdir/bin/python3 ]; then
    PYTHONHOME=$(dirname $(dirname $(readlink $thisdir/bin/python3))); export PYTHONHOME
fi


if [ -f $thisdir/lib/libQt5Gui.so ]; then
    export QT_PLUGIN_PATH=$(dirname $(dirname $(readlink $thisdir/lib/libQt5Gui.so )))/plugins
    export QT_XKB_CONFIG_ROOT=/usr/share/X11/xkb
fi

if [ -f $thisdir/etc/fonts/fonts.conf ]; then
    export FONTCONFIG_PATH=$thisdir/etc/fonts
fi

#---then PKG_CONFIG_PATH
if [ -z "$PKG_CONFIG_PATH" ]; then
    export PKG_CONFIG_PATH="$thisdir/lib/pkgconfig"
else
    export PKG_CONFIG_PATH="$thisdir/lib/pkgconfig:$PKG_CONFIG_PATH"
fi
if [ -d ${thisdir}/lib64/pkgconfig ]; then
    PKG_CONFIG_PATH=${thisdir}/lib64/pkgconfig:$PKG_CONFIG_PATH; export PKG_CONFIG_PATH
fi



