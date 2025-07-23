#!/bin/bash
#
# ==== WARNING: NOT INTENDED TO BE RUN AS ROOT; OTHERWISE, WILL INSTALL OPENMM INTO YOUR SYSTEM, POSSIBLY OVERWRITING EXISTING FILES ====
#
# Here is what we're using :
uname -a 
echo
gcc -v
echo
cmake --version
echo
cat /proc/cpuinfo |head -n10
ncores=`cat /proc/cpuinfo|grep -i "cpu cores" | tail -n1 | awk '{print $NF}'`
echo "Will use $ncores cores"
ROOT=${PWD} ;# this directory
#
ommver="7.3.0" # other versions may require plugin updates; patches may not be needed
dynamo_opencl=1
OMM_HOME=${ROOT}/omm-${ommver}
PYTHON=`which python`
#
# 1 build openmm
#
if  [[ ! -d openmm ]]; then
 git clone https://github.com/pandegroup/openmm.git
 pushd openmm
 git checkout tags/$ommver
 swig_version=`swig -version | grep -i version | awk -F["\ ."] '{print $3}'`
 if [[ $swig_version -ge 4 ]]; then
  cp -f ../patch7.2i/generateWrappers.py wrappers ;# dirty but good enough-- update for swig 4
 fi
 rm -fr build
 mkdir -p build
 pushd build
 cmake ../ -DCMAKE_INSTALL_PREFIX=$OMM_HOME -DPYTHON_EXECUTABLE=$PYTHON
 make -j$ncores
 make install
 make PythonInstall # this will build & attempt to install into root filesystem, and hopefully fail
 popd
 popd
 pushd $OMM_HOME/include/swig
 patch < ${ROOT}/patch7.2i/OpenMMSwigHeaders.diff
 popd
fi
#
# 2 get plugin source
#
PLUGINDIR=${ROOT}/vo-projects
#PLUGINBRANCH=protforce # there are two branches, master and protforce ; protforce is often slightly faster for reasons that are not fully clear ...
PLUGINBRANCH=master
if  [[ ! -d "$PLUGINDIR" ]]; then
 git clone https://github.com/ovchinnv/vo-projects $PLUGINDIR
fi
pushd $PLUGINDIR
# make clean # optional rebuild
git checkout $PLUGINBRANCH
#
# 3 build fortran library :
#
unset MAKE_COMMAND # must compile in serial mode
# first, compile and run with the mindist method (only in watershell v.1):
sed -i 's/version=.*$/version=1/' watershell/version.inc # specify version in config file:
make --silent -j1 watershell # in case we are rebuilding only the watershell
make --silent -j1 plugin_master # parallel compile not supported
#
# 4 build openmm wrapper
#
pushd plugin_master
pushd openmm-dynamo
mkdir -p build
pushd build
if [[ $dynamo_opencl -eq 0 ]] ; then
  EXTRAFLAGS+="-DDYNAMO_BUILD_OPENCL_LIB=OFF"
fi
cmake ../ -DOPENMM_DIR=$OMM_HOME -DCMAKE_INSTALL_PREFIX=$OMM_HOME -DDYNAMO_INCLUDE_DIR=$PLUGINDIR/include \
 -DDYNAMO_LIBRARY_DIR=$PLUGINDIR/lib -DPYTHON_EXECUTABLE=$PYTHON $EXTRAFLAGS
make -j$ncores
make install
make PythonInstall # will build but fail to install as non-root
popd
popd
popd
popd
#fi
#exit
#run benchmarks:
./run-dhfr
./run-b12
#exit
# rebuild plugin to provide watershell v2 and rerun :
pushd $PLUGINDIR
# make clean # optional rebuild
sed -i 's/version=.*$/version=2/' watershell/version.inc # specify version in config file:
make --silent -j1 watershell
make --silent -j1 plugin_master # parallel compile not supported
popd
./run-dhfr # b12 test case does not support watershell 2.0
