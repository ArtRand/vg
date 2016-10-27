# added for mac os builds, taken from README
export PATH="/usr/local/opt/coreutils/libexec/gnubin:/usr/local/bin:$PATH"

# Use glibtool/ize
export LIBTOOL=glibtool
export LIBTOOLIZE=glibtoolize

export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH;
export LIBRARY_PATH=$LD_LIBRARY_PATH;

export LIBRARY_PATH=`pwd`/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=`pwd`/lib:$LD_LIBRARY_PATH
export LD_INCLUDE_PATH=`pwd`/include:$LD_INCLUDE_PATH
export C_INCLUDE_PATH=`pwd`/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=`pwd`/include:$CPLUS_INCLUDE_PATH
export INCLUDE_PATH=`pwd`/include:$INCLUDE_PATH
export PATH=`pwd`/bin:`pwd`/scripts:$PATH
export CC=$(which gcc)
export CXX=$(which g++)

#
#  disable until file arguments work as in normal bash :(
#
# add bash autocompletion
#if test -n "$BASH_VERSION"
#then
#
#	 . ./autocomp.bash
#fi
