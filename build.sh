#!/bin/bash

make -j3 1>/dev/null
make uninstall 1>/dev/null
make install 1>/dev/null
if [ `uname` == "Darwin" ]
then
	dsymutil latan/.libs/liblatan.0.dylib --out=${HOME}/local/lib/liblatan.0.dylib.dSYM
fi

