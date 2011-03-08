#!/bin/bash

make -j3 1>/dev/null
sudo make uninstall 1>/dev/null
sudo make install 1>/dev/null
if [ `uname` == "Darwin" ]
then
	sudo dsymutil latan/.libs/liblatan.0.dylib --out=/usr/local/lib/liblatan.0.dylib.dSYM
fi

