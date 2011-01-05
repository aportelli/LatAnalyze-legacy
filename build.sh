#!/bin/bash

make -j3 1>/dev/null
sudo make uninstall 1>/dev/null
sudo make install 1>/dev/null
