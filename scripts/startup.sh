#!/bin/bash
#!/bin/bash

error() {
    echo "************************"
    echo "   $1 not found         "
    echo "   Installing now       "
    echo "************************"
}

title() {
    echo "************************"
    echo "   checking for $1      "
    echo "************************"
}

install_command() {
  title $1
    if ! command -v $1 &> /dev/null; then
        error $1
        sudo apt install $1 -y
    fi  
}

install_pkg() {
  title $1
    if ! dpkg -l | grep -q $1; then
        error $1
        sudo apt install $1 -y
    fi  
}


echo "Downloading required repos"

#sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y
#sudo apt-get update -y

install_command gcc

install_pkg g++-10

install_pkg make

echo "Setting env vars"
source tbb2020/bin/tbbvars.sh intel64
echo "Starting build, this is expected to take several minutes...."
sleep 2
make -j CC=g++-10

echo "Install complete"