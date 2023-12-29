#!/bin/bash

if [ 0 -eq $# ]
then
    echo Invalid argument number, must provide one macro to use. Run \'./macros.sh \-h\' for more info.
    exit 0
fi

# Environment Variables (Set these up before using!)
# doxygen_exe=/mnt/g/Program\ Files/doxygen/bin/doxygen.exe

print_help () {
    echo -e "macros.sh [-h] macro"
    echo -e "Supported macros:"
    echo -e '\t'build_docs
    echo -e '\t\t'generates documentation for the library in the \'docs\' folder.
    echo -e '\t'build [options] [install]
    echo -e '\t\t'auto-builds the HypOptLib library.
    echo -e '\t\t'options:
    echo -e '\t\t' - [none] defaults to all
    echo -e '\t\t' - all: will build cmake files if none present, otherwise just run make file.
    echo -e '\t\t' - clean: deletes all build files and objects and rebuilds fresh.
    echo -e '\t\t' install: installs the HypOptLib library as a python module. Can then be accessed system-wide by just importing HypoptLib.
    echo -e '\t'setup -- TODO! Not yet implemented.
    echo -e '\t\t'Sets up the build environment. This can be done manually if this fails. See HypOptLib documentation for more details.
    echo -e '\t\t\t'1. Installs the latest PETSc version if not present
    echo -e '\t\t\t'2. Installs the latest mpi version if not present
    echo -e '\t\t\t'3. Installs the latest pybind11 version if not present
    echo -e 
}

install_dependencies () {
    os=ubuntu
    mpi=mpich
    petscdir=../

    for var in "$@"
    do
        if [ "-mpi=openmpi" =  $var ]
        then
            mpi=openmpi
        elif [ "-mpi=mpich" =  $var ]
        then
            mpi=mpich
        fi
    done

    echo $mpi

    exit 0

    # Install HypOptLib dependencies
    sudo apt update

    sudo apt install cmake
    sudo apt install make

    if [ "mpich" = $mpi ]
    then
        sudo apt install mpich
    elif [ "openmpi" = $mpi ]
    then
        sudo apt install openmpi-bin
    fi

    sudo apt install python3-pip
    pip3 install "pybind11[global]"
    sudo apt install libhdf5-serial-dev

    # Install PetSc Dependencies
    sudo apt-get install libblas-dev liblapack-dev

    # Install PetSc
    cd $petscdir
    git clone -b release https://gitlab.com/petsc/petsc.git petsc
    cd petsc
    git pull
    ./configure --download-hdf5
    make all check
}

build_docs () {
    if [ 1 -lt $# ]
    then
        echo Invalid argument number for macro \'build_docs\'. Run \'./macros.sh \-h\' for more info.
        exit 0
    fi

    cd ./docs
    rm -r build
    rm -r doxyxml
    doxygen Doxyfile
    sphinx-build -M html source build
}

build () {
    if [ 3 -lt $# ]
    then
        echo Invalid argument number for macro \'build\'. Run \'./macros.sh \-h\' for more info.
        exit 0
    elif [ 1 -eq $# ]
    then
        option="all"
    else
        option="$2"
    fi

    if [ 3 -eq $# ]
    then
        if [ "install" = $3 ]
        then
            install=True
        else
            echo Invalid argument "$3" for macro \'build "$2"\'. Run \'./macros.sh \-h\' for more info.
        fi
    else
        install=False
    fi

    if [ "clean" = "$option" ]
    then
        cd run
        rm ./HypOptLib.cpython*
        cd ../HypOptLib/
        rm build -r
        mkdir build
        cd build
        cmake ..
        if [ "True" = "$install" ]
        then
            make install
        else
            make
        fi
    elif [ "all" = "$option" ]
    then
        cd ./HypOptLib
        if ! test -f ./build/Makefile
        then
            rm build -r
            mkdir build
            cd build
            cmake ..
        else
            cd build
        fi
        if [ "True" = "$install" ]
        then
            make install
        else
            make
        fi
    else
        echo Invalid argument "$2" for macro \'build\'. Run \'./macros.sh \-h\' for more info.
    fi
}

if [ "-h" = "$1" ]
then
    print_help
elif [ "build_docs" = "$1" ]
then
    build_docs "$@"
elif [ "build" = "$1" ]
then
    build "$@"
elif [ "setup" = "$1" ]
then
    install_dependencies "$@"
else
    echo Invalid macro. Run \'./macros.sh \-h\' for more info.
fi
