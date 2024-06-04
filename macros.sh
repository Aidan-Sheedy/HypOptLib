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
    echo -e '\t'build [clean] [install \<--PySysInstall\> \<--InstallPath=path\>]
    echo -e '\t\t'auto-builds the HypOptLib library.
    echo -e '\t\t'options:
    echo -e '\t\t' - clean: optional flag to delete make and cmake cache and files
    echo -e '\t\t' - install: optional flag to install the library on the system.
    echo -e '\t\t\t' --PySysInstall additional flag to use the system python dist-packages directory.
    echo -e '\t\t\t' --InstallPath=path additional flag to use the provided path as install sub-directory.
    echo -e '\t'setup [--mpi=\<option\>]
    echo -e '\t\t'Sets up the build environment. This can be done manually if this fails. See HypOptLib documentation for more details.
    echo -e '\t\t\t'1. Installs the latest PETSc version if not present
    echo -e '\t\t\t'2. Installs the latest mpi version if not present
    echo -e '\t\t\t'3. Installs the latest pybind11 version if not present\n
    echo -e '\t\t'-mpi: specifies which version of mpi to use, either openmpi or mpich. eg -mpi=openmpi.
    echo -e 
}

install_dependencies () {
    os=ubuntu
    mpi=mpich
    petscdir=../
    petscCfgArgs='--download-hdf5'

    for var in "$@"
    do
        if [ "--mpi=openmpi" =  $var ]
        then
            mpi=openmpi
        elif [ "--mpi=mpich" =  $var ]
        then
            mpi=mpich
        fi
    done

    if [ "openmpi" =  $mpi ]
    then
        petscCfgArgs+=' --download-openmpi'
    elif [ "mpich" =  $mpi ]
    then
        petscCfgArgs+=' --download-mpich'
    fi

    # Install HypOptLib dependencies
    sudo apt update

    sudo apt install cmake
    sudo apt install make

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
    ./configure $petscCfgArgs -j
    make all check -j

    sudo cp arch-linux-c-debug/bin/mpicc /usr/local/bin
    sudo cp arch-linux-c-debug/bin/mpicxx /usr/local/bin
    sudo cp arch-linux-c-debug/bin/mpiexec /usr/local/bin
    sudo cp arch-linux-c-debug/bin/mpirun /usr/local/bin
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

    makeArgs=""
    cmakeArgs=""
    scriptDir="$PWD"

    for var in "$@"
    do
        if [ "install" = $var ]
        then
            makeArgs="$makeArgs install"
        elif [ "clean" = $var ]
        then
            clean=True
        elif [ "--PySysInstall" = $var ]
        then
            if [ $customInstall ]
            then
                echo "Can't have both custom and python installs!"
                exit 0
            fi
            pyInstall=True
            cmakeArgs="$cmakeArgs -DInstallPythonSysPath=ON"
        elif [[ $var == *"--InstallPath="* ]]
        then
            if [ $pyInstall ]
            then
                echo "Can't have both custom and python installs!"
                exit 0
            fi
            customInstall=True
            cmakeArgs="$cmakeArgs -DCustomInstallPath=${var#*=}"
        elif [ "--DefaultPetscLocation" = $var ]
        then
            cmakeArgs="$cmakeArgs -DPETSC_DIR=$script_dir/../petsc -DPETSC_ARCH=arch-linux-c-debug"
        else
            echo "Invalid argument: $var"
            exit 0
        fi
    done

    if [ $clean ]
    then
        echo "-- Cleaning!!!"
        cd run
        rm ./HypOptLib.cpython*
        cd ../HypOptLib/
        rm build -r
    fi

    cd ./HypOptLib
    if ! test -f ./build/Makefile
    then
        mkdir build
    fi

    cd build
    echo "cmake $cmakeArgs .."
    cmake $cmakeArgs ..
    echo "make $makeArgs -j"
    make $makeArgs -j
}

if [ "-h" = "$1" ]
then
    print_help
elif [ "build_docs" = "$1" ]
then
    build_docs "$@"
elif [ "build" = "$1" ]
then
    shift
    build "$@"
elif [ "setup" = "$1" ]
then
    install_dependencies "$@"
else
    echo Invalid macro. Run \'./macros.sh \-h\' for more info.
fi
