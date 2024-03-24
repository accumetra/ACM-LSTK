#!/bin/bash
#### Description : Setup and build ACM-LSTK and the LungNoduleSegmenter algorithm
#### Note : this script is run by AWS with elevated privileges when the aws instance stands up

#TODO: initial checks to ensure running as root (for chmod)

export DEV_DIR=/home/ubuntu/Dev

function configureInitialPrerequisites() {
  #
  # Update the instance and upgrade the installed libraries
  #
  apt-get update
  apt-get upgrade -y

  #
  # Install development libraries and build tools
  #
  apt-get -y install g++ build-essential cmake-curses-gui python-minimal

  #
  # Install cmake
  #
  apt remove cmake -y
  apt purge --auto-remove cmake -y
}

function buildAndConfigureCmake() {

  mkdir -p $DEV_DIR/cmake
  cd $DEV_DIR/cmake
  wget https://cmake.org/files/v3.12/cmake-3.12.0-Linux-x86_64.sh 

  sh cmake-3.12.0-Linux-x86_64.sh --skip-license
  ln -s $DEV_DIR/cmake/bin/cmake /usr/local/bin/cmake
  ln -s $DEV_DIR/cmake/bin/ccmake /usr/local/bin/ccmake
  ln -s $DEV_DIR/cmake/bin/cmake-gui /usr/local/bin/cmake-gui
}

function installAndConfigureJava() {
  # Install OpenJDK 8.0
  add-apt-repository ppa:openjdk-r/ppa -y # Only ubuntu 17.4 or earlier
  apt-get update
  apt-get -y install openjdk-8-jdk 

  # Add JAVA_HOME to /etc/environment
  echo "JAVA_HOME=/usr/lib/jvm/java-8-openjdk-amd64/jre" >> /etc/environment
}

configureAndBuildVtk() {
  # Now we will prepare to build VTK
  # First we need to install some libraries
  #
  # Install graphics libraries for visualization tasks
  #
  apt-get install -y mesa-common-dev freeglut3 freeglut3-dev libglew1.5-dev libglm-dev libgles2-mesa libgles2-mesa-dev 
  # possibly add sudo apt install ocl-icd-opencl-dev for ITK GPU setting

  #Build VTK
  cd $DEV_DIR
  git clone https://gitlab.kitware.com/vtk/vtk.git VTK
  cd VTK
  git checkout v8.1.0

  # Build VTK
  mkdir $DEV_DIR/VTK-build
  cd  $DEV_DIR/VTK-build
  #cmake ../VTK -DBUILD_SHARED_LIBS:BOOL=OFF -DBUILD_TESTING:BOOL=OFF -DCMAKE_BUILD_TYPE:STRING=Release && make -j 8
  cmake \
  -DBUILD_SHARED_LIBS=OFF \
  -DVTK_WRAP_PYTHON=OFF \
  -DVTK_ENABLE_VTKPYTHON=OFF \
  -DModule_vtkWrappingPythonCore=OFF \
  -DVTK_OPENGL_HAS_OSMESA=ON \
  -DVTK_DEFAULT_RENDER_WINDOW_OFFSCREEN:BOOL=ON \
  -DVTK_USE_X=OFF \
  -DVTK_USE_GL2PS:BOOL=ON \
  -DOSMESA_INCLUDE_DIR=$DEV_DIR/mesa-17.2.8/include \
  -DOSMESA_LIBRARY=$DEV_DIR/mesa-17.2.8/lib/gallium/libOSMesa.so \
  -DOPENGL_INCLUDE_DIR=$DEV_DIR/mesa-17.2.8/include \
  -DOPENGL_gl_LIBRARY=$DEV_DIR/mesa-17.2.8/lib/libglapi.so \
  -DOPENGL_glu_LIBRARY=$DEV_DIR/glu-9.0.0/.libs/libGLU.so \
  ../VTK && make  -j 8

  if [ $? -eq 0 ]
  then
    echo "Successfully configured and built VTK"
  else
    echo "Failed to configure & build VTK" >&2
    exit 2
  fi

}

function configureAndBuildItk() {
  # Build ITK
  cd $DEV_DIR
  #git clone git://itk.org/ITK.git
  git clone https://github.com/InsightSoftwareConsortium/ITK
  cd $DEV_DIR/ITK
  git checkout v4.13.0
  mkdir $DEV_DIR/ITK-build
  cd $DEV_DIR/ITK-build
  cmake ../ITK -DBUILD_TESTING:BOOL=OFF -DModule_ITKVtkGlue:BOOL=ON -DModule_LesionSizingToolkit:BOOL=ON -DVTK_DIR:PATH=$DEV_DIR/VTK-build  && make -j 8
  if [ $? -eq 0 ]
  then
    echo "Successfully configured and built ITK"
  else
    echo "Failed to configure & build ITK" >&2
    exit 3
  fi
}

function buildlungNoduleSegmentationAlgo() {
  # Now lets build the LungNoduleSegmenter application (relies on ITK)
  cd $DEV_DIR/ACM-LSTK/src
  mkdir $DEV_DIR/ACM-LSTK/src/LungNoduleSegmenter-build
  cd $DEV_DIR/ACM-LSTK/src/LungNoduleSegmenter-build
  cmake ../LungNoduleSegmenter/ -DCMAKE_BUILD_TYPE:STRING=Release -DITK_DIR:PATH=$DEV_DIR/ITK-build  && make
  if [ $? -eq 0 ]
  then
    echo "Successfully configured and built LungSegmenter algorithm"
  else
    echo "Failed to configure & build LungNoduleSegmenter algorithm" >&2
    exit 5
  fi
}

function buildAndInstallMesa() {
  sudo add-apt-repository ppa:xorg-edgers/ppa -y
  sudo apt-get update && sudo apt-get dist-upgrade -y

  # Need the following for autoconf to work below
  sudo apt-get install -y autoconf autogen automake pkg-config libgtk-3-dev libtool llvm-dev

  cd $DEV_DIR
  wget https://archive.mesa3d.org/older-versions/17.x/mesa-17.2.8.tar.xz
  tar xf mesa-17.2.8.tar.xz
  cd $DEV_DIR/mesa-17.2.8
  export MESA_INSTALL=$(pwd)


  ./configure \
    --disable-xvmc \
    --disable-glx \
    --enable-dri \
    --with-dri-drivers= \
    --with-gallium-drivers=swrast \
    --enable-texture-float \
    --disable-egl \
    --with-platforms= \
    --enable-gallium-osmesa \
    --enable-llvm=yes



  make -j 8
  make install


  # This is to allow Mesa to run multi-threaded

  export GALLIUM_DRIVER=llvmpipe
  export LP_NUM_THREADS=4
}


function buildAndConfigureGlu() {
  ##### Make GLU
  cd $DEV_DIR
  wget ftp://ftp.freedesktop.org/pub/mesa/glu/glu-9.0.0.tar.bz2
  tar xvf glu-9.0.0.tar.bz2
  cd $DEV_DIR/glu-9.0.0
  ./configure \
    --prefix=$DEV_DIR/mesa-17.2.8 \
    --enable-osmesa \
    PKG_CONFIG_PATH=$DEV_DIR/mesa-17.2.8/lib/pkgconfig

  make -j 8
  make install
}

function finalizeNodeConfiguration() {
  # When it's all done, make sure to change the permissions
  chown -R ubuntu.ubuntu $DEV_DIR
}




#configureInitialPrerequisites
#buildAndConfigureCmake
#installAndConfigureJava

#buildAndInstallMesa
#buildAndConfigureGlu
#configureAndBuildVtk
#configureAndBuildItk

buildlungNoduleSegmentationAlgo

# Done 
echo "The ACM-LSTK build is complete."
