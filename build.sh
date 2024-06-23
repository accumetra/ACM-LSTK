#!/bin/bash
#### Description : Setup and build ACM-LSTK and the LungNoduleSegmenter algorithm
#### Note : this script is run by AWS with elevated privileges when the aws instance stands up

#TODO: initial checks to ensure running as root (for chmod)

export DEV_DIR=/home/rick/Dev/accumetra

function configureInitialPrerequisites() {
  #
  # Update the instance and upgrade the installed libraries
  #
  apt-get update
  apt-get upgrade -y

  #
  # Install development libraries and build tools
  #
  apt-get -y install g++ build-essential cmake-curses-gui 

  # Required for MESA https://docs.mesa3d.org/install.html
  # Mesa has moved to the Meson build system so we need to install that 
  # https://mesonbuild.com/Quick-guide.html
  sudo apt-get install python3 python3-pip python3-setuptools python3-wheel ninja-build
  pip3 install mako
  pip3 install meson

  #
  # Install cmake
  #
  apt remove cmake -y
  apt purge --auto-remove cmake -y
}

function buildAndConfigureCmake() {

  mkdir -p $DEV_DIR/cmake
  cd $DEV_DIR/cmake
  wget https://github.com/Kitware/CMake/releases/download/v3.30.0-rc2/cmake-3.30.0-rc2-linux-x86_64.sh

  sh cmake-3.30.0-rc2-linux-x86_64.sh --skip-license
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
  apt-get install -y mesa-common-dev mesa-utils freeglut3-dev ninja-build libosmesa6-dev freeglut3 freeglut3-dev libglew-dev libglm-dev libgles2-mesa libgles2-mesa-dev
  # possibly add sudo apt install ocl-icd-opencl-dev for ITK GPU setting

  #Build VTK
  cd $DEV_DIR
  wget https://www.vtk.org/files/release/9.3/VTK-9.3.0.tar.gz
  tar -xvzf VTK-9.3.0.tar.gz 

  # Build VTK
  mkdir $DEV_DIR/VTK-build
  cd  $DEV_DIR/VTK-build

  cmake \
  -DBUILD_SHARED_LIBS=OFF \
  -DVTK_WRAP_PYTHON=OFF \
  -DVTK_ENABLE_VTKPYTHON=OFF \
  -DModule_vtkWrappingPythonCore=OFF \
  ../VTK-9.3.0 && make  -j 8
  # -DVTK_OPENGL_HAS_OSMESA=ON \
  # -DVTK_DEFAULT_RENDER_WINDOW_OFFSCREEN:BOOL=ON \
  # -DVTK_USE_X=OFF \
  # -DVTK_USE_GL2PS:BOOL=ON \
  # -DOSMESA_INCLUDE_DIR=$DEV_DIR/mesa-24.1.1/include \
  # -DOSMESA_LIBRARY=$DEV_DIR/mesa-24.1.1-install/lib/x86_64-linux-gnu/libOSMesa.so \
  # -DOPENGL_INCLUDE_DIR=$DEV_DIR/mesa-24.1.1/include \
  # -DOPENGL_gl_LIBRARY=$DEV_DIR/mesa-24.1.1-install/lib/x86_64-linux-gnu/libglapi.so \
  # -DOPENGL_glu_LIBRARY=$DEV_DIR/glu-9.0.3-install/lib/x86_64-linux-gnu/libGLU.so \

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
  git checkout v5.3.0 
  mkdir $DEV_DIR/ITK-build
  cd $DEV_DIR/ITK-build

  cmake ../ITK \
  -DBUILD_TESTING:BOOL=OFF \
  -DModule_ITKVtkGlue:BOOL=ON \
  -DModule_LesionSizingToolkit:BOOL=ON -DVTK_DIR:PATH="$DEV_DIR/VTK-build"  && make -j 8

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

  # Install the dependencies for Mesa
  sudo apt-get build-dep mesa -y
  sudo apt install cargo -y
  cargo install --force cbindgen
  export PATH=$PATH:$HOME/.cargo/bin

  cd $DEV_DIR
  wget https://archive.mesa3d.org/mesa-24.1.1.tar.xz
  tar xf mesa-24.1.1.tar.xz
  mkdir $DEV_DIR/mesa-24.1.1-install
  export MESA_INSTALL=$DEV_DIR/mesa-24.1.1-install

  cd $DEV_DIR/mesa-24.1.1

  # Build Mesa and install
  meson setup "../mesa-24.1.1-build/" \
    -Dprefix=$MESA_INSTALL \
    -Dgallium-drivers=swrast \
    -Dvulkan-drivers=intel \
    -Dbuildtype=release \
    -Dllvm=enabled \
    -Dosmesa=true \
    -Dvulkan-drivers=[] \
    -Degl=false \

  
  meson install -C "../mesa-24.1.1-build/"

  # This is to allow Mesa to run multi-threaded

  export GALLIUM_DRIVER=llvmpipe
  export LP_NUM_THREADS=4
}


function buildAndConfigureGlu() {
  ##### Make GLU
  cd $DEV_DIR
  wget "https://archive.mesa3d.org/glu/glu-9.0.3.tar.xz"
  tar xf glu-9.0.3.tar.xz
  mkdir $DEV_DIR/glu-9.0.3-build
  cd $DEV_DIR/glu-9.0.3

  export GLU_INSTALL=$DEV_DIR/glu-9.0.3-install
  mkdir $GLU_INSTALL

  meson setup "../glu-9.0.3-build/" \
    -Dprefix=$GLU_INSTALL \
    -Dbuildtype=release

  meson install -C "../glu-9.0.3-build/"
  
}

function finalizeNodeConfiguration() {
  # When it's all done, make sure to change the permissions
  chown -R ubuntu.ubuntu $DEV_DIR
}




# configureInitialPrerequisites
# buildAndConfigureCmake
# installAndConfigureJava

# buildAndInstallMesa
# buildAndConfigureGlu
# configureAndBuildVtk
# configureAndBuildItk

buildlungNoduleSegmentationAlgo

# Done 
echo "The ACM-LSTK build is complete."
