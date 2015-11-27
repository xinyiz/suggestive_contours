echo "*** 2015 Fall CS348A Programming Assignments ***"
echo ""

current_dir=$(pwd)

if [ ! -d lib ]; then
  mkdir lib
fi
cd lib

# OpenMesh
echo "Install OpenMesh-4.1..."
rm -rf OpenMesh
wget http://www.openmesh.org/media/Releases/4.1/OpenMesh-4.1.tar.gz --no-check-certificate
tar xvf OpenMesh-4.1.tar.gz
rm OpenMesh-4.1.tar.gz -rf
mv OpenMesh-* OpenMesh
cd OpenMesh
mkdir build && cd build
cmake ..
make
cd ../../
echo "Done."

# Eigen
echo "Install Eigen-3.2.6..."
rm -rf eigen-3.2.6
wget https://bitbucket.org/eigen/eigen/get/3.2.6.tar.gz --no-check-certificate
tar xvf 3.2.6.tar.gz
rm -rf 3.2.6.tar.gz
mv eigen-eigen-* Eigen
echo "Done."

cd $current_dir
