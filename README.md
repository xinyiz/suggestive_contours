;:q
jku2015 Fall CS348A Homework4
================
Last modified:
October 2015  
Minhyuk Sung  
mhsung@cs.stanford.edu  


## Problem description
In this homework, you will
1) decimate meshes to a smaller number of triangles,
2) render suggestive contours to help illustrate their shape, and
3) implement an additional feature of your choosing.

There are several slots in need of your code, which are marked clearly and explained in the homework document.

To run the viewer, type:

	./bin/viewer ./models/horse.off 
	
at the command line, where './models/horse.off' is the model you wish to use.


## Setup
### Disclaimer
**The setup process is the exactly same with the previous homework. You can reuse the library files you previously installed. Create `lib` directory and move all library files to there, or change the paths written on `config.cmake.in` file.**

You can work on either [corn](https://web.stanford.edu/group/farmshare/cgi-bin/wiki/index.php/Main_Page#corn_info) machines, or on your personal machines. If you use your own machines, c++ complier and [cmake](https://cmake.org/) must be installed. The code is tested on the corn machines, which use the following versions:

* g++ (ver. 4.8.2)
* cmake (ver. 2.8.12)

### Install Libraries
The following libraries are required:

* [Freeglut](http://freeglut.sourceforge.net/): For rendering. Tested with ver. 2.8.1.
* [OpenMesh](http://www.openmesh.org/): For general mesh processing. Tested with ver. 4.1.
* [Eigen](http://eigen.tuxfamily.org/): For numerical computation using vectors/matrices. Tested with ver. 3.2.6.

GLUT is already installed on *corn* machines and Mac OS X system, so you can skip installing *Freeglut* if you use them.

For *Ubuntu* users, we recommend to install *Freeglut* using *apt-get*:

	sudo apt-get install freeglut3-dev
	
*Windows* users need to download and install *Freeglut* manually from the link above.


For the rest of libraries, if you use *Unix (Linux or Mac OS X)*, run the following script file.

	./install_libraries.sh

This will create `lib` directory and install libraries in this directory.

If you use *Windows* or the script file above has some problems, you can also manually download and install the libraries. In this case, you may need to change the library directories written in `config.cmake.in` file.

Please contact to CA if you have any problem installing libraries.


### Compile code
If you use *Unix (Linux or Mac OS X)*, run the following script file.

	./install.sh

This will create `build` directory and compile code. The output binary file will be created in `bin` directory. After the initial installation, you don't need to run `install.sh` again, but just run `make` in `build` directory.

You can also compile code using `CMake` yourself. Please take a look at `CMakeLists.txt` and `install.sh` files for details.


## Submit your code
First create your own README file that gives the names of your group members and an index to your program files and functions. Please compress the entire directory **except library files** into one file, name it as `cs348a-hw4-(group_member_names)`, and send it to the following email address.

	mhsung@cs.stanford.edu

Please include [cs348a-hw4] in the subject of your email, and make a clear note about collaboration in README.

For any question, please contact to CA.
