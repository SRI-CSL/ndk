Installation on Linux and MacOSX
================================

NDK was originally developed in Python 2.7 under Ubuntu Linux and Mac
OS X.  Current and future development will rely on Python 3, but will
be backward compatible with Python 2.7.

NDK is available at the git repository https://github.com/SRI-CSL/ndk

The ndk can be installed by cd'ing into the ndk top-level directory
and doing the following:

    pip install .

NDK's setup.py file includes all of the needed dependencies, so they
should be autoloaded when you use pip.

here are the module dependencies:

	neo
	numpy
	scipy
	sqlite3
	sklearn
	pyOpenGL
	quantities
	matplotlib

	Theano (optional)

The 'mkspikes2' script uses Theano for fast convolution.  If you don't
have theano installed, 'mkspikes' will do the job more slowly.

Some features (those that depend on raw data, like 'vcoefs') will not
work unless you set your "NDKDATA" environment variable to point to a
colon-separated list of directories that contain .smr or .wav files,
e.g.

     export NDKDATA=/opt/spikes/FI/dist:/opt/spikes/new_data:


The ndk module has two submodules, ds, and es, that provide access to
raw data stores and event stores respectively (hence the names).

Matlab notes: with version R2014b and higher, python functionality is
easily accessed within Matlab.  The ndk subdirectory 'matlab' contains
.m files that define the NDK Matlab API.



Installation on Windows 7
=========================

This installation procedure has been tested on Windows 7, but not on
other versions of Windows.  On Windows, NDK requires some other
packages to be installed before you can install NDK.  You will need:

a. Python 2.7 (Anaconda) from https://www.continuum.io/downloads#windows

b. PyOpenGL packages from http://www.lfd.uci.edu/~gohlke/pythonlibs/#pyopengl

   * PyOpenGL-3.1.1-cp27-cp27m-win32.whl
   * PyOpenGL_accelerate-3.1.1-cp27-cp27m-win32.whl`

both 32 and 64 bit versions have been tested on Windows 7.


Here is the full recipe:

1. Install the 32-bit Anaconda Python Distribution for Python version 2.7.

2. Go to http://www.lfd.uci.edu/~gohlke/pythonlibs/#pyopengl and download the appropriate PyOpenGL and PyOpenGL_accelerate packages.

3. Open the Anaconda Command Shell (preferably as Administrator) and use pip to install PyOpenGL as follows (examples given for 32-bit Python 2.7 versions):

   * "pip install PyOpenGL-3.1.1-cp27-cp27m-win32.whl"
   * "pip install PyOpenGL_accelerate-3.1.1-cp27-cp27m-win32.whl"

4. *If you don't already have an NDK distribution ZIP file:*

   #. Install git (windows version) https://git-scm.com/downloads
   #. Clone the NDK repo: "git clone https://github.com/SRI-CSL/ndk.git"
   #. cd to the 'ndk' directory

5. *If you DO already have a zip file:*

   #. Unzip the NDK distribution.  This will create a directory called 'ndk-1.XXX' where XXX is the NDK sub-version number.
   #. cd to the 'ndk-1.XXX' directory

6.  Use pip to install NDK:  "pip install ."



(Note, some of these steps are thanks to stackoverflow answers. See the
thread in http://stackoverflow.com/questions/26700719/pyopengl-glutinit-nullfunctionerror)

