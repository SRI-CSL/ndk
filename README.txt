The Neurogram Deconvolution Kit:

    A mostly-python library for analysis and visualization of local
    field potential recordings of peripheral nerves.  Although the
    primary focus is peripheral nerve activity, this library will also
    be exercised on Cardiac, EEG and possibly other biological
    time-series data.

Coauthors and contributors:

	  Chris Connolly
	  Maneesh Yadav
	  Andy Poggio
	  Pat Lincoln


The high-level view: Dense, raw electrical recordings are stored in
Apache Cassandra, or any format supported by the 'neo' python package.
While Cassandra was an initial data store substrate, we also use an
internal binary format called 'nbf'.  The 'nbf' format consists of
memory-mappable raw data files, along with metadata stored in pickled
python form.  Look at nbf.py for details. 

The ndk package extracts sparser discrete features (action potential
spikes, for example) from dense time series data and stores them in a
separate data store that allows such features to be classified and
visualized.  These data stores are currently restricted to sqlite3 .db
files.  These .db files always point to the raw data source from which
they're derived.

In the interests of standardization of view and efficiency of data
access, NDK also supports a memory-mapped raw data format called
'nbf'.

NDK works using Python 2.7 and Python 3.5 on Ubuntu Linux, Mac OSX,
and Windows.  See below for Windows installation instructions.  Note
that the Windows build has only been tested under Windows 7.



INSTALLATION ON LINUX AND MACOSX
================================

NDK was originally developed in Python 2.7 under Ubuntu Linux and Mac
OS X.  Current and future development will rely on Python 3, but will
be backward compatible with Python 2.7.

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
have theano installed, 'mkspikes' will do the job somewhat more slowly.

The ndk module has two submodules, ds, and es, that provide access to
raw data stores and event stores respectively (hence the names).

Matlab notes: with version R2014b and higher, python functionality is
easily accessed within Matlab.  The ndk subdirectory 'matlab' contains
.m files that define the NDK Matlab API.



NOTES FOR INSTALLING ON WINDOWS:
================================


You can do a pip install of NDK on Windows, but there's a chance that
the 'vscatter' and 'vrasters' utilities won't work.  This is because
these utilities require OpenGL.  This can be a little tricky,
primarily because PyOpenGL.GLUT depends on the correct version and
placement of the GLUT DLL.  The following recipe has been tested on a
bare Windows 7 VM using:

   o 32-bit Python 2.7 (Anaconda)
     32-bit PyOpenGL‑3.1.1‑cp27‑cp27m‑win32.whl
     32-bit PyOpenGL_accelerate‑3.1.1‑cp27‑cp27m‑win32.whl



1) Install the 32-bit Anaconda Python Distribution for Python version
   2.7 from https://www.continuum.io/downloads#windows

[Note, the following two steps are thanks to stackoverflow answers. See the thread in
  http://stackoverflow.com/questions/26700719/pyopengl-glutinit-nullfunctionerror  ]


2) Go to http://www.lfd.uci.edu/~gohlke/pythonlibs/#pyopengl and 
download these two files:
         PyOpenGL‑3.1.1‑cp27‑cp27m‑win32.whl
	 PyOpenGL_accelerate‑3.1.1‑cp27‑cp27m‑win32.whl

3) Open the Anaconda Command Shell (preferably as Administrator) and use pip to install them:
         pip install PyOpenGL‑3.1.1‑cp27‑cp27m‑win32.whl
	 pip install PyOpenGL_accelerate‑3.1.1‑cp27‑cp27m‑win32.whl

4) Install git (windows version) https://git-scm.com/downloads

5) Clone the NDK repo: "git clone https://github.com/SRI-CSL/ndk.git"

6) cd to the 'ndk' directory and install:
     pip install .




More to come!


