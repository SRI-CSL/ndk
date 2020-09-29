Tips
====

In case you run into trouble, here are some accumulated tips:

* 'virtualenv' can be your friend.  It can help isolate problems and will create an environment you can trust (with contents you can trust).

* On Macs, there are a few different ways to install Python, and this can be a pain.  I found this out the hard way when mixing a "native" Mac Python distro with "brew".  Ultimately, I was able to get everything to work within a 'virtualenv' of a 'brew' python.

* Be aware that 'matplotlib' might (probably?) needs guidance before it will run properly within a 'brew'ed python environment.  Put this line:
        backend: TkAgg
  into your ~/.matplotlib/matplotlibrc file.  If you do not, you will see complaints from matplotlib, and OpenGL could hang.
