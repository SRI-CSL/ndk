audio		 - Generates .wav files from .db files

cluster-dpgmm	 - Uses Dirichlet-process Gaussian mixture models to cluster .db files.

cluster-gmm	 - Uses plain Gaussian mixture models to cluster .db files.

cluster-kmeans	 - Uses K-means to cluster .db files.

createdb	 - Shell script to populate a .db file from .smr or .wav files.
		   If you have raw recordings, this script is an easy way to create .db files.

fix-recname	 - Fixes the .db recording file designation to be a simple filename.

mkcoefs		 - Populates the waveform coefficients in a .db file.

mkdb		 - Initializes a .db file using a specified raw data recording.

mklabels	 - Initializes the cluster labels in a .db file.

mkspikes	 - Populates a .db file with spike timestamps and descriptors.

mkspikes2	 - Same as 'mkspikes', but uses Theano for fast convolutions.

smr2wav		 - Creates a .wav file from one channel of a .smr file.

vcoefs		 - Plots the coefficient space and waveform reconstructions
		   for spikes in a .db file.

vrasters	 - Renders a .db file as a scrollable raster plot, separating
		   spikes by cluster labels.

vscatter	 - Renders a .db file as a point cloud using the channel and
		   coefficient spaces.  Good for checking cluster quality.
