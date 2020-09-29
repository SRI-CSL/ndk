Matlab API
==========

The NDK supports Matlab through the Matlab-Python interface introduced
in Matlab R2014b.  Certain functions are provided for convenience to
avoid having to import the Python functions or classes.  Note that
data store and event store objects have methods (previously described
in the module documentation pages) that can be called directly from
Matlab.

In Matlab, the following functions are defined in
the file 'ndk.m':

.. function:: ndk() 

   Returns the NDK version number as a string.

.. function:: ndk_open(filename)

   Returns a database object corresponding to the specified filename.

.. function:: ndk_amp_min_max(dbobj, channel)

   Returns the min and max *reconstructed* amplitude values on the given channel

.. function:: ndk_dataset_interval(dbobj, channel)

   Returns a list [tmin,tmax] of sample numbers for the dataset on the
   given channel.  That is, the dataset starts at sample number
   *tmin*, and ends at sample number *tmax*.  *(Note: the API doesn't
   yet provide access to the sampling rate, but that would be used to
   convert sample numbers to seconds or milliseconds)*

.. function:: ndk_get_waveform(dbobj, spike_id, channel)

   Returns an array of values for the waveform of spike *spike_id* on
   the specified *channel*.  This is a reconstructed waveform that
   approximates the original signal.

.. function:: ndk_isi(dbobj, channel)

   Returns an array of double floats that represent the interspike
   intervals on the given channel.

.. function:: ndk_waveform_minmax(dbobj, channel)

   Returns a list [wmin,wmax] containing the waveform minimum and
   maximum values for the given *channel*.
