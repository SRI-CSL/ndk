This directory defines the Matlab API that fits the current version of
NDK.  The convention used here simply defines functions that generally
begin with "ndk_".


NOTES:

Example Matlab session that illustrates database opening and waveform
extraction:

>> db = py.ndk.es.db3.db3('/opt/spikes/Feinstein/db/set1/CervicalVagusRecording7.db')

db = 

  Python db3 with properties:

    sample_rate: 32000
        nspikes: 0
            cur: [1×1 py.sqlite3.Cursor]
      recording: [1×38 py.unicode]
       modified: 0
            min: [1×1 py.NoneType]
          bases: [1×1 py.dict]
            max: [1×1 py.NoneType]
      nchannels: 3
           conn: [1×1 py.sqlite3.Connection]

    <ndk.es.db3.db3 instance at 0x7f325249ee18>

>> db.get_spike_waveform(20,0)
[(20, 0, -0.00062885670708489, -0.00028810989416646704, -0.0008190827760929635, -2.6526526165779972e-05, -0.0005336186007610461)]

ans = 

  Python ndarray with properties:

           T: [1×1 py.numpy.ndarray]
        base: [1×1 py.NoneType]
      ctypes: [1×1 py.numpy.core._internal._ctypes]
        data: [1×512 py.buffer]
       dtype: [1×1 py.numpy.dtype]
       flags: [1×1 py.numpy.flagsobj]
        flat: [1×1 py.numpy.flatiter]
        imag: [1×1 py.numpy.ndarray]
    itemsize: 8
      nbytes: 512
        ndim: 1
        real: [1×1 py.numpy.ndarray]
       shape: [1×1 py.tuple]
        size: 64
     strides: [1×1 py.tuple]

    [ -2.84501158e-06  -3.54879072e-06  -6.38829565e-06  -8.67386986e-06
      -1.37364867e-05  -2.03203853e-05  -2.93256390e-05  -4.12726447e-05
      -5.66242545e-05  -7.56966417e-05  -9.85541320e-05  -1.24902463e-04
      -1.54000012e-04  -1.84612377e-04  -2.15039492e-04  -2.43229976e-04
      -2.66986657e-04  -2.84257262e-04  -2.93459696e-04  -2.93779166e-04
      -2.85389116e-04  -2.69518641e-04  -2.48304868e-04  -2.24464101e-04
      -2.00829572e-04  -1.79812781e-04  -1.62936162e-04  -1.50555468e-04
      -1.41839367e-04  -1.35044530e-04  -1.28003234e-04  -1.18743176e-04
      -1.06115725e-04  -9.01829091e-05  -7.22900677e-05  -5.48813042e-05
      -4.10134676e-05  -3.36503098e-05  -3.49941330e-05  -4.59836169e-05
      -6.60372750e-05  -9.31214718e-05  -1.24097437e-04  -1.55270412e-04
      -1.83012730e-04  -2.04303795e-04  -2.17122069e-04  -2.20637434e-04
      -2.15168091e-04  -2.01964334e-04  -1.82900139e-04  -1.60125243e-04
      -1.35745233e-04  -1.11584957e-04  -8.90451553e-05  -6.90507021e-05
      -5.20771960e-05  -3.82269229e-05  -2.73282420e-05  -1.90378224e-05
      -1.21865730e-05  -9.55910361e-06  -5.44102020e-06  -4.40565442e-06]

>> x = db.get_spike_waveform(21,0)
[(21, 0, -0.0007475351115639622, 0.0003338823623537952, -0.00041230394792025703, 0.0002051689843594639, -0.0005326511448824403)]

x = 

  Python ndarray with properties:

           T: [1×1 py.numpy.ndarray]
        base: [1×1 py.NoneType]
      ctypes: [1×1 py.numpy.core._internal._ctypes]
        data: [1×512 py.buffer]
       dtype: [1×1 py.numpy.dtype]
       flags: [1×1 py.numpy.flagsobj]
        flat: [1×1 py.numpy.flatiter]
        imag: [1×1 py.numpy.ndarray]
    itemsize: 8
      nbytes: 512
        ndim: 1
        real: [1×1 py.numpy.ndarray]
       shape: [1×1 py.tuple]
        size: 64
     strides: [1×1 py.tuple]

    [ -2.26137545e-06  -2.91186009e-06  -5.04008450e-06  -6.57848007e-06
      -1.03852588e-05  -1.50816595e-05  -2.13201364e-05  -2.93163891e-05
      -3.91766501e-05  -5.08259937e-05  -6.39356088e-05  -7.78655255e-05
      -9.16418525e-05  -1.03990015e-04  -1.13444744e-04  -1.18538048e-04
      -1.18053180e-04  -1.11321184e-04  -9.84930702e-05  -8.07178768e-05
      -6.01811361e-05  -3.99392246e-05  -2.35180428e-05  -1.43508532e-05
      -1.51423943e-05  -2.72542299e-05  -5.02877598e-05  -8.19836247e-05
      -1.18486997e-04  -1.54975889e-04  -1.86511284e-04  -2.08974301e-04
      -2.19910127e-04  -2.18993521e-04  -2.08046922e-04  -1.90685118e-04
      -1.71572786e-04  -1.55436733e-04  -1.46143893e-04  -1.46004990e-04
      -1.55409367e-04  -1.72872775e-04  -1.95435889e-04  -2.19313576e-04
      -2.40634610e-04  -2.56089190e-04  -2.63404947e-04  -2.61586967e-04
      -2.50887816e-04  -2.32579585e-04  -2.08617968e-04  -1.81260479e-04
      -1.52716750e-04  -1.24890711e-04  -9.92256656e-05  -7.66507706e-05
      -5.76127536e-05  -4.21608487e-05  -3.00563719e-05  -2.08843303e-05
      -1.34079449e-05  -1.04690645e-05  -5.84172349e-06  -4.74544605e-06]

>> x(2)
Array formation and parentheses-style indexing with objects of class 'py.numpy.ndarray' is not allowed.  Use objects of class
'py.numpy.ndarray' only as scalars or use a cell array.
 
>> x[2]
 x[2]
  ↑
Error: Unbalanced or unexpected parenthesis or bracket.
 
>> x

x = 

  Python ndarray with properties:

           T: [1×1 py.numpy.ndarray]
        base: [1×1 py.NoneType]
      ctypes: [1×1 py.numpy.core._internal._ctypes]
        data: [1×512 py.buffer]
       dtype: [1×1 py.numpy.dtype]
       flags: [1×1 py.numpy.flagsobj]
        flat: [1×1 py.numpy.flatiter]
        imag: [1×1 py.numpy.ndarray]
    itemsize: 8
      nbytes: 512
        ndim: 1
        real: [1×1 py.numpy.ndarray]
       shape: [1×1 py.tuple]
        size: 64
     strides: [1×1 py.tuple]

    [ -2.26137545e-06  -2.91186009e-06  -5.04008450e-06  -6.57848007e-06
      -1.03852588e-05  -1.50816595e-05  -2.13201364e-05  -2.93163891e-05
      -3.91766501e-05  -5.08259937e-05  -6.39356088e-05  -7.78655255e-05
      -9.16418525e-05  -1.03990015e-04  -1.13444744e-04  -1.18538048e-04
      -1.18053180e-04  -1.11321184e-04  -9.84930702e-05  -8.07178768e-05
      -6.01811361e-05  -3.99392246e-05  -2.35180428e-05  -1.43508532e-05
      -1.51423943e-05  -2.72542299e-05  -5.02877598e-05  -8.19836247e-05
      -1.18486997e-04  -1.54975889e-04  -1.86511284e-04  -2.08974301e-04
      -2.19910127e-04  -2.18993521e-04  -2.08046922e-04  -1.90685118e-04
      -1.71572786e-04  -1.55436733e-04  -1.46143893e-04  -1.46004990e-04
      -1.55409367e-04  -1.72872775e-04  -1.95435889e-04  -2.19313576e-04
      -2.40634610e-04  -2.56089190e-04  -2.63404947e-04  -2.61586967e-04
      -2.50887816e-04  -2.32579585e-04  -2.08617968e-04  -1.81260479e-04
      -1.52716750e-04  -1.24890711e-04  -9.92256656e-05  -7.66507706e-05
      -5.76127536e-05  -4.21608487e-05  -3.00563719e-05  -2.08843303e-05
      -1.34079449e-05  -1.04690645e-05  -5.84172349e-06  -4.74544605e-06]

>> 

This incantation in Matlab will convert the array into something
Matlab wants to see - looks like Python is doing much of the work, but
the end result is a Matlab double array:

 y = double(py.array.array('d', py.numpy.nditer(x, pyargs('order','F'))))

There are other quirks - be aware that numerical args to Python
functions often seem to end up as floats if Matlab is calling the
function.  This is particularly troublesome when the arg is used as an
index.  Having the function coerce these to ints is probably the
easiest way to deal with it.
