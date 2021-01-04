from distutils.core import setup
setup(name='ndk',
      version='2.155',
      description='Neurogram Deconvolution Kit',
      author='Chris Connolly',
      author_email='chris@neurome.com',
      url='http://www.neurome.com/',
      packages = [ 'ndk', 'ndk.es', 'ndk.ds' ],
      install_requires = [
          'neo',
          'cqlsh',
          'numpy',
          'scipy',
		  'pandas',
          'seaborn',
          'joblib',  # wtf?
#          'sqlite3',
          'PyMySQL',
          'sklearn',
          'uritools',
          'pyOpenGL',
          'quantities',
          'matplotlib'
      ],
      scripts = [ 'bin/audio',
                  'bin/addbasis',
                  'bin/isi',
                  'bin/mkdb',
                  'bin/mkspikes',
                  'bin/mkcoefs',
                  'bin/mklabels',
                  'bin/mkevents',
                  'bin/cluster-kmeans',
                  'bin/cluster-gmm',
                  'bin/cluster-dpgmm',
                  'bin/cluster-dbscan',
                  'bin/profile', 
                  'bin/profile2', 
                  'bin/rmspikes',
                  'bin/rebasis', 
                  'bin/reorg', 
                  'bin/smrevents',
                  'bin/spikepca',
                  'bin/vcoefs',
                  'bin/vsignal',
                  'bin/cplot',
                  'bin/vrasters',
                  'bin/vscatter',
                  'bin/vwave',
                  'bin/ndk2wav',
                  'bin/neo2nbf'
                  ]
)
