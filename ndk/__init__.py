#from pkg_resources import get_distribution
import importlib.metadata

__all__ = [ 'es',  'ds', 'features', 'ui', 'cluster' ]
__version__ = importlib.metadata.version('ndk')

def version():
    return __version__
