from pkg_resources import get_distribution
__all__ = [ 'es',  'ds', 'features', 'ui', 'cluster' ]
__version__ = get_distribution('ndk').version

def version():
    return __version__
