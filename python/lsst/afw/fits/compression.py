__all__ = ["imageCompressionDisabled"]

from contextlib import contextmanager
from ._fits import getAllowImageCompression, setAllowImageCompression


@contextmanager
def imageCompressionDisabled():
    old = getAllowImageCompression()
    setAllowImageCompression(False)
    yield
    setAllowImageCompression(old)
