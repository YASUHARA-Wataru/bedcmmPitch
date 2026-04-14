# bedcmmPitch/_config.py
import importlib.util

# Cython版が存在するかチェック
_has_cython = importlib.util.find_spec(".cy_impl", package="bedcmmPitch") is not None
implementation = "Cython" if _has_cython else "Python"