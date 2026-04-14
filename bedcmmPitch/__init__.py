# bedcmmPitch/__init__.py
from ._config import implementation
import warnings
from .py_impl import calc_Pitch,calc_bedcmm
import sys

if implementation != 'Cython':
    python_exe = sys.executable
    warnings.warn(
        f"\n{'='*60}\n"
        f"[bedcmm] 高速版(Cython)がビルドされていないか、読み込めません。\n"
        f"現在は低速なPython版で動作しています。\n\n"
        f"高速化を有効にするには、以下のコマンドを実行してビルドしてください：\n"
        f" {python_exe} setup.py build_ext --inplace\n"
        f"{'='*60}\n",
        RuntimeWarning
    )

__all__ = ['calc_Pitch','calc_bedcmm']
