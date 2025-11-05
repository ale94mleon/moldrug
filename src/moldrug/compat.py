"""
This module provides backward compatibility with moldrug <= 4.0.

It allows pickle objects created with older versions of moldrug
to be successfully loaded after internal refactoring (e.g. moved classes).
"""
import dill
import importlib
from typing import Dict, Tuple

# Mapping from old (module, class) → new (module, class)
_CLASS_RENAMES: Dict[Tuple[str, str], Tuple[str, str]] = {
    ("moldrug.utils", "GA"): ("moldrug.opt", "GA"),
    ("moldrug.utils", "Local"): ("moldrug.opt", "Local"),
    # Add more mappings here as needed
}


class BackCompatUnpickler(dill.Unpickler):
    """Custom unpickler that redirects renamed or moved classes."""

    def find_class(self, module: str, name: str):
        key = (module, name)
        if key in _CLASS_RENAMES:
            new_module, new_name = _CLASS_RENAMES[key]
            mod = importlib.import_module(new_module)
            return getattr(mod, new_name)
        return super().find_class(module, name)
