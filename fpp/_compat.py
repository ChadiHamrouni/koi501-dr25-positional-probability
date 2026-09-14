"""Compatibility shim for triceratops 1.0.20. Import before triceratops.

triceratops pins pytransit 2.2, which imports names that numpy and scipy have
since removed: numpy.int, numpy.float, numpy.bool (removed in numpy 1.24);
numpy.NaN, numpy.Inf, numpy.float_ (numpy 2.0); scipy.integrate.trapz (scipy
1.14). Each was a rename, so restoring the old name changes no behaviour.
"""

import numpy as np
import scipy.integrate

_NUMPY_ALIASES = {
    "int": int, "float": float, "bool": bool, "object": object,
    "str": str, "complex": complex,
    "NaN": np.nan, "NAN": np.nan, "Inf": np.inf, "infty": np.inf,
    "float_": np.float64, "complex_": np.complex128, "unicode_": np.str_,
    "string_": np.bytes_, "bool8": np.bool_,
}
_SCIPY_INTEGRATE_ALIASES = {
    "trapz": "trapezoid", "cumtrapz": "cumulative_trapezoid", "simps": "simpson",
}

for _name, _value in _NUMPY_ALIASES.items():
    if not hasattr(np, _name):
        setattr(np, _name, _value)
for _old, _new in _SCIPY_INTEGRATE_ALIASES.items():
    if not hasattr(scipy.integrate, _old) and hasattr(scipy.integrate, _new):
        setattr(scipy.integrate, _old, getattr(scipy.integrate, _new))
