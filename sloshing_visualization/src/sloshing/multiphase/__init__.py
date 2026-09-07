"""Independent nonlinear diffuse-interface branch. No automatic tank runs.

Heavy DOLFINx imports live in solver modules, so pure constitutive tests also run
in the historical project's environment. This package never calls SloshingSolver.
"""
