**Added:**

* <news item>

**Changed:**

* Change the characteristic functions in ``diffpy.srfit.pdf.characteristicfunctions`` to emit a ``RuntimeWarning`` when a non-physical shape parameter makes them return zero, since a zero return flattens the fit residual and stalls a refinement without any other sign to the user. The warning is issued once per process for each distinct problem so a refinement loop does not repeat it.
* Change ``lognormal_spherical_particle`` in ``diffpy.srfit.pdf.characteristicfunctions`` to return zero for a negative ``particle_diameter_sigma`` instead of silently returning ``spherical_particle``. A ``particle_diameter_sigma`` of zero is still the sphere limit.
* Change ``sheet_particle`` in ``diffpy.srfit.pdf.characteristicfunctions`` to return an array of zeros for a non-positive ``sheet_thickness`` when ``r`` is an array, instead of a scalar zero.

**Deprecated:**

* <news item>

**Removed:**

* <news item>

**Fixed:**

* Fix ``shell_particle`` in ``diffpy.srfit.pdf.characteristicfunctions`` returning a smoothly varying branch of unphysical values, some below negative one, for a negative ``thickness`` or ``radius``, which a refinement could mistake for a real solution. It now returns zero.
* Fix ``shell_particle`` in ``diffpy.srfit.pdf.characteristicfunctions`` returning one at every ``r`` for a zero ``thickness``, where the normalizing denominator vanishes. It now returns zero.
* Fix ``spheroidal_particle`` in ``diffpy.srfit.pdf.characteristicfunctions`` raising ``ZeroDivisionError`` for a zero ``equatorial_radius``. It now returns zero.

**Security:**

* <news item>
