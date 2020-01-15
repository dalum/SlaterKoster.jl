##################################################
# Helper functions
##################################################

"""
    spline_catmullrom(u, y0, y1, y2, y3)

Smooth interpolation between `y1, y2`, using the Catmull-Rom spline.
The function uses `y0` and `y3` to ensure that the derivative at `y1`
and `y2` is continuous.  `u` is a parameter between `0` and `1`.

"""
spline_catmullrom(u, y0, y1, y2, y3) = y1 + 0.5u*(-y0 + y2 + u*(2y0 - 5y1 + 4y2 - y3 + u*(-y0 + 3y1 - 3y2 + y3)))

_orbitalrank(orbital_symbol) =
    startswith(string(orbital_symbol), "s") ? 1 :
    startswith(string(orbital_symbol), "p") ? 2 :
    startswith(string(orbital_symbol), "d") ? 3 :
    startswith(string(orbital_symbol), "f") ? 4 :
    error("invalid orbital symbol: $orbital_symbol")

##################################################
