"""Resolution conventions, not a claim that a wall film has been predicted."""
def classify_film(thickness, epsilon, wall_spacing):
    if min(epsilon, wall_spacing) <= 0 or thickness < 0:
        raise ValueError("Invalid film geometry")
    if thickness < 2 * epsilon or thickness < 4 * wall_spacing:
        return "unresolved"
    if thickness < 4 * epsilon or thickness < 6 * wall_spacing:
        return "marginal"
    return "resolved_candidate"


def nusselt_flux(h, g, nu, slip_length=0.0):
    """Single-liquid, stress-free surface, downward volumetric flux per width.

    Navier slip is explicit. The user no-slip formula is recovered at Ls=0.
    This does not silently account for gas shear, curvature or film end effects.
    """
    if min(h, g, slip_length) < 0 or nu <= 0:
        raise ValueError("Invalid falling-film parameters")
    return g / nu * (h**3 / 3 + slip_length * h**2)
