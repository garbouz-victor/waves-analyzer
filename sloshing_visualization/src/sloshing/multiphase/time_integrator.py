"""Fixed-step BE/BDF2 coefficients; validation never silently changes dt."""
def derivative_coefficients(step, dt, scheme="bdf2"):
    if step < 1 or dt <= 0 or scheme not in ("be", "bdf2"):
        raise ValueError("Invalid time integration arguments")
    if step == 1 or scheme == "be":
        return 1 / dt, -1 / dt, 0.0
    return 1.5 / dt, -2 / dt, 0.5 / dt
