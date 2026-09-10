import numpy as np
import pytest
from janpa.gto.polynom3d import Polynom3D
from janpa.gto.polynom_rho_z import Polynom_rho_z

def test_polynom3d_creation():
    p = Polynom3D()
    assert p.n_terms == 1
    assert p.coefs[0] == 0.0

def test_polynom3d_multiply():
    # x * y = xy
    p1 = Polynom3D(np.array([1.0]), np.array([[1, 0, 0]]))
    p2 = Polynom3D(np.array([1.0]), np.array([[0, 1, 0]]))
    p1.multiply_by(p2)
    assert p1.n_terms == 1
    assert p1.powers[0, 0] == 1
    assert p1.powers[0, 1] == 1
    assert p1.powers[0, 2] == 0

def test_polynom3d_evaluate():
    # 2*x^2 + 3*y at (2, 3, 0) -> 2*4 + 3*3 = 8 + 9 = 17
    p = Polynom3D(np.array([2.0, 3.0]), np.array([[2, 0, 0], [0, 1, 0]]))
    val = p.evaluate_at_point([2.0, 3.0, 0.0])
    assert np.isclose(val, 17.0)

def test_polynom_rho_z_stack():
    pz = Polynom_rho_z(2, 2)
    pz.push_to_stack()
    pz.mul_lin(1.0) # multiply by (z - 1.0)
    assert pz.last_active_z_power == 1
    pz.restore_from_stack(True)
    assert pz.last_active_z_power == 0