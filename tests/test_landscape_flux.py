import numpy as np

from hpcs_model.landscape_flux.density import estimate_density
from hpcs_model.landscape_flux.potential import compute_potential


def test_potential_no_nans():
    Z = np.random.randn(50, 2)
    cfg = {"landscape_flux": {"density_k": 5, "eps": 1e-8}}
    density, _ = estimate_density(Z, cfg)
    U = compute_potential(density, cfg)
    assert np.isfinite(U).all()
