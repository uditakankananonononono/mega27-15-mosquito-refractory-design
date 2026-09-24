import numpy as np

from mosqdesign.drive_sim import homing_recurse, simulate, wright_fisher


def test_allele_freqs_sum_to_one():
    w, d, r = homing_recurse(0.9, 0.1, 0.0, h=0.99, e=0.5, c_hom=1.0)
    assert abs(w + d + r - 1.0) < 1e-9


def test_drive_invades_from_rare():
    traj = simulate(0.99, 0.01, 0.0, h=0.99, e=0.01, c_hom=0.2, generations=40)
    assert traj[:, 1].max() > 0.9  # drive sweeps before resistance erodes it


def test_no_homing_no_spread():
    traj = simulate(0.9, 0.1, 0.0, h=0.0, e=0.0, c_hom=1.0, s_het=0.5)
    assert traj[-1, 1] < 0.1


def test_high_resistance_rate_blocks_drive():
    traj = simulate(0.99, 0.01, 0.0, h=0.9, e=0.99, c_hom=1.0, generations=60)
    assert traj[-1, 2] > 0.5


def test_full_sterile_female_drive_suppresses_like_kyrou():
    traj = simulate(0.99, 0.01, 0.0, h=1.0, e=0.0, c_hom=1.0, s_het=0.0, generations=20)
    assert traj[-1, 1] > 0.99


def test_wright_fisher_shape_and_seed():
    a = wright_fisher(1000, 0.99, 0.01, 0.0, h=0.99, e=0.01, c_hom=0.2, generations=10, seed=1)
    b = wright_fisher(1000, 0.99, 0.01, 0.0, h=0.99, e=0.01, c_hom=0.2, generations=10, seed=1)
    assert a.shape == (11, 3)
    assert np.allclose(a, b)
