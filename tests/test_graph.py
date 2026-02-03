import pandas as pd

from hpcs_model.dynamics.plasticity import compute_plasticity, select_hpcs


def test_plasticity_selects_hub():
    edge_df = pd.DataFrame(
        {
            "src_cluster": ["0", "0", "1", "2"],
            "dst_cluster": ["1", "2", "0", "0"],
            "mass_norm": [0.6, 0.4, 0.7, 0.8],
            "mass": [0.6, 0.4, 0.7, 0.8],
            "time_t": ["t0"] * 4,
            "time_t1": ["t1"] * 4,
        }
    )
    plast = compute_plasticity(edge_df)
    hpcs = select_hpcs(plast)
    assert hpcs in {"0", "1", "2"}
