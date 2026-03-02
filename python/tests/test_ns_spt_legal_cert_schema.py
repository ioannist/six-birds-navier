import numpy as np

from nswave.spt_legal_cert import (
    certificate_verdict,
    compute_spt_legal_metrics,
    make_spt_legal_certificate,
)


def test_spt_legal_cert_schema_on_plane_wave():
    N = 16
    L = 2 * np.pi
    u_hat = np.zeros((3, N, N, N), dtype=complex)
    amp = 1.0
    # Cosine mode at |k|=2 (shell j=1), globally delocalized in space.
    u_hat[0, 2, 0, 0] = amp
    u_hat[0, N - 2, 0, 0] = amp
    u_hat_hist = np.expand_dims(u_hat, axis=0)

    metrics = compute_spt_legal_metrics(
        u_hat_hist,
        N=N,
        L=L,
        window_c=1.0,
        energy_threshold=0.95,
    )
    thresholds = {
        "kappa_threshold": 10.0,
        "burst_threshold": 50.0,
        "require_m95": False,
        "m95_threshold": None,
    }
    config = {
        "script": "test_ns_spt_legal_cert_schema.py",
        "N": N,
        "seed": 0,
        "dt": 0.0,
        "t_max": 0.0,
        "nu": 0.0,
        "mu": 0.0,
        "alpha": 2.0,
    }
    cert = make_spt_legal_certificate(config=config, metrics=metrics, thresholds=thresholds)
    ok, reasons = certificate_verdict(cert)

    assert ok, reasons
    assert cert["verdict"] == "PASS"
    assert cert["schema_version"] == "ns_spt_legal_cert.v1"
    assert "metrics" in cert
    assert "summary" in cert
    assert "thresholds" in cert
    assert len(metrics["j_values"]) > 0
    assert np.all(np.isfinite(np.asarray(metrics["mean_conc"], dtype=float)))
    assert np.all(np.isfinite(np.asarray(metrics["max_conc"], dtype=float)))
    assert np.all(np.isfinite(np.asarray(metrics["kappa_max_over_time"], dtype=float)))
    assert np.all(np.isfinite(np.asarray(metrics["burst_ratio"], dtype=float)))
