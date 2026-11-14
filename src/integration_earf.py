"""Evidence-Anchored Reliability Fusion (EARF): pick the controller's signal ambiguity set from measured evidence.

Protocol: paper/integration_earf_protocol_2026-11-14.md. Simulation only.
Anchors are fixed numbers copied from result files (cited below); this module reads no patient data.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from trackc_integrated_signal_model import SignalModelParameters, standard_normal_marginal_signature


@dataclass(frozen=True)
class Anchor:
    rho: float
    n: int

    @property
    def z(self) -> float:
        return float(np.arctanh(self.rho))

    @property
    def se(self) -> float:
        return float(1.0 / np.sqrt(self.n - 3))


# output/trackA_hgcc_tmz_102/result.json: rho -0.102, n 83
MOLECULAR_ANCHOR = Anchor(rho=-0.102, n=83)
# output/integration_methoddev/confirmation.json: combined rho 0.114, n 23+26
MRI_ANCHOR = Anchor(rho=0.114, n=49)
MASS = 0.90


def _logpdf(x: float, a: Anchor) -> float:
    z = np.arctanh(np.clip(x, -0.999, 0.999))
    return float(-0.5 * ((z - a.z) / a.se) ** 2)


def evidence_weights(models: Sequence[SignalModelParameters], *, molecular=MOLECULAR_ANCHOR, mri=MRI_ANCHOR) -> np.ndarray:
    """Normalized evidence weight per model (q uniform)."""
    logw = np.array([
        _logpdf(m.a, molecular) + _logpdf(m.b, mri) + _logpdf(standard_normal_marginal_signature(m)[1], mri)
        for m in models
    ])
    w = np.exp(logw - logw.max())
    return w / w.sum()


Z90 = 1.6448536269514722


def credible_box_model_set(q_levels: Sequence[float], *, molecular=MOLECULAR_ANCHOR, mri=MRI_ANCHOR,
                           z: float = Z90) -> tuple:
    """EARF ambiguity set: a, b and the MRI-molecular link rho each at (low, estimate, high) of the
    90% Fisher-z interval; d solved from rho; x q levels. 27 x len(q_levels) models."""
    def levels(an):
        return tuple(float(np.tanh(an.z + k * z * an.se)) for k in (-1, 0, 1))
    out = []
    for q in q_levels:
        for a in levels(molecular):
            for b in levels(mri):
                for rho in levels(mri):
                    d = float(np.clip((rho - a * b) / np.sqrt((1 - a * a) * (1 - b * b)), -0.95, 0.95))
                    out.append(SignalModelParameters(q=float(q), a=a, b=b, d=d))
    return tuple(out)


def evidence_anchored_model_set(models: Sequence[SignalModelParameters], *, mass: float = MASS, **kw) -> tuple:
    """Smallest set of models with weight >= mass; stable grid order breaks ties."""
    models = tuple(models)
    w = evidence_weights(models, **kw)
    order = sorted(range(len(models)), key=lambda i: (-w[i], i))
    keep, total = [], 0.0
    for i in order:
        keep.append(i)
        total += w[i]
        if total >= mass - 1e-12:
            break
    return tuple(models[i] for i in sorted(keep))


def sample_evidence_world(rng: np.random.Generator, q_levels: Sequence[float], *,
                          molecular=MOLECULAR_ANCHOR, mri=MRI_ANCHOR) -> SignalModelParameters:
    """Draw one evidence_anchored world's true (q, a, b, d)."""
    a = float(np.tanh(rng.normal(molecular.z, molecular.se)))
    b = float(np.tanh(rng.normal(mri.z, mri.se)))
    rho = float(np.tanh(rng.normal(mri.z, mri.se)))
    denom = np.sqrt(max(1e-12, (1 - a * a) * (1 - b * b)))
    d = float(np.clip((rho - a * b) / denom, -0.95, 0.95))
    q = float(rng.choice(np.asarray(q_levels, float)))
    return SignalModelParameters(q=q, a=a, b=b, d=d)
