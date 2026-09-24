"""gRNA efficacy CNN (architecture shared with MEGA27-06; trained on Doench 2016)."""
from __future__ import annotations

import numpy as np
import torch
from torch import nn

BASES = "ACGT"
BASE_TO_IDX = {b: i for i, b in enumerate(BASES)}


def one_hot(seq: str, length: int = 30) -> np.ndarray:
    seq = seq.upper()
    arr = np.zeros((4, length), dtype=np.float32)
    for i, b in enumerate(seq[:length]):
        j = BASE_TO_IDX.get(b)
        if j is not None:
            arr[j, i] = 1.0
    return arr


def gc_content(seq: str) -> float:
    s = seq.upper()
    acgt = sum(s.count(b) for b in "ACGT")
    return 0.0 if acgt == 0 else (s.count("G") + s.count("C")) / acgt


def dinucleotide_features(seq: str) -> np.ndarray:
    s = seq.upper()
    feats = np.zeros(16, dtype=np.float32)
    counts = 0
    for i in range(len(s) - 1):
        a, b = BASE_TO_IDX.get(s[i]), BASE_TO_IDX.get(s[i + 1])
        if a is not None and b is not None:
            feats[a * 4 + b] += 1.0
            counts += 1
    if counts:
        feats /= counts
    return feats


def make_aux_features(seq: str) -> np.ndarray:
    return np.concatenate([[gc_content(seq)], dinucleotide_features(seq)]).astype(np.float32)


class GuideEfficacyCNN(nn.Module):
    def __init__(self, seq_len: int = 30, channels=(32, 64, 64), kernel: int = 5,
                 hidden: int = 128, n_aux: int = 17):
        super().__init__()
        layers: list[nn.Module] = []
        in_ch = 4
        for ch in channels:
            layers += [nn.Conv1d(in_ch, ch, kernel_size=kernel, padding=kernel // 2),
                       nn.BatchNorm1d(ch), nn.ReLU()]
            in_ch = ch
        self.conv = nn.Sequential(*layers)
        self.pool = nn.AdaptiveAvgPool1d(1)
        self.head = nn.Sequential(nn.Linear(in_ch + n_aux, hidden), nn.ReLU(),
                                  nn.Dropout(0.2), nn.Linear(hidden, 1))

    def forward(self, x: torch.Tensor, aux: torch.Tensor) -> torch.Tensor:
        h = self.pool(self.conv(x)).squeeze(-1)
        return self.head(torch.cat([h, aux], dim=1)).squeeze(-1)


class EfficacyScorer:
    """Loads a trained GuideEfficacyCNN checkpoint and scores 30-mer contexts."""

    def __init__(self, checkpoint_path: str):
        self.model = GuideEfficacyCNN()
        self.model.load_state_dict(torch.load(checkpoint_path, map_location="cpu"))
        self.model.eval()

    def score(self, contexts: list[str]) -> np.ndarray:
        X = torch.from_numpy(np.stack([one_hot(c, 30) for c in contexts]))
        A = torch.from_numpy(np.stack([make_aux_features(c) for c in contexts]))
        with torch.no_grad():
            return self.model(X, A).numpy()
