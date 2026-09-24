## Per-candidate resistance-rate (e) prediction

Model: e = p_EJ x w, w = fraction of cut-site microhomology deletions that remove PAM-GG or seed sequence. p_EJ swept (not fitted) over 0.0001, 0.001, 0.01, 0.05, 0.2. L = upper-bound single-substitution resistance target size (PAM GG 2x3 + seed 8x3). r0 = standing compromised-allele fraction (Ag1000G).

| guide | L | MH dels (res/total) | w | r0 | e interval | gen to 95% drive (e_low/e_high) | gen to 50% R (e_low/e_high) |
|---|---|---|---|---|---|---|---|
| dsx-v3-1 | 30 | 136/396 | 0.343 | 0.013605442176870748 | [3.43e-05, 6.87e-02] | >200/>200 | 13/12 |
| dsx-v3-2 | 30 | 70/389 | 0.180 | 0.0 | [1.80e-05, 3.60e-02] | 12/>200 | 16/14 |
| kyrou | 30 | 64/303 | 0.211 | 0.0 | [2.11e-05, 4.22e-02] | 12/>200 | 16/14 |
