## SBLIPMS-Impulse-KH-Controllability/

## Kalman–Hewer Controllability Equivalence for Singular Bilinear Lyapunov Periodic Systems with Impulses
DOI: https://doi.org/10.5281/zenodo.20196726 
**Authors:** Sri Venkata Durga Sudarsan Madhyannapu¹ and Pradheep Kumar S.²

¹ Freshmen Engineering Department, Dr. RVR NRI Institute of Technology (Deemed to be University), Pothavarappadu, Agiripalli, Eluru District 521212, Andhra Pradesh, India. Email: msvdsudarsan@gmail.com · ORCID: 0009-0001-2126-6428

² School of Basic Sciences, SRM University AP, Neerukonda, Mangalagiri, Guntur 522240, Andhra Pradesh, India. Email: sravanampradheepkumar@gmail.com

**Target Journal:** Nonlinear Analysis: Hybrid Systems (Elsevier, ISSN 1751-570X) · IF 4.8 · Q1 · SCI/SCIE

**Status:** Submitted to Nonlinear Analysis: Hybrid Systems (April 2026)

---

## Abstract

Impulsive events — thruster firings, drug bolus injections, and gear-tooth impacts — create discontinuous jumps in a system's controllability manifold, a genuinely nonlinear dynamical phenomenon that can instantaneously restore controllability when continuous dynamics fail. No unified framework exists for systems that combine all five structural features: singular descriptor pencil (rank(E) = r < n), Lyapunov bilinear coupling, impulsive resets, periodic coefficients, and Kalman–Hewer controllability equivalence.

This paper introduces the **Singular Bilinear Lyapunov Impulsive Periodic Matrix Differential System (SBLIPMS)** and establishes five results:

1. An impulsive reachability Gramian **W_c^imp** and Kalman controllability rank criterion
2. Kalman–Hewer equivalence proved via a monodromy argument, without Kronecker factorisations
3. An **Impulse Recovery Theorem** with computable index κ* ≤ 0 certifying when impulses restore full rank
4. Stability of controllability rank under timing jitter up to ±8%T (a 432 s window within a 5400 s orbital period) and 10% parameter perturbations
5. A Kronecker-free O((N+q)n³) algorithm achieving **5,400× speedup** at n=16 over the classical O(Nn⁶) approach

A 16-dimensional satellite attitude covariance stress test and a pharmacokinetic compartmental model validate the framework.

---

## Repository Structure

```
SBLIPMS-Impulse-KH-Controllability/
│
├── README.md                          ← This file
│
├── MATLAB_Codes/
│   ├── SBLIPMSMain.m                  ← Main script: all 4 scenarios + Figures 2–5
│   ├── Figure1_SatelliteSV.m          ← Figure 1: Satellite SV trajectory (16-D)
│   ├── algo1_gramian.m                ← Algorithm 1: Kronecker-free Gramian assembly
│   ├── build_STM_eff.m                ← State transition matrix builder (with impulses)
│   ├── compute_Mk.m                   ← Impulsive jump operator M_k computation
│   └── rk4step.m                      ← RK4 stepper (handles singular E via pinv)
│
├── MATLAB_Outputs/
│   └── MATLAB_OUTPUTS_IMPULSE_KH.txt  ← Verified numerical output values
│
├── Figures/
│   ├── Figure1_SatelliteSV.pdf        ← Singular value trajectories (satellite)
│   ├── Figure2_SigmaVsQ.pdf           ← σ_min vs number of impulses q
│   ├── Figure3_RankRecovery.pdf       ← Rank recovery: continuous vs impulsive
│   ├── Figure4_Jitter.pdf             ← Jitter sensitivity (±8%T threshold)
│   └── Figure5_Speedup.pdf            ← Algorithm 1 speedup vs dimension n
│
└── LICENSE                            ← MIT License
```

---

## Numerical Results Summary

### Table 1 — Scenario Comparison (4×4 System, n=4, q=3 impulses)

| Scenario | Input | σ_min(W_c) | Controllable? |
|---|---|---|---|
| A | Strong continuous only | 0.2102 | ✅ Yes |
| B | Strong continuous + impulses | 0.7238 | ✅ Yes (+245%) |
| C | Weak continuous only | 0.0021 ≈ 0 | ❌ No |
| D | Weak continuous + impulses | 0.5032 | ✅ Yes (rank recovered) |

**Scenario D is the key result:** impulses alone restore full controllability when continuous dynamics are deficient (κ* = −2 ≤ 0).

### Table 2 — Algorithm 1 Speedup (Kronecker-free vs Classical)

| Dimension n | Algorithm 1 time | Classical time | Speedup |
|---|---|---|---|
| 4 | 0.31 s | 23.6 s | **76×** |
| 6 | 1.08 s | 337 s | **312×** |
| 8 | 4.23 s | 1,820 s | **430×** |
| 16 | 868 s (14 min) | 4.7×10⁶ s (54 days) | **5,400×** |

### Table 3 — Jitter Sensitivity

| Jitter η | σ_min | Rank preserved? |
|---|---|---|
| 0% | 0.503 | ✅ Yes |
| 4% | 0.481 | ✅ Yes |
| 8% | 0.447 | ✅ Yes |
| 12% | 0.000 | ❌ No (bifurcation) |

**Design safety margin: ±8%T** (432 s per 5400 s orbital period).

---

## How to Reproduce All Results

### Requirements
- MATLAB R2021b or later (R2024b recommended)
- No additional toolboxes required

### Steps

**Step 1: Generate Figures 2–5 and Table 1**
```matlab
run('SBLIPMSMain.m')
```
This produces `Figure2_SigmaVsQ.pdf`, `Figure3_RankRecovery.pdf`, `Figure4_Jitter.pdf`, `Figure5_Speedup.pdf` and prints the scenario table to the console.

**Step 2: Generate Figure 1 (Satellite Singular Value Trajectory)**
```matlab
run('Figure1_SatelliteSV.m')
```
This produces `Figure1_SatelliteSV.pdf` — the 16-dimensional satellite stress test.

Step 3: Verify console output.

Compare the console output with `MATLAB_OUTPUTS_IMPULSE_KH.txt` to confirm reproducibility.
---

## System Definition

The SBLIPMS is governed by:

```
E·Ẋ(t) = A(t)X(t) + X(t)A(t)ᵀ + Σᵢ uᵢ(t)NᵢX(t) + F(t)u(t),   t ≠ tₖ
X(tₖ⁺) = (I + Bₖ)X(tₖ⁻)(I + Cₖᵀ) + Dₖvₖ,                       k = 1,...,q
```

where:
- `E ∈ ℝⁿˣⁿ` is singular with rank(E) = r < n
- `A(t)` is T-periodic
- `tₖ` are impulse instants within [0, T]
- `Dₖvₖ` is the impulsive input injection

**Parameters used in the 4×4 example:**
```matlab
E = diag([1, 1, 0, 0])          % rank-2 singular pencil
A(t) = [-2, cos(t), 0, 0;
         sin(t), -1, 1, 0;
         0, 0, -1, cos(t);
         0, 0, 0, 1]             % T-periodic, T=1
t_imp = [0.25, 0.50, 0.75]      % three impulse instants
B = diag([0.1, -0.05, 0, 0])    % state jump matrix
D = [1,0; 0,1; 0,0; 0,0]        % input distribution matrix
```

---

## Key Theoretical Results

### Theorem 1 — Kalman Controllability Criterion
The SBLIPMS is Kalman controllable if and only if:
```
rank(W_c^imp(0,T)) = r
```
where `W_c^imp = W_c^cts + W_c^jump` is the composite impulsive Gramian.

### Theorem 2 — Kalman–Hewer Equivalence
For SBLIPMS, Kalman controllability ⟺ Hewer controllability, proved via monodromy argument without any Kronecker product decomposition.

### Theorem 3 — Impulse Recovery Theorem
Define the recovery index:
```
κ* = r - rank(W_c^cts) - Σₖ rank(M_k Dₖ projected onto S⊥)
```
If κ* ≤ 0, impulsive inputs restore full controllability. The index κ* is computable via Algorithm 1.

### Algorithm 1 — Kronecker-Free Gramian Assembly
**Complexity:** O((N+q)n³) vs classical O(Nn⁶)  
**Key idea:** Direct RK4 integration of the state transition matrix using `E⁺ = pinv(E)`, without lifting to the n²-dimensional Kronecker product space.

---

## Applications

### Satellite Station-Keeping (16-dimensional)
The 16-dimensional attitude covariance system demonstrates:
- Two thruster firings per orbit restore full controllability
- Classical Kronecker approach: ~54 days; Algorithm 1: ~14 minutes
- Speedup: **5,400×**

### Pharmacokinetic Compartmental Model
Three-compartment model with bolus injections:
- Bolus at t₁ = 0.25T, t₂ = 0.50T, t₃ = 0.75T
- Each bolus feasible with standard infusion pump (finite ‖Dₖvₖ‖)
- Rank recovery confirmed: κ* = −2 ≤ 0

---

## Companion Papers

This repository is part of a research series on Kalman–Hewer equivalence for structured matrix systems:

| Paper | System Class | Journal | Repository |
|---|---|---|---|
| **This paper** | SBLIPMS (Lyapunov bilinear + impulses) | Nonlinear Analysis: Hybrid Systems | — |
| [Bilinear-Matrix-Periodic-Controllability](https://github.com/msvdsudarsan/Bilinear-Matrix-Periodic-Controllability) | Generalized bilinear periodic | MCSS | ✅ |
| [SBMPMS-observability](https://github.com/msvdsudarsan/SBMPMS-observability) | Melnikov observability breakdown | Chaos, Solitons & Fractals | ✅ |

---

## Citation

If you use this code or data in your research, please cite:

```bibtex
@article{Madhyannapu2026sblipms,
  author  = {Madhyannapu, Sri Venkata Durga Sudarsan and {Pradheep Kumar}, S.},
  title   = {{Kalman--Hewer} Controllability Equivalence for Singular
             Bilinear {Lyapunov} Periodic Systems with Impulses},
  journal = {Nonlinear Analysis: Hybrid Systems},
  year    = {2026},
  publisher = {Elsevier},
  note    = {Submitted April 2026}
}
```

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

## Contact

**Sri Venkata Durga Sudarsan Madhyannapu**  
Email: msvdsudarsan@gmail.com  
ORCID: [0009-0001-2126-6428](https://orcid.org/0009-0001-2126-6428)  
Institution: Dr. RVR NRI Institute of Technology (Deemed to be University), Andhra Pradesh, India
