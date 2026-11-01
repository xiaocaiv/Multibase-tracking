# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## How to Run

```bash
# Run the full simulation (MATLAB R2024a):
matlab -batch "main"
```

All files must be on the MATLAB path (they are in a flat directory). Requires Statistics and Machine Learning Toolbox (`chi2inv`, `poissrnd`) and Signal Processing Toolbox (`spectrogram`, `hilbert`, `butter`, `filtfilt` — only when `use_features=true` with wav import).

## Architecture

Dual-platform cooperative multi-target bearing-only tracking system. Two stationary platforms observe bearing angles to moving targets. UKF handles nonlinear bearing measurement, JPDA handles measurement-to-track association, triangulation initializes new tracks.

### Data Flow (per time step)

```
generate_measurements  ──→  z1, z2  (struct arrays with .bearing, .spectrum)
    │
    ├─  UKF predict (all tracks):  F·x,  F·P·Fᵀ + Q
    │
    ├─  Platform 1:  ukf_measurement_predict  →  JPDA  →  UKF update
    ├─  Platform 2:  ukf_measurement_predict  →  JPDA  →  UKF update  (sequential)
    │
    └─  manage_tracks:  delete stale  +  triangulate new from unassociated pairs
```

### Module Dependency Order (implementation / modification order)

| Layer | Files | Role |
|-------|-------|------|
| **Params** | `init_parameters.m` | All simulation, UKF, JPDA, and feature parameters |
| **Ground Truth** | `generate_truth.m` | Deterministic CV trajectories + spectrum templates (via `generate_ship_spectrum`) |
| **Measurements** | `generate_measurements.m` | Struct array: `.bearing` (rad) + `.spectrum` (96×1 normed); includes Poisson clutter |
| **UKF Core** | `generate_sigma_points.m`, `measurement_bearing.m`, `ukf_predict.m`, `ukf_measurement_predict.m` | Pure kinematic — **never modify for feature work** |
| **JPDA Core** | `feasible_events.m`, `jpda_associate.m` | Validation matrix → feasible joint events → marginal β probabilities |
| **Orchestration** | `jpda_process_platform.m` | Per-platform: measurement predict → JPDA → UKF update → spectrum update |
| **Management** | `manage_tracks.m`, `triangulate_init.m` | Track birth (triangulation + spectrum consistency check) and death |
| **Entry** | `main.m` | Simulation loop + history storage |
| **Output** | `plot_results.m` | 6 figures: trajectory, bearings, position RMSE, velocity RMSE, track stats, spectrum templates |

### Key Design Decisions

- **Sequential platform processing**: Platform 1 updated first, Platform 2 operates on already-updated states. This avoids joint multi-sensor JPDA complexity.
- **Measurements are struct arrays**: Each element has `.bearing` (scalar) and `.spectrum` (96×1 vector). Backward-compatible via `[z1.bearing]` to extract a numeric vector.
- **Feature aid is multiplicative in JPDA**: `L_total = L_kinematic × L_feature`. The UKF state vector remains 4D (x, y, vx, vy) — features never enter the filter.
- **Spectrum similarity**: Manual Pearson correlation (no toolbox), sharpened by `exp(gamma × max(0,ρ)²)`. Gamma=5 is conservative; increase for stronger feature reliance.
- **Feature toggle**: Set `params.use_features = false` in `init_parameters.m` to revert to pure bearing-only mode. All feature-related code paths are guarded by this flag.
- **Spectrum template caching**: `generate_truth.m` calls `generate_ship_spectrum()` once and caches in `truth.spectrum{tgt}`. `generate_measurements.m` reads from cache — never re-reads wav files per frame.

### Track Struct Fields

```matlab
track.x             % 4×1 state [x; y; vx; vy]
track.P             % 4×4 covariance
track.id            % unique integer
track.confirmed     % logical (hit_count >= confirm_hits)
track.hit_count     % consecutive associations
track.miss_count    % consecutive misses (deleted at max_misses)
track.spectrum      % 96×1 normalized spectrum (NaN if feature_init==false)
track.feature_init  % logical
```

### Key Tuning Parameters

| Param | Default | Effect |
|-------|---------|--------|
| `Pd` | 0.9 | Detection probability per target per platform |
| `lambda_c` | 1 | Expected clutter per frame per platform |
| `sigma_bearing` | 1° | Bearing measurement noise std |
| `sigma_a` | 0.5 m/s² | Process noise (acceleration std) |
| `feature_gamma` | 5 | Spectrum similarity sharpening; higher = more feature reliance |
| `max_misses` | 3 | Frames before track deletion |
| `confirm_hits` | 3 | Hits before track confirmed |
| `max_new_tracks` | 3 | Max new tracks per frame |

### Spectrum Feature System

`generate_ship_spectrum.m` operates in two modes:
- **Simulation** (`wav_file=''`): Synthesizes LOFAR (random line spectra + Gaussian broadening, 10–500Hz) and DEMON (shaft+blade harmonics, 1–50Hz) → 96D vector.
- **Wav import** (valid `.wav` path): STFT → LOFAR extraction; bandpass → Hilbert envelope → FFT → DEMON extraction. Falls back to simulation on read failure.

New tracks created from triangulation require both platform measurements' spectra to have Pearson ρ ≥ 0.4 (consistency check to suppress ghost tracks from clutter pairs).
