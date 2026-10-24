# gep-ana

Skeleton for the SBS GEp elastic ep analysis, modeled on
[gmn-ana](https://github.com/provakar1994/gmn-ana) (P. Datta): same layout of
namespaces in `include/` + `src/`, one umbrella header, JSON-configured macros,
a global cut evaluated with `TTreeFormula`, per-event cut flags, and an output
tree `Tout` that keeps every event passing the global cut.

The polarization observables (`src/Polarization.cpp`) are **stubs** to be
implemented. Everything up to them (selection, kinematics, FPP angles,
helicity bookkeeping) runs.

## Layout

```
include/gep-ana.h        umbrella header: include this in a macro
include/Constants.h      masses, mu_p
include/JSONManager.h    config reader (needs nlohmann/json.hpp, shipped with ROOT >= 6.24)
include/GEpConfig.h      kinematic setting read from the config (Ebeam, arm angles, ECAL distance)
src/SetROOTVar.cpp       branch helpers; a missing branch warns instead of crashing
src/KinematicVar.cpp     elastic ep from the proton angle, projection to ECAL, coplanarity
src/Cut.cpp              window, elliptical, area/fiducial cuts
src/FPP.cpp              FPP scattering angles (theta, phi) and closest approach (sclose, zclose)
src/Polarization.cpp     phi_fpp histograms per helicity/theta bin + polarization stubs
src/Utilities.cpp        histogram helpers
scripts/gep_elastic_ana.cpp        main macro
config/example/conf_gep_elastic_ana.json   example config (placeholder numbers)
```

## Running

```
cd gep-ana
mkdir -p out
root -l -b -q 'scripts/gep_elastic_ana.cpp("config/gep5/conf_gep5_5969.json","out/gep5")'
```
Set `"max_files_per_run": 5` in the config for a quick test, `-1` for all segments.
`config/gep5/conf_gep5_5969.json` is set up for the GEP5 run 5969 full replay;
`config/example/` is the generic template.
Writes `out/test_<kinematic>.root` with histograms, canvases `c1` (elastic
selection) and `c2` (FPP), the per-helicity `phi_fpp` histograms, and `Tout`.

## What the macro does (per event)

1. **Global cut** from the config string (`TTreeFormula`, leaves refreshed when the chain changes file).
2. **Elastic kinematics** read from the replay's `heep` module: `dpp`, ECAL residuals
   `dxECAL/dyECAL`, coplanarity `dphi`, Q2, epsilon, E', theta_e (requires `heep.datavalid==1`).
3. **Expected e' position** at ECAL = cluster position - residual (for the fiducial cut).
4. **Cut flags**: `dppCut`, `ecalCut` (ellipse in dx-dy), `dphiCut`, `fiduCut`
   (expected e' inside ECAL minus a margin), `elasticCut = dpp && ecal && fidu`.
5. **FPP**: `sbs.gemFPP.track.{theta,phi,sclose,zclose}` from the replay (count =
   `Ndata.sbs.gemFPP.track.theta`, filled only when a front track exists); the FPP track
   with smallest sclose inside the theta/zclose windows is used; `fppCut` on theta, zclose, sclose.
6. **Helicity**: `scalhel.hel` x sign (0 if `scalhel.errcode != 0`), per-run override `helicity_sign_by_run` for IHWP.
7. Fill histograms (dpp with ECAL cut / anti-cut, like gmn-ana's W cut/anti-cut)
   and `Tout` (all kinematics, target and focal-plane variables, FPP angles, helicity, flags).

After the loop the polarization stage is called (route A binned, route B
event-by-event). Nothing is split per run; runs are chained per kinematic.
Beam energy is the nominal value from the config.

## To implement (`src/Polarization.cpp`)

- `HelicityAsymmetry` - A(phi) = (N+ - N-)/(N+ + N-) per theta_fpp bin
- `AnalyzingPower` - A_y(theta_fpp, p) parametrization
- `FitAsymmetry` - fit A(phi) -> focal-plane P_x, P_y
- `SpinTransport` - focal plane -> P_t, P_l at the target
- `FitMaxLikelihood` - event-by-event alternative using `Tout`

Already implemented: `FormFactorRatio` (mu_p GE/GM = -mu_p (Pt/Pl)(E+E')/(2M) tan(theta_e/2))
and `SpinPrecessionAngle` (dipole approximation, chi = gamma kappa_p theta_bend).

## Check before using on real data

- **Branch names** are taken from the GEP5 run 5969 full replay (heep.*, sbs.tr.*, sbs.gemFPP.track.*, earm.ecal.*, scalhel.*).
  The FPP angles/sclose/zclose come from the replay; `src/FPP.cpp` is kept only as a cross-check tool.
- **phi_fpp convention** (zero and sign) must match the analyzing-power / asymmetry
  formula you use.
- **ECAL coordinates**: assumed x vertical (down), y horizontal, origin on the
  central ray at `ecal_dist`; adjust offsets if the replay uses another origin.
- All numbers in the example config are placeholders.
