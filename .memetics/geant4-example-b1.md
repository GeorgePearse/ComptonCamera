---
upstream:
  - repo: Geant4/geant4
    paths:
      - examples/basic/B1/exampleB1.cc
      - examples/basic/B1/src/DetectorConstruction.cc
      - examples/basic/B1/src/PrimaryGeneratorAction.cc
      - examples/basic/B1/src/EventAction.cc
      - examples/basic/B1/src/RunAction.cc
      - examples/basic/B1/src/SteppingAction.cc
      - examples/basic/B1/src/ActionInitialization.cc
    ref: master
    last_seen: b4a16de652ec244f7a0ecc0a5cacfee21930bf75
local:
  paths:
    - exampleB1.cc
    - src/B1DetectorConstruction.cc
    - src/B1PrimaryGeneratorAction.cc
    - src/B1EventAction.cc
    - src/B1RunAction.cc
    - src/B1SteppingAction.cc
    - src/B1ActionInitialization.cc
---

# Geant4 basic example B1 as the simulation skeleton

The Compton camera simulation is the Geant4 `examples/basic/B1` application copied (Geant4 10.x era, when
the classes were still prefixed `B1`) and extended: the main program, `DetectorConstruction`,
`PrimaryGeneratorAction`, `EventAction`, `RunAction`, `SteppingAction` and `ActionInitialization` keep
B1's structure and user-action wiring, while the physics content was replaced with cylindrical
scatterer/absorber detectors, an isotropic 662 keV cone source, coincidence-event file writers and
peak broadening. Copying the example is how Geant4 applications are expected to start; the upstream
files are worth watching because Geant4 periodically reworks the example (namespace `B1`, modern
`G4RunManagerFactory`, analysis manager usage) and those changes show what a maintained application
should look like.

Adaptation notes:

- Class names stay `B1*` (headers in `include/`); do not rename to upstream's unprefixed, namespaced
  classes.
- Geometry is driven by `B1DetectorMessenger` (`/B1/scat`, `/B1/det`, 1:1 .. 3:1 setups) and the source
  by `B1PrimaryGeneratorMessenger`; these messengers are local additions with no upstream equivalent.
- Physics list is `EMStandard_Option4` (not B1's `QBBC`); keep it.
- Data output (coincidence writers in `B1EventAction` / `B1RunAction`, peak broadening flag, output
  folder command) is local and must survive any sync; B1's dose scoring is gone.
