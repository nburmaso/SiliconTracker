## Simple silicon tracker

#### Description:

This is a model of a simple silicon tracker based on Geant4.

Current status:

* Ideal track finder is used
* A class for the Kalman filter method is implemented ``STTrackFitter``
  * The track model is following: {x, y, tx, ty, q/p}
  * Track parameters are propagated in magnetic field with RK4
  * At the moment, energy losses are calculated with Bethe-Bloch formula only
  * Track parameters and covariance matrix are corrected according to energy losses [WIP]
* KF parameters of fitted tracks are stored into the `tracks.root`:

```
+-- reco_tracks
|   +-- z
|   +-- x
|   +-- y
|   +-- tx
|   +-- ty
|   +-- qp
```

* Initial parameters of MC tracks are stored into the `mc_tracks.root` file:

```
+-- mc_tracks
|   +-- mcEventID
|   +-- mcTrackID
|   +-- pdgID
|   +-- motherID
|   +-- vx
|   +-- vy
|   +-- vz
|   +-- px
|   +-- py
|   +-- pz
```

#### Macros

* `draw_pulls.cc` -- can be used to draw pulls and residuals for the fitted parameters. For now one needs to manually
  enable collection of these distributions by setting `fPulls = true` in the `STEventAction::fitTracksKF()`.

#### Requirements:

* ROOT
* Geant4
