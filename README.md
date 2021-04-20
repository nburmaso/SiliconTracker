## Simple silicon tracker

#### Description:
This is a model of a simple silicon tracker based on Geant4.

Current status:
* Ideal track finder is used
* A class for the Kalman filter method is implemented ``STTrackFitter``
    * The track model is following: {x, y, tx, ty, q/p}
    * Track parameters are propagated in magnetic field with RK4
    * At the moment, energy losses is calculated with Bethe-Bloch formula only
    * Track parameters and covariance matrix are corrected according to energy losses [WIP]

####Requirements:
* ROOT
* Geant4
