# Optimal-Estimation
Application of Extended Kalman Filter for state estimation in re-entry vehicle dynamics

we consider the problem of entering a hyper-sonic vehicle into theearth atmosphere.

We assume that we have a dynamic model of the vehicle describingthe position, velocity and a term describing the drag exerted on the vehicle. So we have total of 5 state variables alongside with the set of observations got from radarshowing the the distance from the vehicle and the angle with which the radar observesthe vehicle. 

In this problem we make use of the given dynamic model in parallel tomeasurement samples to get the estimation of the states in the sense of sequential estimation algorithm of Extended Kalman filter(EKF).
