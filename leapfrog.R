# This is a basic leapfrog integrator. It integrates the Hamiltonian
# H(q, p) = 0.5 * p^2 + V(q), for the state (q) and momentum (p)
#
# dt : timestep
# T : number of timesteps to integrate for
# gradV : gradient of V with respect to q
#
# returns : q, p, and time in a list
leapfrog = function(dt, T, gradH) {
  p = rep(0, T)
  q = rep(0, T)
  time = rep(0, T)
  
  p[1] = 0.0
  q[1] = 1.0
  
  for(t in 1:(T - 1)) {
    ph = p[t] - (dt / 2.0) * gradV(q[t])
    q[t + 1] = q[t] + dt * ph
    p[t + 1] = ph - (dt / 2.0) * gradV(q[t + 1])
    time[t + 1] = t * dt
  }
  
  list(p = p, q = q, time = time)
}