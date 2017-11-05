# This is the semi-explicit Sundman transformed generated leapfrog from
# page 243 of Leimkuhler & Reich. It integrates the Hamiltonian
# H(q, p) = 0.5 * p^2 + V(q), for the state (q) and momentum (p)
# with a timescaling dtau/dt = g(q)
#
# dtau : timestep
# T : number of timesteps to integrate for
# gradV : gradient of V with respect to q
# g : timescaling function (one argument)
# dg : derivative of timescaling function wrt argument
#
# returns : q, p, tau, and time in a list
semi_sundman = function(dtau, T, gradV, g, dg) {
  newt = function(q, ph) {
    ln = g(q)
    x = ln
    
    for(n in 1:10) {
      a = q + (dtau / 2.0) * (ln + x) * ph
      f = x - g(a)
      fp = 1 - gp(a) * (dtau / 2.0) * ph
      x = x - f / fp
    }
    
    x
  }

  p = rep(0, T)
  q = rep(0, T)
  time = rep(0, T)
  tau = rep(0, T)
  
  p[1] = 0.0
  q[1] = 1.0
  
  for(t in 1:(T - 1)) {
    ph = p[t] - (dtau / 2.0) * g(q[t]) * gradV(q[t])
    ln = g(q[t])
    lnp = newt(q[t], ph)
    q[t + 1] = q[t] + (dtau / 2.0) * (ln + lnp) * ph
    p[t + 1] = ph - (dtau / 2.0) * lnp * gradV(q[t + 1])
    time[t + 1] = time[t] + g(q[t]) * dtau
    tau[t + 1] = (t - 1) * dtau
  }
  
  list(p = p, q = q, tau = tau, time = time)
}