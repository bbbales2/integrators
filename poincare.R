# This is the Poincare transformed Hamiltonian thing from page 237 of
# Leimkuhler & Reich. It integrates the Hamiltonian
# H(q, p) = g(q) * (0.5 * p^2 + V(q) - E), for the state (q) and momentum (p)
# with a timescaling function g(q) with the integrator on page 308 of Hairer,
# Lubich, and Wanner. 
#
# dtau : timestep
# T : number of timesteps to integrate for
# V : potential energy
# gradV : gradient of V with respect to q
# g : timescaling function (one argument)
# dg : derivative of timescaling function wrt argument
#
# returns : q, p, tau, and time in a list
poincare = function(dtau, T, V, gradV, g, dg) {
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

  E = 0.5 * p[1]^2 + V(q[1])
  
  for(t in 1:(T - 1)) {
    a = p[t] - (dtau / 2.0) * (g(q[t]) * gradV(q[t]) + dg(q[t]) * (V(q[t]) - E))
    binv = 1 / (-(dtau / 4.0) * dg(q[t]))
    #php = (binv + sqrt(binv^2 - 4 * a * binv)) / 2.0
    ph = (binv - sqrt(binv^2 - 4 * a * binv)) / 2.0
    ln = g(q[t])
    lnp = newt(q[t], ph)
    q[t + 1] = q[t] + (dtau / 2.0) * (ln + lnp) * ph
    p[t + 1] = ph - (dtau / 2.0) * (lnp * gradV(q[t + 1]) + dg(q[t + 1]) * (0.5 * ph^2 + V(q[t + 1]) - E))
    time[t + 1] = time[t] + (dtau / 2.0) * (ln + lnp)
    tau[t + 1] = (t - 1) * dtau
  }
  
  list(p = p, q = q, tau = tau, time = time)
}