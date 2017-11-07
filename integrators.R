library(tidyverse)
library(ggplot2)

V = function(q) { exp(2 * q) / 2.0 - q }
energy = function(q, p) { p^2 / 2.0 + V(q) }
gradV = function(q) { exp(2 * q) - 1.0 }

dt = 1e-1
T = 100

# Leapfrog
source("leapfrog.R")
lf = leapfrog(dt, T, gradV) %>% as.tibble %>%
  mutate(energy = energy(q, p))

lf %>% ggplot(aes(q, p)) + geom_path() + geom_point()
lf %>% ggplot(aes(time, energy)) + geom_path() + geom_point()

# Semi-explicit Sundman stuff
g = function(q) { 1 / (1.0 + exp(q)) }
dg = function(q) { -exp(q) / (1.0 + exp(q))^2 }

source("semi_sundman.R")
df = semi_sundman(dt, T, gradV, g, dg) %>% as.tibble %>%
  mutate(energy = energy(q, p), dtdtau = g(q))

df %>% ggplot(aes(q, p)) + geom_path() + geom_point()
df %>% ggplot(aes(tau, time)) + geom_path() + geom_point()
df %>% ggplot(aes(tau, dtdtau)) + geom_path() + geom_point()
df %>% ggplot(aes(time, energy)) + geom_path() + geom_point()

# Semi-explicit Sundman stuff
g = function(q) { 1 / (1.0 + exp(q)) }
dg = function(q) { -exp(q) / (1.0 + exp(q))^2 }

source("poincare.R")
pc = poincare(dt, T, V, gradV, g, dg) %>% as.tibble %>%
  mutate(energy = energy(q, p), dtdtau = g(q))

pc %>% ggplot(aes(q, p)) + geom_path() + geom_point()
pc %>% ggplot(aes(tau, time)) + geom_path() + geom_point()
pc %>% ggplot(aes(tau, dtdtau)) + geom_path() + geom_point()
pc %>% ggplot(aes(time, energy)) + geom_path() + geom_point()
