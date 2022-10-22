% calculate alpha h and beta h
function [alpha_h, beta_h] = h_equations(V, Vrest)
u=Vrest-V;
alpha_h = .07 * exp(u/20);
beta_h = 1 ./ (1+exp(3 + .1*u));
end