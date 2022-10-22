% calculate alpha n and beta n 
function [alpha_n, beta_n] = n_equations(V, Vrest)
u=Vrest-V;
alpha_n = (.1 * u + 1)./(exp(1 + .1 * u) - 1) / 10;
beta_n = .125 * exp(u/80);
end