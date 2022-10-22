
% calculate alpha m and beta m 
function [alpha_m, beta_m] = m_equations(V, Vrest)
u=Vrest-V;
alpha_m = (u+25) / (exp(2.5+.1*u)-1)/10;
beta_m = 4*exp(u/18);
end