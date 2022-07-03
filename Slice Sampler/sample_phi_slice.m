function [phih,uu] = sample_phi_slice(phih,h,muh,sig2h,a0,b0)
%% sample_phi_slice

phi_set = linspace(-1,1,100);    % set of possible values for phi

% target
% func_phi = @(phi) betapdf((phi+1)/2,a0,b0) .* exp( -sum( (h(2:end)-muh-phi.*(h(1:end-1)-muh)).^2 ) / 2*sig2h );
func_phi = @(phi) betapdf((phi+1)/2,a0,b0) .* sum( normpdf(h(2:end), muh+phi.*(h(1:end-1)-muh), sqrt(sig2h)) );

uu = unifrnd(0,func_phi(phih));     % auxiliary slice variable
Ap = (func_phi(phi_set) >= uu);     % set of values
phih = randsample(phi_set(Ap),1);

end