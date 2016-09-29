function [torque,tau,theta,salt] = loop_equilibrium(w,R,nl,Zf,foldtrue,lambda,mu,ref_t,xi_t,F_s)
% [torque,tau,theta,salt] = loop_equilibrium(w,R,nl,Zf,foldtrue,lambda,mu,ref_t,xi_t,F_s)
%
% Compute equilibrium buoyancy torque, required torque tau=w-torque, 
% equilibrium temperature and salinity distributions
%  in a discrete loop (nl discrete levels)
% The equation of state is nonlinear (set by lambda and mu parameters)
% The loop may be folded at level Zf
% Note that the applied torque that would be required for a steady-state is tau=w-torque
%
% loop toolbox, distributed on GitHub: http://github.com/fabien-roquet/loop
% F. Roquet 2016
% GNU General Public License

% loop parameters
dl       = 2*pi/nl;
l        = (1:nl)'*dl;

% compute position of sink and source
z        = cos(l);
n_sink   = 0;
while ( z(n_sink+1) - Zf ) >= 1e-10 && n_sink<nl,
    n_sink = n_sink+1;
end
n_source = nl - n_sink;

% compute mask
mask=ones(nl,1);
if foldtrue == 1
    mask(1:n_sink)=0;
    mask(n_source:nl)=0;
end

% compute geometrical variables
z        = z.*mask + (1-mask).*z(n_sink);
curv     = sin(l) .* mask;

% compute analytically steady-state t/s distrib as a function of eccentric anomaly
theta  = loop_tracer_relax(w,R,nl,n_sink,ref_t,xi_t);
salt   = loop_tracer_fixed(w,R,nl,n_sink,F_s);
sigma  = -theta.*(1+lambda/2*theta-mu*z)+salt;
torque = sum(sigma.*curv,1)'/nl;
tau    = w - torque;


