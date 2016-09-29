function [tracer]=tracer_relax(w, R, nphi, n_sink, tracer_ref, xi)
% [tracer]=tracer_relax(w,R,nphi,n_sink,tracer_ref,xi)
% 
% compute equilibrium distribution of a tracer in a discrete loop (nphi discrete levels)
% forced by a relaxation flux added and removed on both horizontal sides of the loop: 
%  sink applied at   n = n_sink
%  source applied at n = nphi-n_sink
% advective/diffusive equation with loop velocity w and diffusivity R
%
% loop toolbox, distributed on GitHub: http://github.com/fabien-roquet/loop
% F. Roquet 2016
% GNU General Public License

dphi=2*pi/nphi;
phi=(1:nphi)*dphi;

n_source   = nphi - n_sink;
phi_sink   = phi(n_sink);
phi_source = phi(n_source);
phi_diff   = phi_source-phi_sink;

% avoid w=0 case
w(abs(w)<1e-10)=1e-10;

tracer=zeros(nphi,length(w));
for i=1:length(w)
    % Solution for positive velocities
    if w(i) >= 1e-10
        
        delt=R/w(i);
        phi_1=find(phi_sink<=phi & phi<phi_source);
        phi_2=find(phi<phi_sink);
        phi_3=find(phi_source<=phi);
        C=(1-exp((phi_diff-2*pi)/delt))/(1-exp(-(phi_diff)/delt));
        D=(1+exp((phi_diff-2*pi)/delt))/(1-exp((phi_diff-2*pi)/delt)+2*C*exp(-phi_diff/delt));
        t_inf_m=2*pi*(xi/w(i))*tracer_ref/(-(1+D)*(exp((phi_diff-2*pi)/delt)+C)/(1+C*exp(-phi_diff/delt)) ...
            -(2*pi*xi/w(i))*(D-exp((phi_diff-2*pi)/delt)*(1+D)/(1+C*exp(-phi_diff/delt))));
        t_inf_p=-t_inf_m*D;
        tracer(phi_1,i)=t_inf_m ...
            -(C*(t_inf_m-t_inf_p)/(1+C*exp(-phi_diff/delt)))*exp((phi(phi_1)-phi_source)/delt);
        
        tracer(phi_2,i)=t_inf_p ...
            +((t_inf_m-t_inf_p)/(1+C*exp(-phi_diff/delt)))*exp((phi(phi_2)-phi_sink)/delt);
        
        tracer(phi_3,i)=t_inf_p ...
            +((t_inf_m-t_inf_p)/(1+C*exp(-phi_diff/delt)))*exp((phi(phi_3)-phi_sink-2*pi)/delt);
        
        % Solution for negative velocities
    elseif w(i) < -1e-10
        
        
        delt=R/w(i);
        phi_1=find(phi_sink<=phi & phi<phi_source);
        phi_2=find(phi<phi_sink);
        phi_3=find(phi_source<=phi);
        C=(1-exp((2*pi-phi_diff)/delt))/(1-exp((phi_diff)/delt));
        D=C*(1+exp(phi_diff/delt))/(1+exp((2*pi-phi_diff)/delt));
        t_inf_m=2*pi*(xi/w(i))*tracer_ref/((1+D)*(1+C*exp(phi_diff/delt))/(1+C*exp(phi_diff/delt)) ...
            -(2*pi*xi/w(i))*(D-C*exp(phi_diff/delt)*(1+D)/(1+C*exp(phi_diff/delt))));
        t_inf_p=-t_inf_m*D;
        tracer(phi_1,i)=t_inf_p ...
            -(C*(t_inf_p-t_inf_m)/(1+C*exp(phi_diff/delt)))*exp((phi(phi_1)-phi_sink)/delt);
        
        tracer(phi_2,i)=t_inf_m ...
            +((t_inf_p-t_inf_m)/(1+C*exp(phi_diff/delt)))*exp((phi(phi_2)+2*pi-phi_source)/delt);
        
        tracer(phi_3,i)=t_inf_m ...
            +((t_inf_p-t_inf_m)/(1+C*exp(phi_diff/delt)))*exp((phi(phi_3)-phi_source)/delt);
        
    end
end


