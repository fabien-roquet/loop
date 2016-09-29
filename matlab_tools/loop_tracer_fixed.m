function [tracer]=tracer_fixed(w, R, nphi, n_sink, flux)
% [tracer]=tracer_fixed(w,R,nphi,n_sink,flux)
% 
% compute equilibrium distribution of a tracer in a discrete loop (nphi discrete levels)
% forced by a fixed flux added and removed on both horizontal sides of the loop: 
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

tracer=zeros(nphi,length(w));
for i=1:length(w)
    
    % small velocity case
    if abs(w(i))<1e-5, w(i)=1e-5; end
    
    % Solution for positive velocities
    if w(i) >= 1e-10
        
        delt=R/w(i);
        phi_1=find(phi_sink<=phi & phi<phi_source);
        phi_2=find(phi<phi_sink);
        phi_3=find(phi_source<=phi);
        C=(1-exp((phi_diff-2*pi)/delt))/(1-exp(-(phi_diff)/delt));
        tracer(phi_1,i)=flux*(phi_diff-2*pi)/w(i) ...
            +C*(flux*2*pi/w(i))*(exp((phi(phi_1)-phi_source)/delt)/(1+C*exp(-phi_diff/delt)));
        
        tracer(phi_2,i)=flux*(phi_diff)/w(i) ...
            -(flux*2*pi/w(i))*(exp((phi(phi_2)-phi_sink)/delt)/(1+C*exp(-phi_diff/delt)));
        
        tracer(phi_3,i)=flux*(phi_diff)/w(i) ...
            -(flux*2*pi/w(i))*(exp((phi(phi_3)-phi_sink-2*pi)/delt)/(1+C*exp(-phi_diff/delt)));
        
        % Solution for negative velocities
    elseif w(i) < -1e-10
        
        delt=R/w(i);
        phi_1=find(phi_sink<=phi & phi<phi_source);
        phi_2=find(phi<phi_sink);
        phi_3=find(phi_source<=phi);
        C=(1-exp(-(phi_diff-2*pi)/delt))/(1-exp((phi_diff)/delt));
        tracer(phi_1,i)=flux*(phi_diff-2*pi)/w(i) ...
            +C*(flux*2*pi/w(i))*(exp((phi(phi_1)-phi_sink)/delt)/(1+C*exp(phi_diff/delt)));
        
        tracer(phi_2,i)=flux*(phi_diff)/w(i) ...
            -(flux*2*pi/w(i))*(exp((phi(phi_2)-phi_source+2*pi)/delt)/(1+C*exp(phi_diff/delt)));
        
        tracer(phi_3,i)=flux*(phi_diff)/w(i) ...
            -(flux*2*pi/w(i))*(exp((phi(phi_3)-phi_source)/delt)/(1+C*exp(phi_diff/delt)));
        
    end
end

