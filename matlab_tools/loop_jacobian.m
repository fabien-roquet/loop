function [jacobian] = loop_jacobian(w,R,nl,Zf,foldtrue,lambda,mu,ref_t,xi_t,F_s)
% [jacobian] = loop_jacobian(w,R,nl,Zf,foldtrue,lambda,mu,ref_t,xi_t,F_s)
%
% Compute the jacobian associated with the steady-state defined by loop parameters
% Note that the steady-state requires an applied torque such as, tau=w-torque
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

% Compute Jacobian

% init matrix
J = zeros(2*nl);

%temperature
for nn = 1:nl, J(nn    ,nn ) =  - 2*R / dl^2; end
for nn = 2:nl, J(nn-1  ,nn ) = (   w/2 +   R/dl ) / dl ; end
J(nl    ,1       ) = (  w/2 +   R/dl )  / dl ;
for nn = 1:nl-1, J(nn+1  ,nn ) = ( - w/2 +   R/dl ) / dl ; end
J(1       ,nl   ) = ( - w/2 +   R/dl ) / dl ;
J(n_source  ,n_source  ) = J(n_source  ,n_source  ) - ref_t * xi_t * nl   ;
J(n_sink    ,n_sink    ) = J(n_sink    ,n_sink    ) - ref_t * xi_t * nl   ;
for nm = 1:nl,
    J(nm,1       ) = J(nm,1       ) + ...
        curv(nm) .* ( 1 + lambda * theta(nm) - mu * z(nm) ) / nl *...
        ( theta(2) - theta(nl  ) ) / (2 * dl) ;
    J(nm+nl,1  ) = J(nm+nl,1  ) - ...
        curv(nm) / nl * ( theta(2) - theta(nl  ) ) / (2 * dl) ;
    for nn = 2:nl-1,
        J(nm ,nn ) = J(nm ,nn ) + ...
            curv(nm) .* ( 1 + lambda * theta(nm) - mu * z(nm) ) / nl *...
            ( theta(nn+1) - theta(nn-1) ) / (2 * dl) ;
        J(nm+nl ,nn ) = J(nm+nl ,nn ) - ...
            curv(nm) / nl * ( theta(nn+1) - theta(nn-1) ) / (2 * dl) ;
    end
    J(nm,nl    ) = J(nm,nl    ) + ...
        curv(nm) .* ( 1 + lambda * theta(nm) - mu * z(nm) ) / nl *...
        ( theta(1) - theta(nl-1) ) / (2 * dl) ;
    J(nm+nl,nl) = J(nm+nl,nl) - ...
        curv(nm) / nl * ( theta(1) - theta(nl-1) ) / (2 * dl) ;
end

% salinity
for nn = nl+(1:nl), J(nn    ,nn ) =  - 2*R / dl^2; end
for nn = nl+(2:nl), J(nn-1  ,nn ) = (   w/2 +   R/dl ) / dl ; end
J(2*nl  ,nl+1  ) = (  w/2 +   R/dl )  / dl ;
for nn = nl+(1:nl-1), J(nn+1  ,nn ) = ( - w/2 +   R/dl ) / dl ; end
J(nl+1   ,2*nl ) = ( - w/2 +   R/dl ) / dl ;
for nm = 1:nl,
    J(nm,1+nl   ) = J(nm,1+nl   ) + ...
        curv(nm) .* ( 1 + lambda * theta(nm) - mu * z(nm) ) / nl *...
        ( salt(2) - salt(nl  ) ) / (2 * dl) ;
    J(nm+nl,1+nl) = J(nm+nl,1+nl) - ...
        curv(nm) / nl * ( salt(2) - salt(nl  ) ) / (2 * dl) ;
    for nn = 2:nl-1,
        J(nm ,nn+nl ) = J(nm ,nn+nl ) + ...
            curv(nm) .* ( 1 + lambda * theta(nm) - mu * z(nm) ) / nl *...
            ( salt(nn+1) - salt(nn-1) ) / (2 * dl) ;
        J(nm+nl ,nn+nl ) = J(nm+nl ,nn+nl ) - ...
            curv(nm) / nl * ( salt(nn+1) - salt(nn-1) ) / (2 * dl) ;
    end
    J(nm,nl+nl    ) = J(nm,nl+nl    ) + ...
        curv(nm) .* ( 1 + lambda * theta(nm) - mu * z(nm) ) / nl *...
        ( salt(1) - salt(nl-1) ) / (2 * dl) ;
    J(nm+nl,nl+nl) = J(nm+nl,nl+nl) - ...
        curv(nm) / nl * ( salt(1) - salt(nl-1) ) / (2 * dl) ;
end
jacobian = J ;


