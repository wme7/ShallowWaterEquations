function [res] = MUSCL_SWEres1d(q,smax,N,dx,limiter,fluxMethod)
%   MUSCL Monotonic Upstreat Centered Scheme for Conservation Laws
%   Van Leer's MUSCL reconstruction scheme using piece wise linear
%   reconstruction
%  
%   where: limiter='MC'; fluxMethod='HLL';
%
%   Flux at j+1/2
% 
%     j+1/2         Cell's grid:
%   | wL|   |
%   |  /|wR |           1   2   3   4        N-2 N-1  N
%   | / |\  |   {x=0} |-o-|-o-|-o-|-o-| ... |-o-|-o-|-o-| {x=L}
%   |/  | \ |         1   2   3   4   5        N-1  N  N+1
%   |   |  \|
%   |   |   |       NC: Here cells 1 and N are ghost cells%
%     j  j+1
%
% Written by Manuel Diaz, NTU, 2015.12.17.

    % Compute and limit slopes
    res=zeros(2,N); dq=zeros(2,N); 
    flux=zeros(2,N-1); qL=zeros(2,N-1); qR=zeros(2,N-1);
    for i = 1:2
        for j = 2:N-1 % for all internal faces
            switch limiter
                case 'MC' % MC limiter
                    % Find dq_j = minmod{fwd diff, bwd diff, cntrl diff}
                    dqR = 2*(q(i,j+1) - q(i,j));%/dx;
                    dqL = 2*(q(i,j) - q(i,j-1));%/dx;
                    dqC = (q(i,j+1) - q(i,j-1))/2;%/(2*dx);
                    dq(i,j) = minmod([dqR,dqL,dqC]);
                case 'MM' % Minmod limiter
                    % Find dq_j = minmod{fwd diff, bwd diff}
                    dqR = (q(i,j+1) - q(i,j));%/dx;
                    dqL = (q(i,j) - q(i,j-1));%/dx;
                    dq(i,j) = minmod([dqR,dqL]);
                case 'VA' % Van Albada limiter
                    dqR = (q(i,j+1) - q(i,j));%/dx;
                    dqL = (q(i,j) - q(i,j-1));%/dx;
                    dq(i,j) = vanalbada(dqR,dqL,dx);
            end
        end
    end

    % Left and Right extrapolated q-values at the boundary j+1/2
    for j = 2:N-2 % for the domain cells
        qL(:,j) = q(:, j ) + dq(:, j )/2;%*dx	% q_{j+1/2}^{-} from j
        qR(:,j) = q(:,j+1) - dq(:,j+1)/2;%*dx	% q_{j+1/2}^{+} from j+1
    end

    % Flux contribution to the residual of every cell
    for j = 2:N-2 % for all faces the domain cells
        % compute flux at j+1/2
        switch fluxMethod
            case 'LF' % Lax-Friedrichs
                flux(:,j) = LFflux(qL(:,j),qR(:,j),smax);
            case 'RUS' % Rusanov
                flux(:,j) = RUSflux(qL(:,j),qR(:,j));
            case 'HLL' % HLL
                flux(:,j) = HLLflux(qL(:,j),qR(:,j));
        end
        res(:, j ) = res(:, j ) + flux(:,j)/dx;
        res(:,j+1) = res(:,j+1) - flux(:,j)/dx;
    end

    % Flux contribution of the LEFT MOST FACE: left face of cell j=2.
    qR(:,1)=q(:,2)-dq(:,2)*dx/2;    qL(:,1) = qR(:,1);
    % compute: flux(:,1) = AUSMflux(qL(:,1),qR(:,1),gamma);
    switch fluxMethod
        case 'LF' % Lax-Friedrichs
            flux(:,1) = LFflux(qL(:,1),qR(:,1),smax);
        case 'RUS' % Rusanov
            flux(:,1) = RUSflux(qL(:,1),qR(:,1));
        case 'HLL' % HLL
            flux(:,1) = HLLflux(qL(:,1),qR(:,1));
    end
    res(:,2) = res(:,2) - flux(:,1)/dx;

    % Flux contribution of the RIGTH MOST FACE: right face of cell j=N-1.
    qL(:,N-1)=q(:,N-1)+dq(:,N-1)*dx/2;      qR(:,N-1) = qL(:,N-1);
    % compute: flux(:,N-1) = Xflux(qL(:,N-1),qR(:,N-1),gamma);
    switch fluxMethod
        case 'LF' % Lax-Friedrichs
            flux(:,N-1) = LFflux(qL(:,N-1),qR(:,N-1),smax);
        case 'RUS' % Rusanov
            flux(:,N-1) = RUSflux(qL(:,N-1),qR(:,N-1));
        case 'HLL' % HLL
            flux(:,N-1) = HLLflux(qL(:,N-1),qR(:,N-1));
    end
    res(:,N-1) = res(:,N-1) + flux(:,N-1)/dx;
end

function mm = minmod(v)
    % Using Harten's generalized definition
    % minmod: zero if opposite sign, otherwise the one of smaller magnitude.
    %m=size(v,1); mm=zeros(size(v,2),1); s=sum(sign(v),2)/m; ids=find(abs(s)==1);
    %if(~isempty(ids)); mm(ids)=s(ids).*min(abs(v(ids,:)),[],2); end
    s = sum(sign(v))/numel(v); 
    if abs(s)==1; mm = s*min(abs(v(:))); else, mm=0; end
end

function va = vanalbada(da,db,h)
    % Van Albada Slope Limiter Function
    % vanAlbada: extend the simetric formulation of the van leer limiter
    eps2=(0.3*h)^3; 
    va=0.5*(sign(da)*sign(db)+1)*((db^2+eps2)*da+(da^2+eps2)*db)/(da^2+db^2+2*eps2);
end

function FL = LFflux(qL,qR,smax)
    % Lax-Friedrichs flux:
    %
    % P. D. Lax, Weak Solutions of Nonlinear Hyperbolic Equations and Their
    % Numerical Computation, Commun. Pure and Applied Mathematics, 7, 159-193, 
    % 1954.
    %
    global g;
    
    % Left state
    hL = qL(1);
    uL = qL(2)./qL(1);
    
    % Right state
    hR = qR(1);
    uR = qR(2)./qR(1);
        
    % Left and Right fluxes
    FL = [hL*uL; hL*uL^2 + g*hL^2/2];
    FR = [hR*uR; hR*uR^2 + g*hR^2/2];
    
    % Lax-Friedrichs Numerical Flux
    FL = 0.5*( FR + FL + smax*(qL-qR) );
end

function Rusanov = RUSflux(qL,qR)
    % Rusanov flux 
    %
    global g;

    % Left state
    hL = qL(1);
    uL = qL(2)./qL(1);
    
    % Right state
    hR = qR(1);
    uR = qR(2)./qR(1);
        
    % Left and Right fluxes
    FL = [hL*uL; hL*uL^2 + g*hL^2/2];
    FR = [hR*uR; hR*uR^2 + g*hR^2/2];
    
    % Rusanov numerical flux
    %smax = max([abs(uL+sqrt(g*hL)),abs(uL-sqrt(g*hL)),... 
    %            abs(uR+sqrt(g*hR)),abs(uR-sqrt(g*hR))]); %Davis approach
	smax = max([abs(uL)+sqrt(g*hL),abs(uR)+sqrt(g*hR)]); %Toro's approach
    Rusanov = 0.5*( FR + FL + smax*(qL-qR) );
end

function HLL = HLLflux(qL,qR)
    % HLL flux:
    %
    global g;
    
    % Left state
    hL = qL(1);
    uL = qL(2)./qL(1);
    
    % Right state
    hR = qR(1);
    uR = qR(2)./qR(1);
    
    % Left and Right fluxes
    FL = [hL*uL; hL*uL^2 + g*hL^2/2];
    FR = [hR*uR; hR*uR^2 + g*hR^2/2];

    % Left and right speeds
    sL = min([ uR-sqrt(g*hR), uL-sqrt(g*hL) ]);
    sR = max([ uR+sqrt(g*hR), uL+sqrt(g*hL) ]);
        
    % HLL algorithm
    if 0<sL
        HLL = FL;
    elseif sR<0
        HLL = FR;
    else
        HLL = (sR*FL-sL*FR + sL*sR*(qR-qL))/(sR-sL);
    end
end
