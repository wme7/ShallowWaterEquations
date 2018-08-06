%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Solving 1-D formulation of the Saint-Venant Equations 
%                [i.e. the shallow water equations (SWE)]
%                      with a 2nd-order MUSCL scheme
%
%                   dq/dt + df(q)/dx = 0, for x in [a,b] 
%
%           coded by Manuel A. Diaz, manuel.ade'at'gmail.com 
%            Institute of Applied Mechanics, NTU, 2015.11.20
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refs: 
% [1] Toro, Eleuterio F., and Eleuterio Toro. Shock-capturing methods for
%     free-surface shallow flows. New York: John Wiley, 2001. pp:128.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: This snippet is incomplete! I'll return to finish it later ;)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; %close all; clc;
global g

%% Parameters
CFL     = 0.5;	% CFL number
tFinal	= 0.1;	% Final time
nE      = 200;  % Number of cells/Elements
flowIC  = 08;	% see details in CommonIC.m
topoIC  = 00;   % see details in TopographiIC.m
limiter ='MC';	% MC, MM, VA.
fluxMth ='RUS';	% LF, RUS, HLL.
plot_fig= true;

% Gravity
g = 9.81; % [m/s]

% Discretize spatial domain
a=-1; b=1; dx=(b-a)/nE; nx=nE+1; xc=linspace(a,b,nx);

% Set IC
[h0,u0]=CommonIC(xc,flowIC); 
[b]=TopographyIC(xc,topoIC);

% Exact solution for Riemann Problems
if flowIC >= 4 && flowIC <= 8 && topoIC == 0
    [qe,xe]=ExactRiemannSWE_Toro2001(xc,tFinal,u0(1),u0(end),h0(1),h0(end));
    he=qe(1,:); ue=qe(2,:); 
end

% Set q-array & adjust grid for ghost cells
nx=nx+2; q0=[h0; h0.*u0]; zero=[0;0]; q0=[zero,q0,zero]; b=[b(1),b,b(end)];

% Boundary Conditions in ghost cells
q0(:,1)=q0(:,2); q0(:,nx)=q0(:,nx-1);   % Natural BCs

% initial time step
h=q0(1,:); u=q0(2,:)./h; lambda0=max([abs(u+sqrt(g*h)),abs(u-sqrt(g*h))]);
dt0=CFL*dx/lambda0;  % using the system's largest eigenvalue

%% Solver 

% Load initial condition
q=q0; it=0; dt=dt0; t=0; lambda=lambda0;

while t<tFinal
    % Compute primary properties
    h=q(1,:); u=q(2,:)./q(1,:); 
    if min(h)<0; error('negative water height found!'); end

    % Update dt
    lambda=max([abs(u+sqrt(g*h)),abs(u-sqrt(g*h))]);
    dt=CFL*dx/lambda; if t+dt>tFinal; dt=tFinal-t; end

    % RK Initial step
    qo = q;

    % 1st stage
    L=MUSCL_SWEres1d(q,lambda,nx,dx,limiter,fluxMth);	q=qo-dt*L;
    q(:,1)=q(:,2); q(:,nx)=q(:,nx-1); % Neumann BCs

    % 2nd Stage
    L=MUSCL_SWEres1d(q,lambda,nx,dx,limiter,fluxMth);	q=0.75*qo+0.25*(q-dt*L);
    q(:,1)=q(:,2); q(:,nx)=q(:,nx-1); % Neumann BCs

    % 3rd stage
    L=MUSCL_SWEres1d(q,lambda,nx,dx,limiter,fluxMth);	q=(qo+2*(q-dt*L))/3;
    q(:,1)=q(:,2); q(:,nx)=q(:,nx-1); % Neumann BCs

    % Update time and iteration counter
    t=t+dt; it=it+1;

    % Plot figure
    if plot_fig==true && rem(it,10)==0
        subplot(2,1,1); plot(xc,h0,xc,h(2:nx-1),'.r',xc,b(2:nx-1),'-k');
        subplot(2,1,2); plot(xc,u0,xc,u(2:nx-1),'.b');
        drawnow
    end
end

%% Final plot
h=q(1,:); u=q(2,:)./q(1,:);
if flowIC >= 4 && flowIC <= 8 && topoIC == 0
    s1=subplot(2,1,1); plot(xc,h0,'--k',xe,he,'-r',xc,h(2:nx-1),'or',xc,b(2:nx-1),'-k');
    s2=subplot(2,1,2); plot(xc,u0,'--k',xe,ue,'-b',xc,u(2:nx-1),'ob');
else
    s1=subplot(2,1,1); plot(xc,h0,'--k',xe,h(2:nx-1),'or',xc,b(2:nx-1),'-k');
    s2=subplot(2,1,2); plot(xc,u0,'--k',xe,u(2:nx-1),'ob');
end
title(s1,['MUSCL-',fluxMth,'-',limiter,' SWE Solver']);