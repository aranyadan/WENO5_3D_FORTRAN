clear; %close all; clc;
global gamma
global R_gas_const
global Suth_const
global Re
global Pr

%% Parameters
CFL     = 0.75;	% CFL number
tFinal	= 0.80;	% Final time
nEx      = 50;  % Number of cells in x
nEy      = 50;  % Number of cells in y
gamma   = 1.4;  % Ratio of specific heats for ideal di-atomic gas
IC      = 3;	% 10 IC cases are available
plot_fig= 1;

R_gas_const = 287;
p_ref = 101325;             % Reference air pressure (N/m^2)
rho_ref= 1.225;             % Reference air density (kg/m^3)
T_ref = p_ref / (rho_ref * R_gas_const);
Cp = gamma * R_gas_const / (gamma-1);
Cv = Cp - gamma;
Re = 1000;
Suth_const = 110.4/T_ref;
Pr=0.73;

% Discretize spatial domain
a=0; b=1; dx=(b-a)/nEx; nx=nEx+1; x=linspace(a,b,nx);
c=0; d=1; dy=(d-c)/nEy; ny=nEy+1; y=linspace(c,d,ny);
[Y,X] = meshgrid(y,x);

% Xact = X;
% factor = 4.5;
% Yact = (exp(factor.*Y)-1)./(exp(factor)-1);
%
% dYdy = (exp(factor)-1)./(factor.*exp(factor.*Y));
% dYdy = repmat(dYdy,[1,1,4]);

dYdy=1;
Xact = X;
Yact = Y;
delta = abs(Yact(1,1) - Yact(1,2));

% Set IC
[rho0,u0,v0,p0,tFinal,CFL] = Euler_IC2d(Xact,Yact,IC);
E0 = p0./((gamma-1)*rho0)+0.5*(u0.^2+v0.^2);  % Total Energy density
a0 = sqrt(gamma*p0./rho0);            % Speed of sound
q0=reshape([rho0 rho0.*u0 rho0.*v0 rho0.*E0],nx,ny,4);        % vec. of conserved properties
q0 = set_boundary(q0);


% Discretize time domain
lambda0 = max(max(abs(u0) + abs(v0) +a0));
dt0=CFL*min([dx,dy,delta])/lambda0;  % using the system's largest eigenvalue

%% Solver Loop
% Load initial condition
q=q0; it=0; dt=dt0; t=0; lambda=lambda0;

%% Solver Loop
% figure;
id=0;
% h=createButton;
while t<tFinal

    % RK Initial step
    qo = q;

    % 1st stage
    [dFx,dFy] = WENOcaller(lambda,q,dx,dy);
    dFy = dFy .* dYdy;
    dF = dFx+dFy;

    q = qo-dt*dF;
    q = set_boundary(q);


    % 2nd Stage
    [dFx,dFy] = WENOcaller(lambda,q,dx,dy);
    dFy = dFy .* dYdy;
    dF = dFx+dFy;

    q = 0.75*qo+0.25*(q-dt*dF);
    q = set_boundary(q);


    % 3rd stage
    [dFx,dFy] = WENOcaller(lambda,q,dx,dy);
    dFy = dFy .* dYdy;
    dF = dFx+dFy;

    q = (qo+2*(q-dt*dF))/3;
    q = set_boundary(q);



    % compute primary properties
    rho=q(:,:,1); u=q(:,:,2)./rho; v=q(:,:,3)./rho; E=q(:,:,4)./rho; p=(gamma-1)*rho.*(E-0.5*(u.^2+v.^2));
    a=sqrt(gamma*p./rho); if min(min(p))<0; error('negative pressure found!'); end;

    % Update dt and time
    lambda = max(max(abs(u)+abs(v)+a));


    dt=CFL*min([dx,dy,delta])/lambda; if t+dt>tFinal; dt=tFinal-t; end

    % Update time and iteration counter
	t=t+dt; it=it+1;

    if(mod(it,100)==0)
        display(['t=',num2str(t),', iterations=',num2str(it)]);
    end

    cla;fdisplay(Xact,Yact,rho);title(['t=' num2str(t)]); id=id+1;
    pause(0.001)
%     print(['./plots/',num2str(id)],'-dpng')

end

% Calculation of flow parameters
% a = sqrt(gamma*p./rho); M = u./a; % Mach number [-]
% s = 1/(gamma-1)*(log(p/p_ref)+gamma*log(rho_ref./rho));
%                             % Entropy w.r.t reference condition
% ss = log(p./rho.^gamma);    % Dimensionless Entropy
% Q = rho.*u;                 % Mass Flow rate per unit area
% e = p./((gamma-1)*rho);     % internal Energy
