%% Discrete Mechanics and Optimal Control -- Simple Pendulum
function pend_dmoc( duration, N, targetangle )

% Discrete Mechanics and Optimal Control (DMOC) approach is used to 
% swing up a simple pendulum from it's rest position to a desired target angle. 
% 
% Author: Mahan Fathi 
%
% SQP Package: IPOPT
% Download IPOPT and many other cool SQP solvers at: 
% <https://www.inverseproblem.co.nz/OPTI/>
%
% Notes:
% 1. Feel free to apply DMOC by using this as a template. 
% 2. SQP methods are used, hence the solution is subjected to local minima.
% 3. The algorithm is lighning fast and solution is reached quickly, 
%    although running the animation may be a little contly. 
%
% Inputs: 
%   1. duration    : Time to get the target angle.
%   2. N           : Numebr of nodes, excluding starting point.
%   3. targetangle : Desired target angle. Pendulum gets there with zero
%                 velocity.
% Outputs:
%   1. An animation of pendulum performing the swing-up is auto-played and
%      a figure containing control policy is auto-dislayed.
%   2. Optimum argument, being states and control actions over time, is
%      saved as a .mat file in project's directory.
%

tic
global model;

% Initializaiton
close all
MaxIterations = 1000;
model.N = N; % Numebr of nodes, excluding starting point
model.Nd = N+1; % Number of total nodes, including starting point
model.targetangle = targetangle;
model.h = duration/N;
model.times = model.h*(0:N)';
model.Ncon = (N-1) + 2 + 2;     % (N-1) dmoc + 2 initial AND final position
                                % /velocity condition.

% Model Properties
model.l = 1;
model.m = 1;
model.g = 9.81;

% Initial Guess / Indices
xg = zeros(model.Nd,1);
ug = zeros(model.Nd,1);
X0 = [ xg; ug ];
model.ix = (1:model.Nd);
model.iu = model.Nd + (1:model.Nd);
model.NX = size(X0,1);

    funcs.objective         = @objfun;
		funcs.gradient          = @objgrad;
		funcs.constraints       = @confun;
		funcs.jacobian          = @conjac;
		funcs.jacobianstructure = @conjacstructure;
%       funcs.hessian           = @objhess;
		options.cl              = zeros(model.Ncon,1);
		options.cu              = zeros(model.Ncon,1);
		options.ipopt.max_iter  = MaxIterations;
		options.ipopt.hessian_approximation = 'limited-memory';


        [X, info] = ipopt(X0,funcs,options);

        filename = 'result.mat';
        save(filename,'X')
        
state = X(model.ix);
input = X(model.iu);
figure(1); clf; title('claimed trajectory');
plot(model.times,state)
x = [0 0]';
dt = 0.0005;
count = duration/dt;
u = interp1( model.times, input, 0:dt:duration );
for i=1:count

    draw((i-1)*dt,x);
    x=x+dt*dynamics(x,u(i));

end

end

% ====================================
%         Constraint Function
% ====================================
function c = confun(X)

global model;
h = model.h;
N = model.N;
Ncon = model.Ncon;
l = model.l;
m = model.m;
g = model.g;
targetangle = model.targetangle;

% quick index finder
q = @(i) i+1;
f = @(i) model.Nd+i+1;

% Allocate space to c
c = zeros( Ncon,1 );

for k=1:N-1

    c(k) =  h^2 * ( X(f(k-1))+2*X(f(k))+X(f(k+1)) ) - ...
            2*m*(  g*l*h^2 *( sin((X(q(k-1))+X(q(k)))/2)+sin((X(q(k+1))+X(q(k)))/2) ) + ...
            2*X(q(k-1)) - 4*X(q(k)) + 2*X(q(k+1)) );

end

    c(k+1) = X(q(0));
    c(k+2) = X(q(N)) - targetangle ;
    c(k+3) = h*(( X(f(0))+X(f(1)) )/4 - 0.5*g*l*m*sin( (X(q(0))+X(q(1)))/2 ) - ...
               m*( X(q(1))-X(q(0)) )/h^2);
    c(k+4) = h*(( X(f(N-1))+X(f(N)) )/4 - 0.5*g*l*m*sin( (X(q(N-1))+X(q(N)))/2 ) + ...
               m*( X(q(N))-X(q(N-1)) )/h^2);

end


% ====================================
%         Jacobian Function
% ====================================
function J = conjac(X)

global model;
h = model.h;
N = model.N;
NX = model.NX;
Ncon = model.Ncon;
l = model.l;
m = model.m;
g = model.g;

% quick index finder
q = @(i) i+1;
f = @(i) model.Nd+i+1;

% Allocate Sparce Space to J
J = spalloc( Ncon, NX, (N-1)*6+8 );

for k=1:N-1

    qpk = X( q(k-1) );
    qk  = X( q(k) );
    qnk = X( q(k+1) );
%     fpk = X( f(k-1) );
%     fk  = X( f(k) );
%     fnk = X( f(k+1) );
    J( k, [ q(k-1) q(k) q(k+1) f(k-1) f(k) f(k+1) ] ) = ...
        [(-1).*m.*(4+g.*h.^2.*l.*cos((1/2).*(qk+qpk))),(-1).*m.*(( ...
        -8)+g.*h.^2.*l.*(cos((1/2).*(qk+qnk))+cos((1/2).*(qk+qpk)))) ...
        ,(-1).*m.*(4+g.*h.^2.*l.*cos((1/2).*(qk+qnk))),h.^2,2.*h.^2, ...
        h.^2];

end

    J( k+1 , q(0) ) = 1;
    J( k+2 , q(N) ) = 1;
    J( k+3 , [ q(0) q(1) f(0) f(1) ] ) = ...
        [h*(m/h^2-1/4*g*l*m*cos( (X(q(0))+X(q(1)))/2 )),...
        h*(-m/h^2-1/4*g*l*m*cos( (X(q(0))+X(q(1)))/2 )),h/4,h/4];
    J( k+4 , [ q(N-1) q(N) f(N-1) f(N) ] ) = ...
        [h*(-m/h^2-1/4*g*l*m*cos( (X(q(N-1))+X(q(N)))/2 )),...
        h*(m/h^2-1/4*g*l*m*cos( (X(q(N-1))+X(q(N)))/2 )),h/4,h/4];

end

% ====================================
%         Jacobian Structure
% ====================================
function J = conjacstructure(X)

global model
Ncon = model.Ncon;
NX   = model.NX;
N    = model.N;

% quick index finder
q = @(i) i+1;
f = @(i) model.Nd+i+1;

% Allocate Sparce Space to J
J = spalloc( Ncon, NX, (N-1)*6+8 );

for k=1:N-1

    J( k, [ q(k-1) q(k) q(k+1) f(k-1) f(k) f(k+1) ] ) = ones(1,6);

end

    J( k+1 , q(0) ) = 1;
    J( k+2 , q(N) ) = 1;
    J( k+3 , [ q(0) q(1) f(0) f(1) ] ) = ones(1,4);
    J( k+4 , [ q(N-1) q(N) f(N-1) f(N) ] ) = ones(1,4);


end
% ====================================
%         Objective Function
% ====================================
function F = objfun(X)

    global model
    iu = model.iu;
	F  = sum(X(iu).^2);

end

% ====================================
%         Gradient Function
% ====================================
function G = objgrad(X)

    global model
    NX = model.NX;
    iu = model.iu;
    h  = model.h;
    G  = zeros(NX,1);
    G(iu) = 2 * X(iu);

end

% ====================================
%           Draw Function
% ====================================
function status = draw(t,x)
persistent hFig base a1 ac1 raarm;

if (isempty(hFig))
    hFig = figure(25);
    set(hFig,'DoubleBuffer', 'on');

    a1 = 0.75;  ac1 = 0.415;
    av = pi*[0:.05:1];
    rb = .03; hb=.07;
    aw = .01;
    base = rb*[1 cos(av) -1 1; -hb/rb sin(av) -hb/rb -hb/rb]';
    arm = [aw*cos(av-pi/2) -a1+aw*cos(av+pi/2)
        aw*sin(av-pi/2) aw*sin(av+pi/2)]';
    raarm = [(arm(:,1).^2+arm(:,2).^2).^.5, atan2(arm(:,2),arm(:,1))];
end

figure(hFig); cla; hold on; view(0,90);
patch(base(:,1), base(:,2),1+0*base(:,1),'b','FaceColor',[.3 .6 .4])
patch(raarm(:,1).*sin(raarm(:,2)+x(1)-pi),...
    -raarm(:,1).*cos(raarm(:,2)+x(1)-pi), ...
    0*raarm(:,1),'r','FaceColor',[.9 .1 0])
plot3(ac1*sin(x(1)), -ac1*cos(x(1)),1, 'ko',...
    'MarkerSize',10,'MarkerFaceColor','b')
plot3(0,0,1.5,'k.')
title(['t = ', num2str(t(1),'%.2f') ' sec']);
set(gca,'XTick',[],'YTick',[])

axis image; axis([-1.0 1.0 -1.0 1.0]);
drawnow;

status = 0;
end

% ====================================
%               Dynamics
% ====================================
function xdot = dynamics(x,u)

% pendulum parameters
global model;
l = model.l;
m = model.m;
g = model.g;
I = m*l^2;

xdot = [x(2,1); u-m*g*l*sin(x(1,1))./I];

end
