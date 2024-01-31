%% RBC model - Carbon Tax 
%% By Ginevra Vimercati Sanseverino, Matteo Meneghello, Barbara Holguin, Anthony Rahme, Pepe Garcia, Eduardo Pulgar - October 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Endogenous variables: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var 
    c       % Consumption  
    l       % Labour
    g       % Government Spending
    t       % Taxes/Transfers
    w       % Wages
    y       % Output 
    rF      % Risk-free rate 
    b       % Bonds supply
    lambda  % The lagrangian with respect to consumption and savings (bonds) - marginal utility of consumption
    e_a     % The shock variable within the law of motion to output
    em      % Emissions
    ctax    % Carbon tax
    miu     % Abatement
    div     % Dividends
;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exogenous variables: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varexo 
    eta_a     % The exogenous output shock
;                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters 
    beta
    alpha
    chi
    phi
    l0 
    b0 
    g0
    sdev_a 
    sig_a 
    rho_a 
    phi_em
    theta1
    theta2
    ctax_p 
;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calibration for all the parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta       = .9959;                 % Time preference 
alpha      = 1/3;                   % Share of capital income
chi        = 12;                    % Level of disutility (such that l=1/3 at the Steady state)
phi        = 1;                     % Frisch elasticity of labour disutility
l0         = 1/3;                   % Hours worked (1/3 of 24h a day)
b0         = 0;                     % Fixed supply of bonds (Net-zero supply)
g0         = .3;                    % Gov spending as % of GDP (30% here)
sdev_a     = .1;                    % Standard deviation of the exogenous shock to income
rho_a      = .9;                    % The shock persistence level
sig_a      = 1;                     % The sign of the shock (1 if positive, -1 is negative)
phi_em     = .3;                    % Emission intensity
theta1     = .05;                   % Abatement cost parameter 1 ∈ [0.05,0.1]
theta2     = 2.7;                   % Abatement cost parameter 2
ctax_p     = .02;                   % Carbon tax ∈ [0.02,0.1]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model FOCs and Equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model;  
%%%% Households Side %%%%
[name='FOC w.r.t to b: Euler Equation']
lambda = rF*beta*lambda(+1);

[name='FOC w.r.t to c: Marginal Utility of Consumption']
c = 1/lambda;

[name='FOC w.r.t to l: Disutility of labour']
w = chi*l^phi/lambda;

[name='Consumer Budget Constraint']
b = rF*b(-1) + w*l + div - t - c; 
% y = c + i + g; % (this is the same equation)

%%%% Firms Side %%%%
[name='Production Function']
y = e_a*l^(1-alpha);

[name='Emission Abatement']
em = (1-miu)*phi*y;

[name='FOC w.r.t l']
w = (1-alpha)*e_a*l^(-alpha);

[name='FOC w.r.t miu']
ctax/phi_em = theta1*theta2*(miu^(theta2-1));

[name='Dividends']
div = y - w*l - theta1*miu^theta2 - ctax*em;

%%%% Government Spending %%%%
[name='Gov Spending']
g = g0*y;

[name='Gov Budget Constraint']
g = t;

%%%% Carbon Tax %%%%
[name='Carbon Tax']
ctax = ctax_p;

%%%% Market Clearing conditions (to close the model) %%%%
[name='Supply of bonds']
b = b0;

%%%% Exogenous shock law of motion %%%%
[name='Law of motion of the output shock']
log(e_a) = rho_a*log(e_a(-1)) +sig_a*eta_a;

end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Steady state 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
steady_state_model;
%% Analytically without a solver
e_a = 1;                                             % Fixing the shock
ctax = ctax_p;                                       % Carbon Tax
b = b0;                                              % Supply of Bonds
l = l0;                                              % Supply of Labor
rF = 1/beta;                                         % Euler Equation
miu = (ctax/(phi_em*theta1*theta2))^(1/(theta2-1));  % FOC w.r.t miu
y = e_a*l^(1-alpha);                                 % Production Function
w = (1-alpha)*l^(-alpha);                            % FOC w.r.t l
em = (1-miu)*phi*y;                                  % Emission Abatement
div = y - w*l - theta1*miu^theta2 - ctax*em;         % Dividends
g = g0*y;                                            % Gov Spending
t = g;                                               % Gov Budget Constraint
c = b*(rF-1)+w*l+div-t;                              % Consumer Budget Constraint
lambda = 1/c;                                        % Marginal Utility of Consumption
chi = w*lambda/l^phi;                                % Disutility of labour
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conducting all the checks (Residual, Blanchard-Khan ...) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
steady;     % Provides the steady state values
resid(1);   % Residual
check;      % Blanchard-Khan conditions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Declaring the shocks 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shocks;
var eta_a; stderr sdev_a; 
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performing the simulation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stoch_simul (order=1, pruning, irf=30) c t y l lambda rF em miu ctax w; 
