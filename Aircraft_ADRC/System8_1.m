function [sys,x0,str,ts] = System8_1(t,x,u,flag,pa)
%SFUNTMPL General M-file S-function template
%   With M-file S-functions, you can define you own ordinary differential
%   equations (ODEs), discrete system equations, and/or just about
%   any type of algorithm to be used within a Simulink block diagram.
%
%   The general form of an M-File S-function syntax is:
%       [SYS,X0,STR,TS] = SFUNC(T,X,U,FLAG,P1,...,Pn)
%
%   What is returned by SFUNC at a given point in time, T, depends on the
%   value of the FLAG, the current state vector, X, and the current
%   input vector, U.
%
%   FLAG   RESULT             DESCRIPTION
%   -----  ------             --------------------------------------------
%   0      [SIZES,X0,STR,TS]  Initialization, return system sizes in SYS,
%                             initial state in X0, state ordering strings
%                             in STR, and sample times in TS.
%   1      DX                 Return continuous state derivatives in SYS.
%   2      DS                 Update discrete states SYS = X(n+1)
%   3      Y                  Return outputs in SYS.
%   4      TNEXT              Return next time hit for variable step sample
%                             time in SYS.
%   5                         Reserved for future (root finding).
%   9      []                 Termination, perform any cleanup SYS=[].
%
%
%   The state vectors, X and X0 consists of continuous states followed
%   by discrete states.
%
%   Optional parameters, P1,...,Pn can be provided to the S-function and
%   used during any FLAG operation.
%
%   When SFUNC is called with FLAG = 0, the following information
%   should be returned:
%
%      SYS(1) = Number of continuous states.
%      SYS(2) = Number of discrete states.
%      SYS(3) = Number of outputs.
%      SYS(4) = Number of inputs.
%               Any of the first four elements in SYS can be specified
%               as -1 indicating that they are dynamically sized. The
%               actual length for all other flags will be equal to the
%               length of the input, U.
%      SYS(5) = Reserved for root finding. Must be zero.
%      SYS(6) = Direct feedthrough flag (1=yes, 0=no). The s-function
%               has direct feedthrough if U is used during the FLAG=3
%               call. Setting this to 0 is akin to making a promise that
%               U will not be used during FLAG=3. If you break the promise
%               then unpredictable results will occur.
%      SYS(7) = Number of sample times. This is the number of rows in TS.
%
%
%      X0     = Initial state conditions or [] if no states.
%
%      STR    = State ordering strings which is generally specified as [].
%
%      TS     = An m-by-2 matrix containing the sample time
%               (period, offset) information. Where m = number of sample
%               times. The ordering of the sample times must be:
%
%               TS = [0      0,      : Continuous sample time.
%                     0      1,      : Continuous, but fixed in minor step
%                                      sample time.
%                     PERIOD OFFSET, : Discrete sample time where
%                                      PERIOD > 0 & OFFSET < PERIOD.
%                     -2     0];     : Variable step discrete sample time
%                                      where FLAG=4 is used to get time of
%                                      next hit.
%
%               There can be more than one sample time providing
%               they are ordered such that they are monotonically
%               increasing. Only the needed sample times should be
%               specified in TS. When specifying than one
%               sample time, you must check for sample hits explicitly by
%               seeing if
%                  abs(round((T-OFFSET)/PERIOD) - (T-OFFSET)/PERIOD)
%               is within a specified tolerance, generally 1e-8. This
%               tolerance is dependent upon your model's sampling times
%               and simulation time.
%
%               You can also specify that the sample time of the S-function
%               is inherited from the driving block. For functions which
%               change during minor steps, this is done by
%               specifying SYS(7) = 1 and TS = [-1 0]. For functions which
%               are held during minor steps, this is done by specifying
%               SYS(7) = 1 and TS = [-1 1].

%   Copyright 1990-2002 The MathWorks, Inc.
%   $Revision: 1.18 $

%
% The following outlines the general structure of an S-function.
%
switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    [sys,x0,str,ts]=mdlInitializeSizes(pa);

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
  case 1,
    sys=mdlDerivatives(t,x,u,pa);

  %%%%%%%%%%
  % Update %
  %%%%%%%%%%
  case 2,
    sys=mdlUpdate(t,x,u);

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3,
    sys=mdlOutputs(t,x,u,pa);

  %%%%%%%%%%%%%%%%%%%%%%%
  % GetTimeOfNextVarHit %
  %%%%%%%%%%%%%%%%%%%%%%%
  case 4,
    sys=mdlGetTimeOfNextVarHit(t,x,u);

  %%%%%%%%%%%%%
  % Terminate %
  %%%%%%%%%%%%%
  case 9,
    sys=mdlTerminate(t,x,u);

  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
  otherwise
    error(['Unhandled flag = ',num2str(flag)]);

end

% end sfuntmpl

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts]=mdlInitializeSizes(pa)

%
% call simsizes for a sizes structure, fill it in and convert it to a
% sizes array.
%
% Note that in this example, the values are hard coded.  This is not a
% recommended practice as the characteristics of the block are typically
% defined by the S-function parameters.
%
sizes = simsizes;

sizes.NumContStates  = 5;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 5;
sizes.NumInputs      = 1;
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% initialize the initial conditions
%
x0  =pa.x0;

%
% str is always an empty matrix
%
str = [];

%
% initialize the array of sample times
%
ts  = [0 0];
% end mdlInitializeSizes

%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(t,x,u,pa)

V=x(1);
gamma=x(2);
h=x(3);
alpha=x(4);
q=x(5);
delta_e=u;

D=0.5*pa.rho*V^2*pa.S*CD(t,alpha,V,delta_e);
r=h+pa.RE;
Myy=0.5*pa.rho*V^2*pa.S*pa.c*( CMalpha(alpha,V)+ CMdeltaE( alpha, delta_e,V) +CMq(q,alpha, V));
L=0.5*pa.rho*V^2*pa.S*CL(alpha,V, delta_e);

T=0.5*pa.rho*V^2*pa.S*CT( pa.delta_T); %发动机平衡点飞行时，视为delta_T 常值；


    dx1=(T*cos(alpha)-D)/pa.m-pa.mu*sin(gamma)/ (r^2);


    dx2=(L+T*sin(alpha))/(pa.m*V)-(pa.mu-r*V^2)*cos(gamma)/(V*r^2);

    dx3=V*sin(gamma);
    dx4=q-dx2;
    dx5=Myy/pa.Jyy;
    
    
sys = [dx1;dx2;dx3;dx4;dx5];


% function y=CMalpha(alpha)
%     y=-0.035*alpha^2+0.036617*alpha+5.3261*1e-6;

function y=CMalpha(alpha,V)
    M=V/301.7;a=rad2deg(alpha);
    y=0.000014*M^3+0.000001*a^3-0.000322*M^2-0.000043*a^2-0.000059*M*a+0.003361*M+0.000107*a-0.020889;

% function y=CMq(q, c ,alpha, V)
%     y=c/(2*V)*q*(-6.796*alpha^2+0.3015*alpha-0.2289);
function y=CMq(q,alpha, V)
    M=V/301.7;a=rad2deg(alpha);
    y=0.0002*M^3+0.00000*a^3-0.0053*M^2-0.0003*a^2-0.0006*M*a+0.0511*M+0.000*a-0.2857;
    y=q*y/2/V;
% function y=CMdeltaE(ce, alpha, delta_e)
%     y=ce*(delta_e- alpha);
function y=CMdeltaE( alpha, delta_e, V)
    a=rad2deg(alpha);da=rad2deg(delta_e);M=V/301.7;
    y=-9.26*10^-5-a*3.19*10^-5-M*(1.22)*10^-6+da*(1.67*10^-3)+a*da*7.12*10^-6-(a*M)*(4.36)*10^-6-M*da*6.17*10^-6+a*M*da*1.72*10^-7;
% function y=CL( alpha, V)
%     y=0.6203*alpha;
    
function y=CL( alpha, V, delta_e)
    M=V/301.7;a=rad2deg(alpha);da=rad2deg(delta_e);
    y=0.0005*M^3+0.000000*a^3-0.0068*M^2+0.0004*a^2-0.00000*M*a+0.0132*M+0.0179*a+0.1319;
    y1=-1.45*10^(-5)+a*1.01*10^(-4)+M*7.10*10^-6-da*4.14*10^-4-a*da*3.51*10^-6+a*M*4.7*10^-6+M*da*8.72*10^-6-a*M*da*1.70*10^-7;
    y=y+y1;
% function y=CD( alpha)
%     y=0.645*alpha^2+0.0043378*alpha+0.003772;    
   
function y=CD(t, alpha, V,delta_e)
    M=V/301.7;a=rad2deg(alpha);da=rad2deg(delta_e);
    y=0.001*M^3+0.000000*a^3-0.0182*M^2+0.0005*a^2-0.0002*M*a+0.1082*M+0.0038*a-0.1763; 
   y1=4.1398*10^-4+a*2.1687*10^-5+M*(-0.9586)*10^-4+da*(-3.5219)*10^-5+a*M*da*(-6.9843)*10^-7+a^2*2.4897*10^-6+M^2*4.3644*10^-6+da^2*7.4629*10^-6+(a*M*da)^2*1.4026*10^-12;
    y=y+y1;
function y=CT( delta_T)
    if(delta_T<1)
        y=0.02576*delta_T;
    else
        y=0.0224+0.00336*delta_T;
    end


% end mdlDerivatives

%
%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%
function sys=mdlUpdate(t,x,u)

sys = [];

% end mdlUpdate

%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u,pa)

sys = x

% end mdlOutputs

%
%=============================================================================
% mdlGetTimeOfNextaVarHit
% Return the time of the next hit for this block.  Note that the result is
% absolute time.  Note that this function is only used when you specify a
% variable discrete-time sample time [-2 0] in the sample time array in
% mdlInitializeSizes.
%=============================================================================
%
function sys=mdlGetTimeOfNextVarHit(t,x,u)

sampleTime = 1;    %  Example, set the next hit to be one second later.
sys = t + sampleTime;

% end mdlGetTimeOfNextVarHit

%
%=============================================================================
% mdlTerminate
% Perform any end of simulation tasks.
%=============================================================================
%
function sys=mdlTerminate(t,x,u)

sys = [];

% end mdlTerminate
