%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Author: MANGESH SONAWANE
% Program Details: An implementation of Finite difference scheme for a 
% non linear simple damped harmonic oscillator which rather than a decaying
% sinusoidal oscillation,generates a downward pitch glide 
% when ? = 0.000000001, and  ?1 = sqrt(2?400) and initial position and 
% velocity are 1 and 0.1 respectively
% In this x3 is future variable i.e. (n+1)
% x2 is present variable i.e. (n)
% x1 is past variable i.e. (n-1)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
clear all;
clc;

%Sample rate
Fs=44100;

%Time step
T=1/Fs;

%Damping Coefficient
alpha = 0.000000001;

% Initial displacement
x0 = 1;

%Initial Velocity
v0 = 0.1;

% Simulation duration
tEnd=10;

%Number of frames
N=tEnd*Fs;

% Angular frequency
w1 =sqrt(2*pi*400);

%-------------------------------------------------------------------------%
               %Finite difference scheme calculations
%-------------------------------------------------------------------------%
% Calculation of coefficient for FD equation
coefficient1 = 4/(((1+alpha/(2*T))));
coefficient2 = (1- alpha/ (2*T))/ ((1+alpha/(2*T)));

%Update value 
x1 = x0 ;
x2 = (x0+T*v0) ;
out1(1) = x1 ;
out1(2) = x2;

for i=3:N-1
    
    %Finite difference equation calculation
    x3= ((x2*coefficient1)./(2+(w1.^4)* (T.^2 )*(x2).^2)) - coefficient2*x1;
   
    % Update value in the output vector
    out1(i) = x3;
    
    % Update values in past and present variables
    x1=x2;
    x2=x3;
end

%Play sound
soundsc(out1,Fs);
plot(out1);
