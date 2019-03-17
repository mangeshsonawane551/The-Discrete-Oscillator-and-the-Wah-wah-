%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Author: Mangesh Sonawane 
%Program Details: An implementation of Finite difference scheme for a 
%simple damped harmonic oscillator
% In this x3 is future variable i.e. (n+1)
% x2 is present variable i.e. (n)
% x1 is past variable i.e. (n-1)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

clear all;
clc;

%Sample rate
Fs = 44100;

%Fundamental frequency in Hz
f0 = 125;

%Simulation duration in sec
tEnd = 0.2;

%Initial displacement
x0 = 0.9;

%Initial Velocity
v0 = 0.2;

%Time Step
T = 1/Fs;

%Angular Frequency
w0 = 2*pi*f0;

%Number of frames
N = floor(tEnd*Fs);

%Value of output at timestep n=1
x1 = x0;

%value of output at timestep n=2
x2 = (x0+T*v0);

%Declare empty vectors for output of 'underdamp' , 'critical damped'  and
%'over damp'
underdamp = zeros(N,1); 
criticaldamp = zeros(N,1);
overdamp = zeros(N,1); 

%Declare damping factor for underdamp, critical damp and overdamp
%respectively

%[ alpha<w0 underdamped ,  alpha=w0 critical , alpha>w0 overdamp]
alpha = [100 w0 3*w0]; 
%-------------------------------------------------------------------------%
               %Finite difference scheme calculations
%-------------------------------------------------------------------------%
%This loop is used to get the alpha values for each condition of damping
%Counter checks for the index of alpha and assigns the condition of damping
%accordingly
for counter = 1:3
    %Stability check
    if T >= ((2/w0.^2)*(-alpha(counter)+sqrt(alpha(counter).^2 + w0.^2)))
        error('This is unstable')
    end
    
    %Calculation of coefficients  for  each state of finite difference
    %equation
    coefficient1 = ((2-(w0.^2) * (T.^2))/(1+T*alpha(counter)/2));
    
    coefficient2 = ((1-T*alpha(counter)/2)/(1+T*alpha(counter)/2));
    
    %Initialise an temporary output vector for calculations
    out1 = zeros(N,1); 
    
    %Update the values of temporary output vector for time step 1 and 2
    out1(1) = x1; 
    out1(2) = x2; 
    
  %This for loop is used for calculation of Finite difference scheme  
for n=3:N
    
    x3 = (coefficient1*x2)-x1*coefficient2;
    %Update value in temporary output vector
    out1(n) = x3; 
    
    %Updates values of the variables for given 'present' and 'past' step for
    %each updation
    x1 = x2;
    x2 = x3;
end

%Checks for the value of damped used and assigns temporary out to permanent
%declared output respectively
if counter == 1
    %Updates value of temporary output to underdamp vector
    underdamp = out1;
    
elseif counter == 2
    %Updates value of temporary output to criticaldamp vector
    criticaldamp = out1;
    
else
    %Updates value of temporary output to overdamp vector
    overdamp=out1;
end

end

%Plot
figure(1)
%Plot underdamp condition
subplot(3,1,1)
plot([0:N-1]*T,underdamp);
xlabel('time');
ylabel('x(t)');
title('Underdamped');
grid on

%Plot Critical damp condition
subplot(3,1,2)
plot([0:N-1]*T,criticaldamp);
xlabel('time');
ylabel('x(t)');
title('Critical');
grid on

%Plot Overdamp condition
subplot(3,1,3)
plot([0:N-1]*T,overdamp);
xlabel('time');
ylabel('x(t)');
title('Overdamped');
grid on

