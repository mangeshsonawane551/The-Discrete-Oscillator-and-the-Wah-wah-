%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Author:MANGESH SONAWANE
% Program Details: An implementation of Finite difference scheme for a 
%  undamped harmonic oscillator with three spring and two masses
%  |-~~~-M1-~~~-M2-~~~-|   where '~~~' represents spring k
% In this x3 is future variable i.e. (n+1)
% x2 is present variable i.e. (n)
% x1 is past variable i.e. (n-1)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
clear all;
clc;
% Sample rate
Fs=44100;

% Time Step 
T=1/Fs;

% Initial displacement value at time step 0
x1 = [0;1];

% Initial velocity
v0 = [0;1];

% Displacememt at time step 1
x2 = [x1(1)+T*v0(1);x1(2)+T*v0(2)];

% Duration of simulation
tEnd = 50;

% Number of frames
N =tEnd*Fs;

% Masses
M1=30 ;     % value of Mass 1
M2=50;      % Value of mass 2

% M matrix initialisation
M=[M1 0; M2 0];

% Spring constants
k1=50;     %value for spring 1
k2=90;     % value for spring 2
k3=10;       % value for spring 3

% Initialising matrix for K
K=[ (-k1-k2)/M1 k2/M1; k2/M2 (-k3-k2)/M2];

%-------------------------------------------------------------------------%
               % Stability check
%-------------------------------------------------------------------------%

w1 = sqrt(k1+k2 /M1);
w2 = sqrt(k1+k3/M2);

if T>2/w1 || T>2/w2
    error('This is not stable')
end

% Calculation of coefficient for x(n)
coefficient1= ((T.^2 *K )+ 2*eye(2));

% Initialise output  for M1 ar time step 0 and 1
out1(1) = x1(1);
out1(2) = x2(1);

% Initialise output for M2 at time step 0 and 1
out2(1) = x1(2);
out2(1) = x2(2);
%-------------------------------------------------------------------------%
               %Finite difference scheme calculations
%-------------------------------------------------------------------------%
for n=3:N
    
    % Finite difference scheme calculation
   x3= (coefficient1*x2 - x1);
   
   %update values in output vector
   out1(n)=x3(1);
   out2(n)=x3(2);
   
   %Update values of past and present variables
   x1=x2;
   x2=x3;
end
   
% Plot
figure(1)
plot([0:N-1]*T, out2,[0:N-1]*T, out1);
legend('Oscillation at Mass 2','Oscillation at Mass 1');
grid on;
xlabel('Time');
ylabel('Amplitude');
axis tight
