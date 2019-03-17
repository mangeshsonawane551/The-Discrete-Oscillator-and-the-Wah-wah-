%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Program Details:  implementation of Finite difference scheme for a 
%simple undamped harmonic oscillator
% In this x3 is future variable i.e. (n+1)
% x2 is present variable i.e. (n)
% x1 is past variable i.e. (n-1)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
clear all;
clc;

%Sample rate
Fs = 44100;

%Fundamental frequency in Hz
f0 = 20;

%Simulation duration
tEnd = 0.2;An

%Initial displacement
x0 = 3;

%Initial velocity
v0 = 0.6;


%Timestep
T = 1/Fs;

%Angular frequency
w0 = 2*pi*f0;
 
%Stability check 
if(T >= 2/w0)
    error('This is unstable');
end
 
%coefficient value of FDTD equation 
coefficient1 = 2-(T.^2)*(w0).^2;

%Number of frames
N = floor(tEnd*Fs);

%Value of output at timestep n=1
x1 = x0;

%value of output at timestep n=2
x2 = x0+T*v0;

%Initialise output vector
out1 = zeros(N,1); 
out1(1) = x1; 
out1(2) = x2; 
%-------------------------------------------------------------------------%
               %Finite difference scheme calculations
%-------------------------------------------------------------------------%

for n = 3:N
    
    %Finite difference scheme equation 
    x3 = coefficient1*x2-x1; 
    
    %Update value of future i.e. x(n+1) to output variable
    out1(n) = x3; 
    
    %Update the values of each step
    x1 = x2;
    x2 = x3;
end
%-------------------------------------------------------------------------%
                %Calculation uding general solution
%-------------------------------------------------------------------------%
% using the equation x(t) = C sin(?0t + ?) where C=sqrt(A?2 + B?2)
A = x0;

B = v0/w0;

C= sqrt((A.^2)+(B.^2));

%Phase shift
phi = atan(-B/A);


n=0:tEnd/Fs:tEnd;

%general equation
y= C*cos(w0*n +phi);



% play sound
%soundsc((y),Fs);

figure(1)

%Plot finite difference scheme method calculations
subplot(3,1,1)
plot([0:N-1]*T,out1);
axis tight
grid on
title('Finite difference scheme method output');

%Plot general solution 
subplot(3,1,2)
plot([0:tEnd/Fs:tEnd],y);
axis tight
grid on
title('General solution output')

%If zoomed at peak,the ?output? angular frequency  is actually slightly 
%higher than would  in the physical system.
subplot(3,1,3)
plot([0:tEnd/Fs:tEnd],y);
hold on
plot([0:N-1]*T,out1);
grid on
xlabel('time'); 
ylabel('x'); 
title('SHO: Scheme Output'); 
axis tight
