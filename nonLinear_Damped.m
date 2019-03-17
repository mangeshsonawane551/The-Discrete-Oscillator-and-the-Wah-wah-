%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Author: Mangesh Sonawane
% Program Details: An implementation of Finite difference scheme for a 
% nonlinear damped cubic harmonic oscillator.
% In this x3 is future variable i.e. (n+1)
% x2 is present variable i.e. (n)
% x1 is past variable i.e. (n-1)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
clear all;
clc;
%-------------------------------------------------------------------------%
                    % Reads input audio file
%-------------------------------------------------------------------------%

%Input audio file is read in variable 'f' which acts as a force in simple 
%harmonic oscillation system
[f,Fs] = audioread('guitardry.wav');

if size(f,2) == 2                            %detects if stereo or mono
    f = sum(f,2) / size(f,2);                %convert to mono
end

%Determines length of input
len = length(f);
%-------------------------------------------------------------------------%

% Sample Rate
Fs=44100;

%Time Step
T = 1/Fs;

%Initial displacement
x0 = 0;

% Initial Velocity
v0 = 0;

% Damping coefficient
alpha = 8000;
if alpha<600 || alpha>11000
    error('Please maintain range of damping coefficiect')
end


% Maximum frequency
wmax = 70000*2*pi;
if wmax>100000*2*pi || wmax<1000*2*pi
    error('Maintain wmax range ')
end

% Minimum frequency
wmin = 40000*2*pi;
if wmin>wmax
    error('Please enter value for wmin less than wmax')
end

if wmin<1*2*pi || wmin>100000*2*pi
    error('Please maintain the wmin range')
end



% depth size
Delta = wmax-wmin;

%The ?rate? (measured in Hz) at which the auto-wah oscillates up and down in frequency
fwah=2;

%Calculation for sinusoidal function
n=0:len;
% In this function , I increased size of depth to explore different effects that
% sounds amazing
w1 = (wmax-fwah/2) *cos(2*pi*fwah*n*T) + (fwah/2 + Delta);

%-------------------------------------------------------------------------%
               %Finite difference scheme calculations
%-------------------------------------------------------------------------%

% Value of output at timestep n=1
x1 = x0 ;
%value of output at timestep n=2
x2 = (x0+T*v0) ;
%Initialise output vector and its updates
out1 = zeros(len,1); 
out1(1) = x1 +f(1); 
out1(2) = x2 + f(2);

for i=3:len
    % Calculation of coefficients of finite difference scheme equation at
    % each step as angular frequency changes as per time
    
    coefficient1 = ( 1 -(alpha*T/2) + ((w1(i).^4 * (x2).^2 * T.^2)/2))/...
        ((1+(alpha*T/2) + ((w1(i).^2 * (x2).^2 * T.^2)/2)));
    
    coefficient2 = 2/...
        ((1+(alpha*T/2) + ((w1(i).^2 * (x2).^2 * T.^2)/2)));
    
    coefficient3 = T.^2 /...
        ((1+(alpha*T/2) + ((w1(i).^2 * (x2).^2 * T.^2)/2)));
    
    % Finite difference equation calculation
    x3 = -coefficient1*x1 + coefficient2 * x2 + f(i)*coefficient3;
   
    %Update value in output vector
    out1(i) = x3;
    
    % Update values in past and present variables
    x1=x2;
    x2=x3;
end

% Normalise output
 out1 = 3*(out1)/max(abs(out1));
 
 % Play sound
soundsc(out1,Fs);

%Plot Output
xaxis = 0:len-1;
plot(xaxis/Fs,out1);
xlabel('Time');
ylabel('Amplitude');
title('Non-Linear Wah-Wah Effect');
grid on;
axis tight;


