%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%Author: Mangesh Sonawane 
%Program Details: An implementation  that simulates a wah-wah effect 
%using a damped simple harmonic oscillator and the finite difference time
%domain method
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

[f,Fs]=audioread('guitardry.wav' );

% Checks if stereo and converts to mono
if size(f,2) == 2                            %detects if stereo or mono
    f = sum(f,2) / size(f,2);                %convert to mono
end
%-------------------------------------------------------------------------%

% Time step
T = 1/Fs;

% Initial Displacement at timestep 0
x0 = 0;

% Initial volocity at timestep 0
v0 = 0.0;

% Peak natural frequency
wmax = (1000*2*pi);

% Minimum natural frequency
wmin = (100*2*pi);

% rate? (measured in Hz) at which the auto-wah oscillates up and down in frequency
fwah=2;

% Depth of Wah Effect
Depth = wmax-wmin;

% Damping coefficient
alpha = 10000;

% Length of input force/audio
len = length(f);

% Initialise lenght for angular frequency
n = 0:len;
% Angular frequency is the sinusoidal which oscillates between a 
%?peak natural frequency? ?max and a ?minimum natural frequency? ?min
w0 = wmax*cos(2*pi*fwah*n*T) + Depth;

%-------------------------------------------------------------------------%
               %Finite difference scheme calculations
%-------------------------------------------------------------------------%
% Calculation of coefficients foe finite difference scheme
coefficient1 = ( (alpha*T/2)-1)/((alpha*T/2)+1);

coefficient2 =((( T.^2 * w0.^2 -2)/((1 + alpha*T/2))));

coefficient3 = T.^2/ ((1 + alpha*T/2));

% Initialising values for output for step 1 
x1 = x0 + f(1);
% Initialising values for output for step 2
x2 = (x0+T*v0) + f(2);

% updating values of time step 1 and 2 in the output vector
wahwah(1) = x1+f(1); 
wahwah(2) = x2+f(2);


for i = 3:len
     if T >= ((2/w0(i).^2)*(-alpha+sqrt(alpha.^2 + w0(i).^2)))
        error('This is unstable')
    end
    % Finite difference scheme equation claculation
    x3 = coefficient1*x1 - x2 .* coefficient2(i) + f(i)*coefficient3;
    
    % Update future value in output vector
    wahwah(i) = x3; 
    
    %Update values 'present' and 'past' variables
    x1 = x2;
    x2 = x3;
end


% Normalise output
wahwah = wahwah/max(abs(wahwah));

%Plot output of wahwah
figure(1)
plot([0:len-1]*T,wahwah);
xlabel('Time');
ylabel('input signal ');
axis tight;
grid on;
title('Wah - Wah effect ');

% Play sound
soundsc(wahwah,Fs)
