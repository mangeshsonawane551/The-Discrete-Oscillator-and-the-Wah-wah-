classdef Wah_Wah_GUI < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                     matlab.ui.Figure
        WahFrequencyKnobLabel        matlab.ui.control.Label
        WahFrequencyKnob             matlab.ui.control.Knob
        Peaknaturalfrequency103KnobLabel  matlab.ui.control.Label
        Peaknaturalfrequency103Knob  matlab.ui.control.Knob
        DepthKnobLabel               matlab.ui.control.Label
        DepthKnob                    matlab.ui.control.Knob
        PlayoriginalsoundButton      matlab.ui.control.Button
        PlaywahwaheffectButton       matlab.ui.control.Button
        WAHWAHLabel                  matlab.ui.control.Label
    end

    methods (Access = private)

        % Value changed function: WahFrequencyKnob
        function WahFrequencyKnobValueChanged(app, event)
            value = app.WahFrequencyKnob.Value;
            
        end

        % Value changed function: Peaknaturalfrequency103Knob
        function Peaknaturalfrequency103KnobValueChanged(app, event)
            value = app.Peaknaturalfrequency103Knob.Value;
            
        end

        % Value changed function: DepthKnob
        function DepthKnobValueChanged(app, event)
            value = app.DepthKnob.Value;
            
        end

        % Button pushed function: PlayoriginalsoundButton
        function PlayoriginalsoundButtonPushed(app, event)
            [f,Fs] = audioread('guitardry.wav');
            soundsc(f,Fs);
        end

        % Button pushed function: PlaywahwaheffectButton
        function PlaywahwaheffectButtonPushed(app, event)

            
            
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

wmax = app.Peaknaturalfrequency103Knob.Value;
% Maximum frequency
wmax = wmax * 10.^3*2*pi;

if wmax>100000*2*pi || wmax<1000*2*pi
    error('Maintain wmax range ')
end
depth=app.DepthKnob.Value;
depth1=depth*10^3 * 2*pi;
wmin=wmax-depth1;

% Minimum frequency

if wmin>wmax
    depth1=depth*10 * 2*pi
    wmin = wmax-depth1;
end

% if wmin<1*2*pi || wmin>100000*2*pi
%     error('Please maintain the wmin range')
% end



% depth size
Delta = wmax-wmin;

%The ?rate? (measured in Hz) at which the auto-wah oscillates up and down in frequency
fwah=app.WahFrequencyKnob.Value;

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
 out1 = (out1)/max(abs(out1));
 
 % Play sound
soundsc(out1,Fs);



        end
    end

    % App initialization and construction
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure
            app.UIFigure = uifigure;
            app.UIFigure.Position = [100 100 640 480];
            app.UIFigure.Name = 'UI Figure';

            % Create WahFrequencyKnobLabel
            app.WahFrequencyKnobLabel = uilabel(app.UIFigure);
            app.WahFrequencyKnobLabel.HorizontalAlignment = 'center';
            app.WahFrequencyKnobLabel.Position = [461 191 89 22];
            app.WahFrequencyKnobLabel.Text = 'Wah Frequency';

            % Create WahFrequencyKnob
            app.WahFrequencyKnob = uiknob(app.UIFigure, 'continuous');
            app.WahFrequencyKnob.Limits = [0.1 20];
            app.WahFrequencyKnob.ValueChangedFcn = createCallbackFcn(app, @WahFrequencyKnobValueChanged, true);
            app.WahFrequencyKnob.Position = [474 247 60 60];
            app.WahFrequencyKnob.Value = 4;

            % Create Peaknaturalfrequency103KnobLabel
            app.Peaknaturalfrequency103KnobLabel = uilabel(app.UIFigure);
            app.Peaknaturalfrequency103KnobLabel.HorizontalAlignment = 'center';
            app.Peaknaturalfrequency103KnobLabel.Position = [9 191 190 22];
            app.Peaknaturalfrequency103KnobLabel.Text = 'Peak natural frequency (*10^3)';

            % Create Peaknaturalfrequency103Knob
            app.Peaknaturalfrequency103Knob = uiknob(app.UIFigure, 'continuous');
            app.Peaknaturalfrequency103Knob.Limits = [2 99];
            app.Peaknaturalfrequency103Knob.ValueChangedFcn = createCallbackFcn(app, @Peaknaturalfrequency103KnobValueChanged, true);
            app.Peaknaturalfrequency103Knob.Position = [74 247 60 60];
            app.Peaknaturalfrequency103Knob.Value = 80;

            % Create DepthKnobLabel
            app.DepthKnobLabel = uilabel(app.UIFigure);
            app.DepthKnobLabel.HorizontalAlignment = 'center';
            app.DepthKnobLabel.Position = [290 191 38 22];
            app.DepthKnobLabel.Text = 'Depth';

            % Create DepthKnob
            app.DepthKnob = uiknob(app.UIFigure, 'continuous');
            app.DepthKnob.Limits = [1 95];
            app.DepthKnob.ValueChangedFcn = createCallbackFcn(app, @DepthKnobValueChanged, true);
            app.DepthKnob.Position = [278 247 60 60];
            app.DepthKnob.Value = 9;

            % Create PlayoriginalsoundButton
            app.PlayoriginalsoundButton = uibutton(app.UIFigure, 'push');
            app.PlayoriginalsoundButton.ButtonPushedFcn = createCallbackFcn(app, @PlayoriginalsoundButtonPushed, true);
            app.PlayoriginalsoundButton.Position = [93 94 149 38];
            app.PlayoriginalsoundButton.Text = 'Play original sound';

            % Create PlaywahwaheffectButton
            app.PlaywahwaheffectButton = uibutton(app.UIFigure, 'push');
            app.PlaywahwaheffectButton.ButtonPushedFcn = createCallbackFcn(app, @PlaywahwaheffectButtonPushed, true);
            app.PlaywahwaheffectButton.Position = [361 94 164 38];
            app.PlaywahwaheffectButton.Text = 'Play wah wah effect';

            % Create WAHWAHLabel
            app.WAHWAHLabel = uilabel(app.UIFigure);
            app.WAHWAHLabel.FontName = 'Georgia';
            app.WAHWAHLabel.FontSize = 34;
            app.WAHWAHLabel.FontWeight = 'bold';
            app.WAHWAHLabel.Position = [202 390 214 42];
            app.WAHWAHLabel.Text = 'WAH WAH ';
        end
    end

    methods (Access = public)

        % Construct app
        function app = S1889125_wahwah_GUI

            % Create and configure components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end
