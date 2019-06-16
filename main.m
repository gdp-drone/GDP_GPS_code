%{
Author: Lukasz Sliwinski luki3141@gmail.com
date: June 2019
project: ADD GPD 2019

This code analyzes the data in the localized<>.txt files
The data has to be fisrt  imported without the third column (altitude) and
then saved as data.mat.

The code plots the GPS signal, calculates the Fourier decomposition and
calculates characteristic quantities of the GPS signal for both static and
circle tests.
The code presented here is only for base accuracy for 2 metres because it
is repetetive - however the while uncommented code can be found in raw_code
folder along with data.mat from all the tests.
The purpose of this file to explain the mode of work in simple way.
%}

clc
clear
close all

load data.mat

% Transformng data to arrays
C2RF = table2array(localizedAcc2mCircle);
C2GL = table2array(localizedAcc2mCircleGlobal);
S2RF = table2array(localizedAcc2mStatic);
S2GL = table2array(localizedAcc2mStaticGlobal);


% Adding miliseconds to the time of measurement
C2RF(:,3)=C2RF(:,3)+C2RF(:,4)/1000000000;
C2GL(:,3)=C2GL(:,3)+C2GL(:,4)/1000000000;
S2RF(:,3)=S2RF(:,3)+S2RF(:,4)/1000000000;
S2GL(:,3)=S2GL(:,3)+S2GL(:,4)/1000000000;


%% Offsetting static data around the mean
S2RF(:,1)=S2RF(:,1)-mean(S2RF(:,1));
S2RF(:,2)=S2RF(:,2)-mean(S2RF(:,2));
S2GL(:,1)=S2GL(:,1)-mean(S2GL(:,1));
S2GL(:,2)=S2GL(:,2)-mean(S2GL(:,2));




%% Plotting circle data
fsize=12;

h= figure
plot(C2RF(:,1),C2RF(:,2),'.')
grid on
grid minor
xlim([-8 8])
ylim([-4 16])
xlabel('x coordinate [m]','FontSize',fsize) 
ylabel('y coordinate [m]','FontSize',fsize) 
set(h, 'WindowStyle', 'Docked');

h = figure
plot(C2GL(:,1),C2GL(:,2),'.')
grid on
grid minor
xlim([-8 8])
ylim([-4 16])
xlabel('x coordinate [m]','FontSize',fsize) 
ylabel('y coordinate [m]','FontSize',fsize) 
set(h, 'WindowStyle', 'Docked');


close all
%% PLotting static

h = figure
plot(S2RF(:,1),S2RF(:,2),'.')
grid on
grid minor
ylim([-1 0.5])
xlabel('x coordinate [m]','FontSize',fsize) 
ylabel('y coordinate [m]','FontSize',fsize)
set(h, 'WindowStyle', 'Docked');
saveas(gcf,'Static2RF.eps','epsc')

h = figure
plot(S2GL(:,1),S2GL(:,2),'.')
grid on
grid minor
ylim([-1 0.5])
xlabel('x coordinate [m]','FontSize',fsize) 
ylabel('y coordinate [m]','FontSize',fsize) 
set(h, 'WindowStyle', 'Docked');
saveas(gcf,'Static2GL.eps','epsc')

close all
%% Spectral analysis using fourier transform according to
% https://www.mathworks.com/help/matlab/ref/fft.html


% Low frequency analysis
% x signal
% Interpolation to resample the signal
T = [S2RF(1,3):0.1:S2RF(end,3)];
vq2 = interp1(S2RF(:,3),S2RF(:,1),T,'spline');
Y = fft(vq2);
Fs = 10; 
L =size(T);
L=L(2)
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
h = figure
plot(f,(P1)) 
grid on
grid minor
xlabel('Frequency (Hz)','FontSize',fsize)
ylabel('Amplitude [m]','FontSize',fsize)
set(h, 'WindowStyle', 'Docked');

%y signal
%Interpolation to resample the signal
T = [S2RF(1,3):0.1:S2RF(end,3)];
vq2 = interp1(S2RF(:,3),S2RF(:,2),T,'spline');
Y = fft(vq2);
Fs = 10; 
L =size(T);
L=L(2)
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
hold on
plot(f,(P1)) 
legend({'x coordinate','y coordinate'},'FontSize',fsize);
xlim([0 0.1])
saveas(gcf,'Fourier2low.eps','epsc')


% High frequency analysis (same bur different limits
% x signal
% Interpolation to resample the domain
T = [S6RF(1,3):0.1:S6RF(end,3)];
vq2 = interp1(S6RF(:,3),S6RF(:,1),T,'spline');
Y = fft(vq2);
Fs = 10; 
L =size(T);
L=L(2)
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
h = figure
plot(f,P1) 
grid on
grid minor
xlabel('Frequency (Hz)','FontSize',fsize)
ylabel('Amplitude [m]','FontSize',fsize)
set(h, 'WindowStyle', 'Docked');

% y signal
% Interpolation to resample the domain
T = [S6RF(1,3):0.1:S6RF(end,3)];
vq2 = interp1(S6RF(:,3),S6RF(:,2),T,'spline');
Y = fft(vq2);
Fs = 10; 
L =size(T);
L=L(2)
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
hold on
plot(f,P1) 
legend({'x coordinate','y coordinate'},'FontSize',fsize);
xlim([0.1 4])
saveas(gcf,'Fourier6tail.eps','epsc')


T = [S2RF(1,3):0.1:S2RF(end,3)];
vq2 = interp1(S2RF(:,3),S2RF(:,1),T,'spline');
Y = fft(vq2);
Fs = 10; 
L =size(T);
L=L(2)
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
h = figure
plot(f,(P1)) 
grid on
grid minor
xlabel('Frequency (Hz)','FontSize',fsize)
ylabel('Amplitude [m]','FontSize',fsize)
set(h, 'WindowStyle', 'Docked');

close all

%% Calculating standard deviations for /raw/fix

stdS2x = std(S2RF(:,1));
stdS2y = std(S2RF(:,2));
stdS2 = sqrt(stdS2x.^2+stdS2y.^2);

%% Reduction of circle data

% Singing-out data lying on a circle
circle = C2RF(350:780,:);
h=figure
plot(C2RF(:,1),C2RF(:,2),'.')
minErr =1000000;
xCircle=0;
% Finding the best fit to minimise the error 
for yc=[7.76:0.0001:7.80]
    for xc=[-0.72:(0.0001):-0.71]
        radius = ((circle(:,1)-xc).^2+(circle(:,2)-yc).^2).^(1/2);
        [xc,yc]
        mRad = mean(radius);
        err = radius-mRad;
        sErr=sum(abs(err));
        if (minErr>sErr)
            minErr=sErr;
            xCircle2 = xc;
            yCircle2= yc;
            radCricle = mRad;
        end
    end
end
radius = ((circle(:,1)-xCircle2).^2+(circle(:,2)-yCircle2).^2).^(1/2);
mRad2 = mean(radius);
err2 = radius-mRad2;
stdC2 = std(err2)
hold on
th = 0:pi/50:2*pi;
xunit = mRad2 * cos(th) + xCircle2;
yunit = mRad2 * sin(th) + yCircle2;
plot(xunit, yunit);
grid on
grid minor
xlabel('x coordinate [m]','FontSize',fsize)
ylabel('y coordinate [m]','FontSize',fsize)
legend({'GPS data','circle fit'},'FontSize',fsize,'Location','northoutside');
set(h, 'WindowStyle', 'Docked');
saveas(gcf,'Circle2.eps','epsc')
hold off


%Plotting error histograms
h = figure
histogram(err2,20)
xlabel('error [m]','FontSize',fsize) 
set(h, 'WindowStyle', 'Docked');
saveas(gcf,'histogram2.eps','epsc')


%% Calculating drift velocity

spac=10; % data points space by 2 seconds (5 HZ rate)
L = size(S2RF,1);
i=[1:1:(L-spac)];
dist = hypot(S2RF(i+spac,1)-S2RF(i,1),S2RF(i+spac,2)-S2RF(i,2));
vel = dist./(S2RF(i+spac,3)-S2RF(i,3));     
meanVelS2 = mean(vel)
maxVelS2 = max(vel)
