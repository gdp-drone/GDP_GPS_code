clc
clear
close all

load data.mat

C2RF = table2array(localizedAcc2mCircle);
C2GL = table2array(localizedAcc2mCircleGlobal);
ARF = table2array(localizedAndroid);
S2RF = table2array(localizedAcc2mStatic);
S2GL = table2array(localizedAcc2mStaticGlobal);
WRF = table2array(localizedAcc2mWlkRF);
WGL = table2array(localizedAcc2mWlkGl);
S6RF = table2array(localizedAcc6mStatic);
S6GL = table2array(localizedAcc6mStaticGlobal);
C6RF = table2array(localizedAcc6mCircle);
C6GL = table2array(localizedAcc6mCircleGlobal);
SRRF = table2array(localizedStaticRoverRF);
SRGL = table2array(localizedStaticRoverGL);


C2RF(:,3)=C2RF(:,3)+C2RF(:,4)/1000000000;
C2GL(:,3)=C2GL(:,3)+C2GL(:,4)/1000000000;
ARF(:,3)=ARF(:,3)+ARF(:,4)/1000000000;
S2RF(:,3)=S2RF(:,3)+S2RF(:,4)/1000000000;
S2GL(:,3)=S2GL(:,3)+S2GL(:,4)/1000000000;
WRF(:,3)=WRF(:,3)+WRF(:,4)/1000000000;
WGL(:,3)=WGL(:,3)+WGL(:,4)/1000000000;
C6RF(:,3)=C6RF(:,3)+C6RF(:,4)/1000000000;
C6GL(:,3)=C6GL(:,3)+C6GL(:,4)/1000000000;
S6RF(:,3)=S6RF(:,3)+S6RF(:,4)/1000000000;
S6GL(:,3)=S6GL(:,3)+S6GL(:,4)/1000000000;
SRRF(:,3)=SRRF(:,3)+SRRF(:,4)/1000000000;
SRGL(:,3)=SRGL(:,3)+SRGL(:,4)/1000000000;

%% Centering about mean
S6RF(:,1)=S6RF(:,1)-mean(S6RF(:,1));
S6RF(:,2)=S6RF(:,2)-mean(S6RF(:,2));
S6GL(:,1)=S6GL(:,1)-mean(S6GL(:,1));
S6GL(:,2)=S6GL(:,2)-mean(S6GL(:,2));
S2RF(:,1)=S2RF(:,1)-mean(S2RF(:,1));
S2RF(:,2)=S2RF(:,2)-mean(S2RF(:,2));
S2GL(:,1)=S2GL(:,1)-mean(S2GL(:,1));
S2GL(:,2)=S2GL(:,2)-mean(S2GL(:,2));
SRRF(:,1)=SRRF(:,1)-mean(SRRF(:,1));
SRRF(:,2)=SRRF(:,2)-mean(SRRF(:,2));
SRGL(:,1)=SRGL(:,1)-mean(SRGL(:,1));
SRGL(:,2)=SRGL(:,2)-mean(SRGL(:,2));




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
h= figure
plot(C6RF(:,1),C6RF(:,2),'.')
grid on
grid minor
xlim([-8 8])
ylim([-4 16])
xlabel('x coordinate [m]','FontSize',fsize) 
ylabel('y coordinate [m]','FontSize',fsize) 
set(h, 'WindowStyle', 'Docked');

h = figure
subplot(1,2,1);
plot(C2GL(:,1),C2GL(:,2),'.')
grid on
grid minor
xlim([-8 8])
ylim([-4 16])
xlabel('x coordinate [m]','FontSize',fsize) 
ylabel('y coordinate [m]','FontSize',fsize) 
set(h, 'WindowStyle', 'Docked');
subplot(1,2,2);
plot(C6GL(:,1),C6GL(:,2),'.')
grid on
grid minor
xlim([-8 8])
ylim([-4 16])
xlabel('x coordinate [m]','FontSize',fsize) 
ylabel('y coordinate [m]','FontSize',fsize) 
set(h, 'WindowStyle', 'Docked');
saveas(gcf,'globalcircles.eps','epsc')

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
plot(S6RF(:,1),S6RF(:,2),'.')
grid on
grid minor
ylim([-2 2])
xlabel('x coordinate [m]','FontSize',fsize) 
ylabel('y coordinate [m]','FontSize',fsize)
set(h, 'WindowStyle', 'Docked');
saveas(gcf,'Static6RF.eps','epsc')

h = figure
plot(SRRF(:,1),SRRF(:,2),'.')
grid on
grid minor
ylim([-2 2])
xlabel('x coordinate [m]','FontSize',fsize) 
ylabel('y coordinate [m]','FontSize',fsize)
set(h, 'WindowStyle', 'Docked');
saveas(gcf,'StaticRRF.eps','epsc')

h = figure
plot(S2GL(:,1),S2GL(:,2),'.')
grid on
grid minor
ylim([-1 0.5])
xlabel('x coordinate [m]','FontSize',fsize) 
ylabel('y coordinate [m]','FontSize',fsize) 
set(h, 'WindowStyle', 'Docked');
saveas(gcf,'Static2GL.eps','epsc')
h = figure
plot(S6GL(:,1),S6GL(:,2),'.')
grid on
grid minor
ylim([-2 2])
xlabel('x coordinate [m]','FontSize',fsize) 
ylabel('y coordinate [m]','FontSize',fsize) 
set(h, 'WindowStyle', 'Docked');
saveas(gcf,'Static6GL.eps','epsc')
h = figure
plot(SRGL(:,1),SRGL(:,2),'.')
grid on
grid minor
ylim([-2 2])
xlabel('x coordinate [m]','FontSize',fsize) 
ylabel('y coordinate [m]','FontSize',fsize) 
set(h, 'WindowStyle', 'Docked');
saveas(gcf,'StaticRGL.eps','epsc')

close all
%% Spectral analysis using fourier transform

%Interpolation
%Sampling domain
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
xlim([0 0.1])
saveas(gcf,'Fourier6low.eps','epsc')



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


T = [SRRF(1,3):0.1:SRRF(end,3)];
vq2 = interp1(SRRF(:,3),SRRF(:,1),T,'spline');
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

T = [SRRF(1,3):0.1:SRRF(end,3)];
vq2 = interp1(SRRF(:,3),SRRF(:,2),T,'spline');
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
saveas(gcf,'FourierRlow.eps','epsc')

%OGOn

%Interpolation
%Sampling domain
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
xlim([0.1 4])
saveas(gcf,'Fourier2tail.eps','epsc')


T = [SRRF(1,3):0.1:SRRF(end,3)];
vq2 = interp1(SRRF(:,3),SRRF(:,1),T,'spline');
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

T = [SRRF(1,3):0.1:SRRF(end,3)];
vq2 = interp1(SRRF(:,3),SRRF(:,2),T,'spline');
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
xlim([0.1 4])
saveas(gcf,'FourierRtail.eps','epsc')


close all
%% Spektrum using pspektrum

ts =timeseries(S6RF(:,1),S6RF(:,3));
% Parameters
timeLimits = seconds([0.5141964 491.6949]); % seconds
frequencyLimits = [0 2.500292]; % Hz
% Index into signal time region of interest
ts_unnamed_ROI = squeeze(ts.Data);
timeValues = ts.Time;
ts_unnamed_ROI = timetable(seconds(timeValues(:)),ts_unnamed_ROI,'VariableNames',{'Data'});
ts_unnamed_ROI = ts_unnamed_ROI(timerange(timeLimits(1),timeLimits(2),'closed'),1);

% Compute spectral estimate
% Run the function call below without output arguments to plot the results
[Pxx,F] = pspectrum(ts_unnamed_ROI, ...
    'FrequencyLimits',frequencyLimits);

h= figure;
plot(F,pow2db(Pxx))
hold on

ts =timeseries(S6RF(:,2),S6RF(:,3));
timeLimits = seconds([0.5141964 491.6949]); % seconds
frequencyLimits = [0 2.500292]; % Hz

% Index into signal time region of interest
ts_unnamed_ROI = squeeze(ts.Data);
timeValues = ts.Time;
ts_unnamed_ROI = timetable(seconds(timeValues(:)),ts_unnamed_ROI,'VariableNames',{'Data'});
ts_unnamed_ROI = ts_unnamed_ROI(timerange(timeLimits(1),timeLimits(2),'closed'),1);

% Compute spectral estimate
% Run the function call below without output arguments to plot the results
[Pxx,F] = pspectrum(ts_unnamed_ROI, ...
    'FrequencyLimits',frequencyLimits);
plot(F,pow2db(Pxx))

xlabel('Frequency (Hz)','FontSize',fsize)
ylabel('Power [dB]','FontSize',fsize)
legend({'x coordinate','y coordinate'},'FontSize',fsize);
xlim([0 2.5])
set(h, 'WindowStyle', 'Docked');



ts =timeseries(S2RF(:,1),S2RF(:,3));

% Parameters
timeLimits = seconds([0.9093992 335.1274]); % seconds
frequencyLimits = [0 2.495169]; % Hz
% Index into signal time region of interest
ts_unnamed_ROI = squeeze(ts.Data);
timeValues = ts.Time;
ts_unnamed_ROI = timetable(seconds(timeValues(:)),ts_unnamed_ROI,'VariableNames',{'Data'});
ts_unnamed_ROI = ts_unnamed_ROI(timerange(timeLimits(1),timeLimits(2),'closed'),1);

% Compute spectral estimate
% Run the function call below without output arguments to plot the results
[Pxx,F] = pspectrum(ts_unnamed_ROI, ...
    'FrequencyLimits',frequencyLimits);

h= figure;
plot(F,pow2db(Pxx))
hold on

ts =timeseries(S6RF(:,2),S6RF(:,3));
timeLimits = seconds([0.5141964 491.6949]); % seconds
frequencyLimits = [0 2.500292]; % Hz

% Index into signal time region of interest
ts_unnamed_ROI = squeeze(ts.Data);
timeValues = ts.Time;
ts_unnamed_ROI = timetable(seconds(timeValues(:)),ts_unnamed_ROI,'VariableNames',{'Data'});
ts_unnamed_ROI = ts_unnamed_ROI(timerange(timeLimits(1),timeLimits(2),'closed'),1);

% Compute spectral estimate
% Run the function call below without output arguments to plot the results
[Pxx,F] = pspectrum(ts_unnamed_ROI, ...
    'FrequencyLimits',frequencyLimits);
plot(F,pow2db(Pxx))

xlabel('Frequency (Hz)','FontSize',fsize)
ylabel('Power [dB]','FontSize',fsize)
legend({'x coordinate','y coordinate'},'FontSize',fsize);
xlim([0 2.5])
set(h, 'WindowStyle', 'Docked');






%% Calculating standard deviations for RF
stdS6x = std(S6RF(:,1));
stdS6y = std(S6RF(:,2));
stdS6 = sqrt(stdS6x.^2+stdS6y.^2);
stdS2x = std(S2RF(:,1));
stdS2y = std(S2RF(:,2));
stdS2 = sqrt(stdS2x.^2+stdS2y.^2);
stdAx = std(ARF(:,1));
stdAy = std(ARF(:,2));
stdA = sqrt(stdAx.^2+stdAy.^2);
stdSRx = std(SRRF(:,1));
stdSRy = std(SRRF(:,2));
stdSR = sqrt(stdSRx.^2+stdSRy.^2);


close all
%% Reduction of circle data
circle = C2RF(350:780,:);
h=figure
plot(C2RF(:,1),C2RF(:,2),'.')
minErr =1000000;
xCircle=0;
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


h=figure
circle = C6RF(370:end-110,:);
plot(C6RF(:,1),C6RF(:,2),'.')
minErr =1000000;
xCircle=0;
for yc=[5.69:0.0001:5.7]
    for xc=[0.79:(0.0001):0.8]
        radius = ((circle(:,1)-xc).^2+(circle(:,2)-yc).^2).^(1/2);
        [xc,yc]
        mRad = mean(radius);
        err = radius-mRad;
        sErr=sum(abs(err));
        if (minErr>sErr)
            minErr=sErr;
            xCircle6 = xc;
            yCircle6 = yc;
            radCricle = mRad;
        end
    end
end
radius = ((circle(:,1)-xCircle6).^2+(circle(:,2)-yCircle6).^2).^(1/2);
mRad6 = mean(radius);
err6 = radius-mRad6;
stdC6 = std(err6)
hold on
th = 0:pi/50:2*pi;
xunit = mRad6 * cos(th) + xCircle6;
yunit = mRad6 * sin(th) + yCircle6;
plot(xunit, yunit);
grid on
grid minor
xlabel('x coordinate [m]','FontSize',fsize)
ylabel('y coordinate [m]','FontSize',fsize)
legend({'GPS data','circle fit'},'FontSize',fsize,'Location','northoutside');
set(h, 'WindowStyle', 'Docked');
saveas(gcf,'Circle6.eps','epsc')
hold off



h = figure
histogram(err2,20)
xlabel('error [m]','FontSize',fsize) 
set(h, 'WindowStyle', 'Docked');
saveas(gcf,'histogram2.eps','epsc')
h = figure
histogram(err6,20)
xlabel('error [m]','FontSize',fsize) 
set(h, 'WindowStyle', 'Docked');
saveas(gcf,'histogram6.eps','epsc')


%% Calculating drift velocity

spac=10;
L = size(S6RF,1);
i=[1:1:(L-spac)];
dist = hypot(S6RF(i+spac,1)-S6RF(i,1),S6RF(i+spac,2)-S6RF(i,2));
vel = dist./(S6RF(i+spac,3)-S6RF(i,3));     
mean(vel)
max(vel)

spac=10;
L = size(S2RF,1);
i=[1:1:(L-spac)];
dist = hypot(S2RF(i+spac,1)-S2RF(i,1),S2RF(i+spac,2)-S2RF(i,2));
vel = dist./(S2RF(i+spac,3)-S2RF(i,3));     
mean(vel)
max(vel)

spac=10;
L = size(ARF,1);
i=[1:1:(L-spac)];
dist = hypot(ARF(i+spac,1)-ARF(i,1),ARF(i+spac,2)-ARF(i,2));
vel = dist./(ARF(i+spac,3)-ARF(i,3));     
mean(vel)
max(vel)


spac=10;
L = size(SRRF,1);
i=[1:1:(L-spac)];
dist = hypot(SRRF(i+spac,1)-SRRF(i,1),SRRF(i+spac,2)-SRRF(i,2));
vel = dist./(SRRF(i+spac,3)-SRRF(i,3));     
meanVelSR = mean(vel)
maxVelSR = max(vel)

