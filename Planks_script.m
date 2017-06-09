clear all; clc; close all
p.freq = [5.18672*10^(14), 5.48996*10^(14), 6.87858*10^(14), 7.40858*10^(14), 8.20264*10^(14)]; %freq = sort(p.freq, 'descend');
p.wl   = [578, 546.074, 435.835, 404.656, 365.483]; wl = sort(p.wl);
p.c    = 1.60217662 * 10^(-19);
p.h    = 6.62607004*10^(-34); 
%//////////////// YELLOW ///////////////////////////////////
yellow.vol = [0.73, 0.73, 0.73, 0.729, 0.728];
yellow.per = [100 80 60 40 20];
yellow.t   = [5.55 5.11 5.47 5.94;
    4.91 4.31 6.14 4.77;
    6.47 6.98 7.37 7.72;
    10.58 11.92 8.91 12.07;
    17.92 16.44 20.04 18.35];
yellow.t_mean = mean(yellow.t,2)';
yellow.t_std = std(yellow.t,1,2);
%Plot of Yellow



%______________________________________________________________

%//////////////// GREEN ///////////////////////////////////
green.vol = [0.858 0.858 0.858 0.856 0.850];
green.per = [100 80 60 40 20];
green.t   = [17.12 18.29 15.35 15.61;
    22.16 20.29 19.55 22.68;
    18.12 17.27 21.45 19.06;
    21.41 25.13 27.79 24.16;
    30.71 38.33 36.1 36.86];
green.t_mean = mean(green.t,2)';
green.t_std = std(green.t,1,2);


%______________________________________________________________

%//////////////// Individual ///////////////////////////////////

SP1.sp_gen    = [2.118 2.105 2.111 2.118 2.118;
            1.782 1.781 1.780 1.790 1.779;
            1.566 1.564 1.555 1.569 1.568;
            .8710 .8706 .8713 .8713 .8710; 
            .7690 .7694 .7692 .7689 .7691];
        
SP1.sp_gen = flipud(SP1.sp_gen);
SP1.mean_gen = mean(SP1.sp_gen,2)';
SP1.std_gen  = std(SP1.sp_gen,1,2)';


SP2.sp_gen    = [2.118 2.110 2.115 2.113 2.118;
            1.782 1.775 1.769 1.772 1.785;
            1.568 1.565 1.565 1.555 1.544;
            1.1710 1.1690 1.1680 1.1730 1.1730; 
            .7690 .7690 .7670 .7674 .7670];
        SP2.sp_gen = flipud(SP2.sp_gen);
SP2.mean_gen = mean(SP2.sp_gen,2)';
SP2.std_gen  = std(SP2.sp_gen,1,2)';


%______________________________________________________________



%//////////////// Analysis Experiment 1 ///////////////////////////////////
figure(1)
hold on
plot(yellow.per,yellow.t')
plot(yellow.per,yellow.t_mean,'k-.','linewidth',3)
errorbar(yellow.per,yellow.t_mean,yellow.t_std,'--k','linewidth',2)
ylabel('Time (s)');xlabel('Filter Percentage (%)');
title(sprintf('Yellow\n Time to reach stopping potential for each filter'))
legend('Run 1','Run 1','Run 1','Run 1','Mean')
hold off


figure(2)
hold on
plot(green.per,green.t')
plot(green.per,green.t_mean,'k-.','linewidth',3)
errorbar(green.per,green.t_mean,green.t_std,'--k','linewidth',2)
ylabel('Time (s)');xlabel('Filter Percentage (%)');
title(sprintf('GREEN\n Time to reach stopping potential for each filter'))
legend('Run 1','Run 1','Run 1','Run 1','Mean')
hold off

%Passing different amounts of the same colored light through the various
%filter types, where 100% filter means 100% of the light passes through and
%the 20% filter means 20% of the light is let through, changing the filters
%does not change the charging time, which is expected since the light will 
%have the same amount of energy regardless of filter. However, the charging
%time does increase. This is due to the overall intensity of the light
%meaning there is less current overall, but still the same amount of
%energy.

%We know different colors of light have different wavelengths and as such
%have different energies. As we move from color to color, the wavelengths
%increase which means the frequency and energy decreases, from Ultraviolet 
%-> Violet -> Blue -> Green -> Yellow and therefore, due to a decrease in 
%energy between subsequent colors, there is a decrease in stopping
%potential as expected.

%This does support the photon model of light. We come to this conclusion by
%studying the stopping potential and effects of intensity of the light.
%Using the analogy of a water wave, we expect that if we vary the flow
%of water ( electricl current ) in a hose and fill a bucket, the bucket
%will always get filled, but depending on the rate of flow, it will fill
%faster or slower. Likewise, decreasing the intensity of the light makes
%increases the overall charging time but doesn't change the energy of teh
%system.



%______________________________________________________________

%//////////////// Analysis Experiment 2 ///////////////////////////////////
[SP1.fit, SP1.S] = polyfit(p.freq,SP1.mean_gen,1);
[SP1.val, SP1.del] = polyval(SP1.fit,p.freq,SP1.S); clc
SP1.planck = SP1.fit * p.c; SP1.WF = SP1.planck(2); SP1.planck = SP1.planck(1);
%SP1.fit_yel = polyfit(p.freq,yellow.vol,1);
%SP1.val_yel = polyval(SP1.fit_yel,p.freq);

[SP2.fit, SP2.S] = polyfit(p.freq,SP2.mean_gen,1);
[SP2.val, SP2.del] = polyval(SP2.fit,p.freq,SP2.S); clc
SP2.planck = SP2.fit * p.c; SP2.WF = SP2.planck(2); SP2.planck = SP2.planck(1);
%SP2.fit_yel = polyfit(p.freq,yellow.vol,1);
%SP2.val_yel = polyval(SP2.fit_yel,p.freq);

%D = @(v1,v2) abs(v1 - v2)/((v1 + v2)/2) * 100;
%diff = [D(p.h,-SP1.planck(1)) D(p.h,-SP2.planck(1))];
%avg  = (SP1.planck(1) + SP2.planck(1))/2 ;

figure(3)
hold on
errorbar(p.freq,SP1.val,SP1.del,'sb')
plot(p.freq,SP1.mean_gen,'or')
plot(p.freq,SP1.val,'b')
errorbar(p.freq,SP1.mean_gen,SP1.std_gen,'xr')
hold off
xlabel('Frequency (Hz)');ylabel('Stopping Potential (V)');
str = sprintf('y = %.004ex + %.003f',SP1.fit(1),SP1.fit(2)); 
legend('Data',str)

figure(4)
hold on
errorbar(p.freq,SP2.val,SP2.del,'sb')
plot(p.freq,SP2.mean_gen,'dr')
plot(p.freq,SP2.val,'b')
errorbar(p.freq,SP2.mean_gen,SP2.std_gen,'xr')
hold off
xlabel('Frequency (Hz)');ylabel('Stopping Potential (V)');
str = sprintf('y = %.004ex + %.003f',SP2.fit(1),SP2.fit(2)); 
legend('Fit w/ error',str)


%______________________________________________________________
bint = SP1.S;
%SP1.herr = sqrt(diag((bint.R)\inv(bint.R'))./bint.normr.^2./bint.df);
SP1.herr = [-1*SP1.S.R(1) * 10^-30*p.c, SP1.S.R(3) * p.c];
bint = SP2.S;
%SP2.herr = p.c*sqrt(diag((bint.R)\inv(bint.R'))./bint.normr.^2./bint.df);
SP2.herr = [-1*SP2.S.R(1) * 10^-30*p.c, SP2.S.R(3) * p.c];


s = SP1.S; 
T1 = p.c*(sqrt(diag(inv(s.R)*inv(s.R'))./s.normr.^2./s.df));

s = SP2.S; 
T2 = p.c*(sqrt(diag(inv(s.R)*inv(s.R'))./s.normr.^2./s.df));

WMp = (SP1.planck*8.7134^(-35)) + (SP2.planck*7.3258^(-35))/(7.3258^(-35) + 8.7134^(-35));
WMw = abs((SP1.WF*SP1.herr(2) + SP2.WF*SP2.herr(2))/(SP1.herr(2) + SP2.herr(2)));
%Output
clc
fprintf('Solutions\n--------------------------------------------------------------------\n')
fprintf('Planck Estimate:%.06E , Work Function Estimate:%.06E\n',SP1.planck,SP1.WF)
fprintf('Planck Estimate:%.06E , Work Function Estimate:%.06E\n',SP2.planck,SP2.WF)
fprintf('Planck Estimate:%.06E , Work Function Estimate:%.06E\n',WMp,WMw)
fprintf('\n')
fprintf('Error(STD)\n--------------------------------------------------------------------\n')
fprintf('                 First Order       Second Order\n')
fprintf('Original Data:+- %.06E , +- %.06E\n',SP1.std_gen,SP2.std_gen) 
fprintf('Plancks  Fit :+- %.06E , +- %.06E\n',SP1.herr(1),SP2.herr(1)) 
fprintf('Work     Fit :+- %.06E , +- %.06E\n',SP1.herr(2),SP2.herr(2))

%%
clear str WW WM1 D avg ans wl WM diff