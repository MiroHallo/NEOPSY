% PLOT_DATA Plot input observed data (SPAC dispersion curves)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Read and plot input observed SPAC data from ASCII files
% NEOPSY: Stand-alone script for data preparation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Miroslav HALLO
% ETH Zürich, Swiss Seismological Service
% E-mail: miroslav.hallo@sed.ethz.ch
% Tested in Matlab R2025b
% Method:
% Hallo, M., Imperatori, W., Panzera, F., Fäh, D. (2021). Joint multizonal
%      transdimensional Bayesian inversion of surface wave dispersion and
%      ellipticity curves for local near-surface imaging, Geophys. J. Int.,
%      226 (1), 627-659. https://doi.org/10.1093/gji/ggab116
%
% Multizonal Transdimensional Inversion (NEOPSY)
% Version 2026/03: Fjord
%
% Copyright (C) 2019-2021 ETH Zurich
% This program is published under the GNU General Public License (GNU GPL).
%
% This program is free software: you can modify it and/or redistribute it
% or any derivative version under the terms of the GNU General Public
% License as published by the Free Software Foundation, either version 3
% of the License, or (at your option) any later version.
%
% This code is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY. We would like to kindly ask you to acknowledge the authors
% and don't remove their names from the code.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INIT:
close all;
clearvars;
addpath('../lib')
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:

% Prefix for data files
pref = '';

% Use modes? (1=YES, 0=NO)
%   R0, R1, R2, R3;
%   L0, L1, L2, L3;
%   ELL, ELA, -, -;
use = [
1 0 0 0
0 0 0 0
1 0 NaN NaN];

%  Data filenames
fnR0 = [pref,'Data_R0.dat'];
fnR1 = [pref,'Data_R1.dat'];
fnR2 = [pref,'Data_R2.dat'];
fnR3 = [pref,'Data_R3.dat'];
fnL0 = [pref,'Data_L0.dat'];
fnL1 = [pref,'Data_L1.dat'];
fnL2 = [pref,'Data_L2.dat'];
fnL3 = [pref,'Data_L3.dat'];
fnELL = [pref,'Data_ELL.dat'];
fnELA = [pref,'Data_ELA.dat'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% -------------------------------------------------------------------
% Read R and L
for w=1:2
    for m=1:4
        if use(w,m)==1
            
            % Select mode
            if w==1
                if m==1
                    fn = fnR0;
                    ftxt = 'R0';
                elseif m==2
                    fn = fnR1;
                    ftxt = 'R1';
                elseif m==3
                    fn = fnR2;
                    ftxt = 'R2';
                else
                    fn = fnR3;
                    ftxt = 'R3';
                end
            else
                if m==1
                    fn = fnL0;
                    ftxt = 'L0';
                elseif m==2
                    fn = fnL1;
                    ftxt = 'L1';
                elseif m==3
                    fn = fnL2;
                    ftxt = 'L2';
                else
                    fn = fnL3;
                    ftxt = 'L3';
                end
            end
            
            % Read data
            fid = fopen(fn,'r');
            tline = fgets(fid);
            count = 0;
            inData = [];
            while ~feof(fid)
                count = count+1;
                tline = fgets(fid);
                inData(count,1:3) = str2num(tline);
            end
            fclose(fid);
    
            % Convert to velocity
            inData_v = zeros(length(inData(:,3)),4);
            inData_v(:,1) = inData(:,1);
            inData_v(:,2) = 1./inData(:,2);
            inData_v(:,3) = (1./inData(:,2)) - (1./(inData(:,2)+inData(:,3)));
            inData_v(:,4) =  - (1./inData(:,2)) + (1./(inData(:,2)-inData(:,3)));
            
            % Find used data
            okf = inData(:,3)>0;

            % Find frequency axis limits and intelligent ticks
            [amin,amax,fticks] = fscale(min(inData(okf,1)),max(inData(okf,1)));
            frange = [amin, amax];
            
            % plot data
            fi = figure('Units','normalized','Color',[1 1 1],'OuterPosition',[0.2,0.1,0.6,0.8]);
            set(fi, 'PaperOrientation', 'portrait');
            set(fi, 'PaperPositionMode', 'manual');
            set(fi, 'PaperUnits', 'centimeters');
            set(fi, 'Render', 'painters');
            Px = 16;
            Py = 12;
            set(fi,'PaperSize', [Px Py]);
            set(fi,'PaperPosition', [0 0 Px Py]);
            
            subplot(2,1,1)
            set(gca,'XColor',[0 0 0],'YColor',[0 0 0])
            semilogx(inData(okf,1),inData(okf,2).*1000,'color','k','LineStyle',':'); hold on
            errorbar(inData(okf,1),inData(okf,2).*1000,inData(okf,3).*1000,'.','color','k','MarkerSize',5,'CapSize',2);
            hold off
            grid on
            box on
            set(gca,'Xtick',fticks)
            xlabel('Frequency (Hz)');
            ylabel('Slowness (ms/m)');
            title([ftxt,' mode'])
            legend('Data line','Data errors','location','best')
            xlim(frange)
            set(gca,'FontSize',8)
            
            subplot(2,1,2)
            set(gca,'XColor',[0 0 0],'YColor',[0 0 0])
            semilogx(inData_v(okf,1),inData_v(okf,2),'color',[0.1 0.1 0.9],'LineStyle',':'); hold on
            errorbar(inData_v(okf,1),inData_v(okf,2),inData_v(okf,3),inData_v(okf,4),'.','color',[0.0 0.0 0.9],'MarkerSize',5,'CapSize',2);
            hold off
            grid on
            box on
            set(gca,'Xtick',fticks)
            xlabel('Frequency (Hz)');
            ylabel('Velocity (m/s)');
            xlim(frange)
            set(gca,'FontSize',8)
            
            % Save to .pdf file
            saveas(fi,[pref,'Data_',ftxt,'.pdf'])
            
        end
    end
end


%% -------------------------------------------------------------------
% Read ELL
if use(3,1)==1
    fn = fnELL;
    ftxt = 'Rayleigh wave ellipticity';
    
    % read data
    fid = fopen(fn,'r');
    tline = fgets(fid);
    count = 0;
    inData = [];
    while ~feof(fid)
        count = count+1;
        tline = fgets(fid);
        inData(count,1:3) = str2num(tline);
    end
    fclose(fid);
    
    % Prepare for log-plot
    inData_log = zeros(length(inData(:,3)),4);
    inData_log(:,1) = inData(:,1);
    inData_log(:,2) = inData(:,2);
    inData_log(:,3) = inData(:,2)-(inData(:,2)./inData(:,3));
    inData_log(:,4) = -inData(:,2)+(inData(:,2).*inData(:,3));
    
    % Find used data
    okf = inData(:,3)>0;

    % Find frequency axis limits and intelligent ticks
    [amin,amax,fticks] = fscale(min(inData(okf,1)),max(inData(okf,1)));
    frange = [amin, amax];
    
    % plot data
    fi = figure('Units','normalized','Color',[1 1 1],'OuterPosition',[0.2,0.1,0.6,0.8]);
    set(fi, 'PaperOrientation', 'portrait');
    set(fi, 'PaperPositionMode', 'manual');
    set(fi, 'PaperUnits', 'centimeters');
    set(fi, 'Render', 'painters');
    Px = 16;
    Py = 12;
    set(fi,'PaperSize', [Px Py]);
    set(fi,'PaperPosition', [0 0 Px Py]);
    
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0])
    loglog(inData_log(okf,1),inData_log(okf,2),'color','k','LineStyle',':'); hold on;
    errorbar(inData_log(okf,1),inData_log(okf,2),inData_log(okf,3),inData_log(okf,4),'.','color','k','MarkerSize',5,'CapSize',2);
    hold off
    grid on
    box on
    set(gca,'Xtick',fticks)
    xlabel('Frequency (Hz)');
    ylabel('Ellipticity');
    title(ftxt)
    legend('Data line','Data errors','location','best')
    xlim(frange)
    set(gca,'FontSize',8)
    
    % Save to .pdf file
    saveas(fi,[pref,'Data_ELL.pdf'])
    
end


%% -------------------------------------------------------------------
% Read ELA
if use(3,2)==1
    fn = fnELA;
    ftxt = 'Rayleigh wave ellipticity angle';
    
    % read data
    fid = fopen(fn,'r');
    tline = fgets(fid);
    count = 0;
    inData = [];
    while ~feof(fid)
        count = count+1;
        tline = fgets(fid);
        inData(count,1:3) = str2num(tline);
    end
    fclose(fid);
    
    % Convert to Rayleigh wave ellipticity
    inData_e = zeros(length(inData(:,3)),4);
    inData_e(:,1) = inData(:,1);
    inData_e(:,2) = abs(tand(inData(:,2)));
    inData_e(:,3) = abs(tand(min(inData(:,2)+inData(:,3),90)));
    inData_e(:,4) = abs(tand(max(inData(:,2)-inData(:,3),-90)));
    for i=1:length(inData_e(:,1))
        inData_e(i,5) = inData_e(i,2) - min(inData_e(i,3),inData_e(i,4));
        inData_e(i,6) = - inData_e(i,2) + max(inData_e(i,3),inData_e(i,4));
    end
    inData_e(isinf(inData_e(:,6)),6) = 100;
    
    
    % Find used data
    okf = inData(:,3)>0;

    % Find frequency axis limits and intelligent ticks
    [amin,amax,fticks] = fscale(min(inData(okf,1)),max(inData(okf,1)));
    frange = [amin, amax];
    
    % plot data
    fi = figure('Units','normalized','Color',[1 1 1],'OuterPosition',[0.2,0.1,0.6,0.8]);
    set(fi, 'PaperOrientation', 'portrait');
    set(fi, 'PaperPositionMode', 'manual');
    set(fi, 'PaperUnits', 'centimeters');
    set(fi, 'Render', 'painters');
    Px = 16;
    Py = 12;
    set(fi,'PaperSize', [Px Py]);
    set(fi,'PaperPosition', [0 0 Px Py]);
    
    subplot(2,1,1)
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0])
    semilogx(inData(okf,1),inData(okf,2),'color','k','LineStyle',':'); hold on
    errorbar(inData(okf,1),inData(okf,2),inData(okf,3),'.','color','k','MarkerSize',5,'CapSize',2);
    hold off
    grid on
    box on
    set(gca,'Xtick',fticks)
    set(gca,'Ytick',-90:30:90)
    ylim([-90 90])
    xlabel('Frequency (Hz)');
    ylabel('Angle (deg)');
    title(ftxt)
    legend('Data line','Data errors','location','best')
    xlim(frange)
    set(gca,'FontSize',8)
    
    subplot(2,1,2)
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0])
    loglog(inData_e(okf,1),inData_e(okf,2),'color',[0.0 0.0 0.9],'LineStyle',':'); hold on
    errorbar(inData_e(okf,1),inData_e(okf,2),inData_e(okf,5),inData_e(okf,6),'.','color',[0.0 0.0 0.9],'MarkerSize',5,'CapSize',2);
    hold off
    grid on
    box on
    set(gca,'Xtick',fticks)
    xlabel('Frequency (Hz)');
    ylabel('Ellipticity');
    xlim(frange)
    set(gca,'FontSize',8)
    
    % Save to .pdf file
    saveas(fi,[pref,'Data_ELA.pdf'])
    
end

