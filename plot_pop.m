% PLOT_POP results of the Multizonal Transdimensional Inversion (NEOPSY)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot results of the Multizonal Transdimensional Inversion by NEOPSY code
% It reads inversion results and plots statistics (marginal histograms)
% from the ensemble of solutions produced by the "tires" code
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
%
% INIT:
close all;
clearvars;
projRoot = fileparts(which(mfilename));
addpath(fullfile(projRoot, 'lib'));
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:

invName = 'inv'; % NEOPSY working directory with inversion results
invPath = fullfile(projRoot, invName);

datafile =  fullfile(invPath, 'in_data.txt');     % Data of DC and ELL curves (file)
MAXprefix = fullfile(invPath, 'out_modelML');     % ML solution (prefix of input files)
MAPprefix = fullfile(invPath, 'out_modelMAP-SP'); % MAP solution (prefix of input files)
POPprefix = fullfile(invPath, 'out_pop');         % Ensemble statistics (prefix of input files)
OUTfolder = fullfile(invPath, 'results');         % Path to save output figures and results

REFfile = fullfile(projRoot, 'data/vs_ref_Swiss.ascii');  % Reference rock velocity model (file)

OUTprefix = 'PDF';                 % Prefix of output files (Site/Station/Code/Version)
Nfits = 300;                       % Number of random data fits from the solution ensemble to plot
plotDepthMAX = 9000;               % Maximal depth to be plotted [m]

syntPlot = 0;                      % Read and plot synthetic model (synthetic test; 1=YES, 0=NO)
syntfile = fullfile(projRoot, 'data/data.model');  % Target model (if syntPlot == 1)

plotInDPD = 1;                     % Plot ML and MAP models over PDF (1=YES, 0=NO)
plotLeg = 1;                       % Plot Legends (1=YES, 0=NO)
FigureRendering = 1;               % Render .pdf files from figures (0=NO, 1=YES, 2=YES and optimize)
useDataThr = 0;                    % Use joint axes for the data fit plots (1=YES, 0=NO)
SlowThr = [0.5 8.0];               % The joint slowness axes thresholds [ms/m] (only if useDataThr = 1)
FreqThr = [1.0 30.0];              % The joint freqency axes thresholds  [Hz] (only if useDataThr = 1)
maxQWLprop = 0.5;                  % Upper threshold of probability for the QWL plots
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create results folder
if ~exist(OUTfolder,'dir')
    mkdir(OUTfolder)
end
OUTprefix = fullfile(OUTfolder, OUTprefix);

%% -------------------------------------------------------------------
% Read POP headers
fid = fopen([POPprefix,'_headers.txt'],'r');
tline = fgets(fid);
tline = fgets(fid);
nBins = str2num(tline);
tline = fgets(fid);
nDepth = str2num(tline);
tline = fgets(fid);
models = str2num(tline);
tline = fgets(fid);
dummy = str2num(tline);
allBins = zeros(nBins,5);
for i=1:5
    tline = fgets(fid);
    allBins(1:nBins,i) = str2num(tline);
end
tline = fgets(fid);
dBins(1:nDepth) = str2num(tline);
tline = fgets(fid);
logdBins(1:nDepth) = str2num(tline);
fclose(fid);

% Prepare P1 vectors
dBinsP1 = [dBins, dBins(end)+(dBins(2)-dBins(1))];
allBinsP1 = zeros(nBins+1,4);
for i=1:4
    allBinsP1(1:nBins+1,i) = [allBins(1:nBins,i);allBins(end,i)+(allBins(2,i)-allBins(1,i))];
end

% N of models
display(['Number of sampling models: ',num2str(models)])


%% -------------------------------------------------------------------
% Read input data (DC and ELL curves)
fid = fopen(datafile,'r');
tline = fgets(fid);
tline = fgets(fid);
f_n = str2num(tline);
inData = zeros(max(f_n),4,length(f_n));
for m=1:length(f_n)
    if f_n(m)<=0
        continue
    end
    tline = fgets(fid);
    for i=1:f_n(m)
        tline = fgets(fid);
        inData(i,1:4,m) = str2num(tline);
    end
end
fclose(fid);

% Retrive data sigmas
inData_s = zeros(length(inData(:,1,1)),length(inData(1,1,:)));
for m=1:length(f_n)
    for i=1:f_n(m)
        inData_s(i,m) = 1/inData(i,3,m);
    end
end

% Find unused data
okf = false(max(f_n),length(f_n));
for m=1:length(f_n)
    okf(1:f_n(m),m) = inData(1:f_n(m),3,m)~=0;
end


%% -------------------------------------------------------------------
% Read MAX model data (DC and ELL curves)
try
    fid = fopen([MAXprefix,'_data.txt'],'r');
    tline = fgets(fid);
    tline = fgets(fid);
    f_n = str2num(tline);
    maxData = zeros(max(f_n),4,length(f_n));
    for m=1:length(f_n)
        if f_n(m)<=0
            continue
        end
        tline = fgets(fid);
        for i=1:f_n(m)
            tline = fgets(fid);
            maxData(i,1:4,m) = str2num(tline);
        end
    end
    fclose(fid);
catch
    % Prior PDF expert mode
    maxData=inData;
end

% Read MAP model data (DC and ELL curves)
try
    fid = fopen([MAPprefix,'_data.txt'],'r');
    tline = fgets(fid);
    tline = fgets(fid);
    f_n = str2num(tline);
    mapData = zeros(max(f_n),4,length(f_n));
    for m=1:length(f_n)
        if f_n(m)<=0
            continue
        end
        tline = fgets(fid);
        for i=1:f_n(m)
            tline = fgets(fid);
            mapData(i,1:4,m) = str2num(tline);
        end
    end
    fclose(fid);
catch
    % Prior PDF expert mode
    mapData=inData;
end


%% -------------------------------------------------------------------
% Read ML velocity model
fid = fopen([MAXprefix,'.txt'],'r');
tline = fgets(fid);
nL = str2num(tline);
tline = fgets(fid);
VR = str2num(tline);
tmpMod = zeros(nL,5);
for i=1:5
    tline = fgets(fid);
    tmpMod(1:nL,i) = str2num(tline);
end
fclose(fid);
NnL_ML = nL;

% Prepare the velocity model for plots
count = 0;
maxMod = zeros(2*nL,5);
for i=1:nL-1
    count = count + 1;
    maxMod(count,1:5) = tmpMod(i,1:5);
    count = count + 1;
    maxMod(count,1:5) = maxMod(count-1,1:5);
    if count==2
        maxMod(count-1,1) = 0;
    else
        maxMod(count-1,1) = maxMod(count-2,1);
        maxMod(count,1) = maxMod(count-1,1) + maxMod(count,1);
    end
end
count = count + 1;
maxMod(count,1:5) = tmpMod(nL,1:5);
count = count + 1;
maxMod(count,1:5) = maxMod(count-1,1:5);
if(count>2)
    maxMod(count-1,1) = maxMod(count-2,1);
end
maxMod(count,1) = plotDepthMAX;
% Swith Vp and Vs to be consistent
tmpFK = maxMod(:,2);
maxMod(:,2) = maxMod(:,3);
maxMod(:,3) = tmpFK;


% Read MAP velocity model
fid = fopen([MAPprefix,'.txt'],'r');
tline = fgets(fid);
nL = str2num(tline);
tline = fgets(fid);
VR2 = str2num(tline);
tmpMod2 = zeros(nL,5);
for i=1:5
    tline = fgets(fid);
    tmpMod2(1:nL,i) = str2num(tline);
end
fclose(fid);
NnL_MAP = nL;

% Prepare the velocity model for plots
count = 0;
mapMod = zeros(2*nL,5);
for i=1:nL-1
    count = count + 1;
    mapMod(count,1:5) = tmpMod2(i,1:5);
    count = count + 1;
    mapMod(count,1:5) = mapMod(count-1,1:5);
    if count==2
        mapMod(count-1,1) = 0;
    else
        mapMod(count-1,1) = mapMod(count-2,1);
        mapMod(count,1) = mapMod(count-1,1) + mapMod(count,1);
    end
end
count = count + 1;
mapMod(count,1:5) = tmpMod2(nL,1:5);
count = count + 1;
mapMod(count,1:5) = mapMod(count-1,1:5);
if(count>2)
    mapMod(count-1,1) = mapMod(count-2,1);
end
mapMod(count,1) = plotDepthMAX;
% Swith Vp and Vs to be consistent
tmpFK = mapMod(:,2);
mapMod(:,2) = mapMod(:,3);
mapMod(:,3) = tmpFK;


%% -------------------------------------------------------------------
% Read POP data binary ensemble
count = 0;
binData = zeros(ceil(models/300),2,sum(f_n(:)));
fid = fopen([POPprefix,'_data.bin'],'r');
while ~feof(fid)
    d_tmp = fread(fid,[1,sum(f_n(:))],'double');
    dw_tmp = fread(fid,[1,sum(f_n(:))],'double');
    if ~isempty(d_tmp) && ~isempty(dw_tmp)
        count = count+1;
        binData(count,1,1:sum(f_n(:))) = d_tmp;
        binData(count,2,1:sum(f_n(:))) = dw_tmp;
    end
end
fclose(fid);

Nfits = max(1,min(Nfits,count));
synData = zeros(Nfits,2,sum(f_n(:)));
j = 0;
for i=1:count
    if rem(i,max(1,floor(count/Nfits)))==0
        j = j+1;
        synData(j,1:2,1:sum(f_n(:))) = binData(i,1:2,1:sum(f_n(:)));
    end
end


%% -------------------------------------------------------------------
% Plot data fit (slowness)

maxMode = length(f_n);
maxMode_a = sum(f_n(:)>0);

if maxMode_a<4
    Px = maxMode_a*7;
    Py = 7;
    Plin = 1;
elseif maxMode_a<7
    Px = 21;
    Py = 14;
    Plin = 2;
else
    Px = 21;
    Py = 21;
    Plin = 3;
end
fi = figure('Units','centimeters','Position',[2, 2, Px, Py],'Color',[1 1 1]);
set(fi, 'PaperOrientation', 'portrait');
set(fi, 'PaperPositionMode', 'manual');
set(fi, 'PaperUnits', 'centimeters');
set(fi, 'Render', 'painters');
set(fi,'PaperSize', [Px Py]);
set(fi,'PaperPosition', [0 0 Px Py]);

m_a = 0;
for m=1:maxMode
    if f_n(m)>0
        m_a = m_a+1;
    else
        continue
    end
    if m==1
        ftxt = 'Rayleigh (fundamental)';
    elseif m==2
        ftxt = 'Rayleigh (1st higher)';
    elseif m==3
        ftxt = 'Rayleigh (2nd higher)';
    elseif m==4
        ftxt = 'Rayleigh (3rd higher)';
    elseif m==5
        ftxt = 'Love (fundamental)';
    elseif m==6
        ftxt = 'Love (1st higher)';
    elseif m==7
        ftxt = 'Love (2nd higher)';
    elseif m==8
        ftxt = 'Love (3rd higher)';
    elseif m==9
        ftxt = 'Rayleigh wave ellipticity';
    elseif m==10
        ftxt = 'Rayleigh wave ellipticity angle';
    else
        continue
    end

    % Find frequency axis limits and intelligent ticks
    [amin,amax,fticks] = fscale(min(inData(okf(1:f_n(m),m),1,m)),max(inData(okf(1:f_n(m),m),1,m)));
    frange = [amin, amax];

    ci = sum(f_n(1:m-1));

    subplot(Plin,min(maxMode_a,3),m_a)
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0])
    if m<=8
        for i=1:Nfits
            lh(5)=semilogx(inData(1:f_n(m),1,m),squeeze(synData(i,1,ci+1:ci+f_n(m)).*1000),'color',[0.7 0.7 0.8]); hold on
        end
        lh(1)=semilogx(inData(okf(1:f_n(m),m),1,m),inData(okf(1:f_n(m),m),2,m).*1000,'color','k','LineStyle',':');
        lh(2)=errorbar(inData(okf(1:f_n(m),m),1,m),inData(okf(1:f_n(m),m),2,m).*1000,inData_s(okf(1:f_n(m),m),m).*1000,'.','color','k','MarkerSize',5,'CapSize',2);
        lh(3)=semilogx(maxData(1:f_n(m),1,m),maxData(1:f_n(m),2,m).*1000,'color','b');
        lh(4)=semilogx(mapData(1:f_n(m),1,m),mapData(1:f_n(m),2,m).*1000,'color','m');
        ylabel('Slowness (ms/m)');
        if useDataThr == 1
            ylim(SlowThr)
            xlim(FreqThr)
        else
            xlim(frange)
            set(gca,'Xtick',fticks)
        end
    elseif m==9
        for i=1:Nfits
            lh(5)=loglog(inData(1:f_n(m),1,m),10.^squeeze(synData(i,1,ci+1:ci+f_n(m))),'color',[0.7 0.7 0.8]); hold on
        end
        lh(1)=loglog(inData(okf(1:f_n(m),m),1,m),10.^inData(okf(1:f_n(m),m),2,m),'color','k','LineStyle',':');
        lh(2)=errorbar(inData(okf(1:f_n(m),m),1,m),10.^inData(okf(1:f_n(m),m),2,m),...
            10.^inData(okf(1:f_n(m),m),2,m)-10.^(inData(okf(1:f_n(m),m),2,m)-inData_s(okf(1:f_n(m),m),m)),...
            -10.^inData(okf(1:f_n(m),m),2,m)+10.^(inData(okf(1:f_n(m),m),2,m)+inData_s(okf(1:f_n(m),m),m)),'.','color','k','MarkerSize',5,'CapSize',2);
        lh(3)=loglog(maxData(1:f_n(m),1,m),10.^maxData(1:f_n(m),2,m),'color','b');
        lh(4)=loglog(mapData(1:f_n(m),1,m),10.^mapData(1:f_n(m),2,m),'color','m');
        ylabel('Ellipticity');
        ylim([0.2 100])
        xlim(frange)
        set(gca,'Xtick',fticks)
    elseif m==10
        for i=1:Nfits
            tmp = squeeze(synData(i,1,ci+1:ci+f_n(m)));
            ind = find(abs(tmp(2:f_n(m)) - tmp(1:f_n(m)-1))>90,1);
            if isempty(ind), ind = 0; end
            lh(5)=semilogx([inData(1:ind,1,m);NaN;inData(ind+1:f_n(m),1,m)],[tmp(1:ind);NaN;tmp(ind+1:f_n(m))],'color',[0.7 0.7 0.8]); hold on
        end
        ind = find(abs(inData(okf(2:f_n(m),m),2,m) - inData(okf(1:f_n(m)-1,m),2,m))>90,1);
        if isempty(ind), ind = 0; end
        lh(1)=semilogx([inData(okf(1:ind,m),1,m);NaN;inData(okf(ind+1:f_n(m),m),1,m)],[inData(okf(1:ind,m),2,m);NaN;inData(okf(ind+1:f_n(m),m),2,m)],'color','k','LineStyle',':');
        lh(2)=errorbar(inData(okf(1:f_n(m),m),1,m),inData(okf(1:f_n(m),m),2,m),inData_s(okf(1:f_n(m),m),m),'.','color','k','MarkerSize',5,'CapSize',2);
        ind = find(abs(maxData(2:f_n(m),2,m) - maxData(1:f_n(m)-1,2,m))>90,1);
        if isempty(ind), ind = 0; end
        lh(3)=semilogx([maxData(1:ind,1,m);NaN;maxData(ind+1:f_n(m),1,m)],[maxData(1:ind,2,m);NaN;maxData(ind+1:f_n(m),2,m)],'color','b');
        ind = find(abs(mapData(2:f_n(m),2,m) - mapData(1:f_n(m)-1,2,m))>90,1);
        if isempty(ind), ind = 0; end
        lh(4)=semilogx([mapData(1:ind,1,m);NaN;mapData(ind+1:f_n(m),1,m)],[mapData(1:ind,2,m);NaN;mapData(ind+1:f_n(m),2,m)],'color','m');
        set(gca,'ytick',-90:30:90)
        ylabel('Angle (deg)');
        ylim([-90 90])
        xlim(frange)
        set(gca,'Xtick',fticks)
    end
    title(ftxt)
    hold off
    grid on
    box on
    xlabel('Frequency (Hz)');
    if m_a==1
        if plotLeg == 1
            legend(lh,{'Data','Data errors','ML model','MAP model','Predictive dist.'},'Location','northwest')
        end
    end
    set(gca,'FontSize',8)
end
set(findall(gcf,'type','text'),'fontSize',8)
drawnow

% Save to .pdf file
if (FigureRendering>0)
    % saveas(fi,[OUTprefix,'_','fit.pdf'])
    exportgraphics(fi,[OUTprefix,'_','fit_slowness.pdf'], 'Resolution', 300);
    exportgraphics(fi,[OUTprefix,'_','fit_slowness.png'], 'Resolution', 300);
end
if (FigureRendering>1)
    close(fi)
end


%% -------------------------------------------------------------------
% Plot data fit (velocity)

maxMode = length(f_n);
maxMode_a = sum(f_n(:)>0);

if maxMode_a<4
    Px = maxMode_a*7;
    Py = 7;
    Plin = 1;
elseif maxMode_a<7
    Px = 21;
    Py = 14;
    Plin = 2;
else
    Px = 21;
    Py = 21;
    Plin = 3;
end
fi = figure('Units','centimeters','Position',[2, 2, Px, Py],'Color',[1 1 1]);
set(fi, 'PaperOrientation', 'portrait');
set(fi, 'PaperPositionMode', 'manual');
set(fi, 'PaperUnits', 'centimeters');
set(fi, 'Render', 'painters');
set(fi,'PaperSize', [Px Py]);
set(fi,'PaperPosition', [0 0 Px Py]);

m_a = 0;
for m=1:maxMode
    if f_n(m)>0
        m_a = m_a+1;
    else
        continue
    end
    if m==1
        ftxt = 'Rayleigh (fundamental)';
    elseif m==2
        ftxt = 'Rayleigh (1st higher)';
    elseif m==3
        ftxt = 'Rayleigh (2nd higher)';
    elseif m==4
        ftxt = 'Rayleigh (3rd higher)';
    elseif m==5
        ftxt = 'Love (fundamental)';
    elseif m==6
        ftxt = 'Love (1st higher)';
    elseif m==7
        ftxt = 'Love (2nd higher)';
    elseif m==8
        ftxt = 'Love (3rd higher)';
    elseif m==9
        ftxt = 'Rayleigh wave ellipticity';
    elseif m==10
        ftxt = 'Rayleigh wave ellipticity angle';
    else
        continue
    end

    % Find frequency axis limits and intelligent ticks
    [amin,amax,fticks] = fscale(min(inData(okf(1:f_n(m),m),1,m)),max(inData(okf(1:f_n(m),m),1,m)));
    frange = [amin, amax];

    ci = sum(f_n(1:m-1));

    subplot(Plin,min(maxMode_a,3),m_a)
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0])
    if m<=8
        for i=1:Nfits
            lh(5)=semilogx(inData(1:f_n(m),1,m),squeeze(1./synData(i,1,ci+1:ci+f_n(m))),'color',[0.7 0.7 0.8]); hold on
        end
        lh(1)=semilogx(inData(okf(1:f_n(m),m),1,m),1./inData(okf(1:f_n(m),m),2,m),'color','k','LineStyle',':');
        inData_vs1 = (1./inData(okf(1:f_n(m),m),2,m)) - (1./(inData(okf(1:f_n(m),m),2,m)+inData_s(okf(1:f_n(m),m),m)));
        inData_vs2 = - (1./inData(okf(1:f_n(m),m),2,m)) + (1./(inData(okf(1:f_n(m),m),2,m)-inData_s(okf(1:f_n(m),m),m)));
        lh(2)=errorbar(inData(okf(1:f_n(m),m),1,m),1./inData(okf(1:f_n(m),m),2,m),inData_vs1,inData_vs2,'.','color','k','MarkerSize',5,'CapSize',2);
        lh(3)=semilogx(maxData(1:f_n(m),1,m),1./maxData(1:f_n(m),2,m),'color','b');
        lh(4)=semilogx(mapData(1:f_n(m),1,m),1./mapData(1:f_n(m),2,m),'color','m');
        ylabel('Velocity (m/s)');
        if useDataThr == 1
            ylim([1000/SlowThr(2),1000/SlowThr(1)])
            xlim(FreqThr)
        else
            xlim(frange)
            set(gca,'Xtick',fticks)
        end
    elseif m==9
        for i=1:Nfits
            lh(5)=loglog(inData(1:f_n(m),1,m),10.^squeeze(synData(i,1,ci+1:ci+f_n(m))),'color',[0.7 0.7 0.8]); hold on
        end
        lh(1)=loglog(inData(okf(1:f_n(m),m),1,m),10.^inData(okf(1:f_n(m),m),2,m),'color','k','LineStyle',':');
        lh(2)=errorbar(inData(okf(1:f_n(m),m),1,m),10.^inData(okf(1:f_n(m),m),2,m),...
            10.^inData(okf(1:f_n(m),m),2,m)-10.^(inData(okf(1:f_n(m),m),2,m)-inData_s(okf(1:f_n(m),m),m)),...
            -10.^inData(okf(1:f_n(m),m),2,m)+10.^(inData(okf(1:f_n(m),m),2,m)+inData_s(okf(1:f_n(m),m),m)),'.','color','k','MarkerSize',5,'CapSize',2);
        lh(3)=loglog(maxData(1:f_n(m),1,m),10.^maxData(1:f_n(m),2,m),'color','b');
        lh(4)=loglog(mapData(1:f_n(m),1,m),10.^mapData(1:f_n(m),2,m),'color','m');
        ylabel('Ellipticity');
        ylim([0.2 100])
        xlim(frange)
        set(gca,'Xtick',fticks)
    elseif m==10
        for i=1:Nfits
            tmp = squeeze(synData(i,1,ci+1:ci+f_n(m)));
            ind = find(abs(tmp(2:f_n(m)) - tmp(1:f_n(m)-1))>90,1);
            if isempty(ind), ind = 0; end
            lh(5)=semilogx([inData(1:ind,1,m);NaN;inData(ind+1:f_n(m),1,m)],[tmp(1:ind);NaN;tmp(ind+1:f_n(m))],'color',[0.7 0.7 0.8]); hold on
        end
        ind = find(abs(inData(okf(2:f_n(m),m),2,m) - inData(okf(1:f_n(m)-1,m),2,m))>90,1);
        if isempty(ind), ind = 0; end
        lh(1)=semilogx([inData(okf(1:ind,m),1,m);NaN;inData(okf(ind+1:f_n(m),m),1,m)],[inData(okf(1:ind,m),2,m);NaN;inData(okf(ind+1:f_n(m),m),2,m)],'color','k','LineStyle',':');
        lh(2)=errorbar(inData(okf(1:f_n(m),m),1,m),inData(okf(1:f_n(m),m),2,m),inData_s(okf(1:f_n(m),m),m),'.','color','k','MarkerSize',5,'CapSize',2);
        ind = find(abs(maxData(2:f_n(m),2,m) - maxData(1:f_n(m)-1,2,m))>90,1);
        if isempty(ind), ind = 0; end
        lh(3)=semilogx([maxData(1:ind,1,m);NaN;maxData(ind+1:f_n(m),1,m)],[maxData(1:ind,2,m);NaN;maxData(ind+1:f_n(m),2,m)],'color','b');
        ind = find(abs(mapData(2:f_n(m),2,m) - mapData(1:f_n(m)-1,2,m))>90,1);
        if isempty(ind), ind = 0; end
        lh(4)=semilogx([mapData(1:ind,1,m);NaN;mapData(ind+1:f_n(m),1,m)],[mapData(1:ind,2,m);NaN;mapData(ind+1:f_n(m),2,m)],'color','m');
        set(gca,'ytick',-90:30:90)
        ylabel('Angle (deg)');
        ylim([-90 90])
        xlim(frange)
        set(gca,'Xtick',fticks)
    end
    title(ftxt)
    hold off
    grid on
    box on
    xlabel('Frequency (Hz)');
    if m_a==1
        if plotLeg == 1
            legend(lh,{'Data','Data errors','ML model','MAP model','Predictive dist.'},'Location','northeast')
        end
    end
    set(gca,'FontSize',8)
end
set(findall(gcf,'type','text'),'fontSize',8)
drawnow

% Save to .pdf file
if (FigureRendering>0)
    % saveas(fi,[OUTprefix,'_','fitV.pdf'])
    exportgraphics(fi,[OUTprefix,'_','fit_velocity.pdf'], 'Resolution', 300);
    exportgraphics(fi,[OUTprefix,'_','fit_velocity.png'], 'Resolution', 300);
end
if (FigureRendering>1)
    close(fi)
end


%% -------------------------------------------------------------------
% Plot Standardized data misfit

maxMode = length(f_n);
maxMode_a = sum(f_n(:)>0);

if maxMode_a<4
    Px = maxMode_a*7;
    Py = 7;
    Plin = 1;
elseif maxMode_a<7
    Px = 21;
    Py = 14;
    Plin = 2;
else
    Px = 21;
    Py = 21;
    Plin = 3;
end
fi = figure('Units','centimeters','Position',[2, 2, Px, Py],'Color',[1 1 1]);
set(fi, 'PaperOrientation', 'portrait');
set(fi, 'PaperPositionMode', 'manual');
set(fi, 'PaperUnits', 'centimeters');
set(fi, 'Render', 'painters');
set(fi,'PaperSize', [Px Py]);
set(fi,'PaperPosition', [0 0 Px Py]);

m_a = 0;
for m=1:maxMode
    if f_n(m)>0
        m_a = m_a+1;
    else
        continue
    end
    if m==1
        ftxt = 'Rayleigh (fundamental)';
    elseif m==2
        ftxt = 'Rayleigh (1st higher)';
    elseif m==3
        ftxt = 'Rayleigh (2nd higher)';
    elseif m==4
        ftxt = 'Rayleigh (3rd higher)';
    elseif m==5
        ftxt = 'Love (fundamental)';
    elseif m==6
        ftxt = 'Love (1st higher)';
    elseif m==7
        ftxt = 'Love (2nd higher)';
    elseif m==8
        ftxt = 'Love (3rd higher)';
    elseif m==9
        ftxt = 'Rayleigh wave ellipticity';
    elseif m==10
        ftxt = 'Rayleigh wave ellipticity angle';
    else
        continue
    end
    % Find frequency axis limits and intelligent ticks
    [amin,amax,fticks] = fscale(min(inData(okf(1:f_n(m),m),1,m)),max(inData(okf(1:f_n(m),m),1,m)));
    frange = [amin, amax];

    ci = sum(f_n(1:m-1));

    axh(m_a) = subplot(Plin,min(maxMode_a,3),m_a);
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0])
    yyaxis left
    norm = max(inData(okf(1:f_n(m),m),3,m));
    for i=1:Nfits
        tmpData = squeeze(synData(i,2,ci+1:ci+f_n(m)));
        lh2(5)=semilogx(inData(okf(1:f_n(m),m),1,m),tmpData(okf(1:f_n(m),m)),'-','color',[0.7 0.7 0.8]); hold on
    end
    yyaxis right
    lh2(2)=semilogx(inData(okf(1:f_n(m),m),1,m),inData(okf(1:f_n(m),m),3,m)/norm,'-','color','g');
    yyaxis left
    lh2(1)=semilogx(inData(okf(1:f_n(m),m),1,m),inData(okf(1:f_n(m),m),4,m),'o','MarkerSize',5,'color','r');
    lh2(3)=semilogx(maxData(okf(1:f_n(m),m),1,m),maxData(okf(1:f_n(m),m),4,m),'.-','color','b','MarkerSize',5);
    lh2(4)=semilogx(mapData(okf(1:f_n(m),m),1,m),mapData(okf(1:f_n(m),m),4,m),'.-','color','m','MarkerSize',5);
    hold off
    grid on
    box on
    xlabel('Frequency (Hz)');
    ylabel('Standardized units');
    xlim(frange)
    set(gca,'Xtick',fticks)
    title(ftxt)
    if m_a==1
        if plotLeg == 1
            legend(lh2,{'Data','Data weights','ML model','MAP model','Predictive dist.'},'Location','southwest')
        end
    end
    yyaxis right
    ylabel('min(\sigma) / \sigma');
    set(gca,'Ycolor','g')
    set(gca,'Ytick',[0 1])
    set(gca,'Ylim',[0 1])
    yyaxis left
    set(gca,'Ycolor','k')
    set(gca,'FontSize',8)
end
linkaxes(axh,'y')
set(findall(gcf,'type','text'),'fontSize',8)
drawnow

% Save to .pdf file
if (FigureRendering>0)
    % saveas(fi,[OUTprefix,'_','fitS.pdf'])
    exportgraphics(fi,[OUTprefix,'_','fit_weight.pdf'], 'Resolution', 300);
    exportgraphics(fi,[OUTprefix,'_','fit_weight.png'], 'Resolution', 300);
end
if (FigureRendering>1)
    close(fi)
end


%% -------------------------------------------------------------------
% Read synthetic model
if syntPlot == 1
    fid = fopen(syntfile,'r');
    tline = fgets(fid);
    tline = fgets(fid);
    nL = str2num(tline);
    tline = fgets(fid);
    count = 0;
    synMod = zeros(2*nL,5);
    for i=1:nL-1
        tline = fgets(fid);
        count = count + 1;
        synMod(count,1:4) = str2num(tline);
        nu = ((synMod(count,2)^2) - 2*(synMod(count,3)^2)) / (2*((synMod(count,2)^2) - (synMod(count,3)^2)));
        synMod(count,5) = nu;
        count = count + 1;
        synMod(count,1:4) = synMod(count-1,1:4);
        if count==2
            synMod(count-1,1) = 0;
        else
            synMod(count-1,1) = synMod(count-2,1);
            synMod(count,1) = synMod(count-1,1) + synMod(count,1);
        end
        nu = ((synMod(count,2)^2) - 2*(synMod(count,3)^2)) / (2*((synMod(count,2)^2) - (synMod(count,3)^2)));
        synMod(count,5) = nu;
    end
    tline = fgets(fid);
    tline = fgets(fid);
    count = count + 1;
    synMod(count,1:4) = str2num(tline);
    nu = ((synMod(count,2)^2) - 2*(synMod(count,3)^2)) / (2*((synMod(count,2)^2) - (synMod(count,3)^2)));
    synMod(count,5) = nu;
    count = count + 1;
    synMod(count,1:4) = synMod(count-1,1:4);
    synMod(count-1,1) = synMod(count-2,1);
    nu = ((synMod(count,2)^2) - 2*(synMod(count,3)^2)) / (2*((synMod(count,2)^2) - (synMod(count,3)^2)));
    synMod(count,5) = nu;
    synMod(count,1) = plotDepthMAX;
    fclose(fid);
end


%% -------------------------------------------------------------------
% POP basic histograms (layer interfaces and VR)

fid = fopen([POPprefix,'_dep1D.txt'],'r');
pop_1D = zeros(nDepth,1);
for i=1:nDepth
    tline = fgets(fid);
    pop_1D(i,1) = str2num(tline);
end
fclose(fid);

fid = fopen([POPprefix,'_lay1D.txt'],'r');
tline = fgets(fid);
pop_lay = str2num(tline);
fclose(fid);

fid = fopen([POPprefix,'_vr1D.txt'],'r');
tline = fgets(fid);
pop_vr = str2num(tline);
fclose(fid);

dBinsW = (dBins(2)-dBins(1));
dBinsC = dBins(1:nDepth)+(dBinsW/2);

logdBinsW = (log(logdBins(3))-log(logdBins(2)));
logdBinsC = log(logdBins(2:nDepth))+(logdBinsW/2);

vrBinsW = (allBins(2,5)-allBins(1,5));
vrBinsC = allBins(1:nBins,5)+(vrBinsW/2);

% Find interfaces
smoo_data = smoothdata(pop_1D(2:end),'sgolay',9);
smoo_am = mean(smoo_data);
[pks,locs] = findpeaks(smoo_data,'MinPeakDistance',9,'MinPeakHeight',1.5*smoo_am);

% Save into text file
fid = fopen([OUTprefix,'_','layers.ascii'],'w');
fprintf(fid,'%s\r\n','# The most probable depths of layer interfaces');
fprintf(fid,'%s\r\n','# Depth[m], Significance(large number = more significant interface)');
for k=1:length(locs)
    fprintf(fid,'%8.2f %6.2f\r\n',exp(logdBinsC(locs(k))),pks(k)/smoo_am);
end
fclose(fid);

% Plot
Px = 21;
Py = 14;
fi = figure('Units','centimeters','Position',[2, 2, Px, Py],'Color',[1 1 1]);
set(fi, 'PaperOrientation', 'portrait');
set(fi, 'PaperPositionMode', 'manual');
set(fi, 'PaperUnits', 'centimeters');
set(fi, 'Render', 'painters');
set(fi,'PaperSize', [Px Py]);
set(fi,'PaperPosition', [0 0 Px Py]);

subplot(2,3,1)
set(gca,'XColor',[0 0 0],'YColor',[0 0 0])
set(gca, 'OuterPosition', [0.02, 0.025, 0.3, 0.97]);
barh(logdBinsC,pop_1D(2:end)./models,1.0, 'FaceColor',[0.4 0.4 0.4], 'EdgeColor', [0.4 0.4 0.4])
hold on
for k=1:length(locs)
    plot(0.05*max(get(gca,'XLim')),logdBinsC(locs(k)),'Color','g','Marker','+')
end
hold off
set(gca,'ydir','reverse');
set(gca,'ylim',[log(logdBins(2)) log(logdBins(nDepth))+logdBinsW])
box on
grid on;
xlabel('Probability');
ylabel('ln(depth)');
title('Interfaces')
set(gca,'FontSize',8)

subplot(2,3,2)
set(gca,'XColor',[0 0 0],'YColor',[0 0 0])
set(gca, 'OuterPosition', [0.33, 0.025, 0.3, 0.97]);
hold on
for i=2:length(pop_1D)-1
    rectangle('position',[0 logdBins(i) pop_1D(i)/models (logdBins(i+1)-logdBins(i))], 'FaceColor',[0.4 0.4 0.4],'EdgeColor',[0.4 0.4 0.4])
end
for k=1:length(locs)
    plot(0.05*max(get(gca,'XLim')),exp(logdBinsC(locs(k))),'Color','g','Marker','+')
end
hold off
set(gca,'ydir','reverse');
set(gca,'ylim',[0 max(dBins)+dBinsW])
box on
grid on;
xlabel('Probability');
ylabel('Depth (m)');
title('Interfaces')
set(gca,'FontSize',8)

subplot(2,3,3)
set(gca,'XColor',[0 0 0],'YColor',[0 0 0])
bar(1:length(pop_lay),pop_lay./models,1.0, 'FaceColor',[1 0.4 0.6])
hold on;
text(length(pop_lay),max(get(gca,'ylim')),{[' ML: ',num2str(NnL_ML),' layers'],[' MAP: ',num2str(NnL_MAP),' layers']},'Color','k','HorizontalAlignment','right','VerticalAlignment','cap','FontSize',8)
hold off;
set(gca,'xlim',[0.5 length(pop_lay)+0.5])
box on
grid on;
xlabel('{\itk} layers')
ylabel('Probability')
title('Number of layers')
set(gca,'FontSize',8)

subplot(2,3,6)
set(gca,'XColor',[0 0 0],'YColor',[0 0 0])
bar(vrBinsC,pop_vr,1.0, 'FaceColor',[0.4 0.6 0.6],'EdgeColor',[0.4 0.6 0.6]); hold on;
yl = get(gca,'ylim');
plot([0,0],yl,'-','color','k');
plot([43.75,43.75],yl,'-','color','k');
plot([75,75],yl,'-','color','k');
plot([93.75,93.75],yl,'-','color','k');
text(1,yl(2),'Fit (within data errors) ','Color','k','Rotation',90,'HorizontalAlignment','right','VerticalAlignment','top','FontSize',8)
text(43.75,yl(2),'Fair fit ','Color','k','Rotation',90,'HorizontalAlignment','right','VerticalAlignment','top','FontSize',8)
text(75,yl(2),'Good fit ','Color','k','Rotation',90,'HorizontalAlignment','right','VerticalAlignment','top','FontSize',8)
text(93,yl(2),'Perfect fit ','Color','k','Rotation',90,'HorizontalAlignment','right','VerticalAlignment','cap','FontSize',8)
text(1,0,{[' ML: ',num2str(round(VR)),'%'],[' MAP: ',num2str(round(VR2)),'%']},'Color','k','HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',8)
hold off;
set(gca,'ylim',yl)
set(gca,'xlim',[min(allBins(1:nBins,5)) max(allBins(1:nBins,5))+vrBinsW])
set(gca,'xTick',[0 43.75 75 93.75])
box on
grid on;
xlabel('(%)')
ylabel('Number of models')
title('Data variance reduction')
set(gca,'FontSize',8)

set(findall(gcf,'type','text'),'fontSize',8)
drawnow

% Save to .pdf file
if (FigureRendering>0)
    % saveas(fi,[OUTprefix,'_lay.pdf'])
    exportgraphics(fi,[OUTprefix,'_layers.pdf'], 'Resolution', 300);
    exportgraphics(fi,[OUTprefix,'_layers.png'], 'Resolution', 300);
end
if (FigureRendering>1)
    close(fi)
end


%% -------------------------------------------------------------------
% clear from previous before ploting PDFs
try
    clear binData synData;
    clear inData mapData maxData;
catch
end


%% -------------------------------------------------------------------
% Pop statistics

nmi = zeros(1,4);

for ptype=1:4
    
    pop_2D = zeros(nDepth,nBins);
    dBins2 = zeros(2*nDepth,1);
    MaxMeaSig = zeros(2*nDepth,4);
    
    % Read pop statistics
    if ptype==1
        fid = fopen([POPprefix,'_vs2D.txt'],'r');
        tmpBins = allBins(:,1) + (allBins(2,1)-allBins(1,1))/2;
        tmpBinsP1 = allBinsP1(:,1);
        textP = '{\itV}_S';
        textU = '(m/s)';
        textOUT = 'vs';
        textAM = 'AM of PDF';
        maxModX = maxMod(:,3);
        mapModX = mapMod(:,3);
        if syntPlot == 1
            synModX = synMod(:,3);
        end
    elseif ptype==2
        fid = fopen([POPprefix,'_vp2D.txt'],'r');
        tmpBins = allBins(:,2) + (allBins(2,2)-allBins(1,2))/2;
        tmpBinsP1 = allBinsP1(:,2);
        textP = '{\itV}_P';
        textU = '(m/s)';
        textOUT = 'vp';
        textAM = 'AM of PDF';
        maxModX = maxMod(:,2);
        mapModX = mapMod(:,2);
        if syntPlot == 1
            synModX = synMod(:,2);
        end
    elseif ptype==3
        fid = fopen([POPprefix,'_nu2D.txt'],'r');
        tmpBins = allBins(:,3) + (allBins(2,3)-allBins(1,3))/2;
        tmpBinsP1 = allBinsP1(:,3);
        textP = '\nu';
        textU = '';
        textOUT = 'nu';
        textAM = 'AM of PDF';
        maxModX = maxMod(:,5);
        mapModX = mapMod(:,5);
        if syntPlot == 1
            synModX = synMod(:,5);
        end
    else
        fid = fopen([POPprefix,'_rho2D.txt'],'r');
        tmpBins = allBins(:,4) + (allBins(2,4)-allBins(1,4))/2;
        tmpBinsP1 = allBinsP1(:,4);
        textP = '\rho';
        textU = '(kg/m^3)';
        textOUT = 'rho';
        textAM = 'AM of PDF';
        maxModX = maxMod(:,4);
        mapModX = mapMod(:,4);
        if syntPlot == 1
            synModX = synMod(:,4);
        end
    end
    for i=1:nDepth
        tline = fgets(fid);
        pop_2D(i,1:nBins) = str2num(tline);
    end
    fclose(fid);
    
    % Skip fixed parameters
    if abs(tmpBins(end)-tmpBins(1))<0.0000000000001
        continue
    end
    
    % velocity to slowness
    if ptype==1 || ptype==2 % Vs or Vp
        tmpBins = 1./tmpBins(:);
    end
    
    % Prepare MAP/mean profiles
    count = 0;
    for i=1:nDepth
        [mV,mI] = max(pop_2D(i,1:nBins));
        count = count+1;
        MaxMeaSig(count,1) = tmpBins(mI);
        MaxMeaSig(count,2) = sum(tmpBins(:) .* pop_2D(i,1:nBins)') / models;
        sigma_tmp = std(tmpBins(:),pop_2D(i,1:nBins));
        MaxMeaSig(count,3) = MaxMeaSig(count,2)-sigma_tmp;
        MaxMeaSig(count,4) = MaxMeaSig(count,2)+sigma_tmp;
        
        dBins2(count) = dBins(i);
        count = count+1;
        MaxMeaSig(count,1:4) = MaxMeaSig(count-1,1:4);
        if i<nDepth
            dBins2(count) = dBins(i+1);
        else
            dBins2(count) = dBins(end) + (dBins(2)-dBins(1));
        end
    end
    
    % slowness to velocity
    if ptype==1 || ptype==2 % Vs or Vp
        MaxMeaSig(:,1:4) = 1./MaxMeaSig(:,1:4); % slowness to velocity
    end
    
    % pop to PDF
    pdf_2D = pop_2D / models;
    
    % Plot the figure with pop statistics %%%%%%
    Px = 18;
    Py = 12;
    fi = figure('Units','centimeters','Position',[2, 2, Px, Py],'Color',[1 1 1]);
    set(fi, 'PaperOrientation', 'portrait');
    set(fi, 'PaperPositionMode', 'manual');
    set(fi, 'PaperUnits', 'centimeters');
    set(fi, 'Render', 'painters');
    set(fi,'PaperSize', [Px Py]);
    set(fi,'PaperPosition', [0 0 Px Py]);
    
    xLimVar = [tmpBinsP1(1),tmpBinsP1(end)];
    yLimVar = [0,min(plotDepthMAX,dBinsP1(end))];
    
    subplot(1,3,1:2)
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0])
    colormap(dusk);
    surf(tmpBinsP1(:), dBinsP1, [pdf_2D,zeros(nDepth,1);zeros(1,nBins+1)]);
    if plotInDPD == 1
        hold on
        if syntPlot == 1
            plot3(synModX,synMod(:,1),ones(1,length(synModX)),'-','Linewidth',0.5,'color',[0.3 0.3 0.3]);
        end
        plot3(maxModX,maxMod(:,1),ones(1,length(maxModX)),'b-','Linewidth',0.5)
        plot3(mapModX,mapMod(:,1),ones(1,length(mapModX)),'m-','Linewidth',0.5)
        hold off
    end
    shading flat;
    view(0,90);
    set(gca,'ydir','reverse');
    set(gca,'xlim',xLimVar)
    set(gca,'ylim',yLimVar)
    h = colorbar('Location','westoutside');
    ylabel(h, 'Probability')
    nmi(ptype) = max(pdf_2D(:));
    if ptype==1 || ptype==2 % Vs or Vp
        clim([0 max(nmi(1:2))])
    else
        delta = max(nmi(1:2)) - max(nmi(3:4));
        clim([0 max(nmi(1:2)) - (delta/2)])
    end
    box on
    set(gca, 'Layer', 'top');
    
    xlabel([textP,' ',textU]);
    ylabel('Depth (m)');
    title(['Posterior marginal PDF of ',textP])
    set(gca,'FontSize',8)
    
    % pop to ML/MAP/AM
    subplot(1,3,3)
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0])
    hold on
    pop_1Dscaled = (pop_1D/max(pop_1D))*((xLimVar(2)-xLimVar(1))/1.5);
    for i=2:length(pop_1D)-1
        rectangle('position',[xLimVar(1) logdBins(i) pop_1Dscaled(i) (logdBins(i+1)-logdBins(i))], 'FaceColor',[0.6 0.6 0.6],'EdgeColor',[0.6 0.6 0.6])
    end
    line(NaN,NaN,'LineWidth',5,'LineStyle','-','Color',[0.6 0.6 0.6]);
    if syntPlot == 1
        plot(synModX,synMod(:,1),'-','Linewidth',0.5,'color',[0.3 0.3 0.3])
    end
    plot(maxModX,maxMod(:,1),'b-','Linewidth',0.5)
    plot(mapModX,mapMod(:,1),'m-','Linewidth',0.5)
    plot(MaxMeaSig(:,1),dBins2,'r-','Linewidth',0.5)
    plot(MaxMeaSig(:,2),dBins2,'g-','Linewidth',0.5)
    plot(MaxMeaSig(:,3),dBins2,'g:','Linewidth',0.5)
    plot(MaxMeaSig(:,4),dBins2,'g:','Linewidth',0.5)
    hold off
    set(gca,'ydir','reverse');
    set(gca,'xlim',xLimVar)
    set(gca,'ylim',yLimVar)
    box on
    grid on;
    if plotLeg == 1
        if syntPlot == 1
            legend({'Interfaces','Target model','ML model','MAP model','MAX of PDF',textAM,'AM ±\sigma'},'Location','southwest')
        else
            legend({'Interfaces','ML model','MAP model','MAX of PDF',textAM,'AM ±\sigma'},'Location','southwest')
        end
    end
    xlabel([textP,' ',textU]);
    ylabel('Depth (m)');
    title('Profiles')
    set(gca,'FontSize',8)
    set(findall(gcf,'type','text'),'fontSize',8)
    drawnow

    % Save to .pdf file
    if (FigureRendering>0)
        saveas(fi,[OUTprefix,'_marginal_',textOUT,'.pdf'])
        exportgraphics(fi,[OUTprefix,'_marginal_',textOUT,'.png'], 'Resolution', 300);
    end
    if (FigureRendering>1)
        close(fi)
    end

    if ptype==1
        MaxMeaSigVs = MaxMeaSig(1:2:end,:);
    elseif ptype==2
        MaxMeaSigVp = MaxMeaSig(1:2:end,:);
    end

    % Save ML/MAX/AM models to .asii files
    if ptype==1
        % ML model
        fid = fopen([OUTprefix,'_','ML_model.ascii'],'w');
        fprintf(fid,'%s\r\n','# MAXIMUM LIKELIHOOD MODEL. Number of layers:');
        fprintf(fid,'%d\r\n',length(tmpMod(:,1)));
        fprintf(fid,'%s\r\n','# Thickness[m],  Vp[m/s],  Vs[m/s]  and density[kg/m3]');
        for i=1:length(tmpMod(:,1))-1
            fprintf(fid,'%7.1f %7.1f %7.1f %7.1f \r\n',tmpMod(i,1),tmpMod(i,3),tmpMod(i,2),tmpMod(i,4));
        end
        fprintf(fid,'%s\r\n','# Last line is the half-space');
        i = length(tmpMod(:,1));
        fprintf(fid,'%7.1f %7.1f %7.1f %7.1f \r\n',tmpMod(i,1),tmpMod(i,3),tmpMod(i,2),tmpMod(i,4));
        fclose(fid);
        
        % MAP model
        fid = fopen([OUTprefix,'_','MAP_model.ascii'],'w');
        fprintf(fid,'%s\r\n','# MAXIMUM A POSTERIORI MODEL. Number of layers:');
        fprintf(fid,'%d\r\n',length(tmpMod2(:,1)));
        fprintf(fid,'%s\r\n','# Thickness[m],  Vp[m/s],  Vs[m/s]  and density[kg/m3]');
        for i=1:length(tmpMod2(:,1))-1
            fprintf(fid,'%7.1f %7.1f %7.1f %7.1f \r\n',tmpMod2(i,1),tmpMod2(i,3),tmpMod2(i,2),tmpMod2(i,4));
        end
        fprintf(fid,'%s\r\n','# Last line is the half-space');
        i = length(tmpMod2(:,1));
        fprintf(fid,'%7.1f %7.1f %7.1f %7.1f \r\n',tmpMod2(i,1),tmpMod2(i,3),tmpMod2(i,2),tmpMod2(i,4));
        fclose(fid);
    end
    if (ptype==1) || (ptype==2)
        % MAX/AM profiles
        fid = fopen([OUTprefix,'_AM_profile_',textOUT,'.ascii'],'w');
        fprintf(fid,'%s\r\n',['# ',textOUT,' profile from the posterior PDF']);
        fprintf(fid,'%s\r\n','# Number of pseudo-layers');
        fprintf(fid,'%d\r\n',nDepth);
        fprintf(fid,'%s\r\n',['# Thickness[m],  MAXPDF(',textOUT,')[m/s],  AM(',textOUT,')[m/s],  AM_low(',textOUT,')[m/s]  and AM_high(',textOUT,')[m/s]']);
        for i=1:nDepth
            if i<nDepth
                thick = dBins(i+1)-dBins(i);
            else
                thick = 0;
            end
            v_tmp = MaxMeaSig((2*i)-1,1:4);
            fprintf(fid,'%7.1f %7.1f %7.1f %7.1f %7.1f \r\n',thick,v_tmp(1),v_tmp(2),v_tmp(4),v_tmp(3));
        end
        fclose(fid);
    end
end

%% Find and save MAP model uncertainty
DepthCenterMAP = cumsum(tmpMod2(:,1)) - (tmpMod2(:,1)/2);
DepthCenterMAP(end) = dBinsP1(end);
DepthCenterMMS = (dBins + ((dBins(2)-dBins(1))/2))';

fid = fopen([OUTprefix,'_','MAP_model_uncertainty.ascii'],'w');
fprintf(fid,'%s\r\n','# MAXIMUM A POSTERIORI MODEL - UNCERTAINTY. Number of layers:');
fprintf(fid,'%d\r\n',length(tmpMod2(:,1)));
fprintf(fid,'%s\r\n','# Thickness[m],  Vp_low[m/s],  Vp_high[m/s], Vs_low[m/s],  Vs_high[m/s]');
for i=1:length(DepthCenterMAP)
    pos = find(DepthCenterMMS>=DepthCenterMAP(i),1,'first');
    if isempty(pos)
        pos = length(DepthCenterMMS);
    end
    correct = MaxMeaSigVs(pos,2)-tmpMod2(i,2);
    thresVs = MaxMeaSigVs(pos,3:4)-correct;
    thresVs(thresVs<0) = 0;
    correct = MaxMeaSigVp(pos,2)-tmpMod2(i,3);
    thresVp = MaxMeaSigVp(pos,3:4)-correct;
    thresVp(thresVp<0) = 0;
    if i==length(DepthCenterMAP)
        fprintf(fid,'%s\r\n','# Last line is the half-space');
    end
    fprintf(fid,'%7.1f %7.1f %7.1f %7.1f %7.1f\r\n',tmpMod2(i,1),thresVp(2),thresVp(1),thresVs(2),thresVs(1));
end
fclose(fid);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    QWL    %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------
% Read quarter-wavelength representation headers

fid = fopen([POPprefix,'_QWL_head.txt'],'r');
tline = fgets(fid);
tline = fgets(fid);
nBins = str2num(tline);
tline = fgets(fid);
nFreq = str2num(tline);
tline = fgets(fid);
models = str2num(tline);
tline = fgets(fid);
Freq(1:nFreq) = str2num(tline);
tline = fgets(fid);
VsVec(1:nBins) = str2num(tline);
tline = fgets(fid);
DepVec(1:nBins) = str2num(tline);
tline = fgets(fid);
ImpVec(1:nBins) = str2num(tline);
tline = fgets(fid);
AmpVec(1:nBins) = str2num(tline);
tline = fgets(fid);
nXBins = str2num(tline);
tline = fgets(fid);
vs30Bins(1:nXBins) = str2num(tline);
fclose(fid);

% QWL bins for surf plots
fBinsW = log10(Freq(2)) - log10(Freq(1));
fBinsP1 = [10.^(log10(Freq)-fBinsW/2), 10^(log10(Freq(end))+0.5*fBinsW)];
qwlBinsP1(1:nBins+1,1) = [VsVec(1:nBins)';VsVec(end)+(VsVec(2)-VsVec(1))];
qwlBinsP1(1:nBins+1,2) = [DepVec(1:nBins)';10^(log10(DepVec(end))+(log10(DepVec(2))-log10(DepVec(1))))];
qwlBinsP1(1:nBins+1,3) = [ImpVec(1:nBins)';ImpVec(end)+(ImpVec(2)-ImpVec(1))];
qwlBinsP1(1:nBins+1,4) = [AmpVec(1:nBins)';10^(log10(AmpVec(end))+(log10(AmpVec(2))-log10(AmpVec(1))))];

% QWL-frequency tics
fTick = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3 4 5 6 7 8 9 10 20 30 40 50];
fTicklabel = {'0.1' '' '0.3' '' '0.5' '' '' '' '' '1' '2' '3' '4' '5' '' '' '' '' '10' '20' '30' '40' '50'};

% QWL-depth tics
dTick = [1 2 3 4 5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000];
dTicklabel = {'1' '' '' '' '' '' '' '' '' '10' '' '' '' '' '' '' '' '' '100' '' '' '' '' '' '' '' '' '1000' '' '' '' '' '' '' '' '' '10000'};

% SH-amplification tics
aTick = [0.5 0.6 0.7 0.8 0.9 1 2 3 4 5 6 7 8 9 10 20 30 40];
aTicklabel = {'0.5' '' '' '' '' '1' '2' '3' '4' '5' '' '' '' '' '10' '20' '30' '40'};

% N of models
%display(['Number of models: ',num2str(models)])


%% -------------------------------------------------------------------
% Read ML model QWL

fid = fopen([MAXprefix,'_QWL.txt'],'r');
tline = fgets(fid);
nFreqM = str2num(tline);
tline = fgets(fid);
maxDataVs30(1:2) = str2num(tline);
tline = fgets(fid);
FreqDataQWL(1:nFreqM) = str2num(tline);
tline = fgets(fid);
maxDataQWL(1:nFreqM,1) = str2num(tline);
tline = fgets(fid);
maxDataQWL(1:nFreqM,2) = str2num(tline);
tline = fgets(fid);
maxDataQWL(1:nFreqM,3) = str2num(tline);
tline = fgets(fid);
maxDataQWL(1:nFreqM,4) = str2num(tline);
tline = fgets(fid);
maxDataSH(1:nFreqM,1) = str2num(tline);
tline = fgets(fid);
maxDataSH(1:nFreqM,2) = str2num(tline);
fclose(fid);

% Read MAP model QWL
fid = fopen([MAPprefix,'_QWL.txt'],'r');
tline = fgets(fid);
nFreqM = str2num(tline);
tline = fgets(fid);
mapDataVs30(1:2) = str2num(tline);
tline = fgets(fid);
FreqDataQWL2(1:nFreqM) = str2num(tline);
tline = fgets(fid);
mapDataQWL(1:nFreqM,1) = str2num(tline);
tline = fgets(fid);
mapDataQWL(1:nFreqM,2) = str2num(tline);
tline = fgets(fid);
mapDataQWL(1:nFreqM,3) = str2num(tline);
tline = fgets(fid);
mapDataQWL(1:nFreqM,4) = str2num(tline);
tline = fgets(fid);
mapDataSH(1:nFreqM,1) = str2num(tline);
tline = fgets(fid);
mapDataSH(1:nFreqM,2) = str2num(tline);
fclose(fid);


%% -------------------------------------------------------------------
%  Read synthetic model and compute QWL
if syntPlot == 1
    fid = fopen(syntfile,'r');
    tline = fgets(fid);
    tline = fgets(fid);
    nL = str2num(tline);
    tline = fgets(fid);
    synModQWL = zeros(nL-1,4);
    for i=1:nL-1
        tline = fgets(fid);
        synModQWL(i,1:4) = str2num(tline);
    end
    tline = fgets(fid);
    tline = fgets(fid);
    synModQWL(nL,1:4) = str2num(tline);
    fclose(fid);
    
    % Thick to depth
    synModQWL(:,1) = [0;cumsum(synModQWL(1:end-1,1))];
    
    %  Compute QWL
    [z_f,Vs_f,Rho_f,A,z2_f,Vs2_f,IC_f] = respQWL(synModQWL(:,1)',synModQWL(:,3)',synModQWL(:,4)',FreqDataQWL);
    
    %  Compute VS30
    [synf30,synVs30] = qwl30(synModQWL(:,1)',synModQWL(:,3)',30);
    
    %  Compute SH
    Fox = ones(1,nL)*0;
    [h1,h2,h3] = respSH(synModQWL(:,1)',synModQWL(:,3)',synModQWL(:,4)',Fox,FreqDataQWL,0.0);
    
    %  Compute referenced ampl. to the Swiss profile
    %  Read ref. model
    fid = fopen(REFfile,'r');
    refMod = zeros(2,1000);
    for f = 1:1000
        tline = fgets(fid);
        refMod(1:2,f) = str2num(tline);
    end
    fclose(fid);
    refMod(1,:) = [0, cumsum(refMod(1,1:end-1))];
    %  Compute reference QWL
    [refz_f,refVs_f,refRho_f,refA,refz2_f,refVs2_f,refIC_f] = respQWL(refMod(1,:),refMod(2,:),ones(1,1000).*2500,FreqDataQWL);
    %  Correction factor
    ampcor = zeros(1,length(FreqDataQWL));
    for f = 1:length(FreqDataQWL)
        ampcor(f) = sqrt( (synModQWL(end,3)*2500)/(refVs_f(f)*refRho_f(f)) );
    end
end


%% -------------------------------------------------------------------
% Pop QWL statistics

Px = 21;
Py = 14;
fi = figure('Units','centimeters','Position',[2, 2, Px, Py],'Color',[1 1 1]);
set(fi, 'PaperOrientation', 'portrait');
set(fi, 'PaperPositionMode', 'manual');
set(fi, 'PaperUnits', 'centimeters');
set(fi, 'Render', 'painters');
set(fi,'PaperSize', [Px Py]);
set(fi,'PaperPosition', [0 0 Px Py]);

for ptype=1:3
    
    pop_qwl2D = zeros(nFreq,nBins);
    
    % Read pop statistics
    if ptype==1
        fid = fopen([POPprefix,'_QWL_dep2D.txt'],'r');
        tmpBins = 10.^(log10(DepVec) + (log10(DepVec(2))-log10(DepVec(1)))/2);
        tmpBinsP1 = qwlBinsP1(:,2);
        textP = '{\itd}^{QWL}';
        textU = '(m)';
        maxModY = maxDataQWL(:,1);
        mapModY = mapDataQWL(:,1);
        if syntPlot == 1
            synModY = z_f;
        end
    elseif ptype==2
        fid = fopen([POPprefix,'_QWL_vs2D.txt'],'r');
        tmpBins = VsVec + (VsVec(2)-VsVec(1))/2;
        tmpBinsP1 = qwlBinsP1(:,1);
        textP = '{\itV}_S^{QWL}';
        textU = '(m/s)';
        maxModY = maxDataQWL(:,2);
        mapModY = mapDataQWL(:,2);
        if syntPlot == 1
            synModY = Vs_f;
        end
    else
        fid = fopen([POPprefix,'_QWL_imp2D.txt'],'r');
        tmpBins = ImpVec + (ImpVec(2)-ImpVec(1))/2;
        tmpBinsP1 = qwlBinsP1(:,3);
        textP = 'Imp.';
        textU = '';
        maxModY = maxDataQWL(:,4);
        mapModY = mapDataQWL(:,4);
        if syntPlot == 1
            synModY = 1./IC_f;
        end
    end
    for i=1:nFreq
        tline = fgets(fid);
        pop_qwl2D(i,1:nBins) = str2num(tline);
    end
    fclose(fid);
    
%     MaxMeaSig = zeros(nFreq,4);
%     
%     % velocity to slowness
%     if ptype==1 % Logatithmic depth
%         tmpBins = log10(tmpBins(:));
%     elseif ptype==2 % Vs velocity to slowness
%         tmpBins = 1./tmpBins(:);
%     end
%     
%     % Prepare MAP/mean profiles
%     for i=1:nFreq
%         [mV,mI] = max(pop_qwl2D(i,1:nBins));
%         MaxMeaSig(i,1) = tmpBins(mI);
%         MaxMeaSig(i,2) = sum(tmpBins(:) .* pop_qwl2D(i,1:nBins)') / models;
%         sigma_tmp = std(tmpBins(:),pop_qwl2D(i,1:nBins));
%         MaxMeaSig(i,3) = MaxMeaSig(i,2)-sigma_tmp;
%         MaxMeaSig(i,4) = MaxMeaSig(i,2)+sigma_tmp;
%     end
%     
%     % slowness to velocity
%     if ptype==1 % Logatithmic depth
%         MaxMeaSig(:,1:4) = 10.^(MaxMeaSig(:,1:4)); % power of log10
%     elseif ptype==2 % slowness to Vs velocity
%         MaxMeaSig(:,1:4) = 1./MaxMeaSig(:,1:4); % slowness to velocity
%     end
    
    % pop to PDF
    pdf_qwl2D = pop_qwl2D / models;
    
    % pop to PDF
    subplot(2,3,ptype)
    colormap(dusk);
    surf(tmpBinsP1(:),fBinsP1, [pdf_qwl2D,zeros(nFreq,1);zeros(1,nBins+1)]);
    hold on
    if plotInDPD == 1
        cihg=0;
        if syntPlot == 1
            cihg=1;
            hnlg(0+cihg)=plot3(synModY,FreqDataQWL,ones(1,length(synModY)),'-','Linewidth',0.5,'color',[0.3 0.3 0.3]);
        end
        hnlg(1+cihg)=plot3(maxModY,FreqDataQWL,ones(1,length(maxModY)),'b-','Linewidth',0.5);
        hnlg(2+cihg)=plot3(mapModY,FreqDataQWL2,ones(1,length(mapModY)),'m-','Linewidth',0.5);
    end
    hold off
    shading flat;
    view(90,90);
    set(gca,'xdir','reverse');
    set(gca,'yscale','log');
    if ptype==1
        set(gca,'xscale','log');
        set(gca,'xtick',dTick);
        set(gca,'xticklabel',dTicklabel);
    end
    set(gca,'ytick',fTick);
    set(gca,'yticklabel',fTicklabel);
    set(gca,'YTickLabelRotation',0);
    xLimVar = [tmpBinsP1(1),tmpBinsP1(end)];
    yLimVar = [fBinsP1(1),fBinsP1(end)];
    set(gca,'xlim',xLimVar)
    set(gca,'ylim',yLimVar)
    box on
    clim([0 maxQWLprop])
    set(gca, 'Layer', 'top');
    if ptype==1
        xlabel(['QWL depth ',textU]);
    elseif ptype==2
        xlabel(['QWL velocity ',textU]);
    elseif ptype==3
        xlabel('QWL impedance');
    else
        xlabel([textP,' ',textU]);
    end
    ylabel('Frequency (Hz)');
    title('Posterior marginal PDF')
    set(gca,'FontSize',8)
    
end

% ---------------------------------
% Read vs30 statistics
fid = fopen([POPprefix,'_vs30.txt'],'r');
tline = fgets(fid);
pop_vs30 = str2num(tline);
fclose(fid);

% Prepare centered histogram values
vsBinsWS30 = (vs30Bins(2)-vs30Bins(1));
vsBinsCS30 = vs30Bins+(vsBinsWS30/2);

% Prepare vs30 MAP/mean
vsBinsTMP = 1./vsBinsCS30;
MaxMeaSig = zeros(1,4);
[mV,mI] = max(pop_vs30(:));
MaxMeaSig(1) = vsBinsTMP(mI);
MaxMeaSig(2) = sum(vsBinsTMP(:) .* pop_vs30(:)) / models;
sigma_tmp = std(vsBinsTMP(:),pop_vs30(:));
MaxMeaSig(3) = MaxMeaSig(2)-sigma_tmp;
MaxMeaSig(4) = MaxMeaSig(2)+sigma_tmp;
MaxMeaSig(5) = MaxMeaSig(2)-2*sigma_tmp;
MaxMeaSig(6) = MaxMeaSig(2)+2*sigma_tmp;
MaxMeaSig = 1./MaxMeaSig;

% Plot
pop_vs30MAX = max(pop_vs30(:)./models);
nXLim = [0, 1.2*pop_vs30MAX];

subplot(2,3,4:5)
hold on
i=0;
if syntPlot == 1
    hch(2) = line([synVs30 synVs30],[0 max(nXLim)/10],'Linewidth',0.5,'color',[0.3 0.3 0.3]); % Target model
    i=1;
end
hch(2+i) = line([maxDataVs30(2) maxDataVs30(2)],[0 max(nXLim)/10],'Color','b','Linewidth',0.5); % ML model
hch(3+i) = line([mapDataVs30(2) mapDataVs30(2)],[0 max(nXLim)/10],'Color','m','Linewidth',0.5); % MAP model
hch(1) = bar(vsBinsCS30, pop_vs30./models, 1.0, 'FaceColor',[0.5 0.5 0.5], 'EdgeColor',[0.5 0.5 0.5]); % PDF
if syntPlot == 1
    hch(4+i) = plot(synVs30,0, 'x','MarkerSize',5,'MarkerEdgeColor',[0.3 0.3 0.3],'MarkerFaceColor',[0.3 0.3 0.3]); % Target
    i=2;
end
hch(4+i) = plot(maxDataVs30(2),0, 'x','MarkerSize',5,'MarkerEdgeColor','b','MarkerFaceColor','b');% ML 
hch(5+i) = plot(mapDataVs30(2),0, 'x','MarkerSize',5,'MarkerEdgeColor','m','MarkerFaceColor','m');% MAP 
hch(6+i) = plot(MaxMeaSig(1),pop_vs30MAX, 'x','MarkerSize',5,'MarkerEdgeColor','r','MarkerFaceColor','r'); % MAX
hch(8+i) = plot([MaxMeaSig(5) MaxMeaSig(6)],[max(nXLim)/10 max(nXLim)/10],'LineStyle',':','Color','g','Marker','none'); % AM (2 sigma)
hch(7+i) = errorbar(MaxMeaSig(2),max(nXLim)/10,MaxMeaSig(2)-MaxMeaSig(4),MaxMeaSig(3)-MaxMeaSig(2),'o','horizontal','Color','g','MarkerSize',5,'CapSize',5); % AM (1 sigma)
set(gca,'YLim',nXLim)
gylim = get(gca,'XLim');
hold off
set(gca,'XLim',[floor(MaxMeaSig(2)-5*(MaxMeaSig(2)-MaxMeaSig(4))) ceil(MaxMeaSig(2)+5*(MaxMeaSig(3)-MaxMeaSig(2)))])
box on
grid on
ylabel('Probability');
xlabel('{\itV}_{S30} (m/s)');           
set(gca,'Xcolor','k')
if syntPlot == 1
    legend(hch,{'Distribution','Target model','ML model','MAP model','{\itV}_{S30}^{ Target}','{\itV}_{S30}^{ ML}',...
        '{\itV}_{S30}^{ MAP}','{\itV}_{S30}^{ MAX}','{\itV}_{S30}^{ AM} (1\sigma)','{\itV}_{S30}^{ AM} (2\sigma)'},'Location','northeastoutside')
else
    legend(hch,{'Distribution','ML model','MAP model','{\itV}_{S30}^{ ML}','{\itV}_{S30}^{ MAP}',...
        '{\itV}_{S30}^{ MAX}','{\itV}_{S30}^{ AM} (1\sigma)','{\itV}_{S30}^{ AM} (2\sigma)'},'Location','northeastoutside')
end
title('Average {\itV}_{S} down to 30 m')
set(gca,'FontSize',8)

% ---------------------------------
% Plot Vs30 text results
subplot(2,3,6)

if syntPlot == 1
    vs30txt =['{\itV}_{S30}^{Target}=',num2str(round(synVs30)),' (m/s)'];
    f30txt =['f_{30}^{Target}=',num2str(0.01*round(100*synVs30/(4*30))),' (Hz)'];
else
    vs30txt ='';
    f30txt ='';
end

text(0,1,{vs30txt,...
    ['{\itV}_{S30}^{ ML}=',num2str(round(maxDataVs30(2))),' m/s'],...
    ['{\itV}_{S30}^{ MAP}=',num2str(round(mapDataVs30(2))),' m/s'],...
    ['{\itV}_{S30}^{ MAX}=',num2str(round(MaxMeaSig(1))),' m/s'],...
    ['{\itV}_{S30}^{ AM}=',num2str(round(MaxMeaSig(2))),' m/s'],...
    ['(',num2str(round(MaxMeaSig(4))),' - ',num2str(round(MaxMeaSig(3))),' m/s)']},...
    'VerticalAlignment','top','HorizontalAlignment','Left','FontSize',8)

text(0.5,1,{f30txt,...
    ['f_{30}^{ ML}=',num2str(0.01*round(100*maxDataVs30(2)/(4*30))),' Hz'],...
    ['f_{30}^{ MAP}=',num2str(0.01*round(100*mapDataVs30(2)/(4*30))),' Hz'],...
    ['f_{30}^{ MAX}=',num2str(0.01*round(100*MaxMeaSig(1)/(4*30))),' Hz'],...
    ['f_{30}^{ AM}=',num2str(0.01*round(100*MaxMeaSig(2)/(4*30))),' Hz'],...
    ['(',num2str(0.1*round(10*MaxMeaSig(4)/(4*30))),' - ',num2str(0.1*round(10*MaxMeaSig(3)/(4*30))),' Hz)']},...
    'VerticalAlignment','top','HorizontalAlignment','Left','FontSize',8)

hold off
set(gca,'XLim',[0 1])
set(gca,'YLim',[0 1.2])
axis off

colormap(dusk);
clim([0 maxQWLprop])
h = colorbar('Location','northoutside');
ylabel(h, 'Probability')
set(h,'FontSize',8)

% ---------------------------------
% Save to .pdf file
if (FigureRendering>0)
    saveas(fi,[OUTprefix,'_QWL.pdf'])
    exportgraphics(fi,[OUTprefix,'_QWL.png'], 'Resolution', 300);
end
if (FigureRendering>1)
    close(fi)
end

% ---------------------------------
% Save Vs30 into text file
fid = fopen([OUTprefix,'_','Vs30.ascii'],'w');
fprintf(fid,'%s\r\n','# Vs30[m/s],  f30[Hz],  Model Type');
fprintf(fid,'%8.2f %8.4f  %s\r\n',maxDataVs30(2),maxDataVs30(2)/(4*30),'ML model');
fprintf(fid,'%8.2f %8.4f  %s\r\n',mapDataVs30(2),mapDataVs30(2)/(4*30),'MAP model');
fprintf(fid,'%8.2f %8.4f  %s\r\n',MaxMeaSig(1),MaxMeaSig(1)/(4*30),'MAX from PDF');
fprintf(fid,'%8.2f %8.4f  %s\r\n',MaxMeaSig(2),MaxMeaSig(2)/(4*30),'AM from PDF');
fprintf(fid,'%8.2f %8.4f  %s\r\n',MaxMeaSig(4),MaxMeaSig(4)/(4*30),'Low threshold from PDF (-1 sigma)');
fprintf(fid,'%8.2f %8.4f  %s\r\n',MaxMeaSig(3),MaxMeaSig(3)/(4*30),'High threshold from PDF (+1 sigma)');
fprintf(fid,'%8.2f %8.4f  %s\r\n',MaxMeaSig(6),MaxMeaSig(6)/(4*30),'Low threshold from PDF (-2 sigma)');
fprintf(fid,'%8.2f %8.4f  %s\r\n',MaxMeaSig(5),MaxMeaSig(5)/(4*30),'High threshold from PDF (+2 sigma)');
fclose(fid);


%% -------------------------------------------------------------------
% Pop SH amplification of the whole profile

Px = 20;
Py = 14;
fi = figure('Units','centimeters','Position',[2, 2, Px, Py],'Color',[1 1 1]);
set(fi, 'PaperOrientation', 'portrait');
set(fi, 'PaperPositionMode', 'manual');
set(fi, 'PaperUnits', 'centimeters');
set(fi, 'Render', 'painters');
set(fi,'PaperSize', [Px Py]);
set(fi,'PaperPosition', [0 0 Px Py]);

%--------------------------------------------------------------------------
% Rozvrhnout subplots

sx2(1) = subplot(2,2,3);
set(sx2(1), 'OuterPosition', [0.02, 0.02, 0.45, 0.38]);
sx2(2) = subplot(2,2,4);
set(sx2(2), 'OuterPosition', [0.50, 0.02, 0.45, 0.38]);

sx1(1) = subplot(2,2,1);
set(sx1(1), 'OuterPosition', [0.02, 0.35, 0.45, 0.6]);
sx1(2) = subplot(2,2,2);
set(sx1(2), 'OuterPosition', [0.50, 0.35, 0.45, 0.6]);

for ptype=1:2
    
    pop_qwl2D = zeros(nFreq,nBins);
    MaxMeaSig = zeros(nFreq,4);
    
    % Read pop statistics
    if ptype==1
        fid = fopen([POPprefix,'_SH_unr2D.txt'],'r');
        tmpBins = 10.^(log10(AmpVec) + (log10(AmpVec(2))-log10(AmpVec(1)))/2);
        tmpBinsP1 = qwlBinsP1(:,4);
        textP = 'Amplification';
        textU = '';
        maxModY = maxDataSH(:,1);
        mapModY = mapDataSH(:,1);
        if syntPlot == 1
            synModY = abs(h3);
        end
    else
        fid = fopen([POPprefix,'_SH_ref2D.txt'],'r');
        tmpBins = 10.^(log10(AmpVec) + (log10(AmpVec(2))-log10(AmpVec(1)))/2);
        tmpBinsP1 = qwlBinsP1(:,4);
        textP = 'Amplification';
        textU = '';
        maxModY = maxDataSH(:,2);
        mapModY = mapDataSH(:,2);
        if syntPlot == 1
            synModY = abs(h3)./ampcor;
        end
    end
    for i=1:nFreq
        tline = fgets(fid);
        pop_qwl2D(i,1:nBins) = str2num(tline);
    end
    fclose(fid);
    
    % Logatithmic depth
    tmpBins = log10(tmpBins(:));
    
    % Prepare MAP/mean profiles
    for i=1:nFreq
        [mV,mI] = max(pop_qwl2D(i,1:nBins));
        MaxMeaSig(i,1) = tmpBins(mI);
        MaxMeaSig(i,2) = sum(tmpBins(:) .* pop_qwl2D(i,1:nBins)') / models;
        sigma_tmp = std(tmpBins(:),pop_qwl2D(i,1:nBins));
        MaxMeaSig(i,3) = MaxMeaSig(i,2)-sigma_tmp;
        MaxMeaSig(i,4) = MaxMeaSig(i,2)+sigma_tmp;
    end
    
    % Logatithmic depth
    MaxMeaSig(:,1:4) = 10.^(MaxMeaSig(:,1:4)); % power of log10
    
    % pop to PDF
    pdf_qwl2D = pop_qwl2D / models;
    
    % plot PDF
    axes(sx1(ptype));
    set(gca,'XColor',[0 0 0],'YColor',[0 0 0])
    colormap(dusk);
    surf(tmpBinsP1(:),fBinsP1, [pdf_qwl2D,zeros(nFreq,1);zeros(1,nBins+1)]);
    shading flat;
    view(90,90);
    set(gca,'xdir','reverse');
    set(gca,'yscale','log');
    set(gca,'xscale','log');
    set(gca,'xtick',aTick);
    set(gca,'xticklabel',aTicklabel);
    set(gca,'ytick',fTick);
    set(gca,'yticklabel',fTicklabel);
    set(gca,'yTickLabelRotation',0);
    xLimVar = [tmpBinsP1(1),tmpBinsP1(end)];
    yLimVar = [fBinsP1(1),fBinsP1(end)];
    set(gca,'xlim',xLimVar)
    set(gca,'ylim',yLimVar)
    box on
    set(gca, 'Layer', 'top');
    h = colorbar('Location','southoutside');
    ylabel(h, 'Probability')
    xlabel([textP,' ',textU]);
    ylabel('Frequency (Hz)');
    if ptype==1
        title({'Posterior marginal PDF of','SH-wave transfer function'})
    else
        title({'Posterior marginal PDF of','Amplification to the ref. profile'})
    end
    set(gca,'FontSize',8)
    
    % plot ML and MAP
    axes(sx2(ptype));
    hold on
    if syntPlot == 1
        plot(FreqDataQWL,synModY,'-','Linewidth',0.5,'color',[0.3 0.3 0.3])
    end
    plot(FreqDataQWL,maxModY,'b-','Linewidth',0.5)
    plot(FreqDataQWL,mapModY,'m-','Linewidth',0.5)
    plot(Freq(1:nFreq),MaxMeaSig(:,2),'g-','Linewidth',0.5)
    plot(Freq(1:nFreq),MaxMeaSig(:,3),'g:','Linewidth',0.5)
    plot(Freq(1:nFreq),MaxMeaSig(:,4),'g:','Linewidth',0.5)
    hold off
    set(gca,'xscale','log');
    set(gca,'xtick',fTick);
    set(gca,'xticklabel',fTicklabel);
    set(gca,'xTickLabelRotation',0);
    set(gca,'yscale','log');
    set(gca,'ytick',aTick);
    set(gca,'yticklabel',aTicklabel);
    set(gca,'ylim',xLimVar)
    set(gca,'xlim',yLimVar)
    box on
    grid on;
    xlabel('Frequency (Hz)');
    ylabel([textP,' ',textU]);
    %title('ML/MAP/AM')
    set(gca,'FontSize',8)
    
    if ptype==1
        if plotLeg == 1
            if syntPlot == 1
                legend({'Target model','ML model','MAP model','AM of PDF','AM ±\sigma'},'Location','northwest')
            else
                legend({'ML model','MAP model','AM of PDF','AM ±\sigma'},'Location','northwest')
            end
        end
    end

    % Save amplification to .asii files
    if ptype==2
        % ML model amplification
        fid = fopen([OUTprefix,'_','ML_amplification.ascii'],'w');
        fprintf(fid,'%s\r\n','# MAXIMUM LIKELIHOOD MODEL: Amplification to the reference velocity profile');
        fprintf(fid,'%s\r\n','# Frequency[Hz],  Amplification[-]');
        for i=1:length(FreqDataQWL(:))
            fprintf(fid,'%9.5f %9.5f\r\n',FreqDataQWL(i),maxModY(i));
        end
        fclose(fid);
        % MAP model amplification
        fid = fopen([OUTprefix,'_','MAP_amplification.ascii'],'w');
        fprintf(fid,'%s\r\n','# MAXIMUM A POSTERIORI MODEL: Amplification to the reference velocity profile');
        fprintf(fid,'%s\r\n','# Frequency[Hz],  Amplification[-]');
        for i=1:length(FreqDataQWL(:))
            fprintf(fid,'%9.5f %9.5f\r\n',FreqDataQWL(i),mapModY(i));
        end
        fclose(fid);
        % AM profiles
        fid = fopen([OUTprefix,'_','AM_profile_amplification.ascii'],'w');
        fprintf(fid,'%s\r\n','# Arithmetic Mean profile: Amplification to the reference velocity profile (elastic)');
        fprintf(fid,'%s\r\n','# Frequency[Hz],  Amplification[-],  Ampl_low[-],  Ampl_high[-]');
        for i=1:nFreq
            fprintf(fid,'%9.5f %9.5f %9.5f %9.5f\r\n',Freq(i),MaxMeaSig(i,2),MaxMeaSig(i,3),MaxMeaSig(i,4));
        end
        fclose(fid);
    end

end


% ---------------------------------
% Save to .pdf file
if (FigureRendering>0)
    saveas(fi,[OUTprefix,'_SH.pdf'])
    exportgraphics(fi,[OUTprefix,'_SH.png'], 'Resolution', 300);
end
if (FigureRendering>1)
    close(fi)
end

disp('All done.')

