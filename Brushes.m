clc;
clear;
close all;
warning('off', 'MATLAB:handle_graphics:exceptions:SceneNode');
tic

dualframeRate = 10; % Hz or fps
movieLength = "5 Minutes";
dualSecPerFrame = 1/dualframeRate;
secPerFrame = dualSecPerFrame;

doExcludeTraces = 1;
includeFastFlops = 1;
plotDwellScatter = 0;
doSliderPlot = 1;
plotIndividualCumSum = 0;
plotTotalTimeHistograms = 0;
cutoffFraction = 0.95;
doPlotHists = 0;
appendFastFlops = 1;
plotBrushesSamples = 0;
plotNonMarkov = 0;
plotSimpleRates = 0;
numBins = 20;
barPlots = 0;
brushBarPlots = 0;
printHypothesisTest = 1;

doPlotTraces = 1; % Plot a Selection of Traces?
plotTraces = {[43],...
    [1],[5],[33]};

% plotTraces = {[41, 20],...
%     [1, 4],[5, 8],[12, 16]};

%% Loading Data
dualData = {'0_0' '0_10' '30_0' '30_10';...
    [], [], [],[]; % Left (Cy5) Path [2]
    [], [], [],[]; % Left (Cy5) Raw "FRET" [3]
    [], [], [],[]; % Right (Cy3) Path [4]
    [], [], [],[]; % Right (Cy3) Raw "FRET" [5]
    [], [], [],[]}; % Three State [6]
numSamples = size(dualData, 2);
dualNames = {'0\_0', '0\_10', '30\_0', '30\_10'};
for i=1:length(dualNames)
    simpleDualDwells(i).name = dualNames(i);    
    threeStateDwells(i).name = dualNames(i);
end

outsideBoxSpecs = [0.68 0.15 0.3 0.16];
insideBoxSpecs = [0.601666666666666 0.13 0.182500000000002 0.0875];
foursampleOutsideBox = [0.814166666666669 0.14 0.181666666666676 0.21875];
foursampleLegendSpecs = [0.821483335655433 0.466562506705523 0.162499995355805 0.459374986588955];
twosampleLegendSpecs = [0.614816668988766 0.239687503967434 0.162499995355805 0.278124992065133];

old0_0exclude = [7, 21, 24, 49, 50, 53, 56];

dualExcludeTraces = {[50, 73, 178, 7, 56, 86, 89, 92, 93, 94, 109, 135, 184, 97, ...
    21, 24, 62, 80, 116, 120, 122, 127, 149, 152, 49, 53, 74, 177, 181], ...
    [20, 29, 39, 50, 56, 57, 33, 10, 35, 40, 41, 49, 16, 26], ...
    [8, 12, 65, 21, 22, 28, 17, 68, 9, 10],...
    [8, 16, 22, 31, 48, 57, 77, 85, 94, 106, 17, 27, 101, 97]};

olddualExcludeTraces = {[11, 14, 15, 18, 26, 37], ...
    [20, 29, 39, 50, 56, 57, 33, 10, 35, 40, 41, 49, 16, 26], ...
    [8, 12, 65, 21, 22, 28, 17, 68, 9, 10],...
    [8, 16, 22, 31, 48, 57, 77, 85, 94, 106, 17, 27, 101, 97]};

% dualExcludeTraces = olddualExcludeTraces;

excludeAllButFastFlops  = {[50, 73, 178, 7, 56, 86, 89, 92, 93, 94, 109, 135, 184, 97, ...
    49, 53, 74, 177, 181, 67, 120, 122, 127], ...
    [20, 29, 39, 50, 56, 57, 10, 35, 40, 41, 49, 16, 26], ...
    [8, 12, 65, 17, 68, 9, 10, 28],...
    [8, 16, 22, 31, 48, 57, 77, 85, 94, 106, 97, 95, 101]};

if includeFastFlops
    dualExcludeTraces = excludeAllButFastFlops;
end

for i=1:size(dualData, 2)
    % Loading in paths and excluding traces
    cy5File = strcat("croppedTraces/vbOutput_cy5_", dualData{1,i}, ".mat");
    HMMleft = load(cy5File, "path");
    HMMleft = HMMleft.path;
    Rawleft = load(cy5File, "FRET");
    Rawleft = Rawleft.FRET;
    cy3File = strcat("croppedTraces/vbOutput_cy3_", dualData{1,i}, ".mat");
    HMMright = load(cy3File, "path");
    HMMright = HMMright.path;
    Rawright = load(cy3File, "FRET");
    Rawright = Rawright.FRET;
    
    if doExcludeTraces
        HMMleft(dualExcludeTraces{i}) = [];
        HMMright(dualExcludeTraces{i}) = [];
        Rawleft(dualExcludeTraces{i}) = [];
        Rawright(dualExcludeTraces{i}) = [];
    end
    
    % Compute Threestate traces
    % i is sample number, k is number of traces in i'th sample, j is frame
    % of k'th trace
    dualQuenchHMM = cell(1,size(HMMleft,2));
    for k = 1:size(HMMleft,2)
        lq = normalize(HMMleft{k}', "range");
        rq = normalize(HMMright{k}', "range");
        threestate = [];
        for j = 1:size(lq, 2)
            if lq(j)==0
                if rq(j)==0
                    threestate(j) = 0;
                else
                    threestate(j) = 1;
                end
            else
                if rq(j)==0
                    threestate(j) = 3;
                else
                    threestate(j) = 2;
                end
            end
        end
        dualQuenchHMM{k} = threestate;
    end
    dualData{2, i} = HMMleft;
    dualData{3, i} = Rawleft;
    dualData{4, i} = HMMright;
    dualData{5, i} = Rawright;
    dualData{6, i} = dualQuenchHMM;
end


% Printing Num Traces
for i = 1:size(dualData,2)
    disp(dualData{1,i});
    fprintf('kept: %d \n', size(dualData{2,i},2));
    fprintf('rejected: %d \n', size(dualExcludeTraces{i},2));
    fprintf('total: %d \n', size(dualData{2,i},2)+size(dualExcludeTraces{i},2));
end
%% Plotting Traces
if doSliderPlot
    for i=1:size(dualData, 2)
        % slider_plot(dualData{3, i}, dualData{2, i}, dualData{5, i}, dualData{4, i}, dualData{6, i}, dualSecPerFrame, dualData{1,i})
        fret_slider_plot(dualData{3, i}, dualData{2, i}, dualData{5, i}, dualData{4, i}, dualData{6, i}, dualSecPerFrame, dualData{1,i})
    end
end

if doPlotTraces
    for sample = 1:numSamples
        for trace = plotTraces{:,sample}
            leftQuenchHMMSubset = dualData{2, sample};
            leftQuenchRawSubset = dualData{3, sample};
            rightQuenchHMMSubset = dualData{4, sample};
            rightQuenchRawSubset = dualData{5, sample};
            dualQuenchHMM = dualData{6, sample};
            numBins = 30;

            % Prepare time axis
            seconds = 0:secPerFrame:(size(leftQuenchHMMSubset{trace}, 1)-1)*secPerFrame;

            % Prepare FRET traces
            % Old method of reading out "FRET" as Left/(Right + Left) 

            % fretRaw = leftQuenchRawSubset{trace} ./ (rightQuenchRawSubset{trace} + leftQuenchRawSubset{trace});
            % fretHMM = leftQuenchHMMSubset{trace} ./ (rightQuenchHMMSubset{trace} + leftQuenchHMMSubset{trace});
            % midHMM = max(leftQuenchHMMSubset{trace}) / (max(rightQuenchHMMSubset{trace}) + max(leftQuenchHMMSubset{trace}));

            % % Alternate "FRET" calculation
            % leftNormRaw = (leftQuenchRawSubset{trace}-min(leftQuenchHMMSubset{trace}))./(max(leftQuenchHMMSubset{trace}) - min(leftQuenchHMMSubset{trace}));
            % leftNormHMM = (leftQuenchHMMSubset{trace}-min(leftQuenchHMMSubset{trace}))./(max(leftQuenchHMMSubset{trace}) - min(leftQuenchHMMSubset{trace}));
            % rightNormRaw = (rightQuenchRawSubset{trace}-min(rightQuenchHMMSubset{trace}))./(max(rightQuenchHMMSubset{trace}) - min(rightQuenchHMMSubset{trace}));
            % rightNormHMM = (rightQuenchHMMSubset{trace}-min(rightQuenchHMMSubset{trace}))./(max(rightQuenchHMMSubset{trace}) - min(rightQuenchHMMSubset{trace}));
            
            % fretRaw = (leftNormRaw - rightNormRaw) ./ (leftNormRaw + rightNormRaw);
            % fretHMM = (leftNormHMM - rightNormHMM) ./ (leftNormHMM + rightNormHMM);
            % midHMM = (max(leftNormHMM) - max(rightNormHMM)) / (max(rightNormHMM) + max(leftNormHMM));

            % Alternate Alternate "FRET" calculation
            % Normalize each trace so that the AVERAGE high and low values are 1 and 0 respectively
            % then take the difference (left - right)
            leftNormRaw = (leftQuenchRawSubset{trace}-min(leftQuenchHMMSubset{trace}))./(max(leftQuenchHMMSubset{trace}) - min(leftQuenchHMMSubset{trace}));
            leftNormHMM = (leftQuenchHMMSubset{trace}-min(leftQuenchHMMSubset{trace}))./(max(leftQuenchHMMSubset{trace}) - min(leftQuenchHMMSubset{trace}));
            rightNormRaw = (rightQuenchRawSubset{trace}-min(rightQuenchHMMSubset{trace}))./(max(rightQuenchHMMSubset{trace}) - min(rightQuenchHMMSubset{trace}));
            rightNormHMM = (rightQuenchHMMSubset{trace}-min(rightQuenchHMMSubset{trace}))./(max(rightQuenchHMMSubset{trace}) - min(rightQuenchHMMSubset{trace}));
            
            fretRaw = (leftNormRaw - rightNormRaw);
            fretHMM = (leftNormHMM - rightNormHMM); 
            midHMM = (max(leftNormHMM) - max(rightNormHMM));


            % Create figure
            plFig = figure('Position', [10 10 1310 900]);

            % Axes positions (in pixels)
            ax2 = axes('Units','pixels','Position',[180 85+2*265 1000 200]);   % Right channel (top)
            ax5 = axes('Units','pixels','Position',[180+1010 85+2*265 80 200]);   % Right hist
            ax3 = axes('Units','pixels','Position',[180 85+265 1000 200]);       % Left channel (middle)
            ax6 = axes('Units','pixels','Position',[180+1010 85+265 80 200]);       % Left hist
            ax1 = axes('Units','pixels','Position',[180 85 1000 200]); % FRET (bottom)
            ax4 = axes('Units','pixels','Position',[180+1010 85 80 200]); % FRET hist

            % Right channel plot (top)
            axes(ax2);
            plot(seconds, rightQuenchHMMSubset{trace}', '-','LineWidth', 2,'Color',[0 0 0.2], 'DisplayName',"Right Quenching HMM")
            hold on
            plot(seconds, rightQuenchRawSubset{trace}', '-', 'Color',[0 0 0.6], 'DisplayName',"Right Quenching Raw")
            % xlabel('Time (s)', 'fontweight','bold','fontsize',14)
            ylabel('Right (Cy3)', 'fontweight','bold','fontsize',14)
            set(gca, 'linewidth',5, 'fontweight','bold', 'fontsize',20)
            % legend

            % Right histogram
            axes(ax5);
            histogram(rightQuenchRawSubset{trace}, numBins, 'Normalization','probability', 'FaceColor',[0 0 0.6])
            axis off
            view([90 -90])

            % Left channel plot (middle)
            axes(ax3);
            plot(seconds, leftQuenchHMMSubset{trace}', '-','LineWidth', 2, 'Color',[0.4 0 0], 'DisplayName',"Left Quenching HMM")
            hold on
            plot(seconds, leftQuenchRawSubset{trace}', '-','Color', [0.8 0 0], 'DisplayName',"Left Quenching Raw")
            % xlabel('Time (s)', 'fontweight','bold','fontsize',14)
            ylabel('Left (Cy5)', 'fontweight','bold','fontsize',14)
            set(gca, 'linewidth', 5, 'fontweight','bold', 'fontsize',20)
            % legend

            % Left histogram
            axes(ax6);
            histogram(leftQuenchRawSubset{trace}, numBins, 'Normalization','probability', 'FaceColor',[0.8 0 0])
            axis off
            view([90 -90])

            % FRET plot (bottom)
            axes(ax1);
            plot(seconds, fretRaw, '-','LineWidth', 1,'Color',[0 0 0], 'DisplayName',"Dual Quenching Raw")
            hold on
            plot(seconds, fretHMM, '-','LineWidth', 2,'Color',[0.60 0 0.60], 'DisplayName',"Dual Quenching HMM")
            xlabel('Time (s)', 'fontweight','bold','fontsize',14)
            tickvals = [min(fretHMM) midHMM max(fretHMM)];
            yticks(tickvals)
            ylim([min(fretRaw) max(fretRaw)])
            names = {'Left Bound'; 'Unbound'; 'Right Bound'};
            set(gca,'ytick',tickvals,'yticklabel',names, 'linewidth', 5, 'fontweight','bold', 'fontsize',20)
            legend off

            % FRET histogram
            axes(ax4);
            histogram(fretRaw, numBins, 'Normalization','probability', 'FaceColor',[0 0 0])
            axis off
            view([90 -90]);

            % Add one title at the top for trace number
            sgtitle(strcat(dualNames{sample},sprintf(" Trace %d", trace)), 'fontweight','bold','fontsize',25)

            % Save figure
            if ~exist('./PlotTraces', 'dir')
            mkdir('./PlotTraces');
            end
            saveas(plFig, strcat("./PlotTraces/",dualData{1,sample}, sprintf(" Trace %d.png", trace)))
        end
    end
end

%% Calculating Dwell Times
for i=1:size(dualData, 2)
    
    [~, leftThigh, leftTlow, leftThighInd, leftTlowInd, ~, ~, ~, ~] = eb_dwelltimes_traceInd_MN(dualData{2, i}', dualSecPerFrame);
    simpleDualDwells(i).leftLow = leftTlow;
    simpleDualDwells(i).leftTlowInd = leftTlowInd;
    simpleDualDwells(i).leftHigh = leftThigh;
    simpleDualDwells(i).leftThighInd = leftThighInd;
    simpleDualDwells(i).leftThighInd = leftThighInd;
    
    disp(strcat('Sample ', threeStateDwells(i).name{1}))
    [leftDwells, rightDwells, middleLeftRebindDwells, middleRightRebindDwells, ...
        middleLeftRightTransitionDwells, middleRightLeftTransitionDwells, fastFlopsLR, fastFlopsRL] = eb_dwelltimes_dualquench_bothbound(dualData{2, i}, dualData{4, i}, dualSecPerFrame);
    threeStateDwells(i).leftDwells = leftDwells;
    threeStateDwells(i).rightDwells = rightDwells;
    threeStateDwells(i).middleLeftRebindDwells = middleLeftRebindDwells;
    threeStateDwells(i).middleRightRebindDwells = middleRightRebindDwells;
    if appendFastFlops
        appendableFastFlopsLR = [dualSecPerFrame.*ones(1,size(fastFlopsLR,2));fastFlopsLR];
        appendableFastFlopsRL = [dualSecPerFrame.*ones(1,size(fastFlopsRL,2));fastFlopsRL];
        threeStateDwells(i).middleLeftRightTransitionDwells = horzcat(middleLeftRightTransitionDwells, appendableFastFlopsLR);
        threeStateDwells(i).middleRightLeftTransitionDwells = horzcat(middleRightLeftTransitionDwells, appendableFastFlopsRL);
    else
        threeStateDwells(i).middleLeftRightTransitionDwells = middleLeftRightTransitionDwells;
        threeStateDwells(i).middleRightLeftTransitionDwells = middleRightLeftTransitionDwells;
    end
    %% Collating all Dwell Times out of the middle state independent of path
    threeStateDwells(i).middleLeftDwells = [middleLeftRebindDwells middleRightLeftTransitionDwells];
    threeStateDwells(i).middleRightDwells = [middleRightRebindDwells middleLeftRightTransitionDwells];
    threeStateDwells(i).fastFlopsLR = fastFlopsLR;
    threeStateDwells(i).fastFlopsRL = fastFlopsRL;
    
    [~, rightThigh, rightTlow, rightThighInd, rightTlowInd, ~, ~, ~, ~] = eb_dwelltimes_traceInd_MN(dualData{4, i}', dualSecPerFrame);
    simpleDualDwells(i).rightLow = rightTlow;
    simpleDualDwells(i).rightTlowInd = rightTlowInd;
    simpleDualDwells(i).rightHigh = rightThigh;
    simpleDualDwells(i).rightThighInd = rightThighInd;
    
    numTraces = max(leftTlowInd);
    for j = 1:numTraces
        simpleDualDwells(i).leftMeanTlow(j) = mean(leftTlow(leftTlowInd == j));
        simpleDualDwells(i).leftMeanThigh(j) = mean(leftThigh(leftThighInd == j));
        simpleDualDwells(i).rightMeanTlow(j) = mean(rightTlow(rightTlowInd == j));
        simpleDualDwells(i).rightMeanThigh(j) = mean(rightThigh(rightThighInd == j));
    end
end

%% Per Molecule Average Dwell times
if plotDwellScatter
    % Simple Dual Avg. Dwell Times
    simpleDualAvgDwells = figure('Position', [10 100 1900 900]);
    for i=1:numSamples
        subplot(2,numSamples,i*2-1)
        scatter(simpleDualDwells(i).leftMeanTlow, simpleDualDwells(i).leftMeanThigh, 'LineWidth',2)
        grid on
        xlabel('Mean Low Dwell Time (s)')
        ylabel('Mean High Dwell Time (s)')
        xlim([0 100])
        ylim([0 210])
        title(strcat(simpleDualDwells(i).name, " Left"), 'Interpreter','none')
        set(gca,'fontweight','bold', 'FontSize', 14)
        
        subplot(2,numSamples,i*2)
        scatter(simpleDualDwells(i).rightMeanTlow, simpleDualDwells(i).rightMeanThigh, 'LineWidth',2)
        grid on
        xlabel('Mean Low Dwell Time (s)')
        ylabel('Mean High Dwell Time (s)')
        xlim([0 100])
        ylim([0 210])
        title(strcat(simpleDualDwells(i).name, " Right"), 'Interpreter','none')
        set(gca,'fontweight','bold', 'FontSize', 14)
    end
    sgtitle("Simple Dual Latch Per Molecule Average Dwell times", 'FontSize', 20, 'Fontweight', 'bold')
    saveas(simpleDualAvgDwells, "Simple Dual Latch Per Molecule Average Dwell times", 'png')
end

%% Plotting Total Time Histograms
if plotTotalTimeHistograms
    for i = 1:size(dualData, 2)
        states(i).name = dualData(1, i);
        states(i).leftLow = [];
        states(i).leftHigh = [];
        states(i).rightLow = [];
        states(i).rightHigh = [];
        for j = 1:size(dualData{2, i}, 2)
            states(i).leftHigh = cat(1, states(i).leftHigh, dualData{3, i}{j}(dualData{2, i}{j}==max(dualData{2, i}{j})));
            states(i).leftLow = cat(1, states(i).leftLow, dualData{3, i}{j}(dualData{2, i}{j}==min(dualData{2, i}{j})));
            states(i).rightHigh = cat(1, states(i).rightHigh, dualData{5, i}{j}(dualData{4, i}{j}==max(dualData{4, i}{j})));
            states(i).rightLow = cat(1, states(i).rightLow, dualData{5, i}{j}(dualData{4, i}{j}==min(dualData{4, i}{j})));
        end
        histosLeft = figure();
        histogram(states(i).leftHigh, 'BinWidth', .01, 'DisplayName', 'HMM High States')
        hold on
        histogram(states(i).leftLow, 'BinWidth', .01, 'DisplayName', 'HMM Low States')
        framesHighFRET = size(states(i).leftHigh,1);
        framesLowFRET = size(states(i).leftLow,1);
        timeAvgFRET = framesHighFRET/(framesHighFRET + framesLowFRET);
        title(strcat(states(i).name{1}, ' brush Left (Cy5) Relative Fluorescence'), 'Interpreter','none')
        ylabel('Frames Spent at Relative Fluorescence')
        xlabel('Fluorescence Intensity (arbitrary units)')
        lgd = legend('Location', 'best');
        title(lgd, sprintf('Normalized Average Quenching: %.2f', 1-timeAvgFRET))
        saveas(histosLeft, strcat(states(i).name{1}, ' brush Left (Cy5) Relative Fluorescence Histogram'), 'png')
        
        histosRight = figure();
        histogram(states(i).rightHigh, 'BinWidth', .01, 'DisplayName', 'HMM High States')
        hold on
        histogram(states(i).rightLow, 'BinWidth', .01, 'DisplayName', 'HMM Low States')
        framesHighFRET = size(states(i).rightHigh,1);
        framesLowFRET = size(states(i).rightLow,1);
        timeAvgFRET = framesHighFRET/(framesHighFRET + framesLowFRET);
        title(strcat(states(i).name{1}, ' brush Right (Cy3) Relative Fluorescence'), 'Interpreter','none')
        ylabel('Frames Spent at Relative Fluorescence')
        xlabel('Fluorescence Intensity (arbitrary units)')
        lgd = legend('Location', 'best');
        title(lgd, sprintf('Normalized Average Quenching: %.2f', 1-timeAvgFRET))
        saveas(histosRight, strcat(states(i).name{1}, ' brush Right (Cy3) Relative Fluorescence Histogram'), 'png')
    end
end

%% Plotting Dwell Time Histograms
if doPlotHists
    % bins = 0:1:max(singleLeftTlow(1, :));
    % Simple Dual Avg. Dwell Times
    simpleDualAvgDwells = figure('Position', [10 100 1900 900]);
    for i=1:numSamples
        bins = 0:1:max(simpleDualDwells(i).leftLow(1, :));
        subplot(2,numSamples,i*2-1)
        histogram(simpleDualDwells(i).leftLow, bins)
        grid on
        xlabel('Mean Low Dwell Time (s)')
        ylabel('Mean High Dwell Time (s)')
        xlim([0 100])
        ylim([0 210])
        title(strcat(simpleDualDwells(i).name, " Left"), 'Interpreter','none')
        set(gca,'fontweight','bold', 'FontSize', 14)
        
        subplot(2,numSamples,i*2)
        scatter(simpleDualDwells(i).rightMeanTlow, simpleDualDwells(i).rightMeanThigh, 'LineWidth',2)
        grid on
        xlabel('Mean Low Dwell Time (s)')
        ylabel('Mean High Dwell Time (s)')
        xlim([0 100])
        ylim([0 210])
        title(strcat(simpleDualDwells(i).name, " Right"), 'Interpreter','none')
        set(gca,'fontweight','bold', 'FontSize', 14)
    end
    sgtitle("Simple Dual Latch Per Molecule Average Dwell times", 'FontSize', 20, 'Fontweight', 'bold')
    saveas(simpleDualAvgDwells, "Simple Dual Latch Per Molecule Average Dwell times", 'png')
    xlabel('Time (s)')
    title(sprintf('Low Dwell Times (N = %d)', numLowDwells))
    histogram(singleLeftThigh, bins)
    xlabel('Time (s)')
    title(sprintf('High Dwell Times (N = %d)', numHighDwells))
    sgtitle(strcat(side, " Side Dwell Times (", bp, ' bp', movieLength, sprintf(", %d Hz)", frameRate)), 'Interpreter','none')
    saveas(f1, strcat(bp, " Histograms.png"))
    
    %% Low dwell time greater than t histogram
    ehist = zeros(1,size(bins,2));
    for i = 1:size(bins,2)
        ehist(i) = sum(singleLeftTlow(1,:) >= bins(i))/size(singleLeftTlow, 2);
    end
end

%% Plotting and Fitting Cumulative Sums
% Default Fitting Parameters
ub = '1,10,10';
lb = '0,0.001,0.001';
guess = '0.5,1.5,.05';
annealTemp=20;
fitParams = {ub, lb, guess, annealTemp};
singleFitParams = {'10', '0', '0.1', annealTemp};

if plotBrushesSamples
    if plotSimpleRates
        % Binding
        cmap = linspecer(4);
        fourSingleFitParams = {singleFitParams; singleFitParams; singleFitParams; singleFitParams};
        fourDoubleFitParams = {fitParams; fitParams; fitParams; fitParams};
        fourMarkers = {'ks', 'ko', 'k*', 'kx'};
        numSamples = size(threeStateDwells, 2);
        plotDwells = cell(numSamples, 1);
        names = cell(numSamples, 1);
        for i = 1:numSamples
            plotDwells{i} = threeStateDwells(i).middleLeftDwells(1, :);
            names{i} = threeStateDwells(i).name{1};
        end
        leftBindRates = plotNmixed(plotDwells, cutoffFraction, fourDoubleFitParams, names, cmap,...
            fourMarkers, foursampleOutsideBox,foursampleLegendSpecs, "Left Time to Bind");
        xlim([0 15])
        ylim([1e-2 1])
        saveas(gcf, strcat('plotRates/', 'Brushes Left Time to Bind', ".png"))
        
        plotDwells = cell(numSamples, 1);
        names = cell(numSamples, 1);
        for i = 1:numSamples
            plotDwells{i} = threeStateDwells(i).middleRightDwells(1, :);
            names{i} = threeStateDwells(i).name{1};
        end
        rightBindRates = plotNmixed(plotDwells, cutoffFraction, fourDoubleFitParams, names, cmap,...
            fourMarkers, foursampleOutsideBox,foursampleLegendSpecs, "Right Time to Bind");
        xlim([0 15])
        ylim([1e-2 1])
        saveas(gcf, strcat('plotRates/', 'Brushes Right Time to Bind', ".png"))
        
        % Dissociation
        for i = 1:numSamples
            plotDwells{i} = threeStateDwells(i).leftDwells(1, :);
            names{i} = threeStateDwells(i).name{1};
        end
        leftDissociateRates = plotNmixed(plotDwells, cutoffFraction, fourDoubleFitParams, names, cmap,...
            fourMarkers, foursampleOutsideBox,foursampleLegendSpecs, "Left Time to Dissociate");
        xlim([0 15])
        ylim([1e-2 1])
        saveas(gcf, strcat('plotRates/', 'Brushes Left Time to Dissociate', ".png"))
        
        plotDwells = cell(numSamples, 1);
        names = cell(numSamples, 1);
        for i = 1:numSamples
            plotDwells{i} = threeStateDwells(i).rightDwells(1, :);
            names{i} = threeStateDwells(i).name{1};
        end
        rightDissociateRates = plotNmixed(plotDwells, cutoffFraction, fourDoubleFitParams, names, cmap,...
            fourMarkers, foursampleOutsideBox,foursampleLegendSpecs, "Right Time to Dissociate");
        xlim([0 15])
        ylim([1e-2 1])
        saveas(gcf, strcat('plotRates/', 'Brushes Right Time to Dissociate', ".png"))
        rateStruct = struct();
        for i = 1:numSamples                
            rateStruct(i).name = rightDissociateRates{i,1};
            rateStruct(i).rightDissociateRates = rightDissociateRates{i,2};
            rateStruct(i).leftDissociateRates = leftDissociateRates{i,2};
            rateStruct(i).rightBindRates = rightBindRates{i,2};
            rateStruct(i).leftBindRates = leftBindRates{i,2};

            rateStruct(i).rightDissociateErrors = rightDissociateRates{i,3};
            rateStruct(i).leftDissociateErrors = leftDissociateRates{i,3};
            rateStruct(i).rightBindErrors = rightBindRates{i,3};
            rateStruct(i).leftBindErrors = leftBindRates{i,3};
        end
        save('simpleRates.mat', 'rateStruct');
    end
    if plotNonMarkov
        for sample = 1:numSamples
            % Plot Left Markov
            cmap = linspecer(2);
            twoSingleFitParams = {singleFitParams; singleFitParams};
            names = {'$K_{LUL}$', '$K_{RUL}$'};
            twoMarkers = {'ks', 'ko'};
            plotDwells = {threeStateDwells(sample).middleLeftRebindDwells(1, :), threeStateDwells(sample).middleRightLeftTransitionDwells(1, :)};
            plotNmixed(plotDwells, cutoffFraction, twoSingleFitParams, names, cmap, twoMarkers, foursampleOutsideBox, foursampleLegendSpecs, strcat(threeStateDwells(sample).name, ' Left Complex Rates'))
            ylim([1e-2 1])
            xlim([0 15])
            saveas(gcf, strcat('plotRates/complexRates/', threeStateDwells(sample).name, 'Left Complex Rates', ".png"))
            
            % Plot Right Markov
            cmap = linspecer(2);
            twoSingleFitParams = {singleFitParams; singleFitParams};
            names = {'$K_{RUR}$', '$K_{LUR}$'};
            twoMarkers = {'ks', 'ko'};
            plotDwells = {threeStateDwells(sample).middleRightRebindDwells(1, :), threeStateDwells(sample).middleLeftRightTransitionDwells(1, :)};
            plotNmixed(plotDwells, cutoffFraction, twoSingleFitParams, names, cmap, twoMarkers, foursampleOutsideBox, foursampleLegendSpecs, strcat(threeStateDwells(sample).name, ' Right Complex Rates'))
            ylim([1e-2 1])
            xlim([0 15])
            saveas(gcf, strcat('plotRates/complexRates/', threeStateDwells(sample).name, 'Right Complex Rates', ".png"))
        end
    end
end

if printHypothesisTest
    for i =1:4
        p1 = hypothesisTesting(threeStateDwells(i).middleLeftDwells(1,:), cutoffFraction);
        if (p1>0.05)
            disp(strcat(threeStateDwells(i).name, " Left Bind Dwells"))
        end
        p2 = hypothesisTesting(threeStateDwells(i).middleRightDwells(1,:), cutoffFraction);
        if (p2>0.05)
            disp(strcat(threeStateDwells(i).name, " Right Bind Dwells"))
        end
        p3 = hypothesisTesting(threeStateDwells(i).leftDwells(1,:), cutoffFraction);
        if (p3>0.05)
            disp(strcat(threeStateDwells(i).name, " Left Dwells"))
        end
        p4 = hypothesisTesting(threeStateDwells(i).rightDwells(1,:), cutoffFraction);
        if (p4>0.05)
            disp(strcat(threeStateDwells(i).name, " Right Dwells"))
        end
    end
end

if barPlots
    %% Make Brushless Rate Bar Plots
    load('rates.mat', 'fittedRates');
    % Extract names, rates, and errors for all entries
    primaryNames = fittedRates(1:4,1);
    primaryRates = cell2mat(fittedRates(1:4,2));
    primaryErrors = cell2mat(fittedRates(1:4,3));

    markovNames = fittedRates(5:end,1);
    markovRates = cell2mat(fittedRates(5:end,2));
    markovErrors = cell2mat(fittedRates(5:end,3));

    % Convert names to categorical for bar plotting
    primaryCatNames = categorical(primaryNames);
    markovCatNames = categorical(markovNames);

    % Plot Primary Rates
    figure;
    bar(primaryCatNames, primaryRates);
    hold on
    errorbar(primaryCatNames, primaryRates, primaryErrors, '.k', 'LineWidth', 2)
    ylabel('Rate (s^{-1})')
    title('Primary Rates')
    set(gca, 'FontWeight', 'bold', 'FontSize', 14)
    saveas(gcf, 'plotRates/Primary_Rates.png')

    % Plot Complex Rates
    figure;
    bar(markovCatNames, markovRates);
    hold on
    errorbar(markovCatNames, markovRates, markovErrors, '.k', 'LineWidth', 2)
    ylabel('Rate (s^{-1})')
    title('Complex Rates')
    set(gca, 'FontWeight', 'bold', 'FontSize', 14)
    saveas(gcf, 'plotRates/Complex_Rates.png')
end
toc

if brushBarPlots
    load('simpleRates.mat', 'rateStruct');
    cmap = linspecer(numSamples);
    plotNames = categorical(dualNames);

    % Left Dissociation Rates
    figure;
    for sample = 1:numSamples
        bar(plotNames(sample), rateStruct(sample).leftDissociateRates, 'FaceColor', cmap(sample, :));
        hold on
        errorbar(plotNames(sample), rateStruct(sample).leftDissociateRates, rateStruct(sample).leftDissociateErrors, '.k', 'LineWidth', 2)
    end
    set(gca, 'FontWeight', 'bold', 'FontSize', 20, 'TickLabelInterpreter', 'tex', 'LineWidth', 3.3);
    ylabel('Rate (s^{-1})')
    title('Left Dissociation Rates')
    saveas(gcf, 'plotRates/Brush_Left_Dissociation_Rates.png')

    % Right Dissociation Rates
    figure;
    for sample = 1:numSamples
        bar(plotNames(sample), rateStruct(sample).rightDissociateRates, 'FaceColor', cmap(sample, :));
        hold on
        errorbar(plotNames(sample), rateStruct(sample).rightDissociateRates, rateStruct(sample).rightDissociateErrors, '.k', 'LineWidth', 2)
    end
    set(gca, 'FontWeight', 'bold', 'FontSize', 20, 'TickLabelInterpreter', 'tex', 'LineWidth', 3.3);
    ylabel('Rate (s^{-1})')
    title('Right Dissociation Rates')
    saveas(gcf, 'plotRates/Brush_Right_Dissociation_Rates.png')

    % Binding Fractions
    figure;
    for sample = 2:numSamples
        bar(plotNames(sample), rateStruct(sample).leftBindRates(1), 'FaceColor', cmap(sample, :));
        hold on
        errorbar(plotNames(sample), rateStruct(sample).leftBindRates(1), rateStruct(sample).leftBindErrors(1), '.k', 'LineWidth', 2)
    end
    set(gca, 'FontWeight', 'bold', 'FontSize', 20, 'TickLabelInterpreter', 'tex', 'LineWidth', 3.3);
    ylabel('Fraction Fast')
    title('Left Binding Fractions')
    ylim([0 1])
    saveas(gcf, 'plotRates/Brush_Left_Binding_Fractions.png')

    figure;
    for sample = 2:numSamples
        bar(plotNames(sample), rateStruct(sample).rightBindRates(1), 'FaceColor', cmap(sample, :));
        hold on
        errorbar(plotNames(sample), rateStruct(sample).rightBindRates(1), rateStruct(sample).rightBindErrors(1), '.k', 'LineWidth', 2)
    end
    set(gca, 'FontWeight', 'bold', 'FontSize', 20, 'TickLabelInterpreter', 'tex', 'LineWidth', 3.3);
    ylabel('Fraction Fast')
    ylim([0 1])
    title('Right Binding Fractions')
    saveas(gcf, 'plotRates/Brush_Right_Binding_Fractions.png')

    %% Left Binding Rates
    % Fast Rates
    % figure;
    % fastSlow  = {'NA','Fast', 'Slow','Fast', 'Slow','Fast', 'Slow'};
    % fastSlow = categorical(fastSlow);
    % bars = {};
    % bars{1} = bar(plotNames(1), rateStruct(1).leftBindRates, 'FaceColor', cmap(1, :));
    % hold on
    % for sample = 2:numSamples
    %     fastSlowNumber = 2*(sample-1);
    %     bars{sample} = bar(plotNames(sample), [rateStruct(sample).leftBindRates(2), rateStruct(sample).leftBindRates(3)], 'FaceColor', cmap(sample, :));
    %     hold on
    % end
    % errorbar(bars{1}.XEndPoints, rateStruct(1).leftBindRates, rateStruct(1).leftBindErrors, '.k', 'LineWidth', 2)
    % for i = 2:numSamples
    %     errorbar(bars{i}(1).XEndPoints, rateStruct(i).leftBindRates(2), rateStruct(i).leftBindErrors(2), '.k', 'LineWidth', 2)
    %     errorbar(bars{i}(2).XEndPoints, rateStruct(i).leftBindRates(3), rateStruct(i).leftBindErrors(3), '.k', 'LineWidth', 2)
    % end
    % set(gca, 'FontWeight', 'bold', 'FontSize', 20, 'TickLabelInterpreter', 'tex', 'LineWidth', 3.3);
    % ylabel('Rate (s^{-1})')
    % title('Left Binding Fast Rates')
    % saveas(gcf, 'plotRates/Brush_Left_Binding_Fast_Rates.png')
    plotRateBar(cmap, plotNames, rateStruct, 'leftBind', 'Left Binding Rates', 'Brush_Left_Binding_Rates')
    plotRateBar(cmap, plotNames, rateStruct, 'rightBind', 'Right Binding Rates', 'Brush_Right_Binding_Rates')


    % % Slow Rates
    % figure;
    % bar(plotNames(1), rateStruct(1).leftBindRates, 'FaceColor', cmap(1, :));
    % hold on
    % errorbar(plotNames(1), rateStruct(1).leftBindRates, rateStruct(1).leftBindErrors, '.k', 'LineWidth', 2)
    % for sample = 2:numSamples
    %     bar(plotNames(sample), rateStruct(sample).leftBindRates(3), 'FaceColor', cmap(sample, :));
    %     hold on
    %     errorbar(plotNames(sample), rateStruct(sample).leftBindRates(3), rateStruct(sample).leftBindErrors(3), '.k', 'LineWidth', 2)
    % end
    % set(gca, 'FontWeight', 'bold', 'FontSize', 20, 'TickLabelInterpreter', 'tex', 'LineWidth', 3.3);
    % ylabel('Rate (s^{-1})')
    % title('Left Binding Slow Rates')
    % saveas(gcf, 'plotRates/Brush_Left_Binding_Slow_Rates.png')

    % %% Right Binding Rates
    % % Fast Rates
    % % Fast Rates
    % figure;
    % bar(plotNames(1), rateStruct(1).rightBindRates, 'FaceColor', cmap(1, :));
    % hold on
    % errorbar(plotNames(1), rateStruct(1).rightBindRates, rateStruct(1).rightBindErrors, '.k', 'LineWidth', 2)
    % for sample = 2:numSamples
    %     bar(plotNames(sample), rateStruct(sample).rightBindRates(2), 'FaceColor', cmap(sample, :));
    %     hold on
    %     errorbar(plotNames(sample), rateStruct(sample).rightBindRates(2), rateStruct(sample).rightBindErrors(2), '.k', 'LineWidth', 2)
    % end
    % set(gca, 'FontWeight', 'bold', 'FontSize', 20, 'TickLabelInterpreter', 'tex', 'LineWidth', 3.3);
    % ylabel('Rate (s^{-1})')
    % title('Right Binding Fast Rates')
    % saveas(gcf, 'plotRates/Brush_Right_Binding_Fast_Rates.png')

    % % Slow Rates
    % figure;
    % bar(plotNames(1), rateStruct(1).rightBindRates, 'FaceColor', cmap(1, :));
    % hold on
    % errorbar(plotNames(1), rateStruct(1).rightBindRates, rateStruct(1).rightBindErrors, '.k', 'LineWidth', 2)
    % for sample = 2:numSamples
    %     bar(plotNames(sample), rateStruct(sample).rightBindRates(3), 'FaceColor', cmap(sample, :));
    %     hold on
    %     errorbar(plotNames(sample), rateStruct(sample).rightBindRates(3), rateStruct(sample).rightBindErrors(3), '.k', 'LineWidth', 2)
    % end
    % set(gca, 'FontWeight', 'bold', 'FontSize', 20, 'TickLabelInterpreter', 'tex', 'LineWidth', 3.3);
    % ylabel('Rate (s^{-1})')
    % title('Right Binding Slow Rates')
    % saveas(gcf, 'plotRates/Brush_Right_Binding_Slow_Rates.png')
    
end

%% Plotting Functions
% Single exp plotting functions
function fittedRates = plotNSingleCutoff(dwells, cutoffFraction, params, names, colors, markers, boxSpecs, legendSpecs, title)
    Fig = figure('Position', [100 100 1200 800]);
    N = length(names);
    fittedRates = cell(N,3);
    str = cell(N);
    for i = 1:N
        plotCumDistsOneMinusCutoff(dwells{i}', cutoffFraction, Fig, markers{i}, names{i});
        [oneMinusFirstParams, error] = plotFitsCutoffOneMinus(Fig, dwells{i}', cutoffFraction, params{i}, colors(i, :),legendSpecs, title, names{i});
        fittedRates{i,1} = names{i};
        fittedRates{i,2} = oneMinusFirstParams;
        fittedRates{i,3} = error;
        str{i} = strcat(names{i}, sprintf(": k = %.3f \\pm %.3f s^{-1}",oneMinusFirstParams, error));
    end
    t = annotation('textbox',[0.495 0.13 0.293333333333334 0.13375],'String',str, 'Interpreter','tex');%,'FitBoxToText','on');
    t.FontSize = 18;
    t.FontWeight = 'bold';
    xlim([0 30])
    ylim([9e-3 1])
    set(gca, 'LineWidth', 3, 'FontSize', 22, 'FontWeight', 'bold')
    set(gca,'OuterPosition', [0 0 0.88 1])
    % saveas(Fig, strcat(title, ".png"))
end

function [params, error] = plotFitsCutoffOneMinus(fig, dwells, cutoffFraction, fitParams, clr, legendSpecs, Title, fitName)
    fig = figure(fig);
    [CumDist, CumDistTimes] = ecdf(dwells);
    CumDistTimes(1) = 0;
    dwellTimeCutoff = min(CumDistTimes(CumDist>=cutoffFraction));
    sortDwells = sort(dwells);
    [params, error] = fitSingleExpBootstrap(dwells(dwells<dwellTimeCutoff), fitParams, clr, fitName);
    lgd = legend('Location', legendSpecs, 'Interpreter','tex');
    % title(lgd, 'Brush Lengths')
    title(Title, 'Interpreter','none')
    xlabel('Time (s)')
    ylabel('1 - Cumulative Probability')
    set(gca, 'linewidth', 2, 'fontweight','bold', 'fontsize',16)
    % params;
end

function fig = plotCumDistsOneMinusCutoff(dwells, cutoff, fig, mk, displayName)
    % justLeftDwells = dwells(1,:);
    % Plot cumulative sum
    [totalCumDist, totalCumDistTimes] = ecdf(dwells);
    totalCumDistTimes(1) = 0;
    dwellTimeCutoff = min(totalCumDistTimes(totalCumDist>=cutoff));
    [CumDist, CumDistTimes] = ecdf(dwells(dwells<dwellTimeCutoff));
    CumDistTimes(1) = 0;
    if fig == "null"
        fig = figure();
    else
        fig = figure(fig);
    end
    semilogy(CumDistTimes, 1-CumDist, mk, 'DisplayName',displayName, 'LineWidth',4, 'MarkerSize',18)
    hold on
end

function [singleExpParams, paramError] = fitSingleExpBootstrap(dwellTimes, fitParams, clr, fitName)
    [userPDF, dataVar, fitVar, ~,~, ~]=PDFList('Single Exp', 'all', 0);
    lb = fitParams{2};
    ub = fitParams{1};
    guess = fitParams{3};
    annealTemp = fitParams{4};
    sortDwell = sort(dwellTimes);
    dwellsToFit = sortDwell(1:round(length(sortDwell)));
    [singleExpParams, logLi]= MEMLETCL(dwellsToFit, userPDF, dataVar, fitVar, lb,ub, guess,annealTemp); % Fit parameters and Log Likelihood output
    [bootstrapParams, bootstrapLogLi]= MEMLETCL(dwellsToFit, userPDF, dataVar, fitVar, lb,ub, guess,annealTemp, 1000); % Fit parameters and Log Likelihood output
    paramError = std(bootstrapParams);
    fitxvals=linspace(0,max(dwellsToFit),10000)'; %create variables for plotting along x
    oneMinus = exppdfoneminus(fitxvals, singleExpParams);
    semilogy(fitxvals,oneMinus, 'Color', clr, 'LineWidth', 6 , 'DisplayName',...
    fitName)
end

% Double Exp plotting functions

function fittedRates = plotNmixed(dwells, cutoffFraction, params, names, colors, markers, boxSpecs, legendSpecs, title)
    Fig = figure('Position', [50 50 1200 1080]);
    singleFitParams = {'10', '0', '0.1', 30};
    N = length(names);
    fittedRates = cell(N,3);
    str = cell(N);
    for i = 1:N
        plotCumDistsOneMinusCutoff(dwells{i}', cutoffFraction, Fig, markers{i}, names{i});
        pval = hypothesisTesting(dwells{i}, cutoffFraction)
        if pval > 0.05
            % Fit Single Exponential
            [oneMinusFirstParams, error] = plotFitsCutoffOneMinus(Fig, dwells{i}', cutoffFraction, singleFitParams, colors(i, :),legendSpecs, title, names{i});
            fittedRates{i,1} = names{i};
            fittedRates{i,2} = oneMinusFirstParams;
            fittedRates{i,3} = error;
            str{i} = strcat(names{i}, sprintf(": k = %.3f \\pm %.3f s^{-1}",oneMinusFirstParams, error));
        else
            % Fit Double Exponential
            [oneMinusFirstParams, error] = plotDoubleFitsCutoff(Fig, dwells{i}', cutoffFraction, params{i}, colors(i, :),legendSpecs, title, names{i});
            fittedRates{i,1} = names{i};
            fittedRates{i,2} = oneMinusFirstParams;
            fittedRates{i,3} = error;
            str{i} = strcat(names{i}, sprintf(": A = %.2f \\pm %.2f, k_1 = %.3f \\pm %.3f s^{-1}, k_2 = %.3f \\pm %.3f s^{-1}", ...
                oneMinusFirstParams(1),error(1), oneMinusFirstParams(2), error(2), oneMinusFirstParams(3), error(3)));
        end
    end
    t = annotation('textbox',[0.120833333333334 0.0126849894291755 0.690833333333332 0.210280126849895],'String',str, 'Interpreter','tex');%,'FitBoxToText','on');
    t.FontSize = 19;
    t.FontWeight = 'bold';
    xlim([0 30])
    ylim([9e-3 1])
    set(gca, 'LineWidth', 3.3, 'FontSize', 28, 'FontWeight', 'bold')
    set(gca,'OuterPosition', [0 0.22 0.88 0.78])
    % saveas(Fig, strcat(title, ".png"))
end

function [params, error] = plotDoubleFitsCutoff(fig, dwells, cutoffFraction, fitParams, clr, legendSpecs, Title, fitName)
    fig = figure(fig);
    [CumDist, CumDistTimes] = ecdf(dwells);
    CumDistTimes(1) = 0;
    dwellTimeCutoff = min(CumDistTimes(CumDist>=cutoffFraction));
    sortDwells = sort(dwells);
    [params, error] = fitDoubleExpBootstrap(dwells(dwells<dwellTimeCutoff), fitParams, clr, fitName);
    lgd = legend('Location', legendSpecs, 'Interpreter','tex');
    % title(lgd, 'Brush Lengths')
    title(Title, 'Interpreter','none')
    xlabel('Time (s)')
    ylabel('1 - Cumulative Probability')
    set(gca, 'linewidth', 2, 'fontweight','bold', 'fontsize',16)
    % params;
end

function [doubleExpParams, paramError] = fitDoubleExpBootstrap(dwellTimes, fitParams, clr, fitName)
    [userPDF, dataVar, fitVar, ~,~, ~]=PDFList('Double Exp (Independent)', 'all', 0);
    lb = fitParams{2};
    ub = fitParams{1};
    guess = fitParams{3};
    annealTemp = fitParams{4};
    sortDwell = sort(dwellTimes);
    dwellsToFit = sortDwell(1:round(length(sortDwell)));
    [doubleExpParams, logLi]= MEMLETCL(dwellsToFit, userPDF, dataVar, fitVar, lb,ub, guess,annealTemp); % Fit parameters and Log Likelihood output
    % Reorder parameters for actual fit so fastest rate is k1
    if doubleExpParams(2)<doubleExpParams(3)
        temp = doubleExpParams(2);
        doubleExpParams(2) = doubleExpParams(3);
        doubleExpParams(3) = temp;
        doubleExpParams(1) = 1-doubleExpParams(1);
    end
    [bootstrapParams, bootstrapLogLi]= MEMLETCL(dwellsToFit, userPDF, dataVar, fitVar, lb,ub, guess,annealTemp, 1000); % Fit parameters and Log Likelihood output
    % Re-order parameters in bootstrap matrix to keep fastest rate as k1
    misorderedRates = bootstrapParams(:,2)<bootstrapParams(:,3);
    orderedBootstrapParams = bootstrapParams;
    orderedBootstrapParams(misorderedRates,2) = bootstrapParams(misorderedRates,3);
    orderedBootstrapParams(misorderedRates,3) = bootstrapParams(misorderedRates,2);
    orderedBootstrapParams(misorderedRates,1) = 1-bootstrapParams(misorderedRates,1);
    paramError = std(bootstrapParams);
    fitxvals=linspace(0,max(dwellsToFit),10000)'; %create variables for plotting along x
    oneMinus = dbexppdfnotminoneminus(fitxvals, doubleExpParams);
    semilogy(fitxvals,oneMinus, 'Color', clr, 'LineWidth', 6 , 'DisplayName',...
        fitName)
end

function pval = hypothesisTesting(dwellTimes, cutoff)
    % Compute single exponential fit
    [userPDF, dataVar, fitVar, ~,~, ~]=PDFList('Single Exp', 'all', 0);
    annealTemp=15;
    fitParams = {'10', '0', '0.1', annealTemp};
    lb = fitParams{2};
    ub = fitParams{1};
    guess = fitParams{3};
    annealTemp = fitParams{4};
    [CumDist, CumDistTimes] = ecdf(dwellTimes);
    CumDistTimes(1) = 0;
    dwellTimeCutoff = min(CumDistTimes(CumDist>=cutoff));
    dwellSubset = dwellTimes(dwellTimes<dwellTimeCutoff);
    [singleExpParams, singleLogL] = MEMLETCL(dwellSubset', userPDF, dataVar, fitVar, lb,ub, guess,annealTemp); % Fit parameters and Log Likelihood output
    % Compute double exponential fit
    [userPDF, dataVar, fitVar, ~,~, ~]=PDFList('Double Exp (Independent)', 'all', 0);
    ub = '1,10,10';
    lb = '0,0.001,0.001';
    guess = '0.5,1,.1';
    [doubleExpParams, doubleLogL] = MEMLETCL(dwellSubset', userPDF, dataVar, fitVar, lb,ub, guess,annealTemp); % Fit parameters and Log Likelihood output
    delDF = 2; %% two constraints on double exp yields single exp, hence 2
    RLL=-2*(singleLogL-doubleLogL); % the log of the ratio of  the likelihoods
    pval=1-chi2cdf(RLL,delDF); %calculate a p-value from the chi2cdf
    % disp(singleLogL)
    % disp(doubleLogL)
    % disp(pval)
    % if plow == 1
    %     plow = '> 1 - 1e-16';
    % elseif plow == 0
    %     plow = '< 1e-16';
    % else
    %     plow = num2str(plow);
    % end
end

function [] = slider_plot(leftRaw, leftHMM, rightRaw, rightHMM, dualHMM, secPerFrame, titlestr)
    % Plot different plots according to slider location.
    S.fh = figure('units','pixels',...
        'position',[50 50 1360 870],...
        'menubar','none',...
        'name','slider_plot',...
        'numbertitle','off',...
        'resize','off');
    set(S.fh, 'Name', titlestr);
    S.ax1 = axes('unit','pix','position',[120 80 1200 200]);
    S.ax2 = axes('unit','pix','position',[120 80+280 1200 200]);
    S.ax3 = axes('unit','pix','position',[120 80+2*280 1200 200]);

    for a =[S.ax1 S.ax2 S.ax3]
        cla(a)
        hold(a, 'on')
    end

    seconds = 0:secPerFrame:(size(leftHMM{1}, 1)-1)*secPerFrame;
    plot(S.ax3, seconds, leftHMM{1}', '-','LineWidth', 2, 'Color',[0.4 0 0], DisplayName="Left Quenching HMM")
    xlabel(S.ax3,'Time (s)')
    ylabel(S.ax3,'Cy5 Emission')
    plot(S.ax3, seconds, leftRaw{1}', '-','Color', [0.8 0 0], DisplayName="Left Quenching Raw")
    title(S.ax3,sprintf("Trace %d", 1))

    plot(S.ax2, seconds, rightHMM{1}, '-','LineWidth', 2,'Color',[0 0 0.2], DisplayName="Right Quenching HMM")
    xlabel(S.ax2,'Time (s)')
    ylabel(S.ax2,'Cy3 Emission')
    plot(S.ax2, seconds, rightRaw{1}', '-', 'Color',[0 0 0.6], DisplayName="Right Quenching Raw")

    % Three State HMM PLot
    plot(S.ax1, seconds, dualHMM{1}', '-','LineWidth', 1,'Color',[0 0 0.2], DisplayName="Dual Quenching HMM")
    xlabel(S.ax1,'Time (s)')
    yticks(S.ax1,[0 1 2 3])
    ylim(S.ax1, [-.2 3.2])
    names = {'Both Bound'; 'Left Bound'; 'Unbound'; 'Right Bound'};
    set(S.ax1,'ytick',[0 1 2 3],'yticklabel',names)
    grid(S.ax1)
    legend(S.ax1, 'off')

    for a =[S.ax2 S.ax3]
        legend(a)
    end
    saveas(S.fh, "trace.png")

    S.sl = uicontrol('style','slide',...
    'unit','pix',...
    'position',[0 10 1200 30],...
    'min',1,'max',size(dualHMM,2),'val',1,...
    'sliderstep',[1/size(dualHMM,2) 1/size(dualHMM,2)],...
    'callback',{@sl_call,S, leftRaw, leftHMM, rightRaw, rightHMM, dualHMM, secPerFrame});
end

function [] = sl_call(varargin)
    % Callback for the slider.
    [h,S] = varargin{[1,3]};  % calling handle and data structure.
    leftRaw = varargin{4};
    leftHMM = varargin{5};
    rightRaw = varargin{6};
    rightHMM = varargin{7};
    duaHMM = varargin{8};
    secPerFrame = varargin{9};
    % cla
    traceNum = round(get(h,'value'));
    seconds = 0:secPerFrame:(size(leftHMM{traceNum}, 1)-1)*secPerFrame;

    for a =[S.ax1 S.ax2 S.ax3]
        cla(a)
        hold(a, 'on')
    end

    plot(S.ax3, seconds, leftHMM{traceNum}', '-','LineWidth', 2, 'Color',[0.4 0 0], DisplayName="Left Quenching HMM")
    xlabel(S.ax3,'Time (s)')
    ylabel(S.ax3,'Cy5 Emission')
    plot(S.ax3, seconds, leftRaw{traceNum}', '-','Color', [0.8 0 0], DisplayName="Left Quenching Raw")
    title(S.ax3,sprintf("Trace %d", traceNum))

    plot(S.ax2, seconds, rightHMM{traceNum}', '-','LineWidth', 2,'Color',[0 0 0.2], DisplayName="Right Quenching HMM")
    xlabel(S.ax2,'Time (s)')
    ylabel(S.ax2,'Cy3 Emission')
    plot(S.ax2, seconds, rightRaw{traceNum}', '-', 'Color',[0 0 0.6], DisplayName="Right Quenching Raw")

    % Three State HMM PLot
    plot(S.ax1, seconds, duaHMM{traceNum}', '-','LineWidth', 1,'Color',[0 0 0.2], DisplayName="Dual Quenching HMM")
    xlabel(S.ax1,'Time (s)')
    yticks(S.ax1,[0 1 2 3])
    ylim(S.ax1, [-.2 3.2])
    names = {'Both Bound'; 'Left Bound'; 'Unbound'; 'Right Bound'};
    set(S.ax1,'ytick',[0 1 2 3],'yticklabel',names)
    grid(S.ax1, 'on')
    legend(S.ax1, 'off')

    for a =[S.ax2 S.ax3]
        legend(a)
    end
end

function [] = fret_slider_plot(leftRaw, leftHMM, rightRaw, rightHMM, dualHMM, secPerFrame, titlestr)
    % Plot different plots according to slider location.
    % Normalizes and computes a composite trace from the Cy3 and Cy5 traces
    S.fh = figure('units','pixels',...
        'position',[50 50 1340 870],...
        'menubar','none',...
        'name','slider_plot',...
        'numbertitle','off',...
        'resize','off');
    set(S.fh, 'Name', titlestr);

    % S.ax1 = subplot(1,8,[1:7]);
    % S.ax2 = subplot(1,8,[1:7]);
    % S.ax3 = subplot(1,8,[1:7]);
    S.ax1 = axes('unit','pix','position',[120 80 1000 200]);
    S.ax2 = axes('unit','pix','position',[120 80+280 1000 200]);
    S.ax3 = axes('unit','pix','position',[120 80+2*280 1000 200]);

    S.ax4 = axes('unit','pix','position',[120+1010 80 180 200]);
    S.ax5 = axes('unit','pix','position',[120+1010 80+280 180 200]);
    S.ax6 = axes('unit','pix','position',[120+1010 80+2*280 180 200]);

    for a =[S.ax1 S.ax2 S.ax3]
        cla(a)
        hold(a, 'on')
    end

    seconds = 0:secPerFrame:(size(leftHMM{1}, 1)-1)*secPerFrame;
    plot(S.ax3, seconds, leftHMM{1}', '-','LineWidth', 2, 'Color',[0.4 0 0], DisplayName="Left Quenching HMM")
    xlabel(S.ax3,'Time (s)')
    ylabel(S.ax3,'Cy5 Emission')
    plot(S.ax3, seconds, leftRaw{1}', '-','Color', [0.8 0 0], DisplayName="Left Quenching Raw")
    title(S.ax3,sprintf("Trace %d", 1))

    histogram(S.ax6, leftRaw{1}, 'Normalization','probability', 'FaceColor',[0.8 0 0])
    axis(S.ax6, 'off')
    view(S.ax6, [90 -90])

    plot(S.ax2, seconds, rightHMM{1}, '-','LineWidth', 2,'Color',[0 0 0.2], DisplayName="Right Quenching HMM")
    xlabel(S.ax2,'Time (s)')
    ylabel(S.ax2,'Cy3 Emission')
    plot(S.ax2, seconds, rightRaw{1}', '-', 'Color',[0 0 0.6], DisplayName="Right Quenching Raw")
    
    histogram(S.ax5, rightRaw{1}, 'Normalization','probability', 'FaceColor',[0 0 0.6])
    axis(S.ax5, 'off')
    view(S.ax5, [90 -90])


    % fretRaw = leftRaw{1} ./ (rightRaw{1} + leftRaw{1});
    % fretHMM = leftHMM{1} ./ (rightHMM{1} + leftHMM{1});
    % midHMM = max(leftHMM{1}) / (max(rightHMM{1}) + max(leftHMM{1}));

    leftNormRaw = (leftRaw{1}-min(leftHMM{1}))./(max(leftHMM{1}) - min(leftHMM{1}));
    leftNormHMM = (leftHMM{1}-min(leftHMM{1}))./(max(leftHMM{1}) - min(leftHMM{1}));
    rightNormRaw = (rightRaw{1}-min(rightHMM{1}))./(max(rightHMM{1}) - min(rightHMM{1}));
    rightNormHMM = (rightHMM{1}-min(rightHMM{1}))./(max(rightHMM{1}) - min(rightHMM{1}));

    fretRaw = (leftNormRaw - rightNormRaw);
    fretHMM = (leftNormHMM - rightNormHMM); 
    midHMM = (max(leftNormHMM) - max(rightNormHMM));


    % Three State Composite Plot
    plot(S.ax1, seconds, fretRaw, '-','LineWidth', 1,'Color',[0 0 0], DisplayName="Dual Quenching Raw")
    hold on
    plot(S.ax1, seconds, fretHMM, '-','LineWidth', 2,'Color',[0.65 0 0.65], DisplayName="Dual Quenching HMM")
    xlabel(S.ax1,'Time (s)')
    tickvals = [min(fretHMM) midHMM max(fretHMM)];
    yticks(S.ax1,tickvals)
    ylim(S.ax1, [min(fretRaw) max(fretRaw)])
    names = {'Left Bound'; 'Unbound'; 'Right Bound'};
    set(S.ax1,'ytick',tickvals,'yticklabel',names)
    grid(S.ax1)
    legend(S.ax1, 'off')

    histogram(S.ax4, fretRaw, 'Normalization','probability', 'FaceColor',[0 0 0])
    axis(S.ax4, 'off')
    view(S.ax4, [90 -90]);


    for a =[S.ax2 S.ax3]
        legend(a)
    end
    saveas(S.fh, "trace.png")

    S.sl = uicontrol('style','slide',...
    'unit','pix',...
    'position',[0 10 1200 30],...
    'min',1,'max',size(dualHMM,2),'val',1,...
    'sliderstep',[1/size(dualHMM,2) 1/size(dualHMM,2)],...
    'callback',{@sl2_call,S, leftRaw, leftHMM, rightRaw, rightHMM, dualHMM, secPerFrame});
end

function [] = sl2_call(varargin)
    % Callback for the slider.
    [h,S] = varargin{[1,3]};  % calling handle and data structure.
    leftRaw = varargin{4};
    leftHMM = varargin{5};
    rightRaw = varargin{6};
    rightHMM = varargin{7};
    duaHMM = varargin{8};
    secPerFrame = varargin{9};
    % cla
    traceNum = round(get(h,'value'));
    seconds = 0:secPerFrame:(size(leftHMM{traceNum}, 1)-1)*secPerFrame;

    for a =[S.ax1 S.ax2 S.ax3 S.ax4 S.ax5 S.ax6]
        cla(a)
        hold(a, 'on')
    end

    plot(S.ax3, seconds, leftHMM{traceNum}', '-','LineWidth', 2, 'Color',[0.4 0 0], DisplayName="Left Quenching HMM")
    xlabel(S.ax3,'Time (s)')
    ylabel(S.ax3,'Cy5 Emission')
    plot(S.ax3, seconds, leftRaw{traceNum}', '-','Color', [0.8 0 0], DisplayName="Left Quenching Raw")
    title(S.ax3,sprintf("Trace %d", traceNum))

    histogram(S.ax6, leftRaw{traceNum}, 'Normalization','probability', 'FaceColor',[0.8 0 0])
    axis(S.ax6, 'off')
    view(S.ax6, [90 -90])

    plot(S.ax2, seconds, rightHMM{traceNum}', '-','LineWidth', 2,'Color',[0 0 0.2], DisplayName="Right Quenching HMM")
    xlabel(S.ax2,'Time (s)')
    ylabel(S.ax2,'Cy3 Emission')
    plot(S.ax2, seconds, rightRaw{traceNum}', '-', 'Color',[0 0 0.6], DisplayName="Right Quenching Raw")

    histogram(S.ax5, rightRaw{traceNum}, 'Normalization','probability', 'FaceColor',[0 0 0.6])
    axis(S.ax5, 'off')
    view(S.ax5, [90 -90])

    % fretRaw = leftRaw{traceNum} ./ (rightRaw{traceNum} + leftRaw{traceNum});
    % fretHMM = leftHMM{traceNum} ./ (rightHMM{traceNum} + leftHMM{traceNum});
    % midHMM = max(leftHMM{traceNum}) / (max(rightHMM{traceNum}) + max(leftHMM{traceNum}));

    leftNormRaw = (leftRaw{traceNum}-min(leftHMM{traceNum}))./(max(leftHMM{traceNum}) - min(leftHMM{traceNum}));
    leftNormHMM = (leftHMM{traceNum}-min(leftHMM{traceNum}))./(max(leftHMM{traceNum}) - min(leftHMM{traceNum}));
    rightNormRaw = (rightRaw{traceNum}-min(rightHMM{traceNum}))./(max(rightHMM{traceNum}) - min(rightHMM{traceNum}));
    rightNormHMM = (rightHMM{traceNum}-min(rightHMM{traceNum}))./(max(rightHMM{traceNum}) - min(rightHMM{traceNum}));
    
    fretRaw = (leftNormRaw - rightNormRaw);
    fretHMM = (leftNormHMM - rightNormHMM); 
    midHMM = (max(leftNormHMM) - max(rightNormHMM));

    % Three State Composite Plot
    plot(S.ax1, seconds, fretRaw, '-','LineWidth', 1,'Color',[0 0 0], DisplayName="Dual Quenching Raw")
    hold on
    plot(S.ax1, seconds, fretHMM, '-','LineWidth', 2,'Color',[0.65 0 0.65], DisplayName="Dual Quenching HMM")
    xlabel(S.ax1,'Time (s)')
    tickvals = [min(fretHMM) midHMM max(fretHMM)];
    yticks(S.ax1,tickvals)
    ylim(S.ax1, [min(fretRaw) max(fretRaw)])
    names = {'Left Bound'; 'Unbound'; 'Right Bound'};
    set(S.ax1,'ytick',tickvals,'yticklabel',names)
    grid(S.ax1)
    legend(S.ax1, 'off')

    histogram(S.ax4, fretRaw, 'Normalization','probability', 'FaceColor',[0 0 0])
    axis(S.ax4, 'off')
    view(S.ax4, [90 -90])

    for a =[S.ax2 S.ax3]
        legend(a)
    end
end

function plotRateBar(cmap, plotNames, rateStruct, field, plotTitle, savename)
    figure;
    bars = {};
    bars{1} = bar(plotNames(1), getfield(rateStruct(1), field+"Rates"), 'FaceColor', cmap(1, :));
    hold on
    for sample = 2:4
        barval = getfield(rateStruct(sample), field+"Rates");
        bars{sample} = bar(plotNames(sample), [barval(2), barval(3)], 'FaceColor', cmap(sample, :));
        hold on
    end
    errorbar(bars{1}.XEndPoints, getfield(rateStruct(1), field+"Rates"), getfield(rateStruct(1), field+"Errors"), '.k', 'LineWidth', 2)
    for i = 2:4
        barval = getfield(rateStruct(i), field+"Rates");
        errorval = getfield(rateStruct(i), field+"Errors");
        errorbar(bars{i}(1).XEndPoints, barval(2), errorval(2), '.k', 'LineWidth', 2)
        errorbar(bars{i}(2).XEndPoints, barval(3), errorval(3), '.k', 'LineWidth', 2)
    end
    set(gca, 'FontWeight', 'bold', 'FontSize', 20, 'TickLabelInterpreter', 'tex', 'LineWidth', 3.3);
    ylabel('Rate (s^{-1})')
    title(plotTitle)
    saveas(gcf, "plotRates/"+savename +".png")
end