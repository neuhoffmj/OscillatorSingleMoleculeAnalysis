clear;
close all;

dualframeRate = 10; % Hz or fps
movieLength = "5 Minutes";
dualSecPerFrame = 1/dualframeRate;

doExcludeTraces = 1;
plotDwellScatter = 0;
doSliderPlot = 0;
plotIndividualCumSum = 0;
plotTotalTimeHistograms = 0;
doPlotHists = 0;
plotFourSample = 1;
plotFourSampleThreeState = 0;
plotNonMarkov = 0;
plotOneSample = 0;
plotOneSampleThreeState = 0;
numBins = 20;

doPlotTraces = 0; % Plot a Selection of Traces?
plotTraces = [3 12 13 20]; %[1 4 10 15]; % Which Traces do you want to plot?

%% Loading Data
dualData = {'0_0' '0_10' '30_0' '30_10';...
    [], [], [],[]; % Left (Cy5) Path [2]
    [], [], [],[]; % Left (Cy5) Raw "FRET" [3]
    [], [], [],[]; % Right (Cy3) Path [4] 
    [], [], [],[]; % Right (Cy3) Raw "FRET" [5] 
    [], [], [],[]}; % Three State [6]
numDualSamples = size(dualData, 2);
dualNames = ["0_0", "0_10", "30_0", "30_10"];
for i=1:length(dualNames)
    simpleDualDwells(i).name = dualNames(i);
    threeStateDwells(i).name = dualNames(i);
end

outsideBoxSpecs = [0.68 0.15 0.3 0.16];

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

excludeAllButFastFlops  = {[7, 49, 50, 53, 56], ...
    [20, 29, 39, 50, 56, 57, 10, 35, 40, 41, 49, 16, 26], ...
    [8, 12, 65, 17, 68, 9, 10],...
    [8, 16, 22, 31, 48, 57, 77, 85, 94, 106, 97]};


for i=1:size(dualData, 2)
    % Loading in paths and excluding traces
    cy5File = strcat("vbOutput_cy5_", dualData{1,i}, ".mat");
    HMMleft = load(cy5File, "path");
    HMMleft = HMMleft.path;
    Rawleft = load(cy5File, "FRET");
    Rawleft = Rawleft.FRET;
    cy3File = strcat("vbOutput_cy3_", dualData{1,i}, ".mat");
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
        slider_plot(dualData{3, i}, dualData{2, i}, dualData{5, i}, dualData{4, i}, dualData{6, i}, dualSecPerFrame, dualData{1,i})
    end
end

if doPlotTraces
    for i = plotTraces
        % Make Raw and HMM Plots
        plFig = figure('Position', [100 100 1300 400]);
        seconds = 0:secPerFrame:(size(quenchHMMSubset{i}, 1)-1)*secPerFrame;
        if strcmp(side, 'Left')
            clrHMM = [0.4 0 0];
            clrRaw = [0.8 0 0];
            fluor = 'Cy5';
        else
            clrHMM = [0 0.2 0];
            clrRaw = [0 0.6 0];
            fluor = 'Cy3';
        end
        plot(seconds, quenchHMMSubset{i}', '-','LineWidth', 2, 'Color',clrHMM, DisplayName=strcat(side, " Quenching HMM"))
        xlabel('Time (s)', 'fontweight','bold','fontsize',16)
        ylabel(strcat(fluor, ' Emission'), 'fontweight','bold','fontsize',16)
        set(gca, 'linewidth', 2, 'fontweight','bold', 'fontsize',16)
        hold on
        plot(seconds, quenchRawSubset{i}', '-','Color', clrRaw, DisplayName=strcat(side, " Quenching Raw"))
        title(sprintf("Trace %d", i))
        legend()
        saveas(plFig,sprintf("Trace %d.png", i))
    end
end
%% Calculating Dwell Times
for i=1:size(dualData, 2)
    % Need to run both dwell time fitting programs
    % [leftThigh, leftTlow, ~, ~, ~, ~] = eb_dwelltimes(dualData{2, i}', dualSecPerFrame);% = HMMleft;
    % simpleDualDwells(i).leftLow = leftTlow;
    % simpleDualDwells(i).leftHigh = leftThigh;
    % 
    % [rightThigh, rightTlow, ~, ~, ~, ~]= eb_dwelltimes(dualData{4, i}', dualSecPerFrame);% = HMMright;
    % simpleDualDwells(i).rightLow = rightTlow;
    % simpleDualDwells(i).rightHigh = rightThigh;

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
    threeStateDwells(i).middleLeftRightTransitionDwells = middleLeftRightTransitionDwells;
    threeStateDwells(i).middleRightLeftTransitionDwells = middleRightLeftTransitionDwells;
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
    for i=1:numDualSamples
        subplot(2,numDualSamples,i*2-1)
        scatter(simpleDualDwells(i).leftMeanTlow, simpleDualDwells(i).leftMeanThigh, 'LineWidth',2)
        grid on
        xlabel('Mean Low Dwell Time (s)')
        ylabel('Mean High Dwell Time (s)')
        xlim([0 100])
        ylim([0 210])
        title(strcat(simpleDualDwells(i).name, " Left"), 'Interpreter','none')
        set(gca,'fontweight','bold', 'FontSize', 14)
    
        subplot(2,numDualSamples,i*2)
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
    for i=1:numDualSamples
        bins = 0:1:max(simpleDualDwells(i).leftLow(1, :));
        subplot(2,numDualSamples,i*2-1)
        histogram(simpleDualDwells(i).leftLow, bins)
        grid on
        xlabel('Mean Low Dwell Time (s)')
        ylabel('Mean High Dwell Time (s)')
        xlim([0 100])
        ylim([0 210])
        title(strcat(simpleDualDwells(i).name, " Left"), 'Interpreter','none')
        set(gca,'fontweight','bold', 'FontSize', 14)
    
        subplot(2,numDualSamples,i*2)
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
annealTemp=30;  
fitParams = {ub, lb, guess, annealTemp};
singleFitParams = {'10', '0', '0.1', annealTemp};


if plotFourSample
    % Binding
    cmap = linspecer(4);
    numSamples = size(simpleDualDwells, 2);
    plotDwells = cell(numSamples, 1);
    names = cell(numSamples, 1);
    for i = 1:numSamples
        plotDwells{i} = simpleDualDwells(i).leftHigh;
        names{i} = simpleDualDwells(i).name{1};
    end
    plotNDoubleExpOneMinus(plotDwells, {fitParams; fitParams; fitParams; fitParams}, names, cmap,...
        {'ks', 'ko', 'k*', 'kx'}, [0.48 0.15 0.35 0.16], "Left Time to Bind")
    plotDwells = cell(numSamples, 1);
    names = cell(numSamples, 1);
    for i = 1:numSamples
        plotDwells{i} = simpleDualDwells(i).rightHigh;
        names{i} = simpleDualDwells(i).name{1};
    end
    plotNDoubleExpOneMinus(plotDwells, {fitParams; fitParams; fitParams; fitParams}, names, cmap,...
        {'ks', 'ko', 'k*', 'kx'}, [0.48 0.15 0.35 0.16], "Right Time to Bind")

    % Dissociation
    for i = 1:numSamples
        plotDwells{i} = simpleDualDwells(i).leftLow;
    end
    plotNSingleExpOneMinus(plotDwells, {singleFitParams; singleFitParams; singleFitParams; singleFitParams}, names, cmap,...
        {'ks', 'ko', 'k*', 'kx'}, [0.7 0.15 0.17 0.16], "Left Time To Dissociate")
    for i = 1:numSamples
        plotDwells{i} = simpleDualDwells(i).rightLow;
    end
    plotNSingleExpOneMinus(plotDwells, {singleFitParams; singleFitParams; singleFitParams; singleFitParams}, names, cmap,...
        {'ks', 'ko', 'k*', 'kx'}, [0.7 0.15 0.17 0.16], "Right Time To Dissociate")
end

if plotFourSampleThreeState
    % Binding
    cmap = linspecer(4);
    numSamples = size(simpleDualDwells, 2);
    plotDwells = cell(numSamples, 1);
    names = cell(numSamples, 1);
    if plotNonMarkov
        for i = 1:numSamples
            plotDwells{i} = threeStateDwells(i).middleLeftRebindDwells(1, :);
            names{i} = simpleDualDwells(i).name{1};
        end
        plotNDoubleExpOneMinus(plotDwells, {fitParams; fitParams; fitParams; fitParams}, names, cmap,...
            {'ks', 'ko', 'k*', 'kx'}, [0.48 0.15 0.35 0.16], "Middle State Dwell (Left Time to Rebind)")
        xlim([0 20])
        % ylim([4e-4 1])
        for i = 1:numSamples
            plotDwells{i} = threeStateDwells(i).middleLeftRightTransitionDwells(1, :);
        end
        plotNDoubleExpOneMinus(plotDwells, {fitParams; fitParams; fitParams; fitParams}, names, cmap,...
            {'ks', 'ko', 'k*', 'kx'}, [0.48 0.15 0.35 0.16], "Middle State Dwell (Left Right Transition)")
        xlim([0 20])
        % ylim([1e-3 1])
        for i = 1:numSamples
            plotDwells{i} = threeStateDwells(i).middleRightLeftTransitionDwells(1, :);
        end
        plotNDoubleExpOneMinus(plotDwells, {fitParams; fitParams; fitParams; fitParams}, names, cmap,...
            {'ks', 'ko', 'k*', 'kx'}, [0.48 0.15 0.35 0.16], "Middle State Dwell (Right Left Transistion)")
        xlim([0 20])
        % ylim([1e-3 1])
        for i = 1:numSamples
            plotDwells{i} = threeStateDwells(i).middleRightRebindDwells(1, :);
        end
        plotNDoubleExpOneMinus(plotDwells, {fitParams; fitParams; fitParams; fitParams}, names, cmap,...
            {'ks', 'ko', 'k*', 'kx'}, [0.48 0.15 0.35 0.16], "Middle State Dwell (Right Time to Rebind)")
        xlim([0 20])
        % ylim([1e-3 1])
    end
    for i = 1:numSamples
        plotDwells{i} = threeStateDwells(i).middleRightDwells(1, :);
        names{i} = simpleDualDwells(i).name{1};
    end
    plotNDoubleExpOneMinus(plotDwells, {fitParams; fitParams; fitParams; fitParams}, names, cmap,...
        {'ks', 'ko', 'k*', 'kx'}, [0.68 0.15 0.3 0.16], "Middle Right Dwell")
    set(gca,'OuterPosition', [0 0 0.88 1])
    xlim([0 20])
    ylim([3e-3 1])
    for i = 1:numSamples
        plotDwells{i} = threeStateDwells(i).middleLeftDwells(1, :);
    end
    plotNDoubleExpOneMinus(plotDwells, {fitParams; fitParams; fitParams; fitParams}, names, cmap,...
        {'ks', 'ko', 'k*', 'kx'}, [0.68 0.15 0.3 0.16], "Middle Left Dwell")
    set(gca,'OuterPosition', [0 0 0.88 1])
    xlim([0 20])
    ylim([3e-3 1])
end


if plotOneSampleThreeState
    % Binding
    cmap = linspecer(4);
    sampleNumber = 1;
    name = simpleDualDwells(sampleNumber).name{1};
    cutoffFraction = 0.99;
    if plotNonMarkov
        plotOneDoubleExpOneMinus(threeStateDwells(sampleNumber).middleLeftRebindDwells(1, :), cutoffFraction,{fitParams}, name, cmap(1, :),...
            {'ks'}, outsideBoxSpecs, sprintf('Brushless Middle State Dwell (Left Time to Rebind) Cutoff = %0.3f', cutoffFraction))
        plotOneDoubleExpOneMinus(threeStateDwells(sampleNumber).middleLeftRightTransitionDwells(1, :),cutoffFraction, {fitParams}, name, cmap(1, :),...
            {'ks'}, outsideBoxSpecs, sprintf('Brushless Middle State Dwell (Left Right Transition) Cutoff = %0.3f', cutoffFraction))
        plotOneDoubleExpOneMinus(threeStateDwells(sampleNumber).middleRightLeftTransitionDwells(1, :), cutoffFraction,{fitParams}, name, cmap(1, :),...
            {'ks'}, outsideBoxSpecs, sprintf('Brushless Middle State Dwell (Right Left Transistion) Cutoff = %0.3f', cutoffFraction))
         plotOneDoubleExpOneMinus(threeStateDwells(sampleNumber).middleRightRebindDwells(1, :),cutoffFraction, {fitParams}, name, cmap(1, :),...
            {'ks'}, outsideBoxSpecs, sprintf('Brushless Middle State Dwell (Right Time to Rebind) Cutoff = %0.3f', cutoffFraction))
    end
    plotOneDoubleExpOneMinus(threeStateDwells(sampleNumber).middleRightDwells(1, :),cutoffFraction, {fitParams}, name, cmap(1, :),...
            {'ks'}, outsideBoxSpecs, sprintf('Brushless Middle Right Dwell Cutoff = %0.3f', cutoffFraction))
    plotOneDoubleExpOneMinus(threeStateDwells(sampleNumber).middleLeftDwells(1, :),cutoffFraction, {fitParams}, name, cmap(1, :),...
            {'ks'}, outsideBoxSpecs, sprintf('Brushless Middle Left Dwell Cutoff = %0.3f', cutoffFraction))
        
end

%% Plotting Functions

function compileTwoStateFrames()
    highStates = [];
    lowStates = [];
    highRawStates = [];
    lowRawStates = [];
    for i = 1:size(viterbiMean, 2)
        highStates = cat(1, highStates, viterbiMean{i}(viterbiMean{i}==max(viterbiMean{i})));
        lowStates = cat(1, lowStates, viterbiMean{i}(viterbiMean{i}==min(viterbiMean{i})));
        highRawStates = cat(1, highRawStates, FRET{i}(viterbiMean{i}==max(viterbiMean{i})));
        lowRawStates = cat(1, lowRawStates, FRET{i}(viterbiMean{i}==min(viterbiMean{i})));
    end
end

function plotAndFitOneMinus(lowDwell, highDwell, lowFitParams, highFitParams, title)
    Fig = figure('Position', [100 100 1000 700]);
    plotCumDistsOneMinus(lowDwell', Fig, 'ks', 'Dissociation');
    oneMinusLowParams = plotSingleFitsOneMinus(Fig, lowDwell', lowFitParams, 'r', title, 'Fit (Dissociation)');
    plotCumDistsOneMinus(highDwell', Fig, 'ko', 'Binding');
    oneMinusHighParams = plotFitsOneMinus(Fig, highDwell', highFitParams, 'b', title, 'Fit (Binding)');
    str = {sprintf('Dissociation: k = %.2f',oneMinusLowParams(1)), ...
        sprintf('Binding: A = %.2f, k1 = %.2f, k2 = %.2f',oneMinusHighParams(1), oneMinusHighParams(2), oneMinusHighParams(3))};
    t = annotation('textbox',[0.4 0.12 0.5 0.12],'String',str);%,'FitBoxToText','on');
    t.FontSize = 16;
    t.FontWeight = 'bold';
    xlim([0 45])
    ylim([5e-3 1])
    saveas(Fig, strcat(title, ".png"))
end

function plotNDoubleExpOneMinus(dwells, params, names, colors, markers, boxSpecs, title)
    Fig = figure('Position', [100 100 1500 800]);
    N = length(names);
    str = cell(N);
    for i = 1:N
        plotCumDistsOneMinus(dwells{i}', Fig, markers{i}, names{i});
        oneMinusFirstParams = plotFitsOneMinus(Fig, dwells{i}', params{i}, colors(i, :), title, names{i});
        str{i} = strcat(names{i}, sprintf(': A = %.2f, k1 = %.2f, k2 = %.2f',oneMinusFirstParams(1), oneMinusFirstParams(2), oneMinusFirstParams(3)));
    end
    t = annotation('textbox',boxSpecs,'String',str, 'Interpreter','none');%,'FitBoxToText','on');
    t.FontSize = 18;
    t.FontWeight = 'bold';
    xlim([0 20])
    ylim([4e-4 1])
    set(gca, 'LineWidth', 3, 'FontSize', 22, 'FontWeight', 'bold')
    set(gca,'OuterPosition', [0 0 0.88 1])
    saveas(Fig, strcat(title, ".png"))
end

function plotOneSingleExpOneMinusCamera(dwells, cutoffFraction, params, names, colors, markers, boxSpecs, title)
    Fig = figure('Position', [100 100 1500 800]);
    N = length(names);
    str = cell(N);
    plotCumDistsOneMinusCutoff(dwells', cutoffFraction, Fig, markers{1}, names);
    oneMinusFirstParams = plotFitsCutoffOneMinus(Fig, dwells', cutoffFraction, params{1}, colors(1, :), title, names);
    % oneMinusFirstParams = plotFitsOneMinus(Fig, dwells', params{1}, colors(1, :), title, names);
    str{1} = strcat(names, sprintf(': A = %.2f, k1 = %.2f, k2 = %.2f',oneMinusFirstParams(1), oneMinusFirstParams(2), oneMinusFirstParams(3)));
    str{2} = sprintf('N = %d Dwell Times',length(dwells));
    str{3} = sprintf('Cutoff Fraction = %0.3f',cutoffFraction);
    t = annotation('textbox',boxSpecs,'String',str, 'Interpreter','none');%,'FitBoxToText','on');
    t.FontSize = 18;
    t.FontWeight = 'bold';
    xlim([0 10])
    ylim([1e-3 1])
    set(gca, 'LineWidth', 3, 'FontSize', 22, 'FontWeight', 'bold')
    set(gca,'OuterPosition', [0 0 0.88 1])
    saveas(Fig, strcat(title, ".png"))
end

function plotOneDoubleExpOneMinus(dwells, cutoffFraction, params, names, colors, markers, boxSpecs, title)
    Fig = figure('Position', [100 100 1500 800]);
    N = length(names);
    str = cell(N);
    plotCumDistsOneMinusCutoff(dwells', cutoffFraction, Fig, markers{1}, names);
    oneMinusFirstParams = plotFitsCutoffOneMinus(Fig, dwells', cutoffFraction, params{1}, colors(1, :), title, names);
    % oneMinusFirstParams = plotFitsOneMinus(Fig, dwells', params{1}, colors(1, :), title, names);
    str{1} = strcat(names, sprintf(': A = %.2f, k1 = %.2f, k2 = %.2f',oneMinusFirstParams(1), oneMinusFirstParams(2), oneMinusFirstParams(3)));
    str{2} = sprintf('N = %d Dwell Times',length(dwells));
    str{3} = sprintf('Cutoff Fraction = %0.3f',cutoffFraction);
    t = annotation('textbox',boxSpecs,'String',str, 'Interpreter','none');%,'FitBoxToText','on');
    t.FontSize = 18;
    t.FontWeight = 'bold';
    xlim([0 10])
    ylim([1e-3 1])
    set(gca, 'LineWidth', 3, 'FontSize', 22, 'FontWeight', 'bold')
    set(gca,'OuterPosition', [0 0 0.88 1])
    saveas(Fig, strcat(title, ".png"))
end

function plotNSingleExpOneMinus(dwells, params, names, colors, markers, boxSpecs, title)
    Fig = figure('Position', [100 100 1200 800]);
    N = length(names);
    str = cell(N);
    for i = 1:N
        plotCumDistsOneMinus(dwells{i}', Fig, markers{i}, names{i});
        oneMinusFirstParams = plotSingleFitsOneMinus(Fig, dwells{i}', params{i}, colors(i, :), title, names{i});
        str{i} = strcat(names{i}, sprintf(': k = %.2f',oneMinusFirstParams));
    end
    t = annotation('textbox',boxSpecs,'String',str, 'Interpreter','none');%,'FitBoxToText','on');
    t.FontSize = 18;
    t.FontWeight = 'bold';
    xlim([0 30])
    ylim([9e-3 1])
    set(gca, 'LineWidth', 3, 'FontSize', 22, 'FontWeight', 'bold')
    set(gca,'OuterPosition', [0 0 0.88 1])
    saveas(Fig, strcat(title, ".png"))
end

function params = plotFitsOneMinus(fig, dwells, fitParams, clr, Title, fitName)
    fig = figure(fig);
    cutoffFraction = 0.95; % Only Fit this percentile of the fasest dwell times
    [CumDist, CumDistTimes] = ecdf(dwells);
    CumDistTimes(1) = 0;
    dwellTimeCutoff = min(CumDistTimes(CumDist>=cutoffFraction));
    sortDwells = sort(dwells);
    params = plotDoubleExpOneMinus(dwells(dwells<dwellTimeCutoff), fitParams, clr, fitName);
    lgd = legend('Location', 'bestoutside', 'Interpreter','none');
    title(lgd, 'Brush Lengths')
    title(Title, 'Interpreter','none')
    xlabel('Time (s)')
    ylabel('CCDF')
    set(gca, 'linewidth', 2, 'fontweight','bold', 'fontsize',16)
    % params;
end

function params = plotFitsCutoffOneMinus(fig, dwells, cutoffFraction, fitParams, clr, Title, fitName)
    fig = figure(fig);
    [CumDist, CumDistTimes] = ecdf(dwells);
    CumDistTimes(1) = 0;
    dwellTimeCutoff = min(CumDistTimes(CumDist>=cutoffFraction));
    sortDwells = sort(dwells);
    params = plotDoubleExpOneMinus(dwells(dwells<dwellTimeCutoff), fitParams, clr, fitName);
    lgd = legend('Location', 'bestoutside', 'Interpreter','none');
    title(lgd, 'Brush Lengths')
    title(Title, 'Interpreter','none')
    xlabel('Time (s)')
    ylabel('CCDF')
    set(gca, 'linewidth', 2, 'fontweight','bold', 'fontsize',16)
    % params;
end

function params = plotFitsCutoffOneMinusSingleCamera(fig, dwells, cutoffFraction, fitParams, clr, Title, fitName)
    fig = figure(fig);
    [CumDist, CumDistTimes] = ecdf(dwells);
    CumDistTimes(1) = 0;
    dwellTimeCutoff = min(CumDistTimes(CumDist>=cutoffFraction));
    % sortDwells = sort(dwells);
    params = plotDoubleExpOneMinus(dwells(dwells<dwellTimeCutoff), fitParams, clr, fitName);
    lgd = legend('Location', 'bestoutside', 'Interpreter','none');
    title(lgd, 'Brush Lengths')
    title(Title, 'Interpreter','none')
    xlabel('Time (s)')
    ylabel('CCDF')
    set(gca, 'linewidth', 2, 'fontweight','bold', 'fontsize',16)
    % params;
end

function mu = plotSingleFitsOneMinus(fig, dwells, fitParams, clr, Title, fitName)
    fig = figure(fig);
    justDwells = dwells(1,:);
    mu = plotSingleExpOneMinus(dwells, fitParams, clr, fitName);
    lgd = legend('Location', 'southwest', 'Interpreter','none');
    title(Title, 'Interpreter','none')
    title(lgd, 'Brush Lengths')
    xlabel('Time (s)')
    ylabel('CCDF')
    set(gca, 'linewidth', 2, 'fontweight','bold', 'fontsize',16)
    % params;
end

function fig = plotCumDists(dwells, fig, mk, displayName)
    justLeftDwells = dwells(1,:);
    % Plot cumulative sum
    [CumDist, CumDistTimes] = ecdf(dwells);
    CumDistTimes(1) = 0;
    if fig == "null"
        fig = figure();
    else
        fig = figure(fig);
    end
    plot(CumDistTimes, CumDist, mk, 'DisplayName',displayName, 'LineWidth',2, 'MarkerSize',7)
    hold on
end

function fig = plotCumDistsOneMinus(dwells, fig, mk, displayName)
    justLeftDwells = dwells(1,:);
    % Plot cumulative sum
    [CumDist, CumDistTimes] = ecdf(dwells);
    CumDistTimes(1) = 0;
    if fig == "null"
        fig = figure();
    else
        fig = figure(fig);
    end
    semilogy(CumDistTimes, 1-CumDist, mk, 'DisplayName',displayName, 'LineWidth',2, 'MarkerSize',10)
    hold on
end

function fig = plotCumDistsOneMinusCutoff(dwells, cutoff, fig, mk, displayName)
    justLeftDwells = dwells(1,:);
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
    semilogy(CumDistTimes, 1-CumDist, mk, 'DisplayName',displayName, 'LineWidth',2, 'MarkerSize',10)
    hold on
end

function [doubleExpParams] = plotDoubleExpOneMinus(dwellTimes, fitParams, clr, fitName)
    [userPDF, dataVar, fitVar, ~,~, ~]=PDFList('Double Exp (Independent)', 'all', 0);
    lb = fitParams{2};
    ub = fitParams{1};
    guess = fitParams{3};
    annealTemp = fitParams{4};
    fractionToFit = 1;
    sortDwell = sort(dwellTimes);
    dwellsToFit = sortDwell(1:round(fractionToFit*length(sortDwell)));
    [doubleExpParams, ~]= MEMLETCL(dwellsToFit, userPDF, dataVar, fitVar, lb,ub, guess,annealTemp); % Fit parameters and Log Likelihood output
    % Ensure first fitted rate is the faster rate to aid interpretation
    if doubleExpParams(2)<doubleExpParams(3)
        rearrangedParams = doubleExpParams;
        rearrangedParams(1) = 1-doubleExpParams(1);
        rearrangedParams(2) = doubleExpParams(3);
        rearrangedParams(3) = doubleExpParams(2);
        doubleExpParams = rearrangedParams;
    end
    fitxvals=linspace(0,max(dwellsToFit),10000)'; %create variables for plotting along x
    % fitteddh=dbexppdfnotmin(fitxvals,doubleExpParams);  
    % fittedCDF=cumtrapz(fitteddh(~isnan(fitteddh))); %take out any NaNs when doing cumulative (maybe at x=0?) 
    % fittedCDFh=fittedCDF/max(fittedCDF); %normalize CDF
    % oneMinus = 1 - fittedCDFh; % allows for log scale linear fitting
    oneMinus = dbexppdfnotminoneminus(fitxvals, doubleExpParams);
    semilogy(fitxvals,oneMinus, 'Color', clr, 'LineWidth', 3 , 'DisplayName',...
        fitName)
    % semilogy(fitxvals(fitxvals<sortDwell(round(0.95*length(dwellTimes)))),oneMinus(fitxvals<sortDwell(round(0.95*length(dwellTimes)))) , 'Color', clr, 'LineWidth', 3 , 'DisplayName',...
    %     fitName)
end

function [singleExpParams] = plotSingleExpOneMinus(dwellTimes, singleFitParams, clr, fitName)
    [userPDF, dataVar, fitVar, ~,~, ~]=PDFList('Single Exp', 'all', 0);
    lb = singleFitParams{2};
    ub = singleFitParams{1};
    guess = singleFitParams{3};
    annealTemp = singleFitParams{4};
    [singleExpParams, ~]= MEMLETCL(dwellTimes, userPDF, dataVar, fitVar, lb,ub, guess,annealTemp); % Fit parameters and Log Likelihood output
    sortDwell = sort(dwellTimes);
    fitxvals=linspace(0,max(dwellTimes),10000)'; %create variables for plotting along x
    fitteddh=exppdfcalcnotmin(fitxvals,singleExpParams);  
    fittedCDF=cumtrapz(fitteddh(~isnan(fitteddh))); %take out any NaNs when doing cumulative (maybe at x=0?) 
    fittedCDFh=fittedCDF/max(fittedCDF); %normalize CDF
    oneMinus = 1 - fittedCDFh; % allows for log scale linear fitting
    % plot(fitxvals(~isnan(fitteddh),oneMinus(1:end-1) , clr, 'LineWidth', 2 , 'DisplayName',...
        % fitName)
    semilogy(fitxvals(fitxvals<sortDwell(round(0.99*length(dwellTimes)))),oneMinus(fitxvals<sortDwell(round(0.99*length(dwellTimes)))) , 'Color', clr, 'LineWidth', 3 , 'DisplayName',...
        fitName)
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

    plot(S.ax2, seconds, rightHMM{1}, '-','LineWidth', 2,'Color',[0 0.2 0], DisplayName="Right Quenching HMM")
    xlabel(S.ax2,'Time (s)')
    ylabel(S.ax2,'Cy3 Emission')
    plot(S.ax2, seconds, rightRaw{1}', '-', 'Color',[0 0.6 0], DisplayName="Right Quenching Raw")

    % Three State HMM PLot
    plot(S.ax1, seconds, dualHMM{1}', '-','LineWidth', 1,'Color',[0 0.2 0], DisplayName="Dual Quenching HMM")
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
  
    S.sl = uicontrol('style','slide',...
                     'unit','pix',...
                     'position',[0 10 1200 30],...
                     'min',1,'max',size(dualHMM,2),'val',1,...
                     'sliderstep',[1/size(dualHMM,2) 1/size(dualHMM,2)],...
                     'callback',{@sl_call,S, leftRaw, leftHMM, rightRaw, rightHMM, dualHMM, secPerFrame});  
end

function [] = slider_plot_oneside(side, raw, HMM, secPerFrame, titlestr)
    % Make different plots according to slider location.
    % top plot, raw side of interest
    % middle plot, just HMM
    % bottom plot, "FRET" (which should be same as side of interest) w/ HMM

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
    
    if strcmpi(side, 'left')
        HMMclr = [0.4 0 0];
        rawclr = [0.8 0 0];
        fluor = 'Cy5';
    elseif strcmpi(side, 'right')
        HMMclr = [0 0.2 0];
        rawclr = [0 0.6 0];
        fluor = 'Cy3';
    else
        disp("Check side variable")
    end

    seconds = 0:secPerFrame:(size(raw{1}, 1)-1)*secPerFrame;
    plot(S.ax3, seconds, raw{1}', '-','LineWidth', 2, 'Color',rawclr, DisplayName=strcat(side, " Quenching raw"))
    xlabel(S.ax3,'Time (s)')
    ylabel(S.ax3,strcat(fluor, ' Emission'))
    title(S.ax3,sprintf("Trace %d", 1))

    plot(S.ax2, seconds, HMM{1}, '-','LineWidth', 2,'Color',HMMclr, DisplayName=strcat(side, " Quenching HMM"))
    xlabel(S.ax2,'Time (s)')
    ylabel(S.ax2,strcat(fluor, ' Emission'))
    
    % Three State HMM PLot
    plot(S.ax1, seconds, raw{1}', '-','LineWidth', 2, 'Color',rawclr, DisplayName=strcat(side, " Quenching raw"))
    xlabel(S.ax1,'Time (s)')
    ylabel(S.ax1,strcat(fluor, ' Emission'))
    
    plot(S.ax1, seconds, HMM{1}, '-','LineWidth', 2,'Color',HMMclr, DisplayName=strcat(side, " Quenching HMM"))
    xlabel(S.ax1,'Time (s)')
    ylabel(S.ax1,strcat(fluor, ' Emission'))

    for a =[S.ax1 S.ax2 S.ax3]
        legend(a)
    end
  
    S.sl2 = uicontrol('style','slide',...
                     'unit','pix',...
                     'position',[0 10 1200 30],...
                     'min',1,'max',size(HMM,2),'val',1,...
                     'sliderstep',[1/size(HMM,2) 1/size(HMM,2)],...
                     'callback',{@sl2_call, S, raw, HMM, secPerFrame, HMMclr, rawclr, side, fluor});  
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

    plot(S.ax2, seconds, rightHMM{traceNum}', '-','LineWidth', 2,'Color',[0 0.2 0], DisplayName="Right Quenching HMM")
    xlabel(S.ax2,'Time (s)')
    ylabel(S.ax2,'Cy3 Emission')
    plot(S.ax2, seconds, rightRaw{traceNum}', '-', 'Color',[0 0.6 0], DisplayName="Right Quenching Raw")

    % Three State HMM PLot
    plot(S.ax1, seconds, duaHMM{traceNum}', '-','LineWidth', 1,'Color',[0 0.2 0], DisplayName="Dual Quenching HMM")
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

function [] = sl2_call(varargin)
% Callback for the slider.
    [h,S] = varargin{[1,3]};  % calling handle and data structure.
    raw = varargin{4};
    HMM = varargin{5};
    secPerFrame = varargin{6};
    HMMclr = varargin{7};
    rawclr = varargin{8};
    side = varargin{9};
    fluor = varargin{10};
    % cla
    traceNum = round(get(h,'value'));
    
    for a =[S.ax1 S.ax2 S.ax3]
        cla(a)
        hold(a, 'on')
    end
    
    seconds = 0:secPerFrame:(size(raw{traceNum}, 1)-1)*secPerFrame;
    plot(S.ax3, seconds, raw{traceNum}', '-','LineWidth', 2, 'Color',rawclr, DisplayName=strcat(side, " Quenching raw"))
    xlabel(S.ax3,'Time (s)')
    ylabel(S.ax3,strcat(fluor, ' Emission'))
    title(S.ax3,sprintf("Trace %d", traceNum))

    plot(S.ax2, seconds, HMM{traceNum}, '-','LineWidth', 2,'Color',HMMclr, DisplayName=strcat(side, " Quenching HMM"))
    xlabel(S.ax2,'Time (s)')
    ylabel(S.ax2,strcat(fluor, ' Emission'))
    
    % Three State HMM PLot
    plot(S.ax1, seconds, raw{traceNum}', '-','LineWidth', 2, 'Color',rawclr, DisplayName=strcat(side, " Quenching raw"))
    xlabel(S.ax1,'Time (s)')
    ylabel(S.ax1,strcat(fluor, ' Emission'))
    
    plot(S.ax1, seconds, HMM{traceNum}, '-','LineWidth', 2,'Color',HMMclr, DisplayName=strcat(side, " Quenching HMM"))
    xlabel(S.ax1,'Time (s)')
    ylabel(S.ax1,strcat(fluor, ' Emission'))

    for a =[S.ax1 S.ax2 S.ax3]
        legend(a)
    end

end