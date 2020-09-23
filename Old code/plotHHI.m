function plotHHI(TrialData)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    %
    
    %Luke Drnach
    %May 31, 2017
    %Neuromechanics Lab
    %Georgia Tech and Emory
   
    %% Initialization
    NTrials = length(TrialData);
    % Keep track of the maximum compressive, tensile, and shear forces and total power/work delivered/absorbed
    maxPowAbs = zeros(1,NTrials);   %Maximum power absorbed
    maxPowDel = zeros(1,NTrials);   %Maximum power delivered
    maxWork = zeros(1,NTrials);     %Maximum Work done
    maxCForce = zeros(1,NTrials);   %Maximum Compressive Force (Vertical)
    maxTForce = zeros(1,NTrials);   %Maximum Tensile Force (Vertical)
    maxCXForce = zeros(1,NTrials);  %Maximum Compressive X Force
    maxTXForce = zeros(1,NTrials);  %Maximum Tensile X Force
    % Keep track of which trials are assist Ground and which are assist
    % beam
    assistBeam = false(1,NTrials);
    assistGround = false(1,NTrials);
    % Check if Results have been processed yet
    if ~isfield(TrialData(1),'Results')
        fprintf('\nResults have not been computed yet\n');
        return;
    end
    %% Main Loop for plotting individual trials
    for n = 1:NTrials
        %First, check if this is an assist ground or assist beam trial, and
        %if the Results exist
        if any(strcmp(TrialData(n).Info.Condition,{'Assist Beam','Assist Ground'})) && ~isempty(TrialData(n).Results)
            if strcmp(TrialData(n).Info.Condition,'Assist Beam')
                assistBeam(n) = true;
            else
                assistGround(n) = true;
            end
            %Now, get the time, force, arm, and work vectors from the Results
            time = TrialData(n).Results.time;
            force = TrialData(n).Results.Force;
            BeamSway = TrialData(n).Results.Clav;
            AssistArm = TrialData(n).Results.AssistArm;
            AssistWork = TrialData(n).Results.AssistWork;
            AssistPow = TrialData(n).Results.AssistPow;
            % Get the aggregate data
            p = max(AssistPow(AssistPow < 0));
            if isempty(p)
                maxPowAbs(n) = NaN;
            else
                maxPowAbs(n) = abs(p);
            end
            p = max(AssistPow(AssistPow > 0));
            if isempty(p)
                maxPowDel(n) = NaN;
            else
                maxPowDel(n) = p;
            end
            p = max(AssistWork);
            if isempty(p)
                maxWork(n) = NaN;
            else
                maxWork(n) = p;
            end
            p = max(force(3,force(3,:)<0));
            if isempty(p)
                maxCForce(n) = NaN;
            else
                maxCForce(n) = abs(p);
            end
            p = max(force(3,force(3,:)>0));
            if isempty(p)
                maxTForce(n) = NaN;
            else
                maxTForce(n) = p;
            end
            p = max(force(1,force(1,:)<0));
            if isempty(p)
                maxCXForce(n) = NaN;
            else
                maxCXForce(n) = abs(p);
            end
            p = max(force(1,force(1,:)>0));
            if isempty(p)
                maxTXForce(n) = NaN;
            else
                maxTXForce(n) = p;
            end
            % Open the figure and make the plots
            label = [TrialData(n).Info.Trial,': ',TrialData(n).Info.Condition];
            figure('Name',label);
            subplot(6,1,1);
            plot(time,BeamSway);
            ylabel('Beamer Lateral Sway (mm)')
            subplot(6,1,2)
            plot(time,force(3,:));
            ylabel('Vertical Force (N)')
            subplot(6,1,3)
            hold on;
            plot(time,force(1,:))
            plot(time,force(2,:))
            legend('X','Y')
            legend boxoff
            ylabel('Horizontal Forces (N)')
            subplot(6,1,4)
            ArmL = sqrt(AssistArm(1,:).^2 + AssistArm(2,:).^2 + AssistArm(3,:).^2);
            ArmL = ArmL - ArmL(1);
            plot(time,ArmL);
            ylabel('Arm Length Change (mm)')
            subplot(6,1,5);
            plot(time,AssistPow)
            ylabel('Assistive Power (W)')
            subplot(6,1,6)
            plot(time,AssistWork);
            ylabel('Cumulative Assistive Work (J)')
            set(get(gcf,'Children'),'FontName','Garamond','FontSize',12,'FontWeight','bold');
        end
    end
    %% Plots of aggregate data
    figure('Name','Specifications');
    %Compressive Vertical Force
    subplot(1,7,1);
    meanCFA = mean(maxCForce(assistGround),'omitnan');
    p = sum(~isnan(maxCForce(assistGround)));
    if p~= 0
        stdCFA = std(maxCForce(assistGround),'omitnan')./sqrt(p);
    else
        stdCFA = NaN;
    end
    meanCFB = mean(maxCForce(assistBeam),'omitnan');
    p = sum(~isnan(maxCForce(assistBeam)));
    if p~=0
        stdCFB = std(maxCForce(assistBeam),'omitnan')./sqrt(p);
    else
        stdCFB = NaN;
    end
    fprintf('\nMax Compressive Vertical Force')
    fprintf('\nAssist Ground: %f%s%f (N)',meanCFA,char(177),stdCFA);
    fprintf('\nAssist Beam: %f%s%f (N)',meanCFB,char(177),stdCFB);
    bar(0.5,meanCFA,'c')
    hold on;
    bar(1,meanCFB,'b')
    legend('Ground','Beam');
    legend('boxoff');
    set(get(gca,'Children'),'BarWidth',0.5);
    set(gca,'box','off','XLim',[0,1.5],'XTickLabel',{''});
    errorbar([0.5,1],[meanCFA,meanCFB],[stdCFA,stdCFB],'k.');
    ylabel('Max Compressive Vertical Force (N)')
    %Tensile Vertical Force
    subplot(1,7,2)
    meanTFA = mean(maxTForce(assistGround),'omitnan');
    meanTFB = mean(maxTForce(assistBeam),'omitnan');
    p = sum(~isnan(maxTForce(assistGround)));
    if p~=0
        stdTFA = std(maxTForce(assistGround),'omitnan')./sqrt(p);
    else
        stdTFA = NaN;
    end
    p = sum(~isnan(maxTForce(assistBeam)));
    if p~=0
        stdTFB = std(maxTForce(assistBeam),'omitnan')./sqrt(p);
    else
        stdTFB = NaN;
    end
    fprintf('\nMax Tensile Vertical Force')
    fprintf('\nAssist Ground: %f%s%f (N)',meanTFA,char(177),stdTFA);
    fprintf('\nAssist Beam: %f%s%f (N)',meanTFB,char(177),stdTFB);
    bar(0.5,meanTFA,'c');
    hold on;
    bar(1,meanTFB,'b');
    set(get(gca,'Children'),'BarWidth',0.5);
    set(gca,'box','off','XLim',[0,1.5],'XTickLabel',{''});
    errorbar([0.5,1],[meanTFA,meanTFB],[stdTFA,stdTFB],'k.');
    ylabel('Max Tensile Vertical Force (N)')
    %Compressive Lat Force
    subplot(1,7,3);
    meanCFA = mean(maxCXForce(assistGround),'omitnan');
    p = sum(~isnan(maxCXForce(assistGround)));
    if p>0
        stdCFA = std(maxCXForce(assistGround),'omitnan')./sqrt(p);
    else
        stdCFA = NaN;
    end
    meanCFB = mean(maxCXForce(assistBeam),'omitnan');
    p = sum(~isnan(maxCXForce(assistBeam)));
    if p>0
        stdCFB = std(maxCXForce(assistBeam),'omitnan')./sqrt(p);
    else
        stdCFB = NaN;
    end
    fprintf('\nMax Compressive Lateral Force')
    fprintf('\nAssist Ground: %f%s%f (N)',meanCFA,char(177),stdCFA);
    fprintf('\nAssist Beam: %f%s%f (N)',meanCFB,char(177),stdCFB);
    bar(0.5,meanCFA,'c')
    hold on;
    bar(1,meanCFB,'b')
    set(get(gca,'Children'),'BarWidth',0.5);
    set(gca,'box','off','XLim',[0,1.5],'XTickLabel',{''});
    errorbar([0.5,1],[meanCFA,meanCFB],[stdCFA,stdCFB],'k.');
    ylabel('Max Compressive Lateral Force (N)')
    %Tensile Vertical Force
    subplot(1,7,4)
    meanTFA = mean(maxTXForce(assistGround),'omitnan');
    meanTFB = mean(maxTXForce(assistBeam),'omitnan');
    p = sum(~isnan(maxTXForce(assistGround)));
    if p>0
        stdTFA = std(maxTXForce(assistGround),'omitnan')./sqrt(p);
    else
        stdTFA = NaN;
    end
    p = sum(~isnan(maxTXForce(assistGround)));
    if p>0
        stdTFB = std(maxTXForce(assistBeam),'omitnan')./sqrt(p);
    else
        stdTFB = NaN;
    end
    
    fprintf('\nMax Tensile Lateral Force')
    fprintf('\nAssist Ground: %f%s%f (N)',meanTFA,char(177),stdTFA);
    fprintf('\nAssist Beam: %f%s%f (N)',meanTFB,char(177),stdTFB);
    bar(0.5,meanTFA,'c');
    hold on;
    bar(1,meanTFB,'b');
    set(get(gca,'Children'),'BarWidth',0.5);
    set(gca,'box','off','XLim',[0,1.5],'XTickLabel',{''});
    errorbar([0.5,1],[meanTFA,meanTFB],[stdTFA,stdTFB],'k.');
    ylabel('Max Tensile Lateral Force (N)')
    % Power Absorbed
    subplot(1,7,5);
    meanPowAbsA = mean(maxPowAbs(assistGround),'omitnan');
    p = sum(~isnan(maxPowAbs(assistGround)));
    if p>0
        stdPowAbsA = std(maxPowAbs(assistGround),'omitnan')./sqrt(p);
    else
        stdPowAbsA = NaN;
    end
    
    meanPowAbsB = mean(maxPowAbs(assistBeam),'omitnan');
    p = sum(~isnan(maxPowAbs(assistBeam)));
    if p>0
        stdPowAbsB = mean(maxPowAbs(assistBeam),'omitnan')./sqrt(p);
    else
        stdPowAbsB = NaN;
    end
    
    fprintf('\nMax Power Absorbed')
    fprintf('\nAssist Ground: %f%s%f (W)',meanPowAbsA,char(177),stdPowAbsA);
    fprintf('\nAssist Beam: %f%s%f (W)',meanPowAbsB,char(177),stdPowAbsB);
    bar(0.5,meanPowAbsA,'c');
    hold on;
    bar(1,meanPowAbsB,'b');
    set(get(gca,'Children'),'BarWidth',0.5);
    set(gca,'box','off','XLim',[0,1.5],'XTickLabel',{''});
    errorbar([0.5,1],[meanPowAbsA,meanPowAbsB],[stdPowAbsA,stdPowAbsB],'k.');
    ylabel('Max Power Absorbed (W)')
    % Power Delivered
    subplot(1,7,6)
    meanPowDelA = mean(maxPowDel(assistGround),'omitnan');
    p = sum(~isnan(maxPowDel(assistGround)));
    if p > 0
        stdPowDelA = std(maxPowDel(assistGround),'omitnan')./sqrt(p);
    else
        stdPowDelA = NaN;
    end
    meanPowDelB = mean(maxPowDel(assistBeam),'omitnan');
    p = sum(~isnan(maxPowDel(assistBeam)));
    if p>0
        stdPowDelB = mean(maxPowDel(assistBeam),'omitnan')./sqrt(p);
    else
        stdPowDelB = NaN;
    end
    fprintf('\nMax Power Delivered')
    fprintf('\nAssist Ground: %f%s%f (W)',meanPowDelA,char(177),stdPowDelA);
    fprintf('\nAssist Beam: %f%s%f (W)',meanPowDelB,char(177),stdPowDelB);
    bar(0.5,meanPowDelA,'c');
    hold on;
    bar(1,meanPowDelB,'b');
    set(get(gca,'Children'),'BarWidth',0.5);
    set(gca,'box','off','XLim',[0,1.5],'XTickLabel',{''});
    errorbar([0.5,1],[meanPowDelA,meanPowDelB],[stdPowDelA,stdPowDelB],'k.');
    ylabel('Max Power Delivered (W)')
    % Work Performed
    subplot(1,7,7)
    meanWorkA = mean(maxWork(assistGround),'omitnan');
    p = sum(~isnan(maxWork(assistGround)));
    if p>0
        stdWorkA = std(maxWork(assistGround),'omitnan')./sqrt(p);
    else
        stdWorkA = NaN;
    end
    meanWorkB = mean(maxWork(assistBeam),'omitnan');
    p = sum(~isnan(maxWork(assistBeam)));
    if p>0
        stdWorkB = mean(maxWork(assistBeam),'omitnan')./sqrt(p);
    else
        stdWorkB = NaN;
    end
    fprintf('\nMax Work Done')
    fprintf('\nAssist Ground: %f%s%f (J)',meanWorkA,char(177),stdWorkA);
    fprintf('\nAssist Beam: %f%s%f (J)',meanWorkB,char(177),stdWorkB);
    bar(0.5,meanWorkA,'c');
    hold on;
    bar(1,meanWorkB,'b');
    set(get(gca,'Children'),'BarWidth',0.5);
    set(gca,'box','off','XLim',[0,1.5],'XTickLabel',{''});
    errorbar([0.5,1],[meanWorkA,meanWorkB],[stdWorkA,stdWorkB],'k.');
    ylabel('Max Work (J)')
    fprintf('\n');
    
    set(get(gcf,'Children'),'FontSize',14,'FontWeight','bold');
    
    
end

