function plotHHI2017_MW(TrialData, types)
    %PLOTHHI2017: Generates figures from processed HHI data recorded after
    %September 2017
    %   Detailed explanation goes here
    
    % first component of marker data is x or lateral (MW)    
    
    % Add in:
    %   plot force magnitude
    %   plot force alignment
    
    
    % Luke Drnach
    % October 30, 2017
    
    %% Initialization
    % Number of trials
    num_trials = length(TrialData);
    if nargin == 1
        types = {'Assist Ground','Assist Beam','Solo Beam','Assist Solo'};
    end
    % Keep track of which trials are assist ground and which are assist
    % beam
    if ~isfield(TrialData(1), 'Results')
        fprintf('\nResults have not been computed yet\n')
        return;
    end
    for n = 1:num_trials
        if any(strcmpi(TrialData(n).Info.Condition, types))
            if isempty(TrialData(n).Results)
                fprintf('\nResults for Trial %d have not been computed yet\n', n);
                continue;
            end
            if any(strcmpi(TrialData(n).Info.Condition, {'Assist Ground','Assist Beam','Solo Beam'}))
                % Unpack Time and Beamer Sway
                time = TrialData(n).Results.time;
                sway = TrialData(n).Results.beamerSway;
            end
            if any(strcmpi(TrialData(n).Info.Condition, {'Assist Ground','Assist Beam'}))
                % Unpack the time, force, arm, work, and power
                force = TrialData(n).Results.Forces;
                arm = TrialData(n).Results.AssistArm;
                work = TrialData(n).Results.AssistWork;
                power = TrialData(n).Results.AssistPower;
                align = TrialData(n).Results.PowerAlignment;
                % Now make individual plots
                label = [TrialData(n).Info.Trial,': ',TrialData(n).Info.Condition];
                figure('Name', label)
                subplot(3,1,1)
                plot(time, sway,'LineWidth',1.5)
                hold on
                plot(time, zeros(1, length(time)), 'k--','LineWidth',1.5)
                ylabel({'PoB Lateral','Sway (mm)'})
                title(label)
                subplot(3,1,2)
                plot(time, force(3,:),'LineWidth',1.5)
                ylabel('Vertical Force (N)')
                subplot(3,1,3)
                hold on
                plot(time, force(1,:),'LineWidth',1.5)
                plot(time, force(2,:),'LineWidth',1.5)
                legend('X','Y')
                legend boxoff
                ylabel({'Horizontal','Forces (N)'})
                xlabel('Time (s)')
                set(get(gcf,'Children'), 'FontName','Garamond','FontSize',12,'FontWeight','bold')
                figure('Name',label)
                subplot(3,1,1)
                plot(time, sqrt(sum(force.^2, 1)),'LineWidth',1.5)
                ylabel({'Net Force','Magnitude (N)'})
                title(label)
                subplot(3,1,2)
                ArmL = diff(arm,1,2)*TrialData(n).Markers.samplerate;
                ArmL = sqrt(sum(ArmL.^2,1));
                plot(time, [0, ArmL],'LineWidth',1.5)
                ylabel({'Arm Length','Velocity (mm/s)'})
                subplot(3,1,3)
                plot(time,[0, align],'LineWidth',1.5)
                hold on
                plot(time, zeros(1, length(time)), 'k--','LineWidth',1.5)
                ylabel('cos(\theta)')
                xlabel('Time (s)')
                set(get(gcf,'Children'), 'FontName','Garamond','FontSize',12,'FontWeight','bold')
                figure()
                subplot(3,1,1)
                plot(time, power,'LineWidth',1.5)
                hold on
                plot(time, zeros(1, length(time)), 'k--','LineWidth',1.5)
                ylabel('Power (W)')
                title(label)
                subplot(3,1,2)
                plot(time, work,'LineWidth',1.5)
                hold on
                plot(time, zeros(1, length(time)), 'k--','LineWidth',1.5)
                ylabel({'Cumulative','Work (J)'})
                xlabel('Time (s)')
                set(get(gcf,'Children'), 'FontName','Garamond','FontSize',12,'FontWeight','bold')
            end
        end
    end
end

