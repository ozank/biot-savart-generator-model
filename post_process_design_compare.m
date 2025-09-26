%---------------------------------------------------
%  NAME:      Post Process Design Compare.m
%  WHAT:      Takes the outputs of several optimization results and calculates their parameters to compare and plot different results
%  AUTHOR:    Ozan Keysan (09/2025)
%----------------------------------------------------



%Requires a design_solution variable to be previously saved containing
%results of optimization,
%one row per design, with the same order of th GA output "solution"
%variable

%Get the size of the design inputs (to be used in for loops)
design_count = size(design_inputs,1); 

% Run each design post process for all the design parameters

for row = design_count:-1:1  %for loop runs backwards for struct allocation
 %https://uk.mathworks.com/matlabcentral/answers/1447869-array-of-structs-how-to-preallocate-and-assign   
     row
%Assign each row to solution variable (compatible with othe post process
%functions)
solution = design_inputs(row,:);

%Run post process tool for each design
post_process_optimization;

%Assign output parameters to a design struct
design_outputs(row).design_input = solution;
design_outputs(row).machine = machine;
design_outputs(row).HTS = HTS;
design_outputs(row).stator = stator;
design_outputs(row).material = material;

end


%Convert to table for easy operation and plotting functions
design_table = struct2table(design_outputs);

% save design tables as mat files for reviewing later on
save("optimized_results", "design_table", "design_inputs")

%% Plots
% Generate Comparison plots in tiled layout

% Create plots
t = tiledlayout(2,2);

%Plot Properties
t.TileSpacing = 'compact';
t.Padding = 'compact';

% Global aesthetics
set(groot, 'defaultAxesFontSize', 12, ...
           'defaultLineLineWidth', 1.5, ...
           'defaultLineMarkerSize', 5);

title(t,'20 MW, 10 rpm, 95% Efficient, Variable Pole Number', ...
    'FontSize', 12)
xlabel(t,'Number of Poles', 'FontSize', 12)

%HTS Length
nexttile
pole = [design_table.machine.Npole];
hts_length = [design_table.HTS.length_total]/1e3;

plot(pole(1:6), hts_length(1:6), ':o')
hold on
plot(pole(7:14), hts_length(7:14), ':^r')
hold on
plot(pole(15:20), hts_length(15:20), ':*k')
ylabel('Total Length of HTS Tape (km)', 'FontSize', 12)

ylim([0 500])

hold off

legend('D = 6m, N stack = 3','D = 6m, N stack = 5','D = 7m, N stack = 3','Location','southeast','Orientation','vertical')

grid on
box on

%Stator Mass
nexttile

hts_turn = [design_table.HTS.N_turns];

plot(pole(1:6), hts_turn(1:6), ':o')
hold on
plot(pole(7:14), hts_turn(7:14), ':^r')
hold on
plot(pole(15:20), hts_turn(15:20), ':*k')
ylabel('Number of HTS turns per coil', 'FontSize', 12)

ylim([0 200])
grid on

%Stator Current Density
nexttile
current_density =  [design_table.stator.current_density];

plot(pole(1:6), current_density(1:6), ':o')
hold on
plot(pole(7:14), current_density(7:14), ':^r')
hold on
plot(pole(15:20), current_density(15:20), ':*k')

ylabel('Stator Current Density (A/mm2)', 'FontSize', 12)
ylim([0 8])
grid on
box on

%Airgap Flux Density
nexttile

flux_density =  [design_table.stator.B_max];

plot(pole(1:6), flux_density(1:6), ':o')
hold on
plot(pole(7:14), flux_density(7:14), ':^r')
hold on
plot(pole(15:20), flux_density(15:20), ':*k')

%plot([design_table.machine.Npole], [design_table.stator.B_max], '-o')
ylabel('Airgap Flux Density (T), top 10%', 'FontSize', 12)
ylim([0 2.5])
grid on
box on

% Set common xticks for all subplots
%set([ax1 ax2 ax3 ax4], 'XTick', [40 80 120 160]);

ax = findall(gcf,'Type','axes');   % get all axes in the current figure
set(ax,'XTick',[40 80 120 160]);   % apply to all
set(ax,'XLim',[0 165]);   % apply to all
% Adjust figure size for publication
%set(gcf, 'Units', 'centimeters', 'Position', [5 5 18 14]); % width=18cm, height=14cm

%% Extra Figure

% Create plots
t = tiledlayout(3,1);

%Plot Properties
t.TileSpacing = 'compact';
t.Padding = 'compact';

% Global aesthetics
set(groot, 'defaultAxesFontSize', 12, ...
           'defaultLineLineWidth', 1.5, ...
           'defaultLineMarkerSize', 5);

title(t,'20 MW, 10 rpm, 95% Efficient, Variable Pole Number', ...
    'FontSize', 12)
xlabel(t,'Number of Poles', 'FontSize', 12)

%HTS Length
nexttile
pole = [design_table.machine.Npole];

HTS_mass = [design_table.HTS.mass]/1e3;

plot(pole(1:6), HTS_mass(1:6), ':o')
hold on
plot(pole(7:14), HTS_mass(7:14), ':^r')
hold on
plot(pole(15:20), HTS_mass(15:20), ':*k')

ylabel('HTS Tape Mass (t)', 'FontSize', 12)
ylim([0 10])
grid on
box on

hold off

legend('D = 6m, N stack = 3','D = 6m, N stack = 5','D = 7m, N stack = 3','Location','southeast','Orientation','vertical')

grid on
box on

%Stator Mass
nexttile
stator_mass = [design_table.stator.mass]/1e3;

plot(pole(1:6), stator_mass(1:6), ':o')
hold on
plot(pole(7:14), stator_mass(7:14), ':^r')
hold on
plot(pole(15:20), stator_mass(15:20), ':*k')

ylabel('Total Mass of Stator (t)', 'FontSize', 12)
ylim([0 30])
grid on
box on



%Stator Mass
nexttile
total_mass = HTS_mass + stator_mass;

plot(pole(1:6), total_mass(1:6), ':o')
hold on
plot(pole(7:14), total_mass(7:14), ':^r')
hold on
plot(pole(15:20), total_mass(15:20), ':*k')

ylabel('Total Active Material Mass (t)', 'FontSize', 12)
%ylim([0 10])
grid on
box on

% %Stator Current Density
% nexttile
% phase_voltage =  [design_table.stator.induced_voltage_per_phase]/1e3;
% 
% plot(pole(1:6), phase_voltage(1:6), ':o')
% hold on
% plot(pole(7:14), phase_voltage(7:14), ':^r')
% hold on
% plot(pole(15:20), phase_voltage(15:20), ':*k')
% 
% ylabel('Induced Phase Voltage (kV)', 'FontSize', 12)
% ylim([0 22])
% grid on
% box on
% 
% %Airgap Flux Density
% nexttile
% 
% phase_current =  [design_table.stator.phase_current];
% 
% plot(pole(1:6), phase_current(1:6), ':o')
% hold on
% plot(pole(7:14), phase_current(7:14), ':^r')
% hold on
% plot(pole(15:20), phase_current(15:20), ':*k')
% 
% %plot([design_table.machine.Npole], [design_table.stator.B_max], '-o')
% ylabel('Phase Current (A) per stack', 'FontSize', 12)
% %ylim([0 2.5])
% grid on
% box on
% 
% % Set common xticks for all subplots
% %set([ax1 ax2 ax3 ax4], 'XTick', [40 80 120 160]);
% 
% ax = findall(gcf,'Type','axes');   % get all axes in the current figure
% set(ax,'XTick',[40 80 120 160]);   % apply to all
% set(ax,'XLim',[0 165]);   % apply to all
% % Adjust figure size for publication
% %set(gcf, 'Units', 'centimeters', 'Position', [5 5 18 14]); % width=18cm, height=14cm