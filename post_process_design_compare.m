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
%save("optimized_results", "design_table", "design_inputs")

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

title(t,'20 MW, 10 rpm, 3 stacks, 6 m diameter, 95% Efficiency with changing Number of Poles', ...
    'FontSize', 12)
xlabel(t,'Number of Poles', 'FontSize', 12)

%HTS Length
nexttile
plot([design_table.machine.Npole], [design_table.HTS.length_total]/1e3, '-*')
ylabel('Total Length of HTS Tape (km)', 'FontSize', 12)
ylim([0 500])
grid on
box on

%Stator Mass
nexttile
plot([design_table.machine.Npole], [design_table.stator.mass]/1e3, '-o')
ylabel('Total Mass of Stator (t)', 'FontSize', 12)
ylim([0 20])
grid on
box on

% %Efficiency
% nexttile
% plot([design_table.machine.Npole], [design_table.machine.efficiency]*100, '-o')
% ylabel('Efficiency (%)')
% ylim([90 100])
% grid on

%Stator Current Density
nexttile
plot([design_table.machine.Npole], [design_table.stator.current_density], '-o')
ylabel('Stator Current Density (A/mm2)', 'FontSize', 12)
ylim([0 8])
grid on
box on

%Airgap Flux Density
nexttile
plot([design_table.machine.Npole], [design_table.stator.B_max], '-o')
ylabel('Airgap Flux Density (T), top 10%', 'FontSize', 12)
ylim([0 2.5])
grid on
box on

