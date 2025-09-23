%---------------------------------------------------
%  NAME:      Post Process Excel Export.m
%  WHAT:      Read machine parameters and export them into an excel file
%
%  AUTHOR:    Ozan Keysan (09/2025)


%Main Parameters

main_parameters = {'Parameter', 'Value', 'Unit';
                   'Power Output', round(machine.P_output/1e6), 'MW';
                   'Speed',  machine.Nrpm, 'rpm';
                   'Mean Diameter', 2*stator.R_mean, 'm';
                   'Number of Poles', machine.Npole, '';
                   'Number of Stator Stacks', machine.Nstacks, '';
                   'Efficiency', round(100*machine.efficiency,1), '%';
                   'Total length of HTS Tape', round(HTS.length_total/1e3,1), 'km';
                   'Total copper mass', round(stator.mass/1e3,1), 't';
                   'Max. Airgap Flux Density', round(stator.B_max,2), 'T';
                   };
%Create File Name
%filename = 'deneme.xls';
filename = round(machine.P_output/1e6) + "_MW_" + machine.Npole + "_pole_" + 2*stator.R_mean +"_m_diameter_" + machine.Nstacks + "_stack_generator" +".xls";

writecell(main_parameters, filename, 'Sheet', 'Main_Parameters');

%Write Machine Parameters Struct
writetable(rows2vars(struct2table(machine)), filename, 'Sheet', 'Machine')

%Write Stator Parameters Struct
writetable(rows2vars(struct2table(stator)), filename, 'Sheet', 'Stator')

%Write HTS Parameters Struct
writetable(rows2vars(struct2table(HTS)), filename, 'Sheet', 'HTS Rotor')