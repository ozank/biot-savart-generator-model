%---------------------------------------------------
%  NAME:      Post Process Report Generatiom.m
%  WHAT:      Read from the optimization solutions to run the main function
% to calculate main machine parameters and then Create a report that
% contains important parameters of the optimised design and export as a
% report document (PDF, Word or Html)
%
%  AUTHOR:    Ozan Keysan (09/2025)
%----------------------------------------------------

%Read solution variable from optimization outputs
%Input arrangements
%Make sure they are aligned with the actual optimization cost function

machine.Npole = 4*solution(1);   %Number of Poles divided by 4  (make sure it is a multiple of 4
%machine.Npole = 4*floor(inputs(1)/4);   %Number of Poles (make sure it is a multiple of 4
stator.current_density = solution(2);     % Current Density (A/mm^2)
HTS.coil_length = solution(3);
HTS.N_turns = solution(4);
stator.N_turns = solution(5);
stator.coil_width_to_coil_pitch_ratio = solution(6);
stator.coil_thickness = solution(7);
machine.Nstacks = solution(8);


%% MAIN Function Call

main;

%% Report Genertation Functions
% https://uk.mathworks.com/help/rptgen/ug/create-a-report-generator.html

import mlreportgen.report.* 
import mlreportgen.dom.* 

%Define report output file type
% https://uk.mathworks.com/help/rptgen/ug/mlreportgen.report.report-class.html
rpt = Report("optimization_output","pdf"); 

%Title Page
tp = TitlePage; 
%tp.Title = "Design Report of Optimised Generator"; 
tp.Subtitle = append(num2str(round(machine.P_output/1e6)), " MW ", num2str(round(machine.Nrpm)), " rpm ", "Generator");
tp.Author = "Ozan Keysan"; 
append(rpt,tp); 



%Main Parameters

tableData = cell2table(main_parameters) 
% Generating Table into the report 
tObj = Table(tableData); 
append(rpt,tObj); 


%Appendix
%All parametes



close(rpt)
rptview(rpt)

