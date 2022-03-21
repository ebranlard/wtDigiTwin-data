%% Documentation   
% Postprocess linearization file. 

%% Initialization
clear all; close all; %clc;
restoredefaultpath;
addpath(genpath('C:/Work/FAST/matlab-toolbox/'))

%% Script Parameters
BladeLen = 63; % taken as tiprad
TowerLen = 67; % Tuned TowerHt -TowerBsHt

%% Derived parameters
FileName = 'Main_Spar_ED.1.lin';
FileNames = {FileName};

%% Performing MBC
[mbc_data, matData, FAST_linData] = fx_mbc3( FileNames );
[CampbellData] = campbell_diagram_data(mbc_data, BladeLen, TowerLen);%, xlsFileName);

%% Outputs to screen
CD=CampbellData;
fprintf(    '%40s , %15s , %15s\n','Mode','Frequencies_[Hz]','Damping_Ratio_[-]');
for i = 1:min(length(CD.NaturalFreq_Hz),15)
    Desc=CD.Modes{i}.DescStates(CD.Modes{i}.StateHasMaxAtThisMode);
    DescCat='';
    DescCatED='';
    if length(Desc)==0
        DescCat = '' ;
        DescCatED = 'NoMax -' ;
        Desc = CD.Modes{i}.DescStates(1:5);
    end
    nBD=0;
    for iD=1:length(Desc)
        s=Desc{iD};
        s = strrep(s,'First time derivative of'     ,'d/dt of');
        s = strrep(s,'fore-aft bending mode DOF, m'    ,'FA'     );
        s = strrep(s,'side-to-side bending mode DOF, m','SS'     );
        s = strrep(s,'bending-mode DOF of blade '    ,''     );
        s = strrep(s,' rotational-flexibility DOF, rad','-rot'   );
        s = strrep(s,'rotational displacement in ','rot'   );
        s = strrep(s,'translational displacement in ','trans'   );
        s = strrep(s,', rad','');
        s = strrep(s,', m','');
        s = strrep(s,'finite element node ','N'   );
        s = strrep(s,'cosine','cos'   );
        s = strrep(s,'sine','sin'   );
        s = strrep(s,'collective','coll.');
        s = strrep(s,'Blade','Bld');
        s = strrep(s,'rotZ','TORS-ROT');
        s = strrep(s,'transX','FLAP-DISP');
        s = strrep(s,'transY','EDGE-DISP');
        s = strrep(s,'rotX','EDGE');
        s = strrep(s,'rotY','FLAP');
        if Desc{iD}(1:2)=='BD'
            nBD=nBD+1;
        elseif Desc{iD}(1:2)=='ED'
            DescCatED = [s ' - ' DescCatED];
        else
            DescCat = [DescCat ' - ' s];
        end
    end
    DescCat=[DescCatED, DescCat];
    if nBD>0
        DescCat = sprintf('BD%d/%d %s',nBD,sum(CD.Modes{i}.StateHasMaxAtThisMode),DescCat);
    end
    fprintf('%8.3f ; %11.8f ; %s\n',CD.NaturalFreq_Hz(i),CD.DampingRatio(i),DescCat(1:min(80,length(DescCat))));
end
