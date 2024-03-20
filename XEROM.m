%%
clear variables; close all; clc;


% %% create sparse CR matrix
% 
% CR_SA1 = zeros(size(D1));
% CR_SA2 = zeros(size(D1));
% %del_SA1 = 8.9E-4; % delta Sigma_a for the fast group
% %del_SA2 = 3.4E-4; % delta Sigma_a for the thermal group
% 
% del_SA1 = 8.9; % delta Sigma_a for the fast group
% del_SA2 = 3.4; % delta Sigma_a for the thermal group
% 
% 
% CR_SA1(12:13,12:13,1:11) = del_SA1; CR_SA2(12:13,12:13,1:11) = del_SA2; % 1st control rod thermal and fast
% CR_SA1(12:13,20:21,1:11) = del_SA1; CR_SA2(12:13,20:21,1:11) = del_SA2; % 2nd control rod thermal and fast
% CR_SA1(20:21,12:13,1:11) = del_SA1; CR_SA2(20:21,12:13,1:11) = del_SA2; % 3rd control rod thermal and fast
% CR_SA1(20:21,20:21,1:11) = del_SA1; CR_SA2(20:21,20:21,1:11) = del_SA2; % 4th control rod thermal and fast
% CR_SA1(16:17,28:29,1:11) = del_SA1; CR_SA2(16:17,28:29,1:11) = del_SA2; % 5th control rod thermal and fast
% CR_SA1(28:29,16:17,1:11) = del_SA1; CR_SA2(28:29,16:17,1:11) = del_SA2; % 6th control rod thermal and fast
% CR_SA1(4:5,16:17,1:11) = del_SA1; CR_SA2(4:5,16:17,1:11) = del_SA2; % 7th control rod thermal and fast
% CR_SA1(16:17,4:5,1:11) = del_SA1; CR_SA2(16:17,4:5,1:11) = del_SA2; % 8th control rod thermal and fast
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                             CORE_SIM_v1.3
%                Copyright (C) 2011, 2020  Christophe Demazière
%                   Chalmers University of Technology
%                     Department of Applied Physics
%                    Division of Nuclear Engineering
%                      SE-412 96 Gothenburg, Sweden
%                        E-mail: demaz@chalmers.se
%                          Tel: +46-31-772 3082
%                          Fax: +46-31-772 3079
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                          MAIN PROGRAM FILE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% For explanatory notes about the use of this tool, see the file
% USERS_GUIDE.PDF in the directory "manuals".
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%
% INITIALIZATION
%%%%%%%%%%%%%%%%%

clc

fprintf('\nCORE SIM v1.2  Copyright (C) 2011  Christophe Demazière\nThis program comes with ABSOLUTELY NO WARRANTY.\nThis is free software, and you are welcome to redistribute it under certain conditions.\nFor details, type "open(''license.txt'')".\n\n')

pause(2)

fprintf('\nINITIALIZATION IN PROGRESS.\n')

clear variables
format compact
RUN_TEST=1;

% Parameters for sparse matrix operations (can be modified in case of
% memory problems)
spparms('default')

% Loading parameters for program execution:
% VIZ, EIG_MET, BYP, neigs, m, n_restart, conv_ERAM, n_iter, conv_POW
if exist('SETTINGS.m','file')==2
    SETTINGS
else
    fprintf('\nERROR: Your SETTINGS.m file is missing.\nEXECUTION TERMINATED\n')
    return
end

VAR={'VIZ','EIG_MET','BYP','neigs','m','n_restart','conv_ERAM','n_iter','conv_POW'};
TEST_VAR=ismember((VAR),who);

for i=1:size(VAR,2)
    if TEST_VAR(1,i)==0
        fprintf('\nERROR: The variable %s is missing from your SETTINGS.m file.\n',VAR{i});
        RUN_TEST=0;
    end
end
if RUN_TEST==0
    fprintf('\nEXECUTION TERMINATED\n')
    return
end
clear VAR TEST_VAR files

if ((EIG_MET~=1) && (EIG_MET~=2))
    fprintf('\nERROR: The value of your EIG_MET variable in your SETTINGS.m file does\nnot correspond to one of the permissible values.\nEXECUTION TERMINATED\n')
    return
end

if ((BYP~=0) && (BYP~=1))
    fprintf('\nERROR: The value of your BYP variable in your SETTINGS.m file does\nnot correspond to one of the permissible values.\nEXECUTION TERMINATED\n')
    return
end

if ((VIZ~=0) && (VIZ~=1))
    fprintf('\nERROR: The value of your VIZ variable in your SETTINGS.m file does\nnot correspond to one of the permissible values.\nEXECUTION TERMINATED\n')
    return
end

if exist('output','dir')~=7
    mkdir('output')
end

% Loading input variables:
% D1 D2 REM ABS1 ABS2 NUFIS1 NUFIS2 DX DY DZ
files={
    'XS_data_100_50.mat'
    'GEOM_data_100_50.mat'
    'additional_data_100_50.mat'
    'RESULTS.mat'
    };
for i=1:size(files,1)
    if exist(sprintf('../For_Christophe/input_100_50/%s',files{i}),'file')==2
        load (sprintf('../For_Christophe/input_100_50/%s',files{i}))
        fprintf("input data /%s loaded \n",files{i})
    else
        fprintf('\nERROR: Missing input file or missing input directory.\nEXECUTION TERMINATED\n')
        return
    end
end



VAR={'D1','D2','REM','ABS1','ABS2','NUFIS1','NUFIS2','DX','DY','DZ','STA_FLX1','STA_FLX2','KF1','KF2'};
TEST_VAR=ismember((VAR),who);

for i=1:size(VAR,2)
    if TEST_VAR(1,i)==0
        fprintf('\nERROR: The variable %s is missing from your input files.\n',VAR{i});
        RUN_TEST=0;
    end
end
if RUN_TEST==0
    fprintf('\nEXECUTION TERMINATED\n')
    return
end
clear VAR TEST_VAR files
RUN_CORESIM = 0;

% Checking the geometry of the core
I_MAX=size(D2,1);
J_MAX=size(D2,2);
K_MAX=size(D2,3);

if I_MAX<3
    fprintf('\nERROR: You need at least 3 rows to properly define this 3-D problem.\nEXECUTION TERMINATED\n')
    return
end

if J_MAX<3
    fprintf('\nERROR: You need at least 3 columns to properly define this 3-D problem.\nEXECUTION TERMINATED\n')
    return
end

if K_MAX<3
    fprintf('\nERROR: You need at least 3 axial planes to properly define this 3-D problem.\nEXECUTION TERMINATED\n')
    return
end

% Creating a variable TYP defining the location of the fuel and reflector nodes
% (TYP=1 for fuel nodes and TYP=-1 for reflector nodes)
TYP=zeros(I_MAX,J_MAX,K_MAX);
for I=1:I_MAX
    for J=1:J_MAX
        for K=1:K_MAX
            if D2(I,J,K)==0
                TYP(I,J,K)=0;
            elseif (D2(I,J,K)~=0 && NUFIS2(I,J,K)==0)
                TYP(I,J,K)=-1;
            elseif (D2(I,J,K)~=0 && NUFIS2(I,J,K)~=0)
                TYP(I,J,K)=1;
            end
        end
    end
end



if RUN_TEST==0
    fprintf('\nEXECUTION TERMINATED\n')
    return
end

% Creating a variable allowing reordering the 3-D variables into column vectors
tmp=0;
CONV=zeros(I_MAX,J_MAX,K_MAX);
SHIFT_XY=zeros(K_MAX,1);
for K=1:K_MAX
    for I=1:I_MAX
        for J=1:J_MAX
            if TYP(I,J,K)~=0
                tmp=tmp+1;
                CONV(I,J,K)=tmp;
            else
                CONV(I,J,K)=0;
            end
        end
    end
    if K>1
        SHIFT_XY(K)=max(max(CONV(:,:,K)))-max(max(CONV(:,:,K-1)));
    else
        SHIFT_XY(K)=max(max(CONV(:,:,K)));
    end
end
clear tmp
SHIFT_XY_test=SHIFT_XY(1);
for K=2:K_MAX
    if SHIFT_XY(K)~=SHIFT_XY_test
        fprintf('\nERROR: Each of the axial planes should contain the same number of nodes.\nEXECUTION TERMINATED\n')
        return
    end
end
clear SHIFT_XY
SHIFT_XY=SHIFT_XY_test;
clear SHIFT_XY_test
SHIFT_XYZ=max(max(max(CONV)));

% Creating variables allowing determining the positions of consecutive
% nodes in the matrices/vectors for a positive and negative, respectively,
% increment of the row index (IP_SHIFT and IM_SHIFT, respectively), column
% index (JP_SHIFT and JM_SHIFT, respectively), and axial plane index
% (KP_SHIFT and KM_SHIFT, respectively)
IP_SHIFT=zeros(I_MAX,J_MAX,K_MAX);
IM_SHIFT=zeros(I_MAX,J_MAX,K_MAX);
JP_SHIFT=zeros(I_MAX,J_MAX,K_MAX);
JM_SHIFT=zeros(I_MAX,J_MAX,K_MAX);
KP_SHIFT=zeros(I_MAX,J_MAX,K_MAX);
KM_SHIFT=zeros(I_MAX,J_MAX,K_MAX);
for I=1:I_MAX
    for J=1:J_MAX
        for K=1:K_MAX
            if TYP(I,J,K)~=0
                if ((I>1) && (TYP(I-1,J,K)~=0))
                    IM_SHIFT(I,J,K)=CONV(I-1,J,K)-CONV(I,J,K);
                end
                if ((I<I_MAX) && (TYP(I+1,J,K)~=0))
                    IP_SHIFT(I,J,K)=CONV(I+1,J,K)-CONV(I,J,K);
                end
                if ((J>1) && (TYP(I,J-1,K)~=0))
                    JM_SHIFT(I,J,K)=CONV(I,J-1,K)-CONV(I,J,K);
                end
                if ((J<J_MAX) && (TYP(I,J+1,K)~=0))
                    JP_SHIFT(I,J,K)=CONV(I,J+1,K)-CONV(I,J,K);
                end
                if ((K>1) && (TYP(I,J,K-1)~=0))
                    KM_SHIFT(I,J,K)=CONV(I,J,K-1)-CONV(I,J,K);
                end
                if ((K<K_MAX) && (TYP(I,J,K+1)~=0))
                    KP_SHIFT(I,J,K)=CONV(I,J,K+1)-CONV(I,J,K);
                end
            end
        end
    end
end

% Calculating the coupling coefficients between nodes according to the
% box-scheme (finite differences)
SHIFT=0;
for K=1:K_MAX
    % AX
    ax11_tmp=zeros(SHIFT_XY,1);
    ax22_tmp=zeros(SHIFT_XY,1);
    for I=1:I_MAX
        for J=1:J_MAX
            if (TYP(I,J,K)~=0)
                if ( ((I>1) && (TYP(I,J,K)~=0) && (TYP((I-1),J,K)==0)) || ((I==1) && (TYP(I,J,K)~=0)) )
                    % Group 1->1
                    ax11_tmp(CONV(I,J,K)-SHIFT,1)=2*D1((I+1),J,K)*D1(I,J,K)/(DX*(D1((I+1),J,K)+D1(I,J,K)))+1/2/(1+DX/(4*D1(I,J,K)));
                    % Group 2->2
                    ax22_tmp(CONV(I,J,K)-SHIFT,1)=2*D2((I+1),J,K)*D2(I,J,K)/(DX*(D2((I+1),J,K)+D2(I,J,K)))+1/2/(1+DX/(4*D2(I,J,K)));
                elseif ( ((I<I_MAX) && (TYP(I,J,K)~=0) && (TYP((I+1),J,K)==0)) || ((I==I_MAX) && (TYP(I,J,K)~=0)) )
                    % Group 1->1
                    ax11_tmp(CONV(I,J,K)-SHIFT,1)=2*D1(I,J,K)*D1((I-1),J,K)/(DX*(D1(I,J,K)+D1((I-1),J,K)))+1/2/(1+DX/(4*D1(I,J,K)));
                    % Group 2->2
                    ax22_tmp(CONV(I,J,K)-SHIFT,1)=2*D2(I,J,K)*D2((I-1),J,K)/(DX*(D2(I,J,K)+D2((I-1),J,K)))+1/2/(1+DX/(4*D2(I,J,K)));
                else
                    % Group 1->1
                    ax11_tmp(CONV(I,J,K)-SHIFT,1)=2*D1((I+1),J,K)*D1(I,J,K)/(DX*(D1((I+1),J,K)+D1(I,J,K)))+2*D1(I,J,K)*D1((I-1),J,K)/(DX*(D1(I,J,K)+D1((I-1),J,K)));
                    % Group 2->2
                    ax22_tmp(CONV(I,J,K)-SHIFT,1)=2*D2((I+1),J,K)*D2(I,J,K)/(DX*(D2((I+1),J,K)+D2(I,J,K)))+2*D2(I,J,K)*D2((I-1),J,K)/(DX*(D2(I,J,K)+D2((I-1),J,K)));
                end
            end
        end
    end
    ax11=spdiags(ax11_tmp,0,SHIFT_XY,SHIFT_XY);
    ax22=spdiags(ax22_tmp,0,SHIFT_XY,SHIFT_XY);
    AX(SHIFT+1:SHIFT+SHIFT_XY,SHIFT+1:SHIFT+SHIFT_XY)=ax11;
    AX(SHIFT_XYZ+SHIFT+1:SHIFT_XYZ+SHIFT+SHIFT_XY,SHIFT_XYZ+SHIFT+1:SHIFT_XYZ+SHIFT+SHIFT_XY)=ax22;
    clear ax11 ax22 ax11_tmp ax22_tmp
 
    
    %%
    % BX
    bx11_tmp=zeros(max(max(max(CONV)))+max(max(abs(IP_SHIFT(:,:,K))))-SHIFT,max(max(abs(IP_SHIFT(:,:,K)))));
    bx22_tmp=zeros(max(max(max(CONV)))+max(max(abs(IP_SHIFT(:,:,K))))-SHIFT,max(max(abs(IP_SHIFT(:,:,K)))));
    for C=1:max(max(abs(IP_SHIFT(:,:,K))))
        for I=1:I_MAX
            for J=1:J_MAX
                if (IP_SHIFT(I,J,K)==C)
                    if ( ((I>1) && (TYP(I,J,K)~=0) && (TYP((I-1),J,K)==0)) || ((I==1) && (TYP(I,J,K)~=0)) )
                        % Group 1->1
                        bx11_tmp(CONV(I,J,K)+C-SHIFT,C)=-2*D1((I+1),J,K)*D1(I,J,K)/(DX*(D1((I+1),J,K)+D1(I,J,K)));
                        % Group 2->2
                        bx22_tmp(CONV(I,J,K)+C-SHIFT,C)=-2*D2((I+1),J,K)*D2(I,J,K)/(DX*(D2((I+1),J,K)+D2(I,J,K)));
                    elseif ( ((I<I_MAX) && (TYP(I,J,K)~=0) && (TYP((I+1),J,K)==0)) || ((I==I_MAX) && (TYP(I,J,K)~=0)) )
                    else
                        % Group 1->1
                        bx11_tmp(CONV(I,J,K)+C-SHIFT,C)=-2*D1((I+1),J,K)*D1(I,J,K)/(DX*(D1((I+1),J,K)+D1(I,J,K)));
                        % Group 2->2
                        bx22_tmp(CONV(I,J,K)+C-SHIFT,C)=-2*D2((I+1),J,K)*D2(I,J,K)/(DX*(D2((I+1),J,K)+D2(I,J,K)));
                    end
                end
            end
        end
    end

    if size(bx11_tmp,1)<SHIFT_XY
        bx11_tmp(SHIFT_XY,1)=0;
    end
    if size(bx22_tmp,1)<SHIFT_XY
        bx22_tmp(SHIFT_XY,1)=0;
    end
    d=zeros(max(max(abs(IP_SHIFT(:,:,K)))),1);
    for C=1:max(max(abs(IP_SHIFT(:,:,K))))
        d(C,1)=C;
    end
    bx11=spdiags(bx11_tmp,d,SHIFT_XY,SHIFT_XY);
    bx22=spdiags(bx22_tmp,d,SHIFT_XY,SHIFT_XY);    
    BX(SHIFT+1:SHIFT+SHIFT_XY,SHIFT+1:SHIFT+SHIFT_XY)=bx11;
    BX(SHIFT_XYZ+SHIFT+1:SHIFT_XYZ+SHIFT+SHIFT_XY,SHIFT_XYZ+SHIFT+1:SHIFT_XYZ+SHIFT+SHIFT_XY)=bx22;
    clear bx11 bx22 bx11_tmp bx22_tmp C d
    
    % CX
    cx11_tmp=zeros(max(max(max(CONV)))-max(max(abs(IM_SHIFT(:,:,K))))-SHIFT,max(max(abs(IP_SHIFT(:,:,K)))));
    cx22_tmp=zeros(max(max(max(CONV)))-max(max(abs(IM_SHIFT(:,:,K))))-SHIFT,max(max(abs(IP_SHIFT(:,:,K)))));
    for C=1:max(max(abs(IM_SHIFT(:,:,K))))
        for I=1:I_MAX
            for J=1:J_MAX
                if (IM_SHIFT(I,J)==-C)
                    if ( ((I>1) && (TYP(I,J,K)~=0) && (TYP((I-1),J,K)==0)) || ((I==1) && (TYP(I,J,K)~=0)) )
                    elseif ( ((I<I_MAX) && (TYP(I,J,K)~=0) && (TYP((I+1),J,K)==0)) || ((I==I_MAX) && (TYP(I,J,K)~=0)) )
                        % Group 1->1
                        cx11_tmp(CONV(I,J,K)-C-SHIFT,C)=-2*D1(I,J,K)*D1((I-1),J,K)/(DX*(D1(I,J,K)+D1((I-1),J,K)));
                        % Group 2->2
                        cx22_tmp(CONV(I,J,K)-C-SHIFT,C)=-2*D2(I,J,K)*D2((I-1),J,K)/(DX*(D2(I,J,K)+D2((I-1),J,K)));
                    else
                        % Group 1->1
                        cx11_tmp(CONV(I,J,K)-C-SHIFT,C)=-2*D1(I,J,K)*D1((I-1),J,K)/(DX*(D1(I,J,K)+D1((I-1),J,K)));
                        % Group 2->2
                        cx22_tmp(CONV(I,J,K)-C-SHIFT,C)=-2*D2(I,J,K)*D2((I-1),J,K)/(DX*(D2(I,J,K)+D2((I-1),J,K)));
                    end
                end
            end
        end
    end

    if size(cx11_tmp,1)<SHIFT_XY
        cx11_tmp(SHIFT_XY,1)=0;
    end
    if size(cx22_tmp,1)<SHIFT_XY
        cx22_tmp(SHIFT_XY,1)=0;
    end
    d=zeros(max(max(abs(IM_SHIFT(:,:,K)))),1);
    for C=1:max(max(abs(IM_SHIFT(:,:,K))))
        d(C,1)=-C;
    end
    cx11=spdiags(cx11_tmp,d,SHIFT_XY,SHIFT_XY);
    cx22=spdiags(cx22_tmp,d,SHIFT_XY,SHIFT_XY);
    CX(SHIFT+1:SHIFT+SHIFT_XY,SHIFT+1:SHIFT+SHIFT_XY)=cx11;
    CX(SHIFT_XYZ+SHIFT+1:SHIFT_XYZ+SHIFT+SHIFT_XY,SHIFT_XYZ+SHIFT+1:SHIFT_XYZ+SHIFT+SHIFT_XY)=cx22;
    clear cx11 cx22 cx11_tmp cx22_tmp C d
    
    % AY
    ay11_tmp=zeros(SHIFT_XY,1);
    ay22_tmp=zeros(SHIFT_XY,1);
    for I=1:I_MAX
        for J=1:J_MAX
            if (TYP(I,J,K)~=0)
                if ( ((J>1) && (TYP(I,J,K)~=0) && (TYP(I,(J-1),K)==0)) || ((J==1) && (TYP(I,J,K)~=0)) )
                    % Group 1->1
                    ay11_tmp(CONV(I,J,K)-SHIFT,1)=2*D1(I,(J+1),K)*D1(I,J,K)/(DY*(D1(I,(J+1),K)+D1(I,J,K)))+1/2/(1+DY/(4*D1(I,J,K)));
                    % Group 2->2
                    ay22_tmp(CONV(I,J,K)-SHIFT,1)=2*D2(I,(J+1),K)*D2(I,J,K)/(DY*(D2(I,(J+1),K)+D2(I,J,K)))+1/2/(1+DY/(4*D2(I,J,K)));
                elseif ( ((J<J_MAX) && (TYP(I,J,K)~=0) && (TYP(I,(J+1),K)==0)) || ((J==J_MAX) && (TYP(I,J,K)~=0)) )
                    % Group 1->1
                    ay11_tmp(CONV(I,J,K)-SHIFT,1)=2*D1(I,J,K)*D1(I,(J-1),K)/(DY*(D1(I,J,K)+D1(I,(J-1),K)))+1/2/(1+DY/(4*D1(I,J,K)));
                    % Group 2->2
                    ay22_tmp(CONV(I,J,K)-SHIFT,1)=2*D2(I,J,K)*D2(I,(J-1),K)/(DY*(D2(I,J,K)+D2(I,(J-1),K)))+1/2/(1+DY/(4*D2(I,J,K)));
                else
                    % Group 1->1
                    ay11_tmp(CONV(I,J,K)-SHIFT,1)=2*D1(I,(J+1),K)*D1(I,J,K)/(DY*(D1(I,(J+1),K)+D1(I,J,K)))+2*D1(I,J,K)*D1(I,(J-1),K)/(DY*(D1(I,J,K)+D1(I,(J-1),K)));
                    % Group 2->2
                    ay22_tmp(CONV(I,J,K)-SHIFT,1)=2*D2(I,(J+1),K)*D2(I,J,K)/(DY*(D2(I,(J+1),K)+D2(I,J,K)))+2*D2(I,J,K)*D2(I,(J-1),K)/(DY*(D2(I,J,K)+D2(I,(J-1),K)));
                end
            end
        end
    end
    ay11=spdiags(ay11_tmp,0,SHIFT_XY,SHIFT_XY);
    ay22=spdiags(ay22_tmp,0,SHIFT_XY,SHIFT_XY);
    AY(SHIFT+1:SHIFT+SHIFT_XY,SHIFT+1:SHIFT+SHIFT_XY)=ay11;
    AY(SHIFT_XYZ+SHIFT+1:SHIFT_XYZ+SHIFT+SHIFT_XY,SHIFT_XYZ+SHIFT+1:SHIFT_XYZ+SHIFT+SHIFT_XY)=ay22;
    clear ay11 ay22 ay11_tmp ay22_tmp

    % BY
    by11_tmp=zeros(SHIFT_XY+1,1);
    by22_tmp=zeros(SHIFT_XY+1,1);
    for I=1:I_MAX
        for J=1:J_MAX
            if (TYP(I,J,K)~=0)
                if ( ((J>1) && (TYP(I,J,K)~=0) && (TYP(I,(J-1),K)==0)) || ((J==1) && (TYP(I,J,K)~=0)) )
                    % Group 1->1
                    by11_tmp(CONV(I,J,K)+1-SHIFT,1)=-2*D1(I,(J+1),K)*D1(I,J,K)/(DY*(D1(I,(J+1),K)+D1(I,J,K)));
                    % Group 2->2
                    by22_tmp(CONV(I,J,K)+1-SHIFT,1)=-2*D2(I,(J+1),K)*D2(I,J,K)/(DY*(D2(I,(J+1),K)+D2(I,J,K)));
                elseif ( ((J<J_MAX) && (TYP(I,J,K)~=0) && (TYP(I,(J+1),K)==0)) || ((J==J_MAX) && (TYP(I,J,K)~=0)) )
                else
                    % Group 1->1
                    by11_tmp(CONV(I,J,K)+1-SHIFT,1)=-2*D1(I,(J+1),K)*D1(I,J,K)/(DY*(D1(I,(J+1),K)+D1(I,J,K)));
                    % Group 2->2
                    by22_tmp(CONV(I,J,K)+1-SHIFT,1)=-2*D2(I,(J+1),K)*D2(I,J,K)/(DY*(D2(I,(J+1),K)+D2(I,J,K)));
                end
            end
        end
    end
    by11=spdiags(by11_tmp,1,SHIFT_XY,SHIFT_XY);
    by22=spdiags(by22_tmp,1,SHIFT_XY,SHIFT_XY);
    BY(SHIFT+1:SHIFT+SHIFT_XY,SHIFT+1:SHIFT+SHIFT_XY)=by11;
    BY(SHIFT_XYZ+SHIFT+1:SHIFT_XYZ+SHIFT+SHIFT_XY,SHIFT_XYZ+SHIFT+1:SHIFT_XYZ+SHIFT+SHIFT_XY)=by22;
    clear by11 by22 by11_tmp by22_tmp
    
    % CY
    cy11_tmp=zeros(SHIFT_XY-1,1);
    cy22_tmp=zeros(SHIFT_XY-1,1);
    for I=1:I_MAX
        for J=1:J_MAX
            if (TYP(I,J,K)~=0)
                if ( ((J>1) && (TYP(I,J,K)~=0) && (TYP(I,(J-1),K)==0)) || ((J==1) && (TYP(I,J,K)~=0)) )
                elseif ( ((J<J_MAX) && (TYP(I,J,K)~=0) && (TYP(I,(J+1),K)==0)) || ((J==J_MAX) && (TYP(I,J,K)~=0)) )
                    % Group 1->1
                    cy11_tmp(CONV(I,J,K)-1-SHIFT,1)=-2*D1(I,J,K)*D1(I,(J-1),K)/(DY*(D1(I,J,K)+D1(I,(J-1),K)));
                    % Group 2->2
                    cy22_tmp(CONV(I,J,K)-1-SHIFT,1)=-2*D2(I,J,K)*D2(I,(J-1),K)/(DY*(D2(I,J,K)+D2(I,(J-1),K)));
                else
                    % Group 1->1
                    cy11_tmp(CONV(I,J,K)-1-SHIFT,1)=-2*D1(I,J,K)*D1(I,(J-1),K)/(DY*(D1(I,J,K)+D1(I,(J-1),K)));
                    % Group 2->2
                    cy22_tmp(CONV(I,J,K)-1-SHIFT,1)=-2*D2(I,J,K)*D2(I,(J-1),K)/(DY*(D2(I,J,K)+D2(I,(J-1),K)));
                end
            end
        end
    end
    cy11=spdiags(cy11_tmp,-1,SHIFT_XY,SHIFT_XY);
    cy22=spdiags(cy22_tmp,-1,SHIFT_XY,SHIFT_XY);
    CY(SHIFT+1:SHIFT+SHIFT_XY,SHIFT+1:SHIFT+SHIFT_XY)=cy11;
    CY(SHIFT_XYZ+SHIFT+1:SHIFT_XYZ+SHIFT+SHIFT_XY,SHIFT_XYZ+SHIFT+1:SHIFT_XYZ+SHIFT+SHIFT_XY)=cy22;
    clear cy11 cy22 cy11_tmp cy22_tmp

    % AZ
    az11_tmp=zeros(SHIFT_XY,1);
    az22_tmp=zeros(SHIFT_XY,1);
    for I=1:I_MAX
        for J=1:J_MAX
            if (TYP(I,J,K)~=0)
                if ( K==1 )
                    % Group 1->1
                    az11_tmp(CONV(I,J,K)-SHIFT,1)=2*D1(I,J,(K+1))*D1(I,J,K)/(DZ*(D1(I,J,(K+1))+D1(I,J,K)))+1/2/(1+DZ/(4*D1(I,J,K)));
                    % Group 2->2
                    az22_tmp(CONV(I,J,K)-SHIFT,1)=2*D2(I,J,(K+1))*D2(I,J,K)/(DZ*(D2(I,J,(K+1))+D2(I,J,K)))+1/2/(1+DZ/(4*D2(I,J,K)));
                elseif ( K==K_MAX )
                    % Group 1->1
                    az11_tmp(CONV(I,J,K)-SHIFT,1)=2*D1(I,J,K)*D1(I,J,(K-1))/(DZ*(D1(I,J,K)+D1(I,J,(K-1))))+1/2/(1+DZ/(4*D1(I,J,K)));
                    % Group 2->2
                    az22_tmp(CONV(I,J,K)-SHIFT,1)=2*D2(I,J,K)*D2(I,J,(K-1))/(DZ*(D2(I,J,K)+D2(I,J,(K-1))))+1/2/(1+DZ/(4*D2(I,J,K)));
                else
                    % Group 1->1
                    az11_tmp(CONV(I,J,K)-SHIFT,1)=2*D1(I,J,(K+1))*D1(I,J,K)/(DZ*(D1(I,J,(K+1))+D1(I,J,K)))+2*D1(I,J,K)*D1(I,J,(K-1))/(DZ*(D1(I,J,K)+D1(I,J,(K-1))));
                    % Group 2->2
                    az22_tmp(CONV(I,J,K)-SHIFT,1)=2*D2(I,J,(K+1))*D2(I,J,K)/(DZ*(D2(I,J,(K+1))+D2(I,J,K)))+2*D2(I,J,K)*D2(I,J,(K-1))/(DZ*(D2(I,J,K)+D2(I,J,(K-1))));
                end
            end
        end
    end
    az11=spdiags(az11_tmp,0,SHIFT_XY,SHIFT_XY);
    az22=spdiags(az22_tmp,0,SHIFT_XY,SHIFT_XY);
    AZ(SHIFT+1:SHIFT+SHIFT_XY,SHIFT+1:SHIFT+SHIFT_XY)=az11;
    AZ(SHIFT_XYZ+SHIFT+1:SHIFT_XYZ+SHIFT+SHIFT_XY,SHIFT_XYZ+SHIFT+1:SHIFT_XYZ+SHIFT+SHIFT_XY)=az22;
    clear az11 az22 az11_tmp az22_tmp
    
    % BZ
    bz11_tmp=zeros(SHIFT_XY,1);
    bz22_tmp=zeros(SHIFT_XY,1);
    for I=1:I_MAX
        for J=1:J_MAX
            if (TYP(I,J,K)~=0)
                if ( K==1 )
                    % Group 1->1
                    bz11_tmp(CONV(I,J,K)-SHIFT,1)=-2*D1(I,J,(K+1))*D1(I,J,K)/(DZ*(D1(I,J,(K+1))+D1(I,J,K)));
                    % Group 2->2
                    bz22_tmp(CONV(I,J,K)-SHIFT,1)=-2*D2(I,J,(K+1))*D2(I,J,K)/(DZ*(D2(I,J,(K+1))+D2(I,J,K)));
                elseif ( K==K_MAX )
                else
                    % Group 1->1
                    bz11_tmp(CONV(I,J,K)-SHIFT,1)=-2*D1(I,J,(K+1))*D1(I,J,K)/(DZ*(D1(I,J,(K+1))+D1(I,J,K)));
                    % Group 2->2
                    bz22_tmp(CONV(I,J,K)-SHIFT,1)=-2*D2(I,J,(K+1))*D2(I,J,K)/(DZ*(D2(I,J,(K+1))+D2(I,J,K)));
                end
            end
        end
    end
    if K<K_MAX
        bz11=spdiags(bz11_tmp,0,SHIFT_XY,SHIFT_XY);
        bz22=spdiags(bz22_tmp,0,SHIFT_XY,SHIFT_XY);
        BZ(SHIFT+1:SHIFT+SHIFT_XY,SHIFT+1+SHIFT_XY:SHIFT+2*SHIFT_XY)=bz11;
        BZ(SHIFT_XYZ+SHIFT+1:SHIFT_XYZ+SHIFT+SHIFT_XY,SHIFT_XYZ+SHIFT+1+SHIFT_XY:SHIFT_XYZ+SHIFT+2*SHIFT_XY)=bz22;
    end
    clear bz11 bz22 bz11_tmp bz22_tmp
    
    % CZ
    cz11_tmp=zeros(SHIFT_XY,1);
    cz22_tmp=zeros(SHIFT_XY,1);
    for I=1:I_MAX
        for J=1:J_MAX
            if (TYP(I,J,K)~=0)
                if ( K==1 )
                elseif ( K==K_MAX )
                    % Group 1->1
                    cz11_tmp(CONV(I,J,K)-SHIFT,1)=-2*D1(I,J,K)*D1(I,J,(K-1))/(DZ*(D1(I,J,K)+D1(I,J,(K-1))));
                    % Group 2->2
                    cz22_tmp(CONV(I,J,K)-SHIFT,1)=-2*D2(I,J,K)*D2(I,J,(K-1))/(DZ*(D2(I,J,K)+D2(I,J,(K-1))));
                else
                    % Group 1->1
                    cz11_tmp(CONV(I,J,K)-SHIFT,1)=-2*D1(I,J,K)*D1(I,J,(K-1))/(DZ*(D1(I,J,K)+D1(I,J,(K-1))));
                    % Group 2->2
                    cz22_tmp(CONV(I,J,K)-SHIFT,1)=-2*D2(I,J,K)*D2(I,J,(K-1))/(DZ*(D2(I,J,K)+D2(I,J,(K-1))));
                end
            end
        end
    end
    if K>1
        cz11=spdiags(cz11_tmp,0,SHIFT_XY,SHIFT_XY);
        cz22=spdiags(cz22_tmp,0,SHIFT_XY,SHIFT_XY);
        CZ(SHIFT+1:SHIFT+SHIFT_XY,SHIFT+1-SHIFT_XY:SHIFT)=cz11;
        CZ(SHIFT_XYZ+SHIFT+1:SHIFT_XYZ+SHIFT+SHIFT_XY,SHIFT_XYZ+SHIFT+1-SHIFT_XY:SHIFT_XYZ+SHIFT)=cz22;
    end
    clear cz11 cz22 cz11_tmp cz22_tmp

    SHIFT=SHIFT+SHIFT_XY;
end

if size(BZ,1)<2*SHIFT_XYZ
    BZ(2*SHIFT_XYZ,2*SHIFT_XYZ)=0;
end
if size(CZ,2)<2*SHIFT_XYZ
    CZ(2*SHIFT_XYZ,2*SHIFT_XYZ)=0;
end
save("temp/temp_input_vars.mat","AX","AY","AZ","BX","BY","BZ","CX","CY","CZ","DX","DY","DZ","D1","D2","NUFIS1","NUFIS2","ABS1","ABS2","REM","SHIFT_XYZ","SHIFT","SHIFT_XY","I_MAX","J_MAX","K_MAX","TYP","CONV","STA_FLX1","STA_FLX2","KF1","KF2")

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATING THE EIGENMODES OF THE SOURCE-FREE PROBLEM (FORWARD PROBLEM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if RUN_CORESIM
    fprintf('\nCALCULATION OF THE EIGENMODES OF THE SOURCE-FREE PROBLEM (FORWARD PROBLEM) IN PROGRESS.\n')
    Dp_FLX_tmp=zeros(2*SHIFT_XYZ,3);
    for I=1:I_MAX
    	for J=1:J_MAX
            for K=1:K_MAX
                if (TYP(I,J,K)~=0)
                    % Group 1->1
                    Dp_FLX_tmp(CONV(I,J,K),1)=-ABS1(I,J,K)-REM(I,J,K);
                    % Group 2->1
                    Dp_FLX_tmp(CONV(I,J,K)+SHIFT_XYZ,2)=0;
                    % Group 1->2
                    Dp_FLX_tmp(CONV(I,J,K),3)=REM(I,J,K);
                    % Group 2->2
                    Dp_FLX_tmp(CONV(I,J,K)+SHIFT_XYZ,1)=-ABS2(I,J,K);
                end
            end
    	end
    end
    Dp_FLX=spdiags(Dp_FLX_tmp,[0;SHIFT_XYZ;-SHIFT_XYZ],2*SHIFT_XYZ,2*SHIFT_XYZ);
    STA=Dp_FLX-AX/DX-AY/DY-AZ/DZ-BX/DX-BY/DY-BZ/DZ-CX/DX-CY/DY-CZ/DZ;
    clear Dp_FLX Dp_FLX_tmp

    F_tmp=zeros(2*SHIFT_XYZ,2);
    for I=1:I_MAX
    	for J=1:J_MAX
            for K=1:K_MAX
                if TYP(I,J,K)~=0
                    % Group 1->1
                    F_tmp(CONV(I,J,K),1)=-NUFIS1(I,J,K);
                    % Group 2->1
                    F_tmp(CONV(I,J,K)+SHIFT_XYZ,2)=-NUFIS2(I,J,K);
                end
            end
    	end
    end
    F=spdiags(F_tmp,[0;SHIFT_XYZ],2*SHIFT_XYZ,2*SHIFT_XYZ);
    clear F_tmp

    count=0;
    q_old=ones(2*SHIFT_XYZ,1);
    q_old=q_old/norm(q_old);

    if EIG_MET==1
        res=Inf;
        [L,U,P,Q]=lu(STA);
        while res>conv_ERAM
            v=q_old;
            H=zeros(m,m);
            for j=1:m
                w=Q*(U\(L\(P*(F*v(:,j)))));
                for i=1:j
                    H(i,j)=w'*v(:,i);
                    w=w-H(i,j)*v(:,i);
                end
                H(j+1,j)=norm(w,2);
                if H(j+1,j)==0
                    break
                else
                    v(:,j+1)=w/norm(w,2);
                end
            end
            Hr=H(1:(size(H,1)-1),1:size(H,2));
            [X_H,D_H]=eig(Hr);
            D_H=diag(D_H);
            Vr=v(1:size(v,1),1:size(v,2)-1);
            X=Vr*X_H;
            [lambda,order]=sort(D_H,'descend');
            X=X(:,order);
            for i=1:(m-1)
                X(:,i)=X(:,i)/norm(X(:,i),2);
            end
            res_X=zeros(neigs,1);
            for i=1:neigs
                res_X(i)=norm((lambda(i)*STA*X(:,i)-F*X(:,i)),2);
            end
            res=norm(res_X,inf);
            if res>conv_ERAM
                count=count+1;
                fprintf('Residuals before restart # %d: %1.3e\n',count,res);
            else
                fprintf('Residuals on the converged eigenvectors: %1.3e\n',res);
            end
            q_new=zeros(size(X,1),1);
            for i=1:neigs
                q_new=q_new+X(:,i).*1;
            end
            q_old=q_new;
            if count>n_restart
                fprintf('\nERROR: The Explicitely Restarted Arnoldi Method (ERAM) is not converging as expected.\nEither change the parameters of the Explicitely Restarted Arnoldi Method in your SETTINGS.m file\n')
                fprintf('or switch to the power iteration method (POW)\n(setting EIG_MET=2 in your SETTINGS.m file).\n')
                if BYP==0
                    fprintf('EXECUTION TERMINATED\n')
                    return
                else
                    fprintf('EXECUTION ALLOWED TO CONTINUE\n')
                    break
                end
            end
        end
        clear L U P Q v w H Hr X_H D_H Vr order q_old q_new count
    else
        [L,U,P,Q]=lu(STA);
        v=q_old;
        H=zeros(m,m);
        for j=1:m
            w=Q*(U\(L\(P*(F*v(:,j)))));
            for i=1:j
                H(i,j)=w'*v(:,i);
                w=w-H(i,j)*v(:,i);
            end
            H(j+1,j)=norm(w,2);
            if H(j+1,j)==0
                break
            else
                v(:,j+1)=w/norm(w,2);
            end
        end
        Hr=H(1:(size(H,1)-1),1:size(H,2));
        [X_H,D_H]=eig(Hr);
        D_H=diag(D_H);
        Vr=v(1:size(v,1),1:size(v,2)-1);
        X_est=Vr*X_H;
        [lambda_est,order]=sort(D_H,'descend');
        X_est=X_est(:,order);
        for i=1:(m-1)
            X_est(:,i)=X_est(:,i)/norm(X_est(:,i),2);
        end
        clear L U P Q v w H Hr X_H D_H Vr order q_old
        cwd=pwd;
        cd(cwd);
        save TMP
        clear variables
        load TMP
        delete TMP.mat
        clear cwd
        res_X=zeros(neigs,1);
        lambda=zeros(neigs,1);
        X=zeros(2*SHIFT_XYZ,neigs);
        for k=1:neigs
            count=0;
            res_X(k)=Inf;
            fprintf('\nComputing eigenmode # %d\n',k);
            q_old=X_est(:,k);
            k_est=lambda_est(k);
            [L,U,P,Q]=lu(STA-1/k_est*F);
            while res_X(k)>conv_POW
                q_new=Q*(U\(L\(P*F*q_old)));
                q_old=q_new/norm(q_new,2);
                lambda(k)=((q_old'*Q*(U\(L\(P*(F*q_old)))))^-1+1/k_est)^-1;
                X(:,k)=q_old;
                res_X(k)=norm((lambda(k)*STA*X(:,k)-F*X(:,k)),2);
                if res_X(k)>conv_POW
                    count=count+1;
                    if floor(count/10)==count/10
                        fprintf('Residuals after iteration # %d: %1.3e\n',count,res_X(k));
                    end
                else
                    fprintf('Residuals on the converged eigenvector: %1.3e\n',res_X(k));
                end
                if count>n_iter
                    fprintf('\nERROR: The power iteration method (POW) is not converging as expected.\nEither change the parameters of the power iteration method in your SETTINGS.m file\n')
                    fprintf('or switch to the Explicitely Restarted Arnoldi Method (ERAM)\n(setting EIG_MET=1 in your SETTINGS.m file).\n')
                    if BYP==0
                        fprintf('EXECUTION TERMINATED\n')
                        return
                    else
                        fprintf('EXECUTION ALLOWED TO CONTINUE\n')
                        break
                    end
                end
            end
            clear L U P Q q_old q_new count
        end
        clear X_est lambda_est
    end

    for k=1:neigs
        X(:,k)=X(:,k)/max(abs(X(:,k)));
    end

    cwd=pwd;
    cd(cwd);
    save TMP
    clear variables
    load TMP
    delete TMP.mat
    clear cwd
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATING THE EIGENMODES OF THE SOURCE-FREE PROBLEM (ADJOINT PROBLEM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if RUN_CORESIM
    fprintf('\nCALCULATION OF THE EIGENMODES OF THE SOURCE-FREE PROBLEM (ADJOINT PROBLEM) IN PROGRESS.\n')
    Dp_FLX_adj_tmp=zeros(2*SHIFT_XYZ,3);
    for I=1:I_MAX
    	for J=1:J_MAX
            for K=1:K_MAX
                if (TYP(I,J,K)~=0)
                    % Group 1->1
                    Dp_FLX_adj_tmp(CONV(I,J,K),1)=-ABS1(I,J,K)-REM(I,J,K);
                    % Group 2->1
                    Dp_FLX_adj_tmp(CONV(I,J,K)+SHIFT_XYZ,2)=+REM(I,J,K);
                    % Group 1->2
                    Dp_FLX_adj_tmp(CONV(I,J,K),3)=0;
                    % Group 2->2
                    Dp_FLX_adj_tmp(CONV(I,J,K)+SHIFT_XYZ,1)=-ABS2(I,J,K);
                end
            end
    	end
    end
    Dp_FLX_adj=spdiags(Dp_FLX_adj_tmp,[0;SHIFT_XYZ;-SHIFT_XYZ],2*SHIFT_XYZ,2*SHIFT_XYZ);
    STA_adj=Dp_FLX_adj-AX/DX-AY/DY-AZ/DZ-BX/DX-BY/DY-BZ/DZ-CX/DX-CY/DY-CZ/DZ;
    clear Dp_FLX_adj Dp_FLX_adj_tmp

    F_adj_tmp=zeros(2*SHIFT_XYZ,2);
    for I=1:I_MAX
    	for J=1:J_MAX
            for K=1:K_MAX
                if TYP(I,J,K)~=0
                    % Group 1->1
                    F_adj_tmp(CONV(I,J,K),1)=-NUFIS1(I,J,K);
                    F_adj_tmp(CONV(I,J,K)+SHIFT_XYZ,1)=0;
                    % Group 1->2
                    F_adj_tmp(CONV(I,J,K),2)=-NUFIS2(I,J,K);
                end
            end
    	end
    end
    F_adj=spdiags(F_adj_tmp,[0;-SHIFT_XYZ],2*SHIFT_XYZ,2*SHIFT_XYZ);
    clear F_adj_tmp

    count=0;
    q_old=ones(2*SHIFT_XYZ,1);
    q_old=q_old/norm(q_old);

    if EIG_MET==1
        res=Inf;
        [L,U,P,Q]=lu(STA_adj);
        while res>conv_ERAM
            v=q_old;
            H=zeros(m,m);
            for j=1:m
                w=Q*(U\(L\(P*(F_adj*v(:,j)))));
                for i=1:j
                    H(i,j)=w'*v(:,i);
                    w=w-H(i,j)*v(:,i);
                end
                H(j+1,j)=norm(w,2);
                if H(j+1,j)==0
                    break
                else
                    v(:,j+1)=w/norm(w,2);
                end
            end
            Hr=H(1:(size(H,1)-1),1:size(H,2));
            [X_H,D_H]=eig(Hr);
            D_H=diag(D_H);
            Vr=v(1:size(v,1),1:size(v,2)-1);
            X_adj=Vr*X_H;
            [lambda_adj,order]=sort(D_H,'descend');
            X_adj=X_adj(:,order);
            for i=1:(m-1)
                X_adj(:,i)=X_adj(:,i)/norm(X_adj(:,i),2);
            end
            res_X_adj=zeros(neigs,1);
            for i=1:neigs
                res_X_adj(i)=norm((lambda_adj(i)*STA_adj*X_adj(:,i)-F_adj*X_adj(:,i)),2);
            end
            res=norm(res_X_adj,inf);
            if res>conv_ERAM
                count=count+1;
                fprintf('Residuals before restart # %d: %1.3e\n',count,res);
            else
                fprintf('Residuals on the converged eigenvectors: %1.3e\n',res);
            end
            q_new=zeros(size(X_adj,1),1);
            for i=1:neigs
                q_new=q_new+X_adj(:,i).*1;
            end
            q_old=q_new;
            if count>n_restart
                fprintf('\nERROR: The Explicitely Restarted Arnoldi Method (ERAM) is not converging as expected.\nEither change the parameters of the Explicitely Restarted Arnoldi Method in your SETTINGS.m file\n')
                fprintf('or switch to the power iteration method (POW)\n(setting EIG_MET=2 in your SETTINGS.m file).\n')
                if BYP==0
                    fprintf('EXECUTION TERMINATED\n')
                    return
                else
                    fprintf('EXECUTION ALLOWED TO CONTINUE\n')
                    break
                end
            end
        end
        clear L U P Q v w H Hr X_H D_H Vr order q_old q_new count
    else
        [L,U,P,Q]=lu(STA_adj);
        v=q_old;
        H=zeros(m,m);
        for j=1:m
            w=Q*(U\(L\(P*(F_adj*v(:,j)))));
            for i=1:j
                H(i,j)=w'*v(:,i);
                w=w-H(i,j)*v(:,i);
            end
            H(j+1,j)=norm(w,2);
            if H(j+1,j)==0
                break
            else
                v(:,j+1)=w/norm(w,2);
            end
        end
        Hr=H(1:(size(H,1)-1),1:size(H,2));
        [X_H,D_H]=eig(Hr);
        D_H=diag(D_H);
        Vr=v(1:size(v,1),1:size(v,2)-1);
        X_est=Vr*X_H;
        [lambda_est,order]=sort(D_H,'descend');
        X_est=X_est(:,order);
        for i=1:(m-1)
            X_est(:,i)=X_est(:,i)/norm(X_est(:,i),2);
        end
        clear L U P Q v w H Hr X_H D_H Vr order q_old
        res_X_adj=zeros(neigs,1);
        lambda_adj=zeros(neigs,1);
        X_adj=zeros(2*SHIFT_XYZ,1);
        for k=1:neigs
            count=0;
            res_X_adj(k)=Inf;
            fprintf('\nComputing eigenmode # %d\n',k);
            q_old=X_est(:,k);
            k_est=lambda_est(k);
            [L,U,P,Q]=lu(STA_adj-1/k_est*F_adj);
            while res_X_adj(k)>conv_POW
                q_new=Q*(U\(L\(P*F_adj*q_old)));
                q_old=q_new/norm(q_new,2);
                lambda_adj(k)=((q_old'*Q*(U\(L\(P*(F_adj*q_old)))))^-1+1/k_est)^-1;
                X_adj(:,k)=q_old;
                res_X_adj(k)=norm((lambda_adj(k)*STA_adj*X_adj(:,k)-F_adj*X_adj(:,k)),2);
                if res_X_adj(k)>conv_POW
                    count=count+1;
                    if floor(count/10)==count/10
                        fprintf('Residuals after iteration # %d: %1.3e\n',count,res_X_adj(k));
                    end
                else
                    fprintf('Residuals on the converged eigenvector: %1.3e\n',res_X_adj(k));
                end
                if count>n_iter
                    fprintf('\nERROR: The power iteration method (POW) is not converging as expected.\nEither change the parameters of the power iteration method in your SETTINGS.m file\n')
                    fprintf('or switch to the Explicitely Restarted Arnoldi Method (ERAM)\n(setting EIG_MET=1 in your SETTINGS.m file).\n')
                    if BYP==0
                        fprintf('EXECUTION TERMINATED\n')
                        return
                    else
                        fprintf('EXECUTION ALLOWED TO CONTINUE\n')
                        break
                    end
                end
            end
            clear L U P Q q_old q_new count
        end
        clear X_est lambda_est
    end

    for k=1:neigs
        X_adj(:,k)=X_adj(:,k)/max(abs(X_adj(:,k)));
    end

    cwd=pwd;
    cd(cwd);
    save TMP
    clear variables
    load TMP
    delete TMP.mat
    clear cwd

    keff=lambda(1);
    FLX=X(:,1);
    fprintf('\nSAVING OF THE CALCULATED RESULTS IN PROGRESS.\n')
    MOD1=zeros(I_MAX,J_MAX,K_MAX,neigs);
    MOD2=zeros(I_MAX,J_MAX,K_MAX,neigs);
    MOD1_adj=zeros(I_MAX,J_MAX,K_MAX,neigs);
    MOD2_adj=zeros(I_MAX,J_MAX,K_MAX,neigs);
    FLX1=zeros(I_MAX,J_MAX,K_MAX,neigs);
    FLX2=zeros(I_MAX,J_MAX,K_MAX,neigs);
    FLX1_adj=zeros(I_MAX,J_MAX,K_MAX,neigs);
    FLX2_adj=zeros(I_MAX,J_MAX,K_MAX,neigs);
    EXT_S=0;
    for I=1:I_MAX
        for J=1:J_MAX
            for K=1:K_MAX
                if TYP(I,J,K)~=0
                    for n=1:neigs
                        MOD1(I,J,K,n)=X(CONV(I,J,K),n);
                        MOD2(I,J,K,n)=X(CONV(I,J,K)+SHIFT_XYZ,n);
                        MOD1_adj(I,J,K,n)=X_adj(CONV(I,J,K),n);
                        MOD2_adj(I,J,K,n)=X_adj(CONV(I,J,K)+SHIFT_XYZ,n);
                        FLX1(I,J,K)=abs(MOD1(I,J,K,1));
                        FLX2(I,J,K)=abs(MOD2(I,J,K,1));
                    end
                end
            end
        end
    end
    lambda_tmp=lambda;
    clear lambda
    lambda=lambda_tmp(1:neigs);
    clear lambda_tmp
    save("input\XEROM_input.mat")
    clearvars
end
%% End of CORE SIM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculating the coupling coefficients between nodes according to the
% box-scheme (finite differences) for the feedback data

% Loading input variables:
% D1 D2 REM ABS1 ABS2 NUFIS1 NUFIS2 DX DY DZ
files={
    'XS_data_100_50.mat'
    'GEOM_data_100_50.mat'
    'additional_data_100_50.mat'
    'RESULTS.mat'
    };
% for i=1:size(files,1)
%     if exist(sprintf('feedback/%s',files{i}),'file')==2
%         load (sprintf('feedback/%s',files{i}))
%         fprintf("feedback data /%s \n",files{i})
%     else
%         fprintf('\nERROR: Missing input file or missing input directory.\nEXECUTION TERMINATED\n')
%         return
%     end
% end
for i=1:size(files,1)
    if exist(sprintf('../For_Christophe/feedback_100_50/%s',files{i}),'file')==2
        load (sprintf('../For_Christophe/feedback_100_50/%s',files{i}))
        fprintf("feedback data /%s \n",files{i})
    else
        fprintf('\nERROR: Missing input file or missing input directory.\nEXECUTION TERMINATED\n')
        return
    end
end


VAR={'D1','D2','REM','ABS1','ABS2','NUFIS1','NUFIS2','DX','DY','DZ','STA_FLX1','STA_FLX2','KF1','KF2'};
TEST_VAR=ismember((VAR),who);

for i=1:size(VAR,2)
    if TEST_VAR(1,i)==0RERUN
        fprintf('\nERROR: The variable %s is missing from your input files.\n',VAR{i});
        RUN_TEST=0;
    end
end

% if RUN_TEST==0
%     fprintf('\nEXECUTION TERMINATED\n')
%     return
% end
clear VAR TEST_VAR files

% Checking the geometry of the core
I_MAX=size(D2,1);
J_MAX=size(D2,2);
K_MAX=size(D2,3);

if I_MAX<3
    fprintf('\nERROR: You need at least 3 rows to properly define this 3-D problem.\nEXECUTION TERMINATED\n')
    return
end

if J_MAX<3
    fprintf('\nERROR: You need at least 3 columns to properly define this 3-D problem.\nEXECUTION TERMINATED\n')
    return
end

if K_MAX<3
    fprintf('\nERROR: You need at least 3 axial planes to properly define this 3-D problem.\nEXECUTION TERMINATED\n')
    return
end

% Creating a variable TYP defining the location of the fuel and reflector nodes
% (TYP=1 for fuel nodes and TYP=-1 for reflector nodes)
TYP=zeros(I_MAX,J_MAX,K_MAX);
for I=1:I_MAX
    for J=1:J_MAX
        for K=1:K_MAX
            if D2(I,J,K)==0
                TYP(I,J,K)=0;
            elseif (D2(I,J,K)~=0 && NUFIS2(I,J,K)==0)
                TYP(I,J,K)=-1;
            elseif (D2(I,J,K)~=0 && NUFIS2(I,J,K)~=0)
                TYP(I,J,K)=1;
            end
        end
    end
end




% if RUN_TEST==0
%     fprintf('\nEXECUTION TERMINATED\n')
%     return
% end

% Creating a variable allowing reordering the 3-D variables into column vectors
tmp=0;
CONV=zeros(I_MAX,J_MAX,K_MAX);
SHIFT_XY=zeros(K_MAX,1);
for K=1:K_MAX
    for I=1:I_MAX
        for J=1:J_MAX
            if TYP(I,J,K)~=0
                tmp=tmp+1;
                CONV(I,J,K)=tmp;
            else
                CONV(I,J,K)=0;
            end
        end
    end
    if K>1
        SHIFT_XY(K)=max(max(CONV(:,:,K)))-max(max(CONV(:,:,K-1)));
    else
        SHIFT_XY(K)=max(max(CONV(:,:,K)));
    end
end
clear tmp
SHIFT_XY_test=SHIFT_XY(1);
for K=2:K_MAX
    if SHIFT_XY(K)~=SHIFT_XY_test
        fprintf('\nERROR: Each of the axial planes should contain the same number of nodes.\nEXECUTION TERMINATED\n')
        return
    end
end
clear SHIFT_XY
SHIFT_XY=SHIFT_XY_test;
clear SHIFT_XY_test
SHIFT_XYZ=max(max(max(CONV)));

% Creating variables allowing determining the positions of consecutive
% nodes in the matrices/vectors for a positive and negative, respectively,
% increment of the row index (IP_SHIFT and IM_SHIFT, respectively), column
% index (JP_SHIFT and JM_SHIFT, respectively), and axial plane index
% (KP_SHIFT and KM_SHIFT, respectively)
IP_SHIFT=zeros(I_MAX,J_MAX,K_MAX);
IM_SHIFT=zeros(I_MAX,J_MAX,K_MAX);
JP_SHIFT=zeros(I_MAX,J_MAX,K_MAX);
JM_SHIFT=zeros(I_MAX,J_MAX,K_MAX);
KP_SHIFT=zeros(I_MAX,J_MAX,K_MAX);
KM_SHIFT=zeros(I_MAX,J_MAX,K_MAX);
for I=1:I_MAX
    for J=1:J_MAX
        for K=1:K_MAX
            if TYP(I,J,K)~=0
                if ((I>1) && (TYP(I-1,J,K)~=0))
                    IM_SHIFT(I,J,K)=CONV(I-1,J,K)-CONV(I,J,K);
                end
                if ((I<I_MAX) && (TYP(I+1,J,K)~=0))
                    IP_SHIFT(I,J,K)=CONV(I+1,J,K)-CONV(I,J,K);
                end
                if ((J>1) && (TYP(I,J-1,K)~=0))
                    JM_SHIFT(I,J,K)=CONV(I,J-1,K)-CONV(I,J,K);
                end
                if ((J<J_MAX) && (TYP(I,J+1,K)~=0))
                    JP_SHIFT(I,J,K)=CONV(I,J+1,K)-CONV(I,J,K);
                end
                if ((K>1) && (TYP(I,J,K-1)~=0))
                    KM_SHIFT(I,J,K)=CONV(I,J,K-1)-CONV(I,J,K);
                end
                if ((K<K_MAX) && (TYP(I,J,K+1)~=0))
                    KP_SHIFT(I,J,K)=CONV(I,J,K+1)-CONV(I,J,K);
                end
            end
        end
    end
end

% Calculating the coupling coefficients between nodes according to the
% box-scheme (finite differences)
SHIFT=0;
for K=1:K_MAX
    % AX
    ax11_tmp=zeros(SHIFT_XY,1);
    ax22_tmp=zeros(SHIFT_XY,1);
    for I=1:I_MAX
        for J=1:J_MAX
            if (TYP(I,J,K)~=0)
                if ( ((I>1) && (TYP(I,J,K)~=0) && (TYP((I-1),J,K)==0)) || ((I==1) && (TYP(I,J,K)~=0)) )
                    % Group 1->1
                    ax11_tmp(CONV(I,J,K)-SHIFT,1)=2*D1((I+1),J,K)*D1(I,J,K)/(DX*(D1((I+1),J,K)+D1(I,J,K)))+1/2/(1+DX/(4*D1(I,J,K)));
                    % Group 2->2
                    ax22_tmp(CONV(I,J,K)-SHIFT,1)=2*D2((I+1),J,K)*D2(I,J,K)/(DX*(D2((I+1),J,K)+D2(I,J,K)))+1/2/(1+DX/(4*D2(I,J,K)));
                elseif ( ((I<I_MAX) && (TYP(I,J,K)~=0) && (TYP((I+1),J,K)==0)) || ((I==I_MAX) && (TYP(I,J,K)~=0)) )
                    % Group 1->1
                    ax11_tmp(CONV(I,J,K)-SHIFT,1)=2*D1(I,J,K)*D1((I-1),J,K)/(DX*(D1(I,J,K)+D1((I-1),J,K)))+1/2/(1+DX/(4*D1(I,J,K)));
                    % Group 2->2
                    ax22_tmp(CONV(I,J,K)-SHIFT,1)=2*D2(I,J,K)*D2((I-1),J,K)/(DX*(D2(I,J,K)+D2((I-1),J,K)))+1/2/(1+DX/(4*D2(I,J,K)));
                else
                    % Group 1->1
                    ax11_tmp(CONV(I,J,K)-SHIFT,1)=2*D1((I+1),J,K)*D1(I,J,K)/(DX*(D1((I+1),J,K)+D1(I,J,K)))+2*D1(I,J,K)*D1((I-1),J,K)/(DX*(D1(I,J,K)+D1((I-1),J,K)));
                    % Group 2->2
                    ax22_tmp(CONV(I,J,K)-SHIFT,1)=2*D2((I+1),J,K)*D2(I,J,K)/(DX*(D2((I+1),J,K)+D2(I,J,K)))+2*D2(I,J,K)*D2((I-1),J,K)/(DX*(D2(I,J,K)+D2((I-1),J,K)));
                end
            end
        end
    end
    ax11=spdiags(ax11_tmp,0,SHIFT_XY,SHIFT_XY);
    ax22=spdiags(ax22_tmp,0,SHIFT_XY,SHIFT_XY);
    AX(SHIFT+1:SHIFT+SHIFT_XY,SHIFT+1:SHIFT+SHIFT_XY)=ax11;
    AX(SHIFT_XYZ+SHIFT+1:SHIFT_XYZ+SHIFT+SHIFT_XY,SHIFT_XYZ+SHIFT+1:SHIFT_XYZ+SHIFT+SHIFT_XY)=ax22;
    clear ax11 ax22 ax11_tmp ax22_tmp
 
    
    %%
    % BX
    bx11_tmp=zeros(max(max(max(CONV)))+max(max(abs(IP_SHIFT(:,:,K))))-SHIFT,max(max(abs(IP_SHIFT(:,:,K)))));
    bx22_tmp=zeros(max(max(max(CONV)))+max(max(abs(IP_SHIFT(:,:,K))))-SHIFT,max(max(abs(IP_SHIFT(:,:,K)))));
    for C=1:max(max(abs(IP_SHIFT(:,:,K))))
        for I=1:I_MAX
            for J=1:J_MAX
                if (IP_SHIFT(I,J,K)==C)
                    if ( ((I>1) && (TYP(I,J,K)~=0) && (TYP((I-1),J,K)==0)) || ((I==1) && (TYP(I,J,K)~=0)) )
                        % Group 1->1
                        bx11_tmp(CONV(I,J,K)+C-SHIFT,C)=-2*D1((I+1),J,K)*D1(I,J,K)/(DX*(D1((I+1),J,K)+D1(I,J,K)));
                        % Group 2->2
                        bx22_tmp(CONV(I,J,K)+C-SHIFT,C)=-2*D2((I+1),J,K)*D2(I,J,K)/(DX*(D2((I+1),J,K)+D2(I,J,K)));
                    elseif ( ((I<I_MAX) && (TYP(I,J,K)~=0) && (TYP((I+1),J,K)==0)) || ((I==I_MAX) && (TYP(I,J,K)~=0)) )
                    else
                        % Group 1->1
                        bx11_tmp(CONV(I,J,K)+C-SHIFT,C)=-2*D1((I+1),J,K)*D1(I,J,K)/(DX*(D1((I+1),J,K)+D1(I,J,K)));
                        % Group 2->2
                        bx22_tmp(CONV(I,J,K)+C-SHIFT,C)=-2*D2((I+1),J,K)*D2(I,J,K)/(DX*(D2((I+1),J,K)+D2(I,J,K)));
                    end
                end
            end
        end
    end

    if size(bx11_tmp,1)<SHIFT_XY
        bx11_tmp(SHIFT_XY,1)=0;
    end
    if size(bx22_tmp,1)<SHIFT_XY
        bx22_tmp(SHIFT_XY,1)=0;
    end
    d=zeros(max(max(abs(IP_SHIFT(:,:,K)))),1);
    for C=1:max(max(abs(IP_SHIFT(:,:,K))))
        d(C,1)=C;
    end
    bx11=spdiags(bx11_tmp,d,SHIFT_XY,SHIFT_XY);
    bx22=spdiags(bx22_tmp,d,SHIFT_XY,SHIFT_XY);    
    BX(SHIFT+1:SHIFT+SHIFT_XY,SHIFT+1:SHIFT+SHIFT_XY)=bx11;
    BX(SHIFT_XYZ+SHIFT+1:SHIFT_XYZ+SHIFT+SHIFT_XY,SHIFT_XYZ+SHIFT+1:SHIFT_XYZ+SHIFT+SHIFT_XY)=bx22;
    clear bx11 bx22 bx11_tmp bx22_tmp C d
    
    % CX
    cx11_tmp=zeros(max(max(max(CONV)))-max(max(abs(IM_SHIFT(:,:,K))))-SHIFT,max(max(abs(IP_SHIFT(:,:,K)))));
    cx22_tmp=zeros(max(max(max(CONV)))-max(max(abs(IM_SHIFT(:,:,K))))-SHIFT,max(max(abs(IP_SHIFT(:,:,K)))));
    for C=1:max(max(abs(IM_SHIFT(:,:,K))))
        for I=1:I_MAX
            for J=1:J_MAX
                if (IM_SHIFT(I,J)==-C)
                    if ( ((I>1) && (TYP(I,J,K)~=0) && (TYP((I-1),J,K)==0)) || ((I==1) && (TYP(I,J,K)~=0)) )
                    elseif ( ((I<I_MAX) && (TYP(I,J,K)~=0) && (TYP((I+1),J,K)==0)) || ((I==I_MAX) && (TYP(I,J,K)~=0)) )
                        % Group 1->1
                        cx11_tmp(CONV(I,J,K)-C-SHIFT,C)=-2*D1(I,J,K)*D1((I-1),J,K)/(DX*(D1(I,J,K)+D1((I-1),J,K)));
                        % Group 2->2
                        cx22_tmp(CONV(I,J,K)-C-SHIFT,C)=-2*D2(I,J,K)*D2((I-1),J,K)/(DX*(D2(I,J,K)+D2((I-1),J,K)));
                    else
                        % Group 1->1
                        cx11_tmp(CONV(I,J,K)-C-SHIFT,C)=-2*D1(I,J,K)*D1((I-1),J,K)/(DX*(D1(I,J,K)+D1((I-1),J,K)));
                        % Group 2->2
                        cx22_tmp(CONV(I,J,K)-C-SHIFT,C)=-2*D2(I,J,K)*D2((I-1),J,K)/(DX*(D2(I,J,K)+D2((I-1),J,K)));
                    end
                end
            end
        end
    end

    if size(cx11_tmp,1)<SHIFT_XY
        cx11_tmp(SHIFT_XY,1)=0;
    end
    if size(cx22_tmp,1)<SHIFT_XY
        cx22_tmp(SHIFT_XY,1)=0;
    end
    d=zeros(max(max(abs(IM_SHIFT(:,:,K)))),1);
    for C=1:max(max(abs(IM_SHIFT(:,:,K))))
        d(C,1)=-C;
    end
    cx11=spdiags(cx11_tmp,d,SHIFT_XY,SHIFT_XY);
    cx22=spdiags(cx22_tmp,d,SHIFT_XY,SHIFT_XY);
    CX(SHIFT+1:SHIFT+SHIFT_XY,SHIFT+1:SHIFT+SHIFT_XY)=cx11;
    CX(SHIFT_XYZ+SHIFT+1:SHIFT_XYZ+SHIFT+SHIFT_XY,SHIFT_XYZ+SHIFT+1:SHIFT_XYZ+SHIFT+SHIFT_XY)=cx22;
    clear cx11 cx22 cx11_tmp cx22_tmp C d
    
    % AY
    ay11_tmp=zeros(SHIFT_XY,1);
    ay22_tmp=zeros(SHIFT_XY,1);
    for I=1:I_MAX
        for J=1:J_MAX
            if (TYP(I,J,K)~=0)
                if ( ((J>1) && (TYP(I,J,K)~=0) && (TYP(I,(J-1),K)==0)) || ((J==1) && (TYP(I,J,K)~=0)) )
                    % Group 1->1
                    ay11_tmp(CONV(I,J,K)-SHIFT,1)=2*D1(I,(J+1),K)*D1(I,J,K)/(DY*(D1(I,(J+1),K)+D1(I,J,K)))+1/2/(1+DY/(4*D1(I,J,K)));
                    % Group 2->2
                    ay22_tmp(CONV(I,J,K)-SHIFT,1)=2*D2(I,(J+1),K)*D2(I,J,K)/(DY*(D2(I,(J+1),K)+D2(I,J,K)))+1/2/(1+DY/(4*D2(I,J,K)));
                elseif ( ((J<J_MAX) && (TYP(I,J,K)~=0) && (TYP(I,(J+1),K)==0)) || ((J==J_MAX) && (TYP(I,J,K)~=0)) )
                    % Group 1->1
                    ay11_tmp(CONV(I,J,K)-SHIFT,1)=2*D1(I,J,K)*D1(I,(J-1),K)/(DY*(D1(I,J,K)+D1(I,(J-1),K)))+1/2/(1+DY/(4*D1(I,J,K)));
                    % Group 2->2
                    ay22_tmp(CONV(I,J,K)-SHIFT,1)=2*D2(I,J,K)*D2(I,(J-1),K)/(DY*(D2(I,J,K)+D2(I,(J-1),K)))+1/2/(1+DY/(4*D2(I,J,K)));
                else
                    % Group 1->1
                    ay11_tmp(CONV(I,J,K)-SHIFT,1)=2*D1(I,(J+1),K)*D1(I,J,K)/(DY*(D1(I,(J+1),K)+D1(I,J,K)))+2*D1(I,J,K)*D1(I,(J-1),K)/(DY*(D1(I,J,K)+D1(I,(J-1),K)));
                    % Group 2->2
                    ay22_tmp(CONV(I,J,K)-SHIFT,1)=2*D2(I,(J+1),K)*D2(I,J,K)/(DY*(D2(I,(J+1),K)+D2(I,J,K)))+2*D2(I,J,K)*D2(I,(J-1),K)/(DY*(D2(I,J,K)+D2(I,(J-1),K)));
                end
            end
        end
    end
    ay11=spdiags(ay11_tmp,0,SHIFT_XY,SHIFT_XY);
    ay22=spdiags(ay22_tmp,0,SHIFT_XY,SHIFT_XY);
    AY(SHIFT+1:SHIFT+SHIFT_XY,SHIFT+1:SHIFT+SHIFT_XY)=ay11;
    AY(SHIFT_XYZ+SHIFT+1:SHIFT_XYZ+SHIFT+SHIFT_XY,SHIFT_XYZ+SHIFT+1:SHIFT_XYZ+SHIFT+SHIFT_XY)=ay22;
    clear ay11 ay22 ay11_tmp ay22_tmp

    % BY
    by11_tmp=zeros(SHIFT_XY+1,1);
    by22_tmp=zeros(SHIFT_XY+1,1);
    for I=1:I_MAX
        for J=1:J_MAX
            if (TYP(I,J,K)~=0)
                if ( ((J>1) && (TYP(I,J,K)~=0) && (TYP(I,(J-1),K)==0)) || ((J==1) && (TYP(I,J,K)~=0)) )
                    % Group 1->1
                    by11_tmp(CONV(I,J,K)+1-SHIFT,1)=-2*D1(I,(J+1),K)*D1(I,J,K)/(DY*(D1(I,(J+1),K)+D1(I,J,K)));
                    % Group 2->2
                    by22_tmp(CONV(I,J,K)+1-SHIFT,1)=-2*D2(I,(J+1),K)*D2(I,J,K)/(DY*(D2(I,(J+1),K)+D2(I,J,K)));
                elseif ( ((J<J_MAX) && (TYP(I,J,K)~=0) && (TYP(I,(J+1),K)==0)) || ((J==J_MAX) && (TYP(I,J,K)~=0)) )
                else
                    % Group 1->1
                    by11_tmp(CONV(I,J,K)+1-SHIFT,1)=-2*D1(I,(J+1),K)*D1(I,J,K)/(DY*(D1(I,(J+1),K)+D1(I,J,K)));
                    % Group 2->2
                    by22_tmp(CONV(I,J,K)+1-SHIFT,1)=-2*D2(I,(J+1),K)*D2(I,J,K)/(DY*(D2(I,(J+1),K)+D2(I,J,K)));
                end
            end
        end
    end
    by11=spdiags(by11_tmp,1,SHIFT_XY,SHIFT_XY);
    by22=spdiags(by22_tmp,1,SHIFT_XY,SHIFT_XY);
    BY(SHIFT+1:SHIFT+SHIFT_XY,SHIFT+1:SHIFT+SHIFT_XY)=by11;
    BY(SHIFT_XYZ+SHIFT+1:SHIFT_XYZ+SHIFT+SHIFT_XY,SHIFT_XYZ+SHIFT+1:SHIFT_XYZ+SHIFT+SHIFT_XY)=by22;
    clear by11 by22 by11_tmp by22_tmp
    
    % CY
    cy11_tmp=zeros(SHIFT_XY-1,1);
    cy22_tmp=zeros(SHIFT_XY-1,1);
    for I=1:I_MAX
        for J=1:J_MAX
            if (TYP(I,J,K)~=0)
                if ( ((J>1) && (TYP(I,J,K)~=0) && (TYP(I,(J-1),K)==0)) || ((J==1) && (TYP(I,J,K)~=0)) )
                elseif ( ((J<J_MAX) && (TYP(I,J,K)~=0) && (TYP(I,(J+1),K)==0)) || ((J==J_MAX) && (TYP(I,J,K)~=0)) )
                    % Group 1->1
                    cy11_tmp(CONV(I,J,K)-1-SHIFT,1)=-2*D1(I,J,K)*D1(I,(J-1),K)/(DY*(D1(I,J,K)+D1(I,(J-1),K)));
                    % Group 2->2
                    cy22_tmp(CONV(I,J,K)-1-SHIFT,1)=-2*D2(I,J,K)*D2(I,(J-1),K)/(DY*(D2(I,J,K)+D2(I,(J-1),K)));
                else
                    % Group 1->1
                    cy11_tmp(CONV(I,J,K)-1-SHIFT,1)=-2*D1(I,J,K)*D1(I,(J-1),K)/(DY*(D1(I,J,K)+D1(I,(J-1),K)));
                    % Group 2->2
                    cy22_tmp(CONV(I,J,K)-1-SHIFT,1)=-2*D2(I,J,K)*D2(I,(J-1),K)/(DY*(D2(I,J,K)+D2(I,(J-1),K)));
                end
            end
        end
    end
    cy11=spdiags(cy11_tmp,-1,SHIFT_XY,SHIFT_XY);
    cy22=spdiags(cy22_tmp,-1,SHIFT_XY,SHIFT_XY);
    CY(SHIFT+1:SHIFT+SHIFT_XY,SHIFT+1:SHIFT+SHIFT_XY)=cy11;
    CY(SHIFT_XYZ+SHIFT+1:SHIFT_XYZ+SHIFT+SHIFT_XY,SHIFT_XYZ+SHIFT+1:SHIFT_XYZ+SHIFT+SHIFT_XY)=cy22;
    clear cy11 cy22 cy11_tmp cy22_tmp

    % AZ
    az11_tmp=zeros(SHIFT_XY,1);
    az22_tmp=zeros(SHIFT_XY,1);
    for I=1:I_MAX
        for J=1:J_MAX
            if (TYP(I,J,K)~=0)
                if ( K==1 )
                    % Group 1->1
                    az11_tmp(CONV(I,J,K)-SHIFT,1)=2*D1(I,J,(K+1))*D1(I,J,K)/(DZ*(D1(I,J,(K+1))+D1(I,J,K)))+1/2/(1+DZ/(4*D1(I,J,K)));
                    % Group 2->2
                    az22_tmp(CONV(I,J,K)-SHIFT,1)=2*D2(I,J,(K+1))*D2(I,J,K)/(DZ*(D2(I,J,(K+1))+D2(I,J,K)))+1/2/(1+DZ/(4*D2(I,J,K)));
                elseif ( K==K_MAX )
                    % Group 1->1
                    az11_tmp(CONV(I,J,K)-SHIFT,1)=2*D1(I,J,K)*D1(I,J,(K-1))/(DZ*(D1(I,J,K)+D1(I,J,(K-1))))+1/2/(1+DZ/(4*D1(I,J,K)));
                    % Group 2->2
                    az22_tmp(CONV(I,J,K)-SHIFT,1)=2*D2(I,J,K)*D2(I,J,(K-1))/(DZ*(D2(I,J,K)+D2(I,J,(K-1))))+1/2/(1+DZ/(4*D2(I,J,K)));
                else
                    % Group 1->1
                    az11_tmp(CONV(I,J,K)-SHIFT,1)=2*D1(I,J,(K+1))*D1(I,J,K)/(DZ*(D1(I,J,(K+1))+D1(I,J,K)))+2*D1(I,J,K)*D1(I,J,(K-1))/(DZ*(D1(I,J,K)+D1(I,J,(K-1))));
                    % Group 2->2
                    az22_tmp(CONV(I,J,K)-SHIFT,1)=2*D2(I,J,(K+1))*D2(I,J,K)/(DZ*(D2(I,J,(K+1))+D2(I,J,K)))+2*D2(I,J,K)*D2(I,J,(K-1))/(DZ*(D2(I,J,K)+D2(I,J,(K-1))));
                end
            end
        end
    end
    az11=spdiags(az11_tmp,0,SHIFT_XY,SHIFT_XY);
    az22=spdiags(az22_tmp,0,SHIFT_XY,SHIFT_XY);
    AZ(SHIFT+1:SHIFT+SHIFT_XY,SHIFT+1:SHIFT+SHIFT_XY)=az11;
    AZ(SHIFT_XYZ+SHIFT+1:SHIFT_XYZ+SHIFT+SHIFT_XY,SHIFT_XYZ+SHIFT+1:SHIFT_XYZ+SHIFT+SHIFT_XY)=az22;
    clear az11 az22 az11_tmp az22_tmp
    
    % BZ
    bz11_tmp=zeros(SHIFT_XY,1);
    bz22_tmp=zeros(SHIFT_XY,1);
    for I=1:I_MAX
        for J=1:J_MAX
            if (TYP(I,J,K)~=0)
                if ( K==1 )
                    % Group 1->1
                    bz11_tmp(CONV(I,J,K)-SHIFT,1)=-2*D1(I,J,(K+1))*D1(I,J,K)/(DZ*(D1(I,J,(K+1))+D1(I,J,K)));
                    % Group 2->2
                    bz22_tmp(CONV(I,J,K)-SHIFT,1)=-2*D2(I,J,(K+1))*D2(I,J,K)/(DZ*(D2(I,J,(K+1))+D2(I,J,K)));
                elseif ( K==K_MAX )
                else
                    % Group 1->1
                    bz11_tmp(CONV(I,J,K)-SHIFT,1)=-2*D1(I,J,(K+1))*D1(I,J,K)/(DZ*(D1(I,J,(K+1))+D1(I,J,K)));
                    % Group 2->2
                    bz22_tmp(CONV(I,J,K)-SHIFT,1)=-2*D2(I,J,(K+1))*D2(I,J,K)/(DZ*(D2(I,J,(K+1))+D2(I,J,K)));
                end
            end
        end
    end
    if K<K_MAX
        bz11=spdiags(bz11_tmp,0,SHIFT_XY,SHIFT_XY);
        bz22=spdiags(bz22_tmp,0,SHIFT_XY,SHIFT_XY);
        BZ(SHIFT+1:SHIFT+SHIFT_XY,SHIFT+1+SHIFT_XY:SHIFT+2*SHIFT_XY)=bz11;
        BZ(SHIFT_XYZ+SHIFT+1:SHIFT_XYZ+SHIFT+SHIFT_XY,SHIFT_XYZ+SHIFT+1+SHIFT_XY:SHIFT_XYZ+SHIFT+2*SHIFT_XY)=bz22;
    end
    clear bz11 bz22 bz11_tmp bz22_tmp
    
    % CZ
    cz11_tmp=zeros(SHIFT_XY,1);
    cz22_tmp=zeros(SHIFT_XY,1);
    for I=1:I_MAX
        for J=1:J_MAX
            if (TYP(I,J,K)~=0)
                if ( K==1 )
                elseif ( K==K_MAX )
                    % Group 1->1
                    cz11_tmp(CONV(I,J,K)-SHIFT,1)=-2*D1(I,J,K)*D1(I,J,(K-1))/(DZ*(D1(I,J,K)+D1(I,J,(K-1))));
                    % Group 2->2
                    cz22_tmp(CONV(I,J,K)-SHIFT,1)=-2*D2(I,J,K)*D2(I,J,(K-1))/(DZ*(D2(I,J,K)+D2(I,J,(K-1))));
                else
                    % Group 1->1
                    cz11_tmp(CONV(I,J,K)-SHIFT,1)=-2*D1(I,J,K)*D1(I,J,(K-1))/(DZ*(D1(I,J,K)+D1(I,J,(K-1))));
                    % Group 2->2
                    cz22_tmp(CONV(I,J,K)-SHIFT,1)=-2*D2(I,J,K)*D2(I,J,(K-1))/(DZ*(D2(I,J,K)+D2(I,J,(K-1))));
                end
            end
        end
    end
    if K>1
        cz11=spdiags(cz11_tmp,0,SHIFT_XY,SHIFT_XY);
        cz22=spdiags(cz22_tmp,0,SHIFT_XY,SHIFT_XY);
        CZ(SHIFT+1:SHIFT+SHIFT_XY,SHIFT+1-SHIFT_XY:SHIFT)=cz11;
        CZ(SHIFT_XYZ+SHIFT+1:SHIFT_XYZ+SHIFT+SHIFT_XY,SHIFT_XYZ+SHIFT+1-SHIFT_XY:SHIFT_XYZ+SHIFT)=cz22;
    end
    clear cz11 cz22 cz11_tmp cz22_tmp

    SHIFT=SHIFT+SHIFT_XY;
end

if size(BZ,1)<2*SHIFT_XYZ
    BZ(2*SHIFT_XYZ,2*SHIFT_XYZ)=0;
end
if size(CZ,2)<2*SHIFT_XYZ
    CZ(2*SHIFT_XYZ,2*SHIFT_XYZ)=0;
end

save("temp/temp_feedback_vars.mat","AX","AY","AZ","BX","BY","BZ","CX","CY","CZ","DX","DY","DZ","D1","D2","NUFIS1","NUFIS2","KF1","KF2","ABS1","ABS2","REM","STA_FLX1","STA_FLX2")


%% load data for delta values
clearvars
input_additional=load("../For_Christophe/input_100_50/additional_data_100_50.mat");
input = load("temp/temp_input_vars.mat");
feedback_additional = load("../For_Christophe/feedback_100_50/additional_data_100_50.mat");
feedback = load("temp/temp_feedback_vars.mat");
input_results = load("../For_Christophe/input_100_50/RESULTS.mat");
feedback_results = load("../For_Christophe/feedback_100_50/RESULTS.mat");
DV = input.DX*input.DY*input.DZ;

SHIFT = input.SHIFT;
SHIFT_XYZ = input.SHIFT_XYZ;
SHIFT_XY = input.SHIFT_XY;
I_MAX = input.I_MAX;
J_MAX = input.J_MAX;
K_MAX = input.K_MAX;
TYP = input.TYP;
CONV = input.CONV;
DX = input.DX;
DY = input.DY;
DZ = input.DZ;

sizex = I_MAX; %number of nodes in the x-direction
sizey = J_MAX; %number of nodes in the y-direction
sizez = K_MAX; %number of nodes in the z-direction

input.KFIS1 = input.KF1;
input.KFIS2 = input.KF2;
feedback.KFIS1 = feedback.KF1;
feedback.KFIS2 = feedback.KF2;
input.MOD = [input_results.MOD1;input_results.MOD2]; % vector of solutions to the forward problem
input.MOD_EQ = abs([input_results.MOD1(:,:,:,1);input_results.MOD2(:,:,:,1)]); % Vector of only the equilibrium neutron flux solution
input.KFISINT =  DV*sum(G2_inner_product([input.KFIS1,input.KFIS2],input.MOD_EQ,"vector","vector"),"all");
input.PS = input_additional.reactor_power*input_results.lambda(1)/input.KFISINT;

feedback.MOD = [feedback_results.MOD1;feedback_results.MOD2]; % vector of solutions to the forward problem
feedback.MOD_EQ = abs([feedback_results.MOD1(:,:,:,1);feedback_results.MOD2(:,:,:,1)]); % Vector of only the equilibrium neutron flux solution
feedback.KFISINT =  DV*sum(G2_inner_product([feedback.KFIS1,feedback.KFIS2],feedback.MOD_EQ,"vector","vector"),"all");
feedback.PS = feedback_additional.reactor_power*input_results.lambda(1)/input.KFISINT;

input.MOD_EQ_scaled= input.PS*input.MOD_EQ; % Vector of only the equilibrium neutron flux solution scaled by the power
input.MOD_EQ_1_scaled= input.MOD_EQ_scaled(1:sizey,:,:);
input.MOD_EQ_2_scaled= input.MOD_EQ_scaled(sizey+1:end,:,:);

feedback.MOD_EQ_scaled= feedback.PS*feedback.MOD_EQ; % Vector of only the equilibrium neutron flux solution scaled by the power
feedback.MOD_EQ_1_scaled= feedback.MOD_EQ_scaled(1:sizey,:,:);
feedback.MOD_EQ_2_scaled= feedback.MOD_EQ_scaled(sizey+1:end,:,:);

input.power = DV * sum(input.KF1.*input.MOD_EQ_1_scaled+input.KF2.*input.MOD_EQ_2_scaled,'all');
feedback.power = DV * sum(feedback.KF1.*feedback.MOD_EQ_1_scaled+feedback.KF2.*feedback.MOD_EQ_2_scaled,'all');


STA_FLX_col=zeros(1,SHIFT_XYZ*2);

for I=1:I_MAX
	for J=1:J_MAX
        for K=1:K_MAX
            if (TYP(I,J,K)~=0)
                % Group 1
                input.STA_FLX_col(CONV(I,J,K)) = input.MOD_EQ_1_scaled(I,J,K);
                feedback.STA_FLX_col(CONV(I,J,K)) = feedback.MOD_EQ_1_scaled(I,J,K);
                % Group 2
                input.STA_FLX_col(CONV(I,J,K)+SHIFT_XYZ) = input.MOD_EQ_2_scaled(I,J,K);
                feedback.STA_FLX_col(CONV(I,J,K)+SHIFT_XYZ) = feedback.MOD_EQ_2_scaled(I,J,K);
            end
        end
	end
end

input.NapJ = (input.AX/input.DX+input.AY/input.DY+input.AZ/input.DZ+input.BX/input.DX+input.BY/input.DY+input.BZ/input.DZ+input.CX/input.DX+input.CY/input.DY+input.CZ/input.DZ) * input.STA_FLX_col';

feedback.NapJ = (feedback.AX/feedback.DX+feedback.AY/feedback.DY+feedback.AZ/feedback.DZ+feedback.BX/feedback.DX+feedback.BY/feedback.DY+feedback.BZ/feedback.DZ+feedback.CX/feedback.DX+feedback.CY/feedback.DY+feedback.CZ/feedback.DZ) * feedback.STA_FLX_col';

for I = 1:I_MAX
    for J = 1:J_MAX
        for K = 1:K_MAX
            input.NJ1(I,J,K) = input.NapJ(CONV(I,J,K)+SHIFT_XY);
            feedback.NJ1(I,J,K) = feedback.NapJ(CONV(I,J,K)+SHIFT_XY);
            input.NJ2(I,J,K) = input.NapJ(CONV(I,J,K)+SHIFT_XYZ);
            feedback.NJ2(I,J,K) = feedback.NapJ(CONV(I,J,K)+SHIFT_XYZ);
        end
    end
end

input.NUFIS1PHI1 = input.NUFIS1.*input.MOD_EQ_1_scaled; 
input.NUFIS2PHI2 = input.NUFIS2.*input.MOD_EQ_2_scaled;
input.ABS1PHI1 = input.ABS1.*input.MOD_EQ_1_scaled;
input.ABS2PHI2 = input.ABS2.*input.MOD_EQ_2_scaled;
input.REMPHI1 = input.REM.*input.MOD_EQ_1_scaled;


feedback.NUFIS1PHI1 = feedback.NUFIS1.*feedback.MOD_EQ_1_scaled;
feedback.NUFIS2PHI2 = feedback.NUFIS2.*feedback.MOD_EQ_2_scaled;
feedback.ABS1PHI1 = feedback.ABS1.*feedback.MOD_EQ_1_scaled;
feedback.ABS2PHI2 = feedback.ABS2.*feedback.MOD_EQ_2_scaled;
feedback.REMPHI1 = feedback.REM.*feedback.MOD_EQ_1_scaled;

DNapJ = feedback.NapJ - input.NapJ;
DNUFIS1PHI1 = feedback.NUFIS1PHI1 - input.NUFIS1PHI1;
DNUFIS2PHI2 = feedback.NUFIS2PHI2 - input.NUFIS2PHI2;  
DABS1PHI1 = feedback.ABS1PHI1 - input.ABS1PHI1;
DABS2PHI2 = feedback.ABS2PHI2 - input.ABS2PHI2;
DREMPHI1 = feedback.REMPHI1 - input.REMPHI1;
DFLX1 = feedback.MOD_EQ_1_scaled - input.MOD_EQ_1_scaled;
DFLX2 = feedback.MOD_EQ_2_scaled - input.MOD_EQ_2_scaled;

%Extract only thermal information and convert back to 3D array
DNJ1 = zeros(I_MAX,J_MAX,K_MAX);
DNJ2 = zeros(I_MAX,J_MAX,K_MAX);

for I = 1:I_MAX
    for J = 1:J_MAX
        for K = 1:K_MAX
            DNJ1(I,J,K) = DNapJ(CONV(I,J,K)+SHIFT_XY);
            DNJ2(I,J,K) = DNapJ(CONV(I,J,K)+SHIFT_XYZ);
        end
    end
end

DABS1PHI1DPHI1 = DABS1PHI1 ./DFLX1;
DNUFIS1PHI1DPHI1 = DNUFIS1PHI1 ./ DFLX1;
DNapJ1DPHI1 = DNJ1 ./ DFLX1;
DREMPHI1DPHI1 = DREMPHI1 ./ DFLX1;

DABS2PHI2DPHI2 = DABS2PHI2./ DFLX2;
DNapJDPHI2 = DNJ2 ./ DFLX2;
DNUFIS2PHI2DPHI2 = DNUFIS2PHI2 ./ DFLX2;

DABS1PHI1DPHI1(isinf(DABS1PHI1DPHI1)) = 0;
DNUFIS1PHI1DPHI1(isinf(DNUFIS1PHI1DPHI1)) = 0;
DNapJ1DPHI1(isinf(DNapJ1DPHI1)) = 0;
DREMPHI1DPHI1(isinf(DREMPHI1DPHI1)) = 0;


DABS2PHI2DPHI2(isinf(DABS2PHI2DPHI2)) = 0;
DNapJDPHI2(isinf(DNapJDPHI2)) = 0;
DNUFIS2PHI2DPHI2(isinf(DNUFIS2PHI2DPHI2)) = 0;

DABS1PHI1DPHI1(isnan(DABS1PHI1DPHI1)) = 0;
DNUFIS1PHI1DPHI1(isnan(DNUFIS1PHI1DPHI1)) = 0;
DNapJ1DPHI1(isnan(DNapJ1DPHI1)) = 0;
DREMPHI1DPHI1(isnan(DREMPHI1DPHI1)) = 0;

DABS2PHI2DPHI2(isnan(DABS2PHI2DPHI2)) = 0;
DNapJDPHI2(isnan(DNapJDPHI2)) = 0;
DNUFIS2PHI2DPHI2(isnan(DNUFIS2PHI2DPHI2)) = 0;

h = DZ:DZ:DZ*sizez;

figure()
hold on 
vec1(:) = mean(input.MOD_EQ_1_scaled,[1,2]);
vec2(:) = mean(feedback.MOD_EQ_1_scaled,[1,2]);
plot(h,vec1)
plot(h,vec2)
xlabel("Height (cm)","FontSize",16)
ylabel("Fast neutron flux (cm^{-2}s^{-1})","FontSize",16)
legend("100% fast flux", "50% fast flux","Location","south")



%clearvars input feedback
delete temp/temp_input_vars.mat
delete temp/temp_feedback_vars.mat
%% Load Xerom_data
load("../For_Christophe/input_100_50/XS_data_100_50.mat")
load("../For_Christophe/input_100_50/additional_data_100_50.mat")
load("../For_Christophe/input_100_50/GEOM_data_100_50.mat")
load("../For_Christophe/input_100_50/RESULTS.mat")
% load("../For_Christophe/feedback_100_50/XS_data_100_50.mat")
% load("../For_Christophe/feedback_100_50/additional_data_100_50.mat")
% load("../For_Christophe/feedback_100_50/GEOM_data_100_50.mat")
% load("../For_Christophe/feedback_100_50/RESULTS.mat")

%% Create vectors and matrices
DV = DX*DY*DZ;
power = input_additional.reactor_power;       
sizex = I_MAX; %number of nodes in the x-direction
sizey = J_MAX; %number of nodes in the y-direction
sizez = K_MAX; %number of nodes in the z-direction
M=10;
v1 = 2E9; % fast neutron velocity in cm/s
v2 = 2.2E5; % thermal neutron velocity in cm/s 
KAPPA = 0.2976e-10; % Guessed value of Kappa (source unknown) same as is used in the 1G homogenous model
%FB = 3.08228E-19; % original Feedback coeficient Calculated in mathematica from the one group homogenous model drdp*(Sigma_a + D*B^2)Kappa*Int(Phi_eq,V)*SigmaF
%FB= 6.9799e-19; % test feedback
%FB = 1.3e-18;
sigmaX = 2.7e-18;
V_inv = [1/v1, 0 ; 0, 1/v2];
K_VALUE = lambda; % K values / eigenvalues of the modes
DV = DX*DY*DZ; % Discreet volume element
%KN = XS.KN; % kappa / nu
%NU = kappa./KN; % space dependent nu based on guess of a non space dependent kappa (APPROXIMATION)
KFIS1 = input.KF1;
KFIS2 = input.KF2;
SIGF1 = KFIS1./KAPPA; % Fast fission cross section
SIGF2 = KFIS2./KAPPA; % Thermal fission cross section
ZERO = zeros(size(NUFIS1)); % zero element matching the size of the reactor
ONE = ones(size(NUFIS1)); % unit element matching the size of the reactor
%CR_SA = [CR_SA1,ZERO;ZERO,CR_SA2];
F = 1/K_VALUE(1).*[input.NUFIS1, input.NUFIS2;ZERO,ZERO]; % Fission matrix
%D = [-D1,ZERO;ZERO,-D2]; % Diffusion coefficient matrix
MOD = [input_results.MOD1;input_results.MOD2]; % vector of solutions to the forward problem
MOD_adj = [input_results.MOD1_adj,input_results.MOD2_adj]; % vector of solutions to the adjoint problem
MOD_EQ = abs([input_results.MOD1(:,:,:,1);input_results.MOD2(:,:,:,1)]); % Vector of only the equilibrium neutron flux solution
% MOD_EQ_INT = DV * sum(MOD1(:,:,:,1)+MOD2(:,:,:,1),'all');
% SIGF_PHI_INT_2G=DV*sum(G2_inner_product([SIGF1,SIGF2],MOD_EQ,"vector","vector"),"all")
bsq =   0.0002 ; % Bsq taken from 2G-HOM model Assumes cylindrical shape of the reactor
KFISINT =  DV*sum(G2_inner_product([KFIS1,KFIS2],MOD_EQ,"vector","vector"),"all");
PS = power*K_VALUE(1)/KFISINT;
MOD_EQ_scaled= PS*MOD_EQ; % Vector of only the equilibrium neutron flux solution scaled by the power
MOD_EQ_1_scaled= MOD_EQ_scaled(1:sizey,:,:);
MOD_EQ_2_scaled= MOD_EQ_scaled(sizey+1:end,:,:);
MOD_eq_MAT =[MOD_EQ_1_scaled,ZERO; ZERO,MOD_EQ_2_scaled]; % Diagonal matrix containing the solutions to the forward problem
MOD_UPPER = [MOD_EQ_2_scaled, ZERO ; ZERO, ZERO]; % Costom matrix used in the equations 
MOD_LOWER = [ZERO, ZERO; MOD_EQ_2_scaled, ZERO]; % Costom matrix used in the equations
%Make matrices of delta values only for termal feedbacks
DeltaXS = [DNUFIS1PHI1DPHI1-DREMPHI1DPHI1-DNapJ1DPHI1-DABS1PHI1DPHI1,DNUFIS2PHI2DPHI2; DREMPHI1DPHI1, -DNapJDPHI2-DABS2PHI2DPHI2 ];
intial_AO = (DV*sum(KFIS1(:,:,sizez/2+1:sizez).*MOD_EQ_1_scaled(:,:,sizez/2+1:sizez)+KFIS2(:,:,sizez/2+1:sizez).*MOD_EQ_2_scaled(:,:,sizez/2+1:sizez),'all')-DV*sum(KFIS1(:,:,1:sizez/2).*MOD_EQ_1_scaled(:,:,1:sizez/2)+KFIS2(:,:,1:sizez/2).*MOD_EQ_2_scaled(:,:,1:sizez/2),'all'))/(DV*sum(KFIS1(:,:,:).*MOD_EQ_1_scaled(:,:,:)+KFIS2(:,:,:).*MOD_EQ_2_scaled(:,:,:),'all'));

mat1(:,:) = mean(input.STA_FLX2,3);
vec1(:) = mean(input.STA_FLX2,[1,2]);
mat2(:,:) = mean(MOD_EQ_2_scaled,3);
vec2(:) = mean(MOD_EQ_2_scaled,[1,2]);
figure(1)
diff_mat = (mat1-mat2)./mat1;
surf(diff_mat)
xlabel("X (nodes)")
ylabel("Y (Nodes)")
title("radial thermal relative error(averaged) ((SIM-CS)/SIM)")
view(2)
colorbar
figure(2)
plot(vec1)
hold on
plot(vec2)
xlabel("Height (nodes)")
ylabel("Thermal neutron flux (cm^{-2}s^{-1})")
legend("SIMULATE","CORE SIM")
figure(3)
plot((vec1-vec2)./vec1)
xlabel("Height (Nodes)")
ylabel("Relative thermal error ((SIM-CS)/SIM)")


mat1(:,:) = mean(input.STA_FLX1,3);
vec1(:) = mean(input.STA_FLX1,[1,2]);
mat2(:,:) = mean(MOD_EQ_1_scaled,3);
vec2(:) = mean(MOD_EQ_1_scaled,[1,2]);
figure(4)
reldiff_thermal = (STA_FLX2-MOD_EQ_2_scaled)./STA_FLX2;
diff_mat = mean(reldiff_thermal,3);
surf(diff_mat)
xlabel("X (nodes)")
ylabel("Y (Nodes)")
title("Radial fast relative error(averaged) ((SIM-CS)/SIM)")
view(2)
colorbar
figure(5)
plot(vec1)
hold on
plot(vec2)
xlabel("Height (nodes)")
ylabel("Fast neutron flux (cm^{-2}s^{-1})")
legend("SIMULATE","CORE SIM")
figure(6)
plot((vec1-vec2)./vec1)
xlabel("Height (Nodes)")
ylabel("Relative fast error ((SIM-CS)/SIM)")

PHI_bottom = zeros(1,M);
PHI_top = zeros(1,M);
for mode = 1:M
    PHI_bottom(mode) = DV*sum(G2_inner_product([SIGF1(:,:,1:sizez/2),SIGF2(:,:,1:sizez/2)],MOD(:,:,1:sizez/2,mode),"vector","vector"),1:3);
    PHI_top(mode) = DV*sum(G2_inner_product([SIGF1(:,:,sizez/2+1:end),SIGF2(:,:,sizez/2+1:end)], MOD(:,:,sizez/2+1:end,mode),"vector","vector"),1:3);
end

%SIG = [ABS1+REM, ZERO ; -REM, ABS2]; % Matrix containing the absorbtion and removal cross sections
GAMMAI = 1/K_VALUE(1).*[gammaI*SIGF1,gammaI*SIGF2 ; ZERO, ZERO ]; % matrix containing the production of Iodine from fission
GAMMAX = 1/K_VALUE(1).*[gammaX*SIGF1,gammaX*SIGF2 ; ZERO, ZERO ]; % matrix containing the production of xenon from fission
%VECX = [ONE;ZERO]; % transformation vector
VECXT = [ONE,ZERO]; % transformation vector
%MATX = [ZERO,ONE;ZERO,ZERO]; %transformation matrix


I0 = 1/lambdaI*G2_inner_product(VECXT, G2_inner_product(GAMMAI,MOD_EQ_scaled,"matrix","vector"),"vector","vector"); % I_0(r) = \hat(X)^T \cdot 1/(keff*lambda_I)Gamma_I\times Phi_0
X0 = (lambdaI*I0 + G2_inner_product(VECXT,G2_inner_product(GAMMAX,MOD_EQ_scaled,"matrix","vector"),"vector","vector"))./(lambdaX+sigmaX*MOD_EQ_scaled(sizex+1:end,:,:)); 
%X0 = 1/K_VALUE(1)*(gammaI+gammaX)*(SIGF1.*MOD_EQ_scaled(1:32,:,:)+SIGF2.*MOD_EQ_scaled(33:end,:,:))./(lambdaX+sigmaX*MOD_EQ_scaled(33:end,:,:));
Xe_UPPER = [ZERO,X0;ZERO,ZERO]; % Custom matrix used in the equations



%% test properties
close all
Volume = DV * sum(MOD1(:,:,:,1)~=0,'all');
active_volume = DV * sum(NUFIS1~=0,'all');
SIGF_PHI = SIGF1.*MOD_EQ_scaled(1:sizex,:,:)+SIGF2.*MOD_EQ_scaled(sizex+1:end,:,:);
SIGF_PHI_average = DV * sum(SIGF_PHI,'all')/active_volume;
SIGF1_average = DV * sum(SIGF1,'all')/active_volume;
SIGF2_average = DV * sum(SIGF2,'all')/active_volume;
MOD_EQ1_average = DV * sum(MOD_EQ_scaled(1:sizex,:,:),'all')/Volume;
MOD_EQ2_average = DV * sum(MOD_EQ_scaled(sizex+1:end,:,:),'all')/Volume;
X0_Average = DV * sum(X0,'all')/active_volume;
I0_Average = DV * sum(I0,'all')/active_volume;
MOD_EQ2_max = max(MOD_EQ_scaled(sizex+1:end,:,:),[],'all');
X0_max = max(X0,[],'all');
I0_max = max(I0,[],'all');
DSA_average = DV*sum(DABS2PHI2DPHI2,'all')/Volume;
DSF_average = DV*sum(DABS2PHI2DPHI2,'all')/active_volume;
DNaJ_average = DV*sum(DNapJDPHI2,'all')/Volume;
delta_NUFIS1 = feedback.NUFIS1-input.NUFIS1;
delta_ABS1 = feedback.ABS1 - input.ABS1;
delta_NUFIS2 = feedback.NUFIS2-input.NUFIS2;
delta_ABS2 = feedback.ABS2 - input.ABS2;
delta_REM = feedback.REM - input.REM;
Deltarho = DV * sum( (delta_NUFIS1-delta_ABS1-delta_REM).*MOD1_adj(:,:,:,1).*MOD1(:,:,:,1) + delta_NUFIS2.*MOD1_adj(:,:,:,1).*MOD2(:,:,:,1)+delta_REM.*MOD2_adj(:,:,:,1).*MOD1(:,:,:,1) - delta_ABS2.*MOD2_adj(:,:,:,1).*MOD2(:,:,:,1),'all') ...
    / (DV * sum(input.NUFIS1.*MOD1_adj(:,:,:,1).*MOD1(:,:,:,1)+input.NUFIS2.*MOD1_adj(:,:,:,1).*MOD2(:,:,:,1),'all'));
Deltarho_reduced = DV * sum(  delta_NUFIS2.*MOD1_adj(:,:,:,1).*MOD2(:,:,:,1) - delta_ABS2.*MOD2_adj(:,:,:,1).*MOD2(:,:,:,1),'all') ...
    / (DV * sum(input.NUFIS1.*MOD1_adj(:,:,:,1).*MOD1(:,:,:,1)+input.NUFIS2.*MOD1_adj(:,:,:,1).*MOD2(:,:,:,1),'all'));
keff2_SIM = 1.01602;
reac1_SIM = 0;
reac2_SIM = (keff2_SIM-1)/keff2_SIM;
x = DX*(1:sizex);
y = DY*(1:sizey);
z = DZ*(1:sizez);
[X,Y] = meshgrid(x,y);
mat_DSA(:,:) = DABS2PHI2DPHI2(:,:,13);
vec_DSA(:) = DABS2PHI2DPHI2(17,17,:);
figure()
surf(X,Y,mat_DSA);
view(2)
colorbar
zlabel('^{\Delta\Sigma_{a,2}\phi_{2}}/_{\Delta\phi_{1}}',"FontSize",20)
xlabel("X (cm)")
ylabel("Y (cm)")
title("Absorption")
figure()
plot(vec_DSA);
xlabel("Height (cm)")
ylabel('^{\Delta\Sigma_{a,2}\phi_{2}}/_{\Delta\phi_{2}}',"FontSize",20)
title("Absorption")
mat_DSF(:,:) = DNUFIS2PHI2DPHI2(:,:,13);
vec_DSF(:) = DNUFIS2PHI2DPHI2(17,17,:);
figure()
surf(mat_DSF)
view(2)
colorbar
zlabel('^{\Delta\Sigma_{f,2}\phi_{2}}/_{\Delta\phi_{2}}',"FontSize",20)
xlabel("X (cm)")
ylabel("Y (cm)")
title("Fission")
figure()
plot(vec_DSF);
xlabel("Height (cm)")
ylabel('^{\Delta\Sigma_{f,2}\phi_{2}}/_{\Delta\phi_{2}}',"FontSize",20)
title("Fission")
mat_DNaJ(:,:) = DNapJDPHI2(:,:,13);
vec_DNaJ(:) = DNapJDPHI2(17,17,:);
figure()
surf(mat_DNaJ)
view(2)
colorbar
zlabel('^{\Delta(\nabla J_{2}\phi_{2})}/_{\Delta\phi_{2}}',"FontSize",20)
xlabel("X (cm)")
ylabel("Y (cm)")
title("Leakage")
figure()
plot(vec_DNaJ);
xlabel("Height (cm)")
ylabel('^{\Delta(\nabla J_{2}\phi_{2})}/_{\Delta\phi_{2}}',"FontSize",20)
title("Leakage")
%% clear variables
%clear NUFIS1 NUFIS2 n_iter n_restart ABS1 ABS2 FLX1 FLX2 EIG_MET D1 D2 MOD1_adj MOD2_adj REM RES_FLX RES_MOD RES_MOD_adj ...
 %   XS KN DX DY DZ conv_ERAM conv_POW lambda lambda_adj gammaI gammaX v1 v2
%% initialise parameters
PHID_F_PHI = zeros(1,M);
PHID_V_PHI = zeros(1,M);
PHID_PHI = zeros(1,M);
PHID_GAMMAI_PHI = zeros(M);
PHID_GAMMAX_PHI = zeros(M);
PHID_FB_PHI = zeros(M);
PHID_PHIUPPER_PHI=zeros(M);
PHID_PHILOWER_PHI=zeros(M);
PHID_X0_PHI=zeros(M);
%PHID_CR_PHI = zeros(M);

temp_GAMMAX_PHI = zeros(sizex*2,sizey,sizez,M);
temp_GAMMAI_PHI = zeros(sizex*2,sizey,sizez,M);
temp_F_PHI = zeros(sizex*2,sizey,sizez,M);
temp_V_PHI = zeros(sizex*2,sizey,sizez,M);
%temp_PHI_eq_mat_PHI=zeros(sizex*2,sizey,sizez,M);
temp_FB_PHI= zeros(sizex*2,sizey,sizez,M);
temp_PHIUPPER_PHI=zeros(sizex*2,sizey,sizez,M);
temp_PHILOWER_PHI=zeros(sizex*2,sizey,sizez,M);
temp_X0_PHI = zeros(sizex*2,sizey,sizez,M);
%temp_CR_PHI = zeros(sizex*2,sizey,sizez,M);
%% Calculate kets
for n = 1:M
    temp_F_PHI(:,:,:,n)=G2_inner_product(F(:,:,:),MOD(:,:,:,n),"matrix","vector"); % | F* Phi_n>
    %temp_PHI_eq_mat_PHI(:,:,:,n) = G2_inner_product(MOD_eq_MAT,MOD(:,:,:,n),"matrix","vector"); %|Phi_0_mat * Phi_n>
    temp_FB_PHI(:,:,:,n) = G2_inner_product(DeltaXS,MOD(:,:,:,n),'matrix','vector'); % | FB *Phi_n>
    temp_GAMMAX_PHI(:,:,:,n) = G2_inner_product(GAMMAX,MOD(:,:,:,n),"matrix","vector"); %  |1/k_0 *Gamma_X * Phi_n>
    temp_GAMMAI_PHI(:,:,:,n) = G2_inner_product(GAMMAI,MOD(:,:,:,n),"matrix","vector"); % | 1/k_0 Gamma_I * Phi_n>
    temp_PHIUPPER_PHI(:,:,:,n) = G2_inner_product(MOD_UPPER,temp_F_PHI(:,:,:,n),"matrix","vector"); % | \bar{X} \Phi_0 \hat{X}^T * F* Phi_n >
    temp_PHILOWER_PHI(:,:,:,n) = G2_inner_product(MOD_LOWER,temp_F_PHI(:,:,:,n),"matrix","vector"); % | \tilde{X} \Phi_0 \hat{X}^T * F* Phi_n >
    temp_X0_PHI(:,:,:,n) = G2_inner_product(Xe_UPPER,MOD(:,:,:,n),"matrix","vector"); % |\bar{X} * X0 * Phi_n >
    %temp_CR_PHI(:,:,:,n) = G2_inner_product(CR_SA,MOD(:,:,:,n),"matrix","vector"); % CP
end

%% Calculate bra-kets 

for m = 1:M
    PHID_F_PHI(m) = DV * sum(G2_inner_product(MOD_adj(:,:,:,m),temp_F_PHI(:,:,:,m),"vector","vector"),'all'); %<Phi^dagger_m| F Phi_m>
    temp_V_PHI(:,:,:,m) = G2_inner_product(V_inv,MOD(:,:,:,m),"scalar_matrix","vector"); % |v^-1 Phi_m>
    PHID_V_PHI(m) = DV * sum(G2_inner_product(MOD_adj(:,:,:,m),temp_V_PHI(:,:,:,m),"vector","vector"),"all"); %<Phi^dagger_m| v^-1 Phi_m>
    PHID_PHI(m) = DV * sum(G2_inner_product(MOD_adj(:,:,:,m),MOD(:,:,:,m),"vector","vector"),"all"); % <Phi^dagger_m|Phi_m>
    for n = 1:M
        PHID_GAMMAX_PHI(m,n) = DV*sum(G2_inner_product(MOD_adj(:,:,:,m),temp_GAMMAX_PHI(:,:,:,n),"vector","vector"),"all"); % <Phi^dagger_m | Gamma_X * Phi_n >
        PHID_GAMMAI_PHI(m,n) =  DV*sum(G2_inner_product(MOD_adj(:,:,:,m),temp_GAMMAI_PHI(:,:,:,n),"vector","vector"),"all"); %<Phi^dagger_m | Gamma_I * Phi_n > 
        %PHID_PHI_eq_mat_PHI(m,n) =  DV*sum(G2_inner_product(MOD_adj(:,:,:,m),temp_PHI_eq_mat_PHI(:,:,:,n),"vector","vector"),"all"); %<Phi^dagger_m | Phi_eq_mat Phi_n >
        PHID_FB_PHI(m,n) = DV * sum(G2_inner_product(MOD_adj(:,:,:,m),temp_FB_PHI(:,:,:,n),"vector","vector"),"all");
        PHID_PHIUPPER_PHI(m,n) = DV*sum(G2_inner_product(MOD_adj(:,:,:,m),temp_PHIUPPER_PHI(:,:,:,n),"vector","vector"),"all"); %<Phi^dagger_m | \bar{X} \Phi_0 \hat{X}^T * F Phi_n >
        PHID_PHILOWER_PHI(m,n) = DV*sum(G2_inner_product(MOD_adj(:,:,:,m),temp_PHILOWER_PHI(:,:,:,n),"vector","vector"),"all"); %<Phi^dagger_m | \tilde{X} \Phi_0 \hat{X}^T * F Phi_n >
        PHID_X0_PHI(m,n) = DV*sum(G2_inner_product(MOD_adj(:,:,:,m),temp_X0_PHI(:,:,:,n),"vector","vector"),"all"); %<Phi^dagger_m | \bar{X} X_0 Phi_n >
       %PHID_CR_PHI(m,n) = DV*sum(G2_inner_product(MOD_adj(:,:,:,m),temp_CR_PHI(:,:,:,n),"vector","vector"),"all"); % <Phi^dagger_m | CR Phi_n >
    end
end


LAMBDA = PHID_V_PHI./ PHID_F_PHI; % <Phi^dagger_m |v^-1 Phi_n>/ <Phi^dagger_m |F Phi_m> 
%clearvars temp*
%% test magnitudes of terms
sprintf("Printing eig separation term")
1./LAMBDA'.*(1./K_VALUE(1)-1./K_VALUE)
sprintf("printing feedback term")
1./LAMBDA.*PHID_FB_PHI./PHID_F_PHI
sprintf("printing diagonal of feedback term")
diag(1./LAMBDA.*PHID_FB_PHI./PHID_F_PHI)
C11 = diag(1./LAMBDA.*PHID_FB_PHI./PHID_F_PHI) + 1./LAMBDA'.*(1./K_VALUE(1)-1./K_VALUE)
sprintf("printing flux xenon term")
C13 = diag(-1./LAMBDA.*sigmaX.*PHID_PHILOWER_PHI./(PHID_F_PHI.^2).*PHID_PHI)   
sprintf("Printing Iodine creation term")
C21 = diag(PHID_GAMMAI_PHI./PHID_PHI)
sprintf("Printing first xenon absorption term")
-sigmaX.*PHID_X0_PHI./PHID_PHI
C31 = diag((PHID_GAMMAX_PHI-sigmaX*PHID_X0_PHI)./PHID_PHI)
sprintf("Printing second xenon absorption term")
C33 = diag(-sigmaX.*PHID_PHIUPPER_PHI./PHID_F_PHI)
FB_spatial = G2_inner_product(MOD_adj(:,:,:,1),temp_FB_PHI(:,:,:,1),"vector","vector");
analyse_3d_var(FB_spatial,"Feedback")
save data/tempFile
%%
[time_2G,state_values_2G] = ode_Nsolve();
%%
%delete data/tempFile.mat
%toc
%%
close all
% figure(1)
% tot_sol_phi = state_values_2G(:,1:3:M*3);
% phi_point_1(1,:) = MOD1(25,25,17,:);
% temp_spatial_point_1 = tot_sol_phi*phi_point_1';
% plot(time_2G(1:end)/3600, temp_spatial_point_1(1:end))
% ylim([-8E6,8E6])
% xlabel("Time [h]",'Fontsize', 14);
% ylabel("Normalized neutron flux [AU]",'Fontsize', 14);

figure(1)
%yyaxis left
plot(time_2G(1:end)/3600,state_values_2G(1:end,([2 1 5]-1)*3+1),"LineWidth",2)
grid on
%yticks(-3e9:0.5e9:3e9)
%xlim([0 150])    
%ylim([-3.1E9,2.1E9])
ylabel("Amplitude [cm^{-2}s^{-1}]",'Fontsize', 22)
xlabel("Time [h]",'Fontsize', 22)
legend("First axial harmonic", "Fundamental mode", "Second axial harmonic","Fontsize",22)
ax2 = gca;
ax2.FontSize = 22;