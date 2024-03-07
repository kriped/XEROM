%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                             CORE_SIM_v1.2
%                Copyright (C) 2011  Christophe Demazière
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
%             FILE DEFINING THE PARAMETERS USED IN CORE_SIM.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% For explanatory notes about the use of this tool, see the file
% USERS_GUIDE.PDF in the directory "manuals".
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% Variable allowing lauching a MatLab Graphical User Interface (GUI) after
% completion of the calculations.
%
% VIZ=1 for launching the GUI (default), otherwise VIZ=0.
VIZ=1;

%%
% Variable allowing choosing whether the Explicitely Restarted Arnoldi
% Method (ERAM) or the power iteration method (POW) is to be used. Please
% note that the power iteration method uses Wielandt's method to calculate
% the different eigenmodes and a first guess of the eigenvalues is
% required. Such a guess of the eigenvalues is provided by an Arnoldi run
% without restart.
%
% In case of convergence problem for ERAM, it is recommended to switch to
% the power iteration method.
%
% EIG_MET=1 for ERAM (default) or EIG_MET=2 for POW.
EIG_MET=1;

%%
% Variable allowing getting the results even if some of the eigenmodes have
% not converged.
%
% BYP=0 (default) if you want to interrupt the program when the eigenmodes
% have not converged after 'n_restart' restarts (for ERAM) or
% after 'n_iter' iterations (for POW). BYP=1 permits the execution of the
% program even if the eigenmodes have not converged.
BYP=0;

%%
% Number of eigenmodes sought for (both for IRAM and POW).
%
neigs=5;% (10 as default)

%%
% Parameters used for the Explicitely Restarted Arnoldi Method (ERAM).
%
% If you are not familiar with ERAM, do not change these settings.
% Nevertheless, if you experience convergence problems during the
% calculation of the eigenmodes, changing the following parameters could
% help resolve such problems.
%
m=300;% dimension of the Krylov subspace (150 as default)
n_restart=20;% number of maximum restarts (20 as default)
conv_ERAM=100*eps;% convergence criteria on the residuals (100*eps as default)

%%
% Parameters used for the power iteration method (POW)
%
n_iter=400;% maximum number of iterations (200 as default)
conv_POW=100*eps;% convergence criteria on the residual (100*eps as default)