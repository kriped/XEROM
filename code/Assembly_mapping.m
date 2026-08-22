clearvars; 
% number of fuel assemblies = 157
% Naming: Void = 0; Fuel assembly = 1:157; Reflector_cylinder = -1; Reflector
% top = -2; reflector bottom = -3

load('../input/R3C40/EOC/Refinement0/XS_data.mat')
load('../input/R3C40/EOC/Refinement0/FUM_data.mat')

% Create a mask marking the position of the fuel assemblies.
Assembly_mask = zeros(17,17,26); % half rods at the periphery in the top half of the core
%Assembly_mask(1,7:11,1:end) = -1; % Reflector at periphery
Assembly_mask(2,8:10,14:end-1) = 1; % Assemblies
%Assembly_mask(2,[7,11],14:end-1) = -1; % Reflector at periphery
%Assembly_mask(2,7:11,1:13) = -1; % Reflector below half rods
%Assembly_mask(2,7:11,1) = -1; % Reflector at bottom
%Assembly_mask(2,7:11,end) = -1; % Reflector at top
Assembly_mask(3,6:12,2:end-1) = 1; % Assemblies
%Assembly_mask(3,[5,13],2:end-1) = -1; % Reflector at periphery
%Assembly_mask(3,5:13,1) = -1; % Reflector at bottom
%Assembly_mask(3,5:13,end) = -1; % Reflector at top
Assembly_mask(4,5:13,2:end-1) = 1; % Assemblies
%Assembly_mask(4,[4,14],2:end-1) = -1; % Reflector at periphery
%Assembly_mask(4,4:14,1) = -1; % Reflector at bottom
%Assembly_mask(4,4:14,end) = -1; % Reflector at top
Assembly_mask(5,4:14,2:end-1) = 1; % Assemblies
%Assembly_mask(5,[3,15],2:end-1) = -1; % Reflector at periphery
%Assembly_mask(5,3:15,1) = -1; % Reflector at bottom
%Assembly_mask(5,3:15,end) = -1; % Reflector at top
Assembly_mask(6:7,3:15,2:end-1) = 1; % Assemblies
%Assembly_mask(6:7,[2,16],2:end-1) = -1; % Reflector at periphery
%Assembly_mask(6:7,2:16,1) = -1; % Reflector at bottom
%Assembly_mask(6:7,2:16,end) = -1; % Reflector at top
Assembly_mask(8:10,2:16,14:end-1) = 1; % half Assemblies at the periphery in the top half of the core
Assembly_mask(8:10,3:15,2:13) = 1; % half Assemblies at the periphery in the bottom half of the core
%Assembly_mask(8:10,[1,17],14:end-1) = -1; % Reflector at periphery with half rods
%Assembly_mask(8:10,[1:2,16:17],1:13) = -1; % Assemblies missing at the periphery in the bottom half of the core
%Assembly_mask(8:10,1:17,1) = -1; % Reflector bottom
%Assembly_mask(8:10,1:17,end) = -1; % Reflector top
Assembly_mask(11:12,3:15,2:end-1) = 1; % Assemblies
%Assembly_mask(11:12,[2,16],2:end-1) = -1; % Reflector at periphery
%Assembly_mask(11:12,2:16,1) = -1; %Reflector bottom
%Assembly_mask(11:12,2:16,end) = -1; %Reflector top
Assembly_mask(13,4:14,2:end-1) = 1; % Assemblies
%Assembly_mask(13,[3,15],2:end-1) = -1; % Reflector at periphery
%Assembly_mask(13,3:15,1) = -1; % Reflector at bottom
%Assembly_mask(13,3:15,end) = -1; % Reflector at top
Assembly_mask(14,5:13,2:end-1) = 1; % Assemblies
%Assembly_mask(14,[4,14],2:end-1) = -1; % Reflector at periphery
%Assembly_mask(14,4:14,1) = -1; % Reflector at bottom
%Assembly_mask(14,4:14,end) = -1; % Reflector at top
Assembly_mask(15,6:12,2:end-1) = 1; % Assemblies
%Assembly_mask(15,[5,13],2:end-1) = -1; % Reflector at periphery
%Assembly_mask(15,5:13,1) = -1; % Reflector at bottom
%Assembly_mask(15,5:13,end) = -1; % Reflector at top
Assembly_mask(16,8:10,14:end-1) = 1; % half rods at the periphery in the top half of the core
%Assembly_mask(16,[7,11],14:end-1) = -1; % Reflector at periphery
%Assembly_mask(16,7:11,1) = -1; % Reflector at bottom
%Assembly_mask(16,7:11,end) = -1; % Reflector at top
%Assembly_mask(16,7:11,1:13) = -1; % Reflector at below the half assemblies
%Assembly_mask(17,7:11,:) = -1; % Reflector at periphery

Reflector_mask = ABS1 ~=0 & NUFIS1 == 0;

% label the each fuel assembly
counter = 0;
Assembly_map = Assembly_mask;
for i = 1:size(Assembly_map, 1)
    for j = 1:size(Assembly_map, 2)
        if Assembly_mask(i, j, 20) == 1 
            counter = counter + 1;
            Assembly_map(i, j, 2:end-1) = counter;
        end
    end
end
clear i j
%Handle half rods
Assembly_map(8:10,2,1:13) = 0;
Assembly_map(8:10,end-1,1:13) = 0;
Assembly_map(2,8:10,1:13) = 0;
Assembly_map(end-1,8:10,1:13) = 0;
Assembly_map = repelem(Assembly_map,2,2,1); % Split all assemblies in 4 parts radially
Assembly_map = Assembly_map + Reflector_mask*(-1);
Assembly_map_refined = repelem(Assembly_map,2,2,2); % Refine by factor of two in all directions


save("Assembly_map.mat","Assembly_map","Assembly_map_refined")