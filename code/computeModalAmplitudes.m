function [a] =  computeModalAmplitudes(A,phi)
    DV = 220.2203; %HARD CODED ONLY FOR TESTING
    V = 5.2952e+07;
    nt = size(A,4);
    M = size(phi,4);
    A_mat = reshape(A, [], nt);        % [N x nt]
    phi_mat = reshape(phi, [], M);     % [N x M]
    a = DV/V*A_mat' * phi_mat;              % [nt x M]
end