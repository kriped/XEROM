function [dSAdPhi, dNFdPhi, dKFdPhi] = Calculate_dSAdphi(type,AX,AY,AZ,BX,BY,BZ,CX,CY,CZ)
    %Calculate the dSA_j/dphi_k = dSAdphi{jk} 
    load("dSAdphi.mat", "Case_1*","Case_3*") 
    
    
    
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
    

    if type == "full"
        
        dSAdPhi = zeros(2*34,34,26);
        dNFdPhi = zeros(2*34,34,26);
        dKFdPhi = zeros(2*34,34,26);
        %core_region = (isnan(Case_1_FLX1)-1)*(-1);
        %dSAdphi11 = (Case_3_SA1-Case_1_SA1)./(Case_3_FLX1-Case_1_FLX1);
        %dSAdphi21 = (Case_3_SA2-Case_1_SA2)./(Case_3_FLX1-Case_1_FLX1);
        dSAdPhi(1:34,:,:) = (Case_3_SA1-Case_1_SA1)./(Case_3_FLX2-Case_1_FLX2);
        dSAdPhi(35:end,:,:) = (Case_3_SA2-Case_1_SA2)./(Case_3_FLX2-Case_1_FLX2);
        dNFdPhi(1:34,:,:) = (Case_3_NF1-Case_1_NF1)./(Case_3_FLX2-Case_1_FLX2);
        dNFdPhi(35:end,:,:) = (Case_3_NF2-Case_1_NF2)./(Case_3_FLX2-Case_1_FLX2);
        dKFdPhi(1:34,:,:) = (Case_3_KF1-Case_1_KF1)./(Case_3_FLX2-Case_1_FLX2);  
        dKFdPhi(35:end,:,:) = (Case_3_KF2-Case_1_KF2)./(Case_3_FLX2-Case_1_FLX2);
    end
    
    % %flxdif1 = (Case_3_FLX1-Case_1_FLX1);
    % flxdif2 = (Case_3_FLX2-Case_1_FLX2);
    % 
    % %radmean1(:) = mean(flxdif1,[1,2],"omitmissing");
    % radmean2(:) = mean(flxdif2,[1,2],"omitmissing");
    % 
    % % figure(1)
    % % plot(radmean1)
    % % title("Taken from intrusive1.out")
    % % xlabel("Height from core bottom (nodes)")
    % % ylabel("Change in fast neutron flux (cm^{-2}s^{-1})")
    % 
    % figure(2)
    if type == "2G"
        % plot(radmean2)
        % title("Taken from intrusive1.out")
        % xlabel("Height from core bottom (nodes)")
        % ylabel("Change in thermal neutron flux (cm^{-2}s^{-1})")
        %Homogenised delta SA
        [delta_Homo_SA, delta_Homo_FLX] = homogenise_delta_XS(Case_1_SA1,Case_1_SA2,Case_1_FLX1,Case_1_FLX2,Case_3_SA1,Case_3_SA2,Case_3_FLX1,Case_3_FLX2,"2G");
        dSAdPhi = delta_Homo_SA./delta_Homo_FLX;

        %Homogenised delta NF
        [delta_Homo_NF,~] = homogenise_delta_XS(Case_1_NF1,Case_1_NF2,Case_1_FLX1,Case_1_FLX2,Case_3_NF1,Case_3_NF2,Case_3_FLX1,Case_3_FLX2,"2G");
        dNFdPhi = delta_Homo_NF./delta_Homo_FLX;

        %Homogenised delta KF
        [delta_Homo_KF,~] = homogenise_delta_XS(Case_1_KF1,Case_1_KF2,Case_1_FLX1,Case_1_FLX2,Case_3_KF1,Case_3_KF2,Case_3_FLX1,Case_3_FLX2,"2G");
        dKFdPhi = delta_Homo_KF./delta_Homo_FLX;        
    end
    % test reaction rate
    
    % hom_test_A = Homogenised_XSA*Homogenised_FLXA*V
    % mean_test_A = sum(Case_1_SA1.*Case_1_FLX1*dV+Case_1_SA2.*Case_1_FLX2*dV,"all","omitmissing")
    % hom_test_B = Homogenised_XSB*Homogenised_FLXB*V
    % mean_test_B = sum(Case_3_SA1.*Case_3_FLX1*dV+Case_3_SA2.*Case_3_FLX2*dV,"all","omitmissing")
    
    %Calculate homogeneous 2-group 
    if type == "1G"
        % calculate homogenous 1 - group feedback coefficient
        %Homogenised delta SA
        [delta_Homo_SA, delta_Homo_FLX] = homogenise_delta_XS(Case_1_SA1,Case_1_SA2,Case_1_FLX1,Case_1_FLX2,Case_3_SA1,Case_3_SA2,Case_3_FLX1,Case_3_FLX2,"1G");
        dSAdPhi = delta_Homo_SA./delta_Homo_FLX;

        %Homogenised delta NF
        [delta_Homo_NF,~] = homogenise_delta_XS(Case_1_NF1,Case_1_NF2,Case_1_FLX1,Case_1_FLX2,Case_3_NF1,Case_3_NF2,Case_3_FLX1,Case_3_FLX2,"1G");
        dNFdPhi = delta_Homo_NF./delta_Homo_FLX;

        %Homogenised delta KF
        [delta_Homo_KF,~] = homogenise_delta_XS(Case_1_KF1,Case_1_KF2,Case_1_FLX1,Case_1_FLX2,Case_3_KF1,Case_3_KF2,Case_3_FLX1,Case_3_FLX2,"1G");
        dKFdPhi = delta_Homo_KF./delta_Homo_FLX;
    end
end