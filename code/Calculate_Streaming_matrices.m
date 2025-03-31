function Calculate_Streaming_matrices(input_dir)
load(input_dir+"XS_data.mat",'D1','D2','REM','ABS1','ABS2','NUFIS1','NUFIS2')
load(input_dir+"GEOM_data.mat",'DX','DY','DZ')
VAR={'D1','D2','REM','ABS1','ABS2','NUFIS1','NUFIS2','DX','DY','DZ'};
TEST_VAR=ismember((VAR),who);
RUN_TEST = 1;
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

save(input_dir + "STREAMING_MATRICES_data.mat","AX","AY","AZ","BX","BY","BZ","CX","CY","CZ","DX")
save(input_dir + "HELPER_data.mat","SHIFT","SHIFT_XY","SHIFT_XYZ","CONV","TYP","I_MAX","J_MAX","K_MAX")
end