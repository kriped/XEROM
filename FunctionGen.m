function functionlist = FunctionGen(mend)

strlist = strings(1,mend);
for m = 1:mend
    %For each m make a new set of 3 equations
    PhiString = sprintf(['1/LAMBDA(%1$i)*(1/K_VALUE(1)-1/K_VALUE(%1$i))*s((%1$i-1)*3+1)'],m);
    IodineString=sprintf(['-lambdaI*s((%1$i-1)*3+2)'],m);
    XenonString = sprintf(['lambdaI*s((%1$i-1)*3+2)' ...
        '-lambdaX*s((%1$i-1)*3+3)'],m);
    for n = 1:mend
        % Add new terms to each phi and xenon eqs for each equation
        newstr_P = sprintf(['-1/LAMBDA(%1$i)*PHID_PHI_eq_mat_PHI(%1$i,%2$i)/PHID_F_PHI(%1$i)*s((%2$i-1)*3+1)' ...
            '-1/LAMBDA(%1$i)*sigmaX*PHID_PHILOWER_PHI(%1$i,%2$i)/(PHID_F_PHI(%1$i)^2)*PHID_PHI(%1$i)*s((%2$i-1)*3+3)'],m,n);
        PhiString = [PhiString newstr_P];
        newstr_I = sprintf(['+PHID_GAMMAI_PHI(%1$i,%2$i)*s((%2$i-1)*3+1)/PHID_PHI(%1$i)'],m,n);
        IodineString = [IodineString newstr_I];
        newstr_X = sprintf(['+(PHID_GAMMAX_PHI(%1$i,%2$i)-sigmaX*PHID_X0_PHI(%1$i,%2$i))/PHID_PHI(%1$i)*s((%1$i-1)*3+1)' ...
            '-sigmaX*PHID_PHIUPPER_PHI(%1$i,%2$i)/PHID_F_PHI(%1$i)*s((%1$i-1)*3+3)'],m,n);
        XenonString = [XenonString newstr_X];
    end
    %edit end of strings
    %PhiString(end-1:end) = [];
    PhiString = [PhiString,';',newline];
    IodineString = [IodineString, ';',newline];
    XenonString = [XenonString,';',newline];
    %XenonString(end-1:end) = [];
    %XenonString(end-2:end) = [');' newline];
    strlist(m)= [PhiString newline IodineString newline XenonString newline]; % Combine for complete function string
end
% make string into vector and edit the end
functionlist = char(join(strlist));
functionlist=functionlist(1:end-1);
%functionString = ['g = @(t,s)[' functionlist ']'];
%eval(functionString);
%func = str2func(functionString); %convert string vector to function
end