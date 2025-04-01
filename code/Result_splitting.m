%Code to split the RESULT.mat file into several files relevant to the
%MScNPP code
function Result_splitting(input_dir)
    load(input_dir + "RESULTS.mat","FLX1","FLX2","keff")
    save(input_dir + "FUM_data.mat","FLX1","FLX2","keff")
end
