function inner_product = G2_inner_product(var1, var2,type1,type2)
    %Calculates the inner product of two spatially dependent matrices
    size1 = size(var1);
    size2 = size(var2);
    
    size_1x = size1(1);
    size_1y = size1(2);
    if type1 ~= "scalar_matrix"
        size_1z = size1(3);
    end
    
    size_2x = size2(1);
    size_2y = size2(2);
    if type2 ~= "scalar_matrix"
       size_2z = size2(3);
    end
    if (length(size1) > 3 || length(size2) > 3)
        error("One of the matrices have more than 3 dimensions")
    end
    
    switch type1
        case "vector"
            if size_1y == 2* size_1x
                var11 = var1(:,1:end/2,:);
                var12 = var1(:,end/2+1:end,:);
            elseif size_1x == 2* size_1y
                permute(var1,[size_1y,size_1x,size_1z])
                 var11 = var1(:,1:end/2,:);
                 var12 = var1(:,end/2+1:end,:);
                 warning("Vector was permuted to match correct dimmensions. Consider entering vector with correct dimensions")
            else
                error("Sizes are incompatible with this function")
            end
        case "matrix" % Matrix notation is [1,2;3,4]
            var11 = var1(1:end/2,1:end/2,:);
            var12 = var1(1:end/2,end/2+1:end,:);
            var13 = var1(end/2+1:end,1:end/2,:);
            var14 = var1(end/2+1:end,end/2+1:end,:);
        case "scalar_matrix" 
            var11 = var1(1,1);
            var12 = var1(1,2);
            var13 = var1(2,1);
            var14 = var1(2,2);
            type1 = "matrix";
    end
    
    switch type2
        case "vector"
            if size_2x == 2* size_2y
                var21 = var2(1:end/2,:,:);
                var22 = var2(end/2+1:end,:,:);
            elseif size_2y == 2* size_2x
                reshape(var2,[size_2y,size_2x,size_2z]);
                var21 = var2(1:end/2,:,:);
                var22 = var2(end/2+1:end,:,:);
                warning("Vector was permuted to match correct dimmensions. Consider entering vector with correct dimensions")
            else
                error("Sizes are incompatible with this function")
            end
        case "matrix" % Matrix notation is [1,2;3,4]
            var21 = var2(1:end/2,1:end/2,:);
            var22 = var2(1:end/2,end/2+1:end,:);
            var23 = var2(end/2+1:end,1:end/2,:);
            var24 = var2(end/2+1:end,end/2+1:end,:);
        case "scalar_matrix" 
            var21 = var2(1,1);
            var22 = var2(1,2);
            var23 = var2(2,1);
            var24 = var2(2,2);
            type2 = "matrix";
    end
    
    switch true
        case (type1=="matrix" && type2 == "matrix")
            inner_product = [var11.*var21+var12.*var23, var11.*var22+var12.*var24; var13.*var21+var14.*var23, var13.*var22+var14.*var24];
        case (type1 == "matrix" && type2 == "vector")
            inner_product = [var11.*var21+var12.*var22; var13.*var21+var14.*var22];
        case (type1=="vector" && type2 == "matrix")
            inner_product = [var11.*var21+var12.*var23, var11.*var22+var12.*var24];
        case (type1 == "vector" && type2 == "vector") 
            inner_product = var11.*var21+var12.*var22;
    end
    
      
end