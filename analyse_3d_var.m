function analyse_3d_var(var,label)
    sizex = size(var,1);
    sizey = size(var,2);
    sizez = size(var,3);
    
    mat(:,:) = var(:,:,ceil(sizez/2));
    vec(:) = var(ceil(sizex/2),ceil(sizey/2),:);
    var(var == 0) = NaN;
    vec_avg(:) = mean(var,[1,2],"omitmissing");

    figure()
    surf(mat)
    zlabel(label,"FontSize",18)
    xlabel("X (cm)")
    ylabel("Y (cm)")
    figure()
    plot(vec)
    xlabel("Height (cm)")
    ylabel(label,"FontSize",18)
    figure()
    plot(vec_avg)
    xlabel("Height (cm)")
    ylabel(label+" Averaged radially","FontSize",18)


end