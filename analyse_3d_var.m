function analyse_3d_var(var,label)
    sizex = size(var,1);
    sizey = size(var,2);
    sizez = size(var,3);
    Assembly_pitch = 21.5036; %cm
    Node_height = 15.2400; % cm
    DX = Assembly_pitch/4;
    DY = Assembly_pitch/4;
    DZ = Node_height/2;
    
    Width = DX*sizex;
    Height = DZ*sizez;
    x = linspace(0,Width,sizex);
    y = linspace(0,Width,sizey);
    z = linspace(0,Height,sizez);
    [X,Y] = meshgrid(x,y);

    mat_middle(:,:) = var(:,:,ceil(sizez/2));
    mat_mean(:,:) = mean(var,3,"omitmissing");
    vec(:) = var(ceil(sizex/2),ceil(sizey/2),:);
    var(var == 0) = NaN;
    vec_avg(:) = mean(var,[1,2],"omitmissing");

    figure()
    surf(X,Y,mat_mean)
    zlabel(label+" averaged axially","FontSize",16)
    xlabel("X (cm)","FontSize",16)
    ylabel("Y (cm)","FontSize",16)
    title(label+" averaged axially",FontSize=16)
    colorbar
    view(2)
    figure()
    surf(X,Y,mat_middle)
    zlabel(label+" middle plane","FontSize",16)
    xlabel("X (cm)","FontSize",16)
    ylabel("Y (cm)","FontSize",16)
    title(label+" middle plane","FontSize",16)
    colorbar
    view(2)
    figure()
    plot(z,vec)
    xlabel("Height (cm)","FontSize",16)
    ylabel(label + " central line","FontSize",16)
    figure()
    plot(z,vec_avg)
    xlabel("Height (cm)","FontSize",16)
    ylabel(label+" averaged radially","FontSize",16)


end