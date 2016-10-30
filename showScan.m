function showScan(x,y,z,depth)
%SHOWSCAN Summary of this function goes here
%…®√Ëœ‘ æµ—ø®∂˚æÿ’Û
figure(2);
hold on;
for i=1:size(z,3)
    himage=pcolor(x(:,:,i),y(:,:,i),depth(:,:,i));
    shading interp
    colormap(jet) 
    colorbar;
    pause(0.05);
    pause;
    delete(himage);
end
close 2;

end

