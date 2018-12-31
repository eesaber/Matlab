function Plotsetting_gray()
    set(gca, 'Color',[1 1 1])
    g = gray(256);
    colormap(g(50:240,:));
end