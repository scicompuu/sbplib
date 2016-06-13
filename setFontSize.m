% setFontSize(fig, pts)
% Sets all fontsizes within a figure to size pts
% The default value for pts is 16.
function setFontSize(fig, pts)
    default_arg('pts', 16);
    set(findall(fig,'-property','FontSize'),'FontSize',pts);
end