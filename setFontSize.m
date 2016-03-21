function setFontSize(fig, pts)
    default_arg('pts', 16);
    set(findall(fig,'-property','FontSize'),'FontSize',pts);
end