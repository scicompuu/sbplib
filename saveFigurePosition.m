function saveFigurePosition()
    defaultPosition = get(0,'defaultfigureposition');
    f = gcf;
    defaultPosition(1:2) = f.Position(1:2);
    set(0,'defaultfigureposition',defaultPosition);
end