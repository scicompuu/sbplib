function saveFigurePosition()
    defaultPosition = get(0,'defaultfigureposition');
    f = gcf;
    defaultPosition(3:4) = f.Position(3:4);
    set(0,'defaultfigureposition',defaultPosition);
end