function centerColorbar(ah)
    old = ah.CLim;

    l = max(abs(old));
    ah.CLim = [-l l];
end