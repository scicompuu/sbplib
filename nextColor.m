function c = nextColor(ah)
    default_arg('ah', gca);

    c = ah.ColorOrder(ah.ColorOrderIndex, :);
end