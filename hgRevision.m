% Returns the short mercurial revision Id.
%  ok is false if there are uncommited changes.
function [revId, ok] = hgRevision()
    [~, s] = system('hg id -i');
    revId = strtrim(s);

    ok = s(end) ~= '+';
end
