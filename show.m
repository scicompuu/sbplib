function show(str)
    assertType(str, {'string', 'char'})
    val = evalin('caller',str);
    fprintf('%s => %s\n\n', str, toString(val));
end