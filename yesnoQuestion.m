function b = yesnoQuestion(question, defaultAnswer)
    default_arg('defaultAnswer','nodefault');

    yesAnswer = {'y','Y','yes','Yes','YES'};
    noAnswer = {'n','N','no','No','NO'};

    switch defaultAnswer
        case 'nodefault'
            optionString = '[y/n]';
        case yesAnswer
            optionString = '[Y/n]';
            yesAnswer{end+1} = '';
        case noAnswer
            optionString = '[y/N]';
            noAnswer{end+1} = '';
        otherwise
            error('Unrecognized default answer: %s', defaultAnswer);
    end

    b = [];
    while isempty(b)
        answer = input([question ' ' optionString ': '],'s');
        switch answer
            case yesAnswer
                b = true;
            case noAnswer
                b = false;
        end
    end
end