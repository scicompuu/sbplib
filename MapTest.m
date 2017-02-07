function tests = MatTest()
    tests = functiontests(localfunctions);
end

function kvp = getKeyValuePairs()
    kvp = {
        {1,3},1;
        struct(), [1; 3; 4];
        [1,2; 4 3], struct();
        'Hej', struct('lol', 6);
        0, 'Nej';
    };
end

function testSetAndGet(testCase)
    keyValuePairs = getKeyValuePairs();

    map = Map();

    % Insert
    for i = 1:length(keyValuePairs)
        map(keyValuePairs{i,1}) = keyValuePairs{i,2};
    end

    % Validate output
    for i = 1:length(keyValuePairs)
        v = map(keyValuePairs{i,1});
        testCase.verifyEqual(v, keyValuePairs{i,2});
    end
end

function map = exampleMap()
    keyValuePairs = getKeyValuePairs();

    map = Map();

    % Insert
    for i = 1:length(keyValuePairs)
        map(keyValuePairs{i,1}) = keyValuePairs{i,2};
    end
end

function testLength(testCase)
    map = Map();
    testCase.verifyEqual(map.length, 0);

    map = exampleMap();
    testCase.verifyEqual(map.length, 5)
end


function testIsKey(testCase)
    map = exampleMap();

    keyValuePairs = getKeyValuePairs();
    keys = keyValuePairs(:,1);

    for i = 1:length(keys)
        testCase.verifyTrue(map.isKey(keys{i}));
    end

    notKeys = {
        'hej',
        [],
        1,
        {2,5},
    };

    for i = 1:length(notKeys)
        testCase.verifyFalse(map.isKey(notKeys{i}));
    end
end


function testRemove(testCase)
    map = exampleMap();

    remove(map, struct());

    testCase.verifyFalse(map.isKey(struct()));
end

% function testValues(testCase)
%     keyValuePairs = getKeyValuePairs();

%     map = exampleMap();

%     testCase.verifyEqual(values(map), keyValuePairs(:,2)');
% end

% function testKeys(testCase)
%     keyValuePairs = getKeyValuePairs();

%     map = exampleMap();

%     testCase.verifyEqual(keys(map), keyValuePairs(:,1)');
% end
