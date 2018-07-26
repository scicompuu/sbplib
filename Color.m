classdef Color
    properties(Constant)
        blue      = [0.000 0.447 0.741];
        red       = [0.850 0.325 0.098];
        yellow    = [0.929 0.694 0.125];
        purple    = [0.494 0.184 0.556];
        green     = [0.466 0.674 0.188];
        lightblue = [0.301 0.745 0.933];
        darkred   = [0.635 0.078 0.184];
        black     = [0.000 0.000 0.000];
        white     = [1.000 1.000 1.000];
        colors = { Color.blue, Color.red, Color.yellow, Color.green, Color.purple, Color.lightblue, Color.darkred, Color.black, Color.white};
        markers = {'+', 'o', '*', '.', 'x', 'square', 'diamond', 'v', '^', '>', '<', 'pentagram', 'hexagram'};
        lineStyles = {'-', '--', ':', '-.'};

        solidMarkers = {'o', 'square', 'diamond', 'v', 'pentagram', '^', '>', '<', 'hexagram'};

        notabilityYellow     = [100.0   99.0    22.0    ]/100;
        notabilityOrange     = [97.0    61.0    15.0    ]/100;
        notabilityRed        = [92.0    22.0    18.0    ]/100;

        notabilityLightGreen = [67.0    99.0    22.0    ]/100;
        notabilityGreen      = [39.0    72.0    30.0    ]/100;
        notabilityDarkGreen  = [5.0     45.0    24.0    ]/100;

        notabilityLightBlue  = [8.0     45.0    98.0    ]/100;
        notabilityBlue       = [11.0    30.0    69.0    ]/100;
        notabilityDarkBlue   = [7.0     5.0     46.0    ]/100;

        notabilityPink       = [99.0    13.0    82.0    ]/100;
        notabilityPurple     = [64.0    11.0    77.0    ]/100;
        notabilityDarkPurple = [38.0    19.0    56.0    ]/100;

        notabilityBrown      = [58.0    40.0    22.0    ]/100;
        notabilityGrey       = [59.0    59.0    59.0    ]/100;
        notabilityBlack      = [0.0     0.0     0.0     ]/100;
    end

    methods(Static)
        function sample()
            % Filled and non-filled markers?


            function showMarkers(x0, y0, lx, ly, color, filled)
                n = length(Color.markers);
                s = ceil(sqrt(n));

                x = linspace(x0, x0 + lx, s);
                y = linspace(y0, y0 + ly, s);

                [X,Y] = meshgrid(x,y);

                for i = 1:n
                    lh = line(X(i),Y(i));
                    lh.Marker = Color.markers{i};
                    lh.MarkerSize = 12;
                    lh.Color = color;

                    if filled
                        lh.MarkerFaceColor = color;
                    end
                end
            end

            function showColors(x0, y0, lx, ly)
                n = length(Color.colors);
                s = ceil(sqrt(n));

                x = linspace(x0, x0 + lx, s);
                y = linspace(y0, y0 + ly, s);

                [X,Y] = meshgrid(x,y);

                for i = 1:n
                    lh = line(X(i),Y(i));
                    lh.Marker = 'o';
                    lh.MarkerFaceColor = Color.colors{i};
                    lh.Color = Color.colors{i};
                    lh.MarkerSize = 12;
                end
            end

            function showLines(y0, ly, A, w)
                n = length(Color.lineStyles);
                x = linspace(0,1,100);
                y = linspace(y0, y0+ ly, n);
                for i = 1:n
                    lh = line(x, y(i) + A*sin(pi*x*w));
                    lh.LineWidth = 2;
                    lh.LineStyle = Color.lineStyles{i};
                end
            end


            fh = figure();

            y0 = 0.1;
            dl = 0.15;
            margin = 0.07;
            showMarkers(     0.1, y0, dl, dl, Color.blue, false);
            showMarkers(  0.9-dl, y0, dl, dl, Color.blue,  true);
            showColors( 0.5-dl/2, y0, dl, dl);

            showLines(0.6, 0.2, 0.1, 5);

            xlim([0 1])
            ylim([0 1])
            axis square

        end
    end
end