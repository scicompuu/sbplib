classdef Discretization < handle
    properties (Abstract)
        name         %Short description
        description  %Longer description
        order        %Order of accuracy
    end

    methods (Abstract)
        % Prints some info about the discretisation
        printInfo()

        % Return the number of DOF
        n = size(obj)

        % Returns a timestepper for integrating the discretisation in time
        %     method is a string that states which timestepping method should be used.
        %          The implementation should switch on the string and deliver
        %          the appropriate timestepper. It should also provide a default value.
        %     time_align is a time that the timesteps should align with so that for some
        %                integer number of timesteps we end up exactly on time_align
        ts = getTimestepper(obj,method,time_align) %% ???

        % Sets up movie recording to a given file.
        %     saveFrame is a function_handle with no inputs that records the current state
        %               as a frame in the moive.
        saveFrame = setupMov(obj, file)

        % Sets up a plot of the discretisation
        %     update is a function_handle accepting a timestepper that updates the plot to the
        %            state of the timestepper
        [update,hand] = setupPlot(obj)

    end

    methods(Abstract,Static)
        % Compare two functions u and v in the discrete l2 norm.
        e = compareSolutions(u, v)

        % Compare the functions u to the analytical function g in the discrete l2 norm.
        e = compareSolutionsAnalytical(u, g)
    end
end