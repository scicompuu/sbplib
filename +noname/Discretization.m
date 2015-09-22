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
        %     k is a desired timestep
        %     cfl is a choses cfl constant used to set the timestep. ignored if k is set.
        ts = getTimestepper(obj,method,k,cfl)

        % Calculates a timestep for the discretization and a given timestepping method.
        % Can take order, differnt types of scaling in h, or other parameters in Discr into
        % account.
        %   method -- time stepping method for which to give a timestep.
        %   cfl    -- [optioanal] a cfl constant to use to calculate the timetep.
        %             if skipped getTimestep should use a precomputed value.
        k = getTimestep(obj, method, cfl)


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