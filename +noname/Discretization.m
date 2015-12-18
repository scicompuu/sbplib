classdef Discretization < handle
    properties (Abstract)
        name         %Short description
        description  %Longer description
        order        %Order of accuracy
        % h            % scalar desciribing the grid spacing.. (IS THIS THE RIGHT PLACE FOR THIS?)
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
        ts = getTimestepper(obj, method, k, cfl)

        % Calculates a timestep for the discretization and a given timestepping method.
        % Can take order, differnt types of scaling in h, or other parameters in Discr into
        % account.  opt is a struct that among other things may contain
        %   method -- time stepping method for which to give a timestep.
        %   cfl    -- [optioanal] a cfl constant to use to calculate the timetep.
        %             if skipped getTimestep should use a precomputed value.
        %   k      -- timestep to use
        k = getTimestep(obj, opt)

        % getTimeSnapshot returns a struct which represents the solution in ts at current time.
        % if ts is empty or 0 a representation of the initial conditions be returned.
        repr = getTimeSnapshot(obj, ts)

        % Sets up a plot of the discretisation
        %     update is a function_handle accepting a timestepper that updates the plot to the
        %            state of the timestepper
        %     type allows for different kinds of plots. Some special values are used by the lib. 'animate' and 'plot' for example
        [update,hand] = setupPlot(obj, type)

    end
end