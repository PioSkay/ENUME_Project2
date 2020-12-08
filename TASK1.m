classdef TASK1
    properties
        %jump between each individual point.
        jump = 0.001;
        %range to plot
        range;
        %plot of a function
        f_x;
        %derivative of a plot
        f_x_d;
        %accuracy of the calculations
        accuracy = 1e-6;
    end
    methods
        function obj = TASK1(jump)
            obj.jump = jump;
            range = -5:obj.jump:10;
            obj.range = range;
            obj.f_x =@(range)(3.1 - 3*range - exp(-range));
            obj.f_x_d = @(range)(exp(-range) - 3);
        end
        function plotTask1(obj)
            if obj.jump < 0.01
                plot(obj.range, obj.f_x(obj.range), "-r");
            else
                plot(obj.range, obj.f_x(obj.range), "-ro");
            end
            grid on;
            line(xlim, [0,0], 'Color', 'k', 'LineWidth', 0.1);
            legend('f(x) = 3.1 - 3*x - e^-^x');
            xlabel('X');
            ylabel('Y');
        end
        %this method finds the cross of a x axis
        %this is rather simple 
        % if f(x) > 0 and f(x-1) < 0 or f(x) < 0 and f(x-1) > 0 does the
        % job
        %there a little change mainly related to the size of range that we
        %find f(x) == 0 elseif part handle that.
        function obj = findZerosA(obj)
            rootID = 1;
            for i = 2:size(obj.range, 2)
                locX = obj.f_x(obj.range(i-1));
                locX1 = obj.f_x(obj.range(i));
                if locX < 0 && locX1 > 0 || locX > 0 && locX1 < 0
                    fprintf('For the root ID: %d, BisectionMethod initial interval: [%f, %f], ',rootID, obj.range(i-1), obj.range(i));
                    fprintf('and Newtons initial guess: %f, we get: \n',(obj.range(i-1)+ obj.range(i))/2);
                    [locRes, iter, obj] = BisectionMethod(obj, obj.range(i-1), obj.range(i));
                    [locRes2, iter2, obj] = NewtonsMethod(obj, (obj.range(i-1)+ obj.range(i))/2);
                    fprintf('Solution, Bisection: %.8f, number of interations: %i\n', locRes, iter);
                    fprintf('Solution, Newtons: %.8f, number of interations: %i\n', locRes2, iter2);
                    rootID = rootID + 1;
                elseif locX == 0
                    disp(obj.range(i - 1));
                elseif locX1 == 0
                    disp(obj.range(i));
                end
            end
        end
        function [result, iter, obj] = BisectionMethod(obj, a, b)
            iter = 0;
            %first of all as we know from the report
            %we need to check if we are able to find a root
            if obj.f_x(a) * obj.f_x(b) > 0
                disp("Could not find the root bacause f(x) with a given range a -> b does not have any.");
                return;
            end
            % now lets precced all of the steps described in the report.
            %step 1 and 2
            result = (a + b)/2;
            err = abs(obj.f_x(result));
            while err > obj.accuracy
                %just to keep a track of a number of iterations
                iter = iter +1;
                %step 4 and 5
                if obj.f_x(a) * obj.f_x(result) < 0
                    b = result;
                else
                    a = result;
                end
                %step 2
                result = (a + b)/2;
                err = abs(obj.f_x(result));
            end
        end
        % @param[in] a : Initial guess.
        function [locRes, iter, obj] = NewtonsMethod(obj, locRes)
            err = abs(obj.f_x(locRes));
            iter = 0;
            while err > obj.accuracy
                %formula from the report
                b = locRes - (obj.f_x(locRes)/obj.f_x_d(locRes));
                err = abs(b - locRes);
                locRes = b;
                iter = iter + 1;
            end
        end
    end
end