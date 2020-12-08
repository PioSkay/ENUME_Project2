classdef TASK3 < TASK2
    properties
        %the order of a polynomial
        n = 4;
    end
    methods
        function obj = TASK3()
            obj@TASK2;
        end
        function [it, x, obj] = Lagurre(obj, x)
            it = 1;
            while abs(obj.f_x(x)) > obj.accuracy
               %the steps here are the same as in the case of MM2
               c = obj.f_x(x);
               b = obj.f_x_d(x);
               a = obj.f_x_dd(x);
               %the only difference here is that we are taking the 
               r = sqrt((obj.n-1)*((obj.n-1)*b^2-obj.n*a*c));
               d1 = b + r;
               d2 = b - r;
               if (abs(d1) > abs(d2))
                   d = d1;
               else
                   d = d2;
               end
               %here we do need to obtain the values for a next iteration
               z = obj.n*c/d;
               x = x - z;
               it = it + 1;
            end
            return
        end
    end
end