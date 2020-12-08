classdef TASK2
   properties
        %plot of a function
        f_x;
        %derivative of a plot
        f_x_d;
        %second derivative of a plot
        f_x_dd;
        %accuracy of the calculations
        accuracy = 1e-6;
   end
   methods
       function obj = TASK2()
           syms x
           obj.f_x = @(x)(-1*x.^4-7*x.^3+7*x.^2+3*x+9);
           obj.f_x_d = matlabFunction(diff(obj.f_x, x));
           obj.f_x_dd = matlabFunction(diff(obj.f_x_d, x));
       end
       function plotTask2(obj)
           range = -8:0.01:4;
           plot(range, obj.f_x(range), "-r");
           grid on;
           line(xlim, [0,0], 'Color', 'k', 'LineWidth', 0.1);
           legend('f(x) = -1*x^4-7*x^3+7*x^2+3*x+9');
           xlabel('X');
           ylabel('Y');
       end
       function [it, x3, obj] = MullerMM1(obj, x0, x1, x2)
           h1 = x1 - x0;
           h2 = x2 - x1;
           [DELTA2, d] = update(obj, x0, x1, x2, h1, h2);
           it = 3;
           %if there is no root it is possible to go into an infinite loop
           %we need to prevent ourselves from that
           while it <= 100
               %calculating output (out1/out2) in the report
                b = DELTA2 + h2*d;
                D = (b^2 - 4*obj.f_x(x2)*d)^(1/2);
                %determining which one is out.min
                d1 = b + D;
                d2 = b - D;
                if abs(d2) < abs(d1)
                    E = d1;
                else
                    E = d2;
                end
                %caluculating actual minout
                minout = -2*obj.f_x(x2)/E;
                %establishing the x3
                x3 = x2 + minout;
                %checking where the new root approximiation sitisfy the
                %accuracy
                if abs(minout) < obj.accuracy
                    return;
                end
                %if not reassign the vaules
                %and redo all steps
                x0 = x1;
                x1 = x2;
                x2 = x3;
                h1 = x1 - x0;
                h2 = x2 - x1;
                [DELTA2, d] = update(obj, x0, x1, x2, h1, h2);
                it = it + 1;
           end
       end
       function [it, x2, obj] = MullerMM2(obj, x2)
           it = 1;
           while it <= 100
               %step 1 in the report
               c = obj.f_x(x2);
               b = obj.f_x_d(x2);
               aa = obj.f_x_dd(x2);
               %step 2 in the report
               out1 = -2*c/(b + (b^2 -2*aa*c)^(1/2));
               out2 = -2*c/(b - (b^2 -2*aa*c)^(1/2));
               if abs(out1) > abs(out2)
                   out = out2;
               else
                   out = out1;
               end
               if abs(out) < obj.accuracy
                   return;
               end
               %step 3 in the report
               x2 = x2 + out;
               it = it + 1;
           end
       end
   end
   methods (Access = private) 
       function [DEL2, d] = update(obj, x0, x1, x2, h1, h2)
           DEL1 = (obj.f_x(x1) - obj.f_x(x0))/h1;
           DEL2 = (obj.f_x(x2) - obj.f_x(x1))/h2;
           d = (DEL2 - DEL1)/(h2 + h1);
       end
   end
end