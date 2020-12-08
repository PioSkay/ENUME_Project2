function compareRootsTask2(~)
    x0 = -.5;
    x1 = .5;
    x2 = 0;
    obj = TASK2;
    obj2 = TASK1(0.01);
    while abs(obj.f_x(x0) - obj.f_x(x1)) > obj.accuracy...
            && abs(obj.f_x(x1) - obj.f_x(x2)) > obj.accuracy...
            && abs(obj.f_x(x0) - obj.f_x(x2)) > obj.accuracy
        [it, x3, obj] = MullerMM2(obj, x2);
        if abs(imag(x3)) < obj.accuracy
            fprintf('Method MM2:\n');
            fprintf('Iterations: %d, %f \n', it, x3);
            fprintf('Newtons method:\n');
            obj2.f_x = obj.f_x;
            obj2.f_x_d = obj.f_x_d;
            [x3_2, it2, obj2] = NewtonsMethod(obj2, x2);
            fprintf('Iterations: %d, %f \n', it2, x3_2);
        end
        syms x;
        locF = @(x)(x - x3);
        obj.f_x = @(x)(obj.f_x(x)./locF(x));
        %calculate new derivatives
        obj.f_x_d = matlabFunction(diff(obj.f_x, x));
        obj.f_x_dd = matlabFunction(diff(obj.f_x_d, x));
    end
end