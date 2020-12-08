function compareRootsTask3(~)
    x0 = -.5;
    x1 = .5;
    x2 = 0;
    obj = TASK3;
    while abs(obj.f_x(x0) - obj.f_x(x1)) > obj.accuracy...
            && abs(obj.f_x(x1) - obj.f_x(x2)) > obj.accuracy...
            && abs(obj.f_x(x0) - obj.f_x(x2)) > obj.accuracy
        fprintf('Method Laguerre:\n');
        [it, x3, obj] = Lagurre(obj, x2);
        if abs(imag(x3)) < obj.accuracy
            fprintf('Iterations: %d, %f \n', it, x3);
        else
            fprintf('Iterations: %d, %f%+fj \n', it, real(x3), imag(x3));
        end
        syms x;
        locF = @(x)(x - x3);
        obj.f_x = @(x)(obj.f_x(x)./locF(x));
        %calculate new derivatives
        obj.f_x_d = matlabFunction(diff(obj.f_x, x));
        obj.f_x_dd = matlabFunction(diff(obj.f_x_d, x));
    end
    x0 = -.5;
    x1 = .5;
    x2 = 0;
    obj2 = TASK2;
    while abs(obj2.f_x(x0) - obj2.f_x(x1)) > obj2.accuracy...
            && abs(obj2.f_x(x1) - obj2.f_x(x2)) > obj2.accuracy...
            && abs(obj2.f_x(x0) - obj2.f_x(x2)) > obj2.accuracy
        fprintf('Method MM2:\n');
        [it, x3, obj2] = MullerMM2(obj2, x2);
        if abs(imag(x3)) < obj2.accuracy
            fprintf('Iterations: %d, %f \n', it, x3);
        else
            fprintf('Iterations: %d, %f%+fj \n', it, real(x3), imag(x3));
        end
        syms x;
        locF = @(x)(x - x3);
        obj2.f_x = @(x)(obj2.f_x(x)./locF(x));
        %calculate new derivatives
        obj2.f_x_d = matlabFunction(diff(obj2.f_x, x));
        obj2.f_x_dd = matlabFunction(diff(obj2.f_x_d, x));
    end
end