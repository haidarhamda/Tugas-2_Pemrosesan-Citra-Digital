classdef filter
    methods(Static)
        function g=mean(f)
            [n,m,c]=size(f);
            h=ones(3,3)*(1/9);
            g=convolution(f,double(h));
        end

        function g=gauss(f,sigma)
            % g=imgaussfilt(f,sigma);
            % filter_size=7;
            % disp(double(fspecial("gaussian",[filter_size,filter_size],sigma)));
            g=convolution(f,double(fspecial("gaussian",[sigma,sigma],sigma)));
        end

        function g=ilpf(f)
            [n,m,c]=size(f);
        end
    end
end