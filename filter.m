classdef filter
    methods(Static)
        function g=mean(f)
            [n,m,c]=size(f);
            h=double(ones(3,3))*(1/9);
            g=uint8(convolution(double(f),double(h)));
        end

        function g=gauss(f,sigma)
            % g=imgaussfilt(f,sigma);
            % filter_size=7;
            % disp(double(fspecial("gaussian",[filter_size,filter_size],sigma)));
            if(sigma<=0)
                sigma=1;
            end
            g=uint8(convolution(double(f),double(fspecial("gaussian",[sigma,sigma],sigma))));
        end

        function g=ilpf(f,d0)
            [n,m,c]=size(f);
            P=2*n;
            Q=2*m;
            g=zeros(n,m,c);
            for channel=1:c
                fp=padarray(im2double(f(:,:,channel)),[n,m],'post');
                % imshow(fp);
                F=fftshift(fft2(fp));
                D0 = d0;
                u = 0:(P-1); 
                v = 0:(Q-1);
                idx = find(u > P/2); 
                u(idx) = u(idx) - P; 
                idy = find(v > Q/2); 
                v(idy) = v(idy) - Q;
                [V, U] = meshgrid(v, u); 
                D = sqrt(U.^2 + V.^2); 
                H = double(D <=D0); 
                H = fftshift(H); 
                G = H.*F; 
                G1 = ifftshift(G);
                G2 = real(ifft2(G1));
                g(:,:,channel) = G2(1:n, 1:m);
            end
        end

        function g=glpf(f,d0)
            [n,m,c]=size(f);
            P=2*n;
            Q=2*m;
            g=zeros(n,m,c);
            for channel=1:c
                fp=padarray(im2double(f(:,:,channel)),[n,m],'post');
                F=fftshift(fft2(fp));
                D0 = d0;
                u = 0:(P-1); 
                v = 0:(Q-1);
                idx = find(u > P/2); 
                u(idx) = u(idx) - P; 
                idy = find(v > Q/2); 
                v(idy) = v(idy) - Q;
                [V, U] = meshgrid(v, u); 
                D = sqrt(U.^2 + V.^2); 
                H = exp(-(D.^2)./(2*(D0^2))); 
                H = fftshift(H); 
                % figure;imshow(H);
                % figure, mesh(H);
                LPF_f = H.*F;
                LPF_f2=real(ifft2(LPF_f));
                g(:,:,channel)=LPF_f2(1:n, 1:m);
            end
        end

        function g=blpf(f,d0,N2)
            [n,m,c]=size(f);
            P=2*n;
            Q=2*m;
            g=zeros(n,m,c);
            for channel=1:c
                fp=padarray(im2double(f(:,:,channel)),[n,m],'post');
                F=fftshift(fft2(fp));
                D0 = d0;
                u = 0:(P-1); 
                v = 0:(Q-1);
                idx = find(u > P/2); 
                u(idx) = u(idx) - P; 
                idy = find(v > Q/2); 
                v(idy) = v(idy) - Q;
                [V, U] = meshgrid(v, u); 
                D = sqrt(U.^2 + V.^2); 
                N = N2;
                H = 1./(1 + (D./D0).^(2*N)); 
                H = fftshift(H); 
                % figure;imshow(H);
                % figure, mesh(H);
                G = H.*F; 
                G1 = ifftshift(G);
                G2 = real(ifft2(G1));
                disp(size(G2));
                % imshow(G2);
                g(:,:,channel) = G2(1:n, 1:m);
            end
        end

        function g=ihpf(f,d0)
            [n,m,c]=size(f);
            P=2*n;
            Q=2*m;
            g=zeros(n,m,c);
            for channel=1:c
                fp=padarray(im2double(f(:,:,channel)),[n,m],'post');
                % imshow(fp);
                F=fftshift(fft2(fp));
                D0 = d0;
                u = 0:(P-1); 
                v = 0:(Q-1);
                idx = find(u > P/2); 
                u(idx) = u(idx) - P; 
                idy = find(v > Q/2); 
                v(idy) = v(idy) - Q;
                [V, U] = meshgrid(v, u); 
                D = sqrt(U.^2 + V.^2); 
                H = double(D >D0); 
                H = fftshift(H); 
                % figure;imshow(H);
                % figure, mesh(H);
                G = H.*F; 
                G1 = ifftshift(G);
                G2 = real(ifft2(G1));
                g(:,:,channel) = G2(1:n, 1:m);
            end
        end

        function g=ghpf(f,d0)
            [n,m,c]=size(f);
            P=2*n;
            Q=2*m;
            g=zeros(n,m,c);
            for channel=1:c
                fp=padarray(im2double(f(:,:,channel)),[n,m],'post');
                F=fftshift(fft2(fp));
                D0 = d0;
                u = 0:(P-1); 
                v = 0:(Q-1);
                idx = find(u > P/2); 
                u(idx) = u(idx) - P; 
                idy = find(v > Q/2); 
                v(idy) = v(idy) - Q;
                [V, U] = meshgrid(v, u); 
                D = sqrt(U.^2 + V.^2); 
                H = exp(-(D.^2)./(2*(D0^2))); 
                H = 1-H;
                H = fftshift(H); 
                % figure;imshow(H);
                % figure, mesh(H);
                LPF_f = H.*F;
                LPF_f2=real(ifft2(LPF_f));
                g(:,:,channel)=LPF_f2(1:n, 1:m);
            end
        end

        function g=bhpf(f,d0,N2)
            [n,m,c]=size(f);
            P=2*n;
            Q=2*m;
            g=zeros(n,m,c);
            for channel=1:c
                fp=padarray(im2double(f(:,:,channel)),[n,m],'post');
                F=fftshift(fft2(fp));
                D0 = d0;
                u = 0:(P-1); 
                v = 0:(Q-1);
                idx = find(u > P/2); 
                u(idx) = u(idx) - P; 
                idy = find(v > Q/2); 
                v(idy) = v(idy) - Q;
                [V, U] = meshgrid(v, u); 
                D = sqrt(U.^2 + V.^2); 
                N = N2;
                H = 1./(1 + (D./D0).^(2*N)); 
                H = 1-H;
                H = fftshift(H); 
                % figure;imshow(H);
                % figure, mesh(H);
                G = H.*F; 
                G1 = ifftshift(G);
                G2 = real(ifft2(G1));
                disp(size(G2));
                % imshow(G2);
                g(:,:,channel) = G2(1:n, 1:m);
            end
        end
    end
end