classdef filter
    methods(Static)
        function g=mean(f)
            [n,m,c]=size(f);
            h=double(ones(3,3))*(1/9);
            g=uint8(convolution(double(f),double(h)));
        end

        function g=gauss(f,sigma)
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

        function result=noise_filter(img, method, Q, alpha)
            filter_size = 3;
            pad = floor(filter_size / 2);

            padded = padarray(img, [pad, pad], 1, 'both');
            [x, y, z] = size(img);
            result = zeros(size(img));

            for i = 1:x
                for j = 1:y
                    for k = 1:z
                        neighborhood_pixels = padded(i:i+filter_size-1, j:j+filter_size-1, k);
                        switch method
                            case 'min'
                                value = min(neighborhood_pixels(:));
                            case 'max'
                                value = max(neighborhood_pixels(:));
                            case 'median'
                                value = median(neighborhood_pixels(:));
                            case 'mean'
                                value = mean(neighborhood_pixels(:));
                            case 'geometric'
                                value = exp(mean(log(neighborhood_pixels(:)) + 1e-6));
                            case 'harmonic'
                                value = filter_size^2 / sum(1 ./ neighborhood_pixels(:));
                            case 'contraharmonic'
                                value = sum(neighborhood_pixels(:).^(Q+1)) / sum(neighborhood_pixels(:).^Q);
                            case 'midpoint'
                                value = (max(neighborhood_pixels(:)) + min(neighborhood_pixels(:))) / 2;
                            case 'alpha-trimmed'
                                sorted = sort(neighborhood_pixels(:));
                                value = mean(sorted(alpha+1:end-alpha));
                        end
                        result(i, j, k) = value;
                    end
                end
            end
        end
    end
end