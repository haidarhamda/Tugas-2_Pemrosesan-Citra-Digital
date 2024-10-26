function [h]=convolution(f,g)
    [n, m, c] = size(f);
    [x,y]=size(g);
    
    padx=floor(x/2);
    pady=floor(y/2);
    h=double(zeros(n,m,c));
    paddedf=padarray(f,[padx,pady]);
    % disp(size(paddedf));
    for channel=1:c
        imgChannel = paddedf(:,:,channel);
        tmpChannel=double(zeros(n,m));
        tmpChannel=padarray(tmpChannel,[padx,pady]);
        % disp(size(tmp));
        for i=1+padx:n+padx
            for j=1+pady:m+pady
                tmp=double(imgChannel(i-padx:i+padx,j-pady:j+pady));
                
                % tmpChannel(i+padx,j+pady)=sum(sum(tmp.*double(g)));
                % disp(size(tmp));
                tmpChannel(i,j)=sum(sum(tmp.*double(g)));
                
            end
        end
        % disp(h(:,:,channel));
        % disp(tmpChannel(4,:));
        % disp(size(tmpChannel))
        % disp(n+padx)

        h(:,:,channel)=tmpChannel(padx+1:n+padx,pady+1:m+pady);
        % disp(h(:,:,channel));
        % disp(h(4,:,channel));
        % disp(size(h(:,:,channel)));
    end
end

% img = imread("img\camera.bmp");
% g=[-1 -1 -1; -1 8 -1; -1 -1 -1];
% h=(coonvolution(im2double(img),g));
% % disp(h)
% imshow(h)
% hm=uint8(conv2(im2double(img), double(g)));
% figure,imshow(hm)