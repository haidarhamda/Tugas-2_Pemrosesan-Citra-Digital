function [h]=convolution(f,g)
    [n, m, c] = size(f);
    [x,y]=size(g);
    h=double(zeros(n,m,c));
    padx=floor(x/2);
    pady=floor(y/2);
    paddedf=padarray(f,[padx,pady]);

    for channel=1:c
        imgChannel = paddedf(:,:,channel);
        tmpChannel=double(zeros(n,m));
        for i=1+padx:n-padx
            for j=1+pady:m-pady
                tmp=double(imgChannel(i:i+x-padx,j:j+y-pady));
                tmpChannel(i+padx,j+pady)=sum(sum(tmp.*double(g)));
            end
        end
        h(:,:,channel)=tmpChannel;
    end
end

% img = imread("img\camera.bmp");
% g=[-1 -1 -1; -1 8 -1; -1 -1 -1];
% h=(coonvolution(im2double(img),g));
% % disp(h)
% imshow(h)
% hm=uint8(conv2(im2double(img), double(g)));
% figure,imshow(hm)