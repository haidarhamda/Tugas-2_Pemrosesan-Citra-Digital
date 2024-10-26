function [noise, result]=saltAndPepper(filename, noise_density, method)
    img = imread(filename);
    img = im2double(img);

    [x, y, z] = size(img);

    pixels = x * y;
    salt = round(noise_density * pixels / 2);
    pepper = round(noise_density * pixels / 2);

    salt_coords = randperm(pixels, salt);
    pepper_coords = randperm(pixels, pepper);

    noise = img;
    for i = 1:z
        noise(salt_coords + (i - 1)*pixels) = 1;
        noise(pepper_coords + (i - 1)*pixels) = 0;
    end

    result = filter.noise_filter(noise, method);
end