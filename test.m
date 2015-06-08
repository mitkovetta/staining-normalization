I1 = imread('example1.tif');
I2 = imread('example2.tif');

[Inorm1 H1 E1] = normalizeStaining(I1);
[Inorm2 H2 E2] = normalizeStaining(I2);

figure('Name', 'Original image 1'), imshow(I1, []);
figure('Name', 'Normalized image 1'), imshow(Inorm1, []);

figure('Name', 'Original image 2'), imshow(I2, []);
figure('Name', 'Normalized image 2'), imshow(Inorm2, []);
