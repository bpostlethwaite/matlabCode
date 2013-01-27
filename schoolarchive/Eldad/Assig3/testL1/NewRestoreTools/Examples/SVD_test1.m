clc
clear
disp(' ')
disp('     Example 2: Grain data, using different boundary conditions ')
disp('--------------------------------------------------------------------')
disp(' ')

load SVD_test1
figure(1)
  imshow(b,[]), title('Blurred image'), colormap(jet)
figure(2)
  imshow(log(PSF+1),[]), title('PSF'), colormap(jet)
figure(3)
  imshow(b, [0, 0])

disp('>> whos')
whos
disp('>> ...')

disp('>> Az = psfMatrix(PSF, ''zero'');')
Az = psfMatrix(PSF, 'zero');
disp('>> Ar = psfMatrix(PSF, ''reflexive'');')
Ar = psfMatrix(PSF, 'reflexive');
disp('>> ...')

disp('>> [Uz, Sz, Vz] = svd(Az);')
[Uz, Sz, Vz] = svd(Az);
disp('>> [Ur, Sr, Vr] = svd(Ar);')
[Ur, Sr, Vr] = svd(Ar);
disp('>> ...')

disp('>> xz = TSVD(Uz, Sz, Vz, b, ''help'');')
xz = TSVD(Uz, Sz, Vz, b, 'help');
disp('>> xr = TSVD(Ur, Sr, Vr, b, ''help'');')
xr = TSVD(Ur, Sr, Vr, b, 'help');

figure(2)
  imshow(xz, [-0.001, 0.01]), colormap(jet), title('TSVD for Zero BC'), drawnow
figure(3)
  imshow(xr, [-0.001, 0.01]), colormap(jet), title('TSVD for Reflexive BC'), drawnow

