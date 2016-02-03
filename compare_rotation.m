folder = 1; 
En = 103;


filedirectory = sprintf('/aline%i/rotdrum%i/o%i/',folder,folder,En);
file_avalanches =sprintf('%sAvalanches_%i.mat',filedirectory,En);
load(file_avalanches);
rotation_center_file = sprintf('%srotation_center_%i.mat',filedirectory,En);
load(rotation_center_file);

angle_aperture = 2.35;
rot_angle = 0.0032;
i_final = find(diff(Rotation_step(1,:))~=0);
i_initial = find(diff(Rotation_step(1,:))~=0)+1;d_R = Rotation_step(1,i_final)-Rotation_step(1,i_initial);
th = d_R*rot_angle;

for ii = 1:length(i_initial)
    before_fn = sprintf('/aline%i/rotdrum%i/o%02d/positions%02d_%05d.mat',folder,folder,En,En,Fn_imafile(i_final(ii)));
    after_fn = sprintf('/aline%i/rotdrum%i/o%02d/positions%02d_%05d.mat',folder,folder,En,En,In_imafile(i_initial(ii)));
    load(before_fn,'pxs','pys','Npf');
    pxb(1:Npf(351),ii) = pxs(1:Npf(351),351);
    pyb(1:Npf(351),ii) = pys(1:Npf(351),351);
    load(after_fn,'pxs','pys','Npf')
    pxa(1:Npf(1),ii) = pxs(1:Npf(1),1);
    pya(1:Npf(1),ii) = pys(1:Npf(1),1);
end

%% Check to see if it matches.
load(rotation_center_file);
na = 70;
d_steps = 30;
alpha = sum(th(na:na+d_steps))/2;
tran_pxb = (pxb(:,na)-xo)*cos(alpha)-(pyb(:,na)-yo)*sin(alpha)+xo;
tran_pyb = (pxb(:,na)-xo)*sin(alpha)+(pyb(:,na)-yo)*cos(alpha)+yo;
tran_pxa = (pxa(:,na+d_steps)-xo)*cos(-alpha)-(pya(:,na+d_steps)-yo)*sin(-alpha)+xo;
tran_pya = (pxa(:,na+d_steps)-xo)*sin(-alpha)+(pya(:,na+d_steps)-yo)*cos(-alpha)+yo;

plot(tran_pxa,tran_pya,'.',tran_pxb,tran_pyb,'.',x_in,y_in,'k');axis('equal');axis('ij');

%% Find overlap points

[~, trivialbondt1,trivialbondt2] = adjacent(tran_pxa,tran_pya,tran_pxb,tran_pyb,4);
x1 = tran_pxa((trivialbondt1))';
x2 = tran_pxb((trivialbondt2))';
x1x2 = [x1 x2];
y1 = tran_pya((trivialbondt1))';
y2 = tran_pyb((trivialbondt2))';
y1y2 = [y1 y2];


[low_x il_x] = min(x1x2);
vx(1) = low_x + D*4;
vy(1) = y1y2(il_x) +D*4;

[low_y il_y] = min(y1y2);
vy(2) = low_y + D*10;
vx(2) = x1x2(il_y);

[max_x im_x] = max(x1x2);
vx(3) = max_x - D*4;
vy(3) = y1y2(im_x) +D*4;

[max_y im_y] = max(y1y2);
vy(4) = max_y;
vx(4) = x1x2(im_y);

low_y = -400;
max_y = yo+R;
ii_b = find(tran_pxb<max_x & tran_pxb>low_x & tran_pyb <max_y &  tran_pyb>low_y);
ii_a = find(tran_pxa <max_x & tran_pxa > low_x & tran_pya < max_y &  tran_pya>low_y);
 [~, trivialbondt1,trivialbondt2] = adjacent(tran_pxa(ii_a),tran_pya(ii_a),tran_pxb(ii_b),tran_pyb(ii_b),3);
 x1 = 
x_c(:,1) = tran_pxa(ii_a(trivialbondt1))';
x_c(:,2) = tran_pxb(ii_b(trivialbondt2))';
y_c(:,1) = tran_pya(ii_a(trivialbondt1))';
y_c(:,2) = tran_pyb(ii_b(trivialbondt2))';
