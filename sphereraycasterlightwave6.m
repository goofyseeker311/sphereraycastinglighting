output_precision = 32; format long; pkg load statistics; pkg load matgeom;

function [k,n] = raysphereintersection(sphere, ray)
  spherecount = size(sphere,1); raycount = size(ray,1); kmatrix = zeros(raycount,2);
  for m=1:spherecount
    a = sum( ray(:,4:6).^2 ,2); sphererow = ones(raycount,1)*sphere(m,1:4); b = 2*sum( ray(:,4:6).*(ray(:,1:3)-sphererow(:,1:3)) ,2);
    c = sum( ray(:,1:3).^2 + sphererow(:,1:3).^2 - 2*ray(:,1:3).*sphererow(:,1:3) ,2) - ones(raycount,1)*(sphererow(1,4).^2);
    kv1 = (-b + sqrt(b.^2 - 4*a.*c)) ./ (2*a); k1 = [((((abs(imag(kv1))==0).*kv1)>0).*kv1) ones(raycount,1)*m];
    kv2 = (-b - sqrt(b.^2 - 4*a.*c)) ./ (2*a); k2 = [((((abs(imag(kv2))==0).*kv2)>0).*kv2) ones(raycount,1)*m];
    kmatrix = ((kmatrix(:,1)<=0).*(k1(:,1)>0).*k1).+(kmatrix(:,1)>0).*((k1(:,1)>0).*(((k1(:,1)<=kmatrix(:,1)).*k1).+((k1(:,1)>kmatrix(:,1)).*kmatrix)).+(k1(:,1)<=0).*kmatrix);
    kmatrix = ((kmatrix(:,1)<=0).*(k2(:,1)>0).*k2).+(kmatrix(:,1)>0).*((k2(:,1)>0).*(((k2(:,1)<=kmatrix(:,1)).*k2).+((k2(:,1)>kmatrix(:,1)).*kmatrix)).+(k2(:,1)<=0).*kmatrix);
  endfor
  k = kmatrix; n = [(kmatrix(:,1).*ray(:,4)).+ray(:,1) (kmatrix(:,1).*ray(:,5)).+ray(:,2) (kmatrix(:,1).*ray(:,6)).+ray(:,3)];
endfunction

pointlightscount = 5; pointlights = [100*(rand(pointlightscount,3)-0.5) 3.5*randi(255,pointlightscount,3)/255];
nspheres = 100; scenespheres = [100*(rand(nspheres,3)-0.5) 6*rand(nspheres,1)+4 randi(255,nspheres,3)/255];
cameravector = createLine3d(pi/2,0); cameraupvector = createLine3d(0,0); camerarightvector = createLine3d(pi/2,-pi/2);
spherebufferwidth = 3840; spherebufferheight = 1080; numbytes = 4; raycount = spherebufferwidth*spherebufferheight; horizontalstep = 360/(spherebufferwidth-1); verticalstep = 180/(spherebufferheight-1);
camerazbuffersphere = zeros(spherebufferwidth, spherebufferheight); camerargbabuffersphere = zeros(spherebufferwidth, spherebufferheight,numbytes);
ratio = pi/180; thetasteps = flip(-180:horizontalstep:180,2)*ratio; phisteps = (0:verticalstep:180)*ratio; rays = zeros(raycount,6); raycounter = 1;
for theta=thetasteps; for phi=phisteps; rays(raycounter,:) = createLine3d(phi,theta); raycounter++; endfor; endfor

[k, kn] = raysphereintersection(scenespheres(:,1:4), rays); kd=k(:,1);kb=k(:,2); m = max(kd);
kc=zeros(raycount,3);kci=find(kb);kbi=kb(kci);kc(kci,:)=scenespheres(kbi,5:7);kx=kn(kci,1);ky=kn(kci,2);kz=kn(kci,3);kcicount=size(kci,1);knc=kn(kci,:);
k1=reshape((1-(kd./max(abs(kd))).^2),spherebufferheight,spherebufferwidth); renderim=zeros(spherebufferheight,spherebufferwidth,3); 
kc1=reshape(kc(:,1),spherebufferheight,spherebufferwidth); kc2=reshape(kc(:,2),spherebufferheight,spherebufferwidth);
kc3=reshape(kc(:,3),spherebufferheight,spherebufferwidth); renderim(:,:,1) = kc1.*k1; renderim(:,:,2) = kc2.*k1; renderim(:,:,3) = kc3.*k1;

[X,Y,Z] = sphere(); sp = scenespheres; figure(1); clf; view(3); axis([-60 60 -60 60 -60 60]); hold on;
plot3(cameravector(1),cameravector(2),cameravector(3),'ro'); plot3(50*cameravector(4),50*cameravector(5),50*cameravector(6),'bx');
plot3([cameravector(1) 50*cameravector(4)],[cameravector(2) 50*cameravector(5)],[cameravector(3) 50*cameravector(6)],'g');
plot3([cameravector(1) 30*cameraupvector(4)],[cameravector(2) 30*cameraupvector(5)],[cameravector(3) 30*cameraupvector(6)],'b');
plot3([cameravector(1) 20*camerarightvector(4)],[cameravector(2) 20*camerarightvector(5)],[cameravector(3) 20*camerarightvector(6)],'r');
plot3(30*cameraupvector(4),30*cameraupvector(5),30*cameraupvector(6),'kx'); plot3(20*camerarightvector(4),20*camerarightvector(5),20*camerarightvector(6),'bx');
for u=1:nspheres; a=sp(u,4); x=sp(u,1); y=sp(u,2); z=sp(u,3); surf(a*X+x,a*Y+y,a*Z+z); endfor; plot3(kx,ky,kz,'y.'); plot3(pointlights(:,1),pointlights(:,2),pointlights(:,3),'rx'); hold off;
figure(2); imshow(renderim, "xdata", [-180 180], "ydata", [-90 90]); axis on; figure(3); hist(nonzeros(kd),30); imwrite(renderim,"testCCCA0.png");
kl=zeros(raycount,3); for v=1:pointlightscount; printf("v/plcount: %i/%i\n",v,pointlightscount); pla=(ones(kcicount,1)*pointlights(v,1:3)); plv=pla.-knc; pld=sqrt(sum(plv.^2,2)); pld3=pld*ones(1,3);
plnv=plv./pld3; lv=[knc.+(0.00000001.*plnv) plnv]; [ck, ckn] = raysphereintersection(scenespheres(:,1:4), lv); cki=find(ck(:,2)==0); ckicount=size(cki,1); ckni=ckn(cki,:); kcni=kci(cki,:); pld3n=pld3(cki,:);
plc=ones(ckicount,1)*pointlights(v,4:6); kl(kcni,:).+=(plc/pointlightscount); figure(1); hold on; plot3(ckni(:,1),ckni(:,2),ckni(:,3),'r.'); hold off; endfor;
kl1=reshape(kl(:,1),spherebufferheight,spherebufferwidth); kl2=reshape(kl(:,2),spherebufferheight,spherebufferwidth); kl3=reshape(kl(:,3),spherebufferheight,spherebufferwidth);
renderim2=zeros(spherebufferheight,spherebufferwidth,3); renderim2(:,:,1) = kc1.*kl1; renderim2(:,:,2) = kc2.*kl2; renderim2(:,:,3) = kc3.*kl3;
figure(4);clf; view(2); imshow(renderim2, "xdata", [-180 180], "ydata", [-90 90]); axis on; imwrite(renderim2,"testCCCB0.png");
