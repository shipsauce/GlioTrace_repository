function [D1,D2,D3,delta]=register_stack_DT(D1,D2,D3)
% This function registers subsequent images of a
% stack using the net shift in the x- and y-axis estimated by the function
% dftregistration.
%
% @authors: Sven Nelander
% @date: 140624

sigmah=8; % typical cell size
hsizeh = sigmah*5;
h = fspecial('log', hsizeh, sigmah);
delta=[0 0];
for i=2:size(D1,3)
    
    correlation1=corr2(D1(:,:,i-1),D1(:,:,i));
    imgA=D1(:,:,i-1);
    imgB=D1(:,:,i);
    imgX=D2(:,:,i-1);
    imgY=D2(:,:,i);
    imgM=D3(:,:,i-1);
    imgN=D3(:,:,i);

    output = dftregistration(fft2(imgA+imgX+imgM),fft2(imgB+imgY+imgN),100);
    dx=output(4);
    dy=output(3);
    delta(i,:)=[dx dy];
    D1(:,:,i)=imtranslate(D1(:,:,i),[dx dy],'Outputview','same');
    D2(:,:,i)=imtranslate(D2(:,:,i),[dx dy],'Outputview','same');
    D3(:,:,i)=imtranslate(D3(:,:,i),[dx dy],'Outputview','same');
    correlation2=corr2(D1(:,:,i-1),D1(:,:,i));
    
    if(mod(i,5)==0)
    fprintf('frame: %g before: %f after %f\n', [i correlation1 correlation2]);
    end
    
end