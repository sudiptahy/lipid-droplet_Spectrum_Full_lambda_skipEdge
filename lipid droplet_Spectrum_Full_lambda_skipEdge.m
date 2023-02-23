nr=10; % number of rows after which write the images
ld=1; % ld number
c=16; % number of wavelength channel
z0=1; % plane 1
z1=37; % plane 13
i_med=7; % kernel size for median filter keep it odd

%% fitting paramaters
% equation  'y = a1+ b1*exp(-((x-c1)^2)/d1)', 'independent', 'x', 'dependent', 'y' );
I_P = [0 0.54 15 sqrt(c)]; % initialization parameter of a1, b1, c1, d1
I_P_L = [-1 0 0 0]; % lower bound for [a1, b1, c1, d1]
I_P_U=  [01 1 c c];% upper bound for [a1, b1, c1, d1]

%% wavelength setting
l0=405; % wavelength of first image
dl=8.9 ; % wavelength resolution (difference between sucessive image wavelength)

for z=z0:z1
im_info = imfinfo('Image 7.tif');% return tiff structure, one element per image
F_name=im_info(1).Filename;
N_I=size(im_info, 1); %  number of images
path_name =  strsplit(F_name,'\');
path_name=flip(path_name);
Im_name= path_name{1}

i0=(z-1)*c+1; % first image index  in the stack 
i1=z*c;  % last image index in the syper stack
I_stack=[];

% correction for edge zero
% for odd i_med, there are (i_med-1)/ 2  zeros in the edges
d=(i_med-1)/ 2; % d is the size of the 0 edge to be replaces

C=im_info(1).Width;
R=im_info(1).Height;

i_skp=1;
index_skip=[];
for ir=1:d  % for left top
    for ic=1:d
 index_skip(i_skp)=sub2ind([R C], ir,ic);
 i_skp=i_skp+1;
    end
end

for ir=R-d+1:R % for Left bottom
    for ic=1:d
 index_skip(i_skp)=sub2ind([R C], ir,ic);
 i_skp=i_skp+1;
    end
end


for ir=1:d % for Right top
    for ic=C-d+1:C
 index_skip(i_skp)=sub2ind([R C], ir,ic);
 i_skp=i_skp+1;
    end
end

for ir=R-d+1:R % for Right bottom
    for ic=C-d+1:C
 index_skip(i_skp)=sub2ind([R C], ir,ic); % list of inear idex for skipping the fitting
 i_skp=i_skp+1;
    end
end
ir=[]; ic=[];

for i=i0:i1  % for reading images to make a stck
I_temp = double(imread(Im_name, i)); % read the c channel of given z

%IedgeLT=I_temp(1:d, 1:d);
% IedgeLB=I_temp(R-d+1:R , 1:d);
% IedgeRT=I_temp(1:d, C-d+1:C);
% IedgeRB=I_temp(R-d+1:R , C-d+1:C);

I_temp=medfilt2(I_temp,[i_med i_med]);

% I_temp(1:d, 1:d)=IedgeLT; % correction for left top edge zero after filtering 
% I_temp(R-d+1:R , 1:d)=IedgeLB;% correction for edge zero after filtering
% I_temp(1:d, C-d+1:C)=IedgeRT;% correction for edge zero after filtering
% I_temp(R-d+1:R , C-d+1:C)=IedgeRB;% correction for edge zero after filtering

 
I_stack =cat(3 , I_stack, I_temp);
end

Itemp=I_stack(:,:, round(c/2)); % select the image with best brightness

imwrite (uint16(Itemp), strcat('Intensity',num2str(ld), '_' ,Im_name));
Itemp=[];

%% find the list of peak position
 
 % %Itemp=wiener2(Itemp,[3 3]);
% Itemp = medfilt2(Itemp,[i_med i_med]);

%             if z==z0
%             imshow(Itemp, [min(min(Itemp)), max(max(Itemp))])
%             [x0,y0] = ginput;  % ginput can be used to interactively find the LD centre
%             y1=double(uint16(y0));
%             x1=double(uint16(x0));
% 
%             prompt = {'Number of pixels from LD centre:'};
%             title = 'Input';
%             dims = [1 35];
%             definput0 = {num2str(N)};
%             answer1 = inputdlg(prompt,title,dims,definput0);
%             N = str2num(answer1{1});
%             end
% x1 and y1 are the coordinates of the centre of the object

% [n_p, x]=size(x1); % n_p is the number of peaks /LDs
% N_theta=32;
% n_per=1;


    
   %%%%%%%%%%%%%%%%%%%%%%%%%


Peri_in= uint8(zeros(R, C)); % blank peri image
I_Peak= double(zeros(R, C));  % image of peak obtained from gaussian fit
I_Width=double(zeros(R, C)); % image of width obtained from gaussian fit
I_GOF= double(zeros(R, C)); % image of GOF goodness of fit

mask2=logical(zeros(R,C)+1);

   
   mask=mask2; % mask= imread ('MASK.tif');
    
   idx1 = find(mask == 1);  % for calculating spectrum on individual pixel
  %--------- maks mask from Peri for spectrum calculation
    Av_spectrum=[];
  Pix_spectrum=[];  
  %%%%%%%%%%%%%%%% if segmention is fine %%%%%%%%%%%

  
          for i=1:c  % for evaluating spectrum using the mask
            Itemp1=[];
            Itemp=I_stack(:,:,i);  %%% here give the image for which the spectrum to be evaluated
             Pix_spectrum=[Pix_spectrum ; Itemp(idx1)' ]; % each coulumn is apectrum of pixel id 'idx1'
          end
          
S = sum(Pix_spectrum,1); % for normalization of spectrum 
[rN, cN] = size (S);  % sums the number along the row
XXX = [1:c];  % x-axis or wavelength for fitting in 'image number' unit
 % XXXm=l0+ (XXX-1)*dl; % x-axis or wavelength for fitting in absolute unit

%% for Normalising the spectrum and fitting

            for i=1:cN  % for Normalising the spectrum and fitting

             % [r,c]=ind2sub([R C],i) for skipping edges
%              index_skip(i_skp)=sub2ind([R C], r,c); % list of inear idex for skipping the fitting
%              

                if (isempty(find (index_skip==i)))
                 Pix_spectrum_N(:,i)=Pix_spectrum(:,i)/S(i); % normalixed spectrum in coulmn i
               YYY=Pix_spectrum_N(:,i); % Yaxis for fitting 
                YYY2=smooth(YYY, 3); % smoothing the spectrum
              
              % equation  'y = a1+ b1*exp(-((x-c1)^2)/d1)', 'independent', 'x', 'dependent', 'y' );
               % I_P = [0 0.54 15 0.15]; % initialization parameter of a1, b1, c1, d1
               I_P(1)=min(YYY2);
               I_P(2)=max(YYY2);
               I_P(3)= XXX(find (YYY2==max(YYY2),1)); % find the image number in lamda unit
             
              [idx]= find (I_P==0);
               I_P(idx)=0.5; % to keep the lower and upper bound different
               %clear idx ;
               
                I_P_L=I_P - 50*I_P; % 10 % lower lower than I_P
                I_P_U=I_P + 50*I_P; % 10 % lower lower than I_P
                I_P_L(2) = 0;
                I_P_L(3) =0;
                I_P_L(4)=0;
                              
               [fitresult, gof] = createFitGaussian(XXX', YYY, I_P, I_P_L, I_P_U);
            I_Peak(idx1(i))= fitresult.c1;  % image of peak obtained from gaussian fit
            I_Width(idx1(i))=(fitresult.d1); % image of width obtained from gaussian fit
            I_GOF(idx1(i))= gof.rsquare; % image of GOF goodness of fit
            
                else
            I_Peak(idx1(i))= 0;
            I_Width(idx1(i))=0;
            I_GOF(idx1(i))= 0;
                
                end % end of if for fitting
            
            %% write image mid way with n coulmn
            if (rem (i, nr*R)==0 );
               Image = uint16 (1000*I_Peak); %%need to put wavelength values%%
                imwrite (Image,strcat('Peak_into1000_',num2str(ld),'_',Im_name));
                
                 Image1 =  (I_Peak); %%need to put wavelength values%%
                 Image1= uint16((l0 + (Image1-1)*dl));
                 imwrite (Image1,strcat('Peak_lambda_',num2str(ld),'_',Im_name));
                
                ImageW = uint16 (I_Width*100); %%what is the unit of width%%
                imwrite (ImageW,strcat('Width_into100_',num2str(ld),'_',Im_name));
                
                 ImageW1 =  (I_Width); %%need to put wavelength values%%
                 ImageW1= uint16(ImageW1*dl);
                 imwrite (ImageW1,strcat('Width_lambda_',num2str(ld),'_',Im_name));
                 
 
                ImageG = uint16 (I_GOF*100); %%GOF in percentage%%
                imwrite (ImageG,strcat('GOF_',num2str(ld),'_',Im_name));
               
            else
                if (i==cN)
                Image = uint16 (1000*I_Peak); %%need to put wavelength values%%
                imwrite (Image,strcat('Peak_into1000_',num2str(ld),'_',Im_name));
                
                 Image1 =  (I_Peak); %%need to put wavelength values%%
                 Image1= uint16((l0 + (Image1-1)*dl));
                 imwrite (Image1,strcat('Peak_lambda_',num2str(ld),'_',Im_name));
                
                ImageW = uint16 (I_Width*100); %%what is the unit of width%%
                imwrite (ImageW,strcat('Width_into100_',num2str(ld),'_',Im_name));
                
                 ImageW1 =  (I_Width); %%need to put wavelength values%%
                 ImageW1= uint16(ImageW1*dl);
                 imwrite (ImageW1,strcat('Width_lambda_',num2str(ld),'_',Im_name));
                 
 
                ImageG = uint16 (I_GOF*100); %%GOF in percentage%%
                imwrite (ImageG,strcat('GOF_',num2str(ld),'_',Im_name));
                    
                end
            end

            end

% figure
% A=(I_Peak(mask==1));
% l=min(A);
% h=max(A);
% imshow (I_Peak, [l, h])
% colormap jet


%%%%%%%%%%%%  write image file to disc
%          if (i==i0)
%          imwrite(I_temp, strcat('Spect_', Im_name),'WriteMode','overwrite');
%          else
%          imwrite(I_temp, strcat('Spect_', Im_name),'WriteMode','append');
%          end 

 
% i1 (end of the spectrum) in the stack


end

% for iSlice=1:1:c
%     imwrite (uint16(I_sub(:,:,iSlice)),strcat('Intensity',num2str(ld),'_',Im_name),'WriteMode','append');
% end
