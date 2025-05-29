function [  x_horizontal_strip,y_lower_horizontal_strip,...
            y_upper_horizontal_strip,...
            x_sector,y_lower_sector,y_upper_sector,...
            x_disk,y_disk,...
            x_alphav,y_alphav,...
            x_betav,y_betav ] = coordinates_D_regions(...
            eigenvalues_poles,alpha_v,beta_v,theta_s,r_d,q_d,w_H,e_P)

%A) Instantiate coordinate vectors for boundary of Horizontal Strip D-region   
x_horizontal_strip=[];
y_lower_horizontal_strip=[];
y_upper_horizontal_strip=[];

%B) Instantiate coordinate vectors for boundary of Sector D-region
x_sector=[] ;
y_lower_sector=[] ;
y_upper_sector=[] ;

%C) Instantiate coordinate vectors for boundary of Disk D-region
x_disk=[] ;
y_disk=[] ;

%D) Instantiate coordinate vectors for boundary of Vertical Strip D-region
x_alphav=[] ;
y_alphav=[] ;
x_betav=[] ;
y_betav=[] ;            
            
%E) Compute the minimum and maximum real part of the eigenvalue-poles
min_real_aut=min(real(eigenvalues_poles)); 
max_real_aut=max(real(eigenvalues_poles));                

%F) Compute the minimum and maximum real part of Disk Region
min_real_disk=-q_d-r_d;
max_real_disk=-q_d+r_d;

%G) Compute the minimum real part 
min_real=min([min_real_aut,min_real_disk,-beta_v]);

%H) Compute the maximum real part
max_real=max([max_real_aut,max_real_disk,-alpha_v]);
max_real=min([max_real,0]);

%I) Compute the maximum imaginary part of the eigenvalue-poles input
max_imag_aut=max(imag(eigenvalues_poles));

%J) Compute the minimum imaginary part of the Sector Region
max_imag_sector=-min_real*tan(theta_s);

%K) Compute the minimum imaginary part 
max_imag=max([max_imag_aut,max_imag_sector,w_H,r_d]);

%M) Construct vectors for boundaries of the D-regions
for i=1:100
    
    xt=min_real+(i-1)*(max_real-min_real)/99;
        
    if ~isempty(w_H)
        x_horizontal_strip(i,1)=xt;
        y_lower_horizontal_strip(i,1)=-w_H;
        y_upper_horizontal_strip(i,1)=w_H;
    end
    
    if ~isempty(theta_s)
        x_sector(i,1)=xt;
        y_upper_sector(i,1)=xt*tan(theta_s);
        y_lower_sector(i,1)=xt*tan(-theta_s);
    end
    
    if ~isempty(r_d)
        if isempty(q_d)
            x_disk(i,1) = r_d*cos(2*pi*i/100);
            y_disk(i,1) = r_d*sin(2*pi*i/100);
        else
            x_disk(i,1) = -q_d+r_d*cos(2*pi*i/100);
            y_disk(i,1) = -q_d+r_d*sin(2*pi*i/100);
        end
    end
    
    if ~isempty(alpha_v)
        x_alphav(i,1) = -alpha_v;
        y_alphav(i,1) = -max_imag+2*max_imag*(i-1)/99;
    end
    
    if ~isempty(beta_v)
        x_betav(i,1) = -beta_v;
        y_betav(i,1) = -max_imag+2*max_imag*(i-1)/99;
    end
    
end

end

