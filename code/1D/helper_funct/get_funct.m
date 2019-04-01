function [x,f,f_hat] = get_funct(N,funct,k,s)

% true function
if strcmp(funct,'sawtooth_m1_1')
    
    xs = -1; xe = 1;
    dx = (xe-xs)/N;
    x = xs + (0:N-1)'*dx;
    
    f = (x+1).*(x<=0) + (x-1).*(x>0);
    
    f_hat = 1./(k.^2).*(1i*(-sin(k)+k));
    if (sum(k==0) == 1)
        f_hat(k==0) = 0;
    end
    
elseif strcmp(funct, 'sawtooth_mpi_pi')
    
    xs = -pi; xe = pi;
    dx = (xe-xs)/N;
    x = xs + (0:N-1)'*dx;
    
    f = (-x-pi)./(2*pi).*(x<=0) + (pi-x)./(2*pi).*(x>0);
    
    f_hat = 1./(2*pi^2*k.^2).*(sin(pi*k)-k.*pi)*1i;
    if (sum(k==0) == 1)
        f_hat(k==0) = 0;
    end
    
elseif strcmp(funct,'step')
    
    xs = -1; xe = 1;
    dx = (xe-xs)/N;
    x = xs + (0:N-1)'*dx;
    
    f = 5.*(x<=0) + -5.*(x>0);
    
    f_hat = (1./k).*(-5*1i*(cos(k)-1));
    if (sum(k==0) == 1)
        f_hat(k==0) = 0;
    end
    
elseif strcmp(funct,'AG')
    
    xs = -1; xe = 1;
    dx = (xe-xs)/N;
    x = xs + (0:N-1)'*dx;
    
    f = cos(pi*x/2).*(x < -1/2) + cos(3*pi*x/2).*(-1/2 <= x & x <1/2) +...
        cos(7*pi*x/2).*(x <= 1 & x >= 1/2);
    
    % integrate[cos(pi*x/2)*exp(I*k*x),{x,-1,-1/2}] + integrate[exp(I*k*x)*cos(3*pi*x/2),{x,-1/2,1/2}] + integrate[cos(7*pi*x/2) * exp(I*k*x),{x,1/2,1}]
    
    f_hat = 1./(2.*(-4.*k.^2+pi^2)) .* (2.*exp(1i.*k).*pi - ...
        sqrt(2).*exp(1i.*k./2).*(2*1i.*k+pi)) + ... 
        1./(2*(-4.*k.^2+49*pi^2)).*(exp(-1i.*k).*(-14*pi+...
        sqrt(2)*exp(1i*k/2).*(2*1i.*k + 7*pi))) + ... 
        1./(-4*k.^2+9*pi^2).*(sqrt(2).*(3*pi*cos(k./2) + 2*k.*sin(k./2))); 
    if (sum(k==-pi/2))
        f_hat(k==-pi/2) = (8-8*1i+3*pi)/(24*pi^2);
    end
    
    if sum(sum(k==pi/2))
        f_hat(k==pi/2) = (8+8*1i+3*pi)/(24*pi);
    end
    
    if sum(sum(k==(-7*pi)/2))
        f_hat(k==(-7*pi)/2) = (184-40*1i+105*pi)/(840*pi); 
    end
    
    if sum(sum(k==(7*pi)/2))
        f_hat(k==(7*pi)/2) = (184+40*1i+105*pi)/(840*pi); 
    end
    
    if sum(sum(k==(-3*pi)/2))
        f_hat(k==(-3*pi)/2) = (-28+18*1i+15*pi)/(60*pi); 
    end
    
    if sum(sum(k==(3*pi)/2))
        f_hat(k==(3*pi)/2) = -(28+18*1i-15*pi)/(60*pi); 
    end
    
    
elseif strcmp(funct,'multijump')
    
    xs = -pi; xe = pi;
    dx = (xe-xs)/N;
    x = xs + (0:N-1)'*dx;
    
    f = 1.5*(x>=-3*pi/4).*(x<-pi/2) + ...
        ( 7/4 - x/2 + sin(x-1/4) ).*(x>=-pi/4).*(x<pi/8) + ...
        ( x*11/4 - 5 ).*(x>=3*pi/8).*(x<3*pi/4);
    
    % Fourier modes computed analytically using Mathematica
    f_hat = exp(3*pi*1i*k/4).*(3./(4*pi*1i*k)) + ...
        exp(1i*pi*k/2).*(-3./(4*pi*1i*k)) + ...
        exp(-3*1i*pi*k/4).*(5./(1i*k) - ...
        (33*pi)./(16*1i*k) -11./(4*1i*1i*k.*...
        k))/(2*pi) + exp(-3*1i*pi*k/8).*...
        (-5./(1i*k) + (33*pi)./(32*1i*k) + ...
        11./(4*1i*1i*k.*k))/(2*pi) + ...
        exp(-1i*pi*k/8).*(-7./(4*1i*k) + ...
        pi./(16*1i*k) + 1./(2*1i*1i*k.*k) ...
        - sin(pi/8 - 1/4)*(k./(1i*(k.*k-1)))...
        + cos(pi/8 - 1/4)*(1./(k.*k-1)) )/(2*pi) + ...
        exp(1i*pi*k/4).*(7./(4*1i*k) + ...
        pi./(8*1i*k) - 1./(2*1i*1i*k.*k) +...
        sin(-pi/4 - 1/4)*(k./(1i*(k.*k-1))) ...
        - cos(-pi/4 -1/4)*(1./(k.*k-1)) )/(2*pi);
    
    % The DC component
    if( sum(k==0) == 1 )
        f_hat(k==0) = ( (-pi/2+3*pi/4)*1.5 - ...
            5*(3*pi/4-3*pi/8) + 11*(9*pi*pi/16 - 9*pi*pi/64)/8 ...
            + 7*(pi/8+pi/4)/4 - .25*(pi*pi/64 - pi*pi/16) - ...
            cos(pi/8 - .25) + cos(-pi/4 -.25) )/(2*pi);
    end
    
    % Components at k=+/-1
    if( sum(k==1) == 1 )
        f_hat(k==1) = -( 48 - 160*(-1)^(1/8) - (56+16i)*...
            (-1)^(3/8) -88*(-1)^(5/8) + (100-60i)*sqrt(2) - 8*...
            sqrt( 28-45i ) + (4-4i)*((1-1i)+sqrt(2))*exp(1i/4) + ...
            (33*(-1)^(1/8)+2*(-1)^(3/8)-(35-35i)*sqrt(2))*pi + ...
            6i*pi*exp(-1i/4) )/(64*pi);
    end
    if( sum(k==-1) == 1 )
        f_hat(k==-1) = exp(-1i/4)*( (-4-4i)*((1+1i)+sqrt(2)) ...
            + 6i*exp(1i/2)*pi + exp(1i/4)*(-8*(6+(2+7i)*(-1)^(1/8) +...
            11*(-1)^(3/8) + 20*(-1)^(7/8) + (8+5i)*sqrt(2)) + ( 2*...
            (-1)^(5/8) +33*(-1)^(7/8) +(35+35i)*sqrt(2))*pi) )/(64*pi);
    end
    
elseif strcmp(funct,'Wolfgang')
    
    x = -pi + (2*pi)/N*(0:N-1)';
    
    % Fourier modes computed analytically using Mathematica
    f_hat = ( ( -1+exp(-pi*1i*(k-1i)/4) )./(-1-k*1i) ...
        + 1i*( exp(5*1i*k/2) - exp(1i*pi*k/4) )./...
        k + ( exp(-11*1i*k/4).*...
        ( cos(55/4)*1i*k - 5*sin(55/4) + ...
        exp(3*1i*k/2).*( -1i*cos(25/4)*k + ...
        5*sin(25/4) ) ) )./(2*(-25+k.^2)) + 2*sin(3*...
        k/4).*exp(-2*1i*k)./k )/(2*pi);
    
    % The DC component
    if( sum(k==0) == 1 )
        f_hat(k==0) = ( -1.5 + pi/4 - cosh(pi/4) + ...
            sinh(pi/4) + (15-sin(25/4)+sin(55/4))/10 )/(2*pi);
    end
    
    % Components at k=+/-5
    if( sum(k==5) == 1 )
        f_hat(k==5) = ( 3/8 - (1/26 - 5i/26)*(-1+(-1)^...
            (.75+1i/4))  + ((-1)^.75 + 1i*exp(25i/2))/5 + ...
            2*(sin(15/4)*exp(-10i))/5 + sin(15/2)*exp(-20i)/20 )/(2*pi);
    end
    if( sum(k==-5) == 1 )
        f_hat(k==-5) = ( (1/26 + 5i/26)*(1+(-1)^...
            (.25+1i/4)) + (-(-1)^.25 - 1i*exp(-25i/2))/5 + ...
            (2*sin(15/2)*exp(20i) + 16*sin(15/4)*exp(10i) +15)/40 )/(2*pi);
    end
    
    % And now, the physical space function
    f = -1*(x>=-2.5).*(x<-pi/4) + ...
        exp(-x).*(x>=0)...
        .*(x<pi/4) + ( 1 + 0.5*cos(5*x) ).*...
        (x>=1.25).*(x<2.75);
    
elseif strcmp(funct,'variation')
    
    x = -pi + (2*pi)/N*(0:N-1)';
    
    % Fourier modes computed analytically using Mathematica
    f_hat = ( 6*(1+exp(5i*k*pi/6))./(-36+k.^2) + ...
        (1-exp((-2-1i*k)*pi)+(-2-1i*k)*pi)...
        ./(pi*(k-2i).^2) )/(2*pi);
    
    % Components at k=+/-6
    if( sum(k==6) == 1 )
        f_hat(k==6) = ( 1/300 + 1i/400 )*( -(6+18i)*pi-...
            (30+40i)*pi^2 )/(pi^2);
    end
    if( sum(k==-6) == 1 )
        f_hat(k==-6) = ( (30+90i)*pi + 250i*pi^2 )/(1200*pi^2);
    end
    
    % And now, the physical space function
    f = sin(6*x).*(x>=-5*pi/6).*(x<0) + ...
        exp(-2*x).*((pi-x)/pi).*(x>=0);
    
    
elseif strcmp(funct,'binary')
    
    x = linspace(-1,1,N);
    ind = sort(randperm(N,s));
    data_height = ones(s,1);
    f = zeros(N,1);
    f(ind) = data_height;
    
elseif strcmp(funct,'uniform')
    
    x = linspace(-1,1,N);
    % indices
    ind = sort(randperm(N,s));
    data_height = rand(s,1);
    f = zeros(N,1);
    f(ind) = data_height;
    
    
elseif strcmp(funct,'gaussian')
    
    x = linspace(-1,1,N);
    % indices
    ind = sort(randperm(N,s));
    data_height = randn(s,1);
    f = zeros(N,1);
    f(ind) = data_height;
    
elseif strcmp('SinCosLin')
    
    x = -pi + (2*pi)/N*(0:N-1)';
    % Fourier modes computed analytically using Mathematica
    f_hat = ( ((-1).^k-1i*k.*exp(1i*pi*k/2))...
        ./(-1+k.^2) + (((-1).^k).*...
        (-4+exp(3i*k*pi/4).*(4-3i*pi*k)))./(4*k.^2) + ...
        (exp(-1i*k*pi/4).*(sqrt(2)*exp(3i*k*pi/4).*...
        (3-2i*k) + 6*cos(pi/8) - 4i*k*sin(pi/8)...
        ))./(-9+4*k.^2) )/(2*pi);
    
    % The DC component
    if( sum(k==0) == 1 )
        f_hat(k==0) = ( 27*pi^2 - 32*(3+sqrt(2)+2*cos(pi/8)) )...
            /(192*pi);
    end
    
    % Components at k=+/-1
    if( sum(k==1) == 1 )
        f_hat(k==1) = ( 30+20*(-1)^(5/8) + 4*(-1)^(7/8) + ...
            (2-22i)*sqrt(2) - 5*pi*(1i+3*(-1)^(1/4)) )/(40*pi);
    end
    if( sum(k==-1) == 1 )
        f_hat(k==-1) = ( 60-8*(-1)^(1/8) - 40*(-1)^(3/8) + ...
            (4+44i)*sqrt(2) + (5+5i)*pi*(1+1i+3i*sqrt(2)) )/(80*pi);
    end
    
    % And now, the physical space function
    f = sin(x).*(x>=-pi).*(x<-pi/2) - ...
        cos(3*x/2).*(x>=-pi/2).*(x<pi/4) + ...
        (pi-x).*(x>=pi/4).*(x<=pi);
    
end


end

