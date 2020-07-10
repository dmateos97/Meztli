%% Method of characteristics Nozzle design code
%%Assumptions: 
%%Adiabatic & isentropic (Reversible process)
%%2-D supersonic axisymmetic flow from M=1 at throat to M=Me at exit
%%Irrotational flow (du/dr=dv/dx)

%%This code uses the prandtl meyer expansion function to find the maximum
%%angle by which a given exhaust flow can turn at the throat for a minimum
%%length nozzle design. The method of characteristics is then used to
%%parametrize the expansion fan and find the axisymmetric flow deflections
%%and new slopes in order to find the optimized nozzle wall contour. This
%%method delivers axial exit flow with uniform mach number for given Pc/Pe
%%and specific heat ratio inputs. This code consists of 5 sections which
%%perform the following functionalities:

%%Section 1: Input section
    %%Creates GUI for modular user input paramaters

%%Section 2: Preliminary paramater calculations
    %%Uses nozzle thermodynamic relations to find the exit mach number for
    %%given combustion chamber conditions, and finds the maximum vacuum
    %%thrust coefficient Cfv for an ideal and reversible vacuum expansion
    %%process. Given user input thrust a vacuum optimized throat radius is
    %%found (Used if rt = 0 in section 1)
   
%%Section 3: Expansion fan
    %%The Prandtl-Meyer function is used to find and plot the expansion fan
    
%%Section 4: Flow deflections
    %%The MOC is used to find and plot the intersections of the left
    %%running and right running characteristics and the flow deflections at
    %%each intersection of the crossing expansion fans.

%%Section 5: Wall intersections
    %%The wall intersetcions are found given the final slopes computed in
    %%section 4 & the wall slopes starting with theta max found as a
    %%function of the prantl-meyer function. The wall slopes iterate n
    %%times until they reach 0 at the nozzle exit to ensure an evenly
    %%distributed flow. The results are plotted to provide the nozzle
    %%contour.
     
%% Input section
clear; clc; clearvars; 
Ru = 8314; %Universal gas constant (J/molK)

%MOC critical inputs (Defaults for LOX-CH4 10atm LRE)
input = inputdlg({'Enter combustion chamber pressure (Pa) ',...
    'Enter ambient pressure (Pa) ','Enter specific heat ratio ',...
    'Enter throat radius (mm) ','Enter number of characteristics '},'inputs');
     
P_c = str2double(cell2mat(input(1))); 
P_amb = str2double(cell2mat(input(2,1))); Pe = P_amb;
gamma = str2double(cell2mat(input(3,1)));
rt = str2double(cell2mat(input(4,1)));
num = str2double(cell2mat(input(5,1)));
% %Defaults:
% P_c = 2e6; %Combustion Chamber pressure
% P_amb = 101325; %Ambient pressure at optimized altitude (sea level)
% Pe = P_amb; %Nozzle exit pressure (Optimized for sea level)
% gamma = 1.1639; %Gamma - Specific heat capacity ratio (cp/cv)
% rt = 0.01647764*1000; %Throat radius (mm)
% num = 44; %number of characteristics/iterations


%Non critical user defined inputs (Default for LOX-CH4 10atm LRE)
T = 618; %Deisgn thrust
T_c = 3172.4465; %Adiabatic flame temperature (e.g RPA software/Test data)
R = 428.1; %Exhaust gas constant (J/kgK)
%hc = 0; he = 0; %ve2 = sqrt(2*(hc-he)); %check with ve to confirm

%% Preliminary thrust chamber parameters
IGR = (Pe/P_c)^((gamma-1)/gamma); %Ideal gas ratio
Te = T_c*IGR; %Nozzle exit temperature
Pc_Pt_RATIO = ((gamma+1)/2)^(gamma/(gamma-1)); %Throat choking condition (M = 1)
Pc_Pe_ratio = P_c/Pe; %Pressure expansion ratio
Pt = P_c/Pc_Pt_RATIO; %Throat pressure
Tt = 2*T_c/(gamma+1); %Throat temperature
Me = sqrt(2*(1-IGR)/((gamma-1)*IGR)); %Mach number at nozzle exit
Beta = (gamma+1)/(2*(gamma-1)); %Constant f(Gamma)
Cfv = ((2/(gamma+1))^Beta)*((gamma*Me+(1/Me))/sqrt(1+((gamma-1)/2)*Me^2)); %Vacuum optimized Cf
vt = sqrt(2*gamma*R*T_c/(gamma+1)); %throat velocity

%%Redundant Mach number check:
ve = sqrt((2*gamma/(gamma-1))*R*T_c*(1-(Pe/P_c)^((gamma-1)/gamma))); %Exit velocity - can use Ru/M instead of R
ae = ve/Me; %Exit speed of sound - f(Me)
ae2 = sqrt(gamma*R*Te); %Ideal exit speed of sound check - Not f(Me)
Me2 = ve/ae2; %Mach number check - f(ae2)

if rt == 0
    At = T/(Cfv*P_c); %Optimized for an ideal vacuum thrust coefficient
    rt = sqrt(At/pi)*1000; %Throat radius to maximize thrust for vacuum expansion
end

%% Method of characteristics - Expansion fan
clf;
A = (gamma+1)/(gamma-1); %Constant f(Gamma)
B = (gamma-1)/(gamma+1);  %Constant f(Gamma)
vPM = @(x)(sqrt(A)*atan(sqrt(B*(x^2-1)))-atan(sqrt(x^2-1))); %Prantl-Meyer Function - rad 
Theta_max = (1/2)*vPM(Me)*180/pi; %Max diverging expansion angle - deg
Mt = 1; M = Mt;
%InvPrandtlMeyer(Theta_max*2); %Function to check Me validity as
%f(Theta_max) (Must dl)

%Initialize vectors:
Mu = zeros(1,num);
Theta = zeros(1,num);
centerline = zeros(1,num);
RRC = zeros(1,num);
LRC = zeros(1,num);
RRC1 = zeros(1,num);
LRC1 = zeros(1,num);
dtheta = zeros(1,num);

Y_int = ones(1,num)*rt;
dt = 1-(Theta_max-fix(Theta_max));
for i=1:num
    dtheta(i) = (2*Theta_max/num)*(i-1); %Change in flow angle - chosen - degrees
    Theta(i) = (2*(Theta_max)-dtheta(i))*pi/180; %rad
    if i==1
        Theta(1) = 2*(Theta_max)*pi/180; %set outermost characteristic flow angle
    end
    SS_condition = [1 1.001*Me]; %Set choked flow (M=1) to Supersonic (M=Me) desired nozzle condition
    ExpansionAngle = @(x)(Theta(i)-vPM(x)); %Flow angle and expansion fan relation
    M(i) = fzero(ExpansionAngle,SS_condition); %Solve for Mach number as f(Theta,vPM)
    Mu(i) = (asin(1/M(i))); %Mach angle - rad
    centerline(i) = rt*tan(Theta(i)); %thrust axis points
    RRC(i) = -rt/centerline(i); %Right running slopes
          
    LRC(i) = -RRC(i); %Left running slopes
    
    if centerline(i) <= 0
        centerline(i) = rt*tan(89*pi/180);
    end
    
    X = [centerline(i) 0];
    Y = [0 Y_int(i)];
    plot(X,Y,'k')
    hold on
    
end


%% Method of Characteristcs - Flow deflections

%Initialize vectors:
ddeflect = zeros(1,num-1);
Add_Slopes = zeros(1,num-1);

x3plot = zeros(1,num-1);
y3plot = zeros(1,num-1);
x3mat = zeros(num-1);
y3mat = zeros(num-1);

new_theta = zeros(1,num-1);
new_slope = zeros(1,num-1);
deflSlopes = zeros(1,num-1);

%Find the flow deflections and new slopes at each intersection:
for j=1:length(centerline)-1
    for k = 1:j
        if k==1
            a = j;
            yint = -rt;
            LRC_new = LRC(j+1);
        else
            a = j-1*(k-1);
            yint = new_y-new_slope*x3;
            LRC_new = new_slope;
        end
        
        ddeflect(k) = (1/2*(RRC(a)+LRC_new))*pi/180; 
        syms x3
        yR = RRC(a)*x3+rt;
        yL = LRC_new*x3+yint;
        x3 = solve(yR==yL,x3); x3 = double(x3);
        x3plot(k) = x3;
        y3 = RRC(a)*x3+rt;
        y3plot(k) = y3;
        new_y = y3;
        x3mat(j,k) = x3;
        y3mat(j,k) = y3;
        new_theta(k) = atan(LRC_new)-ddeflect(k); %rad
        new_slope = tan(new_theta(k));         
    end
deflSlopes(j) = new_slope;
x3plot2 = [centerline(j+1) x3plot(1:j)];
y3plot2 = [0 y3plot(1:j)];
plot(x3plot2,y3plot2,'b')
end


%% Method of Characteristics - Wall intersections

%Find the wall intersections and contour
deflSlopes = deflSlopes(1:num-1);
deflSlopes2 = [LRC(1) deflSlopes(1:num-1)];
bb = zeros(1,num-1);
for m = 1:num-1
    xp = x3mat(m,m);
    yp = y3mat(m,m);
    bb(m) = yp-xp*deflSlopes(m);
end

b = [-rt bb];
WallSlope = tan(Theta_max*pi/180);
NewWallAngle = zeros(1,num);
x_Wallplot = zeros(1,num);
y_Wallplot = zeros(1,num);

X_Me = [centerline(1) diag(x3mat)'];
Y_Me = [0 diag(y3mat)'];
new_X = 0;
for o = 1:num
    NewWallAngle(o) = Theta_max-(Theta_max/(num-1))*(o-1); %deg
    WallSlope(o) = tan(NewWallAngle(o)*pi/180); 
    if o == 1
        new_Y = rt;
        contX = 0;
        contY = rt;
    else
        new_Y = yWall-WallSlope(o)*xWall;
        new_X = xWall;
        contX = new_X;
        contY = yWall;
    end
    syms xWall
    L_wall = WallSlope(o)*xWall+new_Y; 
    L_exp = deflSlopes2(end-(o-1))*xWall+b(end-(o-1)); 
    xWall = solve(L_wall==L_exp,xWall); xWall = double(xWall);
    yWall = WallSlope(o)*xWall+new_Y;
    y_wallCheck = deflSlopes2(end-(o-1))*xWall+b(end-(o-1)); 
    
    x_Wallplot(o) = xWall;
    y_Wallplot(o) = yWall;
    MachX = [X_Me(end-(o-1)) xWall];
    MachY = [Y_Me(end-(o-1)) yWall];
    plot(MachX, MachY,'b')
    contourX = [contX x_Wallplot(o)];
    contourY = [contY y_Wallplot(o)];
    plot(contourX, contourY, 'b')
    title('2-D Minimum length bell nozzle contour')
    xlabel('Thrust axis (mm)')
    ylabel('Radial axis (mm)') 
end

%Nozle wall contour x y points:
clearvars noz
noz(:,1) = [0; x_Wallplot'; xWall];
noz(:,2) = [rt; y_Wallplot'; yWall]; 

epsilon = yWall/rt; %Expansion area ratio
Cf = Cfv-(1/Pc_Pe_ratio)*epsilon;
fprintf('Pressure ratio (Pc/Pe) = %2.4f \tVacuum thrust coeff. = %2.4f \t     Thrust coeff. = %2.4f \n', Pc_Pe_ratio,Cfv,Cf)
fprintf('Area ratio (Ae/At) = %2.4f \t    Throat Temperature = %2.4f K \t Throat velocity = %2.4f m/s \n',epsilon,Tt,vt);
fprintf('Mach exit number = %2.4f \t        Exit Temperature = %2.4f K \t     Exit velocity = %2.4f m/s \n',Me,Te,ve);


%% References
%%Rocket propulsion (Heister, Pourpoint)
%%MIT courseware - Intro to rocket propulsion (Manuel Martinez-Sanchez) 
%%https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-512-rocket-propulsion-fall-2005/lecture-notes/lecture_4_5.pdf
%%http://www.ltas-aea.ulg.ac.be/cms/uploads/Aerothermodynamics05.pdf 
