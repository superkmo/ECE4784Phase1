%Creator: Kevin Mohan
%Date Last Modified: 10/01/14
%Program Functionality: this program was made to integrate various
%components of a cell into an electrical functional model.

%%%%%%%%%%%%%%%%%%%%% PART 1: Defining Constants %%%%%%%%%%%%%%%%%%%%
gK = 36; %[mS/cm^2]
gNa = 120; %[mS/cm^2]
gL = 0.3; %[mS/cm^2]
EK = -12; %[mV]
ENa = 115; %[mV]
EL = 10.6; %[mV]
Vrest = -70; %[mV]
Cm = 1; %[uF/cm^2]

%%%%%%%%%%%%%%%%%% PART 2: Setting Variables %%%%%%%%%%%%%%%%%

%%Simulation Parameters
Vm = 0; %[mV] offset from -70mV after end of simulation
I = 0;
nstep = 100; %number of steps
hstep = 1/nstep; %step size
t = 0:0.1*hstep:0.1; %milliseconds

%%Gating Variables:
alpham = 0.1*((25-Vm)/(exp((25-Vm)/10)-1));
betam = 4*exp(-1*Vm/18);
alphan = 0.01*((10-Vm)/(exp((10-Vm)/10)-1));
betan = 0.125*exp(-1*Vm/80);
alphah = 0.07*exp(-1*Vm/20);
betah = 1/(exp((30-Vm)/10)+1);

mo = alpham/(alpham+betam);
no = alphan/(alphan+betan);
ho = alphah/(alphah+betah);
m(1) = mo;
n(1) = no;
h(1) = ho;
%%Currents
% INa = m.^3*gNa.*h*(Vm-ENa);
% IK = n.^4*gK*(Vm-EK);
% IL = gL.*(Vm-EL);
% Iion = I-IK-INa-IL;
% Iion2 = Vm*Cm;

%%Derivatives
% dVmdt = Iion/Cm;
% dmdt = alpham*(1-m)-betam*m;
% dndt = alphan*(1-n)-betan*n;
% dhdt = alphah*(1-h)-betah*h;

%%%%%%%%%%%%%%%%%% PART 3: Euler's Method %%%%%%%%%%%%%%%%%
%Euler's Method states that
%y'(t) = f(t,y(t)) where y(t_0) = y_0
%y(t_n) = y_0 + h*f(t_n-1,y_n-1)
%n = (sim_duration)/h

for i=1:nstep %euler method, begin iteration
    
    alpham(i+1) = 0.1.*((25-Vm(i))/(exp((25-Vm(i))./10)-1));
    betam(i+1) = 4*exp(-1.*Vm(i)/18);
    alphan(i+1) = 0.01*((10-Vm(i))/(exp((10-Vm(i))/10)-1));
    betan(i+1) = 0.125*exp(-1.*Vm(i)/80);
    alphah(i+1) = 0.07*exp(-1.*Vm(i)/20);
    betah(i+1) = 1/(exp((30-Vm(i))/10)+1);
    
    INa = m(i).^3*gNa.*h(i).*(Vm(i)-ENa);
    IK = n(i).^4*gK.*(Vm(i)-EK);
    IL = gL.*(Vm(i)-EL);
    Iion = I-IK-INa-IL;
    
    dVmdt(i+1) = Iion/Cm;
    dmdt(i+1) = alpham(i).*(1-m(i))-betam(i).*m(i);
    dndt(i+1) = alphan(i).*(1-n(i))-betan(i).*n(i);
    dhdt(i+1) = alphah(i).*(1-h(i))-betah(i).*h(i);
    
    m(i+1) = m(i) + hstep.*dmdt(i);
    n(i+1) = n(i) + hstep.*dndt(i);
    h(i+1) = h(i) + hstep.*dndt(i);
    
    Vm(i+1) = Vm(i) + hstep.*dVmdt(i);
    
end

Vm = Vm - 70; %[mV] scaling back data by 70 mV offset

