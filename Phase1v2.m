%Creator: Kevin Mohan
%Date Last Modified: 10/01/14
%Program Functionality: this program was made to integrate various
%components of a cell into an electrical functional model.

%%%%%%%%%%%%%%%%%%%%% PART 1: Defining Constants %%%%%%%%%%%%%%%%%%%%
%A bit of formatting first
clear;
clc;
gKbar = 36; %[mS/cm^2]
gNabar = 120; %[mS/cm^2]
gLbar = 0.3; %[mS/cm^2]
gK(1) = gKbar;
gNa(1) = gNabar;
gL(1) = gLbar;
EK = -12; %[mV]
ENa = 115; %[mV]
EL = 10.6; %[mV]
Vrest = -70; %[mV]
Cm = 1; %[uF/cm^2]

%%%%%%%%%%%%%%%%%% PART 2: Setting Variables %%%%%%%%%%%%%%%%%

%%Simulation Parameters
Vm(1) = 0; %[mV] offset from -70mV after end of simulation
nstep = input('Please enter the number of iterations of the Euler Method: '); %number of steps
hstep = 1/nstep; %step size
t = 0:0.1*hstep:0.1; %milliseconds
I_mag = input('Please enter the pulse current (uA/cm^2): ');
I_dur = input('Please enter the duration of this pulse in milliseconds: ');
I = zeros(nstep,1);
for r = 1:I_dur/nstep
    I(r) = I_mag;
end
%%Gating Variables:
alpham = 0.1*((25-Vm)/(exp((25-Vm)/10)-1));
betam = 4*exp(-1*Vm/18);
alphan = 0.01*((10-Vm)/(exp((10-Vm)/10)-1));
betan = 0.125*exp(-1*Vm/80);
alphah = 0.07*exp(-1*Vm/20);
betah = 1/(exp((30-Vm)/10)+1);

%First run values
m(1) = alpham/(alpham+betam);
n(1) = alphan/(alphan+betan);
h(1) = alphah/(alphah+betah);

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
    
    alpham(i) = 0.1.*((25-Vm(i))/(exp((25-Vm(i))./10)-1));
    betam(i) = 4*exp(-1.*Vm(i)/18);
    alphan(i) = 0.01*((10-Vm(i))/(exp((10-Vm(i))/10)-1));
    betan(i) = 0.125*exp(-1.*Vm(i)/80);
    alphah(i) = 0.07*exp(-1.*Vm(i)/20);
    betah(i) = 1/(exp((30-Vm(i))/10)+1);
    
    gK(i+1) = gKbar*n(i).^4;
    gNa(i+1) = gNabar*m(i).^3.*h(i);
    gL(i+1) = gLbar;
    
    INa(i) = gNa(i).*(Vm(i)-ENa);
    IK(i) = gK(i).*(Vm(i)-EK);
    IL(i) = gL(i).*(Vm(i)-EL);
    Iion(i) = I(i)-IK(i)-INa(i)-IL(i);
    
    dmdt(i) = alpham(i).*(1-m(i))-betam(i).*m(i);
    dndt(i) = alphan(i).*(1-n(i))-betan(i).*n(i);
    dhdt(i) = alphah(i).*(1-h(i))-betah(i).*h(i);
    
    m(i+1) = m(i) + hstep.*dmdt(i);
    n(i+1) = n(i) + hstep.*dndt(i);
    h(i+1) = h(i) + hstep.*dndt(i);
    
    Vm(i+1) = Vm(i) + hstep.*Iion(i)/Cm;
    
end

Vm = Vm + Vrest; %[mV] scaling back data by 70 mV offset

%Plotting figure for Vm in question 1
figure;
subplot(1,2,1);
hold on;
title('Membrane Potential Vm');
xlabel('Time [s]');
ylabel('Voltage [mV]');
plot(t,Vm,'g*-');
legend('V_m');
axis([0 0.1 -75 25]);
hold off;

%Plotting figure for gK and gNa in question 1
subplot(1,2,2);
hold on;
title('Potassium and Sodium Conductances');
xlabel('Time [ms]');
ylabel('Area-normalized Conductance [mS/cm^2]');
plot(t,gK,'bo-');
plot(t,gNa,'ro-');
legend('g_K','g_Na');

