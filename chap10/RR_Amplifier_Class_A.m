% script RR_Amplifier_Class_A
% Solves the basic equations of a Class A amplifier.
% Renaissance Robotics codebase, Chapter 10, https://github.com/tbewley/RR
% Copyright 2022 by Thomas Bewley, distributed under Modified BSD License.

clear, syms s C0 R1 R2 V_in V_s

% Part i: Solve for V_B as a function of V_in and {V_s,C0,R1,R2}.
% x={I_0,I_1,I_2,V_B}      <-- unknown vector (assumes I_B is negligible)
A  =[ 1   1  -1   0;       % I_0 + I_1 - I_2 \approx 0
      1   0   0  C0*s;     % I_0 + C0*s*V_B = C*s*V_in
      0  R1   0   1;       % R1*I_1 + V_B = V_s
      0   0  R2  -1];      % R2*I_2 - V_B = 0
b  =[ 0; C0*s*V_in; V_s; 0];
x=A\b; V_B=simplify(x(4))

clear, syms s V_B V_d V_s alpha R3 R4 C5 Rspeaker

% Part ii: Solve for V_out as a function of V_B and {alpha,V_d,V_s,R3,R4,C5}
% x={I_C,I_B,I_E,V_C,V_E,I_3,I_5,V_out}    <-- unknown vector
A  =[ 1   1  -1   0   0   0   0    0;      % I_B + I_C - I_E = 0
     -1   0 alpha 0   0   0   0    0;      % alpha*I_E - I_C = 0
      0   0   0   0   1   0   0    0;      % V_E = V_B - V_d
      1   0   0   0   0  -1   1    0;      % I_C + I_5 - I_3 = 0
      0   0   0   1   0  R3   0    0;      % R3*I_3 + V_C = V_s
      0   0  R4   0  -1   0   0    0;      % R4*I_E - V_E = 0
      0   0   0 -C5*s 0   0   1  C5*s;     % I_5 + C5*s*V_out - C5*s*V_C = 0
      0   0   0   0   0 0 Rspeaker -1];    % Rspeaker*I_5 - V_out = 0
b  =[ 0; 0; V_B-V_d; 0; V_s; 0; 0; 0];
x=A\b; V_out=simplify(x(8))

% Representative application
clear, syms R3 R4
Vs=12; Imax=0.05; Rtot=Vs/Imax; Rspeaker=8; 
Icmax=0.05; K=10; omega_ii=10  % Design targets
eqn1 = R3+R4==Rtot;
eqn2 = 0.99*(R3/R4)*(Rspeaker/(R3+Rspeaker)) == 10;
sol=solve(eqn1,eqn2,R3,R4); R3=eval(sol.R3(1)), R4=eval(sol.R4(1))
Rii=R3 + Rspeaker
C5 = (1/omega_ii)/Rii
