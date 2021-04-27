fluid = 'R141b';
P_min = CoolProp.PropsSI('P','T',298.15,'Q',0,fluid);
P_max = 0.95*CoolProp.Props1SI(fluid,'Pcrit');

Final_point = [3.93 4 4000170.719]; %[K_H K_L P_H]