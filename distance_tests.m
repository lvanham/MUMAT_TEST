clc
p = '\\wsl.localhost\Debian\home\lvh\dist_test_solid_2\';
%p = '\\wsl.localhost\Debian\home\lvh\dist_test_solid\';
%cobj = flip({[p 'B_50.dat']; [p 'B_100.dat']; [p 'B_200.dat']; [p 'B_400.dat'] ; [p 'B_600.dat']; [p 'B_800.dat']; [p 'B_1000.dat']});
%cobj = flip({[p 'B_400.dat'] ; [p 'B_600.dat']; [p 'B_800.dat']; [p 'B_1000.dat']});

%cobj = flip({[p 'B_H125.dat']; [p 'B_H250.dat'] ; [p 'B_H500.dat']; [p 'B_H1000.dat']; [p 'B_H1500.dat']});% [p 'B_H2000.dat']});
cobj =  {[p 'B_S1000_1.dat']; [p 'B_S1000_2.dat']; [p 'B_S500_2.dat']; [p 'B_S500_1.dat'];[p 'B_S250_2.dat']; [p 'B_S250_1.dat']; [p 'B_S50_1.dat']};
MUMATplot(2, 'comp_z', [p 'points.dat'], [p 'B_exact.dat'], cobj, 0, [0 0 1])
