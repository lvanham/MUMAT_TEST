clc
p = '\\wsl.localhost\Debian\home\lvh\dist_test_hollow\';
%cobj = flip({[p 'B_50.dat']; [p 'B_100.dat']; [p 'B_200.dat']; [p 'B_400.dat'] ; [p 'B_600.dat']; [p 'B_800.dat']; [p 'B_1000.dat']});
cobj = flip({[p 'B_H125.dat']; [p 'B_H250.dat'] ; [p 'B_H500.dat']; [p 'B_H1000.dat']; [p 'B_H1500.dat']});% [p 'B_H2000.dat']});

MUMATplot(2, 'comp_x', [p 'points.dat'], [p 'B_exact.dat'], cobj, 0)
