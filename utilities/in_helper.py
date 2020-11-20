import channel as ch

# get mesh info
dnsin = ch.read_dnsin('')
mesh = ch.mesh(dnsin)

retau = input('Please enter extimate of friction Re number:   ')
retau = float(retau)

print(dnsin)

print('dx_plus', mesh.dx  * retau)
print('dz_plus', mesh.dz  * retau)
print('dy_wall', mesh.dyw * retau)
print('dy_center', mesh.dyc * retau)

print('simulation time to achieve 150 mixed units', 150 / retau * dnsin['re'])