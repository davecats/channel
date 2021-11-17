import channel as ch
from copy import copy

# get mesh info
dnsin = ch.read_dnsin('dns.in')
mesh = ch.mesh(dnsin)

retau = input('Please enter extimate of friction Re number:   ')
retau = float(retau)

def get_nzd(nzz):
    guess = 3*nzz
    once_more = True
    while once_more:
        prev_guess = copy(guess)
        while (guess%2) == 0:
            guess /= 2
        if guess==1 or guess==3:
            once_more = False
        else:
            guess = copy(prev_guess) + 1

    return prev_guess

print()
print('Parallelisation')
print('---------------')
print('nx+1', mesh.nx+1)
print('nzd (extimate!)', get_nzd(mesh.nz))

print()
print('Resolution')
print('----------')
print('dx_plus', mesh.dx  * retau)
print('dz_plus', mesh.dz  * retau)
print('dy_wall', mesh.dyw * retau)
print('dy_center', mesh.dyc * retau)

print()
print('Time')
print('----')
print('One mixed unit in simulation units',  dnsin['re']/retau)
print('Simulation time to achieve 150 mixed units', 150 / retau * dnsin['re'])
print()