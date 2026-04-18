
Channel is a program that performs direct numerical simulations of plane channel flows.

## Geometry, length and velocity scales

We start by defining the geometry of our channel flow. This consists in two indefinite, parallel plates; the y-coordinate is defined to be orthogonal to both plates. We then specify the y-coordinates of each of the plates; this is done by setting the ymin, ymax values in dns.in, so that - for instance:
```
xxx 0.0 2.0       ! a, ymin, ymax
```
Here, $xxx$ represents the value of the stretching parameter a, which we will discuss later. In this case, ymin=0 and ymax=2, meaning that the bottom wall will be found at y=0 and the top wall at y=2. Notice that all values are non-dimensional in the simulation, including the value of ymin and max that we are setting. What is then the length scale that we are using?

Consider the plane with y=1: this is located midway between the two plates. Because its non-dimensional coordinate is 1, its dimensional coordinate will be equal to the length scale $L^*$ used for non-dimensionalisation. In this case, $L^*$ corresponds to the channel half-height $h^*$ (because y=1 is the midplane); if one were to use ymin=0 and ymax=1, then the length scale would be the full channel height ($L^* = 2h^*$).

Anyway, it is advised to always use ymin=0 and ymax=2, so that any length will be always implicitly made non-dimensional with the channel half height.

As for the velocity scale, its value is defined by the type of forcing used (constant flow rate or constant pressure gradient), or by the velocity at the wall.

### Constant pressure gradient (CPG)

With __constant pressure gradient (CPG)__ we mean that the flow is forced by a mean pressure gradient that is constant in time. The values of the pressure gradient are set by meanpx, meanpz (see dns.in). If the magnitude of the pressure gradient is 1 ($\sqrt{meanpx^2 + meanpz^2} = 1$), then the velocity scale is implicitly set to the friction velocity and all velocity values produced by the simulation will be automatically in viscous units (please notice this statement is not true be wrong if u0 and uN have different values).

A typical Poiseuille flow forced with CPG would be:
```
1.0 0.0 !meanpx, meanpz
0.0 0.0 !meanflowx, meanflowz
0.0 0.0 !u0,uN
```
meaning that there will be a pressure gradient in the x direction only, no CFR forcing is used and the walls are still. The velocity scale is the friciton velocity. Notice that it can be beneficial to set the values of u0, uN to the same negative value (eg `-10 -10 !u0, uN`); a correct choice of this value lowers the maximum absolute value of the velocity that one sees in the flow, which in turn means that larger values of $\Delta t$ fulfill the CFL condition. In other words, this allows to have longer timesteps, which is beneficial for simulations that run for long. Notice that the value of u0, uN needs to be the same for the flow to be a Poiseuille one; the velocity scale remains the friction velocity as long as $u0=uN$.

Do not use CPG in combination with CFR.

### Constant flow rate (CFR)

With __constant flow rate (CFR)__ we mean that the flow is forced so that its flow rate is constant in time. The constant value of the flow rate is set by meanflowx, meanflowz. If the magnitude of the mean flow rate is 1 ($\sqrt{meanflowx^2 + meanflowz^2} = 1$), the velocity scale is implicitly set to be the bulk velocity (XXX or two times the bulk velocity? also, verify this statement if u0, uN are not zero). A typical Poiseuille flow forced with CFR would be:
```
0.0 0.0 !meanpx, meanpz
1.0 0.0 !meanflowx, meanflowz
0.0 0.0 !u0,uN
```
meaning that the flow rate will be forced only in the x-direction, no CPG forcing is used and the walls are still. The velocity scale is the bulk (XXX) velocity.

Do not use CFR in combination with CPG.

### Setting a velocity at the wall

For a Poiseille flow, there is no relative motion between the plates that bound the channel flow; for a Couette flow, instead, the walls move with respect to one another. We then need to set the velocity of the walls; this is done with the u0, uN parameters in dns.in.

A typical Couette flow would be:
```
0.0 0.0 !meanpx, meanpz
0.0 0.0 !meanflowx, meanflowz
-1.0 1.0 !u0,uN
```
meaning that no CPG, no CFR is used and that the walls move in the x direction with the same velocity but opposite sign. The wall velocity is then the velocity scale.

### Reynolds number

Once the velocity and lenght scales are specified, one can enter the Reynolds number (defined indeed with the length and velocity scales implied by previous choices). This is done by setting the value of `ni` field in dns.in.

## Numerical details: box size, resolution

The last things to set are the box size(s) and the resolution. The plates would be infinite in an ideal world, but this is of course not feasible numerically; so, we make a simulation in a periodic, rectangular box of size $L_x$, $L_z$. These two quantities are set in dns.in by alfa0 and beta0, which represent the Fourier resolution in x and z respectively (keep in mind that the code is spectral, meaning that the x and z directions are Fourier transformed). In other words:
$$alfa0 = 2*\pi / L_x$$
$$beta0 = 2*\pi / L_z$$

### Setting the resolution in the x, z directions

Setting the resolution in the x, z directions is rather easy. One increases the number of points nx, nz in dns.in until the resolutions $\Delta x$, $\Delta z$ satisfy:
$$\Delta x^+ \leq 10$$
$$\Delta z^+ \leq 5$$
notice that these constraints are given in plus (viscous) units. Keep in mind that, if $l_g$ is a generic length scaled with $h^*$ (where the asterisk indicates dimensional quantities, so that $l_g = l_g^* / h^*$), the following relationship can be used to switch to viscous units:
$$ l_g^+ = Re_\tau \,l_g$$

Also notice that, owing to the Fourier mesh used in the simulaiton (see the readme.md on github), the resolution is defined like so:
$$\Delta x = L_x / (2 \text{nx})$$
$$\Delta z = L_x / (2 \text{nz})$$
where nx, nz are defined in dns.in. Keep in mind that the spacing is uniform in x, z.

#### Important remarks
- It is wise to set $n_x$ so that $n_x+1$ has only 2 and 3 as prime factors. This will make the choice of parallelisation strategy easier.
- It is wise to set a value of $n_z$ that is a power of 2. This makes the value of $nzd$ smaller and speeds up calculations.

### Setting the resolution in y

This is a bit more complicated. Here, one needs to play with two parameters (ny, a) to fulfill two different constraints:
$$\Delta y_w^+ \leq 1$$
$$\Delta y_c^+ \leq 10$$

Practically, you help yourself with the script in_helper.py.

### Using in_helper.py

In folder `channel/utilities` you can find a python script `in_helper.py` that helps set up the simulations. To use it, simply navigate to the folder containing the `dns.in` file and run:
```
python /path/to/channel/utilities/in_helper.py
```
The script will require an estimate of the friction Reynolds number $Re_\tau$; if unknown, this can estimated using empirical correlations. The script will output the resolution in wall units, some parameters that are relevant for parellelisation and timestepping, like the value of `nzd` and the time interval that is equivalent to 1 mixed unit.

## Timestepping

- Use a spacing between fields equivalent to 1 $h/u_\tau$ (1 mixed unit). The value of this spacing is printed out by in_helper.py; it is set by the value of dt_field in dns.in.
- Usually, good statistics are obtained when the simulations are averaged for a total of 150 $h/u_\tau$ (150 mixed units). One should also consider that the flow undergoes a transient (order of magnitude: 50 $h/u_\tau$) at the beginning of the simulation; this needs to be discarded before averaging. The total simulation time is set by t_max in dns.in; it is actually wise to set it so that every execution of channel terminates when run on the cluster.
