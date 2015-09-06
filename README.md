# sphere-lab

## Matlab code for spherical harmonic computations

Some routines for spherical harmonic transform related work

### Spherical Harmonics

- Spherical harmonics calculations use the Schmidt semi-normalized Legendre functions, which means degrees of the order of l=2000 can be handled.  Using the unnormalized (regular) associated Legendre functions in the obvious way the calculations fail around l=150.
- Evaluations are done on a theta-phi mesh.  When the mesh separates in theta then there are huge reductions in complexity (increases in speed).  The non-separable mesh needs to be done point-by-point
- The spectral representation is a vector indexed by n and not two indices l and m.  Converting between the two respresentations is easy.
- The inverse spherical harmonic transform takes the spectral vector and creates the spatial function on a mesh, returning both the function and the mesh.  The spatial function is a weighted combination of spherical harmonics. If the spectral representation is a delta (a single non-zero weight) then the spherical harmonix is returned.

### Integration on the Sphere

- Inner products and the Spherical Harmonic Transform need a double integral over the sphere in the natural measure.  Numerically this is done using a trapezoid rule and is achieved in two lines of code.
- The region is defined through a mesh.  To do: add a mask to mesh to have a irregular region.

### Indices

- l: l>=0 is the degree
	- when bandlimited the maximum non-zero index is L_max, and L_tot=L_max+1
- m: abs(m)<=l is the order
- n: n>=0 is the single countable index
	- n=l*(l+1)+m
	- l=floor(sqrt(r-1)); m=r-1-l*(l+1);
	- when bandlimited the maximum non-zero index is N_max=(L_max+1)^2-1, and N_tot=(L_max+1)^2
