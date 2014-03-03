import numpy as np
from astropy.io import fits as pyfits
from tvtk.api import tvtk
import sys

directory = sys.argv[1]
B_r_filename     = directory+'/output.sphere__B0r____.fits'
B_theta_filename = directory+'/output.sphere__B0t____.fits'
B_phi_filename   = directory+'/output.sphere__B0p____.fits'


#### Load theta-grid from file ##############################################
print 'Loading theta-grid ... '
fitsfile = pyfits.open('output.theta_grid.fits')
image = fitsfile[0].data
fitsfile.close()
theta_values = np.array(list(image[:,0]))
sinthetas = np.sin(theta_values)
costhetas = np.cos(theta_values)
N_theta = theta_values.size
del image
del fitsfile


#### Set phi-grid ###########################################################
print 'Setting phi-grid ... '
N_phi = N_theta*2
phi_values = np.arange(0,N_phi)*np.pi*2.0/N_phi
sinphis = np.sin(phi_values)
cosphis = np.cos(phi_values)


#### Set r-grid #############################################################
print 'Loading r-grid ... '
f = open('output.model_zgrid_ext','r')
lines = f.readlines()
f.close()
N_r = len(lines)
r_values = np.zeros([N_r])
i = 0
for line in lines:
   r_values[i] = float(line.split()[1])
   i += 1
r_values = r_values / np.max(r_values) 
del lines
del i
del line


#### Load field data ########################################################
print 'Loading field data ... '
fitsfile = pyfits.open(B_r_filename)
B_r = fitsfile[0].data
fitsfile.close()
print '   B_r done (size: '+str(B_r.size)+')'
fitsfile = pyfits.open(B_theta_filename)
B_theta = fitsfile[0].data
fitsfile.close()
print '   B_theta (size: '+str(B_theta.size)+')'
fitsfile = pyfits.open(B_phi_filename)
B_phi = fitsfile[0].data
fitsfile.close()
print '   B_phi done (size: '+str(B_phi.size)+')'
del fitsfile


#### Convert to Cartesian coordinates #######################################

thinning_factor = 8.0
N_r = int(np.floor(N_r / thinning_factor))
N_theta = int(np.floor(N_theta / thinning_factor))
N_phi = int(np.floor(N_phi / thinning_factor))


#### Convert to Cartesian coordinates #######################################
print 'Converting to Cartesian coordinates ... '

rescale = 1.0
factor_rt = rescale**(1+2)
factor_p = rescale**(1+1)
print '>max Br, Bt, Bp>',factor_rt*np.max(np.abs(B_r)),factor_rt*np.max(np.abs(B_theta)),factor_p*np.max(np.abs(B_phi))

dims = (N_r,N_theta,N_phi+1)
pts = np.empty([dims[0]*dims[1]*dims[2],3])
vector = np.empty([dims[0]*dims[1]*dims[2],3])

rsun=6.96e10
n = 1
a = 0.425

i = 0
for phi_i in sum([range(0,N_phi),[0]],[]):
   phi_index = int(thinning_factor*phi_i)
   phival = phi_values[int(phi_index)]
   for theta_i in range(0,N_theta):
      theta_index = int(thinning_factor*theta_i)
      thetaval = theta_values[theta_index]
      for r_i in range(0,N_r):
         r_index = int(thinning_factor*r_i)
         rval = r_values[r_index]

         # Magentic field from file
         Br = B_r[r_index,theta_index,phi_index] * factor_rt
         Bt = B_theta[r_index,theta_index,phi_index] * factor_p
         Bp = B_phi[r_index,theta_index,phi_index] * factor_rt

#         # Low and Lou (1990) solution prescibed here
#         Br = (rval/rsun)**(-(2+n))*dPdtheta[theta_index]/np.sin(thetaval)
#         Bt = (rval/rsun)**(-(2+n))*n*P[theta_index]/np.sin(thetaval)
#         Bp = (rval/rsun)**(-(1+n))*a*P[theta_index]**(1+1.0/n)

         # Convert to Cartesian coordinates
         x = rval * np.sin(thetaval) * np.cos(phival)
         y = rval * np.sin(thetaval) * np.sin(phival)
         z = rval * np.cos(thetaval)
         pts[i,0] = x
         pts[i,1] = y
         pts[i,2] = z
         vector[i,0] = ( np.sin(thetaval) * np.cos(phival) * Br 
                 + np.cos(thetaval) * np.cos(phival) * Bt
                 - np.sin(phival) * Bp )
         vector[i,1] = ( np.sin(thetaval) * np.sin(phival) * Br
                 + np.cos(thetaval) * np.sin(phival) * Bt
                 + np.cos(phival) * Bp )
         vector[i,2] = ( np.cos(thetaval) * Br
                 - np.sin(thetaval) * Bt )

         # Dipole field
#         vector[i,0] = 3.0*x*z / rval**5
#         vector[i,1] = 3.0*y*z / rval**5
#         vector[i,2] = 3.0*z*z / rval**5 - 1.0 / rval**3

         i += 1

sgrid = tvtk.StructuredGrid(dimensions=dims)
sgrid.points = pts
sgrid.point_data.vectors = vector
sgrid.point_data.vectors.name = 'magnetic field vector'

#### Visualize the field ####################################################
from mayavi import mlab
fig = mlab.figure(1, size=(600, 600), bgcolor=(1, 1, 1), fgcolor=(0, 0, 0))

#surf = mlab.pipeline.surface(sgrid, opacity=0.1)
#mlab.pipeline.surface(mlab.pipeline.extract_edges(surf),
#                            color=(0, 0, 0), )


magnitude = mlab.pipeline.extract_vector_norm(sgrid)
#contours = mlab.pipeline.iso_surface(magnitude,
#                                        contours=[0.0, 0.5, 1.0 ],
#                                        transparent=True,
#                                        opacity=0.9,
#                                        colormap='YlGnBu',
#                                        vmin=0, vmax=2)


field_lines_sphere = mlab.pipeline.streamline(magnitude, seedtype='sphere',
                                        integration_direction='both',
                                        colormap='blue-red')
field_lines_sphere.stream_tracer.maximum_propagation = 100.
field_lines_sphere.seed.widget.theta_resolution = 20
field_lines_sphere.seed.widget.phi_resolution = 10
field_lines_sphere.seed.widget.radius = 0.8
field_lines_sphere.seed.widget.enabled = True

#field_lines_plane = mlab.pipeline.streamline(magnitude, seedtype='plane',
#                                        integration_direction='both',
#                                        colormap='blue-red')
#field_lines_plane.seed.widget.point1 = [  1.0, -1.0, 0.1 ]
#field_lines_plane.seed.widget.point2 = [ -1.0,  1.0, 0.1 ]
#field_lines_plane.seed.widget.normal = [  0.0,  0.0, 1.0 ]
#field_lines_plane.seed.widget.resolution = 25
#field_lines_plane.seed.widget.enabled = False

outline = mlab.pipeline.outline(magnitude)
outline.outline_mode = 'cornered'

mlab.view() #42, 73, 104, [79,  75,  76])

mlab.show()

