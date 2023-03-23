from __future__ import print_function
import sys
import numpy as np

lc= 0.001

injector_corner_refine_ratio= 0.3 #0.1
fuel_injector_refine_ratio= 0.15

diameter_chamber= 4.4e-02

# important shape parameters of the injectors
diameter_oxidizer_injector= 2.047e-02 # 0.05,0.03
diameter_fuel_injector= 0.2e-02 #0.066e-02 #0.02*diameter_chamber
# distance_injector= 0.1*diameter_chamber


distance_injector= 0.13e-02 #0.03*diameter_chamber # 0.05
distance_wall= 0.5*diameter_chamber - diameter_fuel_injector -distance_injector - 0.5*diameter_oxidizer_injector

# distance_wall= (diameter_chamber - num_injectors*diameter_injector - distance_injector*(num_injectors-1))/2

diameter_nozzle= 2.08e-02*0.75 #0.5*diameter_chamber  # 0.75,0.5
diameter_exit= 3.9e-02*0.6 #2.5e-02 #0.8*diameter_chamber   # 1.25,0.8

if(distance_injector<0):
	print('Error: negative distance_injector, the value is:',distance_injector)
	sys.exit()


length_oxidizer_injector= 0.04 #0.5*diameter_chamber  #5
length_fuel_injector= 2.79e-02 #0.3*diameter_chamber  #5
length_step= 1.01e-02 #0.02*diameter_chamber

length_chamber= 0.2 #0.1 #0.5*diameter_chamber  #1.0，0.4

length_converge= 1.0e-2*0.5 #0.1*diameter_chamber #0.1
length_diverge= 1.0e-2 #0.1*diameter_chamber  #0.1
diverge_ratio= 0.3

num_cpts_converge= 20 #10
num_cpts_diverge=  20 #10

total_cpts= 14 + num_cpts_converge + num_cpts_diverge


f=open('rocket.geo','w+')

# f.write('Mesh.Algorithm = 8; // Delaunay for quads'+'\n')

f.write('lc='+str(lc)+';'+'\n')
f.write('injector_corner_refine_ratio='+str(injector_corner_refine_ratio)+';'+'\n')
f.write('fuel_injector_refine_ratio='+str(fuel_injector_refine_ratio)+';'+'\n')

# Points
ind=0
z=0

# point 1
ind+=1
x= 0
y= 0.5*diameter_chamber
f.write('Point('+str(ind)+')={'+str(x)+','+str(y)+','+str(z)+',injector_corner_refine_ratio*lc};'+'\n')

# point 2, go down
ind+=1
# x -= length_injector 
#0.5*diameter_chamber-distance_wall - i*(diameter_injector+distance_injector)
y -= distance_wall
f.write('Point('+str(ind)+')={'+str(x)+','+str(y)+','+str(z)+',injector_corner_refine_ratio*lc};'+'\n') ## 4

# point 3, go left
ind+=1
x -= length_step
f.write('Point('+str(ind)+')={'+str(x)+','+str(y)+','+str(z)+',fuel_injector_refine_ratio*lc};'+'\n') ## 4

# point 4, go left
ind+=1
x = -length_fuel_injector
f.write('Point('+str(ind)+')={'+str(x)+','+str(y)+','+str(z)+',fuel_injector_refine_ratio*lc};'+'\n') ## 4

# point 5, go down
ind+=1
# x += length_injector
y -= diameter_fuel_injector 
f.write('Point('+str(ind)+')={'+str(x)+','+str(y)+','+str(z)+',fuel_injector_refine_ratio*lc};'+'\n')

# point 6, go right
ind+=1
x = -length_step
f.write('Point('+str(ind)+')={'+str(x)+','+str(y)+','+str(z)+',fuel_injector_refine_ratio*lc};'+'\n')

# point 7, go down
ind+=1
# x += length_injector
y -= distance_injector
f.write('Point('+str(ind)+')={'+str(x)+','+str(y)+','+str(z)+',injector_corner_refine_ratio*lc};'+'\n')

# point 8, go left
ind+=1
# x += length_injector
x = -length_oxidizer_injector
f.write('Point('+str(ind)+')={'+str(x)+','+str(y)+','+str(z)+',injector_corner_refine_ratio*lc};'+'\n')

# point 9, go down
ind+=1
# x += length_injector
# y -= diameter_oxidizer_injector 
y = 0
f.write('Point('+str(ind)+')={'+str(x)+','+str(y)+','+str(z)+',injector_corner_refine_ratio*lc};'+'\n')

# point 10, go right
ind+=1
x = -length_step
# y -= diameter_oxidizer_injector 
y = 0
f.write('Point('+str(ind)+')={'+str(x)+','+str(y)+','+str(z)+',injector_corner_refine_ratio*lc};'+'\n')

# point 11, go right
ind+=1
# x += length_injector
x = 0.5*length_chamber 
f.write('Point('+str(ind)+')={'+str(x)+','+str(y)+','+str(z)+',injector_corner_refine_ratio*lc};'+'\n')

# point 12, go right
ind+=1
# x += length_injector
x = length_chamber + length_converge + length_diverge
f.write('Point('+str(ind)+')={'+str(x)+','+str(y)+','+str(z)+',injector_corner_refine_ratio*lc};'+'\n')

# point 13, go up
ind+=1
# x += length_injector
# y += 0.5*diameter_exit
dsx= diverge_ratio*length_diverge
dsx= length_diverge - dsx
y= 0.25*(diameter_nozzle + diameter_exit) + 0.25*(diameter_exit - diameter_nozzle)*np.cos( dsx*np.pi/length_diverge)
f.write('Point('+str(ind)+')={'+str(x)+','+str(y)+','+str(z)+',injector_corner_refine_ratio*lc};'+'\n')
# f.write('Point('+str(ind)+')={'+str(x)+','+str(y)+','+str(z)+',lc};'+'\n')



# divergence points
dx= length_diverge/num_cpts_diverge
for i in range(1,num_cpts_diverge+1):
	ind +=1
	# x= length_chamber + length_converge + length_diverge - i*dx
	# y= 0.25*(diameter_nozzle + diameter_exit) + 0.25*(diameter_exit - diameter_nozzle)*np.cos( (i*dx)*np.pi/length_diverge)
	x= length_chamber + length_converge + length_diverge - i*dx
	dsx= diverge_ratio*(length_diverge - i*dx)
	dsx= length_diverge - dsx
	y= 0.25*(diameter_nozzle + diameter_exit) + 0.25*(diameter_exit - diameter_nozzle)*np.cos( dsx*np.pi/length_diverge)
	# f.write('Point('+str(ind)+')={'+str(x)+','+str(y)+','+str(z)+',lc};'+'\n')
	f.write('Point('+str(ind)+')={'+str(x)+','+str(y)+','+str(z)+',injector_corner_refine_ratio*lc};'+'\n')


dx= length_converge/num_cpts_converge
for i in range(1,num_cpts_converge+1):
	ind +=1
	x= length_chamber + length_converge - i*dx
	y= 0.25*(diameter_chamber + diameter_nozzle) + 0.25*(diameter_nozzle - diameter_chamber)* np.cos( (i*dx)*np.pi/length_converge )
	# f.write('Point('+str(ind)+')={'+str(x)+','+str(y)+','+str(z)+',lc};'+'\n')
	f.write('Point('+str(ind)+')={'+str(x)+','+str(y)+','+str(z)+',injector_corner_refine_ratio*lc};'+'\n')

# Middle-point on top of the chamber
ind+=1
# x += length_injector
x = 0.5*length_chamber
y = 0.5*diameter_chamber
f.write('Point('+str(ind)+')={'+str(x)+','+str(y)+','+str(z)+',injector_corner_refine_ratio*lc};'+'\n')

####################
# Lines
####################
for i in range(1,total_cpts+1):
	f.write('Line('+str(i)+')={'+str(i)+','+str(np.mod(i,total_cpts)+1)+'};'+'\n')


f.write('Line('+str(total_cpts+1)+')={'+str(3)+','+str(6)+'};'+'\n')
f.write('Line('+str(total_cpts+2)+')={'+str(7)+','+str(10)+'};'+'\n')
f.write('Line('+str(total_cpts+3)+')={'+str(11)+','+str(total_cpts)+'};'+'\n')


####################
# Curve loops
####################
# Curve loop for fuel injector
f.write('Curve Loop(1)={3,4,5,'+str(-(total_cpts+1))+'};'+'\n')

# Curve loop for oxidizer injector
f.write('Curve Loop(2)={7,8,9,'+str(-(total_cpts+2))+'};'+'\n')

# Curve loop for left part of chamber
f.write('Curve Loop(3)={1,2,'+str(total_cpts+1)+',6,'+str(total_cpts+2)+',10,'+str(total_cpts+3)+','+str(total_cpts)+'};'+'\n')

# Curve loop for right part of chamber
f.write('Curve Loop(4)={11:'+str(total_cpts-1)+','+str(-(total_cpts+3))+'};'+'\n')





####################
# Plane surfaces
####################
f.write('Plane Surface(1)={1};'+'\n')
f.write('Plane Surface(2)={2};'+'\n')
f.write('Plane Surface(3)={3};'+'\n')
f.write('Plane Surface(4)={4};'+'\n')

###########################
# Structured mesh surface
###########################
# f.write('Transfinite Surface{2}={3,4,5,'+str(-(total_cpts+1))+'};'+'\n')
# f.write('Transfinite Line {3,5} = 100 Using Progression 1;'+'\n')
# f.write('Transfinite Line {4,'+str(-(total_cpts+1))+'} = 20 Using Progression 1;'+'\n')
# f.write('Recombine Surface {2};'+'\n')

####################
# Physical curves
####################
fuel_inlet= 'Physical Curve(1)={4};'
f.write(fuel_inlet+'\n')

oxidizer_inlet= 'Physical Curve(2)={8};'
f.write(oxidizer_inlet+'\n')

# solid_wall= 'Physical Curve(3)={1:3,5:7,13:'+str(total_cpts)+'};'
# f.write(solid_wall+'\n')
solid_wall= 'Physical Curve(3)={1:3,5:7,'+str(total_cpts-1)+':'+str(total_cpts)+'};'
f.write(solid_wall+'\n')

nozzle_exit= 'Physical Curve(4)={12};'
f.write(nozzle_exit+'\n')

symmetry= 'Physical Curve(5)={9:11};'
f.write(symmetry+'\n')

solid_wall= 'Physical Curve(6)={13:'+str(total_cpts-2)+'};'
f.write(solid_wall+'\n')

# # wall boundary
# wall_curve= 'Physical Curve('+str(2)+')={'
# for i in range(num_injectors):
# 	wall_curve += str(4*i+1)+','
# 	wall_curve += str(4*i+3)+','
# 	wall_curve += str(4*i+4)+','

# wall_curve += str(4*num_injectors+1)+':'+str(4*num_injectors+bottom_wall_cpts-1)+','
# wall_curve += str(4*num_injectors+bottom_wall_cpts+1)+':'+str(total_cpts)

# wall_curve += '};'
# f.write(wall_curve+'\n')

# # outle boundary
# outlet_curve= 'Physical Curve('+str(3)+')={'+str(4*num_injectors+bottom_wall_cpts)+'};'
# f.write(outlet_curve+'\n')

####################
# Physical surface
####################
f.write('Physical Surface("My surface")={1:4};'+'\n')

f.close()
sys.exit()

####### Boundary Layer Filed ######################
#boundary layer parameters

# Boundary layer for fuel injector
hwall_n= 1.0e-3*diameter_fuel_injector  #2.0e-4,1.0e-5
ratio= 1.1
thickness= 1.0e-1*diameter_fuel_injector

EdgesList='3,5,7'

NodesList = str(3)+':'+str(8)

f.write('Field[10]= BoundaryLayer;'+'\n')
f.write('Field[10].EdgesList= {'+EdgesList+'};'+'\n')
f.write('Field[10].NodesList= {'+NodesList+'};'+'\n')
f.write('Field[10].IntersectMetrics=1;'+'\n')
f.write('Field[10].Quads= 1;'+'\n')
f.write('Field[10].hwall_n= '+str(hwall_n)+';'+'\n')
f.write('Field[10].ratio='+str(ratio)+';'+'\n')
f.write('Field[10].thickness='+str(thickness)+ ';'+'\n')
f.write('BoundaryLayer Field = {10};'+'\n')


# # Boundary layer for oxidizer injector
# hwall_n= 1.0e-4*diameter_oxidizer_injector
# ratio= 1.1
# thickness= 1.0e-1*diameter_oxidizer_injector

# EdgesList='7'

# NodesList = '7,8'

# f.write('Field[11]= BoundaryLayer;'+'\n')
# f.write('Field[11].EdgesList= {'+EdgesList+'};'+'\n')
# f.write('Field[11].NodesList= {'+NodesList+'};'+'\n')
# f.write('Field[11].IntersectMetrics=1;'+'\n')
# f.write('Field[11].Quads= 1;'+'\n')
# f.write('Field[11].hwall_n= '+str(hwall_n)+';'+'\n')
# f.write('Field[11].ratio='+str(ratio)+';'+'\n')
# f.write('Field[10].thickness='+str(thickness)+ ';'+'\n')

# f.write('BoundaryLayer Field = {11};'+'\n')
f.close()



