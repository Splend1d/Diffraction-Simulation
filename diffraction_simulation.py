#this code finds the possible diffraction points of a crystal
#with 1 -1 1 pointing towards the beam.
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
'''
calculates the distance between two points:
if only one point is given, the second point is assumed to be (0,0,0)
2D distance is supported
'''
def dist(*args):
	if len(args) == 2:
		p1 = args[0]
		p2 = args[1]
	else:
		p1 = args[0]
		p2 = [0,0,0]

	d_squared = [0,0,0]
	d_squared[0] = (p1[0]-p2[0]) ** 2
	d_squared[1] = (p1[1]-p2[1]) ** 2
	try:
		d_squared[2] = (p1[2]-p2[2]) ** 2
	except:
		d_squared[2] = 0

	return np.sqrt(sum(d_squared))

# calculates the cross product of two 3D vectors
def cross(a, b):
    c = [a[1]*b[2] - a[2]*b[1],
         a[2]*b[0] - a[0]*b[2],
         a[0]*b[1] - a[1]*b[0]]
    return c

# decomposes the target vector into linear combination of
# the sources, specified by *args
def decompose(target,*args):
	component = []
	for x in range(len(args[0])):
		for y in range(len(args)):
			try:
				component[x].append(args[y][x])
			except:
				component.append([args[y][x]])
	a = np.array(component)
	b = np.array(target)

	return np.linalg.solve(a,b)
#intensity of a plane
def intensity(plane):
	l_p = 1/sin(2*theta)



class Crystal():
	def __init__(self, lattice_constant, structure, orientation_towards, orientation_up):
		self.lc = lattice_constant
		self.structure = structure
		self.o_towards = orientation_towards
		self.o_upwards = orientation_up
		self.o_right = cross(orientation_towards,orientation_up)
		self.rotation = [0,0,0]

	def __repr__(self):
		return '{}.{}(lattice constant={!r}, crystal sturcture={!r})'.format(self.__class__.__module__,
				self.__class__.__name__,self.lc,self.structure)

	def __str__(self):
		return u"crystal is %s structure with lattice constant %f " % (self.structure, self.lc)

	# determine whether a plane satisfies structure
	def structure_factor(self, plane):
		if self.structure == 'diamond':
			if plane[0] % 2 == plane[1] % 2 and plane[1] % 2 == plane[2] % 2:
				if sum(plane) % 4 != 2:
					return True
			return False
		elif self.structure == 'FCC':
			if plane[0] % 2 == plane[1] % 2 and plane[1] % 2 == plane[2] % 2:
				return True
			return False
		elif self.structure == 'BCC':
			if (plane[0] + plane[1] + plane[2])%2 == 0:
				return True
			return False

class Device():
	def __init__(self, lambda_min, lambda_MAX, camera_length, film_size):
		self.ld_min = lambda_min
		self.ld_MAX = lambda_MAX
		self.length = camera_length
		self.size = film_size
		self.crystal_diffraction_points = []
		self.film_diffraction_points = []

	def __repr__(self):
		return '{}.{}(scanning wavelength min = {!r}, Max = {!r}, film size = {!r}, distance = {!r})'.format(self.__class__.__module__,
				self.__class__.__name__, self.ld_min, self.ld_MAX, self.length, self.size)

	def __str__(self):
		return u"scanning wavelength from %f to %f, film size is %f at %f away\n%s" % (self.ld_min, self.ld_MAX, self.size, self.length, self.crystal)

	def diffract(self,crystal,mode):

		#clear previous data
		self.crystal_diffraction_points = []
		self.film_diffraction_points = []

		#determine whether plane lies between two ewald spheres
		#part1: determine seatch region, method: smallest cube that surrounds the sphere
		_axis_towards_film_dist =  dist(crystal.o_towards) * 1/crystal.lc
		_axis_upwards_dist = dist(crystal.o_upwards) * 1/crystal.lc
		_axis_sidewards_dist = dist(crystal.o_right) * 1/crystal.lc

		r_MAX = 1/self.ld_min
		origin_of_serach = [ r_MAX/_axis_towards_film_dist * x for x in crystal.o_towards ]
		origin_of_serach_int = [ int(r_MAX/_axis_towards_film_dist) * x for x in crystal.o_towards ]
		radius_of_serach = int(dist(origin_of_serach_int)) + 1

		#part2: if diffraction point is smaller than larger sphere, it is a possible diffraction point
		diffraction_points = []
		for h in range(origin_of_serach_int[0] - radius_of_serach, origin_of_serach_int[0] + radius_of_serach + 1):
			for k in range(origin_of_serach_int[1] - radius_of_serach, origin_of_serach_int[1] + radius_of_serach + 1):
				for l in range(origin_of_serach_int[2] - radius_of_serach, origin_of_serach_int[2] + radius_of_serach + 1):
					if dist([h,k,l],origin_of_serach)/crystal.lc < r_MAX:
						diffraction_points.append([h,k,l])


		r_min = 1/self.ld_MAX
		origin_of_serach = [ r_min/_axis_towards_film_dist * x for x in crystal.o_towards ]
		origin_of_serach_int = [int(r_min/_axis_towards_film_dist)*x for x in crystal.o_towards]
		radius_of_serach = int(dist(origin_of_serach_int)) + 1

		#part3: if diffraction point is smaller than the smaller sphere, it is *not* a diffraction point
		for h in range(origin_of_serach_int[0] - radius_of_serach, origin_of_serach_int[0] + radius_of_serach + 1):
			for k in range(origin_of_serach_int[1] - radius_of_serach, origin_of_serach_int[1] + radius_of_serach + 1):
				for l in range(origin_of_serach_int[2] - radius_of_serach, origin_of_serach_int[2] + radius_of_serach + 1):
					if dist([h,k,l],origin_of_serach_int)/crystal.lc < r_min:
						try:
							diffraction_points.remove([h,k,l])
						except:
							continue

		# planes after this step are all diffraction planes
		# we then filter the ones that truly landed on the film

		# find diffraction points that are available in this structure factor
		for plane in diffraction_points:
			if not crystal.structure_factor(plane):
				continue
			else:
				# decompose plane into unit directions aligning with the crystal orientation
				s = decompose(plane, crystal.o_towards, crystal.o_upwards, crystal.o_right)

				#crystal position in reciprocal relative to O
				kx = s[0]*_axis_towards_film_dist
				ky = s[1]*_axis_upwards_dist
				kz = s[2]*_axis_sidewards_dist

				# small angle rotation
				phi = crystal.rotation[0] * 2 * np.pi / 360
				theta = crystal.rotation[1] * 2 * np.pi / 360
				kz2 = kz * np.cos(phi) + -kx * np.sin(phi)
				kx2 = kz * np.sin(phi) + kx * np.cos(phi)
				kx = kx2
				kz = kz2
				kx2 = kx * np.cos(theta) + -ky * np.sin(theta)
				ky2 = kx * np.sin(theta) + ky * np.cos(theta)
				kx = kx2
				ky = ky2

				# diffraction beam direction
				beam_side = np.sqrt(kz ** 2 + ky ** 2)
				beam_x = (kx ** 2 - beam_side * beam_side)/(2 * kx)
				beam_y = ky
				beam_z = kz

				# "back" scattering condition
				if mode == 'back' and beam_x < 0:
					continue
				elif mode == 'front' and beam_x > 0:
					continue

				# real space diffraction point position
				rx = self.length
				ry = beam_y * self.length / beam_x
				rz = beam_z * self.length / beam_x
				if np.abs(ry) < 5 and np.abs(rz) < 5:
					self.crystal_diffraction_points.append(plane)
					self.film_diffraction_points.append([rz,ry])
		return

	def showfilm(self):
		#specify plot properties
		plt.figure()
		#film size eliminate diffraction points that cannot be detected by film
		plt.xlim(-self.size/2,self.size/2)
		plt.ylim(-self.size/2,self.size/2)
		plt.scatter(0,0,color='gray',marker = 'x')
		for point in self.film_diffraction_points:
			plt.scatter(point[0],point[1],color = 'black')
		plt.show()
		return

	def showreciprocal(self,lc):
		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')
		for plane in self.crystal_diffraction_points:
			ax.scatter(plane[0],plane[1],plane[2],color = 'black',s = 5)
		u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
		r = 1 / self.ld_min * lc
		x = r/np.sqrt(3) + r*np.cos(u)*np.sin(v)
		y = -r/np.sqrt(3) + r*np.sin(u)*np.sin(v)
		z = r/np.sqrt(3)+r*np.cos(v)
		ax.plot_wireframe(x, y, z, color="r")
		r = 1 / self.ld_MAX * lc
		x = r/np.sqrt(3) + r*np.cos(u)*np.sin(v)
		y = -r/np.sqrt(3) + r*np.sin(u)*np.sin(v)
		z = r/np.sqrt(3)+r*np.cos(v)
		ax.plot_wireframe(x, y, z, color="blue")
		plt.show()
		return

# sample input
c = Crystal(5.431,'diamond',[1,-1,1],[-1,1,2])
d = Device(0.6,1.2,9,10)
d.diffract(c,'back')
d.showreciprocal(5.431)
print(len(d.crystal_diffraction_points))
d.showfilm()
c.rotation = [5, 3]
d.diffract(c,'back')
d.showfilm()

