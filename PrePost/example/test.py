from create_particles import *


'''
Set up the impact geometry
'''


eta = 1.2
col = ['X','Y','Z','VX','VY','VZ','PART','MASS','RHO']


# set parts geometry
def set_parts_geometry():
  # trimesh target: m-kg-s
  tar_cell_r = 25
  tar_rho = 2700
  tar_center = np.array([0,0,0])
  tar_radius = 1000
  tar_vel = np.array([0,0,0])

  # create target particles
  tar = part()
  box_pos = generate_box(box_size=tar_radius*2*np.array([1,1,1]), box_center=tar_center, r=tar_cell_r, lattice='hcp')
  def tar_geometry(x,y,z): 
    return (x-tar_center[0])**2+(y-tar_center[1])**2+(z-tar_center[2])**2 <= tar_radius**2
  tar.create_particles(box_pos,tar_geometry)
  tar.points.PART = 0
  tar.points.MASS = tar_rho*8*tar_cell_r**3
  tar.points.RHO  = tar_rho
  tar.points.VX   = tar_vel[0]
  tar.points.VY   = tar_vel[1]
  tar.points.VZ   = tar_vel[2]
  tar.check()

  # sphere impactor: m-kg-s
  imp_cell_r = 25
  imp_radius = 100
  imp_rho = 2700
  imp_vmod  = 5000
  imp_angle  = np.deg2rad(45) # from surface normal
  imp_center = np.array([tar_center[0]+tar_radius+imp_radius+8*imp_cell_r*np.cos(imp_angle),8*imp_cell_r*np.sin(imp_angle),0])
  imp_vel    = imp_vmod*np.array([-np.cos(imp_angle),-np.sin(imp_angle),0])

  # create impactor particles
  imp = part()
  box_pos = generate_box(box_size=imp_radius*2*np.array([1,1,1]), box_center=imp_center, r=imp_cell_r, lattice='grid')
  box_pos = rand_rotate(box_pos, imp_center)
  def imp_geometry(x,y,z): 
    return (x-imp_center[0])**2+(y-imp_center[1])**2+(z-imp_center[2])**2 <= imp_radius**2
  imp.create_particles(box_pos,imp_geometry)
  imp.points.PART = 1
  imp.points.MASS = imp_rho*8*imp_cell_r**3
  imp.points.RHO  = imp_rho
  imp.points.VX   = imp_vel[0]
  imp.points.VY   = imp_vel[1]
  imp.points.VZ   = imp_vel[2]
  imp.check()


  print("Initial hsml = %e" % (2.0*eta*tar_cell_r))
  print("Initial dt = %e" % (0.2*2.0*eta*tar_cell_r/5000))

  # assemble parts
  parts = {}
  parts[0] = tar
  parts[1] = imp

  return parts


def main():
  print("Start")

  parts = set_parts_geometry()
  visualize_particles(parts)
  export_particles(parts, "../../Input/0-Particles.csv", columns=col)

  print("Done")


if __name__ == "__main__":
  main()
