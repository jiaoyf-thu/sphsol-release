from create_particles import *


'''
Set up the impact geometry
'''


eta = 1.2
col = ['X','Y','Z','VX','VY','VZ','PART','MASS','RHO','D']


# set parts geometry
def set_parts_geometry():
  # trimesh target: m-kg-s
  tar_cell_r = 216 # 270
  tar_rho = 1300
  # tar_center = np.array([0,0,0])
  tar_vel = np.array([0,0,0])

  # create target particles
  tar = part()
  # Load the input mesh as a list of triplets (ie. triangles) of 3d vertices using open3d
  mesh = o3d.io.read_triangle_mesh("MathildeSH.ply")
  trimesh_mesh = trimesh.Trimesh(vertices=np.asarray(mesh.vertices)*1000.0, faces=np.asarray(mesh.triangles))
  # Compute uniform distribution within the axis-aligned bound box for the mesh
  min_corner = np.amin(trimesh_mesh.vertices, axis = 0)
  max_corner = np.amax(trimesh_mesh.vertices, axis = 0)
  mesh_center = 0.5*(min_corner+max_corner)
  box_size = max(max_corner-min_corner)
  box_pos = generate_box(box_size=np.array([1,1,1])*box_size, box_center=mesh_center, r=tar_cell_r, lattice='psd')
  box_pos = rand_rotate(box_pos, mesh_center)
  # Check if the points are inside the mesh
  box_pos = box_pos[trimesh_mesh.contains(box_pos[["x","y","z"]].values)].copy().reset_index(drop=True)
  tar.create_particles(box_pos)
  tar.points.PART = 0
  tar.points.MASS = tar_rho*8*tar_cell_r**3
  tar.points.RHO  = tar_rho
  tar.points.VX   = tar_vel[0]
  tar.points.VY   = tar_vel[1]
  tar.points.VZ   = tar_vel[2]
  tar.points.D    = 1.0
  tar.check()

  # sphere impactor: m-kg-s
  imp_cell_r = 216 # 270
  imp_radius = 1080
  imp_rho = 1300
  imp_vmod  = 5000
  # imp_center = np.array([-1.689e3, -19.308e3, 13.571e3])
  # imp_center = np.array([0., 10.266e3, 22.029e3])
  imp_center = np.array([0.729e3, 7.047e3, 23.085e3]) # calculate from lat. lon.
  imp_center *= 1.0 + (imp_radius + 10*imp_cell_r)/np.linalg.norm(imp_center)
  imp_vel = imp_vmod * np.array([-imp_center[0],-imp_center[1],-imp_center[2]]) / np.linalg.norm(imp_center)

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
  imp.points.D    = 1.0
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
