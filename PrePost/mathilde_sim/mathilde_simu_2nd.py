from create_particles import *


'''
Set up the impact geometry for a second time simulation
'''


eta = 1.2
col = ['X','Y','Z','VX','VY','VZ','PART','MASS','RHO', 'D','P', 'ALPHA']


home_path = 'D:/jiaoyf/Seafile/SPH/sphsol-para/Data/mathilde_simu/paramatrix1/'
# root_path = ['c000-f04/','c000-f07/','c000-f10/','c1e3-f04/','c1e3-f07/','c1e3-f10/',
#              'c1e4-f04/','c1e4-f07/','c1e4-f10/','c1e5-f04/','c1e5-f07/','c1e5-f10/',
#              'c1e6-f04/','c1e6-f07/','c1e6-f10/','no-strength/']
# root_path = ['c1e4-f10-ps1e7/','c1e4-f10-ps1e9/']
root_path = ['c1e4-f10-hres/']

init_setting = 'D:/jiaoyf/Seafile/SPH/sphsol-para/Data/mathilde_simu/paramatrix1/1-Settings.txt'


# set parts geometry
def set_parts_geometry(path):
  # previous result as the target: m-kg-s
  tar_cell_r = 216 # 270 # 360
  tar_particles = pd.read_csv(path)

  # set a domain for the target
  tar_particles = tar_particles[tar_particles.X**2+tar_particles.Y**2+tar_particles.Z**2 < 40000**2].copy().reset_index(drop=True)
  tar_particles.VX = 0
  tar_particles.VY = 0
  tar_particles.VZ = 0
  # find the largest remanent
  clustering = dbscan.DBSCAN(eps=3.0*tar_cell_r,min_samples=3).fit(tar_particles[['X','Y','Z']].values)
  core_label = np.argmax(np.bincount(clustering.labels_[clustering.core_sample_indices_]))
  tar_particles = tar_particles.loc[np.where(clustering.labels_==core_label)[0]].copy().reset_index(drop=True) # largest remanent
  # if 'PART' not in tar_particles.columns: tar_particles['PART'] = 0
  tar_particles['PART'] = 0
  if 'D' not in tar_particles.columns: tar_particles['D'] = 1.0

  # copy target particles
  tar = part()
  tar.points[col] = tar_particles[col].copy()
  tar.check()

  # sphere impactor: m-kg-s
  imp_cell_r = 216 # 270 # 360
  imp_radius = 1350 # 1200
  imp_rho    = 1300
  # imp_angle  = np.deg2rad(70)
  imp_vmod   = 5000
  # imp_center = np.array([0., 10.266e3, 22.029e3])
  # imp_center = np.array([-1.689e3, -19.308e3, 13.571e3])
  imp_center = np.array([-2.832e3, -19.352e3, 12.98e3]) # calculate from lat. lon.
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
  imp.points.P    = 0.0
  imp.points.ALPHA = 2.0
  imp.check()

  # print("Initial hsml = %e" % (2.0*eta*tar_cell_r))
  # print("Initial dt = %e" % (0.2*2.0*eta*tar_cell_r/5000))

  # assemble parts
  parts = {}
  parts[0] = tar
  parts[1] = imp

  return parts


def main():
  print("Start")

  for root in root_path:
    print(root)
    path = home_path + root + 'Output/rot_Particles0030.csv'
    parts = set_parts_geometry(path)
    # visualize_particles(parts)
    pnum_ = export_particles(parts, home_path + root + '0-Particles.csv', columns=col)

    # replace the setting file
    revised_setting = home_path + root + '1-Settings.txt'
    # Process file with line numbers
    replacement_lines = {27: 'PARTICLES_NUM           [I] : %d\n' % pnum_}
    with open(init_setting, 'r') as infile, open(revised_setting, 'w') as outfile:
      for line_number, line in enumerate(infile, start=1):
        if line_number in replacement_lines:
          outfile.write(replacement_lines[line_number])
        else:
          outfile.write(line)

  print("Done")


if __name__ == "__main__":
  main()
