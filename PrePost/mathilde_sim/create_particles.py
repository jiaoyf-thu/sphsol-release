import numpy as np
import pandas as pd
import open3d as o3d
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
import trimesh
import sklearn.cluster._dbscan as dbscan
from scipy.spatial import KDTree
from scipy.spatial.distance import cdist,euclidean
from numpy.random import default_rng


'''
To create particles for sphsol input.
'''


# Hexagonal close packing (HCP) lattice
def hcp(n):
  k, j, i = [v.flatten() for v in np.meshgrid([range(n[2])],[range(n[1])],[range(n[0])], indexing='ij')]
  df = pd.DataFrame({
    'x': 2 * i + (j + k) % 2,
    'y': np.sqrt(3) * (j + 1/3 * (k % 2)),
    'z': 2 * np.sqrt(6) / 3 * k,})
  return df


# Grid lattice
def grid(n):
  k, j, i = [v.flatten() for v in np.meshgrid([range(n[2])],[range(n[1])],[range(n[0])], indexing='ij')]
  df = pd.DataFrame({
    'x': i,
    'y': j,
    'z': k,})
  return df


class ParametrizedSpiralingDistribution:
  def __init__(self, seed: int):
      self.seed = seed

  def generate(self, particle_radius: float, center: np.ndarray, radius: float) -> np.ndarray:
    """
    Generate a spiraling particle distribution efficiently on a single CPU.
    
    Parameters:
        particle_radius (float): Radius of each particle.
        center (np.ndarray): 3D coordinates of the center [x, y, z].
        radius (float): Maximum radial distance from the center.

    Returns:
        np.ndarray: Array of particle positions, shape (n_particles, 3).
    """
    # Volume of each particle: 8 * r^3
    particle_volume = 8.0 * particle_radius**3

    # Volume of the sphere
    sphere_volume = (4. / 3.) * np.pi * radius**3

    # Total number of particles
    n = int(sphere_volume / particle_volume)

    # Interparticle distance (center-to-center spacing)
    h = 2.0 * particle_radius

    # Number of shells
    num_shells = int(np.ceil(radius / h))

    # Shell particle counts normalized by shell weights
    shell_weights = np.array([((i + 1) * h) ** 2 for i in range(num_shells)])
    shell_weights /= np.sum(shell_weights)  # Normalize
    shell_counts = (n * shell_weights).astype(int)

    rng = default_rng(self.seed)
    positions = []

    for shell_idx, r in enumerate(np.arange(h, radius + h, h)):
      m = shell_counts[shell_idx]
      if m == 0:
          continue  # Skip empty shells

      # Generate angles in batches
      k_vals = np.arange(1, m)
      hk = -1.0 + 2.0 * k_vals / m
      theta = np.arccos(hk)
      phi = np.cumsum(3.8 / np.sqrt(m * (1.0 - hk**2)))

      # Convert to Cartesian coordinates
      x, y, z = self.spherical_to_cartesian_batch(r, theta, phi)

      # Apply random rotations to the shell
      rotation_matrix = self.random_rotation_matrix(rng)
      rotated_positions = rotation_matrix @ np.vstack((x, y, z))

      # Shift positions to the center
      shell_positions = rotated_positions.T + center
      positions.append(shell_positions)

    # Combine all positions into a single array
    return pd.DataFrame(np.vstack(positions), columns=['x', 'y', 'z'])

  @staticmethod
  def spherical_to_cartesian_batch(r, theta, phi):
    """
    Vectorized conversion of spherical to Cartesian coordinates.
    """
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x, y, z

  @staticmethod
  def random_rotation_matrix(rng):
    """
    Generate a random 3D rotation matrix.
    """
    axis = rng.normal(size=3)
    axis /= np.linalg.norm(axis)
    angle = 2.0 * np.pi * rng.random()
    return ParametrizedSpiralingDistribution.rotate_axis(axis, angle)

  @staticmethod
  def rotate_axis(axis, angle):
    """
    Generate a rotation matrix for rotating around a given axis.
    """
    x, y, z = axis
    c = np.cos(angle)
    s = np.sin(angle)
    t = 1 - c
    return np.array([
        [t * x * x + c,     t * x * y - s * z, t * x * z + s * y],
        [t * x * y + s * z, t * y * y + c,     t * y * z - s * x],
        [t * x * z - s * y, t * y * z + s * x, t * z * z + c]
    ])



# generate points in a given box
def generate_box(box_size, box_center, r, lattice):
  if lattice == "hcp":
    r = r * pow(2,1/6)
    nsize = [int(np.around(box_size[0]/(r*2.))),int(np.around(box_size[1]/(r*np.sqrt(3)))),int(np.around(box_size[2]/(r*2.*np.sqrt(6)/3.)))]
    df = hcp(nsize)*r
    center_ = np.array([(np.max(df["x"])+np.min(df["x"]))/2.0,(np.max(df["y"])+np.min(df["y"]))/2.0,(np.max(df["z"])+np.min(df["z"]))/2.0])
  elif lattice == "grid":
    nsize = [int(np.around(box_size[0]/(r*2.))),int(np.around(box_size[1]/(r*2.))),int(np.around(box_size[2]/(r*2.)))]
    df = grid(nsize)*r*2
    center_ = np.array([(np.max(df["x"])+np.min(df["x"]))/2.0,(np.max(df["y"])+np.min(df["y"]))/2.0,(np.max(df["z"])+np.min(df["z"]))/2.0])
  elif lattice == "psd":
    df = ParametrizedSpiralingDistribution(seed=0).generate(r, box_center, max(box_size)/2.0)
    center_ = np.array([0.0,0.0,0.0])
  df["x"] += box_center[0]-center_[0]
  df["y"] += box_center[1]-center_[1]
  df["z"] += box_center[2]-center_[2]
  return df



# random rotate particles (dataframe) around center
def rand_rotate(xyz, center):
  xyz = xyz.copy()
  xyz["x"] -= center[0]
  xyz["y"] -= center[1]
  xyz["z"] -= center[2]
  rot = Rotation.random().as_matrix()
  xyz[["x","y","z"]] = np.dot(xyz[["x","y","z"]], rot.T)
  xyz["x"] += center[0]
  xyz["y"] += center[1]
  xyz["z"] += center[2]
  return xyz

# part class
class part():
  def __init__(self):
    self.points = pd.DataFrame(columns=['X','Y','Z','VX','VY','VZ','PART','MASS','RHO','IFLAG','D','P','ALPHA'])
    self.pnum   = 0

  # create particles within given geometry, which returns a boolean array
  def create_particles(self, box_pos, geometry=None, if_print=True):
    if geometry != None:
      pos = box_pos[geometry(box_pos["x"], box_pos["y"], box_pos["z"])].copy()
    else:
      pos = box_pos.copy()
    pos = pos.reset_index(drop=True)
    self.points.X = pos.x
    self.points.Y = pos.y
    self.points.Z = pos.z

  # check part para
  def check(self, if_print=True):
    self.pnum = len(self.points)
    # set all columns to float
    self.points = self.points.astype(float)
    # set part to int
    self.points.PART = self.points.PART.astype(int)
    if not self.points.IFLAG.isnull().all():
      self.points.IFLAG = self.points.IFLAG.astype(int)
    if if_print:
      print(" Part %d includes %8d sph particles, total mass = %.4e" % (self.points.PART[0], self.pnum, self.points.MASS.sum()))

  # add another part
  def add_part(self, _part, if_print=True):
    self.points = pd.concat([self.points, _part.points], ignore_index=True)
    self.check(if_print)

  # export part to csv
  def export(self, path, mode, columns=None, start_idx=0):
    self.points.index = range(start_idx, start_idx+self.pnum)
    if mode == "w":
      self.points.to_csv(path, mode=mode, columns=columns, float_format='%.6e')
    if mode == "a":
      self.points.to_csv(path, mode=mode, columns=columns, float_format='%.6e', header=False)


# visualize particles
def custom_draw_geometry_with_custom_fov(pcd, fov_step):
    vis = o3d.visualization.Visualizer()
    vis.create_window(width=800, height=600)
    for key in pcd.keys():
      vis.add_geometry(pcd[key])
    ctr = vis.get_view_control()
    print("Field of view (before changing) %.2f" % ctr.get_field_of_view())
    ctr.change_field_of_view(step=fov_step)
    print("Field of view (after changing) %.2f" % ctr.get_field_of_view())
    vis.run()
    vis.destroy_window()

def visualize_particles(parts):
  pcd = {}
  colors = plt.cm.viridis(np.linspace(0, 0.9, len(parts.keys())))
  for i, key in enumerate(parts.keys()):
      pcd[key] = o3d.geometry.PointCloud()
      pcd[key].points = o3d.utility.Vector3dVector(parts[key].points[['X','Y','Z']].values)
      pcd[key].paint_uniform_color(colors[i, :3])
  custom_draw_geometry_with_custom_fov(pcd, -90.0)


# export particles
def export_particles(parts, path, columns):
  start_idx = 0
  for key in parts.keys():
    if key == 0:
      parts[key].export(path, "w", columns, start_idx)
    else:
      start_idx += parts[key-1].pnum
      parts[key].export(path, "a", columns, start_idx)

  pnum_ = sum([parts[key].pnum for key in parts.keys()])
  print("Total particles num = %d" % (pnum_))

  return pnum_

  # with open(path) as f:
  #   f.write("ID,X,Y,Z,VX,VY,VZ,PART,MASS,DENSITY\n")
  #   pid = 0
  #   for key in parts.keys():
  #     for i in range(parts[key].pnum):
  #       f.write("%d,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%d,%.6e,%.6e\n" % \
  #         (pid, parts[key].pos.loc[i,"x"],parts[key].pos.loc[i,"y"],parts[key].pos.loc[i,"z"], \
  #         parts[key].vel[0],parts[key].vel[1],parts[key].vel[2], parts[key].part_id, \
  #         parts[key].mass, parts[key].density))
  #       pid += 1
  # print("Total particles num = %d" % (pid))



'''
SPH cubic spline kernel
'''
def kernel(rij, h):
  h1 = 1. / h
  q = rij * h1
  # get the kernel normalizing factor
  fac = 1.0 / np.pi * h1 * h1 * h1
  tmp2 = 2. - q
  if q > 2.0:
    val = 0.0
  elif q > 1.0:
    val = 0.25 * tmp2 * tmp2 * tmp2
  else:
    val = 1 - 1.5 * q * q * (1 - 0.5 * q)
  return val * fac