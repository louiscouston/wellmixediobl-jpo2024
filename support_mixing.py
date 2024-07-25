import os
import numpy as np

##############################
# Constant physical parameters
##############################
Li = 3.35*1e5
cp = 3974
Pt = 1000 # dbar ~ meters of seawater
lambda1 = -5.73*1e-2 # K 
lambda2 = 8.32*1e-2 # K
lambda3 = -7.53*1e-4 # K/dbar
S0 = 35 # ppt
alpha = 3.87*1e-5 # 1/K
beta = 7.86*1e-4 # 1/ppt
H = 1 # m
nu = 1.8*1e-6 # m^2/s

#print('St*1e-5 is',Li/(cp*1e-5)*1e-5,Li/(cp*10)*1e-5)
#print('gamma is',lambda1*S0/1e-5,lambda1*S0/10)

##############################
# Variable physical parameters
##############################
def varphyparams(Re, At, Pr, Le, TD):
  # Dimensional parameters
  g = 9.81/At # Gravity is reduced!
  kapT = nu/Pr # requires doubling of resolution
  kapS = kapT/Le # also requires doubling of resolution
  # Dimensionless parameters
  Sc = Pr*Le
  RaS = g*beta*S0*H**3/(nu*kapS)
  Rrho = alpha*TD/(beta*S0)
  St = Li/(cp*TD)
  sigma = lambda2/TD+lambda3*Pt/TD
  gamma = lambda1*S0/TD 
  Ri = RaS/Sc/Re**2
  return Sc, RaS, Rrho, St, sigma, gamma, Ri  

#print(varphyparams(200, 1, 1, 1, 1e-5)[2])
#print(varphyparams(200, 1, 1, 1, 10)[2])


#print(varphyparams(200, 1, 1, 1, 1e-5)[-1])
#print(varphyparams(200, 1, 10, 10, 1e-5)[-1])
#print(varphyparams(200, 1, 1, 100, 1e-5)[-1])
#print(varphyparams(400, 1, 1, 1, 1e-5)[-1])
#print(varphyparams(400, 1, 1, 10, 1e-5)[-1])
#print(varphyparams(400, 1, 10, 10, 1e-5)[-1])
#print(varphyparams(800, 1, 1, 1, 1e-5)[-1])
#print(varphyparams(800, 1, 1, 10, 1e-5)[-1])

###############################
# Variable numerical parameters
###############################
def varnumparams(Re, At, Pr, Le, TD, size):
  # Length and resolution
  Sc = Pr*Le
  if Re==200:
    Lx, Ly, Lz = (2, 1, 1)
    check_dt = 1 # number of friction time scales between checkpoints
    if Sc==1: Nx, Ny, Nz = (64, 64, 48)
    elif Sc==10: Nx, Ny, Nz = (128, 128, 96)
    elif Sc==25 or Sc==30 or Sc==100: Nx, Ny, Nz = (256, 256, 192)
  if Re==400:
    Lx, Ly, Lz = (1, 0.5, 1)
    check_dt = 1 # number of friction time scales between checkpoints
    if Sc==1: Nx, Ny, Nz = (64, 64, 96)
    elif Sc==10: Nx, Ny, Nz = (128, 128, 192) 
    elif Sc==100: Nx, Ny, Nz = (256, 256, 384)
  if Re==800:
    Lx, Ly, Lz = (0.5, 0.25, 1)
    check_dt = 1 # number of friction time scales between checkpoints
    if Sc==1: Nx, Ny, Nz = (64, 64, 128) # 8 Cascade @64
    elif Sc==10: Nx, Ny, Nz = (128, 128, 256)
    elif Sc==100: Nx, Ny, Nz = (256, 256, 512)
  # Parallel mesh distribution -- assuming 2^n
  mesh = (2**(np.log2(size)//2), 2**(np.log2(size)-np.log2(size)//2))
  # Create directory if new
  name_dir = 'Re_' + str(int(Re)) + '_At_' + str(int(At)) + '_Pr_' + str(int(Pr)) + '_Le_' + str(int(Le)) + '_TD100000_' + str(int(TD*100000)) #+ '_np_' + str(int(size)) 
  try:
    os.makedirs(name_dir, exist_ok = True)
  except OSError as error:
    print(error) 
  return Lx, Ly, Lz, Nx, Ny, Nz, check_dt, mesh, name_dir
  
###############################
# Variable damping parameters
###############################
def vardamparams(Re):
  damping_depth = 0.2
  transition_length = 0.025
  damping_time_scale = 1/Re**(3/2) 
  return damping_depth, transition_length, damping_time_scale
  

def createmap(n,cmap):
  mapvec = []
  for j in range(0,n):
    idd = np.rint(j*256/(n-1))
    idd = idd.astype(int)
    mapvec.append(cmap(idd))
  return mapvec  

  
  
