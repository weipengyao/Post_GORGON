import vtk
import vtk.util.numpy_support as VN
import numpy as np
import matplotlib as mpl
# mpl.use('pdf')
from pyexp import normalize, dot, cross, magnitude
from pylab import imshow, colorbar, clf

class PVTI:

	def __init__(self,wrkdir,pfx,tm):
		self.workdir = wrkdir
		self.prefix = pfx
		self.time = tm
		self.filename=None
# 		self.d = None
# 		self.V = None
# 		self.B = None
# 		self.Te = None
# 		self.Ti = None
# 		self.Z = None
# 		self.J = None
# 		self.E = None
	
	def source(self):
		"""Display data source full path"""
		print("Data Source: ", self.filename)

	def CellSize(self):
		"""Returns the cell size (It only returns on value as it assumes
		the cells are cubic."""
		return self.mesh.GetSpacing()[0]


# 	def GetData(self, density=True, V=True, B=True, Te=True, Ti=True, Z=True, J=True, E=True, coord=True):
# 		"""Get data from pvti file"""
# 		self.filename = self.workdir + self.prefix + "_3D-all."+ str(self.time) + ".pvti"
# 		reader = vtk.vtkXMLPImageDataReader()
# 		reader.SetFileName(self.filename)
# 		reader.Update()
# 		self.mesh = reader.GetOutput()
# 		# Get cell data from the mesh
# 		self.celldata = self.mesh.GetCellData()
# 		#Get Cell Centres
# 		cc = vtk.vtkCellCenters()
# 		cc.SetInput(self.mesh)
# 		cc.Update()
# 		self.mycc = cc.GetOutput()
# 		if(density):
# 		  self.d = VN.vtk_to_numpy(self.celldata.GetArray("rho"))
# 		if(V):
# 		  self.V = VN.vtk_to_numpy(self.celldata.GetArray("V"))
# 		if(B):
# 		  self.B = VN.vtk_to_numpy(self.celldata.GetArray("B"))
# 		if(Te):
# 		  self.Te = VN.vtk_to_numpy(self.celldata.GetArray("Te"))
# 		if(Ti):
# 		  self.Ti = VN.vtk_to_numpy(self.celldata.GetArray("Ti"))
# 		if(Z):
# 		  self.Z = VN.vtk_to_numpy(self.celldata.GetArray("Zeff"))
# 		if(J):
# 		  self.J = VN.vtk_to_numpy(self.celldata.GetArray("J"))
# 		if(E):
# 		  self.E = VN.vtk_to_numpy(self.celldata.GetArray("E"))
# 		self.coord = VN.vtk_to_numpy(self.mycc.GetPoints().GetData())

	def GetData(self,density=True, V=False, B=False, Te=False, Ti=False, Z=False, J=False, E=False, material1=False, material2=False, coord=False, emission=False, tracer=False, tracer2=False, laser_tracer=False, energy_tracer=False):
         """Get data from pvti file"""
         self.filename = self.workdir+self.prefix+"_3D-all."+str(self.time)+".pvti"#POUR LES FICHIER PVTI
#		self.filename=self.workdir + self.prefix + str(self.time) + ".vti" #POUR LES FICHIER VTI
         reader = vtk.vtkXMLPImageDataReader() #POUR LIRE DES FICHIERS PVTI
#		reader = vtk.vtkXMLImageDataReader()  #POUR LIRE LES FICHIERS VTI
         reader.SetFileName(self.filename)
         reader.Update()
#
         if(not density):
             reader.SetCellArrayStatus('rho',0)
         if(not V):
             reader.SetCellArrayStatus('V',0)
         if(not B):
             reader.SetCellArrayStatus('B',0)
         if(not Te):
             reader.SetCellArrayStatus('Te',0)
         if(not Ti):
             reader.SetCellArrayStatus('Ti',0)
         if(not Z):
             reader.SetCellArrayStatus('Zeff',0)
         if(not J):
             reader.SetCellArrayStatus('J',0)
         if(not E):
             reader.SetCellArrayStatus('E',0)
#
         reader.Update()
         self.mesh = reader.GetOutput()
         # Get cell data from the mesh
         self.celldata = self.mesh.GetCellData()
         #Get Cell Centres
         cc = vtk.vtkCellCenters()
         cc.SetInputData(self.mesh)
         cc.Update()
         self.mycc = cc.GetOutput()
         if(density):
             self.d = VN.vtk_to_numpy(self.celldata.GetArray("rho"))
         if(V):
             self.V = VN.vtk_to_numpy(self.celldata.GetArray("V"))
         if(B):
             self.B = VN.vtk_to_numpy(self.celldata.GetArray("B"))
         if(Te):
             self.Te = VN.vtk_to_numpy(self.celldata.GetArray("Te"))
         if(Ti):
             self.Ti = VN.vtk_to_numpy(self.celldata.GetArray("Ti"))
         if(Z):
             self.Z = VN.vtk_to_numpy(self.celldata.GetArray("Zeff"))
         if(J):
             self.J = VN.vtk_to_numpy(self.celldata.GetArray("J"))
         if(E):
             self.E = VN.vtk_to_numpy(self.celldata.GetArray("E"))
#
         self.coord = VN.vtk_to_numpy(self.mycc.GetPoints().GetData())
#
         if(material1):
             tmpfile = self.workdir + self.prefix + "_3D-carbon_rho."+ str(self.time) + ".pvti"
             material1_reader = vtk.vtkXMLPImageDataReader()
             material1_reader.SetFileName(tmpfile)
             material1_reader.Update()
             self.carbon_rho = VN.vtk_to_numpy(material1_reader.GetOutput().GetCellData().GetArray("carbon_rho"))
         if(material2):
             tmpfile = self.workdir + self.prefix + "_3D-gas_band_rho."+ str(self.time) + ".pvti"
             material2_reader = vtk.vtkXMLPImageDataReader()
             material2_reader.SetFileName(tmpfile)
             material2_reader.Update()
             self.gas_band_rho = VN.vtk_to_numpy(material2_reader.GetOutput().GetCellData().GetArray("gas_band_rho"))

         if(emission):
             tmpfile = self.workdir + self.prefix + "_3D-em."+ str(self.time) + ".pvti"
             emreader = vtk.vtkXMLPImageDataReader()
             emreader.SetFileName(tmpfile)
             emreader.Update()
             self.em = VN.vtk_to_numpy(emreader.GetOutput().GetCellData().GetArray("em"))
#
         if(tracer):
             tmpfile = self.workdir + self.prefix + "_3D-trac."+ str(self.time) + ".pvti"
             tracreader = vtk.vtkXMLPImageDataReader()
             tracreader.SetFileName(tmpfile)
             tracreader.Update()
             self.trac = VN.vtk_to_numpy(tracreader.GetOutput().GetCellData().GetArray("trac"))
#
         if(tracer2):
             tmpfile = self.workdir + self.prefix + "_3D-trac2."+ str(self.time) + ".pvti"
             tracreader = vtk.vtkXMLPImageDataReader()
             tracreader.SetFileName(tmpfile)
             tracreader.Update()
             self.trac2 = VN.vtk_to_numpy(tracreader.GetOutput().GetCellData().GetArray("trac2"))
#
         if(laser_tracer):
             tmpfile = self.workdir + self.prefix + "_3D-laser_trac."+ str(self.time) + ".pvti"
             tracreader = vtk.vtkXMLPImageDataReader()
             tracreader.SetFileName(tmpfile)
             tracreader.Update()
             self.laser_tracer = VN.vtk_to_numpy(tracreader.GetOutput().GetCellData().GetArray("laser_trac"))
#
         if(energy_tracer):
             tmpfile = self.workdir + self.prefix + "_3D-energy_e."+ str(self.time) + ".pvti"
             tracreader = vtk.vtkXMLPImageDataReader()
             tracreader.SetFileName(tmpfile)
             tracreader.Update()
             self.energy_e = VN.vtk_to_numpy(tracreader.GetOutput().GetCellData().GetArray("energy_e"))   


#C>============================================================================

	def GetDataMultiMaterial(self, mat1_rho, mat1_Z, mat2_rho, mat2_Z, mat3_rho, mat3_Z, material1=False, material2=False, material3=False):
         """Function to get the data from different materials used in a GORGON simulation.
         Default number of materials is 3, thus it is needed to specify even if the simulation has run with less --> Needed to find a mean to fix that
         The names used when calling (e.g. carbon_rho & carbon_zbar) have to be the ones used by the GORGON code."""

#C> - Section needed to get domain informations used to reshape 1D arrays into 3D cube sets one can deal with.
#C>-----------------------------------------------------------------------------------------------------------
         self.filename = self.workdir+self.prefix+"_3D-all."+str(self.time)+".pvti"     #C> - PVTI file name for box infos.

         reader = vtk.vtkXMLPImageDataReader()                                          #C> - Reader operator for PVTI data files.
         reader.SetFileName(self.filename)
         reader.Update()

         self.mesh = reader.GetOutput()                                        #C> - Domain meshes info
         self.celldata = self.mesh.GetCellData()                               #C> - Domain cell centres info
         cc = vtk.vtkCellCenters()
         cc.SetInputData(self.mesh)
         cc.Update()
         self.mycc = cc.GetOutput()
         self.coord = VN.vtk_to_numpy(self.mycc.GetPoints().GetData())         #C> Domain grid points info


#C> - Materials mass density and ionization rate extraction
#C> -------------------------------------------------------
         if(material1):
             tmpfile = self.workdir + self.prefix + "_3D-" + mat1_rho + "." + str(self.time) + ".pvti"
             material1_reader = vtk.vtkXMLPImageDataReader()
             material1_reader.SetFileName(tmpfile)
             material1_reader.Update()
             self.mat1_rho = VN.vtk_to_numpy(material1_reader.GetOutput().GetCellData().GetArray(mat1_rho))

             tmpfile = self.workdir + self.prefix + "_3D-" + mat1_Z + "." + str(self.time) + ".pvti"
             material1_reader = vtk.vtkXMLPImageDataReader()
             material1_reader.SetFileName(tmpfile)
             material1_reader.Update()
             self.mat1_Z = VN.vtk_to_numpy(material1_reader.GetOutput().GetCellData().GetArray(mat1_Z))
#C> -------------------
         if(material2):
             tmpfile = self.workdir + self.prefix + "_3D-" + mat2_rho + "." + str(self.time) + ".pvti"
             material2_reader = vtk.vtkXMLPImageDataReader()
             material2_reader.SetFileName(tmpfile)
             material2_reader.Update()
             self.mat2_rho = VN.vtk_to_numpy(material2_reader.GetOutput().GetCellData().GetArray(mat2_rho))

             tmpfile = self.workdir + self.prefix + "_3D-" + mat2_Z + "." + str(self.time) + ".pvti"
             material2_reader = vtk.vtkXMLPImageDataReader()
             material2_reader.SetFileName(tmpfile)
             material2_reader.Update()
             self.mat2_Z = VN.vtk_to_numpy(material2_reader.GetOutput().GetCellData().GetArray(mat2_Z))
#C> -------------------
         if(material3):
             tmpfile = self.workdir + self.prefix + "_3D-" + mat3_rho + "." + str(self.time) + ".pvti"
             material3_reader = vtk.vtkXMLPImageDataReader()
             material3_reader.SetFileName(tmpfile)
             material3_reader.Update()
             self.mat3_rho = VN.vtk_to_numpy(material3_reader.GetOutput().GetCellData().GetArray(mat3_rho))

             tmpfile = self.workdir + self.prefix + "_3D-" + mat3_Z + "." + str(self.time) + ".pvti"
             material3_reader = vtk.vtkXMLPImageDataReader()
             material3_reader.SetFileName(tmpfile)
             material3_reader.Update()
             self.mat3_Z = VN.vtk_to_numpy(material3_reader.GetOutput().GetCellData().GetArray(mat3_Z))


#C>============================================================================

	def cube(self,variable):
		"""Returns a 3D data cube from the 1D array variable"""
		(x1,x2,y1,y2,z1,z2) = self.mesh.GetExtent()
		return np.reshape(variable,(z2,y2,x2)).transpose((2,1,0))
		
	def cart2cyl(self,B):
		"""Returns the cylindrical components of the vector B"""
		R,Phi,Z = self.CylindricalTriad()
		Br = dot(B,R)
		Bphi = dot(B,Phi)
		Bz = dot(B,Z)
		return Br,Bphi,Bz
	
	def GetCentredCoord(self):
	  """Returns an array with the centred coordinates"""
	  xc = self.mycc.GetCenter()[0]
	  yc = self.mycc.GetCenter()[1]
	  return  self.coord - (xc,yc,0.0)


	def CylindricalTriad(self):
		"""Returns the cylindrical triad at each point on the grid"""
		cc = self.GetCentredCoord()
		npoints = np.shape(self.coord)[0]
		
		R = np.zeros([npoints,3])
		R[:,0] = cc[:,0]
		R[:,1] = cc[:,1]
		R = normalize(R)
		
		Z = np.zeros([npoints,3])
		Z[:,2] = 1.0
		
		Phi = cross(Z,R)
		
		return R,Phi,Z
		
	def SphericalR(self):
		""""Returns the spherical radial UNIT vector"""
		return normalize(self.GetCentredCoord())

	def CylindricalCoordinates(self):
		"""Returns the cylindrical coordinates at each point on the grid"""
		#First centre the cartesian grid.
		cc = self.GetCentredCoord()
		varpi = np.sqrt(cc[:,0]**2 + cc[:,1]**2)
		phi = np.pi+np.arctan2(cc[:,1],cc[:,0])
		
		return varpi,phi,cc[:,2]

	def MaskDensity(self,value):
		"""Returns a bool array with true where density > value"""
		return self.d > value

	#ORIGINAL FUNCTION, FULLY WORKING ON 22-03-12
# 	def slice(self,variable,point,direction = "X", low = True, up = True, \
# 	  linlog = False, colbar = True, interp = "nearest", clear = True, scale=1.0):
# 		"""Slices either a 3D cube or displays an existing 2D matrix (e.g areal density, see direction directive).
# 		"""
# 		if (clear):
# 		  clf()
		
# 		bounds = self.mycc.GetBounds()
# 		xmin = bounds[0] * scale
# 		xmax = bounds[1] * scale
# 		ymin = bounds[2] * scale
# 		ymax = bounds[3] * scale
# 		zmin = bounds[4] * scale
# 		zmax = bounds[5] * scale

# 		if(linlog):
# 			variable = np.log10(variable)

# 		if(low != True):	
# 			variable = np.ma.masked_where(variable < low, variable)
# 		if(up != True):
# 			variable = np.ma.masked_where(variable > up, variable)

# 		if (direction == "X"):
# 		  extent = [ymin-ymax/2.0,ymax-ymax/2.0,zmin,zmax]
# 		  imshow(np.rot90(self.cube(variable)[point,:,:],1),extent=extent, interpolation = interp)
# 		elif (direction == "Y"):
# 		  extent = [xmin-xmax/2.0,xmax-xmax/2.0,zmin,zmax]
# 		  imshow(np.rot90(self.cube(variable)[:,point,:],1),extent=extent, interpolation = interp)
# 		elif (direction == "Z"):
# 		  extent = [xmin-xmax/2.0,xmax-xmax/2.0,ymin-ymax/2.0,ymax-ymax/2.0]
# 		  imshow(self.cube(variable)[:,:,point],extent=extent, interpolation = interp)
# 		elif (direction == "2Dside"):
# 		  extent = [ymin-ymax/2.0,ymax-ymax/2.0,zmin,zmax]
# 		  imshow(np.rot90(variable),extent=extent, interpolation = interp)
# 		elif (direction == "2Dtop"):
# 		  extent = [xmin-xmax/2.0,xmax-xmax/2.0,ymin-ymax/2.0,ymax-ymax/2.0]
# 		  imshow(variable,extent=extent, interpolation = interp)
		
		
# 		if(colbar):
# 		  colorbar()
	
	#Improved version. Doesn't need any masked arrays. Uses imshow functionalities
	def slice(self,variable,point,direction = "X", low = None, up = None, \
	  linlog = False, colbar = True, colmap = mpl.cm.jet, interp = "nearest", clear = True, scale=1.0, collow='w', colup='k'):
		"""Slices either a 3D cube or displays an existing 2D matrix (e.g areal density, see direction directive).
		"""
		if (clear):
		  clf()
		
		bounds = self.mesh.GetBounds()
		xmin = bounds[0] * scale
		xmax = bounds[1] * scale
		ymin = bounds[2] * scale
		ymax = bounds[3] * scale
		zmin = bounds[4] * scale
		zmax = bounds[5] * scale

		if(linlog):
			variable = np.log10(variable)
		
		cmap = colmap #mpl.cm.jet
		cmap.set_over(colup)
		cmap.set_under(collow)
		
		if(variable.min() > up):
		  	up = None
		 
		if(variable.max() < low):
		  	low = None
		
		h = mpl.colors.Normalize(vmin=low,vmax=up,clip=False)
		
		if (direction == "X"):
		  extent = [ymin-ymax/2.0,ymax-ymax/2.0,zmin,zmax]
		  imshow(np.rot90(self.cube(variable)[point,:,:],1),extent=extent, interpolation = interp, norm=h,cmap=cmap, aspect='equal')
		elif (direction == "Y"):
		  extent = [xmin-xmax/2.0,xmax-xmax/2.0,zmin,zmax]
		  imshow(np.rot90(self.cube(variable)[:,point,:],1),extent=extent, interpolation = interp, norm=h,cmap=cmap, aspect='equal')
		elif (direction == "Z"):
		  extent = [xmin-xmax/2.0,xmax-xmax/2.0,ymin-ymax/2.0,ymax-ymax/2.0]
		  imshow((self.cube(variable)[:,:,point])[::-1,],extent=extent, interpolation = interp, norm=h,cmap=cmap, aspect='equal')
		elif (direction == "2Dside"):
		  extent = [ymin-ymax/2.0,ymax-ymax/2.0,zmin,zmax]
		  imshow(np.rot90(variable),extent=extent, interpolation = interp, norm=h,cmap=cmap, aspect='equal')
		elif (direction == "2Dtop"):
		  extent = [xmin-xmax/2.0,xmax-xmax/2.0,ymin-ymax/2.0,ymax-ymax/2.0]
		  imshow(variable,extent=extent, interpolation = interp, norm=h,cmap=cmap, aspect='equal')
		
		
		if(colbar):
		  acolbay = colorbar()
		  #plt.colorbar(heatmap)
		  #acolbay.set_label(r'$log_{10}(\rho [g/cm^{3}])$',size=20,color="black", rotation=270)
          #plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    
		  #colorbar.ColorbarBase.set_label(r'$log_{10}(\rho [g/cm^{3}])$')
		  #colorbar.set_label(r'$log_{10}(\rho [g/cm^{3}])$')
		
		return extent
		  
	def WriteVTKvector(self,filename,nx,ny,nz,xmin,ymin,zmin,dx,bx,by,bz):
         """Writes a VTK files given by the full path "filename" for the field bx,by,bz. NOTE: It writes first a dummy scalar, otherwise Paraview will not read it."""
#
         #Write VTK file
         #For some odd reason we need to write a dummy scalar otherwise with
         #just a vector it does not work
         f = open(filename,"wt")
         f.write('# vtk DataFile Version 2.0\n')
         f.write('Structured Grid Dataset\n')
         f.write('ASCII\n')
         f.write('DATASET STRUCTURED_POINTS\n')
         f.write('DIMENSIONS %d %d %d\n' % (nx, ny, nz))
         f.write('ORIGIN %g %g %g\n' %(xmin,ymin,zmin))
         f.write('SPACING %g %g %g\n' % (dx, dx, dx))
         f.write('POINT_DATA %d\n' % ((nx)*(ny)*(nz)))
#
         f.write('SCALARS Density float\n')
         f.write('LOOKUP_TABLE default\n')
#
#
         for k in range(0,nz):
             for j in range(0,ny):
                 for i in range(0,nx):
                     f.write("%g\n" % (1.0))
#
         f.write('VECTORS magnetic float\n')
         for k in range(0,nz):
            for j in range(0,ny):
                for i in range(0,nx):
                    f.write("%g %g %g\n" % (bx[i,j,k], by[i,j,k], bz[i,j,k]))
         f.close()
		
	def WriteVTKscalar(self,filename,nx,ny,nz,xmin,ymin,zmin,dx,var):
		"""Write scalar to a VTK file. The variable to be written is defined over a cube"""
		#Write VTK file
		f = open(filename,"wt")
		f.write('# vtk DataFile Version 2.0\n')
		f.write('Structured Grid Dataset\n')
		f.write('ASCII\n')
		f.write('DATASET STRUCTURED_POINTS\n')
		f.write('DIMENSIONS %d %d %d\n' % (nx, ny, nz))
		f.write('ORIGIN %g %g %g\n' %(xmin,ymin,zmin))
		f.write('SPACING %g %g %g\n' % (dx, dx, dx))
		f.write('POINT_DATA %d\n' % ((nx)*(ny)*(nz)))
		
		f.write('SCALARS var float\n')
		f.write('LOOKUP_TABLE default\n')
	
		for k in range(0,nz):
		  for j in range(0,ny):
		    for i in range(0,nx):
		      f.write("%g\n" % (var[i,j,k]))
		f.close()


	def ContourScalar(self,variable,variable2,xmin,ymin,point,scale,direction = "X", low = None, up = None):


          dx=self.CellSize()
          dx=dx*scale
          dy=dx
          xmin=xmin*scale
          ymin=ymin*scale
          variable_cube=self.cube(variable)
          if (direction == "X"):
              xlen=len(variable_cube[0,:,0])
              ylen=len(variable_cube[0,0,:])
          elif (direction == "Y"):
              xlen=len(variable_cube[:,0,0])
              ylen=len(variable_cube[0,0,:])
          elif (direction == "Z"):
              xlen=len(variable_cube[:,0,0])
              ylen=len(variable_cube[0,:,0])

          x = np.arange(xmin, xmin+xlen*dx, dx)
          y = np.arange(ymin, ymin+ylen*dy, dy)
          X, Y = np.meshgrid(x, y)
          X=np.rot90(X)
          Y=np.rot90(Y)
          
          if (direction == "X"):
              mpl.pyplot.figure(figsize=(6,10))              
          elif (direction == "Y"):
              mpl.pyplot.figure(figsize=(6,10))
          elif (direction == "Z"):
              mpl.pyplot.figure()
              
          self.slice(variable,point,direction,low=low,up=up,scale=scale,colmap=mpl.cm.YlOrRd)
          
          VarCont=self.cube(variable2)
          if (direction == "X"):
              VarCont=VarCont[point,:,:]
          elif (direction == "Y"):
              VarCont=VarCont[:,point,:]
          elif (direction == "Z"):
              VarCont=VarCont[:,:,point]

          CS = mpl.pyplot.contour(X, Y, VarCont, 6,colors='k')
          mpl.pyplot.clabel(CS, inline=1, fontsize=8)


	def PlotVectors(self, variable, variableX, variableY, xmin, ymin, point, scale, resolution, direction="Y", logscale=False, low=None, up=None):

          dx=self.CellSize()
          dx=dx*scale
          dy=dx
          xmin=xmin
          ymin=ymin
          variable_cube=self.cube(variable[:,0])   #Here to use the x-component of the vector field is not a problem since one just wants to know the box size
          if (direction == "X"):
              xlen=len(variable_cube[0,:,0])
              ylen=len(variable_cube[0,0,:])
          elif (direction == "Y"):
              xlen=len(variable_cube[:,0,0])
              ylen=len(variable_cube[0,0,:])
          elif (direction == "Z"):
              xlen=len(variable_cube[:,0,0])
              ylen=len(variable_cube[0,:,0])

          x = np.arange(xmin, xmin+xlen*dx, dx)
          y = np.arange(ymin, ymin+ylen*dy, dy)
          
          print(np.shape(x), np.shape(y), x[0], y[0], x[-1], y[-1])
          
          X, Y = np.meshgrid(x, y)
          X=np.rot90(X)
          Y=np.rot90(Y)
                        
          self.slice(magnitude(variable), point, direction, low=low, up=up, linlog=logscale, scale=scale, colmap=mpl.cm.YlOrRd)

          VarVectX=self.cube(variableX)
          VarVectY=self.cube(variableY)
          if (direction == "X"):
              VarVectX=VarVectX[point,:,:]
              VarVectY=VarVectY[point,:,:]
          elif (direction == "Y"):
              VarVectX=VarVectX[:,point,:]
              VarVectY=VarVectY[:,point,:]
          elif (direction == "Z"):
              VarVectX=VarVectX[:,:,point]
              VarVectY=VarVectY[:,:,point]
              
          QP = mpl.pyplot.quiver(X[::resolution,::resolution],Y[::resolution,::resolution],VarVectX[::resolution,::resolution],VarVectY[::resolution,::resolution],pivot='mid', color='g')
