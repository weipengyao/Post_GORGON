import vtk
import vtk.util.numpy_support as VN
import os
import glob
from numpy import zeros, sqrt, sort
from constants import *


def GetFileList(filepath,nstep=1):
	""" Get list of filenames from filepath, from first to last with step nstep.
		The default for nstep = 1
		Returns the filename list
	"""
	filename = sort(glob.glob(filepath + "*.vtk"))
	filename = filename[range(0,len(filename),nstep)]
	return filename

def GetAllDataFromRamses(filename):
	""" Read filename (unstructured vtk Ramses file) and put data in numpy arrays.
		It returns the following arrays and objects:
		return (time, coord, r, p, V, B, vol, mass, mesh, mycc)
		WARNING: It converts the magnetic field to Gauss
	"""
	#Read data from Unstructured Grid
	reader = vtk.vtkUnstructuredGridReader()
	reader.SetFileName(filename)
	reader.Update()
	mesh = reader.GetOutput()
	reader.CloseVTKFile()

	# Get cell data from the mesh
	c = mesh.GetCellData()
	# Get field data --> TIME
	f = mesh.GetFieldData()

	#Get Cell Centres
	cc = vtk.vtkCellCenters()
	cc.SetInput(mesh)
	cc.Update()
	mycc = cc.GetOutput()
	
	#Get time
	time = VN.vtk_to_numpy(f.GetArray("TIME"))
	
	#Put data into numpy arrays and convert magnetic field to Gauss
	r = VN.vtk_to_numpy(c.GetArray("density"))
	p = VN.vtk_to_numpy(c.GetArray("pressure"))
	V = VN.vtk_to_numpy(c.GetArray("velocity"))
	B = VN.vtk_to_numpy(c.GetArray("magnetic"))*1e-6
	
	#Get Coordinates for each cell
	coord = VN.vtk_to_numpy(mycc.GetPoints().GetData())

	#Calculate volume and mass for each cell
	#nc is the total number of cells
	nc = mesh.GetNumberOfCells()
	vol = zeros(nc)
	for i in range(nc):
		vol[i] = mesh.GetCell(i).GetLength2()
	
	vol = sqrt(vol/3.0)**3
	
	mass = (r*mp)*(vol*pc**3)/solar #mass in solar units

	return (time, coord, r, p, V, B, vol, mass, mesh, mycc)
	

