#!/usr/bin/env python
# Version 2.0.1 Created by Z.Cao 8/Feb/2019
import os, sys, getopt, math
import numpy as np
from decimal import Decimal
import pandas as pd
newchgcar = []
Mode = ''
atoms = []

def readfiles(name):
	try:
		files = open(str(os.getcwd()) + '/' + name)
	except:
		print ("Cannot Open " + name)
		sys.exit(1)
	return files


def getfootercontent(skiprows, readrows, times):
	CHGCAR_path = str(os.getcwd()) + '/CHGCAR'
	foot = ''
	count = 1
	footerlist = open(CHGCAR_path).readlines()[skiprows+readrows : len(open(CHGCAR_path).readlines())]
	temp1 = []
	temp2 = []
	atomtemp = []
	footerindex = []
	for i, val in enumerate(footerlist):
		if val.find('augmentation') == 0:
				footerindex.append(i)
	footerindex.append(len(footerlist))

	for i in range(int(len(footerindex)-1)):
		temp1.append (''.join(footerlist[footerindex[i]+1:footerindex[i+1]]))
		temp2.append (footerlist[footerindex[i]].split()[3])

	# for i in atoms:
	# 	atomtemp.append(int(int(i)/times - 0.0001) + 1)
	compa = []
	for k in range(int(atoms[-1])):
		for i in range(times):
			compa.append(k+1)
	for i in atoms:
		atomtemp.append(compa[int(i)-1])

	# DEBUG
	# print ('-------------------------------')
	# print (atoms)
	# print (atomtemp)

	for i in atomtemp:
		foot += '{0} {1:>3} {2:>3}'.format('augmentation occupancies',count, int(temp2[i-1]))+'\n'
		foot += temp1[i-1]
		count += 1
	return foot


def finaloutput(Nx, Ny, Nz, S1, S2, S3, skiprows, readrows):
	print ('Finish Calculating, Writing CHGCAR.new ...')
	fout = open(str(os.getcwd()) + '/CHGCAR.new','w')
	header = open(str(os.getcwd()) + '/POSCAR.new','r')
	content = open(str(os.getcwd()) + '/CHGCAR.new.temp','r')
	fout.writelines (header)
	fout.writelines ('\n  '+str(S1*Nx)+'  '+str(S2*Ny)+'  '+str(S3*Nz)+'\n')
	for i in content:
		fout.writelines (cut(i))
	fout.writelines (getfootercontent(skiprows, readrows, S1*S2*S3))
	fout.close()
	os.remove(str(os.getcwd())+'/CHGCAR.new.temp')

def outputarray (array, Nx, Ny, Nz, S1, S2, S3):
	print ('Calculating ...')
	fout = open(str(os.getcwd()) + '/CHGCAR.new.temp','w')
	header = open(str(os.getcwd()) + '/POSCAR.new','r')
	outputpotential = ''
	flag = 1
	for i in range(len(array)):
		if flag <5:
			outputpotential += ' {:.11E}'.format(array[i])
			flag += 1
		else:
			outputpotential += ' {:.11E}\n'.format(array[i])
			flag = 1
	# fout.writelines (header)
	# fout.writelines ('\n  '+str(S1*Nx)+'  '+str(S2*Ny)+'  '+str(S3*Nz)+'\n')
	fout.writelines (outputpotential)
	fout.close()


def extendarray (array, Nx, Ny, Nz, S1, S2, S3):
	# print (str(S1)+str(S2)+str(S3))
	newarray = np.zeros(S1*S2*S3*Nx*Ny*Nz)
	print ('Extending Potential ...')
	Narray = np.reshape(array,(Nz,Ny,Nx))
	Nnewarray = np.reshape(newarray, (Nz*S3,Ny*S2,Nx*S1))
	for k in range(S3*Nz):
		for j in range (S2*Ny):
			for i in range (S1*Nx):
				Nnewarray[k,j,i] = Narray [(k+Nz)%Nz , (j+Ny)%Ny , (i+Nx)%Nx]
	newarray = Nnewarray.reshape(S1*S2*S3*Nx*Ny*Nz)
	return newarray


def BaderMergeMode(para):
	flag = 0
	atomtemp = []
	print ('Reading Parameter ...')
	Upotential, Nx, Ny, Nz, lattice, skiprows, readrows = read_vasp_density(str(os.getcwd()) + '/CHGCAR')
	times = para[0]*para[1]*para[2]
	# for i in atoms:
	# 	atomtemp.append(int(int(i)/times - 0.0001) + 1)
	compa = []
	for k in range(int(atoms[-1])):
		for i in range(times):
			compa.append(k+1)
	for i in atoms:
		atomtemp.append(compa[int(i)-1])

	for i in atomtemp:
		print ('Processing Atom #' + str(i))
		if flag == 0:
			potential = only_read_vasp_density(str(os.getcwd()) + '/BvAt' + str(i).zfill(4) + '.dat')
			flag = 1
		else:
			potential += only_read_vasp_density(str(os.getcwd()) + '/BvAt' + str(i).zfill(4) + '.dat')

	newpotential = extendarray (potential, Nx, Ny, Nz, para[0], para[1], para[2])
	selectposcar()
	outputarray (newpotential, Nx, Ny, Nz, para[0], para[1], para[2])
	finaloutput(Nx, Ny, Nz, para[0], para[1], para[2], skiprows, readrows)

def CHGCARextendMode(para):
	potential , Nx, Ny, Nz, lattice, skiprows, readrows = read_vasp_density(str(os.getcwd()) + '/CHGCAR')
	newpotential = extendarray (potential, Nx, Ny, Nz, para[0], para[1], para[2])
	fullselectposcar()
	outputarray (newpotential, Nx, Ny, Nz, para[0], para[1], para[2])
	finaloutput(Nx, Ny, Nz, para[0], para[1], para[2], skiprows, readrows)

def getatoms(axis, minnum, maxnum):
	atomcounts = [int(i) for i in readfiles('POSCAR.new.temp').readlines()[6].split()]
	rows = sum (atomcounts)
	matrixs = readfiles('POSCAR.new.temp').readlines()[8:8+rows]
	if axis == 'x' or axis == 'X':
		for i in range(rows):
			if Decimal(matrixs[i].split()[0]) >= minnum and Decimal(matrixs[i].split()[0]) < maxnum:
				atoms.append(str(i+1).zfill(4))
	elif axis == 'y' or axis == 'Y':
		for i in range(rows):
			if Decimal(matrixs[i].split()[1]) >= minnum and Decimal(matrixs[i].split()[1]) < maxnum:
				atoms.append(str(i+1).zfill(4))
	elif axis == 'z' or axis == 'Z':
		for i in range(rows):
			if Decimal(matrixs[i].split()[2]) >= minnum and Decimal(matrixs[i].split()[2]) < maxnum:
				atoms.append(str(i+1).zfill(4))
	

def getparameters():
	parameter = [int(i) for i in readfiles('parameter.txt').readlines()[0].split()]
	getextendposcar(parameter)

	Mode = readfiles('parameter.txt').readlines()[1].split()[0]
	print (Mode + ' Mode, Selected Atoms:')
	if Mode == 'Range' or Mode == 'RANGE' or Mode == 'range':
		getatoms(str(readfiles('parameter.txt').readlines()[2].split()[0]), Decimal(readfiles('parameter.txt').readlines()[3].split()[0]), Decimal(readfiles('parameter.txt').readlines()[4].split()[0]))
		print (atoms)
		BaderMergeMode(parameter)
	elif Mode == 'All' or Mode == 'all' or Mode == 'ALL':
		print (atoms)
		getatoms('z',0.0,1.5)
		CHGCARextendMode(parameter)
	elif Mode == 'list' or Mode == 'LIST' or Mode == 'List':
		for i in range(len(readfiles('parameter.txt').readlines())-2):
			atoms.append(readfiles('parameter.txt').readlines()[i+2].split()[0].zfill(4))
		print (atoms)
		BaderMergeMode(parameter)


	return parameter

def only_read_vasp_density(FILE, use_pandas=None):
	if use_pandas:
		from pandas import read_table as pandas_read_table
	elif use_pandas is None:
		try:
			from pandas import read_table as pandas_read_table
			use_pandas = True
		except ImportError:
			use_pandas = False

	# print("Reading File Header")
	with open(FILE, "r") as f:
		_ = f.readline()
		scale_factor = float(f.readline())

		lattice = np.zeros(shape=(3,3))
		for row in range(3):
			lattice[row] = [float(x) for x in f.readline().split()]
			lattice = lattice * scale_factor

		num_species = len(f.readline().split())
		num_type = [int(x) for x in f.readline().split()]
		num_atoms = sum(num_type)
		coord_type = f.readline().strip()

		coordinates = np.zeros(shape=(num_atoms, 3))
		for atom_i in range(num_atoms):
			coordinates[atom_i] = [float(x) for x in f.readline().split()]

		# Skip blank line
		_ = f.readline()

		NGX, NGY, NGZ = [int(x) for x in f.readline().split()]

		if use_pandas:
			# print("Reading Potential with pandas")
			skiprows = 10 + num_atoms
			readrows = int(math.ceil(NGX * NGY * NGZ / 5))

			dat = pandas_read_table(FILE, delim_whitespace=True, skiprows=skiprows, header=None, nrows=readrows)
			Potential = dat.iloc[:readrows, :5].values.flatten()
			remainder = (NGX * NGY * NGZ) % 5
			if remainder > 0:
				Potential = Potential[:(-5 + remainder)]

		else:
			# print("Reading CHGCAR without pandas")
			Potential = (f.readline().split() for i in range(int(math.ceil(NGX * NGY * NGZ / 5))))
			Potential = np.fromiter(chain.from_iterable(Potential), float)

	# print("Average of the potential = ", np.average(Potential))

	return Potential

def read_vasp_density(FILE, use_pandas=None):
	"""Generic reading of CHGCAR LOCPOT etc files from VASP
    Args:
        FILE (str): Path to density file
        use_pandas (bool): Use Pandas library for faster file reading. If set
            to None, Pandas will be used when available.
    Returns:
        Potential (array), NGX (int), NGY (int), NGZ (int), lattice (array)
        where Potential is a 1-D flattened array of density data with original
        dimensions NGX x NGY x NGZ and lattice is the 3x3 unit-cell matrix.
    @copyright https://github.com/WMD-group/MacroDensity/blob/master/macrodensity/density_tools.py
	"""
# Get Header information by reading a line at a time

	if use_pandas:
		from pandas import read_table as pandas_read_table
	elif use_pandas is None:
		try:
			from pandas import read_table as pandas_read_table
			use_pandas = True
		except ImportError:
			use_pandas = False

	# print("Reading File Header")
	with open(FILE, "r") as f:
		_ = f.readline()
		scale_factor = float(f.readline())

		lattice = np.zeros(shape=(3,3))
		for row in range(3):
			lattice[row] = [float(x) for x in f.readline().split()]
			lattice = lattice * scale_factor

		num_species = len(f.readline().split())
		num_type = [int(x) for x in f.readline().split()]
		num_atoms = sum(num_type)
		coord_type = f.readline().strip()

		coordinates = np.zeros(shape=(num_atoms, 3))
		for atom_i in range(num_atoms):
			coordinates[atom_i] = [float(x) for x in f.readline().split()]

		# Skip blank line
		_ = f.readline()

		NGX, NGY, NGZ = [int(x) for x in f.readline().split()]

		if use_pandas:
			print("Reading Potential with pandas")
			skiprows = 10 + num_atoms
			readrows = int(math.ceil(NGX * NGY * NGZ / 5))

			dat = pandas_read_table(FILE, delim_whitespace=True, skiprows=skiprows, header=None, nrows=readrows)
			Potential = dat.iloc[:readrows, :5].values.flatten()
			remainder = (NGX * NGY * NGZ) % 5
			if remainder > 0:
				Potential = Potential[:(-5 + remainder)]

		else:
			# print("Reading CHGCAR without pandas")
			Potential = (f.readline().split() for i in range(int(math.ceil(NGX * NGY * NGZ / 5))))
			Potential = np.fromiter(chain.from_iterable(Potential), float)

	# print("Average of the potential = ", np.average(Potential))

	return Potential, NGX, NGY, NGZ, lattice, skiprows, readrows


def getextendposcar(parameter, outname = 'POSCAR.new.temp'):
	print ('Writing POSCAR.new ...')
	outputpath = str(os.getcwd())+"/"+str(outname)
	Na = int(parameter[0])
	Nb = int(parameter[1])
	Nc = int(parameter[2])
	index = []
	x = []
	y = []
	z = []
	gridx = []
	gridy = []
	gridz = []
	poscar = open(str(os.getcwd())+"/CHGCAR").readlines()
	rows = poscar[6].split()

	atoms = len(rows)
	for i in range(atoms):
		index.append(int(rows[i]))
	totoalatoms =  sum(index)
	newindex = [str(i*Na*Nb*Nc) for i in index]
	line7 = "     " + "     ".join(newindex)+'\n'
	for i in range(totoalatoms):
		x.append(Decimal(poscar[i+8].split()[0]))
		y.append(Decimal(poscar[i+8].split()[1]))
		z.append(Decimal(poscar[i+8].split()[2]))
	newx = []
	newy = []
	newz = []
	for i in range(len(x)):
		x[i] = Decimal(x[i])/Decimal(Na)
		y[i] = Decimal(y[i])/Decimal(Nb)
		z[i] = Decimal(z[i])/Decimal(Nc)
		for ia in range(Na):
			for ib in range(Nb):
				for ic in range(Nc):
					# print ("X:"+str(x[i])+"Y:"+str(y[i])+"Z:"+str(z[i])+'\n')
					newx.append(x[i] + Decimal(ia/Na))
					newy.append(y[i] + Decimal(ib/Nb))
					newz.append(z[i] + Decimal(ic/Nc))
					# print ("NX:"+str(newx[i])+"NY: "+str(newy[i])+"NZ:"+str(newz[i])+'\n')
	print (newindex)
	
	fout = open(outputpath,'w')
	for i in range (2):
		fout.writelines(poscar[i])
	for i in range (3):
		gridx.append(poscar[2].split()[i])
		gridy.append(poscar[3].split()[i])
		gridz.append(poscar[4].split()[i])
	line3 = ' '
	for strn in [Decimal(gridx[0])*Na,Decimal(gridx[1])*Na,Decimal(gridx[2])*Na]:
		line3 += '{0:>22.16f}'.format(strn)
	line3 += '\n '
	for strn in [Decimal(gridy[0])*Nb,Decimal(gridy[1])*Nb,Decimal(gridy[2])*Nb]:
		line3 += '{0:>22.16f}'.format(strn)
	line3 += '\n '
	for strn in [Decimal(gridz[0])*Nc,Decimal(gridz[1])*Nc,Decimal(gridz[2])*Nc]:
		line3 += '{0:>22.16f}'.format(strn)
	line3 += '\n '
	fout.writelines(line3[:-1])
	fout.writelines(poscar[5])
	fout.writelines(line7)
	fout.writelines('Direct\n')
	for i in range(len(newx)):
		fout.writelines(' '+"{:.16f}".format(newx[i])+' '+"{:.16f}".format(newy[i])+' '+"{:.16f}".format(newz[i])+'\n')
	# fout.writelines('     '.join(str(int(rows)*Na*Nb*Nc)))
	fout.close()

def selectposcar():
	poscar = open(str(os.getcwd())+"/POSCAR.new.temp")
	fout = open(str(os.getcwd())+"/POSCAR.new", 'w')
	for i in range(6):
		fout.writelines (readfiles('POSCAR.new.temp').readlines()[i])
	sumatom = sum([int(i) for i in readfiles('POSCAR.new.temp').readlines()[6].split()])
	factor = int(sumatom/len(atoms))

	line7 = [str(int(int(i)/factor+0.5)) for i in readfiles('POSCAR.new.temp').readlines()[6].split()]

	fout.writelines('  ' + '  '.join(line7)+'\n')
	fout.writelines('Direct\n')
	for i in range(len(atoms)):
		fout.writelines(readfiles('POSCAR.new.temp').readlines()[int(atoms[i]) + 7])
	fout.close()
	os.remove(str(os.getcwd())+'/POSCAR.new.temp')

def fullselectposcar():
	poscar = open(str(os.getcwd())+"/POSCAR.new.temp")
	fout = open(str(os.getcwd())+"/POSCAR.new", 'w')
	for i in readfiles('POSCAR.new.temp').readlines():
		fout.writelines(i)
	fout.close()
	os.remove(str(os.getcwd())+'/POSCAR.new.temp')



def cut(line):
	newline = []
	for i in line.split():
		if i[0] == '-':
			newstr = i[0]+i[2]+i[1]
			if i[-3:] == '-01':
				newstr = i[0]+i[2]+i[1]+i[3:13]+'E'+'+00'
			elif i[-3] == '+':
				newstr = i[0]+i[2]+i[1]+i[3:13]+'E+'+str(int(i[-2:])+1).zfill(2)
			elif i[-3] == '-':
				newstr = i[0]+i[2]+i[1]+i[3:13]+'E-'+str(int(i[-2:])-1).zfill(2)
			newline.append(newstr)
		else:
			newline.append(i)
	final = ' ' + ' '.join(newline)+'\n'
	return final

def main(argv):
	outputfile_name = ''

	try:
		opts, args = getopt.getopt(argv, "ci:o:",["inputfile=","outputfile="])

	except getopt.GetoptError:
		print ("Error parameter: BaderMerge.py -i <AtomSelection> -o <Outputfile_name>")
		sys.exit(2)

	for opt, arg in opts:
		if opt == '-c':
			print ('Copyright CaoZheng @ 2018')
			sys.exit()

		elif opt in ("-i","--inputfile"): #Atom-selection number
			# Atom_selection.append(arg)
			pass

		elif opt in ("-o","--outputfile"):
			outputfile_name = arg

	# Atom_selection = makeuplist()

	#Atom_selection should be list ['0001','0002',etc.]
	# if Atom_selection == []:
	# 	print ("Please add -i atom-selection!")
	# else:
	# 	print ("Atom selection:" + str(Atom_selection))
	# 	if outputfile_name != '':
	# 		check(Atom_selection)
	# 		mergy(Atom_selection,outputfile_name)
	# 	else:
	# 		check(Atom_selection)
	# 		mergy(Atom_selection)
	para = getparameters()


if __name__ == "__main__":
	main(sys.argv[1:])