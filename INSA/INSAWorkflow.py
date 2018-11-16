import sys
sys.path.extend(['/home/nitram/Documents/work/MUPIF/commit'])
from mupif import *
import Pyro4
import argparse
#import win32com.client as win32

from INSAServerConfig import serverConfig
mode = 2
Cfg = serverConfig(mode)
import logging
log = logging.getLogger()
import time as timeT
import mupif.Physics.PhysicalQuantities as PQ
import copy

class INSAWorkflow(Workflow.Workflow):
        # # # # # # # # # # # # # # # # # # # # # # # # #
	def __init__ (self, targetTime=PQ.PhysicalQuantity('1 s')):
		"""
		Initializes the workflow. As the workflow is non-stationary, we allocate individual 
		applications and store them within a class.
		"""
		super(INSAWorkflow, self).__init__(file='', workdir='', targetTime=targetTime)
                #list of recognized input porperty IDs
		self.myInputPropIDs = [PropertyID.PID_ExtensionalInPlaneStiffness,  PropertyID.PID_ExtensionalOutOfPlaneStiffness, PropertyID.PID_ShearInPlaneStiffness, PropertyID.PID_ShearOutOfPlaneStiffness,PropertyID.PID_LocalBendingStiffness]     
                # list of compulsory IDs
		self.myCompulsoryPropIDs = self.myInputPropIDs
                #list of recognized output property IDs
		self.myOutPropIDs =  [PropertyID.PID_CriticalLoadLevel]
                #dictionary of input properties (values)
		self.myInputProps = {}
                #dictionary of output properties (values)
		self.myOutProps = {}
		self.mesh = None
		self.displacement = None
		self.fibreOrientation = None
		self.domainNumber = None
		#locate nameserver
		ns = PyroUtil.connectNameServer(nshost=Cfg.nshost, nsport=Cfg.nsport, hkey=Cfg.hkey)
		self.JobMan = PyroUtil.connectJobManager(ns, Cfg.jobManName, Cfg.hkey)
		self.INSASolver = None
		#allocate the application instances
		try:
			self.INSASolver = PyroUtil.allocateApplicationWithJobManager( ns, self.JobMan, None, Cfg.hkey, sshContext=None)
			log.info('Created INSA job')
		except Exception as e:
			log.exception(e)
		else:
			if (self.INSASolver is not None):
				INSASolverSignature=self.INSASolver.getApplicationSignature()
				log.info("Working INSA solver on server " + INSASolverSignature)
			else:
				log.debug("Connection to server failed, exiting")

	# # # # # # # # # # # # # # # # # # # # # # # # #
	def setProperty(self, property, objectID=0):
		"""
		Register given property in the application
		"""
		propID = property.getPropertyID()
		if (propID in self.myInputPropIDs):
                        self.myInputProps[propID]=property
		else:
                        raise APIError.APIError('Unknown property ID')
		#self.INSASolver.getProperty(property,val, objectID=0)

	# # # # # # # # # # # # # # # # # # # # # # # # #
	def getProperty(self, propID, objectID=0):
                if (propID in self.myOutPropIDs):
                        return self.myOutProps[propID]
                else:
                        raise APIError.APIError ('Unknown property ID', propID) 
		#return self.INSASolver.getProperty(property, objectID=0)

	# # # # # # # # # # # # # # # # # # # # # # # # #
	def solveStep(self, istep, stageID=0, runInBackground=False):
		"""
		Solves the problem for all the imposed displacements (no time involved).
		"""
		for cID in self.myCompulsoryPropIDs:
                        if cID not in self.myInputProps:
                                raise APIError.APIError (self.getApplicationSignature(), ' Missing compulsory property ', cID)
		try:
                        self.INSASolver.setProperty(self.myInputProps[PropertyID.PID_ExtensionalInPlaneStiffness])
                        self.INSASolver.setProperty(self.myInputProps[PropertyID.PID_ExtensionalOutOfPlaneStiffness])
                        self.INSASolver.setProperty(self.myInputProps[PropertyID.PID_ShearInPlaneStiffness])
                        self.INSASolver.setProperty(self.myInputProps[PropertyID.PID_ShearOutOfPlaneStiffness])
                        self.INSASolver.setProperty(self.myInputProps[PropertyID.PID_LocalBendingStiffness])
		except Exception as err:
			print ("Setting INSA params failed: " + repr(err));
			self.terminate()
		try:
                        # solve INSA part
			log.info("Running INSA solver")
			self.INSASolver.solveStep(None)
			self.mesh = self.INSASolver.getMesh(istep)
			self.displacement = self.INSASolver.getField(FieldID.FID_Displacement,0,0)
			self.fibreOrientation0 = self.INSASolver.getField(FieldID.FID_FibreOrientation,0,1)
			self.fibreOrientation90 = self.INSASolver.getField(FieldID.FID_FibreOrientation,0,2)
			self.fibreOrientation45 = self.INSASolver.getField(FieldID.FID_FibreOrientation,0,3)
			self.fibreOrientation_45 = self.INSASolver.getField(FieldID.FID_FibreOrientation,0,4)
			self.domainNumber = self.INSASolver.getField(FieldID.FID_DomainNumber,0,0)

		except Exception as err:
                        print ("Error:" + repr(err))
                        self.terminate()        


	def terminate(self):
		#self.thermalAppRec.terminateAll()
		self.INSASolver.terminate()
		super(INSAWorkflow, self).terminate()
	# # # # # # # # # # # # # # # # # # # # # # # # # 
	def getMesh (self):
		"""
		Returns the computational mesh.
		:return: Returns the representation of mesh
		:rtype: Mesh
		"""
		#return self.INSASolver.getMesh()
		return self.mesh

	# # # # # # # # # # # # # # # # # # # # # # # # # 
	def getField(self, fieldID,time, objectID) :
		"""
		Returns the requested field at given time. Field is	identified by fieldID.
		"""
		if fieldID == FieldID.FID_Displacement:
			return self.displacement
		elif fieldID == FieldID.FID_FibreOrientation:
			if objectID == 0:
				return self.fibreOrientation0
			elif objectID == 90:
				return self.fibreOrientation90
			elif objectID == 45:
				return self.fibreOrientation45
			elif objectID == -45:
				return self.fibreOrientation_45
			else:
				print("Unknown layer orientation:" +objectID)
		elif fieldID == FieldID.FID_DomainNumber:
			return self.domainNumber
		
		#return self.INSASolver.getField(fieldID,time)

	# # # # # # # # # # # # # # # # # # # # # # # # # 
	def getAPIVersion(self):
		"""
		:return: Returns the supported API version
		:rtype: str, int
		"""
		return "1.0"


	# # # # # # # # # # # # # # # # # # # # # # # # # 
	def getApplicationSignature(self):
		"""
		Get application signature.

		:return: Returns the application identification
		:rtype: str
		"""
		return "INSA workflow 1.0"
	#return self.INSASolver.getApplicationSignature()


        # # # # # # # # # # # # # # # # # # # # # # # # # 
	def getCriticalTimeStep(self):
		"""
		Returns true or false depending whether solve has completed.
		"""
		return PQ.PhysicalQuantity(1.0,'s')


def filterField (sourceField, filterField, tresholdValue):
	# filter composite part
	# sourceField: source field from which subFIELD IS CREATED
	# filterField: Field defining quantintity based on which the filtering is done
	# tresholdValue: Only parts of source field with filterField = tresholdValue will be returned
	# return: field defined on submesh of sourceField, where filterField==tresholdValue
	fmesh = sourceField.getMesh()
	targetMesh = Mesh.UnstructuredMesh()
	targetCells = []
	targetVertices = []
	nodeMap = [-1]*fmesh.getNumberOfVertices()

	fvalues = []
	for iv in fmesh.vertices():
		n = iv.getNumber()
		if (abs(filterField.getVertexValue(n).getValue()[0] - tresholdValue) <1.e-3):
			nn = copy.deepcopy(iv)
			nn.number = len(targetVertices)
			nodeMap[n]=nn.number
			targetVertices.append(nn)
			fvalues.append(sourceField.getVertexValue(n).getValue())
		
	for icell in fmesh.cells():
		# determine if icell belongs to composite domain
		if (nodeMap[icell.getVertices()[0].getNumber()]>=0):
			# cell belonging to composite
			c = icell.copy()
			# append all cell vertices
			cvertices = []
			for i in range(len(c.vertices)):
				#inum = c.vertices[i].getNumber()
				inum = c.vertices[i]
				cvertices.append(nodeMap[inum])
			c.vertices = cvertices
			targetCells.append(c)

	targetMesh.setup(targetVertices, targetCells)
	targetField = Field.Field(targetMesh, FieldID.FID_FibreOrientation, sourceField.getValueType(), sourceField.getUnits(), sourceField.getTime(), values=fvalues, fieldType=sourceField.getFieldType(), objectID=sourceField.getObjectID())
	return targetField


if __name__ == "__main__":
	workflow = INSAWorkflow(targetTime=PQ.PhysicalQuantity(1.0,'s'))
	workflow.setProperty(Property.ConstantProperty(100, PropertyID.PID_ExtensionalInPlaneStiffness, ValueType.Scalar, 'MPa'))
	workflow.setProperty(Property.ConstantProperty(100, PropertyID.PID_ExtensionalOutOfPlaneStiffness, ValueType.Scalar, 'MPa'))
	workflow.setProperty(Property.ConstantProperty(100, PropertyID.PID_ShearInPlaneStiffness, ValueType.Scalar, 'MPa'))
	workflow.setProperty(Property.ConstantProperty(100, PropertyID.PID_ShearOutOfPlaneStiffness, ValueType.Scalar, PQ.getDimensionlessUnit()))
	workflow.setProperty(Property.ConstantProperty(100, PropertyID.PID_LocalBendingStiffness, ValueType.Scalar, PQ.getDimensionlessUnit()))
	workflow.solve()
	workflow.getMesh()
	displ = workflow.getField(FieldID.FID_Displacement,0,0)
	dn = workflow.getField(FieldID.FID_DomainNumber,0,0)
	dn.fieldType == Field.FieldType.FT_vertexBased
	fibre0 = workflow.getField(FieldID.FID_FibreOrientation,0,0)
	fibre90 = workflow.getField(FieldID.FID_FibreOrientation,0,90)
	fibre45 = workflow.getField(FieldID.FID_FibreOrientation,0,45)
	fibre_45 = workflow.getField(FieldID.FID_FibreOrientation,0,-45)
	displ.field2VTKData(displ.fieldID.name).tofile('Airbus_'+displ.fieldID.name)


	filDispl = filterField(displ, dn, 1.0)
	filDispl.field2VTKData(filDispl.fieldID.name).tofile('Airbus_'+filDispl.fieldID.name)

	filFibre0 = filterField(fibre0, dn, 1.0)
	filFibre0.field2VTKData(filFibre0.fieldID.name).tofile('Airbus_'+filFibre0.fieldID.name+'_0')

	filFibre90 = filterField(fibre90, dn, 1.0)
	filFibre90.field2VTKData(filFibre90.fieldID.name).tofile('Airbus_'+filFibre90.fieldID.name+'_90')

	filFibre45 = filterField(fibre45, dn, 1.0)
	filFibre45.field2VTKData(filFibre45.fieldID.name).tofile('Airbus_'+filFibre45.fieldID.name+'_45')

	filFibre_45 = filterField(fibre_45, dn, 1.0)
	filFibre_45.field2VTKData(filFibre_45.fieldID.name).tofile('Airbus_'+filFibre_45.fieldID.name+'_-45')

