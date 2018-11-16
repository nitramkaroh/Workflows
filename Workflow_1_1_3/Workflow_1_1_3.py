import sys
sys.path.extend(['/home/nitram/Documents/work/MUPIF/commit'])
from mupif import *
import Pyro4
import logging
log = logging.getLogger()
import time as timeT
import mupif.Physics.PhysicalQuantities as PQ
import LAMMPS_v3 as lammps

nshost = '172.30.0.1'
nsport = 9090
hkey = 'mupif-secret-key'
digimatJobManName='eX_DigimatMF_JobManager'
mul2JobManName='MUL2.JobManager@Demo1'

class Workflow_1_1_3(Workflow.Workflow):
   
    def __init__(self, targetTime=PQ.PhysicalQuantity('0 s')):
        """
        Initializes the workflow. As the workflow is non-stationary, we allocate individual 
        applications and store them within a class.
        """
        super(Workflow_1_1_3, self).__init__(file='', workdir='', targetTime=targetTime)

        #list of recognized input porperty IDs
        self.myInputPropIDs = [PropertyID.PID_SMILE_MOLECULAR_STRUCTURE,PropertyID.PID_MOLECULAR_WEIGHT, PropertyID.PID_CROSSLINKER_TYPE,PropertyID.PID_FILLER_DESIGNATION, PropertyID.PID_CROSSLINKONG_DENSITY,PropertyID.PID_FILLER_CONCENTRATION, PropertyID.PID_TEMPERATURE, PropertyID.PID_PRESSURE, PropertyID.PID_POLYDISPERSITY_INDEX,PropertyID.PID_SMILE_MODIFIER_MOLECULAR_STRUCTURE,PropertyID.PID_SMILE_FILLER_MOLECULAR_STRUCTURE,PropertyID.PID_DENSITY_OF_FUNCTIONALIZATION, PropertyID.PID_InclusionYoung, PropertyID.PID_InclusionPoisson, PropertyID.PID_InclusionVolumeFraction, PropertyID.PID_InclusionAspectRatio]
        # list of compulsory IDs
        self.myCompulsoryPropIDs = self.myInputPropIDs

        #list of recognized output property IDs
        self.myOutPropIDs =  [PropertyID.PID_CriticalLoadLevel]

        #dictionary of input properties (values)
        self.myInputProps = {}
        #dictionary of output properties (values)
        self.myOutProps = {}

        #locate nameserver
        ns = PyroUtil.connectNameServer(nshost, nsport, hkey)
        #connect to digimat JobManager running on (remote) server
        self.digimatJobMan = PyroUtil.connectJobManager(ns, digimatJobManName,hkey)
        #connect to mult2 JobManager running on (remote) server
        self.mul2JobMan = PyroUtil.connectJobManager(ns, mul2JobManName,hkey)

        # solvers
        self.lammpsSolver = lammps.emailAPI(None) # local LAMMPS API instances
        self.digimatSolver = None
        self.mul2Solver = None
        #allocate the Digimat remote instance
        try:
            self.digimatSolver = PyroUtil.allocateApplicationWithJobManager( ns, self.digimatJobMan, None, hkey)
            log.info('Created digimat job')
            self.mul2Solver = PyroUtil.allocateApplicationWithJobManager( ns, self.mul2JobMan, None, hkey)
            log.info('Created mul2 job')            
        except Exception as e:
            log.exception(e)
        else:
            if ((self.lammpsSolver is not None) and (self.digimatSolver is not None) and (self.mul2Solver is not None)):
                lammpsSolverSignature=self.lammpsSolver.getApplicationSignature()
                log.info("Working lammps solver on server " + lammpsSolverSignature)
                digimatSolverSignature=self.digimatSolver.getApplicationSignature()
                log.info("Working digimat solver on server " + digimatSolverSignature)
                mul2SolverSignature=self.mul2Solver.getApplicationSignature()
                log.info("Working mul2 solver on server " + mul2SolverSignature)
            else:
                log.debug("Connection to server failed, exiting")


    def setProperty(self, property, objectID=0):
        propID = property.getPropertyID()
        if (propID in self.myInputPropIDs):
            self.myInputProps[propID]=property
        else:
            raise APIError.APIError('Unknown property ID')

    def getProperty(self, propID, time, objectID=0):
        if (propID in self.myOutPropIDs):
            return self.myOutProps[propID]
        else:
            raise APIError.APIError ('Unknown property ID', propID)   

    def solveStep(self, istep, stageID=0, runInBackground=False):

        for cID in self.myCompulsoryPropIDs:
            if cID not in self.myInputProps:
                raise APIError.APIError (self.getApplicationSignature(), ' Missing compulsory property ', cID)   
            
        try:
            # lammps 
            self.lammpsSolver.setProperty(self.myInputProps[PropertyID.PID_SMILE_MOLECULAR_STRUCTURE])
            self.lammpsSolver.setProperty(self.myInputProps[PropertyID.PID_MOLECULAR_WEIGHT])
            self.lammpsSolver.setProperty(self.myInputProps[PropertyID.PID_CROSSLINKER_TYPE])
            self.lammpsSolver.setProperty(self.myInputProps[PropertyID.PID_FILLER_DESIGNATION])
            self.lammpsSolver.setProperty(self.myInputProps[PropertyID.PID_CROSSLINKONG_DENSITY])
            self.lammpsSolver.setProperty(self.myInputProps[PropertyID.PID_FILLER_CONCENTRATION])
            self.lammpsSolver.setProperty(self.myInputProps[PropertyID.PID_TEMPERATURE])
            self.lammpsSolver.setProperty(self.myInputProps[PropertyID.PID_PRESSURE])
            self.lammpsSolver.setProperty(self.myInputProps[PropertyID.PID_POLYDISPERSITY_INDEX])
            self.lammpsSolver.setProperty(self.myInputProps[PropertyID.PID_SMILE_MODIFIER_MOLECULAR_STRUCTURE])
            self.lammpsSolver.setProperty(self.myInputProps[PropertyID.PID_SMILE_FILLER_MOLECULAR_STRUCTURE])
            self.lammpsSolver.setProperty(self.myInputProps[PropertyID.PID_DENSITY_OF_FUNCTIONALIZATION])
	    # solve (involves operator interaction)
            # first set useCaseID and execID using metadata
            self.lammpsSolver.setMetadata(MetadataKeys.UseCaseID, self.getMetadata(MetadataKeys.UseCaseID))
            self.lammpsSolver.setMetadata(MetadataKeys.ExecID, self.getMetadata(MetadataKeys.ExecID))
            # solve lammps, useCaseID and execID are not inputs!
            self.lammpsSolver.solveStep (TimeStep.TimeStep(0.0, 0.1, 1, 's'))
            # get result of the simulation
            matrixYoung = self.lammpsSolver.getProperty(PropertyID.PID_EModulus, 0.0)
            matrixPoisson = self.lammpsSolver.getProperty(PropertyID.PID_PoissonRatio, 0.0)

        except Exception as err:
            print ("Error:" + repr(err))
            self.terminate()


        # digimat
        try:
            # map properties from lammps to properties of Digimat
            matrixYoung.propID = PropertyID.PID_MatrixYoung
            matrixPoisson.propID = PropertyID.PID_MatrixPoisson
            self.digimatSolver.setProperty(matrixYoung)
            self.digimatSolver.setProperty(matrixPoisson)

            # fixed properties - taken form the database
            self.digimatSolver.setProperty(self.myInputProps[PropertyID.PID_InclusionYoung])
            self.digimatSolver.setProperty(self.myInputProps[PropertyID.PID_InclusionPoisson])
            self.digimatSolver.setProperty(self.myInputProps[PropertyID.PID_InclusionVolumeFraction])
            self.digimatSolver.setProperty(self.myInputProps[PropertyID.PID_InclusionAspectRatio])
        except Exception as err:
            print ("Setting Digimat params failed: " + repr(err));
            self.terminate()

        try:
            # solve digimat part
            log.info("Running digimat")
            self.digimatSolver.solveStep(None)
            ## get the desired properties
            compositeAxialYoung = self.digimatSolver.getProperty(PropertyID.PID_CompositeAxialYoung)
            compositeInPlaneYoung = self.digimatSolver.getProperty(PropertyID.PID_CompositeInPlaneYoung)
            compositeInPlaneShear = self.digimatSolver.getProperty(PropertyID.PID_CompositeInPlaneShear)
            compositeTransverseShear = self.digimatSolver.getProperty(PropertyID.PID_CompositeTransverseShear)
            compositeInPlanePoisson = self.digimatSolver.getProperty(PropertyID.PID_CompositeInPlanePoisson)
            compositeTransversePoisson = self.digimatSolver.getProperty(PropertyID.PID_CompositeTransversePoisson)
            
        except Exception as err:
            print ("Error:" + repr(err))
            self.terminate()


        try:
            # map properties from Digimat to properties of MUL2
            # Young modulus
            compositeAxialYoung.propID = PropertyID.PID_YoungModulus1
            compositeInPlaneYoung1 = compositeInPlaneYoung
            compositeInPlaneYoung1.propID = PropertyID.PID_YoungModulus2
            compositeInPlaneYoung2 = compositeInPlaneYoung
            compositeInPlaneYoung2.propID = PropertyID.PID_YoungModulus3
            # Shear modulus
            compositeInPlaneShear.propID = PropertyID.PID_ShearModulus12
            compositeTransverseShear1 = compositeTransverseShear
            compositeTransverseShear1.propID = PropertyID.PID_ShearModulus13
            compositeTransverseShear2 = compositeTransverseShear
            compositeTransverseShear2.propID = PropertyID.PID_ShearModulus23
            # Poisson ratio
            compositeInPlanePoisson.propID =  PropertyID.PID_PoissonRatio12
            compositeTransversePoisson1 =  compositeTransversePoisson
            compositeTransversePoisson1.propID = PropertyID.PID_PoissonRatio13
            compositeTransversePoisson2 =  compositeTransversePoisson
            compositeTransversePoisson2.propID = PropertyID.PID_PoissonRatio23
            
            self.mul2Solver.setProperty(compositeAxialYoung)
            self.mul2Solver.setProperty(compositeInPlaneYoung1)
            self.mul2Solver.setProperty(compositeInPlaneYoung2)
            
            self.mul2Solver.setProperty(compositeInPlaneShear)          
            self.mul2Solver.setProperty(compositeTransverseShear1)
            self.mul2Solver.setProperty(compositeTransverseShear2)
            
            self.mul2Solver.setProperty(compositeInPlanePoisson)          
            self.mul2Solver.setProperty(compositeTransversePoisson1)
            self.mul2Solver.setProperty(compositeTransversePoisson2)
            
            
        except Exception as err:
            print ("Setting MUL2 params failed: " + repr(err));
            self.terminate()
            
        try:
            # solve digimat part
            log.info("Running mul2")
            self.mul2Solver.solveStep(None)
            ## get the desired properties
            self.myOutProps[PropertyID.PID_CriticalLoadLevel] = self.mul2Solver.getProperty(PropertyID.PID_CriticalLoadLevel,0)
        except Exception as err:
            print ("Error:" + repr(err))
            self.terminate()



    def getCriticalTimeStep(self):
        # determine critical time step
        return PQ.PhysicalQuantity(1.0, 's')

    def terminate(self):
        #self.thermalAppRec.terminateAll()
        self.lammpsSolver.terminate()
        self.digimatSolver.terminate()
        self.mul2Solver.terminate()
        super(Workflow_1_1_3, self).terminate()

    def getApplicationSignature(self):
        return "Composelector workflow 1.0"

    def getAPIVersion(self):
        return "1.0"

    
if __name__=='__main__':

    useCaseID = 1
    execID = 1
    
    workflow = Workflow_1_1_3(targetTime=PQ.PhysicalQuantity(1.,'s'))
    # create and set lammps material properties
    workflow.setProperty(Property.ConstantProperty(1000, PropertyID.PID_SMILE_MOLECULAR_STRUCTURE, ValueType.Scalar, PQ.getDimensionlessUnit(), None, 0))
    workflow.setProperty(Property.ConstantProperty(100, PropertyID.PID_MOLECULAR_WEIGHT, ValueType.Scalar, 'mol', None, 0))
    workflow.setProperty(Property.ConstantProperty(0.2, PropertyID.PID_CROSSLINKER_TYPE, ValueType.Scalar, PQ.getDimensionlessUnit(), None, 0))
    workflow.setProperty(Property.ConstantProperty(12, PropertyID.PID_FILLER_DESIGNATION, ValueType.Scalar, PQ.getDimensionlessUnit(), None, 0))
    workflow.setProperty(Property.ConstantProperty(12, PropertyID.PID_CROSSLINKONG_DENSITY, ValueType.Scalar, PQ.getDimensionlessUnit(), None, 0))
    workflow.setProperty(Property.ConstantProperty(12, PropertyID.PID_FILLER_CONCENTRATION, ValueType.Scalar, PQ.getDimensionlessUnit(), None, 0))
    workflow.setProperty(Property.ConstantProperty(2000, PropertyID.PID_TEMPERATURE, ValueType.Scalar, 'degC', None, 0))
    workflow.setProperty(Property.ConstantProperty(500, PropertyID.PID_PRESSURE, ValueType.Scalar, 'atm', None, 0))
    workflow.setProperty(Property.ConstantProperty(1, PropertyID.PID_POLYDISPERSITY_INDEX, ValueType.Scalar, PQ.getDimensionlessUnit(), None, 0))
    workflow.setProperty(Property.ConstantProperty(1000, PropertyID.PID_SMILE_MODIFIER_MOLECULAR_STRUCTURE, ValueType.Scalar, PQ.getDimensionlessUnit(), None, 0))
    workflow.setProperty(Property.ConstantProperty(123, PropertyID.PID_SMILE_FILLER_MOLECULAR_STRUCTURE, ValueType.Scalar, PQ.getDimensionlessUnit(), None, 0))
    workflow.setProperty(Property.ConstantProperty(800, PropertyID.PID_DENSITY_OF_FUNCTIONALIZATION, ValueType.Scalar, PQ.getDimensionlessUnit(), None, 0))
    workflow.setProperty(Property.ConstantProperty(100, PropertyID.PID_InclusionYoung,  ValueType.Scalar, 'MPa', None, 0))
    workflow.setProperty(Property.ConstantProperty(0.2, PropertyID.PID_InclusionPoisson, ValueType.Scalar, PQ.getDimensionlessUnit(), None, 0))
    workflow.setProperty(Property.ConstantProperty(0.5, PropertyID.PID_InclusionVolumeFraction, ValueType.Scalar, PQ.getDimensionlessUnit(), None, 0))
    workflow.setProperty(Property.ConstantProperty(0.2, PropertyID.PID_InclusionAspectRatio, ValueType.Scalar, PQ.getDimensionlessUnit(), None, 0))

    # set useCaseID and ExecID 
    workflow.setMetadata(MetadataKeys.UseCaseID, useCaseID)
    workflow.setMetadata(MetadataKeys.ExecID, execID)
    workflow.solve()
    time = PQ.PhysicalQuantity(1.0, 's')
    bucklingLoad = workflow.getProperty(PropertyID.PID_CriticalLoadLevel, time).inUnitsOf('N').getValue()

    print ('OK')
