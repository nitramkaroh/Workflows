import sys
sys.path.extend(['/home/nitram/Documents/work/MUPIF/commit'])
from mupif import *
import Pyro4
import logging
log = logging.getLogger()
import time as timeT
import mupif.Physics.PhysicalQuantities as PQ
import LAMMPS_v3 as lammps
import DigimatConfig
import DigimatPropertyID
digimat_cfg = DigimatConfig.digimatConfig()


class Workflow1(Workflow.Workflow):
   
    def __init__(self, targetTime=PQ.PhysicalQuantity('0 s')):
        """
        Initializes the workflow. As the workflow is non-stationary, we allocate individual 
        applications and store them within a class.
        """
        super(Workflow1, self).__init__(file='', workdir='', targetTime=targetTime)

        #list of recognized input porperty IDs
        self.myInputPropIDs = [PropertyID.PID_SMILE_MOLECULAR_STRUCTURE,PropertyID.PID_MOLECULAR_WEIGHT, PropertyID.PID_CROSSLINKER_TYPE,PropertyID.PID_FILLER_DESIGNATION, PropertyID.PID_CROSSLINKONG_DENSITY,PropertyID.PID_FILLER_CONCENTRATION, PropertyID.PID_TEMPERATURE, PropertyID.PID_PRESSURE, PropertyID.PID_POLYDISPERSITY_INDEX,PropertyID.PID_SMILE_MODIFIER_MOLECULAR_STRUCTURE,PropertyID.PID_SMILE_FILLER_MOLECULAR_STRUCTURE,PropertyID.PID_DENSITY_OF_FUNCTIONALIZATION, DigimatPropertyID.PID_InclusionYoung, DigimatPropertyID.PID_InclusionPoisson, DigimatPropertyID.PID_InclusionVolumeFraction, DigimatPropertyID.PID_InclusionAspectRatio]
        # list of compulsory IDs
        self.myCompulsoryPropIDs = self.myInputPropIDs

        #list of recognized output property IDs
        self.myOutPropIDs =  [DigimatPropertyID.PID_CompositeAxialYoung,DigimatPropertyID.PID_CompositeInPlaneYoung, DigimatPropertyID.PID_CompositeInPlaneShear, DigimatPropertyID.PID_CompositeTransverseShear, DigimatPropertyID.PID_CompositeInPlanePoisson, DigimatPropertyID.PID_CompositeTransversePoisson]

        #dictionary of input properties (values)
        self.myInputProps = {}
        #dictionary of output properties (values)
        self.myOutProps = {}

        #locate nameserver
        ns = PyroUtil.connectNameServer(digimat_cfg.nshost, digimat_cfg.nsport, digimat_cfg.hkey)
        #connect to JobManager running on (remote) server
        self.digimatJobMan = PyroUtil.connectJobManager(ns, digimat_cfg.jobManName,digimat_cfg.hkey)
        
        self.lammpsSolver = lammps.emailAPI(None) # local LAMMPS API instances
        self.digimatSolver = None
        #allocate the Digimat remote instance
        try:
            self.digimatSolver = PyroUtil.allocateApplicationWithJobManager( ns, self.digimatJobMan, digimat_cfg.jobNatPorts[0], digimat_cfg.hkey)

            
            log.info('Created digimat job')
        except Exception as e:
            log.exception(e)
        else:
            if ((self.lammpsSolver is not None) and (self.digimatSolver is not None)):
                lammpsSolverSignature=self.lammpsSolver.getApplicationSignature()
                log.info("Working lammps solver on server " + lammpsSolverSignature)
                digimatSolverSignature=self.digimatSolver.getApplicationSignature()
                log.info("Working digimat solver on server " + digimatSolverSignature)
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
            self.lammpsSolver.solveStep (useCaseID, execID, TimeStep.TimeStep(0.0, 0.1, 1, 's'))
            # get result of the simulation
            matrixYoung = self.lammpsSolver.getProperty(PropertyID.PID_EModulus, 0.0)
            matrixPoisson = self.lammpsSolver.getProperty(PropertyID.PID_PoissonRatio, 0.0)

        except Exception as err:
            print ("Error:" + repr(err))
            self.terminate()


        # digimat
        try:
            # map properties from lammps to properties of Digimat
            matrixYoung.propID = DigimatPropertyID.PID_MatrixYoung
            matrixPoisson.propID = DigimatPropertyID.PID_MatrixPoisson
            self.digimatSolver.setProperty(matrixYoung)
            self.digimatSolver.setProperty(matrixPoisson)

            # fixed properties - taken form the database
            self.digimatSolver.setProperty(self.myInputProps[DigimatPropertyID.PID_InclusionYoung])
            self.digimatSolver.setProperty(self.myInputProps[DigimatPropertyID.PID_InclusionPoisson])
            self.digimatSolver.setProperty(self.myInputProps[DigimatPropertyID.PID_InclusionVolumeFraction])
            self.digimatSolver.setProperty(self.myInputProps[DigimatPropertyID.PID_InclusionAspectRatio])
        except Exception as err:
            print ("Setting Digimat params failed: " + repr(err));
            self.terminate()

        try:
            # solve digimat part
            log.info("Running digimat")
            self.digimatSolver.solveStep(None)
            ## set the desired properties
            self.myOutProps[PropertyID.PID_CompositeAxialYoung] = self.digimatSolver.getProperty(DigimatPropertyID.PID_CompositeAxialYoung)
            self.myOutProps[PropertyID.PID_CompositeInPlaneYoung] = self.digimatSolver.getProperty(DigimatPropertyID.PID_CompositeInPlaneYoung)
            self.myOutProps[PropertyID.PID_CompositeInPlaneShear] = self.digimatSolver.getProperty(DigimatPropertyID.PID_CompositeInPlaneShear)
            self.myOutProps[PropertyID.PID_CompositeTransverseShear] = self.digimatSolver.getProperty(DigimatPropertyID.PID_CompositeTransverseShear)
            self.myOutProps[PropertyID.PID_CompositeInPlanePoisson] = self.digimatSolver.getProperty(DigimatPropertyID.PID_CompositeInPlanePoisson)
            self.myOutProps[PropertyID.PID_CompositeTransversePoisson] = self.digimatSolver.getProperty(DigimatPropertyID.PID_CompositeTransversePoisson)
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
        super(Workflow1, self).terminate()

    def getApplicationSignature(self):
        return "Composelector workflow 1.0"

    def getAPIVersion(self):
        return "1.0"

    
if __name__=='__main__':

    useCaseID = 1
    execID = 1
    
    workflow1 = Workflow1(targetTime=PQ.PhysicalQuantity(1.,'s'))
    # create and set lammps material properties
    workflow1.setProperty(Property.ConstantProperty(1000, PropertyID.PID_SMILE_MOLECULAR_STRUCTURE, ValueType.Scalar, PQ.getDimensionlessUnit(), None, 0))
    workflow1.setProperty(Property.ConstantProperty(100, PropertyID.PID_MOLECULAR_WEIGHT, ValueType.Scalar, 'mol', None, 0))
    workflow1.setProperty(Property.ConstantProperty(0.2, PropertyID.PID_CROSSLINKER_TYPE, ValueType.Scalar, PQ.getDimensionlessUnit(), None, 0))
    workflow1.setProperty(Property.ConstantProperty(12, PropertyID.PID_FILLER_DESIGNATION, ValueType.Scalar, PQ.getDimensionlessUnit(), None, 0))
    workflow1.setProperty(Property.ConstantProperty(12, PropertyID.PID_CROSSLINKONG_DENSITY, ValueType.Scalar, PQ.getDimensionlessUnit(), None, 0))
    workflow1.setProperty(Property.ConstantProperty(12, PropertyID.PID_FILLER_CONCENTRATION, ValueType.Scalar, PQ.getDimensionlessUnit(), None, 0))
    workflow1.setProperty(Property.ConstantProperty(2000, PropertyID.PID_TEMPERATURE, ValueType.Scalar, 'degC', None, 0))
    workflow1.setProperty(Property.ConstantProperty(500, PropertyID.PID_PRESSURE, ValueType.Scalar, 'atm', None, 0))
    workflow1.setProperty(Property.ConstantProperty(1, PropertyID.PID_POLYDISPERSITY_INDEX, ValueType.Scalar, PQ.getDimensionlessUnit(), None, 0))
    workflow1.setProperty(Property.ConstantProperty(1000, PropertyID.PID_SMILE_MODIFIER_MOLECULAR_STRUCTURE, ValueType.Scalar, PQ.getDimensionlessUnit(), None, 0))
    workflow1.setProperty(Property.ConstantProperty(123, PropertyID.PID_SMILE_FILLER_MOLECULAR_STRUCTURE, ValueType.Scalar, PQ.getDimensionlessUnit(), None, 0))
    workflow1.setProperty(Property.ConstantProperty(800, PropertyID.PID_DENSITY_OF_FUNCTIONALIZATION, ValueType.Scalar, PQ.getDimensionlessUnit(), None, 0))   
    workflow1.setProperty(Property.ConstantProperty(100, DigimatPropertyID.PID_InclusionYoung,  ValueType.Scalar, 'MPa', None, 0))
    workflow1.setProperty(Property.ConstantProperty(0.2, DigimatPropertyID.PID_InclusionPoisson, ValueType.Scalar, PQ.getDimensionlessUnit(), None, 0))
    workflow1.setProperty(Property.ConstantProperty(0.5, DigimatPropertyID.PID_InclusionVolumeFraction, ValueType.Scalar, PQ.getDimensionlessUnit(), None, 0))
    workflow1.setProperty(Property.ConstantProperty(0.2, DigimatPropertyID.PID_InclusionAspectRatio, ValueType.Scalar, PQ.getDimensionlessUnit(), None, 0))

    
    workflow1.solve()
    time = PQ.PhysicalQuantity(1.0, 's')
    #print (workflow1.myOutProps)

    compositeAxialYoung = workflow1.getProperty(PropertyID.PID_CompositeAxialYoung, time).inUnitsOf('MPa').getValue()
    compositeInPlaneYoung = workflow1.getProperty(PropertyID.PID_CompositeInPlaneYoung, time).inUnitsOf('MPa').getValue()
    compositeInPlaneShear = workflow1.getProperty(PropertyID.PID_CompositeInPlaneShear, time).inUnitsOf('MPa').getValue()
    compositeTransverseShear = workflow1.getProperty(PropertyID.PID_CompositeTransverseShear, time).inUnitsOf('MPa').getValue()
    compositeInPlanePoisson = workflow1.getProperty(PropertyID.PID_CompositeInPlanePoisson,time).getValue()
    compositeTransversePoisson = workflow1.getProperty(PropertyID.PID_CompositeTransversePoisson, time).getValue()

    print ('OK')
