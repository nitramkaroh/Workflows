import sys
sys.path.extend(['..', '../../..'])
from mupif import *
import logging
log = logging.getLogger()
log.setLevel(logging.INFO)

import ApplicationsConfigs.DigimatConfig as DigimatConfig
import ApplicationsConfigs.DigimatPropertyID as DigimatPropertyID

cfg = DigimatConfig.digimatConfig()

import time as timeT

import ComposelectorSimulationTools.MIUtilities as miu

start = timeT.time()

class DemoDigimat(Workflow.Workflow):
   
    def __init__ (self):
        super(DemoDigimat, self).__init__(file='', workdir='')
        
        #locate nameserver
        ns = PyroUtil.connectNameServer(nshost=cfg.nshost, nsport=cfg.nsport, hkey=cfg.hkey)
        #connect to JobManager running on (remote) server and create a tunnel to it
        self.jobMan = PyroUtil.connectJobManager(ns, cfg.jobManName, cfg.hkey)
        log.info('Connected to JobManager')
        self.app1 = None
        log.info('Try to get application')
        try:
            self.app1 = PyroUtil.allocateApplicationWithJobManager( ns, self.jobMan, cfg.jobNatPorts[0], cfg.hkey, PyroUtil.SSHContext(sshClient=cfg.sshClient, options=cfg.options, sshHost=cfg.sshHost ) )
            log.info(self.app1)
        except Exception as e:
            log.exception(e)

            appsig=self.app1.getApplicationSignature()
            log.info("Working application 1 on server " + appsig)

    def solve(self):

        self.app1.solveStep(None)

        print( self.app1.getProperty(DigimatPropertyID.PID_CompositeAxialYoung       ).inUnitsOf('MPa').getValue() )
        print( self.app1.getProperty(DigimatPropertyID.PID_CompositeInPlaneYoung     ).inUnitsOf('MPa').getValue() )
        print( self.app1.getProperty(DigimatPropertyID.PID_CompositeInPlaneShear     ).inUnitsOf('MPa').getValue() )
        print( self.app1.getProperty(DigimatPropertyID.PID_CompositeTransverseShear  ).inUnitsOf('MPa').getValue() )
        print( self.app1.getProperty(DigimatPropertyID.PID_CompositeInPlanePoisson   ).getValue() )
        print( self.app1.getProperty(DigimatPropertyID.PID_CompositeTransversePoisson).getValue() )
        log.info("Sucessfully")
        
    def terminate(self):    
        self.app1.terminate()
        self.jobMan.terminate()
        super(DemoDigimat, self).terminate()
        log.info("Time elapsed %f s" % (timeT.time()-start))
    def getApplicationSignature(self):
        return "DemoDigimat workflow"

    def getAPIVersion(self):
        return "1.0"

def workflow(inputGUID, execGUID):
    # Define properties:units to export
    propsToExp = {"Matrix Young modulus": "MPa",
                  "Matrix Poisson ratio": "",
                  "Inclusion Young modulus": "MPa",
                  "Inclusion Poisson ratio": "",
                  "Inclusion volume fraction": "%",
                  "Inclusion aspect ratio": ""
                  }

    # Export data from database
    ExportedData = miu.ExportData("MI_Composelector", "Inputs-Outputs", inputGUID, propsToExp, miu.unitSystems.METRIC)

    matrixYoung = ExportedData["Matrix Young modulus"]
    matrixPoisson = ExportedData["Matrix Poisson ratio"]
    inclusionYoung = ExportedData["Inclusion Young modulus"]
    inclusionPoisson = ExportedData["Inclusion Poisson ratio"]
    inclusionVolumeFraction = ExportedData["Inclusion volume fraction"]
    inclusionAspectRatio = ExportedData["Inclusion aspect ratio"]

    try:
        myapp = DemoDigimat()
        myapp.app1.setProperty(Property.ConstantProperty(matrixYoung, DigimatPropertyID.PID_MatrixYoung,               ValueType.Scalar, "MPa"))
        myapp.app1.setProperty(Property.ConstantProperty(inclusionYoung, DigimatPropertyID.PID_InclusionYoung,            ValueType.Scalar, "MPa"))
        myapp.app1.setProperty(Property.ConstantProperty(matrixPoisson, DigimatPropertyID.PID_MatrixPoisson,             ValueType.Scalar, "none"))
        myapp.app1.setProperty(Property.ConstantProperty(inclusionPoisson, DigimatPropertyID.PID_InclusionPoisson,          ValueType.Scalar, "none"))
        myapp.app1.setProperty(Property.ConstantProperty(inclusionVolumeFraction, DigimatPropertyID.PID_InclusionVolumeFraction,   ValueType.Scalar, "none"))
        myapp.app1.setProperty(Property.ConstantProperty(inclusionAspectRatio, DigimatPropertyID.PID_InclusionAspectRatio,      ValueType.Scalar, "none"))

        myapp.solve()

        compositeAxialYoung = myapp.app1.getProperty(DigimatPropertyID.PID_CompositeAxialYoung).inUnitsOf('MPa').getValue()
        compositeInPlaneYoung = myapp.app1.getProperty(DigimatPropertyID.PID_CompositeInPlaneYoung).inUnitsOf('MPa').getValue()
        compositeInPlaneShear = myapp.app1.getProperty(DigimatPropertyID.PID_CompositeInPlaneShear).inUnitsOf('MPa').getValue()
        compositeTransverseShear = myapp.app1.getProperty(DigimatPropertyID.PID_CompositeTransverseShear).inUnitsOf('MPa').getValue()
        compositeInPlanePoisson = myapp.app1.getProperty(DigimatPropertyID.PID_CompositeInPlanePoisson).getValue()
        compositeTransversePoisson = myapp.app1.getProperty(DigimatPropertyID.PID_CompositeTransversePoisson).getValue()

        myapp.terminate()

        ImportHelper = miu.Importer("MI_Composelector", "Inputs-Outputs", ["Inputs/Outputs"])
        ImportHelper.CreateAttribute("Axial Young modulus", compositeAxialYoung, "MPa")
        ImportHelper.CreateAttribute("In-plane Young modulus", compositeInPlaneYoung, "MPa")
        ImportHelper.CreateAttribute("In-plane shear modulus", compositeInPlaneShear, "MPa")
        ImportHelper.CreateAttribute("Tranverse shear modulus", compositeTransverseShear, "MPa")
        ImportHelper.CreateAttribute("In-plane Poisson ratio", compositeInPlanePoisson, "")
        ImportHelper.CreateAttribute("Transverse Poisson ratio", compositeTransversePoisson, "")
        return ImportHelper

    except APIError.APIError as err:
        print ("Mupif API for DIGIMAT-MF error: "+ repr(err))
    except Exception as err:
        print ("Error: " + repr(err))
    except:
        print ("Unkown error.")

# if __name__=='__main__':
#     workflow('exampleGUID')


