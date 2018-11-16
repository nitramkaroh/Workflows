#Common configuration for running examples in local, ssh or VPN mode
import sys, os, os.path
import Pyro4
import logging
log = logging.getLogger()

class digimatConfig(object): 
    def __init__(self):
        
        Pyro4.config.SERIALIZER="pickle"
        Pyro4.config.PICKLE_PROTOCOL_VERSION=2 #to work with python 2.x and 3.x
        Pyro4.config.SERIALIZERS_ACCEPTED={'pickle'}
        Pyro4.config.SERVERTYPE="multiplex"

        #Absolute path to mupif directory - used in JobMan2cmd
        mupif_dir = os.path.abspath(os.path.join(os.getcwd(), "../.."))
        sys.path.append(mupif_dir)
        mupif_dir = os.path.abspath(os.path.join(os.getcwd(), ".."))
        sys.path.append(mupif_dir)
        
        #commmon attributes
        #Password for accessing nameServer and applications
        self.hkey = 'mupif-secret-key'
        #Name of job manager
        self.jobManName='Mupif.JobManager@Digimat'
        #Name of first application
        self.appName = 'MuPIFServer'

        #Jobs in JobManager
        #Range of ports to be assigned on the server to jobs
        self.portsForJobs=( 9000, 9100 )
        #NAT client ports used to establish ssh connections
        self.jobNatPorts = list(range(6000, 6100))

        #Maximum number of jobs
        self.maxJobs=20
        #Auxiliary port used to communicate with application daemons on a local computer
        self.socketApps=10000
        #Main directory for transmitting files
        self.jobManWorkDir='.'
        #Path to JobMan2cmd.py 
        self.jobMan2CmdPath = "/../tools/JobMan2cmd.py"
        
        #NAME SERVER
        #IP/name of a name server
        self.nshost = '172.30.0.1'
        #Port of name server
        self.nsport = 9090
        
        #SERVER for a single job or for JobManager
        #IP/name of a server's daemon
        self.server = '172.30.0.112'
        #Port of server's daemon
        self.serverPort = 44388
        #Nat IP/name
        self.serverNathost = None
        #Nat port
        self.serverNatport = None
        self.jobNatPorts = [None]

        self.sshHost = ''
        self.sshClient='manual'
        self.options = ''
        
        
        self.jobManName='eX_DigimatMF_JobManager'#Name of job manager
