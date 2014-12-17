import sys,time
from os import environ
from os import getcwd
import string
sys.path.append(environ["PYTHON_MODULE_PATH"])
sys.path.append(environ["SWIG_LIB_INSTALL_DIR"])

def configureSimulation(sim):
   import CompuCellSetup
   from XMLUtils import ElementCC3D
   cc3d=ElementCC3D("CompuCell3D")
   
   potts=cc3d.ElementCC3D("Potts")
   potts.ElementCC3D("Dimensions",{"x":70,"y":150,"z":70})
   potts.ElementCC3D("Steps",{},10000)
   potts.ElementCC3D("Temperature",{},15)#original 30   
   potts.ElementCC3D("NeighborOrder",{},3)

   cellType=cc3d.ElementCC3D("Plugin",{"Name":"CellType"})
   cellType.ElementCC3D("CellType", {"TypeName":"Medium","TypeId":"0"})
   cellType.ElementCC3D("CellType", {"TypeName":"TypeA" ,"TypeId":"1"}) # Epithelial cells
   cellType.ElementCC3D("CellType", {"TypeName":"BM" ,"TypeId":"2"})    # Meschymal cells

# assign Energe
   contact=cc3d.ElementCC3D("Plugin",{"Name":"Contact"})
   contact.ElementCC3D("Energy", {"Type1":"Medium", "Type2":"Medium"},0)
   contact.ElementCC3D("Energy", {"Type1":"Medium", "Type2":"TypeA"},8)
   contact.ElementCC3D("Energy", {"Type1":"Medium", "Type2":"BM"},10)
   contact.ElementCC3D("Energy", {"Type1":"TypeA",  "Type2":"TypeA"},2)
   contact.ElementCC3D("Energy", {"Type1":"TypeA",  "Type2":"BM"},1)
   contact.ElementCC3D("Energy", {"Type1":"BM", "Type2":"BM"},0)
   contact.ElementCC3D("NeighborOrder",{},5)

   volume = cc3d.ElementCC3D("Plugin",{"Name":"VolumeLocalFlex"})
   ntp = cc3d.ElementCC3D("Plugin",{"Name":"NeighborTracker"})
   epb = cc3d.ElementCC3D("Plugin",{"Name":"ExternalPotential"})
   ctm = cc3d.ElementCC3D("Plugin",{"Name":"CenterOfMass"})
   pit = cc3d.ElementCC3D("Plugin",{"Name":"PixelTracker"})
   BPT = cc3d.ElementCC3D("Plugin",{"Name":"BoundaryPixelTracker"})
 
   uipd = cc3d.ElementCC3D("Steppable",{"Type":"UniformInitializer"})
   region = uipd.ElementCC3D("Region")
   region.ElementCC3D("BoxMin",{"x":6,  "y":6,  "z":6})
   region.ElementCC3D("BoxMax",{"x":65,  "y":130,  "z":65})
   region.ElementCC3D("Types",{},"TypeA")
   region.ElementCC3D("Width", {},10)
      
   CompuCellSetup.setSimulationXMLDescription(cc3d)

import CompuCellSetup
sim,simthread = CompuCellSetup.getCoreSimulationObjects()
configureSimulation(sim)

import CompuCell
CompuCellSetup.initializeSimulationObjects(sim,simthread)
pyAttributeAdder,dictAdder=CompuCellSetup.attachDictionaryToCells(sim)

#Add Python steppables here
steppableRegistry=CompuCellSetup.getSteppableRegistry()

from DeltaNotchSteppables import InitialCondition
initialCondition=InitialCondition(_simulator=sim,_frequency=1)
steppableRegistry.registerSteppable(initialCondition)

from DeltaNotchSteppables import Growth
Growth=Growth(_simulator=sim,_frequency=1)
steppableRegistry.registerSteppable(Growth)

from DeltaNotchSteppables import DeltaNotchClass
deltaNotchClass=DeltaNotchClass(_simulator=sim,_frequency=1)
steppableRegistry.registerSteppable(deltaNotchClass)

from DeltaNotchSteppables import MitosisSteppable
mitosisSteppable=MitosisSteppable(_simulator=sim,_frequency=1)
steppableRegistry.registerSteppable(mitosisSteppable)

#Create extra player fields here or add attributes
from DeltaNotchSteppables import ExtraFields
extraFields=ExtraFields(_simulator=sim,_frequency=1)
steppableRegistry.registerSteppable(extraFields)

from DeltaNotchSteppables import Plotting1
Plotting1=Plotting1(_simulator=sim,_frequency=1)
steppableRegistry.registerSteppable(Plotting1)

from DeltaNotchSteppables import Plotting2
Plotting2=Plotting2(_simulator=sim,_frequency=1)
steppableRegistry.registerSteppable(Plotting2)

from DeltaNotchSteppables import OutputData
outputData=OutputData(_simulator=sim,_frequency=1)
steppableRegistry.registerSteppable(outputData)

CompuCellSetup.mainLoop(sim,simthread,steppableRegistry)
