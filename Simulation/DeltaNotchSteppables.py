import os
import sys
import CompuCell
import random 
import math
import bionetAPI
import CompuCellSetup  
from PySteppables import *
from PlayerPython import *
from math import *
from PySteppablesExamples import MitosisSteppableBase

class InitialCondition(MitosisSteppableBase):
    def __init__(self,_simulator,_frequency):
        MitosisSteppableBase.__init__(self,_simulator,_frequency)
        # Calcuate Cell Death due to Threshold and BM touching
        self.NoOfDeathThre = 0
        self.NoOfDeathBM   = 0
    def start(self):
        radi = self.dim.x/2
        # Assign Property for Cell ID = 1
        cells_to_die=[]
        for cell in self.cellListByType(1):
            # Assign stochastic initial conditions for stem cell/TA cell volume
            if cell.yCOM < self.dim.y*0.7:
                cell.targetVolume = 1000*random.uniform(1.0,1.5)# Make the initial target Volume of diff cells constant
            else:
                cell.targetVolume = 1000
            cell.lambdaVolume = 5

    # Make sure ExternalPotential plugin is loaded
            # Direct Cell migration
    #         cell.lambdaVecX=-0 # force component pointing along X axis - towards positive X's
    #         cell.lambdaVecY=-0 # force component pointing along Y axis - towards negative Y's
    #         cell.lambdaVecZ=0.0 # force component pointing along Z axis 
    
        # Basal Membrane Generation
        radi = int(self.dim.x/2)
        pt=CompuCell.Point3D(0,0,0)
        Wall = self.potts.createCellG(pt) # the arguments are (type,pt,xSize,ySize,zSize)  
        Wall.type = 2
        totalNumOfSeg = 15
        segm = []   
        cellsToDelete=[]
        
        for pty in range(0,self.dim.y):
            for ptx in range(0,self.dim.x):
                for ptz in range(0,self.dim.z):
                    pt.y=pty
                    pt.x=ptx
                    pt.z=ptz        
                    
                    # Generate Tube structure
                    if pty>= radi:
                        if ((ptx-radi)**2+(ptz-radi)**2)<=radi**2 and ((ptx-radi)**2+(ptz-radi)**2)>=(radi-5)**2:
                            overwrittenCell=self.cellField.get(pt)  

                            self.cellField.set(pt,Wall)
                            self.cleanDeadCells()
                            
                            ang = math.atan2(ptx-radi,ptz-radi)*(180./math.pi)+180
                            segN= int(ang/(360/totalNumOfSeg))
                            if segN ==15:
                                segN = 0
                            segm.append([segN,pt])
                        
                        elif ((ptx-radi)**2+(ptz-radi)**2)>radi**2:
                            self.cellField.set(pt,CompuCell.getMediumCell())
                            self.cleanDeadCells()
                    # Generate Semi-Sphere Sturcture ((Bottom of the crypt)
                 
                    elif pty < radi:
                        if ((ptx-radi)**2+(pty-radi)**2+(ptz-radi)**2)<=radi**2 and ((ptx-radi)**2+(pty-radi)**2+(ptz-radi)**2)>=(radi-5)**2:
                            
                            
                            self.cellField.set(pt,Wall)
                            self.cleanDeadCells()
                            
                        elif ((ptx-radi)**2+(pty-radi)**2+(ptz-radi)**2)>radi**2:
                            
                            overwrittenCell=self.cellField.get(pt)
                            if pt.x!=0 and pt.y!=0 and pt.z!=0: # this line is essential because you do not wan to remove wall cell which sits at (0,0,0)
                                self.cellField.set(pt,CompuCell.getMediumCell())
                                self.cleanDeadCells()
                                
        # now we can overwrite (0,0,0) with medium
        pt.x=0
        pt.y=0
        pt.z=0
        self.cellField.set(pt,CompuCell.getMediumCell())
        
        
        self.divideCellOrientationVectorBased(Wall,0,0,1)
        for divN in range(0,3):        
            cell_to_divide=[]    
            # iterating over cells of type 2        
            for cell in self.cellListByType(2):
                cell_to_divide.append(cell)
            for cell in cell_to_divide:
                DiviVectorX = (cell.xCOM-radi)
                DiviVectorZ = (cell.zCOM-radi) 
                self.divideCellOrientationVectorBased(cell,-DiviVectorZ,0,DiviVectorX) 
        for divN in range(0,4):
            cell_to_divide=[]    
            # iterating over cells of type 2        
            for cell in self.cellListByType(2):
                cell_to_divide.append(cell)
            for cell in cell_to_divide:
                self.divideCellOrientationVectorBased(cell,0,1,0) 
                
        # Assign property for Cell type = 2
        for cell in self.cellListByType(2): 
            cell.targetVolume = cell.volume
            cell.lambdaVolume = 10000000
            
        # Making Scientific Plot
        self.pW=CompuCellSetup.viewManager.plotManager.getNewPlotWindow()
        if not self.pW:
            return
        #Plot Title - properties           
        self.pW.setTitle("No of Cell Death")
        self.pW.setTitleSize(12)
        self.pW.setTitleColor("Green") # you may choose different color - type its name
        
        #plot background
        self.pW.setPlotBackgroundColor("orange") # you may choose different color - type its name
        
        # properties of x axis
        self.pW.setXAxisTitle("No Of MCS")
        self.pW.setXAxisTitleSize(10)      
        self.pW.setXAxisTitleColor("blue")  # you may choose different color - type its name            
        
        # properties of y axis
        self.pW.setYAxisTitle("No of Cell Death")        
        #self.pW.setYAxisLogScale()
        self.pW.setYAxisTitleSize(10)        
        self.pW.setYAxisTitleColor("red")  # you may choose different color - type its name                                
        
        # choices for style are NoCurve,Lines,Sticks,Steps,Dots
        self.pW.addPlot("NoOfCellDeathThre",_style='Dots')
        self.pW.addPlot("NoOfCellDeathBM",_style='Dots')
        #self.pW.addPlot("DATA SERIES 2",_style='Steps') # you may add more than one data series
        
        # plot MCS
        self.pW.changePlotProperty("NoOfCellDeathThre","LineWidth",5)
        self.pW.changePlotProperty("NoOfCellDeathThre","LineColor","red")
        self.pW.changePlotProperty("NoOfCellDeathBM","LineWidth",5)
        self.pW.changePlotProperty("NoOfCellDeathBM","LineColor","black")    
        
        self.pW.addGrid()
        #adding automatically generated legend
        # default position is at the bottom of the plot but here we put it at the top
        self.pW.addAutoLegend("top")
        
        self.clearFlag=False

    def updateAttributes(self):
        childCell = self.mitosisSteppable.childCell
        parentCell = self.mitosisSteppable.parentCell
        childCell.type = parentCell.type
#         parentCell.targetVolume = 1000
#         childCell.targetVolume = 1000
#         childCell.lambdaVolume = parentCell.lambdaVolume;

    def step(self,mcs):
        cells_to_die=[]
        for cell in self.cellList:
            if cell.type == 1:
                # Assume Growth only happens at the bottom of crypt
                if cell.targetVolume:
                    if cell.yCOM < self.dim.y*0.7:
                    # Program Cell Growth
                        cell.targetVolume+= 1    
                    # The diff cells remain unchanged    
                    else:
                        cell.targetVolume = 1000
                    
                # Program Cell Death
                # Set up threshold to kill cells when cells go above the threshold
                if cell.yCOM > self.dim.y-10:
                    cells_to_die.append(cell)
                    self.NoOfDeathThre+=1
                cellNeighborList=self.getCellNeighbors(cell) # generates list of neighbors of cell 'cell'
                wallflag=0
                # Kill cells when the cells not touching BM
                for neighbor in cellNeighborList:
                    # neighborSurfaceData.neighborAddress is a pointer to one of 'cell' neighbors stired in cellNeighborList
                    #IMPORTANT: cell may have Medium (NULL pointer) as a neighbor. therefore before accessing neighbor we first check if it is no Medium
                    if neighbor.neighborAddress: 
                        # Detect the cells touching BM
                        if neighbor.neighborAddress.type == 2:
                            wallflag=1
                if wallflag==0:
                    # Delete the cells without contacting BM
                    cells_to_die.append(cell)
                    if mcs > 5:
                        self.NoOfDeathBM+=1
                        print "!!!!Dieing cell",cell.id,cell.type,cell.volume,cell.targetVolume
                   
        # Cell Killing program
        print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------",len(cells_to_die)
        for cell in cells_to_die:    
            cell.targetVolume = 0
  
            self.deleteCell(cell)
            # The next line araises the dead cell register      
            self.cleanDeadCells()
              
        if mcs > 5:
#             self.NoOfDeathBM -=15   
            #self.pW.eraseAllData() # this is how you erase previous content of the plot
            self.pW.addDataPoint("NoOfCellDeathThre",mcs,self.NoOfDeathThre) # arguments are (name of the data series, x, y)
            self.pW.addDataPoint("NoOfCellDeathBM",mcs,self.NoOfDeathBM) # arguments are (name of the data series, x, y)
            self.pW.showAllPlots()    



class DeltaNotchClass(SteppableBasePy):
    def __init__(self,_simulator,_frequency):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        bionetAPI.initializeBionetworkManager(self.simulator)
        
    def start(self):
        #Loading model
        Name = "DeltaNotch"
        Key  = "DN"
        simulationDir=os.path.dirname (os.path.abspath( __file__ ))
        Path= os.path.join(simulationDir,"SLIMI4R.sbml")
        Path=os.path.abspath(Path) # normalizing path
        
        IntegrationStep = 0.2
        bionetAPI.loadSBMLModel(Name, Path, Key, IntegrationStep)
        
        bionetAPI.addSBMLModelToTemplateLibrary(Name,"TypeA")
        bionetAPI.initializeBionetworks()
        
        #Initial conditions
        for cell in self.cellListByType(1):
            D = random.uniform(0.9,1.0)
            N = random.uniform(0.9,1.0)
            B = random.uniform(0.9,1.0)
            R = random.uniform(0.9,1.0)
            bionetAPI.setBionetworkValue("DN_D",D,cell.id)
            bionetAPI.setBionetworkValue("DN_N",N,cell.id)
            bionetAPI.setBionetworkValue("DN_B",B,cell.id)
            bionetAPI.setBionetworkValue("DN_R",R,cell.id)
            cellDict=CompuCell.getPyAttrib(cell)
            cellDict["D"]=D
            cellDict["N"]=N
            cellDict["B"]=B
            cellDict["R"]=R
    def step(self,mcs):
        for cell in self.cellListByType(1):
            Davg=0.0; Navg=0.0; nn=0
            cellNeighborList=self.getCellNeighbors(cell)
            for neighbor , commonSurfaceArea in self.getCellNeighborDataList(cell):                
                if (neighbor and neighbor.type == 1):
                    #print "~~~~~~~~~~~~~~~~~~neighbor.id=",neighbor.id,",neighbor.type=",neighbor.type,",Dll=",bionetAPI.getBionetworkValue("DN_D",neighbor.id)
                    nn+=1
                    Davg+=bionetAPI.getBionetworkValue("DN_D",neighbor.id)
                    Navg+=bionetAPI.getBionetworkValue("DN_N",neighbor.id)
            if (nn>0):
                Davg=Davg/nn
                Navg=Navg/nn  
            #print "~~~~~~~~~~~~~~~~~~cell.id=",cell.id,"\t# of neighbor cells=",nn,"\tDavg=",Davg,"\tNavg=",Navg            
# Make sure CenterOfMass plugin is loaded
# READ ONLY ACCES
# Wnt Gradient       
            yCOM=((1.0*(self.dim.y-cell.yCOM))/self.dim.y)
            WNT = math.pow(yCOM,3)
            bionetAPI.setBionetworkValue("DN_Davg",Davg,cell.id)
            bionetAPI.setBionetworkValue("DN_Navg",Navg,cell.id)
            bionetAPI.setBionetworkValue("DN_WNT",WNT,cell.id)

            cellDict=CompuCell.getPyAttrib(cell)
            cellDict["D"]=bionetAPI.getBionetworkValue("DN_D",cell.id)
            cellDict["N"]=bionetAPI.getBionetworkValue("DN_N",cell.id)
            cellDict["B"]=bionetAPI.getBionetworkValue("DN_B",cell.id)
            cellDict["R"]=bionetAPI.getBionetworkValue("DN_R",cell.id)
            
        bionetAPI.timestepBionetworks()
class MitosisSteppable(MitosisSteppableBase):
    def __init__(self, _simulator, _frequency=1):
        MitosisSteppableBase.__init__(self, _simulator, _frequency)
    
    def step(self,mcs):
        cells_to_divide=[]
        for cell in self.cellListByType(1):
            if cell.volume > 1800:
                cells_to_divide.append(cell)
                
        for cell in cells_to_divide:
            self.divideCellAlongMinorAxis(cell)

    def updateAttributes(self):
        childCell = self.mitosisSteppable.childCell
        parentCell = self.mitosisSteppable.parentCell
        childCell.type = parentCell.type
                
        parentCell.targetVolume = 1000
        childCell.targetVolume = 1000

        childCell.lambdaVolume = parentCell.lambdaVolume;
        
        bionetAPI.copyBionetworkFromParent( parentCell, childCell )
        childCellDict=CompuCell.getPyAttrib(childCell)
        parentCellDict=CompuCell.getPyAttrib(parentCell)
        childCellDict["D"]=random.uniform(0.9,1.0)*parentCellDict["D"]
        childCellDict["N"]=random.uniform(0.9,1.0)*parentCellDict["N"]
        childCellDict["B"]=random.uniform(0.9,1.0)*parentCellDict["B"]
        childCellDict["R"]=random.uniform(0.9,1.0)*parentCellDict["R"]

class ExtraFields(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    
    def setScalarFields(self,_field1,_field2,_field3,_field4):
        self.scalarField1=_field1
        self.scalarField2=_field2  
        self.scalarField3=_field3
        self.scalarField4=_field4  
    def step(self,mcs):
        clearScalarValueCellLevel(self.scalarField1)
        clearScalarValueCellLevel(self.scalarField2)
        clearScalarValueCellLevel(self.scalarField3)
        clearScalarValueCellLevel(self.scalarField4)
        for cell in self.cellListByType(1):
            cellDict=CompuCell.getPyAttrib(cell)
            fillScalarValueCellLevel(self.scalarField1,cell,cellDict["D"])
            fillScalarValueCellLevel(self.scalarField2,cell,cellDict["N"])
            fillScalarValueCellLevel(self.scalarField3,cell,cellDict["B"])
            fillScalarValueCellLevel(self.scalarField4,cell,cellDict["R"])
