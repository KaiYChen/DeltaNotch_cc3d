import os
import sys
import CompuCell
import random 
import math
import CompuCellSetup  
from PySteppables import *
from PlayerPython import *
from math import *
from PySteppablesExamples import MitosisSteppableBase

tVol = 800
# Calcuate Cell Death due to Threshold and BM touching
NoOfDeathBM = 0
NoOfDeathThre = 0
class InitialCondition(MitosisSteppableBase):
    def __init__(self,_simulator,_frequency):
        MitosisSteppableBase.__init__(self,_simulator,_frequency)
    def start(self):
        radi = self.dim.x/2
        # Assign Property for Cell ID = 1
        cells_to_die=[]
        for cell in self.cellListByType(1):
            cell.targetVolume = tVol
            cell.lambdaVolume = 5
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
        cells_to_die=[]
        for cell in self.cellList:
            if cell.type == 1:                    
                # Program Cell Death
                # Set up threshold to kill cells when cells go above the threshold
                if cell.yCOM > self.dim.y-10:
                    cells_to_die.append(cell)
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
                   
        # Cell Killing program
        for cell in cells_to_die:    
            cell.targetVolume = 0
            self.deleteCell(cell)
            # The next line araises the dead cell register      
            self.cleanDeadCells()

        self.BreakCells()# call new BreakCells() to generate apical and basal cells from Epi cells
        
    def BreakCells(self):
        # iterating over cells of type 1        
        for cell in self.cellListByType(1): # list of  cell types (capitalized)
            # you can access/manipulate cell properties here
            # print "id=",cell.id," type=",cell.type
            Apical=self.potts.createCell()
            Basal=self.potts.createCell()
            Apical.type=3
            Basal.type=4
            # Make sure PixelTracker plugin is loaded
            pixelList=self.getCellPixelList(cell)
            for PTD in pixelList:
                #print "pixel of cell id=",cell.id," type:",cell.type, " = ",PTD.pixel," number of pixels=",pixelList.numberOfPixels()
                R = random.randint(1,4)
                if R == 4:
                    x = PTD.pixel.x
                    y = PTD.pixel.y
                    z = PTD.pixel.z
                    self.cellField[x:x,y:y,z:z]=Basal
                elif R== 3:
                    x = PTD.pixel.x
                    y = PTD.pixel.y
                    z = PTD.pixel.z
                    self.cellField[x:x,y:y,z:z]=Apical
            # assign targetVolume        
            Apical.targetVolume = cell.targetVolume/6
            Basal.targetVolume = cell.targetVolume/6
            cell.targetVolume = cell.targetVolume*2/3
            Apical.lambdaVolume = cell.lambdaVolume
            Basal.lambdaVolume = cell.lambdaVolume
            
            # assign Cell cluster
            reassignIdFlag=self.inventory.reassignClusterId(Apical,cell.clusterId) # changing cluster id to 1536 for cell 'cell'
            reassignIdFlag=self.inventory.reassignClusterId(Basal,cell.clusterId) # changing cluster id to 1536 for cell 'cell'
            
    def updateAttributes(self):
        childCell = self.mitosisSteppable.childCell
        parentCell = self.mitosisSteppable.parentCell
        childCell.type = parentCell.type
#############################
class Growth(MitosisSteppableBase):
    def __init__(self,_simulator,_frequency):
        MitosisSteppableBase.__init__(self,_simulator,_frequency)
    def start(self):
        # Assign Property for Cell ID = 1
        cells_to_die=[]
        for cell in self.cellListByType(1):
            cell.targetVolume = tVol*random.uniform(0.75,1.25)# Make the initial target Volume of diff cells constant       
    def updateAttributes(self):
        childCell = self.mitosisSteppable.childCell
        parentCell = self.mitosisSteppable.parentCell
        childCell.type = parentCell.type

    def step(self,mcs):
        cells_to_die=[]
        global NoOfDeathThre,NoOfDeathBM
        NoOfDeathThre=0
        NoOfDeathBM=0
        for cell in self.cellList:
            if cell.type == 1:
                # Assume Growth only happens at the bottom of crypt
                if cell.targetVolume:
                    if mcs >50 and cell.yCOM < self.dim.y*0.5:
                    # Program Cell Growth
                        # access/modification of a dictionary attached to cell - make sure to decalare in main script that you will use such attribute
                        cellDict=self.getDictionaryAttribute(cell)
#                         if cellDict["G"] == True:
                        cell.targetVolume+= 1*random.uniform(0.5,1.25)
#                         else:
#                             cell.targetVolume = 1000    
                    # The diff cells remain unchanged    
#                     else:
#                         cell.targetVolume = tVol*random.uniform(1,1.25)
# Seperate cell death from cell growth 
                # Program Cell Death
                # Set up threshold to kill cells when cells go above the threshold
                if cell.yCOM > self.dim.y-10:
                    cells_to_die.append(cell)
                    NoOfDeathThre+=1
                    print "       ~~~~~~~~~~~~~~~~~~Dieing cell Thre~~~~~~~~~~~~~~~",cell.id,cell.type,cell.volume,cell.targetVolume
                cellNeighborList=self.getCellNeighbors(cell) # generates list of neighbors of cell 'cell'
                wallflag=0
                
                # Kill cells when the cells not touching BM
                compartmentList=self.inventory.getClusterCells(cell.clusterId)
                for cell2 in compartmentList:
                    cellNeighborList=self.getCellNeighbors(cell2) # generates list of neighbors of cell 'cell'
                    for neighbor in cellNeighborList:
                        # neighborSurfaceData.neighborAddress is a pointer to one of 'cell' neighbors stired in cellNeighborList
                        #IMPORTANT: cell may have Medium (NULL pointer) as a neighbor. therefore before accessing neighbor we first check if it is no Medium
                        if neighbor.neighborAddress: 
                            # Detect the cells touching BM
                            #print "       ~~~~~~~~~~~~~~~~~~NeighborCellType~~~~~~~~~~~~~~~",neighbor.neighborAddress.type
                            if neighbor.neighborAddress.type == 2:
                                wallflag=1
                if wallflag==0:
                    # Delete the cells without contacting BM
                    cells_to_die.append(cell)
                    if mcs > 0:
                        NoOfDeathBM+=1
                        print "       ~~~~~~~~~~~~~~~~~~Dieing cell BM~~~~~~~~~~~~~~~",cell.id,cell.type,cell.volume,cell.targetVolume
        # Cell Killing program
        for cell in cells_to_die:    
            compartmentList=self.inventory.getClusterCells(cell.clusterId)
            for cell2 in compartmentList:
                cell2.targetVolume = 0
                self.deleteCell(cell2)
            # The next line araises the dead cell register      
            self.cleanDeadCells()
#############################
class DeltaNotchClass(SteppableBasePy):
    def __init__(self,_simulator,_frequency):
        SteppableBasePy.__init__(self,_simulator,_frequency)        
    def start(self):
        # adding options that setup SBML solver integrator - these are optional but useful when encounteting integration instabilities              
        Name = "DeltaNotch"
        Key  = "DN"
        modelFile='Simulation/DN_Collier_1.sbml' 
        options={'relative':1e-10,'absolute':1e-12,'steps':10}
        self.setSBMLGlobalOptions(options)
        self.addSBMLToCellTypes(_modelFile=modelFile,_modelName="DN",_types=[self.TYPEA],_stepSize=0.2)  
        #Initial conditions
        state={} #dictionary to store state veriables of the SBML model
        for cell in self.cellListByType(1):
            state['D'] = random.uniform(0.2,1.0)
            state['N'] = random.uniform(0.2,1.0)
            state['B'] = random.uniform(0.9,1.0)
            state['R'] = random.uniform(0.9,1.0)
            self.setSBMLState(_modelName=Key,_cell=cell,_state=state)
            cellDict=self.getDictionaryAttribute(cell)
            cellDict['D']=state['D']
            cellDict['N']=state['N']
            cellDict['B']=state['B']
            cellDict['R']=state['R']
    def step(self,mcs):
        for cell in self.cellListByType(1):
            Davg=0.0; nn=0
            compartmentList=self.inventory.getClusterCells(cell.clusterId)
            neighborList = []
            for cell2 in compartmentList:
                for neighbor , commonSurfaceArea in self.getCellNeighborDataList(cell2):                
                    if neighbor:
                        #print "neighbor.id",neighbor.id,"neighbor.type=",neighbor.type,"clusterId=",neighbor.clusterId
                        neighborList.append(neighbor)
            #print "~~~~~~~~~~neighbor List~~~~~",neighborList
            for neighbor , commonSurfaceArea in self.getCellNeighborDataList(cell):
                if (neighbor and neighbor.type == 1):
                    nn+=1
                    state=self.getSBMLState(_modelName='DN',_cell=neighbor)
                    Davg+=state['D']   
            if (nn>0):
                Davg=Davg/nn
            
            state={}
            yCOM=(cell.yCOM/self.dim.y)
            
            if (cell.yCOM<self.dim.y*0.3):
                GammaB = 1
            else:
                GammaB = 10
            state['Davg']=Davg
            state['GammaB']=GammaB
            #print "~~~~~~~~~~~~~~~~~~~~~~~~~~~Davg:%f~~~~~~~~~~~~" %Davg
            #print "cell ID:%d~~~~~~~~~~~~~~~~~~~~~~~~~~~GammaB:%f~~~~~~~~~~~~" %(cell.id,GammaB)       
            self.setSBMLState(_modelName='DN',_cell=cell,_state=state)
            state=self.getSBMLState(_modelName='DN',_cell=cell)
            testB=self.getSBMLValue(_modelName='DN',_valueName='GammaB',_cell=cell)
            #print "cell ID:%d~~~~~~~~~~~~~~~~~~~~~~~~~~~TESTB:%f~~~~~~~~~~~~" %(cell.id,testB)                   
            cellDict=self.getDictionaryAttribute(cell)
            cellDict['D']=state['D']
            cellDict['N']=state['N']  
            cellDict['B']=state['B']
            cellDict['R']=state['R']
            
            #print "cell ID:%d~~~~~~~~~~~~~~~~~~~~~~~N:%f,\tD:%f,\tB:%f~~~~~~~~~~~~\n" %(cell.id,state['N'],state['D'],state['B'])
        self.timestepSBML()
#############################
class MitosisSteppable(MitosisSteppableBase):
    def __init__(self, _simulator, _frequency=1):
        MitosisSteppableBase.__init__(self, _simulator, _frequency)
    def step(self,mcs):
        global NoOfDivCells
        NoOfDivCells = 0
        cells_to_divide=[]
        
        for cell in self.cellListByType(1):
            cellDict=self.getDictionaryAttribute(cell)
            if mcs>50 and cell.yCOM<self.dim.y*0.3 and cell.volume > tVol*1.5 and cellDict['N']>0.5:
                print "~~~~~~~~~~~~~~~~Stem cell to divide~~~~~~~~~N:%f, y:%f" %(cellDict['N'],cell.yCOM)
                cells_to_divide.append(cell)
                NoOfDivCells+=1
            elif mcs>50 and cell.yCOM<self.dim.y*0.5 and cell.yCOM>self.dim.y*0.3 and cell.volume > tVol*1.75:   
                print "~~~~~~~~~~~~~~~~TA cell to divide~~~~~~~~~N:%f, y:%f" %(cellDict['N'],cell.yCOM)
                cells_to_divide.append(cell)
                NoOfDivCells+=1
                
        for cell in cells_to_divide:
            # break links between apical and basal before cell division
            for fpp in self.getInternalFocalPointPlasticityDataList(cell):
                cell2=fpp.neighborAddress
                self.focalPointPlasticityPlugin.deleteInternalFocalPointPlasticityLink(cell,cell2)
            self.divideCellAlongMajorAxis(cell)
            
    def updateAttributes(self):
        childCell = self.mitosisSteppable.childCell
        parentCell = self.mitosisSteppable.parentCell
        childCell.type = parentCell.type    
        parentCell.targetVolume = tVol
        childCell.targetVolume = tVol
        childCell.lambdaVolume = parentCell.lambdaVolume;
        # inherite properties from parent cells
        self.copySBMLs(_fromCell=parentCell,_toCell=childCell)
        childCellDict=CompuCell.getPyAttrib(childCell)
        parentCellDict=CompuCell.getPyAttrib(parentCell)
        childCellDict["D"]=random.uniform(0.9,1.0)*parentCellDict["D"]
        childCellDict["N"]=random.uniform(0.9,1.0)*parentCellDict["N"]
        childCellDict["B"]=random.uniform(0.9,1.0)*parentCellDict["B"]
        childCellDict["R"]=random.uniform(0.9,1.0)*parentCellDict["R"]
#############################
class ExtraFields(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.scalarFieldD=CompuCellSetup.createScalarFieldCellLevelPy("Delta")
        self.scalarFieldN=CompuCellSetup.createScalarFieldCellLevelPy("Notch")
        self.scalarFieldB=CompuCellSetup.createScalarFieldCellLevelPy("B-cat")
        self.scalarFieldR=CompuCellSetup.createScalarFieldCellLevelPy("NICD")   
    def step(self,mcs):     
        self.scalarFieldD.clear()
        self.scalarFieldN.clear()
        self.scalarFieldB.clear()
        self.scalarFieldR.clear()
        for cell in self.cellListByType(1):
            cellDict=CompuCell.getPyAttrib(cell)
            self.scalarFieldD[cell]=cellDict['D']
            self.scalarFieldN[cell]=cellDict['N']
            self.scalarFieldB[cell]=cellDict['B']
            self.scalarFieldR[cell]=cellDict['R']    
#############################
class Plotting1(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):
        # Making Scientific Plot
        self.pW=CompuCellSetup.viewManager.plotManager.getNewPlotWindow()
        if not self.pW:
            return
        #Plot Title - properties           
        self.pW.setTitle("No of Cell Death")
        self.pW.setTitleSize(12)
        self.pW.setTitleColor("Green") # you may choose different color - type its name  
        #plot background
        self.pW.setPlotBackgroundColor("white") # you may choose different color - type its name
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
        self.pW.addPlot("NoOfCellDeathThre",_style='Steps',_color='green',_size=5)
        self.pW.addPlot("NoOfCellDeathBM",_style='Steps',_color='black',_size=5)    
        self.pW.addGrid()
        #adding automatically generated legend
        # default position is at the bottom of the plot but here we put it at the top
        self.pW.addAutoLegend("top")      
        self.clearFlag=False
    def step(self,mcs):     
        global NoOfDeathBM,NoOfDeathThre
        #self.pW.eraseAllData() # this is how you erase previous content of the plot
        self.pW.addDataPoint("NoOfCellDeathThre",mcs,NoOfDeathThre) # arguments are (name of the data series, x, y)
        self.pW.addDataPoint("NoOfCellDeathBM",mcs,NoOfDeathBM) # arguments are (name of the data series, x, y)
        self.pW.showAllPlots() 
#############################        
class Plotting2(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):
        # Making Scientific Plot
        self.pW=CompuCellSetup.viewManager.plotManager.getNewPlotWindow()
        if not self.pW:
            return
        #Plot Title - properties           
        self.pW.setTitle("No of Dividing Cells")
        self.pW.setTitleSize(12)
        self.pW.setTitleColor("Green") # you may choose different color - type its name  
        #plot background
        self.pW.setPlotBackgroundColor("white") # you may choose different color - type its name
        # properties of x axis
        self.pW.setXAxisTitle("No Of MCS")
        self.pW.setXAxisTitleSize(10)      
        self.pW.setXAxisTitleColor("blue")  # you may choose different color - type its name            
        # properties of y axis
        self.pW.setYAxisTitle("No of Cells")        
        #self.pW.setYAxisLogScale()
        self.pW.setYAxisTitleSize(10)        
        self.pW.setYAxisTitleColor("red")  # you may choose different color - type its name                                
        # choices for style are NoCurve,Lines,Sticks,Steps,Dots
        self.pW.addPlot("NoOfDivCells",_style='Steps',_color='red',_size=5)
#         self.pW.addPlot("NoOfCellDeathBM",_style='Dots')    
        self.pW.addGrid()
        #adding automatically generated legend
        # default position is at the bottom of the plot but here we put it at the top
        self.pW.addAutoLegend("top")      
        self.clearFlag=False
    def step(self,mcs):     
        global NoOfDivCells
        #self.pW.eraseAllData() # this is how you erase previous content of the plot
        self.pW.addDataPoint("NoOfDivCells",mcs,NoOfDivCells) # arguments are (name of the data series, x, y)
#             self.pW.addDataPoint("NoOfCellDeathBM",mcs,NoOfDeathBM) # arguments are (name of the data series, x, y)
        self.pW.showAllPlots() 
#############################
class OutputData(SteppableBasePy):
    def __init__(self,_simulator,_frequency):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):
#         self.createFile(0)
#     def createFile(self):
        # output data
        simulationDir=os.path.dirname (os.path.abspath( __file__ ))
        FilePath= os.path.join(simulationDir,'IntegrationStep.txt')
        self.txt = open(FilePath,"w") # create a new file
        textline = "mcs\tNoOfDeathBM\tNoOfDeathThre\tNoOfDivCells\n"
        ## you can access/manipulate cell properties here 
        self.txt.write(textline) 
        print textline
    def step(self,mcs):
        global NoOfDeathBM,NoOfDeathThre,NoOfDivCells
        textline = "%4d\t%2d\t%2d\t%2d\n" %(mcs,NoOfDeathBM,NoOfDeathThre,NoOfDivCells)
        ## you can access/manipulate cell properties here 
#         self.createFile(1)
        self.txt.write(textline)
        print textline