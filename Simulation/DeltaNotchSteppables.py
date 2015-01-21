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

tVol = 1000
# Calcuate Cell Death due to Threshold and BM touching
NoOfDeathBM = 0
NoOfDeathThre = 0
# threshold for simulation
Thr_mcs = 100
Thr_TA = 0.7
Thr_SC = 0.3
Thr_N  = 0.5
GR=2.0
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
            
        radi = int(self.dim.x/2)
        # central tube generation
        pt=CompuCell.Point3D(0,0,0)
        Tube = self.potts.createCellG(pt) # the arguments are (type,pt,xSize,ySize,zSize)
        Tube.type = 3
        totalNumOfSeg = 1
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
                        if ((ptx-radi)**2+(ptz-radi)**2)<=((radi-5)/2)**2:
                        # central lumen
                            overwrittenCell=self.cellField.get(pt)
                            self.cellField.set(pt,Tube)
                            self.cleanDeadCells()
                    # Generate Semi-Sphere Sturcture ((Bottom of the crypt)
                    elif pty < radi:
                        if ((ptx-radi)**2+(pty-radi)**2+(ptz-radi)**2)<=((radi-5)/2)**2:
                            self.cellField.set(pt,Tube)
                            self.cleanDeadCells()
        # Assign property for Cell type = 3
        for cell in self.cellListByType(3): 
            cell.targetVolume = cell.volume
            cell.lambdaVolume = 15000

    # Basal Membrane Generation
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
                        if ((ptx-radi)**2+(ptz-radi)**2)>=(radi-5)**2:
                            overwrittenCell=self.cellField.get(pt)
                            self.cellField.set(pt,Wall)
                            self.cleanDeadCells()
                            ang = math.atan2(ptx-radi,ptz-radi)*(180./math.pi)+180
                            segN= int(ang/(360/totalNumOfSeg))
                            if segN ==15:
                                segN = 0
                            segm.append([segN,pt])

                    # Generate Semi-Sphere Sturcture ((Bottom of the crypt)
                    elif pty < radi:
                        if ((ptx-radi)**2+(pty-radi)**2+(ptz-radi)**2)>=(radi-5)**2:
                            self.cellField.set(pt,Wall)
                            self.cleanDeadCells()
                                
        # now we can overwrite (0,0,0) with medium
        pt.x=0
        pt.y=0
        pt.z=0
        self.cellField.set(pt,CompuCell.getMediumCell())

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
            # Assign stochastic initial conditions for stem cell/TA cell volume
            cell.targetVolume = tVol*random.uniform(0.95,1.5)# Make the initial target Volume of diff cells constant  
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
                cellDict=CompuCell.getPyAttrib(cell)
                if cell.targetVolume:
                    if mcs >Thr_mcs and cell.yCOM < self.dim.y*Thr_TA and cell.yCOM > self.dim.y*Thr_SC:
                        # Program TA Cell Growth
                        cell.targetVolume+= GR*random.uniform(0.85,1.25)#random growth rate
                    elif mcs >Thr_mcs and cell.yCOM < self.dim.y*Thr_SC and cellDict['N']>Thr_N:
                        # Program Stem Cell Growth
                        cell.targetVolume+= GR*random.uniform(0.85,1.25)#random growth rate   
                    elif cell.volume<=tVol*1.5:
                        cell.targetVolume+= GR*random.uniform(0.85,1.25)#random growth rate   
                    else:
                        # The diff cells remain unchanged 
                        cell.targetVolume = tVol*1.5
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
                    if mcs > 0:
                        NoOfDeathBM+=1
                        print "       ~~~~~~~~~~~~~~~~~~Dieing cell BM~~~~~~~~~~~~~~~",cell.id,cell.type,cell.volume,cell.targetVolume
        # Cell Killing program
        for cell in cells_to_die:    
            cell.targetVolume = 0
            self.deleteCell(cell)
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
        modelFile='Simulation/PFLIMI.sbml' 
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
            
#             state['Davg']=0.5
#             state['Navg']=0.5
#             state['BN0'] = 0.0
#             state['BN'] = 1.0
#             state['p'] = 3.0
#             state['KN'] = 0.1
#             state['AN'] = 1.0
#             state['Kc'] = 1.0
#             state['Kt'] = 1.0
#             state['BD0'] = 0.0
#             state['BD'] = 1.0
#             state['KD'] = 0.01
#             state['h'] = 3.0
#             state['AD'] = 1.0
#             state['BetaB'] = 1.0
#             state['GammaB']=1 
#             state['AR'] = 1.0
            
            self.setSBMLState(_modelName=Key,_cell=cell,_state=state)
            cellDict=self.getDictionaryAttribute(cell)
            cellDict['D']=state['D']
            cellDict['N']=state['N']
            cellDict['B']=state['B']
            cellDict['R']=state['R']
    def step(self,mcs):
        for cell in self.cellListByType(1):
            Davg=0.0; nn=0; Navg=0.0
            compartmentList=self.inventory.getClusterCells(cell.clusterId)
            neighborList = []
            for cell2 in compartmentList:
                for neighbor , commonSurfaceArea in self.getCellNeighborDataList(cell2):                
                    if neighbor:
                        neighborList.append(neighbor)
            for neighbor , commonSurfaceArea in self.getCellNeighborDataList(cell):
                if (neighbor and neighbor.type == 1):
                    nn+=1
                    state=self.getSBMLState(_modelName='DN',_cell=neighbor)
                    Davg+=state['D']
                    Navg+=state['N']
            if (nn>0):
                Davg=Davg/nn
                Navg=Navg/nn
            state={}
            if (cell.yCOM<self.dim.y*Thr_SC):
                GammaB = 1.0
            else:
                GammaB = 10.0
            state['Davg']=Davg
            state['Navg']=Navg 
            state['GammaB']=GammaB 
            self.setSBMLState(_modelName='DN',_cell=cell,_state=state)
            state=self.getSBMLState(_modelName='DN',_cell=cell)
#             testB=self.getSBMLValue(_modelName='DN',_valueName='GammaB',_cell=cell)                
            cellDict=self.getDictionaryAttribute(cell)
            cellDict['D']=state['D']
            cellDict['N']=state['N']  
            cellDict['B']=state['B']
            cellDict['R']=state['R']
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
            if mcs>Thr_mcs and cell.yCOM<self.dim.y*Thr_SC and cell.volume > tVol*1.5 and cellDict['N']>Thr_N:
                # programmed stem cell division
                print "~~~~~~~~~~~~~~~~Stem cell to divide~~~~~~~~~N:%f, y:%f" %(cellDict['N'],cell.yCOM)
                cells_to_divide.append(cell)
                NoOfDivCells+=1
            elif mcs>Thr_mcs and cell.yCOM<self.dim.y*Thr_TA and cell.yCOM>self.dim.y*Thr_SC and cell.volume > tVol*2:
                # programmed TA cell division
                print "~~~~~~~~~~~~~~~~TA cell to divide~~~~~~~~~y:%f" %(cell.yCOM)
                cells_to_divide.append(cell)
                NoOfDivCells+=1
        for cell in cells_to_divide:
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
        childCellDict["D"]=parentCellDict["D"]
        childCellDict["N"]=parentCellDict["N"]
        childCellDict["B"]=parentCellDict["B"]
        childCellDict["R"]=parentCellDict["R"]
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