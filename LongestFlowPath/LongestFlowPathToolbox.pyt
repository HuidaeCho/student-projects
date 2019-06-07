# Copyright (C) 2019, Stefen Gray, University of North Georgia

import arcpy
import math
from arcpy import env


class Toolbox(object):
    def __init__(self):
        '''Define the toolbox (the name of the toolbox is the name of the
        .pyt file).'''
        self.label = 'HelloWorld Python Toolbox'
        self.alias = 'HelloWorldPythonToolbox'

        # List of tool classes associated with this toolbox
        self.tools = [LongestFlowPathProject]


class LongestFlowPathProject(object):
    def __init__(self):
        """
        define the tool (Tool name is the name of the class)
        """
        self.label = "Longest Flow Path"
        self.description = "this is my Longest Flow Path python tool"

    def getParameterInfo(self):
        '''Define parameter definitions'''
        fdir = arcpy.Parameter(
                displayName='Flow Direction',
                name='fdir',
                datatype='DERasterDataset',
                parameterType='Required',
                direction='Input')
        facc = arcpy.Parameter(
                displayName='Flow Accumulation',
                name='facc',
                datatype='DERasterDataset',
                parameterType='Required',
                direction='Input')

        outlets = arcpy.Parameter(
                displayName='Outlets',
                name='outlets',
                datatype='GPFeatureRecordSetLayer',
                parameterType='Required',
                direction='Input')
        paths = arcpy.Parameter(
                displayName='Paths',
                name='paths',
                datatype='DEFeatureClass',
                parameterType='Required',
                direction='Output')

        params = [fdir,facc,outlets,paths]
        return params

    def isLicensed(self):
        '''Set whether tool is licensed to execute.'''
        return True

    def updateParameters(self, parameters):
        '''Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed.'''
        return

    def updateMessages(self, parameters):
        '''Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation.'''
        return

    def execute(self, parameters, messages):
        '''The source code of the tool.'''
        fdir = parameters[0].valueAsText
        facc = parameters[1].valueAsText
        outlets = parameters[2].value
        path = parameters[3].value
        cursor = arcpy.da.SearchCursor(outlets,"Shape")
        outlets = []
        for row in cursor:
            outlets.append(row[0])
        features = []

        accrast = arcpy.Raster(facc)
        dirrast = arcpy.Raster(fdir)
        extent = dirrast.extent
        UR = (extent.XMin,extent.YMax)
        LL = (extent.XMax,extent.YMax)
        cellsize = dirrast.meanCellWidth
        polylines = arcpy.Array()
        spatialRef = dirrast.spatialReference

        D8toarray={1: [16, 8, 32, 4, 64, 2, 128],
                   2: [32, 16, 64, 8, 128, 4, 1],
                   4: [64, 32, 128, 16, 1, 8, 2],
                   8: [128, 64, 1, 32, 2, 16, 4],
                   16: [1, 128, 2, 64, 4, 32, 8],
                   32: [2, 1, 4, 128, 8, 64, 16],
                   64: [4, 2, 8, 1, 16, 128, 32],
                   128: [8, 4, 16, 2, 32, 1, 64]}
        D8inv={1: 16,
               2: 32,
               4: 64,
               8: 128,
               16: 1,
               32: 2,
               64: 4,
               128: 8}
        D8torowcol={1: (0, 1),
                    2: (1, 1),
                    4: (1, 0),
                    8: (1, -1),
                    16: (0, -1),
                    32: (-1, -1),
                    64: (-1, 0),
                    128: (-1, 1)}


        def RowColtoCoords(cellsize,UR,LL,rowcol):
            URx = UR[0] + cellsize/2
            URy = UR[1] - cellsize/2
            LLx = LL[0] - cellsize/2
            LLy = LL[1] + cellsize/2

            rowx = URx + (rowcol[1])*cellsize
            rowy = URy - (rowcol[0])*cellsize
            return (rowx,rowy)

        def CoordstoRowCol(cellsize,UR,LL,coords):
            URx = UR[0] + cellsize
            URy = UR[1] - cellsize
            LLx = LL[0] - cellsize
            LLy = LL[1] + cellsize

            x=coords[0]
            y=coords[1]

            row=math.floor((URy-y) / cellsize)
            col=math.floor((x-URx) / cellsize)
            if (x-URx) % cellsize >= 0.5:
                col = col+1
            if  (URy-y) % cellsize >= 0.5:
                row = row+1

            return row,col


        def PathAlg(startrow,startcol,rast,accrast,pArray):
            finished = 0

            row = startrow
            col = startcol
            extent = rast.extent
            UR = (extent.XMin, extent.YMax)
            LL = (extent.XMax, extent.YMin)
            xy = RowColtoCoords(rast.meanCellWidth, UR, LL, (row, col))
            rastarray = arcpy.RasterToNumPyArray(rast)
            accarray = arcpy.RasterToNumPyArray(accrast)
            pArray.add(arcpy.Point(xy[0],xy[1]))
            while finished == 0:
                done=0
                cellValue = rastarray[row][col]
                searchArray = D8toarray[cellValue]
                directions = []
                for dir in searchArray:
                    searchRowCol = D8torowcol[dir]
                    dr = searchRowCol[0]
                    dc = searchRowCol[1]
                    if col+dc < rast.width and row+dr < rast.height:
                        if rastarray[row+dr][col+dc] == D8inv[dir]:
                            directions.append((row+dr,col+dc))
                if len(directions)==1:
                    done=1
                    xy = RowColtoCoords(cellsize,UR,LL,(directions[0][0],directions[0][1]))
                    pArray.append(arcpy.Point(xy[0],xy[1]))
                    row=directions[0][0]
                    col=directions[0][1]
                else:
                    i=0
                    chosen=10
                    max=0
                    for rowcol in directions:
                        if rowcol[0] < accrast.height and rowcol[1] < accrast.width:
                            if accarray[rowcol[0]][rowcol[1]] >= max:
                                max = accarray[rowcol[0]][rowcol[1]]
                                done=1
                                chosen = i
                        i=i+1
                    if chosen != 10:
                        xy = RowColtoCoords(cellsize,UR,LL,(directions[chosen][0],directions[chosen][1]))
                        row=directions[chosen][0]
                        col=directions[chosen][1]
                        pArray.append(arcpy.Point(xy[0],xy[1]))
                if done != 1:
                    finished = 1
            return pArray

        for outlet in outlets:
            point = arcpy.Point(outlet[0],outlet[1])
            pointArray = arcpy.Array()
            rowcol = CoordstoRowCol(cellsize,UR,LL,(outlet[0],outlet[1]))
            pointsArray = PathAlg(rowcol[0],rowcol[1],dirrast,accrast,pointArray)
            #cursor=arcpy.da.InsertCursor(path,["SHAPE@"])
            features.append(arcpy.Polyline(pointsArray,spatialRef))
            #cursor.insertRow([polyline])
            aprx = arcpy.mp.ArcGISProject("CURRENT")
            aprxMap = aprx.listMaps("Map")[0]
            arcpy.CopyFeatures_management(features,path)
