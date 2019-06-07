# Copyright (C) 2019, Anthony Sonsteng, University of North Georgia

import numpy as np
import math
import statistics

class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Toolbox"
        self.alias = ""

        # List of tool classes associated with this toolbox
        self.tools = [Automatic_greenspace_detector]


class Automatic_greenspace_detector(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Automatic Greenspace detector"
        self.description = ""
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        input_raster = arcpy.Parameter(
            name="facc",
            displayName="input_raster",
            direction="Input",
            datatype="GPRasterLayer",
            parameterType="Required",
        )		
        R_mean = arcpy.Parameter(
            name="R_mean",
            displayName="R_mean",
            direction="Input",
            datatype="GPDouble",
            parameterType="Required",            
        )
        G_mean = arcpy.Parameter(
            name="G_mean",
            displayName="G_mean",
            direction="Input",
            datatype="GPDouble",
            parameterType="Required",            
        )
        B_mean = arcpy.Parameter(
            name="B_mean",
            displayName="B_mean",
            direction="Input",
            datatype="GPDouble",
            parameterType="Required",            
        )
        R_stdev = arcpy.Parameter(
            name="R_stdev",
            displayName="R_stdev",
            direction="Input",
            datatype="GPDouble",
            parameterType="Required",            
        )
        G_stdev = arcpy.Parameter(
            name="G_stdev",
            displayName="G_stdev",
            direction="Input",
            datatype="GPDouble",
            parameterType="Required",            
        )
        B_stdev = arcpy.Parameter(
            name="B_stdev",
            displayName="B_stdev",
            direction="Input",
            datatype="GPDouble",
            parameterType="Required",            
        )
        K = arcpy.Parameter(
            name="K",
            displayName="K",
            direction="Input",
            datatype="GPDouble",
            parameterType="Required",            
        )
        output_filepath_1 = arcpy.Parameter(
            name="output_filepath_1",
            displayName="Output_filepath_1",
            direction="Output",
            datatype="DERasterDataset",
            parameterType="Required",
        )
        output_filepath_2 = arcpy.Parameter(
            name="output_filepath_2",
            displayName="Output_filepath_2",
            direction="Output",
            datatype="DERasterDataset",
            parameterType="Required",
        )
		
        params = [input_raster, R_mean, G_mean, B_mean, R_stdev, G_stdev, B_stdev, K, output_filepath_1, output_filepath_2]
        return params

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        """Modify the values and properties of parameters before internal
        validation is performed.  This method is called whenever a parameter
        has been changed."""
        return

    def updateMessages(self, parameters):
        """Modify the messages created by internal validation for each tool
        parameter.  This method is called after internal validation."""
        return

    def execute(self, parameters, messages):
        """The source code of the tool."""
		
        input_ras = parameters[0].valueAsText
        R_mean = parameters[1].value
        G_mean = parameters[2].value
        B_mean = parameters[3].value
        R_stdev = parameters[4].value
        G_stdev = parameters[5].value
        B_stdev = parameters[6].value
        k = parameters[7].value
        output_filepath_1 = parameters[8].valueAsText
        output_filepath_2 = parameters[9].valueAsText
		
        ras = arcpy.Raster(input_ras)
        ras_a = arcpy.RasterToNumPyArray(ras)
        
        a = [R_mean, G_mean, B_mean]
        sd = [R_stdev, G_stdev, B_stdev]
		
        dist_a, ret_a = find_segments(ras_a, a, sd, k)
        dist = arcpy.NumPyArrayToRaster(dist_a, ras.extent.lowerLeft, ras.meanCellWidth, ras.meanCellHeight)
        ret = arcpy.NumPyArrayToRaster(ret_a, ras.extent.lowerLeft, ras.meanCellWidth, ras.meanCellHeight)
        
        dist.save(output_filepath_1)
        ret.save(output_filepath_2)
        return

def find_segments(ras_a, a, sd, k):
    #get useful info from ras_a.shp
    nbands = 3
    nrows = ras_a.shape[1]
    ncols = ras_a.shape[2]
    
    #calculate radius of target color sphere
    radius = 0
    for b in range(nbands):    
        radius += (k*sd[b])**2
    radius = math.sqrt(radius)
    
    #create a new zero distance array
    dist_a = np.zeros(ras_a.shape[1:3])

    #clone the original raster array
    ret_a = ras_a.copy()
    
    #iterations
    for r in range(nrows):
        for c in range(ncols):
            #calculate color distance for this call aat (r, c)
            for b in range(nbands):
                dist_a[r,c] +=(ras_a[b,r,c]-a[b])**2
                
            #sqrt to get the distance, notsum on squared distance
            dist_a[r,c] = math.sqrt(dist_a[r,c])
            
            #if this cell is outside the sphere
            if dist_a[r,c] > radius:
                #return a negative distance
                dist_a[r,c]= -dist_a[r,c]
                
                #return black
                for b in range(nbands):
                    ret_a[b,r,c] = 0

    return dist_a, ret_a  
