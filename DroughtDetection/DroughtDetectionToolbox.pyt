# Copyright (C) 2019, Caleb Davis, University of North Georgia

import arcpy
import numpy as np
import math


class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "FloodToolbox"
        self.alias = "FloodToolbox"

        # List of tool classes associated with this toolbox
        self.tools = [FloodDetection]


class FloodDetection(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Flood Detection"
        self.description = "Given two rasters and sample points, this tool will find the difference in area."
        self.canRunInBackground = False

    def getParameterInfo(self):
        """Define parameter definitions"""
        bef_ras = arcpy.Parameter(
                name="bef_ras",
                displayName="Before Flooding Raster",
                direction="Input",
                parameterType="Required",
                datatype="GPRasterLayer",
                )
        aft_ras = arcpy.Parameter(
                name="aft_ras",
                displayName="After Flooding Raster",
                direction="Input",
                parameterType="Required",
                datatype="GPRasterLayer",
                )
        bef_r_avg = arcpy.Parameter(
                name="bef_r_avg",
                displayName="Input average of red values from before raster",
                direction="Input",
                parameterType="Required",
                datatype="GPDouble",
                )
        bef_g_avg = arcpy.Parameter(
                name="bef_g_avg",
                displayName="Input average of green values from before raster",
                direction="Input",
                parameterType="Required",
                datatype="GPDouble",
                )
        bef_b_avg = arcpy.Parameter(
                name="bef_b_avg",
                displayName="Input average of blue values from before raster",
                direction="Input",
                parameterType="Required",
                datatype="GPDouble",
                )
        aft_r_avg = arcpy.Parameter(
                name="aft_r_avg",
                displayName="Input average of red values from after raster",
                direction="Input",
                parameterType="Required",
                datatype="GPDouble",
                )
        aft_g_avg = arcpy.Parameter(
                name="aft_g_avg",
                displayName="Input average of green values from after raster",
                direction="Input",
                parameterType="Required",
                datatype="GPDouble",
                )
        aft_b_avg = arcpy.Parameter(
                name="aft_b_avg",
                displayName="Input average of blue values from after raster",
                direction="Input",
                parameterType="Required",
                datatype="GPDouble",
                )
        bef_r_std = arcpy.Parameter(
                name="bef_r_std",
                displayName="Input standard deviation of red values from before raster",
                direction="Input",
                parameterType="Required",
                datatype="GPDouble",
                )
        bef_g_std = arcpy.Parameter(
                name="bef_g_std",
                displayName="Input standard deviation of green values from before raster",
                direction="Input",
                parameterType="Required",
                datatype="GPDouble",
                )
        bef_b_std = arcpy.Parameter(
                name="bef_b_std",
                displayName="Input standard deviation of blue values from before raster",
                direction="Input",
                parameterType="Required",
                datatype="GPDouble",
                )
        aft_r_std = arcpy.Parameter(
                name="aft_r_std",
                displayName="Input standard deviation of red values from after raster",
                direction="Input",
                parameterType="Required",
                datatype="GPDouble",
                )
        aft_g_std = arcpy.Parameter(
                name="aft_g_std",
                displayName="Input standard deviation of green values from after raster",
                direction="Input",
                parameterType="Required",
                datatype="GPDouble",
                )
        aft_b_std = arcpy.Parameter(
                name="aft_b_std",
                displayName="Input standard deviation of blue values from after raster",
                direction="Input",
                parameterType="Required",
                datatype="GPDouble",
                )
        k_val = arcpy.Parameter(
                name="k_val",
                displayName="Input desired k-value, use 1.0 as default",
                direction="Input",
                parameterType="Required",
                datatype="GPDouble",
                )
        output_path = arcpy.Parameter(
                datatype="GPString",
                name="output_path",
                displayName="Output Folder Path",
                direction="Input",
                parameterType="Required",
                )
        output_prefix = arcpy.Parameter(
                datatype="GPString",
                name="output_prefix",
                displayName="Output Prefix",
                direction="Input",
                parameterType="Required",
                )
        params = [bef_ras, aft_ras, bef_r_avg, bef_g_avg, bef_b_avg, aft_r_avg, aft_g_avg, aft_b_avg, 
        bef_r_std, bef_g_std, bef_b_std, aft_r_std, aft_g_std, aft_b_std, k_val, output_path, output_prefix]
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
        input_ras_bef = parameters[0].valueAsText
        input_ras_aft = parameters[1].valueAsText
        bef_r_avg = parameters[2].value
        bef_g_avg = parameters[3].value
        bef_b_avg = parameters[4].value
        aft_r_avg = parameters[5].value
        aft_g_avg = parameters[6].value
        aft_b_avg = parameters[7].value
        bef_r_std = parameters[8].value
        bef_g_std = parameters[9].value
        bef_b_std = parameters[10].value
        aft_r_std = parameters[11].value
        aft_g_std = parameters[12].value
        aft_b_std = parameters[13].value
        k_val = parameters[14].value
        output_path = parameters[15].valueAsText
        output_prefix = parameters[16].valueAsText
        

        ras = arcpy.Raster(input_ras_bef)
        ras_a_bef = arcpy.RasterToNumPyArray(ras)
        
        ras = arcpy.Raster(input_ras_aft)
        ras_a_aft = arcpy.RasterToNumPyArray(ras)
        
        dist_a_bef, ret_a_bef = find_segments(ras_a_bef, [bef_r_avg, bef_g_avg, bef_b_avg], [bef_r_std, bef_g_std, bef_b_std], k_val)
        dist_before = arcpy.NumPyArrayToRaster(dist_a_bef, ras.extent.lowerLeft, ras.meanCellWidth, ras.meanCellHeight)
        ret_before = arcpy.NumPyArrayToRaster(ret_a_bef, ras.extent.lowerLeft, ras.meanCellWidth, ras.meanCellHeight)
        
        dist_a_aft, ret_a_aft = find_segments(ras_a_aft, [aft_r_avg, aft_g_avg, aft_b_avg], [aft_r_std, aft_g_std, aft_b_std], k_val)
        dist_after = arcpy.NumPyArrayToRaster(dist_a_aft, ras.extent.lowerLeft, ras.meanCellWidth, ras.meanCellHeight)
        ret_after = arcpy.NumPyArrayToRaster(ret_a_aft, ras.extent.lowerLeft, ras.meanCellWidth, ras.meanCellHeight)
        
        bin_a_aft = np.where(dist_a_aft>= 0,1,0)
        bin_a_bef = np.where(dist_a_bef>= 0,1,0)
        
        arcpy.AddMessage(bin_a_aft.shape[0])
        arcpy.AddMessage(bin_a_aft.shape[1])
        arcpy.AddMessage(bin_a_bef.shape[0])
        arcpy.AddMessage(bin_a_bef.shape[1])
        
        affected = np.where(np.logical_and(bin_a_aft==1, bin_a_bef==0),1,0)
        affected_r = arcpy.NumPyArrayToRaster(affected, ras.extent.lowerLeft, ras.meanCellWidth, ras.meanCellHeight, 0)
        affected_r.save(output_path + output_prefix + '.tif')
        
        return

def find_segments(ras_a, a, sd, k):

    
    # get useful information from ras_a.shape
    nbands = ras_a.shape[0]
    nbands = 3
    nrows = ras_a.shape[1]
    ncols = ras_a.shape[2]
    
    # calculate the radius of the target color sphere
    radius = 0
    for b in range(nbands):
        radius += (k*sd[b])**2
    radius = math.sqrt(radius)

    # create a new zero distance array
    dist_a = np.zeros(ras_a.shape[1:3])
    
    # clone the original raster array
    ret_a = ras_a.copy()
    
    # iterations
    for r in range(nrows):
        for c in range(ncols):
            # calculate the color distance for this cell at (r, c)
            for b in range(nbands):
                dist_a[r,c] += (ras_a[b,r,c]-a[b])**2
            
            # square root to get the distance, not sum of squared distances
            dist_a[r,c] = math.sqrt(dist_a[r,c])
            
            # if this cell is outside the sphere
            if dist_a[r,c] > radius:
                # return a negative distance
                dist_a[r,c] = -dist_a[r,c]
                
                # return black
                for b in range(nbands):
                    ret_a[b,r,c] = 0
    
    return dist_a, ret_a

  
