/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package Astrocytes_Tools;

import ij.ImagePlus;
import ij.measure.Calibration;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.image3d.ImageHandler;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij2.CLIJ2;

/**
 *
 * @author phm
 */
public class Clij_Proccesing {
    private final CLIJ2 clij2 = CLIJ2.getInstance();
    
    /**
     * Median filter using CLIJ2
     * @param img
     * @param sizeXY
     * @param sizeZ
     * @return 
     */ 
    public ImagePlus median_filter(ImagePlus img, double sizeXY, double sizeZ) {
       ClearCLBuffer imgCL = clij2.push(img);
       ClearCLBuffer imgCLMed = clij2.create(imgCL);
       clij2.median3DBox(imgCL, imgCLMed, sizeXY, sizeXY, sizeZ);
       ImagePlus imgFilter = clij2.pull(imgCLMed);
       clij2.release(imgCL);
       return(imgFilter);
    }  
    
    /**
     * Test if GPU
     * @return 
     */
    public boolean isGPU() {
        String gpuName = clij2.getGPUName();
        return(gpuName.isEmpty());
    }
    
     /**
     * Difference of Gaussians 
     * Using CLIJ2
     * @param img
     * @param size1
     * @param size2
     * @return imgGauss
     */ 
    public ImagePlus DOG(ImagePlus img, double size1, double size2) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLDOG = clij2.create(imgCL);
        clij2.differenceOfGaussian3D(imgCL, imgCLDOG, size1, size1, size1, size2, size2, size2);
        ImagePlus imgFilter = clij2.pull(imgCLDOG);
        return(imgFilter);
    }
    
    /**
     * Threshold 
     * USING CLIJ2
     * @param img
     * @param thMed
     */
    public ImagePlus threshold(ImagePlus img, String thMed) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLBin = clij2.create(imgCL);
        clij2.automaticThreshold(imgCL, imgCLBin, thMed);
        ImagePlus imgBin = clij2.pull(imgCLBin);
        return(imgBin);
    }
    
    
}
