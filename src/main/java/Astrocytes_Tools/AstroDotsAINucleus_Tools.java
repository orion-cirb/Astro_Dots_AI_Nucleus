package Astrocytes_Tools;

import AstroDots_StardistOrion.StarDist2D;
import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.measure.*;
import ij.plugin.Duplicator;
import ij.plugin.RGBStackMerge;
import ij.plugin.Thresholder;
import ij.plugin.ZProjector;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import java.awt.Color;
import java.awt.Font;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import javax.swing.ImageIcon;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom2.BoundingBox;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.Objects3DIntPopulationComputation;
import mcib3d.geom2.measurements.Measure2Distance;
import mcib3d.geom2.measurements.MeasureCentroid;
import mcib3d.geom2.measurements.MeasureIntensity;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.image3d.ImageFloat;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import mcib3d.image3d.distanceMap3d.EDT;
import mcib3d.utils.ThreadUtil;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;


 /*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author phm
 */



public class AstroDotsAINucleus_Tools {
    
public  double meanSEMDotsSize = 0;
public String thMethod = "";
public String segMethod = "";
public double bg = 0;
public double stdBg = 0;
// if pureArn exclude no dots
public boolean pureArn = false;
// ratio astro int / dots int for nucleus selection
public static boolean ratioInt = false; 

public double minNucSize = 50;
public double maxNucSize = 1000;
public double minDotsSize = 0.03;
public double maxDotsSize = 20;

private double minDOG = 2;
private double maxDOG = 4;

// Stardist
private final File modelsPath = new File(IJ.getDirectory("imagej")+File.separator+"models");
private String stardistModel = "StandardFluo.zip";
private double stardistProbThresh = 0.6;
private double stardistOverlayThresh = 0.25;

public Calibration cal = new Calibration();

private final Clij_Proccesing clij2 = new Clij_Proccesing();;
private final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));

 
 /*
    Find starDist models in Fiji models folder
    */
    public String[] findStardistModels() {
        FilenameFilter filter = (dir, name) -> name.endsWith(".zip");
        File[] modelList = modelsPath.listFiles(filter);
        String[] models = new String[modelList.length];
        for (int i = 0; i < modelList.length; i++) {
            models[i] = modelList[i].getName();
        }
        Arrays.sort(models);
        return(models);
    }  

/**
     * Find images extension
     * @param imagesFolder
     * @return 
     */
    public String findImageType(File imagesFolder) {
        String ext = "";
        File[] files = imagesFolder.listFiles();
        for (File file: files) {
            if(file.isFile()) {
                String fileExt = FilenameUtils.getExtension(file.getName());
                switch (fileExt) {
                   case "nd" :
                       ext = fileExt;
                       break;
                   case "nd2" :
                       ext = fileExt;
                       break;
                    case "czi" :
                       ext = fileExt;
                       break;
                    case "lif"  :
                        ext = fileExt;
                        break;
                    case "ics" :
                        ext = fileExt;
                        break;
                    case "ics2" :
                        ext = fileExt;
                        break;
                    case "tif" :
                        ext = fileExt;
                        break;
                    case "tiff" :
                        ext = fileExt;
                        break;
                }
            } else if (file.isDirectory() && !file.getName().equals("Results")) {
                ext = findImageType(file);
                if (! ext.equals(""))
                    break;
            }
        }
        return(ext);
    }
     
    
    
    /**
     * Find images in folder
     * @param imagesFolder
     * @param imageExt
     * @return 
     */
    public ArrayList<String> findImages(String imagesFolder, String imageExt) {
        File inDir = new File(imagesFolder);
        String[] files = inDir.list();
        if (files == null) {
            System.out.println("No image found in "+imagesFolder);
            return null;
        }
        ArrayList<String> images = new ArrayList();
        for (String f : files) {
            // Find images with extension
            String fileExt = FilenameUtils.getExtension(f);
            if (fileExt.equals(imageExt))
                images.add(imagesFolder + File.separator + f);
        }
        Collections.sort(images);
        return(images);
    }
    
    /**
     * Find image calibration
     * @param meta
     * @return 
     */
    public Calibration findImageCalib(IMetadata meta) {
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
        cal.pixelHeight = cal.pixelWidth;
        if (meta.getPixelsPhysicalSizeZ(0) != null)
        cal.pixelDepth = (meta.getPixelsPhysicalSizeZ(0) == null) ? cal.pixelDepth = 1 : meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
        cal.setUnit("microns");
        System.out.println("XY calibration = " + cal.pixelWidth + ", Z calibration = " + cal.pixelDepth);
        return(cal);
    }
    
    
     /**
     * Find channels name
     * @param imageName
     * @param meta
     * @param reader
     * @return 
     * @throws loci.common.services.DependencyException
     * @throws loci.common.services.ServiceException
     * @throws loci.formats.FormatException
     * @throws java.io.IOException
     */
    public String[] findChannels (String imageName, IMetadata meta, ImageProcessorReader reader) throws DependencyException, ServiceException, FormatException, IOException {
        int chs = reader.getSizeC();
        String[] channels = new String[chs];
        String imageExt =  FilenameUtils.getExtension(imageName);
        switch (imageExt) {
            case "nd" :
                for (int n = 0; n < chs; n++) 
                {
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n);
                }
                break;
            case "nd2" :
                for (int n = 0; n < chs; n++) 
                {
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n);
                }
                break;
            case "lif" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null || meta.getChannelName(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n);
                break;
            case "czi" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelFluor(0, n);
                break;
            case "ics" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelExcitationWavelength(0, n).value().toString();
                break;    
            case "ics2" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelExcitationWavelength(0, n).value().toString();
                break;  
            default :
                for (int n = 0; n < chs; n++)
                    channels[n] = Integer.toString(n);
        }
        return(channels);         
    }
    
    /**
     * Dialog ask for channels order and if needed spatial calibration
     * @param channels
     * @param showCal
     * @return ch
     */
    public String[] dialog(String[] channels) {
        String[] models = findStardistModels();
        if (!Arrays.asList(models).contains(stardistModel)) {
            IJ.showMessage("Error", "Missing stardist models");
            return(null);
        }
        String[] channelsName = {"Dapi", "Astro", "Dots"};
        String[] dotsSegMet = {"None","DOG"};
        String[] thMethods = new Thresholder().methods;
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets(0, 120, 0);
        gd.addImage(icon);
        gd.addMessage("Channels", new Font(Font.MONOSPACED , Font.BOLD, 12), Color.blue);
        for (int i =0; i < channelsName.length; i++)
            gd.addChoice(channelsName[i] + " : ", channels, channels[i]);
        gd.addMessage("Nucleus parameters", new Font(Font.MONOSPACED , Font.BOLD, 12), Color.blue);
        gd.addNumericField("Min vol : ", minNucSize, 3);
        gd.addNumericField("Max vol : ", maxNucSize, 3);
        gd.addMessage("Dots parameters", new Font(Font.MONOSPACED , Font.BOLD, 12), Color.blue);
        gd.addChoice("Dots Segmentation Method : ", dotsSegMet, dotsSegMet[0]);
        gd.addChoice("Dots Threshold Method : ", thMethods, thMethods[5]);
        gd.addNumericField("Min vol : ", minDotsSize, 3);
        gd.addNumericField("Max vol : ", maxDotsSize, 3);
        gd.addCheckbox(" pure ARN in astrocyte", pureArn);
        // Calibration
        gd.addMessage("Image calibration", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("XY pixel size : ", cal.pixelWidth, 3);
        gd.addNumericField("Z pixel size : ", cal.pixelDepth, 3);
        
        gd.showDialog();
        String[] ch = new String[channelsName.length];
        for (int i = 0; i < channelsName.length; i++)
            ch[i] = gd.getNextChoice();
        if(gd.wasCanceled())
            ch = null;
        minNucSize = gd.getNextNumber();
        maxNucSize = gd.getNextNumber();
        segMethod = gd.getNextChoice();
        thMethod = gd.getNextChoice();
        minDotsSize = gd.getNextNumber();
        maxDotsSize = gd.getNextNumber();
        pureArn = gd.getNextBoolean();
        cal.pixelWidth = gd.getNextNumber();
        cal.pixelDepth = gd.getNextNumber();
        if(gd.wasCanceled())
            ch = null;
        return(ch);
    }
    
    
    /**
     * Find Z with max intensity in stack
     * @param img
     * @return z
     */
    
    private int find_max(ImagePlus img) {
        double max = 0;
        int zmax = 0;
        for (int z = 1; z <= img.getNSlices(); z++) {
            ImageProcessor ip = img.getStack().getProcessor(z);
            ImageStatistics statistics = new ImageStatistics().getStatistics(ip, ImageStatistics.MEAN, img.getCalibration());
            double meanInt = statistics.mean;
            if (meanInt > max) {
                max = meanInt;
                zmax = z;
            }
        }
        return(zmax);
    }
    
    /**
     * Display detected astrocyte nucleus in green, others in red
     * User can select nucleus with checkboxes
     * @param imgAstro
     * @param imgDots
     * @param nucPop
     * @return selected nucleus
     */
    
    public Object3DInt nucleusSelect(ImagePlus imgAstro, ImagePlus imgDots, Objects3DIntPopulation nucPop) {
        float nucSelected = find_astroNuc(nucPop, imgAstro, imgDots);
        ImagePlus img = imgAstro.duplicate();
        ImageHandler imgObjNuc = ImageHandler.wrap(img).createSameDimensions();
        ImageHandler imgObjNucSelected = ImageHandler.wrap(img).createSameDimensions();
        int nucNb = nucPop.getNbObjects();
        String[] nucIndex = new String[nucNb];
        for (Object3DInt obj : nucPop.getObjects3DInt()) {
            int index = (int)obj.getLabel();
            nucIndex[index - 1] = Integer.toString(index);
            if (index == nucSelected) {
                obj.drawObject(imgObjNucSelected);
                labelObject(obj, imgObjNucSelected.getImagePlus(), 20);
            }
            else {
                obj.drawObject(imgObjNuc);
                labelObject(obj, imgObjNuc.getImagePlus(), 20);
            }
        }
        
        IJ.run(imgObjNuc.getImagePlus(), "Enhance Contrast", "saturated=0.35");
        IJ.run(imgObjNucSelected.getImagePlus(), "Enhance Contrast", "saturated=0.35");
        img.setSlice(img.getNSlices()/2);
        IJ.run(img, "Enhance Contrast", "saturated=0.35");
        img.updateAndDraw();
        ImagePlus[] imgColors = {imgObjNuc.getImagePlus(), imgObjNucSelected.getImagePlus(), null, img, null};
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
        int nucZCenter = new MeasureCentroid(nucPop.getObjectByLabel(nucSelected)).getCentroidAsPoint().getRoundZ();
        imgObjects.setZ(nucZCenter);
        imgObjects.updateAndDraw();
        imgObjects.show("Nucleus");
        
        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Astrocyte nucleus");
        gd.addMessage("Choose astrocyte nucleus");
        gd.addRadioButtonGroup("Nucleus", (String[])nucIndex, 2, 1, nucIndex[(int)nucSelected-1]);
        gd.showDialog();            
        imgObjects.hide();
        flush_close(imgObjects);
        Object3DInt nucSel = nucPop.getObjectByLabelString(gd.getNextRadioButton());
        if(gd.wasCanceled())
            nucSel = null;
        flush_close(img);
        imgObjNuc.closeImagePlus();
        imgObjNucSelected.closeImagePlus();
        return(nucSel);
    }
    
    
    /**
     * Threshold images and fill holes
     * @param img
     * @param thMed
     * @param fill 
     */
    public void threshold(ImagePlus img, String thMed, boolean fill) {
        //  Threshold and binarize
        img.setZ(find_max(img));
        img.updateAndDraw();
        IJ.setAutoThreshold(img, thMed + " dark");
        Prefs.blackBackground = false;
        IJ.run(img, "Convert to Mask", "method="+thMed+" background=Dark");
        if (fill) {
            IJ.run(img,"Fill Holes", "stack");
        }
    }
    
    
    public Objects3DIntPopulation getPopFromImage(ImagePlus img, Calibration cal) {
        // label binary images first
        ImageLabeller labeller = new ImageLabeller();
        ImageInt labels = labeller.getLabels(ImageHandler.wrap(img));
        labels.setCalibration(cal);
        Objects3DIntPopulation pop = new Objects3DIntPopulation(labels);
        return pop;
    }
    
    
     /**
     * Clear out side roi
     * @param img
     * @param roi
     */
    public void clearOutSide(ImagePlus img, Roi roi) {
        PolygonRoi poly = new PolygonRoi(roi.getFloatPolygon(), Roi.FREEROI);
        poly.setLocation(0, 0);
        for (int n = 1; n <= img.getNSlices(); n++) {
            ImageProcessor ip = img.getImageStack().getProcessor(n);
            ip.setRoi(poly);
            ip.setBackgroundValue(0);
            ip.setColor(0);
            ip.fillOutside(poly);
        }
        img.deleteRoi();
        img.updateAndDraw();
    }
    
    
    /**
     * Find mean dot volume
     * return meanVol+SEM
     */
    public double find_mean_dots_volume(Objects3DIntPopulation dotsPop) {
        DescriptiveStatistics dotVolStats = new DescriptiveStatistics();
        for (Object3DInt obj : dotsPop.getObjects3DInt()) {
            dotVolStats.addValue(new MeasureVolume(obj).getVolumeUnit());
        }
        double meanSemVol = dotVolStats.getMean() + dotVolStats.getStandardDeviation() / Math.sqrt(dotVolStats.getN());
        return(meanSemVol);
    }
                            
   /**
     * Remove object with size < min and size > max
     * @param pop
     * @param min
     * @param max
     */
    public void popFilterSize(Objects3DIntPopulation pop, double min, double max) {
        pop.getObjects3DInt().removeIf(p -> (new MeasureVolume(p).getVolumeUnit() < min) || (new MeasureVolume(p).getVolumeUnit() > max));
        pop.resetLabels();
    }
    
    
    
    /**
     * Remove object with only one plan
     * @param pop
     */
    public void popFilterOneZ(Objects3DIntPopulation pop, ImagePlus img) {
        if (img != null) {
            Objects3DIntPopulation excludeBd = new Objects3DIntPopulationComputation(pop).getExcludeBorders(ImageHandler.wrap(img), false);
            popFilterOneZ(excludeBd, null);
        }
        else
            pop.getObjects3DInt().removeIf(p -> (p.getObject3DPlanes().size() == 1));
        pop.resetLabels();
    }
    
    private ImagePlus doZProjection(ImagePlus img, int param) {
        ZProjector zproject = new ZProjector();
        zproject.setMethod(param);
        zproject.setStartSlice(1);
        zproject.setStopSlice(img.getNSlices());
        zproject.setImage(img);
        zproject.doProjection();
       return(zproject.getProjection());
    }
    
    /**
    * Find background image intensity
    * Z project min intensity
    * read mean intensity
    * substract to image
    * @param img 
    */
    public void find_background(ImagePlus img) {
      img.deleteRoi();
      ImagePlus imgProj = doZProjection(img, ZProjector.MIN_METHOD);
      ImageProcessor imp = imgProj.getProcessor();
      bg = imp.getStatistics().mean;
      stdBg = imp.getStatistics().stdDev;
      flush_close(imgProj);
    }
    
    
    /**
     * Find roi volume
     * @param roi
     * @param img
     * @return volume
     */
    private double roi_volume(Roi roi, ImagePlus imgAstro) {
        PolygonRoi poly = new PolygonRoi(roi.getFloatPolygon(), PolygonRoi.FREEROI); 
        poly.setLocation(0, 0);
        ImageProcessor ip = imgAstro.getProcessor();
        ip.setRoi(poly);
        ImageStatistics stats = ImageStatistics.getStatistics(ip, Measurements.AREA, imgAstro.getCalibration());
        double volume = stats.area * imgAstro.getCalibration().pixelDepth * imgAstro.getNSlices();
        return(volume);
    }
    
     /**
     * Apply StarDist 2D slice by slice
     * Label detections in 3D
     * @param img
     * @param factor
     * @param resize
     * @param blockRad
     * @return objects population
     * @throws java.io.IOException
     */
    public Objects3DIntPopulation stardistObjectsPop(ImagePlus img, float factor, boolean resize, int blockRad) throws IOException {
        Object syncObject = new Object();
        double stardistPercentileBottom = 0.2;
        double stardistPercentileTop = 99.8;
        String stardistOutput = "Label Image";
        
        // Resize image to be in a StarDist-friendly scale
        int imgWidth = img.getWidth();
        int imgHeight = img.getHeight();
        String method = (factor > 1) ? "bicubic" : "none";
        ImagePlus imgIn = (resize) ? img.resize((int)(imgWidth*factor), (int)(imgHeight*factor), 1, method) : new Duplicator().run(img);
        
        if (blockRad != 0)
             // Remove outliers
             IJ.run(imgIn, "Remove Outliers...", "radius=20 threshold=1 which=Bright stack");

        // StarDist
        File starDistModelFile = new File(modelsPath+File.separator+stardistModel);
        StarDist2D star = new StarDist2D(syncObject, starDistModelFile);
        star.loadInput(imgIn);
        star.setParams(stardistPercentileBottom, stardistPercentileTop, stardistProbThresh, stardistOverlayThresh, stardistOutput);
        star.run();
        flush_close(imgIn);

        // Label detections in 3D
        ImagePlus imgOut = (resize) ? star.getLabelImagePlus().resize(imgWidth, imgHeight, 1, "none") : star.getLabelImagePlus();       
        ImagePlus imgLabels = star.associateLabels(imgOut);
        imgLabels.setCalibration(cal); 
        flush_close(imgOut);
        Objects3DIntPopulation pop = new Objects3DIntPopulation(ImageHandler.wrap(imgLabels));      
        flush_close(imgLabels);
       return(pop);
    }
    
    
    /**
     * Find astrocyte nucleus keep nucleus 
     * that have the higher/lower ratio between intensity in astrocyte channel
     * and intensity in dots channel
     * @param nucPop
     * @param imgAstro
     * @param imgDots
     * @return nucAstro
     */
    
    public float find_astroNuc( Objects3DIntPopulation nucPop, ImagePlus imgAstro, ImagePlus imgDots) {
        double maxRatioInt = 0;
        float index = 0;      
        if (nucPop.getNbObjects() > 1) {
            for (Object3DInt obj : nucPop.getObjects3DInt()) {
                double astroInt = new MeasureIntensity(obj, ImageHandler.wrap(imgAstro)).getValueMeasurement(MeasureIntensity.INTENSITY_MAX);
                double dotsInt = new MeasureIntensity(obj, ImageHandler.wrap(imgDots)).getValueMeasurement(MeasureIntensity.INTENSITY_AVG);
                double ratio;
                if (ratioInt)
                    ratio = dotsInt/astroInt;
                else
                    ratio = astroInt/dotsInt;
                if (ratio > maxRatioInt) {
                    maxRatioInt = ratio;
                    index = obj.getLabel();
                } 
            }   
        }
        return(index);
    }
    
    
    /**
     * compute local thickness
     * @param imgAstro
     * @return astroMap
    **/
    public ImagePlus localThickness3D (ImagePlus imgAstro) {
        Calibration cal = imgAstro.getCalibration();
        ImagePlus img = imgAstro.duplicate();
        img.setCalibration(cal);
        clij2.median_filter(img, 1, 1);
        clij2.threshold(img, "Li");
        ImageFloat edt = localThickness3D (img, false);
        ImagePlus astroMap = edt.getImagePlus();
        astroMap.setCalibration(img.getCalibration());
        return(astroMap);
    }
    
    
   /**
     * compute local thickness
     * @param img
     * @param inverse
     * @return 
    **/
    public ImageFloat localThickness3D (ImagePlus img, boolean inverse) {
        IJ.showStatus("Computing distance map...");
        img.setCalibration(cal);
        ImageFloat edt = new EDT().run(ImageHandler.wrap(img), 0, inverse, ThreadUtil.getNbCpus());
        return(edt);
    }
    
    /**
     * Find dots population with DOG
     * @param imgDots
     * @param roi
     * @return dotsPop
     */
    
    public Objects3DIntPopulation find_dots_DOG(ImagePlus imgDots, Roi roi) {
        ImagePlus imp = imgDots.duplicate();
        clij2.DOG(imgDots, minDOG, maxDOG);
        threshold(imp, thMethod, false);   
        if (roi != null)
            clearOutSide(imp, roi);
        Objects3DIntPopulation dotsPop = getPopFromImage(imp, cal);
        flush_close(imp);
        return(dotsPop);
    }
   
    
    /**
     * Find dots population
     * @param imgDots
     * @param roi
     * @return dotsPop
     */
    
    public Objects3DIntPopulation find_dots(ImagePlus imgDots, Roi roi) {
        ImagePlus imp = imgDots.duplicate();
        clij2.median_filter(imp, 1, 1);
        threshold(imp, thMethod, false);   
        if (roi != null)
            clearOutSide(imp, roi);
        Objects3DIntPopulation dotsPop = getPopFromImage(imp, cal);
        flush_close(imp);
        return(dotsPop);
    }
    
    /**
     * Label object
     * @param obj
     * @param img 
     * @param fontSize 
     */
    public void labelObject(Object3DInt obj, ImagePlus img, int fontSize) {
        if (IJ.isMacOSX())
            fontSize *= 3;
        
        BoundingBox bb = obj.getBoundingBox();
        int x = bb.xmin;
        int y = bb.ymin;
        img.setSlice(new MeasureCentroid(obj).getCentroidAsPoint().getRoundZ());
        ImageProcessor ip = img.getProcessor();
        ip.setFont(new Font("SansSerif", Font.PLAIN, fontSize));
        ip.setColor(255);
        ip.drawString(Integer.toString((int)obj.getLabel()), x, y);
        img.updateAndDraw();
    }
    
    
   /**
    * draw population
    * nucleus in blue
    * over baseline value = 1 in green
    * under baseline value = 0 in red
    * @param nucObj
    * @param dotsPop
    * @param img
     * @param imgDir
     * @param imgName
     * @param roi
    */
    public void tagsObjects(Object3DInt nucObj, Objects3DIntPopulation dotsPop, ImagePlus img, String imgDir, String imgName, int roi) {  
        ImagePlus imAstro = img.duplicate();
        ImageHandler imgObjNuc = ImageHandler.wrap(imAstro).createSameDimensions();
        nucObj.drawObject(imgObjNuc, 255);                        
        ImageHandler imgObjDotsNotAstro = ImageHandler.wrap(imAstro).createSameDimensions();
        ImageHandler imgObjDotsFineProcess = ImageHandler.wrap(imAstro).createSameDimensions();
        ImageHandler imgObjDotsLargeProcess = ImageHandler.wrap(imAstro).createSameDimensions();
        for (Object3DInt obj : dotsPop.getObjects3DInt()) {
            switch ((int)obj.getIdObject()) {
                case 0:
                    obj.drawObject(imgObjDotsNotAstro, 255);
                    break;
                case 1:
                    obj.drawObject(imgObjDotsFineProcess, 255);
                    break;
                default:
                    obj.drawObject(imgObjDotsLargeProcess, 255);
                    break;
            }
        }
        // save image for objects population
        ImagePlus[] imgColors = {imgObjDotsNotAstro.getImagePlus(), imgObjDotsLargeProcess.getImagePlus(), imgObjNuc.getImagePlus(),
            imAstro, imgObjDotsFineProcess.getImagePlus()};
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
        imgObjects.setCalibration(img.getCalibration());
        for (int i = 1; i <= imgObjects.getNSlices(); i++) {
            imgObjects.setC(i);
            IJ.run(imgObjects, "Enhance Contrast", "saturated=0.35");
            imgObjects.updateAndDraw();
        }
        
        FileSaver ImgObjectsFile = new FileSaver(imgObjects);
        ImgObjectsFile.saveAsTiff(imgDir + imgName + "_Astro-"+(roi+1)+"-Objects.tif"); 
        flush_close(imgObjects);
    }
    
    
   
// Flush and close images
    public void flush_close(ImagePlus img) {
        img.flush();
        img.close();
    } 
         
    /**
    * classify dots
    * @param nucObj nucleus
    * @param dotsPop dots population
    * @param imgAstro read dots intensity
     * @param astroMap
    **/
    public void classify_dots(Object3DInt nucObj, Objects3DIntPopulation dotsPop, ImagePlus imgAstro, ImagePlus astroMap) {
        IJ.showStatus("Classify dots ....");
        find_background(imgAstro);
        // if purearn exclude no dot
        // BgThresholdInt = -1
        double BgThresholdInt = bg + stdBg;
        if (pureArn)
            BgThresholdInt= -1;
        // measure dots mean intensity and volume
        for (Object3DInt obj : dotsPop.getObjects3DInt()) {
            int zmax = obj.getBoundingBox().zmax;
            if (zmax > imgAstro.getNSlices()-1)
                obj.translate(0, 0, -1);
            double meanIntDotAstroimg = new MeasureIntensity(obj,ImageHandler.wrap(imgAstro)).getValueMeasurement(MeasureIntensity.INTENSITY_AVG);
            double astroDiameter = new MeasureIntensity(obj, ImageHandler.wrap(astroMap)).getValueMeasurement(MeasureIntensity.INTENSITY_MAX) * cal.pixelWidth;
            /* classify dots
            * 
            * dots 0  dots <= minThresholdInt and dist > 2 (not in astro (red))
            * dots 1  dots > minThresholdInt  and astroDiameter <= spatial res in Z and dist > 2 (fine processes (magenta))
            * dots 2 dots (larger processes + somatic (green))
            */
            double distNuc = 0;
            if (nucObj.contains(new MeasureCentroid(obj).getCentroidRoundedAsVoxelInt()))
                distNuc = 0;
            else
                distNuc = new Measure2Distance(obj, nucObj).getValue(Measure2Distance.DIST_BB_UNIT);
            // dots 0
            if ((meanIntDotAstroimg <= BgThresholdInt) && (distNuc > 2))
                obj.setIdObject(0);
            // dost 1
            
            else if ((meanIntDotAstroimg > BgThresholdInt) && (astroDiameter <= imgAstro.getCalibration().pixelDepth)  && (distNuc > 2))
                obj.setIdObject(1);
            // dots 2
            else
                obj.setIdObject(2);
        }
    }
    
    /**
     * Compute and write global image parameters
     * @param roi
     * @param imgAstro
     * @param imgAstroMap
     * @param nucObj
     * @param dotsPop
     * @param outDir
     * @param imageName
     * @throws java.io.IOException
     */
    public void compute_Image_parameters(Roi roi, int roiIndex, int totalRoi, ImagePlus imgAstro, ImagePlus imgAstroMap, Object3DInt nucObj, 
            Objects3DIntPopulation dotsPop, BufferedWriter results, String imageName) throws IOException {
        /* write headers
        * largeBrDots (dots in large branches) = dots 2 - dotsinSoma
        * fineBrDots (dots in fine branches) = dots 1
        * dotsNoinAstro (dots outiside astrocyte) = dots 0
        * totalDots (all dots) = dots0 + dots1 + dots2
        * dotsinAstro (dots in astrocyte) = largeBrDots + fineBrDots + dotsinSoma
        * dotsinSoma (dots inside soma) = dots with dist to nucleus < 2
        * dotsinNuc (dots in nucleus) = dots with dist = 0
        * dotsdensinAstro (dots density in astrocyte) = dotsinAstro / astroVolume
        * percDotsNotinAstro (% dots not in astrocyte = dotsNotinAstro / totalDots
        * percDotsinSoma (% dots in soma) = dotsinSoma/totalDotsAstro
        * perDostFineBr (% dots in fine branches) = fineBrDots / totalDotsAstro
        * perDostLargeBr (% dots in large branches) = largeBrDots / totalDotsAstro
        */
        
        // measure roi volume
        double astroVol = roi_volume(roi, imgAstro);
        // Compute nucleus volume
        // add to astro parameters
        double nucVol = new MeasureVolume(nucObj).getVolumeUnit(); 
        float dotsNoinAstro = 0;
        float dotsinFineBr = 0;
        float dotsinLargeBr = 0;
        float dotsinSoma= 0;
        DescriptiveStatistics dotVolStats = new DescriptiveStatistics();
        DescriptiveStatistics dotDiameterStats = new DescriptiveStatistics();
        DescriptiveStatistics dotMeanIntStats = new DescriptiveStatistics();
        String roiName = roi.getName();
        for (Object3DInt dotObj : dotsPop.getObjects3DInt()) {
            double dotVol = new MeasureVolume(dotObj).getVolumeUnit();
            dotVolStats.addValue(dotVol);
            double distNuc = 0;
            if (nucObj.contains(new MeasureCentroid(dotObj).getCentroidRoundedAsVoxelInt()))
                distNuc = 0;
            else {
                distNuc = new Measure2Distance(dotObj, nucObj).getValue(Measure2Distance.DIST_BB_UNIT);
            }
            int Dottype = (int)dotObj.getIdObject();
            double astroDiameter = new MeasureIntensity(dotObj, ImageHandler.wrap(imgAstroMap)).getValueMeasurement(MeasureIntensity.INTENSITY_MAX) * imgAstro.getCalibration().pixelWidth;
            if (astroDiameter != 0)
                dotDiameterStats.addValue(astroDiameter);
            switch (Dottype) {
                case 0 :
                    if (dotVol >= meanSEMDotsSize) {
                        double ratioVol = Math.round(dotVol/meanSEMDotsSize);
                        dotsNoinAstro += ratioVol;
                    }
                    else
                        dotsNoinAstro++;
                    break;
                case 1 :
                    if (dotVol >= meanSEMDotsSize) {
                        double ratioVol = Math.round(dotVol/meanSEMDotsSize);
                        dotsinFineBr += ratioVol;
                    }
                    else
                        dotsinFineBr++;
                    break;
                case 2 : 
                    if (dotVol >= meanSEMDotsSize) {
                        double ratioVol = Math.round(dotVol/meanSEMDotsSize);
                        dotsinLargeBr += ratioVol;
                    }
                    else {
                        dotsinLargeBr++;
                        break;
                    }
            }
            if (distNuc < 2) 
                if (dotVol >= meanSEMDotsSize) {
                    double ratioVol = Math.round(dotVol/meanSEMDotsSize);
                    dotsinSoma+=ratioVol;  
                }   
                else 
                    dotsinSoma++;
            
            double IntdotsinAstro = new MeasureIntensity(dotObj, ImageHandler.wrap(imgAstro)).getValueMeasurement(MeasureIntensity.INTENSITY_AVG);
            dotMeanIntStats.addValue(IntdotsinAstro);
        }
        // exclude dot in soma for large dots
        dotsinLargeBr = dotsinLargeBr - dotsinSoma;
        float totalDots = dotsinLargeBr + dotsinFineBr + dotsNoinAstro + dotsinSoma;
        float dotsinAstro = dotsinLargeBr + dotsinFineBr + dotsinSoma;
        double dotsdensinAstro = dotsinAstro / astroVol;
        double percDotsNotinAstro = 100*(dotsNoinAstro / totalDots);
        double percDotsinSoma = 100*(dotsinSoma / dotsinAstro);
        double perDotsFineBr = 100*(dotsinFineBr / dotsinAstro);
        double perDotsLargeBr = 100*((dotsinLargeBr) / dotsinAstro);
        double meanIntDotsinAstro = dotMeanIntStats.getMean();
        double sdIntDotsinAstro = dotMeanIntStats.getStandardDeviation();
        double meanAstroDiameter = dotDiameterStats.getMean();
        double stdAstroDiameter = dotDiameterStats.getStandardDeviation();
        double medAstroDiameter = dotDiameterStats.getPercentile(50); 
        results.write(imageName + "\t" + roiName + "("+(roiIndex+1)+"/"+totalRoi+")" + "\t" + bg + "\t" + stdBg + "\t" + astroVol + "\t" + dotsdensinAstro +
                "\t" + percDotsNotinAstro + "\t" + percDotsinSoma + "\t" + perDotsFineBr + "\t" + perDotsLargeBr +
                "\t" + meanIntDotsinAstro + "\t" + sdIntDotsinAstro + "\t"+meanAstroDiameter+"\t"+stdAstroDiameter+
                "\t"+medAstroDiameter+"\n");
        results.flush();
    }
    
}
