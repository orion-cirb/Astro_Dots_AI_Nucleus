package Astrocytes;



import Astrocytes_Tools.AstroDotsAINucleus_Tools;
import ij.IJ;
import static ij.IJ.setMinAndMax;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.plugin.Duplicator;
import ij.plugin.PlugIn;
import ij.plugin.frame.RoiManager;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.Region;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.ArrayUtils;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author phm
 */



public class AstroDots_AI_Nucleus implements PlugIn {
    
    private String imageDir;
    private final boolean canceled = false;
    public  String outDirResults = "";
    private BufferedWriter outPutResults;
    
    
   
            
    Astrocytes_Tools.AstroDotsAINucleus_Tools tools = new AstroDotsAINucleus_Tools();
            
    /**
     * 
     * @param arg
     */
    @Override
    public void run(String arg) {
        
       
        try {
            if (canceled) {
                IJ.showMessage(" Pluging canceled");
                return;
            }
            imageDir = IJ.getDirectory("Choose Directory Containing image Files...");
            if (imageDir == null) {
                return;
            }
            // Find images with file_ext extension ics or ics2
            String file_ext = "ics2";
            ArrayList<String> imageFiles = tools.findImages(imageDir, file_ext);
            if (imageFiles.size() == 0) {
                IJ.showStatus("Error", "No images found with " + file_ext + " extension");
                file_ext = "ics";
                IJ.showStatus("Trying to find images with " + file_ext + " extension");
                imageFiles = tools.findImages(imageDir, file_ext);
                if (imageFiles.size() == 0) {
                    IJ.showMessage("Error", "No images found with " + file_ext + " extension");
                    return;
                }
            }

            // create output folder
            outDirResults = imageDir + File.separator+ "Results"+ File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }
            // Write headers for results file
            FileWriter fileResults = new FileWriter(outDirResults + "Astro_results.xls", false);
            outPutResults = new BufferedWriter(fileResults);
            outPutResults.write("ImageName\tRoi Name \tMean background\tStd background\tAstrocyte Volume\tDensity dots in Astro"
                    + "\tPercentage of dots not in astro\tPercentage of dots in soma\tPercentage of dots in fine processes\tPercentage of dots in large processes"
                    + "\tDots mean intensity in Astro\tSD intensity in astro\tMean astro diameter(0 exluded)\tStd astro diameter(0 excluded)"
                    + "\tMed astro diameter(0 excluded)\n");
            outPutResults.flush();
            
            
            // Create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            reader.setId(imageFiles.get(0));
            
            // Find image calibration
            tools.cal = tools.findImageCalib(meta);
            
             // Find channel names
            String[] chsName = tools.findChannels(imageFiles.get(0), meta, reader);
            
            // Channels dialog
            String[] channels = tools.dialog(chsName);
            if (channels == null) {
                IJ.showStatus("Plugin cancelled");
                return;
            } 
            for (String f : imageFiles) {
                String rootName = FilenameUtils.getBaseName(f);
                // Find ROI file
                String roi_file = imageDir+rootName+".zip";
                if (!new File(roi_file).exists()) {
                    IJ.showStatus("No ROI file found !");
                    return;
                }

                // find rois
                RoiManager rm = new RoiManager(false);
                rm.runCommand("Open", roi_file);
                // for each roi open image and crop
                for (Roi roi : rm.getRoisAsArray()) {
                    int r = rm.getRoiIndex(roi);
                    // Find in roi name the desired top and bottom stack 
                    // roi name should be roi_number-ztop-zbottom
                    String[] regExp = roi.getName().split("-");
                    int zStart = (Integer.parseInt(regExp[1]) < 1) ? 1 : Integer.parseInt(regExp[1]);
                    int zStop = (Integer.parseInt(regExp[2]) > reader.getSizeZ()) ? reader.getSizeZ() : Integer.parseInt(regExp[2]);
                    ImporterOptions options = new ImporterOptions();
                    options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                    options.setId(f);
                    options.setSplitChannels(true);
                    options.setCrop(true);
                    Region reg = new Region(roi.getBounds().x, roi.getBounds().y, roi.getBounds().width, roi.getBounds().height);
                    options.setCropRegion(0, reg);
                    options.doCrop();
                    options.setZBegin(0, zStart);
                    options.setZEnd(0, zStop);
                    
                
                    /**
                    * Open channels
                    */
                    
                    // Nucleus channel
                    int indexCh = ArrayUtils.indexOf(chsName, channels[0]);
                    System.out.println("Opening Nucleus channel");
                    ImagePlus imgNuc = BF.openImagePlus(options)[indexCh];
                    if (imgNuc.getBitDepth() == 32) {
                        setMinAndMax(imgNuc,0, 65535);
                        IJ.run(imgNuc, "16-bit", "");
                    }
                        
                    // Find nucleus
                    Objects3DIntPopulation nucPop = tools.stardistObjectsPop(imgNuc, 0, false, 20);
                    tools.popFilterOneZ(nucPop, imgNuc);
                    System.out.println("Roi = "+r+" of " + rm.getCount());
                    System.out.println("Nucleus number = "+nucPop.getNbObjects());
                    tools.popFilterSize(nucPop, tools.minNucSize, tools.maxNucSize);
                    System.out.println("After size filter Nucleus number = "+nucPop.getNbObjects());
                    if (nucPop.getNbObjects() != 0) { 
                        // Astrocyte channel
                        IJ.showStatus("Opening Astrocyte channel");
                        indexCh = ArrayUtils.indexOf(chsName, channels[1]);
                        ImagePlus imgAstro = BF.openImagePlus(options)[indexCh];
                        imgAstro.setTitle(rootName+"_Astro");
                        if (imgAstro.getBitDepth() == 32) {
                            setMinAndMax(imgAstro,0, 65535);
                            IJ.run(imgAstro, "16-bit", "");
                        }
                        
                        // Dots channel
                        IJ.showStatus("Opening Dots channel");
                        indexCh = ArrayUtils.indexOf(chsName, channels[2]);
                        ImagePlus imgDots = BF.openImagePlus(options)[indexCh];
                        imgDots.setTitle(rootName+"_Dots");
                        if (imgDots.getBitDepth() == 32) {
                            setMinAndMax(imgDots,0, 65535);
                            IJ.run(imgDots, "16-bit", "");
                        }
                        
                        // if more than one nucleus find nucleus with intensity in astrocyte channel
                        // and ask to choose
                        Object3DInt nucAstro = null;
                        if (nucPop.getNbObjects() > 1) 
                            nucAstro = tools.nucleusSelect(imgAstro, imgDots, nucPop);                                  
                        else
                            nucAstro = nucPop.getFirstObject();
                                
                        // compute distance map image
                        ImagePlus imgAstroMap = tools.localThickness3D(imgAstro);
                        imgAstroMap.setCalibration(tools.cal);        
                        
                        
                        
                        // Find dots population
                        Objects3DIntPopulation dotsPop = (tools.segMethod.equals("None")) ? tools.find_dots(imgDots, roi) : tools.find_dots_DOG(imgDots, roi);
                        System.out.println("Dots number = "+dotsPop.getNbObjects());
                        // Duplicate dots Population to filter dots without "macro-dots" to calculate the mean dot size
                        Objects3DIntPopulation dotsPopNotMacro = dotsPop;
                        tools.popFilterSize(dotsPopNotMacro, tools.minDotsSize, tools.maxDotsSize);
                        tools.meanSEMDotsSize = tools.find_mean_dots_volume(dotsPopNotMacro);
                        System.out.println("Dots mean size volume + sem = " + tools.meanSEMDotsSize);
                        // min size filter including macro-dots
                        tools.popFilterSize(dotsPop, tools.minDotsSize, Double.MAX_VALUE);
                        System.out.println("After min size filter dots number = "+dotsPop.getNbObjects());

                        // calculate parameters
                        tools.classify_dots(nucAstro, dotsPop, imgAstro, imgAstroMap);                              

                        // draw objects
                        tools.tagsObjects(nucAstro, dotsPop, imgAstro, outDirResults, rootName, r);                               

                        // write a global parameters table by image
                        tools.compute_Image_parameters(roi, r, rm.getCount(), imgAstro, imgAstroMap, nucAstro, dotsPop, outPutResults, rootName);
                        
                        //close images
                        tools.flush_close(imgNuc);
                        tools.flush_close(imgAstro);
                        tools.flush_close(imgAstroMap);
                    }
                    else 
                        System.out.println("No nucleus found !");
                }
            }
            outPutResults.close();
        }
        catch (IOException | FormatException | DependencyException | ServiceException ex) {
            Logger.getLogger(AstroDots_AI_Nucleus.class.getName()).log(Level.SEVERE, null, ex);
        }
        IJ.showStatus("Process done ...");
    }
}
