

setBatchMode(false);
threshold = 12;
dilatation = 0;
min_fiber_area = 10000;
max_fiber_area = 200000;
min_circularity = 0.1;

args = split(getArgument(),";");


//sample = args[0];
//sample = "17-099_HSA-1e12-2.5e12_25-1";
//sample = "17-099_HSA-1e12-2.5e12_15-1";
//sample = "17-099_HSA-1e12-2.5e12_13-1";
sample = "17-099_5e12-1e13-5e13_63-5";
//sample = "17-099_5e12-1e13-5e13_78-2";

input = "C:\\Users\\gcorre\\Desktop\\dev_CSA\\ACD_evelyne";
//input = args[1];
processfile(input,sample);


////// do not change after this line


run("Set Measurements...", "area mean standard modal min centroid center perimeter fit shape feret's median display redirect=None decimal=3");


function processfile(input, sample){
	print("Processing file:"+ sample);

	/// process channels
	for(j=0;j<3;j++){
		for(i=0;i<4;i++){
		run("Bio-Formats Importer", "open="+ input + File.separator +sample+"_z"+i+"_ch0"+j+".tif"+" color_mode=Default rois_import=[ROI manager] view=[Standard ImageJ] stack_order=Default");
		run("Z Project...", "projection=[Max Intensity]");
		close(sample+"_z"+i+"_ch0"+j+".tif");
		rename("z"+i);
		}
		run("Concatenate...", "  title=ch0"+j+" image1=z0 image2=z1 image3=z2 image4=z3");
		run("Z Project...", "projection=[Max Intensity]");
		close("ch0"+j);
	}

	run("Properties...", "channels=1 slices=1 frames=1 pixel_width=1.0000 pixel_height=1.0000 voxel_depth=1.0000 global");
	run("Concatenate...", "  title=RGB image1=MAX_ch00 image2=MAX_ch01 image3=MAX_ch02");

	
	/// identify fibers
	selectWindow("RGB");
	setSlice(2);
	
	run("Duplicate...", "title=membrane");
	run("8-bit");

	selectWindow("membrane");
	run("Enhance Contrast", "saturated=0.35");
	run("Apply LUT");
	run("Gaussian Blur...", "sigma=6");
	run("Find Edges");
	//run("Brightness/Contrast...");
	run("Enhance Contrast", "saturated=0.35");
	run("Apply LUT");
	run("Gaussian Blur...", "sigma=6");

	function segmentFile()
			{

				run("Extended Min & Max", "operation=[Extended Minima] dynamic="+threshold+" connectivity=8");
				run("Impose Min & Max", "original=membrane marker=membrane-emin operation=[Impose Minima] connectivity=8");
				selectWindow("membrane-emin");
				run("Connected Components Labeling", "connectivity=8 type=float");
				run("Marker-controlled Watershed", "input=membrane-imp marker=membrane-emin-lbl mask=None binary calculate use");

			}

	

	
	segmentFile();

	selectWindow("membrane-imp-watershed");
	run("Morphological Filters", "operation=Erosion element=Square radius="+dilatation);
	selectWindow("membrane-imp-watershed-Erosion");
	run("Morphological Filters", "operation=Opening element=Square radius=4");
	selectWindow("membrane-imp-watershed-Erosion-Opening");
	
	setThreshold(-1000000000000000000000000000000, 0);
	run("Convert to Mask");
	run("Invert");
	run("Analyze Particles...", "size="+min_fiber_area+"-"+max_fiber_area+" circularity="+min_circularity+"-1.00 display clear add");
	run("Clear Results");
	close("membrane*");
	
	selectWindow("RGB");
	roiManager("Deselect");
	for(s=1;s<=nSlices();s++){
		setSlice(s);
		roiManager("Measure");
	}	
	roiManager("Deselect");
	roiManager("Save", input+File.separator+sample+"_fibre.roi.zip");
	saveAs("Results", input+File.separator+sample+"_fibre.csv");


	run("Clear Results");
	roiManager("Deselect");
	roiManager("reset");
	run("Select None");
	// identify nuclei
	selectWindow("RGB");
	setSlice(1);
	run("Duplicate...", "title=DAPI");
	run("8-bit");
	run("Enhance Contrast", "saturated=0.35");
	run("Gaussian Blur...", "sigma=2");
	run("Apply LUT");

	setAutoThreshold("Huang dark");
	setOption("BlackBackground", false);
	run("Convert to Mask");
	run("Watershed");
	run("Analyze Particles...", "size=50-10000 circularity=0.10-1.00 display clear add");
	run("Clear Results");
	selectWindow("RGB");
	roiManager("Deselect");
	for(s=1;s<=nSlices();s++){
		setSlice(s);
		roiManager("Measure");
	}
	roiManager("Deselect");
	roiManager("Save", input+File.separator+sample+"_nuclei.roi.zip");
	saveAs("Results", input+File.separator+sample+"_nuclei.csv");

	run("Clear Results");
	roiManager("Deselect");
	roiManager("reset");
	run("Select None");

	// identify ACD spots

	selectWindow("RGB");
	setSlice(3);
	run("Duplicate...", "title=ACD");
	
	run("Variance...", "radius=1");
	run("Gaussian Blur...", "sigma=2");
	run("Find Maxima...", "prominence=50 output=[Single Points]");
	selectWindow("ACD Maxima");
	run("Dilate");
	run("Analyze Particles...", "size=2-50 circularity=0.10-1.00 display clear add");
	run("Clear Results");
	selectWindow("RGB");
	roiManager("Deselect");
	for(s=1;s<=nSlices();s++){
		setSlice(s);
		roiManager("Measure");
	}
	roiManager("Deselect");
	roiManager("Save", input+File.separator+sample+"_spots.roi.zip");
	saveAs("Results", input+File.separator+sample+"_spots.csv");

	selectWindow("RGB");
	saveAs("Tiff", input+File.separator+sample+"_Zmax.tif");
	
	close("*");

}

//eval("script","System.exit(0);");
