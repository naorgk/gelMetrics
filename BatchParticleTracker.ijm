// ===== Batch Particle Tracker Macro =====

inputFolder  = "";    // Folder containing nd2 images

// Particle Tracker parameters
// These parameters worked for our data - user may need to tweak them

radius       = 5;
cutoff       = 0.001;
percentile   = 0.1;
linkrange    = 2;
displacement = 10;

list = getFileList(inputFolder);
for (i = 0; i < list.length; i++) {
    if (endsWith(list[i], ".nd2")) {
        print("Processing: " + list[i]);
	run("Bio-Formats Importer", "open=" + inputFolder + list[i] + " autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	rename(list[i]);

        // Run Mosaic Particle Tracker
	print("Running tracker on " + list[i]);
	run("Particle Tracker 2D/3D", "radius=" + radius + " cutoff=" + cutoff + " per/abs=" + percentile + " link=" + 	linkrange + " displacement=" + displacement + " dynamics=Brownian saveMss");

	selectImage(list[i]);
	close;
    }
}
print("Batch processing complete!");