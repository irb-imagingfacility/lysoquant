inputDir=getDirectory("Choose folder with images");
images=getFileList(inputDir);
for (i = 0; i < lengthOf(images); i++) {
	run("Bio-Formats Importer", "open=[" + inputDir + File.separator() + images[i] + "] autoscale color_mode=Colorized view=Hyperstack stack_order=XYCZT");
	// Run Lysoquant and generate only summary table 
	run("LysoQuant", "lysosome=2 protein=3");
	close("*");
}