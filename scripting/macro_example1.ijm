// Cycle through all open images
for (i = 1; i <= nImages; i++) {
	selectImage(i);
	// Run Lysoquant and generate only summary table 
	run("LysoQuant", "lysosome=2 protein=3");
}