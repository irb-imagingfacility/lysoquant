// LysoQuant from an ImageJ macro: segment and count the lysosomes of the open image,
// without a dialog. Works in the Fiji GUI and in headless Fiji.
//
//   ImageJ-macosx --headless -macro macro_example.ijm
//
// The segmentation runs on the LysoQuant server given in Edit>Options>LysoQuant Settings...
// A macro can set it too (the keys are the ones of the settings; add savePreferences to keep them):

call("ij.Prefs.set", "lysoquant.remotehost", "true");              // false: the server runs on this computer
call("ij.Prefs.set", "lysoquant.server", "gpu-server.example.org");
call("ij.Prefs.set", "lysoquant.port", "22");
call("ij.Prefs.set", "lysoquant.username", "lysoquant");
call("ij.Prefs.set", "lysoquant.rsaKeyFilename", "/home/me/.ssh/lysoquant.key");   // required without a display
call("ij.Prefs.set", "lysoquant.processfolder", "/jobs");
call("ij.Prefs.set", "lysoquant.hide_segmentation", "true");       // do not open the LQ_ image
//call("ij.Prefs.savePreferences");

setBatchMode(true);
open("/path/to/image.tif");        // or use the image that is already open

// Segment and count. The keys are the first words of the fields of the dialog:
//   lysosome, protein, slices, frames, show (single values), measurement (channels),
//   measure (cargo area), threshold
run("LysoQuant", "lysosome=2 protein=3 show measurement=1-4");

// Windows: the Results table has one row per lysosome and measurement channel (with "show"),
// the "LysoQuant" table one row per image, slice and frame.
// Headless: no windows, the Results table (read it with getResult) holds the lysosomes with "show",
// else the totals.
print("lysosomes measured: " + nResults);
print("first one: " + getResultString("Lysosome Type", 0) + ", area " + getResult("Area", 0));

// Only the segmentation: the label image (0 background, 1 empty, 2 loaded lysosomes) is the new image
run("LysoQuant Segment...", "lysosome=2 protein=3");
print("labels: " + getTitle() + ", " + getWidth() + "x" + getHeight());
