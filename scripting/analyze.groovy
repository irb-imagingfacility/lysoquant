// LysoQuant from a Fiji script (Groovy; the same works in Jython): segment and count the lysosomes
// of an image, without a dialog. Works headless too:
//   ImageJ-macosx --headless --run analyze.groovy 'image="/path/to/image.tif"'
//
// The server is the one of Edit>Options>LysoQuant Settings...; name only what differs.

#@ ImagePlus image
#@ CommandService command
#@output ImagePlus labels
#@output ResultsTable totals
#@output ResultsTable lysosomes

inputs = [
    "image", image,
    "lysosomeChannel", 2,
    "proteinChannel", 3,
    "singleValues", true,           // one row per lysosome and measurement channel in "lysosomes"
    // "remote", false,             // the server runs on this computer
    // "server", "gpu-server.example.org", "username", "lysoquant",
    // "keyFile", "/home/me/.ssh/lysoquant.key",     // required without a display
    // "processFolder", "/jobs",
]
module = command.run("ch.irb.lysoquant.AnalyzeCommand", true, inputs as Object[]).get()

labels = module.getOutput("labels")         // 16-bit: 0 background, 1 empty, 2 loaded lysosomes
totals = module.getOutput("totals")
lysosomes = module.getOutput("singles")

// Segment and count separately (the segmentation of a Z/T range is counted with the same range):
//   labels = command.run("ch.irb.lysoquant.SegmentCommand", true, "image", image).get().getOutput("labels")
//   tables = command.run("ch.irb.lysoquant.CountCommand", true, "image", image, "labels", labels).get()
