// Start-up macro for ImageJ
//
// With this macro you can press a button to open the LysoQuant dialog.
// Put this code in the Plugins > Macro > Startup Macros... file and restart
// ImageJ. Then open an image and press "l" to start LysoQuant

macro "Run LysoQuant [l]" {
      run("LysoQuant");
}