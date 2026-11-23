// =====================================================================================
//  C2C12_MyHC2_Quant.ijm      (Fiji / ImageJ 1.53k or newer)
//
//  Semi-automated version of your manual MyHC2 protocol (DAPI + phalloidin + MyHC2 images):
//    1. Calibrate (Image > Properties), split channels, display-only brightness/contrast
//    2. DAPI nuclei: duplicate -> Gaussian -> threshold (YOU adjust) -> Watershed -> Analyze Particles
//    3. MyHC2 mask: duplicate -> (optional Subtract Background) -> Gaussian -> threshold, starts at Default (YOU adjust)
//         -> Convert to Mask -> Open / Close / Fill Holes -> OUTLINE -> Analyze Particles (include holes)
//    4. MyHC2-positive area fraction, number of MyHC2 objects
//    5. Differentiation index and fusion index (nuclei inside MyHC2 objects / total DAPI nuclei)
//    6. Myotube morphology (Feret length, MinFeret width, aspect ratio, area, circularity, angle)
//    7. Alignment: (A) object-based (angles of the fitted ellipses), (B) structure tensor on the MyHC2
//       channel = what OrientationJ computes (coherency + energy thresholds), optional OrientationJ plugin outputs
//
//  The Outline step before Analyze Particles is kept as a fixed part of the workflow (you validated it manually).
//  Analyze Particles runs on the outlined image with "include holes", so each object is measured as the full
//  filled myotube. The macro then CHECKS that the objects really sit on MyHC2-positive pixels (polarity guard).
//
//  YOU stay in control of: brightness/contrast, DAPI threshold and MyHC2 threshold.
//  Results, per image, in  <folder>/C2C12_MyHC2_results/<image>/ :
//     <image>_summary.csv                one row of results/quantifications
//     <image>_MyHC2_objects.csv          one row per MyHC2 object (size, shape, angle, nuclei count, ...)
//     <image>_nuclei_assignment.csv      one row per nucleus (which MyHC2 object it lies in)
//     <image>_DAPI_nuclei_ROIs.zip, <image>_DAPI_mask.tif
//     <image>_MyHC2_mask.tif             filled binary mask (after threshold + Open/Close/Fill Holes)
//     <image>_MyHC2_outline.tif          the outlined image that Analyze Particles runs on
//     <image>_MyHC2_objects_ROIs.zip     one ROI per MyHC2 object
//     <image>_overlay_composite.png      MyHC2 (red) + phalloidin (green) + DAPI (blue)
//     <image>_overlay_MyHC2_nuclei.png   MyHC2 (red) + DAPI (blue) with outlines: yellow = myotube, orange = other MyHC2
//                                         object, green = nucleus inside a MyHC2 object, cyan = other nucleus
//     <image>_settings.txt
//  and one master table  C2C12_MyHC2_results/ALL_images_summary.csv  (one row per image).
//  Brightness/contrast is never applied to pixel data; the macro checks this and stops if "Apply" was clicked.
// =====================================================================================
requires("1.53k");

// ---------- settings (filled by the dialog) ----------
var myhcCh, phalCh, dapiCh, pxSize, pxMode, fovWidth;
var sigmaDAPI, thrMethodDAPI, sizeStrDAPI, sizeInPxDAPI, exclEdgesDAPI;
var sigmaMyhc, thrMethodMyhc, doBGsub, rollRadius, doOpen, doClose, doFill, useOutline, sizeStrMyhc, sizeInPxMyhc, exclEdgesMyhc;
var minNuc, refAngle, minAR, tensorSigmaUm, minCoh, minEnergyPct, runOJ, adjustMode;
var outRoot;
// ---------- settings remembered from the first image (for "reuse" mode) ----------
var haveStored = false;
var stDapiLo = -1, stDapiHi = -1, stMyhcLo = -1, stMyhcHi = -1;
var stMyhcMin = 0, stMyhcMax = 0, stDapiMin = 0, stDapiMax = 0;
// ---------- summary row builder / structure-tensor results ----------
var sumH = "", sumV = "";
var stUsedPct = 0, stR = 0, stDomAng = 0, stSref = 0, stScore10 = 0, stScore15 = 0, stOK = false;

// ============================== MAIN ==============================
setBatchMode(false);
setOption("BlackBackground", true);
run("Colors...", "foreground=white background=black selection=yellow");
run("Conversions...", "scale weighted");
run("Options...", "iterations=1 count=1 black do=Nothing");

// ---- Step 1 of 3: images, calibration, DAPI nuclei ----
Dialog.create("MyHC2 quantification - step 1 of 3: images and DAPI");
Dialog.addChoice("Images to process", newArray("Active image only", "All images in a folder (TIF/JPG/PNG)"), "Active image only");
Dialog.addMessage("Channels:  C1 = MyHC2 (red), C2 = phalloidin (green), C3 = DAPI (blue).\nJPG/PNG images are converted with Image > Color > Make Composite.");
Dialog.addNumber("MyHC2 channel", 1, 0, 3, "");
Dialog.addToSameRow();
Dialog.addNumber("Phalloidin", 2, 0, 3, "");
Dialog.addToSameRow();
Dialog.addNumber("DAPI", 3, 0, 3, "");
Dialog.addNumber("Pixel size (um per pixel)", 1.298, 4, 7, "um");
Dialog.addChoice("Pixel size source", newArray("Fixed pixel size (above)", "From field-of-view width (below)"), "Fixed pixel size (above)");
Dialog.addNumber("Field-of-view width (um) - for JPG/PNG", 1329.2, 1, 7, "um");
Dialog.addNumber("DAPI Gaussian blur sigma (px)", 1);
Dialog.addString("Nucleus size filter", "0-Infinity");
Dialog.addCheckbox("Size filter in PIXEL units", true);
Dialog.addToSameRow();
Dialog.addCheckbox("Exclude nuclei on edges", true);
Dialog.addChoice("Manual adjustments (B&C, thresholds)", newArray("On every image", "On the first image only, reuse for the rest"), "On every image");
Dialog.show();
mode          = Dialog.getChoice();
pxMode        = Dialog.getChoice();
adjustMode    = Dialog.getChoice();
myhcCh        = Dialog.getNumber();
phalCh        = Dialog.getNumber();
dapiCh        = Dialog.getNumber();
pxSize        = Dialog.getNumber();
fovWidth      = Dialog.getNumber();
sigmaDAPI     = Dialog.getNumber();
sizeStrDAPI   = Dialog.getString();
sizeInPxDAPI  = Dialog.getCheckbox();
exclEdgesDAPI = Dialog.getCheckbox();

// ---- Step 2 of 3: MyHC2 mask ----
Dialog.create("MyHC2 quantification - step 2 of 3: MyHC2 mask");
Dialog.addNumber("MyHC2 Gaussian blur sigma (px, 0 = none)", 1);
Dialog.addCheckbox("Subtract Background first", false);
Dialog.addToSameRow();
Dialog.addNumber("rolling ball radius (px)", 50);
Dialog.addCheckbox("Binary Open", true);
Dialog.addToSameRow();
Dialog.addCheckbox("Binary Close", true);
Dialog.addToSameRow();
Dialog.addCheckbox("Fill Holes", true);
Dialog.addCheckbox("OUTLINE before Analyze Particles (keep ticked)", true);
Dialog.addString("MyHC2 object size filter", "150-Infinity");
Dialog.addCheckbox("Size filter in PIXEL units", true);
Dialog.addToSameRow();
Dialog.addCheckbox("Exclude objects on edges", true);
Dialog.addMessage("The threshold always starts at Default; you adjust it by hand.");
Dialog.show();
sigmaMyhc     = Dialog.getNumber();
rollRadius    = Dialog.getNumber();
doBGsub       = Dialog.getCheckbox();
doOpen        = Dialog.getCheckbox();
doClose       = Dialog.getCheckbox();
doFill        = Dialog.getCheckbox();
useOutline    = Dialog.getCheckbox();
sizeStrMyhc   = Dialog.getString();
sizeInPxMyhc  = Dialog.getCheckbox();
exclEdgesMyhc = Dialog.getCheckbox();

// ---- Step 3 of 3: fusion, morphology, alignment ----
Dialog.create("MyHC2 quantification - step 3 of 3: fusion and alignment");
Dialog.addNumber("Myotube = object with at least this many nuclei", 2);
Dialog.addNumber("Reference direction (degrees, 0 = horizontal)", 0);
Dialog.addNumber("Alignment: minimum aspect ratio", 2);
Dialog.addNumber("Structure tensor sigma (um)", 13, 1, 6, "um");
Dialog.addNumber("Structure tensor: minimum coherency (0-1)", 0.2);
Dialog.addNumber("Structure tensor: minimum energy (%)", 0.1, 2, 6, "%");
Dialog.addCheckbox("Also run the OrientationJ plugin (colour map + histogram)", true);
Dialog.show();
minNuc        = Dialog.getNumber();
refAngle      = Dialog.getNumber();
minAR         = Dialog.getNumber();
tensorSigmaUm = Dialog.getNumber();
minCoh        = Dialog.getNumber();
minEnergyPct  = Dialog.getNumber();
runOJ         = Dialog.getCheckbox();

thrMethodDAPI = "Default";      // threshold always starts at Default (then adjusted by hand)
thrMethodMyhc = "Default";

if (myhcCh == dapiCh || myhcCh == phalCh || phalCh == dapiCh) exit("Channel numbers must all be different.");
if (minNuc < 1) exit("A myotube needs at least 1 nucleus (use 2 or 3).");

nDone = 0;
if (mode == "Active image only") {
    if (nImages == 0) exit("Open an image first, then run the macro again.");
    id0 = getImageID();
    title0 = getTitle();
    dir0 = getDirectory("image");
    if (dir0 == "") dir0 = getDirectory("Choose the folder where results should be saved");
    outRoot = dir0 + "C2C12_MyHC2_results" + File.separator;
    File.makeDirectory(outRoot);
    ok = processImage(id0, stripExt(title0));
    if (ok) nDone++;
} else {
    inDir = getDirectory("Choose the folder with your images (TIF / JPG / PNG)");
    outRoot = inDir + "C2C12_MyHC2_results" + File.separator;
    File.makeDirectory(outRoot);
    list = getFileList(inDir);
    for (i = 0; i < list.length; i++) {
        f = list[i];
        lf = toLowerCase(f);
        if (endsWith(f, "/")) continue;
        if (!(endsWith(lf, ".tif") || endsWith(lf, ".tiff") || endsWith(lf, ".jpg") || endsWith(lf, ".jpeg") || endsWith(lf, ".png"))) continue;
        open(inDir + f);
        id = getImageID();
        ok = processImage(id, stripExt(f));
        closeID(id);
        if (ok) nDone++;
    }
}
compileMaster();
roiManager("Reset");
run("Clear Results");
showMessage("C2C12 MyHC2 quantification finished", "Images processed: " + nDone + "\n\nResults folder:\n" + outRoot + "\n\nMaster table: ALL_images_summary.csv");


// ============================== FUNCTIONS ==============================
function stripExt(s) {
    i = lastIndexOf(s, ".");
    if (i > 0) return substring(s, 0, i);
    return s;
}

function closeID(id) {
    if (isOpen(id)) {
        selectImage(id);
        close();
    }
}

// save and close every image / table window the OrientationJ plugin has created
function saveAndCloseOJ(bp, tag) {
    ti = getList("image.titles");
    for (i = 0; i < ti.length; i++) {
        if (startsWith(ti[i], "OJ-")) {
            selectImage(ti[i]);
            saveAs("PNG", bp + "_OrientationJ_" + tag + "_" + i + ".png");
            close();
        }
    }
    tw = getList("window.titles");
    for (i = 0; i < tw.length; i++) {
        if (startsWith(tw[i], "OJ-") || startsWith(tw[i], "OrientationJ")) {
            selectWindow(tw[i]);
            saveAs("Text", bp + "_OrientationJ_" + tag + "_table" + i + ".csv");
            run("Close");
        }
    }
}

// close every image whose title starts with the given prefix
function closeByPrefix(prefix) {
    t = getList("image.titles");
    for (i = 0; i < t.length; i++) {
        if (startsWith(t[i], prefix)) {
            selectImage(t[i]);
            close();
        }
    }
}

function median(arr) {
    a = Array.copy(arr);
    Array.sort(a);
    m = a.length;
    if (m % 2 == 1) return a[(m - 1) / 2];
    return (a[m / 2 - 1] + a[m / 2]) / 2;
}

// angular deviation (0-90 deg) between an axis angle and the reference direction (axes are 180-degree periodic)
function absDev(a, ref) {
    d = (a - ref) % 180;
    if (d < 0) d = d + 180;
    if (d > 90) d = 180 - d;
    return d;
}

// append one column (name, value) to the summary row being built
function addCol(name, v) {
    sumH = sumH + name + ",";
    sumV = sumV + v + ",";
}

// ---------- structure tensor alignment (= what OrientationJ computes), done with ImageJ image maths ----------
// Sets globals: stOK, stUsedPct, stR, stDomAng, stSref, stScore10, stScore15.
// Orientation convention: counter-clockwise from the image x-axis, like the "Angle" of Fit Ellipse.
function structureTensor(srcID, sigPx, minC, minEpct, refDeg) {
    stOK = false;
    selectImage(srcID);
    run("Select None");
    run("Duplicate...", "title=ST_I");
    run("32-bit");
    run("Gaussian Blur...", "sigma=1");
    run("Duplicate...", "title=ST_Ix");
    run("Convolve...", "text1=[-1 0 1\n-2 0 2\n-1 0 1\n] normalize");
    selectImage("ST_I");
    run("Duplicate...", "title=ST_Iy");
    run("Convolve...", "text1=[-1 -2 -1\n0 0 0\n1 2 1\n] normalize");
    imageCalculator("Multiply create 32-bit", "ST_Ix", "ST_Ix"); rename("ST_Jxx");
    imageCalculator("Multiply create 32-bit", "ST_Iy", "ST_Iy"); rename("ST_Jyy");
    imageCalculator("Multiply create 32-bit", "ST_Ix", "ST_Iy"); rename("ST_Jxy");
    nm = newArray("ST_Jxx", "ST_Jyy", "ST_Jxy");
    for (k = 0; k < 3; k++) {
        selectImage(nm[k]);
        run("Gaussian Blur...", "sigma=" + sigPx);
    }
    selectImage("ST_Jxy");
    run("Multiply...", "value=2");                                   // now 2*Jxy
    imageCalculator("Subtract create 32-bit", "ST_Jxx", "ST_Jyy"); rename("ST_D");
    imageCalculator("Add create 32-bit", "ST_Jxx", "ST_Jyy"); rename("ST_E");
    imageCalculator("Multiply create 32-bit", "ST_D", "ST_D"); rename("ST_D2");
    imageCalculator("Multiply create 32-bit", "ST_Jxy", "ST_Jxy"); rename("ST_X2");
    imageCalculator("Add create 32-bit", "ST_D2", "ST_X2"); rename("ST_N");
    run("Square Root");
    run("Add...", "value=0.000000001");
    selectImage("ST_E");
    run("Add...", "value=0.000000001");
    imageCalculator("Divide create 32-bit", "ST_N", "ST_E"); rename("ST_C");        // coherency
    imageCalculator("Divide create 32-bit", "ST_D", "ST_N"); rename("ST_cos");      // cos(2*gradient angle)
    imageCalculator("Divide create 32-bit", "ST_Jxy", "ST_N"); rename("ST_sin");    // sin(2*gradient angle)
    selectImage("ST_cos");
    run("Multiply...", "value=-1");                                  // fibre direction is perpendicular to the gradient:  cos2t = -cos2g
    // masks: coherency >= minC and energy >= minEpct % of max energy
    selectImage("ST_C");
    run("Duplicate...", "title=ST_mC");
    setThreshold(minC, 1000000000);
    run("Convert to Mask");
    run("Divide...", "value=255");
    run("32-bit");
    selectImage("ST_E");
    getStatistics(aa, mE, mnE, maxE);
    run("Duplicate...", "title=ST_mE");
    setThreshold(minEpct / 100 * maxE, 1000000000);
    run("Convert to Mask");
    run("Divide...", "value=255");
    run("32-bit");
    imageCalculator("Multiply create 32-bit", "ST_mC", "ST_mE"); rename("ST_m");
    getStatistics(aa, mUsed);
    stUsedPct = 100 * mUsed;
    imageCalculator("Multiply create 32-bit", "ST_C", "ST_m"); rename("ST_w");     // weight = coherency inside the mask
    getStatistics(aa, wMean);
    if (wMean <= 0) {
        closeByPrefix("ST_");
        return;
    }
    imageCalculator("Multiply create 32-bit", "ST_w", "ST_cos"); rename("ST_wc");
    getStatistics(aa, wcMean);
    imageCalculator("Multiply create 32-bit", "ST_w", "ST_sin"); rename("ST_ws");
    getStatistics(aa, wsMean);
    cm = wcMean / wMean;
    sm = wsMean / wMean;
    stR = sqrt(cm * cm + sm * sm);                                   // nematic order magnitude (0 random ... 1 perfectly aligned)
    a2 = atan2(sm, cm);                                              // = 2 * dominant angle
    stDomAng = 0.5 * a2 * 180 / PI;
    if (stDomAng < 0) stDomAng = stDomAng + 180;
    r2 = 2 * refDeg * PI / 180;
    stSref = cos(r2) * cm + sin(r2) * sm;                            // S = <cos 2(theta - reference)>
    // alignment score: share of weighted signal within +-w degrees of the reference direction
    selectImage("ST_cos");
    run("Duplicate...", "title=ST_Q1");
    run("Multiply...", "value=" + d2s(cos(r2), 8));
    selectImage("ST_sin");
    run("Duplicate...", "title=ST_Q2");
    run("Multiply...", "value=" + d2s(sin(r2), 8));
    imageCalculator("Add create 32-bit", "ST_Q1", "ST_Q2"); rename("ST_Q");          // cos 2(theta - reference) per pixel
    wins = newArray(10, 15);
    scs = newArray(2);
    for (k = 0; k < 2; k++) {
        selectImage("ST_Q");
        run("Duplicate...", "title=ST_win");
        setThreshold(cos(2 * wins[k] * PI / 180), 2);
        run("Convert to Mask");
        run("Divide...", "value=255");
        run("32-bit");
        imageCalculator("Multiply create 32-bit", "ST_w", "ST_win"); rename("ST_ww");
        getStatistics(aa, wwMean);
        scs[k] = 100 * wwMean / wMean;
        selectImage("ST_win"); close();
        selectImage("ST_ww"); close();
    }
    stScore10 = scs[0];
    stScore15 = scs[1];
    stOK = true;
    closeByPrefix("ST_");
}

function processImage(origID, base) {
    outDir = outRoot + base + File.separator;
    File.makeDirectory(outDir);
    bp = outDir + base;
    print("------ " + base);
    sumH = ""; sumV = "";

    // ---------- 1. PREPARE THE IMAGE (work on a copy, the original is never touched) ----------
    selectImage(origID);
    run("Duplicate...", "title=work duplicate");
    workID = getImageID();
    srcType = "multichannel image";
    if (bitDepth() == 24) {     // RGB image (JPG/PNG): Image > Color > Make Composite -> C1 = red, C2 = green, C3 = blue
        run("Make Composite", "display=Composite");
        workID = getImageID();
        srcType = "RGB image converted to a 3-channel composite (8-bit, display-scaled)";
    }
    getDimensions(w, h, c, z, t);
    if (c == 1 && z == 3) {     // 3-slice stack instead of 3 channels
        run("Stack to Hyperstack...", "order=xyczt(default) channels=3 slices=1 frames=1 display=Composite");
        workID = getImageID();
        getDimensions(w, h, c, z, t);
    }
    if (c != 3 || z != 1) {
        print("   SKIPPED: expected a 3-channel single-plane image, found c=" + c + " z=" + z + " t=" + t);
        closeID(workID);
        return false;
    }
    px = pxSize;
    if (pxMode == "From field-of-view width (below)") px = fovWidth / w;
    print("   " + srcType + ", " + w + " x " + h + " px, pixel size " + d2s(px, 4) + " um");
    setVoxelSize(px, px, 1, "micron");                      // = Image > Properties (unit micron, pixel width = height)
    imgArea = w * h * px * px;                              // total image area in um^2
    run("Split Channels");
    selectWindow("C" + myhcCh + "-work"); myhcID = getImageID();
    selectWindow("C" + phalCh + "-work"); phalID = getImageID();
    selectWindow("C" + dapiCh + "-work"); dapiID = getImageID();

    interactive = (adjustMode == "On every image") || !haveStored;

    selectImage(myhcID); getStatistics(aa, myhcMean0);      // raw values BEFORE touching the display
    selectImage(dapiID); getStatistics(aa, dapiMean0);

    allIDs = newArray(myhcID, dapiID, phalID);
    for (k = 0; k < 3; k++) {
        selectImage(allIDs[k]);
        run("Enhance Contrast", "saturated=0.35");          // display only
    }
    if (interactive) {
        selectImage(myhcID);
        run("Brightness/Contrast...");
        waitForUser("MyHC2 display", "Click on the MyHC2 window, then adjust Brightness/Contrast (Auto or sliders).\nDo NOT click Apply.\n\nClick OK here when you are happy.");
        selectImage(dapiID);
        run("Brightness/Contrast...");
        waitForUser("DAPI display", "Click on the DAPI window, then adjust Brightness/Contrast (Auto or sliders).\nDo NOT click Apply.\n\nClick OK here when you are happy.");
    } else {
        selectImage(myhcID); setMinAndMax(stMyhcMin, stMyhcMax);
        selectImage(dapiID); setMinAndMax(stDapiMin, stDapiMax);
    }
    close("B&C");
    selectImage(myhcID); getMinAndMax(myhcMin, myhcMax);
    selectImage(dapiID); getMinAndMax(dapiMin, dapiMax);
    selectImage(phalID); getMinAndMax(phalMin, phalMax);

    selectImage(myhcID); getStatistics(aa, myhcMean1);
    selectImage(dapiID); getStatistics(aa, dapiMean1);
    if (abs(myhcMean1 - myhcMean0) > 0.001 || abs(dapiMean1 - dapiMean0) > 0.001)
        exit("Brightness/Contrast 'Apply' was clicked - the pixel data were changed.\nYour original file is untouched. Run the macro again and do NOT click Apply.");

    // ---------- 2. COUNT DAPI NUCLEI (same as the myogenin macro) ----------
    selectImage(dapiID);
    run("Select None");
    run("Duplicate...", "title=DAPI_analysis");
    anID = getImageID();
    run("Gaussian Blur...", "sigma=" + sigmaDAPI);
    setAutoThreshold(thrMethodDAPI + " dark");
    if (interactive) {
        run("Threshold...");
        waitForUser("DAPI threshold", "Adjust the DAPI threshold in the Threshold window so the red area covers the nuclei\n(start: " + thrMethodDAPI + ", dark background). Do NOT click Apply.\n\nClick OK here when you are happy.");
    } else {
        setThreshold(stDapiLo, stDapiHi);
    }
    selectImage(anID);
    getThreshold(dLo, dHi);
    dLoTxt = "NA";
    if (dLo != -1) {
        dLoTxt = d2s(dLo, 3);
        if (!is("binary")) run("Convert to Mask");
    } else if (!is("binary")) {
        exit("No DAPI threshold was set. Run the macro again and adjust the threshold in the Threshold window.");
    }
    close("Threshold");
    run("Watershed");
    saveAs("Tiff", bp + "_DAPI_mask.tif");

    run("Set Measurements...", "area mean shape feret's display redirect=None decimal=4");
    run("Clear Results");
    roiManager("Reset");
    opt = "size=" + sizeStrDAPI;
    if (sizeInPxDAPI) opt = opt + " pixel";
    opt = opt + " circularity=0.00-1.00 show=Nothing display clear";
    if (exclEdgesDAPI) opt = opt + " exclude";
    opt = opt + " add";
    run("Analyze Particles...", opt);
    n = nResults;
    if (n == 0) {
        print("   No nuclei detected - check the DAPI threshold.");
        closeID(myhcID); closeID(phalID); closeID(dapiID); closeID(anID);
        return false;
    }
    nucArea = newArray(n);
    for (i = 0; i < n; i++) nucArea[i] = getResult("Area", i);
    run("Clear Results");
    if (roiManager("count") != n) exit("ROI Manager count (" + roiManager("count") + ") differs from the nucleus count (" + n + ").");
    roiManager("Deselect");
    nucZip = bp + "_DAPI_nuclei_ROIs.zip";
    roiManager("Save", nucZip);
    print("   DAPI nuclei: " + n);

    // ---------- 3. MyHC2 MASK ----------
    selectImage(myhcID);
    run("Select None");
    run("Duplicate...", "title=MyHC2_analysis");
    mAnID = getImageID();
    if (doBGsub) run("Subtract Background...", "rolling=" + rollRadius);
    if (sigmaMyhc > 0) run("Gaussian Blur...", "sigma=" + sigmaMyhc);
    setAutoThreshold(thrMethodMyhc + " dark");
    if (interactive) {
        run("Threshold...");
        waitForUser("MyHC2 threshold", "Adjust the MyHC2 threshold in the Threshold window.\nMyotubes must be RED and the background black (tick 'Dark background' / click Invert if not).\nDo NOT click Apply.\n\nClick OK here when you are happy.");
    } else {
        setThreshold(stMyhcLo, stMyhcHi);
    }
    selectImage(mAnID);
    getThreshold(mLo, mHi);
    mLoTxt = "NA";
    if (mLo != -1) {
        mLoTxt = d2s(mLo, 3);
        if (!is("binary")) run("Convert to Mask");
    } else if (!is("binary")) {
        exit("No MyHC2 threshold was set. Run the macro again and adjust the threshold in the Threshold window.");
    }
    close("Threshold");
    getStatistics(aa, mThrMean);
    fracThr = 100 * mThrMean / 255;                          // % of image above threshold, before any cleaning

    // polarity guard: if more than half of the image is selected, the background is probably selected
    polQC = "ok";
    if (fracThr > 50) {
        if (interactive) {
            if (getBoolean("More than 50% of the image is selected (" + d2s(fracThr, 1) + "%).\nThat usually means the BACKGROUND was selected instead of the myotubes.\n\nInvert the mask?")) {
                run("Invert");
                getStatistics(aa, mThrMean);
                fracThr = 100 * mThrMean / 255;
                polQC = "inverted by user";
            } else polQC = "CHECK: >50% selected";
        } else {
            polQC = "CHECK: >50% selected";
        }
    }

    // clean the mask: Open (specks), Close (small gaps), Fill Holes
    if (doOpen) run("Open");
    if (doClose) run("Close-");
    if (doFill) run("Fill Holes");
    getStatistics(aa, mCleanMean);
    fracClean = 100 * mCleanMean / 255;
    saveAs("Tiff", bp + "_MyHC2_mask.tif");                 // filled binary mask (image stays open)
    mMaskID = getImageID();

    // OUTLINE, then Analyze Particles on the outlined image (include holes -> full filled objects)
    selectImage(mMaskID);
    run("Select None");
    run("Duplicate...", "title=MyHC2_outline");
    olID = getImageID();
    if (useOutline) run("Outline");
    saveAs("Tiff", bp + "_MyHC2_outline.tif");

    run("Set Measurements...", "area mean perimeter shape feret's display redirect=None decimal=4");
    run("Clear Results");
    roiManager("Reset");
    opt = "size=" + sizeStrMyhc;
    if (sizeInPxMyhc) opt = opt + " pixel";
    opt = opt + " circularity=0.00-1.00 show=[Count Masks] display clear include";
    if (exclEdgesMyhc) opt = opt + " exclude";
    opt = opt + " add";
    run("Analyze Particles...", opt);
    nObj = nResults;
    if (nObj == 0) print("   !! No MyHC2 objects found. After cleaning, the mask covers " + d2s(fracClean, 1) + " % of the image. Check the MyHC2 threshold (myotubes red on a black background) and the object size filter (" + sizeStrMyhc + ").");
    cmID = 0;
    if (nObj > 0) cmID = getImageID();                      // "Count Masks of ..." (object k = grey value k)
    run("Clear Results");
    print("   MyHC2 objects: " + nObj);

    // per-object arrays (measured on the RAW MyHC2 channel, ROI geometry from the outlined image)
    nPer = newArray(1);
    nucObj = newArray(n); nucFrac = newArray(n);
    for (i = 0; i < n; i++) { nucObj[i] = 0; nucFrac[i] = 0; }
    objZip = bp + "_MyHC2_objects_ROIs.zip";
    myoArea = 0;
    inMaskQC = 1;
    nIn = 0; nInMyo = 0; nMyo = 0;
    if (nObj > 0) {
        if (roiManager("count") != nObj) exit("ROI Manager count (" + roiManager("count") + ") differs from the MyHC2 object count (" + nObj + ").");
        roiManager("Deselect");
        roiManager("Save", objZip);
        run("Set Measurements...", "area mean integrated perimeter shape feret's fit display redirect=None decimal=4");
        selectImage(myhcID);
        roiManager("Deselect");
        roiManager("Measure");
        oArea = newArray(nObj); oMean = newArray(nObj); oInt = newArray(nObj); oPer = newArray(nObj); oCirc = newArray(nObj);
        oSol = newArray(nObj); oFeret = newArray(nObj); oMinF = newArray(nObj); oMajor = newArray(nObj); oMinor = newArray(nObj); oAng = newArray(nObj);
        for (i = 0; i < nObj; i++) {
            oArea[i]  = getResult("Area", i);
            oMean[i]  = getResult("Mean", i);
            oInt[i]   = getResult("IntDen", i);
            oPer[i]   = getResult("Perim.", i);
            oCirc[i]  = getResult("Circ.", i);
            oSol[i]   = getResult("Solidity", i);
            oFeret[i] = getResult("Feret", i);
            oMinF[i]  = getResult("MinFeret", i);
            oMajor[i] = getResult("Major", i);
            oMinor[i] = getResult("Minor", i);
            oAng[i]   = getResult("Angle", i);
            myoArea = myoArea + oArea[i];
        }
        run("Clear Results");
        // polarity / outline guard: the objects must sit on MyHC2-positive pixels of the filled mask
        run("Set Measurements...", "mean display redirect=None decimal=4");
        selectImage(mMaskID);
        roiManager("Deselect");
        roiManager("Measure");
        sIn = 0;
        oIn = newArray(nObj);
        for (i = 0; i < nObj; i++) {
            oIn[i] = getResult("Mean", i) / 255;
            sIn = sIn + oIn[i];
        }
        inMaskQC = sIn / nObj;
        run("Clear Results");
        if (inMaskQC < 0.8) {
            polQC = "CHECK: objects not on MyHC2 mask (" + d2s(100 * inMaskQC, 0) + "% inside)";
            print("   !! " + polQC + " - the background may have been analysed. Check the threshold polarity.");
        }

        // ---------- 4. WHICH NUCLEI LIE IN WHICH MyHC2 OBJECT ----------
        roiManager("Reset");
        roiManager("Open", nucZip);
        run("Set Measurements...", "min display redirect=None decimal=4");
        selectImage(cmID);
        roiManager("Deselect");
        roiManager("Measure");
        nucLab = newArray(n);
        for (i = 0; i < n; i++) nucLab[i] = getResult("Max", i);
        run("Clear Results");
        run("Set Measurements...", "mean display redirect=None decimal=4");
        selectImage(mMaskID);
        roiManager("Deselect");
        roiManager("Measure");
        for (i = 0; i < n; i++) nucFrac[i] = getResult("Mean", i) / 255;
        run("Clear Results");
        nPer = newArray(nObj);
        for (k = 0; k < nObj; k++) nPer[k] = 0;
        for (i = 0; i < n; i++) {
            if (nucFrac[i] >= 0.5 && nucLab[i] > 0) {       // nucleus counted if at least half of it lies inside the object
                nucObj[i] = nucLab[i];
                nPer[nucLab[i] - 1] = nPer[nucLab[i] - 1] + 1;
                nIn++;
            }
        }
        for (k = 0; k < nObj; k++) {
            if (nPer[k] >= minNuc) { nInMyo = nInMyo + nPer[k]; nMyo++; }
        }
    }
    diffIdx = 100 * nIn / n;                                // differentiation index (nuclei in MyHC2 objects / all nuclei)
    fusionIdx = 100 * nInMyo / n;                           // fusion index (nuclei in myotubes / all nuclei)
    fusionAmong = "NA";
    if (nIn > 0) fusionAmong = d2s(100 * nInMyo / nIn, 3);  // alternative denominator: nuclei in MyHC2 objects
    areaFrac = 100 * myoArea / imgArea;

    // ---------- 5. MORPHOLOGY OF MYOTUBES + OBJECT-BASED ALIGNMENT ----------
    lenA = newArray(0); widA = newArray(0); arA = newArray(0); areaA = newArray(0); circA = newArray(0); nucA = newArray(0);
    nAl = 0; sumDev = 0; c10 = 0; c15 = 0; cs = 0; sn = 0;
    oAR = newArray(1);
    isMyo = newArray(1);
    if (nObj > 0) {
        oAR = newArray(nObj); isMyo = newArray(nObj);
        for (k = 0; k < nObj; k++) {
            oAR[k] = 0;
            if (oMinF[k] > 0) oAR[k] = oFeret[k] / oMinF[k];
            isMyo[k] = 0;
            if (nPer[k] >= minNuc) {
                isMyo[k] = 1;
                lenA = Array.concat(lenA, oFeret[k]);
                widA = Array.concat(widA, oMinF[k]);
                arA = Array.concat(arA, oAR[k]);
                areaA = Array.concat(areaA, oArea[k]);
                circA = Array.concat(circA, oCirc[k]);
                nucA = Array.concat(nucA, nPer[k]);
                if (oAR[k] >= minAR) {
                    dev = absDev(oAng[k], refAngle);
                    sumDev = sumDev + dev;
                    if (dev <= 10) c10++;
                    if (dev <= 15) c15++;
                    th = 2 * (oAng[k] - refAngle) * PI / 180;
                    cs = cs + cos(th);
                    sn = sn + sin(th);
                    nAl++;
                }
            }
        }
    }

    // ---------- 6. SAVE: per-object CSV and per-nucleus CSV ----------
    qBase = "\"" + base + "\"";
    s = "Image,Object,Area_um2,Perimeter_um,Feret_length_um,MinFeret_width_um,Major_um,Minor_um,AspectRatio_Feret,Circ,Solidity,Angle_deg,AngleDev_from_reference_deg,MyHC2_Mean,MyHC2_IntDen,Fraction_inside_mask,N_nuclei,Is_myotube,Used_for_alignment\n";
    for (k = 0; k < nObj; k++) {
        used = 0;
        if (isMyo[k] == 1 && oAR[k] >= minAR) used = 1;
        s = s + qBase + "," + (k + 1) + "," + d2s(oArea[k], 4) + "," + d2s(oPer[k], 4) + "," + d2s(oFeret[k], 4) + "," + d2s(oMinF[k], 4) + ","
              + d2s(oMajor[k], 4) + "," + d2s(oMinor[k], 4) + "," + d2s(oAR[k], 4) + "," + d2s(oCirc[k], 4) + "," + d2s(oSol[k], 4) + ","
              + d2s(oAng[k], 3) + "," + d2s(absDev(oAng[k], refAngle), 3) + "," + d2s(oMean[k], 4) + "," + d2s(oInt[k], 4) + "," + d2s(oIn[k], 3) + ","
              + nPer[k] + "," + isMyo[k] + "," + used + "\n";
    }
    File.saveString(s, bp + "_MyHC2_objects.csv");
    s = "Image,Nucleus,Area_um2,Fraction_inside_MyHC2_mask,MyHC2_object,Object_is_myotube\n";
    for (i = 0; i < n; i++) {
        inMy = 0;
        if (nucObj[i] > 0) inMy = isMyo[nucObj[i] - 1];
        s = s + qBase + "," + (i + 1) + "," + d2s(nucArea[i], 4) + "," + d2s(nucFrac[i], 4) + "," + nucObj[i] + "," + inMy + "\n";
    }
    File.saveString(s, bp + "_nuclei_assignment.csv");

    // ---------- 7. PIXEL-LEVEL ALIGNMENT (structure tensor on the MyHC2 channel) ----------
    sigPx = tensorSigmaUm / px;
    structureTensor(myhcID, sigPx, minCoh, minEnergyPct, refAngle);

    // ---------- 8. SUMMARY ROW ----------
    getDateAndTime(yr, mo, dw, dm, hr, mi, sc, ms);
    stamp = "" + yr + "-" + IJ.pad(mo + 1, 2) + "-" + IJ.pad(dm, 2) + " " + IJ.pad(hr, 2) + ":" + IJ.pad(mi, 2);
    exclN = 0; if (exclEdgesDAPI) exclN = 1;
    bgs = 0; if (doBGsub) bgs = 1;
    op = 0; if (doOpen) op = 1;
    cl = 0; if (doClose) cl = 1;
    fh = 0; if (doFill) fh = 1;
    ou = 0; if (useOutline) ou = 1;
    addCol("Image", qBase);
    addCol("Method", "Fiji macro");
    addCol("Date", stamp);
    addCol("Pixel_um", d2s(px, 4));
    addCol("Image_area_um2", d2s(imgArea, 1));
    addCol("DAPI_threshold_lower", dLoTxt);
    addCol("Gaussian_sigma_DAPI", sigmaDAPI);
    addCol("Nuclei_size_filter", sizeStrDAPI);
    addCol("Exclude_edges_nuclei", exclN);
    addCol("Total_DAPI_nuclei", n);
    addCol("Median_nuclear_area_um2", d2s(median(nucArea), 3));
    addCol("MyHC2_threshold_method", thrMethodMyhc);
    addCol("MyHC2_threshold_lower", mLoTxt);
    addCol("MyHC2_BG_subtract", bgs);
    addCol("MyHC2_gaussian_sigma", sigmaMyhc);
    addCol("Binary_open", op);
    addCol("Binary_close", cl);
    addCol("Fill_holes", fh);
    addCol("Outline_before_particles", ou);
    addCol("MyHC2_min_object_size", sizeStrMyhc);
    addCol("MyHC2_area_fraction_threshold_only_pct", d2s(fracThr, 3));
    addCol("MyHC2_area_fraction_after_cleaning_pct", d2s(fracClean, 3));
    addCol("N_MyHC2_objects", nObj);
    addCol("MyHC2_area_um2", d2s(myoArea, 1));
    addCol("MyHC2_area_fraction_pct", d2s(areaFrac, 3));
    addCol("Nuclei_in_MyHC2_objects", nIn);
    addCol("Differentiation_index_pct", d2s(diffIdx, 3));
    addCol("Myotube_min_nuclei", minNuc);
    addCol("N_myotubes", nMyo);
    addCol("Myotube_density_per_mm2", d2s(nMyo / (imgArea / 1000000), 3));
    addCol("Nuclei_in_myotubes", nInMyo);
    addCol("Fusion_index_pct", d2s(fusionIdx, 3));
    addCol("Fusion_index_among_MyHC2_nuclei_pct", fusionAmong);
    if (nMyo > 0) {
        Array.getStatistics(nucA, mnA, mxA, avgNuc, sdA);
        Array.getStatistics(lenA, mnA, mxA, avgLen, sdA);
        Array.getStatistics(widA, mnA, mxA, avgWid, sdA);
        Array.getStatistics(arA, mnA, mxA, avgAR, sdA);
        Array.getStatistics(areaA, mnA, mxA, avgAreaM, sdA);
        Array.getStatistics(circA, mnA, mxA, avgCirc, sdA);
        addCol("Mean_nuclei_per_myotube", d2s(avgNuc, 3));
        addCol("Mean_myotube_length_um", d2s(avgLen, 3));
        addCol("Median_myotube_length_um", d2s(median(lenA), 3));
        addCol("Mean_myotube_width_um", d2s(avgWid, 3));
        addCol("Mean_myotube_aspect_ratio", d2s(avgAR, 3));
        addCol("Mean_myotube_area_um2", d2s(avgAreaM, 3));
        addCol("Mean_myotube_circularity", d2s(avgCirc, 4));
    } else {
        addCol("Mean_nuclei_per_myotube", "NA");
        addCol("Mean_myotube_length_um", "NA");
        addCol("Median_myotube_length_um", "NA");
        addCol("Mean_myotube_width_um", "NA");
        addCol("Mean_myotube_aspect_ratio", "NA");
        addCol("Mean_myotube_area_um2", "NA");
        addCol("Mean_myotube_circularity", "NA");
    }
    addCol("Reference_angle_deg", refAngle);
    addCol("Alignment_min_aspect_ratio", minAR);
    addCol("N_myotubes_used_for_alignment", nAl);
    if (nAl > 0) {
        rObj = sqrt(cs * cs + sn * sn) / nAl;
        addCol("Mean_abs_angle_dev_deg", d2s(sumDev / nAl, 3));
        addCol("Pct_within_10deg", d2s(100 * c10 / nAl, 3));
        addCol("Pct_within_15deg", d2s(100 * c15 / nAl, 3));
        addCol("Circular_variance", d2s(1 - rObj, 4));
        addCol("Order_parameter_S_objects", d2s(cs / nAl, 4));
    } else {
        addCol("Mean_abs_angle_dev_deg", "NA");
        addCol("Pct_within_10deg", "NA");
        addCol("Pct_within_15deg", "NA");
        addCol("Circular_variance", "NA");
        addCol("Order_parameter_S_objects", "NA");
    }
    if (stOK && stUsedPct < 5) print("   !! structure tensor used only " + d2s(stUsedPct, 1) + " % of the pixels - lower the minimum energy / coherency, otherwise the alignment values are unreliable.");
    addCol("ST_sigma_um", tensorSigmaUm);
    addCol("ST_min_coherency", minCoh);
    addCol("ST_min_energy_pct", minEnergyPct);
    if (stOK) {
        addCol("ST_pixels_used_pct", d2s(stUsedPct, 2));
        addCol("ST_dominant_angle_deg", d2s(stDomAng, 2));
        addCol("ST_order_parameter_S", d2s(stSref, 4));
        addCol("ST_nematic_magnitude_R", d2s(stR, 4));
        addCol("ST_alignment_score_10deg_pct", d2s(stScore10, 2));
        addCol("ST_alignment_score_15deg_pct", d2s(stScore15, 2));
    } else {
        addCol("ST_pixels_used_pct", "NA");
        addCol("ST_dominant_angle_deg", "NA");
        addCol("ST_order_parameter_S", "NA");
        addCol("ST_nematic_magnitude_R", "NA");
        addCol("ST_alignment_score_10deg_pct", "NA");
        addCol("ST_alignment_score_15deg_pct", "NA");
    }
    addCol("Mask_QC", "\"" + polQC + "\"");
    hdr = substring(sumH, 0, lengthOf(sumH) - 1);
    row = substring(sumV, 0, lengthOf(sumV) - 1);
    File.saveString(hdr + "\n" + row + "\n", bp + "_summary.csv");

    st = "DAPI threshold method (start): " + thrMethodDAPI + "\nDAPI threshold lower (on Gaussian-blurred DAPI): " + dLoTxt
       + "\nMyHC2 threshold method (start): " + thrMethodMyhc + "\nMyHC2 threshold lower/upper: " + mLoTxt + " / " + d2s(mHi, 3)
       + "\nMyHC2 Subtract Background: " + doBGsub + " (radius " + rollRadius + ")\nMyHC2 Gaussian sigma (px): " + sigmaMyhc
       + "\nBinary Open / Close / Fill Holes: " + doOpen + " / " + doClose + " / " + doFill + "\nOutline before Analyze Particles: " + useOutline
       + "\nMyHC2 object size filter: " + sizeStrMyhc + " (pixel units: " + sizeInPxMyhc + ")\nMyotube = object with at least " + minNuc + " nuclei"
       + "\nA nucleus belongs to an object if at least 50% of its area lies inside it"
       + "\nReference direction (deg, counter-clockwise from horizontal): " + refAngle
       + "\nStructure tensor sigma (um / px): " + tensorSigmaUm + " / " + d2s(sigPx, 2) + "   min coherency " + minCoh + "   min energy % " + minEnergyPct
       + "\nMyHC2 display range: " + d2s(myhcMin, 2) + " - " + d2s(myhcMax, 2) + "\nDAPI display range: " + d2s(dapiMin, 2) + " - " + d2s(dapiMax, 2)
       + "\nPixel size (um): " + px + "\nSource: " + srcType + " (" + w + " x " + h + " px)\n";
    File.saveString(st, bp + "_settings.txt");

    // ---------- 9. OVERLAYS ----------
    selectImage(myhcID); run("Select None");
    selectImage(phalID); run("Select None");
    selectImage(dapiID); run("Select None");
    selectImage(myhcID); run("Duplicate...", "title=m8"); setMinAndMax(myhcMin, myhcMax); run("8-bit");
    selectImage(phalID); run("Duplicate...", "title=p8"); setMinAndMax(phalMin, phalMax); run("8-bit");
    selectImage(dapiID); run("Duplicate...", "title=d8"); setMinAndMax(dapiMin, dapiMax); run("8-bit");
    run("Merge Channels...", "c1=m8 c2=p8 c3=d8 keep");     // MyHC2 red, phalloidin green, DAPI blue -> RGB (sources kept)
    compID = getImageID();
    saveAs("PNG", bp + "_overlay_composite.png");
    closeID(compID);
    run("Merge Channels...", "c1=m8 c3=d8");                // MyHC2 red + DAPI blue (sources consumed)
    rgbID = getImageID();
    roiManager("Reset");
    roiManager("Open", nucZip);
    if (nObj > 0) roiManager("Open", objZip);
    grpA = newArray(0); grpB = newArray(0); grpC = newArray(0); grpD = newArray(0);
    for (i = 0; i < n; i++) {
        if (nucObj[i] > 0) grpA = Array.concat(grpA, i);
        else grpB = Array.concat(grpB, i);
    }
    for (k = 0; k < nObj; k++) {
        if (isMyo[k] == 1) grpC = Array.concat(grpC, n + k);
        else grpD = Array.concat(grpD, n + k);
    }
    if (grpA.length > 0) { roiManager("Select", grpA); roiManager("Set Color", "green"); }
    if (grpB.length > 0) { roiManager("Select", grpB); roiManager("Set Color", "cyan"); }
    if (grpC.length > 0) { roiManager("Select", grpC); roiManager("Set Color", "yellow"); }
    if (grpD.length > 0) { roiManager("Select", grpD); roiManager("Set Color", "orange"); }
    roiManager("Deselect");
    selectImage(rgbID);
    run("From ROI Manager");
    run("Flatten");
    saveAs("PNG", bp + "_overlay_MyHC2_nuclei.png");
    flatID = getImageID();

    // ---------- optional: OrientationJ plugin outputs (colour map + histogram), saved as files only ----------
    // Runs LAST, so a plugin problem can never lose the results saved above.
    if (runOJ) {
        List.setCommands;
        if (List.get("OrientationJ Distribution") == "") {
            print("   OrientationJ is not installed (Help > Update... > Manage update sites > OrientationJ) - skipped.");
        } else {
            selectImage(myhcID);
            run("Select None");
            run("Duplicate...", "title=OJ_input");
            run("8-bit");
            ojID = getImageID();
            selectImage(ojID);
            run("OrientationJ Analysis", "tensor=" + d2s(sigPx, 1) + " gradient=0 color-survey=on hue=Orientation sat=Coherency bri=Original-Image ");
            saveAndCloseOJ(bp, "colormap");
            selectImage(ojID);                  // Distribution needs the 8-bit image, not the RGB colour survey left in front by Analysis
            run("OrientationJ Distribution", "tensor=" + d2s(sigPx, 1) + " gradient=0 histogram=on table=on min-coherency=" + d2s(minCoh * 100, 1) + " min-energy=" + d2s(minEnergyPct, 1) + " ");
            saveAndCloseOJ(bp, "distribution");
            closeID(ojID);
        }
    }

    // ---------- remember settings for "reuse" mode ----------
    if (interactive && adjustMode != "On every image") {
        if (dLo == -1 || mLo == -1) {
            print("   (A threshold was applied inside the Threshold window, so it cannot be reused. Every image will ask again.)");
        } else {
            stDapiLo = dLo; stDapiHi = dHi; stMyhcLo = mLo; stMyhcHi = mHi;
            stMyhcMin = myhcMin; stMyhcMax = myhcMax; stDapiMin = dapiMin; stDapiMax = dapiMax;
            haveStored = true;
        }
    }

    print("   MyHC2 area fraction " + d2s(areaFrac, 2) + " %   objects " + nObj + "   myotubes (>= " + minNuc + " nuclei) " + nMyo);
    print("   differentiation index " + d2s(diffIdx, 2) + " %   fusion index " + d2s(fusionIdx, 2) + " %   mask QC: " + polQC);

    // ---------- clean up ----------
    closeID(flatID); closeID(rgbID);
    closeID(myhcID); closeID(phalID); closeID(dapiID); closeID(anID); closeID(mAnID); closeID(mMaskID); closeID(olID); closeID(cmID);
    closeByPrefix("m8"); closeByPrefix("p8"); closeByPrefix("d8");
    closeByPrefix("Count Masks");
    roiManager("Reset");
    run("Clear Results");
    return true;
}

// Combine every <image>_summary.csv into one master table
function compileMaster() {
    list = getFileList(outRoot);
    master = "";
    header = "";
    for (i = 0; i < list.length; i++) {
        if (!endsWith(list[i], "/")) continue;
        sub = outRoot + list[i];
        files = getFileList(sub);
        for (j = 0; j < files.length; j++) {
            if (endsWith(files[j], "_summary.csv")) {
                txt = File.openAsString(sub + files[j]);
                lines = split(txt, "\n\r");
                if (lines.length < 2) continue;
                if (header == "") { header = lines[0]; master = header + "\n"; }
                master = master + lines[1] + "\n";
            }
        }
    }
    if (master != "") File.saveString(master, outRoot + "ALL_images_summary.csv");
}
