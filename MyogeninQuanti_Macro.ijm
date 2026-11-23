// =====================================================================================
//  C2C12_Myogenin_Quant.ijm      (Fiji / ImageJ 1.53k or newer)
//
//  Semi-automated version of your manual protocol:
//    1. Calibrate (Image > Properties), split channels, display-only brightness/contrast
//    2. DAPI nuclei: duplicate -> Gaussian -> threshold (YOU adjust) -> Watershed -> Analyze Particles
//    3. Myogenin per nucleus (Mean, Median, IntDen, RawIntDen) measured on the RAW channel
//    4. Cell-free background ROIs (YOU draw them) -> threshold = mean_bg + k * SD_bg
//    5. Myogenin-positive fraction = 100 * positive nuclei / total DAPI nuclei
//
//  YOU stay in control of: brightness/contrast, DAPI threshold, background ROIs.
//  The macro does everything else and saves, for every image, in  <folder>/C2C12_results/<image>/ :
//     <image>_summary.csv                    one row of results/quantifications
//     <image>_myogenin_per_nucleus.csv       one row per nucleus
//     <image>_background_ROIs.csv / .zip     the background ROIs you drew
//     <image>_DAPI_nuclei_ROIs.zip           nuclei ROIs (open with ROI Manager)
//     <image>_DAPI_mask.tif                  binary mask after threshold + watershed
//     <image>_DAPI_outlines.png              Analyze Particles outline drawing
//     <image>_overlay_DAPI_myogenin.png      DAPI (blue) + myogenin (green) overlay
//     <image>_overlay_nuclei_calls.png       same + nuclei outlines (yellow = myogenin+, cyan = negative, red = bg ROIs)
//     <image>_settings.txt                   thresholds / display ranges used
//  and one master table  C2C12_results/ALL_images_summary.csv  (all images, one row each).
//
//  Brightness/contrast is never applied to pixel data. The macro checks this and stops if "Apply" was clicked.
// =====================================================================================
requires("1.53k");

// ---------- settings (filled by the dialog) ----------
var myogCh, phalCh, dapiCh, pxSize, pxMode, fovWidth, overlayRule, sigmaDAPI, thrMethod, sizeStr, sizeInPx, exclEdges, nBG, sdMult, adjustMode;
var outRoot;
// ---------- settings remembered from the first image (for "reuse" mode) ----------
var haveStored = false;
var stDapiLo = -1, stDapiHi = -1, stMyogMin = 0, stMyogMax = 0, stDapiMin = 0, stDapiMax = 0;
var stBgMean = 0, stBgSD = 0, stBgN = 0;

// ============================== MAIN ==============================
setBatchMode(false);
setOption("BlackBackground", true);
run("Colors...", "foreground=white background=black selection=yellow");
run("Conversions...", "scale weighted");

Dialog.create("C2C12 DAPI + myogenin quantification");
Dialog.addChoice("Images to process", newArray("Active image only", "All images in a folder (TIF/JPG/PNG)"), "Active image only");
Dialog.addMessage("Channel order:  C1 = myogenin (red), C2 = phalloidin (green), C3 = DAPI (blue)\nRGB images (JPG/PNG) are converted with Image > Color > Make Composite automatically.");
Dialog.addNumber("Myogenin channel (C#)", 1);
Dialog.addNumber("Phalloidin channel (C#)", 2);
Dialog.addNumber("DAPI channel (C#)", 3);
Dialog.addNumber("Pixel size (um per pixel)", 1.298, 4, 7, "um");
Dialog.addChoice("Pixel size source", newArray("Fixed pixel size (above)", "From field-of-view width (below)"), "Fixed pixel size (above)");
Dialog.addNumber("Field-of-view width (um) - use for JPG/PNG exports", 1329.2, 1, 7, "um");
Dialog.addNumber("Gaussian blur sigma for DAPI (px)", 1);
Dialog.addChoice("Starting DAPI threshold method", newArray("Default", "Otsu", "Huang", "Li", "Triangle", "IsoData", "Mean"), "Default");
Dialog.addString("Nucleus size filter (um^2)", "0-Infinity");
Dialog.addCheckbox("Size filter is in PIXEL units instead", false);
Dialog.addCheckbox("Exclude nuclei touching the image edge", false);
Dialog.addNumber("Number of background ROIs", 5);
Dialog.addNumber("SD multiplier for the myogenin threshold", 2);
Dialog.addChoice("Colour the nuclei-calls overlay by", newArray("Background rule (your protocol)", "Otsu on nuclear means"), "Background rule (your protocol)");
Dialog.addChoice("Manual adjustments (B&C, DAPI threshold, bg ROIs)", newArray("On every image", "On the first image only, reuse for the rest"), "On every image");
Dialog.show();

mode       = Dialog.getChoice();
myogCh     = Dialog.getNumber();
phalCh     = Dialog.getNumber();
dapiCh     = Dialog.getNumber();
pxSize     = Dialog.getNumber();
fovWidth   = Dialog.getNumber();
pxMode     = Dialog.getChoice();
sigmaDAPI  = Dialog.getNumber();
thrMethod  = Dialog.getChoice();
sizeStr    = Dialog.getString();
sizeInPx   = Dialog.getCheckbox();
exclEdges  = Dialog.getCheckbox();
nBG        = Dialog.getNumber();
sdMult     = Dialog.getNumber();
overlayRule = Dialog.getChoice();
adjustMode = Dialog.getChoice();

if (nBG < 2) exit("Use at least 2 background ROIs (a standard deviation needs 2 or more).");
if (myogCh == dapiCh || myogCh == phalCh || phalCh == dapiCh) exit("Channel numbers must all be different.");

nDone = 0;
if (mode == "Active image only") {
    if (nImages == 0) exit("Open an image first, then run the macro again.");
    id0 = getImageID();
    title0 = getTitle();
    dir0 = getDirectory("image");
    if (dir0 == "") dir0 = getDirectory("Choose the folder where results should be saved");
    outRoot = dir0 + "C2C12_results" + File.separator;
    File.makeDirectory(outRoot);
    ok = processImage(id0, stripExt(title0));
    if (ok) nDone++;
} else {
    inDir = getDirectory("Choose the folder with your images (TIF / JPG / PNG)");
    outRoot = inDir + "C2C12_results" + File.separator;
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
showMessage("C2C12 quantification finished", "Images processed: " + nDone + "\n\nResults folder:\n" + outRoot + "\n\nMaster table: ALL_images_summary.csv");


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

// Otsu threshold of a list of numbers (256-bin histogram)
function otsuThreshold(arr) {
    Array.getStatistics(arr, mn, mx, avg, sdv);
    rng = mx - mn;
    if (rng <= 0) return mn;
    nb = 256;
    hist = newArray(nb);
    for (i = 0; i < arr.length; i++) {
        b = floor((arr[i] - mn) / rng * (nb - 1));
        hist[b] = hist[b] + 1;
    }
    total = arr.length;
    sumAll = 0;
    for (i = 0; i < nb; i++) sumAll = sumAll + i * hist[i];
    wB = 0; sumB = 0; best = -1; bestT = 0;
    for (i = 0; i < nb; i++) {
        wB = wB + hist[i];
        if (wB == 0) continue;
        wF = total - wB;
        if (wF == 0) break;
        sumB = sumB + i * hist[i];
        mB = sumB / wB;
        mF = (sumAll - sumB) / wF;
        between = wB * wF * (mB - mF) * (mB - mF);
        if (between > best) { best = between; bestT = i; }
    }
    return mn + (bestT + 1) / (nb - 1) * rng;
}

function median(arr) {
    a = Array.copy(arr);
    Array.sort(a);
    m = a.length;
    if (m % 2 == 1) return a[(m - 1) / 2];
    return (a[m / 2 - 1] + a[m / 2]) / 2;
}

function processImage(origID, base) {
    outDir = outRoot + base + File.separator;
    File.makeDirectory(outDir);
    bp = outDir + base;
    print("------ " + base);

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
    run("Split Channels");
    selectWindow("C" + myogCh + "-work"); myogID = getImageID();
    selectWindow("C" + phalCh + "-work"); phalID = getImageID();
    selectWindow("C" + dapiCh + "-work"); dapiID = getImageID();

    interactive = (adjustMode == "On every image") || !haveStored;

    selectImage(myogID); getStatistics(aa, myogMean0);       // raw values BEFORE touching the display
    selectImage(dapiID); getStatistics(aa, dapiMean0);

    // display-only auto contrast on every channel
    allIDs = newArray(myogID, dapiID, phalID);
    for (k = 0; k < 3; k++) {
        selectImage(allIDs[k]);
        run("Enhance Contrast", "saturated=0.35");
    }
    if (interactive) {
        selectImage(myogID);
        run("Brightness/Contrast...");
        waitForUser("MYOGENIN display", "Click on the MYOGENIN window, then adjust Brightness/Contrast (Auto or sliders).\nDo NOT click Apply.\n\nClick OK here when you are happy.");
        selectImage(dapiID);
        run("Brightness/Contrast...");
        waitForUser("DAPI display", "Click on the DAPI window, then adjust Brightness/Contrast (Auto or sliders).\nDo NOT click Apply.\n\nClick OK here when you are happy.");
    } else {
        selectImage(myogID); setMinAndMax(stMyogMin, stMyogMax);
        selectImage(dapiID); setMinAndMax(stDapiMin, stDapiMax);
    }
    close("B&C");
    selectImage(myogID); getMinAndMax(myogMin, myogMax);
    selectImage(dapiID); getMinAndMax(dapiMin, dapiMax);

    // safety: raw pixel values must be unchanged ("Apply" would change them)
    selectImage(myogID); getStatistics(aa, myogMean1);
    selectImage(dapiID); getStatistics(aa, dapiMean1);
    if (abs(myogMean1 - myogMean0) > 0.001 || abs(dapiMean1 - dapiMean0) > 0.001)
        exit("Brightness/Contrast 'Apply' was clicked - the pixel data were changed.\nYour original file is untouched. Run the macro again and do NOT click Apply.");

    // ---------- 2. COUNT DAPI NUCLEI ----------
    selectImage(dapiID);
    run("Duplicate...", "title=DAPI_analysis");
    anID = getImageID();
    run("Gaussian Blur...", "sigma=" + sigmaDAPI);
    setAutoThreshold(thrMethod + " dark");
    if (interactive) {
        run("Threshold...");
        waitForUser("DAPI threshold", "Adjust the DAPI threshold in the Threshold window so the red area covers the nuclei\n(start: " + thrMethod + ", dark background). Do NOT click Apply.\n\nClick OK here when you are happy.");
    } else {
        setThreshold(stDapiLo, stDapiHi);
    }
    selectImage(anID);
    getThreshold(thrLo, thrHi);
    thrLoTxt = "NA"; thrHiTxt = "NA";
    if (thrLo != -1) {
        thrLoTxt = d2s(thrLo, 3); thrHiTxt = d2s(thrHi, 3);
        if (!is("binary")) run("Convert to Mask");
    } else if (!is("binary")) {
        exit("No DAPI threshold was set. Run the macro again and adjust the threshold in the Threshold window.");
    }
    close("Threshold");
    run("Watershed");
    saveAs("Tiff", bp + "_DAPI_mask.tif");                 // binary mask (nuclei = 255), image stays open

    run("Set Measurements...", "area mean shape feret's display redirect=None decimal=4");
    run("Clear Results");
    roiManager("Reset");
    opt = "size=" + sizeStr;
    if (sizeInPx) opt = opt + " pixel";
    opt = opt + " circularity=0.00-1.00 show=[Overlay Outlines] display clear";
    if (exclEdges) opt = opt + " exclude";
    opt = opt + " add";
    run("Analyze Particles...", opt);
    n = nResults;
    if (startsWith(getTitle(), "Drawing of")) {
        saveAs("PNG", bp + "_DAPI_outlines.png");
        close();
    }
    if (n == 0) {
        print("   No nuclei detected - check the DAPI threshold.");
        closeID(myogID); closeID(phalID); closeID(dapiID); closeID(anID);
        return false;
    }
    area = newArray(n); circ = newArray(n); feret = newArray(n);
    for (i = 0; i < n; i++) {
        area[i]  = getResult("Area", i);
        circ[i]  = getResult("Circ.", i);
        feret[i] = getResult("Feret", i);
    }
    run("Clear Results");
    if (roiManager("count") != n) exit("ROI Manager count (" + roiManager("count") + ") differs from the particle count (" + n + ").");
    roiManager("Deselect");
    roiManager("Save", bp + "_DAPI_nuclei_ROIs.zip");
    print("   DAPI nuclei: " + n);

    // ---------- 3. MYOGENIN (and DAPI) INTENSITY IN EVERY DAPI NUCLEUS ----------
    run("Set Measurements...", "area mean integrated median display redirect=None decimal=4");
    selectImage(dapiID);
    roiManager("Deselect");
    roiManager("Measure");
    dMean = newArray(n);
    for (i = 0; i < n; i++) dMean[i] = getResult("Mean", i);
    run("Clear Results");

    selectImage(myogID);
    roiManager("Deselect");
    roiManager("Measure");
    myMean = newArray(n); myMed = newArray(n); myInt = newArray(n); myRaw = newArray(n);
    for (i = 0; i < n; i++) {
        myMean[i] = getResult("Mean", i);
        myMed[i]  = getResult("Median", i);
        myInt[i]  = getResult("IntDen", i);
        myRaw[i]  = getResult("RawIntDen", i);
    }
    run("Clear Results");

    // ---------- 4. BACKGROUND CORRECTION / THRESHOLD ----------
    bgM = newArray(nBG);
    bgSrc = "reused from first image";
    drewBG = false;
    if (interactive) {
        drewBG = true;
        bgSrc = "drawn on this image (" + nBG + " ROIs)";
        selectImage(myogID);
        roiManager("Show All without labels");
        k = 0;
        while (k < nBG) {
            selectImage(myogID);
            roiManager("Show All without labels");
            setTool("rectangle");
            run("Select None");
            waitForUser("Background ROI " + (k + 1) + " of " + nBG, "On the MYOGENIN image draw a small rectangle/oval in an area WITHOUT cells,\nthen click OK here.\n(Nuclei outlines are shown so you can avoid them.)");
            selectImage(myogID);
            if (selectionType() == -1) {
                showMessage("No selection found", "Draw the background ROI on the myogenin image first.");
                continue;
            }
            getStatistics(aa, mm);                          // raw myogenin mean inside the ROI
            if (aa > 4 * median(area)) {
                if (getBoolean("This ROI is " + d2s(aa / median(area), 1) + "x the size of a typical nucleus.\nLarge ROIs average out the noise, so the SD becomes too small and the threshold too low.\nDraw a nucleus-sized ROI instead?")) continue;
            }
            roiManager("Add");                              // index n + k
            selectImage(anID);
            roiManager("Select", n + k);                    // apply the same ROI to the nuclei mask
            getStatistics(aMask, mMask);                    // mask mean > 0 means the ROI touches a nucleus
            run("Select None");
            if (mMask > 0) {
                if (getBoolean("This ROI overlaps a nucleus. Draw it again?")) {
                    roiManager("Select", n + k);
                    roiManager("Delete");
                    selectImage(myogID);
                    run("Select None");
                    continue;
                }
            }
            roiManager("Select", n + k);
            roiManager("Rename", "BG_" + (k + 1));
            roiManager("Deselect");
            bgM[k] = mm;
            selectImage(myogID);
            run("Select None");
            k++;
        }
        selectImage(myogID);
        idx = newArray(nBG);
        for (k = 0; k < nBG; k++) idx[k] = n + k;
        roiManager("Select", idx);
        roiManager("Save Selected", bp + "_background_ROIs.zip");
        roiManager("Deselect");
        // background table (X,Y,Width,Height in pixels; Mean = raw myogenin)
        bgTxt = "BG_ROI,X,Y,Width,Height,Mean\n";
        for (k = 0; k < nBG; k++) {
            roiManager("Select", n + k);
            getSelectionBounds(sx, sy, sw, sh);
            bgTxt = bgTxt + (k + 1) + "," + sx + "," + sy + "," + sw + "," + sh + "," + d2s(bgM[k], 4) + "\n";
        }
        roiManager("Deselect");
        File.saveString(bgTxt, bp + "_background_ROIs.csv");
        Array.getStatistics(bgM, bgMin, bgMax, bgMean, bgSD);   // SD = sample SD (n-1)
        nBGused = nBG;
    } else {
        bgMean = stBgMean; bgSD = stBgSD; nBGused = stBgN;
    }
    thr = bgMean + sdMult * bgSD;

    // ---------- 5. CLASSIFY NUCLEI AND CALCULATE ----------
    pos = newArray(n); sub = newArray(n);
    npos = 0;
    for (i = 0; i < n; i++) {
        sub[i] = myMean[i] - bgMean;
        if (myMean[i] > thr) { pos[i] = 1; npos++; }
    }
    pct = 100 * npos / n;

    // supplementary (not in the original protocol): DAPI bleed-through into the myogenin channel
    Fit.doFit("Straight Line", dMean, myMean);
    slope = Fit.p(1);
    r2 = Fit.rSquared;
    Array.getStatistics(dMean, dMin, dMax, dAvg, dSd);
    corr = newArray(n); posC = newArray(n); nposC = 0;
    for (i = 0; i < n; i++) {
        corr[i] = myMean[i] - slope * (dMean[i] - dMin);
        if (corr[i] > thr) { posC[i] = 1; nposC++; }
    }
    pctC = 100 * nposC / n;

    // supplementary rule: Otsu split of the per-nucleus myogenin means (use only if the distribution is clearly two-humped;
    // it is the sensible choice for 8-bit JPG exports whose background is clipped to ~0)
    otsuThr = otsuThreshold(myMean);
    posO = newArray(n); nposO = 0;
    for (i = 0; i < n; i++) {
        if (myMean[i] > otsuThr) { posO[i] = 1; nposO++; }
    }
    pctO = 100 * nposO / n;
    regQC = "ok";
    if (r2 < 0.02) regQC = "CHECK";     // myogenin barely follows DAPI at all: possible channel shift (or no bleed-through)

    Array.getStatistics(myMean, mn1, mx1, myAvg, sd1);
    Array.getStatistics(sub, mn2, mx2, subAvg, sd2);

    // ---------- SAVE: per-nucleus CSV ----------
    qBase = "\"" + base + "\"";          // text in quotes for CSV
    qSrc = "\"" + bgSrc + "\"";
    s = "Image,Nucleus,Area_um2,Circ,Feret_um,DAPI_Mean,Myog_Mean,Myog_Median,Myog_IntDen,Myog_RawIntDen,Myog_Mean_BGsub,Myog_Positive,Myog_Mean_bleedcorr,Myog_Positive_bleedcorr,Myog_Positive_otsu\n";
    for (i = 0; i < n; i++) {
        s = s + qBase + "," + (i + 1) + "," + d2s(area[i], 4) + "," + d2s(circ[i], 4) + "," + d2s(feret[i], 4) + "," + d2s(dMean[i], 4) + ","
              + d2s(myMean[i], 4) + "," + d2s(myMed[i], 4) + "," + d2s(myInt[i], 4) + "," + d2s(myRaw[i], 4) + "," + d2s(sub[i], 4) + "," + pos[i] + ","
              + d2s(corr[i], 4) + "," + posC[i] + "," + posO[i] + "\n";
    }
    File.saveString(s, bp + "_myogenin_per_nucleus.csv");

    // ---------- SAVE: summary CSV (one row) ----------
    getDateAndTime(yr, mo, dw, dm, hr, mi, sc, ms);
    stamp = "" + yr + "-" + IJ.pad(mo + 1, 2) + "-" + IJ.pad(dm, 2) + " " + IJ.pad(hr, 2) + ":" + IJ.pad(mi, 2);
    excl = 0; if (exclEdges) excl = 1;
    hdr = "Image,Method,Date,Pixel_um,DAPI_threshold_method,DAPI_threshold_lower,Gaussian_sigma,Size_filter,Exclude_edges,Total_DAPI_nuclei,Median_nuclear_area_um2,Myog_positive_nuclei,Myog_positive_fraction_pct,N_BG_ROIs,BG_mean,BG_SD,SD_multiplier,Threshold,BG_source,Myog_mean_raw_mean,Myog_mean_raw_median,Myog_BGsub_mean,Myog_BGsub_median,DAPI_bleed_slope,Myog_vs_DAPI_R2,Myog_DAPI_shift_dy_dx_px,Registration_QC,Myog_positive_fraction_bleedcorr_pct,Myog_threshold_otsu,Myog_positive_fraction_otsu_pct\n";
    row = qBase + ",Fiji macro," + stamp + "," + d2s(px, 4) + "," + thrMethod + "," + thrLoTxt + "," + sigmaDAPI + "," + sizeStr + "," + excl + "," + n + "," + d2s(median(area), 3) + ","
        + npos + "," + d2s(pct, 3) + "," + nBGused + "," + d2s(bgMean, 3) + "," + d2s(bgSD, 3) + "," + sdMult + "," + d2s(thr, 3) + "," + qSrc + ","
        + d2s(myAvg, 3) + "," + d2s(median(myMean), 3) + "," + d2s(subAvg, 3) + "," + d2s(median(sub), 3) + "," + d2s(slope, 6) + "," + d2s(r2, 4) + ",NA," + regQC + "," + d2s(pctC, 3) + "," + d2s(otsuThr, 3) + "," + d2s(pctO, 3) + "\n";
    File.saveString(hdr + row, bp + "_summary.csv");

    // ---------- SAVE: settings used ----------
    st = "DAPI threshold method (start): " + thrMethod + "\nDAPI threshold lower/upper (on Gaussian-blurred DAPI): " + thrLoTxt + " / " + thrHiTxt
       + "\nGaussian sigma (px): " + sigmaDAPI + "\nParticle size: " + sizeStr + " (pixel units: " + sizeInPx + ")\nExclude edges: " + exclEdges
       + "\nMyogenin display range: " + d2s(myogMin, 2) + " - " + d2s(myogMax, 2) + "\nDAPI display range: " + d2s(dapiMin, 2) + " - " + d2s(dapiMax, 2)
       + "\nPixel size (um): " + px + "\nSource: " + srcType + " (" + w + " x " + h + " px)" + "\nBackground: " + bgSrc + "\nSD multiplier: " + sdMult + "\nOverlay coloured by: " + overlayRule + "\nSD type: sample SD (n-1)\n";
    File.saveString(st, bp + "_settings.txt");

    // ---------- SAVE: overlays ----------
    // clear any leftover selection first - Duplicate would otherwise crop to the last ROI
    selectImage(myogID); run("Select None");
    selectImage(dapiID); run("Select None");
    selectImage(myogID); run("Duplicate...", "title=m8"); setMinAndMax(myogMin, myogMax); run("8-bit");
    selectImage(dapiID); run("Duplicate...", "title=d8"); setMinAndMax(dapiMin, dapiMax); run("8-bit");
    selectImage("m8"); getDimensions(mw, mh, mc, mz, mt);
    selectImage("d8"); getDimensions(dw2, dh2, dc2, dz2, dt2);
    if (mw != dw2 || mh != dh2) exit("Overlay images differ in size (" + mw + "x" + mh + " vs " + dw2 + "x" + dh2 + ").");
    run("Merge Channels...", "c2=m8 c3=d8");                // myogenin green, DAPI blue -> RGB
    rgbID = getImageID();
    saveAs("PNG", bp + "_overlay_DAPI_myogenin.png");
    selectImage(rgbID);

    roiManager("Reset");
    roiManager("Open", bp + "_DAPI_nuclei_ROIs.zip");
    posIdx = newArray(0); negIdx = newArray(0);
    for (i = 0; i < n; i++) {
        flag = pos[i];
        if (overlayRule == "Otsu on nuclear means") flag = posO[i];
        if (flag == 1) posIdx = Array.concat(posIdx, i);
        else negIdx = Array.concat(negIdx, i);
    }
    if (posIdx.length > 0) { roiManager("Select", posIdx); roiManager("Set Color", "yellow"); }
    if (negIdx.length > 0) { roiManager("Select", negIdx); roiManager("Set Color", "cyan"); }
    if (drewBG) {
        roiManager("Open", bp + "_background_ROIs.zip");
        bgIdx = newArray(nBG);
        for (k = 0; k < nBG; k++) bgIdx[k] = n + k;
        roiManager("Select", bgIdx);
        roiManager("Set Color", "red");
    }
    roiManager("Deselect");
    selectImage(rgbID);
    run("From ROI Manager");
    run("Flatten");
    saveAs("PNG", bp + "_overlay_nuclei_calls.png");
    flatID = getImageID();

    // ---------- remember settings for "reuse" mode ----------
    if (interactive && adjustMode != "On every image") {
        if (thrLo == -1) {
            print("   (Threshold was applied inside the Threshold window, so it cannot be reused. Every image will ask again.)");
        } else {
            stDapiLo = thrLo; stDapiHi = thrHi;
            stMyogMin = myogMin; stMyogMax = myogMax; stDapiMin = dapiMin; stDapiMax = dapiMax;
            stBgMean = bgMean; stBgSD = bgSD; stBgN = nBG;
            haveStored = true;
        }
    }

    print("   background " + d2s(bgMean, 2) + " +/- " + d2s(bgSD, 2) + "  -> threshold " + d2s(thr, 2));
    print("   myogenin-positive: " + npos + " / " + n + " = " + d2s(pct, 2) + " %   [bleed-corrected: " + d2s(pctC, 2) + " %]   [Otsu on nuclear means: " + d2s(pctO, 2) + " %, threshold " + d2s(otsuThr, 2) + "]   registration QC: " + regQC);
    if (regQC == "CHECK") print("   !! CHECK: myogenin intensity does not follow the DAPI nuclei at all. Look at the overlay PNG - the channels may be shifted.");

    // ---------- clean up ----------
    closeID(flatID); closeID(rgbID);
    closeID(myogID); closeID(phalID); closeID(dapiID); closeID(anID);
    if (isOpen("m8")) { selectImage("m8"); close(); }
    if (isOpen("d8")) { selectImage("d8"); close(); }
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
