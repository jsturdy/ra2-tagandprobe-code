import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


##                      _              _       
##   ___ ___  _ __  ___| |_ __ _ _ __ | |_ ___ 
##  / __/ _ \| '_ \/ __| __/ _` | '_ \| __/ __|
## | (_| (_) | | | \__ \ || (_| | | | | |_\__ \
##  \___\___/|_| |_|___/\__\__,_|_| |_|\__|___/
##                                              
################################################

datasets = [
    'test'
]


#isMC = False
isMC = True
InputFileDir = "/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/2012/tnp/"
#BaseFileName = datasets[0]+"_tagProbeTree_noHLT"
BaseFileName = "Photon_tagProbeTreePAT.root"
OutputFilePrefix = "efficiency-data-"



################################################
HLTDef       = "probe_passingHLT"
PDFNameSB    = "pdfSignalPlusBackground"
PDFNameG     = "pdfSimpleGaussian"
PDFName      = "pdfSimpleGaussian"
PDFNameGPL   = "gaussPlusLinear"
PDFNameGPLPB = "gaussPlusLinearBkgFloat"
#PDFName = "pdfVBTF"
PDFNameSBHighEt    = "pdfHighEtSignalPlusBackground"
PDFNameGHighEt     = "pdfHighEtSimpleGaussian"
PDFNameGPLHighEt   = "pdfHighEtgaussPlusLinear"
PDFNameGPLPBHighEt = "pdfHighEtgaussPlusLinearBkgFloat"

if isMC:
    #PDFName = ""
    OutputFilePrefix = "efficiency-mc-"
    BaseFileName = "Photon_tagProbeTreePAT_MC.root"
################################################

#specifies the binning of parameters
EfficiencyBinsNPV = cms.PSet(
    probe_npv   = cms.vdouble(0, 1, 5, 9, 11, 13, 15, 19, 23, 25, 27, 29, 31,33,35,37,41,45, 50, 55,60 )
)
EfficiencyBinsEt = cms.PSet(
    probe_et    = cms.vdouble( 20, 80, 100, 120, 500 )
)
EfficiencyBinsEta = cms.PSet(
    probe_eta   = cms.vdouble( -2.5, -2.1, -1.566, -1.4442, -0.9, 0.0, 0.9, 1.4442, 1.566, 2.1, 2.5 )
)
EfficiencyBinsEtEta = cms.PSet(
    probe_et    = cms.vdouble( 75, 80, 100, 120, 200, 300, 1000),
    probe_eta   = cms.vdouble( -2.5, -2.1, -1.566, -1.4442, -0.9, 0.0, 0.9, 1.4442, 1.566, 2.1, 2.5 )
)
EfficiencyBinsEtNPV = cms.PSet(
    probe_npv   = cms.vdouble( 0, 1, 5, 9, 11, 13, 15, 19, 23, 25, 27, 29, 31, 33, 35, 37, 41, 45, 50, 55,60),
    probe_et    = cms.vdouble( 20, 80, 100, 120, 500 )
)
EfficiencyBinsEtaNPV = cms.PSet(
    probe_npv   = cms.vdouble( 0, 1, 5, 9, 11, 13, 15, 19, 23, 25, 27, 29, 31, 33, 35, 37, 41, 45, 50, 55,60),
    probe_eta   = cms.vdouble( -2.5, -2.1, -1.566, -1.4442, -0.9, 0.0, 0.9, 1.4442, 1.566, 2.1, 2.5 )
)
EfficiencyBinsEtEtaNPV = cms.PSet(
    probe_npv   = cms.vdouble( 0, 1, 5, 9, 11, 13, 15, 19, 23, 25, 27, 29, 31, 33, 35, 37, 41, 45, 50, 55,60),
    probe_et    = cms.vdouble( 75, 80, 100, 120, 200, 300, 1000),
    probe_eta   = cms.vdouble( -2.5, -2.1, -1.566, -1.4442, -0.9, 0.0, 0.9, 1.4442, 1.566, 2.1, 2.5 )
)
## for super clusters
EfficiencyBinsSC = cms.PSet(
    probe_sc_et    = cms.vdouble( 75, 80, 100, 120, 200, 300, 1000),
    probe_sc_eta   = cms.vdouble( -2.5, -2.1, -1.566, -1.4442, -0.9, 0.0, 0.9, 1.4442, 1.566, 2.1, 2.5 )
)

#### For data: except for HLT step
EfficiencyBinningSpecificationEt = cms.PSet(
    #specifies what unbinned variables to include in the dataset, the mass is needed for the fit
    UnbinnedVariables = cms.vstring("mass"),
    #specifies the binning of parameters
    BinnedVariables = cms.PSet(EfficiencyBinsEt),
    #first string is the default followed by binRegExp - PDFname pairs
    BinToPDFmap = cms.vstring(PDFName)
)

EfficiencyBinningSpecificationEta = cms.PSet(
    UnbinnedVariables = cms.vstring("mass"),
    BinnedVariables = cms.PSet(EfficiencyBinsEta),
    BinToPDFmap = cms.vstring(PDFName)
)

EfficiencyBinningSpecificationNPV = cms.PSet(
    UnbinnedVariables = cms.vstring("mass"),
    BinnedVariables = cms.PSet(EfficiencyBinsNPV),
    BinToPDFmap = cms.vstring(PDFName)
)

EfficiencyBinningSpecificationEtEta = cms.PSet(
    UnbinnedVariables = cms.vstring("mass"),
    BinnedVariables = cms.PSet(EfficiencyBinsEtEta),
    BinToPDFmap = cms.vstring(PDFName)
)

EfficiencyBinningSpecificationEtNPV = cms.PSet(
    UnbinnedVariables = cms.vstring("mass"),
    BinnedVariables = cms.PSet(EfficiencyBinsEtNPV),
    BinToPDFmap = cms.vstring(PDFName)
)

EfficiencyBinningSpecificationEtaNPV = cms.PSet(
    UnbinnedVariables = cms.vstring("mass"),
    BinnedVariables = cms.PSet(EfficiencyBinsEtaNPV),
    BinToPDFmap = cms.vstring(PDFName)
)

EfficiencyBinningSpecificationEtEtaNPV = cms.PSet(
    UnbinnedVariables = cms.vstring("mass"),
    BinnedVariables = cms.PSet(EfficiencyBinsEtEtaNPV),
    BinToPDFmap = cms.vstring(PDFName)
)
#### For super clusters
EfficiencyBinningSpecificationSC = cms.PSet(
    UnbinnedVariables = cms.vstring("mass"),
    BinnedVariables = cms.PSet(EfficiencyBinsSC),
    BinToPDFmap = cms.vstring(PDFName)
)
EfficiencyBinningSpecificationSCMC = cms.PSet(
    UnbinnedVariables = cms.vstring("mass"),
    BinnedVariables = cms.PSet(EfficiencyBinsSC,mcTrue = cms.vstring("true")),
    BinToPDFmap = cms.vstring()  
)


#### For MC truth: do truth matching
EfficiencyBinningSpecificationMC = cms.PSet(
    UnbinnedVariables = cms.vstring("mass"),
    BinnedVariables = cms.PSet(
    #probe_npv = cms.vdouble( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 30, 50 ),
    #probe_npv   = cms.vdouble( 0, 1, 5, 9, 11, 13, 15, 19, 23, 25, 27, 29, 31, 33, 35, 37, 41, 45, 50, 55,60),
    #probe_et   = cms.vdouble( 25, 30, 35, 40, 45, 50, 60, 75, 200 ),
    probe_eta  = cms.vdouble( -2.5, -2.1, -1.8, -1.566, -1.4442, -1.0, -0.5, 0.0, 0.5, 1.0, 1.4442, 1.566, 1.8, 2.1, 2.5 ),
    mcTrue     = cms.vstring("true")
    ),
    BinToPDFmap = cms.vstring()  
)

##### For HLT step: just do cut & count
#EfficiencyBinningSpecificationHLT = cms.PSet(
#    UnbinnedVariables = cms.vstring("mass"),
#    BinnedVariables = cms.PSet(EfficiencyBins),
#    BinToPDFmap = cms.vstring()  
#)

##########################################################################################
############################################################################################
mcTruthModules = cms.PSet()
##########################################################################################
##########################################################################################


        

############################################################################################
############################################################################################
####### GsfElectron->Id / selection efficiency 
############################################################################################
############################################################################################

process.PhotonToIsoId = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
    # IO parameters:
    InputFileNames = cms.vstring(BaseFileName),
    InputDirectoryName = cms.string("PhotonToIsoId"),
    InputTreeName = cms.string("fitter_tree"),
#    OutputFileName = cms.string(InputFileDir+OutputFilePrefix+BaseFileName+"-PhotonToIsoId_gaussPlusLinear_npvBinned.root"),
    #OutputFileName = cms.string(InputFileDir+OutputFilePrefix+BaseFileName+"-PhotonToIsoId_SimpleGaussian_npvBinned.root"),
    #OutputFileName = cms.string("PhotonToIsoId_SimpleGaussian_npvBinned.root"),
    OutputFileName = cms.string("PhotonToIsoId_SimpleGaussian_npvBinned_MC.root"),
    #numbrer of CPUs to use for fitting
    NumCPU = cms.uint32(8),
    #NumCPU = cms.uint32(1),
    # specifies wether to save the RooWorkspace containing the data for each bin and
    # the pdf object with the initial and final state snapshots
    SaveWorkspace = cms.bool(True),
    floatShapeParameters = cms.bool(True),
    #fixVars = cms.vstring("mean"),
                                                 
    # defines all the real variables of the probes available in the input tree and intended for use in the efficiencies
    Variables = cms.PSet(
        mass         = cms.vstring("Tag-Probe Mass", "60.0", "120.0", "GeV/c^{2}"),
        probe_nJets  = cms.vstring("# Jets", "0.", "10.", "" ),
        probe_npv    = cms.vstring("# PV", "0.", "60.", "" ),
        probe_et     = cms.vstring("Probe E_{T}", "0", "1000", "GeV/c"),
        probe_eta    = cms.vstring("Probe #eta", "-2.5", "2.5", ""),
        probe_sc_et  = cms.vstring("Probe SC E_{T}", "0", "1000", "GeV/c"),
        probe_sc_eta = cms.vstring("Probe SC #eta", "-2.5", "2.5", ""),
        probe_r9     = cms.vstring("Probe R9", "0.", "1.05", ""),
        probe_pfChargedPURel = cms.vstring("Probe pfChargedPURel", "0.", "1.05", ""),
        probe_pfNeutralPURel = cms.vstring("Probe pfNeutralPURel", "0.", "1.05", ""),
        probe_pfGammaPURel   = cms.vstring("Probe pfGammaPURel", "0.", "1.05", ""),
        probe_hadTowOverEm   = cms.vstring("Probe H/E", "0", "0.25", ""),
        probe_sigmaIetaIeta  = cms.vstring("Probe #sigma_{i#eta i#eta}", "0.", "0.05.", "" )
    ),

    # defines all the discrete variables of the probes available in the input tree and intended for use in the efficiency calculations
    Categories = cms.PSet(
        mcTrue = cms.vstring("MC true", "dummy[true=1,false=0]"),
#        probe_passing_pixel_match                               = cms.vstring("probe_passing_pixel_match", "dummy[pass=1,fail=0]"), 
        probe_passingId     = cms.vstring("probe_passingId",     "dummy[pass=1,fail=0]"), 
        probe_passingHLT    = cms.vstring("probe_passingHLT",    "dummy[pass=1,fail=0]"), 
        probe_passingIso    = cms.vstring("probe_passingIso",    "dummy[pass=1,fail=0]"), 
        probe_passingId_iso = cms.vstring("probe_passingId_iso", "dummy[pass=1,fail=0]"), 
    ),
    # defines all the PDFs that will be available for the efficiency calculations; uses RooFit's "factory" syntax;
    # each pdf needs to define "signal", "backgroundPass", "backgroundFail" pdfs, "efficiency[0.9,0,1]" and "signalFractionInPassing[0.9]" are used for initial values  
    PDFs = cms.PSet(
        pdfSignalPlusBackground = cms.vstring(
            #"CBExGaussShape::signalRes(mass, mean[2.0946e-01], sigma[8.5695e-04],alpha[3.8296e-04], n[6.7489e+00], sigma_2[2.5849e+00], frac[6.5704e-01])",  ### the signal function goes here
            ### signal resolution for "pass" sample
            "CBExGaussShape::signalResPass(mass, meanP[0., -50., 50.],         sigmaP[8.5695e-04, 0., 3.], alphaP[3.8296e-04, -5., 5.], nP[6.7489e+00, 0., 10000000.], sigmaP_2[2.5849e+00, -5., 5.], fracP[6.5704e-01, 0., 1.])",
            ### signal resolution for "fail" sample
            "CBExGaussShape::signalResFail(mass, meanF[2.0946e-01, -50., 50.], sigmaF[8.5695e-04, 0., 5.], alphaF[3.8296e-04, -5., 5.], nF[6.7489e+00, 0., 10000000.], sigmaF_2[2.5849e+00, -5., 5.], fracF[6.5704e-01, 0., 1.])",
            ### NLO line shape
            "ZGeneratorLineShape::signalPhy(mass)",
            "RooCMSShape::backgroundPass(mass, alphaPass[60.,50.,70.], betaPass[0.001, 0.,0.1], betaPass, peakPass[91.1876, 91.1855, 91.1897])",
            "RooCMSShape::backgroundFail(mass, alphaFail[60.,50.,70.], betaFail[0.001, 0.,0.1], betaFail, peakFail[91.1876, 91.1855, 91.1897])",
            "FCONV::signalPass(mass, signalPhy, signalResPass)",
            "FCONV::signalFail(mass, signalPhy, signalResFail)",     
            "efficiency[0.9, 0., 1.]",
            "signalFractionInPassing[0.9, 0., 1.]"
        ),
        pdfSimpleGaussian = cms.vstring(
            #91.2, 89.0, 93.0
            "Gaussian::signal(mass, mean[91.1876, 91.1855, 91.1897], sigma[2.3, 0.5, 10.0])",
            "RooExponential::backgroundPass(mass, cPass[-0.02,-5,0])",
            "RooExponential::backgroundFail(mass, cFail[-0.02,-5,0])",
            "efficiency[0.9, 0., 1.]",
            "signalFractionInPassing[0.9, 0., 1.]"
        ),
        gaussPlusLinear = cms.vstring(
            "ZGeneratorLineShape::signalPhy(mass)",
            "RooCBShape::signalResPass(mass, meanPass[0.0, -1.0, 1.0], widthPass[1.8, 1.0, 3.0], alphaPass[1.4], nPass[138.6])",
            "RooCBShape::signalResFail(mass, meanFail[0.0, -2.0, 2.0], widthFail[1.8, 1.0, 3.0], alphaFail[1.4], nFail[138.6])",
            "FCONV::signalPass(mass, signalPhy, signalResPass)",
            "FCONV::signalFail(mass, signalPhy, signalResFail)",
            "RooCMSShape::backgroundPass(mass, aPass[79.8], bPass[0.081], gPass[0.02], peakPass[91.1876, 91.1855, 91.1897])",
            "RooCMSShape::backgroundFail(mass, aFail[79.8], bFail[0.081], gFail[0.02], peakFail[91.1876, 91.1855, 91.1897])",
            "efficiency[0.9, 0., 1.]",
            "signalFractionInPassing[0.9, 0., 1.]"
        ),
        gaussPlusLinearBkgFloat = cms.vstring(
            "ZGeneratorLineShape::signalPhy(mass)",
            "RooCBShape::signalResPass(mass, meanPass[0.0, -1.0, 1.0], widthPass[1.8, 1.0, 3.0], alphaPass[1.4], nPass[138.6])",
            "RooCBShape::signalResFail(mass, meanFail[0.0, -1.0, 1.0], widthFail[1.8, 1.0, 3.0], alphaFail[1.4], nFail[138.6])",
            "FCONV::signalPass(mass, signalPhy, signalResPass)",
            "FCONV::signalFail(mass, signalPhy, signalResFail)",
            ##     "RooCMSShape::backgroundPass(mass, aPass[79.8, 50.0, 90.0], bPass[0.17, 0.16, 0.18], gPass[0.2, 0.1, 0.5], peakPass[91.2, 80.0, 100.0])",
            ##     "RooCMSShape::backgroundFail(mass, aFail[79.8, 50.0, 90.0], bFail[0.081, 0.05, 0.5], gFail[0.02, 0.002, 0.2], peakFail[91.2, 80.0, 100.0])",
            "RooCMSShape::backgroundPass(mass, aPass[79.8], bPass[0.081], gPass[0.02], peakPass[91.2, 50.0, 120.0])",
            "RooCMSShape::backgroundFail(mass, aFail[79.8], bFail[0.081], gFail[0.02], peakFail[91.2, 50.0, 120.0])",
            "efficiency[0.9, 0., 1.]",
            "signalFractionInPassing[0.9, 0., 1.]"
        ),
        pdfVBTF = cms.vstring(
            ### signal resolution for "pass" sample
            "CBExGaussShape::signalResPass(mass, meanPass[2.0946e-01, -5., 5.], sigmaPass[8.5695e-04, 0., 5.],alphaPass[3.8296e-04], nPass[6.7489e+00], sigma_2Pass[2.5849e+00,0.0,5.0], fracPass[6.5704e-01])",
            ### signal resolution for "fail" sample
            "CBExGaussShape::signalResFail(mass, meanFail[2.0946e-01, -5., 5.], sigmaFail[8.5695e-04, 0., 5.],alphaFail[3.8296e-04], nFail[6.7489e+00], sigma_2Fail[2.5849e+00,0.0,5.0], fracFail[6.5704e-01])",
            "ZGeneratorLineShape::signalPhy(mass)",
            "RooExponential::backgroundPass(mass, cPass[-0.02, -5, 0])",
            "RooExponential::backgroundFail(mass, cFail[-0.02, -5, 0])",
            "FCONV::signalPass(mass, signalPhy, signalResPass)",
            "FCONV::signalFail(mass, signalPhy, signalResFail)",
            "efficiency[0.9, 0., 1.]",
            "signalFractionInPassing[1.0, 0., 1.]"
        ),
        pdfHighEtSignalPlusBackground = cms.vstring(
            ### signal resolution for "pass" sample
            "CBExGaussShape::signalResPass(mass, meanPass[2.0946e-01, -50., 50.], sigmaPass[8.5695e-04, 0., 5.],alphaPass[3.8296e-04, 0., 5.], nPass[1000., 0., 1000000.], sigma_2Pass[2.5849e+00, 0., 5.], fracPass[6.5704e-01, 0., 1.])",
            ### signal resolution for "fail" sample
            "CBExGaussShape::signalResFail(mass, meanFail[2.0946e-01, -50., 50.], sigmaFail[8.5695e-04, 0., 5.],alphaFail[3.8296e-04, 0., 5.], nFail[10., 0., 1000000.],   sigma_2Fail[2.5849e+00, 0., 5.], fracFail[6.5704e-01, 0., 1.])",
            "ZGeneratorLineShape::signalPhy(mass)",
            "RooExponential::backgroundPass(mass, cPass[-0.02, -5, 0])",
            "RooExponential::backgroundFail(mass, cFail[-0.02, -5, 0])",
            "FCONV::signalPass(mass, signalPhy, signalResPass)",
            "FCONV::signalFail(mass, signalPhy, signalResFail)",
            "efficiency[0.9, 0., 1.]",
            "signalFractionInPassing[1.0, 0., 1.]"
        ),
        pdfHighEtSimpleGaussian = cms.vstring(
            "Gaussian::signal(mass, mean[91.1876, 91.1855, 91.1897], sigma[2.3, 0.5, 10.0])",
            "RooExponential::backgroundPass(mass, cPass[-0.02,-5,0])",
            "RooExponential::backgroundFail(mass, cFail[-0.02,-5,0])",
            "efficiency[0.9, 0., 1.]",
            "signalFractionInPassing[0.9, 0., 1.]"
        ),
        pdfHighEtgaussPlusLinear = cms.vstring(
            "ZGeneratorLineShape::signalPhy(mass)",
            "RooCBShape::signalResPass(mass, meanPass[0.0, -1.0, 1.0], widthPass[1.8, 1.0, 3.0], alphaPass[1.4], nPass[138.6, 0., 10000000])",
            "RooCBShape::signalResFail(mass, meanFail[0.0, -2.0, 2.0], widthFail[1.8, 1.0, 3.0], alphaFail[1.4], nFail[138.6, 0., 10000000])",
            "FCONV::signalPass(mass, signalPhy, signalResPass)",
            "FCONV::signalFail(mass, signalPhy, signalResFail)",
            "RooCMSShape::backgroundPass(mass, aPass[79.8], bPass[0.081], gPass[0.02], peakPass[91.1876, 91.1855, 91.1897])",
            "RooCMSShape::backgroundFail(mass, aFail[79.8], bFail[0.081], gFail[0.02], peakFail[91.1876, 91.1855, 91.1897])",
            "efficiency[0.9, 0., 1.]",
            "signalFractionInPassing[0.9, 0., 1.]"
        ),
        pdfHighEtgaussPlusLinearBkgFloat = cms.vstring(
            "ZGeneratorLineShape::signalPhy(mass)",
            "RooCBShape::signalResPass(mass, meanPass[0.0, -1.0, 1.0], widthPass[1.8, 1.0, 3.0], alphaPass[1.4], nPass[138.6])",
            "RooCBShape::signalResFail(mass, meanFail[0.0, -1.0, 1.0], widthFail[1.8, 1.0, 3.0], alphaFail[1.4], nFail[138.6])",
            "FCONV::signalPass(mass, signalPhy, signalResPass)",
            "FCONV::signalFail(mass, signalPhy, signalResFail)",
            ##     "RooCMSShape::backgroundPass(mass, aPass[79.8, 50.0, 90.0], bPass[0.17, 0.16, 0.18], gPass[0.2, 0.1, 0.5], peakPass[91.2, 80.0, 100.0])",
            ##     "RooCMSShape::backgroundFail(mass, aFail[79.8, 50.0, 90.0], bFail[0.081, 0.05, 0.5], gFail[0.02, 0.002, 0.2], peakFail[91.2, 80.0, 100.0])",
            "RooCMSShape::backgroundPass(mass, aPass[79.8], bPass[0.081], gPass[0.02], peakPass[91.1876, 50.0, 120.0])",
            "RooCMSShape::backgroundFail(mass, aFail[79.8], bFail[0.081], gFail[0.02], peakFail[91.1876, 50.0, 120.0])",
            "efficiency[0.9, 0., 1.]",
            "signalFractionInPassing[0.9, 0., 1.]"
        ),
     ),

    # defines a set of efficiency calculations, what PDF to use for fitting and how to bin the data;
    # there will be a separate output directory for each calculation that includes a simultaneous fit, side band subtraction and counting. 
    Efficiencies = cms.PSet(
        RA2IdIso_vs_npv = cms.PSet(
            EfficiencyBinningSpecificationNPV,    
            EfficiencyCategoryAndState = cms.vstring("probe_passingId_iso","pass")
        ),
        RA2IdIso_vs_et = cms.PSet(
            EfficiencyBinningSpecificationEt,    
            EfficiencyCategoryAndState = cms.vstring("probe_passingId_iso","pass")
        ),
        RA2IdIso_vs_eta = cms.PSet(
            EfficiencyBinningSpecificationEta,    
            EfficiencyCategoryAndState = cms.vstring("probe_passingId_iso","pass")
        ),
        RA2IdIso_vs_etanpv = cms.PSet(
            EfficiencyBinningSpecificationEtaNPV,    
            EfficiencyCategoryAndState = cms.vstring("probe_passingId_iso","pass")
        ),
        RA2IdIso_vs_et_eta = cms.PSet(
            EfficiencyBinningSpecificationEtEta,    
            EfficiencyCategoryAndState = cms.vstring("probe_passingId_iso","pass")
        ),
        RA2IdIso_vs_et_npv = cms.PSet(
            EfficiencyBinningSpecificationEtNPV,    
            EfficiencyCategoryAndState = cms.vstring("probe_passingId_iso","pass")
        ),
        RA2IdIso_vs_et_eta_npv = cms.PSet(
            EfficiencyBinningSpecificationEtEtaNPV,    
            EfficiencyCategoryAndState = cms.vstring("probe_passingId_iso","pass")
        ),
############################################################################################
############################################################################################
    )
)


#process.PhotonToIsoId.Efficiencies.RA2PixelMatch_vs_npv.BinToPDFmap = cms.vstring(PDFNameGPL)

#process.PhotonToIsoId.Efficiencies.RA2IdIso_vs_et.BinToPDFmap = cms.vstring(PDFNameGPL)
#process.PhotonToIsoId.Efficiencies.RA2IdIso_vs_npv.BinToPDFmap = cms.vstring(PDFNameGPL)
#process.PhotonToIsoId.Efficiencies.RA2IdIso_vs_eta.BinToPDFmap = cms.vstring(PDFNameGPL)

#process.PhotonToIsoId.Efficiencies.RA2IdIso_vs_etanpv.BinToPDFmap     = cms.vstring(PDFNameGPL)
#process.PhotonToIsoId.Efficiencies.RA2IdIso_vs_et_npv.BinToPDFmap     = cms.vstring(PDFNameGPL)
#process.PhotonToIsoId.Efficiencies.RA2IdIso_vs_et_eta.BinToPDFmap     = cms.vstring(PDFNameGPL)
#process.PhotonToIsoId.Efficiencies.RA2IdIso_vs_et_eta_npv.BinToPDFmap = cms.vstring(PDFNameGPL)

###########Signal + Bkgd PDF
process.PhotonToIsoIdSigPlusBkg = process.PhotonToIsoId.clone()
process.PhotonToIsoIdSigPlusBkg.OutputFileName     = cms.string("PhotonToIsoId_gaussPlusLinPlusBkg_unBinned.root")
if isMC:
    process.PhotonToIsoIdSigPlusBkg.OutputFileName     = cms.string("PhotonToIsoId_gaussPlusLinPlusBkg_unBinned_MC.root")

process.PhotonToIsoIdSigPlusBkg.Efficiencies.RA2IdIso_vs_et.BinToPDFmap = cms.vstring(PDFNameGPL)
process.PhotonToIsoIdSigPlusBkg.Efficiencies.RA2IdIso_vs_npv.BinToPDFmap = cms.vstring(PDFNameGPL)
process.PhotonToIsoIdSigPlusBkg.Efficiencies.RA2IdIso_vs_eta.BinToPDFmap = cms.vstring(PDFNameGPL)

##process.PhotonToIsoIdSigPlusBkg.Efficiencies.RA2IdIso_vs_etanpv.BinToPDFmap     = cms.vstring(PDFNameGPL)
##process.PhotonToIsoIdSigPlusBkg.Efficiencies.RA2IdIso_vs_et_npv.BinToPDFmap     = cms.vstring(PDFNameGPL)
##process.PhotonToIsoIdSigPlusBkg.Efficiencies.RA2IdIso_vs_et_eta.BinToPDFmap     = cms.vstring(PDFNameGPL)
##process.PhotonToIsoIdSigPlusBkg.Efficiencies.RA2IdIso_vs_et_eta_npv.BinToPDFmap = cms.vstring(PDFNameGPL)
#
#
##process.PhotonToIsoIdHighEt                    = process.PhotonToIsoId.clone()
##process.PhotonToIsoIdHighEt.InputDirectoryName = cms.string("PhotonToIsoIdHighEt")
##process.PhotonToIsoIdHighEt.OutputFileName     = cms.string(InputFileDir+OutputFilePrefix+BaseFileName+"-PhotonToIsoIdHighEt_gaussPlusLinear_unBinned.root")
##process.PhotonToIsoIdHighEt.Variables.mass     = cms.vstring("Tag-Probe Mass", "140.0", "500.0", "GeV/c^{2}")
###process.PhotonToIsoIdHighEt.Efficiencies.RA2IdIso.BinToPDFmap   = cms.vstring(PDFNameGPLHighEt)


process.fit = cms.Path(
    process.PhotonToIsoId
#    +process.PhotonToIsoIdSigPlusBkg
    #+process.PhotonToIsoIdHighEt 
    #+process.PhotonToIsoIdHighEtSigPlusBkg 
    )

####-- Dump config ------------------------------------------------------------
#file = open('tnp_analyzer_debugging.py','w')
#file.write(str(process.dumpPython()))
#file.close()

