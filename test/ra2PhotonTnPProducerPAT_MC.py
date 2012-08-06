import FWCore.ParameterSet.Config as cms

##                      _              _       
##   ___ ___  _ __  ___| |_ __ _ _ __ | |_ ___ 
##  / __/ _ \| '_ \/ __| __/ _` | '_ \| __/ __|
## | (_| (_) | | | \__ \ || (_| | | | | |_\__ \
##  \___\___/|_| |_|___/\__\__,_|_| |_|\__|___/
##                                              
########################
MC_FLAG = False

GLOBAL_TAG = 'START52_V11C::All'
OUTPUT_FILE_NAME = "Photon_tagProbeTreePAT_MC.root"

HLTPath1 = "HLT_Photon150_v3"
HLTPath2 = "HLT_Photon20_CaloIdVL_IsoL_v15"
HLTPath3 = "HLT_Photon90_CaloIdVL_IsoL_v14"
#InputTagProcess = "REDIGI36X"
InputTagProcess = "HLT"
RECOProcess = "RECO"
JET_COLL = "ak5PFJets"
JET_CUTS = "abs(eta)<2.6 && chargedHadronEnergyFraction>0 && electronEnergyFraction<0.1 && nConstituents>1 && neutralHadronEnergyFraction<0.99 && neutralEmEnergyFraction<0.99 && pt>15.0" 
ELECTRON_ET_CUT_MIN = 20.0
#ELECTRON_COLL = "gsfElectrons"
ELECTRON_COLL = "patElectronsUserData"
ELECTRON_CUTS = "ecalDrivenSeed==1 && (abs(superCluster.eta)<2.5) && !(1.4442<abs(superCluster.eta)<1.566) && (ecalEnergy*sin(superClusterPosition.theta)>" + str(ELECTRON_ET_CUT_MIN) + ")"
ELECTRON_MED_EB = """(isEB &&
                     (abs(deltaEtaSuperClusterTrackAtVtx) < 0.004) &&
                     (abs(deltaPhiSuperClusterTrackAtVtx) < 0.06) &&
                     (sigmaIetaIeta < 0.01) &&
                     (hcalOverEcalBc < 0.12) &&
                     (gsfTrack.dxy < 0.02) &&
                     (gsfTrack.dz < 0.1) &&
                     (abs((1/ecalEnergy) - (1/(ecalEnergy/eSuperClusterOverP))) < 0.05) &&
                     (userFloat('pfIsoPURel')/pt < 0.15) &&
                     (userInt('passConvVeto') > 0) &&
                     (gsfTrack.trackerExpectedHitsInner.numberOfHits <= 1))"""
#                     Conversion rejection: vertex fit probability < 1e-6 &&

ELECTRON_MED_EE = """(isEE &&
                     (abs(deltaEtaSuperClusterTrackAtVtx) < 0.007) &&
                     (abs(deltaPhiSuperClusterTrackAtVtx) < 0.03) &&
                     (sigmaIetaIeta < 0.03) &&
                     (hcalOverEcalBc < 0.10) &&
                     (gsfTrack.dxy < 0.02) &&
                     (gsfTrack.dz < 0.1) &&
                     (abs((1/ecalEnergy) - (1/(ecalEnergy/eSuperClusterOverP))) < 0.05) &&
                     (userFloat('pfIsoPURel')/pt < 0.15) &&
                     (userInt('passConvVeto') > 0) &&
                     (gsfTrack.trackerExpectedHitsInner.numberOfHits <= 1))"""
#                     Conversion rejection: vertex fit probability < 1e-6 &&


##    ___            _           _      
##   |_ _|_ __   ___| |_   _  __| | ___ 
##    | || '_ \ / __| | | | |/ _` |/ _ \
##    | || | | | (__| | |_| | (_| |  __/
##   |___|_| |_|\___|_|\__,_|\__,_|\___|

process = cms.Process("TagProbe")
#stuff needed for prescales
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = GLOBAL_TAG
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

#process.load("Configuration.StandardSequences.Reconstruction_cff")
#process.load("CondCore.DBCommon.CondDBSetup_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

##   ____             _ ____                           
##  |  _ \ ___   ___ | / ___|  ___  _   _ _ __ ___ ___ 
##  | |_) / _ \ / _ \| \___ \ / _ \| | | | '__/ __/ _ \
##  |  __/ (_) | (_) | |___) | (_) | |_| | | | (_|  __/
##  |_|   \___/ \___/|_|____/ \___/ \__,_|_|  \___\___|
##  

#readFiles = cms.untracked.vstring('/store/data/Run2012B/PhotonHad/AOD/PromptReco-v1/000/194/314/3C3C88B7-12A2-E111-93FE-001D09F290BF.root')
#readFiles = cms.untracked.vstring('/store/data/Run2012B/DoubleElectron/AOD/PromptReco-v1/000/194/314/A4EEF8F2-0AA2-E111-AD66-00237DDBE0E2.root')
#readFiles = cms.untracked.vstring('/store/data/Run2012B/SingleElectron/AOD/PromptReco-v1/000/194/314/5828CF45-F9A1-E111-A38E-003048D2C01E.root')
readFiles = cms.untracked.vstring('/store/mc/Summer12/DYToEE_M_120_TuneZ2star_8TeV_pythia6/AODSIM/PU_S7_START52_V9-v1/0000/8200EF9B-0AA0-E111-9E58-003048FFCB6A.root')
if MC_FLAG:
    readFiles = cms.untracked.vstring('/store/mc/Summer12/DYToEE_M_120_TuneZ2star_8TeV_pythia6/AODSIM/PU_S7_START52_V9-v1/0000/8200EF9B-0AA0-E111-9E58-003048FFCB6A.root')

process.source = cms.Source("PoolSource", 
                            fileNames = readFiles
                            )
readFiles.extend([
    #FILES
    ])

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )    
process.source.inputCommands = cms.untracked.vstring("keep *","drop *_MEtoEDMConverter_*_*")

##PAT steps for electrons and photons

from RecoJets.JetProducers.kt4PFJets_cfi import *
process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFPhotonIso
###electrons
process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
process.electronPFIsolationDepositsSequencePFIso = cms.Sequence(process.elPFIsoDepositChargedPFIso+
                                                                process.elPFIsoDepositChargedAllPFIso+
                                                                process.elPFIsoDepositGammaPFIso+
                                                                process.elPFIsoDepositNeutralPFIso+
                                                                process.elPFIsoDepositPUPFIso)
process.pfElectronIsolationSequencePFIso = cms.Sequence(process.electronPFIsolationDepositsSequencePFIso+
                                                        process.elPFIsoValueCharged03PFIdPFIso+
                                                        process.elPFIsoValueChargedAll03PFIdPFIso+
                                                        process.elPFIsoValueGamma03PFIdPFIso+
                                                        process.elPFIsoValueNeutral03PFIdPFIso+
                                                        process.elPFIsoValuePU03PFIdPFIso+
                                                        process.elPFIsoValueCharged03NoPFIdPFIso+
                                                        process.elPFIsoValueChargedAll03NoPFIdPFIso+
                                                        process.elPFIsoValueGamma03NoPFIdPFIso+
                                                        process.elPFIsoValueNeutral03NoPFIdPFIso+
                                                        process.elPFIsoValuePU03NoPFIdPFIso)

process.load('PhysicsTools.PatAlgos.producersLayer1.electronProducer_cfi')
process.patElectrons.userIsolation = cms.PSet(
    user = cms.VPSet(
        cms.PSet( src = cms.InputTag("elPFIsoValueCharged03PFIdPFIso")),
        cms.PSet( src = cms.InputTag("elPFIsoValueChargedAll03PFIdPFIso")),
        cms.PSet( src = cms.InputTag("elPFIsoValueNeutral03PFIdPFIso")),
        cms.PSet( src = cms.InputTag("elPFIsoValueGamma03PFIdPFIso")),
        cms.PSet( src = cms.InputTag("elPFIsoValuePU03PFIdPFIso"))
    )
)
#process.patElectrons.userData.userFloats = cms.PSet(
#    src = cms.VInputTag(cms.InputTag("kt6PFJetsForIsolation","rho"))
#)
process.patElectrons.addGenMatch = cms.bool(False)
process.patElectrons.embedGenMatch = cms.bool(False)
#process.patElectrons.genParticleMatch = cms.InputTag("")
process.patElectrons.isoDeposits = cms.PSet(
    pfNeutralHadrons   = cms.InputTag("elPFIsoDepositNeutralPFIso"),
    pfChargedAll       = cms.InputTag("elPFIsoDepositChargedAllPFIso"),
    pfPUChargedHadrons = cms.InputTag("elPFIsoDepositPUPFIso"),
    pfPhotons          = cms.InputTag("elPFIsoDepositGammaPFIso"),
    pfChargedHadrons   = cms.InputTag("elPFIsoDepositChargedPFIso")
    )
process.patElectrons.isolationValues = cms.PSet(
    pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU03PFIdPFIso"),
    pfChargedHadrons   = cms.InputTag("elPFIsoValueCharged03PFIdPFIso"),
    pfNeutralHadrons   = cms.InputTag("elPFIsoValueNeutral03PFIdPFIso"),
    pfPhotons          = cms.InputTag("elPFIsoValueGamma03PFIdPFIso"),
    pfChargedAll       = cms.InputTag("elPFIsoValueChargedAll03PFIdPFIso")
    )
process.patElectrons.isolationValuesNoPFId = cms.PSet(
    pfPUChargedHadrons = cms.InputTag("elPFIsoValuePU03NoPFIdPFIso"),
    pfChargedHadrons   = cms.InputTag("elPFIsoValueCharged03NoPFIdPFIso"),
    pfNeutralHadrons   = cms.InputTag("elPFIsoValueNeutral03NoPFIdPFIso"),
    pfPhotons          = cms.InputTag("elPFIsoValueGamma03NoPFIdPFIso"),
    pfChargedAll       = cms.InputTag("elPFIsoValueChargedAll03NoPFIdPFIso")
    )

process.eleIsoSequence = cms.Sequence(process.pfElectronIsolationSequencePFIso)
process.makePatElectrons = cms.Sequence(process.eleIsoSequence+
                                        process.patElectrons)

###photons
process.phoIsoSequence = setupPFPhotonIso(process, 'photons')
process.photonPFIsolationDepositsSequencePFIso = cms.Sequence(process.phPFIsoDepositChargedPFIso+
                                                              process.phPFIsoDepositChargedAllPFIso+
                                                              process.phPFIsoDepositGammaPFIso+
                                                              process.phPFIsoDepositNeutralPFIso+
                                                              process.phPFIsoDepositPUPFIso)
process.pfPhotonIsolationSequencePFIso = cms.Sequence(process.photonPFIsolationDepositsSequencePFIso+
                                                      process.phPFIsoValueCharged03PFIdPFIso+
                                                      process.phPFIsoValueChargedAll03PFIdPFIso+
                                                      process.phPFIsoValueGamma03PFIdPFIso+
                                                      process.phPFIsoValueNeutral03PFIdPFIso+
                                                      process.phPFIsoValuePU03PFIdPFIso)

process.load('PhysicsTools.PatAlgos.producersLayer1.photonProducer_cfi')
process.patPhotons.addGenMatch = cms.bool(False)
process.patPhotons.embedGenMatch = cms.bool(False)
process.patPhotons.isoDeposits = cms.PSet(
    user = cms.VInputTag(
        cms.InputTag("phPFIsoDepositChargedPFIso"),
        cms.InputTag("phPFIsoDepositChargedAllPFIso"),
        cms.InputTag("phPFIsoDepositNeutralPFIso"),
        cms.InputTag("phPFIsoDepositGammaPFIso"),
        cms.InputTag("phPFIsoDepositPUPFIso")
    ),
)
process.patPhotons.userIsolation = cms.PSet(
    user = cms.VPSet(
        cms.PSet( src = cms.InputTag("phPFIsoValueCharged03PFIdPFIso")),
        cms.PSet( src = cms.InputTag("phPFIsoValueChargedAll03PFIdPFIso")),
        cms.PSet( src = cms.InputTag("phPFIsoValueNeutral03PFIdPFIso")),
        cms.PSet( src = cms.InputTag("phPFIsoValueGamma03PFIdPFIso")),
        cms.PSet( src = cms.InputTag("phPFIsoValuePU03PFIdPFIso"))
    )
)
#process.patPhotons.userData.userFloats = cms.PSet(
#    src = cms.VInputTag(cms.InputTag("kt6PFJetsForIsolation","rho"))
#)
process.phoIsoSequence = cms.Sequence(process.pfPhotonIsolationSequencePFIso)
process.makePatPhotons = cms.Sequence(process.phoIsoSequence+
                                      process.patPhotons)


process.load('ZInvisibleBkgds.Photons.addphotonuserdata_cfi')
process.load('ZInvisibleBkgds.Photons.addelectronuserdata_cfi')

from ZInvisibleBkgds.Photons.addelectronuserdata_cfi import *
process.patElectronsUser1 = addelectronuserdata1.clone()
process.patElectronsUser1.electronLabel = cms.InputTag("patElectrons")
process.patElectronsUserData = addelectronuserdata2.clone()
process.patElectronsUserData.electronLabel = cms.InputTag("patElectronsUser1")

from ZInvisibleBkgds.Photons.addphotonuserdata_cfi import *
process.patPhotonsUser1 = addphotonuserdata1.clone()
process.patPhotonsUser1.photonLabel = cms.InputTag("patPhotons")
process.patPhotonsUserData = addphotonuserdata2.clone()
process.patPhotonsUserData.photonLabel = cms.InputTag("patPhotonsUser1")

process.patObjectSequence = cms.Sequence(
    process.kt6PFJetsForIsolation+
    process.pfParticleSelectionSequence+
    process.makePatElectrons+
    process.patElectronsUser1+
    process.patElectronsUserData+
    process.makePatPhotons+
    process.patPhotonsUser1+
    process.patPhotonsUserData
)
##  ____  _           _                                _               
## |  _ \| |__   ___ | |_ ___  _ __    _ __  _ __ ___ | |__   ___  ___ 
## | |_) | '_ \ / _ \| __/ _ \| '_ \  | '_ \| '__/ _ \| '_ \ / _ \/ __|
## |  __/| | | | (_) | || (_) | | | | | |_) | | | (_) | |_) |  __/\__ \
## |_|   |_| |_|\___/ \__\___/|_| |_| | .__/|_|  \___/|_.__/ \___||___/
##                                    |_|                              

#basic probe photon selection
#keep EB and EE efficiencies separate
#loose track match requirement cuts down on non-Z-electron background
process.probePhotons = cms.EDProducer("TrackMatchedPATPhotonProducer",
    #srcObject = cms.InputTag("photons", "", RECOProcess),
    #srcObject = cms.InputTag("patPhotonsUserData", "", RECOProcess),
    srcObject = cms.InputTag("patPhotonsUserData"),
    srcObjectsToMatch = cms.VInputTag(cms.InputTag("generalTracks")),
    srcObjectSelection = cms.string("et>20.0 && abs(eta)<2.5"),
    srcObjectsToMatchSelection = cms.string('pt > 20.0 && quality("highPurity")'),  
    deltaRMax = cms.double(0.3)
    )

##     ___           _       _   _             
##    |_ _|___  ___ | | __ _| |_(_) ___  _ __  
##     | |/ __|/ _ \| |/ _` | __| |/ _ \| '_ \ 
##     | |\__ \ (_) | | (_| | |_| | (_) | | | |
##    |___|___/\___/|_|\__,_|\__|_|\___/|_| |_|

                                         
#  Isolation ################
#Particle based ID
process.photonIsolation = cms.EDFilter("PATPhotonRefSelector",
    src = cms.InputTag("probePhotons"),
    cut = cms.string(
    "(userFloat('pfChargedPURel') < userFloat('pfChargedMediumCut'))"
    " && (userFloat('pfNeutralPURel') < userFloat('pfNeutralMediumCut'))"
    " && (userFloat('pfGammaPURel')   < userFloat('pfGammaMediumCut'))"
    )
)

##  ____  _           _                ___    _ 
## |  _ \| |__   ___ | |_ ___  _ __   |_ _|__| |
## | |_) | '_ \ / _ \| __/ _ \| '_ \   | |/ _` |
## |  __/| | | | (_) | || (_) | | | |  | | (_| |
## |_|   |_| |_|\___/ \__\___/|_| |_| |___\__,_|
##        
process.photonId = process.photonIsolation.clone()
process.photonId.cut = cms.string(#"userInt('passElectronConvVeto') > 0"
    " hadTowOverEm < userFloat('hadTowOverEmMediumCut') && hadronicOverEm < 0.5"
    " && sigmaIetaIeta > 0.001 && sigmaIetaIeta < userFloat('showerShapeMediumCut')"
    )
process.photonIDiso = process.photonIsolation.clone()
process.photonIDiso.cut = cms.string(#"userInt('passElectronConvVeto') > 0"
    " hadTowOverEm < userFloat('hadTowOverEmMediumCut') && hadronicOverEm < 0.5"
    " && (userFloat('pfChargedPURel') < userFloat('pfChargedMediumCut'))"
    " && (userFloat('pfNeutralPURel') < userFloat('pfNeutralMediumCut'))"
    " && (userFloat('pfGammaPURel')   < userFloat('pfGammaMediumCut'))"
    " && sigmaIetaIeta > 0.001 && sigmaIetaIeta < userFloat('showerShapeMediumCut')"
    )

##    _____     _                         __  __       _       _     _             
##   |_   _| __(_) __ _  __ _  ___ _ __  |  \/  | __ _| |_ ___| |__ (_)_ __   __ _ 
##     | || '__| |/ _` |/ _` |/ _ \ '__| | |\/| |/ _` | __/ __| '_ \| | '_ \ / _` |
##     | || |  | | (_| | (_| |  __/ |    | |  | | (_| | || (__| | | | | | | | (_| |
##     |_||_|  |_|\__, |\__, |\___|_|    |_|  |_|\__,_|\__\___|_| |_|_|_| |_|\__, |
##                |___/ |___/                                                |___/ 
##   
process.probePhotonsPassingHLT = cms.EDProducer(
    "trgMatchedPhotonProducer",                     
    InputProducer = cms.InputTag("probePhotons"),
    hltTags = cms.VInputTag(
    cms.InputTag(HLTPath1,"",InputTagProcess),
    cms.InputTag(HLTPath2,"",InputTagProcess),
    cms.InputTag(HLTPath3,"",InputTagProcess),
    ),
    triggerEventTag = cms.untracked.InputTag("hltTriggerSummaryAOD","",InputTagProcess),
    triggerResultsTag = cms.untracked.InputTag("TriggerResults", "", InputTagProcess)
    )


##    _____      _                        _  __     __             
##   | ____|_  _| |_ ___ _ __ _ __   __ _| | \ \   / /_ _ _ __ ___ 
##   |  _| \ \/ / __/ _ \ '__| '_ \ / _` | |  \ \ / / _` | '__/ __|
##   | |___ >  <| ||  __/ |  | | | | (_| | |   \ V / (_| | |  \__ \
##   |_____/_/\_\\__\___|_|  |_| |_|\__,_|_|    \_/ \__,_|_|  |___/
##   

## Here we show how to use a module to compute an external variable
#producer of dR < 0.5 photon-cleaned jets
process.cleanJets = cms.EDProducer("JetViewCleaner",
    srcObject = cms.InputTag(JET_COLL, "", "RECO"),
    srcObjectSelection = cms.string(JET_CUTS),
    srcObjectsToRemove = cms.VInputTag( cms.InputTag("photons", "", RECOProcess)),
    deltaRMin = cms.double(0.5)  
    )


#produce dR(photon, nearest IDed uncorrected jet passing cuts on corrected eta and pT)
process.photonDRToNearestJet = cms.EDProducer("DeltaRNearestJetComputer",
    probes = cms.InputTag("probePhotons"),
       # ^^--- NOTA BENE: if probes are defined by ref, as in this case, 
       #       this must be the full collection, not the subset by refs.
    objects = cms.InputTag("cleanJets"),
    objectSelection = cms.string(JET_CUTS)
)


#count jets passing cuts
process.JetMultiplicity = cms.EDProducer("CandMultiplicityCounter",
    probes = cms.InputTag("probePhotons"),
    objects = cms.InputTag("cleanJets"),
    objectSelection = cms.string(JET_CUTS),
    )


process.ext_ToNearestJet_sequence = cms.Sequence(
    process.cleanJets + 
    process.photonDRToNearestJet +
    process.JetMultiplicity
    )


##    _____             ____        __ _       _ _   _             
##   |_   _|_ _  __ _  |  _ \  ___ / _(_)_ __ (_) |_(_) ___  _ __  
##     | |/ _` |/ _` | | | | |/ _ \ |_| | '_ \| | __| |/ _ \| '_ \ 
##     | | (_| | (_| | | |_| |  __/  _| | | | | | |_| | (_) | | | |
##     |_|\__,_|\__, | |____/ \___|_| |_|_| |_|_|\__|_|\___/|_| |_|
##              |___/                                              

## tag should be a well reconstructed electron. We use VBTF WP80.
#process.ElectronPassingWP80 = cms.EDFilter("GsfElectronRefSelector",
process.ElectronPassingWP80 = cms.EDFilter("PATElectronRefSelector",
    src = cms.InputTag( ELECTRON_COLL ),
    cut = cms.string( "("+ELECTRON_CUTS+")"+" && ("+ELECTRON_MED_EB+" || "+ELECTRON_MED_EE+")"
    ) 
)

process.Tag = process.ElectronPassingWP80.clone()
process.photon_sequence = cms.Sequence(
    process.probePhotons +
    process.photonId +
    process.photonIsolation +
    process.photonIDiso +
#    process.probePhotonsPassingHLT + 
    process.ElectronPassingWP80 +
    process.Tag
    )


##    _____ ___   ____    ____       _          
##   |_   _( _ ) |  _ \  |  _ \ __ _(_)_ __ ___ 
##     | | / _ \/\ |_) | | |_) / _` | | '__/ __|
##     | || (_>  <  __/  |  __/ (_| | | |  \__ \
##     |_| \___/\/_|     |_|   \__,_|_|_|  |___/
##                                              
##   
#  Tag & probe selection ######
process.tagPhoton = cms.EDProducer("CandViewShallowCloneCombiner",
                                   decay = cms.string("Tag probePhotons"),
                                   checkCharge = cms.bool(False),
                                   cut = cms.string("60 < mass < 120")
                                   )
process.tagphotonIDiso = process.tagPhoton.clone()
process.tagphotonIDiso.decay = cms.string("Tag photonIDiso")

process.allTagsAndProbes = cms.Sequence(
    process.tagPhoton +
    process.tagphotonIDiso
)


##    __  __  ____   __  __       _       _               
##   |  \/  |/ ___| |  \/  | __ _| |_ ___| |__   ___  ___ 
##   | |\/| | |     | |\/| |/ _` | __/ __| '_ \ / _ \/ __|
##   | |  | | |___  | |  | | (_| | || (__| | | |  __/\__ \
##   |_|  |_|\____| |_|  |_|\__,_|\__\___|_| |_|\___||___/
##                                                        

#process.McMatchTag = cms.EDFilter("MCTruthDeltaRMatcherNew",
process.McMatchTag = cms.EDProducer("MCTruthDeltaRMatcherNew",
    matchPDGId = cms.vint32(11),
    src = cms.InputTag("Tag"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles"),
    checkCharge = cms.bool(True)
)
process.McMatchPhoton = process.McMatchTag.clone()
process.McMatchPhoton.src = cms.InputTag("probePhotons")
process.McMatchIdIso = process.McMatchTag.clone()
process.McMatchIdIso.src = cms.InputTag("photonIDiso")
process.McMatchHLT = process.McMatchTag.clone()
process.McMatchHLT.src = cms.InputTag("probePhotonsPassingHLT")

process.mc_sequence = cms.Sequence(
    process.McMatchTag
   +process.McMatchPhoton
   +process.McMatchIdIso
  # +process.McMatchHLT
)

############################################################################
##    _____           _       _ ____            _            _   _  ____  ##
##   |_   _|_ _  __ _( )_ __ ( )  _ \ _ __ ___ | |__   ___  | \ | |/ ___| ##
##     | |/ _` |/ _` |/| '_ \|/| |_) | '__/ _ \| '_ \ / _ \ |  \| | |  _  ##
##     | | (_| | (_| | | | | | |  __/| | | (_) | |_) |  __/ | |\  | |_| | ##
##     |_|\__,_|\__, | |_| |_| |_|   |_|  \___/|_.__/ \___| |_| \_|\____| ##
##              |___/                                                     ##
##                                                                        ##
############################################################################
##    ____                      _     _           
##   |  _ \ ___ _   _ ___  __ _| |__ | | ___  ___ 
##   | |_) / _ \ | | / __|/ _` | '_ \| |/ _ \/ __|
##   |  _ <  __/ |_| \__ \ (_| | |_) | |  __/\__ \
##   |_| \_\___|\__,_|___/\__,_|_.__/|_|\___||___/


## I define some common variables for re-use later.
## This will save us repeating the same code for each efficiency category

ZVariablesToStore = cms.PSet(
    eta = cms.string("eta"),
    pt  = cms.string("pt"),
    phi  = cms.string("phi"),
    et  = cms.string("et"),
    e  = cms.string("energy"),
    p  = cms.string("p"),
    px  = cms.string("px"),
    py  = cms.string("py"),
    pz  = cms.string("pz"),
    theta  = cms.string("theta"),    
    vx     = cms.string("vx"),
    vy     = cms.string("vy"),
    vz     = cms.string("vz"),
    rapidity  = cms.string("rapidity"),
    mass  = cms.string("mass"),
    mt  = cms.string("mt"),    
)   


TagPhotonVariablesToStore = cms.PSet(
    eta = cms.string("eta"),
    pt  = cms.string("pt"),
    phi  = cms.string("phi"),
    px  = cms.string("px"),
    py  = cms.string("py"),
    pz  = cms.string("pz"),
    passConvVeto  = cms.string("userInt('passConvVeto')"),     
    ## super cluster quantities
    sc_energy = cms.string("superCluster.energy"),
    sc_et     = cms.string("superCluster.energy*sin(superCluster.position.theta)"),    
    sc_eta    = cms.string("superCluster.eta"),
    sc_phi    = cms.string("superCluster.phi"),
)


ProbePhotonVariablesToStore = cms.PSet(
        probe_eta = cms.string("eta"),
        probe_phi  = cms.string("phi"),
        probe_et  = cms.string("et"),
        probe_px  = cms.string("px"),
        probe_py  = cms.string("py"),
        probe_pz  = cms.string("pz"),
        ## isolation 
        probe_pfChargedPURel  = cms.string("userFloat('pfChargedPURel')"),
        probe_pfNeutralPURel  = cms.string("userFloat('pfNeutralPURel')"),
        probe_pfGammaPURel    = cms.string("userFloat('pfGammaPURel')"),
        ## booleans
        probe_isPhoton  = cms.string("isPhoton"),     
        probe_passElectronConvVeto  = cms.string("userInt('passElectronConvVeto')"),     

        ## Hcal energy over Ecal Energy
        probe_hadronicOverEm = cms.string("hadTowOverEm"),
        ## Cluster shape information
        probe_sigmaIetaIeta = cms.string("sigmaIetaIeta"),
        ## Pixel seed
        probe_hasPixelSeed = cms.string("hasPixelSeed")
)


CommonStuffForPhotonProbe = cms.PSet(
   variables = cms.PSet(ProbePhotonVariablesToStore),
   ignoreExceptions =  cms.bool (False),
   #fillTagTree      =  cms.bool (True),
   addRunLumiInfo   =  cms.bool (True),
   addEventVariablesInfo   =  cms.bool (True),
   pairVariables =  cms.PSet(ZVariablesToStore),
   pairFlags     =  cms.PSet(
          mass60to120 = cms.string("60 < mass < 120")
    ),
    tagVariables   =  cms.PSet(TagPhotonVariablesToStore),
    tagFlags     =  cms.PSet(
          flag = cms.string("pt>0")
    ),    
)



if MC_FLAG:
    mcTruthCommonStuff = cms.PSet(
        isMC = cms.bool(MC_FLAG),
        tagMatches = cms.InputTag("McMatchTag"),
        motherPdgId = cms.vint32(22,23),
        makeMCUnbiasTree = cms.bool(MC_FLAG),
        checkMotherInUnbiasEff = cms.bool(MC_FLAG),
        mcVariables = cms.PSet(
        probe_eta = cms.string("eta"),
        probe_pt  = cms.string("pt"),
        probe_phi  = cms.string("phi"),
        probe_et  = cms.string("et"),
        probe_e  = cms.string("energy"),
        probe_p  = cms.string("p"),
        probe_px  = cms.string("px"),
        probe_py  = cms.string("py"),
        probe_pz  = cms.string("pz"),
        probe_theta  = cms.string("theta"),    
        probe_vx     = cms.string("vx"),
        probe_vy     = cms.string("vy"),
        probe_vz     = cms.string("vz"),   
        probe_charge = cms.string("charge"),
        probe_rapidity  = cms.string("rapidity"),    
        probe_mass  = cms.string("mass"),
        probe_mt  = cms.string("mt"),    
        ),
        mcFlags     =  cms.PSet(
        probe_flag = cms.string("pt>0")
        ),      
        )
else:
     mcTruthCommonStuff = cms.PSet(
         isMC = cms.bool(False)
         )


##    ___                 ___    _ 
##  |_ _|___  ___       |_ _|__| |
##   | |/ __|/ _ \       | |/ _` |
##   | |\__ \ (_) |  _   | | (_| |
##   |___|___/\___/  ( ) |___\__,_|
##                   |/            
##  Photon --> isolation, id  etc.

process.PhotonToIsoId = cms.EDAnalyzer("TagProbeFitTreeProducer",
    ## pick the defaults
    mcTruthCommonStuff,
    CommonStuffForPhotonProbe,
    # choice of tag and probe pairs, and arbitration                 
    tagProbePairs = cms.InputTag("tagPhoton"),
    arbitration   = cms.string("None"),                      
    flags = cms.PSet(
        probe_passingId  = cms.InputTag("photonId"),
        probe_passingIso = cms.InputTag("photonIsolation"),
#        probe_passingHLT = cms.InputTag("probePhotonsPassingHLT"),
        probe_passingId_iso = cms.InputTag("photonIDiso"),
    ),
    probeMatches  = cms.InputTag("McMatchPhoton"),
    allProbes     = cms.InputTag("probePhotons")
)
process.PhotonToIsoId.variables.probe_dRjet = cms.InputTag("photonDRToNearestJet")
process.PhotonToIsoId.variables.probe_nJets = cms.InputTag("JetMultiplicity")



##    ___    _       __    _   _ _   _____ 
##   |_ _|__| |      \ \  | | | | | |_   _|
##    | |/ _` |  _____\ \ | |_| | |   | |  
##    | | (_| | |_____/ / |  _  | |___| |  
##   |___\__,_|      /_/  |_| |_|_____|_|  

##  offline selection --> HLT. First specify which quantities to store in the TP tree. 
if MC_FLAG:
    HLTmcTruthCommonStuff = cms.PSet(
        isMC = cms.bool(MC_FLAG),
        tagMatches = cms.InputTag("McMatchTag"),
        motherPdgId = cms.vint32(22,23),
        makeMCUnbiasTree = cms.bool(MC_FLAG),
        checkMotherInUnbiasEff = cms.bool(MC_FLAG),
        mcVariables = cms.PSet(
          probe_eta = cms.string("eta"),
          probe_phi  = cms.string("phi"),
          probe_et  = cms.string("et"),
          probe_charge = cms.string("charge"),
        ),
        mcFlags     =  cms.PSet(
          probe_flag = cms.string("pt>0")
        ),      
        )
else:
     HLTmcTruthCommonStuff = cms.PSet(
         isMC = cms.bool(False)
         )


process.photonIDisoToHLT = cms.EDAnalyzer("TagProbeFitTreeProducer",
    HLTmcTruthCommonStuff,                                
    variables = cms.PSet(
      probe_eta = cms.string("eta"),
      probe_phi  = cms.string("phi"),
      probe_et  = cms.string("et"),
    ),
    ignoreExceptions =  cms.bool (False),
    addRunLumiInfo   =  cms.bool (False),
    addEventVariablesInfo   =  cms.bool (False),                                                        
    tagProbePairs = cms.InputTag("tagphotonIDiso"),
    arbitration   = cms.string("None"),
    flags = cms.PSet( 
        probe_passingHLT = cms.InputTag("probePhotonsPassingHLT")        
    ),
    probeMatches  = cms.InputTag("McMatchIdIso"),
    allProbes     = cms.InputTag("photonIDiso")
)

process.tree_sequence = cms.Sequence(
    process.PhotonToIsoId
    #+process.photonIDisoToHLT
)    

##    ____       _   _     
##   |  _ \ __ _| |_| |__  
##   | |_) / _` | __| '_ \ 
##   |  __/ (_| | |_| | | |
##   |_|   \__,_|\__|_| |_|
##

if MC_FLAG:
    process.tagAndProbe = cms.Path(
        process.patObjectSequence +
        process.photon_sequence +
        process.ext_ToNearestJet_sequence + 
        process.allTagsAndProbes +
        process.mc_sequence + 
        process.tree_sequence
        )
else:
    process.tagAndProbe = cms.Path(
        process.patObjectSequence +
        process.photon_sequence +
        process.ext_ToNearestJet_sequence + 
        process.allTagsAndProbes +
        process.tree_sequence
        )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(OUTPUT_FILE_NAME)
                                   )
##file = open('tnp_cfg_mc.py','w')
##file.write(str(process.dumpPython()))
##file.close()
