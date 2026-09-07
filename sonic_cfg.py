import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

# -------------- Argumentos
options = VarParsing('analysis')
options.register('outputCSV', 'output.csv',         
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 'Archivo CSV de salida')
options.register('isMC', True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 'True para MC, False para datos')
options.setDefault('inputFiles',
    'file:/code/Signal-1000/288B3B53-E2F7-BC4F-843D-E1AD2E1A2335.root')
options.parseArguments()

# ------- Proceso 
process = cms.Process("MYANALYSIS")

# ----- Servicios 
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

# ------- GlobalTag CORRECTO para UL2016 MiniAODv2
# Dataset: RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag,
                               '106X_mcRun2_asymptotic_v17', '')

# ------ Logger ------
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.cerr.threshold = 'INFO'

#----- Eventos ------ 
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# ------ Source ----------------------------
process.source = cms.Source("PoolSource",
    fileNames    = cms.untracked.vstring(options.inputFiles),
    skipBadFiles = cms.untracked.bool(True),
)

# -------- EDAnalyzer ---todos los parametros requeridos por MyAnalyze -----
process.dumper = cms.EDAnalyzer('MyAnalyzer',

    # Colecciones de jets AK8 Puppi (MiniAOD UL2016)
    jets             = cms.InputTag("slimmedJetsAK8"),

    # Vertices primario --- NPV
    primaryVertices  = cms.InputTag("offlineSlimmedPrimaryVertices"),


    genJetsAK8         = cms.InputTag("slimmedGenJetsAK8"),

    genSubJetsAK8      = cms.InputTag("slimmedGenJetsAK8SoftDropSubJets"),

    jecUncFile         = cms.string("Summer19UL16_V7_MC_Uncertainty_AK8PFPuppi.txt"),

    jecUncFileSubjet   = cms.string("Summer19UL16_V7_MC_Uncertainty_AK4PFPuppi.txt"),

    jerResFile         = cms.string("Summer20UL16_JRV3_MC_PtResolution_AK8PFPuppi.txt"),

    jerSFFile          = cms.string("Summer20UL16_JRV3_MC_SF_AK8PFPuppi.txt"),

    jecL1File          = cms.string("Summer19UL16_V7_MC_L1FastJet_AK4PFPuppi.txt"),

    jecL2File          = cms.string("Summer19UL16_V7_MC_L2Relative_AK4PFPuppi.txt"),     

    jecL3File          = cms.string("Summer19UL16_V7_MC_L3Absolute_AK4PFPuppi.txt"), 

    # Pileup summary ---  Pileup_nTrueInt
    pileupInfo       = cms.InputTag("slimmedAddPileupInfo"),

    # Generador --- GenWeight, PtHat
    genEventInfo     = cms.InputTag("generator"),

    # Particulas generadas ---  gen-jet matching, is_signal
    genParticles     = cms.InputTag("prunedGenParticles"),

    # Densidad de energia--- Rho (para JEC L1 y pileup study)
    rho              = cms.InputTag("fixedGridRhoFastjetAll"),

    # True = MC (activa lectura de PU, genInfo, genParticles)
    isMC             = cms.bool(options.isMC),

    # Nombre del CSV de salida
    outputCSV        = cms.string(options.outputCSV),
)

process.p = cms.Path(process.dumper)
