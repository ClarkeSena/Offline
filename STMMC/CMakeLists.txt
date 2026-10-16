cet_build_plugin(ConcatenateDigitizedWaveforms art::module
    REG_SOURCE src/ConcatenateDigitizedWaveforms_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::RecoDataProducts
)

cet_build_plugin(ConcatenationFilter art::module
    REG_SOURCE src/ConcatenationFilter_module.cc
    LIBRARIES REG
    Offline::RecoDataProducts
)

cet_build_plugin(CountMixedEvents art::module
    REG_SOURCE src/CountMixedEvents_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
)

cet_build_plugin(EmptyStepPointMCs art::module
    REG_SOURCE src/EmptyStepPointMCs_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
)

cet_build_plugin(DetectorEfficiency art::module
    REG_SOURCE src/DetectorEfficiency_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::MCDataProducts
)

cet_build_plugin(HPGeTree art::module
    REG_SOURCE src/HPGeTree_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::MCDataProducts
)

cet_build_plugin(HPGeWaveformsFromStepPointMCs art::module
    REG_SOURCE src/HPGeWaveformsFromStepPointMCs_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::DataProducts
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::ProditionsService
      Offline::RecoDataProducts
      Offline::STMConditions
      Offline::GlobalConstantsService
)

cet_build_plugin(LaBrTree art::module
    REG_SOURCE src/LaBrTree_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::MCDataProducts
)

cet_build_plugin(PHTree art::module
    REG_SOURCE src/PHTree_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::RecoDataProducts
)

cet_build_plugin(ShiftVirtualDetectorStepPointMCs art::module
    REG_SOURCE src/ShiftVirtualDetectorStepPointMCs_module.cc
    LIBRARIES REG
      Offline::GlobalConstantsService
      Offline::MCDataProducts
)

cet_build_plugin(STMResamplingProducer art::module
    REG_SOURCE src/STMResamplingProducer_module.cc
    LIBRARIES REG
        Offline::MCDataProducts
)

cet_build_plugin(STMResamplingFilter art::module
    REG_SOURCE src/STMResamplingFilter_module.cc
    LIBRARIES REG
        Offline::MCDataProducts
)

cet_build_plugin(ValidateConsecutiveEventIDs art::module
    REG_SOURCE src/ValidateConsecutiveEventIDs_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
)

cet_build_plugin(VirtualDetectorTree art::module
    REG_SOURCE src/VirtualDetectorTree_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::GlobalConstantsService
      Offline::MCDataProducts
)

cet_build_plugin(VDResamplerConfigure art::module
    REG_SOURCE src/VDResamplerConfigure_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::GeneralUtilities
      Offline::GlobalConstantsService
      Offline::MCDataProducts
)

cet_build_plugin(VDResamplerReconfigure art::module
    REG_SOURCE src/VDResamplerReconfigure_module.cc
    LIBRARIES REG
      Offline::GeneralUtilities
      Offline::MCDataProducts
)

cet_build_plugin(VDResamplerTrain art::module
    REG_SOURCE src/VDResamplerTrain_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::GeneralUtilities
      Offline::GlobalConstantsService
      Offline::MachineLearningTools
      Offline::MCDataProducts
      Offline::SeedService
)

cet_build_plugin(VDResamplerTrainFromRoot art::module
    REG_SOURCE src/VDResamplerTrainFromRoot_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::GeneralUtilities
      Offline::GlobalConstantsService
      Offline::MachineLearningTools
      Offline::MCDataProducts
      Offline::SeedService
)

cet_build_plugin(VDResamplerGenerateFromModel art::module
    REG_SOURCE src/VDResamplerGenerateFromModel_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::DataProducts
      Offline::GeneralUtilities
      Offline::GlobalConstantsService
      Offline::MachineLearningTools
      Offline::MCDataProducts
      Offline::SeedService
)

cet_build_plugin(VDResamplerGenerateMix art::module
    REG_SOURCE src/VDResamplerGenerateMix_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::DataProducts
      Offline::GeneralUtilities
      Offline::GlobalConstantsService
      Offline::MachineLearningTools
      Offline::MCDataProducts
      Offline::SeedService
)

cet_build_plugin(VDResamplerRegeneratePlots art::module
    REG_SOURCE src/VDResamplerRegeneratePlots_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::GeneralUtilities
      Offline::GlobalConstantsService
      Offline::MachineLearningTools
      Offline::MCDataProducts
)

cet_build_plugin(LaBrTruthSpectrum art::module
    REG_SOURCE src/LaBrTruthSpectrum_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::DataProducts
      Offline::MCDataProducts
)

cet_build_plugin(LaBrTruthPulses art::module
    REG_SOURCE src/LaBrTruthPulses_module.cc
    LIBRARIES REG
      Offline::DataProducts
      Offline::MCDataProducts
)

cet_build_plugin(TagPileupPulses art::module
    REG_SOURCE src/TagPileupPulses_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::MCDataProducts
)

cet_build_plugin(LaBrWaveformsFromStepPointMCsTest art::module
    REG_SOURCE src/LaBrWaveformsFromStepPointMCsTest_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::DataProducts
      Offline::MCDataProducts
      Offline::RecoDataProducts
)

cet_build_plugin(La138Source art::module
    REG_SOURCE src/La138Source_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
      Offline::DataProducts
      Offline::SeedService
)

cet_build_plugin(LaBrIntrinsicSource art::module
    REG_SOURCE src/LaBrIntrinsicSource_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
      Offline::DataProducts
      Offline::SeedService
)

cet_build_plugin(CalibrationSource art::module
    REG_SOURCE src/CalibrationSource_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
      Offline::DataProducts
      Offline::SeedService
)

cet_build_plugin(GenParticleMerger art::module
    REG_SOURCE src/GenParticleMerger_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
)

cet_build_plugin(BinaryWaveformWriter art::module
    REG_SOURCE src/BinaryWaveformWriter_module.cc
    LIBRARIES REG
      Offline::RecoDataProducts
)

cet_build_plugin(PoissonPhotonGun art::module
    REG_SOURCE src/PoissonPhotonGun_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
      Offline::DataProducts
      Offline::SeedService
)

cet_build_plugin(STMGeometryDump art::module
    REG_SOURCE src/STMGeometryDump_module.cc
    LIBRARIES REG
      Offline::GeometryService
      Offline::STMGeom
)

cet_build_plugin(STMParticleSource art::module
    REG_SOURCE src/STMParticleSource_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
      Offline::DataProducts
      Offline::GlobalConstantsService
      Offline::SeedService
      ROOT::Hist
      ROOT::Tree
)

cet_build_plugin(STMParticleDump art::module
    REG_SOURCE src/STMParticleDump_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::GlobalConstantsService
      Offline::MCDataProducts
)

install_source(SUBDIRS src)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/STMMC/fcl)
