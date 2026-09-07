/**\class MyAnalyzer MyAnalyzer.cc Analisis/MyAnalyzer/plugins/MyAnalyzer.cc

 Description: Analyzer para extraer datos de FatJets y sus constituyentes (Point Cloud).
              Guarda CSV con todas las variables para los estudios de publicacion:
                - Pileup robustness  (NPV, Pileup_nTrueInt, Rho, PF_fromPV, PF_puppiWeight)
                - JEC/JER systematics (JetPt_raw, JEC_unc, JetPt_jec_up/down)
                - JMS/JMR systematics  (subjets raw + gen soft-drop mass)
                - Substructura completa (Tau3, Tau32, N2b1, N2b2, D2)
                - Gen-level / senal (GenJetAK8, DeltaR_reco_gen, Parton_flavor, PtHat)
                - Verdad del W (GenW_*, isMergedW, nQuarksInJet)
                - Constituyentes PF extendidos (fromPV, puppiWeight, charge,
                                                Dxy_err, Dz_err, DeltaR_jet)

 Dataset: QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8
          RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17
          DOI: 10.7483/OPENDATA.CMS.791Q.WIX6

 Referencia JEC:  CMS-JME-18-001
 Referencia JER:  CMS-JME-13-004
 Referencia PU:   CMS-JME-13-005
 Referencia JMS:  CMS-JME-16-003

-----------------------------
 CAMBIOS v3 -> v4  [ACCION 1-5]
-----------------------------

 [ACCION 1] jecUncFile ahora es parametro REQUERIDO (no untracked con default vacio).
            El analyzer lanza excepcion explicita si el .txt no se carga.
            Eliminado el fallback analitico no-oficial.

 [ACCION 2] SDMass_jec_up/down CORREGIDOS: se propagan las JEC del SUBJET,
            no las del jet AK8 padre. Nuevo bloque "subjet JEC unc" con
            JetCorrectionUncertainty separado para AK4PFPuppi subjets.
            Nuevo parametro: jecUncFileSubjet (ruta al .txt de AK4PFPuppi).
            Nuevas columnas CSV:
              SDMass_sjJEC_up, SDMass_sjJEC_down  <- CORRECTAS para publicacion
              SDMass_jec_up, SDMass_jec_down      <- MANTENIDAS como aproximacion
                                                      con flag JEC_unc_isApproxMSD=1

 [ACCION 3] JER smearing hibrido implementado en-CMSSW usando:
              JetResolution    (AK8PFPuppi_pt_resolution.txt)
              JetResolutionSF  (AK8PFPuppi_SF.txt)
            Metodo: scaling si hay gen-match dentro de dR<R/2 y |pT_reco-pT_gen|<3*sigma,
                    smearing estocastico (Gaussian) en caso contrario.
            Nuevas columnas CSV:
              JetPt_jer_nom, JetPt_jer_up, JetPt_jer_down
              SDMass_jer_nom, SDMass_jer_up, SDMass_jer_down
            Nuevos parametros: jerResFile, jerSFFile

 [ACCION 4] genJetsAK8 y genSubJetsAK8 son ahora parametros REQUERIDOS cuando isMC=True.
            Constructor valida esto y lanza excepcion si faltan.

 [ACCION 5] Los ficheros JEC/JER se descargan con get_jec_files.sh (ver script adjunto).
            La ruta por defecto apunta a Summer19UL16_V7_MC_*.txt descargados localmente.

 [COMPAT] Columnas nuevas insertas despues de las existentes; no se rompe el esquema v3
          para las columnas anteriores a la seccion JEC. Las columnas JEC/JER nuevas
          se agregan al bloque [4] y se documenta en el header.

-----------------------------
 CAMBIOS v4 -> v5  [L1+L2+L3 SUBJET]
-----------------------------

 [v5-A] Cadena JEC nominal completa para subjets (L1FastJet + L2Relative + L3Absolute).
        Motivacion: la JetCorrectionUncertainty esta parametrizada como fraccion del
        pT YA CORREGIDO (nominal), no del pT raw. En v4 se evaluaba unc(pT_raw),
        introduciendo un sesgo de 1-3% en subjets con pT < 150 GeV.
        El procedimiento correcto es:
          pT_raw -> L1 -> L2 -> L3 = pT_nominal_sj
          unc = jecUncSubjet(pT_nominal_sj, eta)   <- evaluado en nominal
          pT_sj_up   = pT_nominal_sj * (1 + unc)
          pT_sj_down = pT_nominal_sj * (1 - unc)
          SDMass_sjJEC_up/down = m(p4_sj1_up + p4_sj2_up/down)

        Nuevos parametros requeridos:
          jecL1File    : Summer19UL16_V7_MC_L1FastJet_AK4PFPuppi.txt
          jecL2File    : Summer19UL16_V7_MC_L2Relative_AK4PFPuppi.txt
          jecL3File    : Summer19UL16_V7_MC_L3Absolute_AK4PFPuppi.txt

        Nota: L2L3Residuals NO se incluye (es identidad en MC puro).

 [v5-B] Nuevo objeto FactorizedJetCorrector para aplicar L1+L2+L3 en secuencia.
        Requiere incluir JetMETCorrections/Objects/interface/JetCorrector.h y
        JetMETCorrections/Objects/interface/FactorizedJetCorrector.h.
        Se inicializa en beginJob() y se invoca por subjet en analyze().

 [v5-C] Nueva columna CSV: SubJetN_ptNom (pT nominal = raw * L1*L2*L3).
        Util para validar la cadena de correccion a posteriori.
*/

// --------------- Includes estandar  ------------------------
#include <memory>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <random>

// -----------  Framework CMSSW ---------------------------
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// ----------- Objetos de datos CMS --------------------
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

// ------ Gen jets -----------------------------------------
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

// ------ JEC oficial  ------------------------------------
// Requiere en BuildFile.xml:  <use name="CondFormats/JetMETObjects"/>
//                             <use name="JetMETCorrections/Objects"/>
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
// [v5-A] Cadena nominal L1+L2+L3 para obtener pT_nominal del subjet
//#include "JetMETCorrections/Objects/interface/JetCorrector.h"
//#include "JetMETCorrections/Objects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

// ------ [ACCION 3] JER oficial  -------------------------
// Requiere en BuildFile.xml:  <use name="JetMETCorrections/Modules"/>
//                             <use name="CondFormats/JetMETObjects"/>
#include "JetMETCorrections/Modules/interface/JetResolution.h"

//---- Utilidades de geometria -------------------------
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/LorentzVector.h"


// =============================================================================
// Clase MyAnalyzer
// =============================================================================
class MyAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
    explicit MyAnalyzer(const edm::ParameterSet&);
    ~MyAnalyzer() override;

private:
    void beginJob() override;
    void analyze(const edm::Event&, const edm::EventSetup&) override;
    void endJob() override;

    // ----- Tokens -------------------------------------------------------
    edm::EDGetTokenT<std::vector<pat::Jet>>             jetToken_;
    edm::EDGetTokenT<edm::TriggerResults>               trigToken_;
    edm::EDGetTokenT<reco::VertexCollection>            pvToken_;
    edm::EDGetTokenT<std::vector<PileupSummaryInfo>>    puToken_;
    edm::EDGetTokenT<GenEventInfoProduct>               genInfoToken_;
    edm::EDGetTokenT<reco::GenParticleCollection>       genPartsToken_;
    edm::EDGetTokenT<double>                            rhoToken_;
    edm::EDGetTokenT<reco::GenJetCollection>            genJetsAK8Token_;     // [ACCION 4]
    edm::EDGetTokenT<reco::GenJetCollection>            genSubJetsAK8Token_;  // [ACCION 4]

    // ---- Parametros configurables -------------------------------------
    bool          isMC_;
    std::string   csvFileName_;
    std::string   subjetLabel_;       // "SoftDropPuppi" en UL106X
    std::string   jecUncFile_;        // [ACCION 1] REQUERIDO: AK8PFPuppi Uncertainty .txt
    std::string   jecUncFileSubjet_;  // [ACCION 2] REQUERIDO: AK4PFPuppi Uncertainty .txt
    std::string   jerResFile_;        // [ACCION 3] REQUERIDO: AK8PFPuppi pt_resolution .txt
    std::string   jerSFFile_;         // [ACCION 3] REQUERIDO: AK8PFPuppi SF .txt
    // [v5-A] Cadena JEC nominal L1+L2+L3 para subjets AK4PFPuppi
    std::string   jecL1File_;         // REQUERIDO: L1FastJet_AK4PFPuppi.txt
    std::string   jecL2File_;         // REQUERIDO: L2Relative_AK4PFPuppi.txt
    std::string   jecL3File_;         // REQUERIDO: L3Absolute_AK4PFPuppi.txt
    bool          dumpLabels_;
    std::ofstream csvFile_;

    // ------ Estado interno -------------------------------------------
    std::unique_ptr<JetCorrectionUncertainty>  jecUnc_;       // AK8 JEC unc
    std::unique_ptr<JetCorrectionUncertainty>  jecUncSubjet_; // [ACCION 2] AK4 subjet JEC unc
    // [v5-B] Cadena nominal L1+L2+L3 para subjet: pT_raw -> pT_nominal
    std::unique_ptr<FactorizedJetCorrector>    subjetJEC_;     // L1+L2+L3 AK4PFPuppi
    // [ACCION 3] JER: se leen desde fichero en beginJob()
    std::unique_ptr<JME::JetResolution>         jerResolution_;
    std::unique_ptr<JME::JetResolutionScaleFactor> jerSF_;
    std::mt19937  rng_;              // [ACCION 3] Generador para smearing estocastico

    size_t nHeaderFields_    = 0;
    bool   firstRowChecked_  = false;
    bool   labelsDumped_     = false;
    bool   warnedNoSubjets_  = false;

    // ------ Helpers internos -----------------------------------------
    std::string vecF(const std::vector<float>& v) const {
        if (v.empty()) return "";
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(6);
        oss << v[0];
        for (size_t k = 1; k < v.size(); ++k) oss << ";" << v[k];
        return oss.str();
    }

    std::string vecI(const std::vector<int>& v) const {
        if (v.empty()) return "";
        std::ostringstream oss;
        oss << v[0];
        for (size_t k = 1; k < v.size(); ++k) oss << ";" << v[k];
        return oss.str();
    }

    int countGoodPV(const reco::VertexCollection& pvs) const {
        int n = 0;
        for (const auto& v : pvs)
            if (!v.isFake() && v.isValid()
                && v.ndof()           > 4.0
                && std::abs(v.z())    < 24.0
                && v.position().rho() < 2.0)
                ++n;
        return n;
    }

    float getTrueNPU(const std::vector<PileupSummaryInfo>& puInfo) const {
        for (const auto& pu : puInfo)
            if (pu.getBunchCrossing() == 0)
                return static_cast<float>(pu.getTrueNumInteractions());
        return -1.0f;
    }

    // JEC unc oficial para jet AK8 padre
    float jecUncertainty(float pt, float eta) {
        jecUnc_->setJetPt(pt);
        jecUnc_->setJetEta(eta);
        return static_cast<float>(jecUnc_->getUncertainty(true));
    }

    // [ACCION 2] JEC unc oficial para subjets (AK4PFPuppi).
    // IMPORTANTE: pt debe ser el pT NOMINAL (ya corregido L1+L2+L3),
    // no el pT raw. La parametrizacion de JetCorrectionUncertainty
    // esta definida sobre el jet ya corregido (CMS-JME-18-001).
    float jecUncertaintySubjet(float pt_nominal, float eta) {
        jecUncSubjet_->setJetPt(pt_nominal);
        jecUncSubjet_->setJetEta(eta);
        return static_cast<float>(jecUncSubjet_->getUncertainty(true));
    }

    // [v5-B] Aplica la cadena L1+L2+L3 al subjet y devuelve pT_nominal.
    // El FactorizedJetCorrector requiere setJetPt(raw), setJetEta, setJetA, setRho
    // antes de llamar a getCorrection(). El factor multiplicativo devuelto se
    // aplica al p4 raw completo (pT y masa escalan igual a primer orden).
    float subjetNominalFactor(float pt_raw, float eta,
                               float area,  float rho) {
        if (!subjetJEC_ || pt_raw <= 0.0f) return 1.0f;
        subjetJEC_->setJetPt(pt_raw);
        subjetJEC_->setJetEta(eta);
        subjetJEC_->setJetA(area);
        subjetJEC_->setRho(rho);
        float factor = static_cast<float>(subjetJEC_->getCorrection());
        // Proteccion: el factor nunca debe ser negativo ni enormemente grande
        return std::max(0.01f, std::min(factor, 5.0f));
    }

    // [ACCION 3] JER smearing hibrido (CMS-JME-13-004)
    // Devuelve {factor_nom, factor_up, factor_down}
    // factor se multiplica directamente al pT o masa del jet reco.
    struct JERFactors { float nom, up, dn; };
    JERFactors computeJERfactors(float pt_reco, float eta, float rho,
                                  float pt_gen, bool hasGenMatch) const
    {
        // Parametros de resolucion del jet reco
        JME::JetParameters jparams;
        jparams.setJetPt(pt_reco);
        jparams.setJetEta(eta);
        jparams.setRho(rho);

        float sigma = jerResolution_->getResolution(jparams);  // sigma/pT relativa
        float sf_nom = jerSF_->getScaleFactor(jparams, Variation::NOMINAL);
        float sf_up  = jerSF_->getScaleFactor(jparams, Variation::UP);
        float sf_dn  = jerSF_->getScaleFactor(jparams, Variation::DOWN);

        auto smear = [&](float sf) -> float {
            if (sf <= 1.0f) return 1.0f; // sin smearing si SF < 1
            if (hasGenMatch) {
                // Metodo scaling (determinista)
                float c = 1.0f + (sf - 1.0f) * (pt_reco - pt_gen) / pt_reco;
                return std::max(0.0f, c);
            } else {
                // Metodo estocastico: Gaussian extra con sigma = sqrt(sf^2 - 1)*sigma_reco
                float sigma_extra = std::sqrt(sf * sf - 1.0f) * sigma;
                std::normal_distribution<float> gauss(0.0f, sigma_extra);
                float delta = gauss(const_cast<std::mt19937&>(rng_));
                return std::max(0.0f, 1.0f + delta);
            }
        };

        return { smear(sf_nom), smear(sf_up), smear(sf_dn) };
    }

    float getUF(const pat::Jet& jet, const std::string& key) const {
        return jet.hasUserFloat(key) ? jet.userFloat(key) : -999.0f;
    }

    float getBD(const pat::Jet& j, const std::string& key) const {
        float v = j.bDiscriminator(key);
        return (v < -900.0f) ? -999.0f : v;
    }

    static size_t countFields(const std::string& s) {
        return 1 + static_cast<size_t>(std::count(s.begin(), s.end(), ','));
    }
};


// =============================================================================
// Constructor
// =============================================================================
MyAnalyzer::MyAnalyzer(const edm::ParameterSet& iConfig)
    : jetToken_     (consumes<std::vector<pat::Jet>>(
                         iConfig.getParameter<edm::InputTag>("jets")))
    , trigToken_    (consumes<edm::TriggerResults>(
                         edm::InputTag("TriggerResults", "", "HLT")))
    , pvToken_      (consumes<reco::VertexCollection>(
                         iConfig.getParameter<edm::InputTag>("primaryVertices")))
    , puToken_      (consumes<std::vector<PileupSummaryInfo>>(
                         iConfig.getParameter<edm::InputTag>("pileupInfo")))
    , genInfoToken_ (consumes<GenEventInfoProduct>(
                         iConfig.getParameter<edm::InputTag>("genEventInfo")))
    , genPartsToken_(consumes<reco::GenParticleCollection>(
                         iConfig.getParameter<edm::InputTag>("genParticles")))
    , rhoToken_     (consumes<double>(
                         iConfig.getParameter<edm::InputTag>("rho")))
    // [ACCION 4] genJetsAK8 y genSubJetsAK8 ahora son REQUERIDOS cuando isMC=True.
    // Se usan mayConsume para que compilen tambien en modo datos (isMC=False).
    , genJetsAK8Token_   (mayConsume<reco::GenJetCollection>(
                         iConfig.getParameter<edm::InputTag>("genJetsAK8")))
    , genSubJetsAK8Token_(mayConsume<reco::GenJetCollection>(
                         iConfig.getParameter<edm::InputTag>("genSubJetsAK8")))
    , isMC_         (iConfig.getParameter<bool>("isMC"))
    , csvFileName_  (iConfig.getParameter<std::string>("outputCSV"))
    , subjetLabel_  (iConfig.getUntrackedParameter<std::string>("subjetLabel",
                                                                "SoftDropPuppi"))
    // [ACCION 1] Parametros JEC/JER requeridos
    , jecUncFile_       (iConfig.getParameter<std::string>("jecUncFile"))
    , jecUncFileSubjet_ (iConfig.getParameter<std::string>("jecUncFileSubjet"))
    , jerResFile_       (iConfig.getParameter<std::string>("jerResFile"))
    , jerSFFile_        (iConfig.getParameter<std::string>("jerSFFile"))
    // [v5-A] Cadena nominal L1+L2+L3 para subjets
    , jecL1File_        (iConfig.getParameter<std::string>("jecL1File"))
    , jecL2File_        (iConfig.getParameter<std::string>("jecL2File"))
    , jecL3File_        (iConfig.getParameter<std::string>("jecL3File"))
    , dumpLabels_   (iConfig.getUntrackedParameter<bool>("dumpLabels", true))
    // [ACCION 3] Semilla fija para reproducibilidad del smearing estocastico
    , rng_          (12345U)
{
    // [ACCION 4] Validacion temprana: en modo MC los gen-jets son obligatorios
    if (isMC_) {
        edm::InputTag gjTag  = iConfig.getParameter<edm::InputTag>("genJetsAK8");
        edm::InputTag gsjTag = iConfig.getParameter<edm::InputTag>("genSubJetsAK8");
        if (gjTag.label().empty() || gsjTag.label().empty())
            throw cms::Exception("MyAnalyzer")
                << "isMC=True requiere genJetsAK8 y genSubJetsAK8 en la configuracion.";
    }
}

MyAnalyzer::~MyAnalyzer() {}


// =============================================================================
// beginJob --- abrir CSV, cargar JEC/JER, escribir cabecera
// =============================================================================
void MyAnalyzer::beginJob()
{

    // --- DEBUG: imprimir rutas antes de cualquier carga ---
    fprintf(stderr, "[DEBUG] jecUncFile_       = '%s'\n", jecUncFile_.c_str());
    fprintf(stderr, "[DEBUG] jecUncFileSubjet_ = '%s'\n", jecUncFileSubjet_.c_str());
    fprintf(stderr, "[DEBUG] jecL1File_        = '%s'\n", jecL1File_.c_str());
    fprintf(stderr, "[DEBUG] jecL2File_        = '%s'\n", jecL2File_.c_str());
    fprintf(stderr, "[DEBUG] jecL3File_        = '%s'\n", jecL3File_.c_str());
    fprintf(stderr, "[DEBUG] jerResFile_       = '%s'\n", jerResFile_.c_str());
    fprintf(stderr, "[DEBUG] jerSFFile_        = '%s'\n", jerSFFile_.c_str());
    fflush(stderr);

    csvFile_.open(csvFileName_);
    if (!csvFile_.is_open())
        throw std::runtime_error("No se pudo abrir el CSV: " + csvFileName_);
    csvFile_ << std::fixed << std::setprecision(6);

   // ------------------------------------------------------------------
    // Helper lambda: verifica que un fichero .txt existe y no esta vacio
    // ANTES de pasarlo a JetCorrectorParameters, que llama abort() en
    // lugar de lanzar excepcion cuando el fichero no tiene definiciones.
    // ------------------------------------------------------------------
    auto checkTxt = [](const std::string& path, const std::string& label) {
        if (path.empty())
            throw cms::Exception("MyAnalyzer")
                << label << ": ruta vacia. Verifica el cfg.";
        std::ifstream f(path);
        if (!f.good())
            throw cms::Exception("MyAnalyzer")
                << label << ": fichero no encontrado:\n  " << path
                << "\n  Descarga con: bash get_jec_files.sh";
        // Leer primera linea no vacia: debe empezar con '{' (formato JEC)
        std::string line;
        bool found = false;
        while (std::getline(f, line)) {
            if (!line.empty() && line[0] != '#') { found = true; break; }
        }
        if (!found)
            throw cms::Exception("MyAnalyzer")
                << label << ": fichero vacio o sin datos:\n  " << path;
    };

    // ------------------------------------------------------------------
    // [ACCION 1] JEC uncertainty AK8 padre
    // Fichero: Summer19UL16_V7_MC_Uncertainty_AK8PFPuppi.txt
    // ------------------------------------------------------------------
    checkTxt(jecUncFile_, "jecUncFile");
    try {
	fprintf(stderr, "[DEBUG] Cargando AK8 uncertainty...\n"); fflush(stderr);
	JetCorrectorParameters parAK8(jecUncFile_);
	jecUnc_ = std::make_unique<JetCorrectionUncertainty>(parAK8);
        edm::LogWarning("MyAnalyzer") << "[JEC-AK8] OK: " << jecUncFile_;
    } catch (const std::exception& e) {
        throw cms::Exception("MyAnalyzer")
            << "[ACCION 1] Error al cargar jecUncFile con tag 'Total'.\n"
            << "  Fichero: " << jecUncFile_ << "\n"
            << "  Verifica que el fichero contenga el tag {... Total ...}\n"
            << "  Error: " << e.what();
    }

    // ------------------------------------------------------------------
    // [ACCION 2] JEC uncertainty subjets AK4PFPuppi
    // Fichero: Summer19UL16_V7_MC_Uncertainty_AK4PFPuppi.txt
    // ------------------------------------------------------------------
    checkTxt(jecUncFileSubjet_, "jecUncFileSubjet");
    try {
        fprintf(stderr, "[DEBUG] AK8 OK. Cargando AK4 subjet uncertainty...\n"); fflush(stderr);
	JetCorrectorParameters parAK4(jecUncFileSubjet_);
	jecUncSubjet_ = std::make_unique<JetCorrectionUncertainty>(parAK4);
        edm::LogWarning("MyAnalyzer") << "[JEC-Subjet] OK: " << jecUncFileSubjet_;
    } catch (const std::exception& e) {
        throw cms::Exception("MyAnalyzer")
            << "[ACCION 2] Error al cargar jecUncFileSubjet con tag 'Total'.\n"
            << "  Fichero: " << jecUncFileSubjet_ << "\n"
            << "  Error: " << e.what();
    }

    // ------------------------------------------------------------------
    // [v5] Cadena nominal L1+L2+L3 para subjets AK4PFPuppi
    // Los ficheros L1/L2/L3 NO tienen tags internos: se cargan sin tag.
    // Solo los ficheros *_Uncertainty_*.txt usan el tag "Total".
    // ------------------------------------------------------------------
    checkTxt(jecL1File_, "jecL1File");
    checkTxt(jecL2File_, "jecL2File");
    checkTxt(jecL3File_, "jecL3File");
    try {
        JetCorrectorParameters l1(jecL1File_);
        JetCorrectorParameters l2(jecL2File_);
        JetCorrectorParameters l3(jecL3File_);
        std::vector<JetCorrectorParameters> vParSJ;
        vParSJ.push_back(l1);
        vParSJ.push_back(l2);
        vParSJ.push_back(l3);
        subjetJEC_ = std::make_unique<FactorizedJetCorrector>(vParSJ);
        edm::LogWarning("MyAnalyzer")
            << "[v5] Cadena subjet JEC OK:"
            << " L1=" << jecL1File_
            << " L2=" << jecL2File_
            << " L3=" << jecL3File_;
    } catch (const std::exception& e) {
        throw cms::Exception("MyAnalyzer")
            << "[v5] Error al construir FactorizedJetCorrector L1+L2+L3.\n"
            << "  L1: " << jecL1File_ << "\n"
            << "  L2: " << jecL2File_ << "\n"
            << "  L3: " << jecL3File_ << "\n"
            << "  Error: " << e.what();
    }

    // ----------------------------------------------
    // ------------------------------------------------------------------
    // [ACCION 3] JER resolution y Scale Factor --- REQUERIDOS para MC
    // Ficheros: AK8PFPuppi_pt_resolution.txt
    //           AK8PFPuppi_SF.txt
    // Metodo hibrido (CMS-JME-13-004):
    //   - Scaling si hay gen-match dentro de dR<0.4 y |dPt|<3*sigma*pT
    //   - Smearing estocastico en caso contrario
    // ------------------------------------------------------------------
    if (isMC_) {
        try {
            jerResolution_ = std::make_unique<JME::JetResolution>(
                JME::JetResolution(jerResFile_));
            jerSF_ = std::make_unique<JME::JetResolutionScaleFactor>(
                JME::JetResolutionScaleFactor(jerSFFile_));
            edm::LogInfo("MyAnalyzer")
                << "[JER] Resolution: " << jerResFile_
                << "\n[JER] SF:         " << jerSFFile_;
        } catch (const std::exception& e) {
            throw cms::Exception("MyAnalyzer")
                << "[ACCION 3] No se pudieron cargar ficheros JER: "
                << jerResFile_ << " / " << jerSFFile_
                << "\n  Error: " << e.what()
                << "\n  Descarga con: bash get_jec_files.sh";
        }
    }

    // ------------------------------------------------------------------
    // Cabecera CSV
    // Bloque [4] expandido con columnas JER y mSD-subjet-JEC correctas.
    // ------------------------------------------------------------------
    const std::string header =
        // [1] Nivel evento (7)
        "EventID,RunNumber,NPV,Pileup_nTrueInt,Rho,GenWeight,PtHat,"
        // [2] Jet cinematica (8)
        "JetPt,JetPt_raw,JetArea,JetEta,JetPhi,JetMass,JetMass_raw,SDMass,"
        // [3] Substructura (9)
        "Tau1,Tau2,Tau3,Tau21,Tau32,N2b1,N2b2,D2,NSubjets,"
        // [4] JEC uncertainty AK8 (4) --- pT systematics
        "JEC_unc,JetPt_jec_up,JetPt_jec_down,"
        // [4b] SDMass con JEC de SUBJET [ACCION 2] --- CORRECTAS para publicacion
        "SDMass_sjJEC_up,SDMass_sjJEC_down,"
        // [4c] SDMass con JEC de AK8 --- APROXIMACION (solo referencia, no usar en publicacion)
        "SDMass_jec_approx_up,SDMass_jec_approx_down,"
        // [4d] JER systematics [ACCION 3]
        "JetPt_jer_nom,JetPt_jer_up,JetPt_jer_down,"
        "SDMass_jer_nom,SDMass_jer_up,SDMass_jer_down,"
        "JER_hasGenMatch,"
        // [5] Subjets reco crudos (17) --- +SubJet1_ptNom, SubJet2_ptNom [v5-C]
        "SDMass_raw,"
        "SubJet1_ptRaw,SubJet1_ptNom,SubJet1_eta,SubJet1_phi,SubJet1_massRaw,SubJet1_area,SubJet1_deepCSV,"
        "SubJet2_ptRaw,SubJet2_ptNom,SubJet2_eta,SubJet2_phi,SubJet2_massRaw,SubJet2_area,SubJet2_deepCSV,"
        "SubJet_dR12,SubJet_ptAsym,"
        // [6] Gen jet AK8 + gen soft-drop (15)
        "GenJetAK8_pt,GenJetAK8_eta,GenJetAK8_phi,GenJetAK8_mass,"
        "GenJetAK8_SDMass,GenJetAK8_nSubJets,DeltaR_recoGenJet,"
        "GenSubJet1_pt,GenSubJet1_eta,GenSubJet1_phi,GenSubJet1_mass,"
        "GenSubJet2_pt,GenSubJet2_eta,GenSubJet2_phi,GenSubJet2_mass,"
        // [7] Gen particle / flavor (5)
        "GenPart_pt,DeltaR_reco_gen,Parton_flavor,Had_flavor,is_signal,"
        // [8] Verdad del W (7)
        "GenW_pt,GenW_eta,GenW_phi,GenW_mass,DeltaR_jet_W,nQuarksInJet,isMergedW,"
        // [9] Constituyentes PF (14)
        "N_constituents,"
        "PF_Pt,PF_Eta,PF_Phi,PF_E,PF_PdgId,"
        "PF_Dxy,PF_Dz,PF_Dxy_err,PF_Dz_err,"
        "PF_fromPV,PF_puppiWeight,PF_charge,PF_DeltaR_jet";

    nHeaderFields_ = countFields(header);
    csvFile_ << header << "\n";
    edm::LogInfo("MyAnalyzer") << "Cabecera CSV con " << nHeaderFields_ << " columnas.";
}


// =============================================================================
// endJob
// =============================================================================
void MyAnalyzer::endJob()
{
    if (csvFile_.is_open()) csvFile_.close();
}


// =============================================================================
// analyze --- bucle principal
// =============================================================================
void MyAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;

    // ---- Trigger --------------------------------------------------------
    edm::Handle<edm::TriggerResults> h_trig;
    iEvent.getByToken(trigToken_, h_trig);
    if (!h_trig.isValid()) {
        edm::LogWarning("MyAnalyzer") << "TriggerResults no valido --- saltando evento";
        return;
    }
    const edm::TriggerNames& names = iEvent.triggerNames(*h_trig);
    bool passHLT = false;
    for (unsigned int i = 0; i < h_trig->size(); ++i) {
        if (names.triggerName(i).find("HLT_AK8PFJet450") != std::string::npos
            && h_trig->accept(i)) {
            passHLT = true;
            break;
        }
    }
    if (!passHLT) return;

    // ---- Informacion de nivel evento ------------------------------------
    unsigned long long eventID = iEvent.id().event();
    unsigned int       runNum  = iEvent.id().run();

    edm::Handle<reco::VertexCollection> h_pvs;
    iEvent.getByToken(pvToken_, h_pvs);
    int npv = (h_pvs.isValid()) ? countGoodPV(*h_pvs) : -1;

    edm::Handle<double> h_rho;
    iEvent.getByToken(rhoToken_, h_rho);
    float rho = (h_rho.isValid()) ? static_cast<float>(*h_rho) : -999.0f;

    float pileup_nTrueInt = -1.0f;
    float genWeight       =  1.0f;
    float ptHat           = -999.0f;

    if (isMC_) {
        edm::Handle<std::vector<PileupSummaryInfo>> h_pu;
        iEvent.getByToken(puToken_, h_pu);
        if (h_pu.isValid()) pileup_nTrueInt = getTrueNPU(*h_pu);

        edm::Handle<GenEventInfoProduct> h_geninfo;
        iEvent.getByToken(genInfoToken_, h_geninfo);
        if (h_geninfo.isValid()) {
            genWeight = static_cast<float>(h_geninfo->weight());
            ptHat     = static_cast<float>(h_geninfo->qScale());
        }
    }

    edm::Handle<reco::GenParticleCollection> h_genparts;
    if (isMC_) iEvent.getByToken(genPartsToken_, h_genparts);

    // ---- [ACCION 4] Gen jets AK8 --- REQUERIDOS en MC -------------------
    edm::Handle<reco::GenJetCollection> h_genjetsAK8;
    edm::Handle<reco::GenJetCollection> h_gensubjetsAK8;
    if (isMC_) {
        iEvent.getByToken(genJetsAK8Token_,    h_genjetsAK8);
        iEvent.getByToken(genSubJetsAK8Token_, h_gensubjetsAK8);
        if (!h_genjetsAK8.isValid())
            edm::LogError("MyAnalyzer")
                << "[ACCION 4] slimmedGenJetsAK8 no valido. "
                << "Verifica que 'genJetsAK8' este en el cfg.";
        if (!h_gensubjetsAK8.isValid())
            edm::LogError("MyAnalyzer")
                << "[ACCION 4] slimmedGenJetsAK8SoftDropSubJets no valido. "
                << "Verifica que 'genSubJetsAK8' este en el cfg.";
    }

    // ---- Jets -----------------------------------------------------------
    edm::Handle<std::vector<pat::Jet>> h_jets;
    iEvent.getByToken(jetToken_, h_jets);
    if (!h_jets.isValid()) {
        edm::LogWarning("MyAnalyzer") << "Jets no validos -- saltando evento";
        return;
    }

    for (const auto& jet : *h_jets) {

        // ---- Cortes cinematicos -----------------------------------------
        if (jet.pt()            < 550.0) continue;
        if (std::abs(jet.eta()) >   2.4) continue;

        // Volcado unico de nombres disponibles
        if (dumpLabels_ && !labelsDumped_) {
            labelsDumped_ = true;
            std::ostringstream dump;
            dump << "\n===== userFloatNames =====\n";
            for (const auto& n : jet.userFloatNames()) dump << "  " << n << "\n";
            dump << "===== subjetCollectionNames =====\n";
            for (const auto& n : jet.subjetCollectionNames()) dump << "  " << n << "\n";
            dump << "===== bDiscriminators (jet AK8) =====\n";
            for (const auto& p : jet.getPairDiscri()) dump << "  " << p.first << "\n";
            edm::LogWarning("MyAnalyzer") << dump.str();
        }

        // ---- Cinematica del jet -----------------------------------------
        float pt       = jet.pt();
        auto  rawP4    = jet.correctedP4(0);
        float pt_raw   = rawP4.pt();
        float mass_raw = rawP4.mass();
        float area     = jet.jetArea();
        float eta      = jet.eta();
        float phi      = jet.phi();
        float mass     = jet.mass();
        float sd_mass  = getUF(jet, "ak8PFJetsPuppiSoftDropMass");

        // ---- Substructura -----------------------------------------------
        float t1 = getUF(jet, "NjettinessAK8Puppi:tau1");
        float t2 = getUF(jet, "NjettinessAK8Puppi:tau2");
        float t3 = getUF(jet, "NjettinessAK8Puppi:tau3");

        float tau21 = (t1 > 0.0f) ? t2 / t1 : -999.0f;
        float tau32 = (t2 > 0.0f) ? t3 / t2 : -999.0f;

        float n2b1 = getUF(jet, "ak8PFJetsPuppiSoftDropValueMap:nb1AK8PuppiSoftDropN2");
        if (n2b1 < -998.0f) n2b1 = getUF(jet, "nb1AK8PuppiSoftDropN2");

        float n2b2 = getUF(jet, "ak8PFJetsPuppiSoftDropValueMap:nb2AK8PuppiSoftDropN2");
        if (n2b2 < -998.0f) n2b2 = getUF(jet, "nb2AK8PuppiSoftDropN2");

        float d2 = getUF(jet, "ak8PFJetsPuppiSoftDropValueMap:nb1AK8PuppiSoftDropD2");
        if (d2 < -998.0f) d2 = getUF(jet, "nb1AK8PuppiSoftDropD2");

        // ---- Subjets soft-drop reco ------------------------------------
        std::vector<const pat::Jet*> sjs;
        {
            auto ptrs = jet.subjets(subjetLabel_);
            if (ptrs.empty() && !warnedNoSubjets_) {
                warnedNoSubjets_ = true;
                edm::LogWarning("MyAnalyzer")
                    << "jet.subjets(\"" << subjetLabel_ << "\") vacio. "
                    << "Revisa subjetCollectionNames en el volcado de arriba.";
            }
            for (const auto& p : ptrs) if (p.isNonnull()) sjs.push_back(p.get());
            std::sort(sjs.begin(), sjs.end(),
                      [](const pat::Jet* a, const pat::Jet* b){ return a->pt() > b->pt(); });
        }
        int nsubjets = static_cast<int>(sjs.size());

        float sj_ptRaw[2]   = {-999.0f, -999.0f};
        float sj_ptNom[2]   = {-999.0f, -999.0f};  // [v5-C] pT nominal = raw*L1*L2*L3
        float sj_eta[2]     = {-999.0f, -999.0f};
        float sj_phi[2]     = {-999.0f, -999.0f};
        float sj_massRaw[2] = {-999.0f, -999.0f};
        float sj_area[2]    = {-999.0f, -999.0f};
        float sj_deepCSV[2] = {-999.0f, -999.0f};

        math::XYZTLorentzVector sdSumRaw(0, 0, 0, 0);
        for (int k = 0; k < std::min(nsubjets, 2); ++k) {
            const pat::Jet* sj = sjs[k];
            math::XYZTLorentzVector raw =
                sj->jecSetsAvailable() ? sj->correctedP4(0) : sj->p4();
            sj_ptRaw[k]   = static_cast<float>(raw.pt());
            sj_eta[k]     = static_cast<float>(raw.eta());
            sj_phi[k]     = static_cast<float>(raw.phi());
            sj_massRaw[k] = static_cast<float>(raw.mass());
            sj_area[k]    = sj->jetArea();
            float pb  = getBD(*sj, "pfDeepCSVJetTags:probb");
            float pbb = getBD(*sj, "pfDeepCSVJetTags:probbb");
            sj_deepCSV[k] = (pb > -900.0f && pbb > -900.0f) ? (pb + pbb) : -999.0f;
            sdSumRaw += raw;
            // [v5-C] pT nominal: aplica L1+L2+L3 sobre el pT raw del subjet
            if (sj_ptRaw[k] > 0.0f && std::abs(sj_eta[k]) < 5.0f) {
                float jec_f = subjetNominalFactor(sj_ptRaw[k], sj_eta[k],
                                                   sj_area[k], rho);
                sj_ptNom[k] = sj_ptRaw[k] * jec_f;
            }
        }
        float sd_mass_raw = (nsubjets >= 1) ? static_cast<float>(sdSumRaw.M()) : -999.0f;

        float sj_dR12   = -999.0f;
        float sj_ptAsym = -999.0f;
        if (nsubjets >= 2) {
            sj_dR12 = reco::deltaR(sj_eta[0], sj_phi[0], sj_eta[1], sj_phi[1]);
            float den = sj_ptRaw[0] + sj_ptRaw[1];
            if (den > 0.0f) sj_ptAsym = (sj_ptRaw[0] - sj_ptRaw[1]) / den;
        }

        // ----------------------------------------------------------------
        // [ACCION 1] JEC uncertainty AK8 para pT del jet padre
        // ----------------------------------------------------------------
        float jec_unc    = jecUncertainty(pt, eta);
        float pt_jec_up  = pt * (1.0f + jec_unc);
        float pt_jec_down= pt * (1.0f - jec_unc);

        // ----------------------------------------------------------------
        // [v5] SDMass con JEC de SUBJET --- procedimiento RIGUROSO
        // Referencia: CMS-JME-16-003, CMS-JME-18-001
        //
        // Procedimiento correcto (v4 evaluaba unc en pT_raw — incorrecto):
        //   pT_raw  ->  L1+L2+L3  ->  pT_nominal_sj
        //   unc = jecUncSubjet(pT_nominal_sj, eta)   <- evaluado en nominal
        //   p4_sj_up   = p4_raw * jec_factor * (1 + unc)
        //   p4_sj_down = p4_raw * jec_factor * (1 - unc)
        //   SDMass_sjJEC_up/down = m(p4_sj1 + p4_sj2) variado
        //
        // Nota: el factor L1+L2+L3 se aplica al p4 raw completo.
        //       pT y masa escalan con el mismo factor (aproximacion
        //       valida; la masa del subjet es pequena respecto a pT).
        // ----------------------------------------------------------------
        float msd_sjJEC_up   = -999.0f;
        float msd_sjJEC_down = -999.0f;

        if (nsubjets >= 1) {
            math::XYZTLorentzVector sdUp(0,0,0,0), sdDn(0,0,0,0);

            for (int k = 0; k < std::min(nsubjets, 2); ++k) {
                const pat::Jet* sj = sjs[k];
                math::XYZTLorentzVector raw =
                    sj->jecSetsAvailable() ? sj->correctedP4(0) : sj->p4();

                float sjPt_raw  = static_cast<float>(raw.pt());
                float sjEta     = static_cast<float>(raw.eta());
                float sjArea    = sj->jetArea();

                if (sjPt_raw <= 0.0f || std::abs(sjEta) >= 5.0f) {
                    // Subjet fuera de rango: contribuye sin variacion
                    sdUp += raw;
                    sdDn += raw;
                    continue;
                }

                // Paso 1: factor nominal L1+L2+L3 (pT_raw -> pT_nominal)
                float jec_factor = subjetNominalFactor(sjPt_raw, sjEta,
                                                        sjArea, rho);
                float sjPt_nom   = sjPt_raw * jec_factor;

                // Paso 2: incertidumbre evaluada en pT_nominal [v5 — CORRECTO]
                float unc_sj = jecUncertaintySubjet(sjPt_nom, sjEta);

                // Paso 3: p4 variados
                // factor_up   = jec_factor * (1 + unc)
                // factor_down = jec_factor * (1 - unc)
                sdUp += raw * (jec_factor * (1.0f + unc_sj));
                sdDn += raw * (jec_factor * (1.0f - unc_sj));
            }

            msd_sjJEC_up   = static_cast<float>(sdUp.M());
            msd_sjJEC_down = static_cast<float>(sdDn.M());
        }

        // SDMass con JEC del AK8 (APROXIMACION --- mantenida como referencia)
        float msd_jec_approx_up   = (sd_mass > 0.0f) ? sd_mass * (1.0f + jec_unc) : -999.0f;
        float msd_jec_approx_down = (sd_mass > 0.0f) ? sd_mass * (1.0f - jec_unc) : -999.0f;

        // ----------------------------------------------------------------
        // [ACCION 3+4] Gen jet AK8 + masa soft-drop generada
        // Necesario para JER matching (Bloque A v3 + ACCION 4)
        // ----------------------------------------------------------------
        float gj_pt = -999.0f, gj_eta = -999.0f, gj_phi = -999.0f, gj_mass = -999.0f;
        float gj_sdmass = -999.0f, gj_dr = -999.0f;
        int   gj_nsub   = 0;
        float gsj_pt[2] = {-999.0f, -999.0f}, gsj_eta[2] = {-999.0f, -999.0f};
        float gsj_phi[2]= {-999.0f, -999.0f}, gsj_mass[2]= {-999.0f, -999.0f};

        // Variables para JER matching
        float pt_gen_forJER = -999.0f;
        bool  hasGenMatch   = false;

        if (isMC_ && h_genjetsAK8.isValid()) {
            int   bestIdx = -1;
            float bestDR  = 999.0f;
            for (size_t g = 0; g < h_genjetsAK8->size(); ++g) {
                const reco::GenJet& gj = (*h_genjetsAK8)[g];
                float dR = reco::deltaR(eta, phi, gj.eta(), gj.phi());
                if (dR < bestDR) { bestDR = dR; bestIdx = static_cast<int>(g); }
            }

            // [ACCION 3] Criterio de matching JER: dR < R/2 = 0.4 y |dPt|<3*sigma*pT
            if (bestIdx >= 0 && bestDR < 0.4f) {
                const reco::GenJet& gj = (*h_genjetsAK8)[bestIdx];
                float pt_gen_cand = gj.pt();
                // Calcular sigma para el jet reco
                float sigma_jer = -1.0f;
                if (jerResolution_) {
                    JME::JetParameters jparams;
                    jparams.setJetPt(pt);
                    jparams.setJetEta(eta);
                    jparams.setRho(rho);
                    sigma_jer = jerResolution_->getResolution(jparams);
                }
                bool dPtOK = (sigma_jer > 0.0f)
                           ? (std::abs(pt - pt_gen_cand) < 3.0f * sigma_jer * pt)
                           : true;  // si no hay JER, acepta el match geometrico
                if (dPtOK) {
                    hasGenMatch     = true;
                    pt_gen_forJER   = pt_gen_cand;
                    gj_pt = gj.pt(); gj_eta = gj.eta(); gj_phi = gj.phi();
                    gj_mass = gj.mass(); gj_dr = bestDR;

                    if (h_gensubjetsAK8.isValid()) {
                        std::vector<const reco::GenJet*> matched;
                        for (const auto& gsj : *h_gensubjetsAK8)
                            if (reco::deltaR(gj.eta(), gj.phi(), gsj.eta(), gsj.phi()) < 0.8f)
                                matched.push_back(&gsj);
                        std::sort(matched.begin(), matched.end(),
                                  [](const reco::GenJet* a, const reco::GenJet* b){
                                      return a->pt() > b->pt(); });
                        gj_nsub = static_cast<int>(matched.size());
                        math::XYZTLorentzVector gsdSum(0,0,0,0);
                        for (int k = 0; k < std::min(gj_nsub, 2); ++k) {
                            gsj_pt[k]   = matched[k]->pt();
                            gsj_eta[k]  = matched[k]->eta();
                            gsj_phi[k]  = matched[k]->phi();
                            gsj_mass[k] = matched[k]->mass();
                            gsdSum += matched[k]->p4();
                        }
                        if (gj_nsub >= 1) gj_sdmass = static_cast<float>(gsdSum.M());
                    }
                }
            }
        }

        // ----------------------------------------------------------------
        // [ACCION 3] JER smearing hibrido para jet AK8
        // Referencia: CMS-JME-13-004
        // ----------------------------------------------------------------
        float pt_jer_nom  = pt;
        float pt_jer_up   = pt;
        float pt_jer_down = pt;
        float msd_jer_nom  = sd_mass;
        float msd_jer_up   = sd_mass;
        float msd_jer_down = sd_mass;

        if (isMC_ && jerResolution_ && jerSF_) {
            JERFactors jf = computeJERfactors(pt, eta, rho, pt_gen_forJER, hasGenMatch);
            pt_jer_nom   = pt * jf.nom;
            pt_jer_up    = pt * jf.up;
            pt_jer_down  = pt * jf.dn;
            // Propagar el factor JER a m_SD: la masa soft-drop escala con el mismo
            // factor de momentum del jet (aproximacion valida a primer orden).
            // Para una correccion mas precisa usar los subjets gen (gj_sdmass).
            if (sd_mass > 0.0f) {
                msd_jer_nom  = sd_mass * jf.nom;
                msd_jer_up   = sd_mass * jf.up;
                msd_jer_down = sd_mass * jf.dn;
            }
        }

        // ---- Gen-level matching a particula (W/Z/H/top) -----------------
        float gen_pt    = -999.0f;
        float gen_dr    = -999.0f;
        int   is_signal = 0;
        int   parton_fl = jet.partonFlavour();
        int   had_fl    = jet.hadronFlavour();

        if (isMC_ && h_genparts.isValid()) {
            const reco::GenParticle* best = nullptr;
            float minDR = 999.0f;
            for (const auto& gp : *h_genparts) {
                int absId = std::abs(gp.pdgId());
                if (absId != 23 && absId != 24 && absId != 25 && absId != 6) continue;
                if (!gp.isLastCopy()) continue;
                float dR = reco::deltaR(eta, phi, gp.eta(), gp.phi());
                if (dR < minDR) { minDR = dR; best = &gp; }
            }
            if (best && minDR < 0.8f) {
                is_signal = 1;
                gen_dr    = minDR;
                gen_pt    = best->pt();
            }
        }

        // ---- Verdad del W -----------------------------------------------
        float gw_pt = -999.0f, gw_eta = -999.0f, gw_phi = -999.0f, gw_mass = -999.0f;
        float gw_dr = -999.0f;
        int   nQuarksInJet = 0;
        int   isMergedW    = 0;

        if (isMC_ && h_genparts.isValid()) {
            const reco::GenParticle* bestW = nullptr;
            float minDRW = 999.0f;
            for (const auto& gp : *h_genparts) {
                if (std::abs(gp.pdgId()) != 24) continue;
                if (!gp.isLastCopy()) continue;
                float dR = reco::deltaR(eta, phi, gp.eta(), gp.phi());
                if (dR < minDRW) { minDRW = dR; bestW = &gp; }
            }
            if (bestW && minDRW < 0.8f) {
                gw_pt = bestW->pt(); gw_eta = bestW->eta();
                gw_phi = bestW->phi(); gw_mass = bestW->mass();
                gw_dr = minDRW;
                for (size_t d = 0; d < bestW->numberOfDaughters(); ++d) {
                    const reco::Candidate* dau = bestW->daughter(d);
                    if (!dau) continue;
                    if (std::abs(dau->pdgId()) > 5) continue;
                    if (reco::deltaR(eta, phi, dau->eta(), dau->phi()) < 0.8f)
                        ++nQuarksInJet;
                }
                isMergedW = (nQuarksInJet >= 2) ? 1 : 0;
            }
        }

        // ---- Constituyentes PF ------------------------------------------
        struct Cand { const pat::PackedCandidate* c; float wpt; };
        std::vector<Cand> cands;
        cands.reserve(jet.numberOfDaughters());
        for (unsigned int i = 0; i < jet.numberOfDaughters(); ++i) {
            const pat::PackedCandidate* c =
                dynamic_cast<const pat::PackedCandidate*>(jet.daughter(i));
            if (!c) continue;
            cands.push_back({c, static_cast<float>(c->pt() * c->puppiWeight())});
        }
        std::sort(cands.begin(), cands.end(),
                  [](const Cand& a, const Cand& b){ return a.wpt > b.wpt; });
        const size_t NMAX = 150;
        if (cands.size() > NMAX) cands.resize(NMAX);

        std::vector<float> vec_pt, vec_eta, vec_phi, vec_e;
        std::vector<float> vec_dxy, vec_dz, vec_dxy_err, vec_dz_err;
        std::vector<float> vec_puppi, vec_dR_jet;
        std::vector<int>   vec_pdgid, vec_fromPV, vec_charge;

        const size_t NR = cands.size();
        vec_pt.reserve(NR);      vec_eta.reserve(NR);
        vec_phi.reserve(NR);     vec_e.reserve(NR);
        vec_dxy.reserve(NR);     vec_dz.reserve(NR);
        vec_dxy_err.reserve(NR); vec_dz_err.reserve(NR);
        vec_puppi.reserve(NR);   vec_dR_jet.reserve(NR);
        vec_pdgid.reserve(NR);   vec_fromPV.reserve(NR);
        vec_charge.reserve(NR);

        for (const auto& cd : cands) {
            const pat::PackedCandidate* c = cd.c;
            float deta  = c->eta() - eta;
            float dphi  = reco::deltaPhi(c->phi(), phi);
            float dR_cj = std::sqrt(deta*deta + dphi*dphi);

            vec_pt.push_back(c->pt());
            vec_eta.push_back(c->eta());
            vec_phi.push_back(c->phi());
            vec_e.push_back(c->energy());
            vec_pdgid.push_back(c->pdgId());
            vec_charge.push_back(c->charge());
            vec_dR_jet.push_back(dR_cj);

            if (c->hasTrackDetails()) {
                vec_dxy.push_back(c->dxy());
                vec_dz.push_back(c->dz());
                vec_dxy_err.push_back(c->dxyError());
                vec_dz_err.push_back(c->dzError());
            } else {
                vec_dxy.push_back(0.0f);   vec_dz.push_back(0.0f);
                vec_dxy_err.push_back(0.0f); vec_dz_err.push_back(0.0f);
            }
            vec_fromPV.push_back(static_cast<int>(c->fromPV()));
            vec_puppi.push_back(c->puppiWeight());
        }
        int n_const = static_cast<int>(vec_pt.size());

        // ---- Construir fila CSV -----------------------------------------
        std::ostringstream row;
        row << std::fixed << std::setprecision(6);

        row // [1] Nivel evento
            << eventID          << ","
            << runNum           << ","
            << npv              << ","
            << pileup_nTrueInt  << ","
            << rho              << ","
            << genWeight        << ","
            << ptHat            << ","
            // [2] Cinematica
            << pt               << ","
            << pt_raw           << ","
            << area             << ","
            << eta              << ","
            << phi              << ","
            << mass             << ","
            << mass_raw         << ","
            << sd_mass          << ","
            // [3] Substructura
            << t1               << ","
            << t2               << ","
            << t3               << ","
            << tau21            << ","
            << tau32            << ","
            << n2b1             << ","
            << n2b2             << ","
            << d2               << ","
            << nsubjets         << ","
            // [4] JEC AK8 (pT)
            << jec_unc          << ","
            << pt_jec_up        << ","
            << pt_jec_down      << ","
            // [4b] SDMass con JEC de SUBJET --- CORRECTA [ACCION 2]
            << msd_sjJEC_up     << ","
            << msd_sjJEC_down   << ","
            // [4c] SDMass con JEC AK8 --- APROXIMACION
            << msd_jec_approx_up   << ","
            << msd_jec_approx_down << ","
            // [4d] JER [ACCION 3]
            << pt_jer_nom       << ","
            << pt_jer_up        << ","
            << pt_jer_down      << ","
            << msd_jer_nom      << ","
            << msd_jer_up       << ","
            << msd_jer_down     << ","
            << (hasGenMatch ? 1 : 0) << ","
            // [5] Subjets reco crudos (+ptNom [v5-C])
            << sd_mass_raw      << ","
            << sj_ptRaw[0]      << "," << sj_ptNom[0] << ","
            << sj_eta[0]        << "," << sj_phi[0]   << ","
            << sj_massRaw[0]    << "," << sj_area[0]  << "," << sj_deepCSV[0] << ","
            << sj_ptRaw[1]      << "," << sj_ptNom[1] << ","
            << sj_eta[1]        << "," << sj_phi[1]   << ","
            << sj_massRaw[1]    << "," << sj_area[1]  << "," << sj_deepCSV[1] << ","
            << sj_dR12          << "," << sj_ptAsym   << ","
            // [6] Gen jet AK8
            << gj_pt            << "," << gj_eta      << "," << gj_phi        << ","
            << gj_mass          << "," << gj_sdmass   << "," << gj_nsub       << ","
            << gj_dr            << ","
            << gsj_pt[0]        << "," << gsj_eta[0]  << "," << gsj_phi[0]    << ","
            << gsj_mass[0]      << ","
            << gsj_pt[1]        << "," << gsj_eta[1]  << "," << gsj_phi[1]    << ","
            << gsj_mass[1]      << ","
            // [7] Gen particle / flavor
            << gen_pt           << ","
            << gen_dr           << ","
            << parton_fl        << ","
            << had_fl           << ","
            << is_signal        << ","
            // [8] Verdad del W
            << gw_pt            << "," << gw_eta      << "," << gw_phi        << ","
            << gw_mass          << "," << gw_dr        << ","
            << nQuarksInJet     << "," << isMergedW    << ","
            // [9] Constituyentes PF
            << n_const          << ","
            << "\"" << vecF(vec_pt)      << "\","
            << "\"" << vecF(vec_eta)     << "\","
            << "\"" << vecF(vec_phi)     << "\","
            << "\"" << vecF(vec_e)       << "\","
            << "\"" << vecI(vec_pdgid)   << "\","
            << "\"" << vecF(vec_dxy)     << "\","
            << "\"" << vecF(vec_dz)      << "\","
            << "\"" << vecF(vec_dxy_err) << "\","
            << "\"" << vecF(vec_dz_err)  << "\","
            << "\"" << vecI(vec_fromPV)  << "\","
            << "\"" << vecF(vec_puppi)   << "\","
            << "\"" << vecI(vec_charge)  << "\","
            << "\"" << vecF(vec_dR_jet)  << "\"";

        // Chequeo unico de consistencia cabecera <-> fila
        if (!firstRowChecked_) {
            firstRowChecked_ = true;
            size_t nRow = countFields(row.str());
            if (nRow != nHeaderFields_) {
                edm::LogError("MyAnalyzer")
                    << "DESAJUSTE CSV: cabecera=" << nHeaderFields_
                    << " campos, fila=" << nRow << " campos.";
            } else {
                edm::LogInfo("MyAnalyzer")
                    << "Consistencia CSV verificada: " << nRow << " campos.";
            }
        }

        csvFile_ << row.str() << "\n";

    } // fin for jets
}

DEFINE_FWK_MODULE(MyAnalyzer);
