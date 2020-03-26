/**
 * @file BhFitSandboxProcessor.cxx
 * @brief Processor used to throw toys to study BH bkg model complexity
 * @author Cameron Bravo, SLAC National Accelerator Laboratory
 */

#include "BhFitSandboxProcessor.h"
#include "TString.h"

BhFitSandboxProcessor::BhFitSandboxProcessor(const std::string& name, Process& process)
    : Processor(name, process) { 
    }

BhFitSandboxProcessor::~BhFitSandboxProcessor() { }

void BhFitSandboxProcessor::configure(const ParameterSet& parameters) {
    std::cout << "Configuring BhFitSandboxProcessor" << std::endl;
    try {
        debug_               = parameters.getInteger("debug");
        massSpectrum_        = parameters.getString("massSpectrum");
        mass_hypo_           = parameters.getDouble("mass_hypo");
        win_factor_          = parameters.getInteger("win_factor");
        poly_order_          = parameters.getInteger("poly_order");
        seed_                = parameters.getInteger("seed");
        nToys_               = parameters.getInteger("nToys");
        toy_sig_samples_     = parameters.getInteger("toy_sig_samples");
        bkg_mult_            = parameters.getInteger("toy_bkg_mult");
        res_scale_           = parameters.getDouble("res_scale");
        signal_shape_h_name_ = parameters.getString("signal_shape_h_name", "");
        signal_shape_h_file_ = parameters.getString("signal_shape_h_file", "");
        //asymptotic_limit_ = parameters.getBoolean("asymptoticLimit");
    } catch(std::runtime_error& error) {
        std::cout << error.what() << std::endl;
    }
}

void BhFitSandboxProcessor::initialize(std::string inFilename, std::string outFilename) {
    // Init Files
    outF_ = new TFile(outFilename.c_str(),"RECREATE");
    outF_->cd();

    // Get mass spectrum from file
    mass_spec_h = new TH1F("mass_spec_h", "mass_spec_h", 6000, 0.0, 0.3);
    mass_spec_h->Sumw2();
    for(int iBin = 0; iBin < 6000; iBin++)
    {
        for(int ff = 0; ff < 100; ff++) mass_spec_h->Fill(mass_spec_h->GetBinCenter(iBin+1));
    }
    mass_spec_h->Write();

    // Initialize the signal histogram, if a file name and histogram name are provided.
    std::cout << "Signal Shape File :: " << signal_shape_h_file_ << std::endl;
    std::cout << "Signal Shape Hist :: " << signal_shape_h_name_ << std::endl;
    if(signal_shape_h_file_ != "" && signal_shape_h_name_ != "") {
        TFile *file = new TFile(signal_shape_h_file_.c_str());
        signal_shape_h_ = (TH1*) file->Get(signal_shape_h_name_.c_str());
    } else if(signal_shape_h_file_ != "" && signal_shape_h_name_ == "") {
        std::cout << "[BumpHunter] :: !! WARNING !! Signal injection file, but no histogram, specified! Defaulting to Gaussian.";
    } else if(signal_shape_h_file_ == "" && signal_shape_h_name_ != "") {
        std::cout << "[BumpHunter] :: !! WARNING !! Signal injection histogram, but no file, specified! Defaulting to Gaussian.";
    }

    // Init bump hunter manager
    bump_hunter_ = new BumpHunter(bkg_model_, poly_order_, win_factor_, res_scale_, asymptotic_limit_);
    bump_hunter_->setBounds(mass_spec_h->GetXaxis()->GetBinUpEdge(mass_spec_h->FindFirstBinAbove()),
            mass_spec_h->GetXaxis()->GetBinLowEdge(mass_spec_h->FindLastBinAbove()));
    if(debug_ > 0) bump_hunter_->enableDebug();

}

bool BhFitSandboxProcessor::process() {

    std::vector<TH1*> toy_hists;
    std::vector<HpsFitResult*> toy_results;
    TRandom3 rand;
    rand.SetSeed(0);
    for(int iToy = 0; iToy < nToys_; iToy++)
    {
        std::cout << "Generating Toy " << iToy << std::endl;
        toy_hists.push_back(new TH1F(Form("toy%i_h", iToy), Form("toy%i_h", iToy), 6000, 0.0, 0.3));
        toy_hists[iToy]->Sumw2();
        for(int iBin = 0; iBin < 6000; iBin++)
        {
            int nFills = rand.Poisson(rand.Uniform(100000));
            for(int ff = 0; ff < nFills; ff++) toy_hists[iToy]->Fill(toy_hists[iToy]->GetBinCenter(iBin+1));
        }

    }

    int toyFitN = 0;
    for(TH1* hist : toy_hists) {
        std::cout << "Fitting Toy " << toyFitN << std::endl;
        toy_results.push_back(bump_hunter_->performSearch(hist, mass_hypo_, false, false));
        hist->Write();
        toyFitN++;
    }

    return true;
}

void BhFitSandboxProcessor::finalize() { 
    outF_->Close();
    delete outF_;
    delete bump_hunter_;
}

DECLARE_PROCESSOR(BhFitSandboxProcessor); 
