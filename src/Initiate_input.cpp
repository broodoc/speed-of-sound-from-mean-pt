#include "StParameters.h"
#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>
#include <thread>
#include <chrono>
#include "TSystem.h"

void Initiate_input() {
    auto start_time = std::chrono::steady_clock::now();
    
    std::streambuf* cout_buf = std::cout.rdbuf();
    std::streambuf* cerr_buf = std::cerr.rdbuf();
    
    StParameters *params = new StParameters();
    
    if (!params->LoadParamTable()) {
        std::cout << "ERROR: Unable to load input table!" << std::endl;
        delete params;
        return;
    }

    // Get directory name from paramMap
    if (params->paramMap["Directory"].empty()) {
        std::cout << "ERROR: Directory name not specified in paramMap!" << std::endl;
        delete params;
        return;
    }

    std::string subdir = params->paramMap["PrepDir"][0];
    std::string fullPath = "/Users/calebbroodo/Desktop/sos_analysis/prepped_files/" + subdir;
    
    namespace fs = std::filesystem;

    if (fs::exists(fullPath)) {
            std::cout << "⚠️  Directory already exists: " << fullPath << std::endl;
            std::cout << "To prevent overwriting press Ctrl + C within the next 3 seconds";
        for (int i = 0; i < 6; ++i) {
            std::this_thread::sleep_for(std::chrono::seconds(1));
            if (i == 0) std::cout << "...[3]" << std::flush;
            if (i == 1) std::cout << "...[2]" << std::flush;
            if (i == 2) std::cout << "...[1]" << std::flush;
            if (i == 3) std::cout << "...[0]" << std::flush;
            if (i == 4) std::cout << "\nproceeding..." << std::endl;
        }
            
    }
    
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);

    std::ostringstream oss;
    oss << fullPath << "/log_" << std::put_time(&tm, "%Y-%b-%d_%I-%M%p") << ".txt";
    std::string logFileName = oss.str();

    std::ofstream logFile(logFileName);
    if (!logFile.is_open()) {
        std::cerr << "❌ Error opening log file for writing: " << logFileName << std::endl;
    } else {
        std::cout.rdbuf(logFile.rdbuf()); // Redirect std::cout to log file
        std::cerr.rdbuf(logFile.rdbuf()); // Redirect std::cerr to log file
        std::cout << "📄 Logging output to: " << logFileName << std::endl;
    }

    try {
        fs::create_directories(fullPath);
        std::cout << "✅ Created directory: " << fullPath << std::endl;
        
        std::string inputPath = "/Users/calebbroodo/Desktop/sos_analysis/input_table.txt";
        std::string destinationPath = fullPath + "/input_table.txt";
        
        fs::copy_file(inputPath, destinationPath, fs::copy_options::overwrite_existing);
        std::cout << "📄 Copied input_table.txt to: " << destinationPath << std::endl;
        
        std::string inputrunPath = "/Users/calebbroodo/Desktop/sos_analysis/run_over_input_table.sh";
        std::string destinationrunPath = fullPath + "/run_over_input_table.sh";
        
        fs::copy_file(inputrunPath, destinationrunPath, fs::copy_options::overwrite_existing);
        std::cout << "📄 Copied run_over_input_table.sh to: " << destinationPath << std::endl;
        
        std::string dir = params->paramMap["Directory"][0];
        
        if (fs::exists(dir)) {
            std::string referenceFile = fullPath + "/RCF_DST_reference.txt";
            
            std::ofstream out(referenceFile);
            if (out.is_open()) {
                out << dir << std::endl;
                out.close();
                std::cout << "📁 RCF_DST_reference.txt written in: " << referenceFile << std::endl;
            } else {
                std::cerr << "❌ Failed to write to: " << referenceFile << std::endl;
            }
        }
        
        std::string refMult = params->paramMap["Multname"][0];
        std::string multfromdata = params->paramMap["MultFromData"][0];
        std::string pTMult = params->paramMap["pTMultname"][0];
        std::string InEx = params->paramMap["Inclusive/Exclusive(I/E)"][0];
        std::string DebugVerbose = params->paramMap["DebugVerbose(Y/N)"][0];
        std::string overwrite = params->paramMap["Overwrite(Y/N)"][0];
        std::string Glauber = params->paramMap["Glauber(Y/N)"][0];
        std::string DisplayFitRanges = params->paramMap["DisplayFitRanges(Y/N)"][0];
        std::string GenerateFigures = params->paramMap["GenerateFigures(Y/N)"][0];
        std::string MomentAna = params->paramMap["MomentAna(Y/N)"][0];
        
        // Before first script
        auto start_clean = std::chrono::steady_clock::now();
        std::string command = "./clean_and_merge_DSTs.sh " + dir + " "+ refMult + " " + pTMult + " " + Glauber + " " + multfromdata + " 2>&1 | tee -a " + logFileName;
        int exitCode_cl = gSystem->Exec(command.c_str());
        auto end_clean = std::chrono::steady_clock::now();
        auto duration_clean = std::chrono::duration_cast<std::chrono::seconds>(end_clean - start_clean);
        std::cout << "⏱️  cleaning and merging runtime: " << duration_clean.count() << " sec" << std::endl;

        // Before second script
        auto start_cent = std::chrono::steady_clock::now();
        std::string command_cd = "./determine_centrality.sh " + dir + " " + fullPath + " " + refMult + " " + overwrite + " 2>&1 | tee -a " + logFileName;
        int exitCode_cd = gSystem->Exec(command_cd.c_str());
        auto end_cent = std::chrono::steady_clock::now();
        auto duration_cent = std::chrono::duration_cast<std::chrono::seconds>(end_cent - start_cent);
        std::cout << "⏱️  centrality determination runtime: " << duration_cent.count() << " sec" << std::endl;
        
        // Before third script
        auto start_stats = std::chrono::steady_clock::now();
        std::string command_sb = "./determine_statistics.sh " + dir + " " + fullPath + " " + multfromdata + " " + pTMult + " " + InEx + " " + DebugVerbose + " " + " 2>&1 | tee -a " + logFileName;
//        std::string command_sb = "./determine_statistics.sh " + dir + " " + fullPath + " " + refMult + " " + pTMult + " " + InEx + " " + DebugVerbose + " " + " 2>&1 | tee -a " + logFileName;
        int exitCode_sb = gSystem->Exec(command_sb.c_str());
        auto end_stats = std::chrono::steady_clock::now();
        auto duration_stats = std::chrono::duration_cast<std::chrono::seconds>(end_stats - start_stats);
        std::cout << "⏱️ statistical bootstrapping runtime: " << duration_stats.count() << " sec" << std::endl;

        
        // Before fourth script
        auto start_ps = std::chrono::steady_clock::now();
        std::string command_ps = "./determine_pdf_shading.sh " + dir + " " + fullPath + " " + refMult + " 2>&1 | tee -a " + logFileName;
        int exitCode_ps = gSystem->Exec(command_ps.c_str());
        auto end_ps = std::chrono::steady_clock::now();
        auto duration_ps = std::chrono::duration_cast<std::chrono::seconds>(end_ps - start_ps);
        std::cout << "⏱️ pdf shading runtime: " << duration_stats.count() << " sec" << std::endl;
        
        
        // Before fifth script
        auto start_sos = std::chrono::steady_clock::now();
        std::string command_sos = "./determine_sos.sh " + dir + " " + fullPath + " " + refMult + " " + pTMult + " " + InEx + " " + DebugVerbose + " " + DisplayFitRanges + " " + GenerateFigures + " " + MomentAna + " " +" 2>&1 | tee -a " + logFileName;
        int exitCode_sos = gSystem->Exec(command_sos.c_str());
        auto end_sos = std::chrono::steady_clock::now();
        auto duration_sos = std::chrono::duration_cast<std::chrono::seconds>(end_sos - start_sos);
        std::cout << "⏱️  speed of sound analysis runtime: " << duration_sos.count() << " sec" << std::endl;

    } catch (const fs::filesystem_error& e) {
            std::cerr << "❌ Filesystem error: " << e.what() << std::endl;
        }
    
    auto end_time = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    
    logFile.close();
    std::cout.rdbuf(cout_buf);
    std::cerr.rdbuf(cerr_buf);

    int minutes = duration.count() / 60;
    int seconds = duration.count() % 60;

    std::cout << "⏱️  Total runtime: ";
    if (minutes > 0)
        std::cout << minutes << " min ";
    std::cout << seconds << " sec" << std::endl;

    

    delete params;
}
