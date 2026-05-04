#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TString.h>
#include <TH1.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm> 

const int numBatches = 10;

struct FileInfo {
    std::string fileName;
    int numEntries;
};

bool compareByEntries(const FileInfo& a, const FileInfo& b) {
    return a.numEntries > b.numEntries; // Sort in descending order of entries
}

void SubSampler(const char* directoryPath) {
    TSystemDirectory dir(directoryPath, directoryPath);
    TList* files = dir.GetListOfFiles();

    if (!files) {
        std::cerr << "Error: Directory not found or no root files found in the directory." << std::endl;
        return;  // Return to avoid further processing if no files are found
    }

    std::vector<FileInfo> fileInfos;

    for (int i = 0; i < files->GetSize(); ++i) {
        const char* fileName = files->At(i)->GetName();
        std::string filePath = TString::Format("%s/%s", directoryPath, fileName).Data();

        TFile file(filePath.c_str(), "READ");
        if (file.IsZombie()) {
            std::cerr << "Error opening file: " << filePath << std::endl;
            continue;
        }
        
        // Navigate to the Correlator directory within the file
        TDirectory* correlatorDir = dynamic_cast<TDirectory*>(file.Get("Correlator"));
        if (!correlatorDir) {
            std::cerr << "Error accessing Correlator directory in file: " << filePath << std::endl;
            file.Close();
            continue;
        }

        TH1* histogram = dynamic_cast<TH1*>(correlatorDir->Get("Cm2_MBin0_Eta_0")); // Replace with your histogram name
        if (!histogram) {
            std::cerr << "Error accessing histogram in file: " << filePath << std::endl;
            file.Close();
            continue;
        }
       
        FileInfo fileInfo;
        fileInfo.fileName = filePath;
        fileInfo.numEntries = static_cast<int>(histogram->GetEntries());

        fileInfos.push_back(fileInfo);

        // Close the file
        file.Close();
    }

    // Sort files based on the number of entries
    std::sort(fileInfos.begin(), fileInfos.end(), compareByEntries);
  
    // Calculate total entries
    int totalEntries = 0;
    for (size_t i = 0; i < fileInfos.size(); ++i) {
        totalEntries += fileInfos[i].numEntries;
    }

    std::cout << "Total entries = " << totalEntries << std::endl;

    // Calculate entries per batch
    int entriesPerBatch = totalEntries / numBatches;
    std::cout << "Entries per batch: " << entriesPerBatch << std::endl;

    // Create batches
    std::vector<std::vector<FileInfo>> batches(numBatches);
    int currentBatch = 0;
    int currentBatchEntries = 0;

    for (size_t i = 0; i < fileInfos.size(); ++i) {
        if (currentBatchEntries + fileInfos[i].numEntries <= entriesPerBatch) {
            batches[currentBatch].push_back(fileInfos[i]);
            currentBatchEntries += fileInfos[i].numEntries;
        } else {
            if (currentBatch + 1 < numBatches) {
                currentBatch++;
                currentBatchEntries = fileInfos[i].numEntries;
                batches[currentBatch].push_back(fileInfos[i]);
            } else {
                // Skip files to prevent overflow
            }
        }
    }

    std::vector<double> means(numBatches, 0.0);
    std::vector<double> means1(numBatches, 0.0);

    for (int i = 0; i < numBatches; ++i) {
        std::cout << "Batch " << i + 1 << " (" << batches[i].size() << " files, "
                  << "Total Entries: " << entriesPerBatch << "):" << std::endl;

        for (size_t j = 0; j < batches[i].size(); ++j) {
            TFile file(batches[i][j].fileName.c_str(), "READ");
            TH1* histogram = dynamic_cast<TH1*>(file.Get("Cm1_MBin0_Eta_0")); // Replace with your histogram name
            TH1* histogram1 = dynamic_cast<TH1*>(file.Get("Cm2_MBin0_Eta_0"));
            if (histogram) {
                means[i] += histogram->GetMean();
                means1[i] += histogram1->GetMean();
            }
            // Close the file
            file.Close();
        }

        // Calculate mean for the batch
        if (batches[i].size() > 0) {
            means[i] /= batches[i].size();
            means1[i] /= batches[i].size();
            std::cout << "Batch Mean1: " << means[i] << std::endl;
            std::cout << "Batch Mean2: " << means1[i] << std::endl;
            std::cout << "Batch size = " << batches[i].size() << std::endl;
        }

        std::cout << std::endl;
    }

    // Display means array
    std::cout << "Means array: [";
    for (size_t i = 0; i < means.size(); ++i) {
        std::cout << sqrt(means[i] - pow(means1[i], 2)) / means1[i] * 1;
        if (i < means.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "]" << std::endl;

    // Clean up
    delete files;
}
