#include "PlotsMaker.h"
#include "PlotsMaker.h"
#include "TApplication.h"
#include <TCanvas.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TGraph.h>
#include "PlotsMaker.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TF2.h"
#include "TMath.h"
#include <TAxis.h>
#include <TFormula.h>

#include <iostream>
#include <vector>
#include <string>
#include <TCanvas.h>
#include <TGraph.h>
#include <TF1.h>
#include <TLegend.h>
#include <TError.h>
#include <TString.h>
#include <TStyle.h>
#include <TApplication.h>

void PlotsMaker::MakeHistogram(std::vector<int> info, std::string filename){
    std::string destination = "../res/" + filename;
    // Create a ROOT canvas
    TCanvas *c1 = new TCanvas("c1", "Histogram", 800, 600);
    // Create a histogram
    TH1F *hist = new TH1F("hist", filename.c_str(), info.size(), 0, info.size());
    // Fill the histogram with the info
    for (int i = 0; i < info.size(); ++i) {
        hist->SetBinContent(i+1, info[i]);
    }
    // Set histogram style
    hist->SetFillColor(38);
    // Draw the histogram
    hist->Draw();
    // Update the canvas to display the histogram
    c1->Update();
    c1->SaveAs(destination.c_str());
    // Wait for user input to exit the application
    c1->WaitPrimitive();
    delete hist;
    delete c1;
}

void PlotsMaker::MakeHistogram(std::vector<double> info, std::string filename, std::string title, std::string xAxisLabel, std::string yAxisLabel) {
    std::string destination = "../res/" + filename;
    
    // Create a ROOT canvas
    TCanvas *c1 = new TCanvas("c1", title.c_str(), 800, 600);
    
    // Create a histogram
    TH1F *hist = new TH1F("hist", title.c_str(), info.size(), 0, info.size());
    
    // Fill the histogram with the info
    for (int i = 0; i < info.size(); ++i) {
        hist->SetBinContent(i+1, info[i]);
    }
    
    // Set histogram style
    hist->SetFillColor(38);

    hist->SetStats(0);
    
    // Draw the histogram
    hist->Draw();
    
    // Set title and axis labels
    hist->SetTitle(title.c_str());
    hist->GetXaxis()->SetTitle(xAxisLabel.c_str());
    hist->GetYaxis()->SetTitle(yAxisLabel.c_str());
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();
    hist->GetYaxis()->SetTitleOffset(1.5);
    hist->GetYaxis()->SetTitleOffset(1.5);
    
    // Add a grid to the background
    gPad->SetGrid();
    
    // Update the canvas to display the histogram
    c1->Update();
    c1->SaveAs(destination.c_str());
    
    // Wait for user input to exit the application
    c1->WaitPrimitive();
    
    // Clean up
    delete hist;
    delete c1;
}

void PlotsMaker::MakeHistogram2D(std::vector<std::vector<int>> info, std::string filename){
    std::string destination = "../res/" + filename;
    
    // Determine the dimensions of the histogram
    int nx = info.size(); // Number of bins along x-axis
    int ny = info[0].size(); // Number of bins along y-axis

    // Create a ROOT canvas
    TCanvas *c1 = new TCanvas("c1", "2D Histogram", nx, ny);
    
    // Set grayscale color palette
    gStyle->SetPalette(kGreyScale);
    
    // Create a 2D histogram
    TH2F *hist2D = new TH2F("hist2D", filename.c_str(), nx, 0, nx, ny, 0, ny);
    hist2D->SetMinimum(-0.1);
    
    // Fill the histogram with the info
    for (int col = 0; col < nx; col++) {
        for (int row = 0; row < ny; row++) {
            hist2D->Fill(col+0.5, ny-row-0.5, info[col][row]);
        }
    }

    int total = nx*ny;
    int count = 0;
    for (int col = 0; col < nx; col++) {
        for (int row = 0; row < ny; row++) {
            if(info[col][row] == hist2D->GetBinContent(col+0.5, ny-row-0.5)){count += 1;}
        }
    }
    
    std::cout << "TOTAL EQUALITY: " << (double)count/(double)total << std::endl;


    // Draw the 2D histogram
    hist2D->Draw("COLZ");
    
    // Update the canvas to display the histogram
    c1->Update();
    
    // Save the canvas as an image file
    c1->SaveAs(destination.c_str());
    
    // Wait for user input to exit the application
    c1->WaitPrimitive();
    
    // Clean up
    delete hist2D;
    delete c1;
}

void PlotsMaker::MakePointsPlot(std::vector<int>info , std::string filename){
    std::cout << "Empty" << std::endl;
}

void PlotsMaker::MakePlotWithStringX(std::map<std::string, std::pair<int, double>> data, std::string graphTitle, std::string xAxisLabel, std::string yAxisLabel, std::string filename){
    TCanvas *canvas = new TCanvas("c1", "Graph", 800, 600);
    TGraph *graph = new TGraph(data.size());

    // Add data points to the graph
    for (const auto& entry : data) {
        std::cout << "Adding point at index "<< ": X = " << entry.second.first << ", Y = " << entry.second.second << std::endl;
        graph->SetPoint(graph->GetN(), entry.second.first, entry.second.second); // SetPoint(index, x, y)
    }

    std::cout << "Number of points added: " << graph->GetN() << std::endl;
    
    graph->SetTitle(graphTitle.c_str());
    graph->GetXaxis()->SetTitle(xAxisLabel.c_str());
    graph->GetYaxis()->SetTitle(yAxisLabel.c_str());
    canvas->cd();
    graph->Draw("AP");
    graph->SetMarkerSize(0.5); // Set marker size (adjust as desired)
    canvas->SetGrid();
    canvas->Update();

    // Save canvas as an image file
    std::string destination = "../res/" + filename;
    canvas->SaveAs(destination.c_str());
    canvas->Draw();

    // Clean up
    delete graph;
    delete canvas;
}

void PlotsMaker::Draw(const std::vector<double>& x, const std::vector<double>& y, std::string s,
                        std::string title, std::string xAxisLabel, std::string yAxisLabel) {
    // Create a canvas
    TCanvas* canvas = new TCanvas("canvas", "Plot", 800, 600);

    // Create a TGraph object with the given x and y vectors
    TGraph* graph = new TGraph(x.size(), &x[0], &y[0]);

    // Set the title if provided
    if (!s.empty()) {
        graph->SetTitle(s.c_str());
    }

    // Draw the graph with points
    graph->SetMarkerStyle(20); // Set marker style
    graph->SetMarkerSize(0.8); // Set marker size
    graph->Draw("AP"); // Draw the graph with axes and points
    graph->SetTitle(title.c_str());
    graph->GetXaxis()->SetTitle(xAxisLabel.c_str());
    graph->GetYaxis()->SetTitle(yAxisLabel.c_str());

    graph->GetXaxis()->CenterTitle();
    graph->GetYaxis()->CenterTitle();
    graph->GetYaxis()->SetTitleOffset(1.5);
    graph->GetYaxis()->SetTitleOffset(1.5);
    

    gPad->SetGrid();

    // Update the canvas
    canvas->Update();

    // Save the plot as a PDF
    std::string dest = "../res/" + s;
    canvas->SaveAs(dest.c_str());

    // Display the canvas
    canvas->Draw();
    delete canvas;
    delete graph;
}

void PlotsMaker::MakeLinearFitPlot(const std::vector<double>& x, const std::vector<double>& error_x, 
                                   const std::vector<double>& y, const std::vector<double>& error_y,
                                   const std::string& filename, const std::string& title,
                                   const std::string& xAxisLabel, const std::string& yAxisLabel) {
    std::string destination = "../res/" + filename;

    // Create a canvas
    TCanvas* canvas = new TCanvas("canvas", "Linear Fit Plot", 800, 600);

    // Create TGraphErrors object with the given x, y, error_x, and error_y arrays
    TGraphErrors* graph = new TGraphErrors(x.size(), &x[0], &y[0], &error_x[0], &error_y[0]);

    // Perform linear fit
    TF1* linearFit = new TF1("linearFit", "pol1", x.front(), x.back());
    graph->Fit("linearFit", "Q");

    // Set titles and labels
    graph->SetTitle(title.c_str());
    graph->GetXaxis()->SetTitle(xAxisLabel.c_str());
    graph->GetYaxis()->SetTitle(yAxisLabel.c_str());

    // Center the axis labels
    graph->GetXaxis()->CenterTitle();
    graph->GetYaxis()->CenterTitle();

    // Set title offsets
    graph->GetXaxis()->SetTitleOffset(1.2);
    graph->GetYaxis()->SetTitleOffset(1.2);

    // Set grid
    gPad->SetGrid();

    // Draw the graph with points and error bars
    graph->SetMarkerStyle(20); // Set marker style
    graph->SetMarkerSize(1.0); // Set marker size
    graph->Draw("AP"); // Draw the graph with axes, points, and error bars

    // Draw the linear fit curve
    linearFit->Draw("SAME");

    // Update the canvas
    canvas->Update();

    // Save the plot as a PDF
    canvas->SaveAs(destination.c_str());

    // Display the canvas
    canvas->Draw();

    // Print fit results
    double slope = linearFit->GetParameter(1);
    double slopeError = linearFit->GetParError(1);
    double intercept = linearFit->GetParameter(0);
    double interceptError = linearFit->GetParError(0);

    std::cout << "Linear Fit Results for: " << filename << std::endl;
    std::cout << "Slope: " << slope << " ± " << slopeError << std::endl;
    std::cout << "Intercept: " << intercept << " ± " << interceptError << std::endl;
    std::cout <<"===================="<< std::endl;

    // Clean up
    delete canvas;
    delete graph;
}

void PlotsMaker::MakeModelFitPlot(const std::vector<double>& v_x, const std::vector<double>& v_y,
                                   const std::string& filename, 
                                   const std::string& title, const std::string& xAxisLabel, 
                                   const std::string& yAxisLabel) {
    std::string destination = "../res/" + filename;

    // Create a canvas
    TCanvas* canvas = new TCanvas("canvas", "Model Fit Plot", 800, 600);
    
    // Create a TGraph from the v_x and v_y data
    TGraph* graph = new TGraph(v_x.size(), &v_x[0], &v_y[0]);

    // Define the fit function for the sum of three Gaussians
    TF1* modelFit = new TF1("modelFit", 
        "[0]*exp(-0.5*((x-[1])^2)/[2]^2) + [3]*exp(-0.5*((x-[4])^2)/[5]^2) + [6]*exp(-0.5*((x-[7])^2)/[8]^2) + [9]", 
        v_x.front(), v_x.back());

    // Parameters: Amplitude1, Mean1, Stddev1, Amplitude2, Mean2, Stddev2, Amplitude3, Mean3, Stddev3, Y-Translation
    modelFit->SetParameters(5.5, 3593*10e5, 6.36117e+06,  // Gaussian 1: Amplitude, Mean, Stddev
                            3.7, 3613*10e5, 9.41019e+06,  // Gaussian 2: Amplitude, Mean, Stddev
                            1.0, 3600*10e5, 5e6,          // Gaussian 3: Amplitude, Mean, Stddev
                            -26);                         // Y-Translation

    // Fix parameters [2] and [5]
    modelFit->FixParameter(2, 6.36117e+06); // Fix Stddev1
    modelFit->FixParameter(5, 9.41019e+06); // Fix Stddev2

    //Set limits for each parameter
    modelFit->SetParLimits(0, 2, 6);             // Amplitude for Gaussian 1
    modelFit->SetParLimits(1, 3580*10e5, 3600*10e5); // Mean for Gaussian 1

    modelFit->SetParLimits(3, 3, 5);             // Amplitude for Gaussian 2
    modelFit->SetParLimits(4, 3600*10e5, 3620*10e5); // Mean for Gaussian 2

    modelFit->SetParLimits(6, -10, 15);             // Amplitude for Gaussian 3
    modelFit->SetParLimits(7, 3560*10e5, 3650*10e5); // Mean for Gaussian 3
    modelFit->SetParLimits(8, 8e5, 7e7);         // Stddev for Gaussian 3
    modelFit->SetParLimits(9, -45, -20);         // Y-Translation

    // Fit the model to the data
    int fitStatus = graph->Fit("modelFit", "SR"); // "S" saves status, "R" restricts range

    // Check fit status
    if (fitStatus != 0) {
        std::cerr << "Fit failed for: " << filename << " (status: " << fitStatus << ")" << std::endl;
        delete canvas;
        delete graph;
        delete modelFit;
        return;
    }

    // Set titles and labels
    graph->SetTitle(title.c_str());
    graph->GetXaxis()->SetTitle(xAxisLabel.c_str());
    graph->GetYaxis()->SetTitle(yAxisLabel.c_str());

    // Center the axis labels
    graph->GetXaxis()->CenterTitle();
    graph->GetYaxis()->CenterTitle();

    // Set grid for better visualization
    gPad->SetGrid();

    // Draw the graph with points
    graph->SetMarkerStyle(20);  // Set marker style
    graph->SetMarkerSize(1.0);  // Set marker size
    graph->Draw("AP");          // Draw the graph with axes and points

    // Update the canvas
    canvas->Update();

    // Save the plot as a PDF
    canvas->SaveAs(destination.c_str());

    // Print the fit parameters with errors
    for (int i = 0; i < 9; i += 3) {
        double amplitude = modelFit->GetParameter(i);
        double mean = modelFit->GetParameter(i + 1);
        double stddev = modelFit->GetParameter(i + 2);
        double amplitudeError = modelFit->GetParError(i);
        double meanError = modelFit->GetParError(i + 1);
        double stddevError = (i == 2 || i == 5) ? 0 : modelFit->GetParError(i + 2); // Error is 0 for fixed params

        std::cout << "Gaussian " << (i / 3 + 1) << " Parameters:" << std::endl;
        std::cout << "Amplitude: " << amplitude << " ± " << amplitudeError << std::endl;
        std::cout << "Mean: " << mean << " ± " << meanError << std::endl;
        std::cout << "Stddev: " << stddev << " ± " << stddevError << std::endl;
        std::cout << "====================" << std::endl;
    }
    double Yparam = modelFit->GetParameter(9);
    double Yerror = modelFit->GetParError(9);
    std::cout << "Y-Translation: " << Yparam << " ± " << Yerror << std::endl;

    // Clean up
    delete canvas;
    delete graph;
    delete modelFit;
}

// void PlotsMaker::MakeModelFitPlot(const std::vector<double>& v_x, const std::vector<double>& v_y,
//                                    const std::vector<double>& initialParams, 
//                                    const std::string& filename, 
//                                    const std::string& title, const std::string& xAxisLabel, 
//                                    const std::string& yAxisLabel) {
//     std::string destination = "../res/" + filename;

//     // Create a canvas
//     TCanvas* canvas = new TCanvas("canvas", "Model Fit Plot", 800, 600);

//     // Filter out points where x == 0 to avoid the singularity
//     std::vector<double> filtered_x, filtered_y;
//     for (size_t i = 0; i < v_x.size(); ++i) {
//         if (v_x[i] != 0) {  // Only include points where x != 0
//             filtered_x.push_back(v_x[i]);
//             filtered_y.push_back(v_y[i]);
//         }
//     }

//     // Create a TGraph from the filtered data
//     TGraph* graph = new TGraph(filtered_x.size(), &filtered_x[0], &filtered_y[0]);

//     // Define the fit function with the specified model
//     TF1* modelFit = new TF1("modelFit", "[0] * pow((sin([1] * x) / ([1] * x)) * (sin([2] * [3] * x) / sin([3] * x)), 2)", filtered_x.front(), filtered_x.back());

//     // Set initial parameters, using defaults if not enough values are provided
//     double param0 = (initialParams.size() > 0) ? initialParams[0] : 3500;  // Default for a
//     double param1 = (initialParams.size() > 1) ? initialParams[1] : 0.3;   // Default for b
//     double param2 = (initialParams.size() > 2) ? initialParams[2] : 1.0;   // Default for m
//     double param3 = (initialParams.size() > 3) ? initialParams[3] : 0.5;   // Default for d

//     modelFit->SetParameter(0, param0);
//     modelFit->SetParameter(1, param1);
//     modelFit->SetParameter(2, param2);
//     modelFit->SetParameter(3, param3);

//     // Perform the fit with status check
//     int fitStatus = graph->Fit("modelFit", "SR");
//     if (fitStatus != 0) {
//         std::cerr << "Fit failed for: " << filename << " (status: " << fitStatus << ")" << std::endl;
//         delete canvas;
//         delete graph;
//         delete modelFit;
//         return;
//     }

//     modelFit->SetNpx(100000);

//     // Set titles and labels
//     graph->SetTitle(title.c_str());
//     graph->GetXaxis()->SetTitle(xAxisLabel.c_str());
//     graph->GetYaxis()->SetTitle(yAxisLabel.c_str());
//     graph->GetXaxis()->CenterTitle();
//     graph->GetYaxis()->CenterTitle();
//     gPad->SetGrid();

//     // Adjust the axes dynamically based on filtered data
//     double xMin = *std::min_element(filtered_x.begin(), filtered_x.end());
//     double xMax = *std::max_element(filtered_x.begin(), filtered_x.end());
//     double yMin = *std::min_element(filtered_y.begin(), filtered_y.end());
//     double yMax = *std::max_element(filtered_y.begin(), filtered_y.end());

//     // Add padding to Y-axis
//     graph->GetXaxis()->SetLimits(xMin, xMax);
//     graph->GetYaxis()->SetRangeUser(yMin * 0.95, yMax * 1.05);

//     // Draw the graph with points
//     graph->SetMarkerStyle(20);
//     graph->SetMarkerSize(1.0);
//     graph->Draw("AP");

//     // Plot the fit function as a smooth line
//     modelFit->SetLineColor(kRed);
//     modelFit->SetLineWidth(2);
//     modelFit->Draw("SAME");

//     // Update the canvas and save
//     canvas->Update();
//     canvas->SaveAs(destination.c_str());
//     canvas->Draw();

//     // Print the fit results
//     std::cout << "Model Fit Results for: " << filename << std::endl;
//     std::cout << "Parameter a: " << modelFit->GetParameter(0) << " ± " << modelFit->GetParError(0) << std::endl;
//     std::cout << "Parameter b: " << modelFit->GetParameter(1) << " ± " << modelFit->GetParError(1) << std::endl;
//     std::cout << "Parameter m: " << modelFit->GetParameter(2) << " ± " << modelFit->GetParError(2) << std::endl;
//     std::cout << "Parameter d: " << modelFit->GetParameter(3) << " ± " << modelFit->GetParError(3) << std::endl;
//     std::cout << "====================" << std::endl;

//     // Clean up
//     delete canvas;
//     delete graph;
//     delete modelFit;
// }

//USED TO GET DECENT APROXIMATED FITS
void PlotsMaker::MakeModelFitPlot(const std::vector<double>& v_x, const std::vector<double>& v_y,
                                   const std::vector<double>& initialParams, 
                                   const std::string& filename, 
                                   const std::string& title, const std::string& xAxisLabel, 
                                   const std::string& yAxisLabel) {
    std::string destination = "../res/" + filename;

    // Create a canvas
    TCanvas* canvas = new TCanvas("canvas", "Model Fit Plot", 800, 600);

    // Filter out points where x == 0 to avoid the singularity
    std::vector<double> filtered_x, filtered_y;
    for (size_t i = 0; i < v_x.size(); ++i) {
        if (v_x[i] != 0) {  // Only include points where x != 0
            filtered_x.push_back(v_x[i]);
            filtered_y.push_back(v_y[i]);
        }
    }

    // Create a TGraph from the filtered data
    TGraph* graph = new TGraph(filtered_x.size(), &filtered_x[0], &filtered_y[0]);

    // Define the fit function with the specified model
    TF1* modelFit = new TF1("modelFit", "[0] * pow(sin([1] * (x-[2])) / abs([1] * (x-[2])), 2)", filtered_x.front(), filtered_x.back());

    // Set initial parameters, using defaults if not enough values are provided
    double param0 = (initialParams.size() > 0) ? initialParams[0] : 3500;  // Default for E
    double param1 = (initialParams.size() > 1) ? initialParams[1] : 0.3;   // Default for a
    double param2 = (initialParams.size() > 2) ? initialParams[2] : 0.0;   // Default offset

    modelFit->SetParameter(0, param0);
    modelFit->SetParameter(1, param1);
    modelFit->SetParameter(2, param2);

    // Perform the fit with status check
    int fitStatus = graph->Fit("modelFit", "SR");
    if (fitStatus != 0) {
        std::cerr << "Fit failed for: " << filename << " (status: " << fitStatus << ")" << std::endl;
        delete canvas;
        delete graph;
        delete modelFit;
        return;
    }

    modelFit->SetNpx(10000);

    // Set titles and labels
    graph->SetTitle(title.c_str());
    graph->GetXaxis()->SetTitle(xAxisLabel.c_str());
    graph->GetYaxis()->SetTitle(yAxisLabel.c_str());
    graph->GetXaxis()->CenterTitle();
    graph->GetYaxis()->CenterTitle();
    gPad->SetGrid();

    // Adjust the axes dynamically based on filtered data
    double xMin = *std::min_element(filtered_x.begin(), filtered_x.end());
    double xMax = *std::max_element(filtered_x.begin(), filtered_x.end());
    double yMin = *std::min_element(filtered_y.begin(), filtered_y.end());
    double yMax = *std::max_element(filtered_y.begin(), filtered_y.end());

    // Add padding to Y-axis
    graph->GetXaxis()->SetLimits(xMin, xMax);
    graph->GetYaxis()->SetRangeUser(yMin * 0.95, yMax * 1.05);

    // Draw the graph with points
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.0);
    graph->Draw("AP");

    // Plot the fit function as a smooth line
    modelFit->SetLineColor(kRed);
    modelFit->SetLineWidth(2);
    modelFit->Draw("SAME");

    // Update the canvas and save
    canvas->Update();
    canvas->SaveAs(destination.c_str());
    canvas->Draw();

    // Print the fit results
    std::cout << "Model Fit Results for: " << filename << std::endl;
    std::cout << "Parameter E: " << modelFit->GetParameter(0) << " ± " << modelFit->GetParError(0) << std::endl;
    std::cout << "Parameter a: " << modelFit->GetParameter(1) << " ± " << modelFit->GetParError(1) << std::endl;
    std::cout << "Parameter offset: " << modelFit->GetParameter(2) << " ± " << modelFit->GetParError(2) << std::endl;
    std::cout << "====================" << std::endl;

    // Clean up
    delete canvas;
    delete graph;
    delete modelFit;
}

void PlotsMaker::MakeFitPlot(const std::vector<double>& v_x, const std::vector<double>& v_y,
                             const std::string& filename, 
                             const std::string& title, const std::string& xAxisLabel, 
                             const std::string& yAxisLabel) {
    std::string destination = "../res/" + filename;

    // Create a canvas
    TCanvas* canvas = new TCanvas("canvas", title.c_str(), 800, 600);

    // Create a TGraph from the filtered v_x and v_y data
    TGraph* graph = new TGraph(v_x.size(), &v_x[0], &v_y[0]);

    // Define the fit function with the specified model
    TF1* modelFit = new TF1("modelFit", "[0] * pow(sin([1] * x) / ([1] * x), 2) + [2]", v_x.front(), v_x.back());

    // Set initial parameters
    double initialE = 0.5 * (*std::max_element(v_y.begin(), v_y.end())); // Example initial guess for E
    double initialA = 0.23; // Example initial guess for a
    double initialOffset = 0.0; // Example initial guess for offset

    // Set parameters
    modelFit->SetParameter(0, initialE); // Initial guess for E
    modelFit->SetParameter(1, initialA); // Initial guess for a
    modelFit->SetParameter(2, initialOffset); // Initial guess for offset

    // Set limits for parameters
    modelFit->SetParLimits(0, initialE - 200, initialE + 200); // Limits for parameter E
    modelFit->SetParLimits(1, 0.0, 4.0); // Limits for parameter a
    modelFit->SetParLimits(2, -100, 100); // Limits for offset

    // Increase the number of points used to draw the fit function curve for smoother plot
    modelFit->SetNpx(1000);  // Set the number of points for smoother fit curve

    // Perform the fit
    int fitStatus = graph->Fit("modelFit", "S");

    // Check for fitting issues
    if (fitStatus != 0) {
        std::cerr << "Fit failed with status: " << fitStatus << std::endl;
        // Handle fitting failure if needed
    }

    // Set titles and labels
    graph->SetTitle(title.c_str());
    graph->GetXaxis()->SetTitle(xAxisLabel.c_str());
    graph->GetYaxis()->SetTitle(yAxisLabel.c_str());

    // Center the axis labels
    graph->GetXaxis()->CenterTitle();
    graph->GetYaxis()->CenterTitle();

    // Adjust the axes dynamically based on filtered data
    double xMin = *std::min_element(v_x.begin(), v_x.end());
    double xMax = *std::max_element(v_x.begin(), v_x.end());
    double yMin = *std::min_element(v_y.begin(), v_y.end());
    double yMaxFinal = modelFit->GetParameter(0); // Use the best E value as yMax

    // Add some padding to the axis limits
    graph->GetXaxis()->SetLimits(xMin, xMax);
    graph->GetYaxis()->SetRangeUser(yMin * 0.95, yMaxFinal * 1.05); // Adjusting Y-axis with padding

    // Set grid
    gPad->SetGrid();

    // Draw the graph with points
    graph->SetMarkerStyle(20);  // Set marker style
    graph->SetMarkerSize(0.5);  // Set marker size
    graph->Draw("AP");  // Draw the graph with axes and points

    // Draw the best fit curve
    modelFit->SetLineColor(kRed);
    modelFit->SetLineWidth(1.5);  // Set the line thickness to 1.5
    modelFit->Draw("SAME");  // Ensure to draw the fit on the same canvas

    // Force an update to the canvas to reflect the drawn objects
    canvas->Update();

    // Save the plot as a PDF
    canvas->SaveAs(destination.c_str());

    // Display the canvas
    canvas->Draw();

    // Print the fit results
    std::cout << "Fit Results for: " << filename << std::endl;
    std::cout << "Fitted E: " << modelFit->GetParameter(0) << std::endl;
    std::cout << "Fitted a: " << modelFit->GetParameter(1) << std::endl;
    std::cout << "Fitted offset: " << modelFit->GetParameter(2) << std::endl;
    std::cout << "Chi2/NDF: " << modelFit->GetChisquare() / modelFit->GetNDF() << std::endl;
    std::cout << "====================" << std::endl;

    // Clean up
    delete canvas;
    delete graph;
    delete modelFit;  // Clean up the fit model
}

void PlotsMaker::OptimizeAndPlot(const std::vector<double>& v_x, const std::vector<double>& v_y,
                                   const std::string& funcStr, const std::vector<double>& initialParams,
                                   const std::string& filename,
                                   const std::string& title, const std::string& xAxisLabel, 
                                   const std::string& yAxisLabel, int numPoints) {
    std::string destination = "../res/" + filename;

    // Initialize best parameters and best chi-squared
    std::vector<double> bestParams = initialParams;
    double bestChi2 = std::numeric_limits<double>::max();

    // Define ranges around the initial parameters
    std::vector<double> parameterRanges = {0.05, 0.05, 0.05}; // Define step size for each parameter
    
    // Iterate over a small range of values around the initial parameters
    for (double e_offset = -1; e_offset <= 1; e_offset += parameterRanges[0]) {
        for (double a_offset = -0.1; a_offset <= 0.1; a_offset += parameterRanges[1]) {
            for (double offset_offset = -10; offset_offset <= 10; offset_offset += parameterRanges[2]) {
                // Set current parameters
                std::vector<double> currentParams = {
                    initialParams[0] + e_offset,
                    initialParams[1] + a_offset,
                    initialParams[2] + offset_offset
                };

                // Create a TGraph from the v_x and v_y data
                TGraph* graph = new TGraph(v_x.size(), &v_x[0], &v_y[0]);

                // Define the function based on the user-defined string
                TF1* userFunc = new TF1("userFunc", funcStr.c_str(), v_x.front(), v_x.back());

                // Set the parameters for the function
                for (size_t i = 0; i < currentParams.size(); ++i) {
                    userFunc->SetParameter(i, currentParams[i]);
                }

                // Perform a fit
                int fitStatus = graph->Fit("userFunc", "S");

                // Check for fitting success
                if (fitStatus == 0) {
                    double chi2 = userFunc->GetChisquare() / userFunc->GetNDF();
                    // If this fit is better than the previous best
                    if (chi2 < bestChi2) {
                        bestChi2 = chi2;
                        bestParams = currentParams;  // Update best parameters
                    }
                }

                // Clean up
                delete userFunc;
                delete graph;
            }
        }
    }

    // Create a final plot with the best parameters
    TCanvas* canvas = new TCanvas("canvas", title.c_str(), 800, 600);
    TGraph* finalGraph = new TGraph(v_x.size(), &v_x[0], &v_y[0]);
    finalGraph->SetMarkerStyle(20);
    finalGraph->SetMarkerSize(1.0);
    finalGraph->Draw("AP");  // Draw data points

    // Define the function with the best parameters
    TF1* finalFunc = new TF1("finalFunc", funcStr.c_str(), v_x.front(), v_x.back());
    for (size_t i = 0; i < bestParams.size(); ++i) {
        finalFunc->SetParameter(i, bestParams[i]);
    }

    // Generate fit points
    double xFit[numPoints];
    double yFit[numPoints];
    for (int i = 0; i < numPoints; ++i) {
        xFit[i] = v_x.front() + i * (v_x.back() - v_x.front()) / (numPoints - 1);
        yFit[i] = finalFunc->Eval(xFit[i]);
    }

    // Create a TGraph for the fitted model
    TGraph* fitGraph = new TGraph(numPoints, xFit, yFit);
    fitGraph->SetLineColor(kRed);
    fitGraph->SetLineWidth(2);
    fitGraph->Draw("SAME");

    // Set titles and labels
    finalGraph->SetTitle(title.c_str());
    finalGraph->GetXaxis()->SetTitle(xAxisLabel.c_str());
    finalGraph->GetYaxis()->SetTitle(yAxisLabel.c_str());
    gPad->SetGrid();

    // Update the canvas
    canvas->Update();
    canvas->SaveAs(destination.c_str());
    canvas->Draw();

    // Print best parameters
    std::cout << "Best Fit Parameters for: " << filename << std::endl;
    for (size_t i = 0; i < bestParams.size(); ++i) {
        std::cout << "Parameter [" << i << "]: " << bestParams[i] << std::endl;
    }
    std::cout << "Chi2/NDF: " << bestChi2 << std::endl;
    std::cout << "====================" << std::endl;

    // Clean up
    delete canvas;
    delete finalGraph;
    delete finalFunc;
    delete fitGraph;
}

void PlotsMaker::PlotFunctionWithData(const std::vector<double>& v_x, const std::vector<double>& v_y,
                                       const std::string& funcStr, const std::vector<double>& params,
                                       const std::string& filename,
                                       const std::string& title, const std::string& xAxisLabel, 
                                       const std::string& yAxisLabel, int numPoints) {
    std::string destination = "../res/" + filename;

    // Create a canvas
    TCanvas* canvas = new TCanvas("canvas", title.c_str(), 800, 600);
    
    // Create a TGraph from the v_x and v_y data
    TGraph* graph = new TGraph(v_x.size(), &v_x[0], &v_y[0]);

    // Define the function based on the user-defined string
    TF1* userFunc = new TF1("userFunc", funcStr.c_str(), v_x.front(), v_x.back());

    // Set the parameters for the function
    for (size_t i = 0; i < params.size(); ++i) {
        userFunc->SetParameter(i, params[i]);
    }

    // Set titles and labels
    graph->SetTitle(title.c_str());
    graph->GetXaxis()->SetTitle(xAxisLabel.c_str());
    graph->GetYaxis()->SetTitle(yAxisLabel.c_str());

    // Center the axis labels
    graph->GetXaxis()->CenterTitle();
    graph->GetYaxis()->CenterTitle();

    // Set grid
    gPad->SetGrid();

    // Draw the graph with points
    graph->SetMarkerStyle(20);  // Set marker style
    graph->SetMarkerSize(1.0);  // Set marker size
    graph->Draw("AP");  // Draw the graph with axes and points

    // Create an array for the points for the user-defined function
    double xFit[numPoints];
    double yFit[numPoints];

    // Check if v_x contains only one point
    if (v_x.size() == 1) {
        // If there is only one point, generate points around it
        double centerX = v_x[0];
        double centerY = v_y[0];
        double spacing = 0.000001;  // Set spacing around the center point

        // Fill xFit with points around the centerX
        for (int i = 0; i < numPoints; ++i) {
            xFit[i] = centerX + (i - numPoints / 2) * spacing;  // Centering around the single point
            yFit[i] = userFunc->Eval(xFit[i]);  // Evaluate the user-defined function
        }

        // Set x and y axis ranges to center the point
        graph->GetXaxis()->SetRangeUser(centerX - 0.001, centerX + 0.001);  // Adjust range to center the point
        graph->GetYaxis()->SetRangeUser(centerY - 0.001, centerY + 0.001);  // Adjust range to center the point
    } else {
        // If more than one point, use the existing range
        for (int i = 0; i < numPoints; ++i) {
            xFit[i] = v_x.front() + i * (v_x.back() - v_x.front()) / (numPoints - 1);
            yFit[i] = userFunc->Eval(xFit[i]);  // Evaluate the user-defined function
        }
        
        // Set the axis limits to the min and max of the data
        graph->GetXaxis()->SetRangeUser(0, v_x.back());
        graph->GetYaxis()->SetRangeUser(*std::min_element(v_y.begin(), v_y.end()), *std::max_element(v_y.begin(), v_y.end()));
    }

    // Create a TGraph for the fitted model
    TGraph* fitGraph = new TGraph(numPoints, xFit, yFit);
    fitGraph->SetLineColor(kRed);  // Set the fit curve color
    fitGraph->SetLineWidth(2);      // Set the fit curve width

    // Draw the fitted model on the same canvas
    fitGraph->Draw("SAME");  // Draw the fit graph on top of the data

    // Update the canvas
    canvas->Update();

    // Save the plot as a PDF
    canvas->SaveAs(destination.c_str());

    // Display the canvas
    canvas->Draw();

    // Print function parameters
    std::cout << "User Function Parameters for: " << filename << std::endl;
    for (size_t i = 0; i < params.size(); ++i) {
        std::cout << "Parameter [" << i << "]: " << userFunc->GetParameter(i) 
                  << " ± " << userFunc->GetParError(i) << std::endl;
    }
    std::cout << "====================" << std::endl;

    // Clean up
    delete canvas;
    delete graph;
    delete userFunc;  // Clean up the user-defined function
    delete fitGraph;  // Clean up the fit graph
}

void PlotsMaker::PlotFunction(const std::string& funcStr, const std::vector<double>& parameters, 
                              double xMin, double xMax, 
                              const std::string& filename, 
                              const std::string& title, const std::string& xAxisLabel, 
                              const std::string& yAxisLabel, int nPoints = 1000) {
    // Create a canvas for plotting
    TCanvas* canvas = new TCanvas("canvas", "Function Plot", 800, 600);

    // Define the TF1 using the provided function string
    TF1* func = new TF1("func", funcStr.c_str(), xMin, xMax);

    // Set the number of points to plot (increase for more precision)
    func->SetNpx(nPoints);

    // Check if the number of parameters matches
    if (parameters.size() != func->GetNpar()) {
        std::cerr << "Mismatch between the number of parameters in the function and the provided parameters." << std::endl;
        delete canvas;
        delete func;
        return;
    }

    // Set the parameters for the function
    for (size_t i = 0; i < parameters.size(); ++i) {
        func->SetParameter(i, parameters[i]);
    }

    // Set titles and labels
    func->SetTitle(title.c_str());
    func->GetXaxis()->SetTitle(xAxisLabel.c_str());
    func->GetYaxis()->SetTitle(yAxisLabel.c_str());

    // Center the axis labels
    func->GetXaxis()->CenterTitle();
    func->GetYaxis()->CenterTitle();

    // Set grid
    gPad->SetGrid();

    // Draw the function
    func->SetLineColor(kBlue);
    func->Draw();

    // Update the canvas
    canvas->Update();

    // Save the plot as a PDF
    std::string destination = "../res/" + filename;
    canvas->SaveAs(destination.c_str());

    // Display the canvas
    canvas->Draw();

    // Clean up
    delete canvas;
    delete func;
}

void PlotsMaker::MakeLinearPlot(const std::vector<double>& x, const std::vector<double>& error_x, 
                                const std::vector<double>& y, const std::vector<double>& error_y,
                                const std::string& filename, const std::string& title,
                                const std::string& xAxisLabel, const std::string& yAxisLabel,
                                double slope, double intercept) {
    std::string destination = "../res/" + filename;

    // Create a canvas
    TCanvas* canvas = new TCanvas("canvas", "Linear Plot", 800, 600);

    // Create TGraphErrors object with the given x, y, error_x, and error_y arrays
    TGraphErrors* graph = new TGraphErrors(x.size(), &x[0], &y[0], &error_x[0], &error_y[0]);

    // Ensure there are at least 4 points
    if (x.size() < 4) {
        std::cerr << "Error: Not enough data points to start the line from the 4th point." << std::endl;
        return;
    }

    // Create the linear function with the given slope and intercept
    // The range of the function is set to start from the x value of the 4th point
    TF1* linearFunction = new TF1("linearFunction", "[0] + [1]*x", x[3], x.back());
    linearFunction->SetParameter(0, intercept);
    linearFunction->SetParameter(1, slope);

    // Set titles and labels
    graph->SetTitle(title.c_str());
    graph->GetXaxis()->SetTitle(xAxisLabel.c_str());
    graph->GetYaxis()->SetTitle(yAxisLabel.c_str());

    // Center the axis labels
    graph->GetXaxis()->CenterTitle();
    graph->GetYaxis()->CenterTitle();

    // Set title offsets
    graph->GetXaxis()->SetTitleOffset(1.2);
    graph->GetYaxis()->SetTitleOffset(1.2);

    // Set grid
    gPad->SetGrid();

    // Draw the graph with points and error bars
    graph->SetMarkerStyle(20); // Set marker style
    graph->SetMarkerSize(1.0); // Set marker size
    graph->Draw("AP"); // Draw the graph with axes, points, and error bars

    // Draw the linear function
    linearFunction->Draw("SAME");

    // Update the canvas
    canvas->Update();

    // Save the plot as a PDF
    canvas->SaveAs(destination.c_str());

    // Display the canvas
    canvas->Draw();

    // Print the provided parameters
    std::cout << "Linear Plot Results for: " << filename << std::endl;
    std::cout << "Slope: " << slope << std::endl;
    std::cout << "Intercept: " << intercept << std::endl;
    std::cout << "====================" << std::endl;

    // Clean up
    delete canvas;
    delete graph;
    delete linearFunction;
}

void PlotsMaker::MakeOverlayHistograms(std::vector<double> info1, std::vector<double> info2, std::string filename,
                                        std::string title, std::string xAxisLabel, std::string yAxisLabel) {
    // Check if the sizes of the vectors match
    if (info1.size() != info2.size()) {
        std::cerr << "Error: Vector sizes do not match!" << std::endl;
        return;
    }

    std::string destination = "../res/" + filename;

    // Create a ROOT canvas
    TCanvas *c1 = new TCanvas("c1", "Overlay Histograms", 800, 600);

    // Create the first histogram
    TH1F *hist1 = new TH1F("hist1", (filename + "_1").c_str(), info1.size(), 0, info1.size());
    // Fill the first histogram with the data from info1
    for (int i = 0; i < info1.size(); ++i) {
        hist1->SetBinContent(i+1, info1[i]);
    }
    hist1->SetFillColorAlpha(kBlue, 0.5); // Set the fill color with some opacity

    // Draw the first histogram
    hist1->Draw("HIST");

    // Create the second histogram
    TH1F *hist2 = new TH1F("hist2", (filename + "_2").c_str(), info2.size(), 0, info2.size());
    // Fill the second histogram with the data from info2
    for (int i = 0; i < info2.size(); ++i) {
        hist2->SetBinContent(i+1, info2[i]);
    }
    hist2->SetFillColor(kRed); // Set the fill color without opacity

    // Draw the second histogram on top of the first one
    hist2->Draw("HIST SAME");

    hist2->SetTitle(title.c_str());
    hist2->GetXaxis()->SetTitle(xAxisLabel.c_str());
    hist2->GetYaxis()->SetTitle(yAxisLabel.c_str());

    gPad->SetGrid();

    // Update the canvas to display the histograms
    c1->Update();
    c1->SaveAs(destination.c_str());

    // Wait for user input to exit the application
    c1->WaitPrimitive();

    // Clean up
    delete hist1;
    delete hist2;
    delete c1;
}

void PlotsMaker::MakeTwoGraphsPlot(const std::vector<double>& x1, const std::vector<double>& y1, 
                                   const std::vector<double>& x2, const std::vector<double>& y2, 
                                   const std::string& filename, std::string title, 
                                   std::string xAxisLabel, std::string yAxisLabel,
                                   std::string legend1, std::string legend2) {
    std::string destination = "../res/" + filename;
    
    // Create a ROOT canvas
    TCanvas *c1 = new TCanvas("c1", "Two Graphs Plot", 1000, 800);
    
    // Create graphs
    TGraph *graph1 = new TGraph(x1.size(), &x1[0], &y1[0]);
    TGraph *graph2 = new TGraph(x2.size(), &x2[0], &y2[0]);
    
    // Set colors and marker styles for the graphs
    graph1->SetLineColor(kBlue);
    graph1->SetMarkerColor(kBlue);
    graph1->SetMarkerStyle(8); // Use circle markers
    graph1->SetMarkerSize(0.5); // Set marker size to 0.5
    
    graph2->SetLineColor(kRed);
    graph2->SetMarkerColor(kRed);
    graph2->SetMarkerStyle(8); // Use circle markers
    graph2->SetMarkerSize(0.5); // Set marker size to 0.5
    
    // Set titles for the graphs
    graph1->SetTitle(title.c_str());
    graph2->SetTitle(title.c_str());
    
    // Set labels for the x and y axes
    graph1->GetXaxis()->SetTitle(xAxisLabel.c_str());
    graph1->GetYaxis()->SetTitle(yAxisLabel.c_str());
    graph2->GetXaxis()->SetTitle(xAxisLabel.c_str());
    graph2->GetYaxis()->SetTitle(yAxisLabel.c_str());

    // Set offsets for the axis titles
    graph1->GetXaxis()->SetTitleOffset(0);
    graph1->GetYaxis()->SetTitleOffset(0);
    graph2->GetXaxis()->SetTitleOffset(0);
    graph2->GetYaxis()->SetTitleOffset(0);
    
    // Center the axis labels
    graph1->GetXaxis()->CenterTitle();
    graph1->GetYaxis()->CenterTitle();
    graph2->GetXaxis()->CenterTitle();
    graph2->GetYaxis()->CenterTitle();

    // Adjust the position of the y-axis label
    graph1->GetYaxis()->SetTitleOffset(1.5);
    graph2->GetYaxis()->SetTitleOffset(1.5);

    gPad->SetGrid();
    
    // Draw the first graph
    graph1->Draw("AP"); // "AP" option to draw with axes and points
    
    // Draw the second graph on top of the first one
    graph2->Draw("P SAME"); // "P SAME" option to overlay the second graph on top of the first one
    
    // Add a legend
    TLegend *legend = new TLegend(0.4, 0.7, 0.9, 0.9);
    legend->SetTextSize(0.03); // Adjust text size
    legend->AddEntry(graph1, legend1.c_str(), "lp");
    legend->AddEntry(graph2, legend2.c_str(), "lp");
    legend->Draw();
    
    // Update the canvas to display the graphs and legend
    c1->Update();
    
    // Save the canvas as an image file
    c1->SaveAs(destination.c_str());
    
    // Wait for user input to exit the application
    c1->WaitPrimitive();
    
    // Clean up
    delete graph1;
    delete graph2;
    delete legend;
    delete c1;
}

void PlotsMaker::MakeFiveGraphsPlot(const std::vector<double>& x1, const std::vector<double>& y1,
                                    const std::vector<double>& x2, const std::vector<double>& y2,
                                    const std::vector<double>& x3, const std::vector<double>& y3,
                                    const std::vector<double>& x4, const std::vector<double>& y4,
                                    const std::vector<double>& x5, const std::vector<double>& y5,
                                    const std::string& filename, std::string title,
                                    std::string xAxisLabel, std::string yAxisLabel,
                                    std::string legend1, std::string legend2,
                                    std::string legend3, std::string legend4,
                                    std::string legend5) {
    std::string destination = "../res/" + filename;
    
    // Create a ROOT canvas
    TCanvas *c1 = new TCanvas("c1", "Five Graphs Plot", 1000, 800);
    
    // Create graphs
    TGraph *graph1 = new TGraph(x1.size(), &x1[0], &y1[0]);
    TGraph *graph2 = new TGraph(x2.size(), &x2[0], &y2[0]);
    TGraph *graph3 = new TGraph(x3.size(), &x3[0], &y3[0]);
    TGraph *graph4 = new TGraph(x4.size(), &x4[0], &y4[0]);
    TGraph *graph5 = new TGraph(x5.size(), &x5[0], &y5[0]);
    
    // Set colors and marker styles for the graphs
    graph1->SetLineColor(kBlue);
    graph1->SetMarkerColor(kBlue);
    graph1->SetMarkerStyle(8); // Use circle markers
    graph1->SetMarkerSize(0.7); // Set marker size to 0.5
    
    graph2->SetLineColor(kRed);
    graph2->SetMarkerColor(kRed);
    graph2->SetMarkerStyle(8); // Use circle markers
    graph2->SetMarkerSize(0.7); // Set marker size to 0.5
    
    graph3->SetLineColor(kGreen);
    graph3->SetMarkerColor(kGreen);
    graph3->SetMarkerStyle(8); // Use circle markers
    graph3->SetMarkerSize(0.7); // Set marker size to 0.5
    
    graph4->SetLineColor(kMagenta);
    graph4->SetMarkerColor(kMagenta);
    graph4->SetMarkerStyle(8); // Use circle markers
    graph4->SetMarkerSize(0.7); // Set marker size to 0.5
    
    graph5->SetLineColor(kBlack);
    graph5->SetMarkerColor(kBlack);
    graph5->SetMarkerStyle(8); // Use circle markers
    graph5->SetMarkerSize(0.7); // Set marker size to 0.5
    
    // Set titles for the graphs
    graph1->SetTitle(title.c_str());
    graph2->SetTitle(title.c_str());
    graph3->SetTitle(title.c_str());
    graph4->SetTitle(title.c_str());
    graph5->SetTitle(title.c_str());
    
    // Set labels for the x and y axes
    graph1->GetXaxis()->SetTitle(xAxisLabel.c_str());
    graph1->GetYaxis()->SetTitle(yAxisLabel.c_str());
    graph2->GetXaxis()->SetTitle(xAxisLabel.c_str());
    graph2->GetYaxis()->SetTitle(yAxisLabel.c_str());
    graph3->GetXaxis()->SetTitle(xAxisLabel.c_str());
    graph3->GetYaxis()->SetTitle(yAxisLabel.c_str());
    graph4->GetXaxis()->SetTitle(xAxisLabel.c_str());
    graph4->GetYaxis()->SetTitle(yAxisLabel.c_str());
    graph5->GetXaxis()->SetTitle(xAxisLabel.c_str());
    graph5->GetYaxis()->SetTitle(yAxisLabel.c_str());

    // Set offsets for the axis titles
    graph1->GetXaxis()->SetTitleOffset(1.2);
    graph1->GetYaxis()->SetTitleOffset(1.2);
    graph2->GetXaxis()->SetTitleOffset(1.2);
    graph2->GetYaxis()->SetTitleOffset(1.2);
    graph3->GetXaxis()->SetTitleOffset(1.2);
    graph3->GetYaxis()->SetTitleOffset(1.2);
    graph4->GetXaxis()->SetTitleOffset(1.2);
    graph4->GetYaxis()->SetTitleOffset(1.2);
    graph5->GetXaxis()->SetTitleOffset(1.2);
    graph5->GetYaxis()->SetTitleOffset(1.2);
    
    // Center the axis labels
    graph1->GetXaxis()->CenterTitle();
    graph1->GetYaxis()->CenterTitle();
    graph2->GetXaxis()->CenterTitle();
    graph2->GetYaxis()->CenterTitle();
    graph3->GetXaxis()->CenterTitle();
    graph3->GetYaxis()->CenterTitle();
    graph4->GetXaxis()->CenterTitle();
    graph4->GetYaxis()->CenterTitle();
    graph5->GetXaxis()->CenterTitle();
    graph5->GetYaxis()->CenterTitle();

    // Adjust the position of the y-axis label
    graph1->GetYaxis()->SetTitleOffset(1.5);
    graph2->GetYaxis()->SetTitleOffset(1.5);
    graph3->GetYaxis()->SetTitleOffset(1.5);
    graph4->GetYaxis()->SetTitleOffset(1.5);
    graph5->GetYaxis()->SetTitleOffset(1.5);

    gPad->SetGrid();
    
    // Draw the graphs
    graph1->Draw("AP"); // "AP" option to draw with axes and points
    graph2->Draw("P SAME"); // Overlay the second graph on top of the first one
    graph3->Draw("P SAME"); // Overlay the third graph
    graph4->Draw("P SAME"); // Overlay the fourth graph
    graph5->Draw("P SAME"); // Overlay the fifth graph
    
    // Add a legend
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->SetTextSize(0.02); // Adjust text size
    legend->AddEntry(graph1, legend1.c_str(), "lp");
    legend->AddEntry(graph2, legend2.c_str(), "lp");
    legend->AddEntry(graph3, legend3.c_str(), "lp");
    legend->AddEntry(graph4, legend4.c_str(), "lp");
    legend->AddEntry(graph5, legend5.c_str(), "lp");
    legend->Draw();
    
    // Update the canvas to display the graphs and legend
    c1->Update();
    
    // Save the canvas as an image file
    c1->SaveAs(destination.c_str());
    
    // Wait for user input to exit the application
    c1->WaitPrimitive();
    
    // Clean up
    delete graph1;
    delete graph2;
    delete graph3;
    delete graph4;
    delete graph5;
    delete legend;
    delete c1;
}

void PlotsMaker::MakeMultipleGraphsPlot(const std::vector<std::vector<double>>& xVectors, 
                                        const std::vector<std::vector<double>>& yVectors, 
                                        const std::string& filename, 
                                        const std::string& title, 
                                        const std::string& xAxisLabel, 
                                        const std::string& yAxisLabel, 
                                        const std::vector<std::string>& legends) {
    // Check if the input vectors are valid
    if (xVectors.size() != yVectors.size() || xVectors.empty() || legends.size() != xVectors.size()) {
        std::cerr << "Error: Mismatch in the number of x, y vectors, or legends." << std::endl;
        return;
    }

    std::string destination = "../res/" + filename;

    // Create a ROOT canvas
    TCanvas *c1 = new TCanvas("c1", "Multiple Graphs Plot", 1000, 800);

    // Create a legend to display graph names
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->SetTextSize(0.03);

    // Enable grid
    gPad->SetGrid();

    // List of colors to use for each graph
    std::vector<int> colors = {kBlue, kRed, kGreen + 2, kMagenta, kCyan, kOrange + 7};
    int colorIndex = 0;

    // Loop over the vectors of vectors to create and plot each graph
    std::vector<TGraph*> graphs;
    for (size_t i = 0; i < xVectors.size(); ++i) {
        if (xVectors[i].size() != yVectors[i].size()) {
            std::cerr << "Error: Mismatch in the size of x and y vectors for graph " << i + 1 << std::endl;
            continue;
        }

        // Create a new TGraph for each dataset
        TGraph *graph = new TGraph(xVectors[i].size(), &xVectors[i][0], &yVectors[i][0]);
        graph->SetLineColor(colors[colorIndex % colors.size()]);
        graph->SetMarkerColor(colors[colorIndex % colors.size()]);
        graph->SetMarkerStyle(8);  // Use circle markers
        graph->SetMarkerSize(0.5);
        graph->SetLineWidth(2);    // Set line width for better visibility of connecting lines

        // Set graph titles and axis labels only for the first graph
        if (i == 0) {
            graph->SetTitle(title.c_str());
            graph->GetXaxis()->SetTitle(xAxisLabel.c_str());
            graph->GetYaxis()->SetTitle(yAxisLabel.c_str());

            // Set axis label offsets and center titles
            graph->GetXaxis()->SetTitleOffset(1.2);
            graph->GetYaxis()->SetTitleOffset(1.5);
            graph->GetXaxis()->CenterTitle();
            graph->GetYaxis()->CenterTitle();
        }

        // Draw the first graph with axes, points, and lines
        if (i == 0) {
            graph->Draw("APL"); // Axes, points, and lines for the first graph
        } else {
            graph->Draw("PL SAME"); // Points and lines for subsequent graphs
        }

        // Add the graph to the legend
        legend->AddEntry(graph, legends[i].c_str(), "lp");

        // Save the graph to the list for later cleanup
        graphs.push_back(graph);

        // Increment color index for the next graph
        colorIndex++;
    }

    // Draw the legend
    legend->Draw();

    // Update the canvas to display the graphs and legend
    c1->Update();

    // Save the canvas as an image file
    c1->SaveAs(destination.c_str());

    // Wait for user input to exit the application
    c1->WaitPrimitive();

    // Clean up
    for (auto graph : graphs) {
        delete graph;
    }
    delete legend;
    delete c1;
}