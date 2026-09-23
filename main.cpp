#include <iostream>
#include <string>
#include <memory>
#include <vector>
#include <map>
#include <sstream>

// VTK includes
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkKdTreePointLocator.h>

#include "AbaqusInpReader.h"
#include "NpyReader.h"
#include "MappingMethods.h"

class TensorMapper {
public:
    enum class TransferMethod {
        DIRECT,
        PROXIMITY,
        N_POINT,
        GAUSSIAN
    };
    
    struct Config {
        TransferMethod method = TransferMethod::DIRECT;
        double searchRadius = 0.0;
        int numNeighbors = 5;
        double gaussianSigma = 1.0;
        bool applyGlobalSmoothing = false;
        int smoothingIterations = 1;
        bool useKdTree = true;
        bool autoDetectDataLocation = true;
        std::vector<double> origin = {0.0, 0.0, 0.0};
        std::vector<double> spacing = {1.0, 1.0, 1.0};
        std::map<int, std::string> materialMap;
        double relaxationFactor = 0.5;
    };
    
    TensorMapper(const Config& config = Config()) : config_(config) {}
    
    // Load source mesh and data from file
    bool loadSource(const std::string& filename, const std::string& dataName = "") {
        sourceFile_ = filename;
        sourceMesh_ = readMesh(filename, config_.origin, config_.spacing);
        if (!sourceMesh_) {
            std::cerr << "Failed to load source mesh: " << filename << std::endl;
            return false;
        }

        std::cout << "Available point data arrays:" << std::endl;
        for (int i = 0; i < sourceMesh_->GetPointData()->GetNumberOfArrays(); i++) {
            std::cout << "  - " << sourceMesh_->GetPointData()->GetArray(i)->GetName() << std::endl;
        }

        std::cout << "Available cell data arrays:" << std::endl;
        for (int i = 0; i < sourceMesh_->GetCellData()->GetNumberOfArrays(); i++) {
            std::cout << "  - " << sourceMesh_->GetCellData()->GetArray(i)->GetName() << std::endl;
        }
        
        // If data name specified, use it, otherwise use first available array
        if (!dataName.empty()) {
            sourceDataArray_ = findDataArray(sourceMesh_, dataName);
        } else {
            sourceDataArray_ = findFirstDataArray(sourceMesh_);
        }
        
        if (!sourceDataArray_) {
            std::cerr << "Failed to find data array in source mesh" << std::endl;
            return false;
        }
        
        // Detect data location if auto-detection is enabled
        if (config_.autoDetectDataLocation) {
            sourceDataLocation_ = detectDataLocation(sourceDataArray_, sourceMesh_);
            std::cout << "Source data location: " 
                      << (sourceDataLocation_ == DataLocation::NODAL ? "Nodal" : "Elemental") 
                      << std::endl;
        }
        
        // Build KD-tree for source mesh immediately after loading
        if (config_.useKdTree) {
            buildSourceKdTree();
        }
        
        return true;
    }
    
    // Load target mesh from file
    bool loadTarget(const std::string& filename) {
        targetFile_ = filename;
        targetMesh_ = readMesh(filename);
        if (!targetMesh_) {
            std::cerr << "Failed to load target mesh: " << filename << std::endl;
            return false;
        }
        return true;
    }
    
    // Perform mapping and return result
    bool mapData(const std::string& outputName = "", DataLocation targetLocation = DataLocation::NODAL) {
        if (!sourceMesh_ || !targetMesh_ || !sourceDataArray_) {
            std::cerr << "Source or target not properly initialized" << std::endl;
            return false;
        }
        
        // Set target data location
        targetDataLocation_ = targetLocation;
        std::cout << "Target data location: " 
                  << (targetDataLocation_ == DataLocation::NODAL ? "Nodal" : "Elemental") 
                  << std::endl;
        
        // Use the output name provided or use source name
        std::string dataName = outputName.empty() ? 
            sourceDataArray_->GetName() : outputName;
        
        // Create result array with same type and components as source
        vtkIdType numTargetItems = (targetDataLocation_ == DataLocation::NODAL) ?
            targetMesh_->GetNumberOfPoints() : targetMesh_->GetNumberOfCells();
        
        mappedDataArray_ = createSimilarArray(sourceDataArray_, dataName, numTargetItems);
        
        // Choose mapping method based on config
        switch (config_.method) {
            case TransferMethod::DIRECT:
                MappingMethods::mapDirect(sourceDataArray_, mappedDataArray_, targetMesh_, kdTree_, sourceDataLocation_, targetDataLocation_);
                break;
            case TransferMethod::PROXIMITY:
                MappingMethods::mapProximity(sourceDataArray_, mappedDataArray_, targetMesh_, kdTree_, sourceDataLocation_, targetDataLocation_, config_.searchRadius);
                break;
            case TransferMethod::N_POINT:
                MappingMethods::mapNPoint(sourceDataArray_, mappedDataArray_, targetMesh_, kdTree_, sourceDataLocation_, targetDataLocation_, config_.numNeighbors);
                break;
            case TransferMethod::GAUSSIAN:
                MappingMethods::mapGaussian(sourceDataArray_, mappedDataArray_, targetMesh_, kdTree_, sourceDataLocation_, targetDataLocation_, 
                                        config_.searchRadius, config_.gaussianSigma);
                break;
        }
        
        // Apply global smoothing if requested
        if (config_.applyGlobalSmoothing) {
            MappingMethods::applyGlobalSmoothing(mappedDataArray_, targetMesh_, targetDataLocation_, config_.smoothingIterations, config_.relaxationFactor);
        }
        
        // Add result to target mesh in the appropriate location
        if (targetDataLocation_ == DataLocation::NODAL) {
            targetMesh_->GetPointData()->AddArray(mappedDataArray_);
        } else {
            targetMesh_->GetCellData()->AddArray(mappedDataArray_);
        }
        
        return true;
    }
    
    // Save result to file
    bool saveResult(const std::string& filename) {
        if (!targetMesh_ || !mappedDataArray_) {
            std::cerr << "No result to save" << std::endl;
            return false;
        }
        
        return writeMesh(filename, targetMesh_);
    }
    
    // Get mapped result
    vtkSmartPointer<vtkUnstructuredGrid> getResultMesh() {
        return targetMesh_;
    }
    
protected:
    Config config_;
    vtkSmartPointer<vtkUnstructuredGrid> sourceMesh_;
    vtkSmartPointer<vtkUnstructuredGrid> targetMesh_;
    vtkSmartPointer<vtkDataArray> sourceDataArray_;
    vtkSmartPointer<vtkDataArray> mappedDataArray_;
    vtkSmartPointer<vtkKdTreePointLocator> kdTree_; // KD-tree for accelerated search
    DataLocation sourceDataLocation_ = DataLocation::NODAL;
    DataLocation targetDataLocation_ = DataLocation::NODAL;
    std::string sourceFile_;
    std::string targetFile_;
    
    // Detect if data is nodal or elemental
    DataLocation detectDataLocation(vtkSmartPointer<vtkDataArray> dataArray, 
                                  vtkSmartPointer<vtkUnstructuredGrid> mesh) {
        // Check if array is in point data
        auto pointData = mesh->GetPointData();
        for (int i = 0; i < pointData->GetNumberOfArrays(); i++) {
            if (pointData->GetArray(i) == dataArray) {
                return DataLocation::NODAL;
            }
        }
        
        // Check if array is in cell data
        auto cellData = mesh->GetCellData();
        for (int i = 0; i < cellData->GetNumberOfArrays(); i++) {
            if (cellData->GetArray(i) == dataArray) {
                return DataLocation::ELEMENTAL;
            }
        }
        
        // Default to nodal if not found (though this should not happen)
        return DataLocation::NODAL;
    }
    
    // Build KD-tree for source mesh
    void buildSourceKdTree() {
        kdTree_ = vtkSmartPointer<vtkKdTreePointLocator>::New();
        
        if (sourceDataLocation_ == DataLocation::NODAL) {
            // For nodal data, use the mesh points directly
            kdTree_->SetDataSet(sourceMesh_);
        } else {
            // For elemental data, create points at element centroids
            auto centroids = MappingMethods::calculateElementCentroids(sourceMesh_);
            
            // Create dataset with centroids
            auto dataset = vtkSmartPointer<vtkUnstructuredGrid>::New();
            dataset->SetPoints(centroids);
            
            // Create a mapping array to preserve original element IDs
            auto idMapping = vtkSmartPointer<vtkIdTypeArray>::New();
            idMapping->SetName("OriginalIds");
            idMapping->SetNumberOfComponents(1);
            idMapping->SetNumberOfTuples(centroids->GetNumberOfPoints());
            for (vtkIdType i = 0; i < centroids->GetNumberOfPoints(); i++) {
                idMapping->SetValue(i, i);
            }
            dataset->GetPointData()->AddArray(idMapping);
            
            kdTree_->SetDataSet(dataset);
        }
        
        // Build the tree
        kdTree_->BuildLocator();
        
        // Debug: Print KD-tree info
        std::cout << "KD-tree built with " << kdTree_->GetDataSet()->GetNumberOfPoints() 
                  << " points for " << (sourceDataLocation_ == DataLocation::NODAL ? "nodal" : "elemental")
                  << " data" << std::endl;
    }
    
    // Read mesh from file (supports various formats)
    vtkSmartPointer<vtkUnstructuredGrid> readMesh(const std::string& filename, 
        const std::vector<double>& origin = {0.0, 0.0, 0.0},
        const std::vector<double>& spacing = {1.0, 1.0, 1.0}) {
        std::string extension = filename.substr(filename.find_last_of(".") + 1);
        std::cout << "Reading file: " << filename << " (." << extension << ")" << std::endl;
        
        if (extension == "vtu") {
            auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
            reader->SetFileName(filename.c_str());
            reader->Update();
            auto output = reader->GetOutput();
            std::cout << "Read " << output->GetNumberOfPoints() << " points and " 
                      << output->GetNumberOfCells() << " cells" << std::endl;
            return output;
        }
        else if (extension == "inp") {
            return AbaqusInpReader::read(filename);
        }
        else if (extension == "npy") {
            // Add a materialMap parameter
            return NpyReader::read(filename, origin, spacing, targetMesh_, config_.materialMap);
        }
        else if (extension == "h5" || extension == "hdf5") {
            std::cerr << "HDF5 format support not available" << std::endl;
            return nullptr;
        }
        
        std::cerr << "Unsupported file format: " << extension << std::endl;
        return nullptr;
    }
    
    // Write mesh to file
    bool writeMesh(const std::string& filename, vtkSmartPointer<vtkUnstructuredGrid> mesh) {
        std::string extension = filename.substr(filename.find_last_of(".") + 1);
        
        if (extension == "vtu") {
            auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
            writer->SetFileName(filename.c_str());
            writer->SetInputData(mesh);
            writer->Write();
            std::cout << "Wrote file: " << filename << std::endl;
            return true;
        } 
        else if (extension == "inp") {
            // Check if target is inp and output is inp
            std::string targetExtension = targetFile_.substr(targetFile_.find_last_of(".") + 1);
            
            if (targetExtension == "inp" && extension == "inp" && 
                targetMesh_ && mappedDataArray_) {
                std::cout << "Appending field data to existing INP file" << std::endl;
                
                // Determine target location (nodal or elemental)
                std::string targetLocation = (targetDataLocation_ == DataLocation::NODAL) ? 
                    "nodal" : "elemental";
                
                // Use the current data array name
                std::string dataName = mappedDataArray_->GetName();
                
                return AbaqusInpReader::AppendFieldData(targetFile_, filename, mesh, dataName, targetLocation);
            }
            
            // Default case - write a new INP file
            return AbaqusInpReader::write(filename, mesh);
        }
        
        std::cerr << "Unsupported output format: " << extension << std::endl;
        return false;
    }
    
    // Find data array by name
    vtkSmartPointer<vtkDataArray> findDataArray(vtkSmartPointer<vtkUnstructuredGrid> mesh, 
                                              const std::string& name) {
        // Try point data
        auto pointData = mesh->GetPointData();
        auto array = pointData->GetArray(name.c_str());
        if (array) {
            std::cout << "Found point data array: " << name << std::endl;
            return array;
        }
        
        // Try cell data
        auto cellData = mesh->GetCellData();
        array = cellData->GetArray(name.c_str());
        if (array) {
            std::cout << "Found cell data array: " << name << std::endl;
            return array;
        }
        
        return nullptr;
    }
    
    // Find the first data array in the mesh
    vtkSmartPointer<vtkDataArray> findFirstDataArray(vtkSmartPointer<vtkUnstructuredGrid> mesh) {
        auto pointData = mesh->GetPointData();
        if (pointData->GetNumberOfArrays() > 0) {
            auto array = pointData->GetArray(0);
            std::cout << "Using first point data array: " << array->GetName() << std::endl;
            return array;
        }
        
        auto cellData = mesh->GetCellData();
        if (cellData->GetNumberOfArrays() > 0) {
            auto array = cellData->GetArray(0);
            std::cout << "Using first cell data array: " << array->GetName() << std::endl;
            return array;
        }
        
        return nullptr;
    }
    
    // Create a similarly-typed array
    vtkSmartPointer<vtkDataArray> createSimilarArray(vtkSmartPointer<vtkDataArray> source,
                                                  const std::string& name,
                                                  vtkIdType numTuples) {
        auto newArray = vtkSmartPointer<vtkDataArray>::Take(
            source->NewInstance());
        newArray->SetName(name.c_str());
        newArray->SetNumberOfComponents(source->GetNumberOfComponents());
        newArray->SetNumberOfTuples(numTuples);
        return newArray;
    }
};

// Main program
int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <source_file> <target_file> <output_file> [options]" << std::endl;
        std::cerr << "Options:" << std::endl;
        std::cerr << "  --gaussian         Use Gaussian weighting" << std::endl;
        std::cerr << "  --direct           Use direct (nearest point) mapping (default)" << std::endl;
        std::cerr << "  --npoint           Use N-point averaging" << std::endl;
        std::cerr << "  --proximity        Use proximity weighting" << std::endl;
        std::cerr << "  --radius <value>   Set search radius" << std::endl;
        std::cerr << "  --neighbors <n>    Set number of neighbor points" << std::endl;
        std::cerr << "  --sigma <value>    Set Gaussian sigma parameter" << std::endl;
        std::cerr << "  --smooth           Apply global smoothing" << std::endl;
        std::cerr << "  --smooth-iterations <n> Number of smoothing iterations" << std::endl;
        std::cerr << "  --data-name <name> Specify data array name to map" << std::endl;
        std::cerr << "  --target-location <nodal|elemental> Target data location (default: auto)" << std::endl;
        std::cerr << "  --origin <x,y,z>   Set origin coordinates for NPY files (default: 0,0,0)" << std::endl;
        std::cerr << "  --spacing <dx,dy,dz> Set spacing between points for NPY files (default: 1,1,1)" << std::endl;
        std::cerr << "  --material <id> <name> Define material ID and name (can be used multiple times)" << std::endl;
        return 1;
    }
    
    std::string sourceFile = argv[1];
    std::string targetFile = argv[2];
    std::string outputFile = argv[3];
    std::string dataName = "";
    std::string targetLocation = "auto";
    
    std::cout << "Source file: " << sourceFile << std::endl;
    std::cout << "Target file: " << targetFile << std::endl;
    std::cout << "Output file: " << outputFile << std::endl;
    
    // Parse options
    TensorMapper::Config config;
    config.useKdTree = true; // Always use KD-tree for efficiency
    
    for (int i = 4; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--gaussian") {
            config.method = TensorMapper::TransferMethod::GAUSSIAN;
            std::cout << "Using Gaussian mapping" << std::endl;
        } else if (arg == "--direct") {
            config.method = TensorMapper::TransferMethod::DIRECT;
            std::cout << "Using Direct mapping" << std::endl;
        } else if (arg == "--npoint") {
            config.method = TensorMapper::TransferMethod::N_POINT;
            std::cout << "Using N-Point mapping" << std::endl;
        } else if (arg == "--proximity") {
            config.method = TensorMapper::TransferMethod::PROXIMITY;
            std::cout << "Using Proximity mapping" << std::endl;
        } else if (arg == "--radius" && i + 1 < argc) {
            config.searchRadius = std::stod(argv[++i]);
            std::cout << "Search radius: " << config.searchRadius << std::endl;
        } else if (arg == "--neighbors" && i + 1 < argc) {
            config.numNeighbors = std::stoi(argv[++i]);
            std::cout << "Number of neighbors: " << config.numNeighbors << std::endl;
        } else if (arg == "--sigma" && i + 1 < argc) {
            config.gaussianSigma = std::stod(argv[++i]);
            std::cout << "Gaussian sigma: " << config.gaussianSigma << std::endl;
        } else if (arg == "--smooth") {
            config.applyGlobalSmoothing = true;
            std::cout << "Global smoothing enabled" << std::endl;
        } else if (arg == "--smooth-iterations" && i + 1 < argc) {
            config.smoothingIterations = std::stoi(argv[++i]);
            std::cout << "Smoothing iterations: " << config.smoothingIterations << std::endl;
        } else if (arg == "--data-name" && i + 1 < argc) {
            dataName = argv[++i];
            std::cout << "Data array name: " << dataName << std::endl;
        } else if (arg == "--target-location" && i + 1 < argc) {
            targetLocation = argv[++i];
            std::cout << "Target location: " << targetLocation << std::endl;
        } else if (arg == "--relax" && i + 1 < argc) {
            config.relaxationFactor = std::stod(argv[++i]);
            std::cout << "Relaxation factor: " << config.relaxationFactor << std::endl;
        }else if (arg == "--material" && i + 2 < argc) {
            int materialId = std::stoi(argv[++i]);
            std::string materialName = argv[++i];
            config.materialMap[materialId] = materialName;
            std::cout << "Material mapping: " << materialId << " = " << materialName << std::endl;
        } else if (arg == "--origin" && i + 1 < argc) {
            std::string originStr = argv[++i];
            std::vector<double> originValues;
            
            // Parse comma-separated values
            std::stringstream ss(originStr);
            std::string item;
            while (std::getline(ss, item, ',')) {
                try {
                    originValues.push_back(std::stod(item));
                } catch(const std::exception& e) {
                    std::cerr << "Warning: Invalid origin value: " << item << std::endl;
                }
            }
            
            if (originValues.size() > 0) {
                config.origin = originValues;
                std::cout << "Origin: [";
                for (size_t j = 0; j < originValues.size(); j++) {
                    std::cout << originValues[j];
                    if (j < originValues.size() - 1) std::cout << ", ";
                }
                std::cout << "]" << std::endl;
            }
        } else if (arg == "--spacing" && i + 1 < argc) {
            std::string spacingStr = argv[++i];
            std::vector<double> spacingValues;
            
            // Parse comma-separated values
            std::stringstream ss(spacingStr);
            std::string item;
            while (std::getline(ss, item, ',')) {
                try {
                    spacingValues.push_back(std::stod(item));
                } catch(const std::exception& e) {
                    std::cerr << "Warning: Invalid spacing value: " << item << std::endl;
                }
            }
            
            if (spacingValues.size() > 0) {
                config.spacing = spacingValues;
                std::cout << "Spacing: [";
                for (size_t j = 0; j < spacingValues.size(); j++) {
                    std::cout << spacingValues[j];
                    if (j < spacingValues.size() - 1) std::cout << ", ";
                }
                std::cout << "]" << std::endl;
            }
        }
    }
    
    // Create mapper
    TensorMapper mapper(config);
    
    // Load source and target meshes
    if (!mapper.loadTarget(targetFile)) {
        std::cerr << "Failed to load target" << std::endl;
        return 1;
    }
    if (!mapper.loadSource(sourceFile, dataName)) {
        std::cerr << "Failed to load source" << std::endl;
        return 1;
    }

    // Determine target data location
    DataLocation targetDataLoc = DataLocation::NODAL;
    if (targetLocation == "elemental") {
        targetDataLoc = DataLocation::ELEMENTAL;
    } else if (targetLocation == "nodal") {
        targetDataLoc = DataLocation::NODAL;
    } else if (targetLocation != "auto") {
        std::cerr << "Warning: Unknown target location '" << targetLocation 
                  << "', using nodal" << std::endl;
    }
    
    // Perform mapping
    if (!mapper.mapData("", targetDataLoc)) {
        std::cerr << "Mapping failed" << std::endl;
        return 1;
    }
    
    // Save result
    if (!mapper.saveResult(outputFile)) {
        std::cerr << "Failed to save result" << std::endl;
        return 1;
    }
    
    std::cout << "Mapping complete. Result saved to " << outputFile << std::endl;
    return 0;
}