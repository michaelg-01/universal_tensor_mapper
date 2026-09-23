#ifndef NPY_READER_H
#define NPY_READER_H

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <cmath>
#include <atomic>
#include <omp.h>
#include <limits>
#include <map>

// VTK includes
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkIntArray.h>
#include <vtkIdList.h>
#include <vtkCell.h>

// Include cnpy for NPY file reading
#include "cnpy.h"

class NpyReader {
public:
    // New method with material map parameter
    static vtkSmartPointer<vtkUnstructuredGrid> read(
        const std::string& filename,
        const std::vector<double>& origin,
        const std::vector<double>& spacing,
        vtkSmartPointer<vtkUnstructuredGrid> targetMesh,
        const std::map<int, std::string>& materialMap) {
        
        try {
            // Target mesh must be provided
            if (!targetMesh) {
                std::cerr << "Error: Target mesh must be provided for material mapping" << std::endl;
                return nullptr;
            }
            
            // Material map must not be empty
            if (materialMap.empty()) {
                std::cerr << "Warning: Material map is empty, no mapping will be performed" << std::endl;
            }
            
            // Load the NPY file
            std::cout << "Reading NPY file: " << filename << std::endl;
            cnpy::NpyArray data_array = cnpy::npy_load(filename);
            
            // Get dimensions
            if (data_array.shape.size() < 3) {
                std::cerr << "Error: NPY must be 3D array, found " << data_array.shape.size() << " dimensions" << std::endl;
                return nullptr;
            }
            
            size_t depth = data_array.shape[0];
            size_t height = data_array.shape[1];
            size_t width = data_array.shape[2];
            std::cout << "NPY dimensions: " << depth << " x " << height << " x " << width << std::endl;
            
            // Verify data type
            if (data_array.word_size != 1) {
                std::cerr << "Warning: Expected uint8 (1-byte) data, found " << data_array.word_size << " bytes" << std::endl;
            }
            
            // Convert NPY data to a simple vector for easier access
            std::vector<uint8_t> material_array;
            if (data_array.word_size == 1) {
                material_array.resize(data_array.num_vals);
                std::memcpy(material_array.data(), data_array.data<uint8_t>(), data_array.num_vals);
            } else {
                // Handle other data types if needed
                std::cerr << "Error: Only uint8 data type is supported" << std::endl;
                return nullptr;
            }
            
            // Ensure origin and spacing are 3D
            std::vector<double> origin_3d = origin;
            std::vector<double> spacing_3d = spacing;
            while (origin_3d.size() < 3) origin_3d.push_back(0.0);
            while (spacing_3d.size() < 3) spacing_3d.push_back(1.0);

            for (auto& val : origin_3d) val *= 1000.0;  // m to mm
            for (auto& val : spacing_3d) val *= 1000.0; // m to mm
            
            // Get mesh info
            vtkIdType numCells = targetMesh->GetNumberOfCells();
            std::cout << "Target mesh has " << numCells << " cells" << std::endl;
            
            // Print the material IDs and names
            std::cout << "Material ID mapping:" << std::endl;
            for (const auto& [id, name] : materialMap) {
                std::cout << "  " << id << " = " << name << std::endl;
            }
            
            // Create arrays to store material counts
            std::vector<vtkSmartPointer<vtkIntArray>> materialCountArrays;
            for (const auto& [materialId, materialName] : materialMap) {
                auto countArray = vtkSmartPointer<vtkIntArray>::New();
                countArray->SetName((materialName + "_Count").c_str());
                countArray->SetNumberOfComponents(1);
                countArray->SetNumberOfTuples(numCells);
                
                // Initialize to zero
                for (vtkIdType i = 0; i < numCells; i++) {
                    countArray->SetValue(i, 0);
                }
                materialCountArrays.push_back(countArray);
            }
            
            // Create dominant material array
            auto dominantMaterialArray = vtkSmartPointer<vtkIntArray>::New();
            dominantMaterialArray->SetName("DominantMaterial");
            dominantMaterialArray->SetNumberOfComponents(1);
            dominantMaterialArray->SetNumberOfTuples(numCells);
            
            // Create dominant material name array
            auto dominantMaterialNameArray = vtkSmartPointer<vtkIntArray>::New();
            dominantMaterialNameArray->SetName("DominantMaterialName");
            dominantMaterialNameArray->SetNumberOfComponents(1);
            dominantMaterialNameArray->SetNumberOfTuples(numCells);
            
            // Process all cells
            std::cout << "Processing all " << numCells << " cells..." << std::endl;
            
            // Progress reporting
            vtkIdType progress_interval = numCells / 20;
            std::atomic<vtkIdType> processed_cells(0);
            
            // Store material name to index mapping
            std::map<int, int> materialIdToArrayIndex;
            int idx = 0;
            for (const auto& [materialId, materialName] : materialMap) {
                materialIdToArrayIndex[materialId] = idx++;
            }
            
            // Cell-based approach with coordinate flipping
            #pragma omp parallel for
            for (int cellId = 0; cellId < numCells; cellId++) {
                try {
                    // Get the cell
                    vtkSmartPointer<vtkCell> cell = targetMesh->GetCell(cellId);
                    
                    // Get cell points
                    vtkSmartPointer<vtkIdList> pointIds = vtkSmartPointer<vtkIdList>::New();
                    targetMesh->GetCellPoints(cellId, pointIds);
                    
                    // Get bounding box of the cell
                    double bounds[6]; // xmin, xmax, ymin, ymax, zmin, zmax
                    bounds[0] = bounds[2] = bounds[4] = std::numeric_limits<double>::max();
                    bounds[1] = bounds[3] = bounds[5] = -std::numeric_limits<double>::max();
                    
                    for (vtkIdType i = 0; i < pointIds->GetNumberOfIds(); i++) {
                        double point[3];
                        targetMesh->GetPoint(pointIds->GetId(i), point);
                        
                        // Update bounds
                        bounds[0] = std::min(bounds[0], point[0]); // xmin
                        bounds[1] = std::max(bounds[1], point[0]); // xmax
                        bounds[2] = std::min(bounds[2], point[1]); // ymin
                        bounds[3] = std::max(bounds[3], point[1]); // ymax
                        bounds[4] = std::min(bounds[4], point[2]); // zmin
                        bounds[5] = std::max(bounds[5], point[2]); // zmax
                    }
                    
                    // Convert bounding box to voxel coordinates (with x-z flipping)
                    int z_min = static_cast<int>((bounds[0] - origin_3d[0]) / spacing_3d[0]);
                    int z_max = static_cast<int>((bounds[1] - origin_3d[0]) / spacing_3d[0]) + 1;
                    int y_min = static_cast<int>((bounds[2] - origin_3d[1]) / spacing_3d[1]);
                    int y_max = static_cast<int>((bounds[3] - origin_3d[1]) / spacing_3d[1]) + 1;
                    int x_min = static_cast<int>((bounds[4] - origin_3d[2]) / spacing_3d[2]);
                    int x_max = static_cast<int>((bounds[5] - origin_3d[2]) / spacing_3d[2]) + 1;
                    
                    // Clamp to array bounds
                    z_min = std::max(0, std::min(z_min, static_cast<int>(width) - 1));
                    z_max = std::max(0, std::min(z_max, static_cast<int>(width)));
                    y_min = std::max(0, std::min(y_min, static_cast<int>(height) - 1));
                    y_max = std::max(0, std::min(y_max, static_cast<int>(height)));
                    x_min = std::max(0, std::min(x_min, static_cast<int>(depth) - 1));
                    x_max = std::max(0, std::min(x_max, static_cast<int>(depth)));
                    
                    // Count materials within cell
                    std::map<int, int> cell_material_counts;
                    
                    for (int x = x_min; x < x_max; x++) {
                        for (int y = y_min; y < y_max; y++) {
                            for (int z = z_min; z < z_max; z++) {
                                // Convert voxel indices to world coordinates (with flipping)
                                double worldPoint[3];
                                worldPoint[0] = origin_3d[0] + z * spacing_3d[0];  // z to x
                                worldPoint[1] = origin_3d[1] + y * spacing_3d[1];  // y to y
                                worldPoint[2] = origin_3d[2] + x * spacing_3d[2];  // x to z
                                
                                // Check if point is inside cell
                                double closestPoint[3], pcoords[3], dist2;
                                int subId;
                                double* weights = new double[cell->GetNumberOfPoints()];
                                
                                if (cell->EvaluatePosition(worldPoint, closestPoint, subId, pcoords, dist2, weights) == 1) {
                                    // Point is inside cell, check material ID
                                    size_t idx = (x * height * width) + (y * width) + z;
                                    if (idx < material_array.size()) {
                                        uint8_t materialId = material_array[idx];
                                        if (materialMap.find(materialId) != materialMap.end()) {
                                            cell_material_counts[materialId]++;
                                        }
                                    }
                                }
                                
                                delete[] weights;
                            }
                        }
                    }
                    
                    // Update material count arrays
                    for (const auto& [materialId, count] : cell_material_counts) {
                        // Skip if material not in our map
                        if (materialIdToArrayIndex.find(materialId) == materialIdToArrayIndex.end()) {
                            continue;
                        }
                        
                        // Update count in appropriate array
                        int arrayIndex = materialIdToArrayIndex[materialId];
                        materialCountArrays[arrayIndex]->SetValue(cellId, count);
                    }
                    
                    // Find dominant material
                    int maxCount = 0;
                    int dominantMaterial = 0;
                    for (const auto& [materialId, count] : cell_material_counts) {
                        if (count > maxCount && materialMap.find(materialId) != materialMap.end()) {
                            maxCount = count;
                            dominantMaterial = materialId;
                        }
                    }
                    dominantMaterialArray->SetValue(cellId, dominantMaterial);
                    
                    // Look up name for dominant material
                    std::string dominantName = "None";
                    if (dominantMaterial > 0 && materialMap.find(dominantMaterial) != materialMap.end()) {
                        dominantName = materialMap.at(dominantMaterial);
                    }
                    
                    // Report progress
                    vtkIdType current = ++processed_cells;
                    if (current % progress_interval == 0 || current == numCells) {
                        #pragma omp critical
                        {
                            double progress = (100.0 * current) / numCells;
                            std::cout << "Progress: " << progress << "% (" << current << "/" << numCells << ")" << std::endl;
                        }
                    }
                }
                catch (const std::exception& e) {
                    #pragma omp critical
                    {
                        std::cerr << "Error processing cell " << cellId << ": " << e.what() << std::endl;
                    }
                }
            }
            
            // Add arrays to target mesh automatically as npy doesnt have any properties to specify on input
            for (auto& countArray : materialCountArrays) {
                targetMesh->GetCellData()->AddArray(countArray);
            }
            targetMesh->GetCellData()->AddArray(dominantMaterialArray);
            
            std::cout << "Material mapping complete" << std::endl;
            return targetMesh;
        }
        catch (const std::exception& e) {
            std::cerr << "Error in NpyReader: " << e.what() << std::endl;
            return nullptr;
        }
    }
    
    // Legacy method for backward compatibility
    static vtkSmartPointer<vtkUnstructuredGrid> read(
        const std::string& filename,
        const std::vector<double>& origin = {0.0, 0.0, 0.0},
        const std::vector<double>& spacing = {1.0, 1.0, 1.0},
        vtkSmartPointer<vtkUnstructuredGrid> targetMesh = nullptr) {
        
        // Default material map with single material ID
        std::map<int, std::string> defaultMap = {{1, "Material_1"}};
        return read(filename, origin, spacing, targetMesh, defaultMap);
    }
};

#endif // NPY_READER_H