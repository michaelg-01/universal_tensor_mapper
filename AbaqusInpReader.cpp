// AbaqusInpReader.cpp
#include "AbaqusInpReader.h"

#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <set>

// VTK includes
#include <vtkTetra.h>
#include <vtkHexahedron.h>
#include <vtkWedge.h>
#include <vtkTriangle.h>
#include <vtkQuad.h>
#include <vtkQuadraticTetra.h>
#include <vtkQuadraticHexahedron.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkCell.h>
#include <vtkNew.h>

// Static method to read an Abaqus INP file
vtkSmartPointer<vtkUnstructuredGrid> AbaqusInpReader::read(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open INP file: " << filename << std::endl;
        return nullptr;
    }
    
    auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    auto points = vtkSmartPointer<vtkPoints>::New();
    std::map<int, vtkIdType> nodeIdMap;
    
    std::string line;
    bool inNodeSection = false;
    bool inElementSection = false;
    bool inElsetSection = false;
    std::string currentElementType;
    std::string currentElementSet;
    std::string currentElsetName;
    
    // Maps for storing element definitions by type
    std::map<std::string, std::vector<std::pair<int, std::vector<int>>>> elementsByType;
    
    // Map to track element sets
    std::map<std::string, std::set<int>> elementSets;
    
    // Store multi-line elements during parsing
    std::string multiLineElement;
    int multiLineElementId = -1;
    
    // First pass: Read nodes and elements
    while (std::getline(file, line)) {
        // Trim whitespace
        line = trimString(line);
        
        // Skip empty lines
        if (line.empty()) continue;
        
        // Handle line continuations for elements
        if (!multiLineElement.empty()) {
            // Continue building the element definition
            multiLineElement += line;
            
            // Check if this is the end of the continuation
            if (line.back() != ',') {
                // Process the complete element definition
                std::vector<int> nodeIds;
                if (parseElementLine(multiLineElement, multiLineElementId, nodeIds)) {
                    elementsByType[currentElementType].push_back(std::make_pair(multiLineElementId, nodeIds));
                    
                    // Add to current set if applicable
                    if (!currentElementSet.empty()) {
                        elementSets[currentElementSet].insert(multiLineElementId);
                    }
                }
                multiLineElement.clear();
                multiLineElementId = -1;
            }
            continue;
        }
        
        // Process keywords (lines starting with *)
        if (line[0] == '*') {
            inNodeSection = false;
            inElementSection = false;
            inElsetSection = false;
            
            std::string keyword = toUpperString(line);
            
            if (keyword.find("*NODE") == 0) {
                inNodeSection = true;
                continue;
            } 
            else if (keyword.find("*ELEMENT") == 0) {
                inElementSection = true;
                
                // Extract element type
                size_t typePos = keyword.find("TYPE=");
                if (typePos != std::string::npos) {
                    size_t start = typePos + 5;
                    size_t end = keyword.find(",", start);
                    if (end == std::string::npos) end = keyword.length();
                    currentElementType = keyword.substr(start, end - start);
                    // Remove any quotes
                    currentElementType.erase(
                        std::remove(currentElementType.begin(), currentElementType.end(), '\"'), 
                        currentElementType.end());
                    currentElementType.erase(
                        std::remove(currentElementType.begin(), currentElementType.end(), '\''), 
                        currentElementType.end());
                }
                
                // Extract element set name if present
                size_t elsetPos = keyword.find("ELSET=");
                if (elsetPos != std::string::npos) {
                    size_t start = elsetPos + 6;
                    size_t end = keyword.find(",", start);
                    if (end == std::string::npos) end = keyword.length();
                    currentElementSet = keyword.substr(start, end - start);
                    // Remove any quotes
                    currentElementSet.erase(
                        std::remove(currentElementSet.begin(), currentElementSet.end(), '\"'), 
                        currentElementSet.end());
                    currentElementSet.erase(
                        std::remove(currentElementSet.begin(), currentElementSet.end(), '\''), 
                        currentElementSet.end());
                } else {
                    currentElementSet = ""; // No element set specified
                }
                
                continue;
            }
            else if (keyword.find("*ELSET") == 0) {
                inElsetSection = true;
                
                // Extract element set name
                size_t namePos = keyword.find("ELSET=");
                if (namePos != std::string::npos) {
                    size_t start = namePos + 6;
                    size_t end = keyword.find(",", start);
                    if (end == std::string::npos) end = keyword.length();
                    currentElsetName = keyword.substr(start, end - start);
                    // Remove any quotes
                    currentElsetName.erase(
                        std::remove(currentElsetName.begin(), currentElsetName.end(), '\"'), 
                        currentElsetName.end());
                    currentElsetName.erase(
                        std::remove(currentElsetName.begin(), currentElsetName.end(), '\''), 
                        currentElsetName.end());
                }
                
                continue;
            }
            
            // Skip other keywords for now
            continue;
        } 
        // Process data lines
        else {
            if (inNodeSection) {
                parseNodeLine(line, points, nodeIdMap);
            } 
            else if (inElementSection) {
                // Handle possible line continuation
                if (line.back() == ',') {
                    // This is the start of a multi-line element
                    std::vector<std::string> tokens = splitString(line, ',');
                    if (!tokens.empty()) {
                        multiLineElementId = std::stoi(trimString(tokens[0]));
                        multiLineElement = line;
                    }
                } else {
                    // Single line element
                    std::vector<int> nodeIds;
                    int elementId;
                    if (parseElementLine(line, elementId, nodeIds)) {
                        elementsByType[currentElementType].push_back(std::make_pair(elementId, nodeIds));
                        
                        // Add to current set if applicable
                        if (!currentElementSet.empty()) {
                            elementSets[currentElementSet].insert(elementId);
                        }
                    }
                }
            }
            else if (inElsetSection) {
                // Process element set line
                std::vector<std::string> tokens = splitString(line, ',');
                for (const auto& token : tokens) {
                    std::string trimmedToken = trimString(token);
                    if (!trimmedToken.empty()) {
                        try {
                            elementSets[currentElsetName].insert(std::stoi(trimmedToken));
                        } catch(const std::exception& e) {
                            std::cerr << "Error parsing element ID in ELSET: " << trimmedToken << std::endl;
                        }
                    }
                }
            }
        }
    }
    
    // Set points in the grid
    grid->SetPoints(points);
    
    // Second pass: Create VTK cells from the elements
    vtkNew<vtkIdTypeArray> originalIds;
    originalIds->SetName("OriginalElementIds");
    originalIds->SetNumberOfComponents(1);
    
    // Create cells and populate originalIds
    createCellsFromElements(grid, elementsByType, nodeIdMap, originalIds);
    
    // Add original element IDs to cell data
    grid->GetCellData()->AddArray(originalIds);
    
    // Create element set arrays
    vtkNew<vtkIntArray> setIdArray;
    setIdArray->SetName("set_id");
    setIdArray->SetNumberOfComponents(1);
    setIdArray->SetNumberOfTuples(grid->GetNumberOfCells());
    setIdArray->Fill(0);
    
    // For each set, create a mask array and update set_id
    int setId = 1;
    for (const auto& setEntry : elementSets) {
        const std::string& setName = setEntry.first;
        const std::set<int>& elemIds = setEntry.second;
        
        vtkNew<vtkIntArray> setMask;
        setMask->SetName(setName.c_str());
        setMask->SetNumberOfComponents(1);
        setMask->SetNumberOfTuples(grid->GetNumberOfCells());
        setMask->Fill(0);
        
        for (vtkIdType i = 0; i < originalIds->GetNumberOfTuples(); i++) {
            vtkIdType abaqusElemId = originalIds->GetValue(i);
            if (elemIds.find(abaqusElemId) != elemIds.end()) {
                setIdArray->SetValue(i, setId);
                setMask->SetValue(i, 1);
            }
        }
        
        grid->GetCellData()->AddArray(setMask);
        setId++;
    }
    
    grid->GetCellData()->AddArray(setIdArray);
    
    // Create a default data field if none exists
    if (grid->GetPointData()->GetNumberOfArrays() == 0) {
        vtkNew<vtkDoubleArray> defaultField;
        defaultField->SetName("DefaultDistanceField");
        defaultField->SetNumberOfComponents(1);
        defaultField->SetNumberOfTuples(grid->GetNumberOfPoints());
        
        for (vtkIdType i = 0; i < grid->GetNumberOfPoints(); i++) {
            double point[3];
            grid->GetPoint(i, point);
            double distance = sqrt(point[0]*point[0] + point[1]*point[1] + point[2]*point[2]);
            defaultField->SetValue(i, distance);
        }
        
        grid->GetPointData()->AddArray(defaultField);
        std::cout << "Created default point data field" << std::endl;
    }
    
    std::cout << "Read " << points->GetNumberOfPoints() << " points and " 
              << grid->GetNumberOfCells() << " cells from INP file" << std::endl;
    std::cout << "Found " << elementSets.size() << " element sets" << std::endl;
              
    return grid;
}

void AbaqusInpReader::parseNodeLine(const std::string& line, vtkSmartPointer<vtkPoints> points, 
                                    std::map<int, vtkIdType>& nodeIdMap) {
    std::vector<std::string> tokens = splitString(line, ',');
    
    if (tokens.size() >= 4) {
        int nodeId = std::stoi(trimString(tokens[0]));
        double x = std::stod(trimString(tokens[1]));
        double y = std::stod(trimString(tokens[2]));
        double z = std::stod(trimString(tokens[3]));
        
        // Add point and map Abaqus ID to VTK ID
        vtkIdType pointId = points->InsertNextPoint(x, y, z);
        nodeIdMap[nodeId] = pointId;
    }
}

bool AbaqusInpReader::parseElementLine(const std::string& line, int& elementId, std::vector<int>& nodeIds) {
    std::vector<std::string> tokens = splitString(line, ',');
    if (tokens.empty()) return false;
    
    elementId = std::stoi(trimString(tokens[0]));
    nodeIds.clear();
    
    for (size_t i = 1; i < tokens.size(); ++i) {
        std::string nodeIdStr = trimString(tokens[i]);
        if (!nodeIdStr.empty()) {
            nodeIds.push_back(std::stoi(nodeIdStr));
        }
    }
    
    return !nodeIds.empty();
}

void AbaqusInpReader::createCellsFromElements(vtkSmartPointer<vtkUnstructuredGrid> grid,
                                             const std::map<std::string, std::vector<std::pair<int, std::vector<int>>>>& elementsByType,
                                             const std::map<int, vtkIdType>& nodeIdMap,
                                             vtkIdTypeArray* originalIds) {
    
    for (const auto& typeEntry : elementsByType) {
        const std::string& elemType = typeEntry.first;
        const auto& elements = typeEntry.second;
        
        for (const auto& element : elements) {
            int elemId = element.first;
            const std::vector<int>& abaqusNodeIds = element.second;
            
            // Convert Abaqus node IDs to VTK point IDs
            vtkNew<vtkIdList> pointIds;
            for (int abaqusId : abaqusNodeIds) {
                if (nodeIdMap.find(abaqusId) != nodeIdMap.end()) {
                    pointIds->InsertNextId(nodeIdMap.at(abaqusId));
                } else {
                    std::cerr << "Warning: Node ID " << abaqusId << " not found in node map" << std::endl;
                }
            }
            
            // Map Abaqus element types to VTK cell types
            int vtkCellType = mapAbaqusElementTypeToVtk(elemType, pointIds->GetNumberOfIds());
            
            if (vtkCellType != VTK_EMPTY_CELL && pointIds->GetNumberOfIds() > 0) {
                grid->InsertNextCell(vtkCellType, pointIds);
                originalIds->InsertNextValue(elemId);
            } else {
                std::cerr << "Warning: Unsupported element type: " << elemType << " with " 
                          << pointIds->GetNumberOfIds() << " nodes" << std::endl;
            }
        }
    }
}

int AbaqusInpReader::mapAbaqusElementTypeToVtk(const std::string& abaqusType, int numNodes) {
    // Convert to uppercase for case-insensitive comparison
    std::string type = toUpperString(abaqusType);
    
    // Standard linear elements
    if (type == "C3D4" || type == "C3D4H") {
        return VTK_TETRA;
    } else if (type == "C3D8" || type == "C3D8R" || type == "C3D8I" || type == "C3D8RH") {
        return VTK_HEXAHEDRON;
    } else if (type == "C3D6" || type == "C3D6H") {
        return VTK_WEDGE;
    } else if (type == "S3" || type == "S3R") {
        return VTK_TRIANGLE;
    } else if (type == "S4" || type == "S4R") {
        return VTK_QUAD;
    }
    // Quadratic elements
    else if (type == "C3D10" || type == "C3D10H" || type == "C3D10M") {
        return VTK_QUADRATIC_TETRA;
    } else if (type == "C3D20" || type == "C3D20R") {
        return VTK_QUADRATIC_HEXAHEDRON;
    }
    
    // For unknown types, try to infer from node count
    if (numNodes == 4) {
        return VTK_TETRA;
    } else if (numNodes == 8) {
        return VTK_HEXAHEDRON;
    } else if (numNodes == 6) {
        return VTK_WEDGE;
    } else if (numNodes == 3) {
        return VTK_TRIANGLE;
    } else if (numNodes == 10) {
        return VTK_QUADRATIC_TETRA;
    } else if (numNodes == 20) {
        return VTK_QUADRATIC_HEXAHEDRON;
    }
    
    return VTK_EMPTY_CELL;
}

bool AbaqusInpReader::write(const std::string& filename, vtkSmartPointer<vtkUnstructuredGrid> mesh) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open output file: " << filename << std::endl;
        return false;
    }
    
    // Write header
    file << "** Generated by TensorMapper" << std::endl;
    file << "**" << std::endl;
    
    // Write nodes
    file << "*NODE" << std::endl;
    for (vtkIdType i = 0; i < mesh->GetNumberOfPoints(); i++) {
        double point[3];
        mesh->GetPoint(i, point);
        file << (i+1) << ", " << point[0] << ", " << point[1] << ", " << point[2] << std::endl;
    }
    
    // Group cells by type
    std::map<int, std::vector<vtkIdType>> cellsByType;
    for (vtkIdType i = 0; i < mesh->GetNumberOfCells(); i++) {
        vtkCell* cell = mesh->GetCell(i);
        cellsByType[cell->GetCellType()].push_back(i);
    }
    
    // Extract element sets if they exist
    std::map<std::string, std::set<vtkIdType>> elementSets;
    vtkDataArray* setIdArray = mesh->GetCellData()->GetArray("set_id");
    
    if (!setIdArray) {
        // Look for individual set arrays
        for (int i = 0; i < mesh->GetCellData()->GetNumberOfArrays(); i++) {
            vtkDataArray* array = mesh->GetCellData()->GetArray(i);
            std::string arrayName = array->GetName();
            
            // Skip known arrays
            if (arrayName == "OriginalElementIds") continue;
            
            // Check if this is a set mask (1 component int array)
            if (array->GetNumberOfComponents() == 1 && 
                (array->GetDataType() == VTK_INT || array->GetDataType() == VTK_ID_TYPE)) {
                // Collect cell IDs for this set
                for (vtkIdType j = 0; j < array->GetNumberOfTuples(); j++) {
                    if (array->GetComponent(j, 0) > 0) {
                        elementSets[arrayName].insert(j);
                    }
                }
            }
        }
    } else {
        // Use set_id array to identify sets
        std::map<int, std::string> setIdToName;
        
        // Find set name arrays
        for (int i = 0; i < mesh->GetCellData()->GetNumberOfArrays(); i++) {
            vtkDataArray* array = mesh->GetCellData()->GetArray(i);
            std::string arrayName = array->GetName();
            
            // Skip known arrays
            if (arrayName == "OriginalElementIds" || arrayName == "set_id") continue;
            
            // Check if this is a set mask
            if (array->GetNumberOfComponents() == 1 && 
                (array->GetDataType() == VTK_INT || array->GetDataType() == VTK_ID_TYPE)) {
                
                // Find the set ID for this array
                for (vtkIdType j = 0; j < array->GetNumberOfTuples(); j++) {
                    if (array->GetComponent(j, 0) > 0) {
                        int setId = static_cast<int>(setIdArray->GetComponent(j, 0));
                        if (setId > 0) {
                            setIdToName[setId] = arrayName;
                            break;
                        }
                    }
                }
            }
        }
        
        // Collect cells for each set
        for (vtkIdType i = 0; i < setIdArray->GetNumberOfTuples(); i++) {
            int setId = static_cast<int>(setIdArray->GetComponent(i, 0));
            if (setId > 0 && setIdToName.find(setId) != setIdToName.end()) {
                elementSets[setIdToName[setId]].insert(i);
            }
        }
    }
    
    // Write elements by type
    for (const auto& entry : cellsByType) {
        int vtkType = entry.first;
        const auto& cells = entry.second;
        
        std::string abaqusType = mapVtkElementTypeToAbaqus(vtkType);
        if (!abaqusType.empty()) {
            file << "*ELEMENT, TYPE=" << abaqusType << std::endl;
            
            for (vtkIdType cellId : cells) {
                vtkCell* cell = mesh->GetCell(cellId);
                file << (cellId+1) << ", ";
                
                // Get maximum points per line (8 is a good number for readability)
                const int maxPointsPerLine = 8;
                int remainingPoints = cell->GetNumberOfPoints();
                int pointsInCurrentLine = 0;
                
                for (vtkIdType j = 0; j < cell->GetNumberOfPoints(); j++) {
                    // Abaqus node IDs are 1-based
                    file << (cell->GetPointId(j) + 1);
                    
                    remainingPoints--;
                    pointsInCurrentLine++;
                    
                    // Add comma if needed
                    if (remainingPoints > 0) {
                        file << ", ";
                        
                        // Start a new line if max points per line is reached
                        if (pointsInCurrentLine >= maxPointsPerLine) {
                            file << std::endl;
                            // Continuation line indentation
                            file << "    ";
                            pointsInCurrentLine = 0;
                        }
                    }
                }
                file << std::endl;
            }
        }
    }
    
    // Write element sets
    for (const auto& entry : elementSets) {
        const std::string& setName = entry.first;
        const std::set<vtkIdType>& cellIds = entry.second;
        
        file << "*ELSET, ELSET=" << setName << std::endl;
        
        // Write cell IDs (16 per line is standard for Abaqus)
        int count = 0;
        for (vtkIdType cellId : cellIds) {
            file << (cellId+1);
            count++;
            
            if (count == cellIds.size()) {
                // Last element, no comma
                file << std::endl;
            } else if (count % 16 == 0) {
                // End of line
                file << std::endl;
            } else {
                // Comma separator
                file << ", ";
            }
        }
        
        // Add newline if needed
        if (!cellIds.empty() && cellIds.size() % 16 != 0) {
            file << std::endl;
        }
    }
    
    // Write field data if available
    writeFieldData(file, mesh, elementSets);
    //do fo points too
    writePointFieldData(file, mesh);
    
    file.close();
    std::cout << "Wrote INP file: " << filename << std::endl;
    return true;
}

std::string AbaqusInpReader::mapVtkElementTypeToAbaqus(int vtkType) {
    switch (vtkType) {
        case VTK_TETRA: return "C3D4";
        case VTK_HEXAHEDRON: return "C3D8";
        case VTK_WEDGE: return "C3D6";
        case VTK_TRIANGLE: return "S3";
        case VTK_QUAD: return "S4";
        case VTK_QUADRATIC_TETRA: return "C3D10";
        case VTK_QUADRATIC_HEXAHEDRON: return "C3D20";
        default: return "";
    }
}

void AbaqusInpReader::writeFieldData(std::ofstream& file, vtkSmartPointer<vtkUnstructuredGrid> mesh, const std::map<std::string, std::set<vtkIdType>>& elementSets) {
    // Get original element IDs if available
    vtkIdTypeArray* originalIds = vtkIdTypeArray::SafeDownCast(
        mesh->GetCellData()->GetArray("OriginalElementIds"));
    
    // Process cell data arrays
    vtkCellData* cellData = mesh->GetCellData();
    for (int i = 0; i < cellData->GetNumberOfArrays(); i++) {
        vtkDataArray* array = cellData->GetArray(i);
        std::string arrayName = array->GetName();
        
        // Skip known non-field arrays
        if (arrayName == "OriginalElementIds" || arrayName == "set_id") {
            continue;
        }
        
        // Check if this array is a set mask (all 0s and 1s)
        bool isSetMask = true;
        for (vtkIdType j = 0; j < array->GetNumberOfTuples() && isSetMask; j++) {
            double value = array->GetComponent(j, 0);
            if (value != 0.0 && value != 1.0) {
                isSetMask = false;
            }
        }
        
        // Skip set masks
        if (isSetMask) continue;
        
        // Create distribution table first (using valid RATIO label)
        file << "*DISTRIBUTION TABLE, NAME=" << arrayName << "_TABLE" << std::endl;
        file << "RATIO," << std::endl;
        
        // Then create the distribution that references the table
        file << "*DISTRIBUTION, NAME=" << arrayName << "_DIST, LOCATION=ELEMENT, TABLE=" << arrayName << "_TABLE" << std::endl;
        
        // Write a single default value line
        double defaultValue = 0.0;
        file << ", " << defaultValue << std::endl;
        
        // Write element-specific values
        for (vtkIdType j = 0; j < array->GetNumberOfTuples(); j++) {
            // Skip elements with zero values to keep file smaller
            if (array->GetComponent(j, 0) == 0.0) continue;
            
            // Get element ID
            vtkIdType elemId = j + 1; // Default to 1-based index
            if (originalIds && j < originalIds->GetNumberOfTuples()) {
                elemId = originalIds->GetValue(j);
            }
            
            // Write element ID and value
            file << elemId << ", " << array->GetComponent(j, 0) << std::endl;
        }
    }
}

void AbaqusInpReader::writePointFieldData(std::ofstream& file, vtkSmartPointer<vtkUnstructuredGrid> mesh) {
    // Process point data arrays
    vtkPointData* pointData = mesh->GetPointData();
    for (int i = 0; i < pointData->GetNumberOfArrays(); i++) {
        vtkDataArray* array = pointData->GetArray(i);
        std::string arrayName = array->GetName();
        
        // Skip DefaultDistanceField
        if (arrayName == "DefaultDistanceField") continue;
        
        // Create distribution table
        file << "*DISTRIBUTION TABLE, NAME=" << arrayName << "_TABLE" << std::endl;
        file << "RATIO," << std::endl;
        
        // Create distribution - use LOCATION=NODE for point data
        file << "*DISTRIBUTION, NAME=" << arrayName << "_DIST, LOCATION=NODE, TABLE=" << arrayName << "_TABLE" << std::endl;
        
        // Write default value
        double defaultValue = 0.0;
        file << ", " << defaultValue << std::endl;
        
        // Write node-specific values
        for (vtkIdType j = 0; j < array->GetNumberOfTuples(); j++) {
            // Skip nodes with zero values
            if (array->GetComponent(j, 0) == 0.0) continue;
            
            // Node IDs are 1-based in Abaqus
            vtkIdType nodeId = j + 1;
            
            // Write node ID and value
            file << nodeId << ", " << array->GetComponent(j, 0) << std::endl;
        }
    }
}

bool AbaqusInpReader::AppendFieldData(const std::string& inputFile, const std::string& outputFile, 
                      vtkSmartPointer<vtkUnstructuredGrid> mesh, 
                      const std::string& dataName, const std::string& targetLocation) {
    // First, copy the original file
    std::ifstream src(inputFile, std::ios::binary);
    std::ofstream dst(outputFile, std::ios::binary);
    
    if (!src.is_open() || !dst.is_open()) {
        std::cerr << "Failed to open input or output file" << std::endl;
        return false;
    }
    
    // Copy original content
    dst << src.rdbuf();
    
    // Close input file
    src.close();
    
    // Add a section marker
    dst << std::endl << "** Field data added by TensorMapper" << std::endl;
    
    // Check target location (elemental or nodal)
    if (targetLocation == "elemental") {
        // Get the data array from cell data
        vtkDataArray* array = mesh->GetCellData()->GetArray(dataName.c_str());
        if (!array) {
            std::cerr << "Cell data array not found: " << dataName << std::endl;
            dst.close();
            return false;
        }
        
        // Get original element IDs if available
        vtkIdTypeArray* originalIds = vtkIdTypeArray::SafeDownCast(
            mesh->GetCellData()->GetArray("OriginalElementIds"));
        
        // Create distribution table
        dst << "*DISTRIBUTION TABLE, NAME=" << dataName << "_TABLE" << std::endl;
        dst << "RATIO," << std::endl;
        
        // Create distribution
        dst << "*DISTRIBUTION, NAME=" << dataName << "_DIST, LOCATION=ELEMENT, TABLE=" << dataName << "_TABLE" << std::endl;
        
        // Write default value
        double defaultValue = 0.0;
        dst << ", " << defaultValue << std::endl;
        
        // Write element-specific values
        for (vtkIdType j = 0; j < array->GetNumberOfTuples(); j++) {
            // Skip elements with zero values
            if (array->GetComponent(j, 0) == 0.0) continue;
            
            // Get element ID
            vtkIdType elemId = j + 1; // Default to 1-based index
            if (originalIds && j < originalIds->GetNumberOfTuples()) {
                elemId = originalIds->GetValue(j);
            }
            
            // Write element ID and value
            dst << elemId << ", " << array->GetComponent(j, 0) << std::endl;
        }
    } 
    else if (targetLocation == "nodal") {
        // Get the data array from point data
        vtkDataArray* array = mesh->GetPointData()->GetArray(dataName.c_str());
        if (!array) {
            std::cerr << "Point data array not found: " << dataName << std::endl;
            dst.close();
            return false;
        }
        
        // Create distribution table
        dst << "*DISTRIBUTION TABLE, NAME=" << dataName << "_TABLE" << std::endl;
        dst << "RATIO," << std::endl;
        
        // Create distribution
        dst << "*DISTRIBUTION, NAME=" << dataName << "_DIST, LOCATION=NODE, TABLE=" << dataName << "_TABLE" << std::endl;
        
        // Write default value
        double defaultValue = 0.0;
        dst << ", " << defaultValue << std::endl;
        
        // Write node-specific values
        for (vtkIdType j = 0; j < array->GetNumberOfTuples(); j++) {
            // Skip nodes with zero values
            if (array->GetComponent(j, 0) == 0.0) continue;
            
            // Node IDs are 1-based in Abaqus
            vtkIdType nodeId = j + 1;
            
            // Write node ID and value
            dst << nodeId << ", " << array->GetComponent(j, 0) << std::endl;
        }
    }
    else {
        std::cerr << "Unknown target location: " << targetLocation << std::endl;
        dst.close();
        return false;
    }
    
    dst.close();
    std::cout << "Added " << dataName << " field data to " << outputFile << std::endl;
    return true;
}

std::vector<std::string> AbaqusInpReader::splitString(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string token;
    
    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    
    return tokens;
}

std::string AbaqusInpReader::trimString(const std::string& str) {
    auto start = str.begin();
    while (start != str.end() && std::isspace(*start)) {
        start++;
    }
    
    auto end = str.end();
    if (start != end) {
        do {
            end--;
        } while (end > start && std::isspace(*end));
    }
    
    return std::string(start, end + 1);
}

std::string AbaqusInpReader::toUpperString(const std::string& str) {
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(), ::toupper);
    return result;
}