// AbaqusInpReader.h
#pragma once

#include <string>
#include <vector>
#include <map>
#include <set>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkIdTypeArray.h>

class AbaqusInpReader {
public:
    // Read an Abaqus INP file and return a VTK unstructured grid
    static vtkSmartPointer<vtkUnstructuredGrid> read(const std::string& filename);
    
    // Write a VTK unstructured grid to an Abaqus INP file
    static bool write(const std::string& filename, vtkSmartPointer<vtkUnstructuredGrid> mesh);

    static bool AppendFieldData(const std::string& inputFile, const std::string& outputFile, 
        vtkSmartPointer<vtkUnstructuredGrid> mesh, 
        const std::string& dataName, const std::string& targetLocation);

private:
    // Helper methods for parsing INP files
    static void parseNodeLine(const std::string& line, vtkSmartPointer<vtkPoints> points, 
                       std::map<int, vtkIdType>& nodeIdMap);
    
    static bool parseElementLine(const std::string& line, int& elementId, std::vector<int>& nodeIds);
    
    static void createCellsFromElements(vtkSmartPointer<vtkUnstructuredGrid> grid,
                             const std::map<std::string, std::vector<std::pair<int, std::vector<int>>>>& elementsByType,
                             const std::map<int, vtkIdType>& nodeIdMap,
                             vtkIdTypeArray* originalIds);
    
    static int mapAbaqusElementTypeToVtk(const std::string& abaqusType, int numNodes);
    
    static std::string mapVtkElementTypeToAbaqus(int vtkType);
    
    static void writeFieldData(std::ofstream& file, 
        vtkSmartPointer<vtkUnstructuredGrid> mesh,
        const std::map<std::string, std::set<vtkIdType>>& elementSets);

    static void writePointFieldData(std::ofstream& file, vtkSmartPointer<vtkUnstructuredGrid> mesh);
    
    // String manipulation utilities
    static std::vector<std::string> splitString(const std::string& str, char delimiter);
    
    static std::string trimString(const std::string& str);
    
    static std::string toUpperString(const std::string& str);
};