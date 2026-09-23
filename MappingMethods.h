#pragma once

#include <iostream>
#include <vector>
#include <set>
#include <cmath>

// VTK includes
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkCell.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkKdTreePointLocator.h>
#include <vtkIdTypeArray.h>
#include <vtkIdList.h>
#include <vtkMath.h>
#include <vtkNew.h>
#include <vtkAttributeSmoothingFilter.h>

// Define data location enum
enum class DataLocation {
    NODAL,
    ELEMENTAL
};

class MappingMethods {
public:
    // Direct transfer method (nearest point)
    static void mapDirect(vtkSmartPointer<vtkDataArray> sourceDataArray,
                       vtkSmartPointer<vtkDataArray> mappedDataArray,
                       vtkSmartPointer<vtkUnstructuredGrid> targetMesh,
                       vtkSmartPointer<vtkKdTreePointLocator> kdTree,
                       DataLocation sourceDataLocation,
                       DataLocation targetDataLocation);
    
    // Proximity-based mapping
    static void mapProximity(vtkSmartPointer<vtkDataArray> sourceDataArray,
                          vtkSmartPointer<vtkDataArray> mappedDataArray,
                          vtkSmartPointer<vtkUnstructuredGrid> targetMesh,
                          vtkSmartPointer<vtkKdTreePointLocator> kdTree,
                          DataLocation sourceDataLocation,
                          DataLocation targetDataLocation,
                          double searchRadius);
    
    // N-point mapping
    static void mapNPoint(vtkSmartPointer<vtkDataArray> sourceDataArray,
                       vtkSmartPointer<vtkDataArray> mappedDataArray,
                       vtkSmartPointer<vtkUnstructuredGrid> targetMesh,
                       vtkSmartPointer<vtkKdTreePointLocator> kdTree,
                       DataLocation sourceDataLocation,
                       DataLocation targetDataLocation,
                       int numNeighbors);
    
    // Gaussian mapping
    static void mapGaussian(vtkSmartPointer<vtkDataArray> sourceDataArray,
                         vtkSmartPointer<vtkDataArray> mappedDataArray,
                         vtkSmartPointer<vtkUnstructuredGrid> targetMesh,
                         vtkSmartPointer<vtkKdTreePointLocator> kdTree,
                         DataLocation sourceDataLocation,
                         DataLocation targetDataLocation,
                         double radius,
                         double sigma);
    
    // Apply global smoothing
   static void applyGlobalSmoothing(vtkSmartPointer<vtkDataArray> mappedDataArray,
                              vtkSmartPointer<vtkUnstructuredGrid> targetMesh,
                              DataLocation targetDataLocation,
                              int smoothingIterations,
                              double relaxationFactor);
    
    // Calculate element centroids
    static vtkSmartPointer<vtkPoints> calculateElementCentroids(vtkSmartPointer<vtkUnstructuredGrid> mesh);
};