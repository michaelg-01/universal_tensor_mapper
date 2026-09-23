#include "MappingMethods.h"

// Direct transfer method (nearest point)
void MappingMethods::mapDirect(vtkSmartPointer<vtkDataArray> sourceDataArray,
                             vtkSmartPointer<vtkDataArray> mappedDataArray,
                             vtkSmartPointer<vtkUnstructuredGrid> targetMesh,
                             vtkSmartPointer<vtkKdTreePointLocator> kdTree,
                             DataLocation sourceDataLocation,
                             DataLocation targetDataLocation) {
    vtkIdType numTargetItems = (targetDataLocation == DataLocation::NODAL) ?
        targetMesh->GetNumberOfPoints() : targetMesh->GetNumberOfCells();
    
    std::cout << "Direct mapping " << sourceDataArray->GetName() 
              << " to " << numTargetItems << " "
              << (targetDataLocation == DataLocation::NODAL ? "points" : "elements")
              << std::endl;
    
    if (targetDataLocation == DataLocation::NODAL) {
        // Mapping to nodal locations
        for (vtkIdType i = 0; i < numTargetItems; i++) {
            double point[3];
            targetMesh->GetPoint(i, point);
            
            // Find closest source point in KD-tree
            vtkIdType id = kdTree->FindClosestPoint(point);
            
            // For nodal-to-nodal, directly copy from source points
            mappedDataArray->SetTuple(i, sourceDataArray->GetTuple(id));
        }
    } else {
        // Mapping to elemental locations (target is elements)
        auto targetCentroids = calculateElementCentroids(targetMesh);
        
        for (vtkIdType i = 0; i < numTargetItems; i++) {
            double centroid[3];
            targetCentroids->GetPoint(i, centroid);
            
            // Find closest source point/element in KD-tree
            vtkIdType id = kdTree->FindClosestPoint(centroid);
            
            // Copy data to target element
            mappedDataArray->SetTuple(i, sourceDataArray->GetTuple(id));
        }
    }
    
    std::cout << "Mapping complete" << std::endl;
}

// Proximity-based mapping
void MappingMethods::mapProximity(vtkSmartPointer<vtkDataArray> sourceDataArray,
                              vtkSmartPointer<vtkDataArray> mappedDataArray,
                              vtkSmartPointer<vtkUnstructuredGrid> targetMesh,
                              vtkSmartPointer<vtkKdTreePointLocator> kdTree,
                              DataLocation sourceDataLocation,
                              DataLocation targetDataLocation,
                              double searchRadius) {
   
   vtkIdType numTargetItems = (targetDataLocation == DataLocation::NODAL) ?
       targetMesh->GetNumberOfPoints() : targetMesh->GetNumberOfCells();
   
   std::cout << "Proximity mapping " << sourceDataArray->GetName() 
             << " to " << numTargetItems << " "
             << (targetDataLocation == DataLocation::NODAL ? "points" : "elements")
             << " using radius " << searchRadius << std::endl;
   
   vtkNew<vtkIdList> idList;
   int numComponents = sourceDataArray->GetNumberOfComponents();
   std::vector<double> tempValues(numComponents);
   
   vtkSmartPointer<vtkPoints> targetPoints;
   if (targetDataLocation == DataLocation::NODAL) {
       targetPoints = targetMesh->GetPoints();
   } else {
       targetPoints = calculateElementCentroids(targetMesh);
   }
   
   int pointsWithNoNeighbors = 0;
   
   for (vtkIdType i = 0; i < numTargetItems; i++) {
       double point[3];
       targetPoints->GetPoint(i, point);
       
       kdTree->FindPointsWithinRadius(searchRadius, point, idList);
       
       if (idList->GetNumberOfIds() == 0) {
           pointsWithNoNeighbors++;
           continue; // Skip to next point (leave unchanged)
       }
       
       // Sort points by distance
       std::vector<std::pair<double, vtkIdType>> distanceIdPairs;
       for (vtkIdType n = 0; n < idList->GetNumberOfIds(); n++) {
           vtkIdType id = idList->GetId(n);
           double sourcePoint[3];
           kdTree->GetDataSet()->GetPoint(id, sourcePoint);
           double dist = sqrt(vtkMath::Distance2BetweenPoints(point, sourcePoint));
           distanceIdPairs.push_back(std::make_pair(dist, id));
       }
       
       // Sort by distance (ascending)
       std::sort(distanceIdPairs.begin(), distanceIdPairs.end());
       
       // Initialize weight calculation
       std::vector<double> weights(distanceIdPairs.size());
       double totalWeight = 0.0;
       
       for (int c = 0; c < numComponents; c++) {
           tempValues[c] = 0.0;
       }
       
       // Calculate weights using inverse square distance
       for (size_t n = 0; n < distanceIdPairs.size(); n++) {
           double dist = distanceIdPairs[n].first;
           
           // Set a minimum distance to prevent domination by exact matches
           if (dist < 0.01) dist = 0.01;
           
           // Use inverse square distance for higher contrast
           weights[n] = 1.0 / (dist * dist);
           
           totalWeight += weights[n];
       }
       
       // Calculate weighted average
       for (size_t n = 0; n < distanceIdPairs.size(); n++) {
           vtkIdType id = distanceIdPairs[n].second;
           double normalized_weight = weights[n] / totalWeight;
           
           for (int c = 0; c < numComponents; c++) {
               tempValues[c] += sourceDataArray->GetComponent(id, c) * normalized_weight;
           }
       }
       
       // Set final result
       mappedDataArray->SetTuple(i, tempValues.data());
   }
   
   if (pointsWithNoNeighbors > 0) {
       std::cout << "Warning: " << pointsWithNoNeighbors << " points had no neighbors within radius "
                 << searchRadius << ". These points were left unchanged." << std::endl;
   }
   
   std::cout << "Proximity mapping complete" << std::endl;
}

// N-point mapping 
void MappingMethods::mapNPoint(vtkSmartPointer<vtkDataArray> sourceDataArray, vtkSmartPointer<vtkDataArray> mappedDataArray,
                             vtkSmartPointer<vtkUnstructuredGrid> targetMesh, vtkSmartPointer<vtkKdTreePointLocator> kdTree,
                             DataLocation sourceDataLocation, DataLocation targetDataLocation, int numNeighbors) {
    
    vtkIdType numTargetItems = (targetDataLocation == DataLocation::NODAL) ?
        targetMesh->GetNumberOfPoints() : targetMesh->GetNumberOfCells();
    
    std::cout << "N-Point mapping with equal weights" << std::endl;
    
    // Create temporary structure for the nearest points search
    vtkNew<vtkIdList> idList;
    idList->SetNumberOfIds(numNeighbors);
    
    // Get number of components in data array
    int numComponents = sourceDataArray->GetNumberOfComponents();
    std::vector<double> tempValues(numComponents);
    
    // Prepare for different mapping types
    vtkSmartPointer<vtkPoints> targetPoints;
    
    if (targetDataLocation == DataLocation::NODAL) {
        targetPoints = targetMesh->GetPoints();
    } else {
        targetPoints = calculateElementCentroids(targetMesh);
    }
    
    // Process each target item
    for (vtkIdType i = 0; i < numTargetItems; i++) {
        double point[3];
        targetPoints->GetPoint(i, point);
        
        // Find N nearest source points in KD-tree
        kdTree->FindClosestNPoints(numNeighbors, point, idList);
        
        // Initialize temp values to zero
        for (int c = 0; c < numComponents; c++) {
            tempValues[c] = 0.0;
        }
        
        // Calculate simple average
        for (vtkIdType n = 0; n < idList->GetNumberOfIds(); n++) {
            vtkIdType id = idList->GetId(n);
            
            for (int c = 0; c < numComponents; c++) {
                tempValues[c] += sourceDataArray->GetComponent(id, c);
            }
        }
        
        for (int c = 0; c < numComponents; c++) {// Divide by number of points
            tempValues[c] /= idList->GetNumberOfIds();
        }
        mappedDataArray->SetTuple(i, tempValues.data()); // Set the result
    }
    
    std::cout << "N-Point mapping complete" << std::endl;
}

// Gaussian mapping
void MappingMethods::mapGaussian(vtkSmartPointer<vtkDataArray> sourceDataArray, 
                                vtkSmartPointer<vtkDataArray> mappedDataArray,
                                vtkSmartPointer<vtkUnstructuredGrid> targetMesh, 
                                vtkSmartPointer<vtkKdTreePointLocator> kdTree,
                                DataLocation sourceDataLocation, 
                                DataLocation targetDataLocation, 
                                double searchRadius, 
                                double sigma) {
    
    vtkIdType numTargetItems = (targetDataLocation == DataLocation::NODAL) ?
        targetMesh->GetNumberOfPoints() : targetMesh->GetNumberOfCells();
    
    std::cout << "Gaussian mapping " << sourceDataArray->GetName() 
              << " to " << numTargetItems << " "
              << (targetDataLocation == DataLocation::NODAL ? "points" : "elements")
              << " using radius " << searchRadius 
              << " and sigma " << sigma << std::endl;
    
    vtkNew<vtkIdList> idList;
    int numComponents = sourceDataArray->GetNumberOfComponents();
    std::vector<double> tempValues(numComponents);
    
    // Get target points (either mesh points or element centroids)
    vtkSmartPointer<vtkPoints> targetPoints;
    if (targetDataLocation == DataLocation::NODAL) {
        targetPoints = targetMesh->GetPoints();
    } else {
        targetPoints = calculateElementCentroids(targetMesh);
    }
    
    int pointsWithNoNeighbors = 0;
    // Pre-calculate sigma squared for efficiency
    double sigma2 = sigma * sigma;
    double twoSigma2 = 2.0 * sigma2;
    
    // If search radius is not specified, use 3*sigma as a reasonable cutoff
    // (at 3*sigma, the Gaussian weight is ~1% of peak)
    if (searchRadius <= 0) {
        searchRadius = 3.0 * sigma;
        std::cout << "Using automatic search radius: " << searchRadius 
                  << " (3 * sigma)" << std::endl;
    }
    
    // Process each target item
    for (vtkIdType i = 0; i < numTargetItems; i++) {
        double point[3];
        targetPoints->GetPoint(i, point);
        
        // Find all points within search radius
        kdTree->FindPointsWithinRadius(searchRadius, point, idList);
        
        if (idList->GetNumberOfIds() == 0) {
            // No neighbors found - try finding at least the closest point
            vtkIdType closestId = kdTree->FindClosestPoint(point);
            
            // Get distance to closest point
            double closestPoint[3];
            kdTree->GetDataSet()->GetPoint(closestId, closestPoint);
            double dist = sqrt(vtkMath::Distance2BetweenPoints(point, closestPoint));
            
            // If closest point is within extended radius (5*sigma), use it
            if (dist < 5.0 * sigma) {
                idList->InsertNextId(closestId);
            } else {
                // Too far away, use default value or keep unchanged
                pointsWithNoNeighbors++;
                // Initialize with zeros or copy from closest point as fallback
                for (int c = 0; c < numComponents; c++) {
                    tempValues[c] = sourceDataArray->GetComponent(closestId, c);
                }
                mappedDataArray->SetTuple(i, tempValues.data());
                continue;
            }
        }
        
        // Calculate Gaussian weights for all neighbors
        std::vector<std::pair<double, vtkIdType>> distanceIdPairs;
        std::vector<double> weights;
        double totalWeight = 0.0;
        
        // Initialize components to zero
        for (int c = 0; c < numComponents; c++) {
            tempValues[c] = 0.0;
        }
        
        // Calculate distances and weights
        for (vtkIdType n = 0; n < idList->GetNumberOfIds(); n++) {
            vtkIdType id = idList->GetId(n);
            double sourcePoint[3];
            kdTree->GetDataSet()->GetPoint(id, sourcePoint);
            
            // Calculate distance
            double dist = sqrt(vtkMath::Distance2BetweenPoints(point, sourcePoint));
            
            // Calculate Gaussian weight: exp(-(dist^2) / (2*sigma^2))
            double weight = exp(-(dist * dist) / twoSigma2);
            
            // Store for weighted averaging
            weights.push_back(weight);
            totalWeight += weight;
            
            // Store for potential debugging
            distanceIdPairs.push_back(std::make_pair(dist, id));
        }
        
        // Normalize weights and calculate weighted average
        if (totalWeight > 0.0) {
            for (size_t n = 0; n < weights.size(); n++) {
                vtkIdType id = idList->GetId(n);
                double normalizedWeight = weights[n] / totalWeight;
                
                // Add weighted contribution of this neighbor
                for (int c = 0; c < numComponents; c++) {
                    tempValues[c] += sourceDataArray->GetComponent(id, c) * normalizedWeight;
                }
            }
        }
        
        // Set the result
        mappedDataArray->SetTuple(i, tempValues.data());
    }
    
    if (pointsWithNoNeighbors > 0) {
        std::cout << "Warning: " << pointsWithNoNeighbors 
                  << " points had no suitable neighbors and used fallback values" << std::endl;
    }
    
    std::cout << "Gaussian mapping complete" << std::endl;
}

void MappingMethods::applyGlobalSmoothing(
    vtkSmartPointer<vtkDataArray> mappedDataArray,
    vtkSmartPointer<vtkUnstructuredGrid> targetMesh,
    DataLocation targetDataLocation,
    int smoothingIterations,
    double relaxationFactor)
{
    if (!mappedDataArray || smoothingIterations <= 0)
        return;
    
    std::cout << "Applying Laplacian smoothing with " << smoothingIterations 
              << " iterations, relaxation factor: " << relaxationFactor << std::endl;
    
    vtkIdType numItems = mappedDataArray->GetNumberOfTuples();
    int numComponents = mappedDataArray->GetNumberOfComponents();
    
    // Create a temporary array for the smoothing process
    auto tempArray = vtkSmartPointer<vtkDataArray>::Take(mappedDataArray->NewInstance());
    tempArray->SetNumberOfComponents(numComponents);
    tempArray->SetNumberOfTuples(numItems);
    
    // Create a list of neighbors for each item (point or cell)
    std::vector<std::vector<vtkIdType>> neighbors(numItems);
    
    // Build neighbor lists
    if (targetDataLocation == DataLocation::NODAL) {
        // For nodal data, find point neighbors via cells
        for (vtkIdType cellId = 0; cellId < targetMesh->GetNumberOfCells(); cellId++) {
            vtkCell* cell = targetMesh->GetCell(cellId);
            vtkIdType numPoints = cell->GetNumberOfPoints();
            
            for (vtkIdType i = 0; i < numPoints; i++) {
                vtkIdType pointId1 = cell->GetPointId(i);
                
                for (vtkIdType j = 0; j < numPoints; j++) {
                    if (i != j) {
                        vtkIdType pointId2 = cell->GetPointId(j);
                        neighbors[pointId1].push_back(pointId2);
                    }
                }
            }
        }
    } else {
        // For cell data, find cells that share at least one node
        std::vector<std::vector<vtkIdType>> pointToCells(targetMesh->GetNumberOfPoints());
        
        for (vtkIdType cellId = 0; cellId < targetMesh->GetNumberOfCells(); cellId++) {
            vtkCell* cell = targetMesh->GetCell(cellId);
            vtkIdType numPoints = cell->GetNumberOfPoints();
            
            for (vtkIdType i = 0; i < numPoints; i++) {
                vtkIdType pointId = cell->GetPointId(i);
                pointToCells[pointId].push_back(cellId);
            }
        }
        
        // Two cells are neighbors if they share at least one point
        for (vtkIdType pointId = 0; pointId < targetMesh->GetNumberOfPoints(); pointId++) {
            const auto& cells = pointToCells[pointId];
            
            for (size_t i = 0; i < cells.size(); i++) {
                for (size_t j = 0; j < cells.size(); j++) {
                    if (i != j) {
                        neighbors[cells[i]].push_back(cells[j]);
                    }
                }
            }
        }
    }
    
    // Remove duplicate neighbors
    for (vtkIdType i = 0; i < numItems; i++) {
        std::sort(neighbors[i].begin(), neighbors[i].end());
        auto last = std::unique(neighbors[i].begin(), neighbors[i].end());
        neighbors[i].erase(last, neighbors[i].end());
    }
    
    // Perform iterations
    for (int iter = 0; iter < smoothingIterations; iter++) {
        // Copy current values to temp array
        tempArray->DeepCopy(mappedDataArray);
        
        // For each item
        for (vtkIdType itemId = 0; itemId < numItems; itemId++) {
            const auto& itemNeighbors = neighbors[itemId];
            
            // Skip if no neighbors
            if (itemNeighbors.empty())
                continue;
            
            // For each component
            for (int comp = 0; comp < numComponents; comp++) {
                double originalValue = tempArray->GetComponent(itemId, comp);
                double sum = 0.0;
                
                // Calculate average of neighbors
                for (vtkIdType neighborId : itemNeighbors) {
                    sum += tempArray->GetComponent(neighborId, comp);
                }
                
                double average = sum / itemNeighbors.size();
                
                // Apply relaxation factor using standard Laplacian formula
                double newValue = originalValue + relaxationFactor * (average - originalValue);
                mappedDataArray->SetComponent(itemId, comp, newValue);
            }
        }
    }
    
    std::cout << "Laplacian smoothing complete" << std::endl;
}

// Calculate element centroids
vtkSmartPointer<vtkPoints> MappingMethods::calculateElementCentroids(vtkSmartPointer<vtkUnstructuredGrid> mesh) {
    vtkIdType numCells = mesh->GetNumberOfCells();
    auto centroids = vtkSmartPointer<vtkPoints>::New();
    centroids->SetNumberOfPoints(numCells);
    
    for (vtkIdType i = 0; i < numCells; i++) {
        double centroid[3] = {0, 0, 0};
        auto cell = mesh->GetCell(i);
        vtkIdType numPoints = cell->GetNumberOfPoints();
        
        // Average all points in the cell
        for (vtkIdType j = 0; j < numPoints; j++) {
            double point[3];
            mesh->GetPoint(cell->GetPointId(j), point);
            centroid[0] += point[0];
            centroid[1] += point[1];
            centroid[2] += point[2];
        }
        
        centroid[0] /= numPoints;
        centroid[1] /= numPoints;
        centroid[2] /= numPoints;
        
        centroids->SetPoint(i, centroid);
    }
    
    return centroids;
}