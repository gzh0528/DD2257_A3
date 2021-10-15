/*********************************************************************
 *  Author  : Anke Friederici
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 **********************************************************************/

#include <inviwo/core/datastructures/geometry/basicmesh.h>
#include <inviwo/core/datastructures/volume/volumeram.h>
#include <labstreamlines/integrator.h>
#include <labutils/scalarvectorfield.h>
#include <labtopo/topology.h>
#include <labtopo/utils/gradients.h>

namespace inviwo {

const vec4 Topology::ColorsCP[6] = {
    vec4(1, 1, 0, 1),    // Saddle - Yellow
    vec4(1, 0, 0, 1),    // AttractingNode - Red
    vec4(0, 0, 1, 1),    // RepellingNode - Blue
    vec4(0.5, 0, 1, 1),  // AttractingFocus - Purple
    vec4(1, 0.5, 0, 1),  // RepellingFocus - Orange
    vec4(0, 1, 0, 1)     // Center - Green
};

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo Topology::processorInfo_{
    "org.inviwo.Topology",    // Class identifier
    "Vector Field Topology",  // Display name
    "KTH Lab",                // Category
    CodeState::Experimental,  // Code state
    Tags::None,               // Tags
};

const ProcessorInfo Topology::getProcessorInfo() const { return processorInfo_; }

Topology::Topology()
    : Processor()
    , inData("inData")
    , outMesh("meshOut")
    , meshBBoxOut("meshBBoxOut")
// TODO: Initialize additional properties
// propertyName("propertyIdentifier", "Display Name of the Propery",
// default value (optional), minimum value (optional), maximum value (optional), increment
// (optional)); propertyIdentifier cannot have spaces
{
    // Register Ports
    addPort(outMesh);
    addPort(inData);
    addPort(meshBBoxOut);

    // TODO: Register additional properties
    // addProperty(propertyName);
}

static void streamline(const VectorField2& vectorField, const dvec2& seed, double stepSize, int numSteps, std::vector<dvec2>& points)
{
    points.clear();

    dvec2 x0 = seed;
    for (int i = 0; i < numSteps; i++) {
        dvec2 x1 = Integrator::RK4_norm(vectorField, x0, stepSize);

        if (!vectorField.isInside(x1))
            break;

        if (vectorField.interpolate(x1) == dvec2(0))
            break;

        points.push_back(x1);

        x0 = x1;
    }
}

void findBoundPoint(const VectorField2& vectorField, double thresholdLength, dvec2 lo, dvec2 hi, vector<dvec2>& BoundaryPoint,double minThreshold,int sign) {


    dvec2 diff = hi-lo;
    dvec2 mid=(hi+lo)/2.0;
    if (max(diff.x,diff.y)<=thresholdLength) {
        if(((sign==1)&&(vectorField.interpolate(lo).x*vectorField.interpolate(hi).x<0))||((sign==0)&&(vectorField.interpolate(lo).y*vectorField.interpolate(hi).y<0)))
        {
            BoundaryPoint.push_back(mid);
        }
        return;
    }
    findBoundPoint(vectorField, thresholdLength,lo, lo+diff/2.0,BoundaryPoint,minThreshold,sign);
    findBoundPoint(vectorField, thresholdLength, lo+diff/2.0, hi,BoundaryPoint,minThreshold,sign);
}

// Linear cell has at most 1 critical point
bool findCriticalPoint(const VectorField2& vectorField, double thresholdArea, dvec2 lo, dvec2 hi, dvec2& criticalPoint) {
    dvec2 values[4];
    values[0] = vectorField.interpolate({lo.x, lo.y});
    values[1] = vectorField.interpolate({lo.x, hi.y});
    values[2] = vectorField.interpolate({hi.x, hi.y});
    values[3] = vectorField.interpolate({hi.x, lo.y});

    for (int dim = 0; dim < 2; dim++) {
        int sgn = std::signbit(values[0][dim]);
        for (int i = 1; i < 4; i++)
            if (std::signbit(values[i][dim]) != sgn)
                goto next_dim;
        return false;
        next_dim:;
    }

    dvec2 diff = hi-lo;
    double area = diff.x*diff.y;
    if (area > thresholdArea) {
        double dx = diff.x/2;
        double dy = diff.y/2;

        dvec2 subcells[4][2] = {
            {lo + dvec2(0, 0), hi - dvec2(dx, dy)},
            {lo + dvec2(0, dy), hi - dvec2(dx, 0)},
            {lo + dvec2(dx, 0), hi - dvec2(0, dy)},
            {lo + dvec2(dx, dy), hi - dvec2(0, 0)}
        };

        bool found = false;
        for (int i = 0; i < 4 && !found; i++)
            found = findCriticalPoint(vectorField, thresholdArea, subcells[i][0], subcells[i][1], criticalPoint);
        return found;
    } else {
        criticalPoint = (lo+hi)/2;
        return true;
    }
}

Topology::TypeCP categorizePoint(const util::EigenResult& eigen) {
    double threshold = 1e-4;
    if (eigen.eigenvaluesIm[0] == 0.0) {
        if (eigen.eigenvaluesRe[0] * eigen.eigenvaluesRe[1] < 0)
            return Topology::TypeCP::Saddle;
        else if (eigen.eigenvaluesRe[0] < 0)
            return Topology::TypeCP::AttractingNode;
        else
            return Topology::TypeCP::RepellingNode;
    } else {
        if (eigen.eigenvaluesRe[0] < -threshold)
            return Topology::TypeCP::AttractingFocus;
        else if (eigen.eigenvaluesRe[0] > threshold)
            return Topology::TypeCP::RepellingFocus;
        else
            return Topology::TypeCP::Center;
    }
}

void Topology::process() {
    // Get input
    if (!inData.hasData()) {
        return;
    }
    auto vol = inData.getData();

    // Retreive data in a form that we can access it
    const VectorField2 vectorField = VectorField2::createFieldFromVolume(vol);

    // Add a bounding box to the mesh
    const dvec2& BBoxMin = vectorField.getBBoxMin();
    const dvec2& BBoxMax = vectorField.getBBoxMax();
    auto bboxMesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> bboxVertices;
    auto indexBufferBBox = bboxMesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    // Bounding Box vertex 0
    vec4 black = vec4(0, 0, 0, 1);
    Integrator::drawNextPointInPolyline(BBoxMin, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMin[0], BBoxMax[1]), black, indexBufferBBox.get(),
                                        bboxVertices);
    Integrator::drawNextPointInPolyline(BBoxMax, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMax[0], BBoxMin[1]), black, indexBufferBBox.get(),
                                        bboxVertices);
    // Connect back to the first point, to make a full rectangle
    indexBufferBBox->add(static_cast<std::uint32_t>(0));
    bboxMesh->addVertices(bboxVertices);
    meshBBoxOut.setData(bboxMesh);

    // Initialize mesh, vertices and index buffers for seperatrices
    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;
    // Either add all line segments to this index buffer (one large buffer, two consecutive points
    // make up one line), or use several index buffers with connectivity type strip.
    auto indexBufferSeparatrices = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
    // auto indexBufferSeparatrices = mesh->addIndexBuffer(DrawType::Lines,
    // ConnectivityType::Strip);

    auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);

    size2_t dims = vectorField.getNumVerticesPerDim();

    dvec2 cellSize = vectorField.getCellSize();
    double thresholdArea = 1e-6;
    double thresholdLength=1e-5;
    double minThreshold=1e-5;
    struct critical_point {
        dvec2 point;
        size2_t vert;
    };
    std::vector<critical_point> criticalPoints;
    //std::vector<std::vector<bool>> occupied(dims.x, std::vector<bool>(dims.y));
    std::vector<dvec2> BoundaryPoint;
    for (size_t j = 0; j < dims[1]; ++j) {
        for (size_t i = 0; i < dims[0]; ++i) {
            dvec2 point;
            dvec2 lo = vectorField.getPositionAtVertex({i, j});
            dvec2 hi = lo + cellSize;
            if (findCriticalPoint(vectorField, thresholdArea, lo, hi, point)) {
                LogProcessorInfo("Found " << point << " value " << vectorField.interpolate(point));
                criticalPoints.push_back({point, {i, j}});
                //occupied[i][j] = true;
            }
            if(i==0)
            {
                dvec2 hi1(lo.x,lo.y+cellSize.y);
                findBoundPoint(vectorField, thresholdLength,lo, hi1,BoundaryPoint,minThreshold,1);
            }
            if(j==0)
            {
                dvec2 hi2(lo.x+cellSize.x,lo.y);
                findBoundPoint(vectorField, thresholdLength,lo, hi2,BoundaryPoint,minThreshold,0);
            }
            if(i==dims[0]-1)
            {
                dvec2 hi3(lo.x,lo.y+cellSize.y);
                findBoundPoint(vectorField, thresholdLength,lo, hi3,BoundaryPoint,minThreshold,1);
            }
            if(j==dims[1]-1)
            {
                dvec2 hi4(lo.x+cellSize.x,lo.y);
                findBoundPoint(vectorField, thresholdLength,lo, hi4,BoundaryPoint,minThreshold,0);
            }
        }
    }

    /*
    std::vector<std::vector<bool>> hasNeighbor(dims.x, std::vector<bool>(dims.y));
    for (size_t j = 0; j < dims.y; j++) {
        bool last = occupied[0][j];
        for (size_t i = 1; i < dims.x; i++) {
            bool curr = occupied[i][j];
            if (curr && last) {
                hasNeighbor[i-1][j] = true;
                hasNeighbor[i][j] = true;
            }
            last = curr;
        }
    }
    for (size_t i = 0; i < dims.x; i++) {
        bool last = occupied[i][0];
        for (size_t j = 1; j < dims.y; j++) {
            bool curr = occupied[i][j];
            if (curr && last) {
                hasNeighbor[i][j-1] = true;
                hasNeighbor[i][j] = true;
            }
            last = curr;
        }
    }
    */

    int maxSteps = 1000;
    dvec2 extent = BBoxMax-BBoxMin;
    double minExtent = std::min(extent.x, extent.y);
    double stepSize = minExtent / maxSteps;
    std::vector<dvec2> tmpBPoints;
    for(auto p:BoundaryPoint)
    {
        vec4 color_(0.5,0.5,0.5,1);
        mat2 jacobian = vectorField.derive(p);
        auto eigen = util::eigenAnalysis(jacobian);
        for (int eigenValue = 0; eigenValue < 2; eigenValue++) {
            double myStepSize = eigen.eigenvaluesRe[eigenValue] > 0 ? stepSize : -stepSize;
            dvec2 step = stepSize*eigen.eigenvectors[eigenValue];

            streamline(vectorField, p + step, myStepSize, maxSteps, tmpBPoints);
            drawSeparatrix(tmpBPoints, vertices, *mesh);

            streamline(vectorField, p - step, myStepSize, maxSteps, tmpBPoints);
            drawSeparatrix(tmpBPoints, vertices, *mesh);
        }
        vertices.push_back({vec3(p[0], p[1], 0), vec3(0, 0, 1), vec3(p[0], p[1], 0), color_});
        Integrator::drawPoint(p, vec4(0.5, 0.5, 0.5, 1), indexBufferPoints.get(), vertices);
    }
    std::vector<dvec2> tmpPoints;
    for (size_t i = 0; i < criticalPoints.size(); i++) {
        dvec2 p = criticalPoints[i].point;
        size2_t vert = criticalPoints[i].vert;
        mat2 jacobian = vectorField.derive(p);

        //bool isIsolated = !hasNeighbor[vert.x][vert.y];
        bool isIsolated = true;
        for (size_t j = 0; j < criticalPoints.size() && isIsolated; j++) {
            if (i == j)
                continue;
            if (glm::length(p - criticalPoints[j].point) < minExtent / 50)
                isIsolated = false;
        }
        bool firstOrder = isIsolated && glm::abs(glm::determinant(jacobian)) > 1e-5;

        if (!firstOrder)
            continue;

        auto eigen = util::eigenAnalysis(jacobian);
        TypeCP type = categorizePoint(eigen);

        if (type == TypeCP::Saddle) {
            for (int eigenValue = 0; eigenValue < 2; eigenValue++) {
                double myStepSize = eigen.eigenvaluesRe[eigenValue] > 0 ? stepSize : -stepSize;
                dvec2 step = stepSize*eigen.eigenvectors[eigenValue];

                streamline(vectorField, p + step, myStepSize, maxSteps, tmpPoints);
                drawSeparatrix(tmpPoints, vertices, *mesh);

                streamline(vectorField, p - step, myStepSize, maxSteps, tmpPoints);
                drawSeparatrix(tmpPoints, vertices, *mesh);
            }
        }

        vec4 color = ColorsCP[static_cast<int>(type)];
        indexBufferPoints->add(static_cast<std::uint32_t>(vertices.size()));
        vertices.push_back({vec3(p[0], p[1], 0), vec3(0, 0, 1), vec3(p[0], p[1], 0), color});
    }

    mesh->addVertices(vertices);
    outMesh.setData(mesh);
}

void Topology::drawLineSegment(const dvec2& v1, const dvec2& v2, const vec4& color,
                               IndexBufferRAM* indexBuffer,
                               std::vector<BasicMesh::Vertex>& vertices) {
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v1[0], v1[1], 0), vec3(0, 0, 1), vec3(v1[0], v1[1], 0), color});
    indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
    vertices.push_back({vec3(v2[0], v2[1], 0), vec3(0, 0, 1), vec3(v2[0], v2[1], 0), color});
}

void Topology::drawSeparatrix(const std::vector<dvec2>& points,
                              std::vector<BasicMesh::Vertex>& vertices,
                              BasicMesh& mesh) {
    vec4 color(1);
    auto indexBuffer = mesh.addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    for (size_t i = 0; i < points.size(); i++) {
        indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
        vertices.push_back({vec3(points[i][0], points[i][1], 0), vec3(0, 0, 1), vec3(points[i][0], points[i][1], 0), color});
    }
}

}  // namespace inviwo
