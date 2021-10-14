/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Tuesday, September 19, 2017 - 15:08:33
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <inviwo/core/interaction/events/mouseevent.h>
#include <inviwo/core/util/utilities.h>
#include <labstreamlines/integrator.h>
#include <labstreamlines/streamlineintegrator.h>
#include <labutils/scalarvectorfield.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming
// scheme
const ProcessorInfo StreamlineIntegrator::processorInfo_{
    "org.inviwo.StreamlineIntegrator",  // Class identifier
    "Streamline Integrator",            // Display name
    "KTH Lab",                          // Category
    CodeState::Experimental,            // Code state
    Tags::None,                         // Tags
};

const ProcessorInfo StreamlineIntegrator::getProcessorInfo() const { return processorInfo_; }

StreamlineIntegrator::StreamlineIntegrator()
    : Processor()
    , inData("volIn")
    , meshOut("meshOut")
    , meshBBoxOut("meshBBoxOut")
    , propStartPoint("startPoint", "Start Point", vec2(0.5, 0.5), vec2(-1), vec2(1), vec2(0.1))
    , propSeedMode("seedMode2", "Seeds")
    , propRandomCount("randomCount", "Number of seed points", 10, 0)
    , propGridPoints("gridPoints", "Number of grid points", ivec2(4, 4))
    , propDisplayPoints("displayPoints", "Display Points", true)
    , propNumStepsTaken("numstepstaken", "Number of actual steps", 0, 0, 100000)
    , mouseMoveStart(
          "mouseMoveStart", "Move Start", [this](Event* e) { eventMoveStart(e); },
          MouseButton::Left, MouseState::Press | MouseState::Move)
// TODO: Initialize additional properties
// propertyName("propertyIdentifier", "Display Name of the Propery",
// default value (optional), minimum value (optional), maximum value (optional),
// increment (optional)); propertyIdentifier cannot have spaces
    , propBackwards("propBackwards", "Integrate backwards")
    , propStepSize("propStepSize", "Step size", 0.1, 0.001, 10)
    , propNormalize("propNormalize", "Normalize")
    , propMaxSteps("propMaxSteps", "Max steps", 100, 0, 10000)
    , propMaxArc("propMaxArc", "Max arc length", 1, 0.1, 100)
    , propMinVelocity("propMinVelocity", "Min velocity", 0, 0, 10)
{
    // Register Ports
    addPort(inData);
    addPort(meshOut);
    addPort(meshBBoxOut);

    // Register Properties
    propSeedMode.addOption("one", "Single Start Point", 0);
    propSeedMode.addOption("random", "Random Seeds", 1);
    propSeedMode.addOption("uniform", "Uniform Seeds", 2);
    addProperty(propSeedMode);
    addProperty(propRandomCount);
    addProperty(propGridPoints);
    addProperty(propStartPoint);
    addProperty(propDisplayPoints);
    addProperty(propNumStepsTaken);
    propNumStepsTaken.setReadOnly(true);
    propNumStepsTaken.setSemantics(PropertySemantics::Text);
    addProperty(mouseMoveStart);

    // TODO: Register additional properties
    // addProperty(propertyName);
    addProperty(propBackwards);
    addProperty(propStepSize);
    addProperty(propNormalize);
    addProperty(propMaxSteps);
    addProperty(propMaxArc);
    addProperty(propMinVelocity);

    // Show properties for a single seed and hide properties for multiple seeds
    // (TODO)
    propSeedMode.onChange([this]() {
        if (propSeedMode.get() == 0) {
            util::show(propStartPoint, mouseMoveStart);
            util::hide(propRandomCount, propGridPoints);
        } else if (propSeedMode.get() == 1) {
            util::hide(propStartPoint, mouseMoveStart, propGridPoints);
            util::show(propRandomCount);
        } else {
            util::hide(propStartPoint, mouseMoveStart, propRandomCount);
            util::show(propGridPoints);
        }
    });
}

void StreamlineIntegrator::eventMoveStart(Event* event) {
    if (!inData.hasData()) return;
    auto mouseEvent = static_cast<MouseEvent*>(event);
    vec2 mousePos = mouseEvent->posNormalized();

    // Map to bounding box range
    mousePos[0] *= static_cast<float>(BBoxMax_[0] - BBoxMin_[0]);
    mousePos[1] *= static_cast<float>(BBoxMax_[1] - BBoxMin_[1]);
    mousePos += static_cast<vec2>(BBoxMin_);

    // Update starting point
    propStartPoint.set(mousePos);
    event->markAsUsed();
}

int StreamlineIntegrator::drawStreamline(
        const VectorField2& vectorField,
        vec2 startPoint,
        std::shared_ptr<IndexBufferRAM> indexBufferPoints,
        std::shared_ptr<BasicMesh> mesh,
        std::vector<BasicMesh::Vertex>& vertices)
{
    // Draw start point
    if (propDisplayPoints.get() != 0)
        Integrator::drawPoint(startPoint, vec4(0, 0, 0, 1), indexBufferPoints.get(), vertices);

    // TODO: Create one stream line from the given start point
    auto indexBufferLine = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    vec4 color(0, 0, 1, 1);
    Integrator::drawNextPointInPolyline(startPoint, color, indexBufferLine.get(), vertices);

    dvec2 x0 = startPoint;
    int stepsTaken = 0;
    double arcLength = 0;
    for (;;) {
        if (stepsTaken >= propMaxSteps.get()) {
            //LogProcessorInfo("Stop: max steps");
            break;
        }

        dvec2 x1 = Integrator::RK4(vectorField, x0, propStepSize.get());

        if (!vectorField.isInside(x1)) {
            //LogProcessorInfo("Stop: domain");
            break;
        }

        if (vectorField.interpolate(x1) == dvec2(0)) {
            //LogProcessorInfo("Stop: zero");
            break;
        }

        /*
        auto jacobian = vectorField.derive(x1);
        LogProcessorInfo("det: " << glm::abs(glm::determinant(jacobian)));
        if (glm::abs(glm::determinant(jacobian)) < propMinVelocity.get()) {
            //LogProcessorInfo("Stop: min velocity");
            break;
        }
        */
        double dist = glm::distance(x0, x1);
        if (dist < propMinVelocity.get()) {
            //LogProcessorInfo("Stop: min velocity");
            break;
        }

        arcLength += dist;
        if (arcLength >= propMaxArc.get()) {
            //LogProcessorInfo("Stop: max arc length");
            break;
        }

        Integrator::drawNextPointInPolyline(x1, color, indexBufferLine.get(), vertices);
        if (propDisplayPoints.get())
            Integrator::drawNextPointInPolyline(x1, color, indexBufferPoints.get(), vertices);
        x0 = x1;
        stepsTaken++;
    }

    return stepsTaken;
}

dvec2 StreamlineIntegrator::randomPoint(const VectorField2& vectorField)
{
    for (;;) {
        dvec2 x = {xdistr(randGenerator), ydistr(randGenerator)};
        double mag = glm::length(vectorField.interpolate(x));
        if (mag >= magdistr(randGenerator))
            return x;
    }
}

void StreamlineIntegrator::process() {
    // Get input
    if (!inData.hasData()) {
        return;
    }
    auto vol = inData.getData();

    // Retreive data in a form that we can access it
    auto vectorField = VectorField2::createFieldFromVolume(vol);
    BBoxMin_ = vectorField.getBBoxMin();
    BBoxMax_ = vectorField.getBBoxMax();

    xdistr = std::uniform_real_distribution(BBoxMin_.x, BBoxMax_.x);
    ydistr = std::uniform_real_distribution(BBoxMin_.y, BBoxMax_.y);
    magdistr = std::uniform_real_distribution(
            glm::length(vectorField.getMinValue()),
            glm::length(vectorField.getMaxValue()));

    auto bboxMesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> bboxVertices;

    // Make bounding box without vertex duplication, instead of line segments which duplicate
    // vertices, create line segments between each added points with connectivity type of the index
    // buffer
    auto indexBufferBBox = bboxMesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
    // Bounding Box vertex 0
    vec4 black = vec4(0, 0, 0, 1);
    Integrator::drawNextPointInPolyline(BBoxMin_, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMin_[0], BBoxMax_[1]), black,
                                        indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(BBoxMax_, black, indexBufferBBox.get(), bboxVertices);
    Integrator::drawNextPointInPolyline(vec2(BBoxMax_[0], BBoxMin_[1]), black,
                                        indexBufferBBox.get(), bboxVertices);
    // Connect back to the first point, to make a full rectangle
    indexBufferBBox->add(static_cast<std::uint32_t>(0));
    bboxMesh->addVertices(bboxVertices);
    meshBBoxOut.setData(bboxMesh);

    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;

    if (propBackwards.get() || propNormalize.get()) {
        const ivec2 nVertPerDim = vectorField.getNumVerticesPerDim();
        VectorField2 newField(nVertPerDim, BBoxMin_, BBoxMax_ - BBoxMin_);

        for (int i = 0; i < nVertPerDim[0]; i++)
            for (int j = 0; j < nVertPerDim[1]; j++) {
                auto v = vectorField.getValueAtVertex({i, j});
                if (propBackwards.get())
                    v = -v;
                if (propNormalize.get() && v != dvec2(0))
                    v = glm::normalize(v);
                newField.setValueAtVertex({i, j}, v);
            }
        vectorField = newField;
    }

    auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
    int stepsTaken = 0;

    if (propSeedMode.get() == 0) {
        vec2 startPoint = propStartPoint.get();
        stepsTaken = drawStreamline(vectorField, startPoint, indexBufferPoints, mesh, vertices);
    } else if (propSeedMode.get() == 1) {

        for (int i = 0; i < propRandomCount.get(); i++) {
            dvec2 startPoint = randomPoint(vectorField);
            stepsTaken += drawStreamline(vectorField, startPoint, indexBufferPoints, mesh, vertices);
        }

    } else {
        int x = propGridPoints.get().x;
        int y = propGridPoints.get().y;
        double xi = (BBoxMax_.x - BBoxMin_.x) / (x+1);
        double yi = (BBoxMax_.y - BBoxMin_.y) / (y+1);
        double x0 = BBoxMin_.x + xi;
        double y0 = BBoxMin_.y + yi;

        for (int i = 0; i < x; i++)
            for (int j = 0; j < y; j++) {
                dvec2 startPoint(x0 + i*xi, y0 + j*yi);
                stepsTaken += drawStreamline(vectorField, startPoint, indexBufferPoints, mesh, vertices);
            }
    }

    propNumStepsTaken.set(stepsTaken);

    mesh->addVertices(vertices);
    meshOut.setData(mesh);
}

}  // namespace inviwo
