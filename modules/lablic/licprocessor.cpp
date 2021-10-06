/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Monday, October 02, 2017 - 13:31:17
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <inviwo/core/datastructures/volume/volumeram.h>
#include <lablic/licprocessor.h>
#include <labstreamlines/integrator.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo LICProcessor::processorInfo_{
    "org.inviwo.LICProcessor",  // Class identifier
    "LICProcessor",             // Display name
    "KTH Labs",                 // Category
    CodeState::Experimental,    // Code state
    Tags::None,                 // Tags
};

const ProcessorInfo LICProcessor::getProcessorInfo() const { return processorInfo_; }

LICProcessor::LICProcessor()
    : Processor()
    , volumeIn_("volIn")
    , noiseTexIn_("noiseTexIn")
    , licOut_("licOut")
// TODO: Register additional properties
{
    // Register ports
    addPort(volumeIn_);
    addPort(noiseTexIn_);
    addPort(licOut_);

    // Register properties
    // TODO: Register additional properties
}

void streamline(const VectorField2& vectorField, const dvec2& seed, double stepSize, int numSteps, std::vector<dvec2>& points)
{
    points.clear();

    dvec2 x0 = seed;
    for (int i = 0; i < numSteps; i++) {
        dvec2 x1 = Integrator::RK4_norm(vectorField, x0, stepSize);
        //LogInfoCustom("", x1);

        if (!vectorField.isInside(x1)) {
            //LogInfoCustom("", "Not inside");
            break;
        }

        if (vectorField.interpolate(x1) == dvec2(0)) {
            //LogInfoCustom("", "Zero");
            break;
        }

        points.push_back(x1);

        x0 = x1;
    }
}

dvec2 LICProcessor::pixelToField(dvec2 pixel)
{
    return vectorFieldMin_ + dvec2(
        (pixel.x / texDims_.x) * (vectorFieldMax_.x - vectorFieldMin_.x),
        (pixel.y / texDims_.y) * (vectorFieldMax_.y - vectorFieldMin_.y));
}

dvec2 LICProcessor::fieldToPixel(dvec2 vec)
{
    return {
        (vec.x - vectorFieldMin_.x) / (vectorFieldMax_.x - vectorFieldMin_.x) * texDims_.x,
        (vec.y - vectorFieldMin_.y) / (vectorFieldMax_.y - vectorFieldMin_.y) * texDims_.y};
}

void LICProcessor::process() {
    // Get input
    if (!volumeIn_.hasData()) {
        return;
    }

    if (!noiseTexIn_.hasData()) {
        return;
    }

    auto vol = volumeIn_.getData();
    const VectorField2 vectorField = VectorField2::createFieldFromVolume(vol);
    vectorFieldDims_ = vol->getDimensions();
    LogProcessorInfo("vectorFieldDims_: " << vectorFieldDims_);

    vectorFieldMin_ = vectorField.getBBoxMin();
    vectorFieldMax_ = vectorField.getBBoxMax();

    LogProcessorInfo("vectorField bbox: " << vectorFieldMin_ << " " << vectorFieldMax_);

    auto tex = noiseTexIn_.getData();
    const RGBAImage texture = RGBAImage::createFromImage(tex);
    texDims_ = tex->getDimensions();
    LogProcessorInfo("texDims_: " << texDims_);

    // Prepare the output, it has the same dimensions as the texture and rgba values in [0,255]
    auto outImage = std::make_shared<Image>(texDims_, DataVec4UInt8::get());
    RGBAImage licImage(outImage);

    std::vector<std::vector<double>> licTexture(texDims_.x, std::vector<double>(texDims_.y, 0.0));

    // Hint: Output an image showing which pixels you have visited for debugging
    std::vector<std::vector<int>> visited(texDims_.x, std::vector<int>(texDims_.y, 0));

    // TODO: Implement LIC and FastLIC

    /*double stepSize = std::min(
            vectorFieldDims_.x / (double)texDims_.x,
            vectorFieldDims_.y / (double)texDims_.y);*/
    auto pixel = pixelToField({1, 1}) - pixelToField({0, 0});
    double stepSize = std::min(pixel.x, pixel.y);
    LogProcessorInfo("stepSize: " << stepSize);

    int numSteps = 100;

    std::vector<dvec2> forward, backward;

    for (size_t j = 0; j < texDims_.y; j++) {
        for (size_t i = 0; i < texDims_.x; i++) {
            dvec2 vecPoint = pixelToField({i+0.5, j+0.5});

            //if (j < 10 && i < 10)
                //LogProcessorInfo("Point " << size2_t(i, j) << " " << vecPoint << " " << fieldToPixel(vecPoint));

            //LogProcessorInfo("Seed inside " << vectorField.isInside(vecPoint));

            streamline(vectorField, vecPoint, stepSize, numSteps, forward);
            streamline(vectorField, vecPoint, -stepSize, numSteps, backward);
            
            int sum = texture.readPixelGrayScale(size2_t(i, j));
            for (size_t k = 0; k < forward.size(); k++)
                sum += texture.sampleGrayScale(fieldToPixel(forward[k]));
            for (size_t k = 0; k < backward.size(); k++)
                sum += texture.sampleGrayScale(fieldToPixel(backward[k]));

            //if (j < 10 && i < 10)
                //LogProcessorInfo("Line: " << forward.size() << " + 1 + " << backward.size());

            int val = sum / (double)(1 + forward.size() + backward.size());

            //int val = int(licTexture[i][j]);
            //int val = texture.readPixelGrayScale(size2_t(i, j));
            licImage.setPixelGrayScale(size2_t(i, j), val);
            //licImage.setPixel(size2_t(i, j), dvec4(val, val, val, 255));
        }
    }

    licOut_.setData(outImage);
}

}  // namespace inviwo
