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
    , propKernelSize("kernelSize", "Kernel size", 10, 1, 1000)
    ,propMu("mu", "Mu", 0.5, 0.0,1.0)
    ,propSigma("propsigma", "prop Sigma", 0.1, 0.0,1.0)
    , propFastLIC("fastLIC", "FastLIC")
    , propContrastEnhance("contrastEnhance","Contrast Enhance")
    , textureColor("texture_Color", "texture color")
    ,propTransferFunc("TransferFunc", "Colors", &volumeIn_)
{
    // Register ports
    addPort(volumeIn_);
    addPort(noiseTexIn_);
    addPort(licOut_);

    // Register properties
    
    addProperty(propKernelSize);
    addProperty(propMu);
    addProperty(propSigma);
    addProperty(propFastLIC);
    addProperty(propContrastEnhance);
    addProperty(textureColor);
      addProperty(propTransferFunc);
    
      // The default transfer function has just two blue points
      propTransferFunc.get().clear();
      //propIsoTransferFunc.get().add(0.0f, vec4(0.0f, 0.0f, 1.0f, 1.0f));
      propTransferFunc.get().add(0.0f, vec4(0.0f, 0.0f, 1.0f, 1.0f));
      //propTransferFunc.get().add(0.5f, vec4(0.0f, 1.0f, 0.0f, 1.0f));
      propTransferFunc.get().add(1.0f, vec4(1.0f, 0.0f, 0.0f, 1.0f));
      propTransferFunc.setCurrentStateAsDefault();

}

void streamline(const VectorField2& vectorField, const dvec2& seed, double stepSize, int numSteps, std::vector<dvec2>& points)
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
    std::vector<std::vector<double>> magnitudeF(texDims_.x,std::vector<double>(texDims_.y));
    double magmax=magnitudeField(vectorField, texDims_, magnitudeF);
    // Prepare the output, it has the same dimensions as the texture and rgba values in [0,255]
    auto outImage = std::make_shared<Image>(texDims_, DataVec4UInt8::get());
    RGBAImage licImage(outImage);

    auto pixel = pixelToField({1, 1}) - pixelToField({0, 0});
    double stepSize = std::min(pixel.x, pixel.y);
    LogProcessorInfo("stepSize: " << stepSize);
    
    std::vector<std::vector<double>> licTex(texDims_.x,std::vector<double>(texDims_.y));
    double Amean=0.0;
        double p2=0.0;
        int nonblack=0;
    bool contrast=propContrastEnhance.get();
    bool colortexturing=textureColor.get();
    if (!propFastLIC.get()) {
        // Half of the steps in either direction
        int numSteps = propKernelSize.get() / 2;

        // Reusable vectors to store streamlines to avoid allocating in loop
        std::vector<dvec2> forward, backward;

        for (size_t j = 0; j < texDims_.y; j++) {
            for (size_t i = 0; i < texDims_.x; i++) {
                // Middle of pixel
                dvec2 vecPoint = pixelToField({i+0.5, j+0.5});

                // Integrate in both directions, store in forward/backward vectors
                streamline(vectorField, vecPoint, stepSize, numSteps, forward);
                streamline(vectorField, vecPoint, -stepSize, numSteps, backward);
                
                // Sample the points and sum
                int sum = texture.sampleGrayScale({i+0.5, j+0.5});
                double sum2=texture.sampleGrayScale({i+0.5, j+0.5});
                for (size_t k = 0; k < forward.size(); k++)
                {
                    sum += texture.sampleGrayScale(fieldToPixel(forward[k]));
                    sum2 += texture.sampleGrayScale(fieldToPixel(forward[k]));
                }
                    
                for (size_t k = 0; k < backward.size(); k++)
                {
                    sum += texture.sampleGrayScale(fieldToPixel(backward[k]));
                    sum2+=texture.sampleGrayScale(fieldToPixel(backward[k]));
                }
                
                // Arithmetic mean
                int val = sum / (double)(1 + forward.size() + backward.size());
                licTex[i][j]=sum2/(double)(1 + forward.size() + backward.size());
                if(contrast&&licTex[i][j]>1.0)
                {
                    Amean+=licTex[i][j];
                    p2+=licTex[i][j]*licTex[i][j];
                    nonblack++;
                }
                licImage.setPixelGrayScale(size2_t(i, j), val);
            }
        }
    } else {
        // Some reasonable max steps (circumference of canvas)
        int maxSteps = 2 * (texDims_.x + texDims_.y);

        // Which pixels have already been taken care of
        auto seen = new bool[texDims_.x * texDims_.y]();

        // Streamlines...
        std::vector<dvec2> forward, backward;

        // ...collected to one and converted to pixel coords and sampled
        std::vector<dvec2> linePix;
        std::vector<int> lineVal;

        for (size_t j = 0; j < texDims_.y; j++) {
            for (size_t i = 0; i < texDims_.x; i++) {
                // Skip already seen pixels
                size_t seen_idx = j*texDims_.x + i;
                if (seen[seen_idx])
                    continue;

                dvec2 vecPoint = pixelToField({i+0.5, j+0.5});

                streamline(vectorField, vecPoint, stepSize, maxSteps, forward);
                streamline(vectorField, vecPoint, -stepSize, maxSteps, backward);

                linePix.clear();
                lineVal.clear();
                for (int k = backward.size()-1; k >= 0; k--) {
                    auto pixel = fieldToPixel(backward[k]);
                    linePix.push_back(pixel);
                    lineVal.push_back(texture.sampleGrayScale(pixel));
                    seen[(size_t)pixel.y*texDims_.x + (size_t)pixel.x] = true;
                }
                linePix.push_back({i+0.5, j+0.5});
                lineVal.push_back(texture.sampleGrayScale({i+0.5, j+0.5}));
                for (size_t k = 0; k < forward.size(); k++) {
                    auto pixel = fieldToPixel(forward[k]);
                    linePix.push_back(pixel);
                    lineVal.push_back(texture.sampleGrayScale(pixel));
                    seen[(size_t)pixel.y*texDims_.x + (size_t)pixel.x] = true;
                }

                // Now we will compute the box filter
                int sum = 0;
                // Left and right of rolling window
                size_t a = 0, b = std::min(linePix.size(), (size_t)propKernelSize.get());
                for (size_t k = 0; k < b; k++)
                    sum += lineVal[k];

                // The first pixels will have the same value since we use the first propKernelSize points for all of them
                int val0 = sum / b;
                for (size_t k = 0; k <= b/2; k++)
                {
                    licImage.setPixelGrayScale(linePix[k], val0);
                    int ii=linePix[k].x;
                    int jj=linePix[k].y;
                    licTex[ii][jj]=double(sum)/b;
                    //auto contrast
                    if(contrast&&licTex[ii][jj]>1.0)
                    {
                        Amean+=licTex[ii][jj];
                        p2+=licTex[ii][jj]*licTex[ii][jj];
                        nonblack++;
                    }
                }
                    

                // Do the rolling average
                size_t k;
                size_t end = linePix.size() - b/2;
                for (k = b/2 + 1; k < end; k++) {
                    sum -= lineVal[a++];
                    sum += lineVal[b++];

                    int val = sum / (b-a);
                    licImage.setPixelGrayScale(linePix[k], val);
                    int ii=linePix[k].x;
                    int jj=linePix[k].y;
                    licTex[ii][jj]=double(sum)/(b-a);
                    //auto contrast
                    if(contrast&&licTex[ii][jj]>1.0)
                    {
                        Amean+=licTex[ii][jj];
                        p2+=licTex[ii][jj]*licTex[ii][jj];
                        nonblack++;
                    }
                }

                // The last ones also have the same value
                int val1 = sum / (b-a);
                while (k < linePix.size()) {
                    int ii=linePix[k].x;
                    int jj=linePix[k].y;
                    licTex[ii][jj]=double(sum)/(b-a);
                    if(contrast&&licTex[ii][jj]>1.0)
                    {
                        Amean+=licTex[ii][jj];
                        p2+=licTex[ii][jj]*licTex[ii][jj];
                        nonblack++;
                    }
                    licImage.setPixelGrayScale(linePix[k++], val1);
                    
                }
            }
        }
        delete[] seen;
    }
    if(contrast)
    {
        contrastAuto(licImage ,vectorField, texDims_, Amean, p2,nonblack, licTex);
    }
    if(colortexturing)
        ColorLic(licImage, texDims_, magnitudeF, magmax,licTex);
    licOut_.setData(outImage);
}

void LICProcessor::contrastAuto(RGBAImage& outImage,const VectorField2& vectorField,const dvec2 textdim,double Amean,double p2,int nonblack,std::vector<std::vector<double>>& licTex)
{
    double mu=Amean/nonblack;
    double mymu=propMu.get()*255;
    double mySigma=propSigma.get()*255;
    for(int i=0;i<textdim.x;i++)
    {
        for(int j=0;j<textdim.y;j++)
        {
            double sigma=sqrt((p2-nonblack*mu*mu)/(nonblack-1));
            
            double f=mySigma/sigma;
            //int pij=outImage.readPixelGrayScale(size2_t(i,j));
            double pij=licTex[i][j];
            int val=(mymu+f*(pij-mu));
            //licTex[i][j]=val;
            outImage.setPixelGrayScale(size2_t(i,j), val);
        }
    }
}



double LICProcessor::magnitudeField(const VectorField2& vectorField,const dvec2 textdim,
                                std::vector<std::vector<double>>& magnitudeF)
{
    double max_mag = 0;

    for(int i=0;i<textdim.x;i++)
    {
        for(int j=0;j<textdim.y;j++)
        {
            dvec2 vec = vectorField.interpolate(pixelToField({i+0.5,j+0.5}));
            if (vec != dvec2(0)) {
                magnitudeF[i][j] = glm::length(vec);
                if(max_mag < magnitudeF[i][j])
                    max_mag =magnitudeF[i][j];
            }
        }
    }
    return max_mag;
}
void LICProcessor::ColorLic(RGBAImage& outImage,const dvec2 & textdim,
                            std::vector<std::vector<double>>& magnitudeF,const double max_mag,std::vector<std::vector<double>>& licTex)
{
    for(int i=0;i<textdim.x;i++)
    {
        for(int j=0;j<textdim.y;j++)
        {
            //int pij=outImage.readPixelGrayScale(size2_t(i,j));
            double pij=licTex[i][j];
            //int val=pij*magnitudeF[i][j] /max_mag;
            double magratio = magnitudeF[i][j] / max_mag;
            //dvec4 color = 255 * propTransferFunc.get().sample(magratio);
            //outImage.setPixel(size2_t(i,j), (color+dvec4(val,val,val,255))*0.5);
            dvec4 transferColor = propTransferFunc.get().sample(magratio);
            dvec4 color = pij * transferColor;
            color.w = transferColor.w;
            outImage.setPixel(size2_t(i,j), color);
        }
    }
}

}  // namespace inviwo
