/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Monday, October 02, 2017 - 13:31:17
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#pragma once

#include <inviwo/core/common/inviwo.h>
#include <inviwo/core/datastructures/image/imageram.h>
#include <inviwo/core/ports/imageport.h>
#include <inviwo/core/ports/volumeport.h>
#include <inviwo/core/processors/processor.h>
#include <inviwo/core/properties/boolproperty.h>
#include <inviwo/core/properties/ordinalproperty.h>
#include <lablic/lablicmoduledefine.h>
#include <labutils/scalarvectorfield.h>
#include <labutils/rgbaimage.h>
#include <inviwo/core/properties/transferfunctionproperty.h>
namespace inviwo {

/** \docpage{org.inviwo.LICProcessor, LICProcessor}
    ![](org.inviwo.LICProcessor.png?classIdentifier=org.inviwo.LICProcessor)

    Line Integral Convolution with a box kernel.

    ### Inports
      * __vectorField__ 2-dimensional vector field (with vectors of
      two components thus two values within each voxel)
      This processor deals with 2-dimensional data only, therefore it is assumed
      the z-dimension will have size 1 otherwise the 0th slice of the volume
      will be processed.
      * __texture__ Texture to be convolved along the streamlines.

    ### Outports
      * __image__ The image resulting from smearing the given texture
      the streamlines of the given vector field.
*/
class IVW_MODULE_LABLIC_API LICProcessor : public Processor {
    // Friends
    // Types
public:
    // Construction / Deconstruction
public:
    LICProcessor();
    virtual ~LICProcessor() = default;

    // Methods
public:
    virtual const ProcessorInfo getProcessorInfo() const override;
    static const ProcessorInfo processorInfo_;

protected:
    /// Our main computation function
    virtual void process() override;

    // (TODO: Helper functions can be defined here and then implemented in the .cpp)
    // e.g. something like a function for standardLIC, fastLIC, autoContrast, ...
    double magnitudeField(const VectorField2& vectorField,const dvec2 textdim,
                                        std::vector<std::vector<double>>& magnitudeF);
    void ColorLic(RGBAImage& outImage,const dvec2 & textdim,
                                std::vector<std::vector<double>>& magnitudeF,const double max_mag,std::vector<std::vector<double>>& licTex);
    void contrastAuto(RGBAImage& outImage,const VectorField2& vectorField,const dvec2 textdim,double Amean,double P2,int nonblack,std::vector<std::vector<double>>& licTex);
private:
    dvec2 pixelToField(dvec2 pixel);
    dvec2 fieldToPixel(dvec2 vec);

    // Ports
public:
    // Input vector field
    VolumeInport volumeIn_;

    // Input texture
    ImageInport noiseTexIn_;

    // Output image
    ImageOutport licOut_;

    // Properties
public:
   
    IntProperty propKernelSize;
    FloatProperty propMu;
    FloatProperty propSigma;
    BoolProperty propFastLIC;
    BoolProperty propContrastEnhance;
    BoolProperty textureColor;
    TransferFunctionProperty propTransferFunc;

    // Attributes
private:
    size3_t vectorFieldDims_;
    size2_t texDims_;
    dvec2 vectorFieldMin_;
    dvec2 vectorFieldMax_;
};

}  // namespace inviwo
