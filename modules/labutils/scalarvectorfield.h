#pragma once

#include <labutils/labutilsmoduledefine.h>
#include <inviwo/core/datastructures/volume/volume.h>
#include <inviwo/core/datastructures/volume/volumeramprecision.h>
#include <type_traits>

namespace inviwo {
/// Sets the value and data range of the data map.
template <int Dim, int VecDim>
class Field {

    // Typedefs
public:
    /** Scalarfield? */
    static constexpr bool IsScalar = (VecDim == 1);

    /** Integer vertex indices. */
    typedef glm::vec<Dim, glm::i32, glm::defaultp> IndexType;
    typedef glm::i32 IndexElementType;

    /** Position in the grid (world or grid coordinates). */
    typedef glm::vec<Dim, glm::f64, glm::defaultp> PositionType;

    /** Data type (double or vector). */
    typedef typename std::conditional<IsScalar, glm::f64,
                                      glm::vec<VecDim, glm::f64, glm::defaultp>>::type VectorType;
    /** Derivative type (gradient or Jacobian). */
    typedef typename std::conditional<IsScalar, glm::vec<Dim, glm::f64, glm::defaultp>,
                                      glm::mat<VecDim, Dim, glm::f64, glm::defaultp>>::type
        DerivativeType;

    // Functions
public:
    /** \brief Create a new field from a given inviwo volume. */
    static const Field<Dim, VecDim> createFieldFromVolume(
        std::shared_ptr<const inviwo::Volume> volume);

    /**
     * \brief Create new empty field.
     *
     * @param size Number of vertices.
     * @param offset Minimum world position of the grid (optional)
     * @param offset Maximum world position distance (optional)
     */
    Field(const IndexType& size, const PositionType& offset = PositionType(0),
          const PositionType& extent = PositionType(1));

    /** Copy constructor */
    Field(const Field<Dim, VecDim>& other);

    /** Copy assignment */
    Field<Dim, VecDim>& operator=(const Field<Dim, VecDim>& other);

    /** Default destructor */
    ~Field() {
        if (ownsData_) delete data_;
    }

    /**
     * \brief Check if world position is inside
     *
     * @param pos World position to be checked
     */
    bool isInside(const PositionType& pos) const;

    /**
     * \brief Clamp a position to the bounding box
     *
     * @param pos World position to be checked
     */
    PositionType clampPositionToBBox(const PositionType& pos) const;

    /**
     * \brief Get the world position of a vertex.
     * Use this when e.g. drawing vertices.
     *
     * @param idx Index of grid vertex.
     */
    PositionType getPositionAtVertex(const IndexType& idx) const;

    /**
     * \brief Get the value directly at a vertex.
     *
     * @param idx A vertex index.
     */
    VectorType getValueAtVertex(const IndexType& idx) const;

    /**
     * \brief Set the value at any vertex in the field.
     *
     * @param idx A vertex index.
     * @param newValue Value to be set at that index.
     */
    void setValueAtVertex(const IndexType& idx, const VectorType& newValue);

    /**
     * \brief Interpolates the data using bi- or tri-linear interpolation.
     *
     * @param pos A location inside the bounding box of the field.
     */
    VectorType interpolate(const PositionType& pos) const;

    /**
     * \brief Get the derivative at any point in the field.
     * Returns the gradient or Jacobian.
     *
     * @param pos A world position, i.e., dependent on where the field is places in space.
     */
    DerivativeType derive(const PositionType& pos) const;

    /** Returns the number of grid vertices in each dimension. */
    const IndexType& getNumVerticesPerDim() const { return size_; }

    /** Get the minimum world coordinate, i.e., the position of vertex 0. */
    const PositionType& getBBoxMin() const { return offset_; }

    /** Get the maximum world coordinate, i.e., the position of vertex size_-1. */
    const PositionType getBBoxMax() const { return offset_ + extent_; }

    /** Get the world extent of a single cell (rectangle or cuboid). */
    PositionType getCellSize() const;

    /** Get minimum value, component-wise. */
    const VectorType& getMinValue() { return minValue_; }

    /** Get maximum value, component-wise. */
    const VectorType& getMaxValue() { return maxValue_; }

protected:
    /** Get the extent of the world coordinate, i.e., the maximum distance between vertices. */
    const PositionType& getExtent() const { return extent_; }

    /**
     * \brief Interpolates the data using bi- or tri-linear interpolation.
     *
     * @param pos A grid position, i.e., a floating value between 0 and size_.
     */
    VectorType interpolateInGridCoords(const PositionType& pos) const;

    /**
     * \brief Get the derivative at any point in the field.
     * Returns the gradient or Jacobian.
     *
     * @param pos A grid position, i.e., a floating value between 0 and size_.
     */
    DerivativeType deriveInGridCoords(const PositionType& pos) const;

    /**
     * \brief Get the derivative directly at a vertex.
     * Returns the gradient or Jacobian.
     *
     * @param idx A vertex index.
     */
    DerivativeType sampleDerivativeAtVertex(const IndexType& idx) const;

    /**
     * \brief Convert a grid position into a world position.
     *
     * @param pos World position to convert.
     * */
    PositionType gridCoordsFromWorldPos(const PositionType& pos) const;

    /**
     * \brief Convert a grid position into a world position.
     *
     * @param pos Grid position to convert.
     * */
    PositionType worldPosFromGridCoords(const PositionType& pos) const;

    // Data
protected:
    Field(std::shared_ptr<const inviwo::Volume> volume);

    IndexType getLowerIndex(const PositionType& pos) const;
    VectorType sample(const size3_t& idx) const;

    VolumeRAM* data_;
    bool ownsData_ = false;
    IndexType size_;
    PositionType offset_, extent_;
    VectorType minValue_, maxValue_;
};

//--> Definitions <--//

template <int Dim, int VecDim>
const Field<Dim, VecDim> Field<Dim, VecDim>::createFieldFromVolume(
    std::shared_ptr<const inviwo::Volume> volume) {
    return Field<Dim, VecDim>(volume);
}

template <int Dim, int VecDim>
Field<Dim, VecDim>::Field(std::shared_ptr<const inviwo::Volume> volume)
    : data_(nullptr), ownsData_(false) {
    IVW_ASSERT(volume.get(), "No valid volume.");
    IVW_ASSERT(volume->getRepresentation<VolumeRAM>(), "No valid volume RAM representation.");
    data_ = const_cast<VolumeRAM*>(volume->getRepresentation<VolumeRAM>());
    size3_t numElements = volume->getDimensions();
    for (int d = 0; d < Dim; ++d) size_[d] = static_cast<IndexElementType>(numElements[d]);

    auto mat = volume->getModelMatrix();
    for (int d = 0; d < Dim; ++d) {
        offset_[d] = mat[3][d];
        extent_[d] = mat[d][d];
    }

    minValue_ = sample({0, 0, 0});
    maxValue_ = minValue_;
    int maxZ = (Dim >= 3) ? size_[2] : 1;
    for (int x = 0; x < size_[0]; ++x)
        for (int y = 0; y < size_[1]; ++y)
            for (int z = 0; z < maxZ; ++z) {
                auto value = sample({x, y, z});
                if constexpr (IsScalar) {
                    minValue_ = std::min(value, minValue_);
                    maxValue_ = std::max(value, maxValue_);
                } else {
                    for (int d = 0; d < Dim; ++d) {
                        minValue_[d] = std::min(value[d], minValue_[d]);
                        maxValue_[d] = std::max(value[d], maxValue_[d]);
                    }
                }
            }
}

template <int Dim, int VecDim>
Field<Dim, VecDim>::Field(const typename Field<Dim, VecDim>::IndexType& size,
                          const typename Field<Dim, VecDim>::PositionType& offset,
                          const typename Field<Dim, VecDim>::PositionType& extent)
    : ownsData_(true), size_(size), offset_(offset), extent_(extent) {

    typedef typename Field<Dim, VecDim>::VectorType VecType;
    size_t numElements = 1;
    ivec3 volSize(1);
    for (int d = 0; d < Dim; ++d) {
        volSize[d] = size_[d];
        numElements *= size_[d];
    }

    // VecType* volData = new VecType[numElements];
    // memset(volData, 0, numElements * sizeof(typename Field<Dim, VecDim>::VectorType));

    data_ = new VolumeRAMPrecision<VecType>(size3_t(volSize.x, volSize.y, volSize.z));
}

template <int Dim, int VecDim>
Field<Dim, VecDim>::Field(const Field<Dim, VecDim>& other)
    : data_(other.data_->clone())
    , ownsData_(true)
    , size_(other.size_)
    , offset_(other.offset_)
    , extent_(other.extent_)
    , minValue_(other.minValue_)
    , maxValue_(other.maxValue_) {}

template <int Dim, int VecDim>
Field<Dim, VecDim>& Field<Dim, VecDim>::operator=(const Field<Dim, VecDim>& other) {
    ownsData_ = true;
    size_ = other.size_;
    offset_ = other.offset_;
    extent_ = other.extent_;
    minValue_ = other.minValue_;
    maxValue_ = other.maxValue_;
    data_ = other.data_->clone();

    return *this;
}

template <int Dim, int VecDim>
bool Field<Dim, VecDim>::isInside(const typename Field<Dim, VecDim>::PositionType& pos) const {
    for (int d = 0; d < Dim; ++d)
        if (pos[d] < offset_[d] || pos[d] > offset_[d] + extent_[d]) return false;
    return true;
}

template <int Dim, int VecDim>
typename Field<Dim, VecDim>::PositionType Field<Dim, VecDim>::clampPositionToBBox(
    const typename Field<Dim, VecDim>::PositionType& pos) const {
    PositionType clampedPos = pos;
    for (int d = 0; d < Dim; ++d) {
        if (pos[d] < offset_[d]) clampedPos[d] = offset_[d];
        if (pos[d] > offset_[d] + extent_[d]) clampedPos[d] = offset_[d] + extent_[d];
    }
    return clampedPos;
}

template <int Dim, int VecDim>
typename Field<Dim, VecDim>::IndexType Field<Dim, VecDim>::getLowerIndex(
    const typename Field<Dim, VecDim>::PositionType& pos) const {
    IndexType idx;
    for (int d = 0; d < Dim; ++d) {
        IVW_ASSERT(extent_[d] >= 0, "0 dimension detected! D" << d);
        idx[d] = static_cast<int>(((pos[d] - offset_[d]) * (size_[d] - 1)) / extent_[d]);
    }
    return idx;
}

template <int Dim, int VecDim>
typename Field<Dim, VecDim>::PositionType Field<Dim, VecDim>::getPositionAtVertex(
    const typename Field<Dim, VecDim>::IndexType& idx) const {
    typename Field<Dim, VecDim>::PositionType pos;
    for (int d = 0; d < Dim; ++d)
        pos[d] = extent_[d] * static_cast<double>(idx[d]) / static_cast<double>(size_[d] - 1) +
                 offset_[d];
    return pos;
}

template <int Dim, int VecDim>
typename Field<Dim, VecDim>::VectorType Field<Dim, VecDim>::sample(const size3_t& idxT) const {
    if constexpr (IsScalar) {
        return data_->getAsDouble(idxT);
    } else {
        glm::dvec4 values = data_->getAsDVec4(idxT);
        typename Field<Dim, VecDim>::VectorType val(0);
        if constexpr (IsScalar)
            val = values[0];
        else
            for (int d = 0; d < Dim; ++d) val[d] = values[d];
        return val;
    }
}

template <int Dim, int VecDim>
typename Field<Dim, VecDim>::VectorType Field<Dim, VecDim>::interpolate(
    const typename Field<Dim, VecDim>::PositionType& pos) const {
    if (!isInside(pos)) {
        return interpolateInGridCoords(gridCoordsFromWorldPos(clampPositionToBBox(pos)));
    }
    return interpolateInGridCoords(gridCoordsFromWorldPos(pos));
}

template <int Dim, int VecDim>
typename Field<Dim, VecDim>::VectorType Field<Dim, VecDim>::interpolateInGridCoords(
    const typename Field<Dim, VecDim>::PositionType& relPos) const {
    typename Field<Dim, VecDim>::IndexType idxMin;
    typename Field<Dim, VecDim>::PositionType fracIdx;
    for (int d = 0; d < Dim; ++d) {
        idxMin[d] = static_cast<int>(relPos[d]);
        fracIdx[d] = relPos[d] - idxMin[d];
    }
    size3_t Interp;
    Interp[0] = (idxMin[0] != size_[0] - 1) ? 2 : 1;
    Interp[1] = (Dim >= 2 && idxMin[1] != size_[1] - 1) ? 2 : 1;
    Interp[2] = (Dim >= 3 && idxMin[2] != size_[2] - 1) ? 2 : 1;

    size3_t offset(0);

    VectorType val(0);
    for (offset[0] = 0; offset[0] < Interp[0]; ++offset[0]) {
        for (offset[1] = 0; offset[1] < Interp[1]; ++offset[1]) {
            for (offset[2] = 0; offset[2] < Interp[2]; ++offset[2]) {
                size3_t cornerPos(0);
                for (int d = 0; d < Dim; ++d) cornerPos[d] = idxMin[d] + offset[d];

                VectorType cornerVec = sample(cornerPos);
                double scale = 1;
                for (int d = 0; d < Dim; ++d)
                    scale *= (offset[d] == 1) ? fracIdx[d] : 1.0 - fracIdx[d];

                if constexpr (IsScalar)
                    val += scale * cornerVec;
                else {
                    for (int d = 0; d < Dim; ++d) val[d] += scale * cornerVec[d];
                }
            }
        }
    }

    return val;
}

template <int Dim, int VecDim>
typename Field<Dim, VecDim>::VectorType Field<Dim, VecDim>::getValueAtVertex(
    const typename Field<Dim, VecDim>::IndexType& idx) const {
    size3_t idxT(0);
    for (int d = 0; d < Dim; ++d) idxT[d] = idx[d];
    return sample(idxT);
}

template <int Dim, int VecDim>
typename Field<Dim, VecDim>::PositionType Field<Dim, VecDim>::gridCoordsFromWorldPos(
    const typename Field<Dim, VecDim>::PositionType& pos) const {

    typename Field<Dim, VecDim>::PositionType relativePos(0);
    for (int d = 0; d < Dim; ++d)
        relativePos[d] = ((pos[d] - offset_[d]) * (size_[d] - 1)) / extent_[d];

    return relativePos;
}

template <int Dim, int VecDim>
typename Field<Dim, VecDim>::PositionType Field<Dim, VecDim>::worldPosFromGridCoords(
    const typename Field<Dim, VecDim>::PositionType& pos) const {

    typename Field<Dim, VecDim>::PositionType worldPos(0);
    for (int d = 0; d < Dim; ++d) worldPos[d] = (pos[d] * extent_[d] / (size_[d] - 1)) + offset_[d];

    return worldPos;
}

template <int Dim, int VecDim>
void Field<Dim, VecDim>::setValueAtVertex(const Field<Dim, VecDim>::IndexType& idx,
                                          const Field<Dim, VecDim>::VectorType& newValue) {
    size3_t idxT(0);
    for (int d = 0; d < Dim; ++d) idxT[d] = idx[d];

    if constexpr (VecDim == 1) data_->setFromDouble(idxT, newValue);
    if constexpr (VecDim == 2) data_->setFromDVec2(idxT, newValue);
    if constexpr (VecDim == 3) data_->setFromDVec3(idxT, newValue);

    if constexpr (IsScalar) {
        minValue_ = std::min(newValue, minValue_);
        maxValue_ = std::max(newValue, maxValue_);
    } else {
        for (int d = 0; d < Dim; ++d) {
            minValue_[d] = std::min(newValue[d], minValue_[d]);
            maxValue_[d] = std::max(newValue[d], maxValue_[d]);
        }
    }
}

template <int Dim, int VecDim>
typename Field<Dim, VecDim>::PositionType Field<Dim, VecDim>::getCellSize() const {
    Field<Dim, VecDim>::PositionType cellSize;
    for (int d = 0; d < Dim; ++d) cellSize[d] = extent_[d] / (size_[d] - 1);
    return cellSize;
}

template <int Dim, int VecDim>
typename Field<Dim, VecDim>::DerivativeType Field<Dim, VecDim>::derive(
    const PositionType& pos) const {
    return deriveInGridCoords(gridCoordsFromWorldPos(pos));
}

template <int Dim, int VecDim>
typename Field<Dim, VecDim>::DerivativeType Field<Dim, VecDim>::deriveInGridCoords(
    const typename Field<Dim, VecDim>::PositionType& relPos) const {
    typename Field<Dim, VecDim>::DerivativeType derivative;
    typename Field<Dim, VecDim>::VectorType partDiff;
    typename Field<Dim, VecDim>::PositionType minPos, maxPos;
    for (int d = 0; d < Dim; ++d) {
        minPos = relPos;
        maxPos = relPos;
        minPos[d] = std::max(minPos[d] - 0.5, 0.0);
        maxPos[d] = std::min(minPos[d] + 0.5, static_cast<double>(size_[d] - 1));

        partDiff = (interpolateInGridCoords(maxPos) - interpolateInGridCoords(minPos)) /
                   (maxPos[d] - minPos[d]);
        derivative[d] = partDiff;
    }
    return derivative;
}

template <int Dim, int VecDim>
typename Field<Dim, VecDim>::DerivativeType Field<Dim, VecDim>::sampleDerivativeAtVertex(
    const typename Field<Dim, VecDim>::IndexType& idx) const {
    typename Field<Dim, VecDim>::DerivativeType derivative;
    typename Field<Dim, VecDim>::VectorType partDiff;
    typename Field<Dim, VecDim>::IndexType minIdx, maxIdx;
    for (int d = 0; d < Dim; ++d) {
        minIdx = idx;
        maxIdx = idx;
        minIdx[d] = std::max(minIdx[d] - 1, 0);
        maxIdx[d] = std::min(minIdx[d] + 1, size_[d] - 1);

        partDiff = (interpolateInGridCoords(maxIdx) - interpolateInGridCoords(minIdx)) /
                   (maxIdx[d] - minIdx[d]);
        derivative[d] = partDiff;
    }
    return derivative;
}

//--> Scalar Fields <--//
typedef Field<2, 1> ScalarField2;
typedef Field<2, 2> VectorField2;
typedef Field<3, 1> ScalarField3;
typedef Field<3, 3> VectorField3;
}  // namespace inviwo
