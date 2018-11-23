#ifndef _MDMA_GRID_H
#define _MDMA_GRID_H

#include <memory>

namespace mdma
{
    template <typename T, size_t SpatialDimensions>
    class Grid
    {
    public:
        using ArrayIndex = std::array<size_t, SpatialDimensions>;
        using EigenIndex = Eigen::Matrix<size_t, 1, SpatialDimensions>;

        explicit Grid(const ArrayIndex& dimensions)
            : Grid(EigenIndex(dimensions.data()))
        { }

        explicit Grid(Eigen::Ref<const EigenIndex> dimensions)
            : data_shape(dimensions)
        {
            size_t size = 1;
            for (size_t c = 0; c < SpatialDimensions; ++c)
            {
                this->data_stride[c] = size;
                size *= dimensions[c];
            }
            this->data = std::vector<T>(size);
        }

        inline size_t size() const { return this->data.size(); }
        inline size_t capacity() const { return this->data.capacity(); }

        auto begin() noexcept
        {
            return this->data.begin();
        }
        auto begin() const noexcept
        {
            return this->data.begin();
        }
        auto end() noexcept
        {
            return this->data.end();
        }
        auto end() const noexcept
        {
            return this->data.end();
        }

        inline T& operator[](size_t index)
        {
            if (index > this->size())
                throw std::runtime_error("out of bounds error!");
            return this->data[index];
        }
        inline const T& operator[](size_t index) const
        {
            if (index > this->size())
                throw std::runtime_error("out of bounds error!");
            return this->data[index];
        }

        // undefined if elements in lookup != elemenets in dim
        // undefined if any element in lookup
        // is greater than or equal to corresponding dim element
        inline T& operator()(const ArrayIndex& index)
        {
            return this->operator[] (this->offset(index));
        }

        inline const T& operator()(const ArrayIndex& index) const
        {
            return this->operator[] (this->offset(index));
        }

        inline T& operator()(Eigen::Ref<const EigenIndex> index)
        {
            return this->operator[] (this->offset(index));
        }

        inline const T& operator()(Eigen::Ref<const EigenIndex> index) const
        {
            return this->operator[] (this->offset(index));
        }

        template <typename ... Indices>
        inline T& operator()(const Indices&&... indices)
        {
            static_assert(sizeof...(Indices) == SpatialDimensions,
                          "Number of args does not match dimensions");
            return this->operator() ({static_cast<unsigned>(indices)...});
        }

        template <typename ... Indices>
        inline const T& operator()(const Indices&&... indices) const
        {
            static_assert(sizeof...(Indices) == SpatialDimensions,
                          "Number of args does not match dimensions");
            return this->operator() ({static_cast<unsigned>(indices)...});
        }

        inline EigenIndex shape() const
        {
            return this->data_shape;
        }
        inline EigenIndex stride() const
        {
            return this->data_stride;
        }

        // for debug purposes
        EigenIndex index_from_offset(size_t offset) const
        {
            EigenIndex index = EigenIndex::Zero();

            for (size_t c = SpatialDimensions-1; c > 0; --c)
            {
                index[c] = offset / this->data_stride[c];
                offset -= index[c] * this->data_stride[c];
            }
            index[0] = offset;

            return index;
        }

        inline size_t offset(Eigen::Ref<const EigenIndex> index) const
        {
            return index.dot(this->data_stride);
        }

        inline size_t offset(const ArrayIndex& index) const
        {
            size_t offset = 0;
            for (size_t c = 0; c < SpatialDimensions; ++c)
                offset += this->data_stride[c] * index[c];
            return offset;
        }

    private:
        //ArrayIndex data_shape, data_stride;
        EigenIndex data_shape, data_stride;

        std::vector<T> data;
    };
}

#endif
