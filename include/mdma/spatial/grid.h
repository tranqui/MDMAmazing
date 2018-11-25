#ifndef _MDMA_GRID_H
#define _MDMA_GRID_H

#include <memory>

namespace mdma
{
    namespace spatial
    {
        namespace details
        {
            template <size_t SpatialDimensions>
            using ArrayIndex = std::array<int, SpatialDimensions>;
            template <size_t SpatialDimensions>
            using EigenIndex = Eigen::Matrix<int, 1, SpatialDimensions>;

            template <typename T, size_t SpatialDimensions>
            class GridBase
            {
            public:
                explicit GridBase(const ArrayIndex<SpatialDimensions>& dimensions)
                    : GridBase(EigenIndex<SpatialDimensions>(dimensions.data()))
                { }

                explicit GridBase(const Eigen::Ref<const EigenIndex<SpatialDimensions>>& dimensions)
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

                inline T& operator()(const ArrayIndex<SpatialDimensions>& index)
                {
                    return this->operator[] (this->offset(index));
                }

                inline const T& operator()(const ArrayIndex<SpatialDimensions>& index) const
                {
                    return this->operator[] (this->offset(index));
                }

                inline T& operator()(const Eigen::Ref< const EigenIndex<SpatialDimensions> >& index)
                {
                    return this->operator[] (this->offset(index));
                }

                inline const T& operator()(const Eigen::Ref< const EigenIndex<SpatialDimensions> >& index) const
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

                inline EigenIndex<SpatialDimensions> shape() const
                {
                    return this->data_shape;
                }
                inline EigenIndex<SpatialDimensions> stride() const
                {
                    return this->data_stride;
                }

                // for debug purposes
                EigenIndex<SpatialDimensions> index_from_offset(size_t offset) const
                {
                    EigenIndex<SpatialDimensions> index = EigenIndex<SpatialDimensions>::Zero();

                    for (size_t c = SpatialDimensions-1; c > 0; --c)
                    {
                        index[c] = offset / this->data_stride[c];
                        offset -= index[c] * this->data_stride[c];
                    }
                    index[0] = offset;

                    return index;
                }

                inline int offset(const Eigen::Ref< const EigenIndex<SpatialDimensions> >& index) const
                {
                    return index.dot(this->data_stride);
                }

                inline int offset(const ArrayIndex<SpatialDimensions>& index) const
                {
                    int offset = 0;
                    for (size_t c = 0; c < SpatialDimensions; ++c)
                        offset += this->data_stride[c] * index[c];
                    return offset;
                }

            protected:
                //ArrayIndex<SpatialDimensions> data_shape, data_stride;
                EigenIndex<SpatialDimensions> data_shape, data_stride;

                std::vector<T> data;
            };
        }

        template <typename T, size_t SpatialDimensions>
        struct Cell
        {
            using EigenIndex = details::EigenIndex<SpatialDimensions>;

            struct Neighbour
            {
                // Rule of three so containers behave as expected.
                Neighbour() = default;
                Neighbour(const Neighbour&) = default;
                Neighbour& operator=(const Neighbour&) = default;
                
                Neighbour(Cell<T,SpatialDimensions>* cell,
                          const Eigen::Ref<const EigenIndex>& image_offset)
                    : cell(cell), image_offset(image_offset)
                { }

                Cell<T,SpatialDimensions>* cell;
                EigenIndex image_offset;
            };

            std::list<T> children;
            std::vector<Neighbour> neighbours;
            std::vector<Neighbour> forward_neighbours;

            auto begin() noexcept
            {
                return this->children.begin();
            }
            auto begin() const noexcept
            {
                return this->children.begin();
            }
            auto end() noexcept
            {
                return this->children.end();
            }
            auto end() const noexcept
            {
                return this->children.end();
            }
        };

        template <typename T, size_t SpatialDimensions>
        class PeriodicGrid : public details::GridBase<Cell<T,SpatialDimensions>, SpatialDimensions>
        {
        public:
            static constexpr size_t d = SpatialDimensions;

            using Base = details::GridBase<Cell<T,d>, d>;
            using ArrayIndex = details::ArrayIndex<d>;
            using EigenIndex = details::EigenIndex<d>;

            // Number of nearest neighbours on a cubic lattice.
            static constexpr size_t total_nearest_neighbours = std::pow(3,d)-1;

            template <typename Index>
            PeriodicGrid(const Index& dimensions) : Base(dimensions)
            {
                for (auto& child : this->data)
                {
                    child.neighbours.reserve(total_nearest_neighbours);
                    child.forward_neighbours.reserve(total_nearest_neighbours/2);
                }

                auto neighbours = this->possible_neighbours();

                for (auto& nb : neighbours)
                {
                    int offset = this->offset(nb);
                    if (not offset) continue;
                    bool forward = offset > 0;

                    for (size_t i = 0; i < this->data.size(); ++i)
                    {
                        EigenIndex index = this->index_from_offset(i) + nb;
                        auto wrapping = this->wrap_index(index);

                        Cell<T,d>* cell = &((*this)(index));

                        auto& child = this->data[i];
                        child.neighbours.emplace_back(cell, wrapping);
                        if (forward) child.forward_neighbours.push_back(child.neighbours.back());
                    }
                }
            }

            auto possible_neighbours() const
            {
                std::array<EigenIndex, 1+total_nearest_neighbours> neighbours;

                EigenIndex index = EigenIndex::Constant(-1);
                for (size_t i = 0; i < neighbours.size(); ++i)
                {
                    neighbours[i] = index;

                    size_t c = 0;
                    while (c < d and ++index[c] > 1)
                    {
                        index[c] = -1;
                        c++;
                    }
                }

                return neighbours;
            }

            EigenIndex wrap_index(EigenIndex& index) const
            {
                EigenIndex wrapping;

                auto shape = this->shape();
                for (size_t c = 0; c < d; ++c)
                {
                    if (index(c) >= shape(c)) wrapping(c) = shape(c);
                    else if (index(c) < 0) wrapping(c) = -shape(c);
                    else wrapping(c) = 0;
                }

                index -= wrapping;
                return wrapping;
            }
        };
    }
}

#endif
