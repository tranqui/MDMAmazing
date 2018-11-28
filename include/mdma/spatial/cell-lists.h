#ifndef _MDMA_CELL_LISTS_H
#define _MDMA_CELL_LISTS_H

#include <list>
#include "grid.h"

namespace mdma
{
    namespace spatial
    {
        namespace details
        {
            static constexpr int AtomsAtCompileTime = Eigen::Dynamic;
            static constexpr int DimensionsAtCompileTime = 3;

            template <typename Scalar>
            using Domain = Eigen::Matrix<Scalar, AtomsAtCompileTime,
                                         DimensionsAtCompileTime, Eigen::RowMajor>;

            template <typename Scalar>
            using Coordinate = Eigen::Matrix<Scalar, 1, DimensionsAtCompileTime>;

            template <typename Scalar>
            using Axis = Eigen::Matrix<Scalar, AtomsAtCompileTime, 1>;
        }

        template <typename Scalar, template<typename,size_t> typename GridType = PeriodicGrid>
        class CellLists
        {
        public:
            using Domain = details::Domain<Scalar>;
            using Coordinate = details::Coordinate<Scalar>;
            using Axis = details::Axis<Scalar>;

            static constexpr size_t d = details::DimensionsAtCompileTime;

            using Grid = GridType<size_t,d>;
            using ArrayIndex = typename Grid::ArrayIndex;
            using EigenIndex = typename Grid::EigenIndex;

            struct Particle
            {
                Particle(int cell, const Coordinate& position)
                    : cell(cell), position(position)
                { }

                int cell;
                const Coordinate& position;
            };

            CellLists(Scalar neighbour_cutoff,
                      const Eigen::Ref<const Coordinate>& box_dimensions,
                      const Eigen::Ref<const Domain>& coordinates)
                : neighbour_cutoff(neighbour_cutoff),
                  box_dimensions(box_dimensions),
                  coordinates(coordinates),
                  grid(optimal_grid(neighbour_cutoff, box_dimensions))
            {
                for (size_t c = 0; c < d; ++c)
                    cell_widths[c] = box_dimensions[c] / this->shape()[c];

                this->particles.reserve(this->n());

                for (size_t atom = 0; atom < this->n(); ++atom)
                {
                    auto&& coordinate = this->coordinates.row(atom);

                    EigenIndex index;
                    for (size_t c = 0; c < d; ++c)
                    {
                        while (coordinate[c] > this->box_dimensions[c])
                            coordinate[c] -= box_dimensions[c];
                        while (coordinate[c] < 0)
                            coordinate[c] += box_dimensions[c];
                        index[c] = std::floor(coordinate[c] / this->cell_widths[c]);
                    }

                    auto offset = this->grid.offset(index);
                    this->grid[offset].children.push_back(atom);
                    this->particles.emplace_back(offset, coordinate);
                }
            }

            /**
             * Determine optimal number of cells in each spatial dimension for determining
             *   neighbours up to a cutoff.
             *
             * Optimal condition is such that (maximum) cell widths are half the cutoff
             *   length. Neighbours are thus guaranteed to be in the current or adjacent
             *   cells.
             *
             * @param neighbour_cutoff: cutoff to optimise for.
             * @param box_dimensions: physical size of bounding box cells should span in each dimension.
             * @return array giving number of cells in each dimension
             */
            inline EigenIndex optimal_grid(Scalar neighbour_cutoff,
                                           const Eigen::Ref<const Coordinate>& box_dimensions)
            {
                EigenIndex ncells;
                for (size_t c = 0; c < d; ++c)
                    ncells[c] = std::floor(box_dimensions[c] / neighbour_cutoff);

                return ncells;
            }

            inline size_t n() const
            {
                return this->coordinates.rows();
            }

            inline EigenIndex shape() const
            {
                return this->grid.shape();
            }

            inline auto get_cells() const
            {
                Eigen::Matrix<int, details::AtomsAtCompileTime, 1> cells(this->n());

                for (size_t atom = 0; atom < this->n(); ++atom)
                    cells(atom) = this->particles[atom].cell;

                return cells;
            }

            inline Coordinate cell_dimensions() const
            {
                return this->cell_widths;
            }

            class Neighbours;
            friend class Neighbours;
            class Neighbourhood;
            friend class Neighbourhood;

            inline Neighbours neighbours() const
            {
                return Neighbours(*this);
            }
            inline Neighbourhood neighbourhood(int atom) const
            {
                return Neighbourhood(*this, atom);
            }

            Domain find_displacements() const
            {
                std::list<Coordinate> displacements;

                auto function = [&](Coordinate&& delta)
                    { displacements.push_back(delta); };
                this->neighbours().for_each(function);

                Domain delta(displacements.size(), d);
                size_t index = 0;
                for (auto dr : displacements)
                {
                    delta.row(index) = dr;
                    ++index;
                }

                return delta;
            }

            auto find_distances() const
            {
                Domain displacements = this->find_displacements();

                Eigen::Matrix<Scalar, details::AtomsAtCompileTime, 1>
                    distances(displacements.rows());

                for (int i = 0; i < displacements.rows(); ++i)
                    distances(i) = displacements.row(i).norm();

                return distances;
            }

        protected:
            Scalar neighbour_cutoff;
            Coordinate box_dimensions;
            Domain coordinates;

            Grid grid;
            Coordinate cell_widths;

            std::vector<Particle> particles;
        };

        template <typename Scalar, template<typename,size_t> typename GridType>
        class CellLists<Scalar, GridType>::Neighbours
        {
        public:
            Neighbours(const CellLists<Scalar, GridType>& cell_lists)
                : cell_lists(cell_lists)
            { }

            inline void for_each(std::function<void (Coordinate&&)> function) const
            {
                for (auto& cell : this->cell_lists.grid)
                    this->for_each(cell, function);
            }

        private:
            const CellLists<Scalar, GridType>& cell_lists;

            inline void for_each(const Cell<size_t,d>& cell,
                                 std::function<void (Coordinate&&)> function) const
            {
                this->for_each_internal(cell, function);

                Coordinate offset;
                for (auto& neighbour : cell.forward_neighbours)
                {
                    for (size_t c = 0; c < d; ++c)
                        offset(c) = neighbour.image_offset(c) *
                            this->cell_lists.cell_widths(c);
                    this->for_each_external(cell, *neighbour.cell, offset, function);
                }
            }

            inline void for_each_internal(const Cell<size_t,d>& cell,
                                          std::function<void (Coordinate&&)> function) const
            {
                auto it = std::begin(cell);
                auto end = std::end(cell);
                if (it == end) return;

                for (; it != std::prev(end); ++it)
                    for (auto it2 = std::next(it); it2 != cell.end(); ++it2)
                        function(this->cell_lists.coordinates.row(*it) -
                                 this->cell_lists.coordinates.row(*it2));
            }

            inline void for_each_external(const Cell<size_t,d>& cell_a,
                                          const Cell<size_t,d>& cell_b,
                                          const Eigen::Ref<const Coordinate>& offset,
                                          std::function<void (Coordinate&&)> function) const
            {
                for (auto i : cell_a)
                    for (auto j : cell_b)
                        function(this->cell_lists.coordinates.row(i) -
                                 this->cell_lists.coordinates.row(j) - offset);
            }
        };

        template <typename Scalar, template<typename,size_t> typename GridType>
        class CellLists<Scalar, GridType>::Neighbourhood
        {
        public:
            Neighbourhood(const CellLists<Scalar, GridType>& cell_lists, size_t particle)
                : cell_lists(cell_lists), particle(particle)
            { }

            inline void for_each(std::function<void (Coordinate&&)> function) const
            {
                auto &cells = this->cell_lists;
                auto& cell = cells.grid[cells.particles[particle].cell];
                auto&& x = cells.coordinates.row(this->particle);

                for (auto neighbour : cell)
                {
                    if (neighbour == particle) continue;
                    function(x - cells.coordinates.row(neighbour));
                }

                Coordinate offset;
                for (auto& adjacent_cell : cell.neighbours)
                {
                    for (size_t c = 0; c < d; ++c)
                        offset(c) = adjacent_cell.image_offset(c) *
                            cells.cell_widths(c);

                    for (auto neighbour : adjacent_cell.cell->children)
                        function(x - cells.coordinates.row(neighbour) - offset);
                }
            }

            inline bool for_each_terminating(std::function<bool (Coordinate&&)> function) const
            {
                auto& cells = this->cell_lists;
                auto& cell = cells.grid[cells.particles[particle].cell];
                auto&& x = cells.coordinates.row(this->particle);
                bool abort;

                for (auto neighbour : cell)
                {
                    if (neighbour == this->particle) continue;
                    abort = function(x - cells.coordinates.row(neighbour));
                    if (abort) return true;
                }

                Coordinate offset;
                for (auto& adjacent_cell : cell.neighbours)
                {
                    for (size_t c = 0; c < d; ++c)
                        offset(c) = adjacent_cell.image_offset(c) *
                            cells.cell_widths(c);

                    for (auto neighbour : adjacent_cell.cell->children)
                    {
                        abort = function(x - cells.coordinates.row(neighbour) - offset);
                        if (abort) return true;
                    }
                }

                return false;
            }

        private:
            const CellLists<Scalar, GridType>& cell_lists;
            const size_t particle;
        };
    }
}

#endif
