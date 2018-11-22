#ifndef _MDMA_CELL_LISTS_H
#define _MDMA_CELL_LISTS_H

#include <list>
#include "grid.h"

namespace mdma
{
    namespace details
    {
        static constexpr size_t SpatialDimensions = 3;
        static constexpr int NumberOfPoints = Eigen::Dynamic;

        template <typename Coordinate>
        struct Pair
        {
            const Coordinate& a;
            const Coordinate& b;
        };

        template <typename Scalar>
        class NeighbourIter
        {
        public:
            /**
             * Main constructor is fairly self-explanatory.
             */
            inline ConstCoordsIter(const Vector& coords, size_t position=0)
                : coords(coords), current_position(position)
            { }

            /**
             * Copy constructor is trivial.
             */
            inline ConstCoordsIter(const ConstCoordsIter&) = default;

            /**
             * Value referencing: returns the row of the matrix, which is the coordinate
             *  of the current particle.
             */
            inline const Pair<Scalar>& operator*() const
            {
                return block<d>(this->coords, this->current_position);
            }

            /**
             * Test for equality, for loop stopping criteria etc.
             */
            inline bool operator==(const ConstCoordsIter& ref) const
            {
                return this->current_position == ref.current_position;
            }

            /**
             * Test for inequality, for loop stopping criteria etc.
             */
            inline bool operator!=(const ConstCoordsIter& ref) const
            {
                return !(*this == ref);
            }

            /**
             * Increment the operator in the prefix case, i.e. ++iter;
             */
            inline ConstCoordsIter& operator++()
            {
                ++this->current_position;
                return *this;
            }

            /**
             * Increment the operator in the postfix case, i.e. iter++;
             * NB: the (int) parameter signals this is a postfix, an awkward C++ syntax.
             */
            inline ConstCoordsIter operator++( int )
            {
                ConstCoordsIter clone(*this);
                ++this->current_position;
                return clone;
            }

        private:
            // The reference coordinate data we iterate through.
            const Vector& coords;
            // The current particle within said coordinates.
            size_t current_position;
        };
    }

    template <typename Scalar>
    class CellListMap :
        public CartesianBasis<Scalar, SpatialDimensions, NumberOfPoints>
    {
    public:
        using Base = CartesianBasis<Scalar, SpatialDimensions, NumberOfPoints>;
        using typename Base::Point;
        using typename Base::Axis;
        using typename Base::Representation;

        using ArrayIndex = std::array<size_t, SpatialDimensions>;
        using EigenIndex = Eigen::Matrix<size_t, 1, SpatialDimensions>;

        static constexpr size_t d = details::SpatialDimensions;

        CellListMap(Scalar neighbour_cutoff,
                    Eigen::Ref<const Point> box_dimensions,
                    Eigen::Ref<const Representation> coordinates)
            : neighbour_cutoff(neighbour_cutoff),
              cells(optimal_cells(neighbour_cutoff, box_dimensions)),
              coordinates(coordinates)
        {
            for (size_t c = 0; c < SpatialDimensions; ++c)
                cell_widths[c] = box_dimensions[c] / cells.shape()[0];

            for (int atom = 0; atom < coordinates.rows(); ++atom)
            {
                const auto& point = coordinates.row(atom);

                EigenIndex cell_index;
                for (size_t c = 0; c < SpatialDimensions; ++c)
                    cell_index[c] = std::floor(point[c] / this->cell_widths[c]);

                std::list<size_t>& list = this->cells(cell_index);
                list.push_back(atom);
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
        EigenIndex
        optimal_cells(Scalar neighbour_cutoff, Eigen::Ref<const Point> box_dimensions)
        {
            EigenIndex ncells;
            for (size_t c = 0; c < SpatialDimensions; ++c)
                ncells[c] = std::floor(box_dimensions[c] / neighbour_cutoff);

            return ncells;
        }

        //std::
        // friend std::ostream& operator<< (std::ostream& stream, const CellList& cells)
        // {
        // }
        // void assignCell(long i);
        // void assignAllCells();
        // void rebuild();
        // void forAllPotentialNeighbours() //long i, Neighbour& operation)
        // {
        //     for (size_t i = 0; i < this->
        // }
        //void updateCell(long i, &cell);

        int neighbour_check(const std::list<size_t>& cell_a,
                            const std::list<size_t>& cell_b) const
        {
            int count = 0;

            for (auto i : cell_a)
            {
                for (auto j : cell_b)
                {
                    Scalar dr = Base::distance(this->coordinates.row(i),
                                               this->coordinates.row(j));
                    if (dr < this->neighbour_cutoff) count++;
                }
            }

            return count;
        }

        int neighbour_check(const std::list<size_t>& cell) const
        {
            int count = 0;

            auto it = std::begin(cell);
            auto end = std::end(cell);
            if (it == end) return count;

            for (; it != std::prev(end); ++it)
            {
                for (auto it2 = std::next(it); it2 != cell.end(); ++it2)
                {
                    Scalar dr = Base::distance(this->coordinates.row(*it),
                                               this->coordinates.row(*it2));
                    if (dr < this->neighbour_cutoff) count++;
                }
            }

            return count;
        }

        void find_neighbours()
        {
            // for (size_t i = 0; i < this->cells.size(); ++i)
            // {
            //     //cell : this->
            // }

            size_t num_neighbours = std::pow(2,SpatialDimensions) - 1;
            std::cout << this->cells.stride() << std::endl;


            std::list<EigenIndex> neighbours;
            char mask = 1;
            for (size_t i = 1; i <= num_neighbours; ++i)
            {
                neighbours.push_back(EigenIndex());
                auto& index = neighbours.back();

                for (size_t c = 0; c < SpatialDimensions; ++c)                    
                    index[c] = (i >> c) & mask;

                std::cout << i << ":" << index << ":" << this->cells.offset(index) << std::endl;
            }

            size_t index = 0;
            for (auto cell : this->cells)
            {
                EigenIndex i1 = this->cells.index_from_offset(index);
                for (auto& nb : neighbours)
                {
                    EigenIndex i2 = i1 + nb;
                    size_t offset = this->cells.offset(i2);
                    if (offset > this->cells.size())
                        offset = offset % this->cells.size();
                    EigenIndex i3 = this->cells.index_from_offset(offset);
                    std::cout << i1 << ":" << i2 << ":" << i3 << std::endl;
                }
                //int n = cell.size();
                //std::cout << this->cells.index_from_offset(index) << ": "
                //          << this->neighbour_check(cell) << ": "
                //          << n << " " << n*(n-1)/2 << "\n";
                // std::cout << index << ":\t"
                //           << this->cells.index_from_offset(index) << "\t:"
                //           << cell.size() << std::endl;
                index += 1;
            }
        }

        Point cell_dimensions() const
        {
            return this->cell_widths;
        }

    protected:
        Scalar neighbour_cutoff;
        Grid<std::list<size_t>, SpatialDimensions> cells;
        Representation coordinates;

        Point cell_widths;
    };


    template <typename Scalar>
    class CellList : public CellListMap<Scalar>
    {
    public:
        using Base = CellListMap<Scalar, SpatialDimensions, NumberOfPoints>;
        using typename Base::Point;
        using typename Base::Axis;
        using typename Base::Representation;

        CellList(Scalar neighbour_cutoff,
                 Eigen::Ref<const Point> box_dimensions,
                 Eigen::Ref<const Representation> coordinates)
            : Base(neighbour_cutoff, box_dimensions, coordinates),
              coordinates(coordinates)
        {
        }

    protected:
        Representation coordinates;
    };
}

#endif
