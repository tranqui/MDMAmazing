#include <iostream>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
namespace py = pybind11;

#include <random>
#include <mdma/spatial/cell-lists.h>
using namespace mdma;

template <typename Scalar>
class HardSphereMonteCarlo : public spatial::CellLists<Scalar>
{
public:
    constexpr static Scalar particle_diameter = 1;
    constexpr static Scalar particle_diameter_squ = particle_diameter*particle_diameter;

    using Base = spatial::CellLists<Scalar>;
    using typename Base::Coordinate;
    using typename Base::Domain;
    using typename Base::EigenIndex;
    using Base::d;

    HardSphereMonteCarlo(Scalar step_size,
                         const Eigen::Ref<const Coordinate>& box_dimensions,
                         const Eigen::Ref<const Domain>& coordinates)
        : Base(particle_diameter, box_dimensions, coordinates),
          random_generator(std::random_device()()),
          random_particle_distribution(0, this->n()-1),
          random_particle( std::bind(random_particle_distribution,
                                     this->random_generator) ),
          random_step_distribution(-step_size, step_size),
          random_step( std::bind(this->random_step_distribution,
                                 this->random_generator) ),
          total_steps(0), successful_steps(0)
    {
    }

    inline Scalar acceptance_ratio() const
    {
        return static_cast<Scalar>(this->successful_steps) / this->total_steps;
    }

    inline void step(int atom)
    {
        auto&& x = this->coordinates.row(atom);
        Coordinate x0(x);
        this->total_steps++;

        // Step in a random direction.
        for (size_t c = 0; c < d; ++c)
        {
            x(c) += this->random_step();
            if (x(c) > this->box_dimensions(c)) x(c) -= this->box_dimensions(c);
            else if (x(c) < 0) x(c) += this->box_dimensions(c);
        }

        // Update the particle cell if it has changed.
        EigenIndex index;
        for (size_t c = 0; c < d; ++c) index[c] = std::floor(x[c] / this->cell_widths[c]);
        auto new_cell = this->grid.offset(index);
        auto old_cell = this->particles[atom].cell;
        if (new_cell != old_cell)
        {
            this->particles[atom].cell = new_cell;
            this->grid[old_cell].children.remove(atom);
            this->grid[new_cell].children.push_back(atom);
        }

        // Check for any overlaps.
        auto does_overlap =
            [&](Coordinate&& delta)
            {
                return (delta.squaredNorm() < particle_diameter_squ);
            };
        bool overlap = this->neighbourhood(atom).for_each_terminating(does_overlap);

        // Metropolis rule for hard spheres: undo any move which causes an overlap.
        if (overlap)
        {
            x = x0;
            if (new_cell != old_cell)
            {
                this->particles[atom].cell = old_cell;
                this->grid[new_cell].children.remove(atom);
                this->grid[old_cell].children.push_back(atom);
            }

            return;
        }

        this->successful_steps++;
    }

    inline void sweep()
    {
        for (size_t i = 0; i < this->n(); ++i)
            this->step( this->random_particle() );
    }

    inline void run(int num_sweeps)
    {
        for (int sweep = 0; sweep < num_sweeps; ++sweep) this->sweep();
    }

    Domain get_coordinates() const
    {
        return this->coordinates;
    }

protected:
    std::default_random_engine random_generator;

    std::uniform_int_distribution<> random_particle_distribution;
    std::function<int (void)> random_particle;

    std::uniform_real_distribution<Scalar> random_step_distribution;
    std::function<Scalar (void)> random_step;

private:
    size_t total_steps, successful_steps;
};

PYBIND11_MODULE(montecarlo, module)
{
    module.doc() = "Hard sphere monte carlo module.";
    py::module::import("mdma.spatial.distance");

    using MonteCarlo = HardSphereMonteCarlo<long double>;
    py::class_<MonteCarlo>(module, "HardSphereMonteCarlo")
        .def(py::init<long double,
             const Eigen::Ref<const MonteCarlo::Coordinate>&,
             const Eigen::Ref<const MonteCarlo::Domain>&>(),
             py::arg("step_size"),
             py::arg("box_dimensions").noconvert(),
             py::arg("coordinates").noconvert())
        .def_property_readonly("n", &MonteCarlo::n)
        .def_property_readonly("d", [](const MonteCarlo&)
                               { return spatial::details::DimensionsAtCompileTime; })
        .def_property_readonly("shape", &MonteCarlo::shape)
        .def_property_readonly("cell_dimensions", &MonteCarlo::cell_dimensions)
        .def_property_readonly("cells", &MonteCarlo::get_cells, py::return_value_policy::move)
        .def("find_displacements", &MonteCarlo::find_displacements, py::return_value_policy::move)
        .def("find_distances", &MonteCarlo::find_distances, py::return_value_policy::move)
        .def("run", &MonteCarlo::run)
        .def_property_readonly("coordinates", &MonteCarlo::get_coordinates, py::return_value_policy::move)
        .def_property_readonly("acceptance_ratio", &MonteCarlo::acceptance_ratio);
}
