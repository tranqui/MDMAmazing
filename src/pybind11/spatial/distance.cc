#include <iostream>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
namespace py = pybind11;

#include <mdma/spatial/cell-lists.h>
using namespace mdma;

// Define matrix types in row major format for compatibility with numpy.
namespace Eigen
{
    template <typename Scalar>
    using RowMatrixX = Matrix<Scalar, Dynamic, Dynamic, RowMajor>;

    template <typename Scalar>
    using RowVectorX = Matrix<Scalar, 1, Dynamic>;

    template <typename Scalar, size_t d>
    using Coordinates = Matrix<Scalar, Dynamic, d, RowMajor>;
}

/**
 *  Calculate the pairwise distances between points.
 *
 *  @param coords is an n by d matrix given n points in a d-dimensional space.
 *
 *  @return a numpy array giving the n*(n-1)/2 distances between points.
 **/
template <typename CoordType>
auto pdist(const Eigen::Ref<const CoordType>& coords)
{
    using Scalar = decltype(std::declval<CoordType>().sum());
    size_t n = coords.rows();

    py::array_t<Scalar> distances(n*(n-1)/2);
    auto data = static_cast<Scalar*>(distances.request().ptr);

    size_t k = 0;
    for (size_t i = 0; i < n-1; ++i)
        for (size_t j = i+1; j < n; ++j)
            data[k++] = (coords.row(i) - coords.row(j)).norm();

    return distances;
}

/**
 *  Calculate the histogram of pairwise distances between points.
 *
 *  @param coords is an n by d matrix given n points in a d-dimensional space.
 *  @param min is the start point of the histogram bin range.
 *  @param max is the final point of the histogram bin range.
 *  @param nbins is the number of histogram bins in the range.
 *
 *  @return a numpy (int) array giving the counts for each bin.
 **/
template <typename CoordType, typename Scalar>
py::array_t<int> pdist_histogram(const Eigen::Ref<const CoordType>& coords,
                                 Scalar min, Scalar max, size_t nbins)
{
    size_t n = coords.rows();
    Scalar min_squ = min*min;
    Scalar max_squ = max*max;
    Scalar bin_width = (max - min) / nbins;

    py::array_t<int> frequencies(nbins);
    auto data = static_cast<int*>(frequencies.request().ptr);
    for (size_t i = 0; i < nbins; ++i) data[i] = 0;

    for (size_t i = 0; i < n-1; ++i)
        for (size_t j = i+1; j < n; ++j)
        {
            Scalar dr_squ = (coords.row(i) - coords.row(j)).squaredNorm();
            if (dr_squ < min_squ or dr_squ > max_squ) continue;

            Scalar dr = std::sqrt(dr_squ);
            int bin = (dr - min)/bin_width;
            data[bin] += 1;
        }


    return frequencies;
}

/**
 *  Calculate the pairwise distances between points under periodic boundary conditions.
 *
 *  @param coords is an n by d matrix given n points in a d-dimensional space.
 *  @param box_dimensions is a d vector giving size of box in each dimensions for periodic wrapping.
 *
 *  @return a numpy array giving the n*(n-1)/2 distances between points.
 **/
template <typename CoordType, typename BoxType>
auto periodic_pdist(const Eigen::Ref<const CoordType>& coords,
                    const Eigen::Ref<const BoxType>& box_dimensions)
{
    using Scalar = decltype(std::declval<CoordType>().sum());
    size_t n = coords.rows();
    size_t d = coords.cols();

    py::array_t<Scalar> distances(n*(n-1)/2);
    auto data = static_cast<Scalar*>(distances.request().ptr);

    size_t k = 0;
    for (size_t i = 0; i < n-1; ++i)
        for (size_t j = i+1; j < n; ++j)
        {
            Eigen::RowVectorX<Scalar> delta = coords.row(i) - coords.row(j);
            for (size_t c = 0; c < d; ++c)
            {
                while (delta(c) >=  0.5*box_dimensions(c)) delta(c) -= box_dimensions(c);
                while (delta(c) <= -0.5*box_dimensions(c)) delta(c) += box_dimensions(c);
            }
            data[k++] = (coords.row(i) - coords.row(j)).norm();
        }

    return distances;
}

/**
 *  Calculate the histogram of pairwise distances between points under periodic boundary conditions.
 *
 *  @param coords is an n by d matrix given n points in a d-dimensional space.
 *  @param min is the start point of the histogram bin range.
 *  @param max is the final point of the histogram bin range.
 *  @param nbins is the number of histogram bins in the range.
 *
 *  @return a numpy (int) array giving the counts for each bin.
 **/
template <typename CoordType, typename BoxType, typename Scalar>
py::array_t<int> periodic_pdist_histogram(const Eigen::Ref<const CoordType>& coords,
                                          const Eigen::Ref<const BoxType>& box_dimensions,
                                          Scalar min, Scalar max, size_t nbins)
{
    size_t n = coords.rows();
    size_t d = coords.cols();
    Scalar min_squ = min*min;
    Scalar max_squ = max*max;
    Scalar bin_width = (max - min) / nbins;

    py::array_t<int> frequencies(nbins);
    auto data = static_cast<int*>(frequencies.request().ptr);
    for (size_t i = 0; i < nbins; ++i) data[i] = 0;

    for (size_t i = 0; i < n-1; ++i)
        for (size_t j = i+1; j < n; ++j)
        {
            Eigen::RowVectorX<Scalar> delta = coords.row(i) - coords.row(j);
            for (size_t c = 0; c < d; ++c)
            {
                while (delta(c) >=  0.5*box_dimensions(c)) delta(c) -= box_dimensions(c);
                while (delta(c) <= -0.5*box_dimensions(c)) delta(c) += box_dimensions(c);
            }
            Scalar dr_squ = delta.squaredNorm();
            if (dr_squ < min_squ or dr_squ > max_squ) continue;

            Scalar dr = std::sqrt(dr_squ);
            int bin = (dr - min)/bin_width;
            data[bin] += 1;
        }


    return frequencies;
}

/**
 *  Define distance functions in the module for a Scalar type specified
 *    by the Scalar template parameter.
 *
 *  @param m is the module object scope to define the functions inside.
 **/
template <typename Scalar, typename Module>
void define_functions(Module m)
{
    {
        using CoordType = Eigen::RowMatrixX<Scalar>;
        using BoxType = Eigen::RowVectorX<Scalar>;
        using FunctionType = py::array_t<Scalar> (*) (const Eigen::Ref<const CoordType>&);
        using PeriodicFunctionType = py::array_t<Scalar> (*)
            (const Eigen::Ref<const CoordType>&, const Eigen::Ref<const BoxType>&);

        m.def("pdist",
              static_cast<FunctionType>( &pdist<CoordType>),
              py::arg("X").noconvert(), py::return_value_policy::move,
              "Pairwise distances between vectors.");
        m.def("pdist_histogram", &pdist_histogram<CoordType,Scalar>,
              py::arg("X").noconvert(),
              py::arg("min"), py::arg("max"), py::arg("nbins"),
              py::return_value_policy::move,
              "Histogram the pairwise distances between vectors.");

        m.def("periodic_pdist",
              static_cast<PeriodicFunctionType>( &periodic_pdist<CoordType,BoxType>),
              py::arg("X").noconvert(),
              py::arg("box_dimensions").noconvert(),
              py::return_value_policy::move,
              "Pairwise distances between vectors with periodic boundary conditions.");
        m.def("periodic_pdist_histogram", &periodic_pdist_histogram<CoordType,BoxType,Scalar>,
              py::arg("X").noconvert(),
              py::arg("box_dimensions").noconvert(),
              py::arg("min"), py::arg("max"), py::arg("nbins"),
              py::return_value_policy::move,
              "Histogram the pairwise distances between vectors with periodic boundary conditions.");
    }
}

PYBIND11_MODULE(distance, module)
{
    module.doc() = "Pairwise distances module";

    define_functions<double>(module);
    define_functions<long double>(module);

    using Cells = spatial::CellLists<long double>;
    py::class_<Cells>(module, "CellLists")
        .def(py::init<long double,
             const Eigen::Ref<const Cells::Coordinate>&,
             const Eigen::Ref<const Cells::Domain>&>(),
             py::arg("neighbour_cutoff"),
             py::arg("box_dimensions").noconvert(),
             py::arg("coordinates").noconvert())
        .def_property_readonly("n", &Cells::n)
        .def_property_readonly("d", [](const Cells&)
                               { return spatial::details::DimensionsAtCompileTime; })
        .def_property_readonly("shape", &Cells::shape)
        .def_property_readonly("cell_dimensions", &Cells::cell_dimensions)
        .def_property_readonly("cells", &Cells::get_cells, py::return_value_policy::move)
        .def("find_displacements", &Cells::find_displacements, py::return_value_policy::move)
        .def("find_distances", &Cells::find_distances, py::return_value_policy::move);
}
