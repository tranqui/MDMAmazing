#include <iostream>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
namespace py = pybind11;

#include <mdma/spatial/cell-lists.h>
using namespace mdma;

/// Define irrational number pi for convenience.
constexpr long double pi = std::atan2(0, -1);

/// Define matrix types in row major format for compatibility with numpy.
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
 *  @param coords is an n by d matrix giving n points in a d-dimensional space.
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
 *  @param coords is an n by d matrix giving n points in a d-dimensional space.
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
 *  Calculate the displacement between points under periodic boundary conditions.
 *
 *  @param ri is the position of the first point.
 *  @param rj is the position of the second point.
 *  @param box_dimensions is a d vector giving size of box in each dimensions for periodic wrapping.
 *
 *  @return a row vector giving the displacement betwwen the points.
 **/
template <typename CoordType, typename BoxType>
inline auto periodic_displacement(CoordType&& ri, CoordType&& rj,
                                  const Eigen::Ref<const BoxType>& box_dimensions)
{
    using Scalar = typename CoordType::Scalar;
    Eigen::RowVectorX<Scalar> delta = ri - rj;
    for (int c = 0; c < box_dimensions.size(); ++c)
    {
        while (delta(c) >=  0.5*box_dimensions(c)) delta(c) -= box_dimensions(c);
        while (delta(c) <= -0.5*box_dimensions(c)) delta(c) += box_dimensions(c);
    }
    return delta;
}

/**
 *  Calculate the pairwise distances between points under periodic boundary conditions.
 *
 *  @param coords is an n by d matrix giving n points in a d-dimensional space.
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

    py::array_t<Scalar> distances(n*(n-1)/2);
    auto data = static_cast<Scalar*>(distances.request().ptr);

    size_t k = 0;
    for (size_t i = 0; i < n-1; ++i)
        for (size_t j = i+1; j < n; ++j)
        {
            auto delta = periodic_displacement(coords.row(i), coords.row(j), box_dimensions);
            data[k++] = delta.norm();
        }

    return distances;
}

/**
 *  Calculate the histogram of pairwise distances between points under periodic boundary conditions.
 *
 *  @param coords is an n by d matrix giving n points in a d-dimensional space.
 *  @param box_dimensions is a d vector giving size of box in each dimensions for periodic wrapping.
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
    Scalar min_squ = min*min;
    Scalar max_squ = max*max;
    Scalar bin_width = (max - min) / nbins;

    py::array_t<int> frequencies(nbins);
    auto data = static_cast<int*>(frequencies.request().ptr);
    for (size_t i = 0; i < nbins; ++i) data[i] = 0;

    for (size_t i = 0; i < n-1; ++i)
        for (size_t j = i+1; j < n; ++j)
        {
            auto delta = periodic_displacement(coords.row(i), coords.row(j), box_dimensions);
            Scalar dr_squ = delta.squaredNorm();
            if (dr_squ < min_squ or dr_squ > max_squ) continue;

            Scalar dr = std::sqrt(dr_squ);
            int bin = (dr - min)/bin_width;
            data[bin] += 1;
        }

    return frequencies;
}

/**
 *  Calculate the two body distribution function for points under periodic boundary conditions.
 *
 *  @param coords is an n by d matrix giving n points in a d-dimensional space.
 *  @param box_dimensions is a d vector giving size of box in each dimensions for periodic wrapping.
 *  @param min is the start point of the histogram bin range.
 *  @param max is the final point of the histogram bin range.
 *  @param nbins is the number of histogram bins in the range.
 *
 *  @return a numpy array giving the pdf for each bin.
 **/
template <typename CoordType, typename BoxType, typename Scalar>
py::array_t<Scalar> periodic_pairwise_correlations(const Eigen::Ref<const CoordType>& coords,
                                                   const Eigen::Ref<const BoxType>& box_dimensions,
                                                   Scalar min, Scalar max, size_t nbins)
{
    size_t n = coords.rows();
    Scalar bin_width = (max - min) / nbins;

    Scalar edges[nbins+1];
    Scalar left_edge = min;
    for (size_t i = 0; i < nbins+1; ++i)
    {
        edges[i] = left_edge;
        left_edge += bin_width;
    }

    Scalar volume = 1;
    for (int c = 0; c < box_dimensions.cols(); ++c)
        volume *= box_dimensions(c);
    Scalar density = n / volume;

    auto histogram = periodic_pdist_histogram(coords, box_dimensions, min, max, nbins);
    auto histogram_data = static_cast<int*>(histogram.request().ptr);
    py::array_t<Scalar> g(histogram.size());
    auto data = static_cast<Scalar*>(g.request().ptr);

    auto cube = [](Scalar x) { return x*x*x; };
    for (size_t i = 0; i < nbins; ++i)
    {
        long double dV = 4*pi/3*(cube(edges[i+1])-cube(edges[i]));
        data[i] = static_cast<Scalar>(2*histogram_data[i]) / (n*density*dV);
    }

    return g;
}

/**
 *  Calculate the histogram of triplets of mutual pairwise distances between points under
 *  periodic boundary conditions.
 *
 *  This is useful for measuring 3-body correlations in atomic systems.
 *
 *  @param coords is an n by d matrix giving n points in a d-dimensional space.
 *  @param box_dimensions is a d vector giving size of box in each dimensions for periodic wrapping.
 *  @param min is the start point of the histogram bin range.
 *  @param max is the final point of the histogram bin range.
 *  @param nbins is the number of histogram bins in the range.
 *
 *  @return a numpy (int) array giving the counts for each bin, in compressed form
 *   (i.e. the upper triangle of the cubic array).
 **/
template <typename CoordType, typename BoxType, typename Scalar>
py::array_t<int> periodic_triplet_dist_histogram(const Eigen::Ref<const CoordType>& coords,
                                                 const Eigen::Ref<const BoxType>& box_dimensions,
                                                 Scalar min, Scalar max, size_t nbins)
{
    size_t n = coords.rows();
    //size_t d = coords.cols();
    //Scalar min_squ = min*min;
    Scalar max_squ = max*max;
    Scalar bin_width = (max - min) / nbins;

    size_t N = nbins * (nbins + 1) * (nbins + 2) / 6;
    py::array_t<int> frequencies(N);
    auto data = static_cast<int*>(frequencies.request().ptr);
    for (size_t i = 0; i < N; ++i) data[i] = 0;

    auto index = [nbins] (int i, int j, int k) {
                     return i*(i-1)*(i-2)/6 + i*nbins*(2-i+nbins)/2 + (nbins-i+1)*(j-i) - (j-i)*(j-i+1)/2 + k - j;
                 };

    for (size_t i = 0; i < n-2; ++i)
    {
        for (size_t j = i+1; j < n-1; ++j)
        {
            auto delta_ij = periodic_displacement(coords.row(i), coords.row(j), box_dimensions);
            Scalar rij_squ = delta_ij.squaredNorm();
            if (rij_squ > max_squ) continue;

            Scalar rij = std::sqrt(rij_squ);
            int bij = 0;
            if (rij > min) bij = (rij - min)/bin_width;

            for (size_t k = j+1; k < n; ++k)
            {
                auto delta_ik = periodic_displacement(coords.row(i), coords.row(k), box_dimensions);
                Scalar rik_squ = delta_ik.squaredNorm();
                if (rik_squ > max_squ) continue;

                auto delta_jk = periodic_displacement(coords.row(j), coords.row(k), box_dimensions);
                Scalar rjk_squ = delta_jk.squaredNorm();
                if (rjk_squ > max_squ) continue;

                Scalar rik = std::sqrt(rik_squ);
                Scalar rjk = std::sqrt(rjk_squ);

                int bik = 0, bjk = 0;
                if (rik > min) bik = (rik - min)/bin_width;
                //else continue;
                if (rjk > min) bjk = (rjk - min)/bin_width;
                //else continue;

                int bij_copy = bij; // save this so the original is not swapped
                if (bij_copy > bik) std::swap(bij_copy,bik);
                if (bik > bjk) std::swap(bik,bjk);
                if (bij_copy > bik) std::swap(bij_copy,bik);

                data[index(bij_copy,bik,bjk)] += 1;
            }
        }
    }

    size_t l = 0;
    for (size_t i = 0; i < nbins; ++i)
        for (size_t j = i; j < nbins; ++j)
            for (size_t k = j; k < nbins; ++k)
            {
                if (i == j and i == k) data[l] *= 6;
                else if (i == j or i == k or j == k) data[l] *= 2;
                l++;
            }

    return frequencies;
}

/**
 *  Calculate the three body distribution function for points under periodic boundary conditions.
 *
 *  @param coords is an n by d matrix giving n points in a d-dimensional space.
 *  @param box_dimensions is a d vector giving size of box in each dimensions for periodic wrapping.
 *  @param min is the start point of the histogram bin range.
 *  @param max is the final point of the histogram bin range.
 *  @param nbins is the number of histogram bins in the range.
 *
 *  @return a numpy array giving the pdf for each bin, in compressed form
 *   (i.e. the upper triangle of the cubic array).
 **/
template <typename CoordType, typename BoxType, typename Scalar>
py::array_t<Scalar> periodic_triplet_correlations(const Eigen::Ref<const CoordType>& coords,
                                                  const Eigen::Ref<const BoxType>& box_dimensions,
                                                  Scalar min, Scalar max, size_t nbins)
{
    size_t n = coords.rows();
    Scalar bin_width = (max - min) / nbins;

    Scalar r[nbins];
    Scalar left_edge = min;
    for (size_t i = 0; i < nbins; ++i)
    {
        r[i] = left_edge + 0.5*bin_width;
        left_edge += bin_width;
    }

    Scalar volume = 1;
    for (int c = 0; c < box_dimensions.cols(); ++c)
        volume *= box_dimensions(c);
    Scalar density = n / volume;

    auto histogram = periodic_triplet_dist_histogram(coords, box_dimensions, min, max, nbins);
    auto histogram_data = static_cast<int*>(histogram.request().ptr);
    py::array_t<Scalar> g(histogram.size());
    auto data = static_cast<Scalar*>(g.request().ptr);

    size_t l = 0;
    for (size_t i = 0; i < nbins; ++i)
        for (size_t j = i; j < nbins; ++j)
            for (size_t k = j; k < nbins; ++k)
            {
                Scalar dV = 8*pi*pi * r[i] * r[j] * r[k] * bin_width*bin_width*bin_width;
                data[l] = static_cast<Scalar>(histogram_data[l]) / (n*density*density*dV);
                l++;
            }

    return g;
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
        m.def("periodic_pairwise_correlations",
              &periodic_pairwise_correlations<CoordType,BoxType,Scalar>,
              py::arg("X").noconvert(),
              py::arg("box_dimensions").noconvert(),
              py::arg("min"), py::arg("max"), py::arg("nbins"),
              py::return_value_policy::move,
              "Two body distribution function for points under periodic boundary conditions.");

        m.def("periodic_triplet_dist_histogram",
              &periodic_triplet_dist_histogram<CoordType,BoxType,Scalar>,
              py::arg("X").noconvert(),
              py::arg("box_dimensions").noconvert(),
              py::arg("min"), py::arg("max"), py::arg("nbins"),
              py::return_value_policy::move,
              "Histogram for triplets of distances between vectors with periodic boundary conditions.");
        m.def("periodic_triplet_correlations",
              &periodic_triplet_correlations<CoordType,BoxType,Scalar>,
              py::arg("X").noconvert(),
              py::arg("box_dimensions").noconvert(),
              py::arg("min"), py::arg("max"), py::arg("nbins"),
              py::return_value_policy::move,
              "Three body distribution function for points under periodic boundary conditions.");
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
