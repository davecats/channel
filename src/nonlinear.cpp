#include "channel/nonlinear.hpp"

#include "channel/memory.hpp"
#include "channel/runtime.hpp"

#include <Kokkos_Profiling_ScopedRegion.hpp>

#include <array>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#include <mpi.h>

namespace channel {
namespace {

KOKKOS_INLINE_FUNCTION constexpr int derivative_index(int y, int order, int offset) {
  return (y * DnsNonlinearVelocityRhsStage::derivative_orders + order) *
             DnsNonlinearVelocityRhsStage::derivative_stencil +
         (offset + 2);
}

void check_span_size(std::size_t got, std::size_t expected, const char* name) {
  if (got != expected) {
    throw std::runtime_error(std::string(name) + " size mismatch");
  }
}

bool debug_line_enabled() {
  const char* value = std::getenv("CHANNEL_DEBUG_TIMESTEP");
  return value != nullptr && std::string(value).find("line") != std::string::npos;
}

void debug_rhs_line(const char* label,
                    const DnsGrid& grid,
                    const DeviceVector<Complex>& eta_rhs,
                    const DeviceVector<Complex>& d2v_rhs) {
  if (!debug_line_enabled()) return;
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  const int z_storage = 14;
  const int x_storage = 5;
  if (z_storage < 0 || z_storage >= grid.nz || x_storage < 0 || x_storage >= grid.nx) return;

  const auto eta = eta_rhs.copy_to_host();
  const auto d2v = d2v_rhs.copy_to_host();
  const int lines = static_cast<int>(grid.line_count());
  const int line = z_storage * grid.nx + x_storage;
  for (int rank_to_print = 0; rank_to_print < size; ++rank_to_print) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank != rank_to_print) continue;
    for (int y = 0; y < grid.ny; ++y) {
      const int global_y = grid.y_first + y;
      const auto p = static_cast<std::size_t>(y * lines + line);
      std::cout << std::scientific << std::setprecision(16)
                << "CPP_RHSLINE " << label
                << " rank=" << rank
                << " y=" << global_y
                << " ix=" << x_storage
                << " iz=" << -3
                << " (" << eta[p].real() << "," << eta[p].imag() << ")"
                << " (" << d2v[p].real() << "," << d2v[p].imag() << ")"
                << '\n';
    }
    std::cout.flush();
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void debug_init_d2v_terms(const DnsGrid& grid,
                          const DnsNonlinearVelocityRhsConfig& cfg,
                          const DnsState& state,
                          const DeviceVector<double>& derivatives,
                          const DeviceVector<double>& k2_values,
                          const DeviceVector<Complex>& old_d2v_rhs) {
  if (!debug_line_enabled()) return;
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank != 0) return;

  const int global_y = 12;
  const int y = global_y - grid.y_first;
  const int z = 14;
  const int x = 5;
  if (y < grid.active_y_storage_first() || y >= grid.active_y_storage_last_exclusive() ||
      z < 0 || z >= grid.nz || x < 0 || x >= grid.nx) {
    return;
  }

  const auto v = state.component_host(cfg.v_component);
  const auto der = derivatives.copy_to_host();
  const auto k2_host = k2_values.copy_to_host();
  const auto old = old_d2v_rhs.copy_to_host();
  const int lines = static_cast<int>(grid.line_count());
  const int line = z * grid.nx + x;
  const int p = y * lines + line;
  const int dy = y - grid.active_y_storage_first();

  Complex d0v(0.0, 0.0);
  Complex d2v_value(0.0, 0.0);
  Complex d4v_value(0.0, 0.0);
  for (int offset = -2; offset <= 2; ++offset) {
    const int yy = y + offset;
    if (yy < 0 || yy >= grid.ny) continue;
    const auto q = static_cast<std::size_t>(yy * lines + line);
    d0v += der[static_cast<std::size_t>(derivative_index(dy, 0, offset))] * v[q];
    d2v_value += der[static_cast<std::size_t>(derivative_index(dy, 2, offset))] * v[q];
    d4v_value += der[static_cast<std::size_t>(derivative_index(dy, 3, offset))] * v[q];
  }

  const double k2_line = k2_host[static_cast<std::size_t>(line)];
  const Complex unkn = d2v_value - k2_line * d0v;
  const Complex implicit =
      cfg.viscosity * (d4v_value - 2.0 * k2_line * d2v_value + k2_line * k2_line * d0v);
  const Complex rhs =
      cfg.implicit_weight / cfg.dt * unkn + implicit - cfg.history_weight * old[static_cast<std::size_t>(p)];
  std::cout << std::scientific << std::setprecision(16)
            << "CPP_INIT_D2V y=" << global_y
            << " ix=" << x
            << " iz=" << -3
            << " (" << unkn.real() << "," << unkn.imag() << ")"
            << " (" << implicit.real() << "," << implicit.imag() << ")"
            << " (" << old[static_cast<std::size_t>(p)].real() << "," << old[static_cast<std::size_t>(p)].imag() << ")"
            << " (" << rhs.real() << "," << rhs.imag() << ")"
            << '\n';
}

} // namespace

void DnsNonlinearProductTransformStage::prepare(const DnsState& state, DnsNonlinearProductTransformConfig cfg) {
  const auto& grid = state.grid();
  grid.validate();
  if (grid.values_per_component() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::runtime_error("DnsNonlinearProductTransformStage component size exceeds int range");
  }
  grid_ = grid;
  cfg_ = std::move(cfg);

  check_component(cfg_.u_component, "u_component");
  check_component(cfg_.v_component, "v_component");
  check_component(cfg_.w_component, "w_component");
  for (const int component : cfg_.product_components) {
    check_component(component, "product_component");
  }
  if ((cfg_.dealiased_physical_x == 0) != (cfg_.dealiased_physical_z == 0)) {
    throw std::runtime_error("dealiased nonlinear transform requires both physical x and z sizes");
  }
  if (cfg_.dealiased_physical_x != 0) {
    if (cfg_.dealiased_physical_x < grid.nx || cfg_.dealiased_physical_z < grid.nz ||
        (cfg_.dealiased_physical_x % 2) != 0) {
      throw std::runtime_error("dealiased nonlinear transform received invalid padded dimensions");
    }
    dealiased_x_half_ = cfg_.dealiased_physical_x / 2 + 1;
    if (cfg_.distributed_dealiased_fft) {
      distributed_fft_ = std::make_unique<DistributedDealiasedFft2DPlan>();
      distributed_fft_->configure(grid_,
                                  {cfg_.distributed_spectral_x_total,
                                   cfg_.distributed_spectral_x_first,
                                   cfg_.dealiased_physical_x,
                                   cfg_.dealiased_physical_z,
                                   cfg_.distributed_npxz,
                                   cfg_.distributed_ipxz,
                                   cfg_.distributed_comm_x});
      const int z_count = distributed_fft_->z_count();
      initialize_memory_pool();
      distributed_u_ = RealView3D("nonlinear_distributed_u", grid_.ny, z_count, cfg_.dealiased_physical_x);
      distributed_v_ = RealView3D("nonlinear_distributed_v", grid_.ny, z_count, cfg_.dealiased_physical_x);
      distributed_w_ = RealView3D("nonlinear_distributed_w", grid_.ny, z_count, cfg_.dealiased_physical_x);
      distributed_product_ =
          RealView3D("nonlinear_distributed_product", grid_.ny, z_count, cfg_.dealiased_physical_x);
      prepared_ = true;
      return;
    }
    initialize_memory_pool();
    dealiased_spectrum_ =
        ComplexView3D("nonlinear_dealiased_spectrum", grid_.ny, cfg_.dealiased_physical_z, dealiased_x_half_);
    dealiased_z_spectrum_ =
        ComplexView3D("nonlinear_dealiased_z_spectrum", grid_.ny, cfg_.dealiased_physical_z, dealiased_x_half_);
    dealiased_u_ = RealView3D("nonlinear_dealiased_u", grid_.ny, cfg_.dealiased_physical_z, cfg_.dealiased_physical_x);
    dealiased_v_ = RealView3D("nonlinear_dealiased_v", grid_.ny, cfg_.dealiased_physical_z, cfg_.dealiased_physical_x);
    dealiased_w_ = RealView3D("nonlinear_dealiased_w", grid_.ny, cfg_.dealiased_physical_z, cfg_.dealiased_physical_x);
    dealiased_product_ =
        RealView3D("nonlinear_dealiased_product", grid_.ny, cfg_.dealiased_physical_z, cfg_.dealiased_physical_x);
    ExecutionSpace exec;
    if (cfg_.use_dealiased_2d_fft) {
      dealiased_spectrum_xzy_ =
          ComplexView3D("nonlinear_dealiased_spectrum_xzy", dealiased_x_half_, cfg_.dealiased_physical_z, grid_.ny);
      dealiased_u_xzy_ = RealView3D("nonlinear_dealiased_u_xzy",
                                    cfg_.dealiased_physical_x,
                                    cfg_.dealiased_physical_z,
                                    grid_.ny);
      dealiased_v_xzy_ = RealView3D("nonlinear_dealiased_v_xzy",
                                    cfg_.dealiased_physical_x,
                                    cfg_.dealiased_physical_z,
                                    grid_.ny);
      dealiased_w_xzy_ = RealView3D("nonlinear_dealiased_w_xzy",
                                    cfg_.dealiased_physical_x,
                                    cfg_.dealiased_physical_z,
                                    grid_.ny);
      dealiased_product_xzy_ = RealView3D("nonlinear_dealiased_product_xzy",
                                          cfg_.dealiased_physical_x,
                                          cfg_.dealiased_physical_z,
                                          grid_.ny);
      const KokkosFFT::axis_type<2> axes = {1, 0};
      dealiased_2d_forward_plan_ = std::make_unique<RealToComplexFft2DPlan>(
          exec, dealiased_product_xzy_, dealiased_spectrum_xzy_, KokkosFFT::Direction::forward, axes);
      dealiased_2d_inverse_plan_ = std::make_unique<ComplexToRealFft2DPlan>(
          exec, dealiased_spectrum_xzy_, dealiased_u_xzy_, KokkosFFT::Direction::backward, axes);
    } else {
      dealiased_z_inverse_plan_ = std::make_unique<ComplexFft1DPlan>(
          exec, dealiased_spectrum_, dealiased_z_spectrum_, KokkosFFT::Direction::backward, 1);
      dealiased_x_forward_plan_ = std::make_unique<RealToComplexFft1DPlan>(
          exec, dealiased_product_, dealiased_z_spectrum_, KokkosFFT::Direction::forward, 2);
      dealiased_x_inverse_plan_ = std::make_unique<ComplexToRealFft1DPlan>(
          exec, dealiased_z_spectrum_, dealiased_u_, KokkosFFT::Direction::backward, 2);
      dealiased_z_forward_plan_ = std::make_unique<ComplexFft1DPlan>(
          exec, dealiased_z_spectrum_, dealiased_spectrum_, KokkosFFT::Direction::forward, 1);
    }
  }

  velocity_inverse_ffts_.clear();
  if (cfg_.enable_velocity_inverse_ffts) {
    for (const int component : {cfg_.u_component, cfg_.v_component, cfg_.w_component}) {
      auto plan = std::make_unique<DnsComponentFftPlan>();
      plan->configure(state, {component, cfg_.velocity_inverse_normalization});
      velocity_inverse_ffts_.push_back(std::move(plan));
    }
  }

  velocity_transposes_.clear();
  velocity_transposes_.reserve(cfg_.velocity_transposes.size());
  for (const auto& transpose_cfg : cfg_.velocity_transposes) {
    check_component(transpose_cfg.component, "velocity_transpose_component");
    auto plan = std::make_unique<DnsComponentTransposePlan>();
    plan->configure(state, transpose_cfg);
    velocity_transposes_.push_back(std::move(plan));
  }

  product_forward_ffts_.clear();
  if (cfg_.enable_product_forward_ffts) {
    for (const int component : cfg_.product_components) {
      auto plan = std::make_unique<DnsComponentFftPlan>();
      plan->configure(state, {component, cfg_.product_forward_normalization});
      product_forward_ffts_.push_back(std::move(plan));
    }
  }

  product_transposes_.clear();
  product_transposes_.reserve(cfg_.product_transposes.size());
  for (const auto& transpose_cfg : cfg_.product_transposes) {
    check_component(transpose_cfg.component, "product_transpose_component");
    auto plan = std::make_unique<DnsComponentTransposePlan>();
    plan->configure(state, transpose_cfg);
    product_transposes_.push_back(std::move(plan));
  }

  prepared_ = true;
}

void DnsNonlinearProductTransformStage::apply(DnsState& state) {
  check_state(state, "DnsNonlinearProductTransformStage::apply");
  if (cfg_.dealiased_physical_x != 0) {
    if (cfg_.distributed_dealiased_fft) {
      apply_distributed_dealiased(state);
      return;
    }
    if (cfg_.use_dealiased_2d_fft) {
      apply_dealiased_2d(state);
      return;
    }
    apply_dealiased(state);
    return;
  }
  {
    Kokkos::Profiling::ScopedRegion region("transform_to_physical inverse_velocity_ffts");
    for (auto& plan : velocity_inverse_ffts_) {
      plan->execute(state, FftDirection::Inverse, "dns_nonlinear_velocity_inverse_fft");
    }
  }
  {
    Kokkos::Profiling::ScopedRegion region("transform_to_physical velocity_transpose");
    for (auto& plan : velocity_transposes_) {
      plan->execute(state);
    }
  }
  {
    Kokkos::Profiling::ScopedRegion region("transform_back build_products");
    build_products(state);
  }
  if (cfg_.enable_velocity_inverse_ffts) {
    Kokkos::Profiling::ScopedRegion region("transform_to_physical restore_velocity_ffts");
    for (auto& plan : velocity_inverse_ffts_) {
      plan->execute(state, FftDirection::Forward, "dns_nonlinear_velocity_restore_fft");
    }
  }
  {
    Kokkos::Profiling::ScopedRegion region("transform_back product_forward_ffts");
    for (auto& plan : product_forward_ffts_) {
      plan->execute(state, FftDirection::Forward, "dns_nonlinear_product_forward_fft");
    }
  }
  {
    Kokkos::Profiling::ScopedRegion region("transform_back product_transpose");
    for (auto& plan : product_transposes_) {
      plan->execute(state);
    }
  }
}

void DnsNonlinearProductTransformStage::check_state(const DnsState& state, const char* caller) const {
  if (!prepared_) {
    throw std::runtime_error(std::string(caller) + " called before prepare");
  }
  const auto& grid = state.grid();
  if (grid.nx != grid_.nx || grid.ny != grid_.ny || grid.nz != grid_.nz || grid.components != grid_.components) {
    throw std::runtime_error(std::string(caller) + " grid mismatch");
  }
}

void DnsNonlinearProductTransformStage::check_component(int component, const char* name) const {
  if (component < 0 || component >= grid_.components) {
    throw std::runtime_error(std::string("DnsNonlinearProductTransformStage ") + name + " out of range");
  }
}

void DnsNonlinearProductTransformStage::build_products(DnsState& state) {
  const auto u = state.component_view_3d(cfg_.u_component);
  const auto v = state.component_view_3d(cfg_.v_component);
  const auto w = state.component_view_3d(cfg_.w_component);
  const auto uu = state.component_view_3d(cfg_.product_components[static_cast<int>(VelocityProduct::UU)]);
  const auto vv = state.component_view_3d(cfg_.product_components[static_cast<int>(VelocityProduct::VV)]);
  const auto ww = state.component_view_3d(cfg_.product_components[static_cast<int>(VelocityProduct::WW)]);
  const auto uv = state.component_view_3d(cfg_.product_components[static_cast<int>(VelocityProduct::UV)]);
  const auto vw = state.component_view_3d(cfg_.product_components[static_cast<int>(VelocityProduct::VW)]);
  const auto uw = state.component_view_3d(cfg_.product_components[static_cast<int>(VelocityProduct::UW)]);
  const double factor = cfg_.product_factor;

  Kokkos::parallel_for(
      "dns_nonlinear_product_transform_build_products",
      Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace>({0, 0, 0}, {grid_.ny, grid_.nz, grid_.nx}),
      KOKKOS_LAMBDA(const int y, const int z, const int x) {
        const double ur = u(y, z, x).real();
        const double vr = v(y, z, x).real();
        const double wr = w(y, z, x).real();
        uu(y, z, x) = Complex(factor * ur * ur, 0.0);
        vv(y, z, x) = Complex(factor * vr * vr, 0.0);
        ww(y, z, x) = Complex(factor * wr * wr, 0.0);
        uv(y, z, x) = Complex(factor * ur * vr, 0.0);
        vw(y, z, x) = Complex(factor * vr * wr, 0.0);
        uw(y, z, x) = Complex(factor * ur * wr, 0.0);
      });
  fence("dns_nonlinear_product_transform_build_products");
}

void DnsNonlinearProductTransformStage::pack_dealiased_spectrum(DnsState& state, int component_id) {
  const auto component = state.component_view_3d(component_id);
  auto spectrum = dealiased_spectrum_;
  Kokkos::deep_copy(spectrum, Complex(0.0, 0.0));
  const int retained_z = (grid_.nz - 1) / 2;
  const int storage_z = grid_.nz;
  const int padded_z = cfg_.dealiased_physical_z;
  Kokkos::parallel_for(
      "dns_nonlinear_pack_dealiased_spectrum",
      Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace>({0, 0, 0}, {grid_.ny, grid_.nz, grid_.nx}),
      KOKKOS_LAMBDA(const int y, const int z, const int x) {
        const int signed_z = z <= retained_z ? z : z - storage_z;
        const int padded_index_z = signed_z >= 0 ? signed_z : signed_z + padded_z;
        spectrum(y, padded_index_z, x) = component(y, z, x);
      });
  fence("dns_nonlinear_pack_dealiased_spectrum");
}

void DnsNonlinearProductTransformStage::unpack_dealiased_product(DnsState& state, int component_id) {
  auto component = state.component_view_3d(component_id);
  auto spectrum = dealiased_spectrum_;
  const int retained_z = (grid_.nz - 1) / 2;
  const int storage_z = grid_.nz;
  const int padded_z = cfg_.dealiased_physical_z;
  Kokkos::parallel_for(
      "dns_nonlinear_unpack_dealiased_product",
      Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace>({0, 0, 0}, {grid_.ny, grid_.nz, grid_.nx}),
      KOKKOS_LAMBDA(const int y, const int z, const int x) {
        const int signed_z = z <= retained_z ? z : z - storage_z;
        const int padded_index_z = signed_z >= 0 ? signed_z : signed_z + padded_z;
        component(y, z, x) = spectrum(y, padded_index_z, x);
      });
  fence("dns_nonlinear_unpack_dealiased_product");
}

void DnsNonlinearProductTransformStage::pack_dealiased_spectrum_xzy(DnsState& state, int component_id) {
  const auto component = state.component_view_3d(component_id);
  auto spectrum = dealiased_spectrum_xzy_;
  Kokkos::deep_copy(spectrum, Complex(0.0, 0.0));
  const int retained_z = (grid_.nz - 1) / 2;
  const int storage_z = grid_.nz;
  const int padded_z = cfg_.dealiased_physical_z;
  Kokkos::parallel_for(
      "dns_nonlinear_pack_dealiased_spectrum_xzy",
      Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace>({0, 0, 0}, {grid_.nx, grid_.nz, grid_.ny}),
      KOKKOS_LAMBDA(const int x, const int z, const int y) {
        const int signed_z = z <= retained_z ? z : z - storage_z;
        const int padded_index_z = signed_z >= 0 ? signed_z : signed_z + padded_z;
        spectrum(x, padded_index_z, y) = component(y, z, x);
      });
  fence("dns_nonlinear_pack_dealiased_spectrum_xzy");
}

void DnsNonlinearProductTransformStage::unpack_dealiased_product_xzy(DnsState& state, int component_id) {
  auto component = state.component_view_3d(component_id);
  auto spectrum = dealiased_spectrum_xzy_;
  const int retained_z = (grid_.nz - 1) / 2;
  const int storage_z = grid_.nz;
  const int padded_z = cfg_.dealiased_physical_z;
  Kokkos::parallel_for(
      "dns_nonlinear_unpack_dealiased_product_xzy",
      Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace>({0, 0, 0}, {grid_.nx, grid_.nz, grid_.ny}),
      KOKKOS_LAMBDA(const int x, const int z, const int y) {
        const int signed_z = z <= retained_z ? z : z - storage_z;
        const int padded_index_z = signed_z >= 0 ? signed_z : signed_z + padded_z;
        component(y, z, x) = spectrum(x, padded_index_z, y);
      });
  fence("dns_nonlinear_unpack_dealiased_product_xzy");
}

void DnsNonlinearProductTransformStage::apply_distributed_dealiased(DnsState& state) {
  if (!distributed_fft_) {
    throw std::runtime_error("distributed nonlinear transform was not configured");
  }
  const int z_count = distributed_fft_->z_count();
  auto u = distributed_u_;
  auto v = distributed_v_;
  auto w = distributed_w_;
  auto product = distributed_product_;

  {
    Kokkos::Profiling::ScopedRegion region("transform_to_physical distributed_inverse_velocity_ffts");
    distributed_fft_->inverse_component_to_physical(state, cfg_.u_component, distributed_u_, "dns_nonlinear_distributed_u");
    distributed_fft_->inverse_component_to_physical(state, cfg_.v_component, distributed_v_, "dns_nonlinear_distributed_v");
    distributed_fft_->inverse_component_to_physical(state, cfg_.w_component, distributed_w_, "dns_nonlinear_distributed_w");
  }

  const double factor = cfg_.product_factor;
  const auto compute_product = [&](VelocityProduct product_id) {
    {
      Kokkos::Profiling::ScopedRegion region("transform_back distributed_build_products");
      Kokkos::parallel_for(
          "dns_nonlinear_distributed_build_product",
          Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace>(
              {0, 0, 0}, {grid_.ny, z_count, cfg_.dealiased_physical_x}),
          KOKKOS_LAMBDA(const int y, const int z, const int x) {
            double value = 0.0;
            if (product_id == VelocityProduct::UU) {
              value = u(y, z, x) * u(y, z, x);
            } else if (product_id == VelocityProduct::VV) {
              value = v(y, z, x) * v(y, z, x);
            } else if (product_id == VelocityProduct::WW) {
              value = w(y, z, x) * w(y, z, x);
            } else if (product_id == VelocityProduct::UV) {
              value = u(y, z, x) * v(y, z, x);
            } else if (product_id == VelocityProduct::VW) {
              value = v(y, z, x) * w(y, z, x);
            } else {
              value = u(y, z, x) * w(y, z, x);
            }
            product(y, z, x) = factor * value;
          });
      fence("dns_nonlinear_distributed_build_product");
    }
    {
      Kokkos::Profiling::ScopedRegion region("transform_back distributed_product_forward_ffts");
      distributed_fft_->forward_physical_to_component(state,
                                                      distributed_product_,
                                                      cfg_.product_components[static_cast<int>(product_id)],
                                                      "dns_nonlinear_distributed_product");
    }
  };

  compute_product(VelocityProduct::UU);
  compute_product(VelocityProduct::VV);
  compute_product(VelocityProduct::WW);
  compute_product(VelocityProduct::UV);
  compute_product(VelocityProduct::VW);
  compute_product(VelocityProduct::UW);
}

void DnsNonlinearProductTransformStage::apply_dealiased_2d(DnsState& state) {
  auto spectrum = dealiased_spectrum_xzy_;
  auto u = dealiased_u_xzy_;
  auto v = dealiased_v_xzy_;
  auto w = dealiased_w_xzy_;
  auto product = dealiased_product_xzy_;

  ExecutionSpace exec;
  const auto inverse_to_physical = [&](RealView3D real_field, const char* label) {
    Kokkos::Profiling::ScopedRegion region("transform_to_physical local_inverse_fft2d");
    Kokkos::Profiling::ScopedRegion f90_region("transform_to_physical RFT2D");
    KokkosFFT::execute(*dealiased_2d_inverse_plan_, spectrum, real_field, KokkosFFT::Normalization::none);
    exec.fence(std::string(label) + "_irfft2");
  };
  const auto forward_to_spectrum = [&](const char* label) {
    Kokkos::Profiling::ScopedRegion region("transform_back local_product_fft2d");
    Kokkos::Profiling::ScopedRegion f90_region("transform_back HFT2D");
    KokkosFFT::execute(*dealiased_2d_forward_plan_, product, spectrum, KokkosFFT::Normalization::none);
    exec.fence(std::string(label) + "_rfft2");
  };

  {
    Kokkos::Profiling::ScopedRegion f90_region("transform_to_physical");
    {
      Kokkos::Profiling::ScopedRegion assemble_region("transform_to_physical assemble_vvdx_2d");
      pack_dealiased_spectrum_xzy(state, cfg_.u_component);
    }
    inverse_to_physical(u, "dns_nonlinear_dealiased_u");
    {
      Kokkos::Profiling::ScopedRegion assemble_region("transform_to_physical assemble_vvdx_2d");
      pack_dealiased_spectrum_xzy(state, cfg_.v_component);
    }
    inverse_to_physical(v, "dns_nonlinear_dealiased_v");
    {
      Kokkos::Profiling::ScopedRegion assemble_region("transform_to_physical assemble_vvdx_2d");
      pack_dealiased_spectrum_xzy(state, cfg_.w_component);
    }
    inverse_to_physical(w, "dns_nonlinear_dealiased_w");
  }

  const double factor = cfg_.product_factor;
  const auto compute_product = [&](VelocityProduct product_id) {
    {
      Kokkos::Profiling::ScopedRegion region("transform_back build_products");
      Kokkos::parallel_for(
          "dns_nonlinear_dealiased_build_product_xzy",
          Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace>(
              {0, 0, 0}, {cfg_.dealiased_physical_x, cfg_.dealiased_physical_z, grid_.ny}),
          KOKKOS_LAMBDA(const int x, const int z, const int y) {
            double value = 0.0;
            if (product_id == VelocityProduct::UU) {
              value = u(x, z, y) * u(x, z, y);
            } else if (product_id == VelocityProduct::VV) {
              value = v(x, z, y) * v(x, z, y);
            } else if (product_id == VelocityProduct::WW) {
              value = w(x, z, y) * w(x, z, y);
            } else if (product_id == VelocityProduct::UV) {
              value = u(x, z, y) * v(x, z, y);
            } else if (product_id == VelocityProduct::VW) {
              value = v(x, z, y) * w(x, z, y);
            } else {
              value = u(x, z, y) * w(x, z, y);
            }
            product(x, z, y) = factor * value;
          });
      fence("dns_nonlinear_dealiased_build_product_xzy");
    }
    forward_to_spectrum("dns_nonlinear_dealiased_product");
    {
      Kokkos::Profiling::ScopedRegion repack_region("transform_back scatter_xTOz_2d");
      unpack_dealiased_product_xzy(state, cfg_.product_components[static_cast<int>(product_id)]);
    }
  };

  compute_product(VelocityProduct::UU);
  compute_product(VelocityProduct::VV);
  compute_product(VelocityProduct::WW);
  compute_product(VelocityProduct::UV);
  compute_product(VelocityProduct::VW);
  compute_product(VelocityProduct::UW);
}

void DnsNonlinearProductTransformStage::apply_dealiased(DnsState& state) {
  auto spectrum = dealiased_spectrum_;
  auto z_spectrum = dealiased_z_spectrum_;
  auto u = dealiased_u_;
  auto v = dealiased_v_;
  auto w = dealiased_w_;
  auto product = dealiased_product_;

  ExecutionSpace exec;
  const auto inverse_to_physical = [&](RealView3D real_field, const char* label) {
    Kokkos::Profiling::ScopedRegion region("transform_to_physical local_inverse_fft");
    {
      Kokkos::Profiling::ScopedRegion f90_region("transform_to_physical IFT");
      KokkosFFT::execute(*dealiased_z_inverse_plan_, spectrum, z_spectrum, KokkosFFT::Normalization::none);
      exec.fence(std::string(label) + "_z_ifft");
    }
    {
      Kokkos::Profiling::ScopedRegion f90_region("transform_to_physical RFT");
      KokkosFFT::execute(*dealiased_x_inverse_plan_, z_spectrum, real_field, KokkosFFT::Normalization::none);
      exec.fence(std::string(label) + "_x_irfft");
    }
  };
  const auto forward_to_spectrum = [&](const char* label) {
    Kokkos::Profiling::ScopedRegion region("transform_back local_product_fft");
    {
      Kokkos::Profiling::ScopedRegion f90_region("transform_back HFT");
      KokkosFFT::execute(*dealiased_x_forward_plan_, product, z_spectrum, KokkosFFT::Normalization::none);
      exec.fence(std::string(label) + "_x_rfft");
    }
    {
      Kokkos::Profiling::ScopedRegion f90_region("transform_back FFT");
      KokkosFFT::execute(*dealiased_z_forward_plan_, z_spectrum, spectrum, KokkosFFT::Normalization::none);
      exec.fence(std::string(label) + "_z_fft");
    }
  };

  {
    Kokkos::Profiling::ScopedRegion f90_region("transform_to_physical");
    {
      Kokkos::Profiling::ScopedRegion assemble_region("transform_to_physical assemble_vvdz");
      pack_dealiased_spectrum(state, cfg_.u_component);
    }
    inverse_to_physical(u, "dns_nonlinear_dealiased_u");
    {
      Kokkos::Profiling::ScopedRegion assemble_region("transform_to_physical assemble_vvdz");
      pack_dealiased_spectrum(state, cfg_.v_component);
    }
    inverse_to_physical(v, "dns_nonlinear_dealiased_v");
    {
      Kokkos::Profiling::ScopedRegion assemble_region("transform_to_physical assemble_vvdz");
      pack_dealiased_spectrum(state, cfg_.w_component);
    }
    inverse_to_physical(w, "dns_nonlinear_dealiased_w");
  }

  const double factor = cfg_.product_factor;
  const auto compute_product = [&](VelocityProduct product_id) {
    {
      Kokkos::Profiling::ScopedRegion region("transform_back build_products");
      Kokkos::parallel_for(
          "dns_nonlinear_dealiased_build_product",
          Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace>(
              {0, 0, 0}, {grid_.ny, cfg_.dealiased_physical_z, cfg_.dealiased_physical_x}),
          KOKKOS_LAMBDA(const int y, const int z, const int x) {
            double value = 0.0;
            if (product_id == VelocityProduct::UU) {
              value = u(y, z, x) * u(y, z, x);
            } else if (product_id == VelocityProduct::VV) {
              value = v(y, z, x) * v(y, z, x);
            } else if (product_id == VelocityProduct::WW) {
              value = w(y, z, x) * w(y, z, x);
            } else if (product_id == VelocityProduct::UV) {
              value = u(y, z, x) * v(y, z, x);
            } else if (product_id == VelocityProduct::VW) {
              value = v(y, z, x) * w(y, z, x);
            } else {
              value = u(y, z, x) * w(y, z, x);
            }
            product(y, z, x) = factor * value;
          });
      fence("dns_nonlinear_dealiased_build_product");
    }
    forward_to_spectrum("dns_nonlinear_dealiased_product");
    {
      Kokkos::Profiling::ScopedRegion repack_region("transform_back repack_xTOz_local");
      unpack_dealiased_product(state, cfg_.product_components[static_cast<int>(product_id)]);
    }
  };

  compute_product(VelocityProduct::UU);
  compute_product(VelocityProduct::VV);
  compute_product(VelocityProduct::WW);
  compute_product(VelocityProduct::UV);
  compute_product(VelocityProduct::VW);
  compute_product(VelocityProduct::UW);
}

void DnsNonlinearVelocityRhsStage::prepare(const DnsState& state, DnsNonlinearVelocityRhsConfig cfg) {
  const auto& grid = state.grid();
  grid.validate();
  if (grid.values_per_component() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::runtime_error("DnsNonlinearVelocityRhsStage component size exceeds int range");
  }
  grid_ = grid;
  cfg_ = std::move(cfg);

  check_component(cfg_.u_component, "u_component");
  check_component(cfg_.v_component, "v_component");
  check_component(cfg_.w_component, "w_component");
  check_component(cfg_.eta_rhs_component, "eta_rhs_component");
  check_component(cfg_.d2v_rhs_component, "d2v_rhs_component");
  for (const int component : cfg_.product_components) {
    check_component(component, "product_component");
  }
  if (cfg_.dt <= 0.0) {
    throw std::runtime_error("DnsNonlinearVelocityRhsStage requires dt > 0");
  }

  const auto line_count = grid_.line_count();
  ialfa_.resize("nonlinear_ialfa", line_count);
  ibeta_.resize("nonlinear_ibeta", line_count);
  k2_.resize("nonlinear_k2", line_count);
  derivatives_.resize("nonlinear_y_derivatives",
                      static_cast<std::size_t>(grid_.active_y_count == 0 ? grid_.ny : grid_.active_y_count) *
                          derivative_orders * derivative_stencil);
  eta_rhs_.resize("nonlinear_eta_rhs", grid_.values_per_component());
  d2v_rhs_.resize("nonlinear_d2v_rhs", grid_.values_per_component());
  old_eta_rhs_.resize("nonlinear_old_eta_rhs", grid_.values_per_component());
  old_d2v_rhs_.resize("nonlinear_old_d2v_rhs", grid_.values_per_component());
  eta_rhs_.fill(Complex(0.0, 0.0));
  d2v_rhs_.fill(Complex(0.0, 0.0));
  reset_history();
  prepared_ = true;
  wavenumbers_ready_ = false;
  derivatives_ready_ = false;
}

void DnsNonlinearVelocityRhsStage::copy_line_wavenumbers_from_host(std::span<const Complex> ialfa,
                                                                   std::span<const Complex> ibeta,
                                                                   std::span<const double> k2) {
  check_prepared("DnsNonlinearVelocityRhsStage::copy_line_wavenumbers_from_host");
  check_span_size(ialfa.size(), grid_.line_count(), "ialfa");
  check_span_size(ibeta.size(), grid_.line_count(), "ibeta");
  check_span_size(k2.size(), grid_.line_count(), "k2");
  ialfa_.copy_from_host(ialfa);
  ibeta_.copy_from_host(ibeta);
  k2_.copy_from_host(k2);
  wavenumbers_ready_ = true;
}

void DnsNonlinearVelocityRhsStage::copy_y_derivatives_from_host(std::span<const double> derivatives) {
  check_prepared("DnsNonlinearVelocityRhsStage::copy_y_derivatives_from_host");
  check_span_size(derivatives.size(), derivatives_.size(), "y derivatives");
  derivatives_.copy_from_host(derivatives);
  derivatives_ready_ = true;
}

void DnsNonlinearVelocityRhsStage::set_time_scheme(double dt,
                                                   double implicit_weight,
                                                   double explicit_weight,
                                                   double history_weight) {
  check_prepared("DnsNonlinearVelocityRhsStage::set_time_scheme");
  if (dt <= 0.0) {
    throw std::runtime_error("DnsNonlinearVelocityRhsStage requires dt > 0");
  }
  cfg_.dt = dt;
  cfg_.implicit_weight = implicit_weight;
  cfg_.explicit_weight = explicit_weight;
  cfg_.history_weight = history_weight;
}

void DnsNonlinearVelocityRhsStage::reset_history() {
  old_eta_rhs_.fill(Complex(0.0, 0.0));
  old_d2v_rhs_.fill(Complex(0.0, 0.0));
}

void DnsNonlinearVelocityRhsStage::build_velocity_products(DnsState& state) {
  check_state(state, "DnsNonlinearVelocityRhsStage::build_velocity_products");
  const auto u = state.component_view_3d(cfg_.u_component);
  const auto v = state.component_view_3d(cfg_.v_component);
  const auto w = state.component_view_3d(cfg_.w_component);
  const auto uu = state.component_view_3d(cfg_.product_components[static_cast<int>(VelocityProduct::UU)]);
  const auto vv = state.component_view_3d(cfg_.product_components[static_cast<int>(VelocityProduct::VV)]);
  const auto ww = state.component_view_3d(cfg_.product_components[static_cast<int>(VelocityProduct::WW)]);
  const auto uv = state.component_view_3d(cfg_.product_components[static_cast<int>(VelocityProduct::UV)]);
  const auto vw = state.component_view_3d(cfg_.product_components[static_cast<int>(VelocityProduct::VW)]);
  const auto uw = state.component_view_3d(cfg_.product_components[static_cast<int>(VelocityProduct::UW)]);
  const double factor = cfg_.product_factor;

  Kokkos::parallel_for(
      "dns_nonlinear_build_velocity_products",
      Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace>({0, 0, 0}, {grid_.ny, grid_.nz, grid_.nx}),
      KOKKOS_LAMBDA(const int y, const int z, const int x) {
        const double ur = u(y, z, x).real();
        const double vr = v(y, z, x).real();
        const double wr = w(y, z, x).real();
        uu(y, z, x) = Complex(factor * ur * ur, 0.0);
        vv(y, z, x) = Complex(factor * vr * vr, 0.0);
        ww(y, z, x) = Complex(factor * wr * wr, 0.0);
        uv(y, z, x) = Complex(factor * ur * vr, 0.0);
        vw(y, z, x) = Complex(factor * vr * wr, 0.0);
        uw(y, z, x) = Complex(factor * ur * wr, 0.0);
      });
  fence("dns_nonlinear_build_velocity_products");
}

void DnsNonlinearVelocityRhsStage::initialize_rhs(DnsState& state) {
  check_state(state, "DnsNonlinearVelocityRhsStage::initialize_rhs");
  if (!wavenumbers_ready_ || !derivatives_ready_) {
    throw std::runtime_error("DnsNonlinearVelocityRhsStage metadata was not uploaded");
  }

  const int ny = grid_.ny;
  const int nx = grid_.nx;
  const int lines = static_cast<int>(grid_.line_count());
  const int active_first = grid_.active_y_storage_first();
  const int active_last = grid_.active_y_storage_last_exclusive();
  const auto u = state.component_view_3d(cfg_.u_component);
  const auto v = state.component_view_3d(cfg_.v_component);
  const auto w = state.component_view_3d(cfg_.w_component);
  const auto eta = eta_rhs_.view();
  const auto d2v = d2v_rhs_.view();
  const auto old_eta = old_eta_rhs_.view();
  const auto old_d2v = old_d2v_rhs_.view();
  const auto ialfa = ialfa_.view();
  const auto ibeta = ibeta_.view();
  const auto k2 = k2_.view();
  const auto der = derivatives_.view();
  const double ni = cfg_.viscosity;
  const double inv_dt_weight = cfg_.implicit_weight / cfg_.dt;
  const double explicit_weight = cfg_.explicit_weight;
  const double history_weight = cfg_.history_weight;
  const Complex mean_pressure = cfg_.mean_pressure;

  Kokkos::parallel_for(
      "dns_nonlinear_initialize_velocity_rhs",
      Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace>({active_first, 0, 0}, {active_last, grid_.nz, grid_.nx}),
      KOKKOS_LAMBDA(const int y, const int z, const int x) {
        const int line = z * nx + x;
        const int p = y * lines + line;
        const int dy = y - active_first;
        Complex d0v(0.0, 0.0);
        Complex d2v_value(0.0, 0.0);
        Complex d4v_value(0.0, 0.0);
        Complex d0u(0.0, 0.0);
        Complex d2u(0.0, 0.0);
        Complex d0w(0.0, 0.0);
        Complex d2w(0.0, 0.0);

        for (int offset = -2; offset <= 2; ++offset) {
          const int yy = y + offset;
          if (yy < 0 || yy >= ny) continue;
          const double d0 = der(derivative_index(dy, 0, offset));
          const double d2 = der(derivative_index(dy, 2, offset));
          const double d4 = der(derivative_index(dy, 3, offset));
          d0v += d0 * v(yy, z, x);
          d2v_value += d2 * v(yy, z, x);
          d4v_value += d4 * v(yy, z, x);
          d0u += d0 * u(yy, z, x);
          d2u += d2 * u(yy, z, x);
          d0w += d0 * w(yy, z, x);
          d2w += d2 * w(yy, z, x);
        }

        const double k2_line = k2(line);
        const Complex d2v_unknown = d2v_value - k2_line * d0v;
        const Complex d2v_implicit =
            ni * (d4v_value - 2.0 * k2_line * d2v_value + k2_line * k2_line * d0v);
        d2v(p) = inv_dt_weight * d2v_unknown + d2v_implicit - history_weight * old_d2v(p);
        old_d2v(p) = Complex(0.0, 0.0);

        Complex eta_unknown(0.0, 0.0);
        Complex eta_implicit(0.0, 0.0);
        if (k2_line == 0.0) {
          eta_unknown = Complex(d0u.real(), d0w.real());
          eta_implicit = ni * Complex(d2u.real(), d2w.real());
        } else {
          eta_unknown = ibeta(line) * d0u - ialfa(line) * d0w;
          eta_implicit = ni * (ibeta(line) * (d2u - k2_line * d0u) -
                               ialfa(line) * (d2w - k2_line * d0w));
        }
        eta(p) = inv_dt_weight * eta_unknown + eta_implicit - history_weight * old_eta(p);
        if (k2_line == 0.0) eta(p) += explicit_weight * mean_pressure;
        old_eta(p) = Complex(0.0, 0.0);
      });
  fence("dns_nonlinear_initialize_velocity_rhs");
}

void DnsNonlinearVelocityRhsStage::accumulate_product(DnsState& state, VelocityProduct product) {
  check_state(state, "DnsNonlinearVelocityRhsStage::accumulate_product");
  if (!wavenumbers_ready_ || !derivatives_ready_) {
    throw std::runtime_error("DnsNonlinearVelocityRhsStage metadata was not uploaded");
  }

  const int product_index = static_cast<int>(product);
  if (product_index < 0 || product_index >= static_cast<int>(cfg_.product_components.size())) {
    throw std::runtime_error("DnsNonlinearVelocityRhsStage product index out of range");
  }

  const int ny = grid_.ny;
  const int nx = grid_.nx;
  const int lines = static_cast<int>(grid_.line_count());
  const int active_first = grid_.active_y_storage_first();
  const int active_last = grid_.active_y_storage_last_exclusive();
  const auto field = state.component_view_3d(cfg_.product_components[static_cast<std::size_t>(product_index)]);
  const auto eta = eta_rhs_.view();
  const auto d2v = d2v_rhs_.view();
  const auto old_eta = old_eta_rhs_.view();
  const auto old_d2v = old_d2v_rhs_.view();
  const auto ialfa = ialfa_.view();
  const auto ibeta = ibeta_.view();
  const auto k2 = k2_.view();
  const auto der = derivatives_.view();
  const double explicit_weight = cfg_.explicit_weight;

  Kokkos::parallel_for(
      "dns_nonlinear_accumulate_velocity_product",
      Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace>({active_first, 0, 0}, {active_last, grid_.nz, grid_.nx}),
      KOKKOS_LAMBDA(const int y, const int z, const int x) {
        const int line = z * nx + x;
        const int p = y * lines + line;
        const int dy = y - active_first;
        Complex dd0(0.0, 0.0);
        Complex dd1(0.0, 0.0);
        Complex dd2(0.0, 0.0);
        for (int offset = -2; offset <= 2; ++offset) {
          const int yy = y + offset;
          if (yy < 0 || yy >= ny) continue;
          dd0 += der(derivative_index(dy, 0, offset)) * field(yy, z, x);
          dd1 += der(derivative_index(dy, 1, offset)) * field(yy, z, x);
          dd2 += der(derivative_index(dy, 2, offset)) * field(yy, z, x);
        }

        const double k2_line = k2(line);
        Complex rhsu(0.0, 0.0);
        Complex rhsw(0.0, 0.0);
        Complex d2v_expl(0.0, 0.0);
        if (product_index == static_cast<int>(VelocityProduct::UU)) {
          rhsu = -ialfa(line) * dd0;
          d2v_expl = ialfa(line) * ialfa(line) * dd1;
        } else if (product_index == static_cast<int>(VelocityProduct::VV)) {
          d2v_expl = k2_line * dd1;
        } else if (product_index == static_cast<int>(VelocityProduct::WW)) {
          rhsw = -ibeta(line) * dd0;
          d2v_expl = ibeta(line) * ibeta(line) * dd1;
        } else if (product_index == static_cast<int>(VelocityProduct::UV)) {
          rhsu = -dd1;
          d2v_expl = ialfa(line) * dd2 + ialfa(line) * k2_line * dd0;
        } else if (product_index == static_cast<int>(VelocityProduct::VW)) {
          rhsw = -dd1;
          d2v_expl = ibeta(line) * dd2 + ibeta(line) * k2_line * dd0;
        } else {
          rhsu = -ibeta(line) * dd0;
          rhsw = -ialfa(line) * dd0;
          d2v_expl = 2.0 * ialfa(line) * ibeta(line) * dd1;
        }

        d2v(p) += explicit_weight * d2v_expl;
        old_d2v(p) += d2v_expl;

        Complex eta_expl = ibeta(line) * rhsu - ialfa(line) * rhsw;
        if (k2_line == 0.0) eta_expl = Complex(rhsu.real(), rhsw.real());
        eta(p) += explicit_weight * eta_expl;
        old_eta(p) += eta_expl;
      });
  fence("dns_nonlinear_accumulate_velocity_product");
}

void DnsNonlinearVelocityRhsStage::apply(DnsState& state) {
  check_state(state, "DnsNonlinearVelocityRhsStage::apply");
  if (cfg_.build_products_from_velocity) {
    Kokkos::Profiling::ScopedRegion region("buildrhs build_velocity_products");
    build_velocity_products(state);
  }
  debug_init_d2v_terms(grid_, cfg_, state, derivatives_, k2_, old_d2v_rhs_);
  {
    Kokkos::Profiling::ScopedRegion region("buildrhs_prepare");
    initialize_rhs(state);
  }
  debug_rhs_line("after_initialize_rhs", grid_, eta_rhs_, d2v_rhs_);
  {
    Kokkos::Profiling::ScopedRegion region("buildrhs accumulate_products");
    for (int product = 0; product < static_cast<int>(cfg_.product_components.size()); ++product) {
      Kokkos::Profiling::ScopedRegion f90_buildrhs_region("transform_back buildrhs");
      accumulate_product(state, static_cast<VelocityProduct>(product));
      const std::array<const char*, 6> labels = {
          "after_accumulate_uu",
          "after_accumulate_vv",
          "after_accumulate_ww",
          "after_accumulate_uv",
          "after_accumulate_vw",
          "after_accumulate_uw",
      };
      debug_rhs_line(labels[static_cast<std::size_t>(product)], grid_, eta_rhs_, d2v_rhs_);
    }
  }

  const int nx = grid_.nx;
  const int lines = static_cast<int>(grid_.line_count());
  const int active_first = grid_.active_y_storage_first();
  const int active_last = grid_.active_y_storage_last_exclusive();
  const auto eta_state = state.component_view_3d(cfg_.eta_rhs_component);
  const auto d2v_state = state.component_view_3d(cfg_.d2v_rhs_component);
  const auto eta = eta_rhs_.view();
  const auto d2v = d2v_rhs_.view();
  {
    Kokkos::Profiling::ScopedRegion region("buildrhs commit_velocity_rhs");
    Kokkos::parallel_for(
        "dns_nonlinear_commit_velocity_rhs",
        Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace>({active_first, 0, 0}, {active_last, grid_.nz, grid_.nx}),
        KOKKOS_LAMBDA(const int y, const int z, const int x) {
          const int line = z * nx + x;
          const int p = y * lines + line;
          eta_state(y, z, x) = eta(p);
          d2v_state(y, z, x) = d2v(p);
        });
    fence("dns_nonlinear_commit_velocity_rhs");
  }
}

void DnsNonlinearVelocityRhsStage::check_state(const DnsState& state, const char* caller) const {
  check_prepared(caller);
  const auto& grid = state.grid();
  if (grid.nx != grid_.nx || grid.ny != grid_.ny || grid.nz != grid_.nz || grid.components != grid_.components) {
    throw std::runtime_error(std::string(caller) + " grid mismatch");
  }
}

void DnsNonlinearVelocityRhsStage::check_prepared(const char* caller) const {
  if (!prepared_) {
    throw std::runtime_error(std::string(caller) + " called before prepare");
  }
}

void DnsNonlinearVelocityRhsStage::check_component(int component, const char* name) const {
  if (component < 0 || component >= grid_.components) {
    throw std::runtime_error(std::string("DnsNonlinearVelocityRhsStage ") + name + " out of range");
  }
}

} // namespace channel
